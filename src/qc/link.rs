//! `qc` step 2d: linkage of minority alleles across candidate SNP sites.
//!
//! A second pass over the alignment file collects **sparse** observations: for
//! every counted record covering at least one candidate SNP site, its state
//! (`ref`/`alt`/`other`) at each such site. Observations are consolidated per
//! read name across primary and supplementary records (one state per read/site,
//! `other` on disagreement), then pair counts are accumulated **per read**.
//!
//! Pairs of sites within `--link-window` are tested genome-wide per contig
//! (no regions, no caps): a pair is linked when it passes the concordance,
//! phasing, one-sided Fisher exact p and per-strand checks. Haplotype groups are
//! formed strongest-first with a conflict guard: a merge is rejected only for a
//! *decisive incompatibility* (`pair_conflict`). Group reads are written
//! to `P.haplotypes.tsv` and each group span becomes a `linked_haplotype` BED
//! region. Member sites get VCF INFO `LINKED`; tested minority sites in no group
//! get `UNLINKED` (or are dropped, `--unlinked-action`).
//!
//! Determinism: chunks are merged per contig in coordinate order, observations
//! are sorted before consolidation, and pair tables are integer sums, so the
//! output does not depend on the thread count.

use std::collections::HashMap;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{Context, Result};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};

use super::{Candidate, ContigCtx, QcParams, Region, Tier, UnlinkedAction, VarType};
use crate::bamio::{
    FIELDS_LINK, OP_D, OP_EQ, OP_I, OP_M, OP_N, OP_S, OP_X, keep, op_parts, open_indexed,
    terminal_clips,
};
use crate::qc::stats;
use crate::scan::{ContigCounts, FxMap, HIQ, ScanResult};

/// Minimum shared reads for a tested pair to count as a *conflict*.
pub(crate) const MIN_SHARED: u32 = 10;
/// Phasing: reference reads at one site must not carry the alt at the other.
pub(crate) const PHASE_C: f64 = 0.1;
/// Phasing: alt reads at one site must not be reference at the other.
pub(crate) const PHASE_B: f64 = 0.2;
/// Both-strands requirement: minimum alt/alt reads on each strand.
pub(crate) const MIN_STRAND_A: u32 = 2;
/// Share of group reads that must fall inside the SA window / clip side.
pub(crate) const ELSEWHERE_FRAC: f64 = 0.30;
/// SA targets on the same contig within this distance form one locus.
pub(crate) const ELSEWHERE_WIN: u32 = 1_000;

/// Read allele at one candidate site.
const REF: u8 = 1;
const ALT: u8 = 2;
const OTHER: u8 = 3;

const NT16_ASCII: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// One candidate SNP site for the linkage pass (0-based position).
#[derive(Clone, Debug)]
pub(crate) struct Site {
    pub pos: u32,
    pub refb: u8,
    pub alt: u8,
    pub af: f64,
    /// Index into the (sorted) candidate list when this site is a SNP candidate.
    pub cand: Option<usize>,
}

/// A 2×2 table `a = alt/alt, b = alt/ref, c = ref/alt, d = ref/ref`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct StrandTab {
    pub a: u32,
    pub b: u32,
    pub c: u32,
    pub d: u32,
}

/// A pair table with the all-reads, forward and reverse strand splits.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct PairTab {
    pub all: StrandTab,
    pub fwd: StrandTab,
    pub rev: StrandTab,
}

/// One tested pair with everything the grouping needs.
#[derive(Clone, Copy, Debug)]
pub(crate) struct PairRec {
    pub i: u32,
    pub j: u32,
    pub shared: u32,
    pub linked: bool,
    pub fisher: f64,
    pub conc: f64,
    /// alt/alt reads.
    pub a: u32,
    /// a + b + c: reads carrying the alt at X or at Y.
    pub abc: u32,
}

/// Is a tested pair a *decisive incompatibility*? A pair conflicts only when
/// enough reads carry the alt at X or Y (`abc >= link_conflict_min_reads`) and
/// the alt/alt share among them is low (`a/abc < link_conflict_max_share`):
/// the two minority alleles mostly sit on different reads. Underpowered but
/// compatible pairs (e.g. a=3, b=c=0) are neutral and never block a merge.
pub(crate) fn pair_conflict(a: u32, abc: u32, p: &QcParams) -> bool {
    abc >= p.link_conflict_min_reads && (a as f64 / abc as f64) < p.link_conflict_max_share
}

/// One-sided (upper tail) Fisher exact p-value for positive association.
/// `ln`-space hypergeometric, reusing the crate `lgamma`.
pub(crate) fn fisher_upper(a: u32, b: u32, c: u32, d: u32) -> f64 {
    let n = (a + b + c + d) as f64;
    let row1 = (a + b) as f64; // alt at X
    let col1 = (a + c) as f64; // alt at Y
    let row2 = (c + d) as f64;
    if row1 == 0.0 || col1 == 0.0 {
        return 1.0;
    }
    let ln_den = stats::ln_gamma(n + 1.0) - stats::ln_gamma(col1 + 1.0) - stats::ln_gamma(n - col1 + 1.0);
    let hi = (row1.min(col1)) as u32;
    let mut terms: Vec<f64> = Vec::new();
    let mut max_ln = f64::NEG_INFINITY;
    for k in a..=hi {
        let kf = k as f64;
        let ln = stats::ln_gamma(row1 + 1.0)
            - stats::ln_gamma(kf + 1.0)
            - stats::ln_gamma(row1 - kf + 1.0)
            + stats::ln_gamma(row2 + 1.0)
            - stats::ln_gamma(col1 - kf + 1.0)
            - stats::ln_gamma(row2 - col1 + kf + 1.0)
            - ln_den;
        if ln > max_ln {
            max_ln = ln;
        }
        terms.push(ln);
    }
    if terms.is_empty() {
        return 1.0;
    }
    let mut sum = 0.0f64;
    for ln in terms {
        sum += (ln - max_ln).exp();
    }
    (max_ln + sum.ln()).exp().min(1.0)
}

/// Does one 2×2 table pass concordance, phasing and Fisher at `fisher_p`?
fn strand_linked(t: &StrandTab, p: &QcParams, min_a: u32, fisher_p: f64) -> bool {
    if t.a < min_a {
        return false;
    }
    let ab = t.a + t.b;
    let ac = t.a + t.c;
    let cd = t.c + t.d;
    if ab == 0 || ac == 0 {
        return false;
    }
    if (t.a as f64 / ab as f64) < p.link_min_conc || (t.a as f64 / ac as f64) < p.link_min_conc {
        return false;
    }
    if (t.c as f64 / cd as f64) > PHASE_C || (t.b as f64 / ab as f64) > PHASE_B {
        return false;
    }
    fisher_upper(t.a, t.b, t.c, t.d) <= fisher_p
}

/// The linkage rule for one pair (see the spec): the pooled table needs
/// `a >= link_min_reads`, concordance `a/(a+b)` and `a/(a+c) >= link_min_conc`,
/// phasing `c/(c+d) <= 0.1` and `b/(a+b) <= 0.2`, and one-sided Fisher
/// p <= `link_p` (1e-6). The forward and reverse tables each need `a >= 2`, the
/// same concordance/phasing, and Fisher p <= `link_strand_p` (0.01): the strand
/// check only guards against strand-specific errors, significance comes from
/// the pooled test.
pub(crate) fn pair_linked(t: &PairTab, p: &QcParams) -> bool {
    let pooled = StrandTab { a: t.all.a, b: t.all.b, c: t.all.c, d: t.all.d };
    if !strand_linked(&pooled, p, p.link_min_reads, p.link_p) {
        return false;
    }
    strand_linked(&t.fwd, p, MIN_STRAND_A, p.link_strand_p)
        && strand_linked(&t.rev, p, MIN_STRAND_A, p.link_strand_p)
}

/// Union-find with member lists (union by size).
struct Dsu {
    parent: Vec<usize>,
    members: Vec<Vec<usize>>,
}

impl Dsu {
    fn new(n: usize) -> Self {
        Dsu { parent: (0..n).collect(), members: (0..n).map(|i| vec![i]).collect() }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        let mut c = x;
        while self.parent[c] != r {
            let nx = self.parent[c];
            self.parent[c] = r;
            c = nx;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) {
        if a == b {
            return;
        }
        let (mut big, mut small) = (a, b);
        if self.members[big].len() < self.members[small].len() {
            std::mem::swap(&mut big, &mut small);
        }
        self.parent[small] = big;
        let mb = std::mem::take(&mut self.members[small]);
        self.members[big].extend(mb);
    }
}

/// Form haplotype groups: connected components of the linked pairs, but a merge
/// is rejected only when a pair connecting the two components is a *decisive
/// incompatibility* (`pair_conflict`: `a+b+c >= link_conflict_min_reads` and
/// `a/(a+b+c) < link_conflict_max_share`). Linked pairs are processed strongest
/// first (Fisher p, then concordance, then positions) for determinism. Returns
/// components of size >= 2, each sorted, ordered by first member.
pub(crate) fn form_groups(n: usize, pairs: &[PairRec], p: &QcParams) -> Vec<Vec<usize>> {
    let mut linked: Vec<&PairRec> = pairs.iter().filter(|x| x.linked).collect();
    linked.sort_by(|a, b| {
        a.fisher
            .partial_cmp(&b.fisher)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(b.conc.partial_cmp(&a.conc).unwrap_or(std::cmp::Ordering::Equal))
            .then((a.i, a.j).cmp(&(b.i, b.j)))
    });
    // Conflict adjacency: decisive incompatibilities (see [`pair_conflict`]).
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
    for x in pairs {
        if !x.linked && pair_conflict(x.a, x.abc, p) {
            adj[x.i as usize].push(x.j);
            adj[x.j as usize].push(x.i);
        }
    }
    let mut dsu = Dsu::new(n);
    for x in linked {
        let (ri, rj) = (dsu.find(x.i as usize), dsu.find(x.j as usize));
        if ri == rj {
            continue;
        }
        let (small, big) =
            if dsu.members[ri].len() <= dsu.members[rj].len() { (ri, rj) } else { (rj, ri) };
        let members = dsu.members[small].clone();
        let mut conflict = false;
        'outer: for s in members {
            for &t in &adj[s] {
                if dsu.find(t as usize) == big {
                    conflict = true;
                    break 'outer;
                }
            }
        }
        if !conflict {
            dsu.union(ri, rj);
        }
    }
    let mut out: Vec<Vec<usize>> = Vec::new();
    let mut seen = vec![false; n];
    for i in 0..n {
        if seen[i] {
            continue;
        }
        let r = dsu.find(i);
        let mut c: Vec<usize> = Vec::new();
        for k in i..n {
            if !seen[k] && dsu.find(k) == r {
                seen[k] = true;
                c.push(k);
            }
        }
        if c.len() >= 2 {
            out.push(c);
        }
    }
    out
}

/// One consolidated read: one state per site, plus the primary record's strand,
/// SA targets and clips.
#[derive(Clone, Debug)]
pub(crate) struct ConsolidatedRead {
    pub name: String,
    pub rev: bool,
    /// `(site index, state)`, sorted by site index.
    pub sites: Vec<(u32, u8)>,
    pub sa: Vec<(u32, u32)>,
    pub lclip: u32,
    pub rclip: u32,
}

/// Indices of the reads carrying a haplotype group: `alt` at >= 2 of the
/// component's sites and `ref` at none of them.
pub(crate) fn group_read_indices(comp: &[usize], reads: &[ConsolidatedRead]) -> Vec<usize> {
    let comp_sites: Vec<u32> = comp.iter().map(|&s| s as u32).collect();
    let mut out = Vec::new();
    for (ri, r) in reads.iter().enumerate() {
        let (mut nalt, mut nref) = (0u32, 0u32);
        for &(site, st) in &r.sites {
            if comp_sites.binary_search(&site).is_ok() {
                match st {
                    ALT => nalt += 1,
                    REF => nref += 1,
                    _ => {}
                }
            }
        }
        if nalt >= 2 && nref == 0 {
            out.push(ri);
        }
    }
    out
}

/// Region-local indices of candidate minority sites (AF < `mixed_af`) with at
/// least two partners each sharing >= `MIN_SHARED` reads that are in no group.
pub(crate) fn unlinked_sites(
    afs: &[f64],
    cand: &[bool],
    partners: &[u32],
    in_group: &[bool],
    mixed_af: f64,
) -> Vec<usize> {
    (0..afs.len())
        .filter(|&i| cand[i] && afs[i] < mixed_af && partners[i] >= 2 && !in_group[i])
        .collect()
}

/// "Belongs elsewhere" evidence for a group.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Elsewhere {
    /// `(tid, 0-based window start, reads inside the window)`.
    pub sa: Option<(u32, u32, u32)>,
    pub clip_left: Option<u32>,
    pub clip_right: Option<u32>,
}

/// Shared SA locus and shared clip side over the group reads. The SA locus is
/// the 1 kb window over the sorted targets holding the most distinct reads
/// (ties -> leftmost); reported only when it holds >= 30 % of the group reads.
/// A read counts once even with duplicate SA entries. A clip side counts reads
/// whose terminal clip sums to >= `clip_min_len`.
pub(crate) fn elsewhere_of(group: &[usize], reads: &[ConsolidatedRead], clip_min_len: u32) -> Elsewhere {
    let n = group.len() as f64;
    if n == 0.0 {
        return Elsewhere::default();
    }
    let need = (ELSEWHERE_FRAC * n).ceil() as u32;
    let mut out = Elsewhere::default();
    // (tid, pos, group index) over every SA target.
    let mut t: Vec<(u32, u32, u32)> = Vec::new();
    for (gi, &ri) in group.iter().enumerate() {
        for &(tid, pos) in &reads[ri].sa {
            t.push((tid, pos, gi as u32));
        }
    }
    if !t.is_empty() {
        t.sort_unstable();
        let mut best: Option<(u32, u32, u32)> = None;
        for l in 0..t.len() {
            let (tid, start) = (t[l].0, t[l].1);
            let end = start.saturating_add(ELSEWHERE_WIN);
            let mut r = l;
            let mut reads_in: Vec<u32> = Vec::new();
            while r < t.len() && t[r].0 == tid && t[r].1 < end {
                reads_in.push(t[r].2);
                r += 1;
            }
            reads_in.sort_unstable();
            reads_in.dedup();
            let cnt = reads_in.len() as u32;
            if best.map(|b| cnt > b.2).unwrap_or(true) {
                best = Some((tid, start, cnt));
            }
        }
        if let Some((tid, start, cnt)) = best {
            if cnt >= need {
                out.sa = Some((tid, start, cnt));
            }
        }
    }
    let mut left = 0u32;
    let mut right = 0u32;
    for &ri in group {
        if reads[ri].lclip >= clip_min_len {
            left += 1;
        }
        if reads[ri].rclip >= clip_min_len {
            right += 1;
        }
    }
    if left >= need {
        out.clip_left = Some(left);
    }
    if right >= need {
        out.clip_right = Some(right);
    }
    out
}

/// One row of `P.haplotypes.tsv`: a (group, read) pair.
#[derive(Clone, Debug, serde::Serialize)]
pub struct HapRow {
    pub group: u32,
    pub contig: String,
    pub span: String,
    pub read: String,
    pub strand: String,
    /// Comma list of the group-site positions where the read carries alt.
    pub alleles: String,
    pub sa: String,
}

/// The `linkage` JSON block.
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct LinkageStats {
    /// Contigs with candidate sites processed (the pass is genome-wide per
    /// contig; there are no sub-regions).
    pub regions: u64,
    pub pairs_tested: u64,
    pub groups: u64,
    pub linked_sites: u64,
    pub unlinked_sites: u64,
    /// `snp_cluster` regions overlapping a `linked_haplotype` region.
    pub clusters_coincident: u64,
    pub timing_s: f64,
}

/// Everything the linkage pass produces.
pub(crate) struct LinkResult {
    pub regions: Vec<Region>,
    pub stats: LinkageStats,
    pub haplotypes: Vec<HapRow>,
}

/// Parse an `SA` tag into `(tid, 0-based pos)` targets; unknown contigs skipped.
fn parse_sa(s: &str, name_to_tid: &HashMap<String, u32>) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    for entry in s.split(';') {
        if entry.is_empty() {
            continue;
        }
        let mut f = entry.split(',');
        let name = f.next().unwrap_or("");
        let pos: i64 = match f.next().and_then(|x| x.parse().ok()) {
            Some(p) => p,
            None => continue,
        };
        if let Some(&tid) = name_to_tid.get(name) {
            if pos >= 1 {
                out.push((tid, (pos - 1) as u32));
            }
        }
    }
    out
}

/// Per-read metadata collected from the counted records.
#[derive(Clone, Debug, Default)]
struct ReadMeta {
    init: bool,
    has_primary: bool,
    rev: bool,
    sa: Vec<(u32, u32)>,
    lclip: u32,
    rclip: u32,
}

impl ReadMeta {
    fn merge(&mut self, rec_rev: bool, is_primary: bool, sa: Vec<(u32, u32)>, lclip: u32, rclip: u32) {
        if !self.init {
            self.init = true;
            self.has_primary = is_primary;
            self.rev = rec_rev;
            self.sa = sa;
            self.lclip = lclip;
            self.rclip = rclip;
        } else if is_primary && !self.has_primary {
            self.has_primary = true;
            self.rev = rec_rev;
            self.sa = sa;
            self.lclip = lclip;
            self.rclip = rclip;
        }
    }
}

/// One sparse observation of a chunk (read = local id within the chunk).
#[derive(Clone, Copy)]
struct Obs {
    read: u32,
    site: u32,
    state: u8,
}

/// One chunk's observations, with its own read-name interner.
struct Chunk {
    tid: usize,
    start: usize,
    names: Vec<Vec<u8>>,
    interner: HashMap<Vec<u8>, u32>,
    meta: Vec<ReadMeta>,
    obs: Vec<Obs>,
}

impl Chunk {
    fn new(tid: usize, start: usize) -> Self {
        Chunk {
            tid,
            start,
            names: Vec::new(),
            interner: HashMap::new(),
            meta: Vec::new(),
            obs: Vec::new(),
        }
    }

    fn rid(&mut self, name: &[u8]) -> u32 {
        if let Some(&id) = self.interner.get(name) {
            return id;
        }
        let id = self.names.len() as u32;
        self.names.push(name.to_vec());
        self.interner.insert(name.to_vec(), id);
        self.meta.push(ReadMeta::default());
        id
    }
}

/// Emit the sparse observations of one counted record against `sites`.
fn emit_record(
    rec: &bam::Record,
    sites: &[Site],
    name_to_tid: &HashMap<String, u32>,
    chunk: &mut Chunk,
) {
    let pos = rec.pos();
    if pos < 0 {
        return;
    }
    let cig = rec.raw_cigar();
    let enc = rec.seq().encoded;
    if rec.seq_len() == 0 || enc.len() * 2 < rec.seq_len() {
        return;
    }
    let mut span = 0usize;
    for &c in cig {
        let (op, l) = op_parts(c);
        if matches!(op, OP_M | OP_D | OP_N | OP_EQ | OP_X) {
            span += l;
        }
    }
    let start = pos as usize;
    let end = start + span;
    let lo = sites.partition_point(|s| (s.pos as usize) < start);
    let hi = sites.partition_point(|s| (s.pos as usize) < end);
    if lo == hi {
        return;
    }
    let rid = chunk.rid(rec.qname());
    let flags = rec.flags();
    let is_primary = flags & 0x800 == 0;
    let rev = flags & 0x10 != 0;
    let (lclip, rclip) = terminal_clips(cig);
    let sa = match rec.aux(b"SA") {
        Ok(Aux::String(s)) => parse_sa(s, name_to_tid),
        _ => Vec::new(),
    };
    chunk.meta[rid as usize].merge(rev, is_primary, sa, lclip, rclip);

    let mut r = start;
    let mut q = 0usize;
    let mut ptr = lo;
    for &c in cig {
        let (op, l) = op_parts(c);
        match op {
            OP_M | OP_EQ | OP_X => {
                let rend = r + l;
                while ptr < hi && (sites[ptr].pos as usize) < rend {
                    let off = sites[ptr].pos as usize - r;
                    let qq = q + off;
                    let nib = (enc[qq >> 1] >> ((!qq & 1) << 2)) & 0xf;
                    let b = NT16_ASCII[nib as usize];
                    let st = if b == sites[ptr].alt {
                        ALT
                    } else if b == sites[ptr].refb {
                        REF
                    } else {
                        OTHER
                    };
                    chunk.obs.push(Obs { read: rid, site: ptr as u32, state: st });
                    ptr += 1;
                }
                r = rend;
                q += l;
            }
            OP_D | OP_N => {
                let rend = r + l;
                while ptr < hi && (sites[ptr].pos as usize) < rend {
                    chunk.obs.push(Obs { read: rid, site: ptr as u32, state: OTHER });
                    ptr += 1;
                }
                r = rend;
            }
            OP_I | OP_S => q += l,
            _ => {}
        }
    }
}

/// Build the candidate SNP sites of one contig: alt >= 2 reads and
/// AF >= `--link-min-af` (the candidate set plus sub-threshold sites).
///
/// Only SNPs are linkage candidates (indels are never built here). Gray calls
/// are INS/DEL, so they can never be linked group members or `UNLINKED`; the
/// gray rule needs no change here.
fn build_sites(
    tid: usize,
    c: &ContigCounts,
    refseq: &[u8],
    ctx: &ContigCtx,
    cand_at: &HashMap<(usize, u32), usize>,
    link_min_af: f64,
) -> Vec<Site> {
    let mut out = Vec::new();
    for i in 0..c.len {
        if ctx.na[i] != 0 {
            continue;
        }
        let r = match refseq[i] {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => continue,
        };
        let b = &c.bases[i];
        let mut best = usize::MAX;
        let mut bx = 0u32;
        for k in 0..4 {
            if k != r && b[k] + b[4 + k] > bx {
                bx = b[k] + b[4 + k];
                best = k;
            }
        }
        if best == usize::MAX || bx < 2 {
            continue;
        }
        let dp = ctx.dp[i];
        let af = if dp == 0 { 0.0 } else { bx as f64 / dp as f64 };
        if af < link_min_af {
            continue;
        }
        out.push(Site {
            pos: i as u32,
            refb: refseq[i],
            alt: b"ACGT"[best],
            af,
            cand: cand_at.get(&(tid, i as u32)).copied(),
        });
    }
    out
}

/// Mean depth and high-MAPQ share over `[s, e)`.
fn span_depth(c: &ContigCounts, ctx: &ContigCtx, s: u64, e: u64) -> (f64, f64) {
    if e <= s {
        return (0.0, 1.0);
    }
    let (mut sd, mut sh) = (0u64, 0u64);
    for x in s..e {
        let i = x as usize;
        sd += ctx.dp[i] as u64;
        sh += c.cov[i][HIQ] as u64;
    }
    let mean = sd as f64 / (e - s) as f64;
    let hiq = if sd > 0 { sh as f64 / sd as f64 } else { 1.0 };
    (mean, hiq)
}

/// One contig's consolidated reads plus the per-read pair counts.
#[derive(Default)]
struct ContigData {
    reads: Vec<ConsolidatedRead>,
    pairmap: FxMap<(u32, u32), PairTab>,
}

/// Run the linkage second pass. Marks `LINKED`/`UNLINKED` on the candidates and
/// returns the `linked_haplotype` regions, the `linkage` block and the
/// haplotypes rows.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_linkage(
    input: &Path,
    ref_path: &Path,
    scan: &ScanResult,
    ctxs: &[ContigCtx],
    refs: &[&[u8]],
    cands: &mut [Candidate],
    existing_regions: &[Region],
    names: &[&str],
    p: &QcParams,
    low_abs: u32,
) -> Result<LinkResult> {
    let t0 = std::time::Instant::now();
    let mut stats = LinkageStats::default();

    // Candidate SNP index by (tid, pos).
    let mut cand_at: HashMap<(usize, u32), usize> = HashMap::new();
    for (i, c) in cands.iter().enumerate() {
        if c.vtype == VarType::Snp {
            cand_at.insert((c.tid, c.pos0), i);
        }
    }
    let name_to_tid: HashMap<String, u32> =
        names.iter().enumerate().map(|(i, n)| (n.to_string(), i as u32)).collect();

    // Candidate sites per contig.
    let sites: Vec<Vec<Site>> = scan
        .contigs
        .iter()
        .enumerate()
        .map(|(tid, c)| build_sites(tid, c, refs[tid], &ctxs[tid], &cand_at, p.link_min_af))
        .collect();
    let active: Vec<usize> = (0..sites.len()).filter(|&t| sites[t].len() >= 2).collect();
    stats.regions = active.len() as u64;
    if active.is_empty() {
        stats.timing_s = t0.elapsed().as_secs_f64();
        return Ok(LinkResult { regions: Vec::new(), stats, haplotypes: Vec::new() });
    }

    // Chunk tasks over the active contigs.
    let total: u64 = active.iter().map(|&t| scan.contigs[t].len as u64).sum();
    let cs = p
        .chunk_size
        .unwrap_or_else(|| (total / (4 * p.threads.max(1) as u64)).clamp(200_000, 1_000_000))
        .max(1) as usize;
    let mut tasks: Vec<(usize, usize, usize)> = Vec::new();
    for &tid in &active {
        let len = scan.contigs[tid].len;
        let mut a = 0;
        while a < len {
            let b = (a + cs).min(len);
            tasks.push((tid, a, b));
            a = b;
        }
    }

    // Second pass: one reader per thread, one chunk per task.
    let next = AtomicUsize::new(0);
    let nthreads = p.threads.max(1).min(tasks.len());
    let chunks: Vec<Result<Vec<Chunk>>> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..nthreads)
            .map(|_| {
                sc.spawn(|| -> Result<Vec<Chunk>> {
                    let mut reader = open_indexed(input, ref_path, FIELDS_LINK)?;
                    let mut rec = bam::Record::new();
                    let mut local: Vec<Chunk> = Vec::new();
                    loop {
                        let k = next.fetch_add(1, Ordering::Relaxed);
                        let Some(&(tid, start, end)) = tasks.get(k) else { break };
                        reader
                            .fetch((tid as i32, start as i64, end as i64))
                            .with_context(|| format!("fetching {}:{start}-{end}", names[tid]))?;
                        let mut chunk = Chunk::new(tid, start);
                        while let Some(res) = reader.read(&mut rec) {
                            res.with_context(|| format!("reading {}", input.display()))?;
                            if !keep(&rec) {
                                continue;
                            }
                            let rp = rec.pos();
                            if rp < start as i64 || rp >= end as i64 {
                                continue;
                            }
                            emit_record(&rec, &sites[tid], &name_to_tid, &mut chunk);
                        }
                        local.push(chunk);
                    }
                    Ok(local)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap_or_else(|_| Err(anyhow::anyhow!("linkage worker panicked"))))
            .collect()
    });

    let mut all: Vec<Chunk> = Vec::new();
    for r in chunks {
        all.extend(r?);
    }
    // Deterministic merge order: by contig, then chunk start.
    all.sort_by_key(|c| (c.tid, c.start));

    // Merge per contig and consolidate.
    let mut per_contig: Vec<ContigData> = Vec::new();
    per_contig.resize_with(sites.len(), || ContigData { reads: Vec::new(), pairmap: FxMap::default() });
    let mut i = 0;
    while i < all.len() {
        let tid = all[i].tid;
        let mut j = i;
        let mut global_index: HashMap<Vec<u8>, u32> = HashMap::new();
        let mut gnames: Vec<Vec<u8>> = Vec::new();
        let mut gmeta: Vec<ReadMeta> = Vec::new();
        let mut gobs: Vec<Obs> = Vec::new();
        while j < all.len() && all[j].tid == tid {
            let c = &all[j];
            let mut remap: Vec<u32> = Vec::with_capacity(c.names.len());
            for (lid, name) in c.names.iter().enumerate() {
                let gid = match global_index.get(name) {
                    Some(&g) => g,
                    None => {
                        let g = gnames.len() as u32;
                        gnames.push(name.clone());
                        gmeta.push(ReadMeta::default());
                        global_index.insert(name.clone(), g);
                        g
                    }
                };
                let m = &c.meta[lid];
                gmeta[gid as usize].merge(m.rev, m.has_primary, m.sa.clone(), m.lclip, m.rclip);
                remap.push(gid);
            }
            for o in &c.obs {
                gobs.push(Obs { read: remap[o.read as usize], site: o.site, state: o.state });
            }
            j += 1;
        }
        // Consolidate: one state per (read, site); disagreement -> other.
        gobs.sort_unstable_by_key(|o| (o.read, o.site));
        let mut reads: Vec<ConsolidatedRead> = Vec::new();
        let mut cur: Vec<(u32, u8)> = Vec::new();
        let mut cur_read = u32::MAX;
        let mut k = 0;
        while k < gobs.len() {
            let read = gobs[k].read;
            if read != cur_read {
                if cur_read != u32::MAX && cur.len() >= 2 {
                    let meta = &gmeta[cur_read as usize];
                    reads.push(ConsolidatedRead {
                        name: String::from_utf8_lossy(&gnames[cur_read as usize]).into_owned(),
                        rev: meta.rev,
                        sites: std::mem::take(&mut cur),
                        sa: meta.sa.clone(),
                        lclip: meta.lclip,
                        rclip: meta.rclip,
                    });
                }
                cur.clear();
                cur_read = read;
            }
            let site = gobs[k].site;
            let mut st = gobs[k].state;
            let mut m = k + 1;
            while m < gobs.len() && gobs[m].read == read && gobs[m].site == site {
                if gobs[m].state != st {
                    st = OTHER;
                }
                m += 1;
            }
            cur.push((site, st));
            k = m;
        }
        if cur_read != u32::MAX && cur.len() >= 2 {
            let meta = &gmeta[cur_read as usize];
            reads.push(ConsolidatedRead {
                name: String::from_utf8_lossy(&gnames[cur_read as usize]).into_owned(),
                rev: meta.rev,
                sites: std::mem::take(&mut cur),
                sa: meta.sa.clone(),
                lclip: meta.lclip,
                rclip: meta.rclip,
            });
        }
        // Pair counts per read, genome-wide on this contig.
        let pos: Vec<u32> = sites[tid].iter().map(|s| s.pos).collect();
        let mut pairmap: FxMap<(u32, u32), PairTab> = FxMap::default();
        for r in &reads {
            let s = &r.sites;
            for a in 0..s.len() {
                let (si, st) = s[a];
                if st != REF && st != ALT {
                    continue;
                }
                for b in (a + 1)..s.len() {
                    let (sj, stj) = s[b];
                    if pos[sj as usize] - pos[si as usize] > p.link_window {
                        break;
                    }
                    if stj != REF && stj != ALT {
                        continue;
                    }
                    let slot = match (st == ALT, stj == ALT) {
                        (true, true) => 0,
                        (true, false) => 1,
                        (false, true) => 2,
                        (false, false) => 3,
                    };
                    let tab = pairmap.entry((si, sj)).or_default();
                    {
                        let t = &mut tab.all;
                        match slot {
                            0 => t.a += 1,
                            1 => t.b += 1,
                            2 => t.c += 1,
                            _ => t.d += 1,
                        }
                    }
                    let t = if r.rev { &mut tab.rev } else { &mut tab.fwd };
                    match slot {
                        0 => t.a += 1,
                        1 => t.b += 1,
                        2 => t.c += 1,
                        _ => t.d += 1,
                    }
                }
            }
        }
        per_contig[tid] = ContigData { reads, pairmap };
        i = j;
    }

    // Groups and output, per contig.
    let mut out_regions: Vec<Region> = Vec::new();
    let mut haplotypes: Vec<HapRow> = Vec::new();
    let mut group_id: u32 = 0;
    for tid in 0..sites.len() {
        let data = std::mem::take(&mut per_contig[tid]);
        let n = sites[tid].len();
        if n < 2 {
            continue;
        }
        let afs: Vec<f64> = sites[tid].iter().map(|s| s.af).collect();
        let is_cand: Vec<bool> = sites[tid].iter().map(|s| s.cand.is_some()).collect();
        let pos: Vec<u32> = sites[tid].iter().map(|s| s.pos).collect();
        // Build the tested pairs (deterministic order).
        let mut pairs: Vec<PairRec> = Vec::with_capacity(data.pairmap.len());
        for (&(i, j), tab) in &data.pairmap {
            let shared = tab.all.a + tab.all.b + tab.all.c + tab.all.d;
            let ab = (tab.all.a + tab.all.b) as f64;
            let ac = (tab.all.a + tab.all.c) as f64;
            let conc = if ab > 0.0 && ac > 0.0 {
                (tab.all.a as f64 / ab).min(tab.all.a as f64 / ac)
            } else {
                0.0
            };
            pairs.push(PairRec {
                i,
                j,
                shared,
                linked: pair_linked(tab, p),
                fisher: fisher_upper(tab.all.a, tab.all.b, tab.all.c, tab.all.d),
                conc,
                a: tab.all.a,
                abc: tab.all.a + tab.all.b + tab.all.c,
            });
        }
        pairs.sort_by_key(|x| (x.i, x.j));
        stats.pairs_tested += pairs.len() as u64;

        let groups = form_groups(n, &pairs, p);
        let mut in_group = vec![false; n];
        for g in &groups {
            for &s in g {
                in_group[s] = true;
            }
        }
        let mut partners = vec![0u32; n];
        for x in &pairs {
            if x.shared >= MIN_SHARED {
                partners[x.i as usize] += 1;
                partners[x.j as usize] += 1;
            }
        }

        // Regions and haplotypes.
        for comp in &groups {
            group_id += 1;
            let s = pos[comp[0]] as u64;
            let e = pos[*comp.last().unwrap()] as u64 + 1;
            let group_reads = group_read_indices(comp, &data.reads);
            let fwd = group_reads.iter().filter(|&&r| !data.reads[r].rev).count() as u32;
            let rev = group_reads.len() as u32 - fwd;
            let mean_af = comp.iter().map(|&si| sites[tid][si].af).sum::<f64>() / comp.len() as f64;
            let (mean, hiq) = span_depth(&scan.contigs[tid], &ctxs[tid], s, e);
            let ratio = if ctxs[tid].mode > 0 { Some(mean / ctxs[tid].mode as f64) } else { None };
            let lowcov = mean < low_abs as f64;
            let own_red = ratio.is_some_and(|r| r >= p.depth_high) && hiq >= p.depth_high_min_mapq;
            // A linked group is the high bar (its Fisher gate); tier by the
            // disagreement fraction f = max(group mean AF, own depth fraction
            // (r-1)/r). Red only when f > 0.5.
            let depth_f = ratio.filter(|&r| r > 1.0).map_or(0.0, |r| (r - 1.0) / r);
            let f = mean_af.max(depth_f);
            let tier = if f > 0.5 { "red" } else { "orange" };
            for &si in comp {
                if let Some(ci) = sites[tid][si].cand {
                    cands[ci].linked_group = Some(group_id);
                    stats.linked_sites += 1;
                }
            }
            let elsewhere = elsewhere_of(&group_reads, &data.reads, p.clip_min_len);
            let mut extra = format!(
                "group={group_id};sites={};reads={};fwd={fwd};rev={rev};mean_af={mean_af:.3}",
                comp.len(),
                group_reads.len()
            );
            if let Some((etid, epos, cnt)) = elsewhere.sa {
                extra.push_str(&format!(";elsewhere={}:{};elsewhere_reads={cnt}", names[etid as usize], epos));
            }
            if let Some(l) = elsewhere.clip_left {
                extra.push_str(&format!(";clipped=left:{l}"));
            }
            if let Some(r) = elsewhere.clip_right {
                extra.push_str(&format!(";clipped=right:{r}"));
            }
            if lowcov {
                extra.push_str(";lowcov=1");
            }
            if own_red {
                // The group is red from its own depth ratio / high-MAPQ share:
                // give the finding the same evidence as a `depth_high` region.
                extra.push_str(";depth_high=1");
            }
            out_regions.push(Region {
                tid,
                start: s,
                end: e,
                label: "linked_haplotype",
                score: group_reads.len().min(1000) as u32,
                tier,
                reason: "linked".to_string(),
                depth_ratio: ratio,
                hiq_share: Some(hiq),
                extra,
                lowcov,
                link: None,
                clip_left: None,
                clip_right: None,
            });
            for &r in &group_reads {
                let rd = &data.reads[r];
                let mut alts: Vec<u32> = comp
                    .iter()
                    .filter(|&&si| rd.sites.binary_search_by_key(&(si as u32), |x| x.0).map(|k| rd.sites[k].1 == ALT).unwrap_or(false))
                    .map(|&si| pos[si])
                    .collect();
                alts.sort_unstable();
                let sa = rd
                    .sa
                    .first()
                    .map(|&(stid, spos)| format!("{}:{}", names[stid as usize], spos))
                    .unwrap_or_else(|| ".".to_string());
                haplotypes.push(HapRow {
                    group: group_id,
                    contig: names[tid].to_string(),
                    span: format!("{s}-{e}"),
                    read: rd.name.clone(),
                    strand: if rd.rev { "-".to_string() } else { "+".to_string() },
                    alleles: alts.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
                    sa,
                });
            }
        }

        // Tested-but-unlinked minority sites (once per contig).
        for si in unlinked_sites(&afs, &is_cand, &partners, &in_group, p.mixed_min_af) {
            if let Some(ci) = sites[tid][si].cand {
                stats.unlinked_sites += 1;
                match p.unlinked_action {
                    UnlinkedAction::Flag => cands[ci].unlinked = true,
                    UnlinkedAction::Drop => {
                        cands[ci].tier =
                            super::TierOut { tier: Tier::None, reason: "unlinked", mixed: false, strand: false };
                    }
                }
            }
        }
    }

    haplotypes.sort_by(|a, b| (a.group, &a.read).cmp(&(b.group, &b.read)));
    stats.groups = group_id as u64;
    stats.clusters_coincident = existing_regions
        .iter()
        .filter(|r| r.label == "snp_cluster")
        .filter(|cl| out_regions.iter().any(|l| l.tid == cl.tid && l.start < cl.end && cl.start < l.end))
        .count() as u64;
    stats.timing_s = t0.elapsed().as_secs_f64();
    Ok(LinkResult { regions: out_regions, stats, haplotypes })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p() -> QcParams {
        QcParams::default()
    }

    fn st(a: u32, b: u32, c: u32, d: u32) -> StrandTab {
        StrandTab { a, b, c, d }
    }

    fn read(name: &str, rev: bool, sites: Vec<(u32, u8)>, sa: Vec<(u32, u32)>, lclip: u32, rclip: u32) -> ConsolidatedRead {
        ConsolidatedRead { name: name.to_string(), rev, sites, sa, lclip, rclip }
    }

    #[test]
    fn fisher_upper_basic() {
        assert!(fisher_upper(10, 0, 0, 10) < 1e-4);
        assert!(fisher_upper(5, 5, 5, 5) > 0.4);
        assert!(fisher_upper(5, 5, 5, 5) < 0.7);
        assert_eq!(fisher_upper(0, 0, 0, 0), 1.0);
        let a = fisher_upper(8, 2, 1, 9);
        let b = fisher_upper(8, 1, 2, 9);
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    }

    #[test]
    fn pair_rule_concordance_phasing_fisher_strands() {
        let pp = p();
        let good = PairTab {
            all: st(18, 0, 0, 42),
            fwd: st(9, 0, 0, 21),
            rev: st(9, 0, 0, 21),
        };
        assert!(pair_linked(&good, &pp));
        let few = PairTab { all: st(2, 0, 0, 40), ..good };
        assert!(!pair_linked(&few, &pp));
        let conc = PairTab {
            all: st(5, 15, 0, 40),
            fwd: st(3, 8, 0, 20),
            rev: st(2, 7, 0, 20),
        };
        assert!(!pair_linked(&conc, &pp));
        let phase = PairTab {
            all: st(10, 0, 12, 28),
            fwd: st(5, 0, 6, 14),
            rev: st(5, 0, 6, 14),
        };
        assert!(!pair_linked(&phase, &pp));
        let noisy = PairTab {
            all: st(3, 0, 0, 3),
            fwd: st(2, 0, 0, 2),
            rev: st(1, 0, 0, 1),
        };
        assert!(!pair_linked(&noisy, &pp));
        let strand = PairTab {
            all: st(18, 0, 0, 42),
            fwd: st(18, 0, 0, 21),
            rev: st(0, 0, 0, 21),
        };
        assert!(!pair_linked(&strand, &pp));
    }

    /// Pooled table passes, but a strand table fails its (looser) Fisher test
    /// even though it has concordance 1.0 and a >= 2 (P2-4).
    #[test]
    fn pooled_passes_strand_fails_fisher() {
        let pp = p();
        // rev: a=2, d=5 -> Fisher 1/21 = 0.048 > link_strand_p (0.01).
        // fwd: a=10, d=20 -> Fisher 3.3e-8 <= 0.01.
        // pooled: a=12, d=25 -> Fisher far below 1e-6.
        let t = PairTab {
            all: st(12, 0, 0, 25),
            fwd: st(10, 0, 0, 20),
            rev: st(2, 0, 0, 5),
        };
        assert!(fisher_upper(12, 0, 0, 25) <= pp.link_p);
        assert!(fisher_upper(10, 0, 0, 20) <= pp.link_strand_p);
        assert!(fisher_upper(2, 0, 0, 5) > pp.link_strand_p);
        assert!(!pair_linked(&t, &pp));
        // The strand threshold is looser than the pooled one: a strand that
        // would fail `link_p` still passes at `link_strand_p`.
        let t = PairTab {
            all: st(10, 0, 0, 16),
            fwd: st(5, 0, 0, 8),
            rev: st(5, 0, 0, 8),
        };
        assert!(fisher_upper(5, 0, 0, 8) > pp.link_p);
        assert!(fisher_upper(5, 0, 0, 8) <= pp.link_strand_p);
        assert!(fisher_upper(10, 0, 0, 16) <= pp.link_p);
        assert!(pair_linked(&t, &pp));
    }

    #[test]
    fn bad_read_with_errors_everywhere_is_not_linkage() {
        // One noisy read alt at every site cannot satisfy a >= 3.
        let mut reads = Vec::new();
        reads.push(read("bad", false, vec![(0, ALT), (1, ALT)], vec![], 0, 0));
        for i in 0..20 {
            reads.push(read(&format!("r{i}"), i % 2 == 1, vec![(0, REF), (1, REF)], vec![], 0, 0));
        }
        let mut t = PairTab::default();
        for r in &reads {
            let slot = match (r.sites[0].1 == ALT, r.sites[1].1 == ALT) {
                (true, true) => 0,
                (true, false) => 1,
                (false, true) => 2,
                _ => 3,
            };
            for x in [&mut t.all, if r.rev { &mut t.rev } else { &mut t.fwd }] {
                match slot {
                    0 => x.a += 1,
                    1 => x.b += 1,
                    2 => x.c += 1,
                    _ => x.d += 1,
                }
            }
        }
        assert!(!pair_linked(&t, &p()));
    }

    fn pair(i: u32, j: u32, a: u32, d: u32) -> PairRec {
        let pt = PairTab { all: st(a, 0, 0, d), fwd: st(a, 0, 0, d), rev: st(a, 0, 0, d) };
        let lk = pair_linked(&pt, &p());
        PairRec {
            i,
            j,
            shared: a + d,
            linked: lk,
            fisher: fisher_upper(a, 0, 0, d),
            conc: 1.0,
            a,
            abc: a,
        }
    }

    #[test]
    fn groups_from_linked_pairs() {
        // 0-1-2 linked, 3-4 linked, 5 isolated.
        let pairs = vec![pair(0, 1, 18, 42), pair(1, 2, 18, 42), pair(3, 4, 18, 42)];
        let g = form_groups(6, &pairs, &p());
        assert_eq!(g, vec![vec![0, 1, 2], vec![3, 4]]);
    }

    /// Haplotype A at 1-3, B at 5-7, site 4 weakly linked to both: the
    /// tested-unlinked conflict between 3 and 5 keeps the two haplotypes apart.
    #[test]
    fn conflict_guard_splits_two_haplotypes() {
        let pairs = vec![
            // A: 1-2-3, B: 5-6-7.
            pair(1, 2, 18, 42),
            pair(2, 3, 18, 42),
            pair(5, 6, 18, 42),
            pair(6, 7, 18, 42),
            // site 4 weakly linked to 3 and to 5.
            pair(3, 4, 18, 42),
            pair(4, 5, 18, 42),
            // 3 and 5 are tested and NOT linked, and their minority alleles sit
            // on different reads: a decisive incompatibility.
            PairRec {
                i: 3,
                j: 5,
                shared: 80,
                linked: false,
                fisher: 1.0,
                conc: 0.0,
                a: 0,
                abc: 40,
            },
        ];
        let g = form_groups(8, &pairs, &p());
        // Site 4 links to both; whichever side wins it, the 3-5 conflict keeps
        // the two haplotypes separate.
        assert!(g.contains(&vec![1, 2, 3, 4]) || g.contains(&vec![1, 2, 3]) || g.contains(&vec![3, 4]), "{g:?}");
        assert!(g.contains(&vec![4, 5, 6, 7]) || g.contains(&vec![5, 6, 7]) || g.contains(&vec![4, 5]), "{g:?}");
        assert!(!g.iter().any(|c| c.contains(&1) && c.contains(&5)), "1 and 5 in one group: {g:?}");
    }

    #[test]
    fn conflict_rule_boundaries() {
        let p = p();
        // abc exactly at the minimum (6) with share < 0.3 -> conflict.
        assert!(pair_conflict(1, 6, &p));
        // abc below the minimum -> never a conflict, even with a = 0.
        assert!(!pair_conflict(0, 5, &p));
        // share exactly 0.3 -> NOT a conflict (strictly below).
        assert!(!pair_conflict(3, 10, &p));
        assert!(pair_conflict(2, 10, &p));
        // share above the threshold -> not a conflict.
        assert!(!pair_conflict(4, 10, &p));
        // a compatible underpowered pair (a=3, b=c=0) is neutral.
        assert!(!pair_conflict(3, 3, &p));
    }

    /// A neutral, underpowered pair must not block a merge.
    #[test]
    fn neutral_pair_does_not_block_merge() {
        let pairs = vec![
            pair(0, 1, 18, 42),
            pair(1, 2, 18, 42),
            // (0,2): a=3, b=c=0 -> not linked but abc=3 < 6, so neutral.
            PairRec { i: 0, j: 2, shared: 10, linked: false, fisher: 0.05, conc: 1.0, a: 3, abc: 3 },
        ];
        assert_eq!(form_groups(3, &pairs, &p()), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn group_reads_alt_two_and_ref_none() {
        let reads = vec![
            read("a", false, vec![(0, ALT), (1, ALT), (2, ALT)], vec![], 0, 0),
            read("b", true, vec![(0, ALT), (1, ALT), (2, REF)], vec![], 0, 0),
            read("c", false, vec![(0, ALT)], vec![], 0, 0),
            read("d", true, vec![(0, ALT), (1, OTHER), (2, ALT)], vec![], 0, 0),
        ];
        let g = group_read_indices(&[0, 1, 2], &reads);
        assert_eq!(g, vec![0, 3]);
    }

    #[test]
    fn unlinked_classification_once() {
        // Site 0 is linked to 1 (in a group) and shares >= 10 with 2 (not
        // linked): it must not be classified unlinked. Site 2 is in no group
        // and has >= 2 partners -> unlinked.
        let afs = vec![0.2, 0.2, 0.2];
        let cand = vec![true; 3];
        let partners = vec![2u32, 2, 2];
        let in_group = vec![true, true, false];
        let u = unlinked_sites(&afs, &cand, &partners, &in_group, 0.4);
        assert_eq!(u, vec![2]);
        // Above the minority threshold: never classified.
        let afs = vec![0.5, 0.5, 0.5];
        assert!(unlinked_sites(&afs, &cand, &partners, &[false; 3], 0.4).is_empty());
    }

    #[test]
    fn elsewhere_window_and_duplicates() {
        // One SA target per read, spaced 900 bp: no 1 kb window holds 30 % of 7
        // reads (need 3), so no locus - a single-linkage chain (span 5400) would
        // wrongly report one.
        let mut reads = Vec::new();
        for i in 0..7 {
            reads.push(read(&format!("r{i}"), false, vec![(0, ALT), (1, ALT)], vec![(0u32, 900 * i as u32)], 0, 0));
        }
        let group: Vec<usize> = (0..7).collect();
        assert_eq!(elsewhere_of(&group, &reads, 50).sa, None);

        // Two reads share 0, one is at 900, and read 0 has a duplicate entry:
        // the window [0, 1000) holds 3 distinct reads (not 4).
        let mut reads = Vec::new();
        reads.push(read("r0", false, vec![(0, ALT), (1, ALT)], vec![(0u32, 0u32), (0, 0)], 100, 0));
        reads.push(read("r1", false, vec![(0, ALT), (1, ALT)], vec![(0u32, 0u32)], 100, 0));
        reads.push(read("r2", false, vec![(0, ALT), (1, ALT)], vec![(0u32, 900u32)], 100, 0));
        for i in 3..7 {
            reads.push(read(&format!("r{i}"), false, vec![(0, ALT), (1, ALT)], vec![(0u32, 5000 + 1000 * i as u32)], 0, 0));
        }
        let group: Vec<usize> = (0..7).collect();
        let e = elsewhere_of(&group, &reads, 50);
        assert_eq!(e.sa, Some((0, 0, 3)));
        assert_eq!(e.clip_left, Some(3));
    }

    #[test]
    fn elsewhere_exactly_thirty_percent() {
        // 10 group reads, 3 with an SA target in one window -> exactly 30 %.
        let mut reads = Vec::new();
        for i in 0..10 {
            let sa = if i < 3 { vec![(0u32, 5000u32)] } else { Vec::new() };
            reads.push(read(&format!("r{i}"), false, vec![(0, ALT), (1, ALT)], sa, 0, 0));
        }
        let group: Vec<usize> = (0..10).collect();
        let e = elsewhere_of(&group, &reads, 50);
        assert_eq!(e.sa, Some((0, 5000, 3)));
    }
}
