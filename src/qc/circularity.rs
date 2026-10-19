//! `qc`: per-contig circularity evidence (JSON only).
//!
//! Sniffles2 reports a whole-contig DUP when reads span the end→start junction
//! of a contig. Here a second pass over the alignment file collects **all** of
//! each read's alignments on a contig, so a read's full alignment chain is known
//! (including supplementaries elsewhere on the same contig). Per contig it
//! derives: how many reads reach each end, how many of those are clipped past
//! the end, how many **span the junction** and the junction offset (clean /
//! overlap / gap). A read spans when two alignments of one read, on the same
//! strand, are **consecutive in the read's query chain**, one reaching the
//! contig end and the next the contig start, with the end part before the start
//! part in read orientation. The FASTA header's own `circular=` flag is parsed
//! and compared.
//!
//! Only alignments on the contig being evaluated are used (matching the
//! independent per-contig re-derivation used to review this feature); low-MAPQ
//! supplementary mappings to other repeat contigs do not break a junction chain.
//! Results are grouped by read name and sorted, so the output does not depend
//! on the thread count.

use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read};

use super::QcParams;
use super::junctions::{contig_junctions, xlink_events};
pub(crate) use super::junctions::{JunctionEvent, XlinkEvent};
use crate::bamio::{
    FIELDS_BASE, OP_D, OP_EQ, OP_I, OP_M, OP_N, OP_X, keep, op_parts, open_indexed, terminal_clips_u64,
};
use crate::fasta::Fasta;
use crate::scan::ScanResult;

/// An alignment reaches a contig end when it comes within this of the boundary
/// (the spec's clip-side threshold for `clipped_at_end`/`clipped_at_start`).
const CLIP_SIDE: u64 = 50;

/// Per-contig circularity verdict.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Verdict {
    Circular,
    Linear,
    Unclear,
}

impl Default for Verdict {
    fn default() -> Self {
        Verdict::Unclear
    }
}

/// Per-contig circularity evidence (JSON only).
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct ContigCircularity {
    pub verdict: Verdict,
    /// The FASTA header's `circular=` flag (`None` when absent/unknown).
    pub header_circular: Option<bool>,
    /// Reads with two same-strand alignments spanning the end→start junction.
    pub spanning_reads: u64,
    /// Reads with an alignment reaching within `--circ-end-tol` of the contig end.
    pub end_reads: u64,
    /// Reads with an alignment reaching within `--circ-end-tol` of position 0.
    pub start_reads: u64,
    /// `end_reads` whose end-reaching alignment has an outer clip >= 50 bp.
    pub clipped_at_end: u64,
    /// `start_reads` whose start-reaching alignment has an outer clip >= 50 bp.
    pub clipped_at_start: u64,
    /// `clean`, `overlap_<n>` or `gap_<n>` for circular contigs; `None` otherwise.
    pub junction: Option<String>,
    /// Median junction offset over the spanning reads (bp).
    pub junction_offset_median: Option<i64>,
    /// Share of spanning reads within ±`--circ-end-tol` of the median.
    pub junction_offset_consistent: Option<f64>,
}

/// Top-level circularity counts (contigs per verdict).
#[derive(Clone, Copy, Debug, Default, serde::Serialize)]
pub struct CircCounts {
    pub circular: u64,
    pub linear: u64,
    pub unclear: u64,
}

/// One aggregated circularity warning (mapped to a `Warning` by the writer).
#[derive(Clone, Debug)]
pub struct CircWarning {
    pub code: &'static str,
    pub contigs: Vec<String>,
    pub bp: u64,
    pub message: String,
}

/// Everything the circularity pass produces.
pub(crate) struct CircResult {
    /// Indexed by scan/BAM `tid`.
    pub per_contig: Vec<ContigCircularity>,
    pub counts: CircCounts,
    pub warnings: Vec<CircWarning>,
    /// Split-read junction events per contig (indexed by scan/BAM `tid`); the
    /// `qc` pass clusters and gates them into `split_junction` regions.
    pub junctions: Vec<Vec<JunctionEvent>>,
    /// Cross-contig clip events: consecutive alignments of one read on different
    /// contigs within `--xlink-max-gap` in read coordinates. The `qc` pass
    /// matches their endpoints against `clip_hotspot` regions to link the two
    /// sides of a cross-contig join.
    pub xlinks: Vec<XlinkEvent>,
    pub timing_s: f64,
}

/// One alignment of a read on the contig.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Align {
    pub(crate) pos: u64,
    pub(crate) ref_end: u64,
    pub(crate) rev: bool,
    /// Leading clip (S+H) bases.
    pub(crate) lead: u64,
    /// Trailing clip (S+H) bases.
    pub(crate) trail: u64,
    /// Aligned query bases (M/I/=/X).
    pub(crate) aln: u64,
}

impl Align {
    #[inline]
    fn reaches_end(&self, len: u64, tol: u64) -> bool {
        self.ref_end + tol >= len
    }
    #[inline]
    fn reaches_start(&self, tol: u64) -> bool {
        self.pos <= tol
    }
    /// Query interval `[start, end)` in **read orientation** (clips counted).
    #[inline]
    pub(crate) fn q_span(&self) -> (u64, u64) {
        if self.rev {
            (self.trail, self.trail + self.aln)
        } else {
            (self.lead, self.lead + self.aln)
        }
    }
}

/// Parse one `key=value` boolean flag from a FASTA header (case-insensitive;
/// surrounding brackets/punctuation such as `[topology=circular]` are
/// tolerated). `key` includes the trailing `=`.
fn parse_bool_flag(header: &str, key: &str) -> Option<bool> {
    for tok in header.split_whitespace() {
        let t = tok.trim_matches(|c: char| matches!(c, '[' | ']' | '(' | ')' | '{' | '}' | '<' | '>' | ',' | ';'));
        let low = t.to_ascii_lowercase();
        if let Some(v) = low.strip_prefix(key) {
            return match v {
                "true" | "yes" | "y" | "circular" => Some(true),
                "false" | "no" | "n" | "linear" => Some(false),
                _ => None,
            };
        }
    }
    None
}

/// Parse the FASTA header's `circular=true|false|yes|no|Y|N` flag (case-insensitive).
pub(crate) fn parse_circular_flag(header: &str) -> Option<bool> {
    parse_bool_flag(header, "circular=")
}

/// Parse the FASTA header's `topology=circular|linear` flag (autocycler /
/// dnaapler, often inside `[...]`; case-insensitive).
pub(crate) fn parse_topology_flag(header: &str) -> Option<bool> {
    parse_bool_flag(header, "topology=")
}

/// The header's circularity flag from both `circular=` and `topology=`: the one
/// present, or `None` when both are absent, unknown, or disagree. The `bool` is
/// `true` on a disagreement (the caller warns).
pub(crate) fn header_circular_flag(header: &str) -> (Option<bool>, bool) {
    match (parse_circular_flag(header), parse_topology_flag(header)) {
        (Some(a), Some(b)) if a != b => (None, true),
        (Some(a), _) => (Some(a), false),
        (None, Some(b)) => (Some(b), false),
        (None, None) => (None, false),
    }
}

/// Distil a record into an `Align`, or `None` if it has no aligned bases.
fn align_of(rec: &bam::Record) -> Option<Align> {
    let pos = rec.pos();
    if pos < 0 {
        return None;
    }
    let pos = pos as u64;
    let cig = rec.raw_cigar();
    let (lead, trail) = terminal_clips_u64(cig);
    let mut span = 0u64;
    let mut aln = 0u64;
    for &c in cig {
        let (op, l) = op_parts(c);
        let l = l as u64;
        match op {
            OP_M | OP_D | OP_N | OP_EQ | OP_X => {
                span += l;
                if matches!(op, OP_M | OP_EQ | OP_X) {
                    aln += l;
                }
            }
            OP_I => aln += l,
            _ => {}
        }
    }
    if span == 0 || aln == 0 {
        return None;
    }
    Some(Align { pos, ref_end: pos + span, rev: rec.flags() & 0x10 != 0, lead, trail, aln })
}

/// Junction offset (bp) for a spanning pair: `e` reaches the contig end, `s`
/// reaches the contig start. In read orientation the end part comes first on the
/// forward strand and the start part first on the reverse strand.
fn junction_offset(e: &Align, s: &Align, len: u64) -> i64 {
    // Reference distance traversed between the parts, corrected for how far
    // each stops short of its contig boundary.
    let d = (len as i64 - e.ref_end as i64) + s.pos as i64;
    let (_, qe_e) = e.q_span();
    let (qs_s, qe_s) = s.q_span();
    let (qs_e, _) = e.q_span();
    let gap = if e.rev { qs_e as i64 - qe_s as i64 } else { qs_s as i64 - qe_e as i64 };
    gap - d
}

/// Both query intervals monotonic and the end part before the start part in
/// read orientation: `first` starts and ends strictly before `second`. This
/// rejects reversed, identical and contained intervals (and a read that visits
/// the start before the end on the forward strand).
fn order_ok(e: &Align, s: &Align) -> bool {
    let (first, second) = if e.rev { (s, e) } else { (e, s) };
    let (fs, fe) = first.q_span();
    let (ss, se) = second.q_span();
    fs < ss && fe < se
}

/// Median of `v` (lower middle for even counts); `None` if empty.
fn median(v: &mut [i64]) -> Option<i64> {
    if v.is_empty() {
        return None;
    }
    v.sort_unstable();
    Some(v[(v.len() - 1) / 2])
}

/// Offset of the read's spanning pair, or `None`. The pair must be **consecutive
/// in the read's same-strand alignment chain** (sorted by query interval, ties
/// included), same strand, and in the right read order. The first valid pair
/// found is used; a pair is never chosen for being closer to a clean offset.
fn spanning_offset(local: &[Align], len: u64, tol: u64) -> Option<i64> {
    let mut order: Vec<usize> = (0..local.len()).collect();
    order.sort_by(|&a, &b| {
        (local[a].q_span(), local[a].pos, local[a].ref_end, local[a].rev).cmp(&(
            local[b].q_span(),
            local[b].pos,
            local[b].ref_end,
            local[b].rev,
        ))
    });
    // Rank of each alignment within its strand's sorted chain.
    let mut rank_fwd = vec![usize::MAX; local.len()];
    let mut rank_rev = vec![usize::MAX; local.len()];
    {
        let (mut f, mut r) = (0usize, 0usize);
        for &idx in &order {
            if local[idx].rev {
                rank_rev[idx] = r;
                r += 1;
            } else {
                rank_fwd[idx] = f;
                f += 1;
            }
        }
    }
    for &i in &order {
        if !local[i].reaches_end(len, tol) {
            continue;
        }
        for &j in &order {
            if i == j {
                continue;
            }
            let (e, s) = (&local[i], &local[j]);
            if !s.reaches_start(tol) || e.rev != s.rev || !order_ok(e, s) {
                continue;
            }
            // Consecutive in the sorted same-strand chain.
            let rank = if e.rev { &rank_rev } else { &rank_fwd };
            if rank[i].abs_diff(rank[j]) != 1 {
                continue;
            }
            return Some(junction_offset(e, s, len));
        }
    }
    None
}

/// Evaluate one contig from all of its (kept) alignments.
fn evaluate(
    len: u64,
    recs: &[(Vec<u8>, Align)],
    tol: u64,
    p: &QcParams,
) -> (ContigCircularity, Option<i64>, Option<i64>, bool) {
    let mut sorted: Vec<&(Vec<u8>, Align)> = recs.iter().collect();
    sorted.sort_by(|a, b| {
        (a.0.as_slice(), a.1.rev, a.1.pos, a.1.ref_end, a.1.lead, a.1.trail).cmp(&(
            b.0.as_slice(),
            b.1.rev,
            b.1.pos,
            b.1.ref_end,
            b.1.lead,
            b.1.trail,
        ))
    });

    let mut end_reads = 0u64;
    let mut start_reads = 0u64;
    let mut clipped_at_end = 0u64;
    let mut clipped_at_start = 0u64;
    let mut offsets: Vec<i64> = Vec::new();

    let mut i = 0;
    while i < sorted.len() {
        let name = &sorted[i].0;
        let mut j = i;
        while j < sorted.len() && sorted[j].0 == *name {
            j += 1;
        }
        let group: Vec<Align> = sorted[i..j].iter().map(|(_, a)| *a).collect();

        let end = group.iter().any(|r| r.reaches_end(len, tol));
        let start = group.iter().any(|r| r.reaches_start(tol));
        if end {
            end_reads += 1;
            if group.iter().any(|r| r.reaches_end(len, tol) && r.trail >= CLIP_SIDE) {
                clipped_at_end += 1;
            }
        }
        if start {
            start_reads += 1;
            if group.iter().any(|r| r.reaches_start(tol) && r.lead >= CLIP_SIDE) {
                clipped_at_start += 1;
            }
        }
        if end && start {
            if let Some(off) = spanning_offset(&group, len, tol) {
                offsets.push(off);
            }
        }
        i = j;
    }

    let spanning_reads = offsets.len() as u64;
    let min_ends = end_reads.min(start_reads);
    let frac = if min_ends == 0 { 0.0 } else { spanning_reads as f64 / min_ends as f64 };
    let verdict = if spanning_reads >= p.circ_min_reads as u64 && frac >= p.circ_min_frac {
        Verdict::Circular
    } else if end_reads >= 5
        && start_reads >= 5
        && spanning_reads == 0
        && (clipped_at_end as f64) / (end_reads as f64) < 0.1
        && (clipped_at_start as f64) / (start_reads as f64) < 0.1
    {
        Verdict::Linear
    } else {
        Verdict::Unclear
    };

    let mut junction = None;
    let mut junction_offset_median = None;
    let mut junction_offset_consistent = None;
    let mut overlap = None;
    let mut gap = None;
    if verdict == Verdict::Circular {
        let med = median(&mut offsets);
        junction_offset_median = med;
        if let Some(m) = med {
            let within = offsets.iter().filter(|&&o| (o - m).abs() <= tol as i64).count() as f64;
            junction_offset_consistent = Some((within / spanning_reads as f64 * 1000.0).round() / 1000.0);
            junction = Some(if m.abs() <= 10 {
                "clean".to_string()
            } else if m < 0 {
                overlap = Some(-m);
                format!("overlap_{}", -m)
            } else {
                gap = Some(m);
                format!("gap_{m}")
            });
        }
    }

    let unclear = verdict == Verdict::Unclear
        && len >= 10_000
        && spanning_reads == 0
        && (clipped_at_end + clipped_at_start) > 0;

    let c = ContigCircularity {
        verdict,
        header_circular: None,
        spanning_reads,
        end_reads,
        start_reads,
        clipped_at_end,
        clipped_at_start,
        junction,
        junction_offset_median,
        junction_offset_consistent,
    };
    (c, overlap, gap, unclear)
}

/// Run the circularity pass. `scan` supplies contig names/lengths, `fasta` the
/// headers.
pub(crate) fn run_circularity(
    input: &Path,
    ref_path: &Path,
    scan: &ScanResult,
    fasta: &Fasta,
    p: &QcParams,
) -> Result<CircResult> {
    let t0 = Instant::now();
    let names: Vec<String> = scan.contigs.iter().map(|c| c.name.clone()).collect();
    let lens: Vec<u64> = scan.contigs.iter().map(|c| c.len as u64).collect();
    let tol = p.circ_end_tol as u64;

    // A read's whole alignment chain on a contig is needed for adjacency, so the
    // contigs are scanned in full (chunked for parallelism; a record belongs to
    // the chunk containing its start, so nothing is collected twice).
    let total: u64 = lens.iter().sum();
    let cs = p
        .chunk_size
        .unwrap_or_else(|| (total / (4 * p.threads.max(1) as u64)).clamp(200_000, 1_000_000))
        .max(1) as usize;
    let mut tasks: Vec<(usize, usize, usize)> = Vec::new();
    for (tid, &len) in lens.iter().enumerate() {
        let mut a = 0usize;
        while (a as u64) < len {
            let b = ((a + cs) as u64).min(len) as usize;
            tasks.push((tid, a, b));
            a = b;
        }
    }

    let next = AtomicUsize::new(0);
    let nthreads = p.threads.max(1).min(tasks.len().max(1));
    let collected: Vec<Result<Vec<(usize, Vec<u8>, Align)>>> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..nthreads)
            .map(|_| {
                sc.spawn(|| -> Result<Vec<(usize, Vec<u8>, Align)>> {
                    let mut reader = open_indexed(input, ref_path, FIELDS_BASE)?;
                    let mut rec = bam::Record::new();
                    let mut local = Vec::new();
                    loop {
                        let k = next.fetch_add(1, Ordering::Relaxed);
                        let Some(&(tid, start, end)) = tasks.get(k) else { break };
                        if start >= end {
                            continue;
                        }
                        reader
                            .fetch((tid as i32, start as i64, end as i64))
                            .with_context(|| format!("fetching {}:{start}-{end}", names[tid]))?;
                        while let Some(res) = reader.read(&mut rec) {
                            res.with_context(|| format!("reading {}", input.display()))?;
                            let rp = rec.pos();
                            if !keep(&rec) || rec.tid() != tid as i32 || rp < start as i64 || rp >= end as i64 {
                                continue;
                            }
                            if let Some(a) = align_of(&rec) {
                                local.push((tid, rec.qname().to_vec(), a));
                            }
                        }
                    }
                    Ok(local)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap_or_else(|_| Err(anyhow::anyhow!("circularity worker panicked"))))
            .collect()
    });

    let mut by_contig: Vec<Vec<(Vec<u8>, Align)>> = (0..lens.len()).map(|_| Vec::new()).collect();
    for r in collected {
        for (tid, name, a) in r? {
            by_contig[tid].push((name, a));
        }
    }

    let mut per_contig: Vec<ContigCircularity> = Vec::with_capacity(lens.len());
    let mut counts = CircCounts::default();
    let mut mismatch_names: Vec<String> = Vec::new();
    let mut mismatch_msgs: Vec<String> = Vec::new();
    let mut conflict_names: Vec<String> = Vec::new();
    let mut conflict_msgs: Vec<String> = Vec::new();
    let mut overlap_names: Vec<String> = Vec::new();
    let mut overlap_msgs: Vec<String> = Vec::new();
    let mut overlap_bp = 0u64;
    let mut gap_names: Vec<String> = Vec::new();
    let mut gap_msgs: Vec<String> = Vec::new();
    let mut gap_bp = 0u64;
    let mut unclear_names: Vec<String> = Vec::new();

    for tid in 0..lens.len() {
        let name = scan.contigs[tid].name.clone();
        let (header_circular, header_conflict) = fasta
            .header_of(&name)
            .map(header_circular_flag)
            .unwrap_or((None, false));
        let (mut c, overlap, gap, unclear) = evaluate(lens[tid], &by_contig[tid], tol, p);
        c.header_circular = header_circular;
        match c.verdict {
            Verdict::Circular => counts.circular += 1,
            Verdict::Linear => counts.linear += 1,
            Verdict::Unclear => counts.unclear += 1,
        }
        let header_says_circular = header_circular == Some(true);
        if header_conflict {
            conflict_names.push(name.clone());
            conflict_msgs.push(format!("{name} (circular= and topology= disagree)"));
        }
        if (header_says_circular && c.verdict == Verdict::Linear)
            || (!header_says_circular && c.verdict == Verdict::Circular)
        {
            let h = match header_circular {
                Some(true) => "circular",
                Some(false) => "linear",
                None => "absent",
            };
            mismatch_names.push(name.clone());
            mismatch_msgs.push(format!("{name} (header={h}, verdict={})", verdict_name(c.verdict)));
        }
        if let Some(n) = overlap {
            overlap_bp += n as u64;
            overlap_names.push(name.clone());
            overlap_msgs.push(format!("{name} {n} bp"));
        }
        if let Some(n) = gap {
            gap_bp += n as u64;
            gap_names.push(name.clone());
            gap_msgs.push(format!("{name} {n} bp"));
        }
        if unclear {
            unclear_names.push(name.clone());
        }
        per_contig.push(c);
    }

    // Split-read junction events, from the same collected alignment chains. The
    // wrap decision (verdict/header/flag, see [`wrap_decision`]) is computed
    // from the verdicts above, so the circularity pass has no circular
    // dependency on the wrap flag.
    let mut junctions: Vec<Vec<JunctionEvent>> = Vec::with_capacity(lens.len());
    for tid in 0..lens.len() {
        let flag = p.linear.iter().any(|l| l == &names[tid]);
        let (wrap, _by) = wrap_decision(&per_contig[tid], flag);
        junctions.push(contig_junctions(&by_contig[tid], lens[tid], wrap, p));
    }

    // Cross-contig clip events, from the same alignment chains (the read names
    // are already in memory; nothing is decoded twice).
    let xlinks = xlink_events(&by_contig, p.xlink_max_gap as u64);

    let mut warnings = Vec::new();
    if !conflict_names.is_empty() {
        warnings.push(CircWarning {
            code: "circularity_header_conflict",
            contigs: conflict_names,
            bp: 0,
            message: format!("header circular/topology flags disagree: {}", conflict_msgs.join(", ")),
        });
    }
    if !mismatch_names.is_empty() {
        warnings.push(CircWarning {
            code: "circularity_mismatch",
            contigs: mismatch_names,
            bp: 0,
            message: format!("header/verdict mismatch: {}", mismatch_msgs.join(", ")),
        });
    }
    if !overlap_names.is_empty() {
        warnings.push(CircWarning {
            code: "junction_overlap",
            contigs: overlap_names,
            bp: overlap_bp,
            message: format!(
                "circular contig ends overlap (duplicated sequence): {}",
                overlap_msgs.join(", ")
            ),
        });
    }
    if !gap_names.is_empty() {
        warnings.push(CircWarning {
            code: "junction_gap",
            contigs: gap_names,
            bp: gap_bp,
            message: format!("circular contig junction misses sequence: {}", gap_msgs.join(", ")),
        });
    }
    if !unclear_names.is_empty() {
        warnings.push(CircWarning {
            code: "circularity_unclear",
            contigs: unclear_names.clone(),
            bp: 0,
            message: format!(
                "clipped reads at the ends but no spanning pairs (>= 10 kb): {}",
                unclear_names.join(", ")
            ),
        });
    }

    Ok(CircResult { per_contig, counts, warnings, junctions, xlinks, timing_s: t0.elapsed().as_secs_f64() })
}

fn verdict_name(v: Verdict) -> &'static str {
    match v {
        Verdict::Circular => "circular",
        Verdict::Linear => "linear",
        Verdict::Unclear => "unclear",
    }
}

/// Whether a contig wraps, and which source decided it:
///
/// 1. `--linear` (the `flag`) forces **linear**, regardless of verdict/header;
/// 2. otherwise a FASTA header `circular=true` forces **circular**;
/// 3. otherwise the verdict decides: `circular` or `unclear` wrap, `linear`
///    does not.
pub(crate) fn wrap_decision(c: &ContigCircularity, linear_flag: bool) -> (bool, &'static str) {
    if linear_flag {
        (false, "flag")
    } else if c.header_circular == Some(true) {
        (true, "header")
    } else {
        (c.verdict != Verdict::Linear, "verdict")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::qc::test_util::*;
    use rust_htslib::bam::record::{Cigar, CigarString};


    fn eval(len: u64, reads: Vec<(Vec<u8>, Align)>, tol: u64, p: &QcParams) -> ContigCircularity {
        evaluate(len, &reads, tol, p).0
    }

    /// `n` clean forward spanning reads (end part then start part).
    fn spanning(n: u64) -> Vec<(Vec<u8>, Align)> {
        let len = 30_000u64;
        let mut reads: Vec<(Vec<u8>, Align)> = Vec::new();
        for i in 0..n {
            reads.push((format!("span{i}").into_bytes(), rec(len - 100, len, false, 0, 0, 100)));
            reads.push((format!("span{i}").into_bytes(), rec(0, 200, false, 100, 0, 200)));
        }
        reads
    }

    /// Hard clips count as query: `lead`/`trail` are the sum of S and H ops.
    #[test]
    fn align_of_counts_hard_clips_as_query() {
        let mut r = bam::Record::new();
        let cigar = CigarString(vec![
            Cigar::HardClip(30),
            Cigar::SoftClip(20),
            Cigar::Match(100),
            Cigar::HardClip(10),
        ]);
        r.set(b"r", Some(&cigar), &vec![b'A'; 120], &vec![30u8; 120]);
        r.set_tid(0);
        r.set_pos(0);
        let er = align_of(&r).unwrap();
        assert_eq!((er.lead, er.trail, er.aln, er.ref_end), (50, 10, 100, 100));
        assert!(er.reaches_start(50));
    }

    #[test]
    fn circular_flag_parsing() {
        assert_eq!(parse_circular_flag("chr circular=true length=100"), Some(true));
        assert_eq!(parse_circular_flag("chr circular=TRUE"), Some(true));
        assert_eq!(parse_circular_flag("chr circular=yes"), Some(true));
        assert_eq!(parse_circular_flag("chr circular=Y"), Some(true));
        assert_eq!(parse_circular_flag("chr circular=false"), Some(false));
        assert_eq!(parse_circular_flag("chr circular=NO"), Some(false));
        assert_eq!(parse_circular_flag("chr circular=n"), Some(false));
        assert_eq!(parse_circular_flag("chr length=100"), None);
        assert_eq!(parse_circular_flag("chr circular=maybe"), None);
        assert_eq!(parse_circular_flag("circular_thing=true"), None);
    }

    /// `topology=circular|linear` (autocycler/dnaapler, often inside `[...]`)
    /// is recognised too, and the combined flag is `None` (with a conflict) when
    /// the two disagree.
    #[test]
    fn topology_flag_parsing() {
        assert_eq!(parse_topology_flag("[topology=circular]"), Some(true));
        assert_eq!(parse_topology_flag("chr topology=circular"), Some(true));
        assert_eq!(parse_topology_flag("[topology=linear]"), Some(false));
        assert_eq!(parse_topology_flag("chr TOPOLOGY=Linear"), Some(false));
        assert_eq!(parse_topology_flag("chr topology=maybe"), None);
        assert_eq!(parse_topology_flag("chr length=100"), None);

        assert_eq!(header_circular_flag("[topology=circular]"), (Some(true), false));
        assert_eq!(header_circular_flag("chr topology=linear"), (Some(false), false));
        assert_eq!(header_circular_flag("chr circular=true"), (Some(true), false));
        assert_eq!(header_circular_flag("chr circular=true topology=circular"), (Some(true), false));
        assert_eq!(header_circular_flag("chr circular=false topology=linear"), (Some(false), false));
        assert_eq!(header_circular_flag("chr circular=true topology=linear"), (None, true));
        assert_eq!(header_circular_flag("chr circular=false topology=circular"), (None, true));
        assert_eq!(header_circular_flag("chr"), (None, false));
    }

    /// The wrap decision and its source: `--linear` forces linear; else a
    /// header `circular=true` forces circular; else the verdict decides
    /// (circular/unclear wrap, linear does not).
    #[test]
    fn wrap_decision_sources() {
        let mk = |verdict, header| ContigCircularity { verdict, header_circular: header, ..Default::default() };
        // `--linear` overrides verdict and header.
        assert_eq!(wrap_decision(&mk(Verdict::Circular, Some(true)), true), (false, "flag"));
        assert_eq!(wrap_decision(&mk(Verdict::Linear, None), true), (false, "flag"));
        // A header `circular=true` overrides a linear verdict.
        assert_eq!(wrap_decision(&mk(Verdict::Linear, Some(true)), false), (true, "header"));
        assert_eq!(wrap_decision(&mk(Verdict::Unclear, Some(true)), false), (true, "header"));
        // Otherwise the verdict decides: linear does not wrap.
        assert_eq!(wrap_decision(&mk(Verdict::Linear, None), false), (false, "verdict"));
        assert_eq!(wrap_decision(&mk(Verdict::Linear, Some(false)), false), (false, "verdict"));
        // Circular and unclear wrap (even with a header circular=false).
        assert_eq!(wrap_decision(&mk(Verdict::Circular, None), false), (true, "verdict"));
        assert_eq!(wrap_decision(&mk(Verdict::Circular, Some(false)), false), (true, "verdict"));
        assert_eq!(wrap_decision(&mk(Verdict::Unclear, None), false), (true, "verdict"));
        assert_eq!(wrap_decision(&mk(Verdict::Unclear, Some(false)), false), (true, "verdict"));
    }

    /// Clean forward junction: end part [len-100, len), start part [0, 200);
    /// read orientation end part first (lead 0), start part after 100 bases.
    #[test]
    fn offset_clean_forward() {
        let len = 1_000u64;
        let e = rec(len - 100, len, false, 0, 0, 100);
        let s = rec(0, 200, false, 100, 0, 200);
        assert_eq!(junction_offset(&e, &s, len), 0);
        // Extra 30 read bases between the parts -> gap 30.
        let s2 = rec(0, 200, false, 130, 0, 200);
        assert_eq!(junction_offset(&e, &s2, len), 30);
        // Overlap of 30 read bases -> overlap 30.
        let s3 = rec(0, 200, false, 70, 0, 200);
        assert_eq!(junction_offset(&e, &s3, len), -30);
    }

    /// Clean reverse junction: start part [0, 200) first (trail 0 in read
    /// orientation), end part [len-100, len) second.
    #[test]
    fn offset_clean_reverse() {
        let len = 1_000u64;
        let s = rec(0, 200, true, 0, 0, 200);
        let e = rec(len - 100, len, true, 0, 200, 100);
        assert_eq!(junction_offset(&e, &s, len), 0);
        let e2 = rec(len - 100, len, true, 0, 230, 100);
        assert_eq!(junction_offset(&e2, &s, len), 30);
        let e3 = rec(len - 100, len, true, 0, 170, 100);
        assert_eq!(junction_offset(&e3, &s, len), -30);
    }

    /// The parts stopping short of the contig boundaries is corrected for.
    #[test]
    fn offset_corrects_short_of_boundary() {
        let len = 1_000u64;
        // End part stops 20 bp short of the end, start part starts at 10:
        // the read must have 20 + 10 = 30 bases between the aligned parts.
        let e = rec(len - 120, len - 20, false, 0, 0, 100);
        let s = rec(10, 210, false, 130, 0, 200);
        assert_eq!(junction_offset(&e, &s, len), 0);
    }

    /// A read visiting start -> end (forward) is not a spanning read.
    #[test]
    fn order_rejects_start_before_end() {
        let len = 1_000u64;
        let tol = 50u64;
        // start part query [0,200); end part query [200,400).
        let start = rec(0, 200, false, 0, 0, 200);
        let end = rec(len - 200, len, false, 200, 0, 200);
        assert!(!order_ok(&end, &start));
        let reads = vec![("r".as_bytes().to_vec(), start), ("r".as_bytes().to_vec(), end)];
        let c = eval(len, reads, tol, &QcParams::default());
        assert_eq!(c.spanning_reads, 0);
    }

    /// Reversed / identical / contained query intervals are not spanning.
    #[test]
    fn order_rejects_contained_and_identical() {
        let e = rec(900, 1_000, false, 0, 0, 400); // [0, 400)
        let s_same = rec(0, 200, false, 0, 0, 200); // [0, 200) contained / same start
        assert!(!order_ok(&e, &s_same));
        let s_short = rec(0, 200, false, 10, 0, 100); // [10, 110) contained
        assert!(!order_ok(&e, &s_short));
        // Reversed: first part starts later than the second.
        let e_rev = rec(900, 1_000, false, 200, 0, 200);
        let s_early = rec(0, 200, false, 0, 0, 200);
        assert!(!order_ok(&e_rev, &s_early));
    }

    /// A read with a third alignment between the end and start in query space
    /// does not span.
    #[test]
    fn adjacency_rejects_interior_alignment() {
        let len = 30_000u64;
        let tol = 50u64;
        // end [0,100), interior [100,200), start [200,300).
        let end = rec(len - 100, len, false, 0, 0, 100);
        let mid = rec(15_000, 15_100, false, 100, 0, 100);
        let start = rec(0, 100, false, 200, 0, 100);
        assert_eq!(spanning_offset(&[end, start], len, tol), Some(100));
        assert_eq!(spanning_offset(&[end, mid, start], len, tol), None);
    }

    /// A tie in query start does not let a non-consecutive pair through: end
    /// `[0,100)` -> interior `[0,300)` -> start `[100,200)` must not span.
    #[test]
    fn adjacency_rejects_tied_query_start() {
        let len = 30_000u64;
        let tol = 50u64;
        let end = rec(len - 100, len, false, 0, 0, 100); // (0, 100)
        let interior = rec(15_000, 15_300, false, 0, 0, 300); // (0, 300)
        let start = rec(0, 100, false, 100, 0, 100); // (100, 200)
        assert_eq!(spanning_offset(&[end, interior, start], len, tol), None);
        assert!(spanning_offset(&[end, start], len, tol).is_some());
    }

    /// Verdict boundaries.
    #[test]
    fn verdict_rules() {
        let p = QcParams::default();
        let len = 30_000u64;
        let tol = p.circ_end_tol as u64;

        // 5 clean spanning reads, plenty at both ends -> circular.
        let mut recs = spanning(5);
        for i in 0..5u64 {
            recs.push((format!("end{i}").into_bytes(), rec(len - 90, len, false, 0, 0, 100)));
            recs.push((format!("start{i}").into_bytes(), rec(0, 100, false, 0, 0, 100)));
        }
        let c = eval(len, recs, tol, &p);
        assert_eq!(c.verdict, Verdict::Circular);
        assert_eq!(c.spanning_reads, 5);
        assert_eq!(c.junction.as_deref(), Some("clean"));
        assert_eq!(c.junction_offset_median, Some(0));

        // 10 reads stop at both ends with no clips and no spanning -> linear.
        let mut recs3: Vec<(Vec<u8>, Align)> = Vec::new();
        for i in 0..10u64 {
            recs3.push((format!("end{i}").into_bytes(), rec(len - 90, len, false, 0, 0, 100)));
            recs3.push((format!("start{i}").into_bytes(), rec(0, 100, false, 0, 0, 100)));
        }
        let c3 = eval(len, recs3, tol, &p);
        assert_eq!(c3.verdict, Verdict::Linear);
        assert_eq!(c3.junction, None);

        // Many reads clipped at the ends but no spanning pairs -> unclear.
        let mut recs4: Vec<(Vec<u8>, Align)> = Vec::new();
        for i in 0..10u64 {
            recs4.push((format!("end{i}").into_bytes(), rec(len - 90, len, false, 0, 80, 100)));
            recs4.push((format!("start{i}").into_bytes(), rec(0, 100, false, 80, 0, 100)));
        }
        let c4 = eval(len, recs4, tol, &p);
        assert_eq!(c4.verdict, Verdict::Unclear);
    }

    /// `--circ-min-frac`: 3 spanning / 10 end reads is below 0.3 -> unclear;
    /// adding ends to bring the share to 0.3 makes it circular.
    #[test]
    fn verdict_min_frac_boundary() {
        let p = QcParams::default();
        let len = 30_000u64;
        let tol = p.circ_end_tol as u64;
        let mut recs = spanning(3);
        // 7 extra end-only reads -> end_reads = 10, start_reads = 10, frac = 0.3.
        for i in 0..7u64 {
            recs.push((format!("end{i}").into_bytes(), rec(len - 90, len, false, 0, 0, 100)));
            recs.push((format!("start{i}").into_bytes(), rec(0, 100, false, 0, 0, 100)));
        }
        let c = eval(len, recs, tol, &p);
        assert_eq!(c.verdict, Verdict::Circular);
        assert_eq!(c.spanning_reads, 3);
        // 8 extra end-only reads -> frac = 3/11 < 0.3 -> unclear.
        let mut recs2 = spanning(3);
        for i in 0..8u64 {
            recs2.push((format!("end{i}").into_bytes(), rec(len - 90, len, false, 0, 0, 100)));
            recs2.push((format!("start{i}").into_bytes(), rec(0, 100, false, 0, 0, 100)));
        }
        let c2 = eval(len, recs2, tol, &p);
        assert_eq!(c2.verdict, Verdict::Unclear);
    }

    /// `--circ-min-reads`: fewer than 3 spanning reads is not circular.
    #[test]
    fn verdict_min_reads_boundary() {
        let p = QcParams::default();
        let len = 30_000u64;
        let tol = p.circ_end_tol as u64;
        let c = eval(len, spanning(2), tol, &p);
        assert_eq!(c.verdict, Verdict::Unclear);
        let c2 = eval(len, spanning(3), tol, &p);
        assert_eq!(c2.verdict, Verdict::Circular);
    }

    /// `--circ-min-frac` boundary is inclusive.
    #[test]
    fn verdict_min_frac_exact_boundary() {
        let p = QcParams::default();
        let len = 30_000u64;
        let tol = p.circ_end_tol as u64;
        // 3 spanning between 10 end reads = exactly 0.3.
        let mut recs = spanning(3);
        for i in 0..7u64 {
            recs.push((format!("e{i}").into_bytes(), rec(len - 90, len, false, 0, 0, 100)));
        }
        for i in 0..7u64 {
            recs.push((format!("s{i}").into_bytes(), rec(0, 100, false, 0, 0, 100)));
        }
        let c = eval(len, recs, tol, &p);
        assert_eq!(c.verdict, Verdict::Circular);
    }

    /// The linear rule's clipped-share boundary: 1/10 = 0.1 is not `< 0.1`.
    #[test]
    fn linear_clipped_share_boundary() {
        let p = QcParams::default();
        let len = 30_000u64;
        let tol = p.circ_end_tol as u64;
        let mut recs: Vec<(Vec<u8>, Align)> = Vec::new();
        for i in 0..10u64 {
            recs.push((format!("end{i}").into_bytes(), rec(len - 90, len, false, 0, 0, 100)));
            recs.push((format!("start{i}").into_bytes(), rec(0, 100, false, 0, 0, 100)));
        }
        // One *of the ten* end reads clipped at the end -> 0.1, not < 0.1 -> not linear.
        recs.push(("end0".as_bytes().to_vec(), rec(len - 90, len, false, 0, 80, 100)));
        let c = eval(len, recs, tol, &p);
        assert_eq!(c.verdict, Verdict::Unclear);
    }

    /// `--circ-end-tol` boundary: an alignment ending exactly `tol` short reaches
    /// the end, one base further does not.
    #[test]
    fn end_tol_boundary() {
        let len = 1_000u64;
        let tol = 50u64;
        assert!(rec(len - 50 - 100, len - 50, false, 0, 0, 100).reaches_end(len, tol));
        assert!(!rec(len - 51 - 100, len - 51, false, 0, 0, 100).reaches_end(len, tol));
        assert!(rec(50, 150, false, 0, 0, 100).reaches_start(tol));
        assert!(!rec(51, 151, false, 0, 0, 100).reaches_start(tol));
    }

    /// A single record reaching both ends (a read longer than the contig) is
    /// not two alignments and must not count as spanning.
    #[test]
    fn single_record_spanning_whole_contig() {
        let p = QcParams::default();
        let len = 2_000u64;
        let tol = p.circ_end_tol as u64;
        let reads = vec![("r".as_bytes().to_vec(), rec(0, len, false, 0, 0, len))];
        let c = eval(len, reads, tol, &p);
        assert_eq!(c.spanning_reads, 0);
        assert_eq!(c.verdict, Verdict::Unclear);
    }








}
