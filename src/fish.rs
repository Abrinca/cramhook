//! `fish`: extract the reads that disagree with the assembly.
//!
//! Two passes over the alignment, fish's own per-read checks and thresholds, one
//! FASTQ plus a reasons table. The first pass selects reads and builds the
//! reasons table; the second decodes the FASTQ payload of the selected reads
//! only. The optional `--bait` category adds a `lift` session on assembly + bait.
//!
//! Counted records skip secondary (0x100), duplicate (0x400) and QC-fail
//! (0x200); primary and supplementary are both checked. Every counted record
//! must be free of hard clips (`fish` requires a `minimap2 -Y` alignment, so
//! supplementary records carry SEQ and STATE is soft-clipped): any hard clip is
//! an error. The per-read predicates reuse the shared `bamio`/`scan` helpers
//! (`terminal_clips`, `raw_indels`, `merge_large`), so "clipped" and "large
//! indel" mean the same as in `qc`.
//!
//! A read is selected when any of its counted records qualifies for any enabled
//! category: `unmapped`, `clip`, `large_ins` / `large_del`, `linked`
//! (`--haplotypes`), `locus` (`--findings`) or `bait` (`--bait`). In the
//! `findings` scope, `clip` and `large_*` only count when the event position
//! lies in a finding of tier >= `--min-tier` ± `--finding-pad`; the `genome`
//! scope counts them anywhere. `unmapped` is genome-wide in both scopes.
//!
//! `bait` (when `--bait` is given) is answered by the `lift` session: the input
//! is lifted onto `P.bait_ref.fa` = assembly + bait contigs and every read with a
//! counted record on a bait contig in the resulting `P.bait.cram` is selected.
//! `--bait-records primary` narrows this to reads whose primary record is on a
//! bait contig (the default `any` also keeps reads that reach the bait only
//! through a supplementary alignment). A read whose original primary record is
//! duplicate/QC-fail is skipped: `lift` remaps such a primary, so its remapped
//! record would otherwise look counted. All other categories are still evaluated
//! on the original `--in` alignment.
//!
//! Parallelism matches `scan`: contigs are split into chunks, a record belongs
//! to the chunk holding its start, and one task handles coordinate-less records.
//! Workers return per-read aggregates that are merged deterministically, so the
//! outputs are byte-identical across thread counts (the JSON's timing is the
//! only difference). The FASTQ payload of each selected read is the smallest
//! `(rank, tid, pos)` donor (counted primary first), captured by a second pass
//! that decodes only the selected reads.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use anyhow::{Context, Result, anyhow, bail};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};

use crate::bamio::{
    FIELDS_FISH_FULL, FIELDS_FISH_MIN, OP_D, OP_EQ, OP_H, OP_M, OP_N, OP_X, aux_int, op_parts,
    open_indexed, terminal_clips,
};
use crate::fasta::Fasta;
use crate::lift::{LiftOptions, lift};
use crate::scan::{RawIndel, merge_large, raw_indels};
use crate::translate::{TranslateConfig, revcomp};
use crate::version::versioned;

/// Reason name for unmapped reads.
const R_UNMAPPED: &str = "unmapped";
/// Reason name for a terminal clip.
const R_CLIP: &str = "clip";
/// Reason name for a merged large insertion.
const R_INS: &str = "large_ins";
/// Reason name for a merged large deletion.
const R_DEL: &str = "large_del";
/// Reason name for a read of a selected haplotype group.
const R_LINKED: &str = "linked";
/// Reason name for any record overlapping a finding.
const R_LOCUS: &str = "locus";
/// Reason name for a read placed on a bait contig by the lift session.
const R_BAIT: &str = "bait";

/// Scope of `clip` / `large_*` selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FishScope {
    /// Only within a finding of tier >= `--min-tier` ± `--finding-pad`.
    Findings,
    /// Anywhere.
    Genome,
}

/// Minimum finding tier for the `findings` scope, `--locus` and linked selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum FishTier {
    Yellow,
    Orange,
    Red,
}

impl FishTier {
    fn name(self) -> &'static str {
        match self {
            FishTier::Yellow => "yellow",
            FishTier::Orange => "orange",
            FishTier::Red => "red",
        }
    }
}

/// Which records of a read select it for the `bait` category.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BaitRecords {
    /// Any counted record on a bait contig: primary or supplementary (default,
    /// the historical behaviour).
    Any,
    /// Only a read whose **primary** record is on a bait contig; a read that
    /// reaches the bait solely through a supplementary alignment is not selected
    /// (its primary record lies elsewhere).
    Primary,
}

impl BaitRecords {
    fn name(self) -> &'static str {
        match self {
            BaitRecords::Any => "any",
            BaitRecords::Primary => "primary",
        }
    }
}

/// `fish` parameters.
#[derive(Clone, Debug)]
pub struct FishParams {
    /// Worker threads.
    pub threads: usize,
    /// Chunk size override (bp) for the contig pass; `None` picks one from the
    /// genome size and threads. Test-only.
    pub chunk_size: Option<u64>,
    /// Linear contigs; all others are circular (a clip at a contig end is an
    /// origin crossing and does not count).
    pub linear: Vec<String>,
    /// `P.findings.bed`; required by `--locus` and `--scope findings` and the
    /// default scope becomes `findings` when it is given.
    pub findings: Option<PathBuf>,
    /// `P.haplotypes.tsv`; linked selection requires it (there is no positive
    /// `--linked` flag).
    pub haplotypes: Option<PathBuf>,
    /// Bait FASTA(s) (e.g. a phage or plasmid contig). Repeatable; a
    /// multi-record FASTA is fine. Enables the `bait` category.
    pub bait: Vec<PathBuf>,
    /// Tiny-input fallback threshold in MB (10^6 bytes) for the bait `lift`
    /// session; `None` uses the per-preset default (like `cramhook lift`), `0`
    /// disables the fallback so reads are translated as well as remapped.
    pub bait_min_input_mb: Option<f64>,
    /// Which records select a read for the `bait` category: any counted record
    /// on a bait contig (primary or supplementary, default) or only a primary
    /// record.
    pub bait_records: BaitRecords,
    /// Explicit scope; `None` picks `findings` when a findings file is given.
    pub scope: Option<FishScope>,
    /// Minimum finding tier for the `findings` scope, `--locus` and linked selection.
    pub min_tier: FishTier,
    /// bp around each finding for the `findings` scope and `finding` column.
    pub finding_pad: u64,
    pub no_unmapped: bool,
    pub no_clips: bool,
    pub no_indels: bool,
    pub no_linked: bool,
    /// Select any record overlapping a finding ± `locus_flank`.
    pub locus: bool,
    pub locus_flank: u64,
    /// Terminal clip length threshold.
    pub clip_min_len: u32,
    /// Minimum merged large-indel size.
    pub indel_min_len: u32,
    /// Merge same-type indels separated by at most this many ref bp.
    pub indel_merge_gap: u32,
    /// Only indels of at least this size take part in the merge.
    pub indel_piece_min: u32,
}

impl Default for FishParams {
    fn default() -> Self {
        FishParams {
            threads: 1,
            chunk_size: None,
            linear: Vec::new(),
            findings: None,
            haplotypes: None,
            bait: Vec::new(),
            bait_min_input_mb: None,
            bait_records: BaitRecords::Any,
            scope: None,
            min_tier: FishTier::Yellow,
            finding_pad: 1000,
            no_unmapped: false,
            no_clips: false,
            no_indels: false,
            no_linked: false,
            locus: false,
            locus_flank: 5000,
            clip_min_len: 50,
            indel_min_len: 50,
            indel_merge_gap: 50,
            indel_piece_min: 10,
        }
    }
}

/// Reads and bp selected by one reason.
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct ReasonCounts {
    pub reads: u64,
    pub bp: u64,
}

/// Selected reads overall and per reason.
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct FishReads {
    /// Distinct selected reads.
    pub reads: u64,
    /// Summed length of the selected reads.
    pub bp: u64,
    /// Per reason (`unmapped`, `clip`, `large_ins`, `large_del`, `linked`,
    /// `locus`); a read can appear under more than one.
    pub by_reason: BTreeMap<String, ReasonCounts>,
}

/// Input record counts (the whole file).
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct InputCounts {
    pub total: u64,
    pub primary: u64,
    pub supplementary: u64,
    pub secondary: u64,
    pub duplicate: u64,
    pub qcfail: u64,
    pub unmapped: u64,
    pub mapped_bp: u64,
    pub unmapped_bp: u64,
}

impl InputCounts {
    fn add(&mut self, o: &InputCounts) {
        self.total += o.total;
        self.primary += o.primary;
        self.supplementary += o.supplementary;
        self.secondary += o.secondary;
        self.duplicate += o.duplicate;
        self.qcfail += o.qcfail;
        self.unmapped += o.unmapped;
        self.mapped_bp += o.mapped_bp;
        self.unmapped_bp += o.unmapped_bp;
    }
}

/// Categories and thresholds actually used.
#[derive(Clone, Debug, serde::Serialize)]
pub struct FishParamsOut {
    pub scope: String,
    pub min_tier: String,
    pub finding_pad: u64,
    pub unmapped: bool,
    pub clips: bool,
    pub indels: bool,
    pub linked: bool,
    /// `bait` is on when `--bait` was given.
    pub bait: bool,
    pub locus: bool,
    pub locus_flank: u64,
    pub clip_min_len: u32,
    pub indel_min_len: u32,
    pub indel_merge_gap: u32,
    pub indel_piece_min: u32,
}

/// Wall times.
#[derive(Clone, Debug, serde::Serialize)]
pub struct FishTiming {
    pub pass: f64,
    pub total: f64,
}

/// The `P.fish.json` document.
#[derive(Clone, Debug, serde::Serialize)]
pub struct FishSummary {
    pub reads: FishReads,
    pub records: InputCounts,
    pub params: FishParamsOut,
    /// The `bait` step, when `--bait` was given.
    pub bait: Option<BaitSummary>,
    /// Notes about categories turned off or skipped.
    pub notes: Vec<String>,
    pub timing_s: FishTiming,
}

/// One bait contig and the reads selected on it.
#[derive(Clone, Debug, serde::Serialize)]
pub struct BaitContig {
    pub name: String,
    pub length: u64,
    /// Reads whose best bait record is on this contig.
    pub reads: u64,
    /// Summed read length of those reads.
    pub bp: u64,
    /// Of `reads`, those whose best bait record has MAPQ 0.
    pub mapq0_reads: u64,
}

/// The `bait` block of `P.fish.json`.
#[derive(Clone, Debug, serde::Serialize)]
pub struct BaitSummary {
    /// Bait contigs in reference order (assembly contigs are not listed).
    pub contigs: Vec<BaitContig>,
    /// Distinct reads with a counted record on a bait contig.
    pub reads: u64,
    /// Summed read length of those reads.
    pub bp: u64,
    /// Reads whose best bait record has MAPQ 0: the read fits the assembly
    /// (the integrated copy) equally well, so both copies are the same sequence.
    pub mapq0_reads: u64,
    /// The `--bait-records` mode actually used: `any` or `primary`.
    pub mode: String,
    /// Reads whose only counted bait record is supplementary: selected under
    /// `any`, excluded under `primary` (a diagnostic, reported in both modes).
    pub supplementary_only_reads: u64,
    /// `lift` reads translated onto the bait reference.
    pub translated_reads: u64,
    /// `lift` reads remapped against the bait reference.
    pub remapped_reads: u64,
    /// `lift` input-contract check passed.
    pub contract_ok: bool,
    /// `lift` contract-violation reason, when the contract failed.
    pub contract_reason: Option<String>,
    /// `lift` fallback reason (e.g. `"tiny_input"`), when lift remapped all
    /// reads without translating any.
    pub fallback: Option<String>,
    /// Wall time of the whole bait step (reference, lift, reading `P.bait.cram`).
    pub timing_s: f64,
}

// ---------------------------------------------------------------------------
// per-read aggregation

/// The best event of one read for one reason.
#[derive(Clone, Debug)]
struct Event {
    /// Event length (clip length / merged indel size / read length).
    len: u64,
    /// 0-based event position.
    pos: u64,
    /// Event came from the primary record.
    primary: bool,
    strand: char,
    detail: String,
    finding: Option<u32>,
    contig: String,
    sa: String,
}

/// Is `a` the better event than `b`? Largest length, then smallest position,
/// then primary over supplementary; the rest only break exact ties so the choice
/// never depends on the thread count.
fn event_better(a: &Event, b: &Event) -> bool {
    if a.len != b.len {
        return a.len > b.len;
    }
    if a.pos != b.pos {
        return a.pos < b.pos;
    }
    if a.primary != b.primary {
        return a.primary;
    }
    if a.finding != b.finding {
        return match (a.finding, b.finding) {
            (Some(x), Some(y)) => x < y,
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => false,
        };
    }
    if a.detail != b.detail {
        return a.detail < b.detail;
    }
    if a.contig != b.contig {
        return a.contig < b.contig;
    }
    a.strand < b.strand
}

/// The FASTQ payload of one record, with the rank that decided it.
///
/// Rank 0 is a counted primary record (used even when it does not qualify);
/// rank 1 is a qualifying record. The best payload is the smallest
/// `(rank, tid, pos)`, so the choice never depends on the scheduling order of
/// the parallel pass. Every counted record of a read carries the same read in
/// its own orientation, so applying the record's own `0x10` recovers the same
/// original read whichever donor is picked.
#[derive(Clone)]
struct Payload {
    rank: u8,
    tid: i32,
    pos: i64,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// One read's aggregate: the best event per `(reason, finding)`, plus the FASTQ
/// payload (filled by the second pass; see `run_payload_pass`).
#[derive(Default)]
struct ReadAgg {
    events: HashMap<(&'static str, Option<u32>), Event>,
    payload: Option<Payload>,
    /// Read length (SEQ length), used for the per-reason bp.
    read_len: u64,
}

impl ReadAgg {
    fn offer(&mut self, reason: &'static str, ev: Event) {
        let key = (reason, ev.finding);
        match self.events.get(&key) {
            Some(old) if !event_better(&ev, old) => {}
            _ => {
                self.events.insert(key, ev);
            }
        }
    }
}

/// The FASTQ SEQ/QUAL (original orientation, Phred+33) of one record, with the
/// same missing-QUAL rule as `lift` (`write_fastq_entry`'s caller): an absent or
/// all-0xff QUAL becomes `!`.
fn payload_bytes(rec: &bam::Record) -> (Vec<u8>, Vec<u8>) {
    let mut seq = rec.seq().as_bytes();
    let mut qual = rec.qual().to_vec();
    if rec.flags() & 0x10 != 0 {
        seq = revcomp(&seq);
        qual.reverse();
    }
    let missing = !qual.is_empty() && qual.iter().all(|&b| b == 0xff);
    let qual_ascii: Vec<u8> = if qual.is_empty() || missing {
        vec![b'!'; seq.len()]
    } else {
        qual.iter().map(|&b| b.saturating_add(33)).collect()
    };
    (seq, qual_ascii)
}

/// The payload rank of a counted record: 0 for a primary record (used even when
/// it does not qualify), 1 otherwise. The best payload of a read is the smallest
/// `(rank, tid, pos)`, so the choice never depends on job scheduling; the bytes
/// are the same whichever donor is picked (each record stores the read in its
/// own orientation).
fn payload_rank(rec: &bam::Record) -> u8 {
    if rec.flags() & 0x800 == 0 {
        0
    } else {
        1
    }
}

/// Record a payload candidate for a selected read if it ranks ahead of the one
/// already held.
fn note_payload(out: &mut HashMap<String, Payload>, selected: &HashSet<String>, rec: &bam::Record) {
    let flags = rec.flags();
    if flags & 0x100 != 0 || flags & 0x200 != 0 || flags & 0x400 != 0 {
        return;
    }
    let name = String::from_utf8_lossy(rec.qname());
    if !selected.contains(name.as_ref()) {
        return;
    }
    let rank = payload_rank(rec);
    let key = (rank, rec.tid(), rec.pos());
    if let Some(e) = out.get(name.as_ref()) {
        if (e.rank, e.tid, e.pos) <= key {
            return;
        }
    }
    let (seq, qual) = payload_bytes(rec);
    out.insert(name.into_owned(), Payload { rank, tid: rec.tid(), pos: rec.pos(), seq, qual });
}

/// Get (or create) a read's aggregate and record its read length.
fn agg_of<'a>(aggs: &'a mut HashMap<String, ReadAgg>, name: &str, rec: &bam::Record) -> &'a mut ReadAgg {
    let a = aggs.entry(name.to_string()).or_default();
    let n = rec.seq_len() as u64;
    if n > a.read_len {
        a.read_len = n;
    }
    a
}

// ---------------------------------------------------------------------------
// findings / haplotypes

/// One finding interval from `P.findings.bed`.
struct FindingIv {
    /// The id from column 4 (`finding_<n>`), kept as written.
    n: u32,
    tid: usize,
    start: u64,
    end: u64,
    tier: FishTier,
}

/// `finding_<n>` attribution (any tier) and scope membership.
struct Findings {
    all: Vec<FindingIv>,
}

impl Findings {
    fn within(&self, f: &FindingIv, pos: u64, tid: usize, pad: u64) -> bool {
        f.tid == tid && pos >= f.start.saturating_sub(pad) && pos < f.end + pad
    }

    /// Every finding (tier >= `min` when given) containing `pos` ± `pad`,
    /// ascending by id.
    fn all_at(&self, tid: usize, pos: u64, pad: u64, min: Option<FishTier>) -> Vec<u32> {
        let mut v: Vec<u32> = self
            .all
            .iter()
            .filter(|f| min.map_or(true, |m| f.tier >= m) && self.within(f, pos, tid, pad))
            .map(|f| f.n)
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    }

    /// Every finding (tier >= `min` when given) overlapping `[start, end)` ±
    /// `pad`, ascending by id.
    fn all_overlap(
        &self,
        tid: usize,
        start: u64,
        end: u64,
        pad: u64,
        min: Option<FishTier>,
    ) -> Vec<u32> {
        let mut v: Vec<u32> = self
            .all
            .iter()
            .filter(|f| min.map_or(true, |m| f.tier >= m) && f.tid == tid)
            .filter(|f| start < f.end + pad && f.start.saturating_sub(pad) < end)
            .map(|f| f.n)
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    }
}

fn parse_tier(s: &str) -> Option<FishTier> {
    match s {
        "yellow" => Some(FishTier::Yellow),
        "orange" => Some(FishTier::Orange),
        "red" => Some(FishTier::Red),
        _ => None,
    }
}

fn read_findings(path: &Path, name_to_tid: &HashMap<String, usize>) -> Result<Findings> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading findings {}", path.display()))?;
    let mut all = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            continue;
        }
        // Preserve the id written in column 4 (`finding_<n>`); never renumber.
        let Some(n) = f[3].strip_prefix("finding_").and_then(|s| s.parse::<u32>().ok()) else {
            continue;
        };
        let Some(&tid) = name_to_tid.get(f[0]) else { continue };
        let (Ok(start), Ok(end)) = (f[1].parse::<u64>(), f[2].parse::<u64>()) else { continue };
        let Some(tier) = parse_tier(f[5]) else { continue };
        all.push(FindingIv { n, tid, start, end, tier });
    }
    Ok(Findings { all })
}

/// Selected haplotype groups: read name -> group id, and group -> span.
struct Groups {
    /// read -> group id (only selected groups).
    read_group: HashMap<String, u32>,
    /// group id -> (tid, start, end).
    span: HashMap<u32, (usize, u64, u64)>,
}

fn read_groups(
    path: &Path,
    name_to_tid: &HashMap<String, usize>,
    findings: Option<&Findings>,
    min_tier: FishTier,
) -> Result<Groups> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading haplotypes {}", path.display()))?;
    // group -> (tid, start, end)
    let mut span: HashMap<u32, (usize, u64, u64)> = HashMap::new();
    let mut rows: Vec<(u32, String)> = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.starts_with("group\t") || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            continue;
        }
        let Ok(group) = f[0].parse::<u32>() else { continue };
        let Some(&tid) = name_to_tid.get(f[1]) else { continue };
        let mut it = f[2].split('-');
        let (Some(Ok(s)), Some(Ok(e))) = (it.next().map(str::parse::<u64>), it.next().map(str::parse::<u64>))
        else {
            continue;
        };
        span.entry(group).or_insert((tid, s, e));
        rows.push((group, f[3].to_string()));
    }
    // Which groups count? `findings` scope: the span must overlap a tier-qualified
    // finding; genome scope: every group.
    let mut selected: HashSet<u32> = HashSet::new();
    for (&g, &(tid, s, e)) in &span {
        let ok = match findings {
            None => true,
            Some(fd) => fd
                .all
                .iter()
                .any(|f| f.tier >= min_tier && f.tid == tid && s < f.end && f.start < e),
        };
        if ok {
            selected.insert(g);
        }
    }
    let mut read_group: HashMap<String, u32> = HashMap::new();
    for (group, read) in rows {
        if !selected.contains(&group) {
            continue;
        }
        read_group
            .entry(read)
            .and_modify(|g| {
                if group < *g {
                    *g = group;
                }
            })
            .or_insert(group);
    }
    // Keep only the spans of selected groups.
    span.retain(|g, _| selected.contains(g));
    Ok(Groups { read_group, span })
}

// ---------------------------------------------------------------------------
// record processing

/// Shared, immutable per-run context.
struct Ctx<'a> {
    params: &'a FishParams,
    /// `clip`/`large_*` limited to findings (scope `findings`).
    scope_findings: bool,
    contig_name: &'a [String],
    contig_len: &'a [usize],
    /// Per contig: circular (not in `--linear`).
    circular: &'a [bool],
    findings: Option<&'a Findings>,
    groups: Option<&'a Groups>,
}

impl Ctx<'_> {
    /// Findings an event at 0-based `pos` on `tid` is attributed to.
    ///
    /// `Some(ids)` lists every applicable finding (ascending); an empty list
    /// means the event lies in none (attribute `.`). `None` means the event is
    /// out of scope entirely: in the `findings` scope only findings of tier >=
    /// `--min-tier` count, so an event inside no such finding is dropped. In the
    /// `genome` scope every finding counts (any tier) and the event is kept.
    fn event_attr(&self, tid: usize, pos: u64) -> Option<Vec<u32>> {
        let min = if self.scope_findings { Some(self.params.min_tier) } else { None };
        let ids = self
            .findings
            .map(|f| f.all_at(tid, pos, self.params.finding_pad, min))
            .unwrap_or_default();
        if self.scope_findings && ids.is_empty() {
            None
        } else {
            Some(ids)
        }
    }

    /// The tier-qualified findings a `locus` record `[start, end)` overlaps
    /// (± `--locus-flank`).
    fn locus_attr(&self, tid: usize, start: u64, end: u64) -> Vec<u32> {
        self.findings
            .map(|f| f.all_overlap(tid, start, end, self.params.locus_flank, Some(self.params.min_tier)))
            .unwrap_or_default()
    }

    /// The findings a haplotype-group span `[start, end)` overlaps: tier-qualified
    /// ones in the `findings` scope, any in the `genome` scope.
    fn linked_attr(&self, tid: usize, start: u64, end: u64) -> Vec<u32> {
        let min = if self.scope_findings { Some(self.params.min_tier) } else { None };
        self.findings.map(|f| f.all_overlap(tid, start, end, 0, min)).unwrap_or_default()
    }
}

/// The record's `SA` targets as `contig:pos:strand` joined by `,`, `.` if none.
fn sa_string(rec: &bam::Record) -> String {
    match rec.aux(b"SA") {
        Ok(Aux::String(s)) => {
            let mut v: Vec<String> = Vec::new();
            for entry in s.split(';') {
                if entry.is_empty() {
                    continue;
                }
                let mut f = entry.split(',');
                let name = f.next().unwrap_or("");
                let pos = f.next().unwrap_or("");
                let strand = f.next().unwrap_or("");
                if name.is_empty() || pos.is_empty() {
                    continue;
                }
                v.push(format!("{name}:{pos}:{strand}"));
            }
            if v.is_empty() { ".".to_string() } else { v.join(",") }
        }
        _ => ".".to_string(),
    }
}

/// Classify and test one record. `tid` is `usize::MAX` for unplaced unmapped
/// records (the unmapped branch returns before using it).
///
/// `uncounted_primary` collects the names of reads whose *primary* record is
/// duplicate (0x400) or QC-fail (0x200): `lift` streams such primaries to
/// minimap2 anyway, so their remapped record in `P.bait.cram` loses the flag and
/// would look counted. Bait selection drops these reads so it stays consistent
/// with the other categories (and the FASTQ payload pass, which ignores the
/// flags).
fn process_record(
    rec: &bam::Record,
    tid: usize,
    ctx: &Ctx,
    aggs: &mut HashMap<String, ReadAgg>,
    counts: &mut InputCounts,
    uncounted_primary: &mut HashSet<String>,
) -> Result<()> {
    let flags = rec.flags();
    counts.total += 1;
    if flags & 0x100 != 0 {
        counts.secondary += 1;
        return Ok(());
    }
    if flags & 0x200 != 0 {
        counts.qcfail += 1;
        if flags & 0x800 == 0 {
            uncounted_primary.insert(String::from_utf8_lossy(rec.qname()).into_owned());
        }
        return Ok(());
    }
    if flags & 0x400 != 0 {
        counts.duplicate += 1;
        if flags & 0x800 == 0 {
            uncounted_primary.insert(String::from_utf8_lossy(rec.qname()).into_owned());
        }
        return Ok(());
    }
    if rec.raw_cigar().iter().any(|&c| op_parts(c).0 == OP_H) {
        bail!(
            "fish requires a hard-clip-free alignment (run minimap2 with -Y): record {} has a hard clip",
            String::from_utf8_lossy(rec.qname())
        );
    }
    let primary = flags & 0x800 == 0;
    let name = String::from_utf8_lossy(rec.qname()).into_owned();

    if flags & 0x4 != 0 {
        counts.unmapped += 1;
        counts.unmapped_bp += rec.seq_len() as u64;
        if !ctx.params.no_unmapped {
            let a = agg_of(aggs, &name, rec);
            let len = a.read_len;
            a.offer(
                R_UNMAPPED,
                Event {
                    len,
                    pos: 0,
                    primary,
                    strand: '.',
                    detail: ".".to_string(),
                    finding: None,
                    contig: ".".to_string(),
                    sa: ".".to_string(),
                },
            );
        }
        return Ok(());
    }

    if primary {
        counts.primary += 1;
        counts.mapped_bp += rec.seq_len() as u64;
    } else {
        counts.supplementary += 1;
    }

    let rpos = rec.pos();
    if rpos < 0 {
        return Ok(());
    }
    let pos = rpos as usize;
    let cig = rec.raw_cigar();
    let mut span = 0usize;
    for &c in cig {
        let (op, l) = op_parts(c);
        if matches!(op, OP_M | OP_D | OP_N | OP_EQ | OP_X) {
            span += l;
        }
    }
    if span == 0 {
        return Ok(());
    }
    let end = pos + span;
    let contig_len = ctx.contig_len.get(tid).copied().unwrap_or(0);
    let circular = ctx.circular.get(tid).copied().unwrap_or(false);
    let contig = ctx.contig_name.get(tid).cloned().unwrap_or_default();
    let (lead, trail) = terminal_clips(cig);
    let lclip = lead > 0;
    let rclip = cig.len() > 1 && trail > 0;
    let rev = flags & 0x10 != 0;
    let strand = if rev { '-' } else { '+' };
    let sa = sa_string(rec);
    let read_len = rec.seq_len() as u64;

    let mut evs: Vec<(&'static str, Event)> = Vec::new();
    // Offer a base event to every applicable finding (or `.` when none applies).
    let offerings = |reason: &'static str, base: Event, ids: Option<Vec<u32>>, evs: &mut Vec<(&'static str, Event)>| {
        match ids {
            None => {}
            Some(ids) if ids.is_empty() => {
                let mut e = base;
                e.finding = None;
                evs.push((reason, e));
            }
            Some(ids) => {
                for id in ids {
                    let mut e = base.clone();
                    e.finding = Some(id);
                    evs.push((reason, e));
                }
            }
        }
    };

    // clip: a terminal soft clip >= clip_min_len; a contig-end clip is an origin
    // crossing on a circular contig (skipped) but counts on a linear one.
    if !ctx.params.no_clips {
        if lclip && lead >= ctx.params.clip_min_len && (pos != 0 || !circular) {
            let p = pos as u64;
            let base = Event {
                len: lead as u64,
                pos: p,
                primary,
                strand,
                detail: format!("side=left;len={lead}"),
                finding: None,
                contig: contig.clone(),
                sa: sa.clone(),
            };
            offerings(R_CLIP, base, ctx.event_attr(tid, p), &mut evs);
        }
        if rclip && trail >= ctx.params.clip_min_len && (end != contig_len || !circular) {
            let p = (end - 1) as u64;
            let base = Event {
                len: trail as u64,
                pos: p,
                primary,
                strand,
                detail: format!("side=right;len={trail}"),
                finding: None,
                contig: contig.clone(),
                sa: sa.clone(),
            };
            offerings(R_CLIP, base, ctx.event_attr(tid, p), &mut evs);
        }
    }

    // merged large indels
    if !ctx.params.no_indels {
        let mut raw: Vec<RawIndel> = Vec::new();
        raw_indels(cig, pos, ctx.params.indel_piece_min, &mut raw);
        for (ins, s, _e, size) in merge_large(&raw, ctx.params.indel_merge_gap, ctx.params.indel_min_len) {
            let p = s as u64;
            let base = Event {
                len: size as u64,
                pos: p,
                primary,
                strand,
                detail: format!("len={size};at={s}"),
                finding: None,
                contig: contig.clone(),
                sa: sa.clone(),
            };
            offerings(if ins { R_INS } else { R_DEL }, base, ctx.event_attr(tid, p), &mut evs);
        }
    }

    // linked: the read belongs to a selected haplotype group. The row is offered
    // to every finding overlapping the group span (tier-qualified ones in the
    // findings scope, any in the genome scope), or `.` when none overlaps.
    if !ctx.params.no_linked {
        if let Some(g) = ctx.groups.and_then(|gs| gs.read_group.get(&name)).copied() {
            let (gtid, gs, ge) = ctx.groups.unwrap().span[&g];
            let base = Event {
                len: read_len,
                pos: gs,
                primary,
                strand,
                detail: format!("group={g}"),
                finding: None,
                contig: ctx.contig_name.get(gtid).cloned().unwrap_or_default(),
                sa: sa.clone(),
            };
            offerings(R_LINKED, base, Some(ctx.linked_attr(gtid, gs, ge)), &mut evs);
        }
    }

    // locus: any counted record overlapping a finding ± locus_flank
    if ctx.params.locus {
        let ids = ctx.locus_attr(tid, pos as u64, end as u64);
        for id in ids {
            evs.push((
                R_LOCUS,
                Event {
                    len: read_len,
                    pos: pos as u64,
                    primary,
                    strand,
                    detail: ".".to_string(),
                    finding: Some(id),
                    contig: contig.clone(),
                    sa: sa.clone(),
                },
            ));
        }
    }

    if !evs.is_empty() {
        let a = agg_of(aggs, &name, rec);
        for (reason, ev) in evs {
            a.offer(reason, ev);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// bait (lift-based)

/// The best bait record of one read: the `Event` merged into the read's
/// aggregate plus the record's MAPQ (kept for the JSON's MAPQ-0 count).
struct BaitHit {
    event: Event,
    mapq: u8,
}

/// Is `a` the better bait record than `b`? Primary over supplementary, then
/// highest MAPQ, then a deterministic tie-break so the choice never depends on
/// the order the records are read.
fn bait_better(a: &BaitHit, b: &BaitHit) -> bool {
    if a.event.primary != b.event.primary {
        return a.event.primary;
    }
    if a.mapq != b.mapq {
        return a.mapq > b.mapq;
    }
    if a.event.contig != b.event.contig {
        return a.event.contig < b.event.contig;
    }
    if a.event.pos != b.event.pos {
        return a.event.pos < b.event.pos;
    }
    a.event.strand < b.event.strand
}

/// Format an optional integer for the `bait` detail (`.` when absent).
fn fmt_int(v: Option<i64>) -> String {
    v.map_or_else(|| ".".to_string(), |n| n.to_string())
}

/// Load and validate the bait contigs without writing anything.
///
/// Bait contig names must not clash with an assembly contig (or another bait
/// contig): a clash is an error, names are never renamed silently. Returning the
/// assembly text and the bait sequences lets the caller validate before it
/// touches any output file. `Err` here leaves every existing output intact.
fn plan_bait_ref(
    assembly: &Path,
    baits: &[PathBuf],
) -> Result<(String, Vec<(String, Vec<u8>)>)> {
    let text = std::fs::read_to_string(assembly)
        .with_context(|| format!("reading assembly FASTA {}", assembly.display()))?;
    let mut seen: HashSet<String> = HashSet::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            let name = rest.split_whitespace().next().unwrap_or("");
            if name.is_empty() {
                bail!("assembly FASTA {}: empty sequence name", assembly.display());
            }
            if !seen.insert(name.to_string()) {
                bail!("assembly FASTA {}: duplicate sequence name {name:?}", assembly.display());
            }
        }
    }
    let mut bait: Vec<(String, Vec<u8>)> = Vec::new();
    for p in baits {
        let b = Fasta::load(p)?;
        for name in b.order() {
            if !seen.insert(name.clone()) {
                bail!(
                    "bait contig {name:?} ({}) clashes with an assembly or earlier bait contig; \
                     rename it (names are never changed silently)",
                    p.display()
                );
            }
            bait.push((name.clone(), b.seq(name).unwrap().to_vec()));
        }
    }
    Ok((text, bait))
}

/// Write `P.bait_ref.fa` = the assembly verbatim, then the bait contigs.
///
/// The file is staged and only renamed into place on success, so a write failure
/// leaves the previous reference intact. Any `.fai` of the target is removed and
/// recreated for the new content: `lift`'s `ensure_fai` only creates an index
/// when none exists, so a stale index would otherwise be reused on a rerun with
/// different bait. Returns the bait `(name, length)` list in output order.
fn write_bait_ref(
    assembly_text: &str,
    bait: &[(String, Vec<u8>)],
    out: &Path,
) -> Result<Vec<(String, u64)>> {
    let staged = out_path(out, ".tmp");
    let res = (|| -> Result<()> {
        let f = std::fs::File::create(&staged)
            .with_context(|| format!("creating bait reference {}", staged.display()))?;
        let mut w = std::io::BufWriter::new(f);
        for line in assembly_text.lines() {
            w.write_all(line.as_bytes())?;
            w.write_all(b"\n")?;
        }
        for (name, seq) in bait {
            writeln!(w, ">{name}")?;
            w.write_all(seq)?;
            w.write_all(b"\n")?;
        }
        w.flush()?;
        Ok(())
    })();
    if let Err(e) = res {
        let _ = std::fs::remove_file(&staged);
        return Err(e);
    }
    std::fs::rename(&staged, out)
        .with_context(|| format!("renaming {} to {}", staged.display(), out.display()))?;
    // The index is now stale; recreate it for the new content (`lift`'s
    // `ensure_fai` only creates an index when none exists, so it must not see a
    // stale one).
    let fai = out_path(out, ".fai");
    match std::fs::remove_file(&fai) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(anyhow!("removing stale bait index {}: {e}", fai.display())),
    }
    crate::fasta::ensure_fai(out)?;
    Ok(bait.iter().map(|(n, s)| (n.clone(), s.len() as u64)).collect())
}

/// The previous `P.bait_ref.fa` and its `.fai`, moved aside while a new bait
/// reference is built and lifted. A failed lift restores them, so the previous
/// `P.bait.cram` is never left paired with a different reference (its `@SQ` MD5
/// would no longer match and decoding would break).
struct RefBackup {
    ref_path: PathBuf,
    fai_path: PathBuf,
    ref_bak: PathBuf,
    fai_bak: PathBuf,
    has_ref: bool,
    has_fai: bool,
}

impl RefBackup {
    /// Move the existing reference + index aside (if present).
    fn stage(ref_path: &Path) -> Result<RefBackup> {
        let fai_path = out_path(ref_path, ".fai");
        let mut b = RefBackup {
            ref_path: ref_path.to_path_buf(),
            fai_path: fai_path.clone(),
            ref_bak: out_path(ref_path, ".prev"),
            fai_bak: out_path(&fai_path, ".prev"),
            has_ref: false,
            has_fai: false,
        };
        if b.ref_path.is_file() {
            std::fs::rename(&b.ref_path, &b.ref_bak)
                .with_context(|| format!("backing up {}", b.ref_path.display()))?;
            b.has_ref = true;
        }
        if b.fai_path.is_file() {
            if let Err(e) = std::fs::rename(&b.fai_path, &b.fai_bak) {
                if b.has_ref {
                    let _ = std::fs::rename(&b.ref_bak, &b.ref_path);
                }
                return Err(anyhow!("backing up {}: {e}", b.fai_path.display()));
            }
            b.has_fai = true;
        }
        Ok(b)
    }

    /// Remove whatever the failed run wrote and put the previous pair back.
    fn restore(&self) {
        let _ = std::fs::remove_file(&self.ref_path);
        let _ = std::fs::remove_file(&self.fai_path);
        if self.has_ref {
            let _ = std::fs::rename(&self.ref_bak, &self.ref_path);
        }
        if self.has_fai {
            let _ = std::fs::rename(&self.fai_bak, &self.fai_path);
        }
    }

    /// Drop the backups after a successful lift.
    fn discard(&self) {
        if self.has_ref {
            let _ = std::fs::remove_file(&self.ref_bak);
        }
        if self.has_fai {
            let _ = std::fs::remove_file(&self.fai_bak);
        }
    }
}

/// Read `P.bait.cram` and return the best bait record per read. Counted records
/// only (secondary 0x100, QC-fail 0x200 and duplicate 0x400 are skipped);
/// primary and supplementary both count.
fn read_bait_records(
    cram: &Path,
    bait_ref: &Path,
    bait: &[(String, u64)],
) -> Result<HashMap<String, BaitHit>> {
    let mut reader = bam::IndexedReader::from_path(cram)
        .with_context(|| format!("opening {}", cram.display()))?;
    reader
        .set_reference(bait_ref)
        .with_context(|| format!("setting reference {} for {}", bait_ref.display(), cram.display()))?;
    let tid_of: HashMap<String, i32> = {
        let hv = reader.header().clone();
        hv.target_names()
            .iter()
            .enumerate()
            .map(|(i, n)| (String::from_utf8_lossy(n).into_owned(), i as i32))
            .collect()
    };
    let mut best: HashMap<String, BaitHit> = HashMap::new();
    let mut rec = bam::Record::new();
    for (name, _len) in bait {
        let Some(&tid) = tid_of.get(name) else {
            bail!("bait contig {name:?} is missing from {}", cram.display());
        };
        reader
            .fetch(tid)
            .with_context(|| format!("fetching bait contig {name} from {}", cram.display()))?;
        while let Some(res) = reader.read(&mut rec) {
            res.with_context(|| format!("reading {}", cram.display()))?;
            let flags = rec.flags();
            if flags & (0x100 | 0x200 | 0x400) != 0 {
                continue;
            }
            let pos = rec.pos();
            if pos < 0 {
                continue;
            }
            let primary = flags & 0x800 == 0;
            let mapq = rec.mapq();
            let strand = if flags & 0x10 != 0 { '-' } else { '+' };
            let record = if primary { "primary" } else { "supplementary" };
            let detail = format!(
                "contig={name};mapq={mapq};s1={};s2={};record={record}",
                fmt_int(aux_int(&rec, b"s1")),
                fmt_int(aux_int(&rec, b"s2")),
            );
            let qname = String::from_utf8_lossy(rec.qname()).into_owned();
            let cand = BaitHit {
                event: Event {
                    len: rec.seq_len() as u64,
                    pos: pos as u64,
                    primary,
                    strand,
                    detail,
                    finding: None,
                    contig: name.clone(),
                    sa: sa_string(&rec),
                },
                mapq,
            };
            match best.get(&qname) {
                Some(old) if !bait_better(&cand, old) => {}
                _ => {
                    best.insert(qname, cand);
                }
            }
        }
    }
    Ok(best)
}

/// The `bait` step: build `P.bait_ref.fa`, lift the input onto it into
/// `P.bait.cram` (same defaults as `cramhook lift`, except the tiny-input
/// threshold), and collect every read placed on a bait contig. Reads whose
/// original primary is duplicate/QC-fail are skipped (`uncounted`): lift remaps
/// them and their remapped record would otherwise look counted. `records`
/// selects which bait records count: `Any` keeps primary and supplementary
/// records, `Primary` keeps only reads whose primary record is on a bait contig.
/// Returns the per-read hits and the JSON summary.
fn run_bait(
    input: &Path,
    ref_path: &Path,
    out_prefix: &Path,
    baits: &[PathBuf],
    min_input_mb: Option<f64>,
    records: BaitRecords,
    uncounted: &HashSet<String>,
    threads: usize,
) -> Result<(HashMap<String, BaitHit>, BaitSummary)> {
    let t0 = Instant::now();
    // Validate every name before touching any output.
    let (assembly_text, bait_seqs) = plan_bait_ref(ref_path, baits)?;
    let bait_ref = out_path(out_prefix, ".bait_ref.fa");
    // Move the previous reference + index aside so a failed lift can restore the
    // pair that matches the previous `P.bait.cram`.
    let backup = RefBackup::stage(&bait_ref)?;
    let bait_contigs = match write_bait_ref(&assembly_text, &bait_seqs, &bait_ref) {
        Ok(v) => v,
        Err(e) => {
            backup.restore();
            return Err(e);
        }
    };
    let bait_cram = out_path(out_prefix, ".bait.cram");
    let tcfg = TranslateConfig::default();
    let opts = LiftOptions { strict: false, keep_tmp: false, min_input_mb };
    let stats = match lift(input, ref_path, &bait_ref, &bait_cram, None, threads, &tcfg, &opts, None)
    {
        Ok(s) => s,
        Err(e) => {
            backup.restore();
            return Err(e.context(format!(
                "fish --bait: lifting {} onto {} failed",
                input.display(),
                bait_ref.display()
            )));
        }
    };
    backup.discard();
    let mut hits = read_bait_records(&bait_cram, &bait_ref, &bait_contigs)?;
    hits.retain(|name, _| !uncounted.contains(name));

    // A read whose bait support is only a supplementary record is selected under
    // `any` but not under `primary`; the count is reported in both modes.
    let supplementary_only_reads = hits.values().filter(|h| !h.event.primary).count() as u64;
    if records == BaitRecords::Primary {
        hits.retain(|_, h| h.event.primary);
    }

    let index: HashMap<&str, usize> =
        bait_contigs.iter().enumerate().map(|(i, (n, _))| (n.as_str(), i)).collect();
    let mut per: Vec<(u64, u64, u64)> = vec![(0, 0, 0); bait_contigs.len()];
    let mut reads = 0u64;
    let mut bp = 0u64;
    let mut mapq0 = 0u64;
    for hit in hits.values() {
        reads += 1;
        bp += hit.event.len;
        let m0 = u64::from(hit.mapq == 0);
        mapq0 += m0;
        if let Some(&i) = index.get(hit.event.contig.as_str()) {
            per[i].0 += 1;
            per[i].1 += hit.event.len;
            per[i].2 += m0;
        }
    }
    let contigs: Vec<BaitContig> = bait_contigs
        .iter()
        .zip(per)
        .map(|((name, length), (r, b, m0))| BaitContig {
            name: name.clone(),
            length: *length,
            reads: r,
            bp: b,
            mapq0_reads: m0,
        })
        .collect();
    let summary = BaitSummary {
        contigs,
        reads,
        bp,
        mapq0_reads: mapq0,
        mode: records.name().to_string(),
        supplementary_only_reads,
        translated_reads: stats.translated_reads,
        remapped_reads: stats.remapped_reads,
        contract_ok: stats.contract_ok,
        contract_reason: stats.contract_reason.clone(),
        fallback: stats.fallback.as_ref().map(|f| f.reason.clone()),
        timing_s: t0.elapsed().as_secs_f64(),
    };
    Ok((hits, summary))
}

// ---------------------------------------------------------------------------
// I/O

fn out_path(prefix: &Path, suffix: &str) -> PathBuf {
    crate::translate::suffix(prefix, suffix)
}

/// One output row of `P.fish.tsv`.
struct Row {
    read: String,
    len: u64,
    reason: &'static str,
    finding: Option<u32>,
    contig: String,
    pos: u64,
    strand: char,
    record: &'static str,
    detail: String,
    sa: String,
}

enum Task {
    Region { tid: usize, start: usize, end: usize },
    Unplaced,
}

/// Second pass: capture the FASTQ payload of each selected read. One reader per
/// thread over the same chunk tasks; only selected reads are decoded. The donor
/// is the smallest `(rank, tid, pos)` among the read's counted records, so the
/// payload does not depend on the thread count.
fn run_payload_pass(
    input: &Path,
    ref_path: &Path,
    tasks: &[Task],
    selected: &HashSet<String>,
    params: &FishParams,
    names: &[String],
) -> Result<HashMap<String, Payload>> {
    let next = AtomicUsize::new(0);
    let nthreads = params.threads.max(1).min(tasks.len());
    let workers: Vec<Result<HashMap<String, Payload>>> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..nthreads)
            .map(|_| {
                let tasks = tasks;
                let next = &next;
                let names = names;
                sc.spawn(move || -> Result<HashMap<String, Payload>> {
                    let mut reader = open_indexed(input, ref_path, FIELDS_FISH_MIN)?;
                    let mut rec = bam::Record::new();
                    let mut out: HashMap<String, Payload> = HashMap::new();
                    loop {
                        let k = next.fetch_add(1, Ordering::Relaxed);
                        let Some(task) = tasks.get(k) else { break };
                        match *task {
                            Task::Region { tid, start, end } => {
                                reader
                                    .fetch((tid as i32, start as i64, end as i64))
                                    .with_context(|| {
                                        format!("fetching {}:{start}-{end}", names[tid])
                                    })?;
                                while let Some(res) = reader.read(&mut rec) {
                                    res.with_context(|| format!("reading {}", input.display()))?;
                                    let rp = rec.pos();
                                    if rp < start as i64 || rp >= end as i64 {
                                        continue;
                                    }
                                    note_payload(&mut out, selected, &rec);
                                }
                            }
                            Task::Unplaced => {
                                reader
                                    .fetch(bam::FetchDefinition::Unmapped)
                                    .context("fetching coordinate-less records")?;
                                while let Some(res) = reader.read(&mut rec) {
                                    res.with_context(|| format!("reading {}", input.display()))?;
                                    note_payload(&mut out, selected, &rec);
                                }
                            }
                        }
                    }
                    Ok(out)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap_or_else(|_| Err(anyhow!("fish payload worker panicked"))))
            .collect()
    });
    let mut payloads: HashMap<String, Payload> = HashMap::new();
    for w in workers {
        for (n, p) in w? {
            match payloads.get(&n) {
                Some(e) if (e.rank, e.tid, e.pos) <= (p.rank, p.tid, p.pos) => {}
                _ => {
                    payloads.insert(n, p);
                }
            }
        }
    }
    Ok(payloads)
}

/// The record count htslib's index reports (per-contig mapped + unmapped plus
/// the no-coordinate count), or `None` (with a note) when index stats are
/// unavailable. Cheap for a BAM; for a CRAM it does a short sequential scan
/// (a `.crai` has no per-contig stats).
fn completeness_expected(
    input: &Path,
    ref_path: &Path,
    notes: &mut Vec<String>,
) -> Result<Option<u64>> {
    let mut r = open_indexed(input, ref_path, FIELDS_FISH_FULL)?;
    match r.index_stats() {
        Ok(stats) => Ok(Some(stats.iter().map(|&(_, _, m, u)| m + u).sum())),
        Err(e) => {
            notes.push(format!("completeness check skipped: {e}"));
            Ok(None)
        }
    }
}

/// Run `cramhook fish`. Returns the summary (also written to `P.fish.json`).
pub fn run_fish(input: &Path, ref_path: &Path, out_prefix: &Path, params: &FishParams) -> Result<FishSummary> {
    let t0 = Instant::now();
    crate::translate::ensure_parent(out_prefix)?;

    // Resolve the scope and validate the inputs.
    let scope = match params.scope {
        Some(s) => s,
        None => {
            if params.findings.is_some() {
                FishScope::Findings
            } else {
                FishScope::Genome
            }
        }
    };
    if scope == FishScope::Findings && params.findings.is_none() {
        bail!("fish --scope findings requires --findings");
    }
    if params.locus && params.findings.is_none() {
        bail!("fish --locus requires --findings");
    }
    let mut notes: Vec<String> = Vec::new();
    let linked_on = !params.no_linked && params.haplotypes.is_some();
    if !params.no_linked && params.haplotypes.is_none() {
        notes.push("linked: skipped (no --haplotypes)".to_string());
    }

    // Header: contig names and lengths.
    let r0 = open_indexed(input, ref_path, FIELDS_FISH_FULL)?;
    let hv = r0.header().clone();
    drop(r0);
    let contig_name: Vec<String> =
        hv.target_names().iter().map(|n| String::from_utf8_lossy(n).into_owned()).collect();
    let mut contig_len: Vec<usize> = Vec::with_capacity(contig_name.len());
    for (tid, _) in contig_name.iter().enumerate() {
        contig_len.push(hv.target_len(tid as u32).unwrap_or(0) as usize);
    }
    let linear: std::collections::HashSet<&str> = params.linear.iter().map(|s| s.as_str()).collect();
    let circular: Vec<bool> = contig_name.iter().map(|n| !linear.contains(n.as_str())).collect();
    let name_to_tid: HashMap<String, usize> =
        contig_name.iter().enumerate().map(|(i, n)| (n.clone(), i)).collect();

    let findings = match &params.findings {
        Some(p) => Some(read_findings(p, &name_to_tid)?),
        None => None,
    };
    let groups = if linked_on {
        // Only the findings scope filters groups by findings; with `--scope
        // genome` every group counts even when `--findings` is given.
        let filter = if scope == FishScope::Findings { findings.as_ref() } else { None };
        Some(read_groups(
            params.haplotypes.as_ref().unwrap(),
            &name_to_tid,
            filter,
            params.min_tier,
        )?)
    } else {
        None
    };

    let ctx = Ctx {
        params,
        scope_findings: scope == FishScope::Findings,
        contig_name: &contig_name,
        contig_len: &contig_len,
        circular: &circular,
        findings: findings.as_ref(),
        groups: groups.as_ref(),
    };

    // Tasks: contig chunks plus one coordinate-less task.
    let total: u64 = contig_len.iter().map(|&l| l as u64).sum();
    let cs = params
        .chunk_size
        .unwrap_or_else(|| (total / (4 * params.threads.max(1) as u64)).clamp(200_000, 1_000_000))
        .max(1) as usize;
    let mut tasks: Vec<Task> = Vec::new();
    for (tid, &len) in contig_len.iter().enumerate() {
        let mut a = 0;
        while a < len {
            let b = (a + cs).min(len);
            tasks.push(Task::Region { tid, start: a, end: b });
            a = b;
        }
    }
    tasks.push(Task::Unplaced);

    let next = AtomicUsize::new(0);
    let nthreads = params.threads.max(1).min(tasks.len());
    let workers: Vec<Result<(HashMap<String, ReadAgg>, InputCounts, HashSet<String>)>> =
        std::thread::scope(|sc| {
            let handles: Vec<_> = (0..nthreads)
                .map(|_| {
                    let ctx = &ctx;
                    let tasks = &tasks;
                    let next = &next;
                    let names = &contig_name;
                    sc.spawn(move || -> Result<(HashMap<String, ReadAgg>, InputCounts, HashSet<String>)> {
                        let mut reader = open_indexed(input, ref_path, FIELDS_FISH_FULL)?;
                        let mut rec = bam::Record::new();
                        let mut aggs: HashMap<String, ReadAgg> = HashMap::new();
                        let mut counts = InputCounts::default();
                        let mut uncounted: HashSet<String> = HashSet::new();
                        loop {
                            let k = next.fetch_add(1, Ordering::Relaxed);
                            let Some(task) = tasks.get(k) else { break };
                            match *task {
                                Task::Region { tid, start, end } => {
                                    reader
                                        .fetch((tid as i32, start as i64, end as i64))
                                        .with_context(|| {
                                            format!("fetching {}:{start}-{end}", names[tid])
                                        })?;
                                    while let Some(res) = reader.read(&mut rec) {
                                        res.with_context(|| format!("reading {}", input.display()))?;
                                        let rp = rec.pos();
                                        if rp < start as i64 || rp >= end as i64 {
                                            continue;
                                        }
                                        process_record(
                                            &rec,
                                            tid,
                                            ctx,
                                            &mut aggs,
                                            &mut counts,
                                            &mut uncounted,
                                        )?;
                                    }
                                }
                                Task::Unplaced => {
                                    // `FetchDefinition::Unmapped` (`"*"`) returns
                                    // exactly the coordinate-less records (pos < 0),
                                    // including tid >= 0 ones; placed-unmapped records
                                    // (pos >= 0) only come from the region tasks. So
                                    // process every record here and do not skip tid >= 0.
                                    reader
                                        .fetch(bam::FetchDefinition::Unmapped)
                                        .context("fetching coordinate-less records")?;
                                    while let Some(res) = reader.read(&mut rec) {
                                        res.with_context(|| format!("reading {}", input.display()))?;
                                        process_record(
                                            &rec,
                                            usize::MAX,
                                            ctx,
                                            &mut aggs,
                                            &mut counts,
                                            &mut uncounted,
                                        )?;
                                    }
                                }
                            }
                        }
                        Ok((aggs, counts, uncounted))
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap_or_else(|_| Err(anyhow!("fish worker panicked")))).collect()
        });

    // Merge workers deterministically (sum counters, fold per-read events).
    let mut counts = InputCounts::default();
    let mut global: HashMap<String, ReadAgg> = HashMap::new();
    let mut uncounted_primary: HashSet<String> = HashSet::new();
    for w in workers {
        let (aggs, c, uncounted) = w?;
        counts.add(&c);
        uncounted_primary.extend(uncounted);
        for (name, agg) in aggs {
            let entry = global.entry(name).or_default();
            if agg.read_len > entry.read_len {
                entry.read_len = agg.read_len;
            }
            for ((reason, _), ev) in agg.events {
                entry.offer(reason, ev);
            }
        }
    }

    // Bait (optional): lift the input onto assembly + bait and select every read
    // that lands on a bait contig. Independent of the categories above, which
    // were all evaluated on the original `--in` alignment.
    let mut bait_summary = None;
    if !params.bait.is_empty() {
        let (hits, summary) = run_bait(
            input,
            ref_path,
            out_prefix,
            &params.bait,
            params.bait_min_input_mb,
            params.bait_records,
            &uncounted_primary,
            params.threads,
        )?;
        for (name, hit) in hits {
            let a = global.entry(name).or_default();
            if hit.event.len > a.read_len {
                a.read_len = hit.event.len;
            }
            a.offer(R_BAIT, hit.event);
        }
        bait_summary = Some(summary);
    }

    // Completeness: the region fetch and the `"*"` iterator together must see
    // every record the index counts. htslib index stats are cheap for a BAM and
    // a short sequential scan for a CRAM; a mismatch means the file holds records
    // no iterator can return (e.g. flag 4 with a set tid and pos -1), which fish
    // cannot support, so fail loudly instead of dropping them.
    if let Some(exp) = completeness_expected(input, ref_path, &mut notes)? {
        if exp != counts.total {
            bail!(
                "fish: the index reports {exp} records but only {} could be fetched; \
                 the input holds records the coordinate/unmapped iterators do not \
                 return (e.g. flag 4 with a tid set and pos -1). Re-sort the alignment \
                 so every record has a usable coordinate.",
                counts.total
            );
        }
    }

    let mut selected: Vec<(String, ReadAgg)> =
        global.into_iter().filter(|(_, a)| !a.events.is_empty()).collect();
    selected.sort_by(|a, b| a.0.cmp(&b.0));

    // Second pass: capture the FASTQ payload of each selected read. Only the
    // selected reads are decoded, so memory stays proportional to the selection
    // rather than to the whole input.
    if !selected.is_empty() {
        let names: HashSet<String> = selected.iter().map(|(n, _)| n.clone()).collect();
        let payloads = run_payload_pass(input, ref_path, &tasks, &names, params, &contig_name)?;
        for (name, agg) in selected.iter_mut() {
            match payloads.get(name) {
                Some(p) => agg.payload = Some(p.clone()),
                None => bail!(
                    "fish: selected read {name} has no sequence payload (input may lack SEQ/QUAL; run minimap2 with -Y)"
                ),
            }
        }
    }

    // Rows: one per (read, reason, finding), keeping the largest event.
    let mut rows: Vec<Row> = Vec::new();
    let mut reads = FishReads::default();
    for (name, agg) in &selected {
        reads.reads += 1;
        reads.bp += agg.read_len;
        for ((reason, _), ev) in &agg.events {
            rows.push(Row {
                read: name.clone(),
                len: ev.len,
                reason,
                finding: ev.finding,
                contig: ev.contig.clone(),
                pos: ev.pos,
                strand: ev.strand,
                record: ev.primary_or_supp(),
                detail: ev.detail.clone(),
                sa: ev.sa.clone(),
            });
        }
    }
    // Sort read, then reason, then finding (numeric; `.` last).
    rows.sort_by(|a, b| {
        (a.read.as_str(), a.reason, a.finding.is_none(), a.finding.unwrap_or(0)).cmp(&(
            b.read.as_str(),
            b.reason,
            b.finding.is_none(),
            b.finding.unwrap_or(0),
        ))
    });
    // Per-reason counts are per distinct read, not per row.
    for (_, agg) in &selected {
        let mut reasons: Vec<&'static str> = agg.events.keys().map(|(r, _)| *r).collect();
        reasons.sort_unstable();
        reasons.dedup();
        for reason in reasons {
            let e = reads.by_reason.entry(reason.to_string()).or_default();
            e.reads += 1;
            e.bp += agg.read_len;
        }
    }

    // Write the FASTQ (BGZF), sorted by read name. A parallel BGZF thread pool
    // roughly halves the encode tail; the decompressed bytes do not depend on it.
    let fq_path = out_path(out_prefix, ".fastq.gz");
    {
        let pool = if params.threads > 1 {
            rust_htslib::tpool::ThreadPool::new(params.threads as u32).ok()
        } else {
            None
        };
        let mut w = rust_htslib::bgzf::Writer::from_path(&fq_path)
            .with_context(|| format!("creating {}", fq_path.display()))?;
        if let Some(p) = &pool {
            w.set_thread_pool(p).context("attaching the BGZF thread pool")?;
        }
        for (name, agg) in &selected {
            let Some(p) = &agg.payload else {
                bail!("fish: selected read {name} has no alignment record");
            };
            if p.seq.is_empty() {
                bail!(
                    "fish: selected read {name} has no sequence payload (input may lack SEQ; run minimap2 with -Y)"
                );
            }
            if p.seq.len() != p.qual.len() {
                bail!("fish: selected read {name} has inconsistent SEQ/QUAL lengths");
            }
            w.write_all(b"@")?;
            w.write_all(name.as_bytes())?;
            w.write_all(b"\n")?;
            w.write_all(&p.seq)?;
            w.write_all(b"\n+\n")?;
            w.write_all(&p.qual)?;
            w.write_all(b"\n")?;
        }
        w.flush().with_context(|| format!("flushing {}", fq_path.display()))?;
    }

    // Write the TSV.
    let tsv_path = out_path(out_prefix, ".fish.tsv");
    {
        let f = std::fs::File::create(&tsv_path)
            .with_context(|| format!("creating {}", tsv_path.display()))?;
        let mut w = std::io::BufWriter::new(f);
        writeln!(w, "read\tlength\treason\tfinding\tcontig\tpos\tstrand\trecord\tdetail\tsa")?;
        for r in &rows {
            let finding = match r.finding {
                Some(n) => format!("finding_{n}"),
                None => ".".to_string(),
            };
            if r.reason == R_UNMAPPED {
                writeln!(
                    w,
                    "{}\t{}\t{}\t.\t.\t.\t.\t.\t.\t.",
                    r.read, r.len, r.reason
                )?;
            } else {
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    r.read,
                    r.len,
                    r.reason,
                    finding,
                    r.contig,
                    r.pos,
                    r.strand,
                    r.record,
                    r.detail,
                    r.sa
                )?;
            }
        }
        w.flush()?;
    }

    let pass_s = t0.elapsed().as_secs_f64();
    let summary = FishSummary {
        reads,
        records: counts,
        params: FishParamsOut {
            scope: match scope {
                FishScope::Findings => "findings".to_string(),
                FishScope::Genome => "genome".to_string(),
            },
            min_tier: params.min_tier.name().to_string(),
            finding_pad: params.finding_pad,
            unmapped: !params.no_unmapped,
            clips: !params.no_clips,
            indels: !params.no_indels,
            linked: linked_on,
            bait: !params.bait.is_empty(),
            locus: params.locus,
            locus_flank: params.locus_flank,
            clip_min_len: params.clip_min_len,
            indel_min_len: params.indel_min_len,
            indel_merge_gap: params.indel_merge_gap,
            indel_piece_min: params.indel_piece_min,
        },
        bait: bait_summary,
        notes,
        timing_s: FishTiming { pass: pass_s, total: t0.elapsed().as_secs_f64() },
    };
    let json_path = out_path(out_prefix, ".fish.json");
    {
        let f = std::fs::File::create(&json_path)
            .with_context(|| format!("creating {}", json_path.display()))?;
        let mut w = std::io::BufWriter::new(f);
        serde_json::to_writer_pretty(&mut w, &versioned(&summary)).context("writing fish json")?;
        writeln!(w)?;
        w.flush()?;
    }
    Ok(summary)
}

impl Event {
    fn primary_or_supp(&self) -> &'static str {
        if self.primary { "primary" } else { "supplementary" }
    }
}

#[cfg(test)]
mod bait_unit {
    use super::*;

    fn hit(primary: bool, mapq: u8, contig: &str, pos: u64) -> BaitHit {
        BaitHit {
            event: Event {
                len: 100,
                pos,
                primary,
                strand: '+',
                detail: String::new(),
                finding: None,
                contig: contig.to_string(),
                sa: String::new(),
            },
            mapq,
        }
    }

    #[test]
    fn bait_prefers_primary_then_mapq() {
        // A primary wins over a supplementary even with a lower MAPQ.
        let prim_low = hit(true, 0, "b", 10);
        let supp_high = hit(false, 60, "b", 20);
        assert!(bait_better(&prim_low, &supp_high));
        assert!(!bait_better(&supp_high, &prim_low));

        // Same kind: the higher MAPQ wins.
        let mq40 = hit(true, 40, "b", 10);
        let mq60 = hit(true, 60, "b", 20);
        assert!(bait_better(&mq60, &mq40));
        assert!(!bait_better(&mq40, &mq60));

        // Deterministic tie-break: contig, then position, then strand.
        let a = hit(true, 60, "a", 10);
        let b = hit(true, 60, "b", 5);
        assert!(bait_better(&a, &b));
        let p1 = hit(true, 60, "a", 10);
        let p2 = hit(true, 60, "a", 20);
        assert!(bait_better(&p1, &p2));
    }
}
