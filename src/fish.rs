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
//! category: `unmapped`, `clip`, `hairpin`, `large_ins` / `large_del`, `linked`
//! (`--haplotypes`), `locus` (`--findings`) or `bait` (`--bait`). In the
//! `findings` scope, `clip`, `hairpin` and `large_*` only count when the event
//! position lies in a finding of tier >= `--min-tier` ± `--finding-pad`; the
//! `genome` scope counts them anywhere. `unmapped` is genome-wide in both scopes.
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
//!
//! `--clip-cap` (off by default) thins clip-only reads per reference bin after
//! the merge and before the payload pass; see `apply_clip_cap`.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use anyhow::{Context, Result, anyhow, bail};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};

use crate::bamio::{
    FIELDS_FISH_FULL, FIELDS_FISH_MIN, OP_D, OP_EQ, OP_H, OP_I, OP_M, OP_N, OP_S, OP_X, aux_int,
    op_parts, open_indexed, terminal_clips,
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
/// Reason name for a hairpin (fold-back) clip paired with an opposite-strand SA arm.
const R_HAIRPIN: &str = "hairpin";
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
    /// Linear contigs; all others are circular (a clip at a contig end — or
    /// within `circ_end_tol` bp with a same-contig, same-strand `SA` partner
    /// near the opposite end — is an origin crossing and does not count).
    pub linear: Vec<String>,
    /// Circular contigs: a contig-end clip whose alignment stops within this
    /// many bp of the end is an origin crossing (skipped) when the read's `SA`
    /// has a same-contig, same-strand segment within this many bp of the
    /// opposite end (`0` = exact ends only).
    pub circ_end_tol: u64,
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
    /// Do not emit `hairpin` events (and do not fall back to `clip` for them).
    pub no_hairpins: bool,
    pub no_indels: bool,
    pub no_linked: bool,
    /// Select any record overlapping a finding ± `locus_flank`.
    pub locus: bool,
    pub locus_flank: u64,
    /// Terminal clip length threshold.
    pub clip_min_len: u32,
    /// Hairpin: minimum reference overlap as a fraction of the shorter arm's reference span.
    pub hairpin_min_overlap: f64,
    /// Hairpin: minimum absolute tolerance (bp) for query gap/overlap between the two arms.
    pub hairpin_max_gap: u64,
    /// Hairpin: fractional tolerance (× read length) for query gap/overlap between the two arms.
    pub hairpin_max_gap_frac: f64,
    /// Hairpin: minimum query length of each arm as a fraction of read length.
    pub hairpin_min_arm_frac: f64,
    /// Minimum merged large-indel size.
    pub indel_min_len: u32,
    /// Merge same-type indels separated by at most this many ref bp.
    pub indel_merge_gap: u32,
    /// Only indels of at least this size take part in the merge.
    pub indel_piece_min: u32,
    /// Maximum clip-selected reads per `clip_bin`-bp reference bin; `0` (default)
    /// disables the cap. See `apply_clip_cap`.
    pub clip_cap: u32,
    /// Bin size (reference bp) for `clip_cap`; must be positive.
    pub clip_bin: u64,
    /// `--progress` reporter (off by default).
    pub progress: crate::progress::Progress,
}

impl Default for FishParams {
    fn default() -> Self {
        FishParams {
            threads: 1,
            chunk_size: None,
            linear: Vec::new(),
            circ_end_tol: 100,
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
            no_hairpins: false,
            no_indels: false,
            no_linked: false,
            locus: false,
            locus_flank: 5000,
            clip_min_len: 50,
            hairpin_min_overlap: 0.5,
            hairpin_max_gap: 200,
            hairpin_max_gap_frac: 0.02,
            hairpin_min_arm_frac: 0.2,
            indel_min_len: 50,
            indel_merge_gap: 50,
            indel_piece_min: 10,
            clip_cap: 0,
            clip_bin: 2000,
            progress: Default::default(),
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
    /// Per reason (`unmapped`, `clip`, `hairpin`, `large_ins`, `large_del`,
    /// `linked`, `locus`, `bait`); a read can appear under more than one.
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
    pub hairpins: bool,
    pub indels: bool,
    pub linked: bool,
    /// `bait` is on when `--bait` was given.
    pub bait: bool,
    pub locus: bool,
    pub locus_flank: u64,
    pub clip_min_len: u32,
    /// `--circ-end-tol` (bp; 0 = exact ends only).
    pub circ_end_tol: u64,
    pub hairpin_min_overlap: f64,
    pub hairpin_max_gap: u64,
    pub hairpin_max_gap_frac: f64,
    pub hairpin_min_arm_frac: f64,
    pub indel_min_len: u32,
    pub indel_merge_gap: u32,
    pub indel_piece_min: u32,
    /// `--clip-cap` (0 = off).
    pub clip_cap: u32,
    /// `--clip-bin` (bp).
    pub clip_bin: u64,
}

/// Reads dropped by the clip cap.
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct CapSummary {
    /// Clip-only reads dropped.
    pub reads: u64,
    /// Summed length of the dropped reads.
    pub bp: u64,
}

/// Hairpin-converted clip events detected (counted regardless of `--no-hairpins`).
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct HairpinDetected {
    /// Distinct reads with >= 1 hairpin-converted clip event.
    pub reads: u64,
    /// Counted alignment records with >= 1 hairpin-converted clip event.
    pub records: u64,
}

/// Origin-crossing clips skipped on circular contigs (unselected; 0 with `--no-clips`).
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct OriginCrossingSkipped {
    /// Distinct reads with >= 1 skipped origin-crossing clip.
    pub reads: u64,
    /// Counted alignment records with >= 1 skipped origin-crossing clip.
    pub records: u64,
    /// Skipped clip events (`clips_exact + clips_near_end`; a record skipped at
    /// both ends counts 2 clips, 1 record).
    pub clips: u64,
    /// Skipped clips at an exact contig end (`pos == 0` or `end == contig_len`).
    pub clips_exact: u64,
    /// Skipped clips within `--circ-end-tol` with a matching opposite-end `SA` partner.
    pub clips_near_end: u64,
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
    /// Reads dropped by `--clip-cap` (zeros when off); they are not in `reads`.
    pub cap: CapSummary,
    /// Hairpin-converted clip events detected (independent of `--no-hairpins`).
    pub hairpin_detected: HairpinDetected,
    /// Origin-crossing clips skipped on circular contigs (not selected).
    pub origin_crossing_skipped: OriginCrossingSkipped,
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

/// Cap the clip-selected reads per reference bin (`--clip-cap`/`--clip-bin`).
///
/// Every read whose sole reason is `clip` is eligible, whatever its clip length
/// (every other reason keeps its read). An eligible read is binned by
/// `(contig, pos / bin)` of its best clip event (the `event_better` order:
/// largest clip first). In a bin with more than `cap` eligible reads, the `cap`
/// reads with the largest best-clip length are kept and the rest are removed
/// from `selected`; equal lengths are broken by read name ascending. A long clip
/// spans furthest into sequence the assembly does not hold, so the longest clips
/// are the reads worth keeping, and ranking by length first also keeps a
/// hotspot's mixed clip lengths from evicting its long clips. The choice depends
/// only on the read data and names, so it is the same for any thread count or
/// input order. A no-op when `cap` or `bin` is 0. Returns the dropped reads and
/// their summed length.
fn apply_clip_cap(selected: &mut Vec<(String, ReadAgg)>, cap: u32, bin: u64) -> (u64, u64) {
    if cap == 0 || bin == 0 {
        return (0, 0);
    }
    // (bin key) -> eligible reads as (index, best clip length)
    let mut bins: HashMap<(&str, u64), Vec<(usize, u64)>> = HashMap::new();
    for (i, (_, agg)) in selected.iter().enumerate() {
        if !agg.events.keys().all(|(r, _)| *r == R_CLIP) {
            continue;
        }
        let Some(best) = agg.events.values().reduce(|a, b| if event_better(b, a) { b } else { a })
        else {
            continue;
        };
        bins.entry((best.contig.as_str(), best.pos / bin)).or_default().push((i, best.len));
    }
    let mut drop = vec![false; selected.len()];
    for v in bins.values_mut() {
        if v.len() <= cap as usize {
            continue;
        }
        v.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| selected[a.0].0.cmp(&selected[b.0].0)));
        for &(i, _) in &v[cap as usize..] {
            drop[i] = true;
        }
    }
    let mut reads = 0u64;
    let mut bp = 0u64;
    let mut i = 0;
    selected.retain(|(_, a)| {
        let d = drop[i];
        i += 1;
        if d {
            reads += 1;
            bp += a.read_len;
        }
        !d
    });
    (reads, bp)
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

/// Reference span (M/D/N/=/X) of an `SA` CIGAR string, or `None` if invalid.
fn sa_ref_span(cigar: &str) -> Option<u64> {
    let cs = bam::record::CigarString::try_from(cigar).ok()?;
    let mut span = 0u64;
    for &c in &cs.0 {
        let len = match c {
            bam::record::Cigar::Match(l)
            | bam::record::Cigar::Del(l)
            | bam::record::Cigar::RefSkip(l)
            | bam::record::Cigar::Equal(l)
            | bam::record::Cigar::Diff(l) => l as u64,
            _ => 0,
        };
        span = span.saturating_add(len);
    }
    Some(span)
}

/// `(lead_clip, trail_clip, ref_span, query_total)` of an `SA` CIGAR string,
/// or `None` if invalid or empty.
fn sa_cigar_stats(cigar: &str) -> Option<(u64, u64, u64, u64)> {
    let cs = bam::record::CigarString::try_from(cigar).ok()?;
    if cs.0.is_empty() {
        return None;
    }
    let is_clip = |c: &bam::record::Cigar| {
        matches!(c, bam::record::Cigar::SoftClip(_) | bam::record::Cigar::HardClip(_))
    };
    let op_len = |c: &bam::record::Cigar| -> u64 {
        match *c {
            bam::record::Cigar::Match(l)
            | bam::record::Cigar::Ins(l)
            | bam::record::Cigar::Del(l)
            | bam::record::Cigar::RefSkip(l)
            | bam::record::Cigar::SoftClip(l)
            | bam::record::Cigar::HardClip(l)
            | bam::record::Cigar::Equal(l)
            | bam::record::Cigar::Diff(l) => l as u64,
            bam::record::Cigar::Pad(_) => 0,
        }
    };
    let lead: u64 = cs.0.iter().take_while(|c| is_clip(c)).map(op_len).sum();
    let trail: u64 = cs.0.iter().rev().take_while(|c| is_clip(c)).map(op_len).sum();
    let mut ref_span = 0u64;
    let mut query_total = 0u64;
    for &c in &cs.0 {
        match c {
            bam::record::Cigar::Match(l)
            | bam::record::Cigar::Equal(l)
            | bam::record::Cigar::Diff(l) => {
                ref_span += l as u64;
                query_total += l as u64;
            }
            bam::record::Cigar::Del(l) | bam::record::Cigar::RefSkip(l) => {
                ref_span += l as u64;
            }
            bam::record::Cigar::Ins(l)
            | bam::record::Cigar::SoftClip(l)
            | bam::record::Cigar::HardClip(l) => {
                query_total += l as u64;
            }
            bam::record::Cigar::Pad(_) => {}
        }
    }
    Some((lead, trail, ref_span, query_total))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OriginClip {
    None,
    Exact,
    NearEnd,
}

/// Classify whether a terminal clip on a contig is an origin crossing that
/// should be skipped.
///
/// On linear contigs, returns `OriginClip::None`. On circular contigs:
/// - `left = true`: returns `OriginClip::Exact` when `pos == 0`, or
///   `OriginClip::NearEnd` when `pos <= tol` and `SA` has a same-contig,
///   same-strand entry whose reference end is `>= contig_len - tol`.
/// - `left = false`: returns `OriginClip::Exact` when `end == contig_len`, or
///   `OriginClip::NearEnd` when `end >= contig_len - tol` and `SA` has a
///   same-contig, same-strand entry whose 0-based start (`pos - 1`) is `<= tol`.
///
/// `tol` is capped at `contig_len / 2`. When it is 0, only exact contig ends
/// are skipped (and `sa_entries` is not invoked). `sa_entries` is called lazily
/// only for near-end candidates; invalid/absent `SA` payloads (`None`) are
/// treated as having no partner (the clip counts). The partner is taken from the
/// `SA` tag (not checked to be a counted alignment record) and is not checked
/// to cover the clipped read bases: any same-contig, same-strand segment near
/// the opposite end qualifies.
#[allow(clippy::too_many_arguments)]
fn classify_circular_origin_clip<'a>(
    circular: bool,
    left: bool,
    pos: u64,
    end: u64,
    contig: &str,
    contig_len: u64,
    strand: char,
    tol: u64,
    sa_entries: impl FnOnce() -> Option<&'a [crate::translate::SaEntry]>,
) -> OriginClip {
    if !circular {
        return OriginClip::None;
    }
    // On a contig shorter than 2 * tol every position would be "near an end";
    // cap the tolerance at half the contig so the two ends stay distinct.
    let tol = tol.min(contig_len / 2);
    if left {
        if pos == 0 {
            return OriginClip::Exact;
        }
        if tol == 0 || pos > tol {
            return OriginClip::None;
        }
        let threshold = contig_len.saturating_sub(tol);
        let Some(entries) = sa_entries() else { return OriginClip::None };
        if entries.iter().any(|e| {
            e.rname == contig
                && e.strand == strand
                && sa_ref_span(&e.cigar)
                    .map(|span| e.pos.saturating_sub(1).saturating_add(span))
                    .is_some_and(|re| re >= threshold)
        }) {
            OriginClip::NearEnd
        } else {
            OriginClip::None
        }
    } else {
        if end == contig_len {
            return OriginClip::Exact;
        }
        if tol == 0 || end < contig_len.saturating_sub(tol) {
            return OriginClip::None;
        }
        let Some(entries) = sa_entries() else { return OriginClip::None };
        if entries
            .iter()
            .any(|e| e.rname == contig && e.strand == strand && e.pos.saturating_sub(1) <= tol)
        {
            OriginClip::NearEnd
        } else {
            OriginClip::None
        }
    }
}

/// Aligned query interval `[start, end)` in 0-based half-open original read orientation.
fn orig_query_iv(strand: char, lead: u64, trail: u64, read_len: u64) -> Option<(u64, u64)> {
    if lead.saturating_add(trail) >= read_len {
        return None;
    }
    let (qs, qe) = match strand {
        '+' => (lead, read_len - trail),
        '-' => (trail, read_len - lead),
        _ => return None,
    };
    (qs < qe).then_some((qs, qe))
}

/// Check if a qualifying terminal clip on record `R` has an `SA` hairpin partner `E`.
///
/// Returns the formatted `detail` string (`arm1=q<s>-<e>,r<s>-<e>,<strand>;arm2=...`,
/// `arm1` earlier in read, 0-based half-open original-orientation query and reference
/// coordinates) if a matching partner covers the clipped region except a residual
/// `< clip_min_len`. Any pair of counted records qualifies (primary/supplementary or
/// supplementary/supplementary); the partner is taken from the `SA` tag (not checked
/// to be a counted alignment record). Invalid `SA` payloads are treated as no partner.
#[allow(clippy::too_many_arguments)]
fn find_hairpin_partner(
    r_contig: &str,
    r_rstart: u64,
    r_rend: u64,
    r_strand: char,
    r_lead: u64,
    r_trail: u64,
    read_len: u64,
    clip_is_left: bool,
    params: &FishParams,
    sa_entries: Option<&[crate::translate::SaEntry]>,
) -> Option<String> {
    let entries = sa_entries?;
    let (r_qs, r_qe) = orig_query_iv(r_strand, r_lead, r_trail, read_len)?;
    let r_qlen = r_qe - r_qs;
    let r_rspan = r_rend.saturating_sub(r_rstart);
    let clip_min = params.clip_min_len as u64;
    if r_rspan == 0
        || r_qlen < clip_min
        || (r_qlen as f64) < params.hairpin_min_arm_frac * (read_len as f64)
    {
        return None;
    }

    // Read-3' clip (`[r_qe, L)`) wants partner after R in read order;
    // Read-5' clip (`[0, r_qs)`) wants partner before R in read order.
    let want_e_after = clip_is_left ^ (r_strand == '+');
    let clip_len = if want_e_after { read_len - r_qe } else { r_qs };
    if clip_len < clip_min {
        return None;
    }

    let gap_limit = (params.hairpin_max_gap as f64).max(params.hairpin_max_gap_frac * (read_len as f64));
    let min_arm_len = params.hairpin_min_arm_frac * (read_len as f64);

    #[derive(PartialEq, Eq, PartialOrd, Ord)]
    struct PartnerScore(u64, std::cmp::Reverse<u64>, std::cmp::Reverse<u64>);

    // Deterministic tie-break if multiple SA entries qualify: largest ref_overlap,
    // smallest query gap, smallest residual, then SA order.
    let mut best: Option<(PartnerScore, String)> = None;

    for e in entries {
        if e.rname != r_contig || e.strand == r_strand {
            continue;
        }
        let Some((e_lead, e_trail, e_rspan, e_qtotal)) = sa_cigar_stats(&e.cigar) else {
            continue;
        };
        if e_rspan == 0 || e_qtotal != read_len {
            continue;
        }
        let Some((e_qs, e_qe)) = orig_query_iv(e.strand, e_lead, e_trail, read_len) else {
            continue;
        };
        let e_qlen = e_qe - e_qs;
        if e_qlen < clip_min || (e_qlen as f64) < min_arm_len {
            continue;
        }

        let e_rstart = e.pos.saturating_sub(1);
        let e_rend = e_rstart.saturating_add(e_rspan);
        let ref_overlap = r_rend.min(e_rend).saturating_sub(r_rstart.max(e_rstart));
        let shorter_rspan = r_rspan.min(e_rspan);
        if (ref_overlap as f64) < params.hairpin_min_overlap * (shorter_rspan as f64) {
            continue;
        }

        let q_gap = r_qs.max(e_qs).abs_diff(r_qe.min(e_qe));
        if (q_gap as f64) > gap_limit {
            continue;
        }

        if want_e_after {
            if !(e_qs > r_qs && e_qe > r_qe) {
                continue;
            }
        } else if !(e_qs < r_qs && e_qe < r_qe) {
            continue;
        }

        let covered = if want_e_after {
            read_len.min(e_qe).saturating_sub(r_qe.max(e_qs))
        } else {
            r_qs.min(e_qe).saturating_sub(e_qs)
        };
        let residual = clip_len.saturating_sub(covered);
        if residual >= clip_min {
            continue;
        }

        let ((q1s, q1e, r1s, r1e, s1), (q2s, q2e, r2s, r2e, s2)) = if want_e_after {
            ((r_qs, r_qe, r_rstart, r_rend, r_strand), (e_qs, e_qe, e_rstart, e_rend, e.strand))
        } else {
            ((e_qs, e_qe, e_rstart, e_rend, e.strand), (r_qs, r_qe, r_rstart, r_rend, r_strand))
        };
        let score = PartnerScore(ref_overlap, std::cmp::Reverse(q_gap), std::cmp::Reverse(residual));
        if best.as_ref().is_none_or(|(b, _)| score > *b) {
            let detail = format!(
                "arm1=q{q1s}-{q1e},r{r1s}-{r1e},{s1};arm2=q{q2s}-{q2e},r{r2s}-{r2e},{s2}"
            );
            best = Some((score, detail));
        }
    }

    best.map(|(_, d)| d)
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

#[derive(Default)]
struct WorkerStats {
    counts: InputCounts,
    hairpin_records: u64,
    hairpin_reads: HashSet<String>,
    origin_skipped_records: u64,
    origin_skipped_reads: HashSet<String>,
    origin_skipped_exact: u64,
    origin_skipped_near_end: u64,
}

impl WorkerStats {
    fn merge(&mut self, other: WorkerStats) {
        self.counts.add(&other.counts);
        self.hairpin_records += other.hairpin_records;
        self.hairpin_reads.extend(other.hairpin_reads);
        self.origin_skipped_records += other.origin_skipped_records;
        self.origin_skipped_reads.extend(other.origin_skipped_reads);
        self.origin_skipped_exact += other.origin_skipped_exact;
        self.origin_skipped_near_end += other.origin_skipped_near_end;
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
    stats: &mut WorkerStats,
    uncounted_primary: &mut HashSet<String>,
) -> Result<()> {
    let flags = rec.flags();
    stats.counts.total += 1;
    if flags & 0x100 != 0 {
        stats.counts.secondary += 1;
        return Ok(());
    }
    if flags & 0x200 != 0 {
        stats.counts.qcfail += 1;
        if flags & 0x800 == 0 {
            uncounted_primary.insert(String::from_utf8_lossy(rec.qname()).into_owned());
        }
        return Ok(());
    }
    if flags & 0x400 != 0 {
        stats.counts.duplicate += 1;
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
        stats.counts.unmapped += 1;
        stats.counts.unmapped_bp += rec.seq_len() as u64;
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
        stats.counts.primary += 1;
        stats.counts.mapped_bp += rec.seq_len() as u64;
    } else {
        stats.counts.supplementary += 1;
    }

    let rpos = rec.pos();
    if rpos < 0 {
        return Ok(());
    }
    let pos = rpos as usize;
    let cig = rec.raw_cigar();
    let mut span = 0usize;
    let mut cigar_qlen = 0u64;
    for &c in cig {
        let (op, l) = op_parts(c);
        if matches!(op, OP_M | OP_D | OP_N | OP_EQ | OP_X) {
            span += l;
        }
        if matches!(op, OP_M | OP_I | OP_S | OP_H | OP_EQ | OP_X) {
            cigar_qlen += l as u64;
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

    // clip (and hairpin): a terminal soft clip >= clip_min_len; on a circular
    // contig a clip at an exact contig end, or within `circ_end_tol` bp with a
    // same-contig, same-strand SA partner within `circ_end_tol` bp of the
    // opposite end, is an origin crossing (skipped). Otherwise, if an SA partner
    // satisfies the hairpin rule and covers the clip (residual < clip_min_len),
    // the clip becomes a `hairpin` event (suppressed by `--no-hairpins`).
    if !ctx.params.no_clips {
        let tol = ctx.params.circ_end_tol;
        let check_left = lclip && lead >= ctx.params.clip_min_len;
        let check_right = rclip && trail >= ctx.params.clip_min_len;
        if check_left || check_right {
            let sa_cell = std::cell::OnceCell::new();
            let get_sa = || -> Option<&[crate::translate::SaEntry]> {
                sa_cell
                    .get_or_init(|| match rec.aux(b"SA") {
                        Ok(Aux::String(s)) => crate::translate::parse_sa(s).ok(),
                        _ => None,
                    })
                    .as_deref()
            };
            let mut skipped_origin = false;
            let mut had_hairpin = false;

            for (is_left, enabled, clip_len, ev_pos, side_str) in [
                (true, check_left, lead, pos as u64, "left"),
                (false, check_right, trail, (end - 1) as u64, "right"),
            ] {
                if !enabled {
                    continue;
                }
                match classify_circular_origin_clip(
                    circular,
                    is_left,
                    pos as u64,
                    end as u64,
                    &contig,
                    contig_len as u64,
                    strand,
                    tol,
                    get_sa,
                ) {
                    OriginClip::Exact => {
                        skipped_origin = true;
                        stats.origin_skipped_exact += 1;
                    }
                    OriginClip::NearEnd => {
                        skipped_origin = true;
                        stats.origin_skipped_near_end += 1;
                    }
                    OriginClip::None => {
                        if let Some(ids) = ctx.event_attr(tid, ev_pos) {
                            if let Some(detail) = find_hairpin_partner(
                                &contig,
                                pos as u64,
                                end as u64,
                                strand,
                                lead as u64,
                                trail as u64,
                                cigar_qlen,
                                is_left,
                                ctx.params,
                                get_sa(),
                            ) {
                                had_hairpin = true;
                                if !ctx.params.no_hairpins {
                                    let base = Event {
                                        len: clip_len as u64,
                                        pos: ev_pos,
                                        primary,
                                        strand,
                                        detail,
                                        finding: None,
                                        contig: contig.clone(),
                                        sa: sa.clone(),
                                    };
                                    offerings(R_HAIRPIN, base, Some(ids), &mut evs);
                                }
                            } else {
                                let base = Event {
                                    len: clip_len as u64,
                                    pos: ev_pos,
                                    primary,
                                    strand,
                                    detail: format!("side={side_str};len={clip_len}"),
                                    finding: None,
                                    contig: contig.clone(),
                                    sa: sa.clone(),
                                };
                                offerings(R_CLIP, base, Some(ids), &mut evs);
                            }
                        }
                    }
                }
            }

            if skipped_origin {
                stats.origin_skipped_records += 1;
                stats.origin_skipped_reads.insert(name.clone());
            }
            if had_hairpin {
                stats.hairpin_records += 1;
                stats.hairpin_reads.insert(name.clone());
            }
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
    let total: u64 = tasks
        .iter()
        .map(|t| match *t {
            Task::Region { start, end, .. } => (end - start) as u64,
            Task::Unplaced => 0,
        })
        .sum();
    let done_bp = crate::progress::ChunkFrac::new(&params.progress, total);
    let workers: Vec<Result<HashMap<String, Payload>>> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..nthreads)
            .map(|_| {
                let tasks = tasks;
                let next = &next;
                let names = names;
                let done_bp = &done_bp;
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
                                let mut n = 0u64;
                                while let Some(res) = reader.read(&mut rec) {
                                    res.with_context(|| format!("reading {}", input.display()))?;
                                    let rp = rec.pos();
                                    if rp < start as i64 || rp >= end as i64 {
                                        continue;
                                    }
                                    n += 1;
                                    note_payload(&mut out, selected, &rec);
                                }
                                done_bp.add((end - start) as u64, n);
                            }
                            Task::Unplaced => {
                                reader
                                    .fetch(bam::FetchDefinition::Unmapped)
                                    .context("fetching coordinate-less records")?;
                                let mut n = 0u64;
                                while let Some(res) = reader.read(&mut rec) {
                                    res.with_context(|| format!("reading {}", input.display()))?;
                                    n += 1;
                                    note_payload(&mut out, selected, &rec);
                                }
                                done_bp.add(0, n);
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

/// `--progress` phases of `fish`, in order (`bait_lift` only with `--bait`).
pub fn fish_phases(bait: bool) -> Vec<&'static str> {
    let mut v = vec!["pass"];
    if bait {
        v.push("bait_lift");
    }
    v.extend(["payload", "output"]);
    v
}

/// Run `cramhook fish`. Returns the summary (also written to `P.fish.json`).
pub fn run_fish(input: &Path, ref_path: &Path, out_prefix: &Path, params: &FishParams) -> Result<FishSummary> {
    params.progress.plan(&fish_phases(!params.bait.is_empty()));
    params.progress.finish(run_fish_inner(input, ref_path, out_prefix, params))
}

fn run_fish_inner(input: &Path, ref_path: &Path, out_prefix: &Path, params: &FishParams) -> Result<FishSummary> {
    let progress = &params.progress;
    let t0 = Instant::now();
    if params.clip_bin == 0 {
        bail!("--clip-bin must be positive");
    }
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
    let records_total = progress.record_total(input);
    progress.start_records("pass", records_total);
    let done_bp = crate::progress::ChunkFrac::new(progress, total);
    let workers: Vec<Result<(HashMap<String, ReadAgg>, WorkerStats, HashSet<String>)>> =
        std::thread::scope(|sc| {
            let handles: Vec<_> = (0..nthreads)
                .map(|_| {
                    let ctx = &ctx;
                    let tasks = &tasks;
                    let next = &next;
                    let names = &contig_name;
                    let done_bp = &done_bp;
                    sc.spawn(move || -> Result<(HashMap<String, ReadAgg>, WorkerStats, HashSet<String>)> {
                        let mut reader = open_indexed(input, ref_path, FIELDS_FISH_FULL)?;
                        let mut rec = bam::Record::new();
                        let mut aggs: HashMap<String, ReadAgg> = HashMap::new();
                        let mut stats = WorkerStats::default();
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
                                    let mut n = 0u64;
                                    while let Some(res) = reader.read(&mut rec) {
                                        res.with_context(|| format!("reading {}", input.display()))?;
                                        let rp = rec.pos();
                                        if rp < start as i64 || rp >= end as i64 {
                                            continue;
                                        }
                                        n += 1;
                                        process_record(
                                            &rec,
                                            tid,
                                            ctx,
                                            &mut aggs,
                                            &mut stats,
                                            &mut uncounted,
                                        )?;
                                    }
                                    done_bp.add((end - start) as u64, n);
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
                                    let mut n = 0u64;
                                    while let Some(res) = reader.read(&mut rec) {
                                        res.with_context(|| format!("reading {}", input.display()))?;
                                        n += 1;
                                        process_record(
                                            &rec,
                                            usize::MAX,
                                            ctx,
                                            &mut aggs,
                                            &mut stats,
                                            &mut uncounted,
                                        )?;
                                    }
                                    done_bp.add(0, n);
                                }
                            }
                        }
                        Ok((aggs, stats, uncounted))
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap_or_else(|_| Err(anyhow!("fish worker panicked")))).collect()
        });

    // Merge workers deterministically (sum counters, fold per-read events).
    let mut stats = WorkerStats::default();
    let mut global: HashMap<String, ReadAgg> = HashMap::new();
    let mut uncounted_primary: HashSet<String> = HashSet::new();
    for w in workers {
        let (aggs, ws, uncounted) = w?;
        stats.merge(ws);
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
    let counts = stats.counts;
    let hairpin_detected = HairpinDetected {
        reads: stats.hairpin_reads.len() as u64,
        records: stats.hairpin_records,
    };
    let origin_crossing_skipped = OriginCrossingSkipped {
        reads: stats.origin_skipped_reads.len() as u64,
        records: stats.origin_skipped_records,
        clips: stats.origin_skipped_exact + stats.origin_skipped_near_end,
        clips_exact: stats.origin_skipped_exact,
        clips_near_end: stats.origin_skipped_near_end,
    };

    // Bait (optional): lift the input onto assembly + bait and select every read
    // that lands on a bait contig. Independent of the categories above, which
    // were all evaluated on the original `--in` alignment.
    progress.end();
    let mut bait_summary = None;
    if !params.bait.is_empty() {
        progress.start("bait_lift");
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
        progress.end();
    }
    progress.start_records("payload", records_total);

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
    let (cap_reads, cap_bp) = apply_clip_cap(&mut selected, params.clip_cap, params.clip_bin);

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

    progress.end();
    progress.start("output");
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
        cap: CapSummary { reads: cap_reads, bp: cap_bp },
        hairpin_detected,
        origin_crossing_skipped,
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
            hairpins: !params.no_hairpins,
            indels: !params.no_indels,
            linked: linked_on,
            bait: !params.bait.is_empty(),
            locus: params.locus,
            locus_flank: params.locus_flank,
            clip_min_len: params.clip_min_len,
            circ_end_tol: params.circ_end_tol,
            hairpin_min_overlap: params.hairpin_min_overlap,
            hairpin_max_gap: params.hairpin_max_gap,
            hairpin_max_gap_frac: params.hairpin_max_gap_frac,
            hairpin_min_arm_frac: params.hairpin_min_arm_frac,
            indel_min_len: params.indel_min_len,
            indel_merge_gap: params.indel_merge_gap,
            indel_piece_min: params.indel_piece_min,
            clip_cap: params.clip_cap,
            clip_bin: params.clip_bin,
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
    progress.end();
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

#[cfg(test)]
mod cap_unit {
    use super::*;

    fn ev(len: u64, pos: u64, contig: &str) -> Event {
        Event {
            len,
            pos,
            primary: true,
            strand: '+',
            detail: String::new(),
            finding: None,
            contig: contig.to_string(),
            sa: String::new(),
        }
    }

    fn read(name: &str, evs: &[(&'static str, Option<u32>, Event)]) -> (String, ReadAgg) {
        let mut a = ReadAgg { read_len: 1000, ..Default::default() };
        for (r, f, e) in evs {
            let mut e = e.clone();
            e.finding = *f;
            a.offer(*r, e);
        }
        (name.to_string(), a)
    }

    fn clip(name: &str, len: u64, pos: u64, contig: &str) -> (String, ReadAgg) {
        read(name, &[(R_CLIP, None, ev(len, pos, contig))])
    }

    fn names(v: &[(String, ReadAgg)]) -> Vec<String> {
        v.iter().map(|(n, _)| n.clone()).collect()
    }

    fn sorted(mut v: Vec<(String, ReadAgg)>) -> Vec<(String, ReadAgg)> {
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }

    fn bin_of(n: usize, prefix: &str, pos: u64, contig: &str) -> Vec<(String, ReadAgg)> {
        (0..n).map(|i| clip(&format!("{prefix}{i}"), 100, pos + i as u64, contig)).collect()
    }

    #[test]
    fn noop_when_off() {
        let mk = || sorted(bin_of(10, "r", 0, "c"));
        let mut v = mk();
        assert_eq!(apply_clip_cap(&mut v, 0, 2000), (0, 0));
        assert_eq!(v.len(), 10);
        assert_eq!(apply_clip_cap(&mut v, 4, 0), (0, 0));
        assert_eq!(v.len(), 10);
    }

    #[test]
    fn fifth_read_dropped_at_cap_4() {
        let mut v = sorted(bin_of(5, "r", 10, "c"));
        let (n, bp) = apply_clip_cap(&mut v, 4, 2000);
        assert_eq!((n, bp), (1, 1000));
        assert_eq!(v.len(), 4);
        // a bin at exactly the cap is untouched
        let mut v = sorted(bin_of(4, "r", 10, "c"));
        assert_eq!(apply_clip_cap(&mut v, 4, 2000), (0, 0));
    }

    #[test]
    fn non_clip_reasons_never_dropped() {
        let mut v = bin_of(5, "c", 10, "c");
        let others = [R_UNMAPPED, R_HAIRPIN, R_INS, R_DEL, R_LINKED, R_LOCUS, R_BAIT];
        for (r, other) in others.iter().enumerate() {
            v.push(read(
                &format!("x{r}"),
                &[(R_CLIP, None, ev(100, 10, "c")), (other, None, ev(100, 10, "c"))],
            ));
        }
        let mut v = sorted(v);
        let (n, _) = apply_clip_cap(&mut v, 1, 2000);
        assert_eq!(n, 4);
        let kept = names(&v);
        for r in 0..others.len() {
            assert!(kept.contains(&format!("x{r}")));
        }
        assert_eq!(kept.len(), others.len() + 1);
    }

    #[test]
    fn bins_are_independent() {
        let mut v = bin_of(5, "a", 10, "c");
        v.extend(bin_of(5, "b", 2010, "c"));
        v.extend(bin_of(3, "d", 10, "other"));
        let mut v = sorted(v);
        let (n, _) = apply_clip_cap(&mut v, 4, 2000);
        assert_eq!(n, 2);
        let kept = names(&v);
        assert_eq!(kept.iter().filter(|k| k.starts_with('a')).count(), 4);
        assert_eq!(kept.iter().filter(|k| k.starts_with('b')).count(), 4);
        assert_eq!(kept.iter().filter(|k| k.starts_with('d')).count(), 3);
    }

    #[test]
    fn long_clip_is_eligible_and_kept_as_longest() {
        // A 500 bp clip is no longer exempt, but being the longest it is kept.
        let mut v: Vec<_> = (0..6).map(|i| clip(&format!("l{i}"), 500, 10, "c")).collect();
        v.extend((0..6).map(|i| clip(&format!("s{i}"), 499, 10, "c")));
        let mut v = sorted(v);
        let (n, _) = apply_clip_cap(&mut v, 3, 2000);
        assert_eq!(n, 9);
        let kept = names(&v);
        assert_eq!(kept, vec!["l0", "l1", "l2"]);
    }

    #[test]
    fn two_longest_kept_at_cap_2() {
        let mut v = vec![
            clip("a", 120, 10, "c"),
            clip("b", 900, 10, "c"),
            clip("c", 400, 10, "c"),
            clip("d", 900, 10, "c"),
            clip("e", 600, 10, "c"),
        ];
        let (n, bp) = apply_clip_cap(&mut v, 2, 2000);
        assert_eq!((n, bp), (3, 3000));
        // The two longest (b, d tie at 900, so name ascending); c and e follow.
        assert_eq!(names(&v), vec!["b", "d"]);
    }

    #[test]
    fn two_clip_events_use_the_largest() {
        // Findings scope: clip events in two findings; the 300 bp clip at 5000
        // decides the bin (bin 2), not the 80 bp clip at 10 (bin 0).
        let two = |name: &str| {
            read(
                name,
                &[
                    (R_CLIP, Some(1), ev(80, 10, "c")),
                    (R_CLIP, Some(2), ev(300, 5000, "c")),
                ],
            )
        };
        let mut v: Vec<_> = (0..3).map(|i| two(&format!("t{i}"))).collect();
        v.extend(bin_of(3, "z", 10, "c")); // bin 0, three reads
        let mut v = sorted(v);
        // cap 3: bin 0 holds 3 reads, bin 2 holds 3 reads: nothing dropped
        assert_eq!(apply_clip_cap(&mut v, 3, 2000), (0, 0));
        // cap 2: one dropped from each bin; the two-clip reads sit in bin 2
        let (n, _) = apply_clip_cap(&mut v, 2, 2000);
        assert_eq!(n, 2);
        let kept = names(&v);
        assert_eq!(kept.iter().filter(|k| k.starts_with('t')).count(), 2);
        assert_eq!(kept.iter().filter(|k| k.starts_with('z')).count(), 2);
    }

    #[test]
    fn kept_set_independent_of_input_order() {
        // Two bins of 40 reads each with a spread of clip lengths, so the ranking
        // mixes length ordering and name ties.
        let mk = || {
            let mut v: Vec<(String, ReadAgg)> = (0..40)
                .map(|i| clip(&format!("read{i}"), 100 + (i as u64 % 7) * 50, 10 + i as u64, "c"))
                .collect();
            v.extend(
                (0..40).map(|i| clip(&format!("q{i}"), 100 + (i as u64 % 5) * 90, 4010 + i as u64, "c")),
            );
            v
        };
        let mut base = mk();
        apply_clip_cap(&mut base, 7, 2000);
        let mut expect = names(&base);
        expect.sort();
        assert_eq!(expect.len(), 14);
        // Deterministic permutations (reverse, rotate, by length, by name desc).
        for k in 0..4usize {
            let mut v = mk();
            match k {
                0 => v.reverse(),
                1 => v.rotate_left(17),
                2 => v.sort_by(|a, b| {
                    (a.0.len(), &a.0).cmp(&(b.0.len(), &b.0))
                }),
                _ => v.sort_by(|a, b| b.0.cmp(&a.0)),
            }
            apply_clip_cap(&mut v, 7, 2000);
            let mut got = names(&v);
            got.sort();
            assert_eq!(got, expect, "permutation {k}");
        }
    }
}

#[cfg(test)]
mod circ_clip_unit {
    use super::*;

    const LEN: u64 = 2_000;
    const TOL: u64 = 100;

    #[allow(clippy::too_many_arguments)]
    fn circ(
        circular: bool,
        left: bool,
        pos: u64,
        end: u64,
        contig: &str,
        contig_len: u64,
        strand: char,
        tol: u64,
        sa_raw: Option<&str>,
    ) -> OriginClip {
        let cell = std::cell::OnceCell::new();
        classify_circular_origin_clip(circular, left, pos, end, contig, contig_len, strand, tol, || {
            cell.get_or_init(|| sa_raw.and_then(|s| crate::translate::parse_sa(s).ok()))
                .as_deref()
        })
    }

    #[test]
    fn exact_ends_skipped_without_sa_even_at_tol_zero() {
        assert_eq!(circ(true, true, 0, 500, "chr", LEN, '+', TOL, None), OriginClip::Exact);
        assert_eq!(circ(true, true, 0, 500, "chr", LEN, '+', 0, None), OriginClip::Exact);
        assert_eq!(circ(true, false, 1500, LEN, "chr", LEN, '-', TOL, None), OriginClip::Exact);
        assert_eq!(circ(true, false, 1500, LEN, "chr", LEN, '-', 0, None), OriginClip::Exact);
    }

    #[test]
    fn linear_contig_never_skips() {
        assert_eq!(circ(false, true, 0, 500, "chr", LEN, '+', TOL, None), OriginClip::None);
        assert_eq!(
            circ(false, false, 1500, LEN, "chr", LEN, '+', TOL, Some("chr,3,+,120S400M,60,0;")),
            OriginClip::None
        );
    }

    #[test]
    fn near_end_with_matching_partner_skipped() {
        // Right clip stopping at 1990 with SA starting at 1-based 3 (0-based 2 <= 100).
        assert_eq!(
            circ(true, false, 1500, 1990, "chr", LEN, '+', TOL, Some("chr,3,+,120S400M,60,0;")),
            OriginClip::NearEnd
        );
        // Inclusive boundaries: end == LEN - TOL (1900) and SA start == TOL (1-based 101).
        assert_eq!(
            circ(true, false, 1500, 1900, "chr", LEN, '-', TOL, Some("chr,101,-,120S400M,60,0;")),
            OriginClip::NearEnd
        );
        // Left clip at pos 5 with SA ending at 1995 (>= 1900), using D/=/X ops.
        // 1-based 1801 -> 0-based 1800 + (100= + 50D + 45X) = 1995.
        assert_eq!(
            circ(true, true, 5, 505, "chr", LEN, '+', TOL, Some("chr,1801,+,100=50D45X120S,60,0;")),
            OriginClip::NearEnd
        );
        // Multiple SA entries: one non-matching, second matching.
        assert_eq!(
            circ(
                true,
                true,
                100,
                600,
                "chr",
                LEN,
                '-',
                TOL,
                Some("other,1,-,100M,60,0;chr,1701,-,200M,60,0;")
            ),
            OriginClip::NearEnd
        );
    }

    #[test]
    fn near_end_without_or_invalid_sa_not_skipped() {
        assert_eq!(circ(true, false, 1500, 1990, "chr", LEN, '+', TOL, None), OriginClip::None);
        assert_eq!(circ(true, true, 5, 505, "chr", LEN, '+', TOL, None), OriginClip::None);
        assert_eq!(circ(true, false, 1500, 1990, "chr", LEN, '+', TOL, Some("bad")), OriginClip::None);
        assert_eq!(circ(true, true, 5, 505, "chr", LEN, '+', TOL, Some("")), OriginClip::None);
    }

    #[test]
    fn partner_other_contig_not_skipped() {
        assert_eq!(
            circ(true, false, 1500, 1990, "chr", LEN, '+', TOL, Some("other,3,+,120S400M,60,0;")),
            OriginClip::None
        );
        assert_eq!(
            circ(true, true, 5, 505, "chr", LEN, '+', TOL, Some("other,1801,+,200M120S,60,0;")),
            OriginClip::None
        );
    }

    #[test]
    fn partner_wrong_strand_not_skipped() {
        assert_eq!(
            circ(true, false, 1500, 1990, "chr", LEN, '+', TOL, Some("chr,3,-,120S400M,60,0;")),
            OriginClip::None
        );
        assert_eq!(
            circ(true, true, 5, 505, "chr", LEN, '-', TOL, Some("chr,1801,+,200M120S,60,0;")),
            OriginClip::None
        );
    }

    #[test]
    fn too_far_record_or_partner_not_skipped() {
        // Partner too far from start (1-based 102 -> 0-based 101 > 100).
        assert_eq!(
            circ(true, false, 1500, 1990, "chr", LEN, '+', TOL, Some("chr,102,+,120S400M,60,0;")),
            OriginClip::None
        );
        // Record too far from right end (1899 < 1900).
        assert_eq!(
            circ(true, false, 1400, 1899, "chr", LEN, '+', TOL, Some("chr,3,+,120S400M,60,0;")),
            OriginClip::None
        );
        // Partner ends too far from right end (1700 + 199 = 1899 < 1900).
        assert_eq!(
            circ(true, true, 5, 505, "chr", LEN, '+', TOL, Some("chr,1701,+,199M120S,60,0;")),
            OriginClip::None
        );
        // Record too far from left end (101 > 100).
        assert_eq!(
            circ(true, true, 101, 601, "chr", LEN, '+', TOL, Some("chr,1801,+,200M120S,60,0;")),
            OriginClip::None
        );
    }

    #[test]
    fn tol_zero_rejects_near_end_even_with_partner() {
        assert_eq!(
            circ(true, false, 1500, 1990, "chr", LEN, '+', 0, Some("chr,1,+,120S400M,60,0;")),
            OriginClip::None
        );
        assert_eq!(
            circ(true, true, 1, 501, "chr", LEN, '+', 0, Some("chr,1801,+,200M120S,60,0;")),
            OriginClip::None
        );
    }

    #[test]
    fn tol_capped_at_half_contig_len() {
        // len 50, tol 100 -> effective 25: a right clip ending at 20 and a left
        // clip starting at 30 are not near an end.
        assert_eq!(
            circ(true, false, 0, 20, "chr", 50, '+', 100, Some("chr,1,+,10S30M,60,0;")),
            OriginClip::None
        );
        assert_eq!(
            circ(true, true, 30, 45, "chr", 50, '+', 100, Some("chr,21,+,30M10S,60,0;")),
            OriginClip::None
        );
    }

    #[test]
    fn saturates_when_tol_exceeds_contig_len() {
        assert_eq!(
            circ(true, false, 10, 40, "chr", 50, '+', 100, Some("chr,1,+,10S30M,60,0;")),
            OriginClip::NearEnd
        );
        assert_eq!(
            circ(true, true, 10, 40, "chr", 50, '+', 100, Some("chr,21,+,30M10S,60,0;")),
            OriginClip::NearEnd
        );
    }
}

#[cfg(test)]
mod hairpin_unit {
    use super::*;

    fn params_150() -> FishParams {
        FishParams { clip_min_len: 150, ..Default::default() }
    }

    #[allow(clippy::too_many_arguments)]
    fn check(
        contig: &str,
        rstart: u64,
        rend: u64,
        strand: char,
        lead: u64,
        trail: u64,
        read_len: u64,
        clip_is_left: bool,
        p: &FishParams,
        sa_raw: Option<&str>,
    ) -> Option<String> {
        let parsed = sa_raw.and_then(|s| crate::translate::parse_sa(s).ok());
        find_hairpin_partner(
            contig,
            rstart,
            rend,
            strand,
            lead,
            trail,
            read_len,
            clip_is_left,
            p,
            parsed.as_deref(),
        )
    }

    #[test]
    fn example_read_positive_both_directions() {
        // Example read: L=5234; arm1 + 3:3835-6489 (q7-2657, 7S..2577S); arm2 - 3:3865-6493 (q2657-5234, 2629M2657S with 52D or 2577M2657S).
        // Let's use 2629=2629 ref span (52D -> 2577M52D2657S = query 5234, ref 2629 -> 3864..6493).
        let p = params_150();
        let l = 5234u64;
        // R = primary (+), right clip (trail = 2577)
        let d_prim = check(
            "3",
            3834,
            6489,
            '+',
            7,
            2577,
            l,
            false,
            &p,
            Some("3,3865,-,2577M52D2657S,60,0;"),
        );
        assert_eq!(
            d_prim.as_deref(),
            Some("arm1=q7-2657,r3834-6489,+;arm2=q2657-5234,r3864-6493,-")
        );
        // R = primary (+) left clip (lead = 7, even if < clip_min_len, no partner before R)
        assert_eq!(
            check("3", 3834, 6489, '+', 7, 2577, l, true, &p, Some("3,3865,-,2577M52D2657S,60,0;")),
            None
        );

        // R = supplementary (-), right clip (trail = 2657 = read 5' [0, 2657))
        let d_supp = check(
            "3",
            3864,
            6493,
            '-',
            0,
            2657,
            l,
            false,
            &p,
            Some("3,3835,+,7S2650M5D2577S,60,0;"),
        );
        assert_eq!(
            d_supp.as_deref(),
            Some("arm1=q7-2657,r3834-6489,+;arm2=q2657-5234,r3864-6493,-")
        );
    }

    #[test]
    fn same_strand_rejected() {
        let p = params_150();
        assert_eq!(
            check("3", 3834, 6489, '+', 7, 2577, 5234, false, &p, Some("3,3865,+,2657S2577M,60,0;")),
            None
        );
    }

    #[test]
    fn different_contig_rejected() {
        let p = params_150();
        assert_eq!(
            check("3", 3834, 6489, '+', 7, 2577, 5234, false, &p, Some("4,3865,-,2577M2657S,60,0;")),
            None
        );
    }

    #[test]
    fn side_by_side_no_ref_overlap_rejected() {
        let p = params_150();
        // Adjacent on read (q0-2500 and q2500-5000), opposite strand, but disjoint on ref (r1000-3500 vs r3500-6000).
        assert_eq!(
            check("chr", 1000, 3500, '+', 0, 2500, 5000, false, &p, Some("chr,3501,-,2500M2500S,60,0;")),
            None
        );
        // Partial overlap < 0.5 of shorter span (2500 * 0.5 = 1250; overlap = 1249 -> r2251..4751).
        assert_eq!(
            check("chr", 1000, 3500, '+', 0, 2500, 5000, false, &p, Some("chr,2252,-,2500M2500S,60,0;")),
            None
        );
        // Overlap == 1250 -> accepted.
        assert!(
            check("chr", 1000, 3500, '+', 0, 2500, 5000, false, &p, Some("chr,2251,-,2500M2500S,60,0;"))
                .is_some()
        );
    }

    #[test]
    fn query_gap_or_overlap_too_large_rejected() {
        let p = params_150(); // L=5000 -> gap limit = max(200, 0.02 * 5000 = 100) = 200
        // Test positive gap with clip_min_len=250 so residual (< 250) doesn't reject gap=200.
        let p250 = FishParams { clip_min_len: 250, ..p.clone() };
        // R: [0, 2400) (trail=2600). E: [2601, 5000) -> gap = 201 > 200.
        assert_eq!(
            check("chr", 1000, 3400, '+', 0, 2600, 5000, false, &p250, Some("chr,1001,-,2399M2601S,60,0;")),
            None
        );
        // Gap == 200 -> E: [2600, 5000) -> accepted.
        assert!(
            check("chr", 1000, 3400, '+', 0, 2600, 5000, false, &p250, Some("chr,1001,-,2400M2600S,60,0;"))
                .is_some()
        );
        // Overlap == 201 -> R: [0, 2600) (trail=2400), E: [2399, 5000) (e_trail=2399, e_q=2601).
        assert_eq!(
            check("chr", 1000, 3600, '+', 0, 2400, 5000, false, &p, Some("chr,1001,-,2601M2399S,60,0;")),
            None
        );
        // Overlap == 200 -> accepted.
        assert!(
            check("chr", 1000, 3600, '+', 0, 2400, 5000, false, &p, Some("chr,1001,-,2600M2400S,60,0;"))
                .is_some()
        );
        // Long read L=20000 -> gap limit = 0.02 * 20000 = 400. Overlap 400 -> accepted; overlap 401 -> rejected.
        // R: [0, 10400) (trail=9600), E: [10000, 20000) -> overlap 400 -> accepted.
        assert!(
            check("chr", 1000, 11400, '+', 0, 9600, 20000, false, &p, Some("chr,1001,-,10000M10000S,60,0;"))
                .is_some()
        );
        // R: [0, 10401) (trail=9599), E: [10000, 20000) -> overlap 401 -> rejected.
        assert_eq!(
            check("chr", 1000, 11401, '+', 0, 9599, 20000, false, &p, Some("chr,1001,-,10000M10000S,60,0;")),
            None
        );
    }

    #[test]
    fn short_arm_rejected() {
        let p = params_150();
        // L=5000, min arm = 0.2 * 5000 = 1000.
        // Partner arm 999 bp ([4001, 5000)) -> rejected.
        assert_eq!(
            check("chr", 1000, 5001, '+', 0, 999, 5000, false, &p, Some("chr,1001,-,999M4001S,60,0;")),
            None
        );
        // R arm 999 bp ([0, 999)) -> rejected.
        assert_eq!(
            check("chr", 1000, 1999, '+', 0, 4001, 5000, false, &p, Some("chr,1001,-,4001M999S,60,0;")),
            None
        );
        // L=600, 0.2 * 600 = 120 < clip_min_len (150), so 140 bp arm >= 20% but < clip_min_len -> rejected.
        assert_eq!(
            check("chr", 1000, 1460, '+', 0, 140, 600, false, &p, Some("chr,1001,-,140M460S,60,0;")),
            None
        );
    }

    #[test]
    fn reverse_strand_primary_clip_side_mapping() {
        let p = params_150();
        let l = 5000u64;
        // Reverse primary: arm1 is `-` at [0, 2500) (so in `-` BAM CIGAR: 2500S2500M, lead=2500, trail=0).
        // Partner arm2 is `+` at [2500, 5000) (SA: 2500S2500M).
        // Left BAM clip (`clip_is_left = true`, lead = 2500 = read 3' [2500, 5000)) matches!
        assert_eq!(
            check("chr", 1000, 3500, '-', 2500, 0, l, true, &p, Some("chr,1001,+,2500S2500M,60,0;"))
                .as_deref(),
            Some("arm1=q0-2500,r1000-3500,-;arm2=q2500-5000,r1000-3500,+")
        );
        // Reverse primary with an extra 5' clip (trail = 300 in `-` BAM CIGAR: 2500S2200M300S -> [300, 2500)):
        // Left clip (lead = 2500 -> read 3' [2500, 5000)) converts to hairpin;
        // Right clip (trail = 300 -> read 5' [0, 300)) has no partner before R -> stays clip (`None`).
        assert!(
            check("chr", 1000, 3200, '-', 2500, 300, l, true, &p, Some("chr,1001,+,2500S2500M,60,0;"))
                .is_some()
        );
        assert_eq!(
            check("chr", 1000, 3200, '-', 2500, 300, l, false, &p, Some("chr,1001,+,2500S2500M,60,0;")),
            None
        );

        // Reverse primary as arm2 at [2500, 5000) with extra 3' clip (lead = 300 -> read 3' [4700, 5000)):
        // R CIGAR: 300S2200M2500S (r_qs=2500, r_qe=4700). Partner arm1 `+` at [0, 2500) (SA: 2500M2500S).
        // Right clip (trail = 2500 -> read 5' [0, 2500)) converts to hairpin;
        // Left clip (lead = 300 -> read 3' [4700, 5000)) stays clip (`None`).
        assert_eq!(
            check("chr", 1000, 3200, '-', 300, 2500, l, false, &p, Some("chr,1001,+,2500M2500S,60,0;"))
                .as_deref(),
            Some("arm1=q0-2500,r1000-3500,+;arm2=q2500-4700,r1000-3200,-")
        );
        assert_eq!(
            check("chr", 1000, 3200, '-', 300, 2500, l, true, &p, Some("chr,1001,+,2500M2500S,60,0;")),
            None
        );
    }

    #[test]
    fn residual_clip_ge_clip_min_len_stays_clip() {
        let p = params_150();
        let l = 5000u64;
        // R (+): [0, 2300) (trail = 2700, read 3' [2300, 5000)).
        // E (-): [2300, 4850) (e_lead=150, e_trail=2300, e_qlen=2550).
        // Covered = 2550, residual = 2700 - 2550 = 150 >= 150 -> stays clip (`None`)!
        assert_eq!(
            check("chr", 1000, 3300, '+', 0, 2700, l, false, &p, Some("chr,1001,-,150S2550M2300S,60,0;")),
            None
        );
        // E (-): [2300, 4851) -> residual = 149 < 150 -> hairpin!
        assert!(
            check("chr", 1000, 3300, '+', 0, 2700, l, false, &p, Some("chr,1001,-,149S2551M2300S,60,0;"))
                .is_some()
        );
    }

    #[test]
    fn sa_parse_failure_treated_as_no_partner() {
        let p = params_150();
        assert_eq!(check("3", 3834, 6489, '+', 7, 2577, 5234, false, &p, None), None);
        assert_eq!(check("3", 3834, 6489, '+', 7, 2577, 5234, false, &p, Some("")), None);
        assert_eq!(check("3", 3834, 6489, '+', 7, 2577, 5234, false, &p, Some("bad;")), None);
        assert_eq!(
            check("3", 3834, 6489, '+', 7, 2577, 5234, false, &p, Some("3,3865,-,2577M2657S,60,0")),
            None
        );
        // Mismatched query length in SA CIGAR -> treated as no partner.
        assert_eq!(
            check("3", 3834, 6489, '+', 7, 2577, 5234, false, &p, Some("3,3865,-,2577M2000S,60,0;")),
            None
        );
    }

    #[test]
    fn sa_hard_clipped_cigar_still_matches() {
        let p = params_150();
        let d = check(
            "3",
            3834,
            6489,
            '+',
            7,
            2577,
            5234,
            false,
            &p,
            Some("3,3865,-,2577M52D2657H,60,0;"),
        );
        assert_eq!(
            d.as_deref(),
            Some("arm1=q7-2657,r3834-6489,+;arm2=q2657-5234,r3864-6493,-")
        );
    }

    #[test]
    fn two_sa_candidates_tie_break() {
        let p = params_150();
        let l = 5000u64;
        // R (+): [0, 2500) at ref 1000..3500.
        // Candidate A: ref 1500..4000 (overlap 2000), q2500-5000.
        // Candidate B: ref 1000..3500 (overlap 2500), q2500-5000.
        // B must win regardless of order in SA tag (by larger ref_overlap).
        let a = "chr,1501,-,2500M2500S,60,0;";
        let b = "chr,1001,-,2500M2500S,60,0;";
        let want_b = Some("arm1=q0-2500,r1000-3500,+;arm2=q2500-5000,r1000-3500,-");
        assert_eq!(check("chr", 1000, 3500, '+', 0, 2500, l, false, &p, Some(&format!("{a}{b}"))).as_deref(), want_b);
        assert_eq!(check("chr", 1000, 3500, '+', 0, 2500, l, false, &p, Some(&format!("{b}{a}"))).as_deref(), want_b);

        // Equal ref_overlap (2500), different q_gap: C has q_gap=50 (q2550-5000), D has q_gap=0 (q2500-5000). D wins.
        let c = "chr,1001,-,2450M2550S,60,0;";
        let d = "chr,1001,-,2500M2500S,60,0;";
        assert_eq!(check("chr", 1000, 3500, '+', 0, 2500, l, false, &p, Some(&format!("{c}{d}"))).as_deref(), want_b);
        assert_eq!(check("chr", 1000, 3500, '+', 0, 2500, l, false, &p, Some(&format!("{d}{c}"))).as_deref(), want_b);

        // Equal ref_overlap (2400) and equal q_gap (0):
        // E covers [2500, 4900) (residual=100), F covers [2500, 4980) (residual=20). F wins.
        let e = "chr,1001,-,100S2400M2500S,60,0;";
        let f = "chr,1001,-,20S2480M2500S,60,0;";
        let want_f = Some("arm1=q0-2500,r1000-3500,+;arm2=q2500-4980,r1000-3480,-");
        assert_eq!(check("chr", 1000, 3500, '+', 0, 2500, l, false, &p, Some(&format!("{e}{f}"))).as_deref(), want_f);
        assert_eq!(check("chr", 1000, 3500, '+', 0, 2500, l, false, &p, Some(&format!("{f}{e}"))).as_deref(), want_f);
    }

    #[test]
    fn both_clips_are_hairpins() {
        let p = params_150();
        let l = 6000u64;
        // R (+) is the middle arm [2000, 4000) at r1000-3000 (lead=2000, trail=2000).
        // Partner before: (-) [0, 2000) at r1000-3000 -> SA: 4000S2000M.
        // Partner after:  (-) [4000, 6000) at r1000-3000 -> SA: 2000M4000S.
        let sa = "chr,1001,-,4000S2000M,60,0;chr,1001,-,2000M4000S,60,0;";
        assert_eq!(
            check("chr", 1000, 3000, '+', 2000, 2000, l, true, &p, Some(sa)).as_deref(),
            Some("arm1=q0-2000,r1000-3000,-;arm2=q2000-4000,r1000-3000,+")
        );
        assert_eq!(
            check("chr", 1000, 3000, '+', 2000, 2000, l, false, &p, Some(sa)).as_deref(),
            Some("arm1=q2000-4000,r1000-3000,+;arm2=q4000-6000,r1000-3000,-")
        );
    }
}
