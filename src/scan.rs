//! `scan`: one pass over a coordinate-sorted, indexed BAM/CRAM.
//!
//! Each record's CIGAR is walked against the reference directly (no pileup
//! API, no realignment). The result is a set of genome-wide per-position
//! counters (bases per strand, depth, deletions, clip starts, supplementary /
//! MAPQ0 / high-MAPQ depth), sparse left-normalised indel allele counts and
//! per-read large-indel events. `qc` turns these into signals.
//!
//! Parallelism: every contig is split into chunks; worker threads (std scoped
//! threads, one `IndexedReader` each) take chunks from a shared counter. A
//! record belongs to the chunk containing its start position, so nothing is
//! counted twice. Chunk arrays cover chunk start .. max record end and are added
//! into the contig arrays; all counters are integers, so the result does not
//! depend on the number of threads.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::path::Path;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{Context, Result, anyhow, bail};
use rust_htslib::bam::{self, Read};

use crate::bamio::{
    FIELDS_SCAN, OP_D, OP_EQ, OP_H, OP_I, OP_M, OP_N, OP_S, OP_X, op_parts, open_indexed,
    terminal_clips,
};
use crate::fasta::Fasta;

/// Coverage column: reads with a base or deletion, forward strand.
pub const DPF: usize = 0;
/// Coverage column: reads with a base or deletion, reverse strand.
pub const DPR: usize = 1;
/// Coverage column: deletions covering the position, forward strand.
pub const DELF: usize = 2;
/// Coverage column: deletions covering the position, reverse strand.
pub const DELR: usize = 3;
/// Coverage column: supplementary records.
pub const SUPP: usize = 4;
/// Coverage column: MAPQ 0 records.
pub const MQ0: usize = 5;
/// Coverage column: MAPQ >= [`HIGH_MAPQ`] records.
pub const HIQ: usize = 6;
/// Number of coverage columns.
pub const NCOV: usize = 7;
/// Minimum MAPQ counted as "high".
pub const HIGH_MAPQ: u8 = 20;

/// 4-bit BAM base code -> 0..4 for A/C/G/T, 4 for anything else.
const NT16_IDX: [u8; 16] = [4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4];
/// 4-bit BAM base code -> ASCII.
const NT16_ASCII: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Scan parameters.
#[derive(Clone, Debug)]
pub struct ScanParams {
    /// Worker threads.
    pub threads: usize,
    /// Merge same-type indels of one read separated by at most this many ref bp.
    pub large_merge_gap: u32,
    /// Keep merged events of at least this size.
    pub large_min: u32,
    /// Only indels of at least this size take part in the large-indel merge.
    pub large_piece_min: u32,
    /// Only clips of at least this length count towards the hotspot counters.
    pub clip_min_len: u32,
    /// Chunk size override (bp); `None` picks one from the genome size and threads.
    pub chunk_size: Option<u64>,
    /// Progress reporter for the scan phase (finished chunks).
    pub progress: crate::progress::Progress,
}

/// Record counts by category (records whose start lies in a chunk).
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct RecordCounts {
    /// All records seen.
    pub total: u64,
    /// Primary mapped records (not secondary, supplementary, duplicate or QC-fail).
    pub primary_mapped: u64,
    pub supplementary: u64,
    pub secondary: u64,
    pub duplicate: u64,
    pub qcfail: u64,
    /// Unmapped records (placed and unplaced).
    pub unmapped: u64,
    /// Read length (incl. hard clips) of primary mapped records.
    pub mapped_bp: u64,
    /// Read length of unmapped records.
    pub unmapped_bp: u64,
    /// Counted records with a hard clip.
    pub hard_clipped_records: u64,
    /// Sum of aligned reference spans (`M`/`D`/`N`/`=`/`X`) over primary mapped
    /// records (supplementary excluded); mean over `primary_mapped` gives the
    /// read reach used for the contig-end depth taper.
    pub aligned_span: u64,
}

impl RecordCounts {
    fn add(&mut self, o: &RecordCounts) {
        self.total += o.total;
        self.primary_mapped += o.primary_mapped;
        self.supplementary += o.supplementary;
        self.secondary += o.secondary;
        self.duplicate += o.duplicate;
        self.qcfail += o.qcfail;
        self.unmapped += o.unmapped;
        self.mapped_bp += o.mapped_bp;
        self.unmapped_bp += o.unmapped_bp;
        self.hard_clipped_records += o.hard_clipped_records;
        self.aligned_span += o.aligned_span;
    }
}

/// Final per-contig counters.
pub struct ContigCounts {
    pub name: String,
    pub len: usize,
    /// Base counts per position: forward A,C,G,T then reverse A,C,G,T.
    pub bases: Vec<[u32; 8]>,
    /// Coverage columns per position (see [`DPF`] .. [`HIQ`]).
    pub cov: Vec<[u32; NCOV]>,
    /// Clip starts per position: [left (first aligned base), right (last aligned base)].
    pub clips: Vec<[u32; 2]>,
    /// Clipped records whose clip lies at a contig end (left clip at 0, right
    /// clip at the contig end): origin crossings on circular contigs.
    pub edge_clips: u64,
    /// Sum of `M`/`=`/`X` (aligned query) bases over all counted records.
    pub aligned_bases: u64,
}

impl ContigCounts {
    /// Depth at a position: reads with a base or a deletion.
    #[inline]
    pub fn dp(&self, i: usize) -> u32 {
        self.cov[i][DPF] + self.cov[i][DPR]
    }

    /// Aligned-base depth at a position: reads with a base at the position
    /// (`M`/`=`/`X`), i.e. [`dp`](Self::dp) minus the reads whose only
    /// contribution is a deletion. This is the `q.depth.bedgraph` value and
    /// matches IGV's coverage track and the `samtools depth` default (any base
    /// counts, `N` included, deletions do not). qc's own DP — base or deletion
    /// — stays the signal input everywhere else.
    #[inline]
    pub fn aligned_dp(&self, i: usize) -> u32 {
        self.dp(i) - self.cov[i][DELF] - self.cov[i][DELR]
    }
}

/// One left-normalised indel allele with its strand counts.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IndelAllele {
    pub tid: u32,
    /// 0-based anchor (the reference base preceding the indel). For an
    /// `origin` allele this is 0 and the anchor base *follows* the event.
    pub anchor: u32,
    /// The indel sits before reference position 0 (insertion before the
    /// first base, or a deletion starting at 0); VCF uses the right-padded
    /// form (POS = 1, anchor base appended).
    pub origin: bool,
    pub ins: bool,
    pub len: u32,
    /// Inserted sequence (empty for deletions); always the exact sequence.
    pub seq: Vec<u8>,
    pub fwd: u32,
    pub rev: u32,
}

/// One per-read large-indel event (after the within-read merge).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct LargeEvent {
    pub tid: u32,
    pub ins: bool,
    /// 0-based start. DEL: first deleted base; INS: ref position before which
    /// the (first) insertion sits.
    pub start: u32,
    /// DEL: exclusive end of the merged deleted span. INS: position of the last
    /// merged insertion.
    pub end: u32,
    pub size: u32,
    pub reverse: bool,
    /// Hash of the read name (distinct-read support).
    pub qhash: u64,
}

/// Everything `scan` produces.
pub struct ScanResult {
    pub contigs: Vec<ContigCounts>,
    /// Sorted by (tid, anchor, origin, DEL before INS, len, seq).
    pub indels: Vec<IndelAllele>,
    /// Sorted.
    pub large: Vec<LargeEvent>,
    pub records: RecordCounts,
    pub chunks: usize,
    /// Mean aligned reference span (`M`/`D`/`N`/`=`/`X` bp) of a primary mapped
    /// record; 0 when there are none. The depth taper near linear contig ends
    /// extends this far (Flagger §2.5).
    pub mean_aligned_span: f64,
}

// ---------------------------------------------------------------------------
// small fast hasher (no extra crates)

/// Multiply-rotate hasher (FxHash style) for integer keys.
#[derive(Default, Clone, Copy)]
pub(crate) struct FxHasher(u64);

impl FxHasher {
    #[inline]
    fn add(&mut self, x: u64) {
        self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}

impl Hasher for FxHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.add(b as u64);
        }
    }
    #[inline]
    fn write_u8(&mut self, x: u8) {
        self.add(x as u64)
    }
    #[inline]
    fn write_u32(&mut self, x: u32) {
        self.add(x as u64)
    }
    #[inline]
    fn write_u64(&mut self, x: u64) {
        self.add(x)
    }
    #[inline]
    fn write_usize(&mut self, x: usize) {
        self.add(x as u64)
    }
}

pub(crate) type FxMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

/// FNV-1a over the read name.
fn qname_hash(q: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in q {
        h ^= b as u64;
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}

// ---------------------------------------------------------------------------
// pure helpers (unit-tested)

/// Left-normalise a deletion of `len` reference bases starting at 0-based `p`
/// (`p >= 1`). Returns the new start; the VCF anchor is `start - 1`.
pub(crate) fn normalize_del(refseq: &[u8], mut p: usize, len: usize) -> usize {
    while p >= 2 && refseq[p - 1] == refseq[p + len - 1] {
        p -= 1;
    }
    p
}

/// Left-normalise an insertion of `ins` placed before reference position `p`
/// (`p >= 1`). Rotates `ins` in place and returns the new position; the VCF
/// anchor is `p - 1`.
pub(crate) fn normalize_ins(refseq: &[u8], mut p: usize, ins: &mut [u8]) -> usize {
    if ins.is_empty() {
        return p;
    }
    while p >= 2 && refseq[p - 1] == ins[ins.len() - 1] {
        ins.rotate_right(1);
        p -= 1;
    }
    p
}

/// Homopolymer run lengths: `(run, fwd)` where `run[i]` is the length of the
/// run containing `i` and `fwd[i]` the run length starting at `i` (to the
/// right). Capped at 255.
pub fn hp_runs(seq: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let n = seq.len();
    let mut fwd = vec![0u8; n];
    let mut run = vec![0u8; n];
    let mut i = n;
    while i > 0 {
        i -= 1;
        fwd[i] = if i + 1 < n && seq[i + 1] == seq[i] { fwd[i + 1].saturating_add(1) } else { 1 };
    }
    let mut s = 0;
    while s < n {
        let mut e = s + 1;
        while e < n && seq[e] == seq[s] {
            e += 1;
        }
        let l = (e - s).min(255) as u8;
        run[s..e].fill(l);
        s = e;
    }
    (run, fwd)
}

/// Short-tandem-repeat span per position: the length (capped at 65535) of the
/// longest region with period 2..=4 and at least 3 copies covering the
/// position, homopolymers excluded; 0 outside any such region.
pub fn str_spans(seq: &[u8]) -> Vec<u16> {
    let n = seq.len();
    let mut out = vec![0u16; n];
    for k in 2..=4usize {
        if n < 3 * k {
            continue;
        }
        let mut a = 0usize;
        while a + k < n {
            if seq[a] != seq[a + k] {
                a += 1;
                continue;
            }
            let mut b = a;
            while b + k < n && seq[b] == seq[b + k] {
                b += 1;
            }
            // seq[i] == seq[i+k] for i in a..b: region a .. b+k has period k.
            let span = b + k - a;
            let mono = seq[a..a + k].iter().all(|&c| c == seq[a]);
            if span >= 3 * k && !mono {
                let v = span.min(u16::MAX as usize) as u16;
                for x in &mut out[a..b + k] {
                    *x = (*x).max(v);
                }
            }
            a = b + 1;
        }
    }
    out
}

/// One raw CIGAR indel of a read (reference coordinates, not normalised).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RawIndel {
    pub ins: bool,
    pub pos: u32,
    pub len: u32,
}

/// Within-read large-indel merge: consecutive same-type indels separated by at
/// most `gap` reference bp are merged (size = sum of lengths); merged events
/// of at least `min` bp are returned as `(ins, start, end, size)` (see
/// [`LargeEvent`] for `end`). `ev` is in CIGAR order.
pub(crate) fn merge_large(ev: &[RawIndel], gap: u32, min: u32) -> Vec<(bool, u32, u32, u32)> {
    let mut out = Vec::new();
    for ins in [false, true] {
        let mut cur: Option<(u32, u32, u32)> = None;
        for e in ev.iter().filter(|e| e.ins == ins) {
            let e_end = if ins { e.pos } else { e.pos + e.len };
            match &mut cur {
                Some((_, end, size)) if e.pos.saturating_sub(*end) <= gap => {
                    *end = e_end;
                    *size += e.len;
                }
                _ => {
                    if let Some((s, en, sz)) = cur.take() {
                        if sz >= min {
                            out.push((ins, s, en, sz));
                        }
                    }
                    cur = Some((e.pos, e_end, e.len));
                }
            }
        }
        if let Some((s, en, sz)) = cur {
            if sz >= min {
                out.push((ins, s, en, sz));
            }
        }
    }
    out
}

/// Pack an ACGT sequence of at most 32 bases into 2 bits per base.
fn pack_acgt(s: &[u8]) -> Option<u64> {
    if s.len() > 32 {
        return None;
    }
    let mut v = 0u64;
    for &b in s {
        let c = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        v = (v << 2) | c;
    }
    Some(v)
}

fn unpack_acgt(v: u64, len: u32) -> Vec<u8> {
    (0..len).rev().map(|i| b"ACGT"[((v >> (2 * i)) & 3) as usize]).collect()
}

// ---------------------------------------------------------------------------
// accumulators

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct IndelKey {
    tid: u32,
    anchor: u32,
    /// len << 2 | origin << 1 | ins
    len_kind: u32,
    /// 2-bit packed inserted sequence (<= 32 ACGT bases; others use `long_ins`).
    seq: u64,
}

/// One normalised indel of the current read (before the per-read merge).
#[derive(Clone, Copy, Debug)]
struct NormIndel {
    /// Normalised position: first deleted base / base the insertion precedes.
    pos: u32,
    origin: bool,
    ins: bool,
    len: u32,
    /// Offset of the inserted sequence in `Worker::seqs`.
    seq_off: u32,
}

/// Arrays for one chunk (or one whole contig), coverage as difference arrays.
struct Acc {
    start: usize,
    bases: Vec<[u32; 8]>,
    cov: Vec<[i32; NCOV]>,
    clips: Vec<[u32; 2]>,
}

impl Acc {
    fn new(start: usize, len: usize) -> Self {
        Acc {
            start,
            bases: vec![[0; 8]; len],
            cov: vec![[0; NCOV]; len + 1],
            clips: vec![[0; 2]; len],
        }
    }

    fn reset(&mut self, start: usize) {
        self.start = start;
        self.bases.clear();
        self.cov.clear();
        self.clips.clear();
    }

    /// Make sure positions up to `end` (exclusive) are addressable.
    #[inline]
    fn ensure(&mut self, end: usize) {
        let need = end - self.start;
        if self.bases.len() < need {
            self.bases.resize(need, [0; 8]);
            self.clips.resize(need, [0; 2]);
            self.cov.resize(need + 1, [0; NCOV]);
        }
    }

    /// Add `o` into `self` (o.start >= self.start).
    fn add(&mut self, o: &Acc) {
        let off = o.start - self.start;
        for (d, s) in self.bases[off..].iter_mut().zip(&o.bases) {
            for k in 0..8 {
                d[k] += s[k];
            }
        }
        for (d, s) in self.cov[off..].iter_mut().zip(&o.cov) {
            for k in 0..NCOV {
                d[k] += s[k];
            }
        }
        for (d, s) in self.clips[off..].iter_mut().zip(&o.clips) {
            d[0] += s[0];
            d[1] += s[1];
        }
    }
}

/// Per-worker state carried across chunks.
#[derive(Default)]
struct Worker {
    indels: FxMap<IndelKey, [u32; 2]>,
    /// Insertions with > 32 bases or non-ACGT bases, keyed by exact sequence.
    long_ins: FxMap<(IndelKey, Vec<u8>), [u32; 2]>,
    large: Vec<LargeEvent>,
    counts: RecordCounts,
    edge_clips: Vec<u64>,
    /// `M`/`=`/`X` bases per contig (see [`ContigCounts::aligned_bases`]).
    aligned_bases: Vec<u64>,
    /// Inserted sequences of the current read's `norm` events.
    seqs: Vec<u8>,
    norm: Vec<NormIndel>,
    /// Depth-contributing reference segments (M/=/X/D runs between `N`
    /// skips) of the current read, `[start, end)`.
    segs: Vec<(usize, usize)>,
    merged_seq: Vec<u8>,
    raw: Vec<RawIndel>,
}

impl Worker {
    fn add_indel(&mut self, tid: u32, anchor: u32, origin: bool, ins: bool, len: u32, seq: &[u8], rev: bool) {
        let len_kind = (len << 2) | ((origin as u32) << 1) | ins as u32;
        let packed = if ins { pack_acgt(seq) } else { Some(0) };
        match packed {
            Some(v) => {
                let key = IndelKey { tid, anchor, len_kind, seq: v };
                self.indels.entry(key).or_insert([0, 0])[rev as usize] += 1;
            }
            None => {
                let key = IndelKey { tid, anchor, len_kind, seq: 0 };
                self.long_ins.entry((key, seq.to_vec())).or_insert([0, 0])[rev as usize] += 1;
            }
        }
    }

    /// Count the current read's normalised indels: pieces that collapse onto
    /// the same normalised key (e.g. two 1-bp insertions in one homopolymer)
    /// are merged into one event (lengths added, inserted sequences
    /// concatenated in read order), so every read supports a key at most once.
    fn flush_read_indels(&mut self, tid: u32, rev: bool) {
        if self.norm.is_empty() {
            return;
        }
        let mut norm = std::mem::take(&mut self.norm);
        // Support must come from the same reads as DP at the anchor: keep an
        // indel only if its normalised anchor (position 0 for origin alleles)
        // lies in a depth-contributing segment of this read, i.e. not before
        // the alignment start and not inside an `N` skip. Indels that
        // left-normalise out of the read's covered span are dropped (the read
        // does not span the indel's repeat there).
        let segs = &self.segs;
        norm.retain(|e| {
            let a = if e.origin { 0 } else { e.pos as usize - 1 };
            segs.iter().any(|&(s, t)| a >= s && a < t)
        });
        norm.sort_by_key(|e| (e.pos, e.origin, e.ins));
        let mut i = 0;
        while i < norm.len() {
            let e0 = norm[i];
            let mut j = i;
            let mut len = 0u32;
            let mut seq = std::mem::take(&mut self.merged_seq);
            seq.clear();
            while j < norm.len() && (norm[j].pos, norm[j].origin, norm[j].ins) == (e0.pos, e0.origin, e0.ins) {
                len += norm[j].len;
                if e0.ins {
                    let o = norm[j].seq_off as usize;
                    seq.extend_from_slice(&self.seqs[o..o + norm[j].len as usize]);
                }
                j += 1;
            }
            let anchor = if e0.origin { 0 } else { e0.pos - 1 };
            self.add_indel(tid, anchor, e0.origin, e0.ins, len, &seq, rev);
            self.merged_seq = seq;
            i = j;
        }
        norm.clear();
        self.norm = norm;
        self.seqs.clear();
    }
}

/// Raw CIGAR `I`/`D` pieces of at least `piece_min` bp, in CIGAR order: the
/// input to [`merge_large`]. The same rule the `scan` hot path applies: a
/// deletion needs `l > 0`, an insertion is kept on its CIGAR length alone (a
/// read without a sequence still has an insertion of a known length). `out` is
/// cleared first.
pub(crate) fn raw_indels(cig: &[u32], pos: usize, piece_min: u32, out: &mut Vec<RawIndel>) {
    out.clear();
    let mut r = pos;
    for &c in cig {
        let (op, l) = op_parts(c);
        match op {
            OP_M | OP_EQ | OP_X => r += l,
            OP_D => {
                if l > 0 && l as u32 >= piece_min {
                    out.push(RawIndel { ins: false, pos: r as u32, len: l as u32 });
                }
                r += l;
            }
            OP_I => {
                if l as u32 >= piece_min {
                    out.push(RawIndel { ins: true, pos: r as u32, len: l as u32 });
                }
            }
            OP_N => r += l,
            _ => {}
        }
    }
}

/// Walk one mapped, counted record into `acc` / `w`.
fn process_record(
    rec: &bam::Record,
    tid: u32,
    refseq: &[u8],
    acc: &mut Acc,
    w: &mut Worker,
    p: &ScanParams,
) -> Result<()> {
    let flags = rec.flags();
    let pos = rec.pos() as usize;
    let cig = rec.raw_cigar();
    let enc = rec.seq().encoded;
    let rev = flags & 0x10 != 0;
    let s = rev as usize;
    let mapq = rec.mapq();
    let mut span = 0usize;
    for &c in cig {
        let (op, l) = op_parts(c);
        if matches!(op, OP_M | OP_D | OP_N | OP_EQ | OP_X) {
            span += l;
        }
    }
    let end = pos + span;
    if end > refseq.len() {
        bail!(
            "record {} ends at {} beyond contig length {}",
            String::from_utf8_lossy(rec.qname()),
            end,
            refseq.len()
        );
    }
    if span == 0 {
        return Ok(());
    }
    // Primary alignments only: the mean reference reach used for the taper.
    if flags & 0x800 == 0 {
        w.counts.aligned_span += span as u64;
    }
    let has_seq = enc.len() * 2 >= rec.seq_len() && rec.seq_len() > 0;
    acc.ensure(end);
    let off = acc.start;
    let mut cols = [0usize; 4];
    let mut ncols = 0;
    cols[ncols] = DPF + s;
    ncols += 1;
    if flags & 0x800 != 0 {
        cols[ncols] = SUPP;
        ncols += 1;
    }
    if mapq == 0 {
        cols[ncols] = MQ0;
        ncols += 1;
    }
    if mapq >= HIGH_MAPQ {
        cols[ncols] = HIQ;
        ncols += 1;
    }
    let cols = &cols[..ncols];
    let seg = |a: usize, b: usize, acc: &mut Acc| {
        if b > a {
            for &c in cols {
                acc.cov[a - off][c] += 1;
                acc.cov[b - off][c] -= 1;
            }
        }
    };

    // Clips. Sum all terminal S/H ops at each end (`5H60S100M`, `100M60S5H`):
    // an end counts once if its total clip is >= `clip_min_len`. A clip at a
    // contig end (left at 0, right at the last base) is an origin crossing /
    // read overhang: excluded from the counters, still counted in `edge_clips`
    // as before.
    let (lead_clip, trail_clip) = terminal_clips(cig);
    let lclip = lead_clip > 0;
    let rclip = cig.len() > 1 && trail_clip > 0;
    if (lclip && pos == 0) || (rclip && end == refseq.len()) {
        w.edge_clips[tid as usize] += 1;
    }
    if lclip && lead_clip >= p.clip_min_len && pos != 0 {
        acc.clips[pos - off][0] += 1;
    }
    if rclip && trail_clip >= p.clip_min_len && end != refseq.len() {
        acc.clips[end - 1 - off][1] += 1;
    }

    w.raw.clear();
    w.segs.clear();
    let mut r = pos;
    let mut q = 0usize;
    let mut seg_start = pos;
    for &c in cig {
        let (op, l) = op_parts(c);
        match op {
            OP_M | OP_EQ | OP_X => {
                w.aligned_bases[tid as usize] += l as u64;
                if has_seq {
                    let bases = &mut acc.bases[r - off..r - off + l];
                    let so = s * 4;
                    for (k, slot) in bases.iter_mut().enumerate() {
                        let qq = q + k;
                        let nib = (enc[qq >> 1] >> ((!qq & 1) << 2)) & 0xf;
                        let b = NT16_IDX[nib as usize] as usize;
                        if b < 4 {
                            slot[so + b] += 1;
                        }
                    }
                }
                r += l;
                q += l;
            }
            OP_I => {
                if l > 0 && has_seq {
                    let off_s = w.seqs.len();
                    for qq in q..q + l {
                        let nib = (enc[qq >> 1] >> ((!qq & 1) << 2)) & 0xf;
                        w.seqs.push(NT16_ASCII[nib as usize]);
                    }
                    let ins = &mut w.seqs[off_s..];
                    let mut origin = false;
                    let mut np = r;
                    if r == 0 {
                        // Before the first base: shift right once if equivalent
                        // (then anchor 0 is an ordinary left anchor), else keep
                        // it as an origin allele.
                        if refseq[0] == ins[0] {
                            ins.rotate_left(1);
                            np = 1;
                        } else {
                            origin = true;
                        }
                    } else {
                        np = normalize_ins(refseq, r, ins);
                    }
                    // Anchor coverage is checked per read in `flush_read_indels`.
                    w.norm.push(NormIndel { pos: np as u32, origin, ins: true, len: l as u32, seq_off: off_s as u32 });
                }
                q += l;
            }
            OP_D => {
                if l > 0 {
                    acc.cov[r - off][DELF + s] += 1;
                    acc.cov[r + l - off][DELF + s] -= 1;
                    let (np, origin) = if r == 0 {
                        // Deletion of the first bases: shift right once if
                        // equivalent, else an origin allele.
                        // (A deletion of the whole contig has no anchor and is dropped.)
                        if l >= refseq.len() {
                            (usize::MAX, false)
                        } else if refseq[0] == refseq[l] {
                            (1, false)
                        } else {
                            (0, true)
                        }
                    } else {
                        (normalize_del(refseq, r, l), false)
                    };
                    // Anchor coverage is checked per read in `flush_read_indels`.
                    if np != usize::MAX {
                        w.norm.push(NormIndel { pos: np as u32, origin, ins: false, len: l as u32, seq_off: 0 });
                    }
                }
                r += l;
            }
            OP_N => {
                seg(seg_start, r, acc);
                w.segs.push((seg_start, r));
                r += l;
                seg_start = r;
            }
            OP_S => q += l,
            _ => {}
        }
    }
    seg(seg_start, r, acc);
    w.segs.push((seg_start, r));
    w.flush_read_indels(tid, rev);

    raw_indels(cig, pos, p.large_piece_min, &mut w.raw);
    if !w.raw.is_empty() {
        let qh = qname_hash(rec.qname());
        for (ins, start, end, size) in merge_large(&w.raw, p.large_merge_gap, p.large_min) {
            w.large.push(LargeEvent { tid, ins, start, end, size, reverse: rev, qhash: qh });
        }
    }
    Ok(())
}

/// Read length including hard clips.
fn read_len(rec: &bam::Record) -> u64 {
    let mut h = 0u64;
    for &c in rec.raw_cigar() {
        let (op, l) = op_parts(c);
        if op == OP_H {
            h += l as u64;
        }
    }
    rec.seq_len() as u64 + h
}

/// Count a record; returns true if it should be walked.
fn classify(rec: &bam::Record, c: &mut RecordCounts) -> bool {
    let f = rec.flags();
    c.total += 1;
    if f & 0x100 != 0 {
        c.secondary += 1;
        return false;
    }
    if f & 0x200 != 0 {
        c.qcfail += 1;
        return false;
    }
    if f & 0x400 != 0 {
        c.duplicate += 1;
        return false;
    }
    if f & 0x4 != 0 {
        c.unmapped += 1;
        c.unmapped_bp += rec.seq_len() as u64;
        return false;
    }
    if rec.raw_cigar().iter().any(|&x| op_parts(x).0 == OP_H) {
        c.hard_clipped_records += 1;
    }
    if f & 0x800 != 0 {
        c.supplementary += 1;
    } else {
        c.primary_mapped += 1;
        c.mapped_bp += read_len(rec);
    }
    true
}

/// Header text of an alignment file (for the preset fallback).
pub fn read_header_text(input: &Path, ref_path: &Path) -> Result<String> {
    let r = open_indexed(input, ref_path, FIELDS_SCAN)?;
    Ok(String::from_utf8_lossy(r.header().as_bytes()).into_owned())
}

/// Pick a chunk size: about four chunks per thread, between 200 kb and 1 Mb.
fn chunk_size(total: u64, threads: usize) -> u64 {
    (total / (4 * threads.max(1) as u64)).clamp(200_000, 1_000_000)
}

enum Task {
    Region { tid: usize, start: usize, end: usize },
    Unplaced,
}

/// Run the scan.
pub fn scan(input: &Path, ref_path: &Path, fasta: &Fasta, p: &ScanParams) -> Result<ScanResult> {
    let r0 = open_indexed(input, ref_path, FIELDS_SCAN)?;
    let hv = r0.header().clone();
    drop(r0);
    let names: Vec<String> =
        hv.target_names().iter().map(|n| String::from_utf8_lossy(n).into_owned()).collect();
    let mut refs: Vec<&[u8]> = Vec::with_capacity(names.len());
    for (tid, name) in names.iter().enumerate() {
        let seq = fasta
            .seq(name)
            .ok_or_else(|| anyhow!("contig {name:?} of {} is not in the reference", input.display()))?;
        let hl = hv.target_len(tid as u32).unwrap_or(0);
        if hl != seq.len() as u64 {
            bail!("contig {name:?}: header length {hl} != reference length {}", seq.len());
        }
        refs.push(seq);
    }
    let total: u64 = refs.iter().map(|s| s.len() as u64).sum();
    let cs = p.chunk_size.unwrap_or_else(|| chunk_size(total, p.threads)).max(1) as usize;
    let mut tasks = Vec::new();
    for (tid, s) in refs.iter().enumerate() {
        let mut a = 0;
        while a < s.len() {
            let b = (a + cs).min(s.len());
            tasks.push(Task::Region { tid, start: a, end: b });
            a = b;
        }
    }
    tasks.push(Task::Unplaced);
    let n_chunks = tasks.len() - 1;

    let global: Vec<Mutex<Option<Acc>>> = refs.iter().map(|_| Mutex::new(None)).collect();
    let next = AtomicUsize::new(0);
    let nthreads = p.threads.max(1).min(tasks.len());
    let done_bp = crate::progress::ChunkFrac::new(&p.progress, total);

    let workers: Vec<Result<Worker>> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..nthreads)
            .map(|_| {
                sc.spawn(|| -> Result<Worker> {
                    let mut reader = open_indexed(input, ref_path, FIELDS_SCAN)?;
                    let mut w = Worker { edge_clips: vec![0; refs.len()], aligned_bases: vec![0; refs.len()], ..Default::default() };
                    let mut rec = bam::Record::new();
                    let mut acc = Acc::new(0, 0);
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(task) = tasks.get(i) else { break };
                        match *task {
                            Task::Region { tid, start, end } => {
                                acc.reset(start);
                                reader
                                    .fetch((tid as i32, start as i64, end as i64))
                                    .with_context(|| format!("fetching {}:{start}-{end}", names[tid]))?;
                                while let Some(res) = reader.read(&mut rec) {
                                    res.with_context(|| format!("reading {}", input.display()))?;
                                    let rp = rec.pos();
                                    if rp < start as i64 || rp >= end as i64 {
                                        continue;
                                    }
                                    if !classify(&rec, &mut w.counts) {
                                        continue;
                                    }
                                    process_record(&rec, tid as u32, refs[tid], &mut acc, &mut w, p)?;
                                }
                                if !acc.bases.is_empty() {
                                    let mut g = global[tid].lock().unwrap();
                                    let ga = g.get_or_insert_with(|| Acc::new(0, refs[tid].len()));
                                    ga.add(&acc);
                                }
                                done_bp.add((end - start) as u64);
                            }
                            Task::Unplaced => {
                                reader
                                    .fetch(bam::FetchDefinition::Unmapped)
                                    .context("fetching unplaced unmapped reads")?;
                                while let Some(res) = reader.read(&mut rec) {
                                    res.with_context(|| format!("reading {}", input.display()))?;
                                    if rec.tid() >= 0 {
                                        continue;
                                    }
                                    classify(&rec, &mut w.counts);
                                }
                            }
                        }
                    }
                    Ok(w)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap_or_else(|_| Err(anyhow!("scan worker panicked")))).collect()
    });

    let mut indels: FxMap<IndelKey, [u32; 2]> = FxMap::default();
    let mut long_ins: FxMap<(IndelKey, Vec<u8>), [u32; 2]> = FxMap::default();
    let mut large = Vec::new();
    let mut records = RecordCounts::default();
    let mut edge = vec![0u64; refs.len()];
    let mut aligned = vec![0u64; refs.len()];
    for w in workers {
        let w = w?;
        for (k, v) in w.indels {
            let e = indels.entry(k).or_insert([0, 0]);
            e[0] += v[0];
            e[1] += v[1];
        }
        for (k, v) in w.long_ins {
            let e = long_ins.entry(k).or_insert([0, 0]);
            e[0] += v[0];
            e[1] += v[1];
        }
        large.extend(w.large);
        records.add(&w.counts);
        for (e, x) in edge.iter_mut().zip(&w.edge_clips) {
            *e += x;
        }
        for (a, x) in aligned.iter_mut().zip(&w.aligned_bases) {
            *a += x;
        }
    }

    let allele = |k: &IndelKey, seq: Vec<u8>, v: [u32; 2]| IndelAllele {
        tid: k.tid,
        anchor: k.anchor,
        origin: k.len_kind & 2 == 2,
        ins: k.len_kind & 1 == 1,
        len: k.len_kind >> 2,
        seq,
        fwd: v[0],
        rev: v[1],
    };
    let mut out_indels: Vec<IndelAllele> = indels
        .into_iter()
        .map(|(k, v)| {
            let ins = k.len_kind & 1 == 1;
            let seq = if ins { unpack_acgt(k.seq, k.len_kind >> 2) } else { Vec::new() };
            allele(&k, seq, v)
        })
        .chain(long_ins.into_iter().map(|((k, seq), v)| allele(&k, seq, v)))
        .collect();
    out_indels.sort_unstable_by(|a, b| {
        (a.tid, a.anchor, a.origin, a.ins, a.len, &a.seq).cmp(&(b.tid, b.anchor, b.origin, b.ins, b.len, &b.seq))
    });
    large.sort_unstable();

    let mut contigs = Vec::with_capacity(refs.len());
    for (tid, g) in global.into_iter().enumerate() {
        let len = refs[tid].len();
        let acc = g.into_inner().unwrap().unwrap_or_else(|| Acc::new(0, len));
        let Acc { bases, cov, clips, .. } = acc;
        let mut run = [0i32; NCOV];
        let cov: Vec<[u32; NCOV]> = cov
            .into_iter()
            .take(len)
            .map(|d| {
                let mut o = [0u32; NCOV];
                for k in 0..NCOV {
                    run[k] += d[k];
                    o[k] = run[k].max(0) as u32;
                }
                o
            })
            .collect();
        contigs.push(ContigCounts {
            name: names[tid].clone(),
            len,
            bases,
            cov,
            clips,
            edge_clips: edge[tid],
            aligned_bases: aligned[tid],
        });
    }
    let mean_aligned_span =
        if records.primary_mapped > 0 { records.aligned_span as f64 / records.primary_mapped as f64 } else { 0.0 };
    Ok(ScanResult { contigs, indels: out_indels, large, records, chunks: n_chunks, mean_aligned_span })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn del_left_normalisation_homopolymer() {
        // 0123456789
        // GCAAAAT...: deleting any A of the run moves to the first A (p=2).
        let r = b"GCAAAATG";
        assert_eq!(normalize_del(r, 5, 1), 2);
        assert_eq!(normalize_del(r, 3, 1), 2);
        assert_eq!(normalize_del(r, 2, 1), 2);
        // Deleting "AA" at the end of the run.
        assert_eq!(normalize_del(r, 4, 2), 2);
        // Non-repeat: stays.
        assert_eq!(normalize_del(r, 6, 1), 6);
    }

    #[test]
    fn del_left_normalisation_dinucleotide() {
        // GCATATATG: delete "AT" at 6 -> "AT" at 2 (anchor C at 1).
        let r = b"GCATATATG";
        assert_eq!(normalize_del(r, 6, 2), 2);
        assert_eq!(normalize_del(r, 4, 2), 2);
        // Deleting "TA" at 3 is the same event shifted: left-normalises to 2 as "AT".
        assert_eq!(normalize_del(r, 3, 2), 2);
        // Never shifts the anchor below 0.
        assert_eq!(normalize_del(b"AAAA", 3, 1), 1);
    }

    #[test]
    fn ins_left_normalisation() {
        // Insert A inside the run GCAAAAT: before pos 5 -> before 2 (anchor C).
        let r = b"GCAAAATG";
        let mut s = b"A".to_vec();
        assert_eq!(normalize_ins(r, 5, &mut s), 2);
        assert_eq!(s, b"A");
        // Dinucleotide: GCATATG, insert "AT" before 6 (after ...ATAT) -> before 2.
        let r = b"GCATATG";
        let mut s = b"AT".to_vec();
        assert_eq!(normalize_ins(r, 6, &mut s), 2);
        assert_eq!(s, b"AT");
        // Insert "TA" before 5 (between A and T at 4/5): G C A T A | T G
        let mut s = b"TA".to_vec();
        let p = normalize_ins(r, 5, &mut s);
        assert_eq!((p, s.as_slice()), (2, b"AT".as_slice()));
        // No repeat: stays.
        let mut s = b"G".to_vec();
        assert_eq!(normalize_ins(r, 3, &mut s), 3);
    }

    #[test]
    fn homopolymer_runs() {
        let (run, fwd) = hp_runs(b"ACCCGTTA");
        assert_eq!(run, vec![1, 3, 3, 3, 1, 2, 2, 1]);
        assert_eq!(fwd, vec![1, 3, 2, 1, 1, 2, 1, 1]);
    }

    #[test]
    fn str_span_detection() {
        //            0123456789
        let s = b"GCATATATGG";
        let v = str_spans(s);
        // ATATAT (2..8) is 3 copies of AT -> span 6.
        assert_eq!(&v[2..8], &[6, 6, 6, 6, 6, 6]);
        assert_eq!(v[0], 0);
        assert_eq!(v[9], 0);
        // Two copies only: not an STR.
        assert!(str_spans(b"GCATATGG").iter().all(|&x| x == 0));
        // Homopolymer is not an STR.
        assert!(str_spans(b"GAAAAAAAAC").iter().all(|&x| x == 0));
        // Trinucleotide x3.
        let v = str_spans(b"TCAGCAGCAGT");
        assert_eq!(v[1], 9);
        assert_eq!(v[9], 9);
    }

    #[test]
    fn large_merge_within_read() {
        let ev = [
            RawIndel { ins: false, pos: 100, len: 30 },
            RawIndel { ins: true, pos: 140, len: 20 },
            // gap from 130 to 160 = 30 <= 50: merged with the first DEL.
            RawIndel { ins: false, pos: 160, len: 25 },
            // gap from 185 to 300 > 50: separate, too small alone.
            RawIndel { ins: false, pos: 300, len: 40 },
            RawIndel { ins: true, pos: 180, len: 35 },
            RawIndel { ins: false, pos: 1000, len: 60 },
        ];
        let m = merge_large(&ev, 50, 50);
        assert_eq!(m, vec![(false, 100, 185, 55), (false, 1000, 1060, 60), (true, 140, 180, 55)]);
        // Gap just beyond the limit: no merge.
        let ev = [RawIndel { ins: false, pos: 0, len: 30 }, RawIndel { ins: false, pos: 81, len: 30 }];
        assert!(merge_large(&ev, 50, 50).is_empty());
        let ev = [RawIndel { ins: false, pos: 0, len: 30 }, RawIndel { ins: false, pos: 80, len: 30 }];
        assert_eq!(merge_large(&ev, 50, 50), vec![(false, 0, 110, 60)]);
    }

    #[test]
    fn pack_roundtrip() {
        let s = b"ACGTTGCA";
        assert_eq!(unpack_acgt(pack_acgt(s).unwrap(), 8), s);
        assert!(pack_acgt(b"ACN").is_none());
    }

    /// Walk one hand-built record through `process_record` and return the
    /// chunk accumulator and worker.
    fn walk(refseq: &[u8], pos: i64, ops: Vec<rust_htslib::bam::record::Cigar>, min_len: u32) -> (Acc, Worker) {
        let qlen: usize = ops
            .iter()
            .map(|c| match *c {
                rust_htslib::bam::record::Cigar::Match(l)
                | rust_htslib::bam::record::Cigar::SoftClip(l)
                | rust_htslib::bam::record::Cigar::Ins(l) => l as usize,
                _ => 0,
            })
            .sum();
        let mut rec = bam::Record::new();
        let seq = vec![b'A'; qlen];
        let qual = vec![30u8; qlen];
        rec.set(b"r", Some(&rust_htslib::bam::record::CigarString(ops)), &seq, &qual);
        rec.set_tid(0);
        rec.set_pos(pos);
        rec.set_mapq(60);
        rec.set_flags(0);
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        let p = ScanParams {
            threads: 1,
            large_merge_gap: 50,
            large_min: 50,
            large_piece_min: 10,
            clip_min_len: min_len,
            chunk_size: None,
            progress: Default::default(),
        };
        let mut acc = Acc::new(0, refseq.len());
        let mut w = Worker { edge_clips: vec![0; 1], aligned_bases: vec![0; 1], ..Default::default() };
        process_record(&rec, 0, refseq, &mut acc, &mut w, &p).unwrap();
        (acc, w)
    }

    #[test]
    fn clip_length_filter_and_edge_exclusion() {
        use rust_htslib::bam::record::Cigar;
        let refseq = b"ACGT".repeat(50); // 200 bp
        // A 60 bp left clip counts at clip_min_len 50 but not at 61.
        let ops = || vec![Cigar::SoftClip(60), Cigar::Match(100)];
        let (acc, _w) = walk(&refseq, 10, ops(), 50);
        assert_eq!(acc.clips[10][0], 1);
        let (acc, _w) = walk(&refseq, 10, ops(), 61);
        assert_eq!(acc.clips[10][0], 0);
        // A left clip at position 0 is an edge clip: excluded from the counters,
        // still counted in edge_clips.
        let (acc, w) = walk(&refseq, 0, ops(), 50);
        assert_eq!(acc.clips[0][0], 0);
        assert_eq!(w.edge_clips[0], 1);
        // A right clip ending at the last base is an edge clip too.
        let (acc, w) = walk(&refseq, 10, vec![Cigar::Match(190), Cigar::SoftClip(60)], 50);
        assert_eq!(acc.clips[199][1], 0);
        assert_eq!(w.edge_clips[0], 1);
        // An interior right clip counts (last aligned base 189).
        let (acc, w) = walk(&refseq, 0, vec![Cigar::Match(190), Cigar::SoftClip(60)], 50);
        assert_eq!(acc.clips[189][1], 1);
        assert_eq!(w.edge_clips[0], 0);
    }

    /// Turn one chunk accumulator into per-position counters the way `scan`
    /// merges chunks (coverage is stored as difference arrays).
    fn contig_from_acc(acc: &Acc) -> ContigCounts {
        let len = acc.cov.len() - 1;
        let mut run = [0i32; NCOV];
        let cov = acc
            .cov
            .iter()
            .take(len)
            .map(|d| {
                let mut o = [0u32; NCOV];
                for k in 0..NCOV {
                    run[k] += d[k];
                    o[k] = run[k].max(0) as u32;
                }
                o
            })
            .collect();
        ContigCounts {
            name: "chr".to_string(),
            len,
            bases: acc.bases.clone(),
            cov,
            clips: acc.clips.clone(),
            edge_clips: 0,
            aligned_bases: 0,
        }
    }

    #[test]
    fn aligned_dp_excludes_deletions() {
        use rust_htslib::bam::record::Cigar;
        let refseq = b"ACGT".repeat(150); // 600 bp
        // One read `100M 50D 100M`: positions 200..250 are a deletion only.
        let (acc, _w) = walk(&refseq, 100, vec![Cigar::Match(100), Cigar::Del(50), Cigar::Match(100)], 50);
        let c = contig_from_acc(&acc);
        // Inside the deletion: `dp` counts the read, the aligned-base depth does not.
        assert_eq!(c.dp(220), 1);
        assert_eq!(c.aligned_dp(220), 0);
        // A read base: both agree. Uncovered: both zero.
        assert_eq!((c.dp(150), c.aligned_dp(150)), (1, 1));
        assert_eq!((c.dp(50), c.aligned_dp(50)), (0, 0));
    }

    #[test]
    fn clip_length_sums_terminal_ops() {
        use rust_htslib::bam::record::Cigar;
        let refseq = b"ACGT".repeat(50); // 200 bp
        // 5H60S100M: the left end is 65 bp (5 + 60) >= 50.
        let (acc, w) = walk(&refseq, 10, vec![Cigar::HardClip(5), Cigar::SoftClip(60), Cigar::Match(100)], 50);
        assert_eq!(acc.clips[10][0], 1);
        assert_eq!(w.edge_clips[0], 0);
        // 100M60S5H: the right end is 65 bp (60 + 5) at the last aligned base.
        let (acc, w) = walk(&refseq, 10, vec![Cigar::Match(100), Cigar::SoftClip(60), Cigar::HardClip(5)], 50);
        assert_eq!(acc.clips[109][1], 1);
        assert_eq!(w.edge_clips[0], 0);
        // Summed below threshold: 5H20S100M (25 < 50) does not count.
        let (acc, _w) = walk(&refseq, 10, vec![Cigar::HardClip(5), Cigar::SoftClip(20), Cigar::Match(100)], 50);
        assert_eq!(acc.clips[10][0], 0);
        // A terminal hard clip alone counts (as before).
        let (acc, _w) = walk(&refseq, 10, vec![Cigar::HardClip(60), Cigar::Match(100)], 50);
        assert_eq!(acc.clips[10][0], 1);
        let (acc, _w) = walk(&refseq, 10, vec![Cigar::Match(100), Cigar::HardClip(60)], 50);
        assert_eq!(acc.clips[109][1], 1);
    }
}
