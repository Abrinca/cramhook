//! Assembly diff: identity blocks, T1/T2 taint and novel sequence.
//!
//! This is the M2 core. It is written so the later `lift` step can reuse it:
//! the resulting [`Diff`] answers "which identity block contains this old
//! interval?" and "does this old interval overlap taint?" by binary search, and
//! each [`Block`] maps old positions/intervals to new.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result, bail};
use serde::Serialize;

use crate::coords::{NewIv, NewPos, OldIv, OldPos};
use crate::fasta::{Fasta, ensure_fai};
use crate::paf::{CsStep, PafRecord, Strand, parse_paf, walk_cs};
use crate::version::versioned;

// ---------------------------------------------------------------------------
// Tunables (T2 windowing). One place, documented. Initial values: tune on
// S6/S6b (Plan 6, 4.1).
/// Window length for T2 projection.
pub const WIN: u64 = 10_000;
/// Step between window starts.
pub const STEP: u64 = 2_000;
/// Flank added on both sides of every novel interval before windowing.
pub const FLANK: u64 = 10_000;
/// Half-width of a T2 junction interval added at every new adjacency
/// (`--t2-junction-flank`). The junction itself adds no base, so a window is
/// centred on `[j - J, j + J)` to give T2 something to probe and project.
pub const JUNCTION_FLANK: u64 = 50;

/// T2 windowing parameters (CLI-overridable; defaults are
/// [`WIN`]/[`STEP`]/[`FLANK`]/[`JUNCTION_FLANK`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DiffParams {
    pub win: u64,
    pub step: u64,
    pub flank: u64,
    /// Half-width of the junction interval `[j - J, j + J)` added to the
    /// T2-novel set at every new adjacency. `0` disables junction taint.
    pub junction_flank: u64,
}

impl Default for DiffParams {
    fn default() -> Self {
        Self { win: WIN, step: STEP, flank: FLANK, junction_flank: JUNCTION_FLANK }
    }
}

impl DiffParams {
    /// Upper bound for the T2 windowing knobs (keeps `flank`/`win` products far
    /// from overflow).
    pub const MAX_T2: u64 = 1_000_000_000;

    /// Validate the T2 windowing parameters before doing any work.
    pub fn validate(&self) -> Result<()> {
        if self.win == 0 {
            bail!("--t2-win must be > 0");
        }
        if self.step == 0 {
            bail!("--t2-step must be > 0");
        }
        if self.step > self.win {
            bail!("--t2-step ({}) must be <= --t2-win ({})", self.step, self.win);
        }
        if self.win > Self::MAX_T2 {
            bail!("--t2-win ({}) must be <= {}", self.win, Self::MAX_T2);
        }
        if self.flank > Self::MAX_T2 {
            bail!("--t2-flank ({}) must be <= {}", self.flank, Self::MAX_T2);
        }
        if self.junction_flank > Self::MAX_T2 {
            bail!("--t2-junction-flank ({}) must be <= {}", self.junction_flank, Self::MAX_T2);
        }
        Ok(())
    }
}

/// minimap2 preset for reads, or for the assembly diff.
#[derive(Clone, Copy, PartialEq, Eq, Debug, clap::ValueEnum)]
pub enum Preset {
    #[value(name = "map-hifi")]
    MapHifi,
    #[value(name = "lr:hq")]
    LrHq,
}

impl Preset {
    /// The name minimap2 understands (`-x`).
    pub fn minimap2_name(self) -> &'static str {
        match self {
            Preset::MapHifi => "map-hifi",
            Preset::LrHq => "lr:hq",
        }
    }

    /// Preset used for the asm5 diff (always asm5 for M2).
    pub fn diff_asm_name(self) -> &'static str {
        "asm5"
    }
}

/// An identity block: a maximal exact match run, in old coordinates, plus the
/// new contig, lowest new coordinate and strand it maps to.
#[derive(Clone, Debug)]
pub struct Block {
    pub old_contig: String,
    pub old: OldIv,
    pub new_contig: String,
    pub new_start: NewPos,
    pub strand: Strand,
}

impl Block {
    /// Length of the block (same in old and new coordinates).
    pub fn len(&self) -> u64 {
        self.old.len()
    }

    /// The new interval this block maps to.
    pub fn new_iv(&self) -> NewIv {
        NewIv::new(self.new_start.0, self.new_start.0 + self.old.len())
    }

    /// Map a single old position into new coordinates.
    ///
    /// `+`: `new = new_start + (p - old.start)`.
    /// `-`: `new = new_start + (old.end - 1 - p)`.
    pub fn map_pos(&self, p: OldPos) -> Option<NewPos> {
        if !self.old.contains(p) {
            return None;
        }
        let d = p.0 - self.old.start.0;
        Some(match self.strand {
            Strand::Plus => NewPos(self.new_start.0 + d),
            Strand::Minus => NewPos(self.new_start.0 + (self.old.end.0 - 1 - p.0)),
        })
    }

    /// Map an old interval into new coordinates. The result is half-open and
    /// forward-oriented even for `-` blocks.
    pub fn map_iv(&self, iv: OldIv) -> Option<NewIv> {
        if !self.old.contains_iv(iv) {
            return None;
        }
        Some(match self.strand {
            Strand::Plus => NewIv::new(
                self.new_start.0 + iv.start.0 - self.old.start.0,
                self.new_start.0 + iv.end.0 - self.old.start.0,
            ),
            Strand::Minus => NewIv::new(
                self.new_start.0 + self.old.end.0 - iv.end.0,
                self.new_start.0 + self.old.end.0 - iv.start.0,
            ),
        })
    }
}

/// The result of the diff: sorted blocks and merged taint per old contig, plus
/// merged novel intervals per new contig.
#[derive(Default)]
pub struct Diff {
    pub blocks: HashMap<String, Vec<Block>>,
    pub taint_t1: HashMap<String, Vec<OldIv>>,
    pub taint_t2: HashMap<String, Vec<OldIv>>,
    /// Old bases tainted only by a *junction* (T2 junction taint; see
    /// [`compute_junctions`]). Kept apart from [`Diff::taint_t2`] so the
    /// deferral reason and stats can distinguish them.
    pub taint_junction: HashMap<String, Vec<OldIv>>,
    pub novel: HashMap<String, Vec<NewIv>>,
    /// New-assembly junction intervals `[j - J, j + J)` (T2-only; *not* novel).
    pub junctions: HashMap<String, Vec<NewIv>>,
}

impl Diff {
    /// Find the block containing `iv` on `contig`, or `None`.
    pub fn block_containing(&self, contig: &str, iv: OldIv) -> Option<&Block> {
        let v = self.blocks.get(contig)?;
        let i = v.partition_point(|b| b.old.start.0 <= iv.start.0);
        if i == 0 {
            return None;
        }
        let b = &v[i - 1];
        if b.old.contains_iv(iv) { Some(b) } else { None }
    }

    /// True if `iv` overlaps any T1 or T2 taint region on `contig`.
    pub fn overlaps_taint(&self, contig: &str, iv: OldIv) -> bool {
        overlaps_t1(&self.taint_t1, contig, iv) || overlaps_t1(&self.taint_t2, contig, iv)
            || overlaps_t1(&self.taint_junction, contig, iv)
    }

    /// True if `iv` overlaps a T1 taint region.
    pub fn overlaps_t1(&self, contig: &str, iv: OldIv) -> bool {
        overlaps_t1(&self.taint_t1, contig, iv)
    }

    /// True if `iv` overlaps a T2 taint region.
    pub fn overlaps_t2(&self, contig: &str, iv: OldIv) -> bool {
        overlaps_t1(&self.taint_t2, contig, iv)
    }

    /// True if `iv` overlaps a T2 *junction* taint region.
    pub fn overlaps_junction(&self, contig: &str, iv: OldIv) -> bool {
        overlaps_t1(&self.taint_junction, contig, iv)
    }
}

fn overlaps_t1(map: &HashMap<String, Vec<OldIv>>, contig: &str, iv: OldIv) -> bool {
    let Some(v) = map.get(contig) else {
        return false;
    };
    let i = v.partition_point(|x| x.end.0 <= iv.start.0);
    i < v.len() && v[i].start.0 < iv.end.0
}

fn novel_contains(novel: &HashMap<String, Vec<NewIv>>, contig: &str, pos: u64) -> bool {
    let Some(v) = novel.get(contig) else {
        return false;
    };
    let i = v.partition_point(|iv| iv.end.0 <= pos);
    i < v.len() && v[i].start.0 <= pos
}

/// True if `pos` (a new-contig coordinate in `[0, len)`) lies in any junction
/// interval on `contig`.
fn junction_contains(junctions: &HashMap<String, Vec<NewIv>>, contig: &str, pos: u64) -> bool {
    novel_contains(junctions, contig, pos)
}

/// A T2 window: a slice of the (circular) new assembly used to probe the old
/// assembly.
///
/// `new_start` is always in `[0, contig_len)`. The window is `[new_start,
/// new_start + (new_end - new_start))`; `new_end` may exceed `contig_len` when
/// the window crosses the contig end, and the sequence/projection are read
/// modulo `contig_len`. Because the window length is capped at `contig_len`,
/// at most one wrap occurs.
#[derive(Clone, Debug)]
struct Window {
    contig: String,
    new_start: u64,
    new_end: u64,
    contig_len: u64,
}

// ---------------------------------------------------------------------------
// interval helpers

fn merge_old(mut v: Vec<OldIv>) -> Vec<OldIv> {
    v.retain(|x| !x.is_empty());
    v.sort_by_key(|x| (x.start.0, x.end.0));
    let mut out: Vec<OldIv> = Vec::new();
    for iv in v {
        match out.last_mut() {
            Some(last) if iv.start.0 <= last.end.0 => {
                if iv.end.0 > last.end.0 {
                    last.end = iv.end;
                }
            }
            _ => out.push(iv),
        }
    }
    out
}

/// Return the merged regions covered a number of times other than exactly once
/// (0 or >= 2), clipped to `[0, len)`.
fn regions_not_exactly_once(mut ivs: Vec<(u64, u64)>, len: u64) -> Vec<(u64, u64)> {
    ivs.retain(|&(a, b)| a < b);
    let mut events: Vec<(u64, i64)> = Vec::with_capacity(ivs.len() * 2);
    for (a, b) in &ivs {
        events.push((*a, 1));
        events.push((*b, -1));
    }
    events.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::new();
    let mut cov: i64 = 0;
    let mut prev: u64 = 0;
    let mut i = 0usize;
    while i < events.len() {
        let pos = events[i].0;
        if prev < pos && cov != 1 {
            push_region(&mut out, prev, pos);
        }
        while i < events.len() && events[i].0 == pos {
            cov += events[i].1;
            i += 1;
        }
        prev = pos;
    }
    if prev < len && cov != 1 {
        push_region(&mut out, prev, len);
    }
    out
}

fn push_region(out: &mut Vec<(u64, u64)>, s: u64, e: u64) {
    match out.last_mut() {
        Some(last) if last.1 == s => last.1 = e,
        _ => out.push((s, e)),
    }
}

/// Merge possibly-overlapping `[start, end)` ranges.
fn merge_ranges(mut v: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    v.retain(|&(a, b)| a < b);
    v.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::new();
    for (a, b) in v {
        match out.last_mut() {
            Some(last) if a <= last.1 => {
                if b > last.1 {
                    last.1 = b;
                }
            }
            _ => out.push((a, b)),
        }
    }
    out
}

/// Merge possibly-overlapping new-contig intervals.
fn merge_new(mut v: Vec<NewIv>) -> Vec<NewIv> {
    v.retain(|iv| iv.start.0 < iv.end.0);
    v.sort_by_key(|iv| (iv.start.0, iv.end.0));
    let mut out: Vec<NewIv> = Vec::new();
    for iv in v {
        match out.last_mut() {
            Some(last) if iv.start.0 <= last.end.0 => {
                if iv.end.0 > last.end.0 {
                    last.end = iv.end;
                }
            }
            _ => out.push(iv),
        }
    }
    out
}

/// Push the T2 junction interval `[j - flank, j + flank)` **modulo** the
/// circular contig length `len`, wrapped and clipped into `[0, len)`. When the
/// interval covers the whole contig (`2*flank >= len`) a single `[0, len)` is
/// pushed. `flank == 0` adds nothing (junction taint disabled).
fn push_junction(v: &mut Vec<NewIv>, j: u64, flank: u64, len: u64) {
    if flank == 0 || len == 0 {
        return;
    }
    let len = len as i128;
    let j = (j as i128).rem_euclid(len);
    let flank = flank as i128;
    if 2 * flank >= len {
        v.push(NewIv::new(0, len as u64));
        return;
    }
    let s = j - flank;
    let e = j + flank;
    if s < 0 {
        // The interval wraps the contig start.
        v.push(NewIv::new((len + s) as u64, len as u64));
        v.push(NewIv::new(0, e as u64));
    } else if e > len {
        // The interval wraps the contig end.
        v.push(NewIv::new(s as u64, len as u64));
        v.push(NewIv::new(0, (e - len) as u64));
    } else {
        v.push(NewIv::new(s as u64, e as u64));
    }
}

/// Counts for the `summary.json` junction block.
#[derive(Clone, Copy, Serialize, Default)]
pub struct JunctionStats {
    /// Number of new adjacencies whose old positions are not adjacent.
    pub count: usize,
    /// Total bp of the merged junction intervals.
    pub bp: u64,
    /// Of the junctions, how many are a circular wrap of one old contig
    /// (`old.end == len` / `old.start == 0` on the same strand). Informational:
    /// a wrap is a junction like any other — minimap2 indexes contigs linearly,
    /// so a copy split by the old origin inside the index becomes contiguous
    /// after a rotation while a full remap still sees it.
    pub rotation_wraps: usize,
}

/// Classification of the boundary between two blocks that are adjacent in new
/// coordinates (identical `new_start` continuity is checked by the caller).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum JunctionKind {
    /// The old positions are adjacent too: no junction.
    Natural,
    /// The old positions are not adjacent: a genuine new adjacency.
    Edit,
    /// Circular-adjacent on one old contig (`old.end == len` / `old.start == 0`,
    /// same strand). Reported separately, but still a junction.
    RotationWrap,
}

fn junction_kind(b1: &Block, b2: &Block, old_len: u64) -> JunctionKind {
    let same = b1.old_contig == b2.old_contig;
    match (b1.strand, b2.strand) {
        (Strand::Plus, Strand::Plus) if same && b1.old.end.0 == b2.old.start.0 => {
            JunctionKind::Natural
        }
        (Strand::Minus, Strand::Minus) if same && b1.old.start.0 == b2.old.end.0 => {
            JunctionKind::Natural
        }
        (Strand::Plus, Strand::Plus)
            if same && old_len != 0 && b1.old.end.0 == old_len && b2.old.start.0 == 0 =>
        {
            JunctionKind::RotationWrap
        }
        (Strand::Minus, Strand::Minus)
            if same && old_len != 0 && b1.old.start.0 == 0 && b2.old.end.0 == old_len =>
        {
            JunctionKind::RotationWrap
        }
        _ => JunctionKind::Edit,
    }
}

/// Deterministic order over the blocks of one new contig.
fn block_key(b: &Block) -> (u64, String, u64, u64, u8) {
    (
        b.new_start.0,
        b.old_contig.clone(),
        b.old.start.0,
        b.old.end.0,
        u8::from(b.strand == Strand::Minus),
    )
}

/// Deterministic ending block at a block-end position: largest `new_start`, then
/// the smallest old key.
fn pick_end(enders: &[usize], blks: &[&Block]) -> Option<usize> {
    enders
        .iter()
        .copied()
        .min_by_key(|&i| {
            let b = blks[i];
            (std::cmp::Reverse(b.new_start.0), block_key(b))
        })
}

/// Deterministic successor block that starts at `e` (the ending block excluded).
fn pick_start(cands: Option<&Vec<usize>>, skip: usize, blks: &[&Block]) -> Option<usize> {
    cands?
        .iter()
        .copied()
        .filter(|&i| i != skip)
        .min_by_key(|&i| block_key(blks[i]))
}

/// Find every new-assembly junction: an adjacency of the identity-block tiling
/// where the old coordinates are not adjacent (a deletion inside a colinear
/// chain, a contig join, a strand flip or a circular wrap). Such an edit adds
/// no new base, so the new position is *covered exactly once* and
/// [`compute_diff_with`]'s `novel` set is empty there; T2 would never probe it.
/// The junction interval `[j - flank, j + flank)` is added to a T2-only set
/// (never to `novel`, so the clip screen is unaffected).
///
/// Junctions are found from **block-end events**: at every new position where
/// some block ends, the successor is a block that starts there. Nested or
/// overlapping blocks cannot hide an adjacency, because a position strictly
/// inside another block is multiply covered (`novel`) and rejected before any
/// successor lookup. The last-block→first-block seam of a new contig (all new
/// contigs are treated as circular, see M2) is the block-end event at the contig
/// length. Block ends and starts are indexed by position, so the whole pass is
/// `O(B log B)` (plus the per-contig sort) rather than `O(B²)`.
///
/// With `flank == 0` junction taint is disabled and the pass returns early.
fn compute_junctions(
    blocks: &HashMap<String, Vec<Block>>,
    old: &Fasta,
    new: &Fasta,
    novel: &HashMap<String, Vec<NewIv>>,
    flank: u64,
) -> (HashMap<String, Vec<NewIv>>, JunctionStats) {
    if flank == 0 {
        return (HashMap::new(), JunctionStats::default());
    }
    let mut by_new: HashMap<String, Vec<&Block>> = HashMap::new();
    for blks in blocks.values() {
        for b in blks {
            by_new.entry(b.new_contig.clone()).or_default().push(b);
        }
    }
    let mut junctions: HashMap<String, Vec<NewIv>> = HashMap::new();
    let mut stats = JunctionStats::default();
    for (contig, mut blks) in by_new {
        let len = new.len_of(&contig).unwrap_or(0);
        if len == 0 {
            continue;
        }
        blks.sort_by_key(|b| block_key(b));
        // Index the blocks by start and end position (modulo the contig length),
        // so each candidate is an O(1)-ish lookup instead of a full scan.
        let mut by_start: HashMap<u64, Vec<usize>> = HashMap::new();
        let mut by_end: HashMap<u64, Vec<usize>> = HashMap::new();
        for (i, b) in blks.iter().enumerate() {
            by_start.entry(b.new_start.0 % len).or_default().push(i);
            by_end
                .entry((b.new_start.0 + b.old.len()) % len)
                .or_default()
                .push(i);
        }
        // Candidate positions: every block end, modulo the contig length (so the
        // seam of a circular new contig is the end at `len` -> 0).
        let mut ends: Vec<u64> = by_end.keys().copied().collect();
        ends.sort_unstable();
        let mut ivs: Vec<NewIv> = Vec::new();
        for &e in &ends {
            // Reject multiply-covered positions *before* any block scan: a
            // position adjacent to a novel region is not a clean adjacency.
            let before = if e == 0 { len - 1 } else { e - 1 };
            if novel_contains(novel, &contig, before) || novel_contains(novel, &contig, e) {
                continue;
            }
            let Some(ei) = pick_end(&by_end[&e], &blks) else {
                continue;
            };
            let b_end = blks[ei];
            // A successor that merely *contains* `e` would overlap `b_end` at
            // `e-1` (novel), already rejected above, so only a block that starts
            // at `e` can succeed.
            let Some(ni) = pick_start(by_start.get(&e), ei, &blks) else {
                continue;
            };
            let b_next = blks[ni];
            let old_len = old.len_of(&b_end.old_contig).unwrap_or(0);
            match junction_kind(b_end, b_next, old_len) {
                JunctionKind::Natural => continue,
                JunctionKind::Edit => {}
                JunctionKind::RotationWrap => stats.rotation_wraps += 1,
            }
            stats.count += 1;
            push_junction(&mut ivs, e, flank, len);
        }
        if !ivs.is_empty() {
            junctions.insert(contig, merge_new(ivs));
        }
    }
    stats.bp = junctions.values().flatten().map(|iv| iv.len()).sum();
    (junctions, stats)
}

/// Tile `[start, end)` with windows of length `win` and step `step`. The last
/// window ends exactly at `end`; if the interval is shorter than `win` a single
/// whole-interval window is returned.
fn tile_windows(start: u64, end: u64, win: u64, step: u64) -> Vec<(u64, u64)> {
    debug_assert!(start < end);
    if end - start <= win {
        return vec![(start, end)];
    }
    let mut out = Vec::new();
    let mut s = start;
    loop {
        if s + win >= end {
            out.push((end - win, end));
            break;
        }
        out.push((s, s + win));
        s += step;
    }
    out
}

/// Query-forward offset -> original query coordinate for a cs walk.
#[inline]
fn qcoord(qstart: u64, qend: u64, strand: Strand, qoff: u64) -> u64 {
    match strand {
        Strand::Plus => qstart + qoff,
        Strand::Minus => qend - 1 - qoff,
    }
}

/// Reduce an (unwrapped) window offset to a contig coordinate.
#[inline]
fn circ_pos(start: u64, off: u64, contig_len: u64) -> u64 {
    if contig_len == 0 {
        return start.saturating_add(off);
    }
    (start + off) % contig_len
}

/// Read `len` bases of `seq` starting at `start`, wrapping past the end.
/// `seq` is the circular contig; `len` may exceed `seq.len()`.
fn circular_seq(seq: &[u8], start: u64, len: u64) -> Vec<u8> {
    let n = seq.len();
    let mut out = Vec::with_capacity(len as usize);
    if n == 0 {
        return out;
    }
    for i in 0..len {
        out.push(seq[((start + i) % n as u64) as usize]);
    }
    out
}

/// Build the circular T2 windows for one novel interval: the interval plus
/// `flank` on both sides, tiled with `win`/`step`. The arc is capped at the
/// whole contig; when it covers the whole contig the whole circle is tiled from
/// every start `0..contig_len` by `step`, so the origin seam is probed as well
/// (otherwise a wholly novel contig shorter than `win` yields one window with no
/// seam-crossing probe). Windows may cross the contig end. All new contigs are
/// treated as circular: on a linear contig the extra junction is artificial but
/// only ever *adds* taint, so the probe stays conservative (T2 errors fail
/// unsafe).
fn interval_windows(
    contig: &str,
    contig_len: u64,
    iv_start: u64,
    iv_end: u64,
    params: DiffParams,
) -> Vec<Window> {
    if contig_len == 0 || iv_start >= iv_end {
        return Vec::new();
    }
    let arc_len = (iv_end - iv_start)
        .saturating_add(params.flank.saturating_mul(2))
        .min(contig_len);
    // Novel interval start minus flank, modulo the contig length.
    let arc_start =
        (iv_start as i128 - params.flank as i128).rem_euclid(contig_len as i128) as u64;
    let mut out: Vec<Window> =
        tile_windows(arc_start, arc_start.saturating_add(arc_len), params.win, params.step)
            .into_iter()
            .map(|(s, e)| Window {
                contig: contig.to_string(),
                new_start: s % contig_len,
                new_end: s % contig_len + (e - s),
                contig_len,
            })
            .collect();
    // When the arc covers the whole contig the tiling can end exactly at the
    // contig end and never wrap, leaving the origin seam unprobed (e.g. a wholly
    // novel contig with `novel == [0, contig_len)`). Add one wrapping window so
    // matches that exist only across the origin are seen.
    if arc_len == contig_len && !out.iter().any(|w| w.new_start > 0 && w.new_end > contig_len) {
        let win = params.win.min(contig_len).max(1);
        let start = contig_len.saturating_sub(win / 2).clamp(1, contig_len - 1);
        out.push(Window {
            contig: contig.to_string(),
            new_start: start,
            new_end: start + win,
            contig_len,
        });
    }
    out
}

// ---------------------------------------------------------------------------
// public API

/// Everything the CLI needs to write outputs, plus the reusable [`Diff`].
pub struct DiffOutput {
    pub diff: Diff,
    pub diff_paf: String,
    pub t2_paf: String,
    pub summary: Summary,
}

/// Count + total bp of an interval or block set.
#[derive(Clone, Copy, Serialize, Default)]
pub struct CountBp {
    pub count: usize,
    pub bp: u64,
}

/// Taint summary per reason.
#[derive(Clone, Copy, Serialize, Default)]
pub struct TaintSummary {
    pub t1: CountBp,
    pub t2: CountBp,
    /// Old bp tainted only by junction taint (`taint_junction`).
    pub t2_junction: CountBp,
}

/// Wall time (seconds) of each step.
#[derive(Clone, Copy, Serialize, Default)]
pub struct Timings {
    pub diff: f64,
    pub t1: f64,
    pub novel: f64,
    pub t2: f64,
    pub total: f64,
}

/// `P.summary.json` contents.
#[derive(Clone, Serialize)]
pub struct Summary {
    pub preset: String,
    pub threads: usize,
    pub n_primary_alignments: usize,
    pub blocks: CountBp,
    pub taint: TaintSummary,
    pub novel: CountBp,
    /// T2 junction intervals (also `novel`-like windows, but kept out of
    /// `novel` so the clip screen is unaffected).
    pub junctions: CountBp,
    /// Junction sites detected (before merging; each may be a rotation wrap).
    pub junctions_detected: usize,
    /// Of `junctions_detected`, how many are a circular wrap of one old contig
    /// (`old.end == len` / `old.start == 0`, same strand). Informational.
    pub junction_rotation_wraps: usize,
    pub windows: usize,
    /// T2 windowing parameters actually used.
    pub t2_win: u64,
    pub t2_step: u64,
    pub t2_flank: u64,
    pub t2_junction_flank: u64,
    pub timings_s: Timings,
}

/// Run the full M2 diff over two in-memory assemblies. `old_path`/`new_path`
/// are needed for the minimap2 invocations. A temporary directory is created
/// for the T2 windows and removed before returning.
pub fn compute_diff(
    old: &Fasta,
    new: &Fasta,
    old_path: &Path,
    new_path: &Path,
    preset: Preset,
    threads: usize,
) -> Result<DiffOutput> {
    compute_diff_with(old, new, old_path, new_path, preset, threads, DiffParams::default())
}

/// [`compute_diff`] with explicit T2 windowing parameters.
pub fn compute_diff_with(
    old: &Fasta,
    new: &Fasta,
    old_path: &Path,
    new_path: &Path,
    preset: Preset,
    threads: usize,
    params: DiffParams,
) -> Result<DiffOutput> {
    params.validate()?;
    let t_all = Instant::now();

    // --- step 1: asm5 diff, query = old, target = new -----------------------
    let t = Instant::now();
    let diff_paf = run_minimap2(&[
        "-cx".into(),
        preset.diff_asm_name().into(),
        "--cs".into(),
        "--secondary=no".into(),
        "-t".into(),
        threads.to_string(),
        new_path.display().to_string(),
        old_path.display().to_string(),
    ])?;
    let records = parse_paf(&diff_paf)?;
    let primary: Vec<&PafRecord> = records.iter().filter(|r| r.primary).collect();
    let n_primary = primary.len();
    let t_diff = t.elapsed().as_secs_f64();

    // --- step 3: identity blocks (one per `:N` run) -------------------------
    let mut blocks: HashMap<String, Vec<Block>> = HashMap::new();
    for rec in &primary {
        walk_cs(&rec.cs, |step, t, qoff| {
            if let CsStep::Match(n) = step {
                if n == 0 {
                    return;
                }
                let old_iv = query_old_iv(rec, qoff, n);
                blocks.entry(rec.qname.clone()).or_default().push(Block {
                    old_contig: rec.qname.clone(),
                    old: old_iv,
                    new_contig: rec.tname.clone(),
                    new_start: NewPos(rec.tstart + t),
                    strand: rec.strand,
                });
            }
        })?;
    }
    for v in blocks.values_mut() {
        v.sort_by_key(|b| (b.old.start.0, b.old.end.0));
    }

    // --- step 4: T1 taint ----------------------------------------------------
    // Covered 0 or >= 2 times by primary alignments. We sweep the union of the
    // primary query intervals and the identity blocks so that sequence inside an
    // alignment's deletion gap (e.g. a collapsed duplication copy) counts as
    // "covered by no primary alignment (deleted)".
    let t = Instant::now();
    let mut align_ivs: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for rec in &primary {
        align_ivs.entry(rec.qname.clone()).or_default().push((rec.qstart, rec.qend));
    }
    let mut taint_t1: HashMap<String, Vec<OldIv>> = HashMap::new();
    for contig in old.order() {
        let len = old.len_of(contig).unwrap_or(0);
        let block_ivs: Vec<(u64, u64)> = blocks
            .get(contig)
            .map(|v| v.iter().map(|b| (b.old.start.0, b.old.end.0)).collect())
            .unwrap_or_default();
        let mut regions = regions_not_exactly_once(align_ivs.remove(contig).unwrap_or_default(), len);
        regions.extend(regions_not_exactly_once(block_ivs, len));
        let merged = merge_ranges(regions)
            .into_iter()
            .map(|(a, b)| OldIv::new(a, b))
            .collect::<Vec<_>>();
        taint_t1.insert(contig.clone(), merged);
    }
    let t_t1 = t.elapsed().as_secs_f64();

    // --- step 5: novel (new coverage by identity blocks) --------------------
    let t = Instant::now();
    let mut novel: HashMap<String, Vec<NewIv>> = HashMap::new();
    for contig in new.order() {
        let len = new.len_of(contig).unwrap_or(0);
        let mut ivs: Vec<(u64, u64)> = Vec::new();
        for (_, blks) in blocks.iter() {
            for b in blks {
                if &b.new_contig == contig {
                    ivs.push((b.new_start.0, b.new_start.0 + b.old.len()));
                }
            }
        }
        let merged = regions_not_exactly_once(ivs, len)
            .into_iter()
            .map(|(a, b)| NewIv::new(a, b))
            .collect::<Vec<_>>();
        novel.insert(contig.clone(), merged);
    }
    let t_novel = t.elapsed().as_secs_f64();

    // --- step 5b: junctions (new adjacencies with non-adjacent old ends) ----
    let (junctions, junction_stats) = compute_junctions(&blocks, old, new, &novel, params.junction_flank);
    // Identity blocks grouped by new contig, for the junction self-hit test.
    let mut blocks_by_new: HashMap<String, Vec<Block>> = HashMap::new();
    for blks in blocks.values() {
        for b in blks {
            blocks_by_new.entry(b.new_contig.clone()).or_default().push(b.clone());
        }
    }
    for v in blocks_by_new.values_mut() {
        v.sort_by_key(|b| (b.new_start.0, b.old.start.0));
    }

    // --- step 6: T2 taint ---------------------------------------------------
    let t = Instant::now();
    let targets =
        T2Targets { novel: &novel, junctions: &junctions, blocks_by_new: &blocks_by_new };
    let (taint, t2_paf, n_windows) =
        compute_t2(old, new, old_path, &targets, preset, threads, params)?;
    let T2Taint { novel: taint_t2, junction: taint_junction } = taint;
    let t_t2 = t.elapsed().as_secs_f64();

    let diff = Diff { blocks, taint_t1, taint_t2, taint_junction, novel, junctions };

    // --- summary ------------------------------------------------------------
    let blocks_cb = CountBp {
        count: diff.blocks.values().map(|v| v.len()).sum(),
        bp: diff.blocks.values().flatten().map(|b| b.old.len()).sum(),
    };
    let t1_cb = count_bp_old(diff.taint_t1.values().flatten());
    let t2_cb = count_bp_old(diff.taint_t2.values().flatten());
    let t2j_cb = count_bp_old(diff.taint_junction.values().flatten());
    let novel_cb = CountBp {
        count: diff.novel.values().map(|v| v.len()).sum(),
        bp: diff.novel.values().flatten().map(|iv| iv.len()).sum(),
    };
    let junction_cb = CountBp {
        count: diff.junctions.values().map(|v| v.len()).sum(),
        bp: diff.junctions.values().flatten().map(|iv| iv.len()).sum(),
    };
    let total_s = {
        // Make sure all step times sum to something at least as large.
        let _ = t_all;
        Instant::now().duration_since(t_all).as_secs_f64()
    };
    let timings = Timings { diff: t_diff, t1: t_t1, novel: t_novel, t2: t_t2, total: total_s };
    let summary = Summary {
        preset: preset.minimap2_name().to_string(),
        threads,
        n_primary_alignments: n_primary,
        blocks: blocks_cb,
        taint: TaintSummary { t1: t1_cb, t2: t2_cb, t2_junction: t2j_cb },
        novel: novel_cb,
        junctions: junction_cb,
        junctions_detected: junction_stats.count,
        junction_rotation_wraps: junction_stats.rotation_wraps,
        windows: n_windows,
        t2_win: params.win,
        t2_step: params.step,
        t2_flank: params.flank,
        t2_junction_flank: params.junction_flank,
        timings_s: timings,
    };

    Ok(DiffOutput { diff, diff_paf, t2_paf, summary })
}

fn query_old_iv(rec: &PafRecord, qoff: u64, n: u64) -> OldIv {
    match rec.strand {
        Strand::Plus => OldIv::new(rec.qstart + qoff, rec.qstart + qoff + n),
        Strand::Minus => OldIv::new(rec.qend - qoff - n, rec.qend - qoff),
    }
}

fn count_bp_old<'a>(ivs: impl Iterator<Item = &'a OldIv>) -> CountBp {
    let mut cb = CountBp::default();
    for iv in ivs {
        cb.count += 1;
        cb.bp += iv.len();
    }
    cb
}

/// T2 probe targets, all per new contig. Passing them together keeps
/// `compute_t2`'s argument count small.
struct T2Targets<'a> {
    novel: &'a HashMap<String, Vec<NewIv>>,
    junctions: &'a HashMap<String, Vec<NewIv>>,
    /// Identity blocks grouped by new contig, for the junction self-hit test.
    blocks_by_new: &'a HashMap<String, Vec<Block>>,
}

/// The identity block covering new position `p` on `contig`, mapped to the old
/// position and strand it comes from (the "source" of `p`). Deterministic when
/// blocks overlap: the block with the largest `new_start` (innermost), then the
/// smallest `(old_contig, old_start)`.
fn source_of(
    blocks_by_new: &HashMap<String, Vec<Block>>,
    contig: &str,
    p: u64,
) -> Option<(String, u64, Strand)> {
    let v = blocks_by_new.get(contig)?;
    let mut best: Option<&Block> = None;
    for b in v {
        if b.new_start.0 > p || p >= b.new_start.0 + b.old.len() {
            continue;
        }
        let better = match best {
            None => true,
            Some(c) => {
                (b.new_start.0, std::cmp::Reverse((&b.old_contig, b.old.start.0)))
                    > (c.new_start.0, std::cmp::Reverse((&c.old_contig, c.old.start.0)))
            }
        };
        if better {
            best = Some(b);
        }
    }
    let b = best?;
    let old_pos = match b.strand {
        Strand::Plus => b.old.start.0 + (p - b.new_start.0),
        Strand::Minus => b.old.end.0 - 1 - (p - b.new_start.0),
    };
    Some((b.old_contig.clone(), old_pos, b.strand))
}

/// Is the hit `(tname, strand, old_pos)` simply the window aligning back to the
/// source of new position `p`? Such a self-hit taints the junction's own
/// sequence and is ignored (see `decisions.md`: reads on the source next to a
/// junction are covered by T1 for a deletion or by the origin-crossing split
/// handling for a wrap).
fn is_self_hit(
    blocks_by_new: &HashMap<String, Vec<Block>>,
    contig: &str,
    p: u64,
    tname: &str,
    strand: Strand,
    old_pos: u64,
) -> bool {
    match source_of(blocks_by_new, contig, p) {
        Some((oc, op, st)) => oc == tname && op == old_pos && st == strand,
        None => false,
    }
}

/// T2 taint split by provenance: novel-base taint and junction taint.
struct T2Taint {
    novel: HashMap<String, Vec<OldIv>>,
    junction: HashMap<String, Vec<OldIv>>,
}

/// Build T2 windows from novel intervals **and junction intervals**, map them
/// back to `old` and project the hits onto old coordinates.
///
/// **Rationale for projection instead of whole-hit taint:** a window's flanks
/// map back to their own origin, so tainting a whole hit would taint roughly
/// `±WIN` around every polishing edit. Only the parts of a hit that are similar
/// to *novel* new sequence can make a read's position ambiguous, also because
/// T2 errors fail unsafe. So we walk each hit's `cs` and taint an old base only
/// when the window base aligned to it is novel.
///
/// Junction intervals are windows too (a new adjacency adds no base, so there is
/// no `novel` position to probe); their taint is collected separately
/// (`taint_junction`) so the deferral reason and stats can distinguish it.
fn compute_t2(
    old: &Fasta,
    new: &Fasta,
    old_path: &Path,
    targets: &T2Targets<'_>,
    preset: Preset,
    threads: usize,
    params: DiffParams,
) -> Result<(T2Taint, String, usize)> {
    // T2 windows are circular: an extended novel interval that crosses either
    // contig end wraps around instead of being clipped (see `interval_windows`).
    let mut windows: Vec<Window> = Vec::new();
    for contig in new.order() {
        let len = new.len_of(contig).unwrap_or(0);
        for ivs in [targets.novel.get(contig), targets.junctions.get(contig)]
            .into_iter()
            .flatten()
        {
            for iv in ivs {
                windows.extend(interval_windows(contig, len, iv.start.0, iv.end.0, params));
            }
        }
    }

    if windows.is_empty() {
        return Ok((T2Taint { novel: HashMap::new(), junction: HashMap::new() }, String::new(), 0));
    }

    let tmp = unique_tmp_dir()?;
    let result = (|| {
        let win_fa = tmp.join("windows.fa");
        write_windows_fasta(new, &windows, &win_fa)?;
        let paf = run_minimap2(&[
            "-c".into(),
            "--cs".into(),
            "-x".into(),
            preset.minimap2_name().into(),
            "-N".into(),
            "1000000".into(),
            "-p".into(),
            "0.000001".into(),
            "-t".into(),
            threads.to_string(),
            old_path.display().to_string(),
            win_fa.display().to_string(),
        ])?;
        // All hits (primary and secondary) are kept.
        let records = parse_paf(&paf)?;
        let mut raw_novel: Vec<(String, OldIv)> = Vec::new();
        let mut raw_junction: Vec<(String, OldIv)> = Vec::new();
        project_t2_hits(&records, &windows, targets, &mut raw_novel, &mut raw_junction)?;
        let merge = |raw: Vec<(String, OldIv)>| -> HashMap<String, Vec<OldIv>> {
            let mut by_contig: HashMap<String, Vec<OldIv>> = HashMap::new();
            for (c, iv) in raw {
                by_contig.entry(c).or_default().push(iv);
            }
            by_contig.into_iter().map(|(c, v)| (c, merge_old(v))).collect()
        };
        Ok((
            T2Taint { novel: merge(raw_novel), junction: merge(raw_junction) },
            paf,
            windows.len(),
        ))
    })();

    let _ = std::fs::remove_dir_all(&tmp);
    let _ = old; // lengths come from rec.tlen; kept for symmetry/future use
    result
}

fn write_windows_fasta(new: &Fasta, windows: &[Window], path: &Path) -> Result<()> {
    let mut buf = String::new();
    for (i, w) in windows.iter().enumerate() {
        let seq = new
            .seq(&w.contig)
            .with_context(|| format!("window contig {:?} not in new FASTA", w.contig))?;
        if w.contig_len == 0 || w.contig_len as usize != seq.len() {
            bail!(
                "window {i} contig {:?} length {} != FASTA length {}",
                w.contig,
                w.contig_len,
                seq.len()
            );
        }
        if w.new_end <= w.new_start {
            bail!("window {i} [{}..{}) is empty on {:?}", w.new_start, w.new_end, w.contig);
        }
        // Circular read: `new_end - new_start` may extend past the contig end.
        let win = circular_seq(seq, w.new_start, w.new_end - w.new_start);
        buf.push('>');
        buf.push_str(&format!("w{i}\n"));
        // Sequences here are <= WIN bases; one line is fine.
        buf.push_str(std::str::from_utf8(&win).context("non-utf8 window")?);
        buf.push('\n');
    }
    std::fs::write(path, buf).with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

fn project_t2_hits(
    records: &[PafRecord],
    windows: &[Window],
    targets: &T2Targets<'_>,
    raw_novel: &mut Vec<(String, OldIv)>,
    raw_junction: &mut Vec<(String, OldIv)>,
) -> Result<()> {
    for rec in records {
        let widx: usize = rec
            .qname
            .strip_prefix('w')
            .and_then(|s| s.parse().ok())
            .with_context(|| format!("unexpected T2 query name {:?}", rec.qname))?;
        let win = windows
            .get(widx)
            .with_context(|| format!("unknown T2 window index {widx}"))?;
        let qspan = rec.q_span();
        let tstart = rec.tstart;
        // Window-relative query coordinate -> new-contig coordinate (circular).
        let wpos = |qc: u64| circ_pos(win.new_start, qc, win.contig_len);
        // Which set the window base at new position `p` belongs to:
        // `Some(true)` = novel, `Some(false)` = junction, `None` = neither.
        // Novel wins when a position is in both, so `taint_t2` and all existing
        // novel behaviour are unchanged.
        let set_of = |p: u64| -> Option<bool> {
            if novel_contains(targets.novel, &win.contig, p) {
                Some(true)
            } else if junction_contains(targets.junctions, &win.contig, p) {
                Some(false)
            } else {
                None
            }
        };
        // For an aligned base: a junction base counts only when the hit is not
        // the window aligning back to its own source (self-hit). The self-hit
        // taints the junction's own sequence; reads there are covered by T1
        // (deletion) or the origin-crossing split handling (wrap). Novel bases
        // are never excluded.
        let base_kind = |qc: u64, old_pos: u64| -> Option<bool> {
            let p = wpos(qc);
            match set_of(p) {
                Some(true) => Some(true),
                Some(false) => {
                    if is_self_hit(
                        targets.blocks_by_new,
                        &win.contig,
                        p,
                        &rec.tname,
                        rec.strand,
                        old_pos,
                    ) {
                        None
                    } else {
                        Some(false)
                    }
                }
                None => None,
            }
        };
        let mut pushes_novel: Vec<OldIv> = Vec::new();
        let mut pushes_junction: Vec<OldIv> = Vec::new();
        walk_cs(&rec.cs, |step, t, qoff| match step {
            CsStep::Match(n) => {
                for j in 0..n {
                    let qc = qcoord(rec.qstart, rec.qend, rec.strand, qoff + j);
                    match base_kind(qc, tstart + t + j) {
                        Some(true) => {
                            pushes_novel.push(OldIv::new(tstart + t + j, tstart + t + j + 1))
                        }
                        Some(false) => {
                            pushes_junction.push(OldIv::new(tstart + t + j, tstart + t + j + 1))
                        }
                        None => {}
                    }
                }
            }
            CsStep::Sub => {
                let qc = qcoord(rec.qstart, rec.qend, rec.strand, qoff);
                match base_kind(qc, tstart + t) {
                    Some(true) => pushes_novel.push(OldIv::new(tstart + t, tstart + t + 1)),
                    Some(false) => pushes_junction.push(OldIv::new(tstart + t, tstart + t + 1)),
                    None => {}
                }
            }
            CsStep::Ins(n) => {
                let mut nov = false;
                let mut junc = false;
                for j in 0..n {
                    let qc = qcoord(rec.qstart, rec.qend, rec.strand, qoff + j);
                    match set_of(wpos(qc)) {
                        Some(true) => {
                            nov = true;
                            break;
                        }
                        Some(false) => junc = true,
                        None => {}
                    }
                }
                if nov || junc {
                    let a = (tstart + t).saturating_sub(1);
                    let b = (tstart + t + 1).min(rec.tlen);
                    if a < b {
                        let iv = OldIv::new(a, b);
                        if nov {
                            pushes_novel.push(iv);
                        } else {
                            pushes_junction.push(iv);
                        }
                    }
                }
            }
            CsStep::Del(n) => {
                let mut nov = false;
                let mut junc = false;
                if qoff > 0 {
                    match set_of(wpos(qcoord(rec.qstart, rec.qend, rec.strand, qoff - 1))) {
                        Some(true) => nov = true,
                        Some(false) => junc = true,
                        None => {}
                    }
                }
                if !nov && qoff < qspan {
                    match set_of(wpos(qcoord(rec.qstart, rec.qend, rec.strand, qoff))) {
                        Some(true) => nov = true,
                        Some(false) => junc = true,
                        None => {}
                    }
                }
                if nov || junc {
                    let iv = OldIv::new(tstart + t, tstart + t + n);
                    if nov {
                        pushes_novel.push(iv);
                    } else {
                        pushes_junction.push(iv);
                    }
                }
            }
        })?;
        for iv in pushes_novel {
            raw_novel.push((rec.tname.clone(), iv));
        }
        for iv in pushes_junction {
            raw_junction.push((rec.tname.clone(), iv));
        }
    }
    Ok(())
}

fn run_minimap2(args: &[String]) -> Result<String> {
    let out = Command::new("minimap2").args(args).output().context("running minimap2")?;
    if !out.status.success() {
        bail!(
            "minimap2 failed ({}): {}",
            out.status,
            String::from_utf8_lossy(&out.stderr).trim()
        );
    }
    String::from_utf8(out.stdout).context("minimap2 output was not valid UTF-8")
}

fn unique_tmp_dir() -> Result<PathBuf> {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("cramhook-t2-{}-{}", std::process::id(), nanos));
    std::fs::create_dir_all(&dir)
        .with_context(|| format!("creating temp dir {}", dir.display()))?;
    Ok(dir)
}

// ---------------------------------------------------------------------------
// top-level convenience used by the binary and the integration tests

/// Result of [`run_diff`].
pub struct DiffRun {
    pub diff: Diff,
    pub summary: Summary,
}

/// Load both assemblies, run the diff and write all `P.*` outputs.
pub fn run_diff(
    old_path: &Path,
    new_path: &Path,
    preset: Preset,
    threads: usize,
    out_prefix: &Path,
) -> Result<DiffRun> {
    run_diff_with(old_path, new_path, preset, threads, out_prefix, DiffParams::default())
}

/// [`run_diff`] with explicit T2 windowing parameters.
pub fn run_diff_with(
    old_path: &Path,
    new_path: &Path,
    preset: Preset,
    threads: usize,
    out_prefix: &Path,
    params: DiffParams,
) -> Result<DiffRun> {
    ensure_fai(old_path)?;
    ensure_fai(new_path)?;
    let old = Fasta::load(old_path)?;
    let new = Fasta::load(new_path)?;
    let out = compute_diff_with(&old, &new, old_path, new_path, preset, threads, params)?;
    write_outputs(out_prefix, &out)?;
    Ok(DiffRun { diff: out.diff, summary: out.summary })
}

/// Write `P.diff.paf`, `P.t2.paf`, `P.blocks.bed`, `P.taint.bed`,
/// `P.novel.bed`, `P.junctions.bed` and `P.summary.json`.
pub fn write_outputs(prefix: &Path, out: &DiffOutput) -> Result<()> {
    if let Some(parent) = prefix.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("creating {}", parent.display()))?;
        }
    }
    std::fs::write(with_suffix(prefix, ".diff.paf"), &out.diff_paf)
        .with_context(|| format!("writing {}", with_suffix(prefix, ".diff.paf").display()))?;
    std::fs::write(with_suffix(prefix, ".t2.paf"), &out.t2_paf)
        .with_context(|| format!("writing {}", with_suffix(prefix, ".t2.paf").display()))?;

    // blocks.bed: old_contig old_start old_end new_contig new_start strand
    let mut block_rows: Vec<String> = Vec::new();
    let mut contigs: Vec<&String> = out.diff.blocks.keys().collect();
    contigs.sort();
    for contig in contigs {
        for b in &out.diff.blocks[contig] {
            block_rows.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                b.old_contig,
                b.old.start.0,
                b.old.end.0,
                b.new_contig,
                b.new_start.0,
                b.strand.as_char()
            ));
        }
    }
    std::fs::write(with_suffix(prefix, ".blocks.bed"), join_lines(&block_rows))?;

    // taint.bed: old_contig start end reason
    let mut taint_rows: Vec<String> = Vec::new();
    for (reason, map) in [
        ("T1", &out.diff.taint_t1),
        ("T2", &out.diff.taint_t2),
        ("T2J", &out.diff.taint_junction),
    ] {
        let mut cs: Vec<&String> = map.keys().collect();
        cs.sort();
        for contig in cs {
            for iv in &map[contig] {
                taint_rows.push(format!("{}\t{}\t{}\t{}", contig, iv.start.0, iv.end.0, reason));
            }
        }
    }
    std::fs::write(with_suffix(prefix, ".taint.bed"), join_lines(&taint_rows))?;

    // novel.bed: new_contig start end
    let mut novel_rows: Vec<String> = Vec::new();
    let mut ncs: Vec<&String> = out.diff.novel.keys().collect();
    ncs.sort();
    for contig in ncs {
        for iv in &out.diff.novel[contig] {
            novel_rows.push(format!("{}\t{}\t{}", contig, iv.start.0, iv.end.0));
        }
    }
    std::fs::write(with_suffix(prefix, ".novel.bed"), join_lines(&novel_rows))?;

    // junctions.bed: new_contig start end (T2-only junction intervals)
    let mut junction_rows: Vec<String> = Vec::new();
    let mut jcs: Vec<&String> = out.diff.junctions.keys().collect();
    jcs.sort();
    for contig in jcs {
        for iv in &out.diff.junctions[contig] {
            junction_rows.push(format!("{}\t{}\t{}", contig, iv.start.0, iv.end.0));
        }
    }
    std::fs::write(with_suffix(prefix, ".junctions.bed"), join_lines(&junction_rows))?;

    let f = std::fs::File::create(with_suffix(prefix, ".summary.json"))
        .with_context(|| "creating summary.json")?;
    serde_json::to_writer_pretty(std::io::BufWriter::new(f), &versioned(&out.summary))
        .context("writing summary.json")?;
    Ok(())
}

fn join_lines(rows: &[String]) -> String {
    let mut s = rows.join("\n");
    if !s.is_empty() {
        s.push('\n');
    }
    s
}

fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diff_params_validation() {
        assert!(DiffParams::default().validate().is_ok());
        assert!(DiffParams { win: 0, step: 1, flank: 0, junction_flank: 0 }.validate().is_err());
        assert!(DiffParams { win: 10, step: 0, flank: 0, junction_flank: 0 }.validate().is_err());
        assert!(DiffParams { win: 10, step: 11, flank: 0, junction_flank: 0 }.validate().is_err());
        assert!(DiffParams { win: 10, step: 10, flank: 0, junction_flank: 0 }.validate().is_ok());
        // Bound the knobs so `2 * flank` cannot overflow.
        let m = DiffParams::MAX_T2;
        assert!(DiffParams { win: m, step: 1, flank: m, junction_flank: m }.validate().is_ok());
        assert!(DiffParams { win: m + 1, step: 1, flank: 0, junction_flank: 0 }.validate().is_err());
        assert!(DiffParams { win: 1, step: 1, flank: m + 1, junction_flank: 0 }.validate().is_err());
        assert!(DiffParams { win: 1, step: 1, flank: 0, junction_flank: m + 1 }.validate().is_err());
        // junction_flank 0 is allowed (junction taint disabled).
        assert!(DiffParams { win: 1, step: 1, flank: 0, junction_flank: 0 }.validate().is_ok());
    }

    #[test]
    fn whole_contig_windows_cover_the_seam() {
        // A wholly novel contig of exactly one window: without circular seam
        // tiling the single window [0, len) misses matches across the origin.
        let params = DiffParams { win: 10_000, step: 2_000, flank: 10_000, junction_flank: 50 };
        let len = 10_000u64;
        let ws = interval_windows("c", len, 0, len, params);
        assert!(!ws.is_empty());
        let mut covered = vec![false; len as usize];
        for w in &ws {
            for p in w.new_start..w.new_end {
                covered[(p % len) as usize] = true;
            }
        }
        assert!(covered.iter().all(|&b| b), "windows must cover the whole contig");
        assert!(
            ws.iter().any(|w| w.new_start <= len - 100 && w.new_end >= len + 100),
            "a window must span the origin seam ± 100 bp: {ws:?}"
        );

        // A contig shorter than `win` must also get a seam-spanning window.
        let len = 3_000u64;
        let ws = interval_windows("c", len, 0, len, params);
        assert!(
            ws.iter().any(|w| w.new_start <= len - 100 && w.new_end >= len + 100),
            "short wholly novel contig must span the seam: {ws:?}"
        );
    }

    #[test]
    fn t1_sweep_zero_one_two() {
        // empty -> whole contig
        assert_eq!(regions_not_exactly_once(vec![], 100), vec![(0, 100)]);
        // one interval -> the uncovered flanks
        assert_eq!(regions_not_exactly_once(vec![(10, 20)], 100), vec![(0, 10), (20, 100)]);
        // two overlapping -> uncovered flanks and the double-covered middle
        assert_eq!(
            regions_not_exactly_once(vec![(10, 30), (20, 40)], 100),
            vec![(0, 10), (20, 30), (40, 100)]
        );
        // two adjacent -> covered exactly once in [10,30), flanks uncovered
        assert_eq!(
            regions_not_exactly_once(vec![(10, 20), (20, 30)], 100),
            vec![(0, 10), (30, 100)]
        );
        // two identical -> double coverage everywhere combined with the flanks
        assert_eq!(
            regions_not_exactly_once(vec![(10, 20), (10, 20)], 100),
            vec![(0, 100)]
        );
        // full cover -> nothing
        assert_eq!(regions_not_exactly_once(vec![(0, 100)], 100), vec![]);
    }

    #[test]
    fn merge_adjacent_and_overlapping() {
        let m = merge_old(vec![OldIv::new(0, 5), OldIv::new(5, 10), OldIv::new(20, 30)]);
        assert_eq!(m, vec![OldIv::new(0, 10), OldIv::new(20, 30)]);
    }

    #[test]
    fn novel_computation() {
        // one block covering [10,20): everything else is novel.
        let blocks = vec![(10u64, 20u64)];
        assert_eq!(regions_not_exactly_once(blocks, 100), vec![(0, 10), (20, 100)]);
        // two blocks covering the same new region -> whole contig is novel.
        let blocks = vec![(10u64, 20u64), (10u64, 20u64)];
        assert_eq!(regions_not_exactly_once(blocks, 100), vec![(0, 100)]);
    }

    fn blk(oc: &str, os: u64, oe: u64, strand: Strand, ns: u64) -> Block {
        Block {
            old_contig: oc.into(),
            old: OldIv::new(os, oe),
            new_contig: "n".into(),
            new_start: NewPos(ns),
            strand,
        }
    }

    #[test]
    fn junction_kind_classifies() {
        let p = Strand::Plus;
        let m = Strand::Minus;
        // Natural plus adjacency: old end == next old start.
        assert_eq!(junction_kind(&blk("c", 0, 100, p, 0), &blk("c", 100, 200, p, 100), 1000), JunctionKind::Natural);
        // Natural minus adjacency: old start == next old end.
        assert_eq!(junction_kind(&blk("c", 100, 200, m, 0), &blk("c", 0, 100, m, 100), 1000), JunctionKind::Natural);
        // Deletion inside a colinear chain: old gap.
        assert_eq!(junction_kind(&blk("c", 0, 100, p, 0), &blk("c", 200, 300, p, 100), 1000), JunctionKind::Edit);
        // Contig join: different old contigs.
        assert_eq!(junction_kind(&blk("c1", 0, 100, p, 0), &blk("c2", 0, 100, p, 100), 1000), JunctionKind::Edit);
        // Strand flip: inversion junction.
        assert_eq!(junction_kind(&blk("c", 0, 100, p, 0), &blk("c", 100, 200, m, 100), 1000), JunctionKind::Edit);
        // Circular-adjacent wrap: a junction like any other (minimap2 indexes
        // the old contig linearly, so the copy can still become contiguous).
        assert_eq!(junction_kind(&blk("c", 900, 1000, p, 0), &blk("c", 0, 100, p, 100), 1000), JunctionKind::RotationWrap);
        assert_eq!(junction_kind(&blk("c", 0, 100, m, 0), &blk("c", 900, 1000, m, 100), 1000), JunctionKind::RotationWrap);
        // A wrap needs an exact circular adjacency: a gap at the origin is an edit.
        assert_eq!(junction_kind(&blk("c", 900, 1000, p, 0), &blk("c", 5, 100, p, 100), 1000), JunctionKind::Edit);
        assert_eq!(junction_kind(&blk("c", 800, 990, p, 0), &blk("c", 0, 100, p, 100), 1000), JunctionKind::Edit);
    }

    fn nb(oc: &str, os: u64, oe: u64, strand: Strand, new_contig: &str, ns: u64) -> Block {
        Block {
            old_contig: oc.into(),
            old: OldIv::new(os, oe),
            new_contig: new_contig.into(),
            new_start: NewPos(ns),
            strand,
        }
    }

    fn write_fasta_seqs(dir: &std::path::Path, file: &str, seqs: &[(&str, usize)]) -> Fasta {
        let mut s = String::new();
        for (name, len) in seqs {
            s.push('>');
            s.push_str(name);
            s.push('\n');
            s.push_str(&"A".repeat(*len));
            s.push('\n');
        }
        let p = dir.join(file);
        std::fs::write(&p, s).unwrap();
        Fasta::load(&p).unwrap()
    }

    /// A block-end event based junction search must find a circular wrap, must
    /// not be hidden by a nested block, and must find a new-contig seam.
    #[test]
    fn compute_junctions_wrap_nested_and_seam() {
        let dir = std::env::temp_dir().join("cramhook_junction_test");
        std::fs::create_dir_all(&dir).unwrap();
        let p = Strand::Plus;

        // 1. Circular wrap, old contig "c" len 1000 (header not needed).
        let old = write_fasta_seqs(&dir, "old.fa", &[("c", 1000)]);
        let new = write_fasta_seqs(&dir, "new.fa", &[("n", 1000)]);
        let mut blocks: HashMap<String, Vec<Block>> = HashMap::new();
        blocks.insert(
            "c".into(),
            vec![nb("c", 500, 1000, p, "n", 0), nb("c", 0, 500, p, "n", 500)],
        );
        let (j, s) = compute_junctions(&blocks, &old, &new, &HashMap::new(), 50);
        assert_eq!((s.count, s.rotation_wraps), (1, 1));
        assert_eq!(j.get("n").map(Vec::len), Some(1));

        // 2. Nested blocks: [0,20000) old c[20000,40000); nested [1000,2000) old
        //    d[0,1000); [20000,40000) old c[0,20000). The c-wrap at 20000 must be
        //    found; the nested end at 2000 falls in a multiply-covered (novel)
        //    region and is dropped.
        let old = write_fasta_seqs(&dir, "old2.fa", &[("c", 40000), ("d", 1000)]);
        let new = write_fasta_seqs(&dir, "new2.fa", &[("n", 40000)]);
        let mut blocks: HashMap<String, Vec<Block>> = HashMap::new();
        blocks.insert(
            "c".into(),
            vec![
                nb("c", 20000, 40000, p, "n", 0),
                nb("c", 0, 20000, p, "n", 20000),
            ],
        );
        blocks.insert("d".into(), vec![nb("d", 0, 1000, p, "n", 1000)]);
        let mut novel: HashMap<String, Vec<NewIv>> = HashMap::new();
        novel.insert("n".into(), vec![NewIv::new(1000, 2000)]);
        let (j, s) = compute_junctions(&blocks, &old, &new, &novel, 50);
        assert_eq!((s.count, s.rotation_wraps), (1, 1));
        assert_eq!(j["n"], vec![NewIv::new(19950, 20050)]);

        // 3. New-contig seam: old c 40000, new n = old[20000,40000)+old[0,10000)
        //    (len 30000). Junctions: the wrap at 20000 and the seam at 0 (the
        //    deletion old[10000,20000) across the new origin).
        let old = write_fasta_seqs(&dir, "old3.fa", &[("c", 40000)]);
        let new = write_fasta_seqs(&dir, "new3.fa", &[("n", 30000)]);
        let mut blocks: HashMap<String, Vec<Block>> = HashMap::new();
        blocks.insert(
            "c".into(),
            vec![nb("c", 20000, 40000, p, "n", 0), nb("c", 0, 10000, p, "n", 20000)],
        );
        let (j, s) = compute_junctions(&blocks, &old, &new, &HashMap::new(), 50);
        assert_eq!((s.count, s.rotation_wraps), (2, 1));
        // j is the merged set of {0 +- 50} and {20000 +- 50}.
        assert_eq!(
            j["n"],
            vec![NewIv::new(0, 50), NewIv::new(19950, 20050), NewIv::new(29950, 30000)]
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Linear-time junction detection: with every boundary a junction, an
    /// `O(B²)` scan took ~8 s on an 8 Mb contig (reviewer). 80 000 blocks must
    /// finish well inside a generous bound.
    #[test]
    fn compute_junctions_scales_linearly() {
        let b: usize = 80_000;
        let dir = std::env::temp_dir().join("cramhook_junction_perf");
        std::fs::create_dir_all(&dir).unwrap();
        let new_fa = write_fasta_seqs(&dir, "perf-new.fa", &[("n", b * 100)]);
        let old_fa = write_fasta_seqs(&dir, "perf-old.fa", &[("c", b * 150 + 1000)]);
        let mut v: Vec<Block> = Vec::with_capacity(b);
        for i in 0..b as u64 {
            // old blocks are 100 bp with a 50 bp gap (stride 150), new stride is
            // 100 -> every boundary is a non-adjacent (edit) junction.
            v.push(Block {
                old_contig: "c".into(),
                old: OldIv::new(i * 150, i * 150 + 100),
                new_contig: "n".into(),
                new_start: NewPos(i * 100),
                strand: Strand::Plus,
            });
        }
        let mut blocks: HashMap<String, Vec<Block>> = HashMap::new();
        blocks.insert("c".into(), v);
        let t = std::time::Instant::now();
        let (j, s) = compute_junctions(&blocks, &old_fa, &new_fa, &HashMap::new(), 50);
        let dt = t.elapsed();
        assert_eq!(s.count, b, "b-1 internal junctions + the circular seam");
        assert!(!j.is_empty());
        assert!(dt.as_secs_f64() < 30.0, "junction detection too slow: {dt:?}");
        // flank 0 disables the whole pass.
        let (_j0, s0) = compute_junctions(&blocks, &old_fa, &new_fa, &HashMap::new(), 0);
        assert_eq!(s0.count, 0);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn push_junction_normalizes() {
        let mut v = Vec::new();
        push_junction(&mut v, 100, 50, 1000);
        assert_eq!(v, vec![NewIv::new(50, 150)]);
        v.clear();
        // Left underflow wraps modulo len (was clamped to 0 -> missed [960,1000)).
        push_junction(&mut v, 10, 50, 1000);
        assert_eq!(v, vec![NewIv::new(960, 1000), NewIv::new(0, 60)]);
        v.clear();
        // Seam exactly at the contig origin.
        push_junction(&mut v, 0, 50, 1000);
        assert_eq!(v, vec![NewIv::new(950, 1000), NewIv::new(0, 50)]);
        v.clear();
        // Right overflow wraps.
        push_junction(&mut v, 990, 50, 1000);
        assert_eq!(v, vec![NewIv::new(940, 1000), NewIv::new(0, 40)]);
        v.clear();
        // 2*flank >= len covers the whole contig (was an interval > len).
        push_junction(&mut v, 20000, 200000, 40000);
        assert_eq!(v, vec![NewIv::new(0, 40000)]);
        v.clear();
        push_junction(&mut v, 100, 500, 1000);
        assert_eq!(v, vec![NewIv::new(0, 1000)]);
        v.clear();
        push_junction(&mut v, 100, 0, 1000); // flank 0 disables
        assert!(v.is_empty());
    }

    #[test]
    fn map_pos_both_strands() {
        let plus = Block {
            old_contig: "c".into(),
            old: OldIv::new(10, 20),
            new_contig: "c".into(),
            new_start: NewPos(100),
            strand: Strand::Plus,
        };
        assert_eq!(plus.map_pos(OldPos(10)), Some(NewPos(100)));
        assert_eq!(plus.map_pos(OldPos(19)), Some(NewPos(109)));
        assert_eq!(plus.map_pos(OldPos(20)), None);
        assert_eq!(plus.map_iv(OldIv::new(12, 15)), Some(NewIv::new(102, 105)));

        let minus = Block {
            old_contig: "c".into(),
            old: OldIv::new(10, 20),
            new_contig: "c".into(),
            new_start: NewPos(100),
            strand: Strand::Minus,
        };
        assert_eq!(minus.map_pos(OldPos(10)), Some(NewPos(109)));
        assert_eq!(minus.map_pos(OldPos(19)), Some(NewPos(100)));
        assert_eq!(minus.map_iv(OldIv::new(12, 15)), Some(NewIv::new(105, 108)));
        assert_eq!(minus.map_iv(OldIv::new(10, 20)), Some(NewIv::new(100, 110)));
    }

    #[test]
    fn window_tiling_edges() {
        // shorter than a window -> one window = whole interval
        assert_eq!(tile_windows(0, 5_000, WIN, STEP), vec![(0, 5_000)]);
        assert_eq!(tile_windows(0, 10_000, WIN, STEP), vec![(0, 10_000)]);
        // longer: last window ends at the interval end
        let t = tile_windows(0, 25_000, WIN, STEP);
        assert_eq!(t.first(), Some(&(0, 10_000)));
        assert_eq!(t.last(), Some(&(15_000, 25_000)));
        for w in &t {
            assert_eq!(w.1 - w.0, WIN);
        }
        // small overhang
        assert_eq!(tile_windows(0, 10_001, WIN, STEP), vec![(0, 10_000), (1, 10_001)]);
    }

    #[test]
    fn block_containing_and_taint() {
        let mut diff = Diff::default();
        diff.blocks.insert(
            "c".into(),
            vec![
                Block {
                    old_contig: "c".into(),
                    old: OldIv::new(0, 100),
                    new_contig: "c".into(),
                    new_start: NewPos(0),
                    strand: Strand::Plus,
                },
                Block {
                    old_contig: "c".into(),
                    old: OldIv::new(200, 300),
                    new_contig: "c".into(),
                    new_start: NewPos(200),
                    strand: Strand::Plus,
                },
            ],
        );
        assert!(diff.block_containing("c", OldIv::new(10, 20)).is_some());
        assert!(diff.block_containing("c", OldIv::new(50, 250)).is_none());
        assert!(diff.block_containing("c", OldIv::new(150, 160)).is_none());

        diff.taint_t2.insert("c".into(), vec![OldIv::new(40, 60)]);
        assert!(diff.overlaps_taint("c", OldIv::new(59, 61)));
        assert!(!diff.overlaps_taint("c", OldIv::new(60, 70)));
        assert!(!diff.overlaps_taint("d", OldIv::new(0, 10)));
    }
}

#[cfg(test)]
mod t2_tests {
    use super::*;

    fn hit(cs: &str, qstart: u64, qend: u64, strand: Strand, tstart: u64, tend: u64) -> PafRecord {
        PafRecord {
            qname: "w0".into(),
            qlen: 1_000,
            qstart,
            qend,
            strand,
            tname: "o".into(),
            tlen: 10_000,
            tstart,
            tend,
            mapq: 60,
            primary: true,
            cs: cs.into(),
        }
    }

    /// Project a single hand-made hit through [`project_t2_hits`] and return the
    /// merged old intervals it taints. All windows start at absolute coordinate
    /// `win_start` (the `w0` window).
    fn proj(rec: PafRecord, win_start: u64, novel: &[(u64, u64)]) -> Vec<OldIv> {
        let windows = vec![Window {
            contig: "n".into(),
            new_start: win_start,
            new_end: win_start + 1_000,
            contig_len: 10_000,
        }];
        let mut novel_map: HashMap<String, Vec<NewIv>> = HashMap::new();
        novel_map.insert("n".into(), novel.iter().map(|&(a, b)| NewIv::new(a, b)).collect());
        let junctions: HashMap<String, Vec<NewIv>> = HashMap::new();
        let blocks: HashMap<String, Vec<Block>> = HashMap::new();
        let t = T2Targets { novel: &novel_map, junctions: &junctions, blocks_by_new: &blocks };
        let mut raw: Vec<(String, OldIv)> = Vec::new();
        let mut raw_j: Vec<(String, OldIv)> = Vec::new();
        project_t2_hits(&[rec], &windows, &t, &mut raw, &mut raw_j).unwrap();
        assert!(raw_j.is_empty(), "no junctions passed, so no junction taint");
        merge_old(raw.into_iter().map(|(_, iv)| iv).collect())
    }

    /// Like [`proj`] but the taint comes from a *junction* interval.
    fn proj_junction(rec: PafRecord, win_start: u64, junctions: &[(u64, u64)]) -> Vec<OldIv> {
        let windows = vec![Window {
            contig: "n".into(),
            new_start: win_start,
            new_end: win_start + 1_000,
            contig_len: 10_000,
        }];
        let novel_map: HashMap<String, Vec<NewIv>> = HashMap::new();
        let mut jmap: HashMap<String, Vec<NewIv>> = HashMap::new();
        jmap.insert("n".into(), junctions.iter().map(|&(a, b)| NewIv::new(a, b)).collect());
        let blocks: HashMap<String, Vec<Block>> = HashMap::new();
        let t = T2Targets { novel: &novel_map, junctions: &jmap, blocks_by_new: &blocks };
        let mut raw: Vec<(String, OldIv)> = Vec::new();
        let mut raw_j: Vec<(String, OldIv)> = Vec::new();
        project_t2_hits(&[rec], &windows, &t, &mut raw, &mut raw_j).unwrap();
        assert!(raw.is_empty(), "no novel passed, so no novel taint");
        merge_old(raw_j.into_iter().map(|(_, iv)| iv).collect())
    }

    /// Like [`proj_junction`] but with identity blocks so the self-hit exclusion
    /// is active. `blocks_by_new` maps new contig -> blocks.
    fn proj_junction_blocks(
        rec: PafRecord,
        win_start: u64,
        junctions: &[(u64, u64)],
        blocks: &HashMap<String, Vec<Block>>,
    ) -> Vec<OldIv> {
        let windows = vec![Window {
            contig: "n".into(),
            new_start: win_start,
            new_end: win_start + 1_000,
            contig_len: 10_000,
        }];
        let novel_map: HashMap<String, Vec<NewIv>> = HashMap::new();
        let mut jmap: HashMap<String, Vec<NewIv>> = HashMap::new();
        jmap.insert("n".into(), junctions.iter().map(|&(a, b)| NewIv::new(a, b)).collect());
        let t = T2Targets { novel: &novel_map, junctions: &jmap, blocks_by_new: blocks };
        let mut raw: Vec<(String, OldIv)> = Vec::new();
        let mut raw_j: Vec<(String, OldIv)> = Vec::new();
        project_t2_hits(&[rec], &windows, &t, &mut raw, &mut raw_j).unwrap();
        merge_old(raw_j.into_iter().map(|(_, iv)| iv).collect())
    }

    #[test]
    fn plus_match_only_novel_subset() {
        // 10 matches at window[100..110); only [105,110) is novel.
        let r = hit(":10", 0, 10, Strand::Plus, 500, 510);
        assert_eq!(proj(r, 100, &[(105, 110)]), vec![OldIv::new(505, 510)]);
    }

    #[test]
    fn plus_nonzero_query_and_target_start() {
        // qstart=3 shifts the window coordinate; tstart=700 shifts the old coords.
        let r = hit(":10", 3, 13, Strand::Plus, 700, 710);
        assert_eq!(proj(r, 50, &[(55, 60)]), vec![OldIv::new(702, 707)]);
    }

    #[test]
    fn minus_query_goes_backwards() {
        // strand '-': qoff -> qend-1-qoff, then new = win_start + that.
        let r = hit(":10", 2, 12, Strand::Minus, 200, 210);
        assert_eq!(proj(r, 40, &[(45, 50)]), vec![OldIv::new(202, 207)]);
    }

    #[test]
    fn ins_at_novel_boundary_taints_flank() {
        // Insertion of 3 query bases at target 100+5; new base 6 is novel, so the
        // old bases [t-1, t+1) = [104,106) are tainted.
        let r = hit(":5+aaa:5", 0, 13, Strand::Plus, 100, 110);
        assert_eq!(proj(r, 0, &[(6, 7)]), vec![OldIv::new(104, 106)]);
    }

    #[test]
    fn ins_outside_novel_does_not_taint() {
        // Inserted window bases 5,6,7 are not novel; only a trailing match base
        // is. After the insertion the target cursor lags the query offset by 3,
        // so the matching window base 9 is at old 100+5+1 = 106.
        let r = hit(":5+aaa:5", 0, 13, Strand::Plus, 100, 110);
        let merged = proj(r, 0, &[(9, 10)]);
        assert_eq!(merged, vec![OldIv::new(106, 107)]);
        assert!(!merged.iter().any(|iv| iv.overlaps(OldIv::new(104, 106))));
    }

    #[test]
    fn del_next_to_novel_taints_deleted_bases() {
        // '-aaa' deletes target [100+5,100+8); the adjacent next query base maps
        // to novel window 5, so the deleted old bases are tainted, as is the
        // matching base that shares that window position.
        let r = hit(":5-aaa:5", 0, 13, Strand::Plus, 100, 113);
        assert_eq!(proj(r, 0, &[(5, 6)]), vec![OldIv::new(105, 109)]);
    }

    #[test]
    fn del_not_next_to_novel_does_not_taint() {
        // Novel window 2 is an earlier match base; the deletion neighbours are not
        // novel, so no deleted base is tainted.
        let r = hit(":5-aaa:5", 0, 13, Strand::Plus, 100, 113);
        let merged = proj(r, 0, &[(2, 3)]);
        assert_eq!(merged, vec![OldIv::new(102, 103)]);
        assert!(!merged.iter().any(|iv| iv.overlaps(OldIv::new(105, 108))));
    }

    // --- circular wrap math -------------------------------------------------

    #[test]
    fn circular_seq_wraps_past_end() {
        let seq = b"0123456789";
        assert_eq!(circular_seq(seq, 0, 10), b"0123456789");
        assert_eq!(circular_seq(seq, 8, 5), b"89012");
        assert_eq!(circular_seq(seq, 5, 0), b"");
        // A window longer than the contig repeats it.
        assert_eq!(circular_seq(seq, 0, 12), b"012345678901");
    }

    #[test]
    fn interval_windows_wrap_across_origin() {
        let params = DiffParams { win: 400, step: 200, flank: 100, junction_flank: 50 };
        // Novel [900,1000) on a 1000 bp contig: arc [800, 1100) crosses the end.
        let ws = interval_windows("n", 1000, 900, 1000, params);
        assert_eq!(ws.len(), 1);
        assert_eq!(
            (ws[0].new_start, ws[0].new_end, ws[0].contig_len),
            (800, 1100, 1000)
        );
        // The window reads 800..1000 then 0..100 (300 circular bases).
        assert_eq!(ws[0].new_end - ws[0].new_start, 300);
        assert_eq!(circular_seq(b"0123456789", 8, 5), b"89012");
    }

    #[test]
    fn interval_windows_interior_no_wrap() {
        let params = DiffParams { win: 10_000, step: 2_000, flank: 10_000, junction_flank: 50 };
        let ws = interval_windows("n", 100_000, 50_000, 51_000, params);
        // arc [40_000, 61_000): first window starts at the arc start, last ends
        // exactly at the arc end, and nothing wraps.
        assert_eq!(ws.first().map(|w| (w.new_start, w.new_end)), Some((40_000, 50_000)));
        assert_eq!(ws.last().map(|w| (w.new_start, w.new_end)), Some((51_000, 61_000)));
        assert!(ws.iter().all(|w| w.new_end <= 100_000));
    }

    #[test]
    fn interval_windows_short_contig_shorter_than_window() {
        // The contig (300 bp) is shorter than `win`; a single window covers the
        // whole arc, which itself is capped at the contig length.
        let params = DiffParams { win: 10_000, step: 2_000, flank: 100, junction_flank: 50 };
        let ws = interval_windows("n", 300, 200, 300, params);
        assert_eq!(ws.len(), 1);
        assert_eq!(
            (ws[0].new_start, ws[0].new_end, ws[0].contig_len),
            (100, 400, 300)
        );
        assert_eq!(ws[0].new_end - ws[0].new_start, 300);
    }

    #[test]
    fn projection_wraps_novel_spanning_origin() {
        // Window [800, 1100) on a 1000 bp contig; novel is [0,100)+[800,1000),
        // i.e. the window is novel across the whole origin junction. A 300 bp
        // match over the window must taint every old base (projection modulo).
        let rec = hit(":300", 0, 300, Strand::Plus, 0, 300);
        let windows = vec![Window {
            contig: "n".into(),
            new_start: 800,
            new_end: 1100,
            contig_len: 1000,
        }];
        let mut novel_map: HashMap<String, Vec<NewIv>> = HashMap::new();
        novel_map.insert("n".into(), vec![NewIv::new(0, 100), NewIv::new(800, 1000)]);
        let junctions: HashMap<String, Vec<NewIv>> = HashMap::new();
        let blocks: HashMap<String, Vec<Block>> = HashMap::new();
        let t = T2Targets { novel: &novel_map, junctions: &junctions, blocks_by_new: &blocks };
        let mut raw: Vec<(String, OldIv)> = Vec::new();
        let mut raw_j: Vec<(String, OldIv)> = Vec::new();
        project_t2_hits(&[rec], &windows, &t, &mut raw, &mut raw_j).unwrap();
        let merged = merge_old(raw.into_iter().map(|(_, iv)| iv).collect());
        assert_eq!(merged, vec![OldIv::new(0, 300)]);
    }

    #[test]
    fn junction_membership_taints_aligned_old_base() {
        // A window base inside a junction interval taints the old base aligned
        // to it, through the junction (not novel) path.
        let r = hit(":10", 0, 10, Strand::Plus, 500, 510);
        assert_eq!(proj_junction(r, 100, &[(105, 110)]), vec![OldIv::new(505, 510)]);
        // No junction -> no taint.
        let r = hit(":10", 0, 10, Strand::Plus, 500, 510);
        assert!(proj_junction(r, 100, &[]).is_empty());
    }

    #[test]
    fn novel_wins_over_junction_at_the_same_position() {
        // When a window base is both novel and junction, it is counted as novel
        // so the existing `taint_t2` behaviour is preserved.
        let windows = vec![Window {
            contig: "n".into(),
            new_start: 0,
            new_end: 1_000,
            contig_len: 10_000,
        }];
        let mut novel_map: HashMap<String, Vec<NewIv>> = HashMap::new();
        novel_map.insert("n".into(), vec![NewIv::new(5, 6)]);
        let mut jmap: HashMap<String, Vec<NewIv>> = HashMap::new();
        jmap.insert("n".into(), vec![NewIv::new(0, 100)]);
        let rec = hit(":10", 0, 10, Strand::Plus, 500, 510);
        let blocks: HashMap<String, Vec<Block>> = HashMap::new();
        let t = T2Targets { novel: &novel_map, junctions: &jmap, blocks_by_new: &blocks };
        let mut raw: Vec<(String, OldIv)> = Vec::new();
        let mut raw_j: Vec<(String, OldIv)> = Vec::new();
        project_t2_hits(&[rec], &windows, &t, &mut raw, &mut raw_j).unwrap();
        assert_eq!(merge_old(raw.into_iter().map(|(_, iv)| iv).collect()), vec![OldIv::new(505, 506)]);
        assert_eq!(
            merge_old(raw_j.into_iter().map(|(_, iv)| iv).collect()),
            vec![OldIv::new(500, 505), OldIv::new(506, 510)]
        );
    }

    #[test]
    fn junction_self_hit_is_not_tainted() {
        // A junction window whose hit maps back to the *source* of the junction
        // (same old contig, position and strand as p's identity block) is a
        // self-hit and is ignored; a hit to a competitor copy is tainted.
        let mut jmap: HashMap<String, Vec<NewIv>> = HashMap::new();
        jmap.insert("n".into(), vec![NewIv::new(105, 110)]);
        // The block covering window positions [0,1000) maps new[100..110) to
        // old "o"[500..510) on '+': the window's own source.
        let mut blocks: HashMap<String, Vec<Block>> = HashMap::new();
        blocks.insert(
            "n".into(),
            vec![Block {
                old_contig: "o".into(),
                old: OldIv::new(500, 600),
                new_contig: "n".into(),
                new_start: NewPos(100),
                strand: Strand::Plus,
            }],
        );
        // Self-hit: window -> old "o" at 500 (matches the block's source).
        let r = hit(":10", 0, 10, Strand::Plus, 500, 510);
        assert!(
            proj_junction_blocks(r, 100, &[(105, 110)], &blocks).is_empty(),
            "a self-hit must not taint"
        );
        // Competitor: the same window maps to old "p" -> tainted.
        let r = hit(":10", 0, 10, Strand::Plus, 500, 510);
        let r = PafRecord { tname: "p".into(), ..r };
        assert_eq!(
            proj_junction_blocks(r, 100, &[(105, 110)], &blocks),
            vec![OldIv::new(505, 510)]
        );
        // Same position but the block maps on '-': a '+' hit is not the source.
        let mut minus_blocks: HashMap<String, Vec<Block>> = HashMap::new();
        minus_blocks.insert(
            "n".into(),
            vec![Block {
                old_contig: "o".into(),
                old: OldIv::new(500, 600),
                new_contig: "n".into(),
                new_start: NewPos(100),
                strand: Strand::Minus,
            }],
        );
        let r = hit(":10", 0, 10, Strand::Plus, 500, 510);
        assert!(
            !proj_junction_blocks(r, 100, &[(105, 110)], &minus_blocks).is_empty(),
            "a '+' hit is not the source of a '-' block"
        );
    }

    /// `is_self_hit` isolates contig, position **and strand**, and does not
    /// accept a hit shifted by even 1 bp.
    #[test]
    fn is_self_hit_isolates_contig_pos_and_strand() {
        let mut blocks: HashMap<String, Vec<Block>> = HashMap::new();
        blocks.insert(
            "n".into(),
            vec![
                Block {
                    old_contig: "o".into(),
                    old: OldIv::new(500, 600),
                    new_contig: "n".into(),
                    new_start: NewPos(100),
                    strand: Strand::Plus,
                },
                Block {
                    old_contig: "o2".into(),
                    old: OldIv::new(500, 600),
                    new_contig: "n".into(),
                    new_start: NewPos(1000),
                    strand: Strand::Minus,
                },
            ],
        );
        // '+' block: new 105 -> source ("o", 505, +).
        assert!(is_self_hit(&blocks, "n", 105, "o", Strand::Plus, 505));
        assert!(
            !is_self_hit(&blocks, "n", 105, "o", Strand::Minus, 505),
            "same contig/pos but opposite strand must not be a self-hit"
        );
        assert!(!is_self_hit(&blocks, "n", 105, "p", Strand::Plus, 505));
        assert!(!is_self_hit(&blocks, "n", 105, "o", Strand::Plus, 506));
        assert!(!is_self_hit(&blocks, "n", 105, "o", Strand::Plus, 504));
        // '-' block: new 1005 -> source ("o2", 594, -).
        assert!(is_self_hit(&blocks, "n", 1005, "o2", Strand::Minus, 594));
        assert!(
            !is_self_hit(&blocks, "n", 1005, "o2", Strand::Plus, 594),
            "same contig/pos but opposite strand must not be a self-hit"
        );
        assert!(!is_self_hit(&blocks, "n", 1005, "o2", Strand::Minus, 595));
        assert!(!is_self_hit(&blocks, "n", 1005, "o2", Strand::Minus, 593));
        // Outside any block there is no source.
        assert!(!is_self_hit(&blocks, "n", 50, "o", Strand::Plus, 505));
    }
}
