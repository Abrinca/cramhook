//! M3: contract check, eligibility, translation and origin-crossing split reads.
//!
//! A single streaming pass over the input CRAM classifies every read as
//! *translated* (record(s) written to `P.translated.bam`) or *deferred* (one
//! FASTQ entry per read, built from its primary record). The diff computed by
//! M2 supplies identity blocks and taint regions.
//!
//! Coordinates are 0-based half-open, matching the rest of the crate. Strand is
//! the record's SAM orientation (`0x10`). See `docs/plans/Plan-v1.md` sections 4.2/4.3.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use serde::Serialize;

use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, Format, Read};

use crate::bamio::aux_int;
use crate::contract::{
    check_contract, check_sq_against_fai, parse_header, read_fai, reject_old_name_lift,
};
use crate::coords::{OldIv, OldPos};
use crate::diff::{Block, Diff, Preset};
use crate::fasta::{Fasta, ensure_fai};
use crate::paf::Strand;
use crate::version::versioned;

/// Tags a translatable record may carry (Plan section 2). `MD` is handled
/// separately (decode artifact, see `docs/decisions.md`) and `s2` is read by
/// rule 3b but is not valid output; both are stripped from translated records by
/// [`strip_translated_tags`].
///
/// `XL` is cramhook's own provenance tag (`XL:A:T` translated / `XL:A:R`
/// remapped). It is accepted as input so a CRAM produced by `cramhook lift` can
/// be lifted again (a chain); `XL:A:T` also marks a record whose `s2` was
/// deliberately stripped by an earlier lift (see rule 3b in [`check_record`]).
pub const ALLOWED_TAGS: &[&[u8]] =
    &[b"NM", b"ms", b"AS", b"nn", b"tp", b"cm", b"s1", b"s2", b"de", b"rl", b"SA", b"XL"];

/// Divergence (`de`) threshold per preset (Plan 4.2).
pub fn de_threshold(preset: Preset) -> f64 {
    match preset {
        Preset::MapHifi => 0.01,
        Preset::LrHq => 0.03,
    }
}

/// minimap2's default Z-drop threshold (`-z`), used for both read presets.
pub const PRESET_ZDROP: i64 = 400;

/// Default `--max-s2-ratio` (the CLI default and `TranslateConfig::default`).
/// Also the fallback when a `cramhook lift` `@PG` `CL` does not spell the flag,
/// used to decide whether an inherited `XL:A:T` record's stripped `s2` is safe.
pub const DEFAULT_MAX_S2_RATIO: f64 = 0.5;

/// ksw2 scoring parameters (`ref/minimap2/options.c`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Scoring {
    pub a: i64,
    pub b: i64,
    pub q: i64,
    pub e: i64,
    pub q2: i64,
    pub e2: i64,
}

impl Default for Scoring {
    /// The lr:hq / map-ont defaults (also `mm_mapopt_init`).
    fn default() -> Self {
        Self { a: 2, b: 4, q: 4, e: 2, q2: 24, e2: 1 }
    }
}

/// The ksw2 scoring of a read preset (`ref/minimap2/options.c`).
pub fn scoring(preset: Preset) -> Scoring {
    match preset {
        Preset::LrHq => Scoring::default(),
        Preset::MapHifi => Scoring { a: 1, b: 4, q: 6, e: 2, q2: 26, e2: 1 },
    }
}

/// Cost of a single I/D op of length `l` (dual gap penalty), applied once.
#[inline]
fn gap_cost(l: i64, sc: &Scoring) -> i64 {
    (sc.q + l * sc.e).min(sc.q2 + l * sc.e2)
}

/// Read a run of decimal `MD` digits, advancing `*i`. `None` on overflow.
/// A single-digit run (by far the common case, especially `0` between
/// consecutive mismatches) avoids the overflow-checked multiply.
#[inline]
fn md_read_i64(b: &[u8], i: &mut usize) -> Option<i64> {
    let mut n: i64 = (b[*i] - b'0') as i64;
    *i += 1;
    if *i >= b.len() || !b[*i].is_ascii_digit() {
        return Some(n);
    }
    while *i < b.len() && b[*i].is_ascii_digit() {
        n = n.checked_mul(10)?.checked_add((b[*i] - b'0') as i64)?;
        *i += 1;
    }
    Some(n)
}

/// Running state of the ksw2 Z-drop recurrence.
struct ZdState {
    s: i64,
    runmax: i64,
    argmax_diag: i64,
    diag: i64,
    drop_max: i64,
}

impl ZdState {
    fn new() -> Self {
        Self { s: 0, runmax: i64::MIN, argmax_diag: 0, diag: 0, drop_max: 0 }
    }

    /// One alignment column: `score` on the path, `ddiag` on the
    /// query-minus-ref diagonal.
    fn column(&mut self, score: i64, ddiag: i64, e2: i64) {
        self.s += score;
        self.diag += ddiag;
        if self.s >= self.runmax {
            self.runmax = self.s;
            self.argmax_diag = self.diag;
        }
        let drop = self.runmax - self.s - e2 * (self.diag - self.argmax_diag).abs();
        if drop > self.drop_max {
            self.drop_max = drop;
        }
    }

    /// A run of `n >= 1` match columns. `S` rises monotonically and the diagonal
    /// is unchanged, so only the first column can raise `drop_max`: every later
    /// column is closer to the running max, and once `S` reaches it the drop is 0.
    fn matches(&mut self, n: i64, sc: &Scoring) {
        if self.runmax == i64::MIN {
            // No column processed yet: the run becomes the running max, drop 0.
            self.s += n * sc.a;
            self.runmax = self.s;
            self.argmax_diag = self.diag;
            return;
        }
        let s_first = self.s + sc.a;
        if s_first < self.runmax {
            let drop = self.runmax - s_first - sc.e2 * (self.diag - self.argmax_diag).abs();
            if drop > self.drop_max {
                self.drop_max = drop;
            }
        }
        self.s += n * sc.a;
        if self.s >= self.runmax {
            self.runmax = self.s;
            self.argmax_diag = self.diag;
        }
    }
}

/// The path Z-drop statistic of an old alignment, mirroring ksw2's
/// `ksw_apply_zdrop` (`ref/minimap2/ksw2.h`, called with `e2` from
/// `ksw2_extd2_sse.c`). Walks the CIGAR left to right: `M`/`=`/`X` columns score
/// `+a`/`-b` from the `MD` string, each `I`/`D` op scores
/// `-min(q + l·e, q2 + l·e2)` once, and soft clips are skipped. The statistic is
/// the maximum of `Smax - S - e2·|diag - diag_at_max|`; `0` for an empty path.
///
/// Consecutive matches are batched: a match run only raises `S`, so the single
/// column that can move `drop_max` is its first. The cost is therefore
/// O(mismatches + indels + MD runs), not O(aligned bases), and the function
/// allocates nothing.
///
/// `None` means the `MD` tag is missing, malformed (strict grammar
/// `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*`), or inconsistent with the CIGAR (the
/// caller defers fail-safe). An `MD` mismatch involving an `N` is scored `-b`,
/// not minimap2's `-sc_ambi` (a manager decision): this is conservative — it can
/// only increase the statistic, and so only defer more, and only on `-` blocks.
pub fn zdrop_statistic(cigar: &[Cigar], md: Option<&str>, sc: &Scoring) -> Option<i64> {
    let b = md?.as_bytes();
    if b.is_empty() {
        return None;
    }
    let mut i = 0usize;
    // MD grammar state: a numeric run is due (`true`) or a letter/`^` is due.
    let mut expect_num = true;
    // Matches remaining in the current numeric run.
    let mut run: i64 = 0;

    let mut zd = ZdState::new();
    for c in cigar {
        match *c {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                // `is_x`: every column must be an `MD` mismatch (`X`); `is_eq`:
                // every column must be a match (`=`); both false: `M`, free.
                let (is_eq, is_x) = match *c {
                    Cigar::Equal(_) => (true, false),
                    Cigar::Diff(_) => (false, true),
                    _ => (false, false),
                };
                let mut rem = l as i64;
                while rem > 0 {
                    if run > 0 {
                        if is_x {
                            return None; // a match run inside an `X` op
                        }
                        let take = run.min(rem);
                        zd.matches(take, sc);
                        run -= take;
                        rem -= take;
                    } else if expect_num {
                        if i >= b.len() || !b[i].is_ascii_digit() {
                            return None;
                        }
                        let n = md_read_i64(b, &mut i)?;
                        if is_x && n > 0 {
                            return None; // matches inside an `X` op
                        }
                        expect_num = false;
                        // Consume the matches now, so a match run costs one
                        // iteration instead of two.
                        let take = n.min(rem);
                        if take > 0 {
                            zd.matches(take, sc);
                            rem -= take;
                        }
                        run = n - take;
                    } else {
                        if i >= b.len() {
                            return None;
                        }
                        let ch = b[i];
                        if ch == b'^' {
                            return None; // deletion where an aligned column is expected
                        }
                        if !ch.is_ascii_uppercase() {
                            return None;
                        }
                        i += 1;
                        if is_eq {
                            return None; // a mismatch inside an `=` op
                        }
                        zd.column(-sc.b, 0, sc.e2);
                        rem -= 1;
                        // A number must follow a mismatch letter; consume it here
                        // so a run of mismatches costs one iteration each. Its
                        // matches belong to this or a later op, so the `X` check
                        // is left to the `run > 0` branch (which knows the op).
                        if i < b.len() && b[i].is_ascii_digit() {
                            run = md_read_i64(b, &mut i)?;
                            expect_num = false;
                        } else {
                            expect_num = true; // trailing letter; finalize rejects
                        }
                    }
                }
            }
            Cigar::Ins(l) => {
                // An insertion consumes no ref, so a match run simply continues
                // through it and no `MD` token is consumed.
                let l = l as i64;
                zd.column(-gap_cost(l, sc), l, sc.e2);
            }
            Cigar::Del(l) => {
                let l = l as i64;
                // Every match before the deletion must already be consumed, so a
                // pending numeric run must be present and zero.
                if expect_num {
                    if i >= b.len() || !b[i].is_ascii_digit() {
                        return None;
                    }
                    if md_read_i64(b, &mut i)? != 0 {
                        return None; // matches pending before a deletion
                    }
                }
                if run > 0 || i >= b.len() || b[i] != b'^' {
                    return None;
                }
                i += 1;
                let start = i;
                while i < b.len() && b[i].is_ascii_uppercase() {
                    i += 1;
                }
                if (i - start) as i64 != l {
                    return None; // `^` run must match the `D` op length
                }
                expect_num = true;
                zd.column(-gap_cost(l, sc), -l, sc.e2);
            }
            // Clips and ref skips consume no column and carry no score.
            Cigar::SoftClip(_)
            | Cigar::HardClip(_)
            | Cigar::Pad(_)
            | Cigar::RefSkip(_) => {}
        }
    }
    // The `MD` must end in a numeric run; a trailing zero run may remain.
    if expect_num {
        if i >= b.len() || !b[i].is_ascii_digit() {
            return None;
        }
        if md_read_i64(b, &mut i)? != 0 {
            return None; // matches with no aligned column to consume them
        }
    }
    if run > 0 || i != b.len() {
        return None;
    }
    Some(zd.drop_max)
}

/// Print a prominent warning when the clip screen is disabled.
pub(crate) fn warn_if_unsafe_clip_screen(cfg: &TranslateConfig) {
    if cfg.unsafe_no_clip_screen {
        crate::note!(
            "WARNING: clip screen DISABLED (--unsafe-no-clip-screen): chimeric reads may be \
             translated WITHOUT their supplementary alignment; the output is NOT fail-safe."
        );
    }
}

/// CLI options for `cramhook translate`.
#[derive(Clone, Copy, Debug)]
pub struct TranslateConfig {
    pub strict: bool,
    /// Tag translated records `XL:A:T` (and remapped records `XL:A:R`). On by
    /// default; set false with `--no-xl-tag`.
    pub xl_tag: bool,
    pub min_mapq: u8,
    /// Maximum soft-clip; `0` means unlimited.
    pub max_clip: u64,
    pub margin: u64,
    /// Divergence (`de`) threshold override; `None` uses the preset default.
    pub max_de: Option<f64>,
    pub t2_win: u64,
    pub t2_step: u64,
    pub t2_flank: u64,
    /// Half-width of the T2 junction interval added at every new adjacency
    /// (`--t2-junction-flank`); 0 disables junction taint.
    pub t2_junction_flank: u64,
    /// Clip screen: k-mer length, 1..=32.
    pub clip_k: u32,
    /// Clip screen: minimum shared k-mers to call a hit.
    pub clip_min_hits: u32,
    /// Clip screen: minimum clip length to screen.
    pub clip_min_len: u64,
    /// Clip screen: extra bp added around novel intervals / edit sites (default
    /// 0: only actual novel bases and boundary-spanning k-mers count).
    pub clip_flank: u64,
    /// Clip screen (debug): also drop k-mers that occur anywhere in old.fa
    /// (Variant B').
    pub clip_subtract_old: bool,
    /// Disable the clip screen (fail-unsafe; requires `--unsafe-no-clip-screen`).
    pub unsafe_no_clip_screen: bool,
    /// Defer a record when minimap2's competing chain score `s2` is at least this
    /// fraction of the best chain score `s1` (MAPQ instability).
    pub max_s2_ratio: f64,
    /// Runtime state, set by `lift`/`translate` from the contract: true only when
    /// the input is a `cramhook lift` chain and the current lift's
    /// eligibility/safety options are all default ([`TranslateConfig::eligibility_is_default`]).
    /// A missing `s2` on an inherited `XL:A:T` record is then trusted (rule 3b).
    /// Not a CLI option.
    pub inherited_xl_ok: bool,
    /// Defer a record translated through a `-` block whose CIGAR contains an
    /// insertion or deletion of at least this many bp.
    pub max_minus_indel: u64,
    /// Defer a record on a `-` block whose old-alignment path Z-drop statistic
    /// exceeds this (minimap2's `-z`; `0` is a legitimate strict value).
    pub max_zdrop: i64,
    /// ksw2 scoring for the Z-drop statistic. Set from the read preset at
    /// runtime (not a CLI option); defaults to lr:hq.
    pub scoring: Scoring,
}

impl Default for TranslateConfig {
    fn default() -> Self {
        Self {
            strict: false,
            xl_tag: true,
            min_mapq: 60,
            max_clip: 200,
            margin: 500,
            max_de: None,
            t2_win: crate::diff::WIN,
            t2_step: crate::diff::STEP,
            t2_flank: crate::diff::FLANK,
            t2_junction_flank: crate::diff::JUNCTION_FLANK,
            clip_k: 17,
            clip_min_hits: 2,
            clip_min_len: 20,
            clip_flank: 0,
            clip_subtract_old: false,
            unsafe_no_clip_screen: false,
            max_s2_ratio: DEFAULT_MAX_S2_RATIO,
            max_minus_indel: 50,
            max_zdrop: PRESET_ZDROP,
            scoring: Scoring::default(),
            inherited_xl_ok: false,
        }
    }
}

impl TranslateConfig {
    /// The T2 windowing parameters this config requests.
    pub fn diff_params(&self) -> crate::diff::DiffParams {
        crate::diff::DiffParams {
            win: self.t2_win,
            step: self.t2_step,
            flank: self.t2_flank,
            junction_flank: self.t2_junction_flank,
        }
    }

    /// True when every option that affects *which* records are translated or
    /// *how* is at its default. A chained input (a CRAM produced by `cramhook
    /// lift`) is only accepted when this also holds for the current lift, so an
    /// inherited `XL:A:T` record's stripped `s2` and the clips screened by the
    /// earlier, also-all-default lifts can be trusted. `strict`, `xl_tag` and the
    /// runtime `scoring`/`inherited_xl_ok` do not affect translation and are
    /// ignored.
    pub fn eligibility_is_default(&self) -> bool {
        let d = Self::default();
        self.min_mapq == d.min_mapq
            && self.max_clip == d.max_clip
            && self.margin == d.margin
            && self.max_de == d.max_de
            && self.t2_win == d.t2_win
            && self.t2_step == d.t2_step
            && self.t2_flank == d.t2_flank
            && self.t2_junction_flank == d.t2_junction_flank
            && self.clip_k == d.clip_k
            && self.clip_min_hits == d.clip_min_hits
            && self.clip_min_len == d.clip_min_len
            && self.clip_flank == d.clip_flank
            && self.clip_subtract_old == d.clip_subtract_old
            && self.unsafe_no_clip_screen == d.unsafe_no_clip_screen
            && self.max_s2_ratio == d.max_s2_ratio
            && self.max_minus_indel == d.max_minus_indel
            && self.max_zdrop == d.max_zdrop
    }

    /// Validate the tuning knobs. Called before any work by `translate`/`lift`.
    pub fn validate(&self) -> Result<()> {
        if self.clip_k < 1 || self.clip_k > 32 {
            bail!("--clip-k must be in 1..=32 (got {})", self.clip_k);
        }
        if let Some(d) = self.max_de {
            if !d.is_finite() || d < 0.0 {
                bail!("--max-de must be finite and >= 0 (got {d})");
            }
        }
        if self.t2_win == 0 {
            bail!("--t2-win must be > 0");
        }
        if self.t2_step == 0 {
            bail!("--t2-step must be > 0");
        }
        if self.t2_step > self.t2_win {
            bail!("--t2-step ({}) must be <= --t2-win ({})", self.t2_step, self.t2_win);
        }
        if self.t2_win > crate::diff::DiffParams::MAX_T2 {
            bail!("--t2-win ({}) must be <= {}", self.t2_win, crate::diff::DiffParams::MAX_T2);
        }
        if self.t2_flank > crate::diff::DiffParams::MAX_T2 {
            bail!("--t2-flank ({}) must be <= {}", self.t2_flank, crate::diff::DiffParams::MAX_T2);
        }
        if self.t2_junction_flank > crate::diff::DiffParams::MAX_T2 {
            bail!(
                "--t2-junction-flank ({}) must be <= {}",
                self.t2_junction_flank,
                crate::diff::DiffParams::MAX_T2
            );
        }
        if !self.max_s2_ratio.is_finite() || self.max_s2_ratio < 0.0 {
            bail!("--max-s2-ratio must be finite and >= 0 (got {})", self.max_s2_ratio);
        }
        if self.max_zdrop < 0 {
            bail!("--max-zdrop must be >= 0 (got {})", self.max_zdrop);
        }
        Ok(())
    }
}

/// Effective tuning parameters, recorded in `stats.json` under `params`.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct Params {
    pub min_mapq: u8,
    /// Effective `xl_tag` (default true): translated records carry `XL:A:T`,
    /// remapped records `XL:A:R`.
    pub xl_tag: bool,
    pub max_clip: u64,
    pub margin: u64,
    pub max_de: f64,
    pub t2_win: u64,
    pub t2_step: u64,
    pub t2_flank: u64,
    pub t2_junction_flank: u64,
    pub clip_k: u32,
    pub clip_min_hits: u32,
    pub clip_min_len: u64,
    pub clip_flank: u64,
    pub clip_subtract_old: bool,
    pub unsafe_no_clip_screen: bool,
    pub max_s2_ratio: f64,
    pub max_minus_indel: u64,
    pub max_zdrop: i64,
    /// M9: effective tiny-input fallback threshold in MB (10^6 bytes);
    /// 0 disables the fallback. `--min-input-mb` or the read-preset default.
    pub min_input_mb: f64,
    /// M9: effective worker thread count.
    pub threads: usize,
}

impl Params {
    pub fn new(cfg: &TranslateConfig, max_de: f64) -> Params {
        Params {
            min_mapq: cfg.min_mapq,
            xl_tag: cfg.xl_tag,
            max_clip: cfg.max_clip,
            margin: cfg.margin,
            max_de,
            t2_win: cfg.t2_win,
            t2_step: cfg.t2_step,
            t2_flank: cfg.t2_flank,
            t2_junction_flank: cfg.t2_junction_flank,
            clip_k: cfg.clip_k,
            clip_min_hits: cfg.clip_min_hits,
            clip_min_len: cfg.clip_min_len,
            clip_flank: cfg.clip_flank,
            clip_subtract_old: cfg.clip_subtract_old,
            unsafe_no_clip_screen: cfg.unsafe_no_clip_screen,
            max_s2_ratio: cfg.max_s2_ratio,
            max_minus_indel: cfg.max_minus_indel,
            max_zdrop: cfg.max_zdrop,
            min_input_mb: 0.0,
            threads: 0,
        }
    }
}

/// Why a read was deferred. `name()` is stable and used in `stats.json`.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub enum Reason {
    Unmapped,
    Split,
    Mapq,
    Competitor,
    De,
    Clip,
    Tags,
    NoBlock,
    TaintT1,
    TaintT2,
    /// T2 junction taint: the read overlaps an old base tainted only by a new
    /// adjacency (deletion/join), not by a novel base.
    TaintJunction,
    Junction,
    ClipHit,
    MinusIndel,
    Zdrop,
    SplitNotTwoRecords,
    SplitDiffContigs,
    SplitNotOrigin,
    SplitSaInconsistent,
    SplitMapq,
    SplitCompetitor,
    SplitDe,
    SplitClip,
    SplitTags,
    SplitNoBlock,
    SplitTaintT1,
    SplitTaintT2,
    /// Split-unit junction taint (see [`Reason::TaintJunction`]).
    SplitTaintJunction,
    SplitJunction,
    SplitClipHit,
    SplitMinusIndel,
    SplitZdrop,
    Contract,
    /// M9: all reads remapped because the input CRAM is below the tiny-input
    /// threshold. Not a contract violation (`contract_ok` stays true).
    TinyInput,
}

/// Deferral class (Plan 4.2): read-level, diff-level, contract or fallback.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub enum Class {
    Read,
    Diff,
    Contract,
    /// M9: the tiny-input all-remap fallback.
    Fallback,
}

impl Class {
    pub fn name(self) -> &'static str {
        match self {
            Class::Read => "read",
            Class::Diff => "diff",
            Class::Contract => "contract",
            Class::Fallback => "fallback",
        }
    }
}

impl Reason {
    pub fn name(self) -> &'static str {
        match self {
            Reason::Unmapped => "unmapped",
            Reason::Split => "split",
            Reason::Mapq => "mapq",
            Reason::Competitor => "competitor",
            Reason::De => "de",
            Reason::Clip => "clip",
            Reason::Tags => "tags",
            Reason::NoBlock => "no_block",
            Reason::TaintT1 => "taint_t1",
            Reason::TaintT2 => "taint_t2",
            Reason::TaintJunction => "taint_junction",
            Reason::Junction => "junction",
            Reason::ClipHit => "clip_hit",
            Reason::MinusIndel => "minus_indel",
            Reason::Zdrop => "zdrop",
            Reason::SplitNotTwoRecords => "split_not_two_records",
            Reason::SplitDiffContigs => "split_diff_contigs",
            Reason::SplitNotOrigin => "split_not_origin",
            Reason::SplitSaInconsistent => "split_sa_inconsistent",
            Reason::SplitMapq => "split_mapq",
            Reason::SplitCompetitor => "split_competitor",
            Reason::SplitDe => "split_de",
            Reason::SplitClip => "split_clip",
            Reason::SplitTags => "split_tags",
            Reason::SplitNoBlock => "split_no_block",
            Reason::SplitTaintT1 => "split_taint_t1",
            Reason::SplitTaintT2 => "split_taint_t2",
            Reason::SplitTaintJunction => "split_taint_junction",
            Reason::SplitJunction => "split_junction",
            Reason::SplitClipHit => "split_clip_hit",
            Reason::SplitMinusIndel => "split_minus_indel",
            Reason::SplitZdrop => "split_zdrop",
            Reason::Contract => "contract",
            Reason::TinyInput => "tiny_input",
        }
    }

    pub fn class(self) -> Class {
        match self {
            Reason::Unmapped
            | Reason::Split
            | Reason::Mapq
            | Reason::Competitor
            | Reason::De
            | Reason::Clip
            | Reason::Tags
            | Reason::SplitNotTwoRecords
            | Reason::SplitDiffContigs
            | Reason::SplitNotOrigin
            | Reason::SplitSaInconsistent
            | Reason::SplitMapq
            | Reason::SplitCompetitor
            | Reason::SplitDe
            | Reason::SplitClip
            | Reason::SplitTags => Class::Read,
            Reason::NoBlock
            | Reason::TaintT1
            | Reason::TaintT2
            | Reason::TaintJunction
            | Reason::Junction
            | Reason::ClipHit
            | Reason::MinusIndel
            | Reason::Zdrop
            | Reason::SplitNoBlock
            | Reason::SplitTaintT1
            | Reason::SplitTaintT2
            | Reason::SplitTaintJunction
            | Reason::SplitJunction
            | Reason::SplitClipHit
            | Reason::SplitMinusIndel
            | Reason::SplitZdrop => Class::Diff,
            Reason::Contract => Class::Contract,
            Reason::TinyInput => Class::Fallback,
        }
    }
}

// ---------------------------------------------------------------------------
// pure helpers (unit-tested)

#[inline]
fn consumes_query(c: Cigar) -> bool {
    matches!(c, Cigar::Match(_) | Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::Equal(_) | Cigar::Diff(_))
}

#[inline]
fn consumes_ref(c: Cigar) -> bool {
    matches!(c, Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::Equal(_) | Cigar::Diff(_))
}

/// Bases consumed on the query by aligned (non-clip) CIGAR ops.
pub fn query_span(cigar: &[Cigar]) -> u64 {
    cigar.iter().filter(|c| consumes_query(**c) && !matches!(c, Cigar::SoftClip(_))).map(|c| c.len() as u64).sum()
}

/// Bases consumed on the reference by aligned CIGAR ops.
pub fn ref_span(cigar: &[Cigar]) -> u64 {
    cigar.iter().filter(|c| consumes_ref(**c)).map(|c| c.len() as u64).sum()
}

/// Length of the longest insertion or deletion in the CIGAR (0 if none).
pub fn max_indel(cigar: &[Cigar]) -> u64 {
    cigar
        .iter()
        .filter_map(|c| match c {
            Cigar::Ins(l) | Cigar::Del(l) => Some(*l as u64),
            _ => None,
        })
        .max()
        .unwrap_or(0)
}

/// Number of leading soft-clip bases.
pub fn leading_soft(cigar: &[Cigar]) -> u64 {
    let mut n = 0;
    for c in cigar {
        match c {
            Cigar::SoftClip(l) => n += *l as u64,
            _ => break,
        }
    }
    n
}

/// Number of trailing soft-clip bases.
pub fn trailing_soft(cigar: &[Cigar]) -> u64 {
    let mut n = 0;
    for c in cigar.iter().rev() {
        match c {
            Cigar::SoftClip(l) => n += *l as u64,
            _ => break,
        }
    }
    n
}

/// Exclusive reference end of an alignment.
pub fn ref_end(pos: i64, cigar: &[Cigar]) -> i64 {
    pos + ref_span(cigar) as i64
}

/// True if the CIGAR contains any hard clip.
pub fn has_hard_clip(cigar: &[Cigar]) -> bool {
    cigar.iter().any(|c| matches!(c, Cigar::HardClip(_)))
}

/// Reverse the operation order of a CIGAR (used for `-` blocks).
pub fn reverse_cigar(cigar: &[Cigar]) -> Vec<Cigar> {
    cigar.iter().rev().copied().collect()
}

#[inline]
fn is_match_op(c: Cigar) -> bool {
    matches!(c, Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_))
}

#[inline]
fn op_query_len(c: Cigar) -> u32 {
    match c {
        Cigar::Match(l) | Cigar::Ins(l) | Cigar::SoftClip(l) | Cigar::Equal(l) | Cigar::Diff(l) => l,
        _ => 0,
    }
}

#[inline]
fn op_ref_len(c: Cigar) -> u32 {
    match c {
        Cigar::Match(l) | Cigar::Del(l) | Cigar::RefSkip(l) | Cigar::Equal(l) | Cigar::Diff(l) => l,
        _ => 0,
    }
}

fn with_len(c: Cigar, l: u32) -> Cigar {
    match c {
        Cigar::Match(_) => Cigar::Match(l),
        Cigar::Ins(_) => Cigar::Ins(l),
        Cigar::Del(_) => Cigar::Del(l),
        Cigar::RefSkip(_) => Cigar::RefSkip(l),
        Cigar::SoftClip(_) => Cigar::SoftClip(l),
        Cigar::HardClip(_) => Cigar::HardClip(l),
        Cigar::Pad(_) => Cigar::Pad(l),
        Cigar::Equal(_) => Cigar::Equal(l),
        Cigar::Diff(_) => Cigar::Diff(l),
    }
}

fn grow(c: Cigar, n: u32) -> Cigar {
    with_len(c, c.len() + n)
}

/// Remove zero-length ops and merge adjacent ops of the same type.
fn tidy_cigar(ops: Vec<Cigar>) -> Vec<Cigar> {
    let mut out: Vec<Cigar> = Vec::with_capacity(ops.len());
    for c in ops {
        if c.len() == 0 {
            continue;
        }
        if let Some(last) = out.last_mut() {
            if last.char() == c.char() {
                *last = grow(*last, c.len());
                continue;
            }
        }
        out.push(c);
    }
    out
}

/// Normalise every indel to its leftmost equivalent position on the new
/// reference (Plan 4.3). minimap2 left-aligns indels in reference-forward
/// orientation; a CIGAR produced by reversing a `-` alignment has each indel at
/// the right end of its equivalent range, so it must be shifted back.
///
/// A deletion shifts left while `ref[base before the gap] == ref[last base of
/// the gap]`; an insertion shifts left while `read[base before] == read[last
/// inserted base]`. Shifting consumes one base of the immediately preceding
/// match run and grows the following match run, so an indel never crosses
/// another indel or a clip.
pub fn left_align_cigar(ops: &[Cigar], read: &[u8], ref_slice: &[u8]) -> Vec<Cigar> {
    let mut ops: Vec<Cigar> = ops.to_vec();
    let mut q = 0usize;
    let mut r = 0usize;
    let mut i = 0usize;
    while i < ops.len() {
        match ops[i] {
            Cigar::Del(d) => {
                let d = d as usize;
                if i > 0 && i + 1 < ops.len() && is_match_op(ops[i - 1]) && is_match_op(ops[i + 1]) {
                    loop {
                        let prev = ops[i - 1].len() as usize;
                        if prev == 0 || r == 0 || r + d > ref_slice.len() {
                            break;
                        }
                        if ref_slice[r - 1] != ref_slice[r + d - 1] {
                            break;
                        }
                        ops[i - 1] = with_len(ops[i - 1], (prev - 1) as u32);
                        ops[i + 1] = grow(ops[i + 1], 1);
                        r -= 1;
                        q -= 1;
                    }
                }
                r += d;
            }
            Cigar::Ins(ins) => {
                let ins = ins as usize;
                if i > 0 && i + 1 < ops.len() && is_match_op(ops[i - 1]) && is_match_op(ops[i + 1]) {
                    loop {
                        let prev = ops[i - 1].len() as usize;
                        if prev == 0 || q == 0 || q + ins > read.len() {
                            break;
                        }
                        if read[q - 1] != read[q + ins - 1] {
                            break;
                        }
                        ops[i - 1] = with_len(ops[i - 1], (prev - 1) as u32);
                        ops[i + 1] = grow(ops[i + 1], 1);
                        q -= 1;
                        r -= 1;
                    }
                }
                q += ins;
            }
            other => {
                q += op_query_len(other) as usize;
                r += op_ref_len(other) as usize;
            }
        }
        i += 1;
    }
    tidy_cigar(ops)
}

/// Reverse-complement a read, complementing IUPAC codes; unknown bases stay.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| COMP[b as usize]).collect()
}

/// Complement lookup table (branch-free), matching [`revcomp`]'s IUPAC map.
const fn build_comp() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        let b = i as u8;
        t[i] = match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'U' => b'A',
            b'R' => b'Y', // A/G
            b'Y' => b'R', // C/T
            b'S' => b'S', // C/G
            b'W' => b'W', // A/T
            b'K' => b'M', // G/T
            b'M' => b'K', // A/C
            b'B' => b'V', // C/G/T
            b'V' => b'B', // A/C/G
            b'D' => b'H', // A/G/T
            b'H' => b'D', // A/C/T
            b'N' => b'N',
            other => other,
        };
        i += 1;
    }
    t
}

static COMP: [u8; 256] = build_comp();

/// Complement of a BAM 4-bit base, indexed by the nibble. `=ACMGRSVTWYHKDBN`
/// is rust-htslib's `DECODE_BASE`; composing it with `COMP` lets the `-`
/// revcomp decode and complement in a single table lookup.
const fn build_decode_comp() -> [u8; 16] {
    let decode = b"=ACMGRSVTWYHKDBN";
    let mut t = [0u8; 16];
    let mut i = 0;
    while i < 16 {
        t[i] = COMP[decode[i] as usize];
        i += 1;
    }
    t
}
static DECODE_COMP: [u8; 16] = build_decode_comp();

/// For a byte holding two BAM-encoded bases (high nibble = earlier base, low
/// nibble = later base), the two complemented bases in reverse order, packed so
/// the low byte is the complement of the low nibble and the high byte the
/// complement of the high nibble. Reversing a read's encoded bytes and emitting
/// this pair per byte performs the whole revcomp, halving the table lookups and
/// avoiding per-base decoding.
const fn build_pair_comp() -> [u16; 256] {
    let mut t = [0u16; 256];
    let mut b = 0usize;
    while b < 256 {
        let lo = DECODE_COMP[b & 0xf] as u16;
        let hi = DECODE_COMP[(b >> 4) & 0xf] as u16;
        t[b] = lo | (hi << 8);
        b += 1;
    }
    t
}
static PAIR_COMP: [u16; 256] = build_pair_comp();

/// Reverse-complement a record's sequence straight from its 4-bit encoding.
/// Reproduces `revcomp(rec.seq().as_bytes())` exactly, without materialising the
/// forward ASCII sequence or re-encoding it.
fn revcomp_seq(rec: &bam::Record) -> Vec<u8> {
    let sq = rec.seq();
    let encoded = sq.encoded;
    let n = sq.len();
    // Guard against an unexpected packing; the canonical BAM layout is
    // `ceil(n / 2)` bytes, which is what the pair table assumes.
    if encoded.len() != n.div_ceil(2) {
        return revcomp(&sq.as_bytes());
    }
    let mut seq: Vec<u8> = vec![0u8; n];
    if n % 2 == 1 {
        // Odd length: the last byte's low nibble is padding, so its high nibble
        // is the final read base and becomes the first output base.
        seq[0] = DECODE_COMP[(encoded[encoded.len() - 1] >> 4) as usize];
        for (&raw, out) in
            encoded[..encoded.len() - 1].iter().rev().zip(seq[1..].chunks_exact_mut(2))
        {
            let p = PAIR_COMP[raw as usize];
            out[0] = p as u8;
            out[1] = (p >> 8) as u8;
        }
    } else {
        for (&raw, out) in encoded.iter().rev().zip(seq.chunks_exact_mut(2)) {
            let p = PAIR_COMP[raw as usize];
            out[0] = p as u8;
            out[1] = (p >> 8) as u8;
        }
    }
    seq
}

/// Eligibility window in old coordinates (Plan 4.2 rule 7).
pub fn window(pos: i64, ref_end: i64, left_clip: u64, right_clip: u64, margin: u64, contig_len: u64) -> OldIv {
    let s = (pos - left_clip as i64 - margin as i64).max(0) as u64;
    let e = (ref_end + right_clip as i64 + margin as i64).min(contig_len as i64).max(s as i64) as u64;
    OldIv::new(s, e)
}

/// One parsed `SA` entry: `rname,pos,strand,CIGAR,mapq,NM`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SaEntry {
    pub rname: String,
    pub pos: u64,
    pub strand: char,
    pub cigar: String,
    pub mapq: u32,
    pub nm: i64,
}

/// Parse an `SA:Z:` payload (one or more `;`-terminated entries). Every field
/// is validated strictly; a malformed payload is an error (caller defers with
/// reason `split`).
pub fn parse_sa(s: &str) -> Result<Vec<SaEntry>> {
    if s.is_empty() {
        return Ok(Vec::new());
    }
    if !s.ends_with(';') {
        bail!("SA payload is not terminated by ';': {s:?}");
    }
    let parts: Vec<&str> = s.split(';').collect();
    let mut out = Vec::new();
    for part in &parts[..parts.len() - 1] {
        if part.is_empty() {
            bail!("SA payload contains an empty entry");
        }
        let f: Vec<&str> = part.split(',').collect();
        if f.len() != 6 {
            bail!("SA entry {part:?} does not have 6 fields");
        }
        if f[0].is_empty() {
            bail!("SA entry {part:?} has an empty rname");
        }
        let pos: u64 = f[1].parse().with_context(|| format!("SA pos {:?}", f[1]))?;
        if pos == 0 {
            bail!("SA entry {part:?} has a non-positive pos");
        }
        let strand = match f[2] {
            "+" => '+',
            "-" => '-',
            other => bail!("SA entry {part:?} has an invalid strand {other:?}"),
        };
        let cigar_ok = bam::record::CigarString::try_from(f[3])
            .map(|c| !c.0.is_empty())
            .unwrap_or(false);
        if !cigar_ok {
            bail!("SA entry {part:?} has an invalid or empty CIGAR {:?}", f[3]);
        }
        let mapq: u32 = f[4].parse().with_context(|| format!("SA mapq {:?}", f[4]))?;
        if mapq > 255 {
            bail!("SA entry {part:?} has mapq {mapq} > 255");
        }
        let nm: i64 = f[5].parse().with_context(|| format!("SA NM {:?}", f[5]))?;
        if nm < 0 {
            bail!("SA entry {part:?} has a negative NM {nm}");
        }
        out.push(SaEntry {
            rname: f[0].to_string(),
            pos,
            strand,
            cigar: f[3].to_string(),
            mapq,
            nm,
        });
    }
    Ok(out)
}

/// Everything `format_sa` needs about the *other* record, in its own
/// orientation, to reproduce minimap2's `SA` CIGAR simplification
/// (`ref/minimap2/format.c`, the `SA:Z:` writing code): terminal soft clips and
/// a single `M` plus at most one `I`/`D` derived from the aligned spans.
#[derive(Clone, Debug)]
pub struct SaSource {
    pub contig: String,
    pub pos1: u64,
    pub reverse: bool,
    pub cigar: Vec<Cigar>,
    pub l_seq: u64,
    pub mapq: u8,
    pub nm: i64,
}

/// The CIGAR field of one `SA` entry, exactly as minimap2 writes it
/// (`ref/minimap2/format.c`: leading soft clip, a single `M`, at most one
/// `I`/`D` from the aligned spans, trailing soft clip).
pub fn sa_cigar(cigar: &[Cigar]) -> String {
    let clip5 = leading_soft(cigar);
    let clip3 = trailing_soft(cigar);
    let q = query_span(cigar);
    let r = ref_span(cigar);
    let (l_m, l_i, l_d) = if q >= r { (r, q - r, 0) } else { (q, 0, r - q) };
    let mut s = String::new();
    use std::fmt::Write as _;
    if clip5 > 0 {
        let _ = write!(s, "{clip5}S");
    }
    if l_m > 0 {
        let _ = write!(s, "{l_m}M");
    }
    if l_i > 0 {
        let _ = write!(s, "{l_i}I");
    }
    if l_d > 0 {
        let _ = write!(s, "{l_d}D");
    }
    if clip3 > 0 {
        let _ = write!(s, "{clip3}S");
    }
    s
}

/// Format one `SA` entry exactly as minimap2 does, including the trailing `;`.
pub fn format_sa(src: &SaSource) -> String {
    format!(
        "{},{},{},{},{},{};",
        src.contig,
        src.pos1,
        if src.reverse { '-' } else { '+' },
        sa_cigar(&src.cigar),
        src.mapq,
        src.nm
    )
}

// ---------------------------------------------------------------------------
// record view used by the classifier (also constructible in unit tests)

/// A read record distilled from the BAM: only what eligibility needs.
#[derive(Clone, Debug)]
pub struct RecInfo {
    pub qname: Vec<u8>,
    pub flag: u16,
    pub contig: String,
    pub pos: i64,
    pub mapq: u8,
    pub cigar: Vec<Cigar>,
    pub seq_len: usize,
    pub de: Option<f64>,
    pub nm: Option<i64>,
    /// minimap2 best chain score (`s1`), used for the MAPQ-instability rule.
    pub s1: Option<i64>,
    /// minimap2 best competing chain score (`s2`).
    pub s2: Option<i64>,
    /// `SA` tag is present (regardless of whether it parsed).
    pub sa_present: bool,
    /// `SA` is present, non-empty, a string and parses cleanly.
    pub sa_ok: bool,
    pub sa: Vec<SaEntry>,
    pub unexpected_tag: bool,
    /// The record carries `XL:A:T`: it is the output of a previous `cramhook
    /// lift` that translated it. Its `s2` was stripped by that lift, so rule 3b
    /// treats the missing competitor as already vetted instead of deferring.
    pub xl_translated: bool,
    /// The rule-7d path Z-drop statistic, computed once from the decode-generated
    /// `MD` tag and the preset scoring (`None` if `MD` is missing/malformed or
    /// inconsistent with the CIGAR). Only consulted on `-` blocks.
    pub zdrop: Option<i64>,
}

impl RecInfo {
    pub fn is_supplementary(&self) -> bool {
        self.flag & 0x800 != 0
    }
    pub fn is_primary(&self) -> bool {
        self.flag & (0x100 | 0x800) == 0
    }
    pub fn is_unmapped(&self) -> bool {
        self.flag & 0x4 != 0
    }
    pub fn is_reverse(&self) -> bool {
        self.flag & 0x10 != 0
    }
    pub fn ref_end(&self) -> i64 {
        ref_end(self.pos, &self.cigar)
    }
    pub fn outer_clip_reaches_end(&self, contig_len: u64) -> bool {
        self.ref_end() + trailing_soft(&self.cigar) as i64 >= contig_len as i64
    }
    pub fn outer_clip_reaches_start(&self) -> bool {
        self.pos - leading_soft(&self.cigar) as i64 <= 0
    }
}

/// A record that passed rules 3-7, with the block that will translate it.
pub struct Eligible<'a> {
    pub block: &'a Block,
}

/// Map a record read-level/diff-level failure to the matching split reason.
fn split_reason(r: Reason) -> Reason {
    match r {
        Reason::Mapq => Reason::SplitMapq,
        Reason::Competitor => Reason::SplitCompetitor,
        Reason::De => Reason::SplitDe,
        Reason::Clip => Reason::SplitClip,
        Reason::Tags => Reason::SplitTags,
        Reason::NoBlock => Reason::SplitNoBlock,
        Reason::TaintT1 => Reason::SplitTaintT1,
        Reason::TaintT2 => Reason::SplitTaintT2,
        Reason::TaintJunction => Reason::SplitTaintJunction,
        Reason::Junction => Reason::SplitJunction,
        Reason::ClipHit => Reason::SplitClipHit,
        Reason::MinusIndel => Reason::SplitMinusIndel,
        Reason::Zdrop => Reason::SplitZdrop,
        other => other,
    }
}

/// Which contig boundary an unclipped window crosses.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Boundary {
    Start,
    End,
}

/// The origin-crossing junction test, factored out: on `boundary`, the block
/// must reach the old contig boundary and map it to the corresponding new
/// boundary (mirrored for `-`).
pub fn junction_preserved(block: &Block, contig_len: u64, boundary: Boundary, new_len: u64) -> bool {
    if new_len == 0 {
        return false;
    }
    match boundary {
        Boundary::End => {
            if block.old.end.0 != contig_len {
                return false;
            }
            let want = match block.strand {
                Strand::Plus => new_len - 1,
                Strand::Minus => 0,
            };
            block.map_pos(OldPos(contig_len - 1)).map(|p| p.get()) == Some(want)
        }
        Boundary::Start => {
            if block.old.start.0 != 0 {
                return false;
            }
            let want = match block.strand {
                Strand::Plus => 0,
                Strand::Minus => new_len - 1,
            };
            block.map_pos(OldPos(0)).map(|p| p.get()) == Some(want)
        }
    }
}

/// Check rules 3-7 for one record. `clip_to_check` is the total clip for a
/// single record, or the outer (non-junction) clip for an origin-crossing one.
///
/// If the *unclipped* window (`[ref_start − left_clip − margin, ref_end +
/// right_clip + margin)`) crosses a contig boundary, translation is only
/// allowed when the junction is preserved at that boundary; otherwise the read
/// is deferred with reason [`Reason::Junction`] (a read whose clipped bases
/// continue across the contig end must not be lifted onto a rotated/edited
/// junction). A boundary-crossing window that no block covers is
/// [`Reason::NoBlock`], exactly like a non-crossing uncovered window.
pub fn check_record<'a>(
    rec: &RecInfo,
    contig_len: u64,
    diff: &'a Diff,
    thr: f64,
    cfg: &TranslateConfig,
    clip_to_check: u64,
    new_len_of: &dyn Fn(&str) -> Option<u64>,
) -> Result<&'a Block, Reason> {
    if (rec.mapq as u32) < cfg.min_mapq as u32 {
        return Err(Reason::Mapq);
    }
    // Rule 3b: MAPQ instability from near-equal chains. A close competitor
    // (s2/s1 >= max_s2_ratio, or missing/invalid s1/s2) can flip the MAPQ when
    // the competitor's score shifts in the new assembly, so the read is deferred.
    // s1 must be positive and s2 non-negative; a float/string-typed or missing
    // tag was already rejected by `aux_int_strict` (-> `None`).
    //
    // A record carrying `XL:A:T` is cramhook's own translated output: the
    // previous lift stripped its (stale, un-recomputable) `s2`, and the read
    // already passed this rule when that lift translated it. That is only sound
    // for an all-default chain (every earlier lift and the current lift use the
    // default `--max-s2-ratio`), which the contract enforces before
    // `cfg.inherited_xl_ok` is set. `XL:A:R` records keep minimap2's fresh `s2`
    // and are checked normally.
    match (rec.s1, rec.s2) {
        (Some(s1), Some(s2))
            if s1 > 0 && s2 >= 0 && (s2 as f64) < cfg.max_s2_ratio * (s1 as f64) => {}
        (Some(s1), None) if rec.xl_translated && cfg.inherited_xl_ok && s1 > 0 => {}
        _ => return Err(Reason::Competitor),
    }
    let Some(de) = rec.de else {
        return Err(Reason::Tags);
    };
    if !de.is_finite() || de < 0.0 || de > thr {
        return Err(Reason::De);
    }
    if cfg.max_clip != 0 && clip_to_check > cfg.max_clip {
        return Err(Reason::Clip);
    }
    if rec.unexpected_tag {
        return Err(Reason::Tags);
    }
    let left = leading_soft(&rec.cigar);
    let right = trailing_soft(&rec.cigar);
    let raw_start = rec.pos - left as i64 - cfg.margin as i64;
    let raw_end = rec.ref_end() + right as i64 + cfg.margin as i64;
    let crosses_start = raw_start < 0;
    let crosses_end = raw_end > contig_len as i64;

    let win = window(rec.pos, rec.ref_end(), left, right, cfg.margin, contig_len);
    let Some(block) = diff.block_containing(&rec.contig, win) else {
        // No block covers the (clipped) window: this is `no_block`, boundary or
        // not. `Junction` is reserved for a block that covers the window but
        // whose mapping across the contig boundary is not preserved.
        return Err(Reason::NoBlock);
    };
    if crosses_start {
        let new_len = new_len_of(&block.new_contig).ok_or(Reason::Junction)?;
        if !junction_preserved(block, contig_len, Boundary::Start, new_len) {
            return Err(Reason::Junction);
        }
    }
    if crosses_end {
        let new_len = new_len_of(&block.new_contig).ok_or(Reason::Junction)?;
        if !junction_preserved(block, contig_len, Boundary::End, new_len) {
            return Err(Reason::Junction);
        }
    }
    if diff.overlaps_t1(&rec.contig, win) {
        return Err(Reason::TaintT1);
    }
    if diff.overlaps_t2(&rec.contig, win) {
        return Err(Reason::TaintT2);
    }
    // Rule 7 (junction): a new adjacency (deletion/join) can make an old region
    // a better competitor without adding a novel base, so T2's `novel` set is
    // empty there. The junction taint covers it.
    if diff.overlaps_junction(&rec.contig, win) {
        return Err(Reason::TaintJunction);
    }
    // Rule 7b: a large indel on a '-' block is minimap2's own non-mirror
    // direction effect at a structural change; translating it can move the
    // position or drop a supplementary. '+' blocks are unaffected.
    if block.strand == Strand::Minus && max_indel(&rec.cigar) >= cfg.max_minus_indel {
        return Err(Reason::MinusIndel);
    }
    // Rule 7d: on a '-' block, defer when the old alignment's path Z-drop
    // statistic exceeds minimap2's `-z`. ksw2's z-drop can split an alignment in
    // one orientation but not the other, so such a read would stay unsplit after
    // translation while a fresh remap splits it. The statistic was computed once
    // in `to_info`; `None` (missing/malformed/inconsistent `MD`) defers fail-safe.
    if block.strand == Strand::Minus && !matches!(rec.zdrop, Some(z) if z <= cfg.max_zdrop) {
        return Err(Reason::Zdrop);
    }
    Ok(block)
}

/// The two translated records of an origin-crossing split read.
pub struct SplitPlan<'a> {
    /// Index (into the input slice) of the record reaching the contig end.
    pub a: usize,
    /// Index of the record reaching the contig start.
    pub b: usize,
    pub a_block: &'a Block,
    pub b_block: &'a Block,
}

/// Decide whether a buffered group of records is a translatable
/// origin-crossing split read (Plan section 4.2). Returns the two blocks.
pub fn classify_origin_crossing<'a>(
    recs: &[RecInfo],
    contig_len: u64,
    diff: &'a Diff,
    thr: f64,
    cfg: &TranslateConfig,
    new_len_of: &dyn Fn(&str) -> Option<u64>,
) -> Result<SplitPlan<'a>, Reason> {
    if recs.len() != 2 {
        return Err(Reason::SplitNotTwoRecords);
    }
    let n_primary = recs.iter().filter(|r| r.is_primary()).count();
    let n_supp = recs.iter().filter(|r| r.is_supplementary()).count();
    if n_primary != 1 || n_supp != 1 {
        return Err(Reason::SplitNotTwoRecords);
    }
    if recs[0].contig != recs[1].contig {
        return Err(Reason::SplitDiffContigs);
    }

    // SA tags must be consistent with the pair.
    for (i, r) in recs.iter().enumerate() {
        let other = &recs[1 - i];
        if r.sa.len() != 1 {
            return Err(Reason::SplitNotOrigin);
        }
        let e = &r.sa[0];
        let expected_strand = if other.is_reverse() { '-' } else { '+' };
        let other_nm = other.nm.ok_or(Reason::SplitTags)?;
        if e.rname != other.contig
            || e.pos != (other.pos + 1) as u64
            || e.strand != expected_strand
            || e.mapq != other.mapq as u32
            || e.nm != other_nm
        {
            return Err(Reason::SplitNotOrigin);
        }
        // The SA CIGAR must be exactly what minimap2 would write for the other
        // buffered record (same simplification as `sa_cigar`).
        if e.cigar != sa_cigar(&other.cigar) {
            return Err(Reason::SplitSaInconsistent);
        }
    }

    // Decide which record reaches the end (A) and which the start (B).
    let mut assignment = None;
    for a in 0..2usize {
        let b = 1 - a;
        if recs[a].outer_clip_reaches_end(contig_len) && recs[b].outer_clip_reaches_start() {
            assignment = Some((a, b));
            break;
        }
    }
    let Some((a, b)) = assignment else {
        return Err(Reason::SplitNotOrigin);
    };

    // Rules 3-7 individually; rule 5 only on the outer clip.
    let a_outer = leading_soft(&recs[a].cigar);
    let b_outer = trailing_soft(&recs[b].cigar);
    let a_block = check_record(&recs[a], contig_len, diff, thr, cfg, a_outer, new_len_of)
        .map_err(split_reason)?;
    let b_block = check_record(&recs[b], contig_len, diff, thr, cfg, b_outer, new_len_of)
        .map_err(split_reason)?;

    // Junction preserved (old boundaries and new mapping). `check_record`
    // already enforces the boundary each record crosses; this is the shared
    // origin-crossing test, applied to both ends.
    if a_block.new_contig != b_block.new_contig || a_block.strand != b_block.strand {
        return Err(Reason::SplitJunction);
    }
    let new_len = new_len_of(&a_block.new_contig).ok_or(Reason::SplitJunction)?;
    if !junction_preserved(a_block, contig_len, Boundary::End, new_len)
        || !junction_preserved(b_block, contig_len, Boundary::Start, new_len)
    {
        return Err(Reason::SplitJunction);
    }

    Ok(SplitPlan { a, b, a_block, b_block })
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    fn m(n: u32) -> Cigar {
        Cigar::Match(n)
    }
    fn s(n: u32) -> Cigar {
        Cigar::SoftClip(n)
    }
    fn i(n: u32) -> Cigar {
        Cigar::Ins(n)
    }
    fn d(n: u32) -> Cigar {
        Cigar::Del(n)
    }

    #[test]
    fn cigar_reversal() {
        let c = vec![s(10), m(20), i(3), d(5), m(15), s(7)];
        assert_eq!(reverse_cigar(&c), vec![s(7), m(15), d(5), i(3), m(20), s(10)]);
        assert_eq!(reverse_cigar(&[]), Vec::<Cigar>::new());
    }

    #[test]
    fn revcomp_iupac() {
        assert_eq!(revcomp(b"ACGTN"), b"NACGT");
        assert_eq!(revcomp(b"acgt"), b"ACGT");
        assert_eq!(revcomp(b"RYSWKM"), b"KMWSRY");
        assert_eq!(revcomp(b"BDHVN"), b"NBDHV");
        assert_eq!(revcomp(b""), b"");
    }

    #[test]
    fn revcomp_seq_matches_revcomp_of_decoded_bytes() {
        // Decode sequence and reverse QUAL for one record, checking both against
        // the reference implementation: `revcomp(seq().as_bytes())` and the
        // reversed qual.
        fn check(seq: &[u8], qual: &[u8]) {
            assert_eq!(seq.len(), qual.len());
            let mut rec = bam::Record::new();
            let cigar = if seq.is_empty() {
                None
            } else {
                Some(bam::record::CigarString(vec![Cigar::Match(seq.len() as u32)]))
            };
            rec.set(b"r", cigar.as_ref(), seq, qual);
            assert_eq!(
                revcomp_seq(&rec),
                revcomp(&rec.seq().as_bytes()),
                "seq {:?}",
                String::from_utf8_lossy(seq)
            );
            // The `-` branch reverses QUAL; exercise it independently.
            let mut q = rec.qual().to_vec();
            q.reverse();
            let expected: Vec<u8> = qual.iter().rev().copied().collect();
            assert_eq!(q, expected, "qual for seq {:?}", String::from_utf8_lossy(seq));
        }

        // Even, odd, IUPAC, '=' (nibble 0), single-base and empty sequences.
        check(b"ACGTACGTAC", &[b'-'; 10]);
        check(b"ACGTACGTA", &[b'.'; 9]);
        check(b"RYSWKMNBDHV", &[b'|'; 11]);
        check(b"A", b"!");
        check(b"=ACGT=N", &[b'.'; 7]);
        check(b"acgtACGTn", &[b'/'; 9]);
        check(b"", b"");

        // Distinct (non-constant) qualities and the all-missing 0xff qual.
        let seq40: &[u8] = b"ACGTN=ACGTRYSWKMBDHVACGTN=ACGTRYSWKMBDHV";
        let distinct: Vec<u8> = (0..40u16).map(|i| (i * 3 % 94) as u8).collect();
        check(seq40, &distinct);
        check(seq40, &[0xffu8; 40]);

        // Nonzero padding nibble: build an even-length record, then drop one
        // base so l_qseq is odd while the last byte still holds the removed
        // base's (nonzero) nibble. It must be ignored.
        let mut rec = bam::Record::new();
        let cigar = bam::record::CigarString(vec![Cigar::Match(4)]);
        rec.set(b"r", Some(&cigar), b"ACGT", b"\x01\x02\x03\x04");
        rec.inner_mut().core.l_qseq = 3;
        assert_ne!(rec.seq().encoded[1] & 0xf, 0, "padding nibble should be nonzero");
        assert_eq!(revcomp_seq(&rec), revcomp(&rec.seq().as_bytes()));
        let mut q = rec.qual().to_vec();
        q.reverse();
        assert_eq!(q, vec![0x03u8, 0x02, 0x01]);
    }

    #[test]
    fn spans_and_clips() {
        let c = vec![s(100), m(50), i(2), d(4), m(60), s(30)];
        assert_eq!(leading_soft(&c), 100);
        assert_eq!(trailing_soft(&c), 30);
        assert_eq!(query_span(&c), 112);
        assert_eq!(ref_span(&c), 114);
        assert_eq!(ref_end(1000, &c), 1114);
    }

    #[test]
    fn window_clipping() {
        // near a contig start, left clip pushed below 0 -> clipped
        let w = window(10, 200, 100, 50, 500, 1000);
        assert_eq!((w.start.0, w.end.0), (0, 750));
        // near the end, right side clipped
        let w = window(900, 980, 10, 100, 500, 1000);
        assert_eq!((w.start.0, w.end.0), (390, 1000));
        // normal
        let w = window(500, 600, 0, 0, 500, 1000);
        assert_eq!((w.start.0, w.end.0), (0, 1100.min(1000)));
    }

    #[test]
    fn sa_format_matches_minimap2() {
        // r000199 primary: other (supplementary) has 0 leading clip, 13066
        // trailing clip, q span 5444, ref span 5441 -> 5441M3I13066S.
        let src = SaSource {
            contig: "chr".into(),
            pos1: 294560,
            reverse: false,
            cigar: vec![m(349), i(1), m(668), m(505), m(321), m(219), m(308), m(431), m(846), m(928), m(203), m(251), s(13066)],
            l_seq: 18510,
            mapq: 60,
            nm: 16,
        };
        // q span = 349+1+668+505+321+219+308+431+846+928+203+251 = 5030? no, craft a simpler one below.
        // Use a synthetic one with known spans instead for exactness.
        let src2 = SaSource {
            contig: "chr".into(),
            pos1: 294560,
            reverse: false,
            cigar: vec![m(5441), i(3), s(13066)],
            l_seq: 18510,
            mapq: 60,
            nm: 16,
        };
        assert_eq!(format_sa(&src2), "chr,294560,+,5441M3I13066S,60,16;");
        // primary described from the supplementary: 5444S13064M2I
        let src3 = SaSource {
            contig: "chr".into(),
            pos1: 1,
            reverse: false,
            cigar: vec![s(5444), m(13064), i(2)],
            l_seq: 18510,
            mapq: 60,
            nm: 41,
        };
        assert_eq!(format_sa(&src3), "chr,1,+,5444S13064M2I,60,41;");
        // reverse record: clips are leading/trailing of the SAM CIGAR
        let src4 = SaSource {
            contig: "chr".into(),
            pos1: 292176,
            reverse: true,
            cigar: vec![m(7825), i(1), s(9509)],
            l_seq: 17336,
            mapq: 60,
            nm: 19,
        };
        assert_eq!(format_sa(&src4), "chr,292176,-,7825M1I9509S,60,19;");
        let _ = src;
    }

    #[test]
    fn parse_sa_roundtrip() {
        let v = parse_sa("chr,294560,+,5441M3I13066S,60,16;").unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].rname, "chr");
        assert_eq!(v[0].pos, 294560);
        assert_eq!(v[0].strand, '+');
        assert_eq!(v[0].cigar, "5441M3I13066S");
        assert_eq!(v[0].mapq, 60);
        assert_eq!(v[0].nm, 16);
        assert!(parse_sa("bad,entry").is_err());
    }

    #[test]
    fn parse_sa_is_strict() {
        // valid
        assert!(parse_sa("chr,1,+,100M,60,0;").is_ok());
        assert!(parse_sa("chr,1,-,10S5M,255,3;").is_ok());
        // strand garbage
        assert!(parse_sa("chr,1,+garbage,100M,60,0;").is_err());
        assert!(parse_sa("chr,1,x,100M,60,0;").is_err());
        // empty rname
        assert!(parse_sa(",1,+,100M,60,0;").is_err());
        // pos not positive / not an integer
        assert!(parse_sa("chr,0,+,100M,60,0;").is_err());
        assert!(parse_sa("chr,-3,+,100M,60,0;").is_err());
        assert!(parse_sa("chr,abc,+,100M,60,0;").is_err());
        // invalid / empty CIGAR
        assert!(parse_sa("chr,1,+,10Z,60,0;").is_err());
        assert!(parse_sa("chr,1,+,,60,0;").is_err());
        // mapq out of range / non-integer
        assert!(parse_sa("chr,1,+,100M,256,0;").is_err());
        assert!(parse_sa("chr,1,+,100M,-1,0;").is_err());
        // NM negative / non-integer
        assert!(parse_sa("chr,1,+,100M,60,-1;").is_err());
        assert!(parse_sa("chr,1,+,100M,60,x;").is_err());
        // not terminated by ';'
        assert!(parse_sa("chr,1,+,100M,60,0").is_err());
        // interior empty entry
        assert!(parse_sa("chr,1,+,100M,60,0;;chr,2,+,100M,60,0;").is_err());
        // wrong field count
        assert!(parse_sa("chr,1,+,100M,60;").is_err());
        // empty payload is not an error, just no entries
        assert!(parse_sa("").unwrap().is_empty());
    }

    fn diag_block(old_start: u64, old_end: u64, new_start: u64, strand: Strand, new_contig: &str) -> Block {
        Block {
            old_contig: "chr".into(),
            old: OldIv::new(old_start, old_end),
            new_contig: new_contig.into(),
            new_start: crate::coords::NewPos(new_start),
            strand,
        }
    }

    fn diff_with(blocks: Vec<Block>, t1: Vec<OldIv>, t2: Vec<OldIv>) -> Diff {
        let mut d = Diff::default();
        d.blocks.insert("chr".into(), blocks);
        d.taint_t1.insert("chr".into(), t1);
        d.taint_t2.insert("chr".into(), t2);
        d
    }

    /// A minimal `MD` consistent with `cigar` and containing no mismatches, so
    /// rule 7d is a no-op unless a test supplies a specific `MD`.
    fn all_match_md(cigar: &[Cigar]) -> String {
        let mut s = String::new();
        let mut run: u32 = 0;
        for c in cigar {
            match c {
                Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => run += *l,
                Cigar::Del(l) => {
                    s.push_str(&run.to_string());
                    run = 0;
                    s.push('^');
                    for _ in 0..*l {
                        s.push('A');
                    }
                }
                _ => {}
            }
        }
        s.push_str(&run.to_string());
        s
    }

    fn info(qname: &str, flag: u16, pos: i64, mapq: u8, cigar: Vec<Cigar>, de: f64, sa_nm: Option<i64>) -> RecInfo {
        // A minimal all-match MD consistent with the CIGAR, so rule 7d is a
        // no-op (small statistic) unless a test overrides `zdrop`.
        let md = all_match_md(&cigar);
        let zdrop = zdrop_statistic(&cigar, Some(&md), &Scoring::default());
        RecInfo {
            qname: qname.as_bytes().to_vec(),
            flag,
            contig: "chr".into(),
            pos,
            mapq,
            cigar,
            seq_len: 20_000,
            de: Some(de),
            nm: sa_nm,
            // Default to a clean competitor ratio so unrelated tests are unaffected.
            s1: Some(100),
            s2: Some(0),
            sa_present: false,
            sa_ok: false,
            sa: Vec::new(),
            unexpected_tag: false,
            xl_translated: false,
            zdrop,
        }
    }

    #[test]
    fn minus_strand_position_math() {
        // old [10,20) mapped to new_start 100 on '-': the old end maps to the
        // lowest new coordinate, so old 19 -> 100 and old 10 -> 109.
        let b = Block {
            old_contig: "c".into(),
            old: OldIv::new(10, 20),
            new_contig: "c".into(),
            new_start: crate::coords::NewPos(100),
            strand: Strand::Minus,
        };
        assert_eq!(b.map_pos(OldPos(19)).unwrap().get(), 100);
        assert_eq!(b.map_pos(OldPos(10)).unwrap().get(), 109);
        // The translated position of a '-' record is map_pos(ref_end - 1).
        let ref_end = OldPos(20);
        assert_eq!(b.map_pos(OldPos(ref_end.get() - 1)).unwrap().get(), 100);
    }

    #[test]
    fn origin_crossing_plus_translatable() {
        // 1000 bp contig, identical mapping: A at [900,1000), B at [0,100).
        let diff = diff_with(
            vec![diag_block(0, 1000, 0, Strand::Plus, "chr")],
            vec![],
            vec![],
        );
        let cfg = TranslateConfig::default();
        // B reaches start (pos 0, leading clip 0), A reaches end.
        let b = info("r", 0, 0, 60, vec![m(100)], 0.001, Some(5));
        let a = info("r", 0x800, 900, 60, vec![m(100)], 0.001, Some(7));
        // SA tags consistent with each other.
        let mut b = b;
        let mut a = a;
        b.sa = parse_sa("chr,901,+,100M,60,7;").unwrap();
        a.sa = parse_sa("chr,1,+,100M,60,5;").unwrap();
        let plan = classify_origin_crossing(&[b, a], 1000, &diff, 0.01, &cfg, &|_| Some(1000)).unwrap();
        assert_eq!(plan.a, 1);
        assert_eq!(plan.b, 0);
    }

    #[test]
    fn origin_crossing_junction_not_preserved() {
        // Valid blocks, but B maps old 0 to new 500 (not new 0) -> junction broken.
        let diff = diff_with(
            vec![
                diag_block(0, 600, 500, Strand::Plus, "chr"),
                diag_block(400, 1000, 0, Strand::Plus, "chr"),
            ],
            vec![],
            vec![],
        );
        let cfg = TranslateConfig::default();
        let mut b = info("r", 0, 0, 60, vec![m(100)], 0.001, Some(5));
        let mut a = info("r", 0x800, 900, 60, vec![m(100)], 0.001, Some(7));
        b.sa = parse_sa("chr,901,+,100M,60,7;").unwrap();
        a.sa = parse_sa("chr,1,+,100M,60,5;").unwrap();
        assert_eq!(
            classify_origin_crossing(&[b, a], 1000, &diff, 0.01, &cfg, &|_| Some(1000)).err(),
            Some(Reason::SplitJunction)
        );
    }

    #[test]
    fn origin_crossing_minus_translatable() {
        // '-' block: old 0 -> new end, old end -> new 0.
        let diff = diff_with(vec![diag_block(0, 1000, 0, Strand::Minus, "chr")], vec![], vec![]);
        let cfg = TranslateConfig::default();
        let mut b = info("r", 16, 0, 60, vec![m(100)], 0.001, Some(5));
        let mut a = info("r", 16 | 0x800, 900, 60, vec![m(100)], 0.001, Some(7));
        b.sa = parse_sa("chr,901,-,100M,60,7;").unwrap();
        a.sa = parse_sa("chr,1,-,100M,60,5;").unwrap();
        assert!(classify_origin_crossing(&[b, a], 1000, &diff, 0.01, &cfg, &|_| Some(1000)).is_ok());
    }

    #[test]
    fn origin_crossing_not_two_records() {
        let diff = diff_with(vec![diag_block(0, 1000, 0, Strand::Plus, "chr")], vec![], vec![]);
        let cfg = TranslateConfig::default();
        let b = info("r", 0, 0, 60, vec![m(100)], 0.001, None);
        let a = info("r", 0x800, 900, 60, vec![m(100)], 0.001, None);
        assert_eq!(
            classify_origin_crossing(&[b, a], 1000, &diff, 0.01, &cfg, &|_| Some(1000)).err(),
            Some(Reason::SplitNotOrigin)
        );
    }

    #[test]
    fn single_record_taint_reason() {
        let diff = diff_with(
            vec![diag_block(0, 1000, 0, Strand::Plus, "chr")],
            vec![],
            vec![OldIv::new(100, 200)],
        );
        let cfg = TranslateConfig::default();
        let r = info("r", 0, 150, 60, vec![m(50)], 0.001, Some(0));
        assert_eq!(
            check_record(&r, 1000, &diff, 0.01, &cfg, 0, &|_| Some(1000)).err(),
            Some(Reason::TaintT2)
        );
    }

    // --- competitor (s2/s1) rule ---------------------------------------------

    /// A plain eligible single record (clean mapq, de and no taint).
    fn competitor_base() -> (RecInfo, Diff) {
        let diff = diff_with(vec![diag_block(0, 1000, 0, Strand::Plus, "chr")], vec![], vec![]);
        let r = info("r", 0, 500, 60, vec![m(50)], 0.001, Some(0));
        (r, diff)
    }

    #[test]
    fn competitor_ratio_below_at_above() {
        let (base, diff) = competitor_base();
        let cfg = TranslateConfig::default(); // max_s2_ratio = 0.5
        let check = |r: &RecInfo| check_record(r, 1000, &diff, 0.01, &cfg, 0, &|_| Some(1000)).err();

        let mut below = base.clone();
        below.s1 = Some(100);
        below.s2 = Some(49);
        assert_eq!(check(&below), None, "49/100 < 0.5 must be eligible");

        let mut at = base.clone();
        at.s1 = Some(100);
        at.s2 = Some(50);
        assert_eq!(check(&at), Some(Reason::Competitor), "50/100 == 0.5 must defer");

        let mut above = base.clone();
        above.s1 = Some(100);
        above.s2 = Some(60);
        assert_eq!(check(&above), Some(Reason::Competitor), "60/100 > 0.5 must defer");
    }

    #[test]
    fn competitor_missing_or_invalid_tags() {
        let (base, diff) = competitor_base();
        let cfg = TranslateConfig::default();
        let check = |r: &RecInfo| check_record(r, 1000, &diff, 0.01, &cfg, 0, &|_| Some(1000)).err();

        let mut no_s1 = base.clone();
        no_s1.s1 = None;
        no_s1.s2 = Some(0);
        assert_eq!(check(&no_s1), Some(Reason::Competitor));

        let mut no_s2 = base.clone();
        no_s2.s1 = Some(100);
        no_s2.s2 = None;
        assert_eq!(check(&no_s2), Some(Reason::Competitor));

        let mut zero_s1 = base.clone();
        zero_s1.s1 = Some(0);
        zero_s1.s2 = Some(0);
        assert_eq!(check(&zero_s1), Some(Reason::Competitor));

        let mut neg_s1 = base.clone();
        neg_s1.s1 = Some(-5);
        neg_s1.s2 = Some(0);
        assert_eq!(check(&neg_s1), Some(Reason::Competitor));

        // Negative s2 is invalid too (it must be >= 0).
        let mut neg_s2 = base.clone();
        neg_s2.s1 = Some(100);
        neg_s2.s2 = Some(-1);
        assert_eq!(check(&neg_s2), Some(Reason::Competitor));

        // Missing both is also a deferral.
        let mut neither = base.clone();
        neither.s1 = None;
        neither.s2 = None;
        assert_eq!(check(&neither), Some(Reason::Competitor));
    }

    #[test]
    fn competitor_lifted_translated_missing_s2_is_eligible_only_when_chain_ok() {
        // A chained input's `XL:A:T` record: the previous lift stripped `s2`.
        // It is eligible only when the contract accepts the all-default chain
        // (`inherited_xl_ok`); otherwise it defers like any missing-`s2` record.
        let (base, diff) = competitor_base();
        let ok_cfg = TranslateConfig { inherited_xl_ok: true, ..Default::default() };
        let no_cfg = TranslateConfig::default();
        let run = |cfg: &TranslateConfig, r: &RecInfo| {
            check_record(r, 1000, &diff, 0.01, cfg, 0, &|_| Some(1000)).err()
        };

        let mut lifted = base.clone();
        lifted.s1 = Some(100);
        lifted.s2 = None;
        assert_eq!(run(&ok_cfg, &lifted), Some(Reason::Competitor), "no XL flag -> defer");
        assert_eq!(run(&no_cfg, &lifted), Some(Reason::Competitor), "chain not ok -> defer");

        lifted.xl_translated = true;
        assert_eq!(run(&no_cfg, &lifted), Some(Reason::Competitor), "non-default chain -> defer");
        assert_eq!(run(&ok_cfg, &lifted), None, "XL:A:T, chain ok, s1 only -> eligible");

        // A lifted record still needs a positive s1.
        let mut no_s1 = lifted.clone();
        no_s1.s1 = None;
        assert_eq!(run(&ok_cfg, &no_s1), Some(Reason::Competitor));

        let mut zero_s1 = lifted.clone();
        zero_s1.s1 = Some(0);
        assert_eq!(run(&ok_cfg, &zero_s1), Some(Reason::Competitor));

        // A lifted record that DOES carry a close competitor still defers.
        let mut close = lifted.clone();
        close.s2 = Some(60);
        assert_eq!(run(&ok_cfg, &close), Some(Reason::Competitor));
    }

    #[test]
    fn eligibility_is_default_flags_every_eligibility_option() {
        let d = TranslateConfig::default();
        assert!(d.eligibility_is_default());
        // Non-eligibility fields must not count.
        let mut ok = d;
        ok.strict = true;
        ok.xl_tag = false;
        ok.scoring = scoring(Preset::MapHifi);
        ok.inherited_xl_ok = true;
        assert!(ok.eligibility_is_default(), "strict/xl_tag/scoring/inherited_xl_ok are not eligibility");

        // Every other option must make it false.
        let muts: Vec<(&str, Box<dyn Fn(&mut TranslateConfig)>)> = vec![
            ("min_mapq", Box::new(|c| c.min_mapq = 59)),
            ("max_clip", Box::new(|c| c.max_clip = 0)),
            ("margin", Box::new(|c| c.margin = 1)),
            ("max_de", Box::new(|c| c.max_de = Some(0.02))),
            ("t2_win", Box::new(|c| c.t2_win = 1)),
            ("t2_step", Box::new(|c| c.t2_step = 1)),
            ("t2_flank", Box::new(|c| c.t2_flank = 1)),
            ("clip_k", Box::new(|c| c.clip_k = 1)),
            ("clip_min_hits", Box::new(|c| c.clip_min_hits = 1)),
            ("clip_min_len", Box::new(|c| c.clip_min_len = 1)),
            ("clip_flank", Box::new(|c| c.clip_flank = 1)),
            ("clip_subtract_old", Box::new(|c| c.clip_subtract_old = true)),
            ("unsafe_no_clip_screen", Box::new(|c| c.unsafe_no_clip_screen = true)),
            ("max_s2_ratio", Box::new(|c| c.max_s2_ratio = 1.0)),
            ("max_minus_indel", Box::new(|c| c.max_minus_indel = 0)),
            ("max_zdrop", Box::new(|c| c.max_zdrop = 0)),
        ];
        for (name, f) in muts {
            let mut c = d;
            f(&mut c);
            assert!(!c.eligibility_is_default(), "{name} must count as non-default");
        }
    }

    #[test]
    fn scan_aux_accepts_xl_and_flags_foreign_tags() {
        // `XL` is cramhook's own tag and is accepted on input.
        let mut rec = bam::Record::new();
        rec.push_aux(b"XL", Aux::Char(b'T')).unwrap();
        assert_eq!(scan_aux(&rec), (false, None));

        let mut rec2 = bam::Record::new();
        rec2.push_aux(b"XL", Aux::Char(b'R')).unwrap();
        assert_eq!(scan_aux(&rec2), (false, None));

        // A tag outside the allow-list is still flagged.
        let mut rec3 = bam::Record::new();
        rec3.push_aux(b"zz", Aux::I32(1)).unwrap();
        assert!(scan_aux(&rec3).0);
    }

    #[test]
    fn strips_inherited_xl_and_adds_fresh_only_when_tagging() {
        // The emit path is `strip_translated_tags` then `add_xl` iff `xl_tag`.
        // An inherited `XL` (from a previous lift) must never survive; the fresh
        // `XL:A:T` is added only when tagging is on.
        fn emitted(xl_tag: bool, inherited: Option<u8>) -> Vec<u8> {
            let mut rec = bam::Record::new();
            rec.push_aux(b"NM", Aux::I32(3)).unwrap(); // a kept neighbour tag
            if let Some(v) = inherited {
                rec.push_aux(b"XL", Aux::Char(v)).unwrap();
            }
            strip_translated_tags(&mut rec);
            if xl_tag {
                add_xl(&mut rec);
            }
            assert!(rec.aux(b"NM").is_ok(), "NM must be kept");
            let mut out = Vec::new();
            for item in rec.aux_iter() {
                if let Ok((tag, aux)) = item {
                    if tag == b"XL" {
                        out.push(match aux {
                            Aux::Char(c) => c,
                            _ => 0,
                        });
                    }
                }
            }
            out
        }
        assert_eq!(emitted(true, Some(b'R')), vec![b'T'], "inherited XL replaced by T");
        assert_eq!(emitted(true, None), vec![b'T']);
        assert_eq!(emitted(false, Some(b'R')), Vec::<u8>::new(), "inherited XL dropped");
        assert_eq!(emitted(false, None), Vec::<u8>::new());
    }

    #[test]
    fn split_group_strips_inherited_xl_per_record() {
        // `process_group` strips every record of a split pair before emitting.
        let mut group = vec![bam::Record::new(), bam::Record::new()];
        for r in group.iter_mut() {
            r.push_aux(b"XL", Aux::Char(b'R')).unwrap();
            r.push_aux(b"s2", Aux::I32(7)).unwrap();
        }
        for r in group.iter_mut() {
            strip_translated_tags(r);
        }
        for r in group.iter_mut() {
            add_xl(r);
        }
        for r in &group {
            assert!(r.aux(b"s2").is_err(), "s2 must be stripped");
            let mut n = 0usize;
            for item in r.aux_iter() {
                if let Ok((tag, aux)) = item {
                    if tag == b"XL" {
                        n += 1;
                        assert_eq!(aux, Aux::Char(b'T'));
                    }
                }
            }
            assert_eq!(n, 1, "exactly one fresh XL");
        }
        // Tagging off: no XL survives on either record.
        let mut group = vec![bam::Record::new(), bam::Record::new()];
        for r in group.iter_mut() {
            r.push_aux(b"XL", Aux::Char(b'R')).unwrap();
            strip_translated_tags(r);
        }
        for r in &group {
            assert!(r.aux(b"XL").is_err(), "inherited XL survived with tagging off");
        }
    }

    #[test]
    fn competitor_tags_strict_parsing() {
        // Integer aux types are accepted; float/double/string/char are rejected
        // (-> None), so `to_info` makes the record defer.
        let mut rec = bam::Record::new();
        rec.push_aux(b"s1", Aux::I32(100)).unwrap();
        rec.push_aux(b"s2", Aux::Float(0.4)).unwrap();
        assert_eq!(aux_int_strict(&rec, b"s1"), Some(100));
        assert_eq!(aux_int_strict(&rec, b"s2"), None, "float s2 must be rejected");
        assert_eq!(aux_int_strict(&rec, b"zz"), None, "missing tag");

        let mut rec2 = bam::Record::new();
        rec2.push_aux(b"s2", Aux::String("5")).unwrap();
        assert_eq!(aux_int_strict(&rec2, b"s2"), None, "string s2 must be rejected");

        let mut rec3 = bam::Record::new();
        rec3.push_aux(b"s1", Aux::U16(250)).unwrap();
        assert_eq!(aux_int_strict(&rec3, b"s1"), Some(250));
    }

    #[test]
    fn competitor_split_only_supplementary_fails() {
        // Valid origin-crossing pair; only the supplementary record (A) has a
        // close competitor, so the unit is deferred with `split_competitor`.
        let diff = diff_with(vec![diag_block(0, 1000, 0, Strand::Plus, "chr")], vec![], vec![]);
        let cfg = TranslateConfig::default();
        let mut b = info("r", 0, 0, 60, vec![m(100)], 0.001, Some(5));
        let mut a = info("r", 0x800, 900, 60, vec![m(100)], 0.001, Some(7));
        b.sa = parse_sa("chr,901,+,100M,60,7;").unwrap();
        a.sa = parse_sa("chr,1,+,100M,60,5;").unwrap();
        // B is clean; A has a near-equal competitor.
        b.s1 = Some(100);
        b.s2 = Some(0);
        a.s1 = Some(100);
        a.s2 = Some(80);
        assert_eq!(
            classify_origin_crossing(&[b, a], 1000, &diff, 0.01, &cfg, &|_| Some(1000)).err(),
            Some(Reason::SplitCompetitor)
        );
    }

    // --- minus-indel rule (7b) ----------------------------------------------

    fn minus_block() -> Diff {
        diff_with(vec![diag_block(0, 1000, 0, Strand::Minus, "chr")], vec![], vec![])
    }

    fn plus_block() -> Diff {
        diff_with(vec![diag_block(0, 1000, 0, Strand::Plus, "chr")], vec![], vec![])
    }

    #[test]
    fn minus_indel_large_on_plus_is_translated() {
        // '+' blocks are unaffected by the rule.
        let diff = plus_block();
        let cfg = TranslateConfig::default(); // max_minus_indel = 50
        let r = info("r", 0, 500, 60, vec![m(50), d(100), m(50)], 0.001, Some(0));
        assert!(check_record(&r, 1000, &diff, 0.01, &cfg, 0, &|_| Some(1000)).is_ok());
    }

    #[test]
    fn minus_indel_49_vs_50() {
        let diff = minus_block();
        let cfg = TranslateConfig::default(); // max_minus_indel = 50
        let check = |cigar: Vec<Cigar>| {
            let r = info("r", 0, 500, 60, cigar, 0.001, Some(0));
            check_record(&r, 1000, &diff, 0.01, &cfg, 0, &|_| Some(1000)).err()
        };
        // 49 bp deletion is below the threshold -> eligible.
        assert_eq!(check(vec![m(50), d(49), m(50)]), None);
        // 50 bp deletion is at the threshold -> defer.
        assert_eq!(check(vec![m(50), d(50), m(50)]), Some(Reason::MinusIndel));
        // A large insertion is caught too.
        assert_eq!(check(vec![m(50), i(50), m(50)]), Some(Reason::MinusIndel));
        // No indel -> eligible.
        assert_eq!(check(vec![m(100)]), None);
    }

    #[test]
    fn minus_indel_split_minus_indel() {
        // '-' origin-crossing pair; the record reaching the start (B) has a large
        // deletion, so the unit is deferred with `split_minus_indel`.
        let diff = minus_block();
        let cfg = TranslateConfig::default();
        let mut b = info("r", 16, 0, 60, vec![m(50), d(50), m(50)], 0.001, Some(5));
        let mut a = info("r", 16 | 0x800, 900, 60, vec![m(100)], 0.001, Some(7));
        // SA tags describe the other record's (simplified) CIGAR.
        b.sa = parse_sa("chr,901,-,100M,60,7;").unwrap();
        a.sa = parse_sa("chr,1,-,100M50D,60,5;").unwrap();
        assert_eq!(
            classify_origin_crossing(&[b, a], 1000, &diff, 0.01, &cfg, &|_| Some(1000)).err(),
            Some(Reason::SplitMinusIndel)
        );
    }

    // --- zdrop rule (7d) ----------------------------------------------------

    #[test]
    fn zdrop_statistic_hand_computed() {
        let sc = Scoring::default(); // lr:hq: a=2 b=4 q=4 e=2 q2=24 e2=1
        // No indel -> no drop.
        assert_eq!(zdrop_statistic(&[m(100)], Some("100"), &sc), Some(0));
        // 5 bp insertion after 100 matches: -min(4+10, 24+5) = -14, and the
        // diagonal penalty e2*|diag| = 1*5 reduces the drop to 14 - 5 = 9.
        assert_eq!(zdrop_statistic(&[m(100), i(5), m(100)], Some("200"), &sc), Some(9));
        // 3 bp deletion: -min(4+6, 24+3) = -10, diagonal penalty 1*3 -> 7.
        assert_eq!(zdrop_statistic(&[m(100), d(3), m(100)], Some("100^AAA100"), &sc), Some(7));
        // 300 matches then 110 mismatches: peak 600, final 600-440 -> drop 440.
        let md = format!("300{}", "A0".repeat(110));
        assert_eq!(zdrop_statistic(&[m(410)], Some(&md), &sc), Some(440));
    }

    /// Port of `zdrop` from scratch/m7-mismatch/signal.py::stats on a few
    /// hand-built CIGAR/MD examples (see scratch/worker-zdrop/crosscheck.py,
    /// which prints the expected `zdrop` and `MD` for each).
    #[test]
    fn zdrop_statistic_matches_python_reference() {
        let sc = Scoring::default();
        // crosscheck.py case A: 100M5I100M, all match -> 9.
        assert_eq!(zdrop_statistic(&[m(100), i(5), m(100)], Some("200"), &sc), Some(9));
        // case B: 100M3D100M, all match -> 7.
        assert_eq!(zdrop_statistic(&[m(100), d(3), m(100)], Some("100^AAA100"), &sc), Some(7));
        // case C: 100M2D100M with mismatches at read 10/40/150 -> 6.
        assert_eq!(
            zdrop_statistic(&[m(100), d(2), m(100)], Some("10A29A59^CC50A49"), &sc),
            Some(6)
        );
        // case D: 300 matches then 110 mismatches -> 440 (the strict MD form
        // "300A0A0...", as printed by crosscheck.py).
        let md = format!("300{}", "A0".repeat(110));
        assert_eq!(zdrop_statistic(&[m(410)], Some(&md), &sc), Some(440));
    }

    /// The `MD` grammar `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*` is enforced: uppercase
    /// letters, numeric runs at both ends and between tokens.
    #[test]
    fn zdrop_md_grammar_is_strict() {
        let sc = Scoring::default();
        // Valid forms.
        assert!(zdrop_statistic(&[m(11)], Some("10A0"), &sc).is_some());
        assert!(zdrop_statistic(&[m(20)], Some("9A10"), &sc).is_some());
        assert!(zdrop_statistic(&[m(10), d(1), m(1)], Some("10^A1"), &sc).is_some());
        assert!(zdrop_statistic(&[m(10)], Some("10"), &sc).is_some());
        // A trailing numeric run is required.
        assert_eq!(zdrop_statistic(&[m(10)], Some("9A"), &sc), None);
        assert_eq!(zdrop_statistic(&[m(10), d(1), m(1)], Some("10^A"), &sc), None);
        // A leading numeric run is required.
        assert_eq!(zdrop_statistic(&[m(10)], Some("A0"), &sc), None);
        assert_eq!(zdrop_statistic(&[m(10)], Some("^A0"), &sc), None);
        // Uppercase only.
        assert_eq!(zdrop_statistic(&[m(10)], Some("9a0"), &sc), None);
        assert_eq!(zdrop_statistic(&[m(10), d(1), m(1)], Some("10^a1"), &sc), None);
        // No stray bytes, no empty `^` run.
        assert_eq!(zdrop_statistic(&[m(10)], Some("10!"), &sc), None);
        assert_eq!(zdrop_statistic(&[m(10)], Some("10^0"), &sc), None);
        // A `^` run must line up with a `D` op.
        assert_eq!(zdrop_statistic(&[m(10)], Some("10^A0"), &sc), None);
        // A mismatch letter must fall inside an aligned op, not a `D` op.
        assert_eq!(zdrop_statistic(&[m(10), d(1)], Some("10A0"), &sc), None);
    }

    /// `=` columns must be `MD` matches and `X` columns `MD` mismatches.
    #[test]
    fn zdrop_equal_diff_cigar_consistency() {
        let sc = Scoring::default();
        assert_eq!(zdrop_statistic(&[Cigar::Equal(10)], Some("10"), &sc), Some(0));
        assert_eq!(zdrop_statistic(&[Cigar::Equal(10)], Some("9A0"), &sc), None);
        assert_eq!(zdrop_statistic(&[Cigar::Diff(1)], Some("0A0"), &sc), Some(0));
        assert_eq!(zdrop_statistic(&[Cigar::Diff(1)], Some("1"), &sc), None);
        // `M` accepts either.
        assert_eq!(zdrop_statistic(&[m(11)], Some("10A0"), &sc), Some(4));
    }

    #[test]
    fn zdrop_noisy_stretch_minus_defers_plus_translates() {
        let cfg = TranslateConfig::default(); // max_zdrop = 400
        let md = format!("300{}", "A0".repeat(110)); // statistic 440 > 400
        let mut r = info("r", 0, 500, 60, vec![m(410)], 0.001, Some(0));
        r.zdrop = zdrop_statistic(&r.cigar, Some(&md), &Scoring::default());
        assert_eq!(r.zdrop, Some(440));
        // '+' block: rule 7d does not apply.
        assert!(check_record(&r, 1000, &plus_block(), 0.01, &cfg, 0, &|_| Some(1000)).is_ok());
        // '-' block: the same record is deferred.
        assert_eq!(
            check_record(&r, 1000, &minus_block(), 0.01, &cfg, 0, &|_| Some(1000)).err(),
            Some(Reason::Zdrop)
        );
    }

    #[test]
    fn zdrop_threshold_boundary() {
        let mut cfg = TranslateConfig::default();
        cfg.max_zdrop = 8; // 4 * (number of mismatches after the peak)
        let diff = minus_block();
        let mut at = info("r", 0, 500, 60, vec![m(12)], 0.001, Some(0));
        at.zdrop = Some(8); // == threshold
        assert!(check_record(&at, 1000, &diff, 0.01, &cfg, 0, &|_| Some(1000)).is_ok());
        let mut over = info("r", 0, 500, 60, vec![m(13)], 0.001, Some(0));
        over.zdrop = Some(12); // > threshold
        assert_eq!(
            check_record(&over, 1000, &diff, 0.01, &cfg, 0, &|_| Some(1000)).err(),
            Some(Reason::Zdrop)
        );
    }

    #[test]
    fn zdrop_missing_md_fails_safe() {
        let cfg = TranslateConfig::default();
        let mut r = info("r", 0, 500, 60, vec![m(50)], 0.001, Some(0));
        r.zdrop = None; // missing/malformed/inconsistent MD
        // '-' block: a missing statistic defers fail-safe.
        assert_eq!(
            check_record(&r, 1000, &minus_block(), 0.01, &cfg, 0, &|_| Some(1000)).err(),
            Some(Reason::Zdrop)
        );
        // '+' block: unaffected.
        assert!(check_record(&r, 1000, &plus_block(), 0.01, &cfg, 0, &|_| Some(1000)).is_ok());
        // An MD inconsistent with the CIGAR is also `None` (checked directly).
        assert_eq!(zdrop_statistic(&[m(50)], Some("49"), &Scoring::default()), None);
        let mut bad = info("r", 0, 500, 60, vec![m(50)], 0.001, Some(0));
        bad.zdrop = None;
        assert_eq!(
            check_record(&bad, 1000, &minus_block(), 0.01, &cfg, 0, &|_| Some(1000)).err(),
            Some(Reason::Zdrop)
        );
    }

    #[test]
    fn zdrop_split_minus_zdrop() {
        // '-' origin-crossing pair; the record reaching the start (B) has a noisy
        // stretch, so the unit is deferred with `split_zdrop`.
        let diff = minus_block();
        let cfg = TranslateConfig::default();
        let mut b = info("r", 16, 0, 60, vec![m(410)], 0.001, Some(5));
        b.zdrop = Some(440); // statistic 440 > 400
        let mut a = info("r", 16 | 0x800, 900, 60, vec![m(100)], 0.001, Some(7));
        b.sa = parse_sa("chr,901,-,100M,60,7;").unwrap();
        a.sa = parse_sa("chr,1,-,410M,60,5;").unwrap();
        assert_eq!(
            classify_origin_crossing(&[b, a], 1000, &diff, 0.01, &cfg, &|_| Some(1000)).err(),
            Some(Reason::SplitZdrop)
        );
    }

    // --- junction preservation ----------------------------------------------

    #[test]
    fn junction_preserved_both_strands() {
        let plus = diag_block(0, 1000, 0, Strand::Plus, "chr");
        assert!(junction_preserved(&plus, 1000, Boundary::Start, 1000));
        assert!(junction_preserved(&plus, 1000, Boundary::End, 1000));
        // '+' start/end map to new start/end only if the block spans the contig.
        assert!(!junction_preserved(&diag_block(0, 900, 0, Strand::Plus, "chr"), 1000, Boundary::End, 1000));

        let minus = diag_block(0, 1000, 0, Strand::Minus, "chr");
        // '-' maps old 0 to the new end and old end to new 0.
        assert!(junction_preserved(&minus, 1000, Boundary::Start, 1000));
        assert!(junction_preserved(&minus, 1000, Boundary::End, 1000));
        // A block shifted in new coordinates breaks the junction.
        assert!(!junction_preserved(&diag_block(0, 1000, 500, Strand::Plus, "chr"), 1000, Boundary::Start, 1000));
    }

    #[test]
    fn check_record_junction_start() {
        let cfg = TranslateConfig::default();
        // crossing the start (left clip 50 + margin 500) with a block that maps
        // old 0 to new 500 -> junction not preserved.
        let bad = diff_with(vec![diag_block(0, 1000, 500, Strand::Plus, "chr")], vec![], vec![]);
        let r = info("r", 0, 0, 60, vec![s(50), m(100)], 0.001, Some(0));
        assert_eq!(
            check_record(&r, 1000, &bad, 0.01, &cfg, 50, &|_| Some(1000)).err(),
            Some(Reason::Junction)
        );
        let good = diff_with(vec![diag_block(0, 1000, 0, Strand::Plus, "chr")], vec![], vec![]);
        assert!(check_record(&r, 1000, &good, 0.01, &cfg, 50, &|_| Some(1000)).is_ok());
    }

    #[test]
    fn check_record_junction_end() {
        let cfg = TranslateConfig::default();
        let r = info("r", 0, 900, 60, vec![m(100), s(50)], 0.001, Some(0));
        let bad = diff_with(vec![diag_block(0, 1000, 100, Strand::Plus, "chr")], vec![], vec![]);
        assert_eq!(
            check_record(&r, 1000, &bad, 0.01, &cfg, 50, &|_| Some(1000)).err(),
            Some(Reason::Junction)
        );
        let good = diff_with(vec![diag_block(0, 1000, 0, Strand::Plus, "chr")], vec![], vec![]);
        assert!(check_record(&r, 1000, &good, 0.01, &cfg, 50, &|_| Some(1000)).is_ok());
    }

    #[test]
    fn check_record_boundary_crossing_no_block_is_no_block() {
        let cfg = TranslateConfig::default();
        // The unclipped window crosses the contig end, but no block covers it:
        // the reason must be `no_block`, not `junction`.
        let diff = diff_with(vec![], vec![], vec![]);
        let r = info("r", 0, 900, 60, vec![m(100), s(50)], 0.001, Some(0));
        assert_eq!(
            check_record(&r, 1000, &diff, 0.01, &cfg, 50, &|_| Some(1000)).err(),
            Some(Reason::NoBlock)
        );
        // A non-crossing uncovered window is `no_block` too.
        let r2 = info("r", 0, 500, 60, vec![m(50)], 0.001, Some(0));
        assert_eq!(
            check_record(&r2, 1000, &diff, 0.01, &cfg, 0, &|_| Some(1000)).err(),
            Some(Reason::NoBlock)
        );
    }

    #[test]
    fn check_record_boundary_crossing_unpreserved_junction_is_junction() {
        let cfg = TranslateConfig::default();
        // A block covers the window but maps old 0 to new 500: junction broken.
        let diff = diff_with(vec![diag_block(0, 1000, 500, Strand::Plus, "chr")], vec![], vec![]);
        let r = info("r", 0, 900, 60, vec![m(100), s(50)], 0.001, Some(0));
        assert_eq!(
            check_record(&r, 1000, &diff, 0.01, &cfg, 50, &|_| Some(1000)).err(),
            Some(Reason::Junction)
        );
    }

    #[test]
    fn origin_crossing_no_block_is_split_no_block() {
        // Origin-crossing pair with no identity block: `split_no_block`, not
        // `split_junction` (the latter needs a covering block).
        let diff = diff_with(vec![], vec![], vec![]);
        let cfg = TranslateConfig::default();
        let mut b = info("r", 0, 0, 60, vec![m(100)], 0.001, Some(5));
        let mut a = info("r", 0x800, 900, 60, vec![m(100)], 0.001, Some(7));
        b.sa = parse_sa("chr,901,+,100M,60,7;").unwrap();
        a.sa = parse_sa("chr,1,+,100M,60,5;").unwrap();
        assert_eq!(
            classify_origin_crossing(&[b, a], 1000, &diff, 0.01, &cfg, &|_| Some(1000)).err(),
            Some(Reason::SplitNoBlock)
        );
    }

    #[test]
    fn clip_screen_rejects_bad_k() {
        let tmp = std::env::temp_dir().join("cramhook_clip_k_test.fa");
        std::fs::write(&tmp, b">c\nACGTACGTACGTACGTACGT\n").unwrap();
        let fa = Fasta::load(&tmp).unwrap();
        let diff = Diff::default();
        assert!(ClipScreen::build(&fa, &fa, &diff, 0, 2, 20, 0, false).is_err());
        assert!(ClipScreen::build(&fa, &fa, &diff, 33, 2, 20, 0, false).is_err());
        assert!(ClipScreen::build(&fa, &fa, &diff, 1, 2, 20, 0, false).unwrap().is_some());
        assert!(ClipScreen::build(&fa, &fa, &diff, 32, 2, 20, 0, false).unwrap().is_some());
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn config_validation_bounds() {
        let mut cfg = TranslateConfig::default();
        assert!(cfg.validate().is_ok());
        cfg.clip_k = 0;
        assert!(cfg.validate().is_err());
        cfg.clip_k = 33;
        assert!(cfg.validate().is_err());
        cfg.clip_k = 32;
        assert!(cfg.validate().is_ok());
        cfg.clip_k = 17;
        for bad in [Some(-0.1), Some(f64::NAN), Some(f64::INFINITY)] {
            cfg.max_de = bad;
            assert!(cfg.validate().is_err(), "max_de={bad:?}");
        }
        cfg.max_de = Some(0.02);
        assert!(cfg.validate().is_ok());
        cfg.t2_win = 0;
        assert!(cfg.validate().is_err());
        cfg.t2_win = 2000;
        cfg.t2_step = 3000;
        assert!(cfg.validate().is_err());
        cfg.t2_step = 0;
        assert!(cfg.validate().is_err());
        cfg.t2_step = 2000;
        assert!(cfg.validate().is_ok());
        cfg.max_zdrop = -1;
        assert!(cfg.validate().is_err());
        cfg.max_zdrop = 0; // 0 is a legitimate very strict value
        assert!(cfg.validate().is_ok());
    }

    #[test]
    fn clip_screen_basic() {
        let k = 5;
        let mut kmers: KmerSet = KmerSet::default();
        for w in b"AACCGGTT".windows(k) {
            kmers.insert(canonical_kmer(w, k).unwrap());
        }
        let s = ClipScreen { k, min_hits: 2, min_len: 5, kmers };
        assert!(s.hits(b"AACCGG")); // AACCG, ACCGG both shared
        assert!(s.hits(b"CCGGTT"));
        assert!(!s.hits(b"AACCG")); // only one shared k-mer (< 2 hits)
        assert!(!s.hits(b"TTTTTT"));
        assert!(!s.hits(b"NNNNNN")); // non-ACGT k-mers are skipped
        // Canonical encoding makes the test strand-invariant.
        assert_eq!(canonical_kmer(b"AAAC", 4), canonical_kmer(&revcomp(b"AAAC"), 4));
    }

    #[test]
    fn max_clip_zero_is_unlimited() {
        let diff = diff_one_block();
        let r = info("r", 0, 500, 60, vec![s(500), m(50)], 0.001, Some(0));
        // Default cap rejects the 500 bp clip.
        let capped = TranslateConfig::default();
        assert_eq!(
            check_record(&r, 1000, &diff, 0.01, &capped, 500, &|_| Some(1000)).err(),
            Some(Reason::Clip)
        );
        // 0 means unlimited.
        let mut unlimited = TranslateConfig::default();
        unlimited.max_clip = 0;
        assert!(check_record(&r, 1000, &diff, 0.01, &unlimited, 500, &|_| Some(1000)).is_ok());
    }

    // --- indel left-alignment ------------------------------------------------

    #[test]
    fn left_align_homopolymer_del() {
        let refs = b"GGGAAAATTT"; // 10
        let read = b"GGGAAATTT"; // 9
        let ops = vec![m(6), d(1), m(3)];
        assert_eq!(left_align_cigar(&ops, read, refs), vec![m(3), d(1), m(6)]);
    }

    #[test]
    fn left_align_homopolymer_ins() {
        let refs = b"GGGAAATTT"; // 9
        let read = b"GGGAAAATTT"; // 10
        let ops = vec![m(6), i(1), m(3)];
        assert_eq!(left_align_cigar(&ops, read, refs), vec![m(3), i(1), m(6)]);
    }

    #[test]
    fn left_align_dinucleotide_repeat() {
        let refs = b"GGACACACGG"; // 10
        let read = b"GGACACGG"; // 8
        let ops = vec![m(6), d(2), m(2)];
        assert_eq!(left_align_cigar(&ops, read, refs), vec![m(2), d(2), m(6)]);
    }

    #[test]
    fn left_align_adjacent_indels_do_not_cross() {
        let read = b"AAAA";
        let refs = b"AAAA";
        let ops = vec![m(2), i(1), d(1), m(1)];
        assert_eq!(left_align_cigar(&ops, read, refs), ops.clone());
    }

    #[test]
    fn left_align_indel_at_match_boundary() {
        // A deletion after a 1-base match run consumes it and cannot cross the
        // leading clip.
        let read = b"AAAA";
        let refs = b"AAAAA";
        let ops = vec![s(2), m(1), d(1), m(3)];
        assert_eq!(left_align_cigar(&ops, read, refs), vec![s(2), d(1), m(4)]);
    }

    #[test]
    fn left_align_noop() {
        let refs = b"ACGTACGT";
        let read = b"ACGACGT"; // delete the T at ref offset 3
        let ops = vec![m(3), d(1), m(4)];
        assert_eq!(left_align_cigar(&ops, read, refs), ops.clone());
    }

    // --- de validation -------------------------------------------------------

    fn diff_one_block() -> Diff {
        diff_with(vec![diag_block(0, 1000, 0, Strand::Plus, "chr")], vec![], vec![])
    }

    #[test]
    fn de_must_be_finite_and_nonnegative() {
        let diff = diff_one_block();
        let cfg = TranslateConfig::default();
        // pos far from the boundary so the junction check is not involved
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.001] {
            let mut r = info("r", 0, 500, 60, vec![m(50)], 0.001, Some(0));
            r.de = Some(bad);
            assert_eq!(
                check_record(&r, 1000, &diff, 0.01, &cfg, 0, &|_| Some(1000)).err(),
                Some(Reason::De),
                "de={bad} must be rejected"
            );
        }
        let r = info("r", 0, 500, 60, vec![m(50)], 0.001, Some(0));
        assert!(check_record(&r, 1000, &diff, 0.01, &cfg, 0, &|_| Some(1000)).is_ok());
    }

    // --- SA consistency ------------------------------------------------------

    #[test]
    fn origin_crossing_sa_cigar_must_match_other_record() {
        let diff = diff_with(vec![diag_block(0, 1000, 0, Strand::Plus, "chr")], vec![], vec![]);
        let cfg = TranslateConfig::default();
        let mut b = info("r", 0, 0, 60, vec![m(100)], 0.001, Some(5));
        let mut a = info("r", 0x800, 900, 60, vec![m(100)], 0.001, Some(7));
        b.sa = parse_sa("chr,901,+,100M,60,7;").unwrap();
        // Tampered SA CIGAR for A: `1M` instead of `100M`.
        a.sa = parse_sa("chr,1,+,1M,60,5;").unwrap();
        assert_eq!(
            classify_origin_crossing(&[b, a], 1000, &diff, 0.01, &cfg, &|_| Some(1000)).err(),
            Some(Reason::SplitSaInconsistent)
        );
    }

    // --- rule 7d: per-base reference + property test -------------------------

    #[derive(Clone, Copy)]
    enum RefMdItem {
        Match(i64),
        Mismatch,
        Del(i64),
    }

    /// Strict `MD` parser used only by the per-base reference (mirrors the
    /// production grammar).
    fn ref_parse_md(md: &str) -> Option<Vec<RefMdItem>> {
        let b = md.as_bytes();
        if b.is_empty() {
            return None;
        }
        let mut i = 0usize;
        let mut out = Vec::new();
        loop {
            if i >= b.len() || !b[i].is_ascii_digit() {
                return None;
            }
            let mut n: i64 = 0;
            while i < b.len() && b[i].is_ascii_digit() {
                n = n.checked_mul(10)?.checked_add((b[i] - b'0') as i64)?;
                i += 1;
            }
            out.push(RefMdItem::Match(n));
            if i >= b.len() {
                break;
            }
            let c = b[i];
            if c == b'^' {
                i += 1;
                let start = i;
                while i < b.len() && b[i].is_ascii_uppercase() {
                    i += 1;
                }
                if i == start {
                    return None;
                }
                out.push(RefMdItem::Del((i - start) as i64));
            } else if c.is_ascii_uppercase() {
                i += 1;
                out.push(RefMdItem::Mismatch);
            } else {
                return None;
            }
        }
        Some(out)
    }

    /// Next aligned `MD` column: `Some(true)` = mismatch, `Some(false)` = match.
    fn ref_next_column(items: &[RefMdItem], idx: &mut usize, rem: &mut i64) -> Option<bool> {
        loop {
            if *rem > 0 {
                *rem -= 1;
                return Some(false);
            }
            match items.get(*idx) {
                Some(RefMdItem::Match(n)) => {
                    *rem = *n;
                    *idx += 1;
                }
                Some(RefMdItem::Mismatch) => {
                    *idx += 1;
                    return Some(true);
                }
                _ => return None,
            }
        }
    }

    /// Independent per-column recurrence for the reference (its own scalar
    /// state; deliberately does not call production `ZdState`/`gap_cost`).
    struct RefState {
        s: i64,
        smax: i64,
        diag_at_max: i64,
        diag: i64,
        e2: i64,
        drop: i64,
    }

    impl RefState {
        fn new(e2: i64) -> Self {
            Self { s: 0, smax: i64::MIN, diag_at_max: 0, diag: 0, e2, drop: 0 }
        }

        fn col(&mut self, score: i64, ddiag: i64) {
            self.s += score;
            self.diag += ddiag;
            if self.s >= self.smax {
                self.smax = self.s;
                self.diag_at_max = self.diag;
            }
            let drop = self.smax - self.s - self.e2 * (self.diag - self.diag_at_max).abs();
            if drop > self.drop {
                self.drop = drop;
            }
        }
    }

    /// Gap penalty `min(q + l·e, q2 + l·e2)`, written out in the test.
    fn ref_gap_cost(l: i64, sc: &Scoring) -> i64 {
        let (q, e, q2, e2) = (sc.q, sc.e, sc.q2, sc.e2);
        (q + l * e).min(q2 + l * e2)
    }

    /// The original O(aligned bases) implementation, kept as the test oracle.
    fn zdrop_statistic_perbase(cigar: &[Cigar], md: Option<&str>, sc: &Scoring) -> Option<i64> {
        let items = ref_parse_md(md?)?;
        let mut idx = 0usize;
        let mut rem: i64 = 0;
        let mut st = RefState::new(sc.e2);
        for c in cigar {
            match *c {
                Cigar::Match(l) => {
                    for _ in 0..l {
                        let mismatch = ref_next_column(&items, &mut idx, &mut rem)?;
                        st.col(if mismatch { -sc.b } else { sc.a }, 0);
                    }
                }
                Cigar::Equal(l) => {
                    for _ in 0..l {
                        if ref_next_column(&items, &mut idx, &mut rem)? {
                            return None;
                        }
                        st.col(sc.a, 0);
                    }
                }
                Cigar::Diff(l) => {
                    for _ in 0..l {
                        if !ref_next_column(&items, &mut idx, &mut rem)? {
                            return None;
                        }
                        st.col(-sc.b, 0);
                    }
                }
                Cigar::Ins(l) => {
                    let l = l as i64;
                    st.col(-ref_gap_cost(l, sc), l);
                }
                Cigar::Del(l) => {
                    let l = l as i64;
                    if rem > 0 {
                        return None;
                    }
                    while matches!(items.get(idx), Some(RefMdItem::Match(0))) {
                        idx += 1;
                    }
                    match items.get(idx) {
                        Some(RefMdItem::Del(n)) if *n == l => idx += 1,
                        _ => return None,
                    }
                    st.col(-ref_gap_cost(l, sc), -l);
                }
                Cigar::SoftClip(_)
                | Cigar::HardClip(_)
                | Cigar::Pad(_)
                | Cigar::RefSkip(_) => {}
            }
        }
        if rem > 0 {
            return None;
        }
        if !items[idx..].iter().all(|it| matches!(it, RefMdItem::Match(0))) {
            return None;
        }
        Some(st.drop)
    }

    fn xorshift(state: &mut u64) -> u64 {
        let mut x = *state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        *state = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }

    /// A random CIGAR plus a consistent strict `MD` string.
    fn random_md_case(state: &mut u64) -> (Vec<Cigar>, String) {
        const BASES: &[u8] = b"ACGT";
        let mut cigar = Vec::new();
        let mut md = String::new();
        let mut run: i64 = 0;
        let n_ops = 1 + (xorshift(state) % 6) as usize;
        for _ in 0..n_ops {
            match xorshift(state) % 8 {
                0..=3 => {
                    let l = 1 + (xorshift(state) % 6) as u32;
                    for _ in 0..l {
                        if xorshift(state) % 3 == 0 {
                            md.push_str(&run.to_string());
                            run = 0;
                            md.push(BASES[(xorshift(state) % 4) as usize] as char);
                        } else {
                            run += 1;
                        }
                    }
                    cigar.push(Cigar::Match(l));
                }
                4 => {
                    let l = 1 + (xorshift(state) % 6) as u32;
                    run += l as i64;
                    cigar.push(Cigar::Equal(l));
                }
                5 => {
                    let l = 1 + (xorshift(state) % 4) as u32;
                    for _ in 0..l {
                        md.push_str(&run.to_string());
                        run = 0;
                        md.push(BASES[(xorshift(state) % 4) as usize] as char);
                    }
                    cigar.push(Cigar::Diff(l));
                }
                6 => {
                    let l = 1 + (xorshift(state) % 3) as u32;
                    cigar.push(Cigar::Ins(l));
                }
                _ => {
                    let l = 1 + (xorshift(state) % 3) as u32;
                    md.push_str(&run.to_string());
                    run = 0;
                    md.push('^');
                    for _ in 0..l {
                        md.push(BASES[(xorshift(state) % 4) as usize] as char);
                    }
                    cigar.push(Cigar::Del(l));
                }
            }
        }
        md.push_str(&run.to_string());
        (cigar, md)
    }

    /// Corrupt an `MD` string to exercise the rejection paths.
    fn corrupt_md(md: &str, state: &mut u64) -> String {
        let mut b = md.as_bytes().to_vec();
        if b.is_empty() {
            return md.to_string();
        }
        match xorshift(state) % 3 {
            0 => {
                let i = (xorshift(state) as usize) % b.len();
                b[i] = b"ACGT^1a"[(xorshift(state) as usize) % 7];
            }
            1 => {
                let n = (xorshift(state) as usize) % b.len();
                b.truncate(n);
            }
            _ => {
                let i = (xorshift(state) as usize) % b.len();
                b[i] = b'A';
            }
        }
        String::from_utf8(b).unwrap()
    }

    #[test]
    fn zdrop_event_matches_perbase_on_random_cases() {
        let scorings =
            [Scoring::default(), Scoring { a: 1, b: 4, q: 6, e: 2, q2: 26, e2: 1 }];
        let mut state: u64 = 0x9E3779B97F4A7C15;
        let mut n_valid = 0usize;
        for case in 0..2000 {
            let (cigar, md) = random_md_case(&mut state);
            let md =
                if xorshift(&mut state) % 5 == 0 { corrupt_md(&md, &mut state) } else { md };
            for sc in &scorings {
                let fast = zdrop_statistic(&cigar, Some(&md), sc);
                let slow = zdrop_statistic_perbase(&cigar, Some(&md), sc);
                assert_eq!(fast, slow, "case {case} md={md:?} cigar={cigar:?}");
                if fast.is_some() {
                    n_valid += 1;
                }
            }
        }
        assert!(n_valid > 1000, "too few valid cases exercised: {n_valid}");
    }

    /// The raw variable-length data buffer (qname + cigar + seq + qual + aux).
    fn data_bytes(rec: &bam::Record) -> Vec<u8> {
        // SAFETY: `l_data` bytes are initialised by any producer of a record.
        unsafe {
            std::slice::from_raw_parts(rec.inner().data, rec.inner().l_data as usize).to_vec()
        }
    }

    fn core_debug(rec: &bam::Record) -> String {
        format!("{:?}", rec.inner().core)
    }

    fn random_cigar(state: &mut u64) -> Vec<Cigar> {
        let nops = (xorshift(state) % 6) as usize + 1;
        let mut out = Vec::with_capacity(nops);
        for _ in 0..nops {
            let len = (xorshift(state) % 20) as u32 + 1;
            out.push(match xorshift(state) % 6 {
                0 => Cigar::Match(len),
                1 => Cigar::Ins(len),
                2 => Cigar::Del(len),
                3 => Cigar::SoftClip(len),
                4 => Cigar::RefSkip(len),
                _ => Cigar::Equal(len),
            });
        }
        out
    }

    /// `rewrite_in_place` must produce exactly the same bytes as rebuilding the
    /// record through `Record::set`, on random records (odd/even lengths, `=`,
    /// IUPAC, missing 0xff qual, varied CIGARs).
    #[test]
    fn rewrite_in_place_matches_set_on_random_records() {
        let bases = b"=ACMGRSVTWYHKDBNacgtn";
        let mut state = 0x1234_5678_9abc_def0u64;
        for case in 0..3000 {
            let n = (xorshift(&mut state) % 40) as usize;
            let seq: Vec<u8> = (0..n)
                .map(|_| bases[(xorshift(&mut state) as usize) % bases.len()])
                .collect();
            let missing = xorshift(&mut state) % 4 == 0;
            let qual: Vec<u8> = if missing {
                vec![0xff; n]
            } else {
                (0..n).map(|_| (xorshift(&mut state) % 60 + 1) as u8).collect()
            };
            let cigar1 = random_cigar(&mut state);

            let mut rec = bam::Record::new();
            let c1 = bam::record::CigarString(cigar1.clone());
            rec.set(b"read_name", Some(&c1), &seq, &qual);
            rec.push_aux(b"NM", Aux::I32(7)).unwrap();
            rec.push_aux(b"MD", Aux::String("5A3")).unwrap();
            rec.set_tid(3);
            rec.set_pos(1234);

            // Same op count, different lengths: the in-place path (the CIGAR,
            // SEQ and QUAL bytes change; qname and aux must be preserved).
            let nops = cigar1.len();
            let mut cigar2 = random_cigar(&mut state);
            cigar2.truncate(1);
            while cigar2.len() < nops {
                cigar2.push(match xorshift(&mut state) % 3 {
                    0 => Cigar::Match((xorshift(&mut state) % 20) as u32 + 1),
                    1 => Cigar::Ins((xorshift(&mut state) % 20) as u32 + 1),
                    _ => Cigar::Del((xorshift(&mut state) % 20) as u32 + 1),
                });
            }
            cigar2.truncate(nops);
            let seq2: Vec<u8> = (0..n)
                .map(|_| bases[(xorshift(&mut state) as usize) % bases.len()])
                .collect();

            // Reference: the old `set`-based rebuild (qual reversed explicitly).
            let mut want = rec.clone();
            let mut qual_rev = rec.qual().to_vec();
            qual_rev.reverse();
            let c2 = bam::record::CigarString(cigar2.clone());
            want.set(rec.qname(), Some(&c2), &seq2, &qual_rev);

            let mut got = rec.clone();
            assert!(
                rewrite_in_place(&mut got, &cigar2, &seq2),
                "case {case} n={n} nops={nops}"
            );
            assert_eq!(data_bytes(&want), data_bytes(&got), "case {case} n={n}");
            assert_eq!(core_debug(&want), core_debug(&got), "case {case} n={n}");
            // Public accessors agree too (cigar cache rebuild).
            assert_eq!(
                want.cigar().iter().collect::<Vec<_>>(),
                got.cigar().iter().collect::<Vec<_>>()
            );
            assert_eq!(want.seq().as_bytes(), got.seq().as_bytes());
            assert_eq!(want.qual(), got.qual());
            assert_eq!(want.qname(), got.qname());

            // A different op count must fall back (return false) untouched.
            let mut longer = cigar2.clone();
            longer.push(Cigar::Match(1));
            let mut untouched = rec.clone();
            let before = data_bytes(&untouched);
            assert!(!rewrite_in_place(&mut untouched, &longer, &seq2));
            assert_eq!(data_bytes(&untouched), before);
        }
    }

    #[test]
    fn pack_seq_matches_record_set_encoding() {
        let bases = b"=ACMGRSVTWYHKDBNacgtn";
        let mut state = 0xfeed_beef_dead_cafeu64;
        for n in 0..64usize {
            let seq: Vec<u8> = (0..n)
                .map(|_| bases[(xorshift(&mut state) as usize) % bases.len()])
                .collect();
            let qual = vec![30u8; n];
            let mut rec = bam::Record::new();
            let cigar = if n == 0 { None } else { Some(bam::record::CigarString(vec![Cigar::Match(n as u32)])) };
            rec.set(b"r", cigar.as_ref(), &seq, &qual);
            let mut dst = vec![0u8; n.div_ceil(2)];
            pack_seq(&mut dst, &seq);
            assert_eq!(dst.as_slice(), rec.seq().encoded, "n={n}");
        }
        // Every byte value, to pin the full ENCODE table against rust-htslib.
        let all: Vec<u8> = (0u16..256).map(|b| b as u8).collect();
        let mut rec = bam::Record::new();
        let cigar = bam::record::CigarString(vec![Cigar::Match(all.len() as u32)]);
        rec.set(b"r", Some(&cigar), &all, &vec![30u8; all.len()]);
        let mut dst = vec![0u8; all.len().div_ceil(2)];
        pack_seq(&mut dst, &all);
        assert_eq!(dst.as_slice(), rec.seq().encoded);
    }

    /// Count aux entries with the given tag.
    fn count_tag(rec: &bam::Record, tag: &[u8]) -> usize {
        let mut n = 0;
        for item in rec.aux_iter() {
            if let Ok((t, _)) = item {
                if t == tag {
                    n += 1;
                }
            }
        }
        n
    }

    /// M10b: `strip_translated_tags` removes every `MD`/`s2` entry (including
    /// duplicates) in one pass, keeps every other tag and the record bytes, and is
    /// a no-op when neither is present. Both tag orders and duplicates are
    /// exercised.
    #[test]
    fn strip_translated_tags_removes_md_and_s2() {
        fn build(order: &[&[u8]], dup: bool) -> bam::Record {
            let mut rec = bam::Record::new();
            let cigar = bam::record::CigarString(vec![Cigar::Match(4)]);
            rec.set(b"r", Some(&cigar), b"ACGT", &[30, 30, 30, 30]);
            rec.push_aux(b"NM", Aux::I32(3)).unwrap();
            for t in order {
                if *t == b"MD" {
                    rec.push_aux(b"MD", Aux::String("4")).unwrap();
                    if dup {
                        rec.push_aux_unchecked(b"MD", Aux::String("0")).unwrap();
                    }
                } else {
                    rec.push_aux(b"s2", Aux::I32(0)).unwrap();
                    if dup {
                        rec.push_aux_unchecked(b"s2", Aux::I32(7)).unwrap();
                    }
                }
            }
            rec.push_aux(b"AS", Aux::I32(100)).unwrap();
            rec.push_aux(b"SA", Aux::String("chr,1,+,4M,60,0;")).unwrap();
            rec
        }

        let md: &[u8] = b"MD";
        let s2: &[u8] = b"s2";
        for dup in [false, true] {
            for order in [&[md, s2][..], &[s2, md][..]] {
                let mut rec = build(order, dup);
                assert!(count_tag(&rec, b"MD") > 0 && count_tag(&rec, b"s2") > 0);
                strip_translated_tags(&mut rec);
                assert_eq!(count_tag(&rec, b"MD"), 0, "MD left (order {order:?}, dup {dup})");
                assert_eq!(count_tag(&rec, b"s2"), 0, "s2 left (order {order:?}, dup {dup})");
                // Everything else is preserved exactly once.
                assert_eq!(count_tag(&rec, b"NM"), 1);
                assert_eq!(count_tag(&rec, b"AS"), 1);
                assert_eq!(count_tag(&rec, b"SA"), 1);
                assert!(matches!(rec.aux(b"NM"), Ok(Aux::I32(3))));
                assert!(matches!(rec.aux(b"AS"), Ok(Aux::I32(100))));
                assert_eq!(rec.seq().as_bytes(), b"ACGT");
                assert_eq!(rec.qname(), b"r");
            }
        }

        // No-op when neither tag is present: the data bytes are untouched.
        let mut rec = build(&[], false);
        let before = data_bytes(&rec);
        strip_translated_tags(&mut rec);
        assert_eq!(data_bytes(&rec), before);
    }

    /// A record with any duplicated aux tag is invalid SAM and is flagged as
    /// unexpected by `scan_aux` (-> deferred `tags`).
    #[test]
    fn scan_aux_flags_duplicate_tags() {
        fn base() -> bam::Record {
            let mut rec = bam::Record::new();
            let cigar = bam::record::CigarString(vec![Cigar::Match(4)]);
            rec.set(b"r", Some(&cigar), b"ACGT", &[30, 30, 30, 30]);
            rec.push_aux(b"NM", Aux::I32(0)).unwrap();
            rec
        }
        // Unique tags are fine.
        let mut ok = base();
        ok.push_aux(b"s2", Aux::I32(0)).unwrap();
        ok.push_aux(b"MD", Aux::String("4")).unwrap();
        assert_eq!(scan_aux(&ok), (false, Some("4")));

        // A duplicate of an allowed tag, of MD, or of an unknown tag is flagged.
        let dup_tags: [&[u8]; 4] = [b"s2", b"MD", b"NM", b"zz"];
        for tag in dup_tags {
            let mut rec = base();
            if tag == b"MD" {
                rec.push_aux_unchecked(b"MD", Aux::String("4")).unwrap();
                rec.push_aux_unchecked(b"MD", Aux::String("0")).unwrap();
            } else {
                rec.push_aux_unchecked(tag, Aux::I32(0)).unwrap();
                rec.push_aux_unchecked(tag, Aux::I32(1)).unwrap();
            }
            assert!(scan_aux(&rec).0, "duplicate {tag:?} not flagged");
        }
    }

    /// `Params::new` records the effective `xl_tag` (default on) for stats.
    #[test]
    fn params_records_effective_xl_tag() {
        let cfg = TranslateConfig::default();
        assert!(cfg.xl_tag, "xl_tag must default to true");
        assert!(Params::new(&cfg, 0.01).xl_tag);
        let off = TranslateConfig { xl_tag: false, ..Default::default() };
        assert!(!Params::new(&off, 0.01).xl_tag);
    }
}

// ---------------------------------------------------------------------------
// statistics

/// Counts of records read from the input CRAM.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct RecordsIn {
    pub primary: u64,
    pub supplementary: u64,
    pub secondary: u64,
    pub unmapped: u64,
}

/// Step timings in seconds.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct Timings {
    pub diff: f64,
    /// Wall time of the translate/lift streaming pass.
    pub stream: f64,
    /// M4: minimap2 spawn to index ready, parsed from the captured `loaded/built
    /// the index` stderr line; 0 when that line is absent or unparseable.
    pub minimap2_index: f64,
    /// M4: wall time of the translate pass (alias of `stream`, kept explicit).
    pub translate_pass: f64,
    /// M4: wall time from pass end until the minimap2 reader thread drained.
    pub writer_drain: f64,
    /// M4: wall time from writer close to sorted/indexed CRAM ready.
    pub sort_index: f64,
    /// M4: total wall time of the minimap2 child.
    pub minimap2_total: f64,
    /// M4: CPU time of all reaped children (user+sys), via getrusage.
    pub children_cpu: f64,
    /// M5: wall time to build the clip k-mer screen.
    pub clip_screen: f64,
    pub total: f64,
}

/// M9: why the tiny-input all-remap fallback was taken. Serialized as part of
/// `Stats`; `null` when the fallback did not apply.
#[derive(Clone, Debug, Serialize)]
pub struct Fallback {
    /// `"tiny_input"` (the only reason today).
    pub reason: String,
    /// Size of the input CRAM in bytes (from `fs::metadata`, symlinks followed).
    pub input_bytes: u64,
    /// Effective threshold in bytes (`min_input_mb * 10^6`).
    pub threshold_bytes: u64,
}

/// `P.stats.json` contents.
#[derive(Clone, Debug, Default, Serialize)]
pub struct Stats {
    pub contract_ok: bool,
    pub contract_reason: Option<String>,
    pub preset: String,
    pub records_in: RecordsIn,
    pub translated_reads: u64,
    pub translated_records: u64,
    pub deferred_reads: u64,
    /// M4: primary records produced by the remap minimap2.
    pub remapped_reads: u64,
    /// M4: primary records written (translated + remapped).
    pub out_primary_records: u64,
    pub deferred_by_reason: BTreeMap<String, u64>,
    pub deferred_by_class: BTreeMap<String, u64>,
    pub origin_crossing_translated: u64,
    pub origin_crossing_deferred: u64,
    pub qual_missing: u64,
    /// M9: tiny-input fallback decision, or `null` when not applied. A
    /// tiny-input fallback leaves `contract_ok` true and `contract_reason` null.
    pub fallback: Option<Fallback>,
    pub timings_s: Timings,
    /// Effective tuning parameters this run used.
    pub params: Params,
}

impl Stats {
    fn bump_reason(&mut self, reason: Reason) {
        *self.deferred_by_reason.entry(reason.name().to_string()).or_insert(0) += 1;
        *self
            .deferred_by_class
            .entry(reason.class().name().to_string())
            .or_insert(0) += 1;
    }
}

// ---------------------------------------------------------------------------
// aux helpers

fn aux_f64(rec: &bam::Record, tag: &[u8]) -> Option<f64> {
    match rec.aux(tag).ok()? {
        Aux::Float(v) => Some(v as f64),
        Aux::Double(v) => Some(v),
        Aux::I8(v) => Some(v as f64),
        Aux::U8(v) => Some(v as f64),
        Aux::I16(v) => Some(v as f64),
        Aux::U16(v) => Some(v as f64),
        Aux::I32(v) => Some(v as f64),
        Aux::U32(v) => Some(v as f64),
        _ => None,
    }
}

/// Strict integer-only aux parse for the competitor tags: a float, double,
/// string, char or array value is rejected (returns `None`), so the record is
/// deferred rather than trusting a mis-typed `s1`/`s2`.
fn aux_int_strict(rec: &bam::Record, tag: &[u8]) -> Option<i64> {
    match rec.aux(tag).ok()? {
        Aux::I8(v) => Some(v as i64),
        Aux::U8(v) => Some(v as i64),
        Aux::I16(v) => Some(v as i64),
        Aux::U16(v) => Some(v as i64),
        Aux::I32(v) => Some(v as i64),
        Aux::U32(v) => Some(v as i64),
        _ => None,
    }
}

/// Single pass over the record's aux tags: reports whether it carries a tag that
/// is neither allowed nor the decode-generated `MD` (read for rule 7d and
/// stripped from the output), **or any duplicated tag** (invalid SAM; minimap2
/// never writes one), and returns that `MD` string if present. Scanning once
/// avoids a second `MD` lookup per record.
fn scan_aux(rec: &bam::Record) -> (bool, Option<&str>) {
    let mut unexpected = false;
    let mut md = None;
    // Distinct 2-byte tag keys seen so far. At most `ALLOWED_TAGS.len() + 1`
    // (MD) distinct tags are valid, so this small stack array suffices; a full
    // buffer implies a duplicate or an unknown tag, both `unexpected`.
    let mut seen = [[0u8; 2]; ALLOWED_TAGS.len() + 1];
    let mut n_seen = 0usize;
    for item in rec.aux_iter() {
        match item {
            Ok((tag, aux)) => {
                let key = [tag[0], tag[1]];
                if seen[..n_seen].contains(&key) {
                    unexpected = true; // duplicated tag
                } else if n_seen < seen.len() {
                    seen[n_seen] = key;
                    n_seen += 1;
                } else {
                    unexpected = true; // more distinct tags than are valid
                }
                if tag == b"MD".as_slice() {
                    if let Aux::String(s) = aux {
                        md = Some(s);
                    }
                } else if !ALLOWED_TAGS.iter().any(|t| *t == tag) {
                    unexpected = true;
                }
            }
            Err(_) => unexpected = true,
        }
    }
    (unexpected, md)
}

// ---------------------------------------------------------------------------
// streaming

/// Outcome of one translate pass over the CRAM.
pub(crate) enum PassOutcome {
    Completed,
    ContractViolation(String),
}

/// Everything the classifier needs, borrowed for the duration of one pass.
pub(crate) struct Ctx<'a> {
    pub(crate) old: &'a Fasta,
    pub(crate) new: &'a Fasta,
    pub(crate) diff: &'a Diff,
    pub(crate) old_names: &'a [String],
    pub(crate) new_tids: &'a HashMap<String, i32>,
    pub(crate) cfg: &'a TranslateConfig,
    pub(crate) thr: f64,
}

impl<'a> Ctx<'a> {
    fn new_len_of(&self, contig: &str) -> Option<u64> {
        self.new.len_of(contig)
    }
}

impl<'a> Ctx<'a> {
    fn to_info(&self, rec: &bam::Record) -> Result<RecInfo> {
        let tid = rec.tid();
        let contig = if tid >= 0 {
            self.old_names.get(tid as usize).cloned().unwrap_or_default()
        } else {
            String::new()
        };
        let cigar: Vec<Cigar> = rec.cigar().iter().copied().collect();
        // One aux scan: the unexpected-tag check and the `MD` string for rule 7d
        // (computed once here; O(mismatches + indels + MD runs), no allocation).
        let (unexpected_tag, md) = scan_aux(rec);
        let zdrop = match md {
            Some(s) => zdrop_statistic(&cigar, Some(s), &self.cfg.scoring),
            None => None,
        };
        // Track the SA tag's presence separately: an empty, malformed or
        // wrongly typed tag is a split-read deferral, not "no SA".
        let (sa_present, sa, sa_ok) = match rec.aux(b"SA") {
            Err(_) => (false, Vec::new(), false),
            Ok(Aux::String(s)) => {
                if s.is_empty() {
                    (true, Vec::new(), false)
                } else {
                    match parse_sa(s) {
                        Ok(v) if !v.is_empty() => (true, v, true),
                        _ => (true, Vec::new(), false),
                    }
                }
            }
            Ok(_) => (true, Vec::new(), false),
        };
        Ok(RecInfo {
            qname: rec.qname().to_vec(),
            flag: rec.flags(),
            contig,
            pos: rec.pos(),
            mapq: rec.mapq(),
            cigar,
            seq_len: rec.seq_len(),
            de: aux_f64(rec, b"de"),
            nm: aux_int(rec, b"NM"),
            s1: aux_int_strict(rec, b"s1"),
            s2: aux_int_strict(rec, b"s2"),
            sa_present,
            sa_ok,
            sa,
            unexpected_tag,
            xl_translated: matches!(rec.aux(b"XL"), Ok(Aux::Char(b'T'))),
            zdrop,
        })
    }
}

/// Open a CRAM reader with an explicit reference. The reference is mandatory
/// (`REF_PATH=""` in the environment makes network lookup impossible).
pub(crate) fn open_reader(path: &Path, ref_path: &Path, threads: usize) -> Result<bam::Reader> {
    let mut reader = bam::Reader::from_path(path)
        .with_context(|| format!("opening {}", path.display()))?;
    reader
        .set_reference(ref_path)
        .with_context(|| format!("setting reference {} for {}", ref_path.display(), path.display()))?;
    if threads > 1 {
        // Best effort; CRAM decode threading is optional.
        let _ = reader.set_threads(threads);
    }
    Ok(reader)
}

pub(crate) fn fai_path(fa: &Path) -> PathBuf {
    let mut s = fa.as_os_str().to_owned();
    s.push(".fai");
    PathBuf::from(s)
}

pub(crate) fn suffix(prefix: &Path, s: &str) -> PathBuf {
    let mut p = prefix.as_os_str().to_owned();
    p.push(s);
    PathBuf::from(p)
}

pub(crate) fn ensure_parent(prefix: &Path) -> Result<()> {
    if let Some(parent) = prefix.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("creating {}", parent.display()))?;
        }
    }
    Ok(())
}

/// Quote one argv element for the recorded `@PG` `CL`. Left unquoted when it
/// needs no quoting (so ordinary command lines are unchanged), else single-quoted
/// with each embedded `'` written as the `'\''` idiom. `contract::shell_split`
/// accepts exactly these forms.
fn shell_quote(s: &str) -> String {
    let safe = |c: char| {
        c.is_ascii_alphanumeric() || matches!(c, '%' | '+' | ',' | '-' | '.' | '/' | ':' | '=' | '@' | '_')
    };
    if !s.is_empty() && s.chars().all(safe) {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len() + 2);
    out.push('\'');
    for c in s.chars() {
        if c == '\'' {
            out.push_str("'\\''");
        } else {
            out.push(c);
        }
    }
    out.push('\'');
    out
}

/// Join argv into the recorded `@PG` `CL`, one shell word per element.
pub fn shell_join(argv: &[String]) -> String {
    argv.iter().map(|a| shell_quote(a)).collect::<Vec<_>>().join(" ")
}

/// The current process argv, shell-quoted for the output `@PG` `CL`.
pub fn recorded_argv() -> String {
    shell_join(&std::env::args().collect::<Vec<_>>())
}

/// argv-encoding marker written into every `cramhook lift` `@PG`'s `DS`. The
/// contract requires it: a lift `@PG` without it is from an older build whose
/// `CL` was recorded with `.join(" ")` (apostrophes in paths could then be read
/// as quoting and hide an option), so a chained input is refused and fully
/// remapped. Bump the version if the encoding ever changes.
pub const ARGV_ENCODING_MARKER: &str = "argv=sh1";

/// Build the output header: `@HD`, `@SQ` from the new `.fai`, all original
/// `@PG` lines and a trailing `cramhook` `@PG` chained with `PP`. `ds` is
/// appended as a `DS:` field (old/new reference and the remap command). The new
/// `@PG` gets the first unused `cramhook`/`cramhook.N` id (like samtools), so a
/// chain of lifts keeps unique ids and an unambiguous `PP` chain.
pub fn build_output_header(
    header_text: &str,
    new_fai: &[(String, u64)],
    argv: &str,
    ds: Option<&str>,
) -> Vec<u8> {
    let mut s = String::new();
    s.push_str("@HD\tVN:1.6\tSO:unsorted\n");
    for (name, len) in new_fai {
        s.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
    }
    for line in header_text.split('\n') {
        if line.starts_with("@PG\t") {
            s.push_str(line);
            s.push('\n');
        }
    }
    let mut last_pg_id: Option<String> = None;
    let mut used_ids: HashSet<String> = HashSet::new();
    for l in parse_header(header_text) {
        if l.rec_type != "PG" {
            continue;
        }
        if let Some(id) = l.field("ID") {
            last_pg_id = Some(id.to_string());
            used_ids.insert(id.to_string());
        }
    }
    // First unused id: `cramhook`, then `cramhook.1`, `cramhook.2`, ...
    let mut id = "cramhook".to_string();
    let mut n = 1u32;
    while used_ids.contains(&id) {
        id = format!("cramhook.{n}");
        n += 1;
    }
    s.push_str(&format!(
        "@PG\tID:{id}\tPN:cramhook\tVN:{}\tCL:{argv}",
        env!("CARGO_PKG_VERSION")
    ));
    if let Some(id) = last_pg_id {
        s.push_str(&format!("\tPP:{id}"));
    }
    // Always write `DS`, with the argv-encoding marker last so the contract can
    // tell a chained `@PG` written by this build from a legacy unquoted one.
    match ds {
        Some(ds) => s.push_str(&format!("\tDS:{ds} {ARGV_ENCODING_MARKER}")),
        None => s.push_str(&format!("\tDS:{ARGV_ENCODING_MARKER}")),
    }
    s.push('\n');
    s.into_bytes()
}

/// Build the FASTQ SEQ and QUAL (Phred+33) for a deferred primary record,
/// reverse-complementing/reversing when `0x10` is set.
fn fastq_bytes(rec: &bam::Record, stats: &mut Stats) -> (Vec<u8>, Vec<u8>) {
    let mut seq = rec.seq().as_bytes();
    let mut qual = rec.qual().to_vec();
    if rec.flags() & 0x10 != 0 {
        seq = revcomp(&seq);
        qual.reverse();
    }
    let missing = !qual.is_empty() && qual.iter().all(|&b| b == 0xff);
    let qual_ascii: Vec<u8> = if qual.is_empty() || missing {
        if missing {
            stats.qual_missing += 1;
        }
        vec![b'!'; seq.len()]
    } else {
        qual.iter().map(|&b| b.saturating_add(33)).collect()
    };
    (seq, qual_ascii)
}

/// Output sink for a translate/lift pass. `translate` writes files; `lift`
/// streams translated records into `samtools sort` and deferred reads into
/// `minimap2`'s stdin.
pub trait LiftSink {
    fn translated(&mut self, rec: &bam::Record) -> Result<()>;
    fn fastq(&mut self, qname: &[u8], seq: &[u8], qual: &[u8]) -> Result<()>;
    fn deferred(&mut self, qname: &[u8], reason: Reason) -> Result<()>;
    /// Flush/close whatever the pass owns. Called once when a pass completes.
    fn finish(&mut self) -> Result<()> {
        Ok(())
    }
}

/// File-based sink used by the `translate` debug subcommand.
pub struct FileSink {
    bam: bam::Writer,
    fq: BufWriter<File>,
    tsv: BufWriter<File>,
}

impl FileSink {
    pub fn new(prefix: &Path, header_bytes: &[u8]) -> Result<Self> {
        let header_view = bam::HeaderView::from_bytes(header_bytes);
        let header = bam::Header::from_template(&header_view);
        let bam = bam::Writer::from_path(suffix(prefix, ".translated.bam"), &header, Format::Bam)
            .context("creating translated.bam")?;
        let fq = BufWriter::new(
            File::create(suffix(prefix, ".deferred.fq")).context("creating deferred.fq")?,
        );
        let tsv = BufWriter::new(
            File::create(suffix(prefix, ".deferred.tsv")).context("creating deferred.tsv")?,
        );
        Ok(Self { bam, fq, tsv })
    }
}

impl LiftSink for FileSink {
    fn translated(&mut self, rec: &bam::Record) -> Result<()> {
        self.bam.write(rec)?;
        Ok(())
    }
    fn fastq(&mut self, qname: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
        self.fq.write_all(b"@")?;
        self.fq.write_all(qname)?;
        self.fq.write_all(b"\n")?;
        self.fq.write_all(seq)?;
        self.fq.write_all(b"\n+\n")?;
        self.fq.write_all(qual)?;
        self.fq.write_all(b"\n")?;
        Ok(())
    }
    fn deferred(&mut self, qname: &[u8], reason: Reason) -> Result<()> {
        let mut row = Vec::new();
        row.extend_from_slice(qname);
        row.push(b'\t');
        row.extend_from_slice(reason.name().as_bytes());
        row.push(b'\t');
        row.extend_from_slice(reason.class().name().as_bytes());
        row.push(b'\n');
        self.tsv.write_all(&row)?;
        Ok(())
    }
    fn finish(&mut self) -> Result<()> {
        self.fq.flush().context("flushing deferred.fq")?;
        self.tsv.flush().context("flushing deferred.tsv")?;
        Ok(())
    }
}

/// The 4-bit base encoding rust-htslib's `Record::set` uses (`ENCODE_BASE`).
/// For the `=ACMGRSVTWYHKDBN` code characters it is the inverse of
/// `DECODE_BASE`, case-insensitive; the four `'0'`..`'3'` entries replicate a
/// rust-htslib quirk so this is byte-identical for **any** input.
const fn build_encode() -> [u8; 256] {
    let dec = b"=ACMGRSVTWYHKDBN";
    let mut t = [15u8; 256];
    let mut i = 0usize;
    while i < 16 {
        let c = dec[i];
        t[c as usize] = i as u8;
        if c >= b'A' && c <= b'Z' {
            t[(c + 32) as usize] = i as u8;
        }
        i += 1;
    }
    t[b'0' as usize] = 1;
    t[b'1' as usize] = 2;
    t[b'2' as usize] = 4;
    t[b'3' as usize] = 8;
    t
}
static ENCODE: [u8; 256] = build_encode();

/// BAM CIGAR op code for a SAM op character (`Cigar::encode` is private in
/// rust-htslib, so this mirrors it).
#[inline]
fn cigar_op_code(c: char) -> u32 {
    match c {
        'M' => 0,
        'I' => 1,
        'D' => 2,
        'N' => 3,
        'S' => 4,
        'H' => 5,
        'P' => 6,
        '=' => 7,
        'X' => 8,
        _ => 0,
    }
}

/// Pack an ASCII sequence into BAM's 4-bit encoding into `dst` (whose length
/// must be `seq.len().div_ceil(2)`), matching `Record::set` byte for byte.
#[inline]
fn pack_seq(dst: &mut [u8], seq: &[u8]) {
    for (k, chunk) in seq.chunks(2).enumerate() {
        let hi = ENCODE[chunk[0] as usize];
        let lo = if chunk.len() > 1 { ENCODE[chunk[1] as usize] } else { 0 };
        dst[k] = (hi << 4) | lo;
    }
}

/// Rewrite CIGAR, SEQ and QUAL of `rec` in place, reusing the existing data
/// buffer: the CIGAR op count and the read length must be unchanged, so no
/// reallocation, `aux` memmove or qname copy is needed. Returns `false` when
/// the op count differs (caller falls back to `Record::set`).
fn rewrite_in_place(rec: &mut bam::Record, cigar: &[Cigar], seq: &[u8]) -> bool {
    let l_qseq = rec.seq_len();
    if cigar.len() != rec.cigar_len() || seq.len() != l_qseq {
        return false;
    }
    let qname_len = rec.inner().core.l_qname as usize;
    let n_cigar = cigar.len();
    let seq_bytes = l_qseq.div_ceil(2);
    if (rec.inner().l_data as usize) < qname_len + 4 * n_cigar + seq_bytes + l_qseq {
        return false;
    }
    // SAFETY: the offsets are the canonical BAM layout (qname is 4-byte
    // aligned, so the CIGAR array is too) and the bounds check above keeps
    // every write inside `l_data`. `rec`'s cached CIGAR is `None` (only
    // `cache_cigar` sets it and it is never called), so the bytes we write are
    // what `cigar()`/`raw_cigar()` see.
    unsafe {
        let inner = rec.inner_mut();
        let base = inner.data;
        let cg = std::slice::from_raw_parts_mut(base.add(qname_len) as *mut u32, n_cigar);
        for (i, c) in cigar.iter().enumerate() {
            cg[i] = (c.len() << 4) | cigar_op_code(c.char());
        }
        let seq_off = qname_len + 4 * n_cigar;
        let sb = std::slice::from_raw_parts_mut(base.add(seq_off), seq_bytes);
        pack_seq(sb, seq);
        let qb = std::slice::from_raw_parts_mut(base.add(seq_off + seq_bytes), l_qseq);
        qb.reverse();
    }
    true
}

/// Translate one record (Plan 4.3) in place, checking the sequence invariant.
///
/// Callers own `rec` (or a buffered group) and strip the decode-generated `MD`
/// and the old-assembly `s2` first (see [`strip_translated_tags`]), so no clone is
/// needed: `+` records only move tid/pos, `-` records
/// rewrite CIGAR/SEQ/QUAL in the existing data buffer when the op count is
/// unchanged (the common case), falling back to `Record::set` otherwise.
fn translate_record_in_place(
    rec: &mut bam::Record,
    info: &RecInfo,
    block: &Block,
    new_tid: i32,
    old: &Fasta,
    new: &Fasta,
    new_rc: Option<&[u8]>,
) -> Result<()> {
    let ref_start = info.pos as u64;
    let re = info.ref_end() as u64;
    let old_iv = OldIv::new(ref_start, re);
    let new_iv = block.map_iv(old_iv).ok_or_else(|| {
        anyhow::anyhow!(
            "block does not contain the alignment for read {:?}",
            String::from_utf8_lossy(&info.qname)
        )
    })?;

    // Invariant: reference bases must match (reverse-complemented for '-').
    let oseq = old
        .seq(&info.contig)
        .with_context(|| format!("old contig {:?} not in old FASTA", info.contig))?;
    let nseq = new
        .seq(&block.new_contig)
        .with_context(|| format!("new contig {:?} not in new FASTA", block.new_contig))?;
    if re > oseq.len() as u64 || new_iv.end.get() > nseq.len() as u64 {
        bail!(
            "invariant check out of range for read {:?}",
            String::from_utf8_lossy(&info.qname)
        );
    }
    let o = &oseq[ref_start as usize..re as usize];
    let n = &nseq[new_iv.start.get() as usize..new_iv.end.get() as usize];
    // Reference bases must match (reverse-complemented for '-'). For '-' this
    // compares against the precomputed reverse-complement of the whole contig
    // (`new_rc`), which turns the per-record check into a single memcmp instead
    // of materialising `revcomp(n)` on every record. Falls back to `revcomp`
    // when the cache is absent (e.g. `run_all_defer`).
    let ok = match block.strand {
        Strand::Plus => o == n,
        Strand::Minus => match new_rc {
            Some(rc) => &rc[rc.len() - new_iv.end.get() as usize..rc.len() - new_iv.start.get() as usize] == o,
            None => revcomp(n).as_slice() == o,
        },
    };
    if !ok {
        bail!(
            "invariant violation for read {:?}: old {} {}..{} ({} bp) does not match new {} {}..{} (strand {})",
            String::from_utf8_lossy(&info.qname),
            info.contig,
            ref_start,
            re,
            old_iv.len(),
            block.new_contig,
            new_iv.start.get(),
            new_iv.end.get(),
            block.strand.as_char()
        );
    }

    // `rec` has the `MD` and `s2` tags already stripped by the caller (see
    // `emit_single` / `emit_split`); the `-` branch rewrites the data buffer in
    // place.
    if block.strand == Strand::Minus {
        let reversed = reverse_cigar(&info.cigar);
        // Reverse-complement straight from the encoded 4-bit sequence: one
        // table lookup per byte pair instead of per-base decoding.
        let seq = revcomp_seq(rec);
        // minimap2 left-aligns indels on the new reference; reversing the old
        // CIGAR puts each indel at the right end of its equivalent range.
        let cigar = left_align_cigar(&reversed, &seq, n);
        if !rewrite_in_place(rec, &cigar, &seq) {
            // CIGAR op count changed: rebuild through the public API.
            let qname = rec.qname().to_vec();
            let mut qual: Vec<u8> = rec.qual().to_vec();
            qual.reverse();
            let c = bam::record::CigarString(cigar);
            rec.set(&qname, Some(&c), &seq, &qual);
        }
        rec.set_flags(info.flag ^ 0x10);
    }
    rec.set_tid(new_tid);
    rec.set_pos(new_iv.start.get() as i64);
    Ok(())
}

fn add_xl(rec: &mut bam::Record) {
    // The inherited `XL` was removed by `strip_translated_tags` before this, so a
    // plain push cannot create a duplicate. `update_aux` does not support `Char`.
    let _ = rec.push_aux(b"XL", Aux::Char(b'T'));
}

/// Remove the aux tags a translated record must not carry, in a single pass over
/// the aux data:
///
/// * `MD` — the decode artifact (see `docs/decisions.md`); meaningless once the
///   CIGAR has been rewritten.
/// * `s2` — minimap2's best *competing* chain score on the **old** assembly. It
///   drifts with repeat/paralog content and cannot be recomputed by translation
///   (see `scratch/m10b/Report.md`), so it is stripped from translated records.
///   Remapped (`XL:A:R`) records keep minimap2's fresh `s2`.
/// * `XL` — the inherited provenance tag of a chained input. It is removed
///   unconditionally; [`add_xl`] adds the fresh `XL:A:T` afterwards only when
///   tagging is on, so `--no-xl-tag` leaves no inherited `XL:A:R` behind.
///
/// Rule 3b reads `s2` and `XL` from the input *before* this runs (in `to_info`).
/// Every matching entry is removed (duplicates are invalid SAM and normally
/// deferred, but stripping must not leave one behind): `bam_aux_remove` deletes
/// in place and returns a pointer to the next tag, so the single traversal
/// continues from there without rescanning.
fn strip_translated_tags(rec: &mut bam::Record) {
    use rust_htslib::htslib;
    // SAFETY: `b`/`p` point at `rec`'s aux data; the pointers come from htslib's
    // own traversal and are used only through htslib calls that mutate the record
    // in place (`bam_aux_remove` memmoves the tail of the same allocation and
    // returns the next entry). The record is not otherwise mutated here.
    unsafe {
        let b = rec.inner() as *const htslib::bam1_t;
        let p = rec.inner_mut() as *mut htslib::bam1_t;
        let mut s = htslib::bam_aux_first(b);
        while !s.is_null() {
            let tag = (*s.sub(2), *s.sub(1));
            if tag == (b'M', b'D') || tag == (b's', b'2') || tag == (b'X', b'L') {
                // Returns the pointer to the next tag (or NULL if it was last).
                s = htslib::bam_aux_remove(p, s);
            } else {
                s = htslib::bam_aux_next(b, s);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// clip k-mer screen

/// A trivial multiplicative hasher for the u64 k-mer keys. The set is
/// membership-only, so the hash function cannot affect output.
#[derive(Default)]
struct U64Hasher(u64);

impl std::hash::Hasher for U64Hasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 ^ b as u64).wrapping_mul(0x100000001B3);
        }
    }
    fn write_u64(&mut self, i: u64) {
        self.0 = i.wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(31);
    }
}

type KmerSet = HashSet<u64, std::hash::BuildHasherDefault<U64Hasher>>;

/// Fail-safe screen for soft clips (Variant B): a clip of at least `min_len` bp
/// is deferred when it shares at least `min_hits` canonical k-mers with the
/// changed regions of `new.fa`, where a k-mer counts only if it **contains at
/// least one novel base** or **spans an internal identity-block boundary** (new
/// edit site). Unchanged flank sequence around a polishing edit does *not* count,
/// so wrapping/plasmid clips that are merely more of the same plasmid no longer
/// over-defer. `flank > 0` adds extra sequence around novel intervals and edit
/// sites. Canonical k-mers make the test strand-invariant. Built once per pass.
///
/// This is a conservative *superset* of the minimap2 prototype: any clip that
/// minimap2 maps to a changed region must share k-mers with it, and the k-mer
/// test is more permissive than an ungapped alignment.
pub(crate) struct ClipScreen {
    k: usize,
    min_hits: u32,
    min_len: u64,
    kmers: KmerSet,
}

fn base_bits(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Canonical (strand-invariant) 2-bit encoding of a k-mer, or `None` if it
/// contains a non-ACGT base. Exact for `k <= 32` (no collisions).
fn canonical_kmer(w: &[u8], k: usize) -> Option<u64> {
    let shift = 2 * (k as u64 - 1);
    let mut fwd = 0u64;
    let mut rev = 0u64;
    for &b in w {
        let x = base_bits(b)?;
        fwd = (fwd << 2) | x;
        rev = (rev >> 2) | ((3 - x) << shift);
    }
    Some(fwd.min(rev))
}

/// Add canonical k-mers of `seq` whose window start lies in `[lo, hi)`, clipped
/// to the part of the contig where a full k-mer fits.
fn add_kmer_starts(kmers: &mut KmerSet, seq: &[u8], lo: u64, hi: u64, k: usize) {
    let n = seq.len();
    if n < k {
        return;
    }
    let lo = (lo as usize).min(n - k + 1);
    let hi = (hi as usize).min(n - k + 1);
    for i in lo..hi {
        if let Some(c) = canonical_kmer(&seq[i..i + k], k) {
            kmers.insert(c);
        }
    }
}

/// Add every canonical k-mer of `seq` (used to subtract old-sequence k-mers).
fn add_all_kmers(kmers: &mut KmerSet, seq: &[u8], k: usize) {
    if seq.len() < k {
        return;
    }
    for w in seq.windows(k) {
        if let Some(c) = canonical_kmer(w, k) {
            kmers.insert(c);
        }
    }
}

impl ClipScreen {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn build(
        new: &Fasta,
        old: &Fasta,
        diff: &Diff,
        k: u32,
        min_hits: u32,
        min_len: u64,
        flank: u64,
        subtract_old: bool,
    ) -> Result<Option<ClipScreen>> {
        if k == 0 || k > 32 {
            bail!("clip screen k must be in 1..=32 (got {k})");
        }
        let k = k as usize;
        let km1 = k as u64 - 1;
        let mut kmers: KmerSet = KmerSet::default();

        // Variant B: only k-mers that contain at least one novel base. A window
        // [i, i+k) contains a novel base iff it overlaps a novel interval, i.e.
        // i in [s-(k-1), e). `flank` widens the novel interval when > 0.
        for (contig, ivs) in &diff.novel {
            let Some(seq) = new.seq(contig) else { continue };
            let n = seq.len() as u64;
            for iv in ivs {
                let s = iv.start.0.saturating_sub(flank);
                let e = (iv.end.0 + flank).min(n);
                if e > s {
                    add_kmer_starts(&mut kmers, seq, s.saturating_sub(km1), e, k);
                }
            }
        }

        // Variant B: only k-mers that span an internal block boundary (new edit
        // site). Around boundary `b`, window starts in [b-(k-1), b+(k-1)), which
        // includes every window that crosses `b`. Contig-end boundaries are not
        // edit sites and are skipped.
        for blks in diff.blocks.values() {
            for b in blks {
                let Some(seq) = new.seq(&b.new_contig) else { continue };
                let n = seq.len() as u64;
                let end = b.new_start.0 + b.old.len();
                for boundary in [b.new_start.0, end] {
                    if boundary > 0 && boundary < n {
                        let lo = boundary.saturating_sub(km1 + flank);
                        let hi = boundary + km1 + flank;
                        add_kmer_starts(&mut kmers, seq, lo, hi, k);
                    }
                }
            }
        }

        // Variant B': drop every k-mer that also occurs anywhere in old.fa.
        if subtract_old {
            let mut old_kmers: KmerSet = KmerSet::default();
            for contig in old.order() {
                if let Some(seq) = old.seq(contig) {
                    add_all_kmers(&mut old_kmers, seq, k);
                }
            }
            kmers.retain(|c| !old_kmers.contains(c));
        }

        Ok(Some(ClipScreen { k, min_hits, min_len, kmers }))
    }

    /// Does this clip share enough k-mers with the changed-region set?
    ///
    /// Rolling 2-bit canonical encoding: O(clip length), not O(clip × k). Windows
    /// containing a non-ACGT base are skipped exactly as before (the running
    /// encoders are reset on `N`), so the result is identical to the previous
    /// per-window computation.
    pub(crate) fn hits(&self, clip: &[u8]) -> bool {
        let k = self.k;
        if (clip.len() as u64) < self.min_len || clip.len() < k {
            return false;
        }
        let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
        let shift = 2 * (k as u64 - 1);
        let mut fwd = 0u64;
        let mut rev = 0u64;
        let mut valid = 0usize;
        let mut n = 0u32;
        for &b in clip {
            let Some(x) = base_bits(b) else {
                fwd = 0;
                rev = 0;
                valid = 0;
                continue;
            };
            fwd = ((fwd << 2) | x) & mask;
            rev = (rev >> 2) | ((3 - x) << shift);
            valid += 1;
            if valid >= k {
                let c = fwd.min(rev);
                if self.kmers.contains(&c) {
                    n += 1;
                    if n >= self.min_hits {
                        return true;
                    }
                }
            }
        }
        false
    }
}

/// Screen every soft clip of one single record (both free ends). `l`/`t` are the
/// leading/trailing soft-clip lengths; the sequence is only materialised when a
/// clip is long enough to matter.
fn record_clip_hit(screen: &ClipScreen, rec: &bam::Record, l: usize, t: usize) -> bool {
    if l < screen.min_len as usize && t < screen.min_len as usize {
        return false;
    }
    let seq = rec.seq().as_bytes();
    (l > 0 && seq.len() >= l && screen.hits(&seq[..l]))
        || (t > 0 && seq.len() >= t && screen.hits(&seq[seq.len() - t..]))
}

/// Screen only the boundary-reaching clip of an origin-crossing record, exactly
/// as before (the clip selection is what produces the current output; it must
/// not change). The scan itself is optimized and `info` supplies the clip
/// lengths so the record sequence is only materialised when needed.
fn outer_clip_hit(screen: &ClipScreen, rec: &bam::Record, info: &RecInfo, contig_len: u64) -> bool {
    let l = leading_soft(&info.cigar) as usize;
    let t = trailing_soft(&info.cigar) as usize;
    let reaches_end = info.ref_end() + t as i64 >= contig_len as i64;
    let reaches_start = info.pos - l as i64 <= 0;
    let do_t = reaches_end && t > 0;
    let do_l = reaches_start && l > 0;
    if (!do_t || t < screen.min_len as usize) && (!do_l || l < screen.min_len as usize) {
        return false;
    }
    let seq = rec.seq().as_bytes();
    (do_t && seq.len() >= t && screen.hits(&seq[seq.len() - t..]))
        || (do_l && seq.len() >= l && screen.hits(&seq[..l]))
}

/// Bookkeeping and output sink for one pass.
struct Pass<'a, S: LiftSink> {
    ctx: Ctx<'a>,
    sink: S,
    stats: Stats,
    /// The clip screen, if enabled.
    screen: Option<ClipScreen>,
    // ponytail: read names of decided reads; ~tens of bytes each. A few million
    // reads is fine (< ~500 MB); a production run would use a hash of the name.
    seen: HashSet<Vec<u8>>,
    /// Reverse-complemented copy of every contig of `new.fa`, built once. Used to
    /// turn the `-` invariant check into a memcmp. Empty in `run_all_defer`.
    new_rc: HashMap<String, Vec<u8>>,
    /// Records of a true origin-crossing candidate, keyed by read name.
    pending: HashMap<Vec<u8>, Vec<bam::Record>>,
    /// Reads decided to be deferred on their first (supplementary) record; we
    /// wait for the primary to source the FASTQ entry, then drop the rest.
    pending_reason: HashMap<Vec<u8>, Reason>,
}

impl<'a, S: LiftSink> Pass<'a, S> {
    fn new(ctx: Ctx<'a>, sink: S, screen: Option<ClipScreen>) -> Self {
        Self {
            ctx,
            sink,
            stats: Stats::default(),
            screen,
            seen: HashSet::new(),
            new_rc: HashMap::new(),
            pending: HashMap::new(),
            pending_reason: HashMap::new(),
        }
    }

    fn mark_name(&mut self, qname: &[u8]) -> Result<()> {
        if !self.seen.insert(qname.to_vec()) {
            bail!("read {:?} appears more than once", String::from_utf8_lossy(qname));
        }
        Ok(())
    }

    fn defer(&mut self, primary: &bam::Record, reason: Reason) -> Result<()> {
        self.mark_name(primary.qname())?;
        let (seq, qual) = fastq_bytes(primary, &mut self.stats);
        self.sink.fastq(primary.qname(), &seq, &qual)?;
        self.sink.deferred(primary.qname(), reason)?;
        self.stats.deferred_reads += 1;
        self.stats.bump_reason(reason);
        Ok(())
    }

    /// Translate and write one eligible single record, editing it in place.
    fn emit_single(&mut self, rec: &mut bam::Record, info: &RecInfo, block: &Block) -> Result<()> {
        let new_tid = *self
            .ctx
            .new_tids
            .get(&block.new_contig)
            .with_context(|| format!("new contig {:?} not in new .fai", block.new_contig))?;
        translate_record_in_place(rec, info, block, new_tid, self.ctx.old, self.ctx.new, self.new_rc.get(&block.new_contig).map(Vec::as_slice))?;
        if self.ctx.cfg.xl_tag {
            add_xl(rec);
        }
        self.mark_name(&info.qname)?;
        self.sink.translated(rec)?;
        self.stats.translated_reads += 1;
        self.stats.translated_records += 1;
        Ok(())
    }

    /// Translate and write both records of an origin-crossing split read.
    fn emit_split(
        &mut self,
        group: &mut [bam::Record],
        infos: &[RecInfo],
        plan: &SplitPlan<'_>,
    ) -> Result<()> {
        let a_info = &infos[plan.a];
        let b_info = &infos[plan.b];
        let a_tid = *self
            .ctx
            .new_tids
            .get(&plan.a_block.new_contig)
            .with_context(|| format!("new contig {:?} not in new .fai", plan.a_block.new_contig))?;
        let b_tid = *self
            .ctx
            .new_tids
            .get(&plan.b_block.new_contig)
            .with_context(|| format!("new contig {:?} not in new .fai", plan.b_block.new_contig))?;
        translate_record_in_place(
            &mut group[plan.a],
            a_info,
            plan.a_block,
            a_tid,
            self.ctx.old,
            self.ctx.new,
            self.new_rc.get(&plan.a_block.new_contig).map(Vec::as_slice),
        )?;
        translate_record_in_place(
            &mut group[plan.b],
            b_info,
            plan.b_block,
            b_tid,
            self.ctx.old,
            self.ctx.new,
            self.new_rc.get(&plan.b_block.new_contig).map(Vec::as_slice),
        )?;

        let nm_a =
            aux_int(&group[plan.a], b"NM").ok_or_else(|| anyhow::anyhow!("split record without NM"))?;
        let nm_b =
            aux_int(&group[plan.b], b"NM").ok_or_else(|| anyhow::anyhow!("split record without NM"))?;
        let sa_for_a = format_sa(&sa_source(&group[plan.b], &plan.b_block.new_contig, nm_b)?);
        let sa_for_b = format_sa(&sa_source(&group[plan.a], &plan.a_block.new_contig, nm_a)?);
        group[plan.a].update_aux(b"SA", Aux::String(&sa_for_a))?;
        group[plan.b].update_aux(b"SA", Aux::String(&sa_for_b))?;
        if self.ctx.cfg.xl_tag {
            add_xl(&mut group[plan.a]);
            add_xl(&mut group[plan.b]);
        }
        self.mark_name(&a_info.qname)?;
        self.sink.translated(&group[plan.a])?;
        self.sink.translated(&group[plan.b])?;
        self.stats.translated_reads += 1;
        self.stats.translated_records += 2;
        self.stats.origin_crossing_translated += 1;
        Ok(())
    }

    /// Decide what to do with a completed group of records that were buffered
    /// because they carry SA (or are supplementary).
    fn process_group(&mut self, mut group: Vec<bam::Record>) -> Result<()> {
        let infos: Vec<RecInfo> = group.iter().map(|r| self.ctx.to_info(r)).collect::<Result<_>>()?;
        let contig = infos.first().map(|i| i.contig.clone()).unwrap_or_default();
        let contig_len = self.ctx.old.len_of(&contig).unwrap_or(0);
        let outcome = {
            let new_len_of = |c: &str| self.ctx.new_len_of(c);
            classify_origin_crossing(&infos, contig_len, self.ctx.diff, self.ctx.thr, self.ctx.cfg, &new_len_of)
        };
        match outcome {
            Ok(plan) => {
                let clip_hit = self.screen.as_ref().is_some_and(|s| {
                    group
                        .iter()
                        .zip(infos.iter())
                        .any(|(r, info)| outer_clip_hit(s, r, info, contig_len))
                });
                if clip_hit {
                    let primary = group
                        .iter()
                        .find(|r| r.flags() & (0x100 | 0x800) == 0)
                        .unwrap_or(&group[0])
                        .clone();
                    self.stats.origin_crossing_deferred += 1;
                    self.defer(&primary, Reason::SplitClipHit)
                } else {
                    // Strip `MD` and the old-assembly `s2` first: the records
                    // are edited in place, and carrying those tags costs a copy.
                    // Both statistics (`zdrop`, `s2` for rule 3b) are already in
                    // `infos`.
                    for r in group.iter_mut() {
                        strip_translated_tags(r);
                    }
                    self.emit_split(&mut group, &infos, &plan)
                }
            }
            Err(reason) => {
                let primary = group
                    .iter()
                    .find(|r| r.flags() & (0x100 | 0x800) == 0)
                    .unwrap_or(&group[0])
                    .clone();
                self.stats.origin_crossing_deferred += 1;
                self.defer(&primary, reason)
            }
        }
    }
}

fn sa_source(rec: &bam::Record, contig: &str, nm: i64) -> Result<SaSource> {
    let cigar: Vec<Cigar> = rec.cigar().iter().copied().collect();
    Ok(SaSource {
        contig: contig.to_string(),
        pos1: (rec.pos() + 1) as u64,
        reverse: rec.flags() & 0x10 != 0,
        cigar,
        l_seq: rec.seq_len() as u64,
        mapq: rec.mapq(),
        nm,
    })
}

/// Stream in all-defer mode: no translation, every primary read deferred with
/// `defer_reason`.
///
/// With `validate` the record-level contract checks that [`run_translate_pass`]
/// performs are kept: a secondary record or a hard-clipped primary aborts the
/// stream with `PassOutcome::ContractViolation` (used by the tiny-input
/// fallback so that a real contract violation still takes precedence). Contract
/// fallbacks pass `validate = false` because every read is already remapped.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_all_defer<S: LiftSink>(
    reader: &mut bam::Reader,
    sink: S,
    ctx: Ctx<'_>,
    contract_reason: Option<String>,
    defer_reason: Reason,
    validate: bool,
    t_stream: Instant,
    progress: &crate::progress::Progress,
) -> Result<PassOutcomeAndStats> {
    let mut frac = crate::progress::RefFrac::from_header(progress, reader.header());
    let mut pass = Pass::new(ctx, sink, None);
    pass.stats.contract_ok = contract_reason.is_none();
    pass.stats.contract_reason = contract_reason.clone();
    pass.stats.preset = "unknown".to_string();

    let mut rec = bam::Record::new();
    while let Some(res) = reader.read(&mut rec) {
        res?;
        if let Some(f) = frac.as_mut() {
            f.record(rec.tid(), rec.pos());
        }
        let _ = rec.remove_aux(b"MD");
        let flag = rec.flags();
        let secondary = flag & 0x100 != 0;
        let supplementary = flag & 0x800 != 0;
        let unmapped = flag & 0x4 != 0;
        pass.stats.records_in.secondary += secondary as u64;
        pass.stats.records_in.supplementary += supplementary as u64;
        pass.stats.records_in.unmapped += unmapped as u64;
        if !secondary && !supplementary {
            pass.stats.records_in.primary += 1;
        }

        if validate {
            if secondary {
                pass.stats.timings_s.stream = t_stream.elapsed().as_secs_f64();
                let stats = pass.stats.clone();
                return Ok(PassOutcomeAndStats {
                    outcome: PassOutcome::ContractViolation("secondary record (0x100) present".into()),
                    stats,
                });
            }
            let info = pass.ctx.to_info(&rec)?;
            if !supplementary && has_hard_clip(&info.cigar) {
                pass.stats.timings_s.stream = t_stream.elapsed().as_secs_f64();
                let stats = pass.stats.clone();
                return Ok(PassOutcomeAndStats {
                    outcome: PassOutcome::ContractViolation("hard-clipped primary record".into()),
                    stats,
                });
            }
        }

        if !secondary && !supplementary {
            pass.defer(&rec, defer_reason)?;
        }
    }
    pass.sink.finish()?;
    pass.stats.timings_s.stream = t_stream.elapsed().as_secs_f64();
    Ok(PassOutcomeAndStats { outcome: PassOutcome::Completed, stats: pass.stats.clone() })
}

/// Stream once in normal mode. Returns a mid-stream contract violation if one
/// is found (`--strict` bails before that in the caller).
pub(crate) fn run_translate_pass<S: LiftSink>(
    reader: &mut bam::Reader,
    sink: S,
    ctx: Ctx<'_>,
    preset: Preset,
    t_stream: Instant,
    progress: &crate::progress::Progress,
) -> Result<PassOutcomeAndStats> {
    let mut frac = crate::progress::RefFrac::from_header(progress, reader.header());
    let t_screen = Instant::now();
    let screen = if ctx.cfg.unsafe_no_clip_screen {
        None
    } else {
        ClipScreen::build(
            ctx.new,
            ctx.old,
            ctx.diff,
            ctx.cfg.clip_k,
            ctx.cfg.clip_min_hits,
            ctx.cfg.clip_min_len,
            ctx.cfg.clip_flank,
            ctx.cfg.clip_subtract_old,
        )?
    };
    let clip_screen_s = t_screen.elapsed().as_secs_f64();
    let mut pass = Pass::new(ctx, sink, screen);
    // Reverse-complement every new contig once, so the `-` invariant check is a
    // memcmp rather than a per-record `revcomp`.
    let new_fa = pass.ctx.new;
    pass.new_rc = new_fa
        .order()
        .iter()
        .map(|name| (name.clone(), revcomp(new_fa.seq(name).unwrap_or(&[]))))
        .collect();
    pass.stats.contract_ok = true;
    pass.stats.preset = preset.minimap2_name().to_string();
    pass.stats.timings_s.clip_screen = clip_screen_s;

    let mut rec = bam::Record::new();
    while let Some(res) = reader.read(&mut rec) {
        res?;
        if let Some(f) = frac.as_mut() {
            f.record(rec.tid(), rec.pos());
        }
        let flag = rec.flags();
        let secondary = flag & 0x100 != 0;
        let supplementary = flag & 0x800 != 0;
        let unmapped = flag & 0x4 != 0;
        pass.stats.records_in.secondary += secondary as u64;
        pass.stats.records_in.supplementary += supplementary as u64;
        pass.stats.records_in.unmapped += unmapped as u64;
        if !secondary && !supplementary {
            pass.stats.records_in.primary += 1;
        }

        if secondary {
            pass.stats.timings_s.stream = t_stream.elapsed().as_secs_f64();
            let stats = pass.stats.clone();
            return Ok(PassOutcomeAndStats {
                outcome: PassOutcome::ContractViolation("secondary record (0x100) present".into()),
                stats,
            });
        }
        let info = pass.ctx.to_info(&rec)?;
        if !supplementary && has_hard_clip(&info.cigar) {
            pass.stats.timings_s.stream = t_stream.elapsed().as_secs_f64();
            let stats = pass.stats.clone();
            return Ok(PassOutcomeAndStats {
                outcome: PassOutcome::ContractViolation("hard-clipped primary record".into()),
                stats,
            });
        }

        // A read's outcome is decided once; later records are dropped. Uses the
        // borrowed qname (no per-record clone; the name is only cloned when it
        // is stored in `pending`/`pending_reason`).
        if pass.seen.contains(info.qname.as_slice()) {
            continue;
        }
        // Non-origin split read decided on its first (supplementary) record:
        // only the name + reason are stored, never the records.
        if pass.pending_reason.contains_key(info.qname.as_slice()) {
            if !supplementary {
                let reason = pass.pending_reason.remove(info.qname.as_slice()).unwrap();
                pass.stats.origin_crossing_deferred += 1;
                pass.defer(&rec, reason)?;
            }
            continue;
        }
        // Mate of a buffered true origin-crossing candidate.
        if pass.pending.contains_key(info.qname.as_slice()) {
            if supplementary || info.sa_present {
                let v = pass.pending.get_mut(info.qname.as_slice()).unwrap();
                v.push(rec.clone());
                if v.len() == 2 {
                    let group = pass.pending.remove(info.qname.as_slice()).unwrap();
                    pass.process_group(group)?;
                }
            } else {
                // The candidate's mate is a primary without `SA`: the read is
                // inconsistent. Defer it from the arriving primary and drop the
                // buffered records.
                pass.pending.remove(info.qname.as_slice());
                pass.stats.origin_crossing_deferred += 1;
                pass.defer(&rec, Reason::SplitSaInconsistent)?;
            }
            continue;
        }

        let is_split = supplementary || info.sa_present;
        if !is_split {
            if unmapped {
                pass.defer(&rec, Reason::Unmapped)?;
            } else {
                let l = leading_soft(&info.cigar) as usize;
                let t = trailing_soft(&info.cigar) as usize;
                let clip = (l + t) as u64;
                let contig_len = pass.ctx.old.len_of(&info.contig).unwrap_or(0);
                let outcome = {
                    let new_len_of = |c: &str| pass.ctx.new_len_of(c);
                    check_record(
                        &info,
                        contig_len,
                        pass.ctx.diff,
                        pass.ctx.thr,
                        pass.ctx.cfg,
                        clip,
                        &new_len_of,
                    )
                };
                match outcome {
                    Ok(block) => {
                        if pass.screen.as_ref().is_some_and(|s| record_clip_hit(s, &rec, l, t)) {
                            pass.defer(&rec, Reason::ClipHit)?;
                        } else {
                            // Strip `MD` and the old-assembly `s2` before
                            // translating: the record is edited in place and no
                            // longer cloned (rule 3b already read `s2` in `info`).
                            strip_translated_tags(&mut rec);
                            pass.emit_single(&mut rec, &info, block)?;
                        }
                    }
                    Err(reason) => pass.defer(&rec, reason)?,
                }
            }
            continue;
        }

        // First record of a split read: decide from it, so non-origin reads are
        // never buffered. Only true origin-crossing candidates hold records.
        let contig_len = pass.ctx.old.len_of(&info.contig).unwrap_or(0);
        let decision = if !info.sa_ok {
            Err(Reason::Split)
        } else if info.sa.len() != 1 {
            Err(Reason::SplitNotTwoRecords)
        } else if info.sa[0].rname != info.contig {
            Err(Reason::SplitDiffContigs)
        } else if contig_len == 0
            || !(info.outer_clip_reaches_end(contig_len) || info.outer_clip_reaches_start())
        {
            Err(Reason::SplitNotOrigin)
        } else {
            Ok(())
        };
        match decision {
            Ok(()) => {
                let v = pass.pending.entry(info.qname.clone()).or_default();
                v.push(rec.clone());
            }
            Err(reason) if !supplementary => {
                pass.stats.origin_crossing_deferred += 1;
                pass.defer(&rec, reason)?;
            }
            Err(reason) => {
                // Supplementary-first non-origin: remember name + reason only.
                pass.pending_reason.insert(info.qname.clone(), reason);
            }
        }
    }

    // EOF: incomplete split reads are deferred with reason `split`. Only the
    // names are collected, never the buffered groups.
    let names: Vec<Vec<u8>> = pass.pending.keys().cloned().collect();
    for name in names {
        let Some(group) = pass.pending.remove(&name) else { continue };
        // Only a read with a primary record counts; a stray supplementary is
        // dropped (it was never counted in `records_in.primary`).
        let Some(primary) = group.iter().find(|r| r.flags() & (0x100 | 0x800) == 0) else {
            continue;
        };
        let primary = primary.clone();
        pass.stats.origin_crossing_deferred += 1;
        pass.defer(&primary, Reason::Split)?;
    }
    // Non-origin names whose primary never arrived: no primary to source from,
    // and they were never counted in `records_in.primary`.
    pass.pending_reason.clear();

    // Invariant: translated + deferred == primary records in; each name once.
    let total_reads = pass.stats.translated_reads + pass.stats.deferred_reads;
    if total_reads != pass.stats.records_in.primary {
        bail!(
            "read-count invariant violated: translated {} + deferred {} != primary in {}",
            pass.stats.translated_reads,
            pass.stats.deferred_reads,
            pass.stats.records_in.primary
        );
    }
    pass.sink.finish()?;
    pass.stats.timings_s.stream = t_stream.elapsed().as_secs_f64();
    let stats = pass.stats.clone();
    Ok(PassOutcomeAndStats { outcome: PassOutcome::Completed, stats })
}

pub(crate) struct PassOutcomeAndStats {
    pub(crate) outcome: PassOutcome,
    pub(crate) stats: Stats,
}

// ---------------------------------------------------------------------------
// top-level entry point

/// Run `cramhook translate`. Returns the stats that were written.
pub fn translate(
    in_cram: &Path,
    old_ref: &Path,
    new_ref: &Path,
    out_prefix: &Path,
    cfg: &TranslateConfig,
    threads: usize,
) -> Result<Stats> {
    let t0 = Instant::now();
    cfg.validate()?;
    warn_if_unsafe_clip_screen(cfg);
    ensure_parent(out_prefix)?;
    ensure_fai(old_ref)?;
    ensure_fai(new_ref)?;
    let old_fai = read_fai(&fai_path(old_ref))?;
    let new_fai = read_fai(&fai_path(new_ref))?;
    let new_tids: HashMap<String, i32> =
        new_fai.iter().enumerate().map(|(i, (n, _))| (n.clone(), i as i32)).collect();

    let mut reader = open_reader(in_cram, old_ref, threads)?;
    let header_text = String::from_utf8_lossy(reader.header().as_bytes()).to_string();
    let old_names: Vec<String> = reader
        .header()
        .target_names()
        .iter()
        .map(|n| String::from_utf8_lossy(n).to_string())
        .collect();

    // Wrong reference is a hard error, not a contract fallback.
    check_sq_against_fai(&header_text, &old_fai)?;
    // A lift `@PG` under the old name is a hard error, not a remap fallback.
    reject_old_name_lift(&header_text).map_err(anyhow::Error::msg)?;

    let argv: String = recorded_argv();
    let output_header = build_output_header(&header_text, &new_fai, &argv, None);

    let contract = check_contract(&header_text);
    let (preset, mut contract_reason, chained) = match contract {
        Ok(c) => (Some(c.preset), None, c.chained),
        Err(r) => {
            if cfg.strict {
                bail!("input contract violation: {r}");
            }
            (None, Some(r), false)
        }
    };
    // A chained input is only accepted when the current lift's own
    // eligibility/safety options are at their defaults too.
    if contract_reason.is_none() && chained && !cfg.eligibility_is_default() {
        let r = "chained input: the current lift must use default eligibility/safety options"
            .to_string();
        if cfg.strict {
            bail!("input contract violation: {r}");
        }
        contract_reason = Some(r);
    }

    if let Some(reason) = contract_reason {
        let old_fa = Fasta::load(old_ref)?;
        let new_fa = Fasta::load(new_ref)?;
        let empty_diff = Diff::default();
        let mut reader2 = open_reader(in_cram, old_ref, threads)?;
        let ctx = Ctx {
            old: &old_fa,
            new: &new_fa,
            diff: &empty_diff,
            old_names: &old_names,
            new_tids: &new_tids,
            cfg,
            thr: 0.0,
        };
        let sink = FileSink::new(out_prefix, &output_header)?;
        let mut stats =
            run_all_defer(
                &mut reader2,
                sink,
                ctx,
                Some(reason),
                Reason::Contract,
                false,
                Instant::now(),
                &crate::progress::Progress::off(),
            )?
                .stats;
        stats.params = Params::new(cfg, cfg.max_de.unwrap_or(0.0));
        stats.params.threads = threads;
        stats.timings_s.total = t0.elapsed().as_secs_f64();
        write_stats(out_prefix, &stats)?;
        return Ok(stats);
    }

    let preset = preset.expect("preset present when contract ok");
    let thr = cfg.max_de.unwrap_or_else(|| de_threshold(preset));
    // Rule 7d scoring is derived from the resolved read preset; the config field
    // is runtime state, not a CLI option.
    let mut cfg_local = *cfg;
    cfg_local.scoring = scoring(preset);
    cfg_local.inherited_xl_ok = chained;
    let cfg = &cfg_local;
    let old = Fasta::load(old_ref)?;
    let new = Fasta::load(new_ref)?;

    let t_diff = Instant::now();
    let diff =
        crate::diff::compute_diff_with(&old, &new, old_ref, new_ref, preset, threads, cfg.diff_params())?
            .diff;
    let diff_s = t_diff.elapsed().as_secs_f64();

    let ctx = Ctx {
        old: &old,
        new: &new,
        diff: &diff,
        old_names: &old_names,
        new_tids: &new_tids,
        cfg,
        thr,
    };
    let t_stream = Instant::now();
    let sink = FileSink::new(out_prefix, &output_header)?;
    let res =
        run_translate_pass(&mut reader, sink, ctx, preset, t_stream, &crate::progress::Progress::off())?;
    let mut stats = res.stats;
    match res.outcome {
        PassOutcome::Completed => {}
        PassOutcome::ContractViolation(reason) => {
            if cfg.strict {
                bail!("input contract violation (mid-stream): {reason}");
            }
            drop(reader);
            let mut reader2 = open_reader(in_cram, old_ref, threads)?;
            let ctx2 = Ctx {
                old: &old,
                new: &new,
                diff: &diff,
                old_names: &old_names,
                new_tids: &new_tids,
                cfg,
                thr,
            };
            let sink2 = FileSink::new(out_prefix, &output_header)?;
            stats = run_all_defer(
                &mut reader2,
                sink2,
                ctx2,
                Some(reason),
                Reason::Contract,
                false,
                t_stream,
                &crate::progress::Progress::off(),
            )?
            .stats;
            // The contract was valid before the mid-stream violation, so the
            // preset is known even though every read is now deferred.
            stats.preset = preset.minimap2_name().to_string();
        }
    }
    stats.timings_s.diff = diff_s;
    stats.timings_s.total = t0.elapsed().as_secs_f64();
    stats.params = Params::new(cfg, thr);
    stats.params.threads = threads;
    write_stats(out_prefix, &stats)?;
    Ok(stats)
}

fn write_stats(prefix: &Path, stats: &Stats) -> Result<()> {
    let f = File::create(suffix(prefix, ".stats.json")).context("creating stats.json")?;
    let mut w = BufWriter::new(f);
    serde_json::to_writer_pretty(&mut w, &versioned(stats)).context("writing stats.json")?;
    w.flush().context("flushing stats.json")?;
    Ok(())
}
