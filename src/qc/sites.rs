//! `qc` site types and low-coverage reason bits.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use super::*;

/// Low-coverage reason bits for one assessable position: bit 0 = `abs`
/// (DP < `low_abs`), bit 1 = `rel` (DP < `low_rel` × mode); 0 = not low.
pub(crate) fn lowcov_bits(dp: u32, mode: u32, low_abs: u32, low_rel: f64) -> u8 {
    let mut m = 0u8;
    if dp < low_abs {
        m |= 1;
    }
    if (dp as f64) < low_rel * mode as f64 {
        m |= 2;
    }
    m
}

/// Reason label for low-coverage bits (`abs`, `rel` or `abs,rel`).
pub(crate) fn lowcov_reason(bits: u8) -> String {
    match bits {
        1 => "abs".to_string(),
        2 => "rel".to_string(),
        3 => "abs,rel".to_string(),
        _ => String::new(),
    }
}

/// A tested allele (pre-filter passed).
#[derive(Clone, Debug)]
pub(crate) struct Candidate {
    pub tid: usize,
    /// 0-based: SNP position, indel anchor (0 for `origin` indels).
    pub pos0: u32,
    /// Indel before reference position 0 (right-padded VCF form).
    pub origin: bool,
    pub vtype: VarType,
    pub len: u32,
    pub ref_allele: Vec<u8>,
    pub alt_allele: Vec<u8>,
    pub dp: u32,
    pub dpf: u32,
    pub dpr: u32,
    pub altf: u32,
    pub altr: u32,
    pub reff: u32,
    pub refr: u32,
    pub hp: u32,
    pub str_span: u32,
    pub mq0f: f64,
    pub class: usize,
    /// Number of exact indel alleles merged into this candidate (1 when not
    /// grouped); only > 1 for grouped INS alleles.
    pub grouped: u32,
    pub ln_p: f64,
    pub q: f64,
    pub wilson: f64,
    pub strand_ok: bool,
    pub significant: bool,
    /// Low-coverage reason bits at the position (0, `abs`, `rel`, `abs,rel`).
    pub lowcov: u8,
    /// Haplotype group id from the linkage pass (SNP candidates).
    pub linked_group: Option<u32>,
    /// Tested against >= 2 partners and linked to none (flagged or dropped per `--unlinked-action`).
    pub unlinked: bool,
    pub tier: TierOut,
}

impl Candidate {
    pub fn alt(&self) -> u32 {
        self.altf + self.altr
    }
    /// Alt fraction. Support is counted only from reads that cover the
    /// anchor (see `scan`), at most once per read, so alt <= DP holds; the
    /// clamp is defensive.
    pub fn af(&self) -> f64 {
        debug_assert!(self.alt() <= self.dp, "alt {} > dp {} at {}", self.alt(), self.dp, self.pos0);
        if self.dp == 0 { 0.0 } else { (self.alt() as f64 / self.dp as f64).min(1.0) }
    }
}

/// A BED region finding.
#[derive(Clone, Debug)]
pub(crate) struct Region {
    pub tid: usize,
    pub start: u64,
    pub end: u64,
    pub label: &'static str,
    pub score: u32,
    pub tier: &'static str,
    pub reason: String,
    pub depth_ratio: Option<f64>,
    pub hiq_share: Option<f64>,
    pub extra: String,
    /// Large-indel regions with mean DP < `low_abs`: `lowcov=1` in the BED
    /// extra column and counted in `alerts.lowcov`.
    pub lowcov: bool,
    /// Two positions that belong to one event (`split_junction`: the two
    /// breakpoints; `depth_segment`: its two ends). Findings link the groups at
    /// the two ends (see [`collect_findings`]). `None` for every other region.
    pub link: Option<(u64, u64)>,
    /// `clip_hotspot` only: the two sides' candidate clusters (count-weighted
    /// median, count, ratio and aligned-side depth). The cross-contig link uses
    /// the matched side's own `pos`/`depth`, so its `f` equals that side's ratio.
    /// `None` for every other region label.
    pub clip_left: Option<ClipSide>,
    pub clip_right: Option<ClipSide>,
}

/// Background fit per class (for the JSON).
#[derive(Clone, Debug, serde::Serialize)]
pub struct ClassFit {
    pub class: String,
    #[serde(flatten)]
    pub fit: BbFit,
    /// Positions used for the fit (AF below the fit limit).
    pub fit_sites: u64,
    /// Assessable positions of this class (BH m).
    pub m: u64,
    /// True if the pooled INS/DEL fit was used.
    pub fallback: bool,
    pub candidates: u64,
}

/// Per-contig derived arrays.
pub(crate) struct ContigCtx {
    pub hp_run: Vec<u8>,
    pub hp_fwd: Vec<u8>,
    pub str_span: Vec<u16>,
    pub dp: Vec<u32>,
    /// 0 = assessable, 1 = ref N, 2 = low depth, 3 = MAPQ0.
    pub na: Vec<u8>,
    /// Low-coverage reason bits per position (0 = not low, see
    /// [`lowcov_bits`]).
    pub lowcov: Vec<u8>,
    pub circular: bool,
    /// Which source decided the wrap: `verdict`, `header` or `flag`
    /// (JSON `circular_by`).
    pub circular_by: &'static str,
    pub mode: u32,
    pub mean: f64,
    /// Depth-band skip reason for this contig (JSON `depth_skipped`).
    pub depth_skipped: Option<&'static str>,
}

pub(crate) const NA_REASONS: [&str; 4] = ["", "ref_n", "low_depth", "mapq0"];


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn af_never_exceeds_one() {
        let mut c = Candidate {
            tid: 0,
            pos0: 0,
            origin: false,
            vtype: VarType::Ins,
            len: 1,
            ref_allele: vec![],
            alt_allele: vec![],
            dp: 10,
            dpf: 5,
            dpr: 5,
            altf: 5,
            altr: 5,
            reff: 0,
            refr: 0,
            hp: 1,
            str_span: 0,
            mq0f: 0.0,
            class: 1,
            grouped: 1,
            ln_p: 0.0,
            q: 1.0,
            wilson: 0.0,
            strand_ok: true,
            significant: true,
            lowcov: 0,
            linked_group: None,
            unlinked: false,
            tier: TierOut { tier: Tier::None, reason: "", mixed: false, strand: false },
        };
        assert_eq!(c.af(), 1.0);
        c.dp = 0;
        c.altf = 0;
        c.altr = 0;
        assert_eq!(c.af(), 0.0);
    }
    #[test]
    fn lowcov_classification_and_boundaries() {
        // mode 40, low_abs 20, low_rel 0.25 => rel threshold 10x.
        let (a, rel, mode) = (20u32, 0.25f64, 40u32);
        // Neither.
        assert_eq!(lowcov_bits(40, mode, a, rel), 0);
        // abs only (below 20 but at/above the 10x rel line).
        assert_eq!(lowcov_bits(12, mode, a, rel), 1);
        assert_eq!(lowcov_bits(10, mode, a, rel), 1);
        // rel only (at/above 20 but below 10x is impossible at mode 40; use a
        // high abs so the abs bit cannot fire).
        assert_eq!(lowcov_bits(8, mode, 5, rel), 2);
        // Both.
        assert_eq!(lowcov_bits(8, mode, a, rel), 3);
        // Boundaries: exactly at each threshold is *not* low.
        assert_eq!(lowcov_bits(a, mode, a, rel), 0); // dp == low_abs
        assert_eq!(lowcov_bits(10, mode, 5, rel), 0); // dp == low_rel * mode
        assert_eq!(lowcov_bits(a - 1, mode, a, rel), 1);
        assert_eq!(lowcov_bits(9, mode, 5, rel), 2);
    }
}
