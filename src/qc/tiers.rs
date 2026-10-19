//! `qc` site tiers, tier classes and the tier rule.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use super::*;

/// Site tier. `Gray` sits below `Yellow`: it is a reported call (written to the
/// VCF), but a homopolymer-ambiguous indel, not a finding member and not part
/// of the red/orange/yellow alert or tier counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Tier {
    None,
    Gray,
    Yellow,
    Orange,
    Red,
}

impl Tier {
    pub fn name(self) -> &'static str {
        match self {
            Tier::None => "none",
            Tier::Gray => "gray",
            Tier::Yellow => "yellow",
            Tier::Orange => "orange",
            Tier::Red => "red",
        }
    }
}

/// Calls per tier (gray kept apart; it is not one of the red/orange/yellow
/// counts).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize)]
pub struct TierCounts {
    pub red: u64,
    pub orange: u64,
    pub yellow: u64,
    pub gray: u64,
}

impl TierCounts {
    pub(crate) fn bump(&mut self, t: Tier) {
        match t {
            Tier::Red => self.red += 1,
            Tier::Orange => self.orange += 1,
            Tier::Yellow => self.yellow += 1,
            Tier::Gray => self.gray += 1,
            Tier::None => {}
        }
    }
}

/// Indel context buckets: HP run length 1..7, 8+, STR.
pub(crate) const NB: usize = 9;
pub(crate) const BUCKETS: [&str; NB] = ["HP1", "HP2", "HP3", "HP4", "HP5", "HP6", "HP7", "HP8+", "STR"];
/// Classes: SNP, INS buckets, DEL buckets.
pub(crate) const NCLASS: usize = 1 + 2 * NB;

pub(crate) fn class_name(c: usize) -> String {
    if c == 0 {
        "SNP".to_string()
    } else if c <= NB {
        format!("INS_{}", BUCKETS[c - 1])
    } else {
        format!("DEL_{}", BUCKETS[c - 1 - NB])
    }
}

/// Indel bucket from the forward homopolymer run and STR span at anchor + 1.
pub(crate) fn indel_bucket(hp_fwd: u8, str_span: u16) -> usize {
    if hp_fwd >= 2 {
        (hp_fwd.min(8) - 1) as usize
    } else if str_span >= 6 {
        8
    } else {
        0
    }
}

/// Strand-balance gate.
pub(crate) fn strand_ok(altf: u32, altr: u32, dpf: u32, dpr: u32, min_dp: u32, ratio: f64) -> bool {
    if altf == 0 || altr == 0 || dpf < min_dp || dpr < min_dp {
        return false;
    }
    let a = altf as f64 / dpf as f64;
    let b = altr as f64 / dpr as f64;
    a.min(b) >= ratio * a.max(b)
}

/// Inputs to the tier rule.
#[derive(Clone, Copy, Debug)]
pub(crate) struct TierIn {
    pub af: f64,
    pub alt: u32,
    pub wilson: f64,
    pub strand_ok: bool,
    pub significant: bool,
}

/// Tier rule result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct TierOut {
    pub tier: Tier,
    pub reason: &'static str,
    pub mixed: bool,
    pub strand: bool,
}

/// The tier rule table (in order):
/// 1. strand-balanced AF > 0.5 → reported regardless of the q-value: red
///    (`majority`) if Wilson LB > `red_min_wilson` and alt >= `red_min_alt`,
///    else yellow (`weak_majority`);
/// 2. not significant → none;
/// 3. AF >= `mixed_min_af`, strand gate failed → yellow + STRAND (`strand_bias`);
/// 4. `mixed_min_af` <= AF <= 0.5 → orange + MIXED (`mixed`);
/// 5. AF < `mixed_min_af`, strand gate failed → none (`strand_bias`, STRAND
///    flag) if `yellow_strand_gate`, else yellow (`minority`) with STRAND;
/// 6. AF < `mixed_min_af` → yellow (`minority`).
///
/// The rule is "red = high bar (Wilson LB > 0.5, alt >= `red_min_alt`) and
/// f > 0.5; orange = high bar and f <= 0.5; yellow = reported (low bar) but
/// below the high bar", with `f` = AF. The majority reporting gate is
/// independent of the q-value (a strand-balanced majority is never dropped).
/// The mixed band keeps its old orange meaning between `mixed_min_af` and 0.5.
pub(crate) fn tier_of(t: &TierIn, p: &QcParams) -> TierOut {
    let out = |tier, reason, mixed, strand| TierOut { tier, reason, mixed, strand };
    if t.strand_ok && t.af > 0.5 {
        if t.wilson > p.red_min_wilson && t.alt >= p.red_min_alt {
            return out(Tier::Red, "majority", false, false);
        }
        return out(Tier::Yellow, "weak_majority", false, false);
    }
    if !t.significant {
        return out(Tier::None, "not_significant", false, !t.strand_ok);
    }
    if t.af >= p.mixed_min_af && !t.strand_ok {
        return out(Tier::Yellow, "strand_bias", false, true);
    }
    if t.af >= p.mixed_min_af {
        return out(Tier::Orange, "mixed", true, false);
    }
    if !t.strand_ok && p.yellow_strand_gate {
        return out(Tier::None, "strand_bias", false, true);
    }
    out(Tier::Yellow, "minority", false, !t.strand_ok)
}


#[cfg(test)]
mod tests {
    use super::*;

    fn tin(af: f64, alt: u32, wilson: f64, sb: bool, sig: bool) -> TierIn {
        TierIn { af, alt, wilson, strand_ok: sb, significant: sig }
    }
    #[test]
    fn tier_rule_table() {
        let p = QcParams::default();
        let cases: &[(TierIn, Tier, &str, bool, bool)] = &[
            // Majority, strong evidence.
            (tin(0.9, 45, 0.8, true, true), Tier::Red, "majority", false, false),
            // Majority without significance is still red.
            (tin(0.9, 45, 0.8, true, false), Tier::Red, "majority", false, false),
            // Red no longer needs AF > 0.6: Wilson LB > 0.5 and alt >= 5 suffice.
            (tin(0.55, 10, 0.51, true, true), Tier::Red, "majority", false, false),
            // AF > 0.5 but alt < red_min_alt: weak majority -> yellow.
            (tin(1.0, 4, 0.51, true, true), Tier::Yellow, "weak_majority", false, false),
            // AF > 0.5 but Wilson bound too low: weak majority -> yellow.
            (tin(0.65, 13, 0.43, true, true), Tier::Yellow, "weak_majority", false, false),
            // Majority reporting gate is independent of q: a strand-balanced
            // AF > 0.5 site with q not significant is still yellow, not dropped.
            (tin(0.8, 8, 0.49, true, false), Tier::Yellow, "weak_majority", false, false),
            // Majority, strand-biased: yellow + STRAND.
            (tin(0.9, 45, 0.8, false, true), Tier::Yellow, "strand_bias", false, true),
            (tin(0.9, 45, 0.8, false, false), Tier::None, "not_significant", false, true),
            // Mixed band (0.40..0.50) stays orange.
            (tin(0.5, 25, 0.36, true, true), Tier::Orange, "mixed", true, false),
            (tin(0.40, 20, 0.27, true, true), Tier::Orange, "mixed", true, false),
            (tin(0.5, 25, 0.36, false, true), Tier::Yellow, "strand_bias", false, true),
            (tin(0.5, 25, 0.36, true, false), Tier::None, "not_significant", false, false),
            // Minority.
            (tin(0.2, 10, 0.11, true, true), Tier::Yellow, "minority", false, false),
            // Minority failing the strand gate: not called (yellow strand gate).
            (tin(0.2, 10, 0.11, false, true), Tier::None, "strand_bias", false, true),
            (tin(0.39, 20, 0.26, false, true), Tier::None, "strand_bias", false, true),
            (tin(0.2, 10, 0.11, true, false), Tier::None, "not_significant", false, false),
        ];
        for (i, (t, tier, reason, mixed, strand)) in cases.iter().enumerate() {
            let o = tier_of(t, &p);
            assert_eq!(o, TierOut { tier: *tier, reason, mixed: *mixed, strand: *strand }, "case {i}");
        }
        // --no-yellow-strand-gate: old rule, minority stays yellow with STRAND;
        // nothing else changes.
        let old = QcParams { yellow_strand_gate: false, ..QcParams::default() };
        for (i, (t, tier, reason, mixed, strand)) in cases.iter().enumerate() {
            let want = if *reason == "strand_bias" && t.af < old.mixed_min_af {
                TierOut { tier: Tier::Yellow, reason: "minority", mixed: false, strand: true }
            } else {
                TierOut { tier: *tier, reason, mixed: *mixed, strand: *strand }
            };
            assert_eq!(tier_of(t, &old), want, "old rule, case {i}");
        }
    }
    #[test]
    fn strand_gate() {
        assert!(strand_ok(10, 10, 20, 20, 3, 1.0 / 3.0));
        assert!(!strand_ok(10, 0, 20, 20, 3, 1.0 / 3.0));
        assert!(!strand_ok(10, 2, 20, 20, 3, 1.0 / 3.0));
        assert!(strand_ok(9, 3, 20, 20, 3, 1.0 / 3.0));
        assert!(!strand_ok(2, 1, 20, 2, 3, 1.0 / 3.0));
    }
    #[test]
    fn indel_buckets() {
        assert_eq!(indel_bucket(1, 0), 0);
        assert_eq!(indel_bucket(2, 0), 1);
        assert_eq!(indel_bucket(8, 0), 7);
        assert_eq!(indel_bucket(30, 0), 7);
        assert_eq!(indel_bucket(1, 6), 8);
        assert_eq!(indel_bucket(1, 5), 0);
        assert_eq!(class_name(0), "SNP");
        assert_eq!(class_name(1), "INS_HP1");
        assert_eq!(class_name(9), "INS_STR");
        assert_eq!(class_name(10), "DEL_HP1");
        assert_eq!(class_name(17), "DEL_HP8+");
    }
}
