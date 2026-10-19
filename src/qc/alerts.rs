//! `qc` alert counts (per type, per contig).
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use super::*;

/// Alert counts: trusted (`PASS`) calls, plus the `LowCov` calls split out.
/// `gray` is kept separate from the red/orange/yellow counts.
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct AlertCounts {
    pub red: u64,
    pub orange: u64,
    pub yellow: u64,
    pub gray: u64,
    pub lowcov: TierCounts,
}

impl AlertCounts {
    fn bump(&mut self, tier: Tier, lowcov: bool) {
        if lowcov {
            self.lowcov.bump(tier);
            return;
        }
        match tier {
            Tier::Red => self.red += 1,
            Tier::Orange => self.orange += 1,
            Tier::Yellow => self.yellow += 1,
            Tier::Gray => self.gray += 1,
            Tier::None => {}
        }
    }
}

/// Tiered calls by type: site types plus tiered region labels.
pub(crate) const ALERT_TYPES: [&str; 13] = [
    "SNP",
    "INS",
    "DEL",
    "snp_cluster",
    "large_del",
    "large_ins",
    "clip_hotspot",
    "depth_low",
    "depth_half",
    "depth_high",
    "depth_step",
    "linked_haplotype",
    "split_junction",
];

/// Alert counts (per contig or whole run): the `PASS` calls, the `LowCov`
/// calls and the same split per call type. The per-contig counts sum to the
/// top-level ones.
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct Alerts {
    #[serde(flatten)]
    pub counts: AlertCounts,
    pub by_type: std::collections::BTreeMap<String, AlertCounts>,
}

impl Alerts {
    fn bump(&mut self, label: &str, tier: Tier, lowcov: bool) {
        self.counts.bump(tier, lowcov);
        self.by_type.entry(label.to_string()).or_default().bump(tier, lowcov);
    }
}

/// Alert counts per contig and top-level. Tiered calls are VCF site records
/// (by type) and tiered BED regions (`snp_cluster`, `large_*`, by label). A
/// site or large-indel region whose FILTER is `LowCov` is counted separately.
pub(crate) fn collect_alerts(ncontigs: usize, cands: &[Candidate], regions: &[Region]) -> (Vec<Alerts>, Alerts) {
    let fresh = || Alerts {
        by_type: ALERT_TYPES.iter().map(|t| (t.to_string(), AlertCounts::default())).collect(),
        ..Alerts::default()
    };
    let mut per: Vec<Alerts> = (0..ncontigs).map(|_| fresh()).collect();
    let mut top = fresh();
    for c in cands {
        if c.tier.tier == Tier::None {
            continue;
        }
        let lowcov = c.lowcov & 1 != 0;
        per[c.tid].bump(c.vtype.name(), c.tier.tier, lowcov);
        top.bump(c.vtype.name(), c.tier.tier, lowcov);
    }
    for r in regions {
        // `depth_segment` is a derived paired interval (output only), not a call.
        if r.label == "depth_segment" {
            continue;
        }
        let tier = match r.tier {
            "red" => Tier::Red,
            "orange" => Tier::Orange,
            "yellow" => Tier::Yellow,
            _ => continue,
        };
        per[r.tid].bump(r.label, tier, r.lowcov);
        top.bump(r.label, tier, r.lowcov);
    }
    (per, top)
}


#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal candidate for the alert tests.
    fn mk_cand(tid: usize, vtype: VarType, tier: Tier, lowcov: u8) -> Candidate {
        Candidate {
            tid,
            pos0: 0,
            origin: false,
            vtype,
            len: 1,
            ref_allele: vec![b'A'],
            alt_allele: vec![b'C'],
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
            class: 0,
            grouped: 1,
            ln_p: 0.0,
            q: 1.0,
            wilson: 0.0,
            strand_ok: true,
            significant: true,
            lowcov,
            linked_group: None,
            unlinked: false,
            tier: TierOut { tier, reason: "minority", mixed: false, strand: false },
        }
    }
    fn region(tid: usize, label: &'static str, tier: &'static str, lowcov: bool) -> Region {
        Region {
            tid,
            start: 0,
            end: 10,
            label,
            score: 1,
            tier,
            reason: String::new(),
            depth_ratio: None,
            hiq_share: None,
            extra: String::new(),
            lowcov,
            link: None,
            clip_left: None,
            clip_right: None,
        }
    }
    #[test]
    fn top_level_alerts_sum_over_contigs() {
        let cands = vec![
            // Verified majority at low depth: tier RED, FILTER LowCov.
            mk_cand(0, VarType::Snp, Tier::Red, 1),
            // rel-only site: tier ORANGE, FILTER PASS.
            mk_cand(0, VarType::Ins, Tier::Orange, 2),
            mk_cand(1, VarType::Snp, Tier::Yellow, 0),
            mk_cand(1, VarType::Del, Tier::None, 0),
        ];
        let regions = vec![
            region(0, "large_del", "red", false),
            region(1, "large_del", "red", true), // LowCov large-del region
            region(1, "snp_cluster", "orange", false),
            region(1, "not_assessable", "none", false),
        ];
        let (per, top) = collect_alerts(2, &cands, &regions);
        let sum = |f: fn(&Alerts) -> u64| per.iter().map(f).sum::<u64>();
        assert_eq!(top.counts.red, sum(|a| a.counts.red));
        assert_eq!(top.counts.orange, sum(|a| a.counts.orange));
        assert_eq!(top.counts.yellow, sum(|a| a.counts.yellow));
        assert_eq!(top.counts.lowcov.red, sum(|a| a.counts.lowcov.red));
        assert_eq!(top.counts.lowcov.orange, sum(|a| a.counts.lowcov.orange));
        // PASS: 1 red (large_del) + 1 orange (rel-only INS) + 1 orange cluster + 1 yellow.
        assert_eq!((top.counts.red, top.counts.orange, top.counts.yellow), (1, 2, 1));
        // LowCov: the low-depth SNP (red) and the low-depth large-del region (red).
        assert_eq!((top.counts.lowcov.red, top.counts.lowcov.orange), (2, 0));
        // Untiered regions are annotations, not alerts.
        assert_eq!(top.by_type["SNP"].lowcov.red, 1);
        assert_eq!(top.by_type["INS"].orange, 1);
        assert_eq!(top.by_type["large_del"].red, 1);
        assert_eq!(top.by_type["large_del"].lowcov.red, 1);
        assert_eq!(top.by_type["snp_cluster"].orange, 1);
    }
}
