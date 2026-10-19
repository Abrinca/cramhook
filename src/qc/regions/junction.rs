//! `qc` split-read junction regions.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use crate::qc::*;
use crate::scan::*;

/// One split-read junction cluster: canonical breakpoints with their sides,
/// orientation class and distinct-read support (split by the reads' strand).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct JunctionCluster {
    /// Lower breakpoint and which side of it the read continues on.
    pub(crate) p1: u64,
    pub(crate) right1: bool,
    /// Upper breakpoint and its side.
    pub(crate) p2: u64,
    pub(crate) right2: bool,
    pub(crate) reads: u64,
    /// Distinct reads in the canonical traversal direction.
    pub(crate) fwd: u64,
    /// Distinct reads in the reverse-complement traversal direction.
    pub(crate) rev: u64,
}

/// Readable orientation class from the two breakpoint sides: a read leaving the
/// reference right of `p1` and continuing left of `p2` is a **tandem**/circle
/// (end→start) join, the mirror is a **deletion** (start→end) jump, and both
/// segments on the left / on the right are the two **inversion** classes.
fn orient_name(right1: bool, right2: bool) -> &'static str {
    match (right1, right2) {
        (true, false) => "tandem",
        (false, true) => "deletion",
        (false, false) => "inverted_l",
        (true, true) => "inverted_r",
    }
}

/// Leader (single-linkage) clusters of junction events: each event joins the
/// first cluster of the **same orientation class** (both breakpoint sides) whose
/// breakpoints are both within `tol`; support is the number of distinct reads
/// (forward and reverse strands merged). Events are sorted first, so the result
/// does not depend on the order the pass collected them.
pub(crate) fn junction_clusters(events: &[circularity::JunctionEvent], tol: u64) -> Vec<JunctionCluster> {
    let mut evs: Vec<&circularity::JunctionEvent> = events.iter().collect();
    evs.sort_by(|x, y| {
        (x.right1, x.right2, x.p1, x.p2, x.fwd, x.read).cmp(&(y.right1, y.right2, y.p1, y.p2, y.fwd, y.read))
    });
    let mut clusters: Vec<(u64, bool, u64, bool, Vec<(u32, bool)>)> = Vec::new();
    for e in evs {
        let mut placed = false;
        for c in clusters.iter_mut() {
            if c.1 == e.right1
                && c.3 == e.right2
                && c.0.abs_diff(e.p1) <= tol
                && c.2.abs_diff(e.p2) <= tol
            {
                c.4.push((e.read, e.fwd));
                placed = true;
                break;
            }
        }
        if !placed {
            clusters.push((e.p1, e.right1, e.p2, e.right2, vec![(e.read, e.fwd)]));
        }
    }
    let mut out: Vec<JunctionCluster> = clusters
        .into_iter()
        .map(|(p1, right1, p2, right2, mut reads)| {
            // One entry per distinct read (a read has one traversal direction).
            reads.sort_unstable();
            reads.dedup_by_key(|x| x.0);
            let fwd = reads.iter().filter(|x| x.1).count() as u64;
            let rev = reads.len() as u64 - fwd;
            JunctionCluster { p1, right1, p2, right2, reads: reads.len() as u64, fwd, rev }
        })
        .collect();
    out.sort_by(|x, y| (x.p1, x.right1, x.p2, x.right2).cmp(&(y.p1, y.right1, y.p2, y.right2)));
    out
}

/// Junction tier: below the support/ratio gate `None`; otherwise the gate is
/// the high bar (every reported junction is "highly significant"), so the tier
/// is the support fraction `f`: red when `f > 0.5`, orange otherwise.
fn junction_tier(reads: u64, frac: f64, p: &QcParams) -> Option<&'static str> {
    if reads < p.junction_min_reads as u64 || frac < p.junction_yellow {
        return None;
    }
    Some(if frac > 0.5 { "red" } else { "orange" })
}

/// Split-read junction regions (`split_junction`, reason `junction`): clusters
/// of reads joining two distant loci, gated on support and support/depth. Two
/// regions per cluster, one at each breakpoint (±`--junction-tol`). The fraction
/// is the support over the mean depth at the two breakpoints; tiers mirror the
/// clip thresholds.
pub(crate) fn junction_regions(
    tid: usize,
    c: &ContigCounts,
    ctx: &ContigCtx,
    p: &QcParams,
    events: &[circularity::JunctionEvent],
    junc_starts: &[u32],
    name: &str,
    out: &mut Vec<Region>,
) {
    let len = c.len as u64;
    if len == 0 || events.is_empty() {
        return;
    }
    let dp_pref = prefix_u32(&ctx.dp);
    let tol = p.junction_tol as u64;
    for cl in junction_clusters(events, tol) {
        if cl.reads < p.junction_min_reads as u64 {
            continue;
        }
        let last = (len - 1) as usize;
        let dlo = ctx.dp.get((cl.p1 as usize).min(last)).copied().unwrap_or(0) as f64;
        let dhi = ctx.dp.get((cl.p2 as usize).min(last)).copied().unwrap_or(0) as f64;
        let depth = (dlo + dhi) / 2.0;
        if depth <= 0.0 {
            continue;
        }
        let frac = cl.reads as f64 / depth;
        let Some(tier) = junction_tier(cl.reads, frac, p) else { continue };
        // Significance: junction reads vs the local junction rate in the window
        // around the junction midpoint. Only the two breakpoint neighbourhoods
        // are excluded (not the whole span between them), so breakpoints far
        // apart keep a usable background. Depth comes from a per-contig prefix
        // sum and events from a per-contig sorted position list, so this is
        // O(log n) per cluster, not a rescan.
        let center = (cl.p1 + cl.p2) / 2;
        let excls = [
            (cl.p1.saturating_sub(tol), (cl.p1 + tol + 1).min(len)),
            (cl.p2.saturating_sub(tol), (cl.p2 + tol + 1).min(len)),
        ];
        let bg_cov = prefix_window_excl(&dp_pref, len, ctx.circular, center, SIG_SPAN, &excls);
        let bg_ev = count_window_excl(junc_starts, len, ctx.circular, center, SIG_SPAN, &excls);
        let rate = if bg_cov > 0 { bg_ev as f64 / bg_cov as f64 } else { 0.0 };
        let s = poisson_score(cl.reads, rate * depth);
        let span = cl.p2 - cl.p1;
        for (pos, partner) in [(cl.p1, cl.p2), (cl.p2, cl.p1)] {
            let start = pos.saturating_sub(tol).min(len - 1);
            let end = (pos + tol + 1).min(len).max(start + 1);
            let (dr, hq) = region_depth(c, ctx, start, end);
            let extra = format!(
                "partner={name}:{partner};reads={};frac={frac:.4};span={span};orient={};fwd={};rev={}",
                cl.reads,
                orient_name(cl.right1, cl.right2),
                cl.fwd,
                cl.rev
            ) + &sig_suffix(s, frac);
            out.push(Region {
                tid,
                start,
                end,
                label: "split_junction",
                score: cl.reads.min(1000) as u32,
                tier,
                reason: "junction".to_string(),
                depth_ratio: dr,
                hiq_share: hq,
                extra,
                lowcov: false,
                link: Some((cl.p1, cl.p2)),
                clip_left: None,
                clip_right: None,
            });
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::qc::test_util::*;

    #[test]
    fn junction_tier_boundaries() {
        let p = QcParams::default();
        // Below the support gate.
        assert_eq!(junction_tier(4, 0.9, &p), None);
        // Below the fraction gate.
        assert_eq!(junction_tier(5, 0.019, &p), None);
        // Every reported junction is at the high bar: orange up to f = 0.5,
        // red above (no yellow).
        assert_eq!(junction_tier(5, 0.02, &p), Some("orange"));
        assert_eq!(junction_tier(5, 0.5, &p), Some("orange"));
        assert_eq!(junction_tier(5, 0.51, &p), Some("red"));
    }
    #[test]
    fn junction_clusters_merge_strands_and_count() {
        // (p1, right1, p2, right2, fwd) — a tandem join at (2_000, 12_000).
        let ev = |p1, right1, p2, right2, fwd, read: u32| circularity::JunctionEvent {
            p1,
            right1,
            p2,
            right2,
            fwd,
            read,
        };
        let events = vec![
            ev(2_000, true, 12_000, false, false, 1), // reverse-complement traversal
            ev(2_000, true, 12_000, false, true, 2),  // canonical traversal
            ev(2_010, true, 12_040, false, false, 3), // within tol on both breakpoints
            ev(2_200, true, 12_000, false, false, 4), // too far at the first breakpoint
            ev(2_000, true, 12_000, false, false, 1), // same read: no extra support
        ];
        let cl = junction_clusters(&events, 50);
        assert_eq!(cl.len(), 2, "{cl:?}");
        let direct = cl.iter().find(|c| c.reads == 3).unwrap();
        assert_eq!((direct.p1, direct.right1, direct.p2, direct.right2), (2_000, true, 12_000, false));
        assert_eq!((direct.fwd, direct.rev), (1, 2));
        assert!(cl.iter().any(|c| c.reads == 1 && c.p1 == 2_200));
    }
    /// Different junction classes at the same breakpoints stay separate: a
    /// tandem (end->start) vs a deletion (start->end), and the two inversion
    /// classes.
    #[test]
    fn junction_clusters_keep_classes_distinct() {
        let ev = |p1, right1, p2, right2, read: u32| circularity::JunctionEvent {
            p1,
            right1,
            p2,
            right2,
            fwd: false,
            read,
        };
        // Same positions, four classes.
        let events = vec![
            ev(2_000, true, 12_000, false, 1),  // tandem
            ev(2_000, false, 12_000, true, 2),  // deletion
            ev(2_000, false, 12_000, false, 3), // inverted_l
            ev(2_000, true, 12_000, true, 4),   // inverted_r
        ];
        let cl = junction_clusters(&events, 50);
        assert_eq!(cl.len(), 4, "{cl:?}");
        let names: Vec<&str> = cl.iter().map(|c| orient_name(c.right1, c.right2)).collect();
        assert!(names.contains(&"tandem") && names.contains(&"deletion"), "{names:?}");
        assert!(names.contains(&"inverted_l") && names.contains(&"inverted_r"), "{names:?}");
        assert!(cl.iter().all(|c| c.reads == 1));
    }
    /// Overlapping breakpoint neighbourhoods (`--junction-tol 1000`, breakpoints
    /// 1 kb apart): the reviewer repro. Must not underflow; with no background
    /// coverage the score is the lambda-floor value.
    #[test]
    fn junction_overlapping_breakpoint_exclusions() {
        let len = 80_000usize;
        let mut dp = vec![0u32; len];
        for x in dp.iter_mut().take(1000).skip(900) {
            *x = 5;
        }
        for x in dp.iter_mut().take(2100).skip(2000) {
            *x = 5;
        }
        let (c, ctx) = mk_contig(len, vec![[0u32; 2]; len], dp);
        let p = QcParams { junction_tol: 1000, ..QcParams::default() };
        let events: Vec<circularity::JunctionEvent> = (0..5)
            .map(|r| circularity::JunctionEvent { p1: 1000, right1: false, p2: 2000, right2: true, fwd: true, read: r })
            .collect();
        let starts: Vec<u32> = events.iter().map(|e| e.p1 as u32).collect();
        let mut out = Vec::new();
        junction_regions(0, &c, &ctx, &p, &events, &starts, "chr", &mut out);
        assert_eq!(out.len(), 2, "{out:?}");
        let want = poisson_score(5, 0.5);
        for r in &out {
            let s: f64 = r.extra.split(";s=").nth(1).unwrap().split(';').next().unwrap().parse().unwrap();
            assert!((s - want).abs() < 0.01, "s={s} want={want} {}", r.extra);
        }
    }
    /// Junction score with breakpoints > 50 kb apart: the background is the
    /// between-region (not an empty window), so the score reflects it.
    #[test]
    fn junction_score_far_breakpoints() {
        let len = 200_000usize;
        let mut dp = vec![1u32; len];
        dp[10_000] = 1000;
        dp[90_000] = 1000;
        let (c, ctx) = mk_contig(len, vec![[0u32; 2]; len], dp);
        let p = QcParams::default();
        // One real cluster: 20 reads joining 10_000 -> 90_000 (frac 0.02).
        let mut events: Vec<circularity::JunctionEvent> = (0..20)
            .map(|r| circularity::JunctionEvent { p1: 10_000, right1: true, p2: 90_000, right2: false, fwd: true, read: r })
            .collect();
        // 200 unrelated background junctions inside the between window, each its
        // own cluster (spaced > tol apart), 1 read each.
        for i in 0..200u64 {
            let a = 25_000 + i * 250;
            events.push(circularity::JunctionEvent { p1: a, right1: true, p2: a + 40_000, right2: false, fwd: true, read: 100 + i as u32 });
        }
        let mut starts: Vec<u32> = events.iter().map(|e| e.p1 as u32).collect();
        starts.sort_unstable();
        let mut out = Vec::new();
        junction_regions(0, &c, &ctx, &p, &events, &starts, "chr", &mut out);
        assert_eq!(out.len(), 2, "one cluster -> two breakpoint regions: {out:?}");
        // bg_cov = 50 000 * 1 (window [25k,75k], depth 1, no exclusion overlap);
        // bg_ev = 200; rate = 0.004; depth = 1000 -> lambda = 4.
        let want = poisson_score(20, 4.0);
        for r in &out {
            let s: f64 = r.extra.split(";s=").nth(1).unwrap().split(';').next().unwrap().parse().unwrap();
            assert!((s - want).abs() < 0.01, "s={s} want={want} {}", r.extra);
        }
    }
    /// Many junction score lookups must be sub-quadratic (the old code rescanned
    /// every event per cluster).
    #[test]
    fn junction_score_lookups_scale() {
        let n = 200_000u32;
        let len = n as u64 * 4 + 10;
        let starts: Vec<u32> = (0..n).map(|i| i * 4).collect();
        let dp = prefix_u32(&vec![10u32; len as usize]);
        let t = std::time::Instant::now();
        let mut acc = 0u64;
        for i in 0..200_000u64 {
            let center = (i * 7) % (len - 1);
            let excl = [(center.saturating_sub(200), center + 201)];
            acc += count_window_excl(&starts, len, false, center, 25_000, &excl);
            acc += prefix_window_excl(&dp, len, false, center, 25_000, &excl);
        }
        assert!(acc > 0);
        let el = t.elapsed().as_secs_f64();
        assert!(el < 1.0, "junction score lookups took {el:.3}s");
    }
}
