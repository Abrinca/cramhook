//! `qc` clip-hotspot regions.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use crate::qc::*;
use crate::scan::*;

/// One side's clip cluster: total clip count, count-weighted median position
/// and the inclusive span of the clustered clip positions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ClipCluster {
    pub pos: u64,
    pub count: u64,
    pub start: u64,
    pub end: u64,
}

/// Single-linkage chains of clip positions (`counts[i] > 0`) whose consecutive
/// gaps are <= `merge`. The cluster position is the count-weighted median; the
/// span is the first..last clip position (inclusive). One linear pass.
pub(crate) fn clip_clusters(counts: &[u32], merge: u32) -> Vec<ClipCluster> {
    let positions: Vec<(u64, u64)> =
        counts.iter().enumerate().filter(|&(_, &c)| c > 0).map(|(i, &c)| (i as u64, c as u64)).collect();
    let mut out = Vec::new();
    let mut k = 0;
    while k < positions.len() {
        let start = positions[k].0;
        let mut total = positions[k].1;
        let mut m = k + 1;
        while m < positions.len() && positions[m].0 - positions[m - 1].0 <= merge as u64 {
            total += positions[m].1;
            m += 1;
        }
        let end = positions[m - 1].0;
        // Count-weighted median: the first position whose cumulative weight
        // reaches half the cluster total (ties resolved to the smaller pos).
        let mut cum = 0u64;
        let mut pos = start;
        for &(p, c) in &positions[k..m] {
            cum += c;
            if cum * 2 >= total {
                pos = p;
                break;
            }
        }
        out.push(ClipCluster { pos, count: total, start, end });
        k = m;
    }
    out
}

/// Mean depth over the aligned-side window of a clip cluster: left clips use
/// `[pos, pos + 10)` (the aligned side is to the right), right clips
/// `[pos - 9, pos + 1)` (aligned side to the left). Truncated to `[0, len)`;
/// `None` when the window is empty.
pub(crate) fn clip_aligned_depth(dp: &[u32], pos: u64, len: u64, left: bool) -> Option<f64> {
    let (lo, hi) = if left { (pos, pos + 10) } else { (pos.saturating_sub(9), pos + 1) };
    let (a, b) = (lo.min(len), hi.min(len));
    if b <= a {
        return None;
    }
    let s: u64 = (a..b).map(|x| dp[x as usize] as u64).sum();
    Some(s as f64 / (b - a) as f64)
}

/// Clip-hotspot candidate (low-bar) gate: enough reads and a fraction at or
/// above `clip_yellow`. Every candidate is reported (yellow at least).
pub(crate) fn clip_candidate_ok(count: u64, ratio: f64, p: &QcParams) -> bool {
    count >= p.clip_min_reads as u64 && ratio >= p.clip_yellow
}

/// Clip-hotspot tier from the significance score `s` (`-log10 p`, the larger
/// of the two sides) and the disagreement fraction `f`: red at the high bar
/// with `f > 0.5`, orange at the high bar otherwise, yellow below it (but
/// reported). `f` is the single-side ratio, or `(left + right) / depth` for a
/// two-sided hotspot.
pub(crate) fn clip_high_tier(s: f64, f: f64, p: &QcParams) -> &'static str {
    if s >= p.clip_high_s {
        if f > 0.5 { "red" } else { "orange" }
    } else {
        "yellow"
    }
}

/// One candidate clip cluster (candidate gate passed) with its aligned-side
/// depth and ratio.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ClipSide {
    /// Count-weighted median position.
    pub(crate) pos: u64,
    pub(crate) count: u64,
    /// Cluster count / aligned-side mean depth.
    pub(crate) ratio: f64,
    /// Mean depth over the aligned-side window.
    pub(crate) depth: f64,
    /// Inclusive span of the clustered clip positions.
    pub(crate) span: (u64, u64),
}

/// Candidate clip clusters of one side, after dropping contig ends and applying
/// the count/ratio gate.
fn clip_candidates(counts: &[u32], dp: &[u32], len: u64, p: &QcParams, left: bool) -> Vec<ClipSide> {
    let edge = p.clip_edge as u64;
    clip_clusters(counts, p.clip_merge)
        .into_iter()
        .filter(|cl| cl.pos >= edge && cl.pos + edge < len)
        .filter_map(|cl| {
            let depth = clip_aligned_depth(dp, cl.pos, len, left)?;
            let ratio = cl.count as f64 / depth;
            if !clip_candidate_ok(cl.count, ratio, p) {
                return None;
            }
            Some(ClipSide {
                pos: cl.pos,
                count: cl.count,
                ratio,
                depth,
                span: (cl.start, cl.end),
            })
        })
        .collect()
}

/// Emit one clip-hotspot region from a candidate left and/or right cluster.
#[allow(clippy::too_many_arguments)]
fn emit_clip(
    tid: usize,
    c: &ContigCounts,
    ctx: &ContigCtx,
    p: &QcParams,
    low_abs: u32,
    left: Option<&ClipSide>,
    right: Option<&ClipSide>,
    dp_pref: &[u64],
    left_pref: &[u64],
    right_pref: &[u64],
    out: &mut Vec<Region>,
) {
    let len = c.len as u64;
    let (mut start, mut end) = (u64::MAX, 0u64);
    for s in [left, right].into_iter().flatten() {
        start = start.min(s.span.0);
        end = end.max(s.span.1);
    }
    let start = start.min(len - 1);
    let end = (end + 1).clamp(start + 1, len);
    let (lc, rc) = (left.map_or(0, |s| s.count), right.map_or(0, |s| s.count));
    let (lr, rr) = (left.map_or(0.0, |s| s.ratio), right.map_or(0.0, |s| s.ratio));
    let depth = match (left, right) {
        (Some(l), Some(r)) => (l.depth + r.depth) / 2.0,
        (Some(s), None) | (None, Some(s)) => s.depth,
        (None, None) => return,
    };
    let depth_ratio = if ctx.mode > 0 { Some(depth / ctx.mode as f64) } else { None };
    let hiq = region_depth(c, ctx, start, end).1;
    let lowcov = depth < low_abs as f64;
    let reason = match (left.is_some(), right.is_some()) {
        (true, true) => "clip_both",
        (true, false) => "clip_left",
        _ => "clip_right",
    };
    // Significance of the side with the larger Poisson score; its clip ratio
    // is the disagreement fraction `f`. Window and call sums come from prefix
    // sums, so this is O(1) per call. The cluster exclusion is clipped to the
    // window pieces, so a cluster wider than the window cannot underflow;
    // `call_cov` stays the full cluster coverage.
    let side_score = |s: &ClipSide, clip_pref: &[u64]| -> f64 {
        let excl = [(s.span.0, s.span.1 + 1)];
        let ev = prefix_window_excl(clip_pref, len, ctx.circular, s.pos, SIG_SPAN, &excl);
        let cov = prefix_window_excl(dp_pref, len, ctx.circular, s.pos, SIG_SPAN, &excl);
        let rate = if cov > 0 { ev as f64 / cov as f64 } else { 0.0 };
        let call_cov = prefix_range(dp_pref, excl[0].0, excl[0].1);
        poisson_score(s.count, rate * call_cov as f64)
    };
    let mut sig = (f64::NEG_INFINITY, 0.0f64);
    if let Some(s) = left {
        sig = (side_score(s, left_pref), s.ratio);
    }
    if let Some(s) = right {
        let sc = side_score(s, right_pref);
        if sc > sig.0 {
            sig = (sc, s.ratio);
        }
    }
    // Disagreement fraction: for a two-sided hotspot each side holds about half
    // the reads at the break, so f = (left + right clips) / depth; a single
    // side uses its own ratio.
    let f = if left.is_some() && right.is_some() && depth > 0.0 {
        (lc + rc) as f64 / depth
    } else {
        sig.1
    };
    let tier = clip_high_tier(sig.0, f, p);
    let mut extra = format!("left={lc};right={rc};ratio_left={lr:.3};ratio_right={rr:.3}");
    if lowcov {
        extra.push_str(";lowcov=1");
    }
    extra.push_str(&sig_suffix(sig.0, f));
    out.push(Region {
        tid,
        start,
        end,
        label: "clip_hotspot",
        score: lc.max(rc).min(1000) as u32,
        tier,
        reason: reason.to_string(),
        depth_ratio,
        hiq_share: hiq,
        extra,
        lowcov,
        link: None,
        clip_left: left.copied(),
        clip_right: right.copied(),
    });
}

/// Clip hotspots: per-side clusters of clips of at least `clip_min_len`, ratio
/// against the aligned-side depth, tiers and left/right pairing. Clusters near
/// a contig end are skipped (origin crossings / read-end jitter).
pub(crate) fn clip_hotspot_regions(
    tid: usize,
    c: &ContigCounts,
    ctx: &ContigCtx,
    p: &QcParams,
    low_abs: u32,
    out: &mut Vec<Region>,
) {
    let len = c.len as u64;
    if len == 0 {
        return;
    }
    let left_counts: Vec<u32> = c.clips.iter().map(|x| x[0]).collect();
    let right_counts: Vec<u32> = c.clips.iter().map(|x| x[1]).collect();
    let dp_pref = prefix_u32(&ctx.dp);
    let left_pref = prefix_u32(&left_counts);
    let right_pref = prefix_u32(&right_counts);
    let lefts = clip_candidates(&left_counts, &ctx.dp, len, p, true);
    let rights = clip_candidates(&right_counts, &ctx.dp, len, p, false);
    // Pair each left with the nearest unused right within `clip_pair_dist`.
    // Both sides are ascending by position and same-side clusters are separated
    // by more than `clip_merge`, so the binary-searched window holds at most
    // about 2·clip_pair_dist / clip_merge rights: O(lefts · (log rights +
    // clip_pair_dist / clip_merge)), i.e. a few comparisons per left at the
    // defaults (50 / 50), and still linear in lefts for wide windows.
    // Deterministic: lefts ascending; ties go to the smaller right index (=
    // smaller position, since rights are ascending).
    let dist = p.clip_pair_dist as u64;
    let mut used = vec![false; rights.len()];
    for l in &lefts {
        let lo = l.pos.saturating_sub(dist);
        let hi = l.pos + dist;
        let mut i = rights.partition_point(|r| r.pos < lo);
        let mut best: Option<usize> = None;
        let mut best_d = u64::MAX;
        while i < rights.len() && rights[i].pos <= hi {
            if !used[i] {
                let d = l.pos.abs_diff(rights[i].pos);
                if d < best_d {
                    best_d = d;
                    best = Some(i);
                }
            }
            i += 1;
        }
        match best {
            Some(i) => {
                used[i] = true;
                emit_clip(tid, c, ctx, p, low_abs, Some(l), Some(&rights[i]), &dp_pref, &left_pref, &right_pref, out);
            }
            None => emit_clip(tid, c, ctx, p, low_abs, Some(l), None, &dp_pref, &left_pref, &right_pref, out),
        }
    }
    for (i, r) in rights.iter().enumerate() {
        if !used[i] {
            emit_clip(tid, c, ctx, p, low_abs, None, Some(r), &dp_pref, &left_pref, &right_pref, out);
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::qc::test_util::*;

    #[test]
    fn clip_chaining_and_weighted_median() {
        // Gap 50 merges, gap 51 splits (single linkage).
        let mut counts = vec![0u32; 400];
        counts[100] = 1;
        counts[150] = 1;
        counts[201] = 1;
        let cl = clip_clusters(&counts, 50);
        assert_eq!(cl.len(), 2, "{cl:?}");
        assert_eq!(cl[0], ClipCluster { pos: 100, count: 2, start: 100, end: 150 });
        assert_eq!(cl[1], ClipCluster { pos: 201, count: 1, start: 201, end: 201 });
        // A gap of exactly 51 is not merged.
        let mut counts = vec![0u32; 400];
        counts[100] = 1;
        counts[151] = 1;
        assert_eq!(clip_clusters(&counts, 50).len(), 2);
        // Count-weighted median picks the heavy position.
        let mut counts = vec![0u32; 300];
        counts[100] = 1;
        counts[120] = 5;
        let cl = clip_clusters(&counts, 50);
        assert_eq!((cl.len(), cl[0].pos, cl[0].count, cl[0].start, cl[0].end), (1, 120, 6, 100, 120));
        let mut counts = vec![0u32; 300];
        counts[100] = 5;
        counts[120] = 1;
        assert_eq!(clip_clusters(&counts, 50)[0].pos, 100);
    }
    #[test]
    fn clip_aligned_window_at_contig_ends() {
        let dp = vec![10u32; 1_000];
        // Left clip: aligned side to the right, [pos, pos + 10).
        assert_eq!(clip_aligned_depth(&dp, 0, 1_000, true), Some(10.0));
        assert_eq!(clip_aligned_depth(&dp, 995, 1_000, true), Some(10.0)); // truncated to 5 bp
        assert_eq!(clip_aligned_depth(&dp, 1_000, 1_000, true), None); // empty window
        // Right clip: aligned side to the left, [pos - 9, pos + 1).
        assert_eq!(clip_aligned_depth(&dp, 0, 1_000, false), Some(10.0)); // truncated to [0, 1)
        assert_eq!(clip_aligned_depth(&dp, 999, 1_000, false), Some(10.0));
        // Values come from the window positions.
        let mut d = vec![0u32; 1_000];
        for x in &mut d[995..1_000] {
            *x = 4;
        }
        assert_eq!(clip_aligned_depth(&d, 995, 1_000, true), Some(4.0));
        assert_eq!(clip_aligned_depth(&d, 1_000, 1_000, true), None);
    }
    #[test]
    fn clip_tier_boundaries() {
        let p = QcParams::default();
        // Candidate (low-bar) gate: count and ratio floors.
        assert!(clip_candidate_ok(3, 0.15, &p));
        assert!(!clip_candidate_ok(3, 0.149, &p));
        assert!(!clip_candidate_ok(2, 0.9, &p));
        // Below the high bar: reported yellow, whatever the fraction.
        assert_eq!(clip_high_tier(99.9, 0.9, &p), "yellow");
        // At the high bar: orange up to f = 0.5, red above.
        assert_eq!(clip_high_tier(100.0, 0.5, &p), "orange");
        assert_eq!(clip_high_tier(100.0, 0.51, &p), "red");
        assert_eq!(clip_high_tier(1e6, 0.2, &p), "orange");
    }
    /// A two-sided `clip_both` hotspot requires the *combined* fraction: each
    /// side alone can be <= 0.5 while their sum is > 0.5 (each side holds about
    /// half the reads at a break). Exactly 0.5 stays orange.
    #[test]
    fn clip_combined_two_sided_fraction() {
        let p = QcParams { clip_high_s: 0.0, ..QcParams::default() };
        let emit = |lc: u64, rc: u64| -> Region {
            let len = 200_000usize;
            let (c, ctx) = mk_contig(len, vec![[0u32; 2]; len], vec![10u32; len]);
            let dp_pref = prefix_u32(&ctx.dp);
            let zero: Vec<u32> = vec![0; len];
            let lpref = prefix_u32(&zero);
            let rpref = prefix_u32(&zero);
            let mk = |pos: u64, count: u64| ClipSide {
                pos,
                count,
                ratio: count as f64 / 10.0,
                depth: 10.0,
                span: (pos, pos),
            };
            let l = mk(100_000, lc);
            let r = mk(100_010, rc);
            let mut out = Vec::new();
            emit_clip(0, &c, &ctx, &p, 10, Some(&l), Some(&r), &dp_pref, &lpref, &rpref, &mut out);
            out.pop().unwrap()
        };
        // Each side 0.4 <= 0.5, combined 0.8 > 0.5 -> red with the combined f.
        let reg = emit(4, 4);
        assert_eq!(reg.reason, "clip_both", "{reg:?}");
        assert_eq!(reg.tier, "red", "{reg:?}");
        assert!(reg.extra.contains("f=0.8000"), "{}", reg.extra);
        // Combined 0.5 -> orange (f is not > 0.5).
        let reg = emit(3, 2);
        assert_eq!(reg.tier, "orange", "{reg:?}");
        assert!(reg.extra.contains("f=0.5000"), "{}", reg.extra);
    }
    #[test]
    fn clip_pairing_and_lowcov() {
        let p = QcParams::default();
        // Left clip at 400, right at 449: 49 bp apart -> one paired region.
        let mut clips = vec![[0u32; 2]; 1_000];
        clips[400][0] = 5;
        clips[449][1] = 5;
        let (c, ctx) = mk_contig(1_000, clips, vec![10u32; 1_000]);
        let mut regs = Vec::new();
        clip_hotspot_regions(0, &c, &ctx, &p, 10, &mut regs);
        assert_eq!(regs.len(), 1, "{regs:?}");
        assert_eq!(regs[0].reason.as_str(), "clip_both");
        assert_eq!((regs[0].start, regs[0].end), (400, 450));
        // f = (5 + 5) / 10 = 1.0, but s is only ~3.7 (< clip_high_s 100): yellow.
        assert_eq!(regs[0].tier, "yellow");
        assert!(regs[0].extra.contains("left=5;right=5"), "{}", regs[0].extra);

        // 51 bp apart: not paired, two single regions.
        let mut clips = vec![[0u32; 2]; 1_000];
        clips[400][0] = 5;
        clips[451][1] = 5;
        let (c, ctx) = mk_contig(1_000, clips, vec![10u32; 1_000]);
        let mut regs = Vec::new();
        clip_hotspot_regions(0, &c, &ctx, &p, 10, &mut regs);
        let mut reasons: Vec<&str> = regs.iter().map(|r| r.reason.as_str()).collect();
        reasons.sort_unstable();
        assert_eq!(reasons, vec!["clip_left", "clip_right"], "{regs:?}");

        // Low aligned-side depth (< low_abs) sets lowcov=1 for the JSON alerts.
        let mut clips = vec![[0u32; 2]; 1_000];
        clips[400][0] = 5;
        clips[449][1] = 5;
        let (c, ctx) = mk_contig(1_000, clips, vec![5u32; 1_000]);
        let mut regs = Vec::new();
        clip_hotspot_regions(0, &c, &ctx, &p, 10, &mut regs);
        assert_eq!(regs.len(), 1, "{regs:?}");
        assert!(regs[0].lowcov, "{regs:?}");
        assert!(regs[0].extra.contains("lowcov=1"), "{}", regs[0].extra);
        // f = (5 + 5) / 5 = 2.0, but s ~ 3.7 < clip_high_s 100: yellow.
        assert_eq!(regs[0].tier, "yellow");

        // Clusters within clip_edge of a contig end are skipped.
        let mut clips = vec![[0u32; 2]; 1_000];
        clips[50][0] = 5;
        clips[950][1] = 5;
        let (c, ctx) = mk_contig(1_000, clips, vec![10u32; 1_000]);
        let mut regs = Vec::new();
        clip_hotspot_regions(0, &c, &ctx, &p, 10, &mut regs);
        assert!(regs.is_empty(), "{regs:?}");
    }
    /// 10,000 separate clusters: the pairing must not be a full cross product.
    #[test]
    fn clip_pairing_many_clusters_is_fast() {
        let n = 10_000usize;
        let spacing = 60usize;
        let off = 200usize;
        let len = off + spacing * n + 100;
        let mut clips = vec![[0u32; 2]; len];
        for i in 0..n {
            clips[off + spacing * i][0] = 5; // left clip
            clips[off + spacing * i + 30][1] = 5; // right clip, 30 bp away
        }
        let dp = vec![10u32; len];
        let (c, ctx) = mk_contig(len, clips, dp);
        let p = QcParams::default();
        let t = std::time::Instant::now();
        let mut regs = Vec::new();
        clip_hotspot_regions(0, &c, &ctx, &p, 10, &mut regs);
        let el = t.elapsed().as_secs_f64();
        assert_eq!(regs.len(), n, "{} != {n}", regs.len());
        assert!(regs.iter().all(|r| r.reason == "clip_both"), "unpaired clusters");
        assert!(el < 1.0, "pairing many clusters took {el:.3}s");
    }
    /// The reviewer repro: an 80 kb contig, 100M50S reads every 50 bp from 1000
    /// to 61 000, so one right-clip cluster spans ~60 kb — wider than the ±25 kb
    /// score window. The cluster must not underflow the window subtraction.
    #[test]
    fn clip_wide_cluster_score_is_finite() {
        let len = 80_000usize;
        let mut clips = vec![[0u32; 2]; len];
        for pos in (1_099..61_100).step_by(50) {
            clips[pos][1] += 1; // right clip at the last aligned base
        }
        let (c, ctx) = mk_contig(len, clips, vec![100u32; len]);
        let p = QcParams::default();
        let mut out = Vec::new();
        clip_hotspot_regions(0, &c, &ctx, &p, 10, &mut out);
        assert_eq!(out.len(), 1, "{out:?}");
        let s: f64 = out[0].extra.split(";s=").nth(1).unwrap().split(';').next().unwrap().parse().unwrap();
        assert!(s.is_finite() && s >= 0.0, "{:?}", out[0]);
    }
}
