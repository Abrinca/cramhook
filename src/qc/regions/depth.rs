//! `qc` depth-band regions.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use crate::qc::*;
use crate::scan::*;

/// Depth band assigned to one window (Plan-v3 "Depth" row).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum DepthBand {
    Low,
    Half,
    High,
}

impl DepthBand {
    pub(crate) fn label(self) -> &'static str {
        match self {
            DepthBand::Low => "depth_low",
            DepthBand::Half => "depth_half",
            DepthBand::High => "depth_high",
        }
    }
    pub(crate) fn name(self) -> &'static str {
        match self {
            DepthBand::Low => "low",
            DepthBand::Half => "half",
            DepthBand::High => "high",
        }
    }
}

/// Band of one window from its depth ratio and high-MAPQ share.
pub(crate) fn depth_band(ratio: f64, r: f64, p: &QcParams) -> Option<DepthBand> {
    if ratio <= p.depth_low {
        Some(DepthBand::Low)
    } else if ratio >= p.depth_half_min && ratio <= p.depth_half_max && r <= p.depth_half_max_mapq {
        Some(DepthBand::Half)
    } else if ratio >= p.depth_high && r >= p.depth_high_min_mapq {
        Some(DepthBand::High)
    } else {
        None
    }
}

/// One fixed depth window.
#[derive(Clone, Copy, Debug)]
pub(crate) struct DepthWindow {
    /// Window index on the contig.
    pub idx: usize,
    pub start: u64,
    pub end: u64,
    /// Local baseline: median of the window means within
    /// [`QcParams::depth_baseline_span`] (0 when every window is empty).
    pub baseline: f64,
    /// Mean DP ÷ local baseline (0 when the baseline is 0).
    pub ratio: f64,
    /// High-MAPQ depth ÷ DP; 1 in a window with no depth.
    pub r: f64,
    /// Assigned band, or `None` (between bands / a skipped contig end).
    pub band: Option<DepthBand>,
}

/// Median of the positive values (0 for an empty slice).
fn median_positive(v: &mut Vec<f64>) -> f64 {
    v.retain(|&x| x > 0.0);
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = v.len() / 2;
    if v.len() % 2 == 1 { v[m] } else { (v[m - 1] + v[m]) / 2.0 }
}

/// Per-window local baselines: the median of the window means within `span` bp
/// of the window centre (measured on `centers`, bp), ignoring windows with mean
/// 0. Circular contigs wrap across the origin; linear contigs truncate. A
/// contig shorter than `2 * span` uses one whole-contig median for every window
/// (the small-plasmid case).
pub(crate) fn local_baselines(
    means: &[f64],
    centers: &[f64],
    len: f64,
    span: f64,
    circular: bool,
) -> Vec<f64> {
    let n = means.len();
    if n == 0 {
        return Vec::new();
    }
    if len < 2.0 * span {
        let mut v = means.to_vec();
        return vec![median_positive(&mut v); n];
    }
    // Three copies on a circle, so each original window is reached from both
    // sides; the query uses the middle copy. `seen` avoids counting an original
    // window twice when `len == 2 * span` puts its two copies at ±span.
    let (ext_c, ext_m, n2, mid) = if circular {
        let mut c = Vec::with_capacity(3 * n);
        c.extend(centers.iter().map(|x| x - len));
        c.extend_from_slice(centers);
        c.extend(centers.iter().map(|x| x + len));
        let mut m = Vec::with_capacity(3 * n);
        m.extend_from_slice(means);
        m.extend_from_slice(means);
        m.extend_from_slice(means);
        (c, m, 3 * n, n)
    } else {
        (centers.to_vec(), means.to_vec(), n, 0)
    };
    let mut seen = vec![false; n];
    let mut out = vec![0.0f64; n];
    let (mut lo, mut hi) = (0usize, 0usize);
    for i in 0..n {
        let c = ext_c[mid + i];
        while lo < n2 && ext_c[lo] < c - span {
            lo += 1;
        }
        if hi < lo {
            hi = lo;
        }
        while hi < n2 && ext_c[hi] <= c + span {
            hi += 1;
        }
        let mut v = Vec::new();
        for k in lo..hi {
            let orig = k % n;
            if !seen[orig] && ext_m[k] > 0.0 {
                seen[orig] = true;
                v.push(ext_m[k]);
            }
        }
        for k in lo..hi {
            seen[k % n] = false;
        }
        out[i] = median_positive(&mut v);
    }
    out
}

/// Fixed, non-overlapping windows of `depth_window` bp over the per-position
/// `dp` and high-MAPQ `hiq` arrays. The last partial window is merged into the
/// previous one; contigs shorter than four windows give none. The ratio is the
/// window mean DP ÷ its local baseline (see [`local_baselines`]). On linear
/// contigs, windows within `reach` bp of either end get no band (their expected
/// depth tapers; Flagger §2.5). Circular contigs never skip.
pub(crate) fn depth_windows(
    dp: &[u32],
    hiq: &[u32],
    circular: bool,
    reach: f64,
    p: &QcParams,
) -> Vec<DepthWindow> {
    let w = p.depth_window as usize;
    let len = dp.len();
    if w == 0 || len / w < 4 {
        return Vec::new();
    }
    let nwin = len / w;
    let mut means = Vec::with_capacity(nwin);
    let mut centers = Vec::with_capacity(nwin);
    let mut rs = Vec::with_capacity(nwin);
    for i in 0..nwin {
        let start = i * w;
        let end = if i == nwin - 1 { len } else { (i + 1) * w };
        let (mut sd, mut sh) = (0u64, 0u64);
        for x in start..end {
            sd += dp[x] as u64;
            sh += hiq[x] as u64;
        }
        let n = (end - start) as f64;
        means.push(sd as f64 / n);
        centers.push((start as f64 + end as f64) / 2.0);
        rs.push(if sd == 0 { 1.0 } else { sh as f64 / sd as f64 });
    }
    let bl = local_baselines(&means, &centers, len as f64, p.depth_baseline_span as f64, circular);
    let mut out = Vec::with_capacity(nwin);
    for i in 0..nwin {
        let start = (i * w) as u64;
        let end = if i == nwin - 1 { len as u64 } else { ((i + 1) * w) as u64 };
        let baseline = bl[i];
        let ratio = if baseline > 0.0 { means[i] / baseline } else { 0.0 };
        let skipped = !circular && ((start as f64) < reach || ((len as u64 - end) as f64) < reach);
        let band = if skipped { None } else { depth_band(ratio, rs[i], p) };
        out.push(DepthWindow { idx: i, start, end, baseline, ratio, r: rs[i], band });
    }
    out
}

/// Merge same-band windows into `(start, end, band, span_windows, mean_ratio,
/// mean_r, mean_baseline)` runs. Consecutive same-band windows at most `gap + 1`
/// indices apart merge (a gap of `gap` windows of another band or of no band is
/// allowed); runs shorter than `min_len` are dropped.
#[allow(clippy::type_complexity)]
pub(crate) fn depth_region_runs(
    windows: &[DepthWindow],
    gap: usize,
    min_len: u64,
) -> Vec<(u64, u64, DepthBand, usize, f64, f64, f64)> {
    let mut out = Vec::new();
    for band in [DepthBand::Low, DepthBand::Half, DepthBand::High] {
        let ws: Vec<&DepthWindow> = windows.iter().filter(|w| w.band == Some(band)).collect();
        let mut i = 0;
        while i < ws.len() {
            let mut j = i;
            while j + 1 < ws.len() && ws[j + 1].idx - ws[j].idx <= gap + 1 {
                j += 1;
            }
            let (start, end) = (ws[i].start, ws[j].end);
            if end - start >= min_len {
                // Region means span the whole first..last window range, including
                // any allowed gap windows.
                let (lo_i, hi_i) = (ws[i].idx, ws[j].idx);
                let mut ratio = 0.0;
                let mut r = 0.0;
                let mut baseline = 0.0;
                let mut count = 0.0;
                for w in windows.iter().filter(|w| w.idx >= lo_i && w.idx <= hi_i) {
                    ratio += w.ratio;
                    r += w.r;
                    baseline += w.baseline;
                    count += 1.0;
                }
                let (mean_ratio, mean_r, mean_baseline) = (ratio / count, r / count, baseline / count);
                out.push((start, end, band, ws[j].idx - ws[i].idx + 1, mean_ratio, mean_r, mean_baseline));
            }
            i = j + 1;
        }
    }
    out
}

/// Tier of a depth-band run from its band and mean ratio: every run is at the
/// high bar, so tier by the disagreement fraction `f` (`low`/`half`: `1 - r`;
/// `high`: `(r - 1) / r`); red when `f > 0.5`.
fn band_tier(band: DepthBand, ratio: f64) -> &'static str {
    let f = match band {
        DepthBand::Low | DepthBand::Half => 1.0 - ratio,
        DepthBand::High => if ratio > 0.0 { (ratio - 1.0) / ratio } else { 0.0 },
    };
    if f > 0.5 { "red" } else { "orange" }
}

/// Depth-band BED regions (see [`depth_windows`] / [`depth_region_runs`]).
/// Every band run is "highly significant", so the tier is the disagreement
/// fraction `f`: `low`/`half` use `f = 1 - r`, `high` uses `f = (r - 1) / r`.
/// Red needs `f > 0.5` (low always; half when `r < 0.5`; high when `r > 2`),
/// orange otherwise. The old `depth_high` red upgrade from an overlapping
/// `snp_cluster`/`clip_hotspot` is gone: it could make a region red at
/// `f <= 0.5`, which the rule forbids. The contig baseline mode is reported in
/// `extra`; a contig with an unusable baseline (`depth_skipped`) is skipped.
pub(crate) fn depth_regions(
    scan: &ScanResult,
    ctxs: &[ContigCtx],
    p: &QcParams,
    mean_span: f64,
    out: &mut Vec<Region>,
) {
    for (tid, ctx) in ctxs.iter().enumerate() {
        if ctx.depth_skipped.is_some() || ctx.mode == 0 {
            continue;
        }
        let c = &scan.contigs[tid];
        let hiq: Vec<u32> = c.cov.iter().map(|x| x[HIQ]).collect();
        let wins = depth_windows(&ctx.dp, &hiq, ctx.circular, mean_span, p);
        for (start, end, band, span, ratio, r, baseline) in depth_region_runs(&wins, p.depth_gap_windows, p.depth_min_len) {
            let tier = band_tier(band, ratio);
            out.push(Region {
                tid,
                start,
                end,
                label: band.label(),
                score: 0,
                tier,
                reason: band.name().to_string(),
                depth_ratio: Some(ratio),
                hiq_share: Some(r),
                extra: format!("windows={span};baseline={baseline:.1}"),
                lowcov: false,
                link: None,
                clip_left: None,
                clip_right: None,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn depth_params(window: u32) -> QcParams {
        QcParams { depth_window: window, ..QcParams::default() }
    }
    fn dw(idx: usize, w: u64, band: Option<DepthBand>) -> DepthWindow {
        DepthWindow { idx, start: idx as u64 * w, end: (idx as u64 + 1) * w, baseline: 1.0, ratio: 1.0, r: 1.0, band }
    }
    #[test]
    fn depth_window_means_and_r() {
        let p = depth_params(10);
        // 4 windows of 10 bp; window 1 has no depth (r undefined -> 1).
        let mut dp = vec![40u32; 40];
        let mut hiq = vec![40u32; 40];
        for x in 0..10 {
            dp[10 + x] = 0;
            hiq[10 + x] = 0;
            dp[20 + x] = 20;
            hiq[20 + x] = 5;
            dp[30 + x] = 80;
            hiq[30 + x] = 80;
        }
        let ws = depth_windows(&dp, &hiq, false, 0.0, &p);
        assert_eq!(ws.len(), 4);
        assert_eq!((ws[0].ratio, ws[0].r), (1.0, 1.0));
        assert_eq!((ws[1].ratio, ws[1].r), (0.0, 1.0));
        assert_eq!((ws[2].ratio, ws[2].r), (0.5, 0.25));
        assert_eq!((ws[3].ratio, ws[3].r), (2.0, 1.0));
        // Whole-contig median (contig < 2 spans) is 40: 20, 40, 80.
        assert!(ws.iter().all(|w| w.baseline == 40.0));
    }
    #[test]
    fn depth_last_partial_window_merges() {
        let p = depth_params(10);
        // 45 bp -> 4 full windows + a 5 bp partial merged into the last.
        let dp = vec![40u32; 45];
        let hiq = vec![40u32; 45];
        let ws = depth_windows(&dp, &hiq, false, 0.0, &p);
        assert_eq!(ws.len(), 4);
        assert_eq!((ws[3].start, ws[3].end), (30, 45));
        // A contig shorter than four windows gives nothing.
        assert!(depth_windows(&vec![40u32; 39], &vec![40u32; 39], false, 0.0, &p).is_empty());
    }
    #[test]
    fn depth_band_boundaries() {
        let p = QcParams::default();
        assert_eq!(depth_band(0.25, 1.0, &p), Some(DepthBand::Low));
        assert_eq!(depth_band(0.26, 0.0, &p), None);
        assert_eq!(depth_band(0.35, 0.25, &p), Some(DepthBand::Half));
        assert_eq!(depth_band(0.65, 0.25, &p), Some(DepthBand::Half));
        assert_eq!(depth_band(0.65, 0.26, &p), None);
        assert_eq!(depth_band(0.34, 0.0, &p), None);
        assert_eq!(depth_band(1.75, 0.75, &p), Some(DepthBand::High));
        assert_eq!(depth_band(1.75, 0.74, &p), None);
        assert_eq!(depth_band(1.74, 1.0, &p), None);
        assert_eq!(depth_band(2.5, 0.9, &p), Some(DepthBand::High));
    }
    #[test]
    fn depth_region_merge_gap_and_min_len() {
        // Window 1000 bp, min_len 1000: a single window survives.
        let p = depth_params(1000);
        let wins = vec![dw(0, 1000, Some(DepthBand::Low)), dw(1, 1000, None), dw(2, 1000, Some(DepthBand::Low))];
        let r = depth_region_runs(&wins, p.depth_gap_windows, p.depth_min_len);
        assert_eq!(r.len(), 1, "{r:?}");
        assert_eq!((r[0].0, r[0].1, r[0].2, r[0].3), (0, 3000, DepthBand::Low, 3));
        // A gap of two windows (indices 3 apart) does not merge.
        let wins = vec![dw(0, 1000, Some(DepthBand::Low)), dw(3, 1000, Some(DepthBand::Low))];
        let r = depth_region_runs(&wins, p.depth_gap_windows, p.depth_min_len);
        assert_eq!(r.len(), 2, "{r:?}");
        // Below min_len: dropped (500 bp window, 1000 bp minimum).
        let p5 = depth_params(500);
        let wins = vec![dw(0, 500, Some(DepthBand::Low))];
        assert!(depth_region_runs(&wins, p5.depth_gap_windows, p5.depth_min_len).is_empty());
        // Exactly min_len is kept.
        let wins = vec![dw(0, 1000, Some(DepthBand::Low)), dw(1, 1000, Some(DepthBand::Low))];
        let r = depth_region_runs(&wins, p.depth_gap_windows, p.depth_min_len);
        assert_eq!(r.len(), 1, "{r:?}");
    }
    #[test]
    fn depth_region_means_include_gap_windows() {
        let p = depth_params(1000);
        let mut a = dw(0, 1000, Some(DepthBand::High));
        a.ratio = 2.0;
        a.r = 1.0;
        let mut b = dw(1, 1000, None);
        b.ratio = 1.0;
        b.r = 0.0;
        let mut c = dw(2, 1000, Some(DepthBand::High));
        c.ratio = 2.0;
        c.r = 1.0;
        let r = depth_region_runs(&[a, b, c], p.depth_gap_windows, p.depth_min_len);
        assert_eq!(r.len(), 1, "{r:?}");
        assert_eq!((r[0].0, r[0].1, r[0].2, r[0].3), (0, 3000, DepthBand::High, 3));
        assert!((r[0].4 - 5.0 / 3.0).abs() < 1e-9, "ratio {:?}", r[0]);
        assert!((r[0].5 - 2.0 / 3.0).abs() < 1e-9, "r {:?}", r[0]);
    }
    #[test]
    fn depth_linear_tapers_circular_does_not() {
        let p = depth_params(1000);
        // Window 0 is 0.25x the whole-contig median (40): Low.
        let mut dp = vec![40u32; 4000];
        let mut hiq = vec![40u32; 4000];
        dp[..1000].fill(10);
        hiq[..1000].fill(10);
        // Linear with reach 1500: window 0 is skipped (no band).
        let ws = depth_windows(&dp, &hiq, false, 1500.0, &p);
        assert!(ws[0].band.is_none(), "{ws:?}");
        assert!(depth_region_runs(&ws, p.depth_gap_windows, p.depth_min_len).is_empty());
        // Circular: nothing is skipped, window 0 is a Low region.
        let ws = depth_windows(&dp, &hiq, true, 1500.0, &p);
        assert_eq!(ws[0].band, Some(DepthBand::Low));
        let r = depth_region_runs(&ws, p.depth_gap_windows, p.depth_min_len);
        assert_eq!(r.len(), 1, "{r:?}");
        assert_eq!((r[0].0, r[0].1, r[0].2), (0, 1000, DepthBand::Low));
    }
    #[test]
    fn depth_low_baseline_skip() {
        let p = QcParams::default();
        let refseq = b"ACGT".repeat(25); // 100 bp
        let make = |dp: u32| {
            let mut c = ContigCounts {
                name: "chr".to_string(),
                len: 100,
                bases: vec![[0u32; 8]; 100],
                cov: vec![[0u32; crate::scan::NCOV]; 100],
                clips: vec![[0u32; 2]; 100],
                edge_clips: 0,
                aligned_bases: 0,
            };
            for x in c.cov.iter_mut() {
                x[crate::scan::DPF] = dp;
            }
            c
        };
        // Mode 5 < low_abs 10: unusable baseline.
        let ctx = contig_ctx(&make(5), &refseq, false, &p, 10, 0.25);
        assert_eq!(ctx.depth_skipped, Some("low_coverage"));
        // Mode 40 >= 10: usable.
        let ctx = contig_ctx(&make(40), &refseq, false, &p, 10, 0.25);
        assert_eq!(ctx.depth_skipped, None);
    }
    #[test]
    fn depth_band_tier_from_fraction() {
        // Low band (r <= 0.25): f = 1 - r >= 0.75 -> always red.
        assert_eq!(band_tier(DepthBand::Low, 0.25), "red");
        assert_eq!(band_tier(DepthBand::Low, 0.05), "red");
        // Half band: red only when r < 0.5 (f > 0.5).
        assert_eq!(band_tier(DepthBand::Half, 0.35), "red");
        assert_eq!(band_tier(DepthBand::Half, 0.49), "red");
        assert_eq!(band_tier(DepthBand::Half, 0.50), "orange");
        assert_eq!(band_tier(DepthBand::Half, 0.65), "orange");
        // High band: red only when r > 2 (f > 0.5); the old 1.75 upgrade is gone.
        assert_eq!(band_tier(DepthBand::High, 1.75), "orange");
        assert_eq!(band_tier(DepthBand::High, 2.00), "orange");
        assert_eq!(band_tier(DepthBand::High, 2.01), "red");
        assert_eq!(band_tier(DepthBand::High, 3.0), "red");
    }
    #[test]
    fn local_baseline_wraps_and_short_contig() {
        // Reviewer repro: 4 windows, centres 500/1500/2500/3500, len 4000.
        let means = vec![10.0, 100.0, 100.0, 100.0];
        let centers = vec![500.0, 1500.0, 2500.0, 3500.0];
        // Linear: window 0 sees only 0 and 1 -> median(10,100) = 55; window 3
        // sees only 2 and 3 -> 100.
        let lin = local_baselines(&means, &centers, 4000.0, 1500.0, false);
        assert_eq!(lin[0], 55.0, "{lin:?}");
        assert_eq!(lin[3], 100.0, "{lin:?}");
        // Circular: window 0 sees 1 and the wrapped window 3 (at −len) -> 100;
        // window 3 sees 2, 3 and the wrapped window 0 (at +len) -> 100. Both
        // wrap directions are covered.
        let cir = local_baselines(&means, &centers, 4000.0, 1500.0, true);
        assert_eq!(cir[0], 100.0, "{cir:?}");
        assert_eq!(cir[3], 100.0, "{cir:?}");
        // Zero-mean windows are excluded from the median.
        let means = vec![0.0, 40.0, 40.0, 100.0];
        let b = local_baselines(&means, &centers, 4000.0, 1500.0, false);
        assert_eq!(b[3], 70.0, "{b:?}");
        // Contig shorter than 2 spans: one whole-contig median everywhere
        // (means 40, 40, 100 with the zero excluded).
        let b = local_baselines(&means, &centers, 2000.0, 1500.0, false);
        assert!(b.iter().all(|&x| x == 40.0), "{b:?}");
        // len == 2 * span: the two copies of a window at exactly ±span are
        // counted once. Every window's range covers all four here.
        let means = vec![10.0, 20.0, 30.0, 1000.0];
        let cir = local_baselines(&means, &centers, 4000.0, 2000.0, true);
        assert!(cir.iter().all(|&x| x == 25.0), "{cir:?}");
    }
    #[test]
    fn local_baseline_gradient_and_step() {
        // 150 kb in 300 windows: a smooth 40x -> 339x gradient.
        let p = QcParams { depth_window: 500, depth_baseline_span: 50_000, ..QcParams::default() };
        let nwin = 300usize;
        let len = nwin * 500;
        let mut dp = vec![0u32; len];
        let mut hiq = vec![0u32; len];
        for i in 0..nwin {
            let v = (40 + i) as u32;
            for x in i * 500..(i + 1) * 500 {
                dp[x] = v;
                hiq[x] = v;
            }
        }
        // A smooth gradient gives no depth regions with the local baseline.
        let ws = depth_windows(&dp, &hiq, false, 0.0, &p);
        assert!(depth_region_runs(&ws, p.depth_gap_windows, p.depth_min_len).is_empty(), "{ws:?}");
        // A 3 kb step (windows 140..146 at 2x) is still flagged.
        for x in 140 * 500..146 * 500 {
            dp[x] *= 2;
            hiq[x] *= 2;
        }
        let ws = depth_windows(&dp, &hiq, false, 0.0, &p);
        let r = depth_region_runs(&ws, p.depth_gap_windows, p.depth_min_len);
        assert_eq!(r.len(), 1, "{r:?}");
        assert_eq!((r[0].0, r[0].1, r[0].2), (70_000, 73_000, DepthBand::High));
    }
    /// Timing only: 10k windows (≈5 Mb contig at window 500), circular. Run
    /// with `cargo test --release local_baseline_cost -- --ignored --nocapture`.
    #[test]
    #[ignore = "timing measurement"]
    fn local_baseline_cost() {
        let n = 10_000usize;
        let means: Vec<f64> = (0..n).map(|i| 40.0 + (i % 37) as f64).collect();
        let centers: Vec<f64> = (0..n).map(|i| i as f64 * 500.0 + 250.0).collect();
        let len = n as f64 * 500.0;
        let t = std::time::Instant::now();
        let b = local_baselines(&means, &centers, len, 50_000.0, true);
        let el = t.elapsed().as_secs_f64();
        assert_eq!(b.len(), n);
        crate::note!("local_baselines n={n}: {:.1} ms", el * 1000.0);
    }
}
