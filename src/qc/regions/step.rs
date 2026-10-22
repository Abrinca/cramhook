//! `qc` depth-step and depth-segment regions.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use crate::qc::*;

// Depth-step detector (qc spec section 2).
//
// Contrast `Δ_i = mean DP[i, i+w) - mean DP[i-w, i)` from prefix sums, and its
// local score over ±`span` bp on a coarse grid, interpolated:
//
//   z_i = (Δ_i - median) / (1.4826 * max(MAD, sqrt(local_depth/w)))
//
// i.e. the MAD carries the 1.4826 normal-consistency factor and the Poisson
// floor `sqrt(local_depth/w)` is inside that max (the 8/20 bars were calibrated
// with this exact form). Peaks are the local maxima of `|z|` on a `stride` grid,
// non-maximum suppressed ±`nms` bp and refined at 1 bp. Circular contigs wrap
// for the contrast, the neighbourhood and the pairing; linear ones do not (a
// step at a real linear end is real).

/// Stride (bp) of the contrast positions.
const STEP_STRIDE: i64 = 10;
/// Coarse-grid spacing (bp) for the local median/MAD.
const STEP_GRID: i64 = 1000;
/// ±bp around a coarse peak for the 1 bp refinement.
const STEP_REFINE: i64 = 50;
/// Half-width (bp) of a reported `depth_step` interval.
const STEP_REGION_HALF: i64 = 50;
/// ±bp around a `depth_segment` end for the inside/flank ratio.
const STEP_SEG_FLANK: i64 = 1000;

/// A refined depth-step peak (position in `[0, n)`, z at 1 bp, signed contrast).
#[derive(Clone, Copy, Debug)]
struct StepPeak {
    pos: i64,
    z: f64,
    delta: f64,
}

/// Circular-aware prefix sum: `cum(i)` over the contig repeated periodically
/// (circular) or the ordinary prefix sum of the `[0, n]` range (linear).
fn step_cum(cs: &[i64], n: i64, i: i64, circular: bool) -> i64 {
    if circular {
        let q = i.div_euclid(n);
        let r = i.rem_euclid(n);
        q * cs[n as usize] + cs[r as usize]
    } else {
        cs[i as usize]
    }
}

/// `Δ_i = mean DP[i, i+w) - mean DP[i-w, i)` at position `p` (`p - w >= 0` and
/// `p + w <= n` on linear contigs).
fn step_delta(cs: &[i64], n: i64, p: i64, w: i64, circular: bool) -> f64 {
    let left = step_cum(cs, n, p, circular) - step_cum(cs, n, p - w, circular);
    let right = step_cum(cs, n, p + w, circular) - step_cum(cs, n, p, circular);
    (right - left) as f64 / w as f64
}

/// Mean depth over `len` positions starting at `start` (wrapping modulo `n`).
fn step_mean_len(dp: &[u32], start: i64, len: i64, n: i64) -> f64 {
    if len <= 0 {
        return 0.0;
    }
    let mut total = 0i64;
    let mut rem = len;
    let mut i = start.rem_euclid(n);
    while rem > 0 {
        let take = rem.min(n - i);
        for k in i..i + take {
            total += dp[k as usize] as i64;
        }
        rem -= take;
        i = 0;
    }
    total as f64 / len as f64
}

/// Mean depth over `[a, b)`: wraps on circular contigs, clamps on linear ones.
fn step_mean_range(dp: &[u32], n: i64, a: i64, b: i64, circular: bool) -> f64 {
    if b <= a {
        return 0.0;
    }
    if circular {
        step_mean_len(dp, a, b - a, n)
    } else {
        let a = a.max(0);
        let b = b.min(n);
        if b <= a {
            return 0.0;
        }
        let s: i64 = (a..b).map(|i| dp[i as usize] as i64).sum();
        s as f64 / (b - a) as f64
    }
}

/// `f = |Δ| / max(depth left, depth right)` and the two side means.
fn step_frac(dp: &[u32], n: i64, pos: i64, delta: f64, w: i64, circular: bool) -> (f64, f64, f64) {
    let left = step_mean_range(dp, n, pos - w, pos, circular);
    let right = step_mean_range(dp, n, pos, pos + w, circular);
    let maxd = left.max(right);
    let frac = if maxd > 0.0 { delta.abs() / maxd } else { 0.0 };
    (frac, left, right)
}

/// Tier of a depth step: `s = |z|`, `f = |Δ| / max(depth sides)`. Red requires
/// the high bar and `f > 0.5`; orange is the high bar with a small fraction;
/// yellow is the low bar.
fn step_tier(z_abs: f64, frac: f64, p: &QcParams) -> Option<&'static str> {
    if z_abs >= p.step_red {
        Some(if frac > 0.5 { "red" } else { "orange" })
    } else if z_abs >= p.step_yellow {
        Some("yellow")
    } else {
        None
    }
}

/// BED pieces of `[pos - half, pos + half)`, split at the origin on a circular
/// contig, clamped on a linear one.
fn step_interval(pos: i64, half: i64, n: i64, circular: bool) -> Vec<(u64, u64)> {
    if !circular {
        let a = (pos - half).max(0);
        let b = (pos + half).min(n);
        return if b > a { vec![(a as u64, b as u64)] } else { Vec::new() };
    }
    if pos - half < 0 {
        vec![((pos - half + n) as u64, n as u64), (0, (pos + half) as u64)]
    } else if pos + half > n {
        vec![((pos - half) as u64, n as u64), (0, (pos + half - n) as u64)]
    } else {
        vec![((pos - half) as u64, (pos + half) as u64)]
    }
}

/// Whether a `depth_step` peak at `pos` (0 <= pos < n) is *at the origin* of a
/// wrapping contig: its reported 100 bp interval `[pos - half, pos + half)`
/// touches the origin, i.e. a BED piece has `start == 0` or `end == len`. This
/// is the report's origin-exact definition; the literal "pos == 0 (mod len)"
/// misses peaks that refine a few bp off the origin but still split there.
fn step_at_origin(pos: i64, n: i64) -> bool {
    if n <= 0 {
        return false;
    }
    let pos = pos.rem_euclid(n);
    step_interval(pos, STEP_REGION_HALF, n, true).iter().any(|&(s, e)| s == 0 || e == n as u64)
}

/// Median of a mutable slice (reorders it). `n == 0` gives 0.
fn step_median(v: &mut [f64]) -> f64 {
    let n = v.len();
    if n == 0 {
        return 0.0;
    }
    let mid = n / 2;
    v.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if n % 2 == 1 {
        v[mid]
    } else {
        let lo = v[..mid].iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (lo + v[mid]) / 2.0
    }
}

/// Linear interpolation of `(xp, fp)` at `x` (clamped at the ends).
fn step_interp(x: f64, xp: &[i64], fp: &[f64]) -> f64 {
    let n = xp.len();
    if n == 0 {
        return 0.0;
    }
    if x <= xp[0] as f64 {
        return fp[0];
    }
    if x >= xp[n - 1] as f64 {
        return fp[n - 1];
    }
    let i = xp.partition_point(|&p| (p as f64) <= x) - 1;
    let x0 = xp[i] as f64;
    let x1 = xp[i + 1] as f64;
    let t = (x - x0) / (x1 - x0);
    fp[i] + t * (fp[i + 1] - fp[i])
}

/// Push `delta[k]` for every `pos[k]` in `[lo_p, hi_p]`. `pos` is the actual
/// (sorted) contrast position of each `delta` entry, so linear contigs — whose
/// grid starts at `w`, not at 0 — gather the right neighbourhood.
fn step_gather(pos: &[i64], delta: &[f64], lo_p: i64, hi_p: i64, vals: &mut Vec<f64>) {
    if hi_p < lo_p {
        return;
    }
    let k0 = pos.partition_point(|&p| p < lo_p);
    let k1 = pos.partition_point(|&p| p <= hi_p);
    vals.extend_from_slice(&delta[k0..k1]);
}

/// Set `s[k] = true` for every index whose `pos[k]` lies in `[lo, hi]`.
fn step_mark(pos: &[i64], s: &mut [bool], lo: i64, hi: i64) {
    let a = pos.partition_point(|&p| p < lo);
    let b = pos.partition_point(|&p| p <= hi);
    s[a..b].fill(true);
}

/// Exact peak position at stride 1 within ±[`STEP_REFINE`] of `p0`.
#[allow(clippy::too_many_arguments)]
fn step_refine(
    cs: &[i64],
    n: i64,
    p0: i64,
    w: i64,
    gpos: &[i64],
    med_g: &[f64],
    mad_g: &[f64],
    floor_g: &[f64],
    circular: bool,
) -> Option<StepPeak> {
    let ps: Vec<i64> = if circular {
        let mut v: Vec<i64> = (p0 - STEP_REFINE..=p0 + STEP_REFINE).map(|x| x.rem_euclid(n)).collect();
        v.sort_unstable();
        v.dedup();
        v
    } else {
        let lo = w.max(p0 - STEP_REFINE);
        let hi = (n - w).min(p0 + STEP_REFINE);
        if hi <= lo {
            return None;
        }
        (lo..hi).collect()
    };
    let mut best = f64::NEG_INFINITY;
    let mut out = StepPeak { pos: p0, z: 0.0, delta: 0.0 };
    for &q in &ps {
        let d = step_delta(cs, n, q, w, circular);
        let med = step_interp(q as f64, gpos, med_g);
        let mut mad = step_interp(q as f64, gpos, mad_g);
        let flr = step_interp(q as f64, gpos, floor_g);
        mad = mad.max(flr.max(1e-9));
        let z = (d - med) / (1.4826 * mad);
        if z.abs() > best {
            best = z.abs();
            out = StepPeak { pos: q, z, delta: d };
        }
    }
    Some(out)
}

/// Rank of a tier name (red > orange > yellow), for picking the lower of two.
fn step_tier_rank(t: &str) -> u8 {
    match t {
        "red" => 3,
        "orange" => 2,
        _ => 1,
    }
}

/// One maximal same-sign run of position-sorted steps: `left` is its leftmost
/// position (used only to order runs) and `strong` its strongest step (max
/// `|z|`, ties to the smaller position).
#[derive(Clone, Copy, Debug)]
struct StepRun {
    left: i64,
    strong: StepPeak,
}

/// Group the position-sorted steps into maximal same-sign runs, ordered by
/// their leftmost position. On a circular contig the list is rotated to start
/// at a sign change, so a run spanning the origin is returned as one run; a
/// linear contig is traversed from index 0 without rotation.
fn step_runs(peaks: &[StepPeak], circular: bool) -> Vec<StepRun> {
    let m = peaks.len();
    if m == 0 {
        return Vec::new();
    }
    let seq: Vec<StepPeak> = if circular {
        match (0..m).find(|&i| (peaks[i].delta > 0.0) != (peaks[(i + m - 1) % m].delta > 0.0)) {
            Some(s) => (0..m).map(|k| peaks[(s + k) % m]).collect(),
            None => peaks.to_vec(),
        }
    } else {
        peaks.to_vec()
    };
    let mut runs: Vec<StepRun> = Vec::new();
    for p in seq {
        if let Some(r) = runs.last_mut()
            && (p.delta > 0.0) == (r.strong.delta > 0.0)
        {
            r.left = r.left.min(p.pos);
            if p.z.abs() > r.strong.z.abs()
                || (p.z.abs() == r.strong.z.abs() && p.pos < r.strong.pos)
            {
                r.strong = p;
            }
            continue;
        }
        runs.push(StepRun { left: p.pos, strong: p });
    }
    // Order runs by their leftmost step.
    runs.sort_by_key(|r| r.left);
    runs
}

/// Pair runs into `depth_segment` intervals (output only).
///
/// The reported `depth_step` peaks of a contig are sorted by position and
/// grouped into maximal same-sign runs; each run is represented by its
/// **strongest** step (max `|z|`, ties to the smaller position). Runs alternate
/// in sign, so a segment always links an up run and a down run. There are two
/// ways to pair neighbouring runs: starting at run 0 (`(0,1), (2,3), ...`) or
/// starting at run 1 (`(1,2), (3,4), ...`); on a circular contig the last and
/// first runs are neighbours (a matching may pair `(m-1, 0)` across the
/// origin), on a linear contig they are not. The matching that pairs more runs
/// wins (on a linear contig the start-1 matching drops the two ends), then the
/// smaller total forward segment length (`n` = contig length, so a wrapping arc
/// is measured the longer way); ties go to the run-0 matching. An odd number of
/// runs on a circle leaves one run unpaired in each matching, so those are
/// compared on length. A **gain**
/// segment runs from an up run to the following down run, a **loss** (dip) from
/// a down run to the following up run, and a segment whose measured
/// inside/flank ratio contradicts its direction is dropped in
/// [`emit_step_segment`].
///
/// Returns `(first, second)` in genome order (the first is the segment start).
fn step_pair_segments(peaks: &[StepPeak], circular: bool, n: i64) -> Vec<(StepPeak, StepPeak)> {
    let runs = step_runs(peaks, circular);
    let m = runs.len();
    if m < 2 {
        return Vec::new();
    }
    // One alternating matching, starting at run `start`.
    let matching = |start: usize| -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        if circular {
            for t in 0..m / 2 {
                let a = (start + 2 * t) % m;
                let b = (a + 1) % m;
                if (runs[a].strong.delta > 0.0) != (runs[b].strong.delta > 0.0) {
                    out.push((a, b));
                }
            }
        } else {
            let mut i = start;
            while i + 1 < m {
                if (runs[i].strong.delta > 0.0) != (runs[i + 1].strong.delta > 0.0) {
                    out.push((i, i + 1));
                }
                i += 2;
            }
        }
        out
    };
    // Total forward segment length: `n + (e - s)` for an arc that wraps.
    let cost = |pairs: &[(usize, usize)]| -> i64 {
        pairs
            .iter()
            .map(|&(i, j)| {
                let (s, e) = (runs[i].strong.pos, runs[j].strong.pos);
                if e >= s {
                    e - s
                } else if circular {
                    e - s + n
                } else {
                    s - e
                }
            })
            .sum()
    };
    let a = matching(0);
    let b = matching(1);
    // Prefer the matching that pairs more runs (on a linear contig the start-1
    // matching drops the two ends); otherwise the smaller total length; ties go
    // to the run-0 matching.
    let pairs = if a.len() != b.len() {
        if a.len() > b.len() { a } else { b }
    } else if cost(&a) <= cost(&b) {
        a
    } else {
        b
    };
    pairs.into_iter().map(|(i, j)| (runs[i].strong, runs[j].strong)).collect()
}

/// The `depth_segment` inside/flank mean ratio (`None` when the flank mean is 0).
fn step_seg_ratio(dp: &[u32], n: i64, start: i64, end: i64, flank: i64, circular: bool) -> Option<f64> {
    let interior = if end > start {
        step_mean_len(dp, start, end - start, n)
    } else if end == start {
        dp[start as usize] as f64
    } else {
        step_mean_len(dp, start, end - start + n, n)
    };
    let (mut sum, mut cnt) = (0f64, 0i64);
    if circular {
        sum += step_mean_len(dp, start - flank, flank, n) * flank as f64;
        sum += step_mean_len(dp, end, flank, n) * flank as f64;
        cnt = 2 * flank;
    } else {
        for (a, b) in [(start - flank, start), (end, end + flank)] {
            let a = a.max(0);
            let b = b.min(n);
            if b > a {
                let s: i64 = (a..b).map(|i| dp[i as usize] as i64).sum();
                sum += s as f64;
                cnt += b - a;
            }
        }
    }
    if cnt == 0 {
        return None;
    }
    let fm = sum / cnt as f64;
    if fm > 0.0 { Some(interior / fm) } else { None }
}

/// Emit one `depth_segment` region for a paired interval. It is written to the
/// BED (tier = the weaker end) but is never a finding member.
#[allow(clippy::too_many_arguments)]
fn emit_step_segment(
    tid: usize,
    ctx: &ContigCtx,
    n: i64,
    a: &StepPeak,
    b: &StepPeak,
    w: i64,
    low_abs: u32,
    p: &QcParams,
    out: &mut Vec<Region>,
) {
    let dp = &ctx.dp;
    let circular = ctx.circular;
    let start = a.pos.rem_euclid(n);
    let end = b.pos.rem_euclid(n);
    let kind = if a.delta > 0.0 { "gain" } else { "loss" };
    let (fa, la, ra) = step_frac(dp, n, a.pos, a.delta, w, circular);
    let (fb, lb, rb) = step_frac(dp, n, b.pos, b.delta, w, circular);
    // The segment tier is the lower of the two endpoint tiers.
    let ta = match step_tier(a.z.abs(), fa, p) {
        Some(t) => t,
        None => return,
    };
    let tb = match step_tier(b.z.abs(), fb, p) {
        Some(t) => t,
        None => return,
    };
    let tier = if step_tier_rank(ta) <= step_tier_rank(tb) { ta } else { tb };
    let lowcov = [la, ra, lb, rb].iter().any(|&v| v < low_abs as f64);
    let ratio = step_seg_ratio(dp, n, start, end, STEP_SEG_FLANK, circular);
    // Drop a segment whose measured ratio contradicts its direction.
    if let Some(r) = ratio {
        let consistent = if a.delta > 0.0 { r > 1.0 } else { r < 1.0 };
        if !consistent {
            return;
        }
    }
    let extra = match ratio {
        Some(r) => format!("ratio={r:.3}"),
        None => "ratio=.".to_string(),
    };
    let mut pieces: Vec<(u64, u64)> = Vec::new();
    if end > start {
        pieces.push((start as u64, end as u64));
    } else if end < start {
        pieces.push((start as u64, n as u64));
        if end > 0 {
            pieces.push((0, end as u64));
        }
    }
    for (s, e) in pieces {
        out.push(Region {
            tid,
            start: s,
            end: e,
            label: "depth_segment",
            score: 0,
            tier,
            reason: kind.to_string(),
            depth_ratio: None,
            hiq_share: None,
            extra: extra.clone(),
            lowcov,
            link: Some((start as u64, end as u64)),
            clip_left: None,
            clip_right: None,
        });
    }
}

/// `depth_step` / `depth_segment` regions for every usable contig.
///
/// Returns the per-contig count of `depth_step` peaks dropped as origin
/// artefacts (indexed by contig, 0 for contigs without one).
pub(crate) fn depth_step_regions(ctxs: &[ContigCtx], p: &QcParams, low_abs: u32, out: &mut Vec<Region>) -> Vec<u64> {
    let w = p.step_window as i64;
    let stride = STEP_STRIDE;
    let grid = STEP_GRID;
    let half = p.step_span as i64;
    let nms = p.step_nms as i64;
    let mut origin_dropped = vec![0u64; ctxs.len()];
    if w <= 0 || stride <= 0 || half < 0 {
        return origin_dropped;
    }
    for (tid, ctx) in ctxs.iter().enumerate() {
        // Only an empty / all-zero contig is skipped: depth steps are computed
        // on low-coverage contigs too (their regions are flagged `lowcov`).
        if ctx.mode == 0 {
            continue;
        }
        let dp = &ctx.dp;
        let n = dp.len() as i64;
        if n < 2 * w + stride {
            continue;
        }
        let circular = ctx.circular;
        // Prefix sums (ordinary, length n + 1).
        let mut cs: Vec<i64> = Vec::with_capacity(n as usize + 1);
        cs.push(0);
        let mut acc = 0i64;
        for &d in dp {
            acc += d as i64;
            cs.push(acc);
        }
        // Contrast positions: the whole contig (circular) or the interior
        // (linear, ends excluded so pos ± w stays in range).
        let mut pos: Vec<i64> = Vec::new();
        if circular {
            let mut k = 0i64;
            while k < n {
                pos.push(k);
                k += stride;
            }
        } else {
            let mut k = w;
            while k < n - w {
                pos.push(k);
                k += stride;
            }
        }
        let m = pos.len();
        if m == 0 {
            continue;
        }
        let mut delta = vec![0f64; m];
        for (i, &pp) in pos.iter().enumerate() {
            delta[i] = step_delta(&cs, n, pp, w, circular);
        }
        // Coarse grid for the local median/MAD (and its Poisson floor).
        let mut gpos: Vec<i64> = Vec::new();
        let mut g = 0i64;
        while g < n {
            gpos.push(g);
            g += grid;
        }
        let ng = gpos.len();
        let mut med_g = vec![0f64; ng];
        let mut mad_g = vec![0f64; ng];
        let mut floor_g = vec![0f64; ng];
        // A small circular contig uses one whole-contig median/MAD.
        let global = circular && n <= 2 * half;
        let mut vals: Vec<f64> = Vec::new();
        for gi in 0..ng {
            let gg = gpos[gi];
            if global {
                if gi == 0 {
                    let mut v = delta.clone();
                    let med = step_median(&mut v);
                    for x in v.iter_mut() {
                        *x = (*x - med).abs();
                    }
                    med_g[0] = med;
                    mad_g[0] = step_median(&mut v);
                } else {
                    med_g[gi] = med_g[0];
                    mad_g[gi] = mad_g[0];
                }
            } else {
                vals.clear();
                step_gather(&pos, &delta, gg - half, gg + half, &mut vals);
                if circular {
                    if gg - half < 0 {
                        step_gather(&pos, &delta, gg - half + n, n - 1, &mut vals);
                    }
                    if gg + half >= n {
                        step_gather(&pos, &delta, 0, gg + half - n, &mut vals);
                    }
                }
                if vals.is_empty() {
                    med_g[gi] = 0.0;
                    mad_g[gi] = 0.0;
                } else {
                    let med = step_median(&mut vals);
                    for x in vals.iter_mut() {
                        *x = (*x - med).abs();
                    }
                    med_g[gi] = med;
                    mad_g[gi] = step_median(&mut vals);
                }
            }
            let d = if circular {
                (step_cum(&cs, n, gg + w, true) - step_cum(&cs, n, gg - w, true)) as f64 / (2 * w) as f64
            } else {
                let a = (gg - w).max(0);
                let b = (gg + w).min(n);
                if b <= a { 0.0 } else { (cs[b as usize] - cs[a as usize]) as f64 / (b - a) as f64 }
            };
            floor_g[gi] = (d.max(0.0) / w as f64).sqrt();
        }
        // z at every stride position; the MAD is floored at the Poisson
        // expectation `sqrt(local_depth / w)` whenever it binds.
        let mut z = vec![0f64; m];
        for i in 0..m {
            let x = pos[i] as f64;
            let med = step_interp(x, &gpos, &med_g);
            let mut mad = step_interp(x, &gpos, &mad_g);
            let flr = step_interp(x, &gpos, &floor_g);
            mad = mad.max(flr.max(1e-9));
            z[i] = (delta[i] - med) / (1.4826 * mad);
        }
        // Local maxima of |z| (NMS ±nms), taken strongest-first.
        let zmin = p.step_yellow;
        let mut cand: Vec<usize> = (0..m).filter(|&i| z[i].abs() >= zmin).collect();
        cand.sort_by(|&a, &b| {
            z[b].abs()
                .partial_cmp(&z[a].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.cmp(&b))
        });
        // Suppression bitmap (strongest-first): a step is dropped when it falls
        // within ±nms bp of an already-taken stronger step.
        let mut suppressed = vec![false; m];
        let mut taken: Vec<usize> = Vec::new();
        for &i in &cand {
            if suppressed[i] {
                continue;
            }
            taken.push(i);
            let (lo, hi) = (pos[i] - nms, pos[i] + nms);
            step_mark(&pos, &mut suppressed, lo, hi);
            if circular {
                step_mark(&pos, &mut suppressed, lo + n, hi + n);
                step_mark(&pos, &mut suppressed, lo - n, hi - n);
            }
        }
        taken.sort_by_key(|&i| pos[i]);
        // Refine to 1 bp.
        let mut peaks: Vec<StepPeak> = Vec::new();
        for &i in &taken {
            if let Some(pk) = step_refine(&cs, n, pos[i], w, &gpos, &med_g, &mad_g, &floor_g, circular) {
                peaks.push(pk);
            }
        }
        // Drop a `depth_step` whose 100 bp interval touches the origin of a
        // wrapping contig. A read crossing the arbitrarily chosen origin is
        // split into a primary and a supplementary alignment (both soft-clipped
        // at the ends), so the read counts just inside each end differ (median
        // 1.34x); the wrapped step contrast therefore peaks at the origin even
        // though the wrap arithmetic is correct. It is an aligner/representation
        // artefact, not a coverage change. The filter runs *after* non-maximum
        // suppression (above), so a dropped peak cannot unmask a neighbour it
        // suppressed and every other call stays byte-identical; `depth_segment`
        // pairing then uses the remaining steps. See docs/project/decisions.md "Depth
        // steps at the origin".
        if circular {
            let mut kept = Vec::with_capacity(peaks.len());
            for pk in peaks {
                if step_at_origin(pk.pos, n) {
                    origin_dropped[tid] += 1;
                } else {
                    kept.push(pk);
                }
            }
            peaks = kept;
        }
        for pk in &peaks {
            let (frac, left, right) = step_frac(dp, n, pk.pos, pk.delta, w, circular);
            let tier = match step_tier(pk.z.abs(), frac, p) {
                Some(t) => t,
                None => continue,
            };
            let lowcov = left < low_abs as f64 || right < low_abs as f64;
            let reason = if pk.delta > 0.0 { "step_up" } else { "step_down" };
            let mean = (left + right) / 2.0;
            let depth_ratio = if ctx.mode > 0 { Some(mean / ctx.mode as f64) } else { None };
            let mut extra = format!(
                "z={:.2};delta={:.2};depth_left={:.1};depth_right={:.1};frac={:.3}",
                pk.z, pk.delta, left, right, frac
            );
            if lowcov {
                extra.push_str(";lowcov=1");
            }
            for (s, e) in step_interval(pk.pos, STEP_REGION_HALF, n, circular) {
                out.push(Region {
                    tid,
                    start: s,
                    end: e,
                    label: "depth_step",
                    score: 0,
                    tier,
                    reason: reason.to_string(),
                    depth_ratio,
                    hiq_share: None,
                    extra: extra.clone(),
                    lowcov,
                    link: None,
                    clip_left: None,
                    clip_right: None,
                });
            }
        }
        for (a, b) in step_pair_segments(&peaks, circular, n) {
            emit_step_segment(tid, ctx, n, &a, &b, w, low_abs, p, out);
        }
    }
    origin_dropped
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::scan::*;

    fn step_ctx(dp: Vec<u32>, circular: bool) -> ContigCtx {
        let len = dp.len();
        let mode = depth_mode(&dp);
        ContigCtx {
            hp_run: vec![0u8; len],
            hp_fwd: vec![0u8; len],
            str_span: vec![0u16; len],
            dp,
            na: vec![0u8; len],
            lowcov: vec![0u8; len],
            circular,
            circular_by: "verdict",
            mode,
            mean: 0.0,
            depth_skipped: None,
        }
    }
    #[test]
    fn depth_step_wrap_and_segment() {
        let n = 60_000usize;
        // High depth on [10_000, 50_000), low elsewhere: an interior up step and
        // an interior down step, with no origin imbalance (both ends are low).
        let dp: Vec<u32> = (0..n).map(|i| if (10_000..50_000).contains(&i) { 40 } else { 10 }).collect();
        let ctx = step_ctx(dp, true);
        let p = QcParams::default();
        let mut regs = Vec::new();
        depth_step_regions(std::slice::from_ref(&ctx), &p, 10, &mut regs);
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        let ups: Vec<&Region> = steps.iter().copied().filter(|r| r.reason == "step_up").collect();
        let downs: Vec<&Region> = steps.iter().copied().filter(|r| r.reason == "step_down").collect();
        assert_eq!(ups.len(), 1, "{regs:?}");
        assert_eq!(downs.len(), 1, "{regs:?}");
        // Neither peak is at the origin, so neither is dropped.
        assert!(steps.iter().all(|r| r.start != 0 && r.end != n as u64), "{steps:?}");
        let up = ups[0];
        assert_eq!(up.tier, "red", "{up:?}");
        assert_eq!((up.start, up.end), (10_000 - 50, 10_000 + 50), "{up:?}");
        assert!(up.extra.contains("delta=30.00"), "{}", up.extra);
        assert!(up.extra.contains("frac=0.750"), "{}", up.extra);
        let down = downs[0];
        assert_eq!(down.tier, "red", "{down:?}");
        assert_eq!((down.start, down.end), (50_000 - 50, 50_000 + 50), "{down:?}");
        assert!(down.extra.contains("delta=-30.00"), "{}", down.extra);
        assert!(down.extra.contains("frac=0.750"), "{}", down.extra);
        // One `loss` segment from the down step at 50 kb to the up step at 10 kb:
        // the shorter way round wraps the origin (two BED rows) and the inside
        // (the low-depth half) is below the flanks.
        let segs: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_segment").collect();
        assert_eq!(segs.len(), 2, "{segs:?}");
        assert!(segs.iter().all(|r| r.reason == "loss" && r.tier == "red"), "{segs:?}");
        assert!(segs.iter().any(|r| r.start == 50_000 && r.end == n as u64), "{segs:?}");
        assert!(segs.iter().any(|r| r.start == 0 && r.end == 10_000), "{segs:?}");
        assert!(!segs.iter().any(|r| r.start == r.end), "empty segment piece: {segs:?}");

        // A low-coverage side marks the step `lowcov=1`.
        let mut regs = Vec::new();
        depth_step_regions(std::slice::from_ref(&ctx), &p, 20, &mut regs);
        assert!(regs.iter().filter(|r| r.label == "depth_step").all(|r| r.lowcov), "{regs:?}");
    }
    /// Emission-level boundary: a step whose interval starts at `0` (pos 50) or
    /// ends at `len` (pos len-50) is dropped; at `51` / `len-51` it is kept.
    #[test]
    fn depth_step_origin_emission_boundary() {
        let n = 60_000i64;
        let p = QcParams::default();
        let mk = |lo: i64, hi: i64| -> Vec<u32> {
            (0..n as usize).map(|i| if (lo..hi).contains(&(i as i64)) { 40 } else { 10 }).collect()
        };
        let run = |dp: Vec<u32>| {
            let mut regs = Vec::new();
            let dropped = depth_step_regions(std::slice::from_ref(&step_ctx(dp, true)), &p, 10, &mut regs);
            (dropped, regs)
        };
        // Up step at exactly 50: its row `(0, 100)` touches the origin -> dropped;
        // the down step at 30 kb is interior and kept.
        let (d, regs) = run(mk(50, 30_000));
        assert_eq!(d, vec![1u64], "{regs:?}");
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        assert!(steps.iter().any(|r| (r.start, r.end) == (29_950, 30_050)), "{steps:?}");
        assert!(steps.iter().all(|r| r.start != 0 && r.end != n as u64), "{steps:?}");
        // Up step at 51: row `(1, 101)`, interior -> kept.
        let (d, regs) = run(mk(51, 30_000));
        assert_eq!(d, vec![0u64], "{regs:?}");
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        assert!(steps.iter().any(|r| (r.start, r.end) == (1, 101)), "{steps:?}");
        // Down step at exactly len-50: its row `(len-100, len)` -> dropped.
        let (d, regs) = run(mk(30_000, n - 50));
        assert_eq!(d, vec![1u64], "{regs:?}");
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        assert!(steps.iter().any(|r| (r.start, r.end) == (29_950, 30_050)), "{steps:?}");
        assert!(steps.iter().all(|r| r.end != n as u64), "{steps:?}");
        // Down step at len-51: row `(len-101, len-1)` -> kept.
        let (d, regs) = run(mk(30_000, n - 51));
        assert_eq!(d, vec![0u64], "{regs:?}");
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        assert!(steps.iter().any(|r| (r.start, r.end) == (n as u64 - 101, n as u64 - 1)), "{steps:?}");
    }
    /// The origin boundary is defined by the emitted interval, not the refined
    /// position: `[pos - 50, pos + 50)` touching `0` or `len` drops.
    #[test]
    fn depth_step_origin_boundary() {
        let n = 60_000i64;
        // Drop: the interval touches 0 or len (split rows, or a single row
        // starting at 0 / ending at len).
        for pos in [0i64, 1, 4, 49, 50, n - 1, n - 4, n - 49, n - 50] {
            assert!(step_at_origin(pos, n), "pos={pos} must be at the origin");
        }
        assert!(step_at_origin(n, n), "pos == len is the same point as 0");
        // Keep: strictly interior, so the row is `(pos - 50, pos + 50)`.
        for pos in [51i64, 100, n - 100, n - 51] {
            assert!(!step_at_origin(pos, n), "pos={pos} must not be at the origin");
        }
    }
    /// An end/start depth imbalance on a wrapping contig peaks at the origin
    /// (an aligner/representation artefact) and is dropped; an interior step a
    /// few hundred bp away is kept.
    #[test]
    fn depth_step_origin_dropped_interior_kept() {
        let n = 60_000usize;
        let p = QcParams::default();
        // 40x on [0, 600), 10x elsewhere: an up step at the origin and a down
        // step at 600 (>= 50 bp away, so not the origin).
        let mk = || (0..n).map(|i| if i < 600 { 40 } else { 10 }).collect::<Vec<u32>>();
        let ctx = step_ctx(mk(), true);
        let mut regs = Vec::new();
        let dropped = depth_step_regions(std::slice::from_ref(&ctx), &p, 10, &mut regs);
        assert_eq!(dropped, vec![1u64], "one origin step dropped");
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        assert_eq!(steps.len(), 1, "only the interior step survives: {regs:?}");
        assert_eq!(steps[0].reason, "step_down", "{:?}", steps[0]);
        assert_eq!((steps[0].start, steps[0].end), (550, 650), "{:?}", steps[0]);
        assert!(steps.iter().all(|r| r.start != 0 && r.end != n as u64), "{steps:?}");
        // A single surviving step cannot anchor a segment.
        assert!(regs.iter().all(|r| r.label != "depth_segment"), "{regs:?}");

        // A linear (non-wrapping) contig is not gated: nothing is dropped and the
        // near-end step is still emitted (its interval is clamped, never split).
        let ctx = step_ctx(mk(), false);
        let mut regs = Vec::new();
        let dropped = depth_step_regions(std::slice::from_ref(&ctx), &p, 10, &mut regs);
        assert_eq!(dropped, vec![0u64]);
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        assert!(steps.iter().any(|r| r.reason == "step_down"), "{regs:?}");
        assert!(steps.iter().all(|r| r.start != 0 && r.end != n as u64), "{regs:?}");
    }
    #[test]
    fn depth_step_linear_no_wrap() {
        let n = 60_000usize;
        let step = 30_000usize;
        let dp: Vec<u32> = (0..n).map(|i| if i < step { 40 } else { 10 }).collect();
        let ctx = step_ctx(dp, false);
        let p = QcParams::default();
        let mut regs = Vec::new();
        depth_step_regions(std::slice::from_ref(&ctx), &p, 10, &mut regs);
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        assert_eq!(steps.len(), 1, "{regs:?}");
        assert_eq!(steps[0].reason, "step_down", "{:?}", steps[0]);
        // No origin wrap and a single peak -> no paired segment.
        assert!(regs.iter().all(|r| r.label != "depth_segment"), "{regs:?}");
    }
    #[test]
    fn depth_step_flat_high_depth_no_peak() {
        let p = QcParams::default();
        for circular in [true, false] {
            let ctx = step_ctx(vec![40u32; 60_000], circular);
            let mut regs = Vec::new();
            depth_step_regions(std::slice::from_ref(&ctx), &p, 10, &mut regs);
            assert!(regs.is_empty(), "circular={circular}: {regs:?}");
        }
    }
    /// The Poisson floor is inside the `1.4826 *` factor: on a flat contig the
    /// MAD is 0, so `z = Δ / (1.4826 * sqrt(depth/w))` and the yellow bar needs
    /// `Δ >= 8 * 1.4826 * sqrt(depth/w)`.
    #[test]
    fn depth_step_score_is_floor_dominated() {
        let n = 60_000usize;
        let depth = 100_000u32;
        let p = QcParams::default();
        let mk = |bump: u32| {
            let mut dp = vec![depth; n];
            for x in &mut dp[30_000..30_100] {
                *x += bump;
            }
            step_ctx(dp, true)
        };
        // floor = sqrt(100_000 / 100) = sqrt(1000) = 31.62 -> z = Δ / 46.88.
        let ctx = mk(500);
        let mut regs = Vec::new();
        depth_step_regions(std::slice::from_ref(&ctx), &p, 10, &mut regs);
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        assert_eq!(steps.len(), 1, "{regs:?}");
        assert_eq!(steps[0].tier, "yellow", "{:?}", steps[0]);
        let z: f64 = steps[0].extra.split("z=").nth(1).unwrap().split(';').next().unwrap().parse().unwrap();
        assert!((z.abs() - 500.0 / (1.4826 * 1000f64.sqrt())).abs() < 0.5, "z={z}");
        // A smaller bump (Δ = 300 -> z ≈ 6.4) stays below the yellow bar.
        let ctx = mk(300);
        let mut regs = Vec::new();
        depth_step_regions(std::slice::from_ref(&ctx), &p, 10, &mut regs);
        assert!(regs.iter().all(|r| r.label != "depth_step"), "{regs:?}");
    }
    /// `step_gather` collects by the actual positions, not by `k * stride`:
    /// linear grids start at `w`, so the old index arithmetic was off by `w`.
    #[test]
    fn step_gather_uses_actual_positions() {
        let pos: Vec<i64> = (0..10).map(|k| 100 + k * 10).collect(); // 100..190
        let delta: Vec<f64> = pos.iter().map(|&p| p as f64).collect();
        let mut vals = Vec::new();
        step_gather(&pos, &delta, 120, 150, &mut vals);
        assert_eq!(vals, vec![120.0, 130.0, 140.0, 150.0]);
        // A circular grid starts at 0.
        let pos0: Vec<i64> = (0..10).map(|k| k * 10).collect();
        let d0: Vec<f64> = pos0.iter().map(|&p| p as f64).collect();
        let mut vals = Vec::new();
        step_gather(&pos0, &d0, 20, 50, &mut vals);
        assert_eq!(vals, vec![20.0, 30.0, 40.0, 50.0]);
    }
    /// A linear contig whose whole neighbourhood is a depth wave (so the local
    /// MAD is genuinely nonzero, ~10) with one real level change. The step must
    /// be scored against that MAD, not the (much smaller) Poisson floor: a
    /// level change of 180 at depth ~130 gives |z| ~ 12, not ~ 98.
    #[test]
    fn depth_step_linear_nonzero_neighbourhood() {
        let n = 40_000usize;
        let step = 24_000usize;
        let dp: Vec<u32> = (0..n)
            .map(|i| {
                let base = if i < step { 200 } else { 40 };
                base + 20 * ((i / 100) % 2) as u32
            })
            .collect();
        let ctx = step_ctx(dp, false);
        let p = QcParams::default();
        let mut regs = Vec::new();
        depth_step_regions(std::slice::from_ref(&ctx), &p, 10, &mut regs);
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        let down = steps
            .iter()
            .find(|r| r.reason == "step_down")
            .unwrap_or_else(|| panic!("no step_down: {regs:?}"));
        let z: f64 = down.extra.split("z=").nth(1).unwrap().split(';').next().unwrap().parse().unwrap();
        assert!(z.abs() >= 8.0 && z.abs() <= 30.0, "z={z} (MAD not used?): {down:?}");
    }
    /// Runs (all singletons here): up,down,up,down gives two gains left to
    /// right, and a same-sign neighbour is skipped.
    #[test]
    fn depth_step_pairing_two_gains() {
        let peaks = vec![
            StepPeak { pos: 1_000, z: 40.0, delta: 1.0 },
            StepPeak { pos: 2_000, z: -30.0, delta: -1.0 },
            StepPeak { pos: 3_000, z: 40.0, delta: 1.0 },
            StepPeak { pos: 4_000, z: -30.0, delta: -1.0 },
        ];
        let segs = step_pair_segments(&peaks, true, 10_000);
        let got: Vec<(i64, i64)> = segs.iter().map(|(a, b)| (a.pos, b.pos)).collect();
        assert_eq!(got, vec![(1_000, 2_000), (3_000, 4_000)], "{segs:?}");
    }
    /// Each step joins at most one segment: a down step between two up steps
    /// goes to the earlier pair, so up,down,up (linear) gives one gain.
    #[test]
    fn depth_step_pairing_each_step_once() {
        let peaks = vec![
            StepPeak { pos: 1_000, z: 40.0, delta: 1.0 },
            StepPeak { pos: 2_000, z: -30.0, delta: -1.0 },
            StepPeak { pos: 3_000, z: 40.0, delta: 1.0 },
        ];
        let segs = step_pair_segments(&peaks, false, 10_000);
        assert_eq!(segs.len(), 1, "{segs:?}");
        assert!(segs[0].0.pos == 1_000 && segs[0].1.pos == 2_000, "{segs:?}");
    }
    /// A dip/loss is a down step immediately followed by an up step.
    #[test]
    fn depth_step_pairing_dip() {
        let peaks = vec![
            StepPeak { pos: 1_000, z: -40.0, delta: -1.0 },
            StepPeak { pos: 2_000, z: 40.0, delta: 1.0 },
        ];
        let segs = step_pair_segments(&peaks, false, 10_000);
        assert_eq!(segs.len(), 1, "{segs:?}");
        assert!(segs[0].0.delta < 0.0 && segs[0].0.pos == 1_000 && segs[0].1.pos == 2_000, "{segs:?}");
    }
    /// Runs, with each run represented by its strongest step: a gain runs from
    /// the strongest step of an up run to the strongest step of the following
    /// down run (the staircase case).
    #[test]
    fn depth_step_pairing_runs_staircase() {
        let peaks = vec![
            StepPeak { pos: 100, z: 40.0, delta: 1.0 },
            StepPeak { pos: 200, z: 40.0, delta: 1.0 },
            StepPeak { pos: 300, z: -30.0, delta: -1.0 },
            StepPeak { pos: 400, z: -30.0, delta: -1.0 },
            StepPeak { pos: 500, z: 40.0, delta: 1.0 },
            StepPeak { pos: 600, z: 40.0, delta: 1.0 },
        ];
        let segs = step_pair_segments(&peaks, false, 10_000);
        let got: Vec<(i64, i64)> = segs.iter().map(|(a, b)| (a.pos, b.pos)).collect();
        // Strongest of the first up run and of the first down run (ties: 100, 300).
        assert_eq!(got, vec![(100, 300)], "{segs:?}");
    }
    /// A run is represented by its strongest step, not its far ends: a weak
    /// step far from the strong one does not move a segment endpoint.
    #[test]
    fn depth_step_pairing_strongest_in_run() {
        // ONT-FAM29080: weak ups at 0 and 1815990, strong up at 543393 (all one
        // origin-spanning up run), strong down at 563743.
        let peaks = vec![
            StepPeak { pos: 0, z: 11.0, delta: 1.0 },
            StepPeak { pos: 543_393, z: 29.0, delta: 1.0 },
            StepPeak { pos: 563_743, z: -28.0, delta: -1.0 },
            StepPeak { pos: 1_815_990, z: 8.0, delta: 1.0 },
        ];
        let segs = step_pair_segments(&peaks, true, 10_000_000);
        let got: Vec<(i64, i64)> = segs.iter().map(|(a, b)| (a.pos, b.pos)).collect();
        assert_eq!(got, vec![(543_393, 563_743)], "{segs:?}");
    }
    /// The down run's strongest step is the endpoint, not its first: a weak
    /// near step before a strong far step does not shorten the segment
    /// (PB-FAM20853).
    #[test]
    fn depth_step_pairing_strongest_down_run() {
        let peaks = vec![
            StepPeak { pos: 1_862_513, z: 61.0, delta: 1.0 },
            StepPeak { pos: 1_888_131, z: -11.0, delta: -1.0 },
            StepPeak { pos: 1_903_397, z: -59.0, delta: -1.0 },
        ];
        let segs = step_pair_segments(&peaks, true, 10_000_000);
        let got: Vec<(i64, i64)> = segs.iter().map(|(a, b)| (a.pos, b.pos)).collect();
        assert_eq!(got, vec![(1_862_513, 1_903_397)], "{segs:?}");
    }
    /// Circular: a same-sign run spanning the origin is one run. On a linear
    /// contig the same steps form two runs.
    #[test]
    fn depth_step_runs_grouping() {
        // Sorted: up@100, up@200, down@5000, up@9900 -> the up run
        // {9900, 100, 200} spans the origin.
        let peaks = vec![
            StepPeak { pos: 100, z: 40.0, delta: 1.0 },
            StepPeak { pos: 200, z: 40.0, delta: 1.0 },
            StepPeak { pos: 5_000, z: -30.0, delta: -1.0 },
            StepPeak { pos: 9_900, z: 40.0, delta: 1.0 },
        ];
        let circ = step_runs(&peaks, true);
        assert_eq!(circ.len(), 2, "circular runs: {circ:?}");
        assert!(circ[0].strong.delta > 0.0 && circ[0].left == 100 && circ[0].strong.pos == 100);
        assert!(circ[1].strong.delta < 0.0 && circ[1].left == 5_000 && circ[1].strong.pos == 5_000);
        let segs = step_pair_segments(&peaks, true, 10_000);
        let got: Vec<(i64, i64)> = segs.iter().map(|(a, b)| (a.pos, b.pos)).collect();
        assert_eq!(got, vec![(100, 5_000)], "{segs:?}");
        // Linear: the origin up run splits into {100,200} and {9900}.
        let lin = step_runs(&peaks, false);
        assert_eq!(lin.len(), 3, "linear runs: {lin:?}");
    }
    /// Regression for the pbipa_04 wrap bug: up 330801, down 332261, up
    /// 3561590, down 3563048 on a circular contig must give the two 1.46 kb
    /// gains, not two ~1.2 Mb wrapped segments.
    #[test]
    fn depth_step_pairing_pbipa_pattern() {
        let peaks = vec![
            StepPeak { pos: 330_801, z: 50.0, delta: 1.0 },
            StepPeak { pos: 332_261, z: -50.0, delta: -1.0 },
            StepPeak { pos: 3_561_590, z: 40.0, delta: 1.0 },
            StepPeak { pos: 3_563_048, z: -40.0, delta: -1.0 },
        ];
        let segs = step_pair_segments(&peaks, true, 3_600_000);
        let got: Vec<(i64, i64)> = segs.iter().map(|(a, b)| (a.pos, b.pos)).collect();
        assert_eq!(got, vec![(330_801, 332_261), (3_561_590, 3_563_048)], "{segs:?}");
    }
    /// Regression for the FAM20862 pflye_05 bug: the run spanning the origin is
    /// ordered first by its leftmost step, so the old greedy matching paired it
    /// with the wrong neighbour and emitted a ~2.4 Mb wrapped *loss* instead of
    /// the real 2.37 kb gain. With both matchings and the smaller total length
    /// the up run (strong at 2413744) pairs with the following down run (strong
    /// at 2416114).
    #[test]
    fn depth_step_pairing_origin_run_picks_gain() {
        let peaks = vec![
            StepPeak { pos: 0, z: -44.0, delta: -1.0 },
            StepPeak { pos: 49_690, z: 45.0, delta: 1.0 },
            StepPeak { pos: 2_413_744, z: 1_123.0, delta: 1.0 },
            StepPeak { pos: 2_416_114, z: -1_090.0, delta: -1.0 },
            StepPeak { pos: 2_501_639, z: -44.0, delta: -1.0 },
        ];
        let segs = step_pair_segments(&peaks, true, 2_501_689);
        let got: Vec<(i64, i64)> = segs.iter().map(|(a, b)| (a.pos, b.pos)).collect();
        assert_eq!(got, vec![(2_413_744, 2_416_114)], "{segs:?}");
        assert!(segs[0].0.delta > 0.0, "must be a gain: {segs:?}");
    }
    /// Regression for the FAM22000 pflye_04 bug: the down run {2465212, 0} is
    /// strongest at 0 (a |z| tie, smaller position), so the old matching made a
    /// 2.42 Mb *loss* from 0 to 2422594. The run-1 matching pairs the up run
    /// (2422594) with that down run and, wrapping the origin, gives the 43 kb
    /// gain to the contig end.
    #[test]
    fn depth_step_pairing_origin_run_picks_short_gain() {
        let peaks = vec![
            StepPeak { pos: 0, z: -42.88, delta: -1.0 },
            StepPeak { pos: 2_422_594, z: 37.28, delta: 1.0 },
            StepPeak { pos: 2_465_212, z: -42.88, delta: -1.0 },
        ];
        let segs = step_pair_segments(&peaks, true, 2_465_262);
        let got: Vec<(i64, i64)> = segs.iter().map(|(a, b)| (a.pos, b.pos)).collect();
        assert_eq!(got, vec![(2_422_594, 0)], "{segs:?}");
        assert!(segs[0].0.delta > 0.0, "must be a gain: {segs:?}");
    }
    /// A segment whose measured inside/flank ratio contradicts its direction is
    /// dropped: a gain with `ratio <= 1`, a loss with `ratio >= 1`.
    #[test]
    fn depth_step_segment_ratio_filter() {
        let p = QcParams::default();
        let emit = |dp: Vec<u32>, a: StepPeak, b: StepPeak| -> usize {
            let ctx = step_ctx(dp, true);
            let mut regs = Vec::new();
            emit_step_segment(0, &ctx, 60_000, &a, &b, 100, 10, &p, &mut regs);
            regs.len()
        };
        let up_at = |pos: i64| StepPeak { pos, z: 25.0, delta: 1.0 };
        let down_at = |pos: i64| StepPeak { pos, z: -25.0, delta: -1.0 };
        // gain = up@20000 then down@40000; loss = down@20000 then up@40000
        // (both span the same [20_000, 40_000) interior).
        let gain = (up_at(20_000), down_at(40_000));
        let loss = (down_at(20_000), up_at(40_000));
        // Inside [20_000, 40_000) at half the flanks: ratio 0.5.
        let low: Vec<u32> = (0..60_000).map(|i| if (20_000..40_000).contains(&i) { 10 } else { 20 }).collect();
        assert_eq!(emit(low.clone(), gain.0, gain.1), 0, "gain with ratio 0.5 must be dropped");
        assert_eq!(emit(low, loss.0, loss.1), 1, "loss with ratio 0.5 must be kept");
        // Inside at twice the flanks: ratio 2.
        let high: Vec<u32> = (0..60_000).map(|i| if (20_000..40_000).contains(&i) { 40 } else { 20 }).collect();
        assert_eq!(emit(high.clone(), gain.0, gain.1), 1, "gain with ratio 2 must be kept");
        assert_eq!(emit(high, loss.0, loss.1), 0, "loss with ratio 2 must be dropped");
    }
    /// Finding 5: the segment tier is the lower of its two endpoint tiers.
    #[test]
    fn depth_step_segment_tier_is_lower_end() {
        let p = QcParams::default();
        let emit = |dp: Vec<u32>, a: StepPeak, b: StepPeak| -> String {
            let ctx = step_ctx(dp, true);
            let mut regs = Vec::new();
            emit_step_segment(0, &ctx, 60_000, &a, &b, 100, 10, &p, &mut regs);
            assert_eq!(regs.len(), 1, "{regs:?}");
            regs[0].tier.to_string()
        };
        // up 10->40 (red: f = 30/40 > 0.5) then down 40->20 (orange: f = 20/40 = 0.5).
        let dp1: Vec<u32> = (0..60_000).map(|i| if i < 20_000 { 10 } else if i < 40_000 { 40 } else { 20 }).collect();
        let t = emit(dp1, StepPeak { pos: 20_000, z: 25.0, delta: 30.0 }, StepPeak { pos: 40_000, z: -25.0, delta: -20.0 });
        assert_eq!(t, "orange");
        // up 10->20 (orange) then down 20->10 (yellow: z = 10).
        let dp2: Vec<u32> = (0..60_000).map(|i| if i < 20_000 { 10 } else if i < 40_000 { 20 } else { 10 }).collect();
        let t = emit(dp2, StepPeak { pos: 20_000, z: 25.0, delta: 10.0 }, StepPeak { pos: 40_000, z: -10.0, delta: -10.0 });
        assert_eq!(t, "yellow");
    }
    /// Finding 4: a low-coverage contig (`depth_skipped`) still gets steps.
    #[test]
    fn depth_step_lowcov_contig_not_skipped() {
        use crate::scan::{DPF, NCOV};
        let n = 60_000usize;
        let mut c = ContigCounts {
            name: "chr".to_string(),
            len: n,
            bases: vec![[0u32; 8]; n],
            cov: vec![[0u32; NCOV]; n],
            clips: Vec::new(),
            edge_clips: 0,
            aligned_bases: 0,
        };
        for i in 0..n {
            c.cov[i][DPF] = if i < 40_000 { 8 } else { 2 };
        }
        let refseq: Vec<u8> = (0..n).map(|i| b"ACGT"[i % 4]).collect();
        let p = QcParams::default();
        let ctx = contig_ctx(&c, &refseq, true, &p, 10, p.lowcov_rel);
        assert_eq!(ctx.depth_skipped, Some("low_coverage"));
        let mut regs = Vec::new();
        depth_step_regions(std::slice::from_ref(&ctx), &p, 10, &mut regs);
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        assert!(!steps.is_empty(), "low-coverage contig skipped: {regs:?}");
        assert!(steps.iter().all(|r| r.lowcov), "{regs:?}");
    }
    /// Finding 7: non-maximum suppression drops a weaker step within ±nms.
    #[test]
    fn depth_step_nms_suppresses_close_peak() {
        let n = 60_000usize;
        let dp: Vec<u32> = (0..n).map(|i| if (30_000..30_400).contains(&i) { 100 } else { 40 }).collect();
        let ctx = step_ctx(dp, true);
        let p = QcParams::default();
        let mut regs = Vec::new();
        depth_step_regions(std::slice::from_ref(&ctx), &p, 10, &mut regs);
        let steps: Vec<&Region> = regs.iter().filter(|r| r.label == "depth_step").collect();
        assert_eq!(steps.len(), 1, "{regs:?}");
    }
}
