//! `qc` region significance scores (Poisson local background).
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use super::*;

// Region significance scores (qc tiers spec, step 2).
//
// Every region signal reduces to "k reads disagree out of n". The score is
// `s = -log10 P(X >= k)` under a Poisson with mean `lambda` fitted from the
// signal's own local background in a +/-`SIG_SPAN` window around the call,
// excluding the call itself:
//
//   rate = (background events in the window) / (covered bp in the window)
//   lambda = rate * call_cov
//
// where `call_cov` is the covered bp (`Σ dp`) of the call region for clips and
// large indels (`depth * width`) and the depth at the point for junction ends
// and bridging reads. The background mean is floored at `SIG_LAMBDA_FLOOR` (half
// a count per window) so `s` stays finite. `f` is the signal's existing
// disagreement fraction. Scores are heuristic and do not change any tier.

/// Half-width (bp) of the local background window.
pub(crate) const SIG_SPAN: u64 = 25_000;
/// Floor of the Poisson background mean (half a count) so `s` stays finite.
pub(crate) const SIG_LAMBDA_FLOOR: f64 = 0.5;

/// `s = -log10 P(X >= k)` for `X ~ Poisson(max(lambda, floor))`.
pub(crate) fn poisson_score(k: u64, lambda: f64) -> f64 {
    let lam = if lambda.is_finite() { lambda.max(SIG_LAMBDA_FLOOR) } else { SIG_LAMBDA_FLOOR };
    let lnp = stats::poisson_ln_sf(k, lam);
    if lnp <= f64::NEG_INFINITY { f64::INFINITY } else { -lnp / std::f64::consts::LN_10 }
}

/// The `;s=..;f=..` suffix appended to a region's `extra`.
pub(crate) fn sig_suffix(s: f64, f: f64) -> String {
    format!(";s={s:.2};f={f:.4}")
}

/// Up to two non-wrapping `[start, end)` pieces of the window
/// `[center - span, center + span)` on a contig of `len`: a circular window
/// that wraps the origin is split into two pieces, and a length `<= 2*span`
/// gives the whole contig as one piece. Returns the pieces and how many are
/// valid. Coverage, membership and exclusions all use these same pieces.
fn sig_window_pieces(center: u64, span: u64, len: u64, circular: bool) -> ([(u64, u64); 2], usize) {
    let empty = ([(0, 0), (0, 0)], 0);
    if len == 0 || span == 0 {
        return empty;
    }
    if !circular {
        let lo = center.saturating_sub(span);
        let hi = center.saturating_add(span).min(len);
        return if hi > lo { ([(lo, hi), (0, 0)], 1) } else { empty };
    }
    let cnt = span.saturating_mul(2).min(len);
    if cnt == len {
        return ([(0, len), (0, 0)], 1);
    }
    let start = (center as i64 - span as i64).rem_euclid(len as i64) as u64;
    let end = start + cnt;
    if end <= len {
        ([(start, end), (0, 0)], 1)
    } else {
        ([(start, len), (0, end - len)], 2)
    }
}

/// Prefix sums of a per-position array (length `len`), so a window sum is O(1).
pub(crate) fn prefix_u32(v: &[u32]) -> Vec<u64> {
    let mut p = Vec::with_capacity(v.len() + 1);
    let mut s = 0u64;
    p.push(0);
    for &x in v {
        s += x as u64;
        p.push(s);
    }
    p
}

/// Merge overlapping/adjacent exclusion intervals into their union (sorted,
/// non-overlapping). Exclusions are non-wrapping `[start, end)` within `[0, len)`.
fn merge_excls(excls: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut v: Vec<(u64, u64)> = excls.iter().copied().filter(|&(a, b)| b > a).collect();
    v.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::with_capacity(v.len());
    for (a, b) in v {
        match out.last_mut() {
            Some(last) if a <= last.1 => last.1 = last.1.max(b),
            _ => out.push((a, b)),
        }
    }
    out
}

/// Sum of `pref` over the window `[center - span, center + span)`, minus the
/// intersection with the union of the exclusion intervals (exclusions are
/// non-wrapping `[start, end)` within `[0, len)`). Overlapping exclusions are
/// merged first, and each merged interval is clipped to each window piece
/// before subtracting, so neither a call wider than the window nor two
/// overlapping exclusions can underflow.
pub(crate) fn prefix_window_excl(
    pref: &[u64],
    len: u64,
    circular: bool,
    center: u64,
    span: u64,
    excls: &[(u64, u64)],
) -> u64 {
    let (pieces, k) = sig_window_pieces(center, span, len, circular);
    let merged_storage;
    let merged: &[(u64, u64)] = if excls.len() <= 1 {
        excls
    } else {
        merged_storage = merge_excls(excls);
        &merged_storage
    };
    let mut total = 0u64;
    for i in 0..k {
        let (a, b) = pieces[i];
        let mut sum = pref[b as usize] - pref[a as usize];
        for &(e0, e1) in merged {
            let lo = a.max(e0);
            let hi = b.min(e1);
            if hi > lo {
                sum -= pref[hi as usize] - pref[lo as usize];
            }
        }
        total += sum;
    }
    total
}

/// Sum of `pref` over `[s, e)` (`s <= e <= len`, no wrap).
pub(crate) fn prefix_range(pref: &[u64], s: u64, e: u64) -> u64 {
    pref[e as usize] - pref[s as usize]
}

/// Number of sorted `positions` in `[lo, hi)`.
fn count_in_sorted(positions: &[u32], lo: u64, hi: u64) -> u64 {
    if hi <= lo || positions.is_empty() {
        return 0;
    }
    let lo = lo.min(u32::MAX as u64) as u32;
    let hi = hi.min(u32::MAX as u64) as u32;
    let a = positions.partition_point(|&x| x < lo);
    let b = positions.partition_point(|&x| x < hi);
    (b - a) as u64
}

/// Number of sorted `positions` in the window `[center - span, center + span)`,
/// excluding the union of the exclusion intervals (merged, then clipped to each
/// window piece).
pub(crate) fn count_window_excl(
    positions: &[u32],
    len: u64,
    circular: bool,
    center: u64,
    span: u64,
    excls: &[(u64, u64)],
) -> u64 {
    let (pieces, k) = sig_window_pieces(center, span, len, circular);
    let merged_storage;
    let merged: &[(u64, u64)] = if excls.len() <= 1 {
        excls
    } else {
        merged_storage = merge_excls(excls);
        &merged_storage
    };
    let mut n = 0u64;
    for i in 0..k {
        let (a, b) = pieces[i];
        n += count_in_sorted(positions, a, b);
        for &(e0, e1) in merged {
            let lo = a.max(e0);
            let hi = b.min(e1);
            if hi > lo {
                n -= count_in_sorted(positions, lo, hi);
            }
        }
    }
    n
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poisson_score_values() {
        // k = 0 is never significant.
        assert_eq!(poisson_score(0, 3.0), 0.0);
        // k = 1, lambda floored at 0.5: -log10(1 - e^-0.5).
        let want = -(1.0 - (-0.5f64).exp()).log10();
        assert!((poisson_score(1, 0.5) - want).abs() < 1e-9, "{}", poisson_score(1, 0.5));
        // The floor binds: a zero (or negative) background equals the floor.
        assert_eq!(poisson_score(4, 0.0), poisson_score(4, 0.5));
        assert_eq!(poisson_score(4, -2.0), poisson_score(4, 0.5));
        // k >= 5 at the floor is ~3.76.
        let s5 = poisson_score(5, 0.5);
        assert!((3.75..3.78).contains(&s5), "{s5}");
        // Score increases with k at fixed lambda.
        assert!(poisson_score(10, 0.5) > s5);
        // Score decreases as the background mean grows.
        assert!(poisson_score(10, 5.0) < poisson_score(10, 1.0));
    }
    #[test]
    fn sig_window_pieces_and_membership() {
        // Linear: half-open [c-span, c+span), clamped at the ends.
        let (p, k) = sig_window_pieces(10, 5, 100, false);
        assert_eq!((k, p[0]), (1, (5, 15)));
        let (p, k) = sig_window_pieces(2, 5, 100, false);
        assert_eq!((k, p[0]), (1, (0, 7)));
        // len > 2*span: a circular window that wraps splits into two pieces.
        let (p, k) = sig_window_pieces(2, 5, 100, true);
        assert_eq!(k, 2);
        assert_eq!((p[0], p[1]), ((97, 100), (0, 7)));
        // Membership comes from the same half-open pieces (boundaries inclusive
        // on the low side, exclusive on the high side).
        let member = |q| (0..k).any(|i| q >= p[i].0 && q < p[i].1);
        assert!(member(97) && member(99) && member(0) && member(6));
        assert!(!member(7) && !member(96));
        // len == 2*span: the whole contig is one piece.
        let (p, k) = sig_window_pieces(3, 5, 10, true);
        assert_eq!((k, p[0]), (1, (0, 10)));
        // len < 2*span: still one whole-contig piece.
        let (p, k) = sig_window_pieces(3, 5, 6, true);
        assert_eq!((k, p[0]), (1, (0, 6)));
    }
    #[test]
    fn prefix_window_and_range() {
        let v = [1u32, 2, 3, 4, 5, 6];
        let p = prefix_u32(&v);
        assert_eq!(p, vec![0, 1, 3, 6, 10, 15, 21]);
        // Linear window clamps at the ends.
        assert_eq!(prefix_window_excl(&p, 6, false, 3, 2, &[]), 2 + 3 + 4 + 5); // [1,5)
        assert_eq!(prefix_window_excl(&p, 6, false, 0, 2, &[]), 1 + 2); // [0,2)
        // Circular window wraps once.
        assert_eq!(prefix_window_excl(&p, 6, true, 0, 2, &[]), 5 + 6 + 1 + 2); // [4,6)+[0,2)
        // Range sums.
        assert_eq!(prefix_range(&p, 1, 4), 2 + 3 + 4);
        assert_eq!(prefix_range(&p, 3, 3), 0);
    }
    #[test]
    fn prefix_window_excl_clips_to_window() {
        // All-ones depth: prefix sums are the position index.
        let p = prefix_u32(&vec![1u32; 100]);
        // A call wider than the window (reviewer repro): subtracts only the
        // window's own intersection, never underflows.
        assert_eq!(prefix_window_excl(&p, 100, false, 50, 5, &[(0, 100)]), 0);
        assert_eq!(prefix_window_excl(&p, 100, false, 50, 5, &[(40, 60)]), 0);
        // Partial overlap subtracts only the intersection [48,55): 10 - 7 = 3.
        assert_eq!(prefix_window_excl(&p, 100, false, 50, 5, &[(48, 60)]), 3);
        // Circular pieces (97,100)+(0,7): [95,100) removes 3 -> 7.
        assert_eq!(prefix_window_excl(&p, 100, true, 2, 5, &[(95, 100)]), 7);
        // No exclusion equals the plain window sum.
        assert_eq!(prefix_window_excl(&p, 100, false, 50, 5, &[]), 10);
        assert_eq!(prefix_window_excl(&p, 100, true, 2, 5, &[]), 10);
        // Two exclusions (junction breakpoints) each clipped to the pieces.
        assert_eq!(prefix_window_excl(&p, 100, false, 50, 10, &[(40, 45), (55, 60)]), 10);
    }
    #[test]
    fn count_window_excl_uses_sorted_starts() {
        let starts: Vec<u32> = vec![0, 10, 20, 30, 40, 50];
        // Window [15,45) holds 20,30,40; excluding [25,35) removes 30 -> 2.
        assert_eq!(count_window_excl(&starts, 100, false, 30, 15, &[(25, 35)]), 2);
        // Excluding the whole window -> 0.
        assert_eq!(count_window_excl(&starts, 100, false, 30, 15, &[(0, 100)]), 0);
        // Circular window around 0 = (85,100)+(0,15): 0 and 10 -> 2.
        assert_eq!(count_window_excl(&starts, 100, true, 0, 15, &[]), 2);
    }
    #[test]
    fn exclusions_merge_overlapping_and_nested() {
        let p = prefix_u32(&vec![1u32; 100]);
        // Overlapping: [(0,60),(40,100)] -> union [0,100); window [45,55) -> 0.
        assert_eq!(prefix_window_excl(&p, 100, false, 50, 5, &[(0, 60), (40, 100)]), 0);
        // Nested: [(0,100),(40,60)] -> union [0,100) -> 0.
        assert_eq!(prefix_window_excl(&p, 100, false, 50, 5, &[(0, 100), (40, 60)]), 0);
        // Disjoint stay disjoint: [(40,45),(55,60)] removes 5+5 of [40,60) -> 10.
        assert_eq!(prefix_window_excl(&p, 100, false, 50, 10, &[(40, 45), (55, 60)]), 10);
        // Overlap that only partly covers the window: [(40,50),(48,60)] -> [40,60)
        // removes all of [45,55) -> 0.
        assert_eq!(prefix_window_excl(&p, 100, false, 50, 5, &[(40, 50), (48, 60)]), 0);
        // Circular pieces (97,100)+(0,7): [(95,100),(98,102)] -> [95,102) removes
        // [97,100) = 3 -> 7.
        assert_eq!(prefix_window_excl(&p, 100, true, 2, 5, &[(95, 100), (98, 102)]), 7);
        // Circular, exclusions on different pieces: [(95,100),(0,3)] -> 10-6 = 4.
        assert_eq!(prefix_window_excl(&p, 100, true, 2, 5, &[(95, 100), (0, 3)]), 4);
        // Event counts use the same union.
        let starts: Vec<u32> = vec![0, 10, 20, 30, 40, 50];
        // [(25,35),(30,40)] -> [25,40): window [15,45) keeps 20 and 40 -> 2.
        assert_eq!(count_window_excl(&starts, 100, false, 30, 15, &[(25, 35), (30, 40)]), 2);
        // Nested [(0,100),(10,20)] -> [0,100) removes all -> 0.
        assert_eq!(count_window_excl(&starts, 100, false, 30, 15, &[(0, 100), (10, 20)]), 0);
    }
}
