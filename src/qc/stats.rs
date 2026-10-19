//! Statistics for `qc`: log-gamma, beta-binomial moment fit and upper tail,
//! Wilson bound, Benjamini–Hochberg. No external crates.

/// Natural log of the gamma function (Lanczos, g = 7, n = 9), for `x > 0`.
pub(crate) fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection.
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + G + 0.5;
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

fn ln_beta(a: f64, b: f64) -> f64 {
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}

/// Beta-binomial log pmf.
pub(crate) fn bb_ln_pmf(k: u32, n: u32, a: f64, b: f64) -> f64 {
    let (k, n) = (k as f64, n as f64);
    ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0) + ln_beta(k + a, n - k + b)
        - ln_beta(a, b)
}

/// Natural log of the upper tail `P(X >= x)` for `X ~ BetaBinomial(n, a, b)`.
pub(crate) fn bb_ln_sf(x: u32, n: u32, a: f64, b: f64) -> f64 {
    if x == 0 {
        return 0.0;
    }
    if x > n {
        return f64::NEG_INFINITY;
    }
    let mean = n as f64 * a / (a + b);
    if x as f64 > mean {
        // Sum upward relative to the first term.
        let l0 = bb_ln_pmf(x, n, a, b);
        let mut term = 1.0f64;
        let mut sum = 1.0f64;
        let mut k = x;
        while k < n {
            let kf = k as f64;
            let nf = n as f64;
            term *= (nf - kf) / (kf + 1.0) * (kf + a) / (nf - kf - 1.0 + b);
            sum += term;
            k += 1;
            if term < sum * 1e-17 && kf > mean + 1.0 {
                break;
            }
        }
        (l0 + sum.ln()).min(0.0)
    } else {
        // 1 - P(X <= x-1), summing the lower tail downward.
        let l0 = bb_ln_pmf(x - 1, n, a, b);
        let mut term = 1.0f64;
        let mut sum = 1.0f64;
        let mut k = x - 1;
        while k > 0 {
            let kf = k as f64;
            let nf = n as f64;
            term *= kf / (nf - kf + 1.0) * (nf - kf + b) / (kf - 1.0 + a);
            sum += term;
            k -= 1;
            if term < sum * 1e-17 && kf < mean - 1.0 {
                break;
            }
        }
        let lower = (l0 + sum.ln()).exp();
        (1.0 - lower).max(0.0).ln().max(f64::MIN)
    }
}

/// Natural log of the Poisson pmf.
fn poisson_ln_pmf(i: u64, mu: f64) -> f64 {
    -mu + i as f64 * mu.ln() - ln_gamma(i as f64 + 1.0)
}

/// Natural log of the Poisson upper tail `P(N >= k)` for mean `mu`, computed
/// in log space so it neither underflows nor cancels at large `mu`.
pub(crate) fn poisson_ln_sf(k: u64, mu: f64) -> f64 {
    if k == 0 {
        return 0.0;
    }
    if mu <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if (k as f64) > mu {
        // Upper tail: sum upward relative to pmf(k); terms decrease since k > mu.
        let (mut term, mut sum) = (1.0f64, 1.0f64);
        let mut i = k;
        loop {
            i += 1;
            term *= mu / i as f64;
            sum += term;
            if term < sum * 1e-17 {
                break;
            }
        }
        return poisson_ln_pmf(k, mu) + sum.ln();
    }
    // k <= mu: ln P(N <= k-1) from pmf(k-1) summing downward (terms decrease
    // since k-1 < mu), then ln(1 - cdf) = ln(-expm1(ln cdf)).
    let top = k - 1;
    let (mut term, mut sum) = (1.0f64, 1.0f64);
    let mut i = top;
    while i > 0 {
        term *= i as f64 / mu;
        sum += term;
        i -= 1;
        if term < sum * 1e-17 {
            break;
        }
    }
    let ln_cdf = (poisson_ln_pmf(top, mu) + sum.ln()).min(0.0);
    (-ln_cdf.exp_m1()).ln().max(f64::MIN)
}

/// Integer moment sums for a beta-binomial fit.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct MomentSums {
    pub sites: u64,
    pub sx: u128,
    pub sn: u128,
    pub sx2: u128,
    pub snx: u128,
    pub sn2: u128,
}

impl MomentSums {
    #[inline]
    pub fn push(&mut self, x: u32, n: u32) {
        let (x, n) = (x as u128, n as u128);
        self.sites += 1;
        self.sx += x;
        self.sn += n;
        self.sx2 += x * x;
        self.snx += n * x;
        self.sn2 += n * n;
    }

    pub fn add(&mut self, o: &MomentSums) {
        self.sites += o.sites;
        self.sx += o.sx;
        self.sn += o.sn;
        self.sx2 += o.sx2;
        self.snx += o.snx;
        self.sn2 += o.sn2;
    }
}

/// Fitted beta-binomial background.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize)]
pub struct BbFit {
    pub p: f64,
    pub rho: f64,
    pub alpha: f64,
    pub beta: f64,
}

pub(crate) const RHO_MIN: f64 = 1e-6;
pub(crate) const RHO_MAX: f64 = 0.5;
pub(crate) const P_MIN: f64 = 1e-6;

/// Method-of-moments fit: p = Σx/Σn,
/// ρ = Σ[(x−np)² − np(1−p)] / Σ[np(1−p)(n−1)], clamped to [1e-6, 0.5];
/// p is clamped to [1e-6, 0.5] so the distribution stays proper.
pub(crate) fn fit_bb(m: &MomentSums) -> BbFit {
    let sn = m.sn as f64;
    let p = if m.sn == 0 { P_MIN } else { (m.sx as f64 / sn).clamp(P_MIN, 0.5) };
    let q = 1.0 - p;
    let num = m.sx2 as f64 - 2.0 * p * m.snx as f64 + p * p * m.sn2 as f64 - p * q * sn;
    let den = p * q * (m.sn2 as f64 - sn);
    let rho = if den > 0.0 { (num / den).clamp(RHO_MIN, RHO_MAX) } else { RHO_MIN };
    let k = 1.0 / rho - 1.0;
    BbFit { p, rho, alpha: p * k, beta: q * k }
}

/// Wilson score interval lower bound (95 %).
pub(crate) fn wilson_lower(k: u32, n: u32) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let z = 1.959_963_984_540_054;
    let n = n as f64;
    let ph = k as f64 / n;
    let z2 = z * z;
    let c = ph + z2 / (2.0 * n);
    let m = z * (ph * (1.0 - ph) / n + z2 / (4.0 * n * n)).sqrt();
    ((c - m) / (1.0 + z2 / n)).max(0.0)
}

/// Benjamini–Hochberg q-values for `pvals` with `m` total tests (`m` is
/// raised to `pvals.len()` if smaller). Input order is preserved; ties are
/// resolved by index so the result is deterministic.
pub(crate) fn bh(pvals: &[f64], m: usize) -> Vec<f64> {
    let m = m.max(pvals.len()) as f64;
    let mut idx: Vec<usize> = (0..pvals.len()).collect();
    idx.sort_by(|&a, &b| pvals[a].total_cmp(&pvals[b]).then(a.cmp(&b)));
    let mut q = vec![1.0; pvals.len()];
    let mut run = 1.0f64;
    for (r, &i) in idx.iter().enumerate().rev() {
        let v = (pvals[i] * m / (r + 1) as f64).min(1.0);
        run = run.min(v);
        q[i] = run;
    }
    q
}

/// Reusable state for the banded bit-parallel (Myers/Edlib) edit distance.
/// The query profile and block rows are allocated once and reused across
/// comparisons, so grouping does not allocate per comparison.
///
/// This is the block-wise Myers bit-vector algorithm as implemented by Edlib
/// (`calculateBlock`), with the Ukkonen band from the same reference: only the
/// blocks within `band` of the diagonal are touched, so a comparison costs
/// `O(len · band / 64)` word operations. It is exact — the band only limits
/// work, never the result.
#[derive(Debug, Default)]
pub(crate) struct Banded {
    /// Query (pattern) length of the last [`Banded::set_pattern`] call.
    qlen: usize,
    /// Number of 64-bit blocks covering the query.
    max_num_blocks: usize,
    /// Wildcard padding at the end of the last block (`max_num_blocks*64 - qlen`).
    w: usize,
    /// Query profile: `peq[c * max_num_blocks + b]`, bit `r - b*64` set when
    /// query position `r` equals byte `c` (or `r >= qlen`, the wildcard padding).
    peq: Vec<u64>,
    p: Vec<u64>,
    m: Vec<u64>,
    score: Vec<i32>,
}

impl Banded {
    pub(crate) fn new() -> Self {
        Banded::default()
    }

    /// Build the query profile for `pat`. Reused across many `distance_to`
    /// calls with the same representative.
    pub(crate) fn set_pattern(&mut self, pat: &[u8]) {
        let qlen = pat.len();
        let mb = qlen.div_ceil(64).max(1);
        self.qlen = qlen;
        self.max_num_blocks = mb;
        self.w = mb * 64 - qlen;
        self.peq.clear();
        self.peq.resize(256 * mb, 0);
        for (r, &c) in pat.iter().enumerate() {
            self.peq[(c as usize) * mb + r / 64] |= 1u64 << (r % 64);
        }
        // The query is padded at the end with wildcard cells that match any
        // symbol (Edlib `buildPeq`), so a comparison never penalizes them.
        for r in qlen..mb * 64 {
            let bit = 1u64 << (r % 64);
            let off = (r / 64) as usize;
            for s in 0..256 {
                self.peq[s * mb + off] |= bit;
            }
        }
        self.p.clear();
        self.p.resize(mb, 0);
        self.m.clear();
        self.m.resize(mb, 0);
        self.score.clear();
        self.score.resize(mb, 0);
    }

    /// Banded Levenshtein (edit) distance between the current pattern and `b`,
    /// with `band` the maximum allowed distance. Returns `Some(d)` when
    /// `d <= band`, else `None` (the band was exceeded).
    pub(crate) fn distance_to(&mut self, b: &[u8], band: usize) -> Option<usize> {
        let qlen = self.qlen;
        let tlen = b.len();
        if qlen == 0 || tlen == 0 {
            let d = qlen.max(tlen);
            return (d <= band).then_some(d);
        }
        if qlen.abs_diff(tlen) > band {
            return None;
        }
        let mb = self.max_num_blocks;
        let w = self.w;
        let ql = qlen as i32;
        let tl = tlen as i32;
        let mut k = band.min(qlen.max(tlen)) as i32;

        // Initial Ukkonen band (Edlib `myersCalcEditDistanceNW`).
        let inner = k.min((k + ql - tl) / 2) + 1;
        let last0 = (inner as usize).div_ceil(64);
        let mut first_block: i32 = 0;
        let mut last_block: i32 = mb.min(last0) as i32 - 1;
        for blb in 0..=last_block as usize {
            self.p[blb] = !0u64;
            self.m[blb] = 0;
            self.score[blb] = ((blb + 1) * 64) as i32;
        }

        for (c, &ch) in b.iter().enumerate() {
            let base = (ch as usize) * mb;
            // Calculate the column, block by block (Edlib `calculateBlock`).
            let mut hout: i32 = 1;
            {
                let mut bl = first_block as usize;
                let mut bnd = first_block;
                while bnd <= last_block {
                    let (np, nm, nh) = calculate_block(self.p[bl], self.m[bl], self.peq[base + bl], hout);
                    self.p[bl] = np;
                    self.m[bl] = nm;
                    self.score[bl] += nh;
                    hout = nh;
                    bl += 1;
                    bnd += 1;
                }
            }
            let mut bl = last_block as usize;
            let ci = c as i32;

            // Tighten k with a valid lower bound on the remaining path.
            let term = (tl - ci - 1).max(ql - ((1 + last_block) * 64 - 1) - 1);
            let add = if last_block == mb as i32 - 1 { w as i32 } else { 0 };
            k = k.min(self.score[bl] + term + add);

            // Expand the band downward by one block if the border is still
            // within reach (Edlib NW).
            if last_block + 1 < mb as i32
                && !(((last_block + 1) * 64 - 1) > k - self.score[bl] + 2 * 64 - 2 - tl + ci + ql)
            {
                last_block += 1;
                bl += 1;
                self.p[bl] = !0u64;
                self.m[bl] = 0;
                let (np, nm, nh) = calculate_block(self.p[bl], self.m[bl], self.peq[base + bl], hout);
                self.p[bl] = np;
                self.m[bl] = nm;
                self.score[bl] = self.score[bl - 1] - hout + 64 + nh;
            }
            // Shrink the band from the bottom.
            while last_block >= first_block
                && (self.score[bl] >= k + 64
                    || ((last_block + 1) * 64 - 1 > k - self.score[bl] + 2 * 64 - 2 - tl + ci + ql + 1))
            {
                last_block -= 1;
                if last_block < first_block {
                    // Band empty: the distance exceeds the bound (and `bl` would underflow).
                    return None;
                }
                bl -= 1;
            }
            // Shrink the band from the top.
            while first_block <= last_block
                && (self.score[first_block as usize] >= k + 64
                    || ((first_block + 1) * 64 - 1
                        < self.score[first_block as usize] - k - tl + ql + ci))
            {
                first_block += 1;
            }
            if last_block < first_block {
                return None;
            }
        }

        if last_block == mb as i32 - 1 {
            // The bottom-right cell is the last query cell, at block index `w`.
            let (p, m, score) = (self.p[last_block as usize], self.m[last_block as usize], self.score[last_block as usize]);
            let mut sc = score;
            let mut mask = 1u64 << 63;
            let mut best = 0i32;
            for i in 0..64 {
                if i == w {
                    best = sc;
                    break;
                }
                if p & mask != 0 {
                    sc -= 1;
                }
                if m & mask != 0 {
                    sc += 1;
                }
                mask >>= 1;
            }
            if best <= k {
                return Some(best as usize);
            }
        }
        None
    }

    /// Convenience wrapper: build the profile for `a`, then compare with `b`.
    #[cfg(test)]
    pub(crate) fn distance(&mut self, a: &[u8], b: &[u8], band: usize) -> Option<usize> {
        self.set_pattern(a);
        self.distance_to(b, band)
    }
}

/// One 64-bit block of a Myers column (Edlib `calculateBlock`). `hin`/`hout`
/// are the horizontal deltas at the block boundaries (+1, -1 or 0).
#[inline]
fn calculate_block(pv: u64, mv: u64, eq: u64, hin: i32) -> (u64, u64, i32) {
    let hin_is_neg = ((hin >> 2) as u64) & 1; // 1 when hin == -1
    let xv = eq | mv;
    let eq = eq | hin_is_neg;
    let xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq;
    let ph = mv | !(xh | pv);
    let mh = pv & xh;
    let hout = ((ph >> 63) as i32) - ((mh >> 63) as i32);
    let ph = (ph << 1) | ((hin > 0) as u64);
    let mh = (mh << 1) | hin_is_neg;
    let pv_out = mh | !(xv | ph);
    let mv_out = ph & xv;
    (pv_out, mv_out, hout)
}

/// Convenience wrapper allocating fresh rolling rows (see [`Banded`]).
#[cfg(test)]
pub(crate) fn banded_edit_distance(a: &[u8], b: &[u8], band: usize) -> Option<usize> {
    Banded::new().distance(a, b, band)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Plain O(n·m) Levenshtein for the exhaustive band check.
    fn plain_levenshtein(a: &[u8], b: &[u8]) -> usize {
        let (n, m) = (a.len(), b.len());
        let mut prev: Vec<usize> = (0..=m).collect();
        for i in 1..=n {
            let mut cur = vec![0usize; m + 1];
            cur[0] = i;
            for j in 1..=m {
                let cost = usize::from(a[i - 1] != b[j - 1]);
                cur[j] = (prev[j] + 1).min(prev[j - 1] + cost).min(cur[j - 1] + 1);
            }
            prev = cur;
        }
        prev[m]
    }

    #[test]
    fn banded_edit_distance_cases() {
        // Identical.
        assert_eq!(banded_edit_distance(b"ACGTACGT", b"ACGTACGT", 0), Some(0));
        // Substitutions.
        assert_eq!(banded_edit_distance(b"ACGT", b"AGGT", 1), Some(1));
        assert_eq!(banded_edit_distance(b"ACGT", b"TGCA", 4), Some(4));
        // Indels, including leading deletions (P1 regression: ACGT -> CGT).
        assert_eq!(banded_edit_distance(b"ACGT", b"ACGTT", 1), Some(1));
        assert_eq!(banded_edit_distance(b"ACGTT", b"ACGT", 1), Some(1));
        assert_eq!(banded_edit_distance(b"ACGT", b"CGT", 1), Some(1));
        assert_eq!(banded_edit_distance(b"ACGT", b"GT", 2), Some(2));
        assert_eq!(banded_edit_distance(b"ACGT", b"AGT", 1), Some(1));
        // Empty strings (including leading deletions to empty).
        assert_eq!(banded_edit_distance(b"", b"", 0), Some(0));
        assert_eq!(banded_edit_distance(b"", b"ACG", 3), Some(3));
        assert_eq!(banded_edit_distance(b"ACG", b"", 3), Some(3));
        assert_eq!(banded_edit_distance(b"ACG", b"", 2), None);
        // Two differences: within band 2, not within band 1.
        assert_eq!(banded_edit_distance(b"ACGT", b"AGGA", 2), Some(2));
        assert_eq!(banded_edit_distance(b"ACGT", b"AGGA", 1), None);
        // Length difference beyond the band exits early.
        assert_eq!(banded_edit_distance(b"ACGT", b"ACGTTTTT", 2), None);
    }

    #[test]
    fn banded_edit_distance_matches_plain_and_is_symmetric() {
        // All binary strings of length 0..=5.
        let mut strs: Vec<Vec<u8>> = vec![Vec::new()];
        let mut level = vec![Vec::new()];
        for _ in 0..5 {
            let mut next = Vec::new();
            for s in &level {
                for c in [b'A', b'B'] {
                    let mut t = s.clone();
                    t.push(c);
                    next.push(t.clone());
                    strs.push(t);
                }
            }
            level = next;
        }
        let mut ed = Banded::new();
        for a in &strs {
            for b in &strs {
                let d = plain_levenshtein(a, b);
                for band in 0..=5 {
                    let want = (d <= band).then_some(d);
                    assert_eq!(ed.distance(a, b, band), want, "a={a:?} b={b:?} band={band} d={d}");
                    // Symmetry.
                    assert_eq!(ed.distance(b, a, band), want, "sym a={a:?} b={b:?} band={band}");
                }
            }
        }
    }

    /// Deterministic PRNG (LCG) sequence over {A,C,G,T}.
    fn lcg_seq(state: &mut u64, len: usize) -> Vec<u8> {
        (0..len)
            .map(|_| {
                *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[((*state >> 33) % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn banded_edit_distance_empty_band_returns_none() {
        // Regression: shrinking the band past block 0 underflowed the block cursor.
        let mut ed = Banded::new();
        ed.set_pattern(&[b'A'; 64]);
        assert_eq!(ed.distance_to(&[b'C'; 64], 0), None);
        assert_eq!(ed.distance_to(&[b'C'; 64], 63), None);
        assert_eq!(ed.distance_to(&[b'C'; 64], 64), Some(64));
        assert_eq!(ed.distance_to(&[b'A'; 64], 0), Some(0));
    }

    #[test]
    fn banded_edit_distance_matches_plain_multiblock() {
        // Exercise the multi-word (length > 64) Myers blocks, including the
        // Ukkonen band adjustment, against the plain DP: unrelated sequences,
        // mutated near-copies (the successful/merge path), and varied lengths.
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut cases: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        for &(la, lb) in &[
            (65usize, 65usize),
            (64, 70),
            (70, 64),
            (100, 100),
            (127, 128),
            (128, 65),
            (200, 190),
            (500, 500),
            (1000, 1000),
            (64, 1000),
        ] {
            let a = lcg_seq(&mut state, la);
            let b = lcg_seq(&mut state, lb);
            cases.push((a.clone(), b));
            // A near-copy: mutate `a` with substitutions and one indel.
            let mut m = a.clone();
            for _ in 0..(la / 20).max(1) {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                let pos = ((state >> 33) as usize) % m.len();
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                m[pos] = b"ACGT"[((state >> 33) % 4) as usize];
            }
            cases.push((a, m));
        }
        let mut ed = Banded::new();
        for (a, b) in &cases {
            let d = plain_levenshtein(a, b);
            // Probe bands around the true distance and the length difference.
            let mut bands: Vec<usize> = vec![0, 1, 2, d.saturating_sub(1), d, d + 1];
            bands.push(a.len().abs_diff(b.len()));
            bands.sort_unstable();
            bands.dedup();
            for band in bands {
                let want = (d <= band).then_some(d);
                assert_eq!(ed.distance(a, b, band), want, "a.len={} b.len={} d={d} band={band}", a.len(), b.len());
                assert_eq!(ed.distance(b, a, band), want, "sym a.len={} b.len={} d={d} band={band}", a.len(), b.len());
            }
        }
    }

    #[test]
    fn ln_gamma_known_values() {
        assert!((ln_gamma(1.0)).abs() < 1e-12);
        assert!((ln_gamma(2.0)).abs() < 1e-12);
        assert!((ln_gamma(5.0) - 24f64.ln()).abs() < 1e-12);
        assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-12);
        // ln Γ(101) = ln(100!)
        let lf: f64 = (1..=100).map(|i| (i as f64).ln()).sum();
        assert!((ln_gamma(101.0) - lf).abs() < 1e-9);
        assert!((ln_gamma(1e-3) - 6.907_178_885_383_854).abs() < 1e-9);
    }

    #[test]
    fn beta_binomial_reduces_to_binomial() {
        // rho -> 0: alpha, beta -> large; binomial(10, 0.1): P(X >= 3).
        let p = 0.1;
        let k = 1.0 / 1e-9 - 1.0;
        let want = 1.0 - (0.9f64.powi(10) + 10.0 * 0.1 * 0.9f64.powi(9) + 45.0 * 0.01 * 0.9f64.powi(8));
        assert!((want - 0.070_190_826_4).abs() < 1e-9);
        let got = bb_ln_sf(3, 10, p * k, (1.0 - p) * k).exp();
        assert!((got - want).abs() < 1e-6, "got {got} want {want}");
        // Lower-tail branch (x <= mean): P(X >= 1) = 1 - 0.9^10.
        let got = bb_ln_sf(1, 10, p * k, (1.0 - p) * k).exp();
        assert!((got - (1.0 - 0.9f64.powi(10))).abs() < 1e-6);
        // Far tail: P(X >= 10) = 1e-10.
        let got = bb_ln_sf(10, 10, p * k, (1.0 - p) * k).exp();
        assert!((got / 1e-10 - 1.0).abs() < 1e-4);
        assert_eq!(bb_ln_sf(0, 10, 1.0, 1.0), 0.0);
    }

    #[test]
    fn beta_binomial_uniform_case() {
        // a = b = 1: X is uniform on 0..=n, so P(X >= x) = (n - x + 1)/(n + 1).
        for x in 0..=20u32 {
            let got = bb_ln_sf(x, 20, 1.0, 1.0).exp();
            let want = (21 - x) as f64 / 21.0;
            assert!((got - want).abs() < 1e-10, "x={x} got {got} want {want}");
        }
        // Overdispersion makes the tail heavier than the binomial.
        let p = 0.01;
        let bin = bb_ln_sf(5, 100, p * 1e9, (1.0 - p) * 1e9);
        let k = 1.0 / 0.05 - 1.0;
        let bb = bb_ln_sf(5, 100, p * k, (1.0 - p) * k);
        assert!(bb > bin);
    }

    #[test]
    fn poisson_tail() {
        // mu = 1: P(N >= 3) = 1 - 2.5/e.
        let want = 1.0 - 2.5 / std::f64::consts::E;
        assert!((poisson_ln_sf(3, 1.0).exp() - want).abs() < 1e-12);
        // Lower-tail branch: P(N >= 1 | mu = 2) = 1 - e^-2.
        assert!((poisson_ln_sf(1, 2.0).exp() - (1.0 - (-2.0f64).exp())).abs() < 1e-12);
        // Far tail: P(N >= 10 | 0.1) ~ pmf(10) = e^-0.1 1e-10 / 10!.
        let pmf = (-0.1f64).exp() * 1e-10 / 3_628_800.0;
        assert!((poisson_ln_sf(10, 0.1).exp() / pmf - 1.0).abs() < 0.02);
        assert_eq!(poisson_ln_sf(0, 3.0), 0.0);
        // Large mu (references: independent log-sum-exp over math.lgamma terms
        // in Python). These underflowed / cancelled before.
        for (k, mu, want) in [
            (1000, 1000.0, -0.6847718632905715),
            (10_000, 10_000.0, -0.6904910943929234),
            (1100, 1000.0, -6.9458410143655005),
            (900, 1000.0, -0.0006227916790093602),
            (3, 1.0, -2.5219682600314),
            (50, 10.0, -43.13137931408243),
        ] {
            let got = poisson_ln_sf(k, mu);
            assert!((got - want).abs() < 1e-9 * want.abs().max(1.0), "k={k} mu={mu}: {got} vs {want}");
        }
        assert_eq!(poisson_ln_sf(2, 0.0), f64::NEG_INFINITY);
    }

    #[test]
    fn moment_fit() {
        // x = n p exactly: no overdispersion -> rho clamped to the minimum.
        let mut m = MomentSums::default();
        for _ in 0..100 {
            m.push(1, 10);
        }
        let f = fit_bb(&m);
        assert!((f.p - 0.1).abs() < 1e-12);
        assert_eq!(f.rho, RHO_MIN);
        // x in {3, 7}, n = 10: p = 0.5, rho = (4 - 2.5) / (2.5 * 9) = 1/15.
        let mut m = MomentSums::default();
        for _ in 0..50 {
            m.push(3, 10);
            m.push(7, 10);
        }
        let f = fit_bb(&m);
        assert!((f.p - 0.5).abs() < 1e-12);
        assert!((f.rho - 1.0 / 15.0).abs() < 1e-12);
        assert!((f.alpha - 0.5 * 14.0).abs() < 1e-9);
        assert!((f.beta - 0.5 * 14.0).abs() < 1e-9);
        // Extreme: {0, 10} -> rho = 1 -> clamped to 0.5.
        let mut m = MomentSums::default();
        m.push(0, 10);
        m.push(10, 10);
        assert_eq!(fit_bb(&m).rho, RHO_MAX);
        // Empty.
        let f = fit_bb(&MomentSums::default());
        assert_eq!(f.p, P_MIN);
    }

    #[test]
    fn wilson() {
        assert!((wilson_lower(8, 10) - 0.4902).abs() < 1e-3);
        assert!((wilson_lower(50, 100) - 0.4038).abs() < 1e-3);
        assert_eq!(wilson_lower(0, 0), 0.0);
        assert!(wilson_lower(0, 10) == 0.0);
    }

    #[test]
    fn benjamini_hochberg() {
        // Classic example: p = .01 .04 .03 .005, m = 4.
        let q = bh(&[0.01, 0.04, 0.03, 0.005], 4);
        let want = [0.02, 0.04, 0.04, 0.02];
        for (a, b) in q.iter().zip(want) {
            assert!((a - b).abs() < 1e-12, "{q:?}");
        }
        // Larger m (untested positions counted): q scales up.
        let q = bh(&[0.001], 1000);
        assert!((q[0] - 1.0).abs() < 1e-12);
        let q = bh(&[1e-6, 0.5], 100);
        assert!((q[0] - 1e-4).abs() < 1e-15);
        assert!((q[1] - 1.0).abs() < 1e-12);
    }
}
