//! `qc` site-calling rules, indels/background and candidate collection.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use crate::qc::*;
use crate::scan::{self, *};

fn base_idx(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Mode of the per-position depth over positions with depth > 0 (0 if none;
/// ties resolved to the smaller depth). Zero-depth positions are excluded on
/// purpose: a mostly uncovered contig still gets the depth of its covered part
/// as baseline (reported as `depth_mode_excludes_zero` and `covered_fraction`).
pub(crate) fn depth_mode(dp: &[u32]) -> u32 {
    let mut hist: Vec<u64> = Vec::new();
    for &d in dp {
        if d == 0 {
            continue;
        }
        let d = d as usize;
        if d >= hist.len() {
            hist.resize(d + 1, 0);
        }
        hist[d] += 1;
    }
    let mut best = 0usize;
    for (d, &c) in hist.iter().enumerate() {
        if c > hist[best] {
            best = d;
        }
    }
    best as u32
}

/// Sliding-window clusters over sorted, unique site positions (< `len`).
/// A window `[x, x + window)` qualifies if it holds at least `min_sites`
/// sites; the union of all qualifying windows is returned, overlapping (or
/// touching) windows merged. For `min_sites` consecutive sites `s_i ..= s_j`
/// with `s_j - s_i < window`, the windows holding them cover
/// `[s_j - window + 1, s_i + window)`. Linear contigs clip regions to
/// `[0, len)`; on circular contigs windows wrap across the origin and a
/// wrapping region is returned with `start < len < end`. Returns
/// `(start, end, sites in region)`.
pub(crate) fn cluster_windows(
    sites: &[u32],
    len: u64,
    circular: bool,
    window: u64,
    min_sites: usize,
) -> Vec<(u64, u64, usize)> {
    let min_sites = min_sites.max(1);
    if sites.len() < min_sites || len == 0 || window == 0 {
        return Vec::new();
    }
    let (l, w) = (len as i64, window as i64);
    // Sites of one turn plus (circular) the next turn, so runs can wrap.
    let mut ext: Vec<i64> = sites.iter().map(|&s| s as i64).collect();
    if circular {
        ext.extend(sites.iter().map(|&s| s as i64 + l));
    }
    let mut ivs: Vec<(i64, i64)> = Vec::new();
    for i in 0..sites.len() {
        let j = i + min_sites - 1;
        if j < ext.len() && ext[j] - ext[i] < w {
            let (mut a, mut b) = (ext[j] - w + 1, ext[i] + w);
            if circular {
                // Normalise so the start lies in [0, len).
                let sh = a.div_euclid(l) * l;
                a -= sh;
                b -= sh;
            } else {
                a = a.max(0);
                b = b.min(l);
            }
            ivs.push((a, b));
        }
    }
    ivs.sort_unstable();
    let mut regs: Vec<(i64, i64)> = Vec::new();
    for (a, b) in ivs {
        match regs.last_mut() {
            Some(last) if a <= last.1 => last.1 = last.1.max(b),
            _ => regs.push((a, b)),
        }
    }
    if circular {
        // A wrapping last region may overlap the first ones.
        while regs.len() >= 2 {
            let last = *regs.last().unwrap();
            if last.1 > l && regs[0].0 <= last.1 - l {
                let first = regs.remove(0);
                let lr = regs.last_mut().unwrap();
                lr.1 = lr.1.max(first.1 + l);
            } else {
                break;
            }
        }
        for r in &mut regs {
            if r.1 - r.0 >= l {
                *r = (0, l);
            }
        }
    }
    let mut all: Vec<i64> = sites.iter().map(|&s| s as i64).collect();
    if circular {
        all.extend(sites.iter().map(|&s| s as i64 + l));
    }
    regs.into_iter()
        .map(|(a, b)| {
            let n = all.iter().filter(|&&x| x >= a && x < b).count().min(sites.len());
            (a as u64, b as u64, n)
        })
        .collect()
}

/// Does a SNP (alt base with `altf`/`altr` reads at depth `dp`) qualify as a
/// minority-SNP cluster site? Base rule: alt >= cluster_min_alt and AF >=
/// cluster_min_af; plus the switchable refinements A (nominal beta-binomial
/// p against the SNP background) and C (strand).
pub(crate) fn cluster_site_ok(altf: u32, altr: u32, dp: u32, dpf: u32, dpr: u32, fit: &BbFit, p: &QcParams) -> bool {
    let alt = altf + altr;
    if dp == 0 || alt < p.cluster_min_alt || (alt as f64 / dp as f64) < p.cluster_min_af {
        return false;
    }
    let strand = match p.cluster_strand {
        ClusterStrand::Off => true,
        ClusterStrand::Both => altf > 0 && altr > 0,
        ClusterStrand::Ratio => strand_ok(altf, altr, dpf, dpr, p.strand_min_dp, p.strand_ratio),
    };
    if !strand {
        return false;
    }
    p.cluster_site_p >= 1.0 || stats::bb_ln_sf(alt, dp, fit.alpha, fit.beta) <= p.cluster_site_p.ln()
}

/// Is a SNP a cluster site? Sites below the absolute low-coverage threshold
/// (`dp < low_abs`) are excluded; otherwise the usual cluster-site rule.
#[allow(clippy::too_many_arguments)]
pub(crate) fn lowcov_cluster_ok(
    altf: u32,
    altr: u32,
    dp: u32,
    dpf: u32,
    dpr: u32,
    low_abs: u32,
    fit: &BbFit,
    p: &QcParams,
) -> bool {
    dp >= low_abs && cluster_site_ok(altf, altr, dp, dpf, dpr, fit, p)
}

/// Number of sliding-window hypotheses on one contig, as `cluster_windows`
/// tests them: linear `len - W + 1` (at least 1), circular `len` (every start
/// position, wrapping).
pub(crate) fn window_tests(len: u64, circular: bool, window: u64) -> u64 {
    if circular { len.max(1) } else { len.saturating_sub(window).saturating_add(1).max(1) }
}

/// Refinement B, per contig: the smallest site count k such that a window
/// of `window` bp with k sites is significant, P(N >= k | λ·W) <= alpha /
/// n_windows, with λ = `sites / assessable` (this contig) and n_windows =
/// [`window_tests`] for this contig (Bonferroni over the window hypotheses of
/// the same contig λ is estimated on). Returns `usize::MAX` if unattainable
/// (no assessable bases, or more sites needed than the contig has or a window
/// can hold).
pub(crate) fn density_min_sites(
    sites: u64,
    assessable: u64,
    len: u64,
    circular: bool,
    window: u64,
    alpha: f64,
) -> usize {
    if assessable == 0 || window == 0 || sites == 0 {
        return usize::MAX;
    }
    let mu = sites as f64 / assessable as f64 * window as f64;
    let thr = (alpha / window_tests(len, circular, window) as f64).ln();
    let sig = |k: u64| stats::poisson_ln_sf(k, mu) <= thr;
    // The tail is monotone in k: binary search on [1, cap].
    let cap = sites.min(window);
    if !sig(cap) {
        return usize::MAX;
    }
    let (mut lo, mut hi) = (1u64, cap);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if sig(mid) { hi = mid } else { lo = mid + 1 }
    }
    lo as usize
}

pub(crate) fn contig_ctx(c: &ContigCounts, refseq: &[u8], circular: bool, p: &QcParams, low_abs: u32, low_rel: f64) -> ContigCtx {
    let (hp_run, hp_fwd) = scan::hp_runs(refseq);
    let str_span = scan::str_spans(refseq);
    let dp: Vec<u32> = (0..c.len).map(|i| c.dp(i)).collect();
    let na: Vec<u8> = (0..c.len)
        .map(|i| {
            if base_idx(refseq[i]).is_none() {
                1
            } else if dp[i] < p.min_depth || dp[i] == 0 {
                2
            } else if c.cov[i][MQ0] as f64 / dp[i] as f64 > p.max_mq0_frac {
                3
            } else {
                0
            }
        })
        .collect();
    let mode = depth_mode(&dp);
    let lowcov: Vec<u8> = (0..c.len)
        .map(|i| if na[i] == 0 { lowcov_bits(dp[i], mode, low_abs, low_rel) } else { 0 })
        .collect();
    let mean = if c.len == 0 { 0.0 } else { dp.iter().map(|&d| d as u64).sum::<u64>() as f64 / c.len as f64 };
    let depth_skipped = if mode < low_abs { Some("low_coverage") } else { None };
    ContigCtx {
        hp_run,
        hp_fwd,
        str_span,
        dp,
        na,
        lowcov,
        circular,
        // Overridden by the caller with the wrap decision source.
        circular_by: "verdict",
        mode,
        mean,
        depth_skipped,
    }
}

/// Low-coverage annotation regions from the per-position reason bits: runs of
/// low-coverage positions with gaps < `merge` merged; regions shorter than
/// `min_len` dropped. Returns `(start, end, reason bits union)`.
pub(crate) fn lowcov_regions(mask: &[u8], merge: u32, min_len: u32) -> Vec<(u64, u64, u8)> {
    let n = mask.len();
    let mut out = Vec::new();
    let mut i = 0;
    while i < n {
        if mask[i] == 0 {
            i += 1;
            continue;
        }
        let start = i;
        let mut end;
        let mut reason = 0u8;
        let mut j = i;
        loop {
            while j < n && mask[j] != 0 {
                reason |= mask[j];
                j += 1;
            }
            end = j;
            let gap_start = j;
            while j < n && mask[j] == 0 {
                j += 1;
            }
            if j < n && ((j - gap_start) as u32) < merge {
                continue;
            }
            break;
        }
        if (end - start) as u32 >= min_len {
            out.push((start as u64, end as u64, reason));
        }
        i = j;
    }
    out
}

pub(crate) fn low_coverage_regions(tid: usize, c: &ContigCounts, ctx: &ContigCtx, p: &QcParams, out: &mut Vec<Region>) {
    for (s, e, reason) in lowcov_regions(&ctx.lowcov, p.lowcov_merge, p.lowcov_min_len) {
        let (dr, hq) = region_depth(c, ctx, s, e);
        out.push(Region {
            tid,
            start: s,
            end: e,
            label: "low_coverage",
            score: 0,
            tier: "none",
            reason: lowcov_reason(reason),
            depth_ratio: dr,
            hiq_share: hq,
            extra: String::new(),
            lowcov: false,
            link: None,
            clip_left: None,
            clip_right: None,
        });
    }
}

/// The one class-context rule for indels: (hp_fwd, str) where the indel
/// lives — anchor + 1 for ordinary alleles, position 0 for origin alleles
/// (which sit before the first base). Used for background, BH m and tests.
#[inline]
fn indel_ctx(ctx: &ContigCtx, anchor: usize, origin: bool) -> (u8, u16) {
    let p = if origin { 0 } else { anchor + 1 };
    if p < ctx.hp_fwd.len() { (ctx.hp_fwd[p], ctx.str_span[p]) } else { (1, 0) }
}

/// One grouped INS candidate: the representative exact allele plus the support
/// of every similar allele merged into it (same anchor window, length within
/// [`QcParams::ins_group_len_frac`] and sequence within
/// [`QcParams::ins_group_max_diff`] edit distance).
#[derive(Clone, Debug)]
pub(crate) struct InsGroup {
    /// Representative anchor (0 for origin alleles).
    pub anchor: u32,
    pub origin: bool,
    pub len: u32,
    /// Representative inserted sequence (exact).
    pub seq: Vec<u8>,
    /// Total fwd/rev support over the group (summed, per the rule).
    pub fwd: u32,
    pub rev: u32,
    /// Exact alleles merged into the group (1 = not grouped).
    pub n: u32,
}

/// Greedy grouping of similar INS alleles of one contig (see `InsGroup`).
/// Only alleles at the *same anchor* are grouped: they come from distinct reads
/// that all cover the anchor, so the group's alt can never exceed that anchor's
/// DP. Alleles shorter than `ins_group_min_len` (or with `ins_group_min_len ==
/// 0`) are not grouped and keep the exact-sequence behaviour. Grouping is a
/// star: the highest-support ungrouped allele of an anchor becomes a
/// representative and absorbs every other ungrouped allele of that anchor
/// matching it; ties are broken by sequence, for determinism.
pub(crate) fn ins_groups(indels: &[IndelAllele], p: &QcParams) -> Vec<InsGroup> {
    if p.ins_group_min_len == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut ed = stats::Banded::new();
    // Indels are sorted by anchor; group the long INS alleles of each anchor.
    let mut i = 0;
    while i < indels.len() {
        let anchor = indels[i].anchor;
        let mut j = i;
        while j < indels.len() && indels[j].anchor == anchor {
            j += 1;
        }
        let mut idx: Vec<usize> =
            (i..j).filter(|&k| indels[k].ins && indels[k].len >= p.ins_group_min_len).collect();
        if !idx.is_empty() {
            idx.sort_by(|&x, &y| {
                let (a, b) = (&indels[x], &indels[y]);
                (b.fwd + b.rev).cmp(&(a.fwd + a.rev)).then(a.seq.cmp(&b.seq))
            });
            let mut used = vec![false; idx.len()];
            for gi in 0..idx.len() {
                if used[gi] {
                    continue;
                }
                used[gi] = true;
                let rep = &indels[idx[gi]];
                let tol = ((p.ins_group_len_frac * rep.len as f64) as u32).max(2);
                let band = (p.ins_group_max_diff * rep.len as f64) as usize;
                // Build the query profile once; every candidate is compared
                // against this representative.
                ed.set_pattern(&rep.seq);
                let mut g =
                    InsGroup { anchor: rep.anchor, origin: rep.origin, len: rep.len, seq: rep.seq.clone(), fwd: rep.fwd, rev: rep.rev, n: 1 };
                for gj in (gi + 1)..idx.len() {
                    if used[gj] {
                        continue;
                    }
                    let m = &indels[idx[gj]];
                    // Cheap prefilters before the DP: same strand/origin, and
                    // a length difference already beyond the length tolerance.
                    if m.origin != rep.origin || m.len.abs_diff(rep.len) > tol {
                        continue;
                    }
                    if ed.distance_to(&m.seq, band).is_none() {
                        continue;
                    }
                    g.fwd += m.fwd;
                    g.rev += m.rev;
                    g.n += 1;
                    used[gj] = true;
                }
                out.push(g);
            }
        }
        i = j;
    }
    out
}

/// Per-anchor fwd/rev totals over all indel alleles (used for the ref counts
/// `DP − all indel alt` at the anchor).
fn anchor_sums(indels: &[IndelAllele]) -> std::collections::HashMap<u32, (u32, u32)> {
    let mut m: std::collections::HashMap<u32, (u32, u32)> = std::collections::HashMap::new();
    for a in indels {
        let e = m.entry(a.anchor).or_insert((0, 0));
        e.0 += a.fwd;
        e.1 += a.rev;
    }
    m
}

/// Per-position max INS and DEL allele counts (by anchor), ordinary alleles
/// only, plus the max INS and DEL counts of origin alleles (a separate test
/// slot at position 0). Long INS alleles are replaced by their group totals at
/// the representative anchor, so the background fit uses the same statistic
/// that is tested.
fn indel_max(len: usize, indels: &[IndelAllele], groups: &[InsGroup], min_len: u32) -> (Vec<u32>, Vec<u32>, [u32; 2]) {
    let mut ins = vec![0u32; len];
    let mut del = vec![0u32; len];
    let mut origin = [0u32; 2];
    for a in indels {
        if min_len > 0 && a.ins && a.len >= min_len {
            continue; // folded into a group
        }
        let n = a.fwd + a.rev;
        let s = if a.origin {
            &mut origin[!a.ins as usize]
        } else if a.ins {
            &mut ins[a.anchor as usize]
        } else {
            &mut del[a.anchor as usize]
        };
        *s = (*s).max(n);
    }
    for g in groups {
        let n = g.fwd + g.rev;
        let s = if g.origin { &mut origin[0] } else { &mut ins[g.anchor as usize] };
        *s = (*s).max(n);
    }
    (ins, del, origin)
}

pub(crate) fn accumulate_background(
    c: &ContigCounts,
    refseq: &[u8],
    ctx: &ContigCtx,
    indels: &[IndelAllele],
    groups: &[InsGroup],
    p: &QcParams,
    sums: &mut [MomentSums; NCLASS],
    m: &mut [u64; NCLASS],
) {
    let (ins_max, del_max, origin_max) = indel_max(c.len, indels, groups, p.ins_group_min_len);
    let lim = p.fit_max_af;
    for i in 0..c.len {
        if ctx.na[i] != 0 {
            continue;
        }
        let n = ctx.dp[i];
        let r = base_idx(refseq[i]).unwrap();
        let b = &c.bases[i];
        let mut x = 0u32;
        for k in 0..4 {
            if k != r {
                x = x.max(b[k] + b[4 + k]);
            }
        }
        let nf = n as f64;
        m[0] += 1;
        if (x as f64) / nf < lim {
            sums[0].push(x, n);
        }
        // Ordinary slot at every position; at position 0 also the origin slot.
        let slots: &[(bool, u32, u32)] =
            if i == 0 { &[(false, ins_max[i], del_max[i]), (true, origin_max[0], origin_max[1])] } else { &[(false, ins_max[i], del_max[i])] };
        for &(origin, xi, xd) in slots {
            let (hp, st) = indel_ctx(ctx, i, origin);
            let bk = indel_bucket(hp, st);
            let (ci, cdl) = (1 + bk, 1 + NB + bk);
            m[ci] += 1;
            m[cdl] += 1;
            if (xi as f64) / nf < lim {
                sums[ci].push(xi, n);
            }
            if (xd as f64) / nf < lim {
                sums[cdl].push(xd, n);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn collect_candidates(
    tid: usize,
    c: &ContigCounts,
    refseq: &[u8],
    ctx: &ContigCtx,
    indels: &[IndelAllele],
    groups: &[InsGroup],
    p: &QcParams,
    indel_min_af: f64,
    low_abs: u32,
    snp_fit: &BbFit,
    out: &mut Vec<Candidate>,
    cluster_sites: &mut Vec<u32>,
) {
    // SNPs.
    for i in 0..c.len {
        if ctx.na[i] != 0 {
            continue;
        }
        let r = base_idx(refseq[i]).unwrap();
        let b = &c.bases[i];
        let mut best = usize::MAX;
        let mut bx = 0u32;
        for k in 0..4 {
            if k != r && b[k] + b[4 + k] > bx {
                bx = b[k] + b[4 + k];
                best = k;
            }
        }
        if best == usize::MAX {
            continue;
        }
        let dp = ctx.dp[i];
        let af = bx as f64 / dp as f64;
        // Sites below the absolute low-coverage threshold are not cluster sites.
        if lowcov_cluster_ok(
            b[best],
            b[4 + best],
            dp,
            c.cov[i][scan::DPF],
            c.cov[i][scan::DPR],
            low_abs,
            snp_fit,
            p,
        ) {
            cluster_sites.push(i as u32);
        }
        if bx < p.min_alt || af < p.snp_min_af {
            continue;
        }
        let cov = &c.cov[i];
        let (dpf, dpr) = (cov[scan::DPF], cov[scan::DPR]);
        out.push(Candidate {
            tid,
            pos0: i as u32,
            vtype: VarType::Snp,
            len: 1,
            origin: false,
            ref_allele: vec![refseq[i]],
            alt_allele: vec![b"ACGT"[best]],
            dp,
            dpf,
            dpr,
            altf: b[best],
            altr: b[4 + best],
            reff: b[r],
            refr: b[4 + r],
            hp: ctx.hp_run[i] as u32,
            str_span: ctx.str_span[i] as u32,
            mq0f: cov[MQ0] as f64 / dp as f64,
            class: 0,
            grouped: 1,
            ln_p: 0.0,
            q: 1.0,
            wilson: stats::wilson_lower(bx, dp),
            strand_ok: strand_ok(b[best], b[4 + best], dpf, dpr, p.strand_min_dp, p.strand_ratio),
            significant: false,
            lowcov: ctx.lowcov[i],
            linked_group: None,
            unlinked: false,
            tier: TierOut { tier: Tier::None, reason: "", mixed: false, strand: false },
        });
    }
    // Indels, grouped by anchor. Long INS alleles are emitted from their group
    // below; short INS and all DEL keep the per-allele path.
    let sums = anchor_sums(indels);
    let mut i = 0;
    while i < indels.len() {
        let anchor = indels[i].anchor as usize;
        let mut j = i;
        while j < indels.len() && indels[j].anchor as usize == anchor {
            j += 1;
        }
        if ctx.na[anchor] == 0 {
            let (sf, sr) = sums[&(anchor as u32)];
            let dp = ctx.dp[anchor];
            let cov = &c.cov[anchor];
            let (dpf, dpr) = (cov[scan::DPF], cov[scan::DPR]);
            for a in &indels[i..j] {
                if p.ins_group_min_len > 0 && a.ins && a.len >= p.ins_group_min_len {
                    continue; // folded into a group
                }
                let alt = a.fwd + a.rev;
                let af = alt as f64 / dp as f64;
                if alt < p.min_alt || af < indel_min_af {
                    continue;
                }
                // Origin alleles live at position 0 (context there, not at anchor + 1).
                let (hp, st) = indel_ctx(ctx, anchor, a.origin);
                let bk = indel_bucket(hp, st);
                let (vtype, class, ref_allele, alt_allele) = if a.origin {
                    // Right-padded VCF form: POS = 1, anchor base appended.
                    if a.ins {
                        let mut alt_a = a.seq.clone();
                        alt_a.push(refseq[0]);
                        (VarType::Ins, 1 + bk, vec![refseq[0]], alt_a)
                    } else {
                        let l = a.len as usize;
                        (VarType::Del, 1 + NB + bk, refseq[..=l].to_vec(), vec![refseq[l]])
                    }
                } else if a.ins {
                    let mut alt_a = vec![refseq[anchor]];
                    alt_a.extend_from_slice(&a.seq);
                    (VarType::Ins, 1 + bk, vec![refseq[anchor]], alt_a)
                } else {
                    let end = (anchor + 1 + a.len as usize).min(refseq.len());
                    (VarType::Del, 1 + NB + bk, refseq[anchor..end].to_vec(), vec![refseq[anchor]])
                };
                out.push(Candidate {
                    tid,
                    pos0: anchor as u32,
                    origin: a.origin,
                    vtype,
                    len: a.len,
                    ref_allele,
                    alt_allele,
                    dp,
                    dpf,
                    dpr,
                    altf: a.fwd,
                    altr: a.rev,
                    reff: dpf.saturating_sub(sf),
                    refr: dpr.saturating_sub(sr),
                    hp: hp as u32,
                    str_span: st as u32,
                    mq0f: cov[MQ0] as f64 / dp as f64,
                    class,
                    grouped: 1,
                    ln_p: 0.0,
                    q: 1.0,
                    wilson: stats::wilson_lower(alt, dp),
                    strand_ok: strand_ok(a.fwd, a.rev, dpf, dpr, p.strand_min_dp, p.strand_ratio),
                    significant: false,
                    lowcov: ctx.lowcov[anchor],
                    linked_group: None,
                    unlinked: false,
                    tier: TierOut { tier: Tier::None, reason: "", mixed: false, strand: false },
                });
            }
        }
        i = j;
    }
    // Grouped long INS alleles: one candidate at the anchor with the group's
    // summed support. Group members are same-anchor alleles from distinct reads
    // that all cover the anchor, so alt stays within that anchor's DP (REF =
    // DP − all indel alleles at the anchor).
    for g in groups {
        let anchor = g.anchor as usize;
        if ctx.na[anchor] != 0 {
            continue;
        }
        let alt = g.fwd + g.rev;
        let dp = ctx.dp[anchor];
        let af = alt as f64 / dp as f64;
        if alt < p.min_alt || af < indel_min_af {
            continue;
        }
        let cov = &c.cov[anchor];
        let (dpf, dpr) = (cov[scan::DPF], cov[scan::DPR]);
        let (hp, st) = indel_ctx(ctx, anchor, g.origin);
        let bk = indel_bucket(hp, st);
        let (ref_allele, alt_allele) = if g.origin {
            let mut alt_a = g.seq.clone();
            alt_a.push(refseq[0]);
            (vec![refseq[0]], alt_a)
        } else {
            let mut alt_a = vec![refseq[anchor]];
            alt_a.extend_from_slice(&g.seq);
            (vec![refseq[anchor]], alt_a)
        };
        let (sf, sr) = sums[&(anchor as u32)];
        out.push(Candidate {
            tid,
            pos0: anchor as u32,
            origin: g.origin,
            vtype: VarType::Ins,
            len: g.len,
            ref_allele,
            alt_allele,
            dp,
            dpf,
            dpr,
            altf: g.fwd,
            altr: g.rev,
            reff: dpf.saturating_sub(sf),
            refr: dpr.saturating_sub(sr),
            hp: hp as u32,
            str_span: st as u32,
            mq0f: cov[MQ0] as f64 / dp as f64,
            class: 1 + bk,
            grouped: g.n,
            ln_p: 0.0,
            q: 1.0,
            wilson: stats::wilson_lower(alt, dp),
            strand_ok: strand_ok(g.fwd, g.rev, dpf, dpr, p.strand_min_dp, p.strand_ratio),
            significant: false,
            lowcov: ctx.lowcov[anchor],
            linked_group: None,
            unlinked: false,
            tier: TierOut { tier: Tier::None, reason: "", mixed: false, strand: false },
        });
    }
}

/// Mean depth and high-MAPQ share over `[s, e)` (positions taken modulo the
/// contig length, for wrapping regions).
pub(crate) fn region_depth(c: &ContigCounts, ctx: &ContigCtx, s: u64, e: u64) -> (Option<f64>, Option<f64>) {
    let len = c.len as u64;
    if e <= s || len == 0 {
        return (None, None);
    }
    let (mut sd, mut sh) = (0u64, 0u64);
    for x in s..e {
        let i = (x % len) as usize;
        sd += ctx.dp[i] as u64;
        sh += c.cov[i][HIQ] as u64;
    }
    let mean = sd as f64 / (e - s) as f64;
    let ratio = if ctx.mode > 0 { Some(mean / ctx.mode as f64) } else { None };
    let hiq = if sd > 0 { Some(sh as f64 / sd as f64) } else { None };
    (ratio, hiq)
}


#[cfg(test)]
mod tests {
    use super::*;

    fn ins_allele(anchor: u32, seq: &[u8], fwd: u32, rev: u32) -> IndelAllele {
        IndelAllele { tid: 0, anchor, origin: false, ins: true, len: seq.len() as u32, seq: seq.to_vec(), fwd, rev }
    }
    fn flip(seq: &[u8], pos: &[usize]) -> Vec<u8> {
        let mut s = seq.to_vec();
        for &i in pos {
            s[i] = match s[i] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        s
    }
    #[test]
    fn mode_of_depth() {
        assert_eq!(depth_mode(&[0, 0, 0, 5, 5, 6, 7, 7]), 5);
        assert_eq!(depth_mode(&[0, 0]), 0);
    }
    #[test]
    fn sliding_window_clusters() {
        // Five sites within 500 bp: union of the windows holding all five.
        let s = [100, 200, 300, 400, 599];
        assert_eq!(cluster_windows(&s, 10_000, false, 500, 5), vec![(100, 600, 5)]);
        // Span of exactly 500 is too wide (a window holds positions x..x+499).
        let s = [100, 200, 300, 400, 600];
        assert!(cluster_windows(&s, 10_000, false, 500, 5).is_empty());
        // Overlapping windows merge; linear contigs clip at 0.
        let s = [100, 150, 200, 250, 300, 500, 700];
        assert_eq!(cluster_windows(&s, 10_000, false, 500, 5), vec![(0, 650, 6)]);
        // Two separate clusters.
        let s = [0, 10, 20, 30, 40, 5000, 5010, 5020, 5030, 5040];
        assert_eq!(cluster_windows(&s, 10_000, false, 500, 5), vec![(0, 500, 5), (4541, 5500, 5)]);
    }
    #[test]
    fn cluster_windows_merge_windows_not_site_spans() {
        // Reviewer counterexample: the site spans 10..15 and 23..28 do not
        // overlap, but the qualifying windows [5,20) and [18,33) do.
        let s = [10, 11, 12, 13, 14, 23, 24, 25, 26, 27];
        assert_eq!(cluster_windows(&s, 100, false, 10, 5), vec![(5, 33, 10)]);
        // Linear clipping at both contig ends.
        assert_eq!(cluster_windows(&[0, 1, 2, 3, 4], 100, false, 10, 5), vec![(0, 10, 5)]);
        assert_eq!(cluster_windows(&[95, 96, 97, 98, 99], 100, false, 10, 5), vec![(90, 100, 5)]);
        // Circular: the same sites wrap instead of being clipped.
        assert_eq!(cluster_windows(&[95, 96, 97, 98, 99], 100, true, 10, 5), vec![(90, 105, 5)]);
        assert_eq!(cluster_windows(&[0, 1, 2, 3, 4], 100, true, 10, 5), vec![(95, 110, 5)]);
    }
    #[test]
    fn cluster_site_refinements() {
        let fit = BbFit { p: 0.01, rho: 0.01, alpha: 0.01 * 99.0, beta: 0.99 * 99.0 };
        // All refinements off: the old base rule.
        let base = QcParams { cluster_site_p: 1.0, cluster_strand: ClusterStrand::Off, ..QcParams::default() };
        // Base rule only: 2 reads at 12x (AF 0.17) qualify; 1 read or AF < 0.10 do not.
        assert!(cluster_site_ok(1, 1, 12, 6, 6, &fit, &base));
        assert!(!cluster_site_ok(1, 0, 12, 6, 6, &fit, &base));
        assert!(!cluster_site_ok(5, 5, 200, 100, 100, &fit, &base));
        // A: nominal p <= 1e-3 rejects 2/12 against a 1 % background, keeps 20/100.
        let a = QcParams { cluster_site_p: 1e-3, cluster_strand: ClusterStrand::Off, ..QcParams::default() };
        let p2 = stats::bb_ln_sf(2, 12, fit.alpha, fit.beta).exp();
        assert!(p2 > 1e-3, "{p2}");
        assert!(!cluster_site_ok(1, 1, 12, 6, 6, &fit, &a));
        assert!(cluster_site_ok(10, 10, 100, 50, 50, &fit, &a));
        // C both: alt on both strands.
        let cb = QcParams { cluster_strand: ClusterStrand::Both, cluster_site_p: 1.0, ..QcParams::default() };
        assert!(!cluster_site_ok(20, 0, 100, 50, 50, &fit, &cb));
        assert!(cluster_site_ok(19, 1, 100, 50, 50, &fit, &cb));
        // C ratio: the site strand gate (19 vs 1 fails, 12 vs 8 passes).
        let cr = QcParams { cluster_strand: ClusterStrand::Ratio, cluster_site_p: 1.0, ..QcParams::default() };
        // Defaults are A + C(ratio), density off.
        let d = QcParams::default();
        assert_eq!((d.cluster_site_p, d.cluster_strand, d.cluster_density_p), (1e-3, ClusterStrand::Ratio, None));
        assert!(!cluster_site_ok(1, 1, 12, 6, 6, &fit, &d));
        assert!(!cluster_site_ok(20, 0, 100, 50, 50, &fit, &d));
        assert!(cluster_site_ok(12, 8, 100, 50, 50, &fit, &d));
        assert!(!cluster_site_ok(19, 1, 100, 50, 50, &fit, &cr));
        assert!(cluster_site_ok(12, 8, 100, 50, 50, &fit, &cr));
    }
    #[test]
    fn cluster_density_threshold() {
        // Expected k values from an independent Python computation
        // (log-sum-exp of math.lgamma Poisson terms, linear k scan).
        // Window hypotheses: linear len - W + 1, circular len.
        assert_eq!(window_tests(1_000_000, false, 500), 999_501);
        assert_eq!(window_tests(1_000_000, true, 500), 1_000_000);
        assert_eq!(window_tests(300, false, 500), 1);
        // 100 sites on 1 Mb, W = 500: mu = 0.05, threshold 0.01 / 999501:
        // P(N >= 4) = 2.6e-7 is not enough, P(N >= 5) = 3e-9 is -> k = 5.
        assert_eq!(density_min_sites(100, 1_000_000, 1_000_000, false, 500, 0.01), 5);
        assert_eq!(density_min_sites(100, 1_000_000, 1_000_000, true, 500, 0.01), 5);
        // Dense noise (mu = 10).
        let k = density_min_sites(20_000, 1_000_000, 1_000_000, false, 500, 0.01);
        assert_eq!(k, 33);
        assert!(density_min_sites(20_000, 1_000_000, 1_000_000, false, 500, 1e-6) >= k);
        // Linear vs circular window counts change k on a small contig (101 vs 2000 tests).
        assert_eq!(density_min_sites(20, 2_000, 2_000, false, 500, 0.01), 18);
        assert_eq!(density_min_sites(20, 2_000, 2_000, true, 500, 0.01), 19);
        // Large mu (5000): binary search, fast and minimal.
        let t = std::time::Instant::now();
        assert_eq!(density_min_sites(5_000_000, 10_000_000, 10_000_000, false, 10_000, 0.01), 5431);
        assert!(t.elapsed().as_secs_f64() < 1.0);
        // Unattainable: more sites needed than exist.
        assert_eq!(density_min_sites(3, 20_000, 20_000, true, 500, 1e-9), usize::MAX);
        assert_eq!(density_min_sites(5, 0, 1_000, false, 500, 0.01), usize::MAX);
        // Applied through cluster_windows: 5 sites in 500 bp survive at k = 5,
        // not at k = 33.
        let s = [100, 200, 300, 400, 599];
        assert_eq!(cluster_windows(&s, 1_000_000, false, 500, 5).len(), 1);
        assert!(cluster_windows(&s, 1_000_000, false, 500, k).is_empty());
    }
    #[test]
    fn sliding_window_wraps_origin() {
        // 3 sites before the end, 2 after the origin: only a circular window finds them.
        let s = [50, 120, 9_800, 9_900, 9_950];
        assert!(cluster_windows(&s, 10_000, false, 500, 5).is_empty());
        assert_eq!(cluster_windows(&s, 10_000, true, 500, 5), vec![(9_621, 10_300, 5)]);
        // A cluster at the origin plus a wrapping one merge into one region.
        let s = [10, 20, 30, 40, 50, 9_990];
        assert_eq!(cluster_windows(&s, 10_000, true, 500, 5), vec![(9_541, 10_510, 6)]);
    }
    #[test]
    fn lowcov_cluster_sites_excluded() {
        let fit = BbFit { p: 0.01, rho: 0.01, alpha: 0.01 * 99.0, beta: 0.99 * 99.0 };
        let p = QcParams { cluster_site_p: 1.0, cluster_strand: ClusterStrand::Off, ..QcParams::default() };
        // 20/100 qualifies as a cluster site at these settings.
        assert!(lowcov_cluster_ok(10, 10, 100, 50, 50, 20, &fit, &p));
        // But not when the site is below the absolute low-coverage threshold.
        assert!(!lowcov_cluster_ok(10, 10, 19, 10, 9, 20, &fit, &p));
    }
    #[test]
    fn lowcov_region_merge_and_min_len() {
        // [0,100) and [130,230): gap 30 < 50 merged into one 230 bp region,
        // reason bits OR to 3.
        let mut m = vec![0u8; 300];
        for x in m[0..100].iter_mut() {
            *x = 1;
        }
        for x in m[130..230].iter_mut() {
            *x = 2;
        }
        assert_eq!(lowcov_regions(&m, 50, 100), vec![(0, 230, 3)]);
        // A gap of exactly 50 is not merged (< 50 required): the 230 bp piece
        // survives, the 100 bp piece survives too (>= min_len).
        let mut m = vec![0u8; 400];
        for x in m[0..100].iter_mut() {
            *x = 1;
        }
        for x in m[150..250].iter_mut() {
            *x = 1;
        }
        assert_eq!(lowcov_regions(&m, 50, 100), vec![(0, 100, 1), (150, 250, 1)]);
        // Regions shorter than min_len are dropped (gap >= merge keeps them apart).
        let mut m = vec![0u8; 400];
        for x in m[0..60].iter_mut() {
            *x = 1;
        }
        for x in m[200..260].iter_mut() {
            *x = 2;
        }
        assert!(lowcov_regions(&m, 50, 100).is_empty(), "{:?}", lowcov_regions(&m, 50, 100));
        // A clean 100 bp region is kept.
        let mut m = vec![0u8; 100];
        for x in m.iter_mut() {
            *x = 1;
        }
        assert_eq!(lowcov_regions(&m, 50, 100), vec![(0, 100, 1)]);
    }
    #[test]
    fn ins_grouping_merges_noisy_copies() {
        let p = QcParams::default();
        let base: Vec<u8> = (0..50).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        // Two noisy copies (3 and 5 substitutions) and one unrelated sequence.
        let other: Vec<u8> = vec![b'A'; 50];
        let indels = vec![
            ins_allele(100, &flip(&base, &[5, 6, 7, 8, 9]), 3, 0),
            ins_allele(100, &base, 6, 4),
            ins_allele(100, &flip(&base, &[1, 2, 3]), 4, 1),
            ins_allele(100, &other, 2, 2),
        ];
        let groups = ins_groups(&indels, &p);
        assert_eq!(groups.len(), 2, "{groups:?}");
        let g = groups.iter().find(|g| g.seq == base).unwrap();
        assert_eq!(g.n, 3);
        assert_eq!((g.anchor, g.len, g.fwd, g.rev), (100, 50, 13, 5));
        let o = groups.iter().find(|g| g.seq == other).unwrap();
        assert_eq!(o.n, 1);
        // Deterministic: reversed input gives the same groups.
        let mut rev = indels.clone();
        rev.reverse();
        let g2 = ins_groups(&rev, &p);
        assert_eq!(groups.len(), g2.len());
        assert!(g2.iter().any(|g| g.seq == base && g.n == 3));
    }
    #[test]
    fn ins_grouping_limits() {
        let p = QcParams::default();
        let base: Vec<u8> = (0..50).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        // Length outside max(2, 10 % of 50) = 5 bp stays separate.
        let mut long = base.clone();
        long.extend_from_slice(b"ACGTAC");
        let g = ins_groups(&[ins_allele(100, &base, 10, 0), ins_allele(100, &long, 8, 0)], &p);
        assert_eq!(g.len(), 2, "{g:?}");
        assert!(g.iter().all(|g| g.n == 1));
        // Different anchors never merge, even for identical sequences.
        let g = ins_groups(&[ins_allele(100, &base, 10, 0), ins_allele(105, &base, 8, 0)], &p);
        assert_eq!(g.len(), 2, "{g:?}");
        assert!(g.iter().all(|g| g.n == 1));
        // Same anchor, identical sequence: merge.
        let g = ins_groups(&[ins_allele(100, &base, 10, 0), ins_allele(100, &base, 8, 0)], &p);
        assert_eq!(g.len(), 1, "{g:?}");
        assert_eq!(g[0].n, 2);
        // min_len 0 disables grouping.
        // Insertions shorter than min_len are not grouped.
        let short = &base[..5];
        assert!(ins_groups(&[ins_allele(100, short, 10, 0), ins_allele(100, short, 8, 0)], &p).is_empty());
        let p0 = QcParams { ins_group_min_len: 0, ..p };
        assert!(ins_groups(&[ins_allele(100, &base, 10, 0), ins_allele(100, &flip(&base, &[1]), 8, 0)], &p0).is_empty());
    }
    #[test]
    fn grouped_ins_max_uses_group_totals() {
        let p = QcParams::default();
        let base: Vec<u8> = (0..50).map(|i| b"ACGT"[(i * 7 % 4) as usize]).collect();
        let indels = vec![ins_allele(100, &base, 6, 4), ins_allele(100, &flip(&base, &[1]), 3, 0)];
        let groups = ins_groups(&indels, &p);
        assert_eq!(groups.len(), 1);
        // The per-position INS max is the group total (13), not the top allele.
        let (ins, _del, _o) = indel_max(200, &indels, &groups, p.ins_group_min_len);
        assert_eq!(ins[100], 13);
        // With grouping off the exact top allele (10) is used instead.
        let (ins0, _d, _o) = indel_max(200, &indels, &[], 0);
        assert_eq!(ins0[100], 10);
    }
}
