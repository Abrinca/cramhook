//! `qc` large-indel regions.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use crate::qc::*;
use crate::scan::*;

pub(crate) fn large_indel_regions(scan: &ScanResult, ctxs: &[ContigCtx], p: &QcParams, low_abs: u32, out: &mut Vec<Region>) {
    let ev = &scan.large;
    let mut i = 0;
    while i < ev.len() {
        let tid = ev[i].tid as usize;
        // All events of this contig (DEL before INS, each sorted by start).
        let mut t_end = i;
        while t_end < ev.len() && ev[t_end].tid as usize == tid {
            t_end += 1;
        }
        // Per-contig depth prefix sum and sorted event starts per type, so the
        // score is O(log n) per call instead of rescanning depth and all events.
        let dp_pref = prefix_u32(&ctxs[tid].dp);
        let mut ins_starts: Vec<u32> = Vec::new();
        let mut del_starts: Vec<u32> = Vec::new();
        for g in &ev[i..t_end] {
            if g.ins {
                ins_starts.push(g.start);
            } else {
                del_starts.push(g.start);
            }
        }
        ins_starts.sort_unstable();
        del_starts.sort_unstable();
        let mut k = i;
        while k < t_end {
            let ins = ev[k].ins;
            let mut l = k + 1;
            let mut reach = if ins { ev[k].start } else { ev[k].end };
            let mut last_start = ev[k].start;
            while l < t_end && ev[l].ins == ins {
                let g = &ev[l];
                let joins = if ins {
                    g.start - last_start <= p.large_ins_window
                } else {
                    g.start < reach
                };
                if !joins {
                    break;
                }
                last_start = g.start;
                reach = reach.max(g.end);
                l += 1;
            }
            emit_large(
                tid,
                ins,
                &ev[k..l],
                scan,
                ctxs,
                p,
                low_abs,
                &dp_pref,
                if ins { &ins_starts } else { &del_starts },
                out,
            );
            k = l;
        }
        i = t_end;
    }
}

fn emit_large(
    tid: usize,
    ins: bool,
    cl: &[LargeEvent],
    scan: &ScanResult,
    ctxs: &[ContigCtx],
    p: &QcParams,
    low_abs: u32,
    dp_pref: &[u64],
    starts: &[u32],
    out: &mut Vec<Region>,
) {
    let mut reads: Vec<u64> = cl.iter().map(|e| e.qhash).collect();
    reads.sort_unstable();
    reads.dedup();
    let support = reads.len() as u32;
    if support < p.large_min_support {
        return;
    }
    let c = &scan.contigs[tid];
    let ctx = &ctxs[tid];
    let len = c.len as u64;
    // An insertion after the last base has start == len: use the preceding
    // anchor so the interval stays inside the contig.
    let s = (cl.iter().map(|e| e.start).min().unwrap() as u64).min(len.saturating_sub(1));
    let e = if ins {
        cl.iter().map(|e| e.end).max().unwrap() as u64 + 1
    } else {
        cl.iter().map(|e| e.end).max().unwrap() as u64
    };
    let e = e.min(c.len as u64).max(s + 1);
    let sd = prefix_range(dp_pref, s, e);
    let mean = sd as f64 / (e - s) as f64;
    if mean <= 0.0 {
        return;
    }
    let ratio = support as f64 / mean;
    if ratio < p.large_min_ratio {
        return;
    }
    let mut sizes: Vec<u32> = cl.iter().map(|e| e.size).collect();
    sizes.sort_unstable();
    let median = sizes[(sizes.len() - 1) / 2];
    let fwd = cl.iter().filter(|e| !e.reverse).count();
    let rev = cl.len() - fwd;
    let (dr, hq) = region_depth(c, ctx, s, e);
    // The candidate gate is the high bar (every reported large indel is
    // "highly significant"), so the tier is the support fraction `f`: red when
    // `f > 0.5`, orange otherwise.
    let tier = if ratio > 0.5 { "red" } else { "orange" };
    let lowcov = mean < low_abs as f64;
    // Significance: supporting reads vs the local large-indel rate (events per
    // covered bp in the window, times the covered bp of the call region). Depth
    // and events come from per-contig prefix sums / sorted starts.
    let center = (s + e) / 2;
    let excl = [(s, e)];
    let bg_cov = prefix_window_excl(dp_pref, len, ctx.circular, center, SIG_SPAN, &excl);
    let bg_ev = count_window_excl(starts, len, ctx.circular, center, SIG_SPAN, &excl);
    let rate = if bg_cov > 0 { bg_ev as f64 / bg_cov as f64 } else { 0.0 };
    let sig = poisson_score(support as u64, rate * sd as f64);
    let mut extra = format!("size={median};support={support};ratio={ratio:.3};fwd={fwd};rev={rev}");
    if lowcov {
        extra.push_str(";lowcov=1");
    }
    extra.push_str(&sig_suffix(sig, ratio));
    out.push(Region {
        tid,
        start: s,
        end: e,
        label: if ins { "large_ins" } else { "large_del" },
        score: support.min(1000),
        tier,
        reason: "large_indel".to_string(),
        depth_ratio: dr,
        hiq_share: hq,
        extra,
        lowcov,
        link: None,
        clip_left: None,
        clip_right: None,
    });
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::qc::test_util::*;

    /// Many large-indel calls must not rescan all events per call.
    #[test]
    fn large_indel_many_calls_is_fast() {
        let n = 10_000u32;
        let spacing = 100u32;
        let len = (100 + spacing * n + 200) as usize;
        let (c, ctx) = mk_contig(len, vec![[0u32; 2]; len], vec![5u32; len]);
        let mut events = Vec::new();
        for g in 0..n {
            let start = 100 + g * spacing;
            for r in 0..3u64 {
                events.push(LargeEvent {
                    tid: 0,
                    ins: false,
                    start,
                    end: start + 60,
                    size: 60,
                    reverse: r % 2 == 1,
                    qhash: g as u64 * 3 + r,
                });
            }
        }
        events.sort_unstable();
        let scan = ScanResult {
            contigs: vec![c],
            indels: Vec::new(),
            large: events,
            records: Default::default(),
            chunks: 1,
            mean_aligned_span: 0.0,
        };
        let p = QcParams::default();
        let mut out = Vec::new();
        let t = std::time::Instant::now();
        large_indel_regions(&scan, &[ctx], &p, 1, &mut out);
        let el = t.elapsed().as_secs_f64();
        assert_eq!(out.len(), n as usize, "{}", out.len());
        assert!(el < 1.0, "large_indel_regions took {el:.3}s");
    }
}
