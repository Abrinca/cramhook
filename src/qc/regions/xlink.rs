//! `qc` cross-contig clip links.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use crate::qc::*;

// ---- Cross-contig clip links ------------------------------------------------
//
// Starting gates (not tuned; see `docs/project/decisions.md`): the partner cluster needs
// at least this many reads and this share of the hotspot's clipped reads on the
// matched side, and its breakpoints cluster within this many bp (leader
// clustering).
const XLINK_MIN_READS: u64 = 5;
const XLINK_MIN_FRAC: f64 = 0.5;
const XLINK_PARTNER_TOL: u64 = 200;

/// One cross-contig clip link: a `clip_hotspot` region end joined to a cluster
/// of breakpoints on another contig by reads whose other alignment starts there.
#[derive(Clone, Debug)]
pub(crate) struct Xlink {
    /// Index into the `regions` slice (`None` never occurs: links point at a
    /// `clip_hotspot` region).
    pub region: usize,
    /// Hotspot (contig id and 0-based half-open interval).
    pub tid1: usize,
    pub start1: u64,
    pub end1: u64,
    /// Partner breakpoint cluster, 0-based half-open `[min, max + 1)`, anchored
    /// inside the partner contig (a breakpoint at `len` anchors on `len - 1`).
    pub tid2: usize,
    pub start2: u64,
    pub end2: u64,
    pub reads: u64,
    /// `reads / depth` at the hotspot side.
    pub f: f64,
    /// Depth at the hotspot side (the denominator of `f`), noted on the partner
    /// finding as its `partner_depth`.
    pub hotspot_depth: f64,
    /// Depth at the partner breakpoint, noted on the hotspot finding as its
    /// `partner_depth`.
    pub partner_depth: f64,
    /// Finding containing the hotspot region and one containing the partner
    /// breakpoint (`None` -> written as `.`).
    pub finding1: Option<usize>,
    pub finding2: Option<usize>,
}

/// A candidate partner-breakpoint cluster (same contig and side).
struct PartnerCluster {
    tid: usize,
    right: bool,
    /// Representative (first) breakpoint, for leader clustering.
    pos: u64,
    min: u64,
    max: u64,
    reads: u64,
}

/// Anchor a breakpoint inside `[0, len)`: a breakpoint at the contig end (`len`)
/// is the "insertion after the last base" case, anchored on the last base. This
/// mirrors `junction_regions` / `emit_large`.
fn anchor(pos: u64, len: u64) -> u64 {
    if len == 0 { 0 } else { pos.min(len - 1) }
}

/// Whether the half-open window `[s, e)` overlaps the finding. On a circular
/// contig a finding may be stored in the wrapped `[start, end)` convention where
/// `end` exceeds `len`; a linear contig never wraps.
fn finding_overlaps(f: &Finding, s: u64, e: u64, len: u64, circular: bool) -> bool {
    finding_covers(f, s, len, circular) || (e > 0 && finding_covers(f, e - 1, len, circular))
}

/// Whether `pos` lies inside the finding. Only a **circular** contig wraps, and
/// then a wrapped finding has `start + arc > len` (its two BED rows are
/// `[start, len)` and `[0, start + arc - len)`); on a linear contig the test is
/// the plain `start <= pos < end`, so `pos == len` is never inside anything.
fn finding_covers(f: &Finding, pos: u64, len: u64, circular: bool) -> bool {
    if len == 0 {
        return false;
    }
    if !circular {
        return pos >= f.start && pos < f.end;
    }
    let arc = f.end.saturating_sub(f.start);
    if arc == 0 || f.start >= len {
        return false;
    }
    (pos + len - f.start) % len < arc
}

/// Cross-contig clip links: for every `clip_hotspot` region, match cross-contig
/// events with an endpoint at the region (within `--clip-merge`, on the side the
/// clip is on: a left clip joins an aligned segment to its right), cluster the
/// partner breakpoints (same contig and side, within [`XLINK_PARTNER_TOL`]) and
/// keep the largest cluster when it passes [`XLINK_MIN_READS`] and
/// [`XLINK_MIN_FRAC`]. The finding ids are resolved against the already-built
/// findings; the linking itself never moves regions or findings.
#[allow(clippy::too_many_arguments)]
pub(crate) fn collect_xlinks(
    events: &[circularity::XlinkEvent],
    regions: &[Region],
    ctxs: &[ContigCtx],
    lens: &[u64],
    findings: &[Finding],
    clip_merge: u32,
) -> Vec<Xlink> {
    let merge = clip_merge as u64;
    let mut out = Vec::new();
    for (ri, r) in regions.iter().enumerate() {
        if r.label != "clip_hotspot" {
            continue;
        }
        let len = lens.get(r.tid).copied().unwrap_or(0);
        let circular1 = ctxs.get(r.tid).is_some_and(|c| c.circular);
        // The matched side's own cluster (count-weighted median position and
        // aligned-side depth), so the link's `f` equals that side's ratio.
        for (side, match_right) in [(r.clip_left, true), (r.clip_right, false)] {
            let Some(side) = side else { continue };
            if side.count == 0 {
                continue;
            }
            let count = side.count;
            let lo = r.start.saturating_sub(merge);
            let hi = r.end.saturating_add(merge);
            // Partner endpoints of events with an endpoint at the hotspot.
            let mut items: Vec<(usize, bool, u64, u32)> = Vec::new();
            for e in events {
                if e.tid1 as usize == r.tid && e.right1 == match_right && e.p1 >= lo && e.p1 <= hi {
                    items.push((e.tid2 as usize, e.right2, e.p2, e.read));
                }
                if e.tid2 as usize == r.tid && e.right2 == match_right && e.p2 >= lo && e.p2 <= hi {
                    items.push((e.tid1 as usize, e.right1, e.p1, e.read));
                }
            }
            if items.is_empty() {
                continue;
            }
            // Leader (single-linkage) clusters by (contig, side, within tol),
            // over the sorted events so the result is deterministic.
            items.sort_unstable();
            let mut clusters: Vec<PartnerCluster> = Vec::new();
            let mut members: Vec<Vec<u32>> = Vec::new();
            for (tid, right_p, p, read) in items {
                let mut placed = false;
                for (ci, c) in clusters.iter_mut().enumerate() {
                    if c.tid == tid && c.right == right_p && c.pos.abs_diff(p) <= XLINK_PARTNER_TOL {
                        c.min = c.min.min(p);
                        c.max = c.max.max(p);
                        members[ci].push(read);
                        placed = true;
                        break;
                    }
                }
                if !placed {
                    clusters.push(PartnerCluster { tid, right: right_p, pos: p, min: p, max: p, reads: 0 });
                    members.push(vec![read]);
                }
            }
            for (c, ms) in clusters.iter_mut().zip(members.iter_mut()) {
                ms.sort_unstable();
                ms.dedup();
                c.reads = ms.len() as u64;
            }
            // Largest cluster; ties break on (contig, side, breakpoint).
            let best = clusters.iter().min_by(|a, b| {
                (std::cmp::Reverse(a.reads), a.tid, a.right, a.pos)
                    .cmp(&(std::cmp::Reverse(b.reads), b.tid, b.right, b.pos))
            });
            let Some(c) = best else { continue };
            if c.reads < XLINK_MIN_READS || (c.reads as f64) < XLINK_MIN_FRAC * count as f64 {
                continue;
            }
            // `f` uses the matched side's own depth (the clip's aligned side),
            // not the mean of the two ends and not the region's boundary.
            if side.depth <= 0.0 {
                continue;
            }
            // Anchor the partner cluster inside its contig: a breakpoint at the
            // contig end (`len`) is not a valid interval.
            let plen = lens.get(c.tid).copied().unwrap_or(0);
            if plen == 0 {
                continue;
            }
            let start2 = anchor(c.min, plen);
            let end2 = (anchor(c.max, plen) + 1).min(plen);
            let circular2 = ctxs.get(c.tid).is_some_and(|cc| cc.circular);
            let finding1 = findings
                .iter()
                .position(|f| f.tid == r.tid && finding_overlaps(f, r.start, r.end, len, circular1));
            let finding2 = findings.iter().position(|f| {
                f.tid == c.tid && finding_covers(f, start2, plen, circular2)
            });
            let partner_depth = ctxs
                .get(c.tid)
                .and_then(|cc| cc.dp.get(start2 as usize))
                .copied()
                .unwrap_or(0) as f64;
            out.push(Xlink {
                region: ri,
                tid1: r.tid,
                start1: r.start,
                end1: r.end,
                tid2: c.tid,
                start2,
                end2,
                reads: c.reads,
                f: c.reads as f64 / side.depth,
                hotspot_depth: side.depth,
                partner_depth,
                finding1,
                finding2,
            });
        }
    }
    out
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::qc::test_util::*;

    /// Cross-contig clip links: a left-clip hotspot links to the partner cluster
    /// its reads join on another contig; the `--clip-merge` window, the read
    /// gate and the side share gate apply, and the partner breakpoints cluster.
    #[test]
    fn xlink_gate_reads_and_share() {
        let p = QcParams::default();
        let len = 10_000usize;
        let (_, ctx) = mk_contig(len, vec![[0u32; 2]; len], vec![100u32; len]);
        let (_, ctx1) = mk_contig(len, vec![[0u32; 2]; len], vec![100u32; len]);
        let ctxs = vec![ctx, ctx1];
        let lens = vec![len as u64, len as u64];
        let mut r = freg(0, "clip_hotspot", "yellow", 1_000, 1_100, false);
        r.reason = "clip_left".to_string();
        r.extra = "left=30;right=0;ratio_left=0.3;ratio_right=0.0".to_string();
        // The matched side's own cluster: count-weighted median 1_050, depth 100
        // (count / depth = ratio 0.3).
        r.clip_left = Some(ClipSide { pos: 1_050, count: 30, ratio: 0.3, depth: 100.0, span: (1_000, 1_100) });
        // A left clip (aligned segment to the right, `right=true`) needs events
        // with their endpoint at contig 0 with `right2=true`.
        let ev = |pos: u64, read: u32| circularity::XlinkEvent {
            tid1: 1,
            p1: pos,
            right1: true,
            tid2: 0,
            p2: 1_050,
            right2: true,
            read,
        };
        let events: Vec<_> = (0..20).map(|i| ev(5_000, i)).collect();
        let x = collect_xlinks(&events, std::slice::from_ref(&r), &ctxs, &lens, &[], p.clip_merge);
        assert_eq!(x.len(), 1, "{x:?}");
        assert_eq!((x[0].tid2, x[0].start2, x[0].end2, x[0].reads), (1, 5_000, 5_001, 20));
        // f = reads / (matched side depth) = 20 / 100.
        assert!((x[0].f - 0.2).abs() < 1e-9, "{}", x[0].f);
        // partner depth = raw depth at the partner breakpoint (100 everywhere).
        assert!((x[0].partner_depth - 100.0).abs() < 1e-9, "{}", x[0].partner_depth);
        // Below the min-reads gate.
        let few: Vec<_> = (0..4).map(|i| ev(5_000, i)).collect();
        assert!(collect_xlinks(&few, std::slice::from_ref(&r), &ctxs, &lens, &[], p.clip_merge).is_empty());
        // Below 50 % of the side's clipped reads: 20 of 60 fails.
        let mut r60 = r.clone();
        if let Some(s) = &mut r60.clip_left {
            s.count = 60;
        }
        assert!(collect_xlinks(&events, std::slice::from_ref(&r60), &ctxs, &lens, &[], p.clip_merge).is_empty());
        // Partner breakpoints cluster (within 200 bp) and the largest wins.
        let mut split: Vec<_> = (0..20).map(|i| ev(5_000, i)).collect();
        split.extend((0..25).map(|i| ev(8_000, 100 + i)));
        // 8_000 is >200 away from 5_000 -> two clusters; the 25-read one wins.
        let x = collect_xlinks(&split, std::slice::from_ref(&r), &ctxs, &lens, &[], p.clip_merge);
        assert_eq!(
            (x.len(), x[0].start2, x[0].reads),
            (1, 8_000, 25),
            "largest partner cluster: {x:?}"
        );
        // An event outside the region ± clip-merge window is not matched.
        let far: Vec<_> = (0..20).map(|i| ev(5_000, i)).collect();
        let mut r_far = r.clone();
        r_far.start = 5_000;
        r_far.end = 5_100;
        assert!(collect_xlinks(&far, std::slice::from_ref(&r_far), &ctxs, &lens, &[], p.clip_merge).is_empty());
    }
    /// A finding literal for the link-attribution tests.
    fn finding(tid: usize, start: u64, end: u64) -> Finding {
        Finding {
            tid,
            start,
            end,
            tier: "red",
            lowcov: false,
            score: 1,
            likely: "unclear",
            evidence: Evidence::default(),
            members: Vec::new(),
        }
    }
    /// The link's `f` is `reads / (matched side depth)` from the side's own
    /// median, not the depth at the combined region boundary (reviewer issue 1):
    /// here the boundary depth is 100 but the side depth is 10, so `f` is 1.0,
    /// not 0.1.
    #[test]
    fn xlink_f_uses_matched_side_depth() {
        let p = QcParams::default();
        let len = 10_000usize;
        let mut dp = vec![10u32; len];
        for x in dp.iter_mut().take(1_010).skip(1_000) {
            *x = 100;
        }
        let (_, ctx) = mk_contig(len, vec![[0u32; 2]; len], dp);
        let ctxs = vec![ctx];
        let lens = vec![len as u64, len as u64];
        let mut r = freg(0, "clip_hotspot", "yellow", 1_000, 1_010, false);
        r.reason = "clip_left".to_string();
        r.extra = "left=10;right=0".to_string();
        r.clip_left = Some(ClipSide { pos: 1_000, count: 10, ratio: 1.0, depth: 10.0, span: (1_000, 1_010) });
        let events: Vec<_> = (0..10)
            .map(|read| circularity::XlinkEvent {
                tid1: 1,
                p1: 5_000,
                right1: true,
                tid2: 0,
                p2: 1_005,
                right2: true,
                read,
            })
            .collect();
        let x = collect_xlinks(&events, std::slice::from_ref(&r), &ctxs, &lens, &[], p.clip_merge);
        assert_eq!(x.len(), 1, "{x:?}");
        assert!((x[0].f - 1.0).abs() < 1e-9, "f must use the side depth: {}", x[0].f);
        // The partner finding's mirrored `link=` note uses this side's depth as
        // its `partner_depth` (`f`'s denominator).
        assert!((x[0].hotspot_depth - 10.0).abs() < 1e-9, "{}", x[0].hotspot_depth);
    }
    /// A partner breakpoint at the contig end (`len`) is anchored on the last
    /// base, so the BEDPE interval stays inside the contig and the finding is
    /// looked up with the same clamped coordinate (reviewer issue 2). On a
    /// **linear** contig `pos == len` must not wrap onto an origin finding.
    #[test]
    fn xlink_terminal_partner_clamped() {
        let p = QcParams::default();
        let len = 10_000usize;
        let (_, ctxa) = mk_contig(len, vec![[0u32; 2]; len], vec![10u32; len]);
        let (_, ctxb) = mk_contig(len, vec![[0u32; 2]; len], vec![10u32; len]);
        let ctxs = vec![ctxa, ctxb];
        let lens = vec![len as u64, len as u64];
        let mut r = freg(0, "clip_hotspot", "yellow", 5_000, 5_001, false);
        r.reason = "clip_left".to_string();
        r.extra = "left=10;right=0".to_string();
        r.clip_left = Some(ClipSide { pos: 5_000, count: 10, ratio: 1.0, depth: 10.0, span: (5_000, 5_001) });
        // The read leaves contig 1 exactly at its end (ref_end == len).
        let events: Vec<_> = (0..10)
            .map(|read| circularity::XlinkEvent {
                tid1: 1,
                p1: len as u64,
                right1: false,
                tid2: 0,
                p2: 5_000,
                right2: true,
                read,
            })
            .collect();
        // An unrelated origin finding on the partner contig and the hotspot's own
        // finding; the clamped breakpoint (9_999) is in neither.
        let findings = vec![finding(0, 4_000, 6_000), finding(1, 0, 5_000)];
        let x = collect_xlinks(&events, std::slice::from_ref(&r), &ctxs, &lens, &findings, p.clip_merge);
        assert_eq!(x.len(), 1, "{x:?}");
        assert_eq!((x[0].tid2, x[0].start2, x[0].end2), (1, 9_999, 10_000), "{x:?}");
        assert!(x[0].start2 < len as u64 && x[0].end2 <= len as u64, "interval out of contig: {x:?}");
        assert_eq!(x[0].finding1, Some(0), "hotspot finding");
        assert_eq!(x[0].finding2, None, "9_999 must not match the origin finding: {x:?}");
    }
    /// `finding_covers` wraps only on a circular contig; on a linear contig
    /// `pos == len` is never inside a finding, so it cannot attach a terminal
    /// breakpoint to an unrelated origin finding.
    #[test]
    fn finding_covers_linear_and_circular() {
        let origin = finding(0, 0, 5_000);
        assert!(finding_covers(&origin, 100, 10_000, false));
        assert!(!finding_covers(&origin, 5_000, 10_000, false));
        assert!(!finding_covers(&origin, 10_000, 10_000, false), "linear pos == len must not wrap");
        // Circular: a finding that wraps the origin (stored end > len).
        let wrap = finding(0, 9_000, 12_000);
        assert!(finding_covers(&wrap, 9_000, 10_000, true));
        assert!(finding_covers(&wrap, 1_000, 10_000, true), "0..2_000 is inside the arc");
        assert!(!finding_covers(&wrap, 5_000, 10_000, true));
    }
}
