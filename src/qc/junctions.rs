//! `qc` split-read junction and cross-contig clip event extraction.
//!
//! Moved verbatim out of `qc/circularity.rs` (cleanup-split step 2).

use crate::qc::*;
use crate::qc::circularity::Align;

/// One split-read junction event on a contig: a read leaves the reference beside
/// `(p1, right1)` and continues beside `(p2, right2)`, where `right` says whether
/// the aligned segment lies to the **right** of the breakpoint (reference
/// coordinates). The pair is stored so that `(p1, right1) <= (p2, right2)`
/// lexicographically; that canonical form is invariant under the
/// reverse-complement traversal of the same junction (a read on the other strand
/// visits the same two `(pos, side)` pairs in the opposite order), so both
/// strands of one junction share a key while a tandem/end→start join, a
/// deletion/start→end jump and the two inversion classes stay distinct. `fwd` is
/// the **traversal direction**: `true` when the read visits `(p1, right1)` before
/// `(p2, right2)` (the canonical order), `false` for the reverse-complement
/// traversal (which visits them in the opposite order). This is a property of the
/// junction, not the read's strand, so the two strands of an inverted junction
/// split into `fwd`/`rev` correctly. `read` is the index of the read's group in
/// the sorted chain (its identity for distinct-read support).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct JunctionEvent {
    pub p1: u64,
    pub right1: bool,
    pub p2: u64,
    pub right2: bool,
    pub fwd: bool,
    pub read: u32,
}

/// One cross-contig clip event: consecutive alignments of one read (ordered by
/// query interval) on **different** contigs. `tid1`/`(p1, right1)` is the
/// leaving breakpoint (the 3' reference boundary of the first alignment) and
/// `tid2`/`(p2, right2)` the entering breakpoint (the 5' boundary of the next),
/// with the same strand-aware side logic as [`JunctionEvent`]. Unlike a
/// same-contig junction the pair is **not** canonicalised: the two ends live on
/// different contigs, so the leaving/entering orientation is meaningful for the
/// clip-side matching. `read` is a global distinct-read index (assigned over the
/// name-sorted alignment chains), its identity for support counting.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct XlinkEvent {
    pub tid1: u32,
    pub p1: u64,
    pub right1: bool,
    pub tid2: u32,
    pub p2: u64,
    pub right2: bool,
    pub read: u32,
}

/// Split-read junction events of one contig, from the same alignment chains the
/// circularity pass collects. Segments are consecutive in read order (sorted by
/// query interval) and must make forward progress (strictly increasing query
/// start and end); the first leaves its 3' reference boundary and the next enters
/// its 5' boundary, strand-aware. Events with `|p1 - p2|` below
/// `--junction-min-dist` or a read-coordinate gap/overlap above
/// `--junction-max-gap` (extra unaligned sequence is an insertion, not a
/// junction) are dropped, as is a junction with either breakpoint within
/// `--junction-edge-excl` of a contig end on a circular contig (the end->start
/// wrap is the special case of both). Output is sorted, so it does not depend on
/// the thread count.
pub(crate) fn contig_junctions(recs: &[(Vec<u8>, Align)], len: u64, circular: bool, p: &QcParams) -> Vec<JunctionEvent> {
    let min_dist = p.junction_min_dist as u64;
    let max_gap = p.junction_max_gap as u64;
    let edge = p.junction_edge_excl as u64;
    let near_edge = |pos: u64| pos < edge || pos + edge >= len;
    let mut sorted: Vec<&(Vec<u8>, Align)> = recs.iter().collect();
    sorted.sort_by(|a, b| {
        (a.0.as_slice(), a.1.q_span(), a.1.pos, a.1.ref_end, a.1.rev).cmp(&(
            b.0.as_slice(),
            b.1.q_span(),
            b.1.pos,
            b.1.ref_end,
            b.1.rev,
        ))
    });
    let mut out = Vec::new();
    let mut i = 0;
    let mut idx = 0u32;
    while i < sorted.len() {
        let mut j = i;
        while j < sorted.len() && sorted[j].0 == sorted[i].0 {
            j += 1;
        }
        // Within one read the outer sort orders the segments by query interval.
        for w in sorted[i..j].windows(2) {
            let (first, next) = (w[0].1, w[1].1);
            // Consecutive alignments must make forward progress in read
            // coordinates (strictly increasing start and end), so identical,
            // contained and tied intervals at distant loci are not junctions.
            let (qs_first, qe_first) = first.q_span();
            let (qs_next, qe_next) = next.q_span();
            if qs_next <= qs_first || qe_next <= qe_first {
                continue;
            }
            if qs_next.abs_diff(qe_first) > max_gap {
                continue;
            }
            // Breakpoints and the side each segment lies on: the 3' boundary of
            // the first segment and the 5' boundary of the next; a segment maps
            // to the right of its 3' boundary when reverse, and of its 5'
            // boundary when forward.
            let (pf, rf) = if first.rev { (first.pos, true) } else { (first.ref_end, false) };
            let (pn, rn) = if next.rev { (next.ref_end, false) } else { (next.pos, true) };
            if pf.abs_diff(pn) < min_dist {
                continue;
            }
            if circular && (near_edge(pf) || near_edge(pn)) {
                continue;
            }
            // Canonical key: the read's (pos, side) order, or its reverse (the
            // reverse-complement traversal of the same junction). The sides are
            // kept, so tandem/deletion and the two inversion classes differ; the
            // traversal direction (which order the read used) becomes fwd/rev.
            let forward = (pf, rf) <= (pn, rn);
            let (p1, right1, p2, right2) = if forward {
                (pf, rf, pn, rn)
            } else {
                (pn, rn, pf, rf)
            };
            out.push(JunctionEvent { p1, right1, p2, right2, fwd: forward, read: idx });
        }
        i = j;
        idx += 1;
    }
    out.sort();
    out
}

/// Cross-contig clip events from the same collected alignment chains. Every
/// read's alignments on all contigs are ordered by query interval (as
/// [`contig_junctions`] does per contig); each consecutive pair on **different**
/// contigs with forward progress in read coordinates (strictly increasing query
/// start and end) and a read-coordinate gap *or overlap* of at most `max_gap`
/// becomes an event. Overlaps are expected at repeat-mediated joins (the two
/// alignments of one read may share a few hundred bp), so `max_gap` is the
/// read-coordinate distance `|qs_next - qe_first|`. No contig-end exclusion
/// (the end->start wrap is a same-contig phenomenon). Output is sorted, so it
/// does not depend on the thread count.
pub(crate) fn xlink_events(by_contig: &[Vec<(Vec<u8>, Align)>], max_gap: u64) -> Vec<XlinkEvent> {
    // One entry per alignment, borrowing the read name from `by_contig` (nothing
    // is copied), sorted by name then query interval: a read's chain on all
    // contigs becomes contiguous.
    let mut all: Vec<(&[u8], u32, Align)> = Vec::new();
    for (tid, recs) in by_contig.iter().enumerate() {
        for (name, a) in recs {
            all.push((name.as_slice(), tid as u32, *a));
        }
    }
    all.sort_unstable_by(|x, y| {
        (x.0, x.2.q_span(), x.2.pos, x.2.ref_end, x.2.rev, x.1).cmp(&(
            y.0,
            y.2.q_span(),
            y.2.pos,
            y.2.ref_end,
            y.2.rev,
            y.1,
        ))
    });
    let mut out = Vec::new();
    let mut i = 0;
    let mut idx = 0u32;
    while i < all.len() {
        let mut j = i;
        while j < all.len() && all[j].0 == all[i].0 {
            j += 1;
        }
        for w in all[i..j].windows(2) {
            let (tid1, first) = (w[0].1, w[0].2);
            let (tid2, next) = (w[1].1, w[1].2);
            if tid1 == tid2 {
                continue;
            }
            let (qs_first, qe_first) = first.q_span();
            let (qs_next, qe_next) = next.q_span();
            if qs_next <= qs_first || qe_next <= qe_first {
                continue;
            }
            if qs_next.abs_diff(qe_first) > max_gap {
                continue;
            }
            let (p1, right1) = if first.rev { (first.pos, true) } else { (first.ref_end, false) };
            let (p2, right2) = if next.rev { (next.ref_end, false) } else { (next.pos, true) };
            out.push(XlinkEvent { tid1, p1, right1, tid2, p2, right2, read: idx });
        }
        i = j;
        idx += 1;
    }
    out.sort_unstable();
    out
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::qc::test_util::*;

    /// A clean forward junction pair: the first segment ends at `b`, the next
    /// starts at `a`, with `gap` read bases between the two query intervals.
    fn fwd_pair(a: u64, b: u64, gap: u64) -> Vec<(Vec<u8>, Align)> {
        vec![
            ("r".as_bytes().to_vec(), rec(b - 500, b, false, 0, 0, 500)),
            ("r".as_bytes().to_vec(), rec(a, a + 500, false, 500 + gap, 0, 500)),
        ]
    }
    /// Junction events: the clean pair is found, the distance and read-gap gates
    /// apply (gap at the limit is kept), and junctions near a contig end are
    /// excluded on a circular contig but kept on a linear one.
    #[test]
    fn junction_events_forward_and_filters() {
        let p = QcParams::default();
        let len = 30_000u64;
        // A forward end->start (tandem/circle) join at 12_000 -> 2_000: the read
        // leaves the reference left of 12_000 and continues right of 2_000. The
        // canonical order starts at 2_000, so this read is the reversed traversal.
        let ev = contig_junctions(&fwd_pair(2_000, 12_000, 0), len, false, &p);
        assert_eq!(ev.len(), 1, "{ev:?}");
        assert_eq!((ev[0].p1, ev[0].right1, ev[0].p2, ev[0].right2, ev[0].fwd), (2_000, true, 12_000, false, false));
        // Below `--junction-min-dist` (1000).
        assert!(contig_junctions(&fwd_pair(11_500, 12_000, 0), len, false, &p).is_empty());
        // Read-coordinate gap above `--junction-max-gap` (200); at the limit kept.
        assert!(contig_junctions(&fwd_pair(2_000, 12_000, 201), len, false, &p).is_empty());
        assert_eq!(contig_junctions(&fwd_pair(2_000, 12_000, 200), len, false, &p).len(), 1);
        // Origin wrap: excluded when circular, kept when linear.
        let wrap = fwd_pair(0, len, 0);
        assert!(contig_junctions(&wrap, len, true, &p).is_empty());
        assert_eq!(contig_junctions(&wrap, len, false, &p).len(), 1);
        // The same wrap on a reverse read swaps the two sides and must still be
        // excluded.
        let wrap_rev = vec![
            ("r".as_bytes().to_vec(), rec(0, 500, true, 0, 0, 500)),
            ("r".as_bytes().to_vec(), rec(len - 500, len, true, 0, 500, 500)),
        ];
        assert!(contig_junctions(&wrap_rev, len, true, &p).is_empty());
    }
    /// `--junction-edge-excl`: on a circular contig a junction with either
    /// breakpoint within 2000 bp of a contig end is dropped (here the far
    /// breakpoint is interior); the same junction is kept on a linear contig.
    #[test]
    fn junction_edge_exclusion_circular_only() {
        let p = QcParams::default();
        let len = 30_000u64;
        // b = 1_500 (within 2000 of the start), a = 20_000 (interior).
        let near_start = fwd_pair(20_000, 1_500, 0);
        assert!(contig_junctions(&near_start, len, true, &p).is_empty());
        assert_eq!(contig_junctions(&near_start, len, false, &p).len(), 1);
        // a = 29_000 (within 2000 of the end), b = 10_000 (interior).
        let near_end = fwd_pair(29_000, 10_000, 0);
        assert!(contig_junctions(&near_end, len, true, &p).is_empty());
        assert_eq!(contig_junctions(&near_end, len, false, &p).len(), 1);
        // Both breakpoints interior: kept on a circular contig.
        assert_eq!(contig_junctions(&fwd_pair(2_000, 12_000, 0), len, true, &p).len(), 1);
    }
    /// Cross-contig clip events: two alignments of one read on different contigs
    /// with forward progress and a read-coordinate gap *or overlap* (a
    /// repeat-mediated join overlaps) within `--xlink-max-gap` become an event,
    /// with the same leaving/entering side logic as `contig_junctions`; a gap
    /// above the limit, backwards progress, and a same-contig pair are dropped.
    #[test]
    fn xlink_events_across_contigs_and_gap_or_overlap() {
        // Forward alignment with an explicit leading clip (`lead`) so `q_span`
        // is `(lead, lead + aln)`.
        let al = |pos: u64, ref_end: u64, lead: u64, aln: u64| rec(pos, ref_end, false, lead, 0, aln);
        let name = b"r".to_vec();
        // ctg0: q_span (0, 1000), 3' boundary 5000. ctg1: q_span (200, 1200):
        // 800 bp of read overlap with ctg0's segment.
        let mut by: Vec<Vec<(Vec<u8>, Align)>> = vec![Vec::new(), Vec::new()];
        by[0].push((name.clone(), al(4_000, 5_000, 0, 1_000)));
        by[1].push((name.clone(), al(100, 1_100, 200, 1_000)));
        let ev = xlink_events(&by, 1_000);
        assert_eq!(ev.len(), 1, "{ev:?}");
        assert_eq!((ev[0].tid1, ev[0].p1, ev[0].right1), (0, 5_000, false));
        assert_eq!((ev[0].tid2, ev[0].p2, ev[0].right2), (1, 100, true));
        // Read gap above the limit is dropped; at the limit kept.
        let mut gap = by.clone();
        gap[1][0].1 = al(100, 1_100, 2_001, 1_000); // q_span (2001, 3001): gap 1001
        assert!(xlink_events(&gap, 1_000).is_empty());
        gap[1][0].1 = al(100, 1_100, 2_000, 1_000); // gap 1000
        assert_eq!(xlink_events(&gap, 1_000).len(), 1);
        // Backwards progress in read coordinates is not a junction.
        let mut back: Vec<Vec<(Vec<u8>, Align)>> = vec![Vec::new(), Vec::new()];
        back[0].push((name.clone(), al(4_000, 5_000, 0, 1_000)));
        back[1].push((name.clone(), al(100, 1_100, 500, 400))); // q_span (500, 900): contained
        assert!(xlink_events(&back, 1_000).is_empty());
        // Same-contig consecutive alignments are not cross-contig events.
        let mut same: Vec<Vec<(Vec<u8>, Align)>> = vec![Vec::new(), Vec::new()];
        same[0].push((name.clone(), al(4_000, 5_000, 0, 1_000)));
        same[0].push((name.clone(), al(9_000, 10_000, 2_000, 1_000)));
        assert!(xlink_events(&same, 1_000).is_empty());
    }
    /// The same junction read on the forward and on the reverse strand is one
    /// canonical event: the reverse-complement traversal swaps the two sides and
    /// reverses their order, and the canonical key is unchanged.
    #[test]
    fn junction_reverse_complement_pairs_match() {
        let p = QcParams::default();
        let len = 30_000u64;
        let fwd = fwd_pair(2_000, 12_000, 0);
        // Reverse-complement layout: the read's 5' part maps reverse right of
        // 2_000, its 3' part reverse left of 12_000.
        let rev = vec![
            ("r".as_bytes().to_vec(), rec(2_000, 2_500, true, 0, 0, 500)),
            ("r".as_bytes().to_vec(), rec(11_500, 12_000, true, 0, 500, 500)),
        ];
        let ef = contig_junctions(&fwd, len, false, &p);
        let er = contig_junctions(&rev, len, false, &p);
        assert_eq!(ef.len(), 1, "{ef:?}");
        assert_eq!(er.len(), 1, "{er:?}");
        assert_eq!(
            (ef[0].p1, ef[0].right1, ef[0].p2, ef[0].right2),
            (er[0].p1, er[0].right1, er[0].p2, er[0].right2)
        );
        // Opposite traversal directions: the canonical order vs its reverse.
        assert!(ef[0].fwd != er[0].fwd, "{ef:?} {er:?}");
    }
    /// An inverted junction (both breakpoints on the same side) read on the two
    /// strands: the traversals differ, so `fwd`/`rev` are split, and they cluster
    /// into one event.
    #[test]
    fn junction_inverted_reverse_complement_counts() {
        let p = QcParams::default();
        let len = 30_000u64;
        // Read `a`: first forward leaves left of 12_000, next reverse enters left
        // of 2_000 (both sides left -> inverted_l). Read `b` is its
        // reverse-complement traversal: first forward leaves left of 2_000, next
        // reverse enters left of 12_000.
        let recs = vec![
            ("a".as_bytes().to_vec(), rec(11_500, 12_000, false, 0, 0, 500)),
            ("a".as_bytes().to_vec(), rec(1_500, 2_000, true, 0, 500, 500)),
            ("b".as_bytes().to_vec(), rec(1_500, 2_000, false, 0, 0, 500)),
            ("b".as_bytes().to_vec(), rec(11_500, 12_000, true, 0, 500, 500)),
        ];
        let ev = contig_junctions(&recs, len, false, &p);
        assert_eq!(ev.len(), 2, "{ev:?}");
        assert_eq!((ev[0].p1, ev[0].right1, ev[0].p2, ev[0].right2), (2_000, false, 12_000, false));
        assert_eq!(
            (ev[1].p1, ev[1].right1, ev[1].p2, ev[1].right2),
            (2_000, false, 12_000, false)
        );
        assert!(ev[0].fwd != ev[1].fwd, "{ev:?}");
        // One cluster, with the traversals split (not all on one side).
        let cl = crate::qc::junction_clusters(&ev, 50);
        assert_eq!(cl.len(), 1, "{cl:?}");
        assert_eq!((cl[0].reads, cl[0].fwd, cl[0].rev), (2, 1, 1), "{cl:?}");
    }
    /// A tandem (end->start) join and the reciprocal direct (start->end) join at
    /// the same two positions are different junctions: their canonical keys keep
    /// the sides and differ.
    #[test]
    fn junction_reciprocal_direct_joins_distinct() {
        let p = QcParams::default();
        let len = 30_000u64;
        // end(12_000) -> start(2_000): leaves left of 12_000, continues right of 2_000.
        let tandem = contig_junctions(&fwd_pair(2_000, 12_000, 0), len, false, &p);
        // start(2_000) -> end(12_000): leaves left of 2_000, continues right of 12_000.
        let deletion = vec![
            ("r".as_bytes().to_vec(), rec(1_500, 2_000, false, 0, 0, 500)),
            ("r".as_bytes().to_vec(), rec(12_000, 12_500, false, 500, 0, 500)),
        ];
        let ed = contig_junctions(&deletion, len, false, &p);
        assert_eq!(tandem.len(), 1, "{tandem:?}");
        assert_eq!(ed.len(), 1, "{ed:?}");
        assert_eq!((tandem[0].p1, tandem[0].p2), (2_000, 12_000));
        assert_eq!((ed[0].p1, ed[0].p2), (2_000, 12_000));
        assert_ne!((tandem[0].right1, tandem[0].right2), (ed[0].right1, ed[0].right2));
        // tandem: continues right of 2_000, left of 12_000; deletion: the mirror.
        assert_eq!((tandem[0].right1, tandem[0].right2), (true, false));
        assert_eq!((ed[0].right1, ed[0].right2), (false, true));
    }
    /// Adjacency requires strictly increasing query start and end: identical,
    /// contained and tied intervals at distant loci are not junctions, while a
    /// small read-coordinate overlap is accepted.
    #[test]
    fn junction_adjacency_requires_forward_progress() {
        let p = QcParams::default();
        let len = 30_000u64;
        // Identical query spans at distant loci.
        let identical = vec![
            ("r".as_bytes().to_vec(), rec(12_000, 12_500, false, 0, 0, 500)),
            ("r".as_bytes().to_vec(), rec(2_000, 2_500, false, 0, 0, 500)),
        ];
        assert!(contig_junctions(&identical, len, false, &p).is_empty(), "identical spans accepted");
        // Contained span (query [100,300) inside [0,500)).
        let contained = vec![
            ("r".as_bytes().to_vec(), rec(12_000, 12_500, false, 0, 0, 500)),
            ("r".as_bytes().to_vec(), rec(2_000, 2_200, false, 100, 0, 200)),
        ];
        assert!(contig_junctions(&contained, len, false, &p).is_empty(), "contained span accepted");
        // Tied query start.
        let tied = vec![
            ("r".as_bytes().to_vec(), rec(12_000, 12_500, false, 0, 0, 500)),
            ("r".as_bytes().to_vec(), rec(2_000, 2_700, false, 0, 0, 700)),
        ];
        assert!(contig_junctions(&tied, len, false, &p).is_empty(), "tied start accepted");
        // Small overlap: query [450,950) after [0,500).
        let overlap = vec![
            ("r".as_bytes().to_vec(), rec(12_000, 12_500, false, 0, 0, 500)),
            ("r".as_bytes().to_vec(), rec(2_000, 2_500, false, 450, 0, 500)),
        ];
        assert_eq!(contig_junctions(&overlap, len, false, &p).len(), 1, "valid overlap rejected");
    }
}
