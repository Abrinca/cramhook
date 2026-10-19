//! Shared `qc` unit-test helpers.
//!
//! Used by more than one module's `#[cfg(test)] mod tests` (cleanup-split step 2).

use crate::qc::*;
use crate::qc::circularity::Align;
use crate::scan::ContigCounts;

/// A tiered region of `label` over `[start, end)`.
pub(crate) fn freg(tid: usize, label: &'static str, tier: &'static str, start: u64, end: u64, lowcov: bool) -> Region {
    Region {
        tid,
        start,
        end,
        label,
        score: 1,
        tier,
        reason: String::new(),
        depth_ratio: None,
        hiq_share: None,
        extra: String::new(),
        lowcov,
        link: None,
        clip_left: None,
        clip_right: None,
    }
}

pub(crate) fn mk_contig(len: usize, clips: Vec<[u32; 2]>, dp: Vec<u32>) -> (ContigCounts, ContigCtx) {
    let c = ContigCounts {
        name: "chr".to_string(),
        len,
        bases: vec![[0u32; 8]; len],
        cov: vec![[0u32; crate::scan::NCOV]; len],
        clips,
        edge_clips: 0,
        aligned_bases: 0,
    };
    let mode = depth_mode(&dp);
    let ctx = ContigCtx {
        hp_run: vec![0u8; len],
        hp_fwd: vec![0u8; len],
        str_span: vec![0u16; len],
        dp,
        na: vec![0u8; len],
        lowcov: vec![0u8; len],
        circular: false,
        circular_by: "verdict",
        mode,
        mean: 0.0,
        depth_skipped: None,
    };
    (c, ctx)
}


pub(crate) fn rec(pos: u64, ref_end: u64, rev: bool, lead: u64, trail: u64, aln: u64) -> Align {
    Align { pos, ref_end, rev, lead, trail, aln }
}

