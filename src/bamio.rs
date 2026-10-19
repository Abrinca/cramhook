//! `bamio`: shared alignment-file I/O and CIGAR helpers.
//!
//! The `is_cram`/`open_indexed` pair opens an indexed BAM/CRAM for the several
//! passes (`scan`, `qc::circularity`, `qc::link`, `fish`); each caller passes the
//! exact CRAM required-field mask it decodes. The CIGAR helpers (`op_parts`, the
//! `OP_*` codes, `terminal_clips`) are shared by every pass that walks a CIGAR.

use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};
use rust_htslib::htslib;

/// CRAM required fields: identity and position only.
const FIELDS_ID: u32 = htslib::sam_fields_SAM_QNAME
    | htslib::sam_fields_SAM_FLAG
    | htslib::sam_fields_SAM_RNAME
    | htslib::sam_fields_SAM_POS;

/// CRAM required fields: identity, position and mapping, without SEQ/QUAL/AUX
/// (the `qc::circularity` mask).
pub(crate) const FIELDS_BASE: u32 =
    FIELDS_ID | htslib::sam_fields_SAM_MAPQ | htslib::sam_fields_SAM_CIGAR;

/// CRAM required fields for `scan`: mapping columns plus SEQ.
pub(crate) const FIELDS_SCAN: u32 = FIELDS_BASE | htslib::sam_fields_SAM_SEQ;

/// CRAM required fields for `qc::link`: mapping columns, SEQ and AUX.
pub(crate) const FIELDS_LINK: u32 =
    FIELDS_BASE | htslib::sam_fields_SAM_SEQ | htslib::sam_fields_SAM_AUX;

/// CRAM required fields for `fish`'s payload pass: identity, position, SEQ and
/// QUAL (no CIGAR/AUX/MAPQ).
pub(crate) const FIELDS_FISH_MIN: u32 =
    FIELDS_ID | htslib::sam_fields_SAM_SEQ | htslib::sam_fields_SAM_QUAL;

/// CRAM required fields for `fish`'s full pass: mapping columns, SEQ, QUAL, AUX.
pub(crate) const FIELDS_FISH_FULL: u32 = FIELDS_BASE
    | htslib::sam_fields_SAM_SEQ
    | htslib::sam_fields_SAM_QUAL
    | htslib::sam_fields_SAM_AUX;

/// `true` when the path looks like CRAM (case-insensitive `.cram`).
pub(crate) fn is_cram(p: &Path) -> bool {
    p.extension().is_some_and(|e| e.eq_ignore_ascii_case("cram"))
}

/// Open an indexed reader with an explicit reference. For CRAM only `fields`
/// are decoded; each caller passes its own [`FIELDS_*`](FIELDS_SCAN) mask.
pub(crate) fn open_indexed(
    input: &Path,
    ref_path: &Path,
    fields: u32,
) -> Result<bam::IndexedReader> {
    let mut r = bam::IndexedReader::from_path(input)
        .with_context(|| format!("opening {} (an index is required)", input.display()))?;
    r.set_reference(ref_path)
        .with_context(|| format!("setting reference {} for {}", ref_path.display(), input.display()))?;
    if is_cram(input) {
        r.set_cram_options(htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS, fields)
            .context("setting CRAM required fields")?;
    }
    Ok(r)
}

/// Same record filter as `scan`: skip secondary, QC-fail, duplicate and
/// unmapped; keep primary and supplementary.
pub(crate) fn keep(rec: &bam::Record) -> bool {
    let f = rec.flags();
    f & 0x100 == 0 && f & 0x200 == 0 && f & 0x400 == 0 && f & 0x4 == 0
}

/// Read an integer aux tag (`s1`/`s2`/`NM`); a missing tag or a non-numeric
/// type is `None`. Float/double values are truncated to `i64`.
pub(crate) fn aux_int(rec: &bam::Record, tag: &[u8]) -> Option<i64> {
    match rec.aux(tag).ok()? {
        Aux::I8(v) => Some(v as i64),
        Aux::U8(v) => Some(v as i64),
        Aux::I16(v) => Some(v as i64),
        Aux::U16(v) => Some(v as i64),
        Aux::I32(v) => Some(v as i64),
        Aux::U32(v) => Some(v as i64),
        Aux::Float(v) => Some(v as i64),
        Aux::Double(v) => Some(v as i64),
        _ => None,
    }
}

#[inline]
pub(crate) fn op_parts(c: u32) -> (u32, usize) {
    (c & 0xf, (c >> 4) as usize)
}

pub(crate) const OP_M: u32 = 0;
pub(crate) const OP_I: u32 = 1;
pub(crate) const OP_D: u32 = 2;
pub(crate) const OP_N: u32 = 3;
pub(crate) const OP_S: u32 = 4;
pub(crate) const OP_H: u32 = 5;
pub(crate) const OP_EQ: u32 = 7;
pub(crate) const OP_X: u32 = 8;

#[inline]
fn clip_op(c: u32) -> bool {
    matches!(op_parts(c).0, OP_S | OP_H)
}

/// Sum of the terminal `S`/`H` ops at each end of a CIGAR: `(leading, trailing)`.
/// The terminal-clip rule `scan` counts (`5H60S100M` has 65 leading).
pub(crate) fn terminal_clips(cig: &[u32]) -> (u32, u32) {
    let lead: u32 = cig.iter().take_while(|&&c| clip_op(c)).map(|&c| op_parts(c).1 as u32).sum();
    let trail: u32 = cig.iter().rev().take_while(|&&c| clip_op(c)).map(|&c| op_parts(c).1 as u32).sum();
    (lead, trail)
}

/// [`terminal_clips`] summed in `u64` (the circularity pass's arithmetic; a
/// non-canonical CIGAR with many maximal clip ops cannot overflow).
pub(crate) fn terminal_clips_u64(cig: &[u32]) -> (u64, u64) {
    let lead: u64 = cig.iter().take_while(|&&c| clip_op(c)).map(|&c| op_parts(c).1 as u64).sum();
    let trail: u64 = cig.iter().rev().take_while(|&&c| clip_op(c)).map(|&c| op_parts(c).1 as u64).sum();
    (lead, trail)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn terminal_clips_u64_does_not_overflow() {
        let max_h = (((1u32 << 28) - 1) << 4) | OP_H;
        let mut cig = vec![max_h; 17];
        cig.push((1 << 4) | OP_M);
        assert_eq!(terminal_clips_u64(&cig), (17 * ((1u64 << 28) - 1), 0));
        assert_eq!(terminal_clips(&[(5 << 4) | OP_H, (60 << 4) | OP_S, (100 << 4) | OP_M]), (65, 0));
    }
}
