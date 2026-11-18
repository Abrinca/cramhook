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

/// `true` when `reader` is reading a BAM file (not CRAM/SAM).
fn is_bam_reader(reader: &bam::IndexedReader) -> bool {
    // SAFETY: `htsfile()` returns a valid pointer for the lifetime of `reader`.
    unsafe { (*reader.htsfile()).format.format == htslib::htsExactFormat_bam }
}

/// Whether the BAM index holds any chunk for reference `tid` (`None`: unknown).
fn reference_has_data(idx: *const htslib::hts_idx_t, tid: i32) -> Option<bool> {
    const HTS_POS_MAX: i64 = ((i32::MAX as i64) << 32) | i32::MAX as i64;
    // SAFETY: `idx` is a live index; the iterator is destroyed before return.
    unsafe {
        let it = htslib::sam_itr_queryi(idx, tid, 0, HTS_POS_MAX);
        if it.is_null() {
            return None;
        }
        let has = (*it).n_off > 0;
        htslib::hts_itr_destroy(it);
        Some(has)
    }
}

/// Total record count reported by a BAM index (`mapped + unmapped` across all
/// references plus `n_no_coor`), or `None` when not a BAM, the index is
/// missing, or a reference has chunks without statistics bins.
pub(crate) fn bam_index_record_count(reader: &bam::IndexedReader) -> Option<u64> {
    if !is_bam_reader(reader) {
        return None;
    }
    let idx = reader.index().inner_ptr();
    if idx.is_null() {
        return None;
    }
    let mut total = 0u64;
    for tid in 0..reader.header().target_count() {
        let tid = i32::try_from(tid).ok()?;
        let (mut mapped, mut unmapped) = (0u64, 0u64);
        // SAFETY: `idx` is the reader's live index; the out-pointers are valid.
        let rc = unsafe { htslib::hts_idx_get_stat(idx, tid, &mut mapped, &mut unmapped) };
        if rc < 0 {
            if reference_has_data(idx, tid)? {
                return None; // statistics absent
            }
            continue;
        }
        total += mapped + unmapped;
    }
    // SAFETY: as above.
    Some(total + unsafe { htslib::hts_idx_get_n_no_coor(idx) })
}

/// Seek to the coordinate-less (`*` / `FetchDefinition::Unmapped`) records of
/// `reader`. Returns `Ok(true)` if the caller should read records, or
/// `Ok(false)` when `fetch` fails only because `reader` is a BAM whose index
/// reports zero records in total (htslib cannot seek to `*` on an empty BAI).
/// Any other fetch error (CRAM, non-empty BAM, or an index without stats) is
/// propagated.
pub(crate) fn fetch_unplaced(reader: &mut bam::IndexedReader) -> Result<bool> {
    match reader.fetch(bam::FetchDefinition::Unmapped) {
        Ok(()) => Ok(true),
        Err(err) => {
            if bam_index_record_count(reader) == Some(0) {
                Ok(false)
            } else {
                Err(err.into())
            }
        }
    }
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

    fn write_test_bam(path: &Path, recs: &[(i32, i64, u16)]) {
        use rust_htslib::bam::header::HeaderRecord;
        let mut h = bam::Header::new();
        h.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
        h.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", 100));
        {
            let mut w = bam::Writer::from_path(path, &h, bam::Format::Bam).unwrap();
            for (i, &(tid, pos, flags)) in recs.iter().enumerate() {
                let mut rec = bam::Record::new();
                let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(4)]);
                let has_cigar = tid >= 0 && pos >= 0 && flags & 0x4 == 0;
                rec.set(
                    format!("r{i}").as_bytes(),
                    has_cigar.then_some(&cigar),
                    b"ACGT",
                    &[30; 4],
                );
                rec.set_tid(tid);
                rec.set_pos(pos);
                rec.set_mtid(-1);
                rec.set_mpos(-1);
                rec.set_flags(flags);
                w.write(&rec).unwrap();
            }
        }
        bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
    }

    #[test]
    fn fetch_unplaced_handles_empty_mapped_unplaced_and_missing_stats() {
        let dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test-bamio");
        std::fs::create_dir_all(&dir).unwrap();

        // Empty BAM (0 records): returns Ok(false) rather than failing.
        let empty = dir.join("empty.bam");
        write_test_bam(&empty, &[]);
        let mut r = bam::IndexedReader::from_path(&empty).unwrap();
        assert_eq!(bam_index_record_count(&r), Some(0));
        assert!(!fetch_unplaced(&mut r).unwrap());

        // Mapped-only BAM: fetch succeeds and yields 0 unplaced records.
        let mapped = dir.join("mapped.bam");
        write_test_bam(&mapped, &[(0, 10, 0)]);
        let mut r = bam::IndexedReader::from_path(&mapped).unwrap();
        assert_eq!(bam_index_record_count(&r), Some(1));
        assert!(fetch_unplaced(&mut r).unwrap());
        let mut rec = bam::Record::new();
        assert!(r.read(&mut rec).is_none());

        // Unplaced-only BAM: fetch succeeds and yields the unplaced record.
        let unplaced = dir.join("unplaced.bam");
        write_test_bam(&unplaced, &[(-1, -1, 4)]);
        let mut r = bam::IndexedReader::from_path(&unplaced).unwrap();
        assert_eq!(bam_index_record_count(&r), Some(1));
        assert!(fetch_unplaced(&mut r).unwrap());
        assert!(r.read(&mut rec).transpose().unwrap().is_some());
        assert!(r.read(&mut rec).is_none());

        // Non-empty BAM whose BAI lacks bin 37450 (META_BIN): fetch fails and
        // the error must propagate rather than treating the BAM as empty.
        let no_meta = dir.join("no_meta.bam");
        write_test_bam(&no_meta, &[(0, 10, 0)]);
        let bai = dir.join("no_meta.bam.bai");
        let b = std::fs::read(&bai).unwrap();
        let i32_at = |o: usize| i32::from_le_bytes(b[o..o + 4].try_into().unwrap());
        let mut stripped = b[..8].to_vec();
        let mut off = 8;
        for _ in 0..i32_at(4) {
            let n_bin = i32_at(off);
            off += 4;
            let mut bins = Vec::new();
            let mut kept = 0i32;
            for _ in 0..n_bin {
                let start = off;
                let bin = u32::from_le_bytes(b[off..off + 4].try_into().unwrap());
                let n_chunk = i32_at(off + 4) as usize;
                off += 8 + 16 * n_chunk;
                if bin != 37450 {
                    bins.extend_from_slice(&b[start..off]);
                    kept += 1;
                }
            }
            stripped.extend_from_slice(&kept.to_le_bytes());
            stripped.extend_from_slice(&bins);
            let n_intv = i32_at(off) as usize;
            stripped.extend_from_slice(&b[off..off + 4 + 8 * n_intv]);
            off += 4 + 8 * n_intv;
        }
        stripped.extend_from_slice(&b[off..]);
        std::fs::write(&bai, stripped).unwrap();
        let mut r = bam::IndexedReader::from_path(&no_meta).unwrap();
        assert_eq!(bam_index_record_count(&r), None);
        assert!(fetch_unplaced(&mut r).is_err());
    }
}
