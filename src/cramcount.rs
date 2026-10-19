//! Record counts of an alignment file without decoding it, for the
//! `--progress` totals: CRAM container headers, or BAM index statistics.
//!
//! Best effort: any unsupported version, parse error or truncation gives
//! `None`, never an error.

use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use rust_htslib::bam::{self, Read as _};
use rust_htslib::htslib;

/// Every record in the file (primary, supplementary, secondary, unmapped;
/// `samtools view -c`): from the container headers of a CRAM 3.x, or the index
/// statistics of an indexed BAM.
pub fn record_count(path: &Path) -> Option<u64> {
    let mut magic = [0u8; 4];
    File::open(path).ok()?.read_exact(&mut magic).ok()?;
    if &magic == b"CRAM" {
        cram_record_count(path)
    } else if path.extension().is_some_and(|e| e.eq_ignore_ascii_case("bam")) {
        bam_record_count(path)
    } else {
        None
    }
}

/// Records of an indexed BAM: mapped + unmapped per reference (the index's
/// optional statistics pseudo-bin) plus the unplaced unmapped count.
///
/// A BAI may lack the pseudo-bins; `hts_idx_get_stat` then fails, and so it
/// does for a reference without any data. The two are told apart by asking the
/// index for the reference's chunks: data without statistics gives `None`.
/// `hts_idx_get_n_no_coor` cannot signal absence (htslib reads a missing
/// trailing count as 0). Every writer that produces the pseudo-bins (htslib,
/// samtools, picard) also writes that count, so it is trusted once the
/// per-reference statistics are present; should it still be short, the
/// `--progress` guard switches the phase to an estimate when more records are
/// read than counted.
fn bam_record_count(path: &Path) -> Option<u64> {
    let r = bam::IndexedReader::from_path(path).ok()?;
    let idx = r.index().inner_ptr();
    let mut total = 0u64;
    for tid in 0..r.header().target_count() {
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

/// Whether the index holds any chunk for reference `tid` (`None`: unknown).
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

/// Start position of the CRAM EOF container (`"EOF"` as itf8).
const EOF_START: i32 = 4_542_278;

/// Records of a CRAM 3.x file from its container headers, without decoding.
pub fn cram_record_count(path: &Path) -> Option<u64> {
    let f = File::open(path).ok()?;
    let len = f.metadata().ok()?.len();
    let mut r = BufReader::new(f);
    let mut def = [0u8; 26];
    r.read_exact(&mut def).ok()?;
    if &def[..4] != b"CRAM" || def[4] != 3 {
        return None;
    }
    let mut pos = 26u64;
    let mut total = 0u64;
    let mut first = true;
    let mut saw_eof = false;
    while pos < len {
        let mut b4 = [0u8; 4];
        r.read_exact(&mut b4).ok()?;
        let data_len = i32::from_le_bytes(b4);
        if data_len < 0 {
            return None;
        }
        let mut hdr = Counted { r: &mut r, n: 4 };
        let ref_id = itf8(&mut hdr)?;
        let start = itf8(&mut hdr)?;
        let _span = itf8(&mut hdr)?;
        let nrec = itf8(&mut hdr)?;
        ltf8(&mut hdr)?; // record counter
        ltf8(&mut hdr)?; // bases
        itf8(&mut hdr)?; // blocks
        let landmarks = itf8(&mut hdr)?;
        if landmarks < 0 {
            return None;
        }
        for _ in 0..landmarks {
            itf8(&mut hdr)?;
        }
        let mut crc = [0u8; 4];
        hdr.read_exact(&mut crc)?;
        let hdr_len = hdr.n;
        if nrec < 0 {
            return None;
        }
        pos += hdr_len + data_len as u64;
        if pos > len {
            return None; // truncated
        }
        if !first {
            total += nrec as u64;
        }
        first = false;
        saw_eof = nrec == 0 && ref_id == -1 && start == EOF_START;
        r.seek(SeekFrom::Start(pos)).ok()?;
    }
    // The last container must be the EOF container (some writers leave more
    // than one): a file cut at a container boundary would otherwise look
    // complete but count short.
    if !saw_eof {
        return None;
    }
    Some(total)
}

/// Reader that counts the bytes consumed.
struct Counted<'a, R: Read> {
    r: &'a mut R,
    n: u64,
}

impl<R: Read> Counted<'_, R> {
    fn byte(&mut self) -> Option<u8> {
        let mut b = [0u8; 1];
        self.read_exact(&mut b)?;
        Some(b[0])
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> Option<()> {
        self.r.read_exact(buf).ok()?;
        self.n += buf.len() as u64;
        Some(())
    }
}

/// CRAM ITF8 integer (1–5 bytes), as a signed 32-bit value.
fn itf8<R: Read>(r: &mut Counted<'_, R>) -> Option<i32> {
    let b0 = r.byte()? as u32;
    let v = if b0 < 0x80 {
        b0
    } else if b0 < 0xC0 {
        ((b0 & 0x3F) << 8) | r.byte()? as u32
    } else if b0 < 0xE0 {
        let mut b = [0u8; 2];
        r.read_exact(&mut b)?;
        ((b0 & 0x1F) << 16) | (b[0] as u32) << 8 | b[1] as u32
    } else if b0 < 0xF0 {
        let mut b = [0u8; 3];
        r.read_exact(&mut b)?;
        ((b0 & 0x0F) << 24) | (b[0] as u32) << 16 | (b[1] as u32) << 8 | b[2] as u32
    } else {
        let mut b = [0u8; 4];
        r.read_exact(&mut b)?;
        ((b0 & 0x0F) << 28) | (b[0] as u32) << 20 | (b[1] as u32) << 12 | (b[2] as u32) << 4 | (b[3] as u32 & 0x0F)
    };
    Some(v as i32)
}

/// CRAM LTF8 integer (1–9 bytes); only skipped here.
fn ltf8<R: Read>(r: &mut Counted<'_, R>) -> Option<()> {
    let b0 = r.byte()?;
    let extra = b0.leading_ones() as usize;
    let mut b = [0u8; 8];
    r.read_exact(&mut b[..extra])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::process::Command;

    fn root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
    }

    /// `samtools view -c`, or `None` when samtools is not available.
    fn samtools_count(p: &Path) -> Option<u64> {
        let o = Command::new("samtools").args(["view", "-c"]).arg(p).output().ok()?;
        if !o.status.success() {
            return None;
        }
        String::from_utf8_lossy(&o.stdout).trim().parse().ok()
    }

    fn fixtures() -> Vec<PathBuf> {
        let r = root();
        let mut v = vec![
            r.join("tests/fixtures/clip_join/oracle.cram"),
            r.join("data/synthetic/S1/hifi/old.cram"),
            r.join("data/synthetic/S1/ont/old.cram"),
            r.join("testdata/data/cramhook-testdata-v1/ont/cons.cram"),
            r.join("testdata/data/cramhook-testdata-v1/hifi/asm.cram"),
            r.join("data/v2/real/ont_FAM29080__cons-to-polished/old.cram"),
        ];
        v.retain(|p| p.exists());
        v
    }

    #[test]
    fn totals_match_samtools() {
        let files = fixtures();
        assert!(files.iter().any(|p| p.starts_with(root().join("tests/fixtures"))), "{files:?}");
        let mut checked = 0;
        for p in &files {
            let Some(want) = samtools_count(p) else { continue };
            assert_eq!(record_count(p), Some(want), "{}", p.display());
            checked += 1;
        }
        assert!(checked >= 2, "samtools missing? checked {checked} files");
    }

    #[test]
    fn bam_from_index_stats() {
        let dir = root().join("target/test-cramcount");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("clip_join.bam");
        std::fs::copy(root().join("tests/fixtures/clip_join/old.bam"), &bam).unwrap();
        // Without an index there is no total.
        let _ = std::fs::remove_file(dir.join("clip_join.bam.bai"));
        assert_eq!(record_count(&bam), None);
        bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
        let n = record_count(&bam).expect("indexed BAM");
        assert!(n > 0);
        if let Some(want) = samtools_count(&bam) {
            assert_eq!(n, want);
        }
    }

    /// Remove the statistics pseudo-bin (bin 37450) from every reference of a
    /// BAI, keeping everything else (a valid index, as some writers produce).
    fn strip_bai_stats(b: &[u8]) -> Vec<u8> {
        let i32_at = |o: usize| i32::from_le_bytes(b[o..o + 4].try_into().unwrap());
        let mut out = b[..8].to_vec();
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
            out.extend_from_slice(&kept.to_le_bytes());
            out.extend_from_slice(&bins);
            let n_intv = i32_at(off) as usize;
            out.extend_from_slice(&b[off..off + 4 + 8 * n_intv]);
            off += 4 + 8 * n_intv;
        }
        out.extend_from_slice(&b[off..]); // n_no_coor
        out
    }

    #[test]
    fn bai_without_stats_is_none() {
        let dir = root().join("target/test-cramcount");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("no_stats.bam");
        std::fs::copy(root().join("tests/fixtures/clip_join/old.bam"), &bam).unwrap();
        let bai = dir.join("no_stats.bam.bai");
        let _ = std::fs::remove_file(&bai);
        bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
        let full = std::fs::read(&bai).unwrap();
        assert!(record_count(&bam).unwrap() > 0);
        let stripped = strip_bai_stats(&full);
        assert!(stripped.len() < full.len(), "no stats bins found to strip");
        std::fs::write(&bai, &stripped).unwrap();
        // Still a usable index, but no statistics: no total (not Some(0)).
        assert!(bam::IndexedReader::from_path(&bam).is_ok());
        assert_eq!(record_count(&bam), None);
    }

    #[test]
    fn bam_with_an_empty_reference_counts() {
        // A reference without any record has no statistics bin either; that
        // must not turn the count into `None`.
        use rust_htslib::bam::header::HeaderRecord;
        let dir = root().join("target/test-cramcount");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("empty_ref.bam");
        let mut h = bam::Header::new();
        h.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
        h.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "a").push_tag(b"LN", 1000));
        h.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "empty").push_tag(b"LN", 1000));
        {
            let mut w = bam::Writer::from_path(&bam, &h, bam::Format::Bam).unwrap();
            for (i, tid) in [(0, 0), (1, 0), (2, -1)] {
                let mut rec = bam::Record::new();
                let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(4)]);
                let placed = tid >= 0;
                rec.set(format!("r{i}").as_bytes(), placed.then_some(&cigar), b"ACGT", &[30; 4]);
                rec.set_tid(tid);
                rec.set_pos(if placed { 10 * i } else { -1 });
                rec.set_mtid(-1);
                rec.set_mpos(-1);
                rec.set_flags(if placed { 0 } else { 4 });
                w.write(&rec).unwrap();
            }
        }
        bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
        assert_eq!(record_count(&bam), Some(3));
    }

    #[test]
    fn truncated_and_garbage_are_none() {
        let src = root().join("tests/fixtures/clip_join/oracle.cram");
        let bytes = std::fs::read(&src).unwrap();
        let dir = root().join("target/test-cramcount");
        std::fs::create_dir_all(&dir).unwrap();
        // Cut inside a container, and cut exactly before the EOF container.
        let cut = dir.join("cut.cram");
        std::fs::write(&cut, &bytes[..bytes.len() / 2]).unwrap();
        assert_eq!(cram_record_count(&cut), None);
        // Drop every trailing EOF container (this fixture has several).
        let eof = &bytes[bytes.len() - 38..];
        let mut end = bytes.len();
        while bytes[..end].ends_with(eof) {
            end -= 38;
        }
        let no_eof = dir.join("no_eof.cram");
        std::fs::write(&no_eof, &bytes[..end]).unwrap();
        assert_eq!(cram_record_count(&no_eof), None);
        let junk = dir.join("junk.cram");
        let mut g = b"CRAM\x03\x01".to_vec();
        g.extend(std::iter::repeat_n(0xffu8, 100));
        std::fs::write(&junk, &g).unwrap();
        assert_eq!(cram_record_count(&junk), None);
        let v2 = dir.join("v2.cram");
        let mut b = bytes.clone();
        b[4] = 2;
        std::fs::write(&v2, &b).unwrap();
        assert_eq!(cram_record_count(&v2), None);
        let empty = dir.join("empty.cram");
        std::fs::write(&empty, b"").unwrap();
        assert_eq!(record_count(&empty), None);
    }
}
