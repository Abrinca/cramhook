//! `qc --crispr`: annotate CRISPR arrays detected in the reference.
//!
//! Arrays are detected by the published `pilercr-rs` crate (MIT, a Rust port
//! of PILER-CR). This is **annotation only**: the rows never change a TIER,
//! are not finding members and are not counted in `alerts`.
//!
//! Detection is split in two so it overlaps the scan: `start_crispr` runs the
//! FASTA read and `pilercr-rs` on its own thread while the scan and the
//! circularity pass run, and `finish_crispr` joins it and maps the arrays onto
//! the scanned contigs with upstream's `pilercr_rs::locate`.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::thread::JoinHandle;

use anyhow::{Context, Result};
use pilercr_rs::{Array, Config, Sequence, locate};

use super::Region;
use crate::note;
use crate::scan::ScanResult;

/// One array mapped onto a contig, with the values the `crispr` BED row needs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CrisprArray {
    /// Contig id in the scan (`@SQ`) order.
    pub tid: usize,
    /// 0-based, half-open.
    pub start: u64,
    pub end: u64,
    /// Number of repeats (array copies).
    pub copies: usize,
    /// Mean ungapped repeat length.
    pub repeat: usize,
    /// Integer mean spacer length.
    pub spacer: usize,
    /// Consensus letters only (gaps and whitespace removed).
    pub consensus: String,
}

impl CrisprArray {
    /// True when the 0-based position is inside the array (`start <= pos < end`).
    pub fn contains(&self, tid: usize, pos0: u32) -> bool {
        self.tid == tid && (pos0 as u64) >= self.start && (pos0 as u64) < self.end
    }

    /// The untiered `crispr` BED row (`tier none`, score 0).
    pub fn region(&self) -> Region {
        Region {
            tid: self.tid,
            start: self.start,
            end: self.end,
            label: "crispr",
            score: 0,
            tier: "none",
            reason: ".".to_string(),
            depth_ratio: None,
            hiq_share: None,
            extra: format!(
                "copies={};repeat={};spacer={};consensus={}",
                self.copies, self.repeat, self.spacer, self.consensus
            ),
            lowcov: false,
            link: None,
            clip_left: None,
            clip_right: None,
        }
    }
}

/// Start detecting arrays in the `--ref` FASTA on a thread of its own.
///
/// Only the scan-independent part runs there: the FASTA read and
/// `pilercr-rs` detection overlap the scan and the circularity pass. Mapping
/// the arrays onto contig ids needs `scan`, so `finish_crispr` does it after
/// the join.
///
/// Whole-pass failures — the FASTA cannot be read/parsed or `pilercr-rs`
/// detection itself errors — surface as an error from `finish_crispr`, which
/// `qc` treats as fatal (the caller asked for the annotation).
///
/// The FASTA is read by `Sequence::from_reader` rather than rebuilt from the
/// in-memory `Fasta`, so the contig padding rules stay upstream's.
pub(crate) fn start_crispr(ref_path: &Path) -> JoinHandle<Result<(Sequence, Vec<Array>)>> {
    let ref_path = ref_path.to_path_buf();
    std::thread::spawn(move || read_and_detect(&ref_path))
}

/// Parse `ref_path` with `pilercr-rs` and detect arrays (worker-thread body).
fn read_and_detect(ref_path: &Path) -> Result<(Sequence, Vec<Array>)> {
    let sequence = Sequence::from_reader(BufReader::new(
        File::open(ref_path).with_context(|| format!("opening FASTA {}", ref_path.display()))?,
    ))
    .with_context(|| format!("reading FASTA {}", ref_path.display()))?;
    let arrays = pilercr_rs::detect(&sequence, &Config::default())
        .map_err(|e| anyhow::anyhow!("pilercr-rs detection failed: {e}"))?;
    Ok((sequence, arrays))
}

/// Join `start_crispr` and map the detected arrays onto the scanned contigs.
///
/// Arrays that do not map cleanly (outside a contig, crossing a contig
/// boundary, malformed alignment rows, or on a contig absent from the
/// alignment) are skipped with a warning: one bad array must not fail the run.
/// A panicked worker is reported as an error, like any other whole-pass
/// failure.
pub(crate) fn finish_crispr(
    handle: JoinHandle<Result<(Sequence, Vec<Array>)>>,
    scan: &ScanResult,
) -> Result<Vec<CrisprArray>> {
    let (sequence, arrays) = handle
        .join()
        .map_err(|_| anyhow::anyhow!("CRISPR detection thread panicked"))??;
    let mut out = Vec::with_capacity(arrays.len());
    for a in &arrays {
        match map_array(&sequence, a, scan) {
            Ok(mapped) => out.push(mapped),
            Err(e) => note!("warning: crispr: {e}; skipping the array"),
        }
    }
    out.sort_by_key(|a| (a.tid, a.start));
    Ok(out)
}

/// Map one detected array onto a scanned contig.
///
/// Coordinates, copy count and the mean repeat/spacer lengths come from
/// `pilercr_rs::locate`, which validates the array (alignment rows, contig
/// bounds) the same way the crate's report writer does.
fn map_array(seq: &Sequence, a: &Array, scan: &ScanResult) -> Result<CrisprArray, String> {
    let loc = locate(seq, a)?;
    let name = seq.contigs[loc.contig]
        .label
        .split_whitespace()
        .next()
        .unwrap_or("");
    let tid = scan
        .contigs
        .iter()
        .position(|c| c.name == name)
        .ok_or_else(|| format!("array at {} is on {name:?}, which is not in the alignment", a.pos))?;
    Ok(CrisprArray {
        tid,
        start: loc.start as u64,
        end: (loc.start + loc.len) as u64,
        copies: loc.copies,
        repeat: loc.repeat_len,
        spacer: loc.spacer_len,
        consensus: letters(&a.consensus),
    })
}

/// Consensus letters only: gaps and whitespace removed.
fn letters(seq: &[u8]) -> String {
    seq.iter()
        .filter(|c| c.is_ascii_alphabetic())
        .map(|&c| c as char)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use pilercr_rs::ArrayElement;
    use pilercr_rs::seq::Contig;

    /// Two contigs with padding between them: `a` at 0..60, `b` at 64..104.
    fn sequence() -> Sequence {
        Sequence {
            bases: vec![b'A'; 104],
            contigs: vec![
                Contig { label: "a first contig".into(), start: 0, len: 60 },
                Contig { label: "b".into(), start: 64, len: 40 },
            ],
        }
    }

    fn array(pos: usize, repeats: &[&[u8]], spacers: &[&[u8]]) -> Array {
        Array {
            pos,
            consensus: b"ACGT".to_vec(),
            aligned_consensus: repeats.first().map(|r| r.to_vec()).unwrap_or_default(),
            elements: repeats
                .iter()
                .zip(spacers)
                .map(|(repeat, spacer)| ArrayElement {
                    left_flank: vec![b' '; 10],
                    repeat: repeat.to_vec(),
                    spacer: spacer.to_vec(),
                })
                .collect(),
            ..Array::default()
        }
    }

    fn scan(names: &[&str]) -> ScanResult {
        ScanResult {
            contigs: names
                .iter()
                .map(|n| crate::scan::ContigCounts {
                    name: n.to_string(),
                    len: 0,
                    bases: Vec::new(),
                    cov: Vec::new(),
                    clips: Vec::new(),
                    edge_clips: 0,
                    aligned_bases: 0,
                })
                .collect(),
            indels: Vec::new(),
            large: Vec::new(),
            records: Default::default(),
            chunks: 1,
            mean_aligned_span: 0.0,
        }
    }

    #[test]
    fn unmappable_arrays_are_skipped() {
        // Arrays that `locate` rejects (or that sit on a contig absent from the
        // alignment) make `map_array` fail, so detection skips them (with a
        // warning) instead of failing the run. Crossing: local 50 of a 60 bp
        // contig with 3 x 4 bp repeats and 2 x 4 bp spacers = 20 bp. Padding:
        // 60..64 lies between contigs `a` and `b`.
        let crossing = array(50, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        assert!(map_array(&sequence(), &crossing, &scan(&["a", "b"])).is_err());
        let padding = array(61, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        assert!(map_array(&sequence(), &padding, &scan(&["a", "b"])).is_err());
        // A mappable array on a contig absent from the alignment.
        let on_a = array(0, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        assert!(map_array(&sequence(), &on_a, &scan(&[])).is_err());
    }

    #[test]
    fn inconsistent_alignment_rows_are_skipped() {
        // One repeat row narrower than the aligned consensus: `locate` rejects
        // it, so the array is skipped rather than mapped with a wrong length.
        let mut a = array(0, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        a.elements[2].repeat.truncate(3);
        assert!(map_array(&sequence(), &a, &scan(&["a", "b"])).is_err());
    }

    #[test]
    fn consensus_letters_drop_gaps_and_whitespace() {
        assert_eq!(letters(b"AC-GT"), "ACGT");
        assert_eq!(letters(b"A C\tG\n"), "ACG");
    }

    #[test]
    fn global_position_maps_to_the_right_contig() {
        // Second contig, local offset 4: global 68.
        let a = array(68, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        let mapped = map_array(&sequence(), &a, &scan(&["a", "b"])).unwrap();
        assert_eq!((mapped.tid, mapped.start, mapped.end), (1, 4, 24));
    }

    #[test]
    fn mapped_array_reports_mean_spacer_and_copy_values() {
        let a = array(
            64,
            &[b"AC-GT", b"ACGTT", b"ACGTT", b"ACGTT"],
            &[b"AAAA", b"CCCCCC", b"GGGGG", b"TTT"],
        );
        let mapped = map_array(&sequence(), &a, &scan(&["a", "b"])).unwrap();
        assert_eq!(mapped.tid, 1);
        // 4 + 5 + 5 + 5 repeat bases + 4 + 6 + 5 spacer bases = 34.
        assert_eq!((mapped.start, mapped.end), (0, 34));
        assert_eq!(mapped.copies, 4);
        assert_eq!(mapped.repeat, 4); // (4 + 5 + 5 + 5) / 4
        assert_eq!(mapped.spacer, 5); // (4 + 6 + 5) / 3
        assert_eq!(mapped.consensus, "ACGT");
    }

    #[test]
    fn mapped_array_uses_the_first_header_token() {
        let a = array(0, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        let mapped = map_array(&sequence(), &a, &scan(&["a", "b"])).unwrap();
        assert_eq!((mapped.tid, mapped.start, mapped.end), (0, 0, 20));
    }

    #[test]
    fn contains_is_half_open() {
        let a = CrisprArray {
            tid: 1,
            start: 10,
            end: 20,
            copies: 3,
            repeat: 4,
            spacer: 2,
            consensus: "ACGT".into(),
        };
        assert!(!a.contains(0, 15));
        assert!(!a.contains(1, 9));
        assert!(a.contains(1, 10));
        assert!(a.contains(1, 19));
        assert!(!a.contains(1, 20));
        assert_eq!(a.region().tier, "none");
        assert_eq!(a.region().score, 0);
    }
}
