//! `qc --crispr`: annotate CRISPR arrays detected in the reference.
//!
//! Arrays are detected by the published `pilercr-rs` crate (MIT, a Rust port
//! of PILER-CR). This is **annotation only**: the rows never change a TIER,
//! are not finding members and are not counted in `alerts`.
//!
//! Detection is split in two so it overlaps the scan: `start_crispr` runs the
//! FASTA read and `pilercr-rs` on its own thread while the scan and the
//! circularity pass run, and `finish_crispr` joins it and maps the arrays onto
//! the scanned contigs.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::thread::JoinHandle;

use anyhow::{Context, Result};
use pilercr_rs::{Array, Config, Sequence};

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
fn map_array(seq: &Sequence, a: &Array, scan: &ScanResult) -> Result<CrisprArray, String> {
    let (ci, local, len) = array_interval(seq, a)?;
    let name = seq.contigs[ci].label.split_whitespace().next().unwrap_or("");
    let tid = scan
        .contigs
        .iter()
        .position(|c| c.name == name)
        .ok_or_else(|| format!("array at {} is on {name:?}, which is not in the alignment", a.pos))?;
    let n = a.repeats.len();
    let repeat = a.repeats.iter().map(|r| ungapped(r)).sum::<usize>() / n;
    let spacer = a.spacers[..n - 1].iter().map(Vec::len).sum::<usize>() / (n - 1).max(1);
    Ok(CrisprArray {
        tid,
        start: local as u64,
        end: (local + len) as u64,
        copies: n,
        repeat,
        spacer,
        consensus: letters(&a.consensus),
    })
}

/// Map an array's global start to `(contig index, 0-based local start, length)`.
///
/// Mirrors `report::location()` in `pilercr-rs`, which is not public in 0.1.0:
/// the same row validation, and the length is the sum of the non-gap repeat
/// characters plus the first `n - 1` spacers (the terminal right-flank spacer
/// is not part of the array). **Replace this with the upstream helper once it
/// is exposed.**
fn array_interval(seq: &Sequence, a: &Array) -> Result<(usize, usize, usize), String> {
    let ci = seq
        .contigs
        .iter()
        .position(|c| c.start <= a.pos && a.pos < c.start + c.len)
        .ok_or_else(|| format!("array at {} starts outside a contig", a.pos))?;
    let n = a.repeats.len();
    // Same checks as `report::location()` (0.1.0): one repeat row per spacer
    // and left flank, a non-empty consensus, every row as wide as the aligned
    // consensus (the array's own aligned sequences) and no all-gap repeat.
    if n == 0
        || a.spacers.len() != n
        || a.left_flanks.len() != n
        || a.aligned_consensus.is_empty()
        || a.consensus.is_empty()
        || a.repeats.iter().any(|r| r.len() != a.aligned_consensus.len() || ungapped(r) == 0)
    {
        return Err(format!("array at {} has inconsistent alignment rows", a.pos));
    }
    // Every row is taken from the reference, so it must be printable text
    // (spaces pad the left flanks at contig starts).
    for row in a
        .repeats
        .iter()
        .chain(&a.spacers)
        .chain(&a.left_flanks)
        .chain(std::iter::once(&a.consensus))
        .chain(std::iter::once(&a.aligned_consensus))
    {
        if !row.iter().all(|c| c.is_ascii_graphic() || *c == b' ') {
            return Err(format!("array at {} contains invalid sequence text", a.pos));
        }
    }
    let length = a
        .repeats
        .iter()
        .map(|r| ungapped(r))
        .chain(a.spacers[..n - 1].iter().map(Vec::len))
        .try_fold(0usize, |sum, len| sum.checked_add(len))
        .ok_or_else(|| format!("array at {} length overflow", a.pos))?;
    let local = a.pos - seq.contigs[ci].start;
    if local
        .checked_add(length)
        .is_none_or(|end| end > seq.contigs[ci].len)
    {
        return Err(format!("array at {} crosses the end of contig {}", a.pos, ci));
    }
    Ok((ci, local, length))
}

/// Non-gap characters of one aligned repeat row.
fn ungapped(row: &[u8]) -> usize {
    row.iter().filter(|&&c| c != b'-').count()
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
            repeats: repeats.iter().map(|r| r.to_vec()).collect(),
            spacers: spacers.iter().map(|s| s.to_vec()).collect(),
            left_flanks: repeats.iter().map(|_| vec![b' '; 10]).collect(),
            consensus: b"ACGT".to_vec(),
            aligned_consensus: repeats.first().map(|r| r.to_vec()).unwrap_or_default(),
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
    fn global_position_maps_to_the_right_contig() {
        // Second contig, local offset 4: global 68.
        let a = array(68, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        let (ci, local, len) = array_interval(&sequence(), &a).unwrap();
        assert_eq!((ci, local, len), (1, 4, 3 * 4 + 2 * 4));
    }

    #[test]
    fn length_counts_non_gap_repeat_characters_only() {
        // Gapped rows: 3, 2 and 4 ungapped bases; spacers 5 and 7.
        let a = array(0, &[b"AC-G", b"A--G", b"ACGT"], &[b"AAAAA", b"AAAAAAA", b"CC"]);
        let (_, _, len) = array_interval(&sequence(), &a).unwrap();
        assert_eq!(len, 3 + 2 + 4 + 5 + 7);
    }

    #[test]
    fn array_crossing_a_contig_boundary_is_rejected() {
        // Starts at local 50 of a 60 bp contig; 3 x 4 bp repeats + 2 x 4 bp
        // spacers = 20 bp crosses the end.
        let a = array(50, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        let err = array_interval(&sequence(), &a).unwrap_err();
        assert!(err.contains("crosses"), "{err}");
    }

    #[test]
    fn unmappable_arrays_are_skipped() {
        // A crossing array and an array on a contig absent from the alignment
        // both make `map_array` fail, so detection skips them (with a warning)
        // instead of failing the run.
        let crossing = array(50, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        assert!(map_array(&sequence(), &crossing, &scan(&["a", "b"])).is_err());
        let on_a = array(0, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        assert!(map_array(&sequence(), &on_a, &scan(&[])).is_err());
    }

    #[test]
    fn array_starting_in_padding_is_rejected() {
        let a = array(61, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        assert!(array_interval(&sequence(), &a).is_err());
    }

    /// `report::location()` rejects these rows; `array_interval` must too, so
    /// the mapping never derives a length from a broken alignment.
    #[test]
    fn malformed_alignment_rows_are_rejected() {
        let base = array(0, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        let ok = array_interval(&sequence(), &base).unwrap();

        // A repeat row with no ungapped bases has no length to add.
        let mut a = base.clone();
        a.repeats[1] = b"----".to_vec();
        let err = array_interval(&sequence(), &a).unwrap_err();
        assert!(err.contains("inconsistent alignment rows"), "{err}");

        // Every repeat row must be as wide as the aligned consensus.
        let mut a = base.clone();
        a.repeats[2].truncate(3);
        let err = array_interval(&sequence(), &a).unwrap_err();
        assert!(err.contains("inconsistent alignment rows"), "{err}");

        // Empty consensus or aligned consensus.
        let mut a = base.clone();
        a.consensus.clear();
        assert!(array_interval(&sequence(), &a).is_err());
        let mut a = base.clone();
        a.aligned_consensus.clear();
        assert!(array_interval(&sequence(), &a).is_err());

        // One left flank per repeat (and per spacer).
        let mut a = base.clone();
        a.left_flanks.pop();
        let err = array_interval(&sequence(), &a).unwrap_err();
        assert!(err.contains("inconsistent alignment rows"), "{err}");

        // One spacer per repeat.
        let mut a = base.clone();
        a.spacers.pop();
        assert!(array_interval(&sequence(), &a).is_err());

        // No repeats at all.
        let mut a = base.clone();
        a.repeats.clear();
        assert!(array_interval(&sequence(), &a).is_err());

        // Rows come from the reference, so control bytes are invalid text.
        let mut a = base.clone();
        a.repeats[0] = b"AC\nT".to_vec();
        let err = array_interval(&sequence(), &a).unwrap_err();
        assert!(err.contains("invalid sequence text"), "{err}");

        assert_eq!(ok, (0, 0, 3 * 4 + 2 * 4));
    }

    /// The last spacer is the right flank up to the contig end, so an array
    /// ending exactly at the contig end is valid even though that spacer would
    /// extend past it.
    #[test]
    fn terminal_right_flank_at_the_contig_end_is_valid() {
        // Contig `a` is 0..60; start at 40 with 3 x 4 bp repeats and two
        // internal 4 bp spacers = 20 bp, ending exactly at 60.
        let a = array(40, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"AAAA", b"GGGGGGGGGG"]);
        assert_eq!(array_interval(&sequence(), &a).unwrap(), (0, 40, 20));
        // Same array one base further right crosses the end.
        let a = array(41, &[b"ACGT", b"ACGT", b"ACGT"], &[b"AAAA", b"AAAA", b"GGGGGGGGGG"]);
        assert!(array_interval(&sequence(), &a).unwrap_err().contains("crosses"));
    }

    #[test]
    fn consensus_letters_drop_gaps_and_whitespace() {
        assert_eq!(letters(b"AC-GT"), "ACGT");
        assert_eq!(letters(b"A C\tG\n"), "ACG");
        let mut a = array(0, &[b"AC-G", b"ACGT", b"ACGT"], &[b"AAAA", b"CCCC", b"GGGG"]);
        a.aligned_consensus = a.repeats[1].clone();
        assert_eq!(letters(&a.consensus), letters(&a.aligned_consensus));
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
