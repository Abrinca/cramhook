//! Integration test for `qc --crispr`: arrays detected in the reference are
//! annotated only (BED rows, VCF INFO flag, JSON counts); TIERs, findings and
//! alerts are identical with and without the flag.

use std::path::{Path, PathBuf};
use std::process::Command;

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{self, Format, Header, Record, Writer};

use cramhook::{Preset, QcParams, run_qc};

const C1: usize = 3_000;
const C2: usize = 4_000;
const C3: usize = 2_500;
const ARRAY_START: usize = 1_000;
const REPEAT: &[u8] = b"GTACCGATGCTAGCTAGGATCCGATCGATCGTAGCT"; // 36 bp
const SPACERS: [&[u8]; 3] = [
    b"GATACCGTCTAGTGCATACGATGCCTATCG",
    b"CTTGATGACCTAGACTACGGATCATGCTAG",
    b"TAGCTACGATGCTGATCGTAGCTAGCATCG",
];
const SNP_C2: usize = 1_120; // inside a spacer of the array
const SNP_C1: usize = 1_100; // no array on c1
/// 4 x 36 bp repeats + 3 x 30 bp spacers.
const ARRAY_LEN: usize = 4 * 36 + 3 * 30;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

fn workdir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/crispr-impl/test").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn random_ref(len: usize, seed: u64) -> Vec<u8> {
    let mut rng = Lcg(seed);
    (0..len).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect()
}

fn alt_of(b: u8) -> u8 {
    match b {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

/// The three-contig reference; `c2` carries the 4-copy array.
fn reference() -> [Vec<u8>; 3] {
    assert_eq!(REPEAT.len(), 36);
    assert!(SPACERS.iter().all(|s| s.len() == 30));
    let mut c1 = random_ref(C1, 11);
    let mut c2 = random_ref(C2, 12);
    let c3 = random_ref(C3, 13);
    let mut p = ARRAY_START;
    for (i, spacer) in SPACERS.iter().enumerate() {
        c2[p..p + REPEAT.len()].copy_from_slice(REPEAT);
        p += REPEAT.len();
        c2[p..p + spacer.len()].copy_from_slice(spacer);
        p += spacer.len();
        if i == SPACERS.len() - 1 {
            c2[p..p + REPEAT.len()].copy_from_slice(REPEAT);
        }
    }
    c1[SNP_C1] = b'A'; // the reads carry C
    c2[SNP_C2] = b'A'; // the reads carry C
    [c1, c2, c3]
}

fn write_fasta(dir: &Path, refs: &[Vec<u8>; 3]) -> PathBuf {
    write_fasta_styled(dir, "ref.fa", refs, false, true)
}

/// Write `refs` with `\r\n` line endings and/or without a final newline.
fn write_fasta_styled(dir: &Path, file: &str, refs: &[Vec<u8>; 3], crlf: bool, final_newline: bool) -> PathBuf {
    let fa = dir.join(file);
    let eol = if crlf { "\r\n" } else { "\n" };
    let mut text = String::new();
    let mut fai = String::new();
    for (name, seq) in ["c1", "c2", "c3"].iter().zip(refs) {
        let header = format!(">{name}{eol}");
        let offset = text.len() + header.len();
        text.push_str(&header);
        for chunk in seq.chunks(60) {
            text.push_str(std::str::from_utf8(chunk).unwrap());
            text.push_str(eol);
        }
        fai.push_str(&format!("{name}\t{}\t{offset}\t60\t{}\n", seq.len(), 60 + eol.len()));
    }
    if !final_newline {
        text.truncate(text.len() - eol.len());
    }
    std::fs::write(&fa, text).unwrap();
    std::fs::write(dir.join(format!("{file}.fai")), fai).unwrap();
    fa
}

struct Rd {
    tid: usize,
    pos: usize,
    seq: Vec<u8>,
    rev: bool,
    /// `None`: one `M` over the read.
    cigar: Option<CigarString>,
}

/// 12 reads (alternating strands, 500 bp) covering `[900, 1400)`; all carry
/// the alt at `snp`.
fn reads_for(tid: usize, seq: &[u8], snp: usize) -> Vec<Rd> {
    reads_at(tid, seq, &[snp])
}

/// `reads_for` with one alt per 0-based reference position in `sites`.
fn reads_at(tid: usize, seq: &[u8], sites: &[usize]) -> Vec<Rd> {
    const START: usize = 900;
    (0..12)
        .map(|i| {
            let mut read = seq[START..START + 500].to_vec();
            for &site in sites {
                assert!((START..START + 500).contains(&site), "site {site} is outside the read");
                read[site - START] = alt_of(read[site - START]);
            }
            Rd { tid, pos: START, seq: read, rev: i % 2 == 1, cigar: None }
        })
        .collect()
}

fn write_bam(dir: &Path, refs: &[Vec<u8>; 3], mut reads: Vec<Rd>) -> PathBuf {
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    for (name, seq) in ["c1", "c2", "c3"].iter().zip(refs) {
        header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", *name).push_tag(b"LN", seq.len()));
    }
    let bam_path = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
    reads.sort_by_key(|r| (r.tid, r.pos));
    for (i, r) in reads.iter().enumerate() {
        let mut rec = Record::new();
        let qual = vec![30u8; r.seq.len()];
        let cigar = r.cigar.clone().unwrap_or_else(|| CigarString(vec![Cigar::Match(r.seq.len() as u32)]));
        rec.set(format!("r{i}").as_bytes(), Some(&cigar), &r.seq, &qual);
        rec.set_tid(r.tid as i32);
        rec.set_pos(r.pos as i64);
        rec.set_mapq(60);
        rec.set_flags(if r.rev { 0x10 } else { 0 });
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        w.write(&rec).unwrap();
    }
    drop(w);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
    bam_path
}

fn params(crispr: bool) -> QcParams {
    QcParams { preset: Some(Preset::MapHifi), threads: 2, crispr, ..Default::default() }
}

fn vcf_data(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.to_string())
        .collect()
}

fn json_without_volatile(path: &Path) -> serde_json::Value {
    let mut v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let o = v.as_object_mut().unwrap();
    o.remove("timing_s");
    o.remove("threads");
    o.remove("chunks");
    o.remove("inputs");
    o.remove("crispr");
    if let Some(l) = o.get_mut("linkage").and_then(|l| l.as_object_mut()) {
        l.remove("timing_s");
    }
    // The flag itself is the only expected difference in `params`.
    if let Some(p) = o.get_mut("params").and_then(|p| p.as_object_mut()) {
        p.remove("crispr");
    }
    v
}

fn bed_rows(path: &Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').map(|s| s.to_string()).collect())
        .collect()
}

#[test]
fn crispr_annotates_only() {
    let dir = workdir("annotate");
    let refs = reference();
    let fa = write_fasta(&dir, &refs);
    let mut reads = reads_for(0, &refs[0], SNP_C1);
    reads.extend(reads_for(1, &refs[1], SNP_C2));
    let bam = write_bam(&dir, &refs, reads);

    let on = dir.join("on/out");
    let off = dir.join("off/out");
    let s_on = run_qc(&bam, &fa, &on, &params(true)).unwrap();
    let s_off = run_qc(&bam, &fa, &off, &params(false)).unwrap();

    // BED: exactly one `crispr` row, 0-based half-open, expected interval.
    let bed_on = bed_rows(&on.with_extension("bed"));
    let crispr: Vec<&Vec<String>> = bed_on.iter().filter(|r| r[3] == "crispr").collect();
    assert_eq!(crispr.len(), 1, "{bed_on:?}");
    let r = crispr[0];
    assert_eq!(r[0], "c2");
    assert_eq!(r[1], ARRAY_START.to_string());
    assert_eq!(r[2], (ARRAY_START + ARRAY_LEN).to_string());
    assert_eq!((r[4].as_str(), r[5].as_str(), r[6].as_str()), ("0", "none", "."));
    assert_eq!((r[7].as_str(), r[8].as_str()), (".", "."));
    // The published detector's bound refinement reports the first repeat row
    // with a leading alignment gap (36 ungapped bases) and the other three rows
    // with one extra leading base (37 ungapped), so all three internal spacers
    // come out 29 bp instead of 30 and the consensus carries that leading `G`.
    // The span is unchanged: 36 + 3 * 37 + 3 * 29 = 4 * 36 + 3 * 30 = 234.
    // `repeat=36` is the integer mean of 36.75; `spacer=29` is that
    // detector-specific refinement, not a guarantee of the true 30 bp spacer.
    assert_eq!(
        r[9],
        format!("copies=4;repeat=36;spacer=29;consensus=G{}", std::str::from_utf8(REPEAT).unwrap())
    );
    let bed_off = bed_rows(&off.with_extension("bed"));
    assert!(!bed_off.iter().any(|r| r[3] == "crispr"), "{bed_off:?}");
    // Everything else in the BED is identical.
    let others: Vec<&Vec<String>> = bed_on.iter().filter(|r| r[3] != "crispr").collect();
    assert_eq!(others, bed_off.iter().collect::<Vec<_>>());

    // VCF: the site inside the array carries CRISPR; TIER is unchanged, and
    // the records are otherwise identical.
    let vcf_on = vcf_data(&on.with_extension("vcf"));
    let vcf_off = vcf_data(&off.with_extension("vcf"));
    assert_eq!(vcf_on.len(), vcf_off.len(), "{vcf_on:?}");
    let pos_of = |l: &String| l.split('\t').nth(1).and_then(|p| p.parse::<usize>().ok());
    let inside = vcf_on
        .iter()
        .find(|l| pos_of(l) == Some(SNP_C2 + 1))
        .unwrap_or_else(|| panic!("no VCF record at POS {}: {vcf_on:?}", SNP_C2 + 1));
    let f: Vec<&str> = inside.split('\t').collect();
    assert!(f[7].contains(";CRISPR"), "{inside}");
    assert!(f[7].contains("TIER=RED"), "{inside}");
    let outside = vcf_on
        .iter()
        .find(|l| pos_of(l) == Some(SNP_C1 + 1))
        .unwrap_or_else(|| panic!("no VCF record at POS {}: {vcf_on:?}", SNP_C1 + 1));
    assert!(!outside.split('\t').nth(7).unwrap().contains("CRISPR"), "{outside}");
    for (a, b) in vcf_on.iter().zip(&vcf_off) {
        assert_eq!(a.replace(";CRISPR", ""), *b);
    }
    let vcf_text_on = std::fs::read_to_string(on.with_extension("vcf")).unwrap();
    let vcf_text_off = std::fs::read_to_string(off.with_extension("vcf")).unwrap();
    assert!(vcf_text_on.contains(
        "##INFO=<ID=CRISPR,Number=0,Type=Flag,Description=\"Inside a CRISPR array detected by pilercr-rs (annotation only; TIER unchanged)\">"
    ));
    assert!(!vcf_text_off.contains("ID=CRISPR"));

    // JSON: counts with the flag, zeros without.
    assert_eq!(s_on.crispr.arrays, 1);
    assert_eq!(s_on.crispr.bp, ARRAY_LEN as u64);
    assert_eq!((s_off.crispr.arrays, s_off.crispr.bp), (0, 0));
    let j: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(on.with_extension("json")).unwrap()).unwrap();
    assert_eq!(j["crispr"], serde_json::json!({"arrays": 1, "bp": ARRAY_LEN}));
    let j_off: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(off.with_extension("json")).unwrap()).unwrap();
    assert_eq!(j_off["crispr"], serde_json::json!({"arrays": 0, "bp": 0}));

    // Findings and alerts are identical with and without the flag.
    assert_eq!(
        std::fs::read_to_string(on.with_extension("findings.bed")).unwrap(),
        std::fs::read_to_string(off.with_extension("findings.bed")).unwrap()
    );
    assert_eq!(json_without_volatile(&on.with_extension("json")), json_without_volatile(&off.with_extension("json")));
}

/// Without `--crispr` there are no rows, no INFO flag and zero counts.
#[test]
fn crispr_off_is_inert() {
    let dir = workdir("off");
    let refs = reference();
    let fa = write_fasta(&dir, &refs);
    let mut reads = reads_for(0, &refs[0], SNP_C1);
    reads.extend(reads_for(1, &refs[1], SNP_C2));
    let bam = write_bam(&dir, &refs, reads);
    let prefix = dir.join("out");
    let s = run_qc(&bam, &fa, &prefix, &params(false)).unwrap();
    assert_eq!((s.crispr.arrays, s.crispr.bp), (0, 0));
    let bed = std::fs::read_to_string(prefix.with_extension("bed")).unwrap();
    assert!(!bed.contains("\tcrispr\t"), "{bed}");
    let vcf = std::fs::read_to_string(prefix.with_extension("vcf")).unwrap();
    assert!(!vcf.contains("CRISPR"), "{vcf}");
}

/// A real reference with no CRISPR arrays reports zeros.
#[test]
fn real_reference_without_arrays_has_zero_counts() {
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("testdata/data/cramhook-testdata-v1/ont");
    let (fa, cram) = (base.join("cons.fa"), base.join("cons.cram"));
    if !fa.exists() || !cram.exists() {
        eprintln!("skipping: test data not fetched ({})", base.display());
        return;
    }
    let dir = workdir("real");
    let prefix = dir.join("out");
    let p = QcParams { preset: Some(Preset::LrHq), threads: 2, crispr: true, ..Default::default() };
    let s = run_qc(&cram, &fa, &prefix, &p).unwrap();
    assert_eq!((s.crispr.arrays, s.crispr.bp), (0, 0));
    let bed = std::fs::read_to_string(prefix.with_extension("bed")).unwrap();
    assert!(!bed.contains("\tcrispr\t"), "{bed}");
    let j: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(prefix.with_extension("json")).unwrap()).unwrap();
    assert_eq!(j["crispr"], serde_json::json!({"arrays": 0, "bp": 0}));
}

/// The `CRISPR` flag follows the half-open array interval: `POS - 1` equals
/// `start - 1` (out), `start` (in), `end - 1` (in) and `end` (out again).
#[test]
fn crispr_vcf_flag_is_half_open() {
    let dir = workdir("bounds");
    let refs = reference();
    let fa = write_fasta(&dir, &refs);
    let (start, end) = (ARRAY_START, ARRAY_START + ARRAY_LEN);
    let reads = reads_at(1, &refs[1], &[start - 1, start, end - 1, end]);
    let bam = write_bam(&dir, &refs, reads);
    let prefix = dir.join("out");
    let s = run_qc(&bam, &fa, &prefix, &params(true)).unwrap();
    assert_eq!((s.crispr.arrays, s.crispr.bp), (1, ARRAY_LEN as u64));

    // The interval is the same one the annotation test expects.
    let bed = bed_rows(&prefix.with_extension("bed"));
    let rows: Vec<&Vec<String>> = bed.iter().filter(|r| r[3] == "crispr").collect();
    assert_eq!(rows.len(), 1, "{bed:?}");
    assert_eq!(rows[0][1], start.to_string());
    assert_eq!(rows[0][2], end.to_string());

    let vcf = vcf_data(&prefix.with_extension("vcf"));
    let flag = |pos0: usize| -> bool {
        let want = (pos0 + 1).to_string();
        let line = vcf
            .iter()
            .find(|l| l.split('\t').nth(1) == Some(want.as_str()))
            .unwrap_or_else(|| panic!("no VCF record at POS {}: {vcf:?}", pos0 + 1));
        line.split('\t').nth(7).unwrap().split(';').any(|f| f == "CRISPR")
    };
    assert!(!flag(start - 1), "0-based {} is left of the array", start - 1);
    assert!(flag(start), "0-based {start} is the first array base");
    assert!(flag(end - 1), "0-based {} is the last array base", end - 1);
    assert!(!flag(end), "0-based {end} is right of the array");
}

/// A 1 bp insertion before reference base `ARRAY_START + 1` starts with its
/// CIGAR anchor inside the array, but the caller left-normalises it to before
/// `ARRAY_START`: the VCF record follows the normalised anchor and is not
/// flagged.
#[test]
fn crispr_flag_follows_the_normalised_indel_anchor() {
    let dir = workdir("indel");
    let refs = reference();
    let fa = write_fasta(&dir, &refs);
    // The first array base is `G` (the repeat starts with it), so inserting a
    // `G` before the second base normalises left to before the first one.
    assert_eq!(refs[1][ARRAY_START], b'G');
    let reads: Vec<Rd> = (0..12)
        .map(|i| {
            let mut seq = refs[1][900..ARRAY_START + 1].to_vec();
            seq.push(b'G');
            seq.extend_from_slice(&refs[1][ARRAY_START + 1..1400]);
            Rd {
                tid: 1,
                pos: 900,
                seq,
                rev: i % 2 == 1,
                cigar: Some(CigarString(vec![
                    Cigar::Match((ARRAY_START + 1 - 900) as u32),
                    Cigar::Ins(1),
                    Cigar::Match((1400 - ARRAY_START - 1) as u32),
                ])),
            }
        })
        .collect();
    let bam = write_bam(&dir, &refs, reads);
    let prefix = dir.join("out");
    let s = run_qc(&bam, &fa, &prefix, &params(true)).unwrap();
    assert_eq!((s.crispr.arrays, s.crispr.bp), (1, ARRAY_LEN as u64));

    let vcf = vcf_data(&prefix.with_extension("vcf"));
    let ins: Vec<(usize, &String)> = vcf
        .iter()
        .filter(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[3].len() == 1 && f[4].len() == 2
        })
        .map(|l| (l.split('\t').nth(1).unwrap().parse::<usize>().unwrap(), l))
        .collect();
    assert_eq!(ins.len(), 1, "{vcf:?}");
    let (pos, line) = ins[0];
    // Unnormalised the insertion sits before 0-based `ARRAY_START + 1`, so its
    // VCF POS would be `ARRAY_START + 1` with POS - 1 = `ARRAY_START` (inside);
    // the normalised POS is `ARRAY_START` with POS - 1 outside.
    assert!(pos - 1 < ARRAY_START, "the anchor did not normalise out of the array (POS {pos}): {line}");
    assert!(!line.split('\t').nth(7).unwrap().contains("CRISPR"), "{line}");
    assert!(
        !vcf.iter().any(|l| l.split('\t').nth(1) == Some((ARRAY_START + 1).to_string().as_str())),
        "the unnormalised anchor must not be reported: {vcf:?}"
    );
}

/// CRLF line endings and a missing final newline must not stop detection.
#[test]
fn crlf_reference_without_a_final_newline_still_detects() {
    let dir = workdir("crlf");
    let refs = reference();
    let fa = write_fasta_styled(&dir, "crlf.fa", &refs, true, false);
    let mut reads = reads_for(0, &refs[0], SNP_C1);
    reads.extend(reads_for(1, &refs[1], SNP_C2));
    let bam = write_bam(&dir, &refs, reads);
    let prefix = dir.join("out");
    let s = run_qc(&bam, &fa, &prefix, &params(true)).unwrap();
    assert_eq!((s.crispr.arrays, s.crispr.bp), (1, ARRAY_LEN as u64));
    let bed = bed_rows(&prefix.with_extension("bed"));
    let rows: Vec<&Vec<String>> = bed.iter().filter(|r| r[3] == "crispr").collect();
    assert_eq!(rows.len(), 1, "{bed:?}");
    assert_eq!(rows[0][1], ARRAY_START.to_string());
    assert_eq!(rows[0][2], (ARRAY_START + ARRAY_LEN).to_string());
}

/// A whole-pass detection failure fails the run with context; only individual
/// unmappable arrays warn and are skipped.
#[test]
fn crispr_detection_failure_fails_the_run() {
    let dir = workdir("fail");
    let refs = reference();
    let fa = write_fasta(&dir, &refs);
    // A header with no newline and no sequence: qc's own FASTA reader accepts
    // it, but pilercr-rs rejects the unfinished header.
    let mut text = std::fs::read_to_string(&fa).unwrap();
    text.push_str(">truncated");
    std::fs::write(&fa, text).unwrap();
    let bam = write_bam(&dir, &refs, reads_for(0, &refs[0], SNP_C1));

    let prefix = dir.join("out");
    let o = Command::new(env!("CARGO_BIN_EXE_cramhook"))
        .args(["qc", "--in"])
        .arg(&bam)
        .args(["--ref"])
        .arg(&fa)
        .args(["--out-prefix"])
        .arg(&prefix)
        .args(["--preset", "map-hifi", "--threads", "2", "--crispr"])
        .output()
        .expect("running cramhook qc");
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success(), "qc --crispr must fail on a bad reference FASTA: {err}");
    assert!(err.contains("reading FASTA"), "{err}");
    assert!(err.contains("unfinished FASTA header"), "{err}");

    // The same reference without `--crispr` still runs: the failure is the pass.
    run_qc(&bam, &fa, &dir.join("ok/out"), &params(false)).unwrap();
}
