//! Integration tests for the per-contig circularity evidence in `cramhook qc`.
//!
//! Synthetic references and reads exercise: a clean circular junction (split
//! primary + supplementary records with valid `SA` tags), a duplicated junction
//! (the last 200 bp duplicate the first 200 bp -> ends overlap), a linear contig
//! whose reads stop at the ends, hard-clipped reverse pairs, supplementary-only
//! pairs, a read with an interior alignment (must not span), wrong-order reads,
//! and the exact `--circ-min-frac`, clipped-share and `--circ-end-tol` boundaries.

use std::path::{Path, PathBuf};

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{self, Format, Header, Record, Writer};

use cramhook::{Preset, QcParams, run_qc};

const L: usize = 30_000;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

fn workdir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/worker-circ/test").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn random_ref(len: usize, seed: u64) -> Vec<u8> {
    let mut rng = Lcg(seed);
    (0..len).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect()
}

#[derive(Clone)]
struct Rd {
    name: String,
    pos: usize,
    cigar: Vec<Cigar>,
    seq: Vec<u8>,
    rev: bool,
    supp: bool,
    sa: Option<String>,
}

fn strand(rev: bool) -> char {
    if rev {
        '-'
    } else {
        '+'
    }
}

fn write_bam(dir: &Path, refseq: &[u8], header_line: &str, mut reads: Vec<Rd>) -> (PathBuf, PathBuf) {
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">{header_line}\n{}\n", std::str::from_utf8(refseq).unwrap())).unwrap();
    let name = header_line.split_whitespace().next().unwrap();
    let off = 1 + header_line.len() + 1;
    std::fs::write(
        dir.join("ref.fa.fai"),
        format!("{name}\t{}\t{}\t{}\t{}\n", refseq.len(), off, refseq.len(), refseq.len() + 1),
    )
    .unwrap();
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", name).push_tag(b"LN", refseq.len()));
    let bam = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();
    reads.sort_by_key(|r| r.pos);
    for r in &reads {
        let mut rec = Record::new();
        let qual = vec![30u8; r.seq.len()];
        rec.set(r.name.as_bytes(), Some(&CigarString(r.cigar.clone())), &r.seq, &qual);
        rec.set_tid(0);
        rec.set_pos(r.pos as i64);
        rec.set_mapq(60);
        let mut flags = 0u16;
        if r.rev {
            flags |= 0x10;
        }
        if r.supp {
            flags |= 0x800;
        }
        rec.set_flags(flags);
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        if let Some(sa) = &r.sa {
            rec.push_aux(b"SA", Aux::String(sa)).unwrap();
        }
        w.write(&rec).unwrap();
    }
    drop(w);
    bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
    (bam, fa)
}

/// An `SA` entry: `rname,pos1,strand,CIGAR,mapQ,NM;`.
fn sa_entry(pos0: usize, rev: bool, cigar: &str) -> String {
    format!("chr,{},{},{cigar},60,0;", pos0 + 1, strand(rev))
}

/// A clean spanning read: end part query `[0,200)`, start part query `[200,400)`.
fn clean_pair(name: &str, refseq: &[u8], rev: bool) -> Vec<Rd> {
    let mut seq = Vec::with_capacity(400);
    seq.extend_from_slice(&refseq[L - 200..L]);
    seq.extend_from_slice(&refseq[0..200]);
    vec![
        Rd {
            name: name.to_string(),
            pos: L - 200,
            cigar: vec![Cigar::Match(200), Cigar::SoftClip(200)],
            seq: seq.clone(),
            rev,
            supp: false,
            sa: Some(sa_entry(0, rev, "200S200M")),
        },
        Rd {
            name: name.to_string(),
            pos: 0,
            cigar: vec![Cigar::SoftClip(200), Cigar::Match(200)],
            seq,
            rev,
            supp: true,
            sa: Some(sa_entry(L - 200, rev, "200M200S")),
        },
    ]
}

/// A spanning read whose parts overlap by 200 bp (the reference's last 200 bp
/// duplicate its first 200 bp). End part query `[0,400)`, start `[200,600)`.
fn overlap_pair(name: &str, refseq: &[u8], rev: bool) -> Vec<Rd> {
    // SEQ is reference-forward: [L-400, L) then the unique [200, 400) bases.
    let mut seq = Vec::with_capacity(600);
    seq.extend_from_slice(&refseq[L - 400..L]);
    seq.extend_from_slice(&refseq[200..400]);
    vec![
        Rd {
            name: name.to_string(),
            pos: L - 400,
            cigar: vec![Cigar::Match(400), Cigar::SoftClip(200)],
            seq: seq.clone(),
            rev,
            supp: false,
            sa: Some(sa_entry(0, rev, "200S400M")),
        },
        Rd {
            name: name.to_string(),
            pos: 0,
            cigar: vec![Cigar::SoftClip(200), Cigar::Match(400)],
            seq,
            rev,
            supp: true,
            sa: Some(sa_entry(L - 400, rev, "400M200S")),
        },
    ]
}

fn params(threads: usize, chunk: Option<u64>) -> QcParams {
    QcParams { preset: Some(Preset::MapHifi), threads, chunk_size: chunk, ..Default::default() }
}

fn json(path: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn json_without_volatile(path: &Path) -> serde_json::Value {
    let mut v = json(path);
    let o = v.as_object_mut().unwrap();
    o.remove("timing_s");
    o.remove("threads");
    o.remove("chunks");
    o.remove("inputs");
    if let Some(l) = o.get_mut("linkage").and_then(|l| l.as_object_mut()) {
        l.remove("timing_s");
    }
    v
}

/// A single-record read reaching the end (`200M` at `L-200`) or the start.
fn end_only(i: usize, refseq: &[u8], at_end: bool, clip: u64) -> Rd {
    let name = if at_end { format!("e{i}") } else { format!("s{i}") };
    let pos = if at_end { L - 200 } else { 0 };
    let seqend = if at_end { &refseq[L - 200..L] } else { &refseq[0..200] };
    let mut cigar = vec![Cigar::Match(200)];
    let mut seq = seqend.to_vec();
    if clip > 0 {
        cigar.push(Cigar::SoftClip(clip as u32));
        seq.extend(std::iter::repeat(b'A').take(clip as usize));
    }
    Rd { name, pos, cigar, seq, rev: false, supp: false, sa: None }
}

fn circ(j: &serde_json::Value) -> &serde_json::Value {
    &j["contigs"]["chr"]["circularity"]
}

#[test]
fn sa_fixtures_are_valid() {
    let refseq = random_ref(L, 7);
    for rev in [false, true] {
        for r in clean_pair("p", &refseq, rev).into_iter().chain(overlap_pair("q", &refseq, rev)) {
            let sa = r.sa.as_deref().unwrap();
            assert!(sa.ends_with(';'), "{sa}");
            for entry in sa.split(';').filter(|e| !e.is_empty()) {
                let f: Vec<&str> = entry.split(',').collect();
                assert_eq!(f.len(), 6, "SA entry is not comma-separated: {entry:?}");
                assert!(!f[0].is_empty());
                assert!(f[1].parse::<u64>().unwrap() > 0);
            }
        }
    }
    // The crate's strict parser accepts it too.
    let r = clean_pair("p", &refseq, false);
    let entries = cramhook::translate::parse_sa(r[0].sa.as_deref().unwrap()).expect("SA parses");
    assert_eq!(entries.len(), 1);
    assert_eq!((entries[0].rname.as_str(), entries[0].strand, entries[0].pos), ("chr", '+', 1));
}

#[test]
fn circular_clean_junction() {
    let dir = workdir("clean");
    let refseq = random_ref(L, 11);
    let mut reads = Vec::new();
    for i in 0..4 {
        reads.extend(clean_pair(&format!("s{i}"), &refseq, i % 2 == 1));
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let j = json(&p.with_extension("json"));
    let c = circ(&j);
    assert_eq!(c["verdict"], "circular", "{c}");
    assert_eq!(c["header_circular"], true, "{c}");
    assert_eq!(c["spanning_reads"], 4, "{c}");
    assert_eq!(c["junction"], "clean", "{c}");
    assert_eq!(c["junction_offset_median"], 0, "{c}");
    assert_eq!(c["junction_offset_consistent"], 1.0, "{c}");
    assert_eq!(j["circularity"]["circular"], 1, "{j}");
    assert!(!j["warnings"].as_array().unwrap().iter().any(|w| w["code"] == "circularity_mismatch"), "{j}");
}

#[test]
fn circular_overlap_junction() {
    let dir = workdir("overlap");
    let mut refseq = random_ref(L, 22);
    // The reference's last 200 bp duplicate its first 200 bp.
    let dup = refseq[0..200].to_vec();
    refseq[L - 200..L].copy_from_slice(&dup);
    let mut reads = Vec::new();
    for i in 0..4 {
        reads.extend(overlap_pair(&format!("s{i}"), &refseq, i % 2 == 1));
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let j = json(&p.with_extension("json"));
    let c = circ(&j);
    assert_eq!(c["verdict"], "circular", "{c}");
    assert_eq!(c["junction"], "overlap_200", "{c}");
    assert_eq!(c["junction_offset_median"], -200, "{c}");
    let w = j["warnings"].as_array().unwrap().iter().find(|w| w["code"] == "junction_overlap").unwrap_or_else(|| panic!("{j}"));
    assert_eq!(w["contigs"], serde_json::json!(["chr"]), "{w}");
    assert_eq!(w["bp"], 200, "{w}");
    assert!(w["message"].as_str().unwrap().contains("200 bp"), "{w}");
}

#[test]
fn linear_reads_stop_at_ends() {
    let dir = workdir("linear");
    let refseq = random_ref(L, 33);
    let mut reads = Vec::new();
    for i in 0..6 {
        reads.push(end_only(i, &refseq, true, 0));
        reads.push(end_only(i, &refseq, false, 0));
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=false", reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let j = json(&p.with_extension("json"));
    let c = circ(&j);
    assert_eq!(c["verdict"], "linear", "{c}");
    assert_eq!(c["end_reads"], 6, "{c}");
    assert_eq!(c["start_reads"], 6, "{c}");
    assert_eq!(c["spanning_reads"], 0, "{c}");
    assert_eq!(c["junction"], serde_json::Value::Null, "{c}");
    assert_eq!(j["circularity"]["linear"], 1, "{j}");
}

#[test]
fn hard_clipped_reverse_pair() {
    let dir = workdir("hardclip");
    let refseq = random_ref(L, 55);
    let mut reads = Vec::new();
    for i in 0..3 {
        let name = format!("h{i}");
        reads.push(Rd {
            name: name.clone(),
            pos: L - 200,
            cigar: vec![Cigar::Match(200), Cigar::HardClip(200)],
            seq: refseq[L - 200..L].to_vec(),
            rev: true,
            supp: false,
            sa: Some(sa_entry(0, true, "200H200M")),
        });
        reads.push(Rd {
            name: name.clone(),
            pos: 0,
            cigar: vec![Cigar::HardClip(200), Cigar::Match(200)],
            seq: refseq[0..200].to_vec(),
            rev: true,
            supp: true,
            sa: Some(sa_entry(L - 200, true, "200M200H")),
        });
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let jv = json(&p.with_extension("json"));
    let c = circ(&jv);
    assert_eq!(c["verdict"], "circular", "{c}");
    assert_eq!(c["junction"], "clean", "{c}");
    assert_eq!(c["junction_offset_median"], 0, "{c}");
}

#[test]
fn supplementary_only_pair() {
    let dir = workdir("supponly");
    let refseq = random_ref(L, 66);
    let mut reads = Vec::new();
    for i in 0..3 {
        let mut pair = clean_pair(&format!("p{i}"), &refseq, false);
        // Both records supplementary: no primary at all.
        for r in &mut pair {
            r.supp = true;
        }
        reads.extend(pair);
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let jv = json(&p.with_extension("json"));
    let c = circ(&jv);
    assert_eq!(c["verdict"], "circular", "{c}");
    assert_eq!(c["spanning_reads"], 3, "{c}");
    assert_eq!(c["junction"], "clean", "{c}");
}

/// A read with three alignments where an interior alignment lies between the end
/// and the start in query space must not count as spanning.
#[test]
fn interior_record_not_spanning() {
    let dir = workdir("interior");
    let refseq = random_ref(L, 77);
    let seq = refseq[0..300].to_vec();
    let mut reads = Vec::new();
    for i in 0..3 {
        let name = format!("m{i}");
        reads.push(Rd {
            name: name.clone(),
            pos: L - 100,
            cigar: vec![Cigar::Match(100), Cigar::SoftClip(200)],
            seq: seq.clone(),
            rev: false,
            supp: false,
            sa: Some(sa_entry(15_000, false, "100S100M100S") + &sa_entry(0, false, "200S100M")),
        });
        reads.push(Rd {
            name: name.clone(),
            pos: 15_000,
            cigar: vec![Cigar::SoftClip(100), Cigar::Match(100), Cigar::SoftClip(100)],
            seq: seq.clone(),
            rev: false,
            supp: true,
            sa: Some(sa_entry(L - 100, false, "100M200S") + &sa_entry(0, false, "200S100M")),
        });
        reads.push(Rd {
            name: name.clone(),
            pos: 0,
            cigar: vec![Cigar::SoftClip(200), Cigar::Match(100)],
            seq: seq.clone(),
            rev: false,
            supp: true,
            sa: Some(sa_entry(L - 100, false, "100M200S") + &sa_entry(15_000, false, "100S100M100S")),
        });
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let jv = json(&p.with_extension("json"));
    let c = circ(&jv);
    assert_eq!(c["spanning_reads"], 0, "{c}");
    assert_ne!(c["verdict"], "circular", "{c}");
}

/// A read visiting start -> end (forward) is not a spanning read.
#[test]
fn wrong_order_not_spanning() {
    let dir = workdir("wrongorder");
    let refseq = random_ref(L, 88);
    let mut seq = Vec::with_capacity(400);
    seq.extend_from_slice(&refseq[0..200]);
    seq.extend_from_slice(&refseq[L - 200..L]);
    let mut reads = Vec::new();
    for i in 0..3 {
        let name = format!("w{i}");
        reads.push(Rd {
            name: name.clone(),
            pos: L - 200,
            cigar: vec![Cigar::SoftClip(200), Cigar::Match(200)],
            seq: seq.clone(),
            rev: false,
            supp: false,
            sa: Some(sa_entry(0, false, "200M200S")),
        });
        reads.push(Rd {
            name: name.clone(),
            pos: 0,
            cigar: vec![Cigar::Match(200), Cigar::SoftClip(200)],
            seq: seq.clone(),
            rev: false,
            supp: true,
            sa: Some(sa_entry(L - 200, false, "200S200M")),
        });
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let jv = json(&p.with_extension("json"));
    let c = circ(&jv);
    assert_eq!(c["spanning_reads"], 0, "{c}");
    assert_ne!(c["verdict"], "circular", "{c}");
}

/// `--circ-min-frac` exact boundary: 3/10 = 0.3 is circular, 0.31 is not.
#[test]
fn min_frac_boundary() {
    let dir = workdir("minfrac");
    let refseq = random_ref(L, 99);
    let mut reads = Vec::new();
    for i in 0..3 {
        reads.extend(clean_pair(&format!("sp{i}"), &refseq, i % 2 == 1));
    }
    for i in 0..7 {
        reads.push(end_only(i, &refseq, true, 0));
        reads.push(end_only(i, &refseq, false, 0));
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let c = json(&p.with_extension("json"));
    assert_eq!(circ(&c)["spanning_reads"], 3, "{c}");
    assert_eq!(circ(&c)["end_reads"], 10, "{c}");
    assert_eq!(circ(&c)["verdict"], "circular", "frac exactly 0.3 is circular");

    // Raising the fraction above 0.3 makes the same data unclear.
    let p2 = dir.join("out2");
    let strict = QcParams { circ_min_frac: 0.31, ..params(1, None) };
    run_qc(&bam, &fa, &p2, &strict).unwrap();
    assert_eq!(circ(&json(&p2.with_extension("json")))["verdict"], "unclear");
}

/// The linear rule's clipped-share boundary: 0/10 < 0.1 is linear, 1/10 is not.
#[test]
fn clipped_share_boundary() {
    let dir = workdir("clipshare");
    let refseq = random_ref(L, 111);
    let mut reads = Vec::new();
    for i in 0..10 {
        reads.push(end_only(i, &refseq, true, 0));
        reads.push(end_only(i, &refseq, false, 0));
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=false", reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let c = json(&p.with_extension("json"));
    assert_eq!(circ(&c)["verdict"], "linear", "{c}");
    assert_eq!(circ(&c)["clipped_at_end"], 0, "{c}");

    // One of the ten end reads is clipped past the end: 0.1 is not < 0.1.
    let dir2 = workdir("clipshare1");
    let mut reads2 = Vec::new();
    for i in 0..10 {
        reads2.push(end_only(i, &refseq, true, if i == 0 { 80 } else { 0 }));
        reads2.push(end_only(i, &refseq, false, 0));
    }
    let (bam2, fa2) = write_bam(&dir2, &refseq, "chr circular=false", reads2);
    let p2 = dir2.join("out");
    run_qc(&bam2, &fa2, &p2, &params(1, None)).unwrap();
    let c2 = json(&p2.with_extension("json"));
    assert_eq!(circ(&c2)["clipped_at_end"], 1, "{c2}");
    assert_eq!(circ(&c2)["end_reads"], 10, "{c2}");
    assert_eq!(circ(&c2)["verdict"], "unclear", "clipped share 0.1 is not < 0.1");
}

/// `--circ-end-tol` boundary: an alignment ending exactly `tol` short of the end
/// reaches it, one base further does not (likewise at the start).
#[test]
fn end_tol_boundary() {
    let dir = workdir("endtol");
    let refseq = random_ref(L, 222);
    let reads = vec![
        // ref_end = L-50 -> reaches with tol 50.
        Rd { name: "hit".into(), pos: L - 250, cigar: vec![Cigar::Match(200)], seq: refseq[L - 250..L - 50].to_vec(), rev: false, supp: false, sa: None },
        // ref_end = L-51 -> does not.
        Rd { name: "miss".into(), pos: L - 251, cigar: vec![Cigar::Match(200)], seq: refseq[L - 251..L - 51].to_vec(), rev: false, supp: false, sa: None },
        // pos 50 -> start reached; pos 51 -> not.
        Rd { name: "shit".into(), pos: 50, cigar: vec![Cigar::Match(200)], seq: refseq[50..250].to_vec(), rev: false, supp: false, sa: None },
        Rd { name: "smiss".into(), pos: 51, cigar: vec![Cigar::Match(200)], seq: refseq[51..251].to_vec(), rev: false, supp: false, sa: None },
    ];
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=false", reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let c = json(&p.with_extension("json"));
    assert_eq!(circ(&c)["end_reads"], 1, "{c}");
    assert_eq!(circ(&c)["start_reads"], 1, "{c}");
}

/// The wrap decision and its source in the JSON: `--linear` forces linear; a
/// header `circular=true` forces circular; otherwise the verdict decides
/// (circular/unclear wrap, linear does not).
#[test]
fn wrap_status_and_source() {
    let refseq = random_ref(L, 101);
    let mut circ_reads = Vec::new();
    for i in 0..4 {
        circ_reads.extend(clean_pair(&format!("s{i}"), &refseq, i % 2 == 1));
    }

    // Header circular=true -> circular, by header.
    let dir = workdir("wrap_header");
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", circ_reads.clone());
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let j = json(&p.with_extension("json"));
    assert_eq!(j["contigs"]["chr"]["circularity"]["verdict"], "circular", "{j}");
    assert_eq!(j["contigs"]["chr"]["circular"], true, "{j}");
    assert_eq!(j["contigs"]["chr"]["circular_by"], "header", "{j}");

    // Same reads, no header flag -> verdict circular, by verdict.
    let dir = workdir("wrap_verdict_circ");
    let (bam, fa) = write_bam(&dir, &refseq, "chr", circ_reads.clone());
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let j = json(&p.with_extension("json"));
    assert_eq!(j["contigs"]["chr"]["circularity"]["verdict"], "circular", "{j}");
    assert_eq!(j["contigs"]["chr"]["circular"], true, "{j}");
    assert_eq!(j["contigs"]["chr"]["circular_by"], "verdict", "{j}");

    // `--linear` overrides the circular verdict/header -> by flag.
    let dir = workdir("wrap_flag");
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", circ_reads);
    let p = dir.join("out");
    let qp = QcParams { linear: vec!["chr".to_string()], ..params(1, None) };
    run_qc(&bam, &fa, &p, &qp).unwrap();
    let j = json(&p.with_extension("json"));
    assert_eq!(j["contigs"]["chr"]["circularity"]["verdict"], "circular", "{j}");
    assert_eq!(j["contigs"]["chr"]["circular"], false, "{j}");
    assert_eq!(j["contigs"]["chr"]["circular_by"], "flag", "{j}");

    // A linear verdict does not wrap (by verdict).
    let dir = workdir("wrap_linear");
    let refseq2 = random_ref(L, 202);
    let mut reads2 = Vec::new();
    for i in 0..6 {
        reads2.push(end_only(i, &refseq2, true, 0));
        reads2.push(end_only(i, &refseq2, false, 0));
    }
    let (bam, fa) = write_bam(&dir, &refseq2, "chr", reads2);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let j = json(&p.with_extension("json"));
    assert_eq!(j["contigs"]["chr"]["circularity"]["verdict"], "linear", "{j}");
    assert_eq!(j["contigs"]["chr"]["circular"], false, "{j}");
    assert_eq!(j["contigs"]["chr"]["circular_by"], "verdict", "{j}");

    // An unclear verdict wraps (by verdict).
    let dir = workdir("wrap_unclear");
    let refseq3 = random_ref(L, 303);
    let mut reads3 = Vec::new();
    for i in 0..10 {
        reads3.push(end_only(i, &refseq3, true, 80));
        reads3.push(end_only(i, &refseq3, false, 80));
    }
    let (bam, fa) = write_bam(&dir, &refseq3, "chr", reads3);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let j = json(&p.with_extension("json"));
    assert_eq!(j["contigs"]["chr"]["circularity"]["verdict"], "unclear", "{j}");
    assert_eq!(j["contigs"]["chr"]["circular"], true, "{j}");
    assert_eq!(j["contigs"]["chr"]["circular_by"], "verdict", "{j}");
}

#[test]
fn threads_identical() {
    let dir = workdir("threads");
    let refseq = random_ref(L, 44);
    let mut reads = Vec::new();
    for i in 0..4 {
        reads.extend(clean_pair(&format!("s{i}"), &refseq, i % 2 == 1));
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", reads);
    let p1 = dir.join("t1/out");
    let p8 = dir.join("t8/out");
    run_qc(&bam, &fa, &p1, &params(1, None)).unwrap();
    run_qc(&bam, &fa, &p8, &params(8, Some(1_500))).unwrap();
    assert_eq!(json_without_volatile(&p1.with_extension("json")), json_without_volatile(&p8.with_extension("json")));
}

/// A read joining two distant loci. `a`/`b` are the two breakpoints; the read's
/// first segment (read order) ends at one and the next starts at the other. The
/// forward and reverse layouts are the two strands of the same junction, so both
/// canonicalise to the same `(lo, hi)` pair.
fn junction_pair(name: &str, refseq: &[u8], a: usize, b: usize, rev: bool) -> Vec<Rd> {
    if rev {
        // Read 5' part maps reverse at [a, a+500) (aligned query [0,500)); 3'
        // part at [b-500, b) (aligned query [500,1000)). SEQ is
        // reference-oriented, so the aligned parts match the reference.
        let mut seq = vec![b'A'; 500];
        seq.extend_from_slice(&refseq[b - 500..b]);
        seq.extend_from_slice(&refseq[a..a + 500]);
        vec![
            Rd {
                name: name.to_string(),
                pos: a,
                cigar: vec![Cigar::SoftClip(1_000), Cigar::Match(500)],
                seq: seq.clone(),
                rev: true,
                supp: false,
                sa: None,
            },
            Rd {
                name: name.to_string(),
                pos: b - 500,
                cigar: vec![Cigar::SoftClip(500), Cigar::Match(500), Cigar::SoftClip(500)],
                seq,
                rev: true,
                supp: true,
                sa: None,
            },
        ]
    } else {
        let mut seq = refseq[b - 500..b].to_vec();
        seq.extend_from_slice(&refseq[a..a + 500]);
        seq.extend(std::iter::repeat(b'A').take(500));
        vec![
            Rd {
                name: name.to_string(),
                pos: b - 500,
                cigar: vec![Cigar::Match(500), Cigar::SoftClip(1_000)],
                seq: seq.clone(),
                rev: false,
                supp: false,
                sa: None,
            },
            Rd {
                name: name.to_string(),
                pos: a,
                cigar: vec![Cigar::SoftClip(500), Cigar::Match(500), Cigar::SoftClip(500)],
                seq,
                rev: false,
                supp: true,
                sa: None,
            },
        ]
    }
}

/// Split-read junction signal. A uniform ~100x reference carries a junction at
/// ~5 % of depth (3 forward + 2 reverse reads -> one cluster; the canonical
/// traversal is the reverse-strand one, so `fwd=2;rev=3`; flagged
/// `split_junction`, yellow), one below `--junction-min-reads` (not flagged), a
/// junction within `--junction-edge-excl` of a contig end and the end->start wrap
/// on the circular contig (both edge-excluded). The junction regions enter the
/// findings as clip-style `misjoin` evidence, and `--threads 1` and 8 give
/// identical output.
#[test]
fn split_junction_signal() {
    let dir = workdir("split_junction");
    let refseq = random_ref(L, 31);
    let (a1, b1) = (3_000usize, 12_000usize);
    let (a2, b2) = (6_000usize, 15_000usize);
    let (a3, b3) = (500usize, 10_000usize); // within --junction-edge-excl (2000)
    let mut reads: Vec<Rd> = Vec::new();
    // Uniform ~100x background: 2 kb reads every 20 bp (reference-forward).
    for i in 0..((L - 2_000) / 20 + 1) {
        let pos = i * 20;
        reads.push(Rd {
            name: format!("bg{i}"),
            pos,
            cigar: vec![Cigar::Match(2_000)],
            seq: refseq[pos..pos + 2_000].to_vec(),
            rev: false,
            supp: false,
            sa: None,
        });
    }
    // Locus 1: 3 forward + 2 reverse junction reads -> one merged cluster at ~5 %.
    for i in 0..3 {
        reads.extend(junction_pair(&format!("j1f_{i}"), &refseq, a1, b1, false));
    }
    for i in 0..2 {
        reads.extend(junction_pair(&format!("j1r_{i}"), &refseq, a1, b1, true));
    }
    // Locus 2: 4 junction reads -> below `--junction-min-reads` (5).
    for i in 0..4 {
        reads.extend(junction_pair(&format!("j2_{i}"), &refseq, a2, b2, false));
    }
    // Near a contig end: 8 reads, but edge-excluded on the circular contig.
    for i in 0..8 {
        reads.extend(junction_pair(&format!("j3_{i}"), &refseq, a3, b3, false));
    }
    // Origin crossing: first segment ends at the contig end, next starts at 0.
    for i in 0..10 {
        reads.extend(junction_pair(&format!("o{i}"), &refseq, 0, L, false));
    }

    let junctions = |bed: &str| -> Vec<Vec<String>> {
        bed.lines()
            .filter(|l| l.contains("\tsplit_junction\t"))
            .map(|l| l.split('\t').map(|s| s.to_string()).collect())
            .collect()
    };

    let (bam, fa) = write_bam(&dir, &refseq, "chr", reads.clone());
    let p1 = dir.join("out");
    let qp = params(1, None);
    let tol = qp.junction_tol as usize;
    run_qc(&bam, &fa, &p1, &qp).unwrap();
    let bed = std::fs::read_to_string(p1.with_extension("bed")).unwrap();
    let js = junctions(&bed);
    // Two regions (one per breakpoint) for locus 1 only.
    assert_eq!(js.len(), 2, "expected only the locus-1 junction:\n{bed}");
    let left = js.iter().find(|r| r[1] == (a1 - tol).to_string()).unwrap_or_else(|| panic!("{bed}"));
    let right = js.iter().find(|r| r[1] == (b1 - tol).to_string()).unwrap_or_else(|| panic!("{bed}"));
    assert_eq!(
        (left[2].as_str(), left[5].as_str(), left[6].as_str()),
        ((a1 + tol + 1).to_string().as_str(), "orange", "junction"),
        "{left:?}"
    );
    assert_eq!(
        (right[2].as_str(), right[5].as_str(), right[6].as_str()),
        ((b1 + tol + 1).to_string().as_str(), "orange", "junction"),
        "{right:?}"
    );
    assert!(
        left[9].contains(&format!("partner=chr:{b1};reads=5;")) && left[9].contains("span=9000")
            && left[9].contains("orient=tandem;fwd=2;rev=3"),
        "{left:?}"
    );
    assert!(right[9].contains(&format!("partner=chr:{a1};reads=5;")), "{right:?}");
    // No region at locus 2 (support 4 < 5), the near-end locus (edge-excluded)
    // or the origin wrap.
    let starts: Vec<&str> = js.iter().map(|r| r[1].as_str()).collect();
    let banned = [
        (a2 - tol).to_string(),
        (b2 - tol).to_string(),
        (a3 - tol).to_string(),
        (b3 - tol).to_string(),
        "0".to_string(),
    ];
    assert!(!starts.iter().any(|s| banned.iter().any(|b| s == b)), "below-k / edge / origin flagged:\n{bed}");

    // The two junction regions are linked (`--finding-link-max` default 100 kb)
    // into one finding spanning both breakpoints, as clip-style `misjoin`
    // evidence.
    let findings = std::fs::read_to_string(p1.with_extension("findings.bed")).unwrap();
    let jf: Vec<&str> = findings.lines().filter(|l| l.contains("members=split_junction:2")).collect();
    assert_eq!(jf.len(), 1, "{findings}");
    assert!(jf[0].contains(&format!("\t{}\t{}\t", a1 - tol, b1 + tol + 1)), "{jf:?}");
    assert!(jf[0].contains("\tmisjoin\t") && jf[0].contains("evidence=clips"), "{jf:?}");

    // JSON: two orange `split_junction` alerts (reported calls are the high bar).
    let j = json(&p1.with_extension("json"));
    assert_eq!(j["alerts"]["by_type"]["split_junction"]["orange"], 2, "{j}");

    // Threads 1 vs 8 (small chunks) agree.
    let dir8 = workdir("split_junction_t8");
    let (bam8, fa8) = write_bam(&dir8, &refseq, "chr", reads);
    let p8 = dir8.join("out");
    run_qc(&bam8, &fa8, &p8, &params(8, Some(1_000))).unwrap();
    let bed8 = std::fs::read_to_string(p8.with_extension("bed")).unwrap();
    assert_eq!(bed, bed8, "BED differs between --threads 1 and 8");
    assert_eq!(json_without_volatile(&p1.with_extension("json")), json_without_volatile(&p8.with_extension("json")));
}

/// The automatic wrap selection (verdict/header/flag) changes a real signal: a
/// split junction within `--junction-edge-excl` of a contig end is dropped when
/// the contig wraps but kept when it is linear; a FASTA header `circular=true`
/// overrides a linear verdict (header wins).
#[test]
fn wrap_selection_changes_split_junction() {
    let refseq = random_ref(L, 404);
    let (a, b) = (500usize, 10_000usize);
    let tol = params(1, None).junction_tol as usize;
    let mut linear_reads: Vec<Rd> = Vec::new();
    for i in 0..8 {
        linear_reads.extend(junction_pair(&format!("j{i}"), &refseq, a, b, false));
    }
    // End-only reads reach both ends without spanning -> verdict linear.
    for i in 0..6 {
        linear_reads.push(end_only(i, &refseq, true, 0));
        linear_reads.push(end_only(i, &refseq, false, 0));
    }
    let junctions = |bed: &str| bed.lines().filter(|l| l.contains("\tsplit_junction\t")).count();

    // (a) verdict linear, header absent -> kept (no edge exclusion).
    let dir = workdir("wrap_junction_linear");
    let (bam, fa) = write_bam(&dir, &refseq, "chr", linear_reads.clone());
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let j = json(&p.with_extension("json"));
    assert_eq!(j["contigs"]["chr"]["circularity"]["verdict"], "linear", "{j}");
    assert_eq!(j["contigs"]["chr"]["circular"], false, "{j}");
    assert_eq!(j["contigs"]["chr"]["circular_by"], "verdict", "{j}");
    let bed = std::fs::read_to_string(p.with_extension("bed")).unwrap();
    assert_eq!(junctions(&bed), 2, "linear: the near-end junction must be kept:\n{bed}");
    assert!(bed.lines().any(|l| l.contains("\tsplit_junction\t") && l.contains(&format!("\t{}\t", a - tol))), "{bed}");

    // (b) same reads + spanning pairs -> verdict circular -> edge-excluded.
    let dir = workdir("wrap_junction_circular");
    let mut circ_reads = linear_reads.clone();
    for i in 0..4 {
        circ_reads.extend(clean_pair(&format!("s{i}"), &refseq, i % 2 == 1));
    }
    let (bam, fa) = write_bam(&dir, &refseq, "chr", circ_reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let j = json(&p.with_extension("json"));
    assert_eq!(j["contigs"]["chr"]["circularity"]["verdict"], "circular", "{j}");
    assert_eq!(j["contigs"]["chr"]["circular"], true, "{j}");
    assert_eq!(j["contigs"]["chr"]["circular_by"], "verdict", "{j}");
    let bed = std::fs::read_to_string(p.with_extension("bed")).unwrap();
    assert_eq!(junctions(&bed), 0, "circular: the near-end junction must be edge-excluded:\n{bed}");

    // (c) linear verdict but header circular=true -> header wins -> excluded.
    let dir = workdir("wrap_junction_header");
    let (bam, fa) = write_bam(&dir, &refseq, "chr circular=true", linear_reads);
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();
    let j = json(&p.with_extension("json"));
    assert_eq!(j["contigs"]["chr"]["circularity"]["verdict"], "linear", "{j}");
    assert_eq!(j["contigs"]["chr"]["circular"], true, "{j}");
    assert_eq!(j["contigs"]["chr"]["circular_by"], "header", "{j}");
    let bed = std::fs::read_to_string(p.with_extension("bed")).unwrap();
    assert_eq!(junctions(&bed), 0, "header circular=true must win and edge-exclude:\n{bed}");
}

/// Cross-contig clip links, on the existing two-contig `clip_join` fixture: a
/// right-clip hotspot on `A` (`A:10999-11001`) joins a left-clip hotspot on `B`
/// (`B:7910-7956`) at 100 % of each side, so both get `clip_to` and both
/// findings get a mirrored `link=` note (each also carrying the partner-end
/// depth); the JSON counts them.
#[test]
fn cross_contig_clip_links() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test-qc-links");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let fx = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/clip_join");
    let (bam, fa) = (fx.join("oracle.cram"), fx.join("new.fa"));
    assert!(bam.exists() && fa.exists(), "missing clip_join fixture");
    let p = dir.join("out");
    run_qc(&bam, &fa, &p, &params(1, None)).unwrap();

    let bed = std::fs::read_to_string(p.with_extension("bed")).unwrap();
    assert!(
        bed.lines().any(|l| l.starts_with("A\t10999\t11001") && l.contains("clip_to=B:7910-7956;clip_to_reads=10")),
        "{bed}"
    );
    assert!(
        bed.lines().any(|l| l.starts_with("B\t7910\t7956") && l.contains("clip_to=A:11000-11002;clip_to_reads=10")),
        "{bed}"
    );

    let findings = std::fs::read_to_string(p.with_extension("findings.bed")).unwrap();
    // Both findings get the mirrored note (each wrapped finding is written twice),
    // with the f of the end they describe and the partner-end depth as a `link=`
    // field, not a BEDPE column.
    assert!(findings.contains("link=finding_2@B:7910-7956,reads=10,f=1.0000,partner_depth=1"), "{findings}");
    assert!(findings.contains("link=finding_1@A:10999-11001,reads=10,f=1.0000,partner_depth=10"), "{findings}");
    assert!(findings.contains("link=finding_2@B:7910-7956,reads=10,f=1.7857,partner_depth=5.6"), "{findings}");
    assert!(findings.contains("link=finding_1@A:11000-11002,reads=10,f=1.7857,partner_depth=2"), "{findings}");
    assert!(!findings.contains("f_partner"), "{findings}");
    assert!(!p.with_extension("links.bedpe").exists(), "links.bedpe must no longer be written");
    let j = json(&p.with_extension("json"));
    assert_eq!(j["xlinks"], 2, "{j}");
}

/// A two-contig read for [`write_bam2`].
struct Rd2 {
    name: String,
    tid: usize,
    pos: usize,
    cigar: Vec<Cigar>,
    seq: Vec<u8>,
    rev: bool,
    supp: bool,
}

/// Write a BAM/FASTA with `refs` contigs (in order) and reads tagged by contig.
fn write_bam2(dir: &Path, refs: &[(&str, Vec<u8>)], mut reads: Vec<Rd2>) -> (PathBuf, PathBuf) {
    let fa = dir.join("ref.fa");
    let mut fa_text = String::new();
    let mut fai = String::new();
    for (name, seq) in refs {
        let off = fa_text.len() + 1 + name.len() + 1;
        fa_text.push_str(&format!(">{name}\n{}\n", std::str::from_utf8(seq).unwrap()));
        fai.push_str(&format!("{name}\t{}\t{}\t{}\t{}\n", seq.len(), off, seq.len(), seq.len() + 1));
    }
    std::fs::write(&fa, fa_text).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), fai).unwrap();
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    for (name, seq) in refs {
        header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", *name).push_tag(b"LN", seq.len()));
    }
    let bam = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();
    reads.sort_by_key(|r| (r.tid, r.pos));
    for r in &reads {
        let mut rec = Record::new();
        let qual = vec![30u8; r.seq.len()];
        rec.set(r.name.as_bytes(), Some(&CigarString(r.cigar.clone())), &r.seq, &qual);
        rec.set_tid(r.tid as i32);
        rec.set_pos(r.pos as i64);
        rec.set_mapq(60);
        let mut flags = 0u16;
        if r.rev {
            flags |= 0x10;
        }
        if r.supp {
            flags |= 0x800;
        }
        rec.set_flags(flags);
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        w.write(&rec).unwrap();
    }
    drop(w);
    bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
    (bam, fa)
}

/// Two linear contigs whose hotspots link to each other. The mirrored `link=`
/// notes are identical on both findings and must be written once each (reviewer
/// issue 3).
#[test]
fn cross_contig_link_notes_are_deduplicated() {
    let dir = workdir("xlink_dupe");
    let a = random_ref(10_000, 71);
    let b = random_ref(10_000, 72);
    let mut reads = Vec::new();
    for i in 0..10 {
        let name = format!("r{i}");
        reads.push(Rd2 {
            name: name.clone(),
            tid: 0,
            pos: 5_000,
            cigar: vec![Cigar::SoftClip(1_000), Cigar::Match(1_000)],
            seq: vec![b'A'; 2_000],
            rev: true,
            supp: false,
        });
        reads.push(Rd2 {
            name,
            tid: 1,
            pos: 5_000,
            cigar: vec![Cigar::SoftClip(1_000), Cigar::Match(1_000)],
            seq: vec![b'A'; 2_000],
            rev: false,
            supp: true,
        });
    }
    let (bam, fa) = write_bam2(&dir, &[("A", a), ("B", b)], reads);
    let p = dir.join("out");
    let qp = QcParams { linear: vec!["A".to_string(), "B".to_string()], ..params(1, None) };
    run_qc(&bam, &fa, &p, &qp).unwrap();

    let findings = std::fs::read_to_string(p.with_extension("findings.bed")).unwrap();
    for line in findings.lines().filter(|l| !l.starts_with('#')) {
        let extra = line.split('\t').nth(7).unwrap_or("");
        let notes: Vec<&str> = extra
            .split(';')
            .flat_map(|part| part.split('|'))
            .filter_map(|n| n.strip_prefix("link="))
            .collect();
        let mut seen = std::collections::HashSet::new();
        for n in &notes {
            assert!(seen.insert(*n), "duplicate link note {n:?} in:\n{findings}");
        }
    }
    // The two findings do carry the mirrored note.
    assert!(findings.contains("link=finding_2@B:5000-5001"), "{findings}");
    assert!(findings.contains("link=finding_1@A:5000-5001"), "{findings}");
}
