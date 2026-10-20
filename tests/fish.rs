//! Integration tests for `cramhook fish`: the per-read extraction, its
//! categories, scopes, flags, FASTQ orientation and thread identity.

use std::io::Read;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Mutex, MutexGuard};

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{self, Format, Header, Record, Writer};

use cramhook::{BaitRecords, FishParams, FishScope, FishTier, run_fish};

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

fn random_seq(len: usize, seed: u64) -> Vec<u8> {
    let mut rng = Lcg(seed);
    (0..len).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect()
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        x => x,
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| comp(b)).collect()
}

fn workdir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/fish/test").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// One low-level BAM record (flags explicit so unmapped/supplementary can be built).
struct RawRec {
    name: String,
    pos: i64,
    cigar: Vec<Cigar>,
    seq: Vec<u8>,
    flags: u16,
}

fn params(threads: usize) -> FishParams {
    FishParams { threads, ..Default::default() }
}

/// Write ref + BAM, run fish with `p`, return the output prefix.
fn run(name: &str, refseq: &[u8], recs: &[RawRec], p: &FishParams) -> PathBuf {
    let dir = workdir(name);
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr\n{}\n", std::str::from_utf8(refseq).unwrap())).unwrap();
    std::fs::write(
        dir.join("ref.fa.fai"),
        format!("chr\t{}\t5\t{}\t{}\n", refseq.len(), refseq.len(), refseq.len() + 1),
    )
    .unwrap();
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", refseq.len()));
    let bam = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();
    let mut sorted: Vec<&RawRec> = recs.iter().collect();
    sorted.sort_by_key(|r| (r.pos < 0, r.pos));
    for r in sorted {
        let qual = vec![30u8; r.seq.len()];
        let mut rec = Record::new();
        rec.set(r.name.as_bytes(), Some(&CigarString(r.cigar.clone())), &r.seq, &qual);
        rec.set_tid(if r.pos < 0 { -1 } else { 0 });
        rec.set_pos(r.pos);
        rec.set_mapq(60);
        rec.set_flags(r.flags);
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        w.write(&rec).unwrap();
    }
    drop(w);
    bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
    let prefix = dir.join("out");
    run_fish(&bam, &fa, &prefix, p).unwrap();
    prefix
}

/// `P.fish.tsv` rows (no header) as field vectors.
fn tsv(prefix: &Path) -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(prefix.with_extension("fish.tsv")).unwrap();
    text.lines().skip(1).map(|l| l.split('\t').map(|s| s.to_string()).collect()).collect()
}

fn reasons_of(rows: &[Vec<String>], reason: &str) -> Vec<String> {
    let mut v: Vec<String> =
        rows.iter().filter(|r| r[2] == reason).map(|r| r[0].clone()).collect();
    v.sort();
    v
}

/// Decompressed FASTQ bytes.
fn fastq(prefix: &Path) -> Vec<u8> {
    let mut r = rust_htslib::bgzf::Reader::from_path(prefix.with_extension("fastq.gz")).unwrap();
    let mut v = Vec::new();
    r.read_to_end(&mut v).unwrap();
    v
}

fn json(prefix: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(prefix.with_extension("fish.json")).unwrap()).unwrap()
}

const L: usize = 9_000;
const READ_LEN: usize = 1_000;
const STEP: usize = 25;
const JUNCTION: usize = 2_000;
const INS_AT: usize = 5_000;
const INS_LEN: usize = 60;
const DEL_AT: usize = 7_000;
const DEL_LEN: usize = 80;

/// Background ~40x reads + a clip junction, a large insertion, a large deletion
/// and one unmapped read. Returns `(reference, reads, clip names, ins names,
/// del names)`.
fn build() -> (Vec<u8>, Vec<RawRec>, Vec<String>, Vec<String>, Vec<String>) {
    let r = random_seq(L, 7);
    let mut recs: Vec<RawRec> = Vec::new();
    let mut s = 0usize;
    while s + READ_LEN <= L {
        recs.push(RawRec {
            name: format!("bg{s}"),
            pos: s as i64,
            cigar: vec![Cigar::Match(READ_LEN as u32)],
            seq: r[s..s + READ_LEN].to_vec(),
            flags: if (s / STEP) % 2 == 1 { 0x10 } else { 0 },
        });
        s += STEP;
    }

    let mut clip = Vec::new();
    for k in 0..20 {
        let mut seq = r[JUNCTION - 500..JUNCTION].to_vec();
        seq.extend(std::iter::repeat_n(b'A', 120));
        recs.push(RawRec {
            name: format!("clipR{k}"),
            pos: (JUNCTION - 500) as i64,
            cigar: vec![Cigar::Match(500), Cigar::SoftClip(120)],
            seq,
            flags: if k % 2 == 0 { 0x10 } else { 0 },
        });
        clip.push(format!("clipR{k}"));
    }
    for k in 0..20 {
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(&r[JUNCTION..JUNCTION + 500]);
        recs.push(RawRec {
            name: format!("clipL{k}"),
            pos: JUNCTION as i64,
            cigar: vec![Cigar::SoftClip(120), Cigar::Match(500)],
            seq,
            flags: if k % 2 == 1 { 0x10 } else { 0 },
        });
        clip.push(format!("clipL{k}"));
    }

    let ins = random_seq(INS_LEN, 99);
    let mut ins_reads = Vec::new();
    for rec in recs.iter_mut() {
        if ins_reads.len() >= 15 {
            break;
        }
        if rec.name.starts_with("bg") && rec.pos >= 0 && (rec.pos as usize) < INS_AT && rec.pos as usize + READ_LEN > INS_AT {
            let p = rec.pos as usize;
            let off = INS_AT - p;
            rec.cigar = vec![
                Cigar::Match(off as u32),
                Cigar::Ins(INS_LEN as u32),
                Cigar::Match((READ_LEN - off) as u32),
            ];
            let mut seq = r[p..INS_AT].to_vec();
            seq.extend_from_slice(&ins);
            seq.extend_from_slice(&r[INS_AT..p + READ_LEN]);
            rec.seq = seq;
            ins_reads.push(rec.name.clone());
        }
    }
    assert_eq!(ins_reads.len(), 15);

    let mut del_reads = Vec::new();
    for rec in recs.iter_mut() {
        if del_reads.len() >= 15 {
            break;
        }
        if rec.name.starts_with("bg")
            && rec.pos >= 0
            && (rec.pos as usize) < DEL_AT
            && rec.pos as usize + READ_LEN >= DEL_AT + DEL_LEN
        {
            let p = rec.pos as usize;
            let off = DEL_AT - p;
            rec.cigar = vec![
                Cigar::Match(off as u32),
                Cigar::Del(DEL_LEN as u32),
                Cigar::Match((READ_LEN - off - DEL_LEN) as u32),
            ];
            let mut seq = r[p..DEL_AT].to_vec();
            seq.extend_from_slice(&r[DEL_AT + DEL_LEN..p + READ_LEN]);
            rec.seq = seq;
            del_reads.push(rec.name.clone());
        }
    }
    assert_eq!(del_reads.len(), 15);

    recs.push(RawRec {
        name: "unplaced1".to_string(),
        pos: -1,
        cigar: Vec::new(),
        seq: random_seq(500, 5),
        flags: 0x4,
    });

    (r, recs, clip, ins_reads, del_reads)
}

#[test]
fn fish_selects_each_category_and_thread_identity() {
    let (r, recs, clip, ins_reads, del_reads) = build();
    let p1 = run("cats_t1", &r, &recs, &params(1));
    let p8 = run("cats_t8", &r, &recs, &params(8));

    let rows = tsv(&p1);
    let mut want_clip = clip.clone();
    want_clip.sort();
    let mut want_ins = ins_reads.clone();
    want_ins.sort();
    let mut want_del = del_reads.clone();
    want_del.sort();

    assert_eq!(reasons_of(&rows, "clip"), want_clip, "{rows:?}");
    assert_eq!(reasons_of(&rows, "large_ins"), want_ins, "{rows:?}");
    assert_eq!(reasons_of(&rows, "large_del"), want_del, "{rows:?}");
    assert_eq!(reasons_of(&rows, "unmapped"), vec!["unplaced1".to_string()], "{rows:?}");

    // A plain read ending exactly at the junction is not a clip.
    assert!(!rows.iter().any(|r| r[2] == "clip" && r[0] == format!("bg{}", JUNCTION - READ_LEN)), "{rows:?}");

    // Clip details and positions.
    let left: Vec<&Vec<String>> = rows.iter().filter(|r| r[8] == "side=left;len=120").collect();
    let right: Vec<&Vec<String>> = rows.iter().filter(|r| r[8] == "side=right;len=120").collect();
    assert_eq!((left.len(), right.len()), (20, 20), "{rows:?}");
    assert!(left.iter().all(|r| r[5] == JUNCTION.to_string()), "{left:?}");
    assert!(right.iter().all(|r| r[5] == (JUNCTION - 1).to_string()), "{right:?}");
    assert!(rows.iter().all(|r| r[2] == "unmapped" || ((r[6] == "+" || r[6] == "-") && r[7] == "primary")), "{rows:?}");

    // Inserts/dels carry their merged length and raw start.
    let ins_rows: Vec<&Vec<String>> = rows.iter().filter(|r| r[2] == "large_ins").collect();
    assert!(ins_rows.iter().all(|r| r[8] == "len=60;at=5000" && r[5] == "5000"), "{ins_rows:?}");
    let del_rows: Vec<&Vec<String>> = rows.iter().filter(|r| r[2] == "large_del").collect();
    assert!(del_rows.iter().all(|r| r[8] == "len=80;at=7000" && r[5] == "7000"), "{del_rows:?}");
    // Unmapped row uses `.` columns.
    let un: &Vec<String> = rows.iter().find(|r| r[2] == "unmapped").unwrap();
    assert_eq!((un[3].as_str(), un[4].as_str(), un[5].as_str(), un[6].as_str(), un[7].as_str()), (".", ".", ".", ".", "."), "{un:?}");

    // JSON reasons and counts.
    let j = json(&p1);
    assert_eq!(j["format_version"], 1, "{j}");
    let raw = std::fs::read_to_string(p1.with_extension("fish.json")).unwrap();
    assert!(raw.starts_with("{\n  \"format_version\": 1,"), "format_version is not the first JSON key:\n{raw}");
    assert_eq!(j["reads"]["by_reason"]["clip"]["reads"], 40, "{j}");
    assert_eq!(j["reads"]["by_reason"]["large_ins"]["reads"], 15, "{j}");
    assert_eq!(j["reads"]["by_reason"]["large_del"]["reads"], 15, "{j}");
    assert_eq!(j["reads"]["by_reason"]["unmapped"]["reads"], 1, "{j}");
    assert_eq!(j["reads"]["reads"], 71, "{j}");

    // Thread identity: identical TSV + decompressed FASTQ; JSON minus timing.
    assert_eq!(std::fs::read(p1.with_extension("fish.tsv")).unwrap(), std::fs::read(p8.with_extension("fish.tsv")).unwrap());
    assert_eq!(fastq(&p1), fastq(&p8), "FASTQ differs between --threads 1 and 8");
    let strip = |path: &Path| -> serde_json::Value {
        let mut v = json(path);
        v.as_object_mut().unwrap().remove("timing_s");
        v
    };
    assert_eq!(strip(&p1), strip(&p8));
}

/// Distinct read -> sorted reasons from the TSV.
fn reasons_by_read(rows: &[Vec<String>]) -> std::collections::BTreeMap<String, Vec<String>> {
    let mut m: std::collections::BTreeMap<String, Vec<String>> = Default::default();
    for r in rows {
        let v = m.entry(r[0].clone()).or_default();
        if !v.contains(&r[2]) {
            v.push(r[2].clone());
        }
    }
    m.values_mut().for_each(|v| v.sort());
    m
}

fn fastq_count(prefix: &Path) -> usize {
    fastq(prefix).iter().filter(|&&b| b == b'\n').count() / 4
}

#[test]
fn fish_clip_cap() {
    let (r, mut recs, _clip, _ins, _del) = build();
    // Reads in the same 100 bp bin as the right-clip junction (pos 1998/1999) that
    // the cap must keep: clip + large insertion (two reasons), and 600 bp clips
    // (eligible now, kept because they are the longest in their bin).
    for k in 0..3 {
        let mut seq = r[1500..1700].to_vec();
        seq.extend_from_slice(&random_seq(INS_LEN, 300 + k));
        seq.extend_from_slice(&r[1700..1940]);
        seq.extend(std::iter::repeat_n(b'A', 120));
        recs.push(RawRec {
            name: format!("mixed{k}"),
            pos: 1500,
            cigar: vec![Cigar::Match(200), Cigar::Ins(INS_LEN as u32), Cigar::Match(240), Cigar::SoftClip(120)],
            seq,
            flags: 0,
        });
        let mut seq = r[1599..1999].to_vec();
        seq.extend(std::iter::repeat_n(b'C', 600));
        recs.push(RawRec {
            name: format!("long{k}"),
            pos: 1599,
            cigar: vec![Cigar::Match(400), Cigar::SoftClip(600)],
            seq,
            flags: 0,
        });
    }

    let off = run("cap_off", &r, &recs, &params(1));
    // A small explicit chunk size keeps the threading comparison meaningful on
    // this 9 kb fixture: at the default automatic chunk size the whole
    // reference is a single task even with 8 workers. The cap equals the number
    // of 600 bp clips in the right-clip bin, so they survive as the longest.
    let capped = |threads: usize, name: &str| {
        run(
            name,
            &r,
            &recs,
            &FishParams {
                clip_cap: 3,
                clip_bin: 100,
                chunk_size: Some(500),
                ..params(threads)
            },
        )
    };
    let on = capped(1, "cap_on_t1");

    let before = reasons_by_read(&tsv(&off));
    let after = reasons_by_read(&tsv(&on));
    let dropped: Vec<&String> = before.keys().filter(|k| !after.contains_key(*k)).collect();
    assert!(!dropped.is_empty());
    // Only clip-only reads disappear; every kept read keeps its reasons.
    for d in &dropped {
        assert_eq!(before[*d], vec!["clip".to_string()], "{d}");
        assert!(d.starts_with("clipR") || d.starts_with("clipL"), "{d}");
    }
    for (k, v) in &after {
        assert_eq!(&before[k], v, "{k}");
    }
    for n in ["mixed0", "mixed1", "mixed2", "long0", "long1", "long2", "unplaced1"] {
        assert!(after.contains_key(n), "{n} must be kept");
    }
    // The 600 bp clips are eligible but outrank the 120 bp clips of the same bin,
    // so all three are kept: every clipR read goes and 17 of the 20 clipL reads.
    let kept_clip: Vec<&str> = after
        .iter()
        .filter(|(k, r)| r.len() == 1 && r[0] == "clip" && k.starts_with("clip"))
        .map(|(k, _)| k.as_str())
        .collect();
    assert_eq!(kept_clip, vec!["clipL0", "clipL1", "clipL10"], "{kept_clip:?}");
    assert_eq!(dropped.iter().filter(|d| d.starts_with("clipR")).count(), 20, "{dropped:?}");
    assert_eq!(dropped.iter().filter(|d| d.starts_with("clipL")).count(), 17, "{dropped:?}");
    assert_eq!(dropped.len(), 37);
    for reason in ["large_ins", "large_del", "unmapped"] {
        assert_eq!(reasons_of(&tsv(&off), reason), reasons_of(&tsv(&on), reason), "{reason}");
    }

    // FASTQ and TSV agree, and the JSON reports the difference.
    assert_eq!(fastq_count(&on), after.len());
    assert_eq!(fastq_count(&off), before.len());
    let (j0, j1) = (json(&off), json(&on));
    assert_eq!(j1["cap"]["reads"].as_u64().unwrap(), (before.len() - after.len()) as u64, "{j1}");
    assert_eq!(j1["cap"]["bp"].as_u64().unwrap(), 37 * 620, "{j1}");
    assert_eq!(j0["cap"]["reads"], 0);
    assert_eq!(j0["cap"]["bp"], 0);
    assert_eq!(j0["params"]["clip_cap"], 0);
    assert_eq!(j0["params"]["clip_bin"], 2000);
    assert_eq!(j1["params"]["clip_cap"], 3);
    assert_eq!(j1["params"]["clip_bin"], 100);
    assert_eq!(j1["reads"]["reads"].as_u64().unwrap(), after.len() as u64);
    assert_eq!(
        j1["reads"]["by_reason"]["clip"]["reads"].as_u64().unwrap(),
        3 + 3 + 3,
        "{j1}"
    );

    // Thread identity with the cap on.
    let on8 = capped(8, "cap_on_t8");
    assert_eq!(
        std::fs::read(on.with_extension("fish.tsv")).unwrap(),
        std::fs::read(on8.with_extension("fish.tsv")).unwrap()
    );
    assert_eq!(fastq(&on), fastq(&on8));
    let strip = |path: &Path| {
        let mut v = json(path);
        v.as_object_mut().unwrap().remove("timing_s");
        v
    };
    assert_eq!(strip(&on), strip(&on8));
}

#[test]
fn fish_clip_bin_zero_is_an_error() {
    let dir = workdir("cap_bin0");
    let r = random_seq(1_000, 3);
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr\n{}\n", std::str::from_utf8(&r).unwrap())).unwrap();
    let err = run_fish(&dir.join("x.bam"), &fa, &dir.join("out"), &FishParams { clip_bin: 0, ..params(1) })
        .unwrap_err()
        .to_string();
    assert!(err.contains("--clip-bin must be positive"), "{err}");
}

#[test]
fn fish_fastq_original_orientation() {
    // A reverse-strand clipped read: SEQ is stored revcomp'd and QUAL reversed;
    // fish must write the read in its original orientation.
    let l = 2_000usize;
    let r = random_seq(l, 17);
    let original = random_seq(560, 31); // 60 clip + 500 match
    let qual: Vec<u8> = (0..560).map(|i| (i % 40) as u8 + 1).collect();
    let seq_rev = revcomp(&original);
    let dir = workdir("rev_orient");
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr\n{}\n", std::str::from_utf8(&r).unwrap())).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), format!("chr\t{l}\t5\t{l}\t{}\n", l + 1)).unwrap();
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", l));
    let bam = dir.join("reads.bam");
    {
        let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();
        let mut rr = Record::new();
        rr.set(b"revclip", Some(&CigarString(vec![Cigar::SoftClip(60), Cigar::Match(500)])), &seq_rev, &qual);
        rr.set_tid(0);
        rr.set_pos(500);
        rr.set_mapq(60);
        rr.set_flags(0x10);
        rr.set_mtid(-1);
        rr.set_mpos(-1);
        w.write(&rr).unwrap();
        drop(w);
    }
    bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
    let prefix = dir.join("out");
    run_fish(&bam, &fa, &prefix, &params(1)).unwrap();

    let rows = tsv(&prefix);
    assert_eq!(reasons_of(&rows, "clip"), vec!["revclip".to_string()], "{rows:?}");
    assert_eq!(rows[0][6], "-", "{rows:?}");

    let fq = String::from_utf8(fastq(&prefix)).unwrap();
    let lines: Vec<&str> = fq.lines().collect();
    assert_eq!(lines[0], "@revclip");
    assert_eq!(lines[1].as_bytes(), original.as_slice(), "SEQ must be the original read");
    let want_q: Vec<u8> = qual.iter().rev().map(|&b| b + 33).collect();
    assert_eq!(lines[3].as_bytes(), want_q.as_slice(), "QUAL must be reversed");
}

#[test]
fn fish_circular_origin_crossing_skipped_linear_counts() {
    // A left clip at pos 0 and a right clip ending at the contig end.
    let l = 2_000usize;
    let r = random_seq(l, 23);
    let mut recs = Vec::new();
    for k in 0..5 {
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(&r[0..500]);
        recs.push(RawRec {
            name: format!("l{k}"),
            pos: 0,
            cigar: vec![Cigar::SoftClip(120), Cigar::Match(500)],
            seq,
            flags: if k % 2 == 0 { 0x10 } else { 0 },
        });
    }
    for k in 0..5 {
        let mut seq = r[1500..2000].to_vec();
        seq.extend(std::iter::repeat_n(b'A', 120));
        recs.push(RawRec {
            name: format!("rr{k}"),
            pos: 1500,
            cigar: vec![Cigar::Match(500), Cigar::SoftClip(120)],
            seq,
            flags: if k % 2 == 1 { 0x10 } else { 0 },
        });
    }

    // Circular (default): both ends are origin crossings -> no clip rows.
    let pc = run("circ", &r, &recs, &params(1));
    assert!(reasons_of(&tsv(&pc), "clip").is_empty(), "{:?}", tsv(&pc));

    // Linear: the contig-end clips count.
    let pl = run("linear", &r, &recs, &FishParams { linear: vec!["chr".to_string()], ..params(1) });
    assert_eq!(reasons_of(&tsv(&pl), "clip").len(), 10, "{:?}", tsv(&pl));
}

#[test]
fn fish_hard_clip_errors() {
    let l = 1_000usize;
    let r = random_seq(l, 29);
    let recs = vec![RawRec {
        name: "hc".to_string(),
        pos: 100,
        cigar: vec![Cigar::HardClip(5), Cigar::Match(500)],
        seq: r[100..600].to_vec(),
        flags: 0,
    }];
    let dir = workdir("hardclip");
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr\n{}\n", std::str::from_utf8(&r).unwrap())).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), format!("chr\t{l}\t5\t{l}\t{}\n", l + 1)).unwrap();
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", l));
    let bam = dir.join("reads.bam");
    {
        let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();
        let mut rec = Record::new();
        rec.set(b"hc", Some(&CigarString(recs[0].cigar.clone())), &recs[0].seq, &vec![30u8; 500]);
        rec.set_tid(0);
        rec.set_pos(100);
        rec.set_mapq(60);
        rec.set_flags(0);
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        w.write(&rec).unwrap();
        drop(w);
    }
    bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
    let prefix = dir.join("out");
    let err = run_fish(&bam, &fa, &prefix, &params(1)).unwrap_err();
    assert!(err.to_string().contains("-Y"), "{err}");
}

#[test]
fn fish_no_flags_disable_categories() {
    let (r, recs, _, _, _) = build();
    let rows = |p: &FishParams, name: &str| tsv(&run(name, &r, &recs, p));

    let none = rows(&FishParams { no_clips: true, no_indels: true, ..params(1) }, "noall");
    assert!(reasons_of(&none, "clip").is_empty());
    assert!(reasons_of(&none, "large_ins").is_empty());
    assert!(reasons_of(&none, "large_del").is_empty());
    assert_eq!(reasons_of(&none, "unmapped").len(), 1);

    let no_unmapped = rows(&FishParams { no_unmapped: true, ..params(1) }, "noun");
    assert!(reasons_of(&no_unmapped, "unmapped").is_empty());
    assert_eq!(reasons_of(&no_unmapped, "clip").len(), 40);
}

/// Write a findings.bed with the clip/ins/del sites for the scope tests.
fn write_findings(dir: &Path, clip_tier: &str) -> PathBuf {
    let p = dir.join("f.findings.bed");
    let text = format!(
        "#contig\tstart\tend\tlabel\tscore\ttier\tlikely\textra\n\
         chr\t1990\t2010\tfinding_1\t1\t{clip_tier}\tmisjoin\tevidence=clip;members=clip_hotspot:1;lowcov=0\n\
         chr\t4990\t5010\tfinding_2\t1\tred\textra_sequence\tevidence=long_ins;members=large_ins:1;lowcov=0\n\
         chr\t6990\t7090\tfinding_3\t1\tred\tmissing_sequence\tevidence=long_del;members=large_del:1;lowcov=0\n"
    );
    std::fs::write(&p, text).unwrap();
    p
}

#[test]
fn fish_findings_scope_and_min_tier() {
    let (r, recs, clip, ins_reads, del_reads) = build();
    let dir = workdir("scope");
    let f = write_findings(&dir, "yellow");

    // Findings scope, min-tier yellow, pad 0: every planted category is inside a
    // finding.
    let p_yellow = FishParams {
        findings: Some(f.clone()),
        min_tier: FishTier::Yellow,
        finding_pad: 0,
        ..params(1)
    };
    let py = run("scope_yellow", &r, &recs, &p_yellow);
    let rows = tsv(&py);
    assert_eq!(reasons_of(&rows, "clip").len(), 40, "{rows:?}");
    assert_eq!(reasons_of(&rows, "large_ins").len(), 15, "{rows:?}");
    assert_eq!(reasons_of(&rows, "large_del").len(), 15, "{rows:?}");
    // The finding column names the finding.
    assert!(rows.iter().any(|r| r[2] == "clip" && r[3] == "finding_1"), "{rows:?}");
    // The JSON reports the resolved scope.
    assert_eq!(json(&py)["params"]["scope"], "findings", "{}", json(&py));

    // min-tier red with pad 0: the yellow clip finding is excluded, the red
    // ins/del findings are kept.
    let p_red = FishParams {
        findings: Some(f),
        min_tier: FishTier::Red,
        finding_pad: 0,
        ..params(1)
    };
    let pr = run("scope_red", &r, &recs, &p_red);
    let rows = tsv(&pr);
    assert!(reasons_of(&rows, "clip").is_empty(), "{rows:?}");
    assert_eq!(reasons_of(&rows, "large_ins").len(), 15, "{rows:?}");
    assert_eq!(reasons_of(&rows, "large_del").len(), 15, "{rows:?}");

    // Genome scope with a real findings file: categories still counted anywhere.
    let f2 = write_findings(&dir, "yellow");
    let p_g = FishParams { findings: Some(f2), scope: Some(FishScope::Genome), ..params(1) };
    let pg = run("scope_genome", &r, &recs, &p_g);
    let rows = tsv(&pg);
    assert_eq!(reasons_of(&rows, "clip").len(), 40, "{rows:?}");
    let mut want_ins = ins_reads;
    want_ins.sort();
    assert_eq!(reasons_of(&rows, "large_ins"), want_ins, "{rows:?}");
    let mut want_del = del_reads;
    want_del.sort();
    assert_eq!(reasons_of(&rows, "large_del"), want_del, "{rows:?}");
    let _ = clip;
}

#[test]
fn fish_linked_and_locus() {
    let l = 4_000usize;
    let r = random_seq(l, 41);
    let mut recs = Vec::new();
    let mut s = 0usize;
    while s + 500 <= l {
        recs.push(RawRec {
            name: format!("bg{s}"),
            pos: s as i64,
            cigar: vec![Cigar::Match(500)],
            seq: r[s..s + 500].to_vec(),
            flags: if (s / 500) % 2 == 1 { 0x10 } else { 0 },
        });
        s += 100;
    }
    let dir = workdir("linked");
    // Haplotype group 1 (reads bg0, bg100, bg200) spans [100, 600); a second
    // group (bg3000) spans [3000, 3500).
    let haps = dir.join("h.haplotypes.tsv");
    std::fs::write(
        &haps,
        "group\tcontig\tspan\tread\tstrand\talleles\tsa\n\
         1\tchr\t100-600\tbg0\t+\t100\t.\n\
         1\tchr\t100-600\tbg100\t-\t100\t.\n\
         1\tchr\t100-600\tbg200\t+\t100\t.\n\
         2\tchr\t3000-3500\tbg3000\t+\t3000\t.\n",
    )
    .unwrap();

    // Genome scope (no findings): every group read is selected.
    let pg = FishParams { haplotypes: Some(haps.clone()), ..params(1) };
    let prefix = run("linked_genome", &r, &recs, &pg);
    let rows = tsv(&prefix);
    let mut got = reasons_of(&rows, "linked");
    got.sort();
    assert_eq!(got, vec!["bg0".to_string(), "bg100".to_string(), "bg200".to_string(), "bg3000".to_string()], "{rows:?}");
    assert!(rows.iter().all(|r| r[2] != "linked" || r[8].starts_with("group=")), "{rows:?}");

    // Findings scope: a red finding over [0, 700) keeps group 1 only.
    let f = dir.join("l.findings.bed");
    std::fs::write(
        &f,
        "chr\t0\t700\tfinding_1\t1\tred\tminority\tevidence=sites;members=SNP:1;lowcov=0\n",
    )
    .unwrap();
    let pf = FishParams { haplotypes: Some(haps), findings: Some(f), ..params(1) };
    let prefix = run("linked_findings", &r, &recs, &pf);
    let rows = tsv(&prefix);
    let mut got = reasons_of(&rows, "linked");
    got.sort();
    assert_eq!(got, vec!["bg0".to_string(), "bg100".to_string(), "bg200".to_string()], "{rows:?}");

    // locus: any record overlapping the finding ± locus_flank.
    let f2 = dir.join("loc.findings.bed");
    std::fs::write(
        &f2,
        "chr\t2000\t2100\tfinding_1\t1\tred\tminority\tevidence=sites;members=SNP:1;lowcov=0\n",
    )
    .unwrap();
    let pl = FishParams {
        findings: Some(f2),
        locus: true,
        locus_flank: 0,
        scope: Some(FishScope::Genome),
        ..params(1)
    };
    let prefix = run("locus", &r, &recs, &pl);
    let rows = tsv(&prefix);
    let locus = reasons_of(&rows, "locus");
    // Records covering [2000, 2100): reads starting every 100 bp whose 500 bp
    // span overlaps it (start 1600..=2000).
    let want: Vec<String> = (16..=20).map(|k| format!("bg{}", k * 100)).collect();
    assert_eq!(locus, want, "{rows:?}");
    assert!(rows.iter().all(|r| r[2] != "locus" || r[3] == "finding_1"), "{rows:?}");

    // --linked needs no haplotypes: skipped, noted in the JSON.
    let pn = FishParams { ..params(1) };
    let prefix = run("linked_skipped", &r, &recs, &pn);
    assert_eq!(json(&prefix)["params"]["linked"], false, "{}", json(&prefix));
    assert!(json(&prefix)["notes"].as_array().unwrap().iter().any(|n| n.as_str().unwrap().contains("haplotypes")), "{}", json(&prefix));
}

// ---------------------------------------------------------------------------
// multi-contig / chapter-boundary fixture and the reviewer-round tests

/// A record for the multi-contig builder.
struct R {
    name: String,
    tid: i32,
    pos: i64,
    cigar: Vec<Cigar>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    flags: u16,
}

fn write_ref_multi(dir: &Path, contigs: &[(String, Vec<u8>)]) -> PathBuf {
    let fa = dir.join("ref.fa");
    let mut s = String::new();
    let mut fai = String::new();
    let mut off = 0usize;
    for (n, seq) in contigs {
        s.push('>');
        s.push_str(n);
        s.push('\n');
        s.push_str(std::str::from_utf8(seq).unwrap());
        s.push('\n');
        fai.push_str(&format!("{n}\t{}\t{}\t{}\t{}\n", seq.len(), off + n.len() + 2, seq.len(), seq.len() + 1));
        off += n.len() + seq.len() + 3;
    }
    std::fs::write(&fa, s).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), fai).unwrap();
    fa
}

/// Write the ref + BAM (or CRAM) and return `(fa, aln)`. The header carries a
/// contract-valid minimap2 `@PG` so `lift` (used by `--bait`) accepts the
/// alignment; `fish` itself ignores `@PG`.
fn write_aln(dir: &Path, contigs: &[(String, Vec<u8>)], recs: &[R], cram: bool) -> (PathBuf, PathBuf) {
    write_aln_cl(dir, contigs, recs, cram, "minimap2 -ax map-hifi -t 1 --secondary=no -Y ref.fa reads.fq")
}

/// Like `write_aln`, with an explicit minimap2 `CL` (to build contract-violation
/// inputs for the bait contract-fallback test).
fn write_aln_cl(
    dir: &Path,
    contigs: &[(String, Vec<u8>)],
    recs: &[R],
    cram: bool,
    cl: &str,
) -> (PathBuf, PathBuf) {
    let fa = write_ref_multi(dir, contigs);
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    for (n, seq) in contigs {
        header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", n).push_tag(b"LN", seq.len()));
    }
    header.push_record(
        HeaderRecord::new(b"PG")
            .push_tag(b"ID", "minimap2")
            .push_tag(b"PN", "minimap2")
            .push_tag(b"VN", "2.31")
            .push_tag(b"CL", cl),
    );
    let aln = dir.join(if cram { "reads.cram" } else { "reads.bam" });
    let fmt = if cram { Format::Cram } else { Format::Bam };
    let mut w = Writer::from_path(&aln, &header, fmt).unwrap();
    if cram {
        w.set_reference(&fa).unwrap();
    }
    let mut sorted: Vec<&R> = recs.iter().collect();
    sorted.sort_by_key(|r| (r.tid < 0 || r.pos < 0, r.tid, r.pos));
    for r in sorted {
        let mut rec = Record::new();
        rec.set(r.name.as_bytes(), Some(&CigarString(r.cigar.clone())), &r.seq, &r.qual);
        rec.set_tid(r.tid);
        rec.set_pos(r.pos);
        rec.set_mapq(60);
        rec.set_flags(r.flags);
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        w.write(&rec).unwrap();
    }
    drop(w);
    if cram {
        let exe = std::env::var("CRAMHOOK_SAMTOOLS").unwrap_or_else(|_| "samtools".to_string());
        let st = Command::new(exe).arg("index").arg(&aln).status().expect("samtools index");
        assert!(st.success(), "samtools index failed");
    } else {
        bam::index::build(&aln, None, bam::index::Type::Bai, 1).unwrap();
    }
    (fa, aln)
}

/// Write a single-contig BAM from `R` records and run fish.
fn run_r(name: &str, contigs: &[(String, Vec<u8>)], recs: &[R], p: &FishParams) -> PathBuf {
    let dir = workdir(name);
    let (fa, aln) = write_aln(&dir, contigs, recs, false);
    let prefix = dir.join("out");
    run_fish(&aln, &fa, &prefix, p).unwrap();
    prefix
}

/// Parse a decompressed FASTQ into name -> (seq, qual).
fn fastq_entries(prefix: &Path) -> std::collections::BTreeMap<String, (String, String)> {
    let text = String::from_utf8(fastq(prefix)).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len() % 4, 0, "FASTQ not a multiple of 4 lines:\n{text}");
    let mut m = std::collections::BTreeMap::new();
    for chunk in lines.chunks(4) {
        assert!(chunk[0].starts_with('@'));
        m.insert(chunk[0][1..].to_string(), (chunk[1].to_string(), chunk[3].to_string()));
    }
    m
}

/// BAM and CRAM of the same records must give the same FASTQ (finding 1):
/// CRAM must decode QUAL, not fall back to `!`.
#[test]
fn fish_bam_cram_fastq_equivalence() {
    let contigs = vec![("chr".to_string(), random_seq(2_000, 61))];
    let r = &contigs[0].1;
    let qual: Vec<u8> = (0..560).map(|i| (i % 40) as u8 + 1).collect();
    let mut recs = Vec::new();
    for (k, rev) in [(0usize, false), (1, true)] {
        let seq = if rev { revcomp(&random_seq(560, 70 + k as u64)) } else { random_seq(560, 70 + k as u64) };
        recs.push(R {
            name: format!("c{k}"),
            tid: 0,
            pos: 500,
            cigar: vec![Cigar::SoftClip(60), Cigar::Match(500)],
            seq,
            qual: qual.clone(),
            flags: if rev { 0x10 } else { 0 },
        });
    }
    recs.push(R {
        name: "plain".to_string(),
        tid: 0,
        pos: 0,
        cigar: vec![Cigar::Match(500)],
        seq: r[0..500].to_vec(),
        qual: vec![30u8; 500],
        flags: 0,
    });

    let bam_dir = workdir("bam_cram_bam");
    let (fa_b, bam) = write_aln(&bam_dir, &contigs, &recs, false);
    let bam_prefix = bam_dir.join("out");
    run_fish(&bam, &fa_b, &bam_prefix, &params(1)).unwrap();

    let cram_dir = workdir("bam_cram_cram");
    let (fa_c, cram) = write_aln(&cram_dir, &contigs, &recs, true);
    let cram_prefix = cram_dir.join("out");
    run_fish(&cram, &fa_c, &cram_prefix, &params(1)).unwrap();

    let bf = fastq_entries(&bam_prefix);
    let cf = fastq_entries(&cram_prefix);
    assert_eq!(bf, cf, "BAM and CRAM FASTQ differ");
    // QUAL is real, not the all-`!` fallback.
    let (_, q) = &cf["c0"];
    assert!(q.chars().any(|c| c != '!'), "CRAM QUAL collapsed to all-`!`: {q:?}");
    assert_eq!(
        std::fs::read(bam_prefix.with_extension("fish.tsv")).unwrap(),
        std::fs::read(cram_prefix.with_extension("fish.tsv")).unwrap()
    );
}

/// Findings ids from column 4 are preserved, not renumbered (finding 4).
#[test]
fn fish_findings_ids_preserved() {
    let l = 3_000usize;
    let r = random_seq(l, 71);
    let mut recs = Vec::new();
    for k in 0..8 {
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(&r[1500..2000]);
        recs.push(R {
            name: format!("c{k}"),
            tid: 0,
            pos: 1500,
            cigar: vec![Cigar::SoftClip(120), Cigar::Match(500)],
            seq,
            qual: vec![30u8; 620],
            flags: if k % 2 == 1 { 0x10 } else { 0 },
        });
    }
    let dir = workdir("find_ids");
    let f = dir.join("f.findings.bed");
    std::fs::write(
        &f,
        "chr\t1400\t1600\tfinding_5\t1\tred\tmisjoin\tevidence=clip;members=clip_hotspot:1;lowcov=0\n\
         chr\t2500\t2600\tfinding_9\t1\tyellow\tminority\tevidence=site_minor;members=DEL:1;lowcov=0\n",
    )
    .unwrap();
    let p = FishParams { findings: Some(f), finding_pad: 0, ..params(1) };
    let prefix = run_r("find_ids_run", &[("chr".to_string(), r.clone())], &recs, &p);
    let rows = tsv(&prefix);
    assert_eq!(reasons_of(&rows, "clip").len(), 8, "{rows:?}");
    assert!(rows.iter().all(|x| x[3] == "finding_5"), "ids must be preserved, not renumbered: {rows:?}");
    assert!(!rows.iter().any(|x| x[3] == "finding_1" || x[3] == "finding_2"), "{rows:?}");
}

/// With `--scope genome`, `--findings` must not filter haplotype groups (finding 3).
#[test]
fn fish_genome_scope_ignores_findings_for_groups() {
    let l = 4_000usize;
    let r = random_seq(l, 81);
    let mut recs = Vec::new();
    let mut s = 0usize;
    while s + 500 <= l {
        recs.push(R {
            name: format!("bg{s}"),
            tid: 0,
            pos: s as i64,
            cigar: vec![Cigar::Match(500)],
            seq: r[s..s + 500].to_vec(),
            qual: vec![30u8; 500],
            flags: if (s / 500) % 2 == 1 { 0x10 } else { 0 },
        });
        s += 500;
    }
    let dir = workdir("genome_vs_findings");
    let haps = dir.join("h.tsv");
    std::fs::write(
        &haps,
        "group\tcontig\tspan\tread\tstrand\talleles\tsa\n\
         1\tchr\t100-600\tbg0\t+\t100\t.\n\
         2\tchr\t3000-3500\tbg3000\t+\t3000\t.\n",
    )
    .unwrap();
    let f = dir.join("f.bed");
    std::fs::write(&f, "chr\t0\t700\tfinding_1\t1\tred\tminority\tevidence=sites;members=SNP:1;lowcov=0\n").unwrap();

    // Findings scope: only the group overlapping the finding.
    let pf = FishParams { haplotypes: Some(haps.clone()), findings: Some(f.clone()), ..params(1) };
    let rows = tsv(&run_r("gvf_findings", &[("chr".to_string(), r.clone())], &recs, &pf));
    assert_eq!(reasons_of(&rows, "linked"), vec!["bg0".to_string()], "{rows:?}");

    // Genome scope: findings given, but every group counts.
    let pg = FishParams {
        haplotypes: Some(haps),
        findings: Some(f),
        scope: Some(FishScope::Genome),
        ..params(1)
    };
    let rows = tsv(&run_r("gvf_genome", &[("chr".to_string(), r.clone())], &recs, &pg));
    assert_eq!(reasons_of(&rows, "linked"), vec!["bg0".to_string(), "bg3000".to_string()], "{rows:?}");
}

/// Coordinate-less records (pos < 0), including tid >= 0 ones, are reached
/// exactly once (finding 5).
#[test]
fn fish_coordinate_less_unmapped_reached() {
    let contigs = vec![("chr".to_string(), random_seq(1_000, 91))];
    let recs = vec![
        R {
            name: "m1".to_string(),
            tid: 0,
            pos: 100,
            cigar: vec![Cigar::Match(500)],
            seq: vec![b'A'; 500],
            qual: vec![30u8; 500],
            flags: 0,
        },
        // flag-4 with a coordinate: returned by the region fetch.
        R {
            name: "placed".to_string(),
            tid: 0,
            pos: 200,
            cigar: Vec::new(),
            seq: vec![b'C'; 100],
            qual: vec![30u8; 100],
            flags: 0x4,
        },
        // coordinate-less: only the `"*"` iterator returns it.
        R {
            name: "unp".to_string(),
            tid: -1,
            pos: -1,
            cigar: Vec::new(),
            seq: vec![b'G'; 100],
            qual: vec![30u8; 100],
            flags: 0x4,
        },
    ];
    let dir = workdir("coordless");
    let (fa, aln) = write_aln(&dir, &contigs, &recs, false);
    let prefix = dir.join("out");
    run_fish(&aln, &fa, &prefix, &params(1)).unwrap();
    let rows = tsv(&prefix);
    assert_eq!(reasons_of(&rows, "unmapped"), vec!["placed".to_string(), "unp".to_string()], "{rows:?}");
    assert_eq!(json(&prefix)["records"]["unmapped"], 2, "{}", json(&prefix));
}

/// `--no-linked` disables linked selection (finding 7).
#[test]
fn fish_no_linked_flag() {
    let l = 1_000usize;
    let r = random_seq(l, 101);
    let recs = vec![R {
        name: "bg0".to_string(),
        tid: 0,
        pos: 0,
        cigar: vec![Cigar::Match(500)],
        seq: r[0..500].to_vec(),
        qual: vec![30u8; 500],
        flags: 0,
    }];
    let dir = workdir("no_linked");
    let haps = dir.join("h.tsv");
    std::fs::write(
        &haps,
        "group\tcontig\tspan\tread\tstrand\talleles\tsa\n1\tchr\t100-600\tbg0\t+\t100\t.\n",
    )
    .unwrap();
    let p = FishParams { haplotypes: Some(haps), no_linked: true, ..params(1) };
    let prefix = run_r("no_linked_run", &[("chr".to_string(), r.clone())], &recs, &p);
    assert!(reasons_of(&tsv(&prefix), "linked").is_empty(), "{:?}", tsv(&prefix));
    assert_eq!(json(&prefix)["params"]["linked"], false, "{}", json(&prefix));
}

/// Absent QUAL becomes `!` (finding 7); SEQ is kept.
#[test]
fn fish_missing_qual_becomes_bang() {
    let contigs = vec![("chr".to_string(), random_seq(2_000, 111))];
    let recs = vec![R {
        name: "nq".to_string(),
        tid: 0,
        pos: 500,
        cigar: vec![Cigar::SoftClip(60), Cigar::Match(500)],
        seq: random_seq(560, 112),
        qual: vec![0xffu8; 560],
        flags: 0,
    }];
    let dir = workdir("missing_qual");
    let (fa, aln) = write_aln(&dir, &contigs, &recs, false);
    let prefix = dir.join("out");
    run_fish(&aln, &fa, &prefix, &params(1)).unwrap();
    let fq = fastq_entries(&prefix);
    let (seq, qual) = &fq["nq"];
    assert_eq!(qual.len(), seq.len());
    assert!(qual.chars().all(|c| c == '!'), "{qual:?}");
}

/// Secondary / QC-fail / duplicate records are not counted (finding 7).
#[test]
fn fish_counted_record_exclusions() {
    let contigs = vec![("chr".to_string(), random_seq(2_000, 121))];
    let mk = |name: &str, flags: u16| R {
        name: name.to_string(),
        tid: 0,
        pos: 500,
        cigar: vec![Cigar::SoftClip(60), Cigar::Match(500)],
        seq: random_seq(560, 122),
        qual: vec![30u8; 560],
        flags,
    };
    let recs = vec![mk("sec", 0x100), mk("qf", 0x200), mk("dup", 0x400)];
    let dir = workdir("exclusions");
    let (fa, aln) = write_aln(&dir, &contigs, &recs, false);
    let prefix = dir.join("out");
    run_fish(&aln, &fa, &prefix, &params(1)).unwrap();
    assert!(tsv(&prefix).is_empty(), "{:?}", tsv(&prefix));
    let j = json(&prefix);
    assert_eq!((j["records"]["secondary"].as_u64(), j["records"]["qcfail"].as_u64(), j["records"]["duplicate"].as_u64()), (Some(1), Some(1), Some(1)), "{j}");
    assert!(String::from_utf8(fastq(&prefix)).unwrap().is_empty(), "no reads selected");
}

/// Thread identity across multiple chunks/contigs, a read split across contigs
/// (primary + supplementary), boundary records and shared names (finding 7).
#[test]
fn fish_thread_identity_multichunk() {
    let chr1 = random_seq(2_000, 131);
    let chr2 = random_seq(2_000, 132);
    let contigs = vec![("chr1".to_string(), chr1.clone()), ("chr2".to_string(), chr2.clone())];
    let mut recs: Vec<R> = Vec::new();
    // Background ~every 250 bp on both contigs (not selected).
    for (tid, contig) in [(0i32, &chr1), (1, &chr2)] {
        let mut s = 0;
        while s + 400 <= contig.len() {
            recs.push(R {
                name: format!("g{tid}_{s}"),
                tid,
                pos: s as i64,
                cigar: vec![Cigar::Match(400)],
                seq: contig[s..s + 400].to_vec(),
                qual: vec![30u8; 400],
                flags: 0,
            });
            s += 250;
        }
    }
    // A read split across the two contigs: primary on chr1, supplementary with a
    // left clip on chr2 -> selected via the supplementary, payload from the primary.
    let shared_seq = chr1[1000..1500].to_vec();
    recs.push(R {
        name: "shared".to_string(),
        tid: 0,
        pos: 1000,
        cigar: vec![Cigar::Match(500)],
        seq: shared_seq.clone(),
        qual: vec![30u8; 500],
        flags: 0,
    });
    recs.push(R {
        name: "shared".to_string(),
        tid: 1,
        pos: 500,
        cigar: vec![Cigar::SoftClip(60), Cigar::Match(440)],
        seq: random_seq(500, 133),
        qual: vec![30u8; 500],
        flags: 0x800,
    });
    // Reads starting exactly on chunk boundaries (chunk_size = 400 below).
    for (i, b) in [400i64, 800, 1200, 1600].iter().enumerate() {
        let mut seq = vec![b'A'; 60];
        seq.extend_from_slice(&chr1[*b as usize..*b as usize + 400]);
        recs.push(R {
            name: format!("b{i}"),
            tid: 0,
            pos: *b,
            cigar: vec![Cigar::SoftClip(60), Cigar::Match(400)],
            seq,
            qual: vec![30u8; 460],
            flags: if i % 2 == 1 { 0x10 } else { 0 },
        });
    }

    let p = |threads: usize| FishParams { threads, chunk_size: Some(400), ..Default::default() };
    let dir = workdir("multichunk");
    let (fa, aln) = write_aln(&dir, &contigs, &recs, false);
    let p1 = dir.join("t1");
    let p4 = dir.join("t4");
    run_fish(&aln, &fa, &p1, &p(1)).unwrap();
    run_fish(&aln, &fa, &p4, &p(4)).unwrap();

    // The environment has >= 4 cores in CI? clamp: only compare if threads take.
    assert_eq!(std::fs::read(p1.with_extension("fish.tsv")).unwrap(), std::fs::read(p4.with_extension("fish.tsv")).unwrap(), "TSV differs");
    assert_eq!(fastq(&p1), fastq(&p4), "FASTQ differs across threads");
    let strip = |path: &Path| {
        let mut v = json(path);
        v.as_object_mut().unwrap().remove("timing_s");
        v
    };
    assert_eq!(strip(&p1), strip(&p4), "JSON differs across threads");

    // The split read: selected once, sequence from the primary record.
    let fq = fastq_entries(&p1);
    let (seq, _q) = fq.get("shared").expect("shared read selected");
    assert_eq!(seq.as_bytes(), shared_seq.as_slice(), "payload must come from the primary record");
    assert_eq!(fq.keys().filter(|k| *k == "shared").count(), 1);
    // Exactly one row for the shared read, reason clip.
    let rows = tsv(&p1);
    assert_eq!(rows.iter().filter(|x| x[0] == "shared").count(), 1, "{rows:?}");
    assert_eq!(rows.iter().find(|x| x[0] == "shared").unwrap()[2], "clip");
}

/// A flag-4 record with a set tid and `pos = -1` cannot be returned by the
/// region or `"*"` iterators. `write_aln` can build+index this layout as a CRAM
/// (a BAM cannot be indexed), and fish must fail the completeness check instead
/// of silently dropping it (finding 1, second round).
#[test]
fn fish_exotic_coordinate_less_layout_errors() {
    let contigs = vec![("chr".to_string(), random_seq(1_000, 141))];
    let recs = vec![
        R {
            name: "m1".to_string(),
            tid: 0,
            pos: 100,
            cigar: vec![Cigar::Match(500)],
            seq: vec![b'A'; 500],
            qual: vec![30u8; 500],
            flags: 0,
        },
        R {
            name: "c1".to_string(),
            tid: 0,
            pos: -1,
            cigar: Vec::new(),
            seq: vec![b'C'; 100],
            qual: vec![30u8; 100],
            flags: 0x4,
        },
    ];
    let dir = workdir("exotic");
    let (fa, aln) = write_aln(&dir, &contigs, &recs, true);
    let prefix = dir.join("out");
    let err = run_fish(&aln, &fa, &prefix, &params(1)).unwrap_err().to_string();
    assert!(err.contains("index reports"), "{err}");
}

/// A linked row is attributed to a finding overlapping the group span: the
/// tier-qualified finding in the findings scope, any finding in the genome scope
/// (finding 2, second round).
#[test]
fn fish_linked_finding_attribution() {
    let l = 4_000usize;
    let r = random_seq(l, 151);
    let recs = vec![R {
        name: "r0".to_string(),
        tid: 0,
        pos: 1_000,
        cigar: vec![Cigar::Match(400)],
        seq: r[1_000..1_400].to_vec(),
        qual: vec![30u8; 400],
        flags: 0,
    }];
    let dir = workdir("linked_finding");
    let haps = dir.join("h.tsv");
    std::fs::write(
        &haps,
        "group\tcontig\tspan\tread\tstrand\talleles\tsa\n1\tchr\t1000-1100\tr0\t+\t1000\t.\n",
    )
    .unwrap();
    let f = dir.join("f.bed");
    // Only finding_17 overlaps; the id from column 4 must be reported verbatim.
    std::fs::write(
        &f,
        "chr\t900\t1200\tfinding_17\t1\tred\tminority\tevidence=sites;members=SNP:1;lowcov=0\n",
    )
    .unwrap();

    // Findings scope: linked-only selection, attributed to finding_17.
    let pf = FishParams { haplotypes: Some(haps.clone()), findings: Some(f.clone()), ..params(1) };
    let rows = tsv(&run_r("linked_finding_findings", &[("chr".to_string(), r.clone())], &recs, &pf));
    assert_eq!(rows.len(), 1, "{rows:?}");
    assert_eq!((rows[0][2].as_str(), rows[0][3].as_str()), ("linked", "finding_17"), "{rows:?}");

    // Genome scope: the overlapping finding is still attributed.
    let pg = FishParams {
        haplotypes: Some(haps),
        findings: Some(f),
        scope: Some(FishScope::Genome),
        ..params(1)
    };
    let rows = tsv(&run_r("linked_finding_genome", &[("chr".to_string(), r)], &recs, &pg));
    assert_eq!(rows.len(), 1, "{rows:?}");
    assert_eq!((rows[0][2].as_str(), rows[0][3].as_str()), ("linked", "finding_17"), "{rows:?}");
}

/// Rows are keyed by (read, reason, finding): a read with a small clip inside a
/// finding and a larger clip outside gets one row per finding id (including `.`)
/// in the genome scope; the findings scope keeps only the inside event.
#[test]
fn fish_genome_scope_lists_every_finding() {
    let l = 4_000usize;
    let r = random_seq(l, 161);
    let recs = vec![
        // Primary: a 50 bp left clip at 1000, inside the finding.
        R {
            name: "r".to_string(),
            tid: 0,
            pos: 1_000,
            cigar: vec![Cigar::SoftClip(50), Cigar::Match(400)],
            seq: random_seq(450, 162),
            qual: vec![30u8; 450],
            flags: 0,
        },
        // Supplementary: a larger 300 bp left clip at 3000, outside every finding.
        R {
            name: "r".to_string(),
            tid: 0,
            pos: 3_000,
            cigar: vec![Cigar::SoftClip(300), Cigar::Match(100)],
            seq: random_seq(400, 163),
            qual: vec![30u8; 400],
            flags: 0x800,
        },
    ];
    let dir = workdir("every_finding");
    let f = dir.join("f.bed");
    std::fs::write(&f, "chr\t950\t1100\tfinding_1\t1\tred\tminority\tevidence=sites;members=SNP:1;lowcov=0\n").unwrap();
    let contigs = [("chr".to_string(), r)];

    // Genome scope: two clip rows for the read, one per finding id (`.` last).
    let pg = FishParams { findings: Some(f.clone()), scope: Some(FishScope::Genome), finding_pad: 0, ..params(1) };
    let pg_prefix = run_r("every_finding_genome", &contigs, &recs, &pg);
    let rows = tsv(&pg_prefix);
    let clip: Vec<&Vec<String>> = rows.iter().filter(|x| x[0] == "r" && x[2] == "clip").collect();
    assert_eq!(clip.len(), 2, "{rows:?}");
    assert_eq!((clip[0][3].as_str(), clip[0][1].as_str()), ("finding_1", "50"), "{rows:?}");
    assert_eq!((clip[1][3].as_str(), clip[1][1].as_str()), (".", "300"), "{rows:?}");
    assert_eq!(json(&pg_prefix)["reads"]["by_reason"]["clip"]["reads"], 1, "distinct reads, not rows");

    // Findings scope: only the event inside the finding.
    let pf = FishParams { findings: Some(f), finding_pad: 0, ..params(1) };
    let rows = tsv(&run_r("every_finding_findings", &contigs, &recs, &pf));
    let clip: Vec<&Vec<String>> = rows.iter().filter(|x| x[0] == "r" && x[2] == "clip").collect();
    assert_eq!(clip.len(), 1, "{rows:?}");
    assert_eq!((clip[0][3].as_str(), clip[0][1].as_str()), ("finding_1", "50"), "{rows:?}");
}

/// Evidence is offered to every applicable finding: two overlapping (padded)
/// findings both get a row, not just the smallest id.
#[test]
fn fish_overlapping_findings_both_rows() {
    let l = 4_000usize;
    let r = random_seq(l, 171);
    let recs = vec![
        R {
            name: "r".to_string(),
            tid: 0,
            pos: 100,
            cigar: vec![Cigar::SoftClip(60), Cigar::Match(200)],
            seq: random_seq(260, 172),
            qual: vec![30u8; 260],
            flags: 0,
        },
        R {
            name: "r".to_string(),
            tid: 0,
            pos: 500,
            cigar: vec![Cigar::SoftClip(100), Cigar::Match(200)],
            seq: random_seq(300, 173),
            qual: vec![30u8; 300],
            flags: 0x800,
        },
    ];
    let dir = workdir("overlapping");
    let f = dir.join("f.bed");
    std::fs::write(
        &f,
        "chr\t90\t120\tfinding_2\t1\tred\tminority\tevidence=clip;members=clip_hotspot:1;lowcov=0\n\
         chr\t490\t520\tfinding_10\t1\tred\tminority\tevidence=clip;members=clip_hotspot:1;lowcov=0\n",
    )
    .unwrap();
    // Default pad (1000): both clips fall in both findings.
    let p = FishParams { findings: Some(f), finding_pad: 1_000, ..params(1) };
    let rows = tsv(&run_r("overlapping_run", &[("chr".to_string(), r)], &recs, &p));
    let clip: Vec<&Vec<String>> = rows.iter().filter(|x| x[0] == "r" && x[2] == "clip").collect();
    assert_eq!(clip.len(), 2, "{rows:?}");
    assert_eq!((clip[0][3].as_str(), clip[1][3].as_str()), ("finding_2", "finding_10"), "{rows:?}");
}

/// One locus record spanning two findings is listed once per finding
/// (`--locus --locus-flank 0 --finding-pad 0`).
#[test]
fn fish_locus_spans_two_findings() {
    let l = 4_000usize;
    let r = random_seq(l, 181);
    let recs = vec![R {
        name: "r".to_string(),
        tid: 0,
        pos: 100,
        cigar: vec![Cigar::Match(600)],
        seq: r[100..700].to_vec(),
        qual: vec![30u8; 600],
        flags: 0,
    }];
    let dir = workdir("locus_two");
    let f = dir.join("f.bed");
    std::fs::write(
        &f,
        "chr\t90\t120\tfinding_2\t1\tred\tminority\tevidence=sites;members=SNP:1;lowcov=0\n\
         chr\t490\t520\tfinding_10\t1\tred\tminority\tevidence=sites;members=SNP:1;lowcov=0\n",
    )
    .unwrap();
    let p = FishParams { findings: Some(f), locus: true, locus_flank: 0, finding_pad: 0, ..params(1) };
    let rows = tsv(&run_r("locus_two_run", &[("chr".to_string(), r)], &recs, &p));
    let locus: Vec<&Vec<String>> = rows.iter().filter(|x| x[2] == "locus").collect();
    assert_eq!(locus.len(), 2, "{rows:?}");
    assert_eq!((locus[0][3].as_str(), locus[1][3].as_str()), ("finding_2", "finding_10"), "{rows:?}");
}

/// A haplotype group spanning two findings is listed once per finding.
#[test]
fn fish_linked_spans_two_findings() {
    let l = 4_000usize;
    let r = random_seq(l, 191);
    let recs = vec![R {
        name: "r".to_string(),
        tid: 0,
        pos: 100,
        cigar: vec![Cigar::Match(200)],
        seq: r[100..300].to_vec(),
        qual: vec![30u8; 200],
        flags: 0,
    }];
    let dir = workdir("linked_two");
    let haps = dir.join("h.tsv");
    std::fs::write(
        &haps,
        "group\tcontig\tspan\tread\tstrand\talleles\tsa\n1\tchr\t100-600\tr\t+\t100\t.\n",
    )
    .unwrap();
    let f = dir.join("f.bed");
    std::fs::write(
        &f,
        "chr\t90\t120\tfinding_2\t1\tred\tminority\tevidence=sites;members=SNP:1;lowcov=0\n\
         chr\t490\t520\tfinding_10\t1\tred\tminority\tevidence=sites;members=SNP:1;lowcov=0\n",
    )
    .unwrap();
    let p = FishParams { haplotypes: Some(haps), findings: Some(f), ..params(1) };
    let rows = tsv(&run_r("linked_two_run", &[("chr".to_string(), r)], &recs, &p));
    let linked: Vec<&Vec<String>> = rows.iter().filter(|x| x[2] == "linked").collect();
    assert_eq!(linked.len(), 2, "{rows:?}");
    assert_eq!((linked[0][3].as_str(), linked[1][3].as_str()), ("finding_2", "finding_10"), "{rows:?}");
}

// ---------------------------------------------------------------------------
// bait (`--bait`, lift-based)

/// Write a bait FASTA into `dir`; returns its path.
fn write_bait_fa(dir: &Path, fname: &str, contigs: &[(String, Vec<u8>)]) -> PathBuf {
    let p = dir.join(fname);
    let mut s = String::new();
    for (n, seq) in contigs {
        s.push('>');
        s.push_str(n);
        s.push('\n');
        s.push_str(std::str::from_utf8(seq).unwrap());
        s.push('\n');
    }
    std::fs::write(&p, s).unwrap();
    p
}

/// Write reads as FASTQ (`name`, SEQ, QUAL as Phred+33).
fn write_fastq(path: &Path, reads: &[(String, Vec<u8>, Vec<u8>)]) {
    let mut s: Vec<u8> = Vec::new();
    for (name, seq, qual) in reads {
        s.extend_from_slice(b"@");
        s.extend_from_slice(name.as_bytes());
        s.push(b'\n');
        s.extend_from_slice(seq);
        s.push(b'\n');
        s.extend_from_slice(b"+\n");
        s.extend_from_slice(qual);
        s.push(b'\n');
    }
    std::fs::write(path, s).unwrap();
}

/// Align `fastq` to `fa` with minimap2 + samtools sort into a sorted, indexed
/// BAM whose `@PG` is contract-valid (so `lift` accepts it).
fn align_reads(fa: &Path, fastq: &Path, bam: &Path) {
    let mut mm = Command::new("minimap2")
        .args(["-ax", "map-hifi", "-t", "1", "--secondary=no", "-Y"])
        .arg(fa)
        .arg(fastq)
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn minimap2");
    let sort = Command::new("samtools")
        .args(["sort", "-O", "BAM", "-o"])
        .arg(bam)
        .stdin(Stdio::from(mm.stdout.take().unwrap()))
        .status()
        .expect("spawn samtools sort");
    assert!(sort.success(), "samtools sort failed");
    assert!(mm.wait().unwrap().success(), "minimap2 failed");
    bam::index::build(bam, None, bam::index::Type::Bai, 1).unwrap();
}

/// Build ref + BAM, write one bait FASTA, run fish with `--bait`; returns the
/// output prefix.
fn run_bait_fish(
    name: &str,
    contigs: &[(String, Vec<u8>)],
    recs: &[R],
    bait: &[(String, Vec<u8>)],
    p: &FishParams,
) -> PathBuf {
    let dir = workdir(name);
    let bait_fa = write_bait_fa(&dir, "bait.fa", bait);
    let (fa, aln) = write_aln(&dir, contigs, recs, false);
    let prefix = dir.join("out");
    let params = FishParams { bait: vec![bait_fa], ..p.clone() };
    run_fish(&aln, &fa, &prefix, &params).unwrap();
    prefix
}

/// JSON with the `timing_s` fields (top level and `bait`) removed, for
/// thread-identity comparisons.
fn strip_timing(path: &Path) -> serde_json::Value {
    let mut v = json(path);
    v.as_object_mut().unwrap().remove("timing_s");
    if let Some(b) = v.get_mut("bait").and_then(|x| x.as_object_mut()) {
        b.remove("timing_s");
    }
    v
}

/// Serializes the bait tests: one of them sets `CRAMHOOK_MINIMAP2`, which is
/// process-global.
fn bait_lock() -> MutexGuard<'static, ()> {
    static LOCK: Mutex<()> = Mutex::new(());
    LOCK.lock().unwrap_or_else(|e| e.into_inner())
}

/// A real minimap2 input: reads from the assembly translate, reads from the bait
/// (absent from the assembly) are remapped onto it and selected with reason
/// `bait`. Their FASTQ payloads are the original reads; outputs are
/// thread-identical.
#[test]
fn fish_bait_translates_and_selects() {
    let _g = bait_lock();
    let a = random_seq(9_000, 301);
    let b = random_seq(3_000, 302);
    let mut reads: Vec<(String, Vec<u8>, Vec<u8>)> = Vec::new();
    let mut want: Vec<String> = Vec::new();
    let mut s = 0usize;
    while s + 500 <= a.len() {
        reads.push((format!("a{s}"), a[s..s + 500].to_vec(), vec![b'I'; 500]));
        s += 250;
    }
    let mut j = 0usize;
    while j + 500 <= b.len() {
        let n = format!("b{j}");
        reads.push((n.clone(), b[j..j + 500].to_vec(), vec![b'I'; 500]));
        want.push(n);
        j += 250;
    }
    want.sort();

    let dir = workdir("bait_translate");
    let contigs = vec![("chr".to_string(), a)];
    let fa = write_ref_multi(&dir, &contigs);
    let fq = dir.join("reads.fq");
    write_fastq(&fq, &reads);
    let aln = dir.join("reads.bam");
    align_reads(&fa, &fq, &aln);
    let bait_fa = write_bait_fa(&dir, "bait.fa", &[("phage".to_string(), b)]);

    let p = |threads: usize| FishParams {
        threads,
        no_unmapped: true,
        bait: vec![bait_fa.clone()],
        // Do not take the tiny-input fallback: the assembly reads must translate.
        bait_min_input_mb: Some(0.0),
        ..Default::default()
    };
    let p1 = dir.join("t1");
    run_fish(&aln, &fa, &p1, &p(1)).unwrap();
    let p4 = dir.join("t4");
    run_fish(&aln, &fa, &p4, &p(4)).unwrap();

    let rows = tsv(&p1);
    let mut got = reasons_of(&rows, "bait");
    got.sort();
    assert_eq!(got, want, "{rows:?}");
    assert!(reasons_of(&rows, "unmapped").is_empty(), "{rows:?}");
    let bait_rows: Vec<&Vec<String>> = rows.iter().filter(|r| r[2] == "bait").collect();
    assert!(bait_rows.iter().all(|r| r[4] == "phage" && r[7] == "primary"), "{bait_rows:?}");
    assert!(bait_rows.iter().all(|r| r[8].starts_with("contig=phage;mapq=")), "{bait_rows:?}");

    let j = json(&p1);
    assert_eq!(j["params"]["bait"], true, "{j}");
    assert_eq!(j["bait"]["reads"], want.len() as u64, "{j}");
    assert_eq!(j["bait"]["contigs"][0]["name"], "phage", "{j}");
    assert_eq!(j["bait"]["contigs"][0]["length"], 3_000, "{j}");
    assert_eq!(j["bait"]["contigs"][0]["reads"], want.len() as u64, "{j}");
    assert!(j["bait"]["translated_reads"].as_u64().unwrap() > 0, "no translation: {j}");
    assert_eq!(j["bait"]["fallback"], serde_json::Value::Null, "{j}");
    assert_eq!(j["bait"]["contract_ok"], true, "{j}");

    // Exact original payloads for the bait-selected (unmapped-in-assembly) reads.
    let fq_out = fastq_entries(&p1);
    for (name, seq, q) in &reads {
        if want.contains(name) {
            let (gs, gq) = fq_out.get(name).unwrap_or_else(|| panic!("{name} not selected"));
            assert_eq!(gs.as_bytes(), seq.as_slice(), "{name} SEQ");
            assert_eq!(gq.as_bytes(), q.as_slice(), "{name} QUAL");
        }
    }

    assert_eq!(
        std::fs::read(p1.with_extension("fish.tsv")).unwrap(),
        std::fs::read(p4.with_extension("fish.tsv")).unwrap(),
        "TSV differs"
    );
    assert_eq!(fastq(&p1), fastq(&p4), "FASTQ differs");
    assert_eq!(strip_timing(&p1), strip_timing(&p4), "JSON differs");
}

/// A bait identical to a region of the assembly (prophage case): only reads of
/// that region may be selected, they fit the assembly equally well (MAPQ 0), the
/// selection is non-empty, and reads far from the region are never selected.
#[test]
fn fish_bait_prophage_region_mapq0() {
    let _g = bait_lock();
    let a = random_seq(9_000, 401);
    let b = a[3_000..6_000].to_vec();
    let mut recs: Vec<R> = Vec::new();
    let mut region: Vec<String> = Vec::new();
    let mut far: Vec<String> = Vec::new();
    let mut s = 0usize;
    while s + 500 <= a.len() {
        let n = format!("r{s}");
        if s >= 3_000 && s + 500 <= 6_000 {
            region.push(n.clone());
        } else {
            far.push(n.clone());
        }
        recs.push(R {
            name: n,
            tid: 0,
            pos: s as i64,
            cigar: vec![Cigar::Match(500)],
            seq: a[s..s + 500].to_vec(),
            qual: vec![30u8; 500],
            flags: if (s / 250) % 2 == 1 { 0x10 } else { 0 },
        });
        s += 250;
    }
    let contigs = vec![("chr".to_string(), a)];
    let bait = vec![("prophage".to_string(), b)];
    let p = FishParams { threads: 1, bait_min_input_mb: Some(0.0), ..Default::default() };
    let prefix = run_bait_fish("bait_prophage", &contigs, &recs, &bait, &p);

    let rows = tsv(&prefix);
    let bait_rows: Vec<&Vec<String>> = rows.iter().filter(|r| r[2] == "bait").collect();
    assert!(!bait_rows.is_empty(), "no prophage read selected: {rows:?}");
    assert!(bait_rows.iter().all(|r| r[4] == "prophage"), "{bait_rows:?}");
    // Only reads of the shared region can be selected...
    assert!(bait_rows.iter().all(|r| region.contains(&r[0])), "{bait_rows:?}");
    // ...and they fit the assembly (the integrated copy) equally well: MAPQ 0.
    assert!(bait_rows.iter().all(|r| r[8].contains(";mapq=0;")), "{bait_rows:?}");
    // Reads far from the shared region are never selected.
    assert!(!far.iter().any(|f| bait_rows.iter().any(|r| &r[0] == f)), "{bait_rows:?}");

    let j = json(&prefix);
    assert!(j["bait"]["reads"].as_u64().unwrap() > 0, "{j}");
    assert_eq!(j["bait"]["mapq0_reads"], j["bait"]["reads"], "{j}");
}

/// The FASTQ payload of a bait-selected read is the original read: an unmapped
/// read (flag 4) and a reverse read (flag 0x10, SEQ stored revcomp'd) both come
/// back in their original orientation.
#[test]
fn fish_bait_payload_orientation() {
    let _g = bait_lock();
    let a = random_seq(4_000, 601);
    let b = random_seq(2_000, 602);
    let orig_f = b[100..600].to_vec();
    let orig_r = b[900..1_400].to_vec();
    let q_f: Vec<u8> = (0..500).map(|i| (i % 40) as u8 + 1).collect();
    let q_r: Vec<u8> = (0..500).map(|i| ((i + 7) % 40) as u8 + 1).collect();
    let mut q_r_rec = q_r.clone();
    q_r_rec.reverse();
    let recs = vec![
        R {
            name: "bg".to_string(),
            tid: 0,
            pos: 0,
            cigar: vec![Cigar::Match(500)],
            seq: a[0..500].to_vec(),
            qual: vec![30u8; 500],
            flags: 0,
        },
        R {
            name: "fwd".to_string(),
            tid: -1,
            pos: -1,
            cigar: Vec::new(),
            seq: orig_f.clone(),
            qual: q_f.clone(),
            flags: 0x4,
        },
        R {
            name: "rev".to_string(),
            tid: -1,
            pos: -1,
            cigar: Vec::new(),
            seq: revcomp(&orig_r),
            qual: q_r_rec,
            flags: 0x4 | 0x10,
        },
    ];
    let contigs = vec![("chr".to_string(), a)];
    let bait = vec![("phage".to_string(), b)];
    let p = FishParams { threads: 1, no_unmapped: true, ..Default::default() };
    let prefix = run_bait_fish("bait_payload", &contigs, &recs, &bait, &p);

    let rows = tsv(&prefix);
    let mut got = reasons_of(&rows, "bait");
    got.sort();
    assert_eq!(got, vec!["fwd".to_string(), "rev".to_string()], "{rows:?}");

    let fq = fastq_entries(&prefix);
    let (fs, fq_q) = &fq["fwd"];
    assert_eq!(fs.as_bytes(), orig_f.as_slice(), "forward SEQ");
    let fq_q_want: Vec<u8> = q_f.iter().map(|&x| x + 33).collect();
    assert_eq!(fq_q.as_bytes(), fq_q_want.as_slice(), "forward QUAL");
    let (rs, rq) = &fq["rev"];
    assert_eq!(rs.as_bytes(), orig_r.as_slice(), "reverse SEQ must be the original read");
    let rq_want: Vec<u8> = q_r.iter().map(|&x| x + 33).collect();
    assert_eq!(rq.as_bytes(), rq_want.as_slice(), "reverse QUAL");
}

/// A read whose original primary is duplicate (0x400) or QC-fail (0x200) is not
/// selected by `bait`: `lift` remaps it and the remapped record loses the flag.
#[test]
fn fish_bait_skips_duplicate_and_qcfail_primaries() {
    let _g = bait_lock();
    let a = random_seq(4_000, 611);
    let b = random_seq(2_000, 612);
    let mk = |name: &str, off: usize, flags: u16| R {
        name: name.to_string(),
        tid: -1,
        pos: -1,
        cigar: Vec::new(),
        seq: b[off..off + 500].to_vec(),
        qual: vec![30u8; 500],
        flags,
    };
    let recs = vec![
        R {
            name: "bg".to_string(),
            tid: 0,
            pos: 0,
            cigar: vec![Cigar::Match(500)],
            seq: a[0..500].to_vec(),
            qual: vec![30u8; 500],
            flags: 0,
        },
        mk("ok", 100, 0x4),
        mk("dup", 600, 0x4 | 0x400),
        mk("qf", 1_100, 0x4 | 0x200),
    ];
    let contigs = vec![("chr".to_string(), a)];
    let bait = vec![("phage".to_string(), b)];
    let p = FishParams { threads: 1, no_unmapped: true, ..Default::default() };
    let prefix = run_bait_fish("bait_dup_qcfail", &contigs, &recs, &bait, &p);

    let rows = tsv(&prefix);
    assert_eq!(reasons_of(&rows, "bait"), vec!["ok".to_string()], "{rows:?}");
    let j = json(&prefix);
    assert_eq!((j["records"]["duplicate"].as_u64(), j["records"]["qcfail"].as_u64()), (Some(1), Some(1)), "{j}");
    assert_eq!(j["bait"]["reads"], 1, "{j}");
}

/// `--bait` is repeatable: reads are selected per bait contig and the JSON lists
/// every bait contig.
#[test]
fn fish_bait_multiple_files() {
    let _g = bait_lock();
    let a = random_seq(4_000, 621);
    let b1 = random_seq(2_000, 622);
    let b2 = random_seq(2_000, 623);
    let recs = vec![
        R {
            name: "bg".to_string(),
            tid: 0,
            pos: 0,
            cigar: vec![Cigar::Match(500)],
            seq: a[0..500].to_vec(),
            qual: vec![30u8; 500],
            flags: 0,
        },
        R {
            name: "r1".to_string(),
            tid: -1,
            pos: -1,
            cigar: Vec::new(),
            seq: b1[100..600].to_vec(),
            qual: vec![30u8; 500],
            flags: 0x4,
        },
        R {
            name: "r2".to_string(),
            tid: -1,
            pos: -1,
            cigar: Vec::new(),
            seq: b2[900..1_400].to_vec(),
            qual: vec![30u8; 500],
            flags: 0x4,
        },
    ];
    let dir = workdir("bait_multi");
    let contigs = vec![("chr".to_string(), a)];
    let (fa, aln) = write_aln(&dir, &contigs, &recs, false);
    let f1 = write_bait_fa(&dir, "b1.fa", &[("plasmid1".to_string(), b1)]);
    let f2 = write_bait_fa(&dir, "b2.fa", &[("plasmid2".to_string(), b2)]);
    let prefix = dir.join("out");
    let p = FishParams {
        threads: 1,
        no_unmapped: true,
        bait: vec![f1, f2],
        ..Default::default()
    };
    run_fish(&aln, &fa, &prefix, &p).unwrap();

    let rows = tsv(&prefix);
    let mut got = reasons_of(&rows, "bait");
    got.sort();
    assert_eq!(got, vec!["r1".to_string(), "r2".to_string()], "{rows:?}");
    let bait_rows: Vec<&Vec<String>> = rows.iter().filter(|r| r[2] == "bait").collect();
    let r1 = bait_rows.iter().find(|r| r[0] == "r1").unwrap();
    let r2 = bait_rows.iter().find(|r| r[0] == "r2").unwrap();
    assert_eq!(r1[4], "plasmid1", "{bait_rows:?}");
    assert_eq!(r2[4], "plasmid2", "{bait_rows:?}");

    let j = json(&prefix);
    let names: Vec<String> =
        j["bait"]["contigs"].as_array().unwrap().iter().map(|c| c["name"].as_str().unwrap().to_string()).collect();
    assert_eq!(names, vec!["plasmid1".to_string(), "plasmid2".to_string()], "{j}");
    assert_eq!(j["bait"]["reads"], 2, "{j}");
}

/// A contract violation in the input header is reported in the bait JSON block;
/// `lift` still remaps every read, so the bait selection is unchanged.
#[test]
fn fish_bait_contract_fallback_reported() {
    let _g = bait_lock();
    let a = random_seq(4_000, 631);
    let b = random_seq(2_000, 632);
    let recs = vec![
        R {
            name: "bg".to_string(),
            tid: 0,
            pos: 0,
            cigar: vec![Cigar::Match(500)],
            seq: a[0..500].to_vec(),
            qual: vec![30u8; 500],
            flags: 0,
        },
        R {
            name: "r".to_string(),
            tid: -1,
            pos: -1,
            cigar: Vec::new(),
            seq: b[100..600].to_vec(),
            qual: vec![30u8; 500],
            flags: 0x4,
        },
    ];
    let dir = workdir("bait_contract");
    let contigs = vec![("chr".to_string(), a)];
    let (fa, aln) = write_aln_cl(
        &dir,
        &contigs,
        &recs,
        false,
        "minimap2 -ax map-hifi -t 1 --secondary=yes -Y ref.fa reads.fq",
    );
    let bait_fa = write_bait_fa(&dir, "bait.fa", &[("phage".to_string(), b)]);
    let prefix = dir.join("out");
    let p = FishParams { threads: 1, no_unmapped: true, bait: vec![bait_fa], ..Default::default() };
    run_fish(&aln, &fa, &prefix, &p).unwrap();

    let rows = tsv(&prefix);
    assert_eq!(reasons_of(&rows, "bait"), vec!["r".to_string()], "{rows:?}");
    let j = json(&prefix);
    assert_eq!(j["bait"]["contract_ok"], false, "{j}");
    assert!(
        j["bait"]["contract_reason"].as_str().unwrap_or("").contains("secondary=yes"),
        "{j}"
    );
    assert_eq!(j["bait"]["fallback"], serde_json::Value::Null, "{j}");
}

/// A same-prefix rerun with a renamed bait must not reuse the stale `.fai`, and
/// a failed validation (name clash) must leave the previous outputs intact.
#[test]
fn fish_bait_rerun_and_failed_validation() {
    let _g = bait_lock();
    let a = random_seq(4_000, 641);
    let b1 = random_seq(2_000, 642);
    let b2 = random_seq(2_000, 643);
    let recs = vec![R {
        name: "bg".to_string(),
        tid: 0,
        pos: 0,
        cigar: vec![Cigar::Match(500)],
        seq: a[0..500].to_vec(),
        qual: vec![30u8; 500],
        flags: 0,
    }];
    let dir = workdir("bait_rerun");
    let contigs = vec![("chr".to_string(), a)];
    let (fa, aln) = write_aln(&dir, &contigs, &recs, false);
    let f1 = write_bait_fa(&dir, "b1.fa", &[("phage1".to_string(), b1)]);
    let f2 = write_bait_fa(&dir, "b2.fa", &[("phage2".to_string(), b2)]);
    let clash = write_bait_fa(&dir, "clash.fa", &[("chr".to_string(), random_seq(300, 644))]);
    let prefix = dir.join("out");
    let run = |bait: Vec<PathBuf>| {
        let p = FishParams { threads: 1, bait, ..Default::default() };
        run_fish(&aln, &fa, &prefix, &p)
    };

    // First run: bait "phage1".
    run(vec![f1]).unwrap();
    let ref1 = std::fs::read(prefix.with_extension("bait_ref.fa")).unwrap();
    assert!(String::from_utf8_lossy(&ref1).contains(">phage1"));

    // Rerun with a renamed bait: must succeed (no stale .fai) and update the ref.
    run(vec![f2]).unwrap();
    let ref2 = std::fs::read(prefix.with_extension("bait_ref.fa")).unwrap();
    assert!(String::from_utf8_lossy(&ref2).contains(">phage2"), "ref not updated");
    assert!(!String::from_utf8_lossy(&ref2).contains(">phage1"));
    let fai = std::fs::read_to_string(prefix.with_extension("bait_ref.fa.fai")).unwrap();
    assert!(fai.contains("phage2\t2000\t"), "stale fai: {fai}");
    assert!(!fai.contains("phage1"), "stale fai: {fai}");
    let cram2 = std::fs::read(prefix.with_extension("bait.cram")).unwrap();
    let tsv2 = std::fs::read(prefix.with_extension("fish.tsv")).unwrap();

    // Failed validation (name clash): the previous outputs must be untouched.
    let err = run(vec![clash]).unwrap_err().to_string();
    assert!(err.contains("clash"), "{err}");
    assert_eq!(std::fs::read(prefix.with_extension("bait_ref.fa")).unwrap(), ref2, "ref clobbered");
    assert_eq!(std::fs::read(prefix.with_extension("bait.cram")).unwrap(), cram2, "cram clobbered");
    assert_eq!(std::fs::read(prefix.with_extension("fish.tsv")).unwrap(), tsv2, "tsv clobbered");
}

/// A failed bait lift (an injected failing minimap2) must restore the previous
/// `P.bait_ref.fa` + `.fai`, so the previous `P.bait.cram` still decodes against
/// it (its `@SQ` MD5 keeps matching).
#[test]
fn fish_bait_failed_lift_restores_reference() {
    let _g = bait_lock();
    let a = random_seq(4_000, 651);
    let b1 = random_seq(2_000, 652);
    let b2 = random_seq(2_000, 653);
    let recs = vec![
        R {
            name: "bg".to_string(),
            tid: 0,
            pos: 0,
            cigar: vec![Cigar::Match(500)],
            seq: a[0..500].to_vec(),
            qual: vec![30u8; 500],
            flags: 0,
        },
        R {
            name: "r".to_string(),
            tid: -1,
            pos: -1,
            cigar: Vec::new(),
            seq: b1[100..600].to_vec(),
            qual: vec![30u8; 500],
            flags: 0x4,
        },
    ];
    let dir = workdir("bait_lift_fail");
    let contigs = vec![("chr".to_string(), a)];
    let (fa, aln) = write_aln(&dir, &contigs, &recs, false);
    let f1 = write_bait_fa(&dir, "b1.fa", &[("phage1".to_string(), b1.clone())]);
    let f2 = write_bait_fa(&dir, "b2.fa", &[("phage2".to_string(), b2)]);
    let prefix = dir.join("out");

    // First run succeeds and produces the reference + CRAM to be preserved.
    let p1 = FishParams { threads: 1, bait: vec![f1], ..Default::default() };
    run_fish(&aln, &fa, &prefix, &p1).unwrap();
    let ref1 = std::fs::read(prefix.with_extension("bait_ref.fa")).unwrap();
    let fai1 = std::fs::read(prefix.with_extension("bait_ref.fa.fai")).unwrap();
    let cram1 = std::fs::read(prefix.with_extension("bait.cram")).unwrap();

    // Inject a minimap2 that fails immediately (the same hook lift's tests use).
    let script = dir.join("mm_fail.sh");
    std::fs::write(&script, "#!/bin/sh\nexit 43\n").unwrap();
    let mut perm = std::fs::metadata(&script).unwrap().permissions();
    perm.set_mode(0o755);
    std::fs::set_permissions(&script, perm).unwrap();

    unsafe { std::env::set_var("CRAMHOOK_MINIMAP2", &script) };
    let p2 = FishParams { threads: 1, bait: vec![f2], ..Default::default() };
    let res = run_fish(&aln, &fa, &prefix, &p2);
    unsafe { std::env::remove_var("CRAMHOOK_MINIMAP2") };
    let err = format!("{:#}", res.expect_err("the injected minimap2 failure must fail fish"));
    assert!(err.contains("minimap2 exited with status 43"), "{err}");

    // The previous reference + index are back, and the previous CRAM still
    // decodes against them (a wrong reference would corrupt the decoded SEQ).
    let ref_now = std::fs::read(prefix.with_extension("bait_ref.fa")).unwrap();
    let fai_now = std::fs::read(prefix.with_extension("bait_ref.fa.fai")).unwrap();
    let cram_now = std::fs::read(prefix.with_extension("bait.cram")).unwrap();
    assert!(ref_now == ref1, "ref not restored ({} vs {} bytes)", ref_now.len(), ref1.len());
    assert!(fai_now == fai1, "fai not restored ({} vs {} bytes)", fai_now.len(), fai1.len());
    assert!(cram_now == cram1, "cram changed ({} vs {} bytes)", cram_now.len(), cram1.len());
    assert!(
        !prefix.with_extension("bait_ref.fa.prev").exists(),
        "left a .prev reference behind"
    );

    let mut reader = bam::Reader::from_path(prefix.with_extension("bait.cram")).unwrap();
    reader.set_reference(prefix.with_extension("bait_ref.fa")).unwrap();
    let want_f = b1[100..600].to_vec();
    let want_r = revcomp(&want_f);
    let mut seen = false;
    for r in rust_htslib::bam::Read::records(&mut reader) {
        let r = r.expect("previous CRAM must decode against the restored reference");
        if r.qname() == b"r" {
            let got = r.seq().as_bytes();
            assert!(got == want_f || got == want_r, "decoded SEQ mismatch: {got:?}");
            seen = true;
        }
    }
    assert!(seen, "the bait-selected read is missing from the previous CRAM");
}

/// A wrapped finding is written as two BED rows with the same id (blocker 1);
/// both rows must attribute reads, so an origin-side read is not lost to a
/// linear comparison. The two clips here sit at the two ends of one wrapped
/// finding.
#[test]
fn fish_wrapped_finding_origin_row_attributed() {
    let l = 12_000usize;
    let r = random_seq(l, 99);
    let mut recs = Vec::new();
    for k in 0..6 {
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(&r[0..500]);
        recs.push(R {
            name: format!("o{k}"),
            tid: 0,
            pos: 0,
            cigar: vec![Cigar::SoftClip(120), Cigar::Match(500)],
            seq,
            qual: vec![30u8; 620],
            flags: if k % 2 == 1 { 0x10 } else { 0 },
        });
    }
    for k in 0..6 {
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(&r[11_500..12_000]);
        recs.push(R {
            name: format!("f{k}"),
            tid: 0,
            pos: 11_500,
            cigar: vec![Cigar::SoftClip(120), Cigar::Match(500)],
            seq,
            qual: vec![30u8; 620],
            flags: 0,
        });
    }
    let dir = workdir("wrapped_finding");
    let f = dir.join("f.findings.bed");
    // One wrapped finding: two rows, same id `finding_1`.
    std::fs::write(
        &f,
        "chr\t11500\t12000\tfinding_1\t1\tred\tmisjoin\tevidence=clip;members=clip_hotspot:1;lowcov=0\n\
         chr\t0\t500\tfinding_1\t1\tred\tmisjoin\tevidence=clip;members=clip_hotspot:1;lowcov=0\n",
    )
    .unwrap();
    let p = FishParams {
        findings: Some(f),
        linear: vec!["chr".to_string()],
        min_tier: FishTier::Yellow,
        finding_pad: 0,
        ..params(1)
    };
    let rows = tsv(&run_r("wrapped_finding_run", &[("chr".to_string(), r)], &recs, &p));
    let clips: Vec<&Vec<String>> = rows.iter().filter(|x| x[2] == "clip").collect();
    assert_eq!(clips.len(), 12, "{rows:?}");
    // The origin-side reads (pos 0) are attributed to the same finding as the
    // far-side ones; without the second row they would be lost.
    assert!(clips.iter().all(|x| x[3] == "finding_1"), "{rows:?}");
}

/// `--bait-records`: a read whose only bait record is supplementary is selected
/// under `any` (the default) but not under `primary`; a read whose primary record
/// is on the bait is selected under both. An absent flag and an explicit `any`
/// behave identically; `P.fish.json` records the mode and the supplementary-only
/// count in every mode.
#[test]
fn fish_bait_records_any_vs_primary() {
    let _g = bait_lock();
    let a = random_seq(4_000, 701);
    let b = random_seq(2_000, 702);
    // A chimeric read: 800 bp of the assembly followed by 300 bp of the bait.
    // Aligned to the assembly alone the bait part is a 300 bp soft clip; remapped
    // onto assembly + bait it yields a primary on `chr` and a *supplementary* on
    // the bait (the 800 bp side scores higher), so it reaches the bait only
    // through the supplementary record.
    let mut chim = a[0..800].to_vec();
    chim.extend_from_slice(&b[0..300]);
    // A read present only in the bait: its primary record is on the bait.
    let plain = b[600..1_100].to_vec();
    let recs = vec![
        R {
            name: "bg".to_string(),
            tid: 0,
            pos: 1_500,
            cigar: vec![Cigar::Match(500)],
            seq: a[1_500..2_000].to_vec(),
            qual: vec![30u8; 500],
            flags: 0,
        },
        R {
            name: "chim".to_string(),
            tid: 0,
            pos: 0,
            cigar: vec![Cigar::Match(800), Cigar::SoftClip(300)],
            seq: chim.clone(),
            qual: vec![30u8; 1_100],
            flags: 0,
        },
        R {
            name: "plain".to_string(),
            tid: -1,
            pos: -1,
            cigar: Vec::new(),
            seq: plain.clone(),
            qual: vec![30u8; 500],
            flags: 0x4,
        },
    ];
    let contigs = vec![("chr".to_string(), a)];
    let bait = vec![("phage".to_string(), b)];
    let base = FishParams { threads: 1, no_unmapped: true, no_clips: true, ..Default::default() };

    // Default (no flag): both reads are selected and the JSON records `any`.
    let default_prefix = run_bait_fish("bait_records_default", &contigs, &recs, &bait, &base);
    let mut got = reasons_of(&tsv(&default_prefix), "bait");
    got.sort();
    assert_eq!(got, vec!["chim".to_string(), "plain".to_string()], "{got:?}");
    let default_rows = tsv(&default_prefix);
    let chim_row = default_rows.iter().find(|r| r[0] == "chim").unwrap();
    assert_eq!(chim_row[7], "supplementary", "{chim_row:?}");
    let jd = json(&default_prefix);
    assert_eq!(jd["bait"]["mode"], "any", "{jd}");
    assert_eq!(jd["bait"]["supplementary_only_reads"], 1, "{jd}");

    // Explicit `any`: behaves identically (identical JSON apart from timing).
    let p_any = FishParams { bait_records: BaitRecords::Any, ..base.clone() };
    let any_prefix = run_bait_fish("bait_records_any", &contigs, &recs, &bait, &p_any);
    assert_eq!(reasons_of(&tsv(&any_prefix), "bait"), got, "explicit any must match default");
    assert_eq!(
        std::fs::read(default_prefix.with_extension("fish.tsv")).unwrap(),
        std::fs::read(any_prefix.with_extension("fish.tsv")).unwrap(),
        "default and explicit any TSV differ"
    );
    assert_eq!(fastq(&default_prefix), fastq(&any_prefix), "default and explicit any FASTQ differ");
    assert_eq!(strip_timing(&default_prefix), strip_timing(&any_prefix), "default and explicit any JSON differ");

    // Explicit `primary`: only the read whose primary record is on the bait.
    let p_primary = FishParams { bait_records: BaitRecords::Primary, ..base.clone() };
    let primary_prefix = run_bait_fish("bait_records_primary", &contigs, &recs, &bait, &p_primary);
    assert_eq!(reasons_of(&tsv(&primary_prefix), "bait"), vec!["plain".to_string()]);
    let jp = json(&primary_prefix);
    assert_eq!(jp["bait"]["mode"], "primary", "{jp}");
    // The diagnostic is the same in both modes: `chim` was supplementary-only.
    assert_eq!(jp["bait"]["supplementary_only_reads"], 1, "{jp}");
    assert_eq!(jp["bait"]["reads"], 1, "{jp}");

    // The dropped read is absent from the primary run's FASTQ.
    let fq = fastq_entries(&primary_prefix);
    assert!(fq.contains_key("plain"), "{:?}", fq.keys().collect::<Vec<_>>());
    assert!(!fq.contains_key("chim"), "chim must not be selected under primary");
}
