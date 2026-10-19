//! Integration test for the `qc` linkage second pass (`src/qc/link.rs`).
//!
//! A synthetic 60× reference with three SNP neighbourhoods:
//!  * a "second strain" carrying 6 SNPs within 3 kb in 30 % of the reads ->
//!    one orange `linked_haplotype` region, the 6 sites marked `LINKED` and the
//!    second-strain reads in `P.haplotypes.tsv`;
//!  * 6 sites with independent 15 % random errors -> no group, sites `UNLINKED`;
//!  * 2 sites with strand-specific errors on forward reads only -> no group.
//! `--threads 1` and 8 must give identical output.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{self, Format, Header, Record, Writer};

use cramhook::{Preset, QcParams, run_qc};

const L: usize = 24_000;
const READ_LEN: usize = 1_200;
/// ~100× coverage. Sites are close enough that every tested pair passes the
/// per-strand check (the conflict guard rejects a merge when an underpowered
/// pair is not linked).
const STEP: usize = 12;
/// Region 1: 6 linked sites within 1 kb, spaced so no 500 bp window holds
/// 5 sites (which would form a `snp_cluster`).
const S1: [usize; 6] = [1_000, 1_130, 1_260, 1_390, 1_520, 1_650];
/// Region 2: 6 sites with independent random errors.
const S2: [usize; 6] = [9_000, 9_300, 9_600, 9_900, 10_200, 10_500];
/// Region 3: 2 sites with strand-specific errors (forward reads only).
const S3: [usize; 2] = [16_000, 16_400];

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

fn alt_of(b: u8) -> u8 {
    match b {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

fn workdir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/qc-linkage/test").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Read index -> true if it is a "second strain" read (carries the 6 SNPs).
fn is_second(i: usize) -> bool {
    i % 10 < 3
}

fn build(dir: &Path) -> (PathBuf, PathBuf) {
    let mut rng = Lcg(11);
    let refseq: Vec<u8> = (0..L).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect();
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr\n{}\n", std::str::from_utf8(&refseq).unwrap())).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), format!("chr\t{L}\t5\t{L}\t{}\n", L + 1)).unwrap();

    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", L));
    let bam = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();

    let mut rng = Lcg(99);
    let n_reads = (L - READ_LEN) / STEP + 1;
    let planted = |p: usize| S1.contains(&p) || S2.contains(&p) || S3.contains(&p);
    for i in 0..n_reads {
        let s = i * STEP;
        let rev = i % 2 == 1;
        let mut seq = refseq[s..s + READ_LEN].to_vec();
        let covers = |p: usize| s <= p && p < s + READ_LEN;
        if is_second(i) {
            for &p in &S1 {
                if covers(p) {
                    seq[p - s] = alt_of(refseq[p]);
                }
            }
        }
        // Region 2: independent (disjoint, strand-balanced) 15 % errors. The
        // residue class per site keeps each site's alt reads on both strands and
        // shares no read between two sites, so no pair can link.
        for (j, &p) in S2.iter().enumerate() {
            if covers(p) && i % 7 == (2 * j) % 7 {
                seq[p - s] = alt_of(refseq[p]);
            }
        }
        if !rev && rng.next() % 2 == 0 {
            for &p in &S3 {
                if covers(p) {
                    seq[p - s] = alt_of(refseq[p]);
                }
            }
        }
        // Background 0.5 % substitution noise everywhere except the planted
        // sites, so the SNP background fit behaves like a real dataset.
        for (k, b) in seq.iter_mut().enumerate() {
            let p = s + k;
            if !planted(p) && rng.next() % 200 == 0 {
                *b = b"ACGT"[(rng.next() % 4) as usize];
            }
        }
        let qual = vec![30u8; seq.len()];
        let mut rec = Record::new();
        rec.set(format!("r{i}").as_bytes(), Some(&CigarString(vec![Cigar::Match(READ_LEN as u32)])), &seq, &qual);
        rec.set_tid(0);
        rec.set_pos(s as i64);
        rec.set_mapq(60);
        rec.set_flags(if rev { 0x10 } else { 0 });
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        w.write(&rec).unwrap();
    }
    drop(w);
    bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
    (bam, fa)
}

fn params(threads: usize, chunk: Option<u64>) -> QcParams {
    QcParams { preset: Some(Preset::MapHifi), threads, candidates: true, chunk_size: chunk, ..Default::default() }
}

fn read_text(p: &Path) -> String {
    std::fs::read_to_string(p).unwrap()
}

#[test]
fn qc_linkage_groups_and_unlinked() {
    let dir = workdir("main");
    let (bam, fa) = build(&dir);
    let p1 = dir.join("t1/out");
    let p8 = dir.join("t8/out");
    // `--unlinked-action flag` (the default is `drop`, tested below) so the
    // unlinked sites stay in the VCF and their UNLINKED flag can be checked.
    let flag = |threads, chunk| QcParams { unlinked_action: cramhook::qc::UnlinkedAction::Flag, ..params(threads, chunk) };
    run_qc(&bam, &fa, &p1, &flag(1, None)).unwrap();
    let s8 = run_qc(&bam, &fa, &p8, &flag(8, Some(1_500))).unwrap();
    assert!(s8.chunks > 10);

    // --- BED: exactly one orange linked_haplotype region over region 1.
    let bed = read_text(&p1.with_extension("bed"));
    let lh: Vec<Vec<&str>> =
        bed.lines().filter(|l| l.contains("\tlinked_haplotype\t")).map(|l| l.split('\t').collect()).collect();
    assert_eq!(lh.len(), 1, "linked_haplotype regions:\n{bed}");
    let r = &lh[0];
    assert_eq!((r[1], r[2], r[5]), ("1000", "1651", "orange"), "{r:?}");
    assert!(r[9].contains("group=1;"), "{r:?}");
    assert!(r[9].contains("sites=6;"), "{r:?}");
    // The group reads are the second-strain reads covering >= 2 sites.
    let want_reads: BTreeSet<String> = (0..(L - READ_LEN) / STEP + 1)
        .filter(|&i| is_second(i) && S1.iter().filter(|&&p| i * STEP <= p && p < i * STEP + READ_LEN).count() >= 2)
        .map(|i| format!("r{i}"))
        .collect();
    assert!(r[9].contains(&format!("reads={};", want_reads.len())), "{r:?} want {}", want_reads.len());

    // --- VCF: region 1 sites LINKED, region 2 sites UNLINKED, region 3 neither.
    let vcf: Vec<Vec<String>> = read_text(&p1.with_extension("vcf"))
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').map(|s| s.to_string()).collect())
        .collect();
    let at = |pos: usize| vcf.iter().find(|v| v[1] == (pos + 1).to_string()).cloned();
    for &p in &S1 {
        let v = at(p).unwrap_or_else(|| panic!("no VCF at {p}"));
        assert!(v[7].contains(";LINKED=1"), "site {p}: {v:?}");
        assert!(!v[7].contains(";UNLINKED"), "site {p}: {v:?}");
    }
    for &p in &S2 {
        let v = at(p).unwrap_or_else(|| panic!("no VCF at {p}"));
        assert!(v[7].contains(";UNLINKED"), "site {p}: {v:?}");
        assert!(!v[7].contains(";LINKED="), "site {p}: {v:?}");
    }
    for &p in &S3 {
        // Strand-specific errors: the site may not even be called (strand gate),
        // but it must never be linked.
        if let Some(v) = at(p) {
            assert!(!v[7].contains(";LINKED=") && !v[7].contains(";UNLINKED"), "site {p}: {v:?}");
        }
    }
    // No linked_haplotype region over region 3.
    assert!(!bed.lines().any(|l| l.contains("\tlinked_haplotype\t") && l.starts_with("chr\t16")));

    // --- haplotypes.tsv: one row per (group 1, second-strain read covering >= 2 sites).
    let hap = read_text(&dir.join("t1/out.haplotypes.tsv"));
    let rows: Vec<Vec<&str>> = hap.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    let names: BTreeSet<String> = rows.iter().filter(|r| r[0] == "1").map(|r| r[3].to_string()).collect();
    assert_eq!(names, want_reads, "haplotypes.tsv read names:\n{hap}");
    assert!(rows.iter().all(|r| r[1] == "chr" && r[2] == "1000-1651"), "{rows:?}");
    // Every row carries at least two site positions.
    for r in &rows {
        assert!(r[5].split(',').count() >= 2, "{r:?}");
    }

    // --- JSON linkage block.
    let j: serde_json::Value = serde_json::from_str(&read_text(&p1.with_extension("json"))).unwrap();
    let lk = &j["linkage"];
    // Background noise creates sub-threshold sites that chain regions together,
    // so the region count is not exact; the group and site counts are.
    assert!(lk["regions"].as_u64().unwrap() >= 1, "{lk}");
    assert_eq!(lk["groups"], 1, "{lk}");
    assert_eq!(lk["linked_sites"], 6, "{lk}");
    assert!(lk["unlinked_sites"].as_u64().unwrap() >= 6, "{lk}");
    assert_eq!(j["alerts"]["by_type"]["linked_haplotype"]["orange"], 1, "{}", j["alerts"]);

    // --- Findings: the linked group is one finding with likely=subpopulation.
    let fbed = read_text(&p1.with_extension("findings.bed"));
    let f: Vec<Vec<&str>> = fbed
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    let linked = f.iter().find(|x| x[7].contains("linked_haplotype:")).unwrap_or_else(|| panic!("{fbed}"));
    assert_eq!(linked[6], "subpopulation", "{linked:?}");
    assert!(linked[7].contains("linked_haplotype:1"), "{linked:?}");
    assert!(linked[7].contains("linked"), "{linked:?}");

    // --- Threads 1 vs 8 identical.
    for suf in ["vcf", "bed", "haplotypes.tsv", "findings.bed", "candidates.tsv"] {
        let a = std::fs::read(dir.join(format!("t1/out.{suf}"))).unwrap();
        let b = std::fs::read(dir.join(format!("t8/out.{suf}"))).unwrap();
        assert!(a == b, "{suf} differs between --threads 1 and 8");
    }
}

/// `--unlinked-action drop`: the tested-but-unlinked sites leave the VCF but
/// stay in candidates.tsv with tier none and reason `unlinked`.
#[test]
fn qc_unlinked_action_drop() {
    let dir = workdir("drop");
    let (bam, fa) = build(&dir);
    let prefix = dir.join("out");
    run_qc(&bam, &fa, &prefix, &QcParams { unlinked_action: cramhook::qc::UnlinkedAction::Drop, ..params(1, None) })
        .unwrap();
    let vcf = read_text(&prefix.with_extension("vcf"));
    for &p in &S2 {
        let pos = (p + 1).to_string();
        assert!(!vcf.lines().any(|l| !l.starts_with('#') && l.split('\t').nth(1) == Some(&pos)), "S2 site {p} still in VCF");
    }
    // The linked sites are unaffected.
    for &p in &S1 {
        let pos = (p + 1).to_string();
        assert!(vcf.lines().any(|l| l.contains(&format!("\t{pos}\t")) && l.contains("LINKED=1")), "S1 site {p} missing");
    }
    let cands = read_text(&dir.join("out.candidates.tsv"));
    let header: Vec<&str> = cands.lines().next().unwrap().split('\t').collect();
    let ti = header.iter().position(|h| *h == "tier").unwrap();
    let ri = header.iter().position(|h| *h == "reason").unwrap();
    for &p in &S2 {
        let row = cands
            .lines()
            .skip(1)
            .find(|l| l.starts_with(&format!("chr\t{p}\t")))
            .unwrap_or_else(|| panic!("S2 candidate {p} missing"));
        let f: Vec<&str> = row.split('\t').collect();
        assert_eq!((f[ti], f[ri]), ("none", "unlinked"), "candidate {p}: {f:?}");
    }
}

// ---------------------------------------------------------------------------
// Realistic CIGARs, supplementary consolidation and SA/clip evidence.

struct Rec {
    name: String,
    pos: usize,
    cigar: Vec<Cigar>,
    /// Reference positions replaced by the alt base.
    alts: Vec<usize>,
    rev: bool,
    supp: bool,
    sa: Option<&'static str>,
}

/// Build a BAM from explicit records over `refseq` (2 sites at 1000/1500). The
/// query is built by walking the CIGAR, so I/D/S/H ops are handled.
fn build_records(dir: &Path, refseq: &[u8], recs: &[Rec]) -> (PathBuf, PathBuf) {
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr\n{}\n", std::str::from_utf8(refseq).unwrap())).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), format!("chr\t{}\t5\t{}\t{}\n", refseq.len(), refseq.len(), refseq.len() + 1))
        .unwrap();
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", refseq.len()));
    let bam = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();
    let mut sorted: Vec<&Rec> = recs.iter().collect();
    sorted.sort_by_key(|r| r.pos);
    for r in sorted {
        let mut seq: Vec<u8> = Vec::new();
        let mut rr = r.pos;
        for c in &r.cigar {
            match *c {
                Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                    for k in 0..l as usize {
                        let p = rr + k;
                        seq.push(if r.alts.contains(&p) { alt_of(refseq[p]) } else { refseq[p] });
                    }
                    rr += l as usize;
                }
                Cigar::Del(l) | Cigar::RefSkip(l) => rr += l as usize,
                Cigar::Ins(l) | Cigar::SoftClip(l) => seq.extend(std::iter::repeat(b'A').take(l as usize)),
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        let qual = vec![30u8; seq.len()];
        let mut rec = Record::new();
        rec.set(r.name.as_bytes(), Some(&CigarString(r.cigar.clone())), &seq, &qual);
        rec.set_tid(0);
        rec.set_pos(r.pos as i64);
        rec.set_mapq(60);
        rec.set_flags((if r.rev { 0x10 } else { 0 }) | (if r.supp { 0x800 } else { 0 }));
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        if let Some(sa) = r.sa {
            rec.push_aux(b"SA", Aux::String(sa)).unwrap();
        }
        w.write(&rec).unwrap();
    }
    drop(w);
    bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
    (bam, fa)
}

#[test]
fn qc_linkage_consolidation_cigars_and_sa() {
    let dir = workdir("consolidation");
    let mut rng = Lcg(5);
    let refseq: Vec<u8> = (0..4_000).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect();
    let mut recs: Vec<Rec> = Vec::new();
    // 80 reference reads over both sites.
    for i in 0..80 {
        recs.push(Rec {
            name: format!("ref{i}"),
            pos: 600,
            cigar: vec![Cigar::Match(1_000)],
            alts: vec![],
            rev: i % 2 == 1,
            supp: false,
            sa: None,
        });
    }
    // 20 second-strain reads (alt at 1000 and 1500), with a 100 bp left clip
    // and an SA target at chr:3000.
    for i in 0..20 {
        recs.push(Rec {
            name: format!("strain{i}"),
            pos: 600,
            cigar: vec![Cigar::SoftClip(100), Cigar::Match(1_000)],
            alts: vec![1_000, 1_500],
            rev: i % 2 == 1,
            supp: false,
            sa: Some("chr,3000,+,1000M,60,0;"),
        });
    }
    // One read with a primary and 3 supplementary records, all alt at both
    // sites: consolidated to ONE read.
    for i in 0..4 {
        recs.push(Rec {
            name: "bad".to_string(),
            pos: 600,
            cigar: vec![Cigar::Match(1_000)],
            alts: vec![1_000, 1_500],
            rev: false,
            supp: i > 0,
            sa: None,
        });
    }
    // A split read: primary covers site 1000 (alt), supplementary covers site
    // 1500 (alt) -> consolidated to alt at both, counted once.
    recs.push(Rec {
        name: "split".to_string(),
        pos: 700,
        cigar: vec![Cigar::Match(500)],
        alts: vec![1_000],
        rev: false,
        supp: false,
        sa: None,
    });
    recs.push(Rec {
        name: "split".to_string(),
        pos: 1_300,
        cigar: vec![Cigar::Match(500)],
        alts: vec![1_500],
        rev: false,
        supp: true,
        sa: None,
    });
    // A read whose deletion covers site 1000 -> state `other`.
    recs.push(Rec {
        name: "delread".to_string(),
        pos: 800,
        cigar: vec![Cigar::Match(200), Cigar::Del(100), Cigar::Match(700)],
        alts: vec![],
        rev: false,
        supp: false,
        sa: None,
    });
    // Two records of one read disagreeing at both sites -> consolidated `other`.
    recs.push(Rec {
        name: "conflict".to_string(),
        pos: 600,
        cigar: vec![Cigar::Match(1_000)],
        alts: vec![1_000, 1_500],
        rev: false,
        supp: false,
        sa: None,
    });
    recs.push(Rec {
        name: "conflict".to_string(),
        pos: 600,
        cigar: vec![Cigar::Match(1_000)],
        alts: vec![],
        rev: false,
        supp: true,
        sa: None,
    });
    let (bam, fa) = build_records(&dir, &refseq, &recs);
    let prefix = dir.join("out");
    run_qc(&bam, &fa, &prefix, &params(1, None)).unwrap();

    let bed = std::fs::read_to_string(prefix.with_extension("bed")).unwrap();
    let lh: Vec<Vec<&str>> =
        bed.lines().filter(|l| l.contains("\tlinked_haplotype\t")).map(|l| l.split('\t').collect()).collect();
    assert_eq!(lh.len(), 1, "{bed}");
    let r = &lh[0];
    assert_eq!((r[1], r[2], r[5]), ("1000", "1501", "orange"), "{r:?}");
    // 20 strain + bad (once) + split (once) = 22 reads.
    assert!(r[9].contains("reads=22;"), "{r:?}");
    // SA evidence from the 10 strain reads, clip evidence from their left clips.
    assert!(r[9].contains("elsewhere=chr:2999;elsewhere_reads=20"), "{r:?}");
    assert!(r[9].contains("clipped=left:20"), "{r:?}");

    // haplotypes.tsv: `bad` and `split` once each, 22 rows.
    let hap = std::fs::read_to_string(dir.join("out.haplotypes.tsv")).unwrap();
    let names: Vec<&str> = hap.lines().skip(1).map(|l| l.split('\t').nth(3).unwrap()).collect();
    assert_eq!(names.len(), 22, "{hap}");
    assert_eq!(names.iter().filter(|n| **n == "bad").count(), 1, "{hap}");
    assert_eq!(names.iter().filter(|n| **n == "split").count(), 1, "{hap}");
    assert!(!names.contains(&"delread"), "deletion read is not a group read: {hap}");
    assert!(!names.contains(&"conflict"), "disagreeing read is not a group read: {hap}");
}

/// CIGAR operations (`I`, `H`, `N`, `=`, `X`, `D`) and the exact carried-allele
/// list: a read with a deletion or an `N` skip at one group site and ALT at two
/// others stays a group read, and the deleted/skipped position is not listed.
#[test]
fn qc_linkage_cigar_ops_and_carried_alleles() {
    let dir = workdir("cigar_ops");
    let mut rng = Lcg(5);
    let refseq: Vec<u8> = (0..4_000).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect();
    let sites = [1_000usize, 1_500, 2_000];
    let mut recs: Vec<Rec> = Vec::new();
    for i in 0..80 {
        recs.push(Rec {
            name: format!("ref{i}"),
            pos: 900,
            cigar: vec![Cigar::Match(1_200)],
            alts: vec![],
            rev: i % 2 == 1,
            supp: false,
            sa: None,
        });
    }
    for i in 0..20 {
        recs.push(Rec {
            name: format!("strain{i}"),
            pos: 900,
            cigar: vec![Cigar::Match(1_200)],
            alts: sites.to_vec(),
            rev: i % 2 == 1,
            supp: false,
            sa: None,
        });
    }
    // Insertion near site 1500; ALT at all three sites.
    recs.push(Rec {
        name: "insread".to_string(),
        pos: 900,
        cigar: vec![Cigar::Match(500), Cigar::Ins(50), Cigar::Match(650)],
        alts: sites.to_vec(),
        rev: false,
        supp: false,
        sa: None,
    });
    // Leading hard clip; ALT at all three sites.
    recs.push(Rec {
        name: "hardread".to_string(),
        pos: 950,
        cigar: vec![Cigar::HardClip(50), Cigar::Match(1_150)],
        alts: sites.to_vec(),
        rev: false,
        supp: false,
        sa: None,
    });
    // `=` and `X` blocks; ALT at all three sites.
    recs.push(Rec {
        name: "eqxread".to_string(),
        pos: 900,
        cigar: vec![Cigar::Equal(500), Cigar::Diff(500), Cigar::Equal(200)],
        alts: sites.to_vec(),
        rev: false,
        supp: false,
        sa: None,
    });
    // `N` skip covers site 1500; ALT at 1000 and 2000 only.
    recs.push(Rec {
        name: "skipread".to_string(),
        pos: 900,
        cigar: vec![Cigar::Match(600), Cigar::RefSkip(100), Cigar::Match(500)],
        alts: vec![1_000, 2_000],
        rev: false,
        supp: false,
        sa: None,
    });
    // Deletion covers site 1500; ALT at 1000 and 2000 only.
    recs.push(Rec {
        name: "delread".to_string(),
        pos: 900,
        cigar: vec![Cigar::Match(600), Cigar::Del(100), Cigar::Match(500)],
        alts: vec![1_000, 2_000],
        rev: false,
        supp: false,
        sa: None,
    });
    let (bam, fa) = build_records(&dir, &refseq, &recs);
    let prefix = dir.join("out");
    run_qc(&bam, &fa, &prefix, &params(1, None)).unwrap();

    // One orange group over the three sites.
    let bed = std::fs::read_to_string(prefix.with_extension("bed")).unwrap();
    let lh: Vec<Vec<&str>> =
        bed.lines().filter(|l| l.contains("\tlinked_haplotype\t")).map(|l| l.split('\t').collect()).collect();
    assert_eq!(lh.len(), 1, "{bed}");
    assert_eq!((lh[0][1], lh[0][2], lh[0][5]), ("1000", "2001", "orange"), "{:?}", lh[0]);
    assert!(lh[0][9].contains("reads=25;"), "{:?}", lh[0]);

    // Exact carried alleles from haplotypes.tsv.
    let hap = std::fs::read_to_string(dir.join("out.haplotypes.tsv")).unwrap();
    let alleles: std::collections::HashMap<&str, &str> =
        hap.lines().skip(1).map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[3], f[5])
        }).collect();
    let want = |name: &str, a: &str| {
        assert_eq!(alleles.get(name).copied(), Some(a), "{name} in {hap}");
    };
    // A deletion / skip at one site, ALT at the other two: still a group read,
    // and the deleted/skipped position is not carried.
    want("delread", "1000,2000");
    want("skipread", "1000,2000");
    // I / H / = / X do not change the carried alleles.
    want("insread", "1000,1500,2000");
    want("hardread", "1000,1500,2000");
    want("eqxread", "1000,1500,2000");
    want("strain0", "1000,1500,2000");
    // The deletion/skip reads really are group reads.
    assert!(alleles.contains_key("delread") && alleles.contains_key("skipread"), "{hap}");
    assert_eq!(alleles.len(), 25, "{hap}");
}

