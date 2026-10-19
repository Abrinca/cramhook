//! Integration test for `cramhook qc`: a synthetic reference and reads with
//! planted events (majority SNP, 50 % SNP, minority SNP, homopolymer deletion
//! placed differently by different reads, strand-biased alt, 80 bp deletion)
//! plus 0.5 % random substitution noise. Checks VCF / BED / JSON and that the
//! output is identical for 1 and 8 threads.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{self, Format, Header, Record, Writer};

use cramhook::{Preset, QcParams, run_qc};

const L: usize = 20_000;
const READ_LEN: usize = 1_000;
const STEP: usize = 25;
const SNP_RED: usize = 3_000;
const SNP_ORANGE: usize = 6_000;
const SNP_YELLOW: usize = 9_000;
const HP_START: usize = 12_000;
const HP_LEN: usize = 6;
const SNP_STRAND: usize = 15_000;
/// Minority alt (AF 0.25) on forward reads only: gated out of the VCF.
const SNP_STRAND_MINOR: usize = 16_000;
const LDEL_START: usize = 17_000;
const LDEL_LEN: usize = 80;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

fn workdir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/qc-scan/test").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn reference() -> Vec<u8> {
    let mut rng = Lcg(7);
    let mut r: Vec<u8> = (0..L).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect();
    // Break accidental homopolymers around the planted run.
    r[HP_START - 2] = b'G';
    r[HP_START - 1] = b'C';
    for x in &mut r[HP_START..HP_START + HP_LEN] {
        *x = b'A';
    }
    r[HP_START + HP_LEN] = b'G';
    // The large deletion must not left-shift (keeps the expected anchor simple).
    if r[LDEL_START - 1] == r[LDEL_START + LDEL_LEN - 1] {
        r[LDEL_START - 1] = alt_of(r[LDEL_START - 1]);
    }
    r
}

fn alt_of(b: u8) -> u8 {
    match b {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

fn build_inputs(dir: &Path) -> (PathBuf, PathBuf) {
    let refseq = reference();
    let fa = dir.join("ref.fa");
    let mut text = String::from(">chr\n");
    for c in refseq.chunks(60) {
        text.push_str(std::str::from_utf8(c).unwrap());
        text.push('\n');
    }
    std::fs::write(&fa, text).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), format!("chr\t{L}\t5\t60\t61\n")).unwrap();

    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", L));
    let bam_path = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
    let mut rng = Lcg(42);
    let planted = [SNP_RED, SNP_ORANGE, SNP_YELLOW, SNP_STRAND, SNP_STRAND_MINOR];
    let mut i = 0usize;
    let mut s = 0usize;
    while s + READ_LEN <= L {
        let e = s + READ_LEN;
        let rev = i % 2 == 1;
        let pair = i / 2;
        let has_hp_del = pair % 5 != 0;
        let hp_off = pair % HP_LEN;
        let has_ldel = pair % 5 != 0 && s < LDEL_START && e > LDEL_START + LDEL_LEN;
        let mut seq = Vec::new();
        let mut ops: Vec<Cigar> = Vec::new();
        let mut m = 0u32;
        let mut p = s;
        while p < e {
            let del_len = if has_hp_del && p == HP_START + hp_off && s < HP_START && e > HP_START + HP_LEN {
                1
            } else if has_ldel && p == LDEL_START {
                LDEL_LEN
            } else {
                0
            };
            if del_len > 0 {
                if m > 0 {
                    ops.push(Cigar::Match(m));
                    m = 0;
                }
                ops.push(Cigar::Del(del_len as u32));
                p += del_len;
                continue;
            }
            let r = refseq[p];
            let b = match p {
                SNP_RED => alt_of(r),
                SNP_ORANGE if pair % 2 == 0 => alt_of(r),
                SNP_YELLOW if pair % 4 == 0 => alt_of(r),
                SNP_STRAND if !rev => alt_of(r),
                SNP_STRAND_MINOR if !rev && pair % 2 == 0 => alt_of(r),
                _ if planted.contains(&p) => r,
                _ if rng.next() % 200 == 0 => b"ACGT"[(rng.next() % 4) as usize],
                _ => r,
            };
            seq.push(b);
            m += 1;
            p += 1;
        }
        if m > 0 {
            ops.push(Cigar::Match(m));
        }
        let qual = vec![30u8; seq.len()];
        let mut rec = Record::new();
        let name = format!("r{i}");
        rec.set(name.as_bytes(), Some(&CigarString(ops)), &seq, &qual);
        rec.set_tid(0);
        rec.set_pos(s as i64);
        rec.set_mapq(60);
        rec.set_flags(if rev { 0x10 } else { 0 });
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        w.write(&rec).unwrap();
        i += 1;
        s += STEP;
    }
    drop(w);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
    (bam_path, fa)
}

fn params(threads: usize, chunk: Option<u64>) -> QcParams {
    QcParams { preset: Some(Preset::MapHifi), threads, candidates: true, chunk_size: chunk, ..Default::default() }
}

/// VCF body records keyed by POS.
fn vcf_records(path: &Path) -> BTreeMap<u64, Vec<String>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out = BTreeMap::new();
    for l in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<String> = l.split('\t').map(|s| s.to_string()).collect();
        out.insert(f[1].parse().unwrap(), f);
    }
    out
}

fn json_without_volatile(path: &Path) -> serde_json::Value {
    let mut v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let o = v.as_object_mut().unwrap();
    o.remove("timing_s");
    o.remove("threads");
    o.remove("chunks");
    o.remove("inputs");
    // The linkage block's own timing is wall-clock; everything else in it is
    // thread-independent.
    if let Some(l) = o.get_mut("linkage").and_then(|l| l.as_object_mut()) {
        l.remove("timing_s");
    }
    v
}

#[test]
fn qc_planted_events_and_thread_identity() {
    let dir = workdir("planted");
    let (bam_path, fa) = build_inputs(&dir);
    let p1 = dir.join("t1/out");
    let p8 = dir.join("t8/out");
    // `--unlinked-action flag`: the lone 25 % minority SNP has only background-noise
    // partners nearby, so under the default `drop` it leaves the VCF (checked below).
    let flag = |threads, chunk| QcParams { unlinked_action: cramhook::qc::UnlinkedAction::Flag, ..params(threads, chunk) };
    run_qc(&bam_path, &fa, &p1, &flag(1, None)).unwrap();
    // Small chunks so records really are spread over many chunks and threads.
    let s8 = run_qc(&bam_path, &fa, &p8, &flag(8, Some(1_500))).unwrap();
    assert!(s8.chunks > 10);
    // Default `drop`: the lone minority SNP is tested against noise partners, linked to
    // none, and dropped (tier none, reason `unlinked`), but stays in candidates.tsv.
    let pd = dir.join("drop/out");
    run_qc(&bam_path, &fa, &pd, &params(1, None)).unwrap();
    assert!(!vcf_records(&pd.with_extension("vcf")).contains_key(&((SNP_YELLOW + 1) as u64)));
    let cand = std::fs::read_to_string(pd.with_extension("candidates.tsv")).unwrap();
    assert!(
        cand.lines().any(|l| l.split('\t').nth(1) == Some(&SNP_YELLOW.to_string()[..]) && l.contains("\tunlinked")),
        "lone minority SNP not marked unlinked in candidates.tsv"
    );

    let vcf = vcf_records(&p1.with_extension("vcf"));
    let get = |pos: usize| vcf.get(&(pos as u64)).unwrap_or_else(|| panic!("no VCF record at POS {pos}: {vcf:?}"));
    // VCF POS is 1-based; indels use the anchor base.
    // FILTER is confidence only (PASS here: full depth); the tier is INFO TIER.
    let red = get(SNP_RED + 1);
    assert_eq!(red[6], "PASS", "{red:?}");
    assert!(red[7].starts_with("TIER=RED;TYPE=SNP;AF=1.0000"));
    let orange = get(SNP_ORANGE + 1);
    assert_eq!(orange[6], "PASS", "{orange:?}");
    assert!(orange[7].contains("TIER=ORANGE") && orange[7].contains(";MIXED"));
    let yellow = get(SNP_YELLOW + 1);
    assert_eq!(yellow[6], "PASS", "{yellow:?}");
    assert!(yellow[7].contains("TIER=YELLOW") && yellow[7].contains("REASON=minority") && !yellow[7].contains("STRAND"));
    // Homopolymer deletion: every read deletes a different A, all normalised to one allele.
    let hp = get(HP_START);
    assert_eq!((hp[3].as_str(), hp[4].as_str(), hp[6].as_str()), ("CA", "C", "PASS"), "{hp:?}");
    assert!(hp[7].contains("TIER=RED") && hp[7].contains("TYPE=DEL") && hp[7].contains("HP=6") && hp[7].contains("CLASS=DEL_HP6"));
    // Strand-biased: flagged, never orange.
    let sb = get(SNP_STRAND + 1);
    assert_eq!(sb[6], "PASS", "{sb:?}");
    assert!(sb[7].contains("TIER=YELLOW") && sb[7].contains(";STRAND") && sb[7].contains("REASON=strand_bias"));
    // Strand-biased minority: not called, but kept as a candidate.
    assert!(!vcf.contains_key(&((SNP_STRAND_MINOR + 1) as u64)));
    let tsv = std::fs::read_to_string(dir.join("t1/out.candidates.tsv")).unwrap();
    let row: Vec<&str> = tsv
        .lines()
        .find(|l| l.starts_with(&format!("chr\t{SNP_STRAND_MINOR}\tSNP\t")))
        .unwrap_or_else(|| panic!("no candidate at {SNP_STRAND_MINOR}"))
        .split('\t')
        .collect();
    assert_eq!((row[14], row[23], row[24], row[25], row[26]), ("0.2500", "1", "none", "strand_bias", "STRAND"), "{row:?}");
    // The large deletion is also a site.
    let ld = get(LDEL_START);
    assert_eq!(ld[6], "PASS", "{ld:?}");
    assert!(ld[7].contains("TIER=RED"));
    assert_eq!(ld[3].len(), LDEL_LEN + 1);
    // Nothing else is called (noise only).
    assert_eq!(vcf.len(), 6, "{:?}", vcf.keys().collect::<Vec<_>>());

    // BED: large deletion region with size and support.
    let bed = std::fs::read_to_string(p1.with_extension("bed")).unwrap();
    let ldel: Vec<Vec<&str>> =
        bed.lines().filter(|l| l.contains("\tlarge_del\t")).map(|l| l.split('\t').collect()).collect();
    assert_eq!(ldel.len(), 1, "{bed}");
    let r = &ldel[0];
    assert_eq!((r[1], r[2], r[5]), ("17000", "17080", "red"));
    assert!(r[9].contains(&format!("size={LDEL_LEN};")), "{r:?}");
    assert!(!bed.contains("snp_cluster"));
    // Contig ends are not assessable (low depth), never clean.
    assert!(bed.lines().any(|l| l.starts_with("chr\t0\t") && l.contains("not_assessable\t0\tnone\tlow_depth")));

    // JSON.
    let json_text = std::fs::read_to_string(p1.with_extension("json")).unwrap();
    assert!(json_text.starts_with("{\n  \"format_version\": 1,"), "format_version is not the first JSON key:\n{json_text}");
    let j: serde_json::Value = serde_json::from_str(&json_text).unwrap();
    assert_eq!(j["format_version"], 1);
    let vcf_text = std::fs::read_to_string(p1.with_extension("vcf")).unwrap();
    assert!(vcf_text.starts_with("##fileformat=VCFv4.2\n##format_version=1\n"), "VCF format_version header missing:\n{vcf_text}");
    // Alerts: VCF site calls by type plus tiered BED regions.
    let bt = &j["alerts"]["by_type"];
    assert_eq!(bt["SNP"]["red"], 1);
    assert_eq!(bt["DEL"]["red"], 2);
    assert_eq!(bt["INS"]["red"], 0);
    assert_eq!(bt["large_del"]["red"], 1);
    assert_eq!(j["alerts"]["red"], 4); // 3 red sites + the large_del region
    assert_eq!(bt["SNP"]["orange"], 1);
    assert_eq!(j["alerts"]["orange"], 1);
    assert_eq!(bt["SNP"]["yellow"], 2);
    assert_eq!(j["alerts"]["yellow"], 2);
    assert_eq!(j["reads"]["records"]["primary"], ((L - READ_LEN) / STEP + 1) as u64);
    let chr = &j["contigs"]["chr"];
    assert_eq!(chr["depth_mode"], 40);
    assert_eq!(chr["depth_mode_excludes_zero"], true);
    assert_eq!(chr["covered_fraction"], 1.0);
    // qc_qv is computed over confident (not low by abs) assessable bases.
    let confident = j["assembly"]["confident_bp"].as_f64().unwrap();
    let want_qv = -10.0 * (5.0 / (confident + 1.0)).log10();
    assert!((j["qc_qv"].as_f64().unwrap() - want_qv).abs() < 1e-9);
    let snp_fit = &j["background"][0];
    assert_eq!(snp_fit["class"], "SNP");
    assert!(snp_fit["p"].as_f64().unwrap() > 1e-4 && snp_fit["p"].as_f64().unwrap() < 0.02, "{snp_fit}");

    // Thread / chunking identity.
    for suf in ["vcf", "bed", "candidates.tsv", "depth.bedgraph", "gc.bedgraph"] {
        let a = std::fs::read(dir.join(format!("t1/out.{suf}"))).unwrap();
        let b = std::fs::read(dir.join(format!("t8/out.{suf}"))).unwrap();
        assert!(a == b, "{suf} differs between --threads 1 and 8");
    }
    assert_eq!(json_without_volatile(&p1.with_extension("json")), json_without_volatile(&p8.with_extension("json")));
    assert_eq!(j["params"]["yellow_strand_gate"], true);

    // --no-yellow-strand-gate restores the old rule: yellow with STRAND.
    let p0 = dir.join("nogate/out");
    run_qc(&bam_path, &fa, &p0, &QcParams { yellow_strand_gate: false, ..flag(1, None) }).unwrap();
    let vcf0 = vcf_records(&p0.with_extension("vcf"));
    let r0 = &vcf0[&((SNP_STRAND_MINOR + 1) as u64)];
    assert_eq!(r0[6], "PASS", "{r0:?}");
    assert!(r0[7].contains("TIER=YELLOW") && r0[7].contains(";STRAND") && r0[7].contains("REASON=minority"));
    assert_eq!(vcf0.len(), 7);
}

/// Low coverage: a reference with a normal-depth (40×) region and a low-depth
/// (≈12×) region, each with a planted majority SNP. With `lr:hq` (low_abs 20)
/// both SNPs are TIER=RED; the low one gets FILTER=LowCov (LOWCOV flag) and is
/// counted in `alerts.lowcov`, not `alerts`. The BED gets a `low_coverage`
/// annotation region and the JSON a `low_coverage` warning on that contig;
/// `--threads 1` and 8 give identical output.
#[test]
fn qc_low_coverage_region_and_filter() {
    // `run_case_p` cleans and writes its own workdirs (`lowcov_t1`/`lowcov_t8`).
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/qc-scan/test");
    let l = 8_000usize;
    let r = random_ref(l, 71);
    let snp_hi = 2_500usize;
    let snp_lo = 6_300usize;
    let mut reads: Vec<Rd> = Vec::new();
    let mut i = 0usize;
    // Normal region: 40× over 0..5000 (1 kb reads starting every 25 bp).
    for s in (0..=4_000).step_by(25) {
        let mut seq = r[s..s + 1_000].to_vec();
        if (s..s + 1_000).contains(&snp_hi) {
            seq[snp_hi - s] = alt_of(r[snp_hi]);
        }
        reads.push(Rd { pos: s, cigar: vec![Cigar::Match(1_000)], seq, rev: i % 2 == 1 });
        i += 1;
    }
    // Low region: ≈12× over 5000..8000 (reads starting every 80 bp).
    for s in (5_000..=7_000).step_by(80) {
        let mut seq = r[s..s + 1_000].to_vec();
        if (s..s + 1_000).contains(&snp_lo) {
            seq[snp_lo - s] = alt_of(r[snp_lo]);
        }
        reads.push(Rd { pos: s, cigar: vec![Cigar::Match(1_000)], seq, rev: i % 2 == 1 });
        i += 1;
    }
    let lr = |threads, chunk| QcParams {
        preset: Some(Preset::LrHq),
        threads,
        candidates: true,
        chunk_size: chunk,
        ..Default::default()
    };
    let (vcf, cands, bed, json) = run_case_p("lowcov_t1", &r, reads.clone(), lr(1, None));
    let (vcf8, _c8, bed8, _j8) = run_case_p("lowcov_t8", &r, reads.clone(), lr(8, Some(1_500)));
    assert_af_le_one(&cands);

    // Normal-depth majority SNP: TIER RED, FILTER PASS, no LOWCOV.
    let hi = vcf.iter().find(|v| v[1] == (snp_hi + 1).to_string()).unwrap_or_else(|| panic!("{vcf:?}"));
    assert_eq!(hi[6], "PASS", "{hi:?}");
    assert!(hi[7].contains("TIER=RED") && !hi[7].contains("LOWCOV"), "{hi:?}");
    // Low-depth majority SNP: TIER RED (uncapped), FILTER LowCov + LOWCOV flag.
    let lo = vcf.iter().find(|v| v[1] == (snp_lo + 1).to_string()).unwrap_or_else(|| panic!("{vcf:?}"));
    assert_eq!(lo[6], "LowCov", "{lo:?}");
    assert!(lo[7].contains("TIER=RED") && lo[7].contains(";LOWCOV"), "{lo:?}");
    let lo_cand = cands
        .iter()
        .find(|c| c["pos0"] == snp_lo.to_string())
        .unwrap_or_else(|| panic!("{cands:?}"));
    assert_eq!((lo_cand["tier"].as_str(), lo_cand["lowcov"].as_str()), ("red", "abs"));
    assert!(!cands.iter().any(|c| c.contains_key("tier_uncapped")), "tier_uncapped column removed");

    // BED: a low_coverage annotation region (untiered).
    let low = bed.lines().filter(|l| l.contains("\tlow_coverage\t")).collect::<Vec<_>>();
    assert!(!low.is_empty(), "{bed}");
    assert!(low.iter().all(|l| l.split('\t').nth(5) == Some("none")), "{bed}");

    // Alerts: the low-depth majority SNP is counted under `lowcov`, not `red`.
    // The 40× SNP stays a PASS red call, so SNP red = 1 and SNP lowcov red = 1.
    let alerts = &json["alerts"];
    let red_pass = vcf.iter().filter(|v| v[7].contains("TIER=RED") && v[6] == "PASS").count() as u64;
    let red_low = vcf.iter().filter(|v| v[7].contains("TIER=RED") && v[6] == "LowCov").count() as u64;
    assert!(red_low >= 1, "{vcf:?}");
    // Top-level red also counts the tiered depth-band BED regions.
    let depth_red = ["depth_low", "depth_half", "depth_high"]
        .iter()
        .map(|l| alerts["by_type"][l]["red"].as_u64().unwrap())
        .sum::<u64>();
    assert_eq!(alerts["red"].as_u64().unwrap(), red_pass + depth_red, "{alerts}");
    assert_eq!(alerts["lowcov"]["red"].as_u64().unwrap(), red_low, "{alerts}");
    assert_eq!(alerts["by_type"]["SNP"]["lowcov"]["red"], 1, "{alerts}");

    // JSON: per-contig and top-level low_coverage warnings with sensible bp;
    // the wording must not claim a cap.
    let cw = json["contigs"]["chr"]["warnings"].as_array().unwrap();
    let w = cw.iter().find(|w| w["code"] == "low_coverage").unwrap_or_else(|| panic!("{json}"));
    // Reason is the union over the region; the ramps also trigger `rel`.
    assert!(w["reason"].as_str().unwrap().contains("abs"), "{w}");
    assert!(w["message"].as_str().unwrap().contains("calls marked LowCov"), "{w}");
    assert!(!w["message"].as_str().unwrap().contains("capped"), "{w}");
    assert!(w["bp"].as_u64().unwrap() >= 1_000, "{w}");
    assert!(w["fraction"].as_f64().unwrap() >= 0.05, "{w}");
    let tw = json["warnings"].as_array().unwrap();
    let tw = tw.iter().find(|w| w["code"] == "low_coverage").unwrap_or_else(|| panic!("{json}"));
    assert_eq!(tw["contigs"], serde_json::json!(["chr"]), "{tw}");
    assert_eq!(tw["bp"], w["bp"], "{tw} vs {w}");

    // Findings BED vs JSON: tier and PASS/lowcov counts must agree.
    let fbed = std::fs::read_to_string(base.join("lowcov_t1/out.findings.bed")).unwrap();
    let mut pass = [0u64; 3]; // red, orange, yellow
    let mut low = [0u64; 3];
    for l in fbed.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = l.split('\t').collect();
        let i = match f[5] {
            "red" => 0,
            "orange" => 1,
            "yellow" => 2,
            _ => continue,
        };
        if f[7].contains("lowcov=1") {
            low[i] += 1
        } else {
            pass[i] += 1
        }
    }
    let jf = &json["findings"];
    assert_eq!(
        (pass[0], pass[1], pass[2]),
        (jf["red"].as_u64().unwrap(), jf["orange"].as_u64().unwrap(), jf["yellow"].as_u64().unwrap()),
        "findings BED PASS counts vs JSON:\n{fbed}"
    );
    assert_eq!(
        (low[0], low[1], low[2]),
        (
            jf["lowcov"]["red"].as_u64().unwrap(),
            jf["lowcov"]["orange"].as_u64().unwrap(),
            jf["lowcov"]["yellow"].as_u64().unwrap()
        ),
        "findings BED lowcov counts vs JSON:\n{fbed}"
    );
    assert!(pass.iter().sum::<u64>() >= 1 && low.iter().sum::<u64>() >= 1, "{fbed}");

    // Thread / chunking identity (the two runs use different workdirs, so the
    // VCF/BED headers embed different reference paths; compare the bodies).
    assert_eq!(vcf, vcf8, "VCF body differs between --threads 1 and 8");
    assert_eq!(bed, bed8, "BED differs between --threads 1 and 8");
    for suf in ["candidates.tsv", "depth.bedgraph", "gc.bedgraph", "findings.bed"] {
        let a = std::fs::read(base.join(format!("lowcov_t1/out.{suf}"))).unwrap();
        let b = std::fs::read(base.join(format!("lowcov_t8/out.{suf}"))).unwrap();
        assert!(a == b, "{suf} differs between --threads 1 and 8");
    }
    assert_eq!(
        json_without_volatile(&base.join("lowcov_t1/out.json")),
        json_without_volatile(&base.join("lowcov_t8/out.json"))
    );
}

/// A strand-balanced majority SNP whose q-value is not significant must still
/// be reported (yellow `weak_majority`), never dropped: the majority reporting
/// gate is independent of the q-value (step-3 regression). FILTER follows the
/// depth: 6x < `--lowcov-abs` (10) gives `LowCov`.
#[test]
fn qc_weak_majority_reported_without_significance() {
    let l = 3_000usize;
    let r = random_ref(l, 909);
    let snp = 1_500usize;
    // 6 reads at 6x, 4 carrying the alt (AF 0.667, Wilson LB < 0.5); strand
    // balanced (3 fwd / 3 rev, alt on both strands).
    let mut reads = Vec::new();
    for i in 0..6usize {
        let mut seq = r[1_000..2_000].to_vec();
        if i < 4 {
            seq[snp - 1_000] = alt_of(r[snp]);
        }
        reads.push(Rd { pos: 1_000, cigar: vec![Cigar::Match(1_000)], seq, rev: i % 2 == 1 });
    }
    // `--fdr 0`: no site can be significant, so the old ordering dropped it.
    let qp = QcParams { preset: Some(Preset::MapHifi), fdr: 0.0, candidates: true, ..Default::default() };
    let (vcf, cands, _bed, _json) = run_case_p("weak_majority_report", &r, reads, qp);
    let v = vcf.iter().find(|v| v[1] == (snp + 1).to_string()).unwrap_or_else(|| panic!("{vcf:?}"));
    assert_eq!(v[6], "LowCov", "{v:?}");
    assert!(v[7].contains("TIER=YELLOW") && v[7].contains("REASON=weak_majority"), "{v:?}");
    assert!(v[7].contains(";LOWCOV"), "{v:?}");
    let c = cands.iter().find(|c| c["pos0"] == snp.to_string()).unwrap_or_else(|| panic!("{cands:?}"));
    assert_eq!(c["tier"], "yellow", "{c:?}");
    assert_eq!(c["reason"], "weak_majority", "{c:?}");
    assert_eq!(c["significant"], "0", "must exercise the non-significant path: {c:?}");
}

/// The JSON `contigs` object follows FASTA order, not the BAM `@SQ` order.
#[test]
fn qc_contigs_json_in_fasta_order() {
    let dir = workdir("fasta_order");
    let a = random_ref(300, 101);
    let b = random_ref(300, 202);
    // FASTA order: chrB first, then chrA.
    let fa = dir.join("ref.fa");
    let mut text = String::new();
    for (name, seq) in [("chrB", &b), ("chrA", &a)] {
        text.push('>');
        text.push_str(name);
        text.push('\n');
        text.push_str(std::str::from_utf8(seq).unwrap());
        text.push('\n');
    }
    std::fs::write(&fa, text).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), "chrB\t300\t6\t300\t301\nchrA\t300\t314\t300\t301\n").unwrap();
    // BAM header order: chrA (tid 0) then chrB (tid 1) - opposite to the FASTA.
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chrA").push_tag(b"LN", 300));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chrB").push_tag(b"LN", 300));
    let bam = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();
    let cig = CigarString(vec![Cigar::Match(100)]);
    for (tid, (name, seq)) in [(0i32, ("chrA", &a)), (1i32, ("chrB", &b))] {
        for i in 0..6 {
            let mut rec = Record::new();
            rec.set(format!("{name}r{i}").as_bytes(), Some(&cig), &seq[100..200], &vec![30u8; 100]);
            rec.set_tid(tid);
            rec.set_pos(100);
            rec.set_mapq(60);
            rec.set_flags(if i % 2 == 1 { 0x10 } else { 0 });
            rec.set_mtid(-1);
            rec.set_mpos(-1);
            w.write(&rec).unwrap();
        }
    }
    drop(w);
    bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();
    let prefix = dir.join("out");
    run_qc(&bam, &fa, &prefix, &params(1, None)).unwrap();
    let text = std::fs::read_to_string(prefix.with_extension("json")).unwrap();
    // The `contigs` object keys are `"<name>": {` (the name also appears in
    // the warnings), so match the key form.
    let ib = text.find("\"chrB\": {").expect("chrB key");
    let ia = text.find("\"chrA\": {").expect("chrA key");
    assert!(ib < ia, "JSON contigs not in FASTA order (chrB before chrA):\n{text}");
}

/// A `rel`-only low-coverage warning: with `lowcov_abs` at the minimum depth
/// the abs bit can never fire, and a long shallow region triggers `rel`.
#[test]
fn qc_rel_only_warning_wording() {
    let l = 10_000usize;
    let r = random_ref(l, 77);
    let mut reads: Vec<Rd> = Vec::new();
    let mut i = 0usize;
    // Background ≈8× everywhere.
    for s in (0..=9_000).step_by(125) {
        reads.push(Rd { pos: s, cigar: vec![Cigar::Match(1_000)], seq: r[s..s + 1_000].to_vec(), rev: i % 2 == 1 });
        i += 1;
    }
    // Extra reads over the first 7 kb → ≈40× there (the depth mode).
    for s in (0..=6_000).step_by(30) {
        reads.push(Rd { pos: s, cigar: vec![Cigar::Match(1_000)], seq: r[s..s + 1_000].to_vec(), rev: i % 2 == 1 });
        i += 1;
    }
    let qp = QcParams { preset: Some(Preset::MapHifi), lowcov_abs: Some(5), threads: 1, candidates: true, ..Default::default() };
    let (_v, _c, _bed, json) = run_case_p("rel_only", &r, reads, qp);
    let cw = json["contigs"]["chr"]["warnings"].as_array().unwrap();
    let w = cw.iter().find(|w| w["code"] == "low_coverage").unwrap_or_else(|| panic!("{json}"));
    assert_eq!(w["reason"], "rel", "{w}");
    let msg = w["message"].as_str().unwrap();
    assert!(msg.contains("of contig depth"), "{w}");
    assert!(!msg.contains("capped"), "{w}");
    // rel-only records are FILTER=PASS (INFO LOWCOV only), so the wording must
    // not claim calls were marked LowCov; it must say they remain PASS.
    assert!(!msg.contains("calls marked LowCov"), "{w}");
    assert!(msg.contains("LOWCOV annotation"), "{w}");
    assert!(msg.contains("calls remain PASS"), "{w}");
}

// ---------------------------------------------------------------------------
// Small regression scenarios (review findings).

#[derive(Clone)]
struct Rd {
    pos: usize,
    cigar: Vec<Cigar>,
    seq: Vec<u8>,
    rev: bool,
}

fn random_ref(len: usize, seed: u64) -> Vec<u8> {
    let mut rng = Lcg(seed);
    (0..len).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect()
}

/// Read sequence for a CIGAR at `pos`: reference bases for M, `ins` chunks for I (in order).
fn read_seq(refseq: &[u8], pos: usize, cigar: &[Cigar], ins: &[&[u8]]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut r = pos;
    let mut k = 0;
    for c in cigar {
        match *c {
            Cigar::Match(l) => {
                out.extend_from_slice(&refseq[r..r + l as usize]);
                r += l as usize;
            }
            Cigar::Ins(_) => {
                out.extend_from_slice(ins[k]);
                k += 1;
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => r += l as usize,
            _ => panic!("unsupported op"),
        }
    }
    out
}

/// One planted homopolymer edit for [`hp_reads`].
#[derive(Clone, Copy, PartialEq, Eq)]
enum Edit {
    None,
    /// Delete one reference base at `p` (normalises into its homopolymer run).
    Del,
    /// Insert one `A` before `p` (normalises into its homopolymer run).
    InsA,
}

/// Reads tiling `r` (1 kb long, starting every 25 bp, 0.5 % substitution
/// noise). For each reference position `p` of a read, `edit(pair, rev, p)` says
/// whether to delete a base (`Del`) or insert an `A` (`InsA`). Deletions and
/// insertions normalise into their homopolymer run (run bases flanked by a
/// non-`A`), so reads that edit different run positions collapse to one allele.
/// Use run anchors that are not multiples of 25 so no read starts inside a run.
fn hp_reads(r: &[u8], edit: impl Fn(usize, bool, usize) -> Edit) -> Vec<Rd> {
    const RL: usize = 1_000;
    const STEP: usize = 25;
    let mut reads = Vec::new();
    let mut rng = Lcg(4242);
    let mut i = 0usize;
    let mut s = 0usize;
    while s + RL <= r.len() {
        let e = s + RL;
        let rev = i % 2 == 1;
        let pair = i / 2;
        let mut seq = Vec::new();
        let mut ops: Vec<Cigar> = Vec::new();
        let mut m = 0u32;
        let mut p = s;
        let mut ins_at: Option<usize> = None;
        while p < e {
            if ins_at != Some(p) {
                match edit(pair, rev, p) {
                    Edit::InsA => {
                        if m > 0 {
                            ops.push(Cigar::Match(m));
                            m = 0;
                        }
                        ops.push(Cigar::Ins(1));
                        seq.push(b'A');
                        ins_at = Some(p);
                        continue;
                    }
                    Edit::Del => {
                        if m > 0 {
                            ops.push(Cigar::Match(m));
                            m = 0;
                        }
                        ops.push(Cigar::Del(1));
                        p += 1;
                        continue;
                    }
                    Edit::None => {}
                }
            }
            let b = if rng.next() % 200 == 0 { b"ACGT"[(rng.next() % 4) as usize] } else { r[p] };
            seq.push(b);
            m += 1;
            p += 1;
        }
        if m > 0 {
            ops.push(Cigar::Match(m));
        }
        reads.push(Rd { pos: s, cigar: ops, seq, rev });
        i += 1;
        s += STEP;
    }
    reads
}

type CaseOut = (Vec<Vec<String>>, Vec<BTreeMap<String, String>>, String, serde_json::Value);

/// Write ref + BAM, run qc (map-hifi, 1 thread), return (VCF records, candidate rows, BED text).
fn run_case(name: &str, refseq: &[u8], reads: Vec<Rd>) -> (Vec<Vec<String>>, Vec<BTreeMap<String, String>>, String) {
    let (v, c, b, _) = run_case_p(name, refseq, reads, params(1, None));
    (v, c, b)
}

/// As `run_case` with explicit parameters; also returns the JSON summary.
fn run_case_p(name: &str, refseq: &[u8], mut reads: Vec<Rd>, qp: QcParams) -> CaseOut {
    let dir = workdir(name);
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr\n{}\n", std::str::from_utf8(refseq).unwrap())).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), format!("chr\t{}\t5\t{}\t{}\n", refseq.len(), refseq.len(), refseq.len() + 1))
        .unwrap();
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", refseq.len()));
    let bam_path = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
    reads.sort_by_key(|r| r.pos);
    for (i, r) in reads.iter().enumerate() {
        let mut rec = Record::new();
        let qual = vec![30u8; r.seq.len()];
        rec.set(format!("r{i}").as_bytes(), Some(&CigarString(r.cigar.clone())), &r.seq, &qual);
        rec.set_tid(0);
        rec.set_pos(r.pos as i64);
        rec.set_mapq(60);
        rec.set_flags(if r.rev { 0x10 } else { 0 });
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        w.write(&rec).unwrap();
    }
    drop(w);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
    let prefix = dir.join("out");
    run_qc(&bam_path, &fa, &prefix, &qp).unwrap();
    let vcf: Vec<Vec<String>> = std::fs::read_to_string(prefix.with_extension("vcf"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').map(|s| s.to_string()).collect())
        .collect();
    let tsv = std::fs::read_to_string(dir.join("out.candidates.tsv")).unwrap();
    let mut lines = tsv.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let cands = lines
        .map(|l| head.iter().map(|h| h.to_string()).zip(l.split('\t').map(|s| s.to_string())).collect())
        .collect();
    let bed = std::fs::read_to_string(prefix.with_extension("bed")).unwrap();
    let json = serde_json::from_str(&std::fs::read_to_string(prefix.with_extension("json")).unwrap()).unwrap();
    (vcf, cands, bed, json)
}

fn assert_af_le_one(cands: &[BTreeMap<String, String>]) {
    for c in cands {
        let af: f64 = c["af"].parse().unwrap();
        let alt: u32 = c["alt_reads"].parse().unwrap();
        let dp: u32 = c["dp"].parse().unwrap();
        assert!(af <= 1.0 && alt <= dp, "{c:?}");
    }
}

/// P1-1: a large insertion after the last base must not panic; the region
/// stays inside the contig (preceding anchor). Its mean DP is below `low_abs`,
/// so it is tagged `lowcov=1` and counted under `alerts.lowcov`; a second
/// large insertion at exactly the absolute threshold stays a main PASS call.
#[test]
fn qc_terminal_large_insertion() {
    let mut r = random_ref(300, 3);
    let ins = random_ref(50, 99);
    // No left shift: the last base differs from the inserted sequence's last base.
    r[299] = alt_of(ins[49]);
    let mut reads: Vec<Rd> = Vec::new();
    // Terminal insertion: 6 reads → region mean DP 6 < low_abs (map-hifi 10).
    let cig = vec![Cigar::Match(5), Cigar::Ins(50)];
    for i in 0..6 {
        reads.push(Rd { pos: 295, seq: read_seq(&r, 295, &cig, &[&ins]), cigar: cig.clone(), rev: i % 2 == 1 });
    }
    // Interior insertion farther than `large_ins_window` (100) so the two
    // events stay separate: 10 reads cover the anchor → mean DP exactly 10,
    // i.e. not below `low_abs`, so it is a PASS call with no lowcov tag.
    let cig2 = vec![Cigar::Match(5), Cigar::Ins(50), Cigar::Match(5)];
    for i in 0..10 {
        reads.push(Rd { pos: 100, seq: read_seq(&r, 100, &cig2, &[&ins]), cigar: cig2.clone(), rev: i % 2 == 1 });
    }
    let (vcf, cands, bed, json) = run_case_p("terminal_ins", &r, reads, params(1, None));

    let reg: Vec<Vec<&str>> =
        bed.lines().filter(|l| l.contains("\tlarge_ins\t")).map(|l| l.split('\t').collect()).collect();
    assert_eq!(reg.len(), 2, "{bed}");
    // Terminal region: mean DP 6 < low_abs → extra carries `lowcov=1`.
    let term = reg.iter().find(|r| r[1] == "299").unwrap_or_else(|| panic!("{bed}"));
    assert_eq!((term[1], term[2], term[5]), ("299", "300", "red"), "{bed}");
    assert!(term[9].contains("lowcov=1"), "{term:?}");
    // Exactly-at-threshold region: mean DP 10 is not < low_abs → stays PASS.
    let at = reg.iter().find(|r| r[1] == "105").unwrap_or_else(|| panic!("{bed}"));
    assert_eq!((at[1], at[2], at[5]), ("105", "106", "red"), "{bed}");
    assert!(!at[9].contains("lowcov=1"), "{at:?}");

    // Alerts: the terminal region is counted under `lowcov`, the threshold one
    // stays in the main PASS red counts.
    let bt = &json["alerts"]["by_type"]["large_ins"];
    assert_eq!(bt["lowcov"]["red"], 1, "{bt}");
    assert_eq!(bt["red"], 1, "{bt}");

    assert_af_le_one(&cands);
    assert!(vcf.iter().any(|v| v[1] == "300" && v[7].contains("TYPE=INS")), "{vcf:?}");
}

/// P1-2: insertions in a homopolymer from reads that start inside the run
/// left-normalise to before their first aligned base; they must not be
/// counted against unrelated depth at the new anchor.
#[test]
fn qc_normalised_before_read_start() {
    // C AAAA...A (1..60) G + random.
    let mut r = random_ref(200, 5);
    r[0] = b'C';
    for x in &mut r[1..60] {
        *x = b'A';
    }
    r[60] = b'G';
    let mut reads = Vec::new();
    // 6 reference reads over 0..10.
    for i in 0..6 {
        let cig = vec![Cigar::Match(10)];
        reads.push(Rd { pos: 0, seq: read_seq(&r, 0, &cig, &[]), cigar: cig, rev: i % 2 == 1 });
    }
    // 12 reads over 50..60 with an extra A inside the run (normalises to anchor 0).
    for i in 0..12 {
        let cig = vec![Cigar::Match(5), Cigar::Ins(1), Cigar::Match(5)];
        reads.push(Rd { pos: 50, seq: read_seq(&r, 50, &cig, &[b"A"]), cigar: cig, rev: i % 2 == 1 });
    }
    // 6 reads spanning the whole run with the insertion: these count.
    for i in 0..6 {
        let cig = vec![Cigar::Match(30), Cigar::Ins(1), Cigar::Match(50)];
        reads.push(Rd { pos: 0, seq: read_seq(&r, 0, &cig, &[b"A"]), cigar: cig, rev: i % 2 == 1 });
    }
    let (_vcf, cands, _bed) = run_case("norm_before_start", &r, reads);
    assert_af_le_one(&cands);
    let ins: Vec<_> = cands.iter().filter(|c| c["type"] == "INS").collect();
    assert_eq!(ins.len(), 1, "{ins:?}");
    assert_eq!((ins[0]["pos0"].as_str(), ins[0]["alt_reads"].as_str(), ins[0]["dp"].as_str()), ("0", "6", "12"));
}

/// P1-3: two 1-bp insertions of one read in the same homopolymer are one
/// 2-bp event, supported once.
#[test]
fn qc_multiple_indels_one_read_count_once() {
    let mut r = random_ref(200, 11);
    r[0] = b'C';
    for x in &mut r[1..21] {
        *x = b'A';
    }
    r[21] = b'G';
    let cig = vec![Cigar::Match(5), Cigar::Ins(1), Cigar::Match(5), Cigar::Ins(1), Cigar::Match(90)];
    let reads = (0..6)
        .map(|i| Rd { pos: 0, seq: read_seq(&r, 0, &cig, &[b"A", b"A"]), cigar: cig.clone(), rev: i % 2 == 1 })
        .collect();
    let (vcf, cands, _bed) = run_case("multi_indel", &r, reads);
    assert_af_le_one(&cands);
    let ins: Vec<_> = cands.iter().filter(|c| c["type"] == "INS").collect();
    assert_eq!(ins.len(), 1, "{ins:?}");
    assert_eq!((ins[0]["len"].as_str(), ins[0]["alt"].as_str(), ins[0]["alt_reads"].as_str()), ("2", "CAA", "6"));
    // dp = 6 < lowcov_abs (map-hifi 10): TIER stays RED, FILTER becomes LowCov.
    assert!(
        vcf.iter().any(|v| v[1] == "1"
            && v[3] == "C"
            && v[4] == "CAA"
            && v[6] == "LowCov"
            && v[7].contains("TIER=RED")
            && v[7].contains(";LOWCOV")),
        "{vcf:?}"
    );
}

/// P2-4: long insertions keep their sequence identity; single-read alleles
/// of equal length are not pooled into another sequence's AD.
#[test]
fn qc_long_insertions_not_pooled() {
    let mut r = random_ref(300, 13);
    let cig = vec![Cigar::Match(100), Cigar::Ins(40), Cigar::Match(100)];
    let shared = random_ref(40, 1000);
    // The insertion (before ref 150) must not left-shift.
    r[149] = alt_of(shared[39]);
    let mut reads = Vec::new();
    for i in 0..6 {
        let own = random_ref(40, 2000 + i);
        reads.push(Rd { pos: 50, seq: read_seq(&r, 50, &cig, &[&own]), cigar: cig.clone(), rev: i % 2 == 1 });
    }
    for i in 0..3 {
        reads.push(Rd { pos: 50, seq: read_seq(&r, 50, &cig, &[&shared]), cigar: cig.clone(), rev: i % 2 == 1 });
    }
    let (_vcf, cands, _bed) = run_case("long_ins", &r, reads);
    let ins: Vec<_> = cands.iter().filter(|c| c["type"] == "INS").collect();
    assert_eq!(ins.len(), 1, "{ins:?}");
    assert_eq!(ins[0]["alt_reads"], "3");
    assert_eq!(ins[0]["alt"].as_bytes()[1..], shared[..]);
}

/// P2-5: indels before the first reference base use the right-padded VCF form.
#[test]
fn qc_origin_indels() {
    let mut r = random_ref(200, 17);
    r[0] = b'C';
    r[3] = b'T';
    // Insertion "GT" before position 0 (ref[0] = C != G: cannot shift right).
    let cig = vec![Cigar::Ins(2), Cigar::Match(100)];
    let reads = (0..6)
        .map(|i| Rd { pos: 0, seq: read_seq(&r, 0, &cig, &[b"GT"]), cigar: cig.clone(), rev: i % 2 == 1 })
        .collect();
    let (vcf, cands, _bed) = run_case("origin_ins", &r, reads);
    assert_af_le_one(&cands);
    let v = vcf.iter().find(|v| v[7].contains("TYPE=INS")).unwrap_or_else(|| panic!("{vcf:?}"));
    assert_eq!((v[1].as_str(), v[3].as_str(), v[4].as_str()), ("1", "C", "GTC"));
    assert!(v[7].contains(";ORIGIN"));

    // Deletion of positions 0..3 (ref[0] = C != ref[3] = T: cannot shift right).
    let cig = vec![Cigar::Del(3), Cigar::Match(100)];
    let reads = (0..6)
        .map(|i| Rd { pos: 0, seq: read_seq(&r, 0, &cig, &[]), cigar: cig.clone(), rev: i % 2 == 1 })
        .collect();
    let (vcf, cands, _bed) = run_case("origin_del", &r, reads);
    assert_af_le_one(&cands);
    let v = vcf.iter().find(|v| v[7].contains("TYPE=DEL")).unwrap_or_else(|| panic!("{vcf:?}"));
    let want_ref = String::from_utf8(r[..4].to_vec()).unwrap();
    assert_eq!((v[1].as_str(), v[3].as_str(), v[4].as_str()), ("1", want_ref.as_str(), "T"));
    // dp = 6 < lowcov_abs (map-hifi 10): TIER stays RED, FILTER LowCov.
    assert!(v[7].contains(";ORIGIN") && v[7].contains(";LOWCOV") && v[7].contains("TIER=RED") && v[6] == "LowCov");
}

/// Re-review P1: an indel that left-normalises into an `N` skip of its read
/// is not supported by that read (the anchor is not in the read's DP).
#[test]
fn qc_normalised_into_refskip() {
    // ... G T | A x 11 | G ...: the run starts at 52, inside the read's skip 51..55.
    let mut r = random_ref(200, 23);
    r[50] = b'G';
    r[51] = b'T';
    for x in &mut r[52..63] {
        *x = b'A';
    }
    r[63] = b'G';
    for (name, op, ins) in [("refskip_ins", Cigar::Ins(1), true), ("refskip_del", Cigar::Del(1), false)] {
        let mut reads = Vec::new();
        for i in 0..6 {
            let cig = vec![Cigar::Match(40)];
            reads.push(Rd { pos: 40, seq: read_seq(&r, 40, &cig, &[]), cigar: cig, rev: i % 2 == 1 });
        }
        for i in 0..12 {
            let cig = vec![Cigar::Match(1), Cigar::RefSkip(4), Cigar::Match(5), op, Cigar::Match(5)];
            let chunks: Vec<&[u8]> = if ins { vec![b"A"] } else { vec![] };
            reads.push(Rd { pos: 50, seq: read_seq(&r, 50, &cig, &chunks), cigar: cig, rev: i % 2 == 1 });
        }
        let (_vcf, cands, _bed) = run_case(name, &r, reads);
        assert_af_le_one(&cands);
        let t = if ins { "INS" } else { "DEL" };
        assert!(cands.iter().all(|c| c["type"] != t), "{name}: {cands:?}");
    }
}

/// Re-review P2: an origin allele is tested, fitted and counted (BH m) in
/// the class of its own context (position 0), not that of anchor + 1.
#[test]
fn qc_origin_allele_class_consistent() {
    // C + A x 20 + G ...: position 0 is HP1, anchor + 1 is HP8+.
    let mut r = random_ref(200, 29);
    r[0] = b'C';
    for x in &mut r[1..21] {
        *x = b'A';
    }
    r[21] = b'G';
    let mut reads = Vec::new();
    for i in 0..20 {
        let (cig, chunks): (Vec<Cigar>, Vec<&[u8]>) = if i < 4 {
            (vec![Cigar::Ins(1), Cigar::Match(100)], vec![b"G"])
        } else {
            (vec![Cigar::Match(100)], vec![])
        };
        reads.push(Rd { pos: 0, seq: read_seq(&r, 0, &cig, &chunks), cigar: cig, rev: i % 2 == 1 });
    }
    let qp = QcParams { min_fit_sites: 0, ..params(1, None) };
    let (_vcf, cands, _bed, json) = run_case_p("origin_class", &r, reads, qp);
    let ins: Vec<_> = cands.iter().filter(|c| c["type"] == "INS").collect();
    assert_eq!(ins.len(), 1, "{ins:?}");
    assert_eq!((ins[0]["class"].as_str(), ins[0]["alt"].as_str(), ins[0]["flags"].as_str()), ("INS_HP1", "GC", "ORIGIN"));
    let bg = |c: &str| json["background"].as_array().unwrap().iter().find(|b| b["class"] == c).unwrap().clone();
    // The origin insertion (AF 0.2 < 0.3, so in the fit) lands in INS_HP1 only.
    assert_eq!(bg("INS_HP8+")["p"], 1e-6, "{}", bg("INS_HP8+"));
    assert!(bg("INS_HP1")["p"].as_f64().unwrap() > 1e-6, "{}", bg("INS_HP1"));
    // BH m: assessable positions 0..100 give 100 ordinary slots + 1 origin slot;
    // anchors 0..=13 have a run of >= 8 after them.
    let m = |c: &str| bg(c)["m"].as_u64().unwrap();
    let total: u64 = (0..9).map(|b| m(&format!("INS_{}", ["HP1", "HP2", "HP3", "HP4", "HP5", "HP6", "HP7", "HP8+", "STR"][b]))).sum();
    assert_eq!(total, 101);
    assert_eq!(m("INS_HP8+"), 13);
}

/// Minority-SNP clusters with the default refinements (A: nominal p <= 1e-3,
/// C: strand ratio): a strand-balanced cluster survives, a forward-only one
/// (strand-biased error pattern) is removed; with `--cluster-site-p 1
/// --cluster-strand off` both are reported.
#[test]
fn qc_snp_clusters_default_refinements() {
    let r = random_ref(6_000, 31);
    let balanced: Vec<usize> = (0..6).map(|k| 1_500 + 60 * k).collect();
    let fwd_only: Vec<usize> = (0..6).map(|k| 4_000 + 60 * k).collect();
    let mut reads = Vec::new();
    let mut i = 0usize;
    let mut s = 0usize;
    while s + 1_000 <= r.len() {
        let rev = i % 2 == 1;
        let pair = i / 2;
        let mut seq = r[s..s + 1_000].to_vec();
        for (k, &p) in balanced.iter().chain(&fwd_only).enumerate() {
            let hit = if k < balanced.len() { pair % 4 == 0 } else { !rev && pair % 2 == 0 };
            if p >= s && p < s + 1_000 && hit {
                seq[p - s] = alt_of(r[p]);
            }
        }
        reads.push(Rd { pos: s, cigar: vec![Cigar::Match(1_000)], seq, rev });
        i += 1;
        s += 25;
    }
    let clusters = |bed: &str| -> Vec<(u64, u64, String)> {
        bed.lines()
            .filter(|l| l.contains("\tsnp_cluster\t"))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[1].parse().unwrap(), f[2].parse().unwrap(), f[5].to_string())
            })
            .collect()
    };
    let rd = || reads.iter().map(|x| Rd { pos: x.pos, cigar: x.cigar.clone(), seq: x.seq.clone(), rev: x.rev }).collect();
    let (_v, _c, bed, json) = run_case_p("clusters_default", &r, rd(), params(1, None));
    let cl = clusters(&bed);
    assert_eq!(cl.len(), 1, "{bed}");
    assert!(cl[0].0 <= 1_500 && cl[0].1 > 1_800 && cl[0].1 < 4_000, "{cl:?}");
    assert_eq!(cl[0].2, "orange");
    assert_eq!(json["params"]["cluster_site_p"], 1e-3);
    assert_eq!(json["params"]["cluster_strand"], "ratio");
    assert!(json["params"]["cluster_density_p"].is_null());

    let off = QcParams {
        cluster_site_p: 1.0,
        cluster_strand: cramhook::qc::ClusterStrand::Off,
        ..params(1, None)
    };
    let (_v, _c, bed, _j) = run_case_p("clusters_off", &r, rd(), off);
    let cl = clusters(&bed);
    assert_eq!(cl.len(), 2, "{bed}");
    assert!(cl[1].0 <= 4_000 && cl[1].1 > 4_300, "{cl:?}");
    // Without the density test the BED extra is the original one (no kmin).
    assert!(bed.lines().filter(|l| l.contains("\tsnp_cluster\t")).all(|l| l.ends_with("\tsites=6;window=500")), "{bed}");

    // Density test on: kmin is reported. (On this 6 kb contig the 6 sites
    // themselves set λ; alpha 0.01 removes the cluster, alpha 1 keeps it.)
    let dens = QcParams { cluster_density_p: Some(0.01), ..params(1, None) };
    let (_v, _c, bed, _j) = run_case_p("clusters_density_strict", &r, rd(), dens);
    assert!(!bed.contains("\tsnp_cluster\t"), "{bed}");
    let dens = QcParams { cluster_density_p: Some(1.0), ..params(1, None) };
    let (_v, _c, bed, _j) = run_case_p("clusters_density", &r, rd(), dens);
    let lines: Vec<&str> = bed.lines().filter(|l| l.contains("\tsnp_cluster\t")).collect();
    assert!(!lines.is_empty() && lines.iter().all(|l| l.contains(";kmin=")), "{bed}");
}

/// Grouping similar long insertions: reads carry a 60 bp insertion with 1-3
/// errors, in three noisy copies. Default (grouping on) merges them into one
/// RED INS with INFO GROUPED; `--ins-group-min-len 0` splits them into multiple
/// low-support alleles and drops the `grouped` column; 1 and 8 threads agree.
#[test]
fn qc_ins_grouping() {
    let mut r = random_ref(1_000, 55);
    let ins = random_ref(60, 77);
    // No left shift: the preceding base differs from the insertion's last base.
    r[499] = alt_of(ins[59]);
    // Three noisy copies: 1, 2 and 3 substitutions inside the insertion.
    let mut variants: Vec<Vec<u8>> = Vec::new();
    for errs in [&[5usize][..], &[10, 20][..], &[30, 40, 50][..]] {
        let mut v = ins.clone();
        for &p in errs {
            v[p] = alt_of(v[p]);
        }
        variants.push(v);
    }
    // Each copy three times, plus one reference read covering the anchor:
    // 9 alt / 10 DP -> AF 0.9, RED.
    let cig = vec![Cigar::Match(100), Cigar::Ins(60), Cigar::Match(400)];
    let mut reads = Vec::new();
    for (k, v) in variants.iter().enumerate() {
        for i in 0..3 {
            reads.push(Rd {
                pos: 400,
                seq: read_seq(&r, 400, &cig, &[v]),
                cigar: cig.clone(),
                rev: (k + i) % 2 == 1,
            });
        }
    }
    let rcig = vec![Cigar::Match(400)];
    reads.push(Rd { pos: 400, seq: read_seq(&r, 400, &rcig, &[]), cigar: rcig, rev: false });

    let (vcf, cands, bed, _json) = run_case_p("ins_group", &r, reads.clone(), params(1, None));
    let ins_c: Vec<_> = cands.iter().filter(|c| c["type"] == "INS").collect();
    assert_eq!(ins_c.len(), 1, "{ins_c:?}");
    assert_eq!((ins_c[0]["alt_reads"].as_str(), ins_c[0]["dp"].as_str()), ("9", "10"));
    assert_eq!(ins_c[0]["grouped"], "3", "{ins_c:?}");
    let v = vcf.iter().find(|v| v[7].contains("TYPE=INS")).unwrap_or_else(|| panic!("{vcf:?}"));
    assert!(v[7].contains("TIER=RED") && v[7].contains("GROUPED=3"), "{v:?}");

    // Threads 1 vs 8 agree (varying chunk size).
    let q8 = QcParams { threads: 8, chunk_size: Some(400), ..params(1, None) };
    let (vcf8, cands8, bed8, _json8) = run_case_p("ins_group_t8", &r, reads.clone(), q8);
    assert_eq!(vcf, vcf8, "VCF differs between threads 1 and 8");
    assert_eq!(bed, bed8, "BED differs between threads 1 and 8");
    assert_eq!(cands, cands8, "candidates differ between threads 1 and 8");

    // Grouping off: multiple low-support alleles, no `grouped` column / INFO.
    let (vcf0, cands0, _bed0, _json0) = run_case_p("ins_group_off", &r, reads, QcParams { ins_group_min_len: 0, ..params(1, None) });
    let ins0: Vec<_> = cands0.iter().filter(|c| c["type"] == "INS").collect();
    assert_eq!(ins0.len(), 3, "{ins0:?}");
    assert!(ins0.iter().all(|c| !c.contains_key("grouped")), "{ins0:?}");
    assert!(ins0.iter().all(|c| c["alt_reads"] == "3"), "{ins0:?}");
    assert!(vcf0.iter().filter(|v| v[7].contains("TYPE=INS")).all(|v| !v[7].contains("GROUPED")), "{vcf0:?}");
}

/// Clip hotspots: a planted junction (a reference whose middle is absent from
/// the reads' source) makes reads clip on both sides of the breakpoint, giving
/// one paired (`clip_both`) RED `clip_hotspot`. A locus with only short 30 bp
/// clips stays below `--clip-min-len` and gives nothing. Same-side clips
/// straddling a chunk boundary (records starting in different chunks) merge
/// into one region, and threads 1 vs 8 produce identical BED and JSON.
#[test]
fn qc_clip_hotspot_junction_and_length_filter() {
    let l = 4_000usize;
    let r = random_ref(l, 91);
    let junction = 2_000usize;
    let mut reads = Vec::new();
    // 10 reads ending at the junction (right clip at junction - 1).
    for i in 0..10 {
        let cig = vec![Cigar::Match(500), Cigar::SoftClip(120)];
        let mut seq = r[junction - 500..junction].to_vec();
        seq.extend_from_slice(&[b'A'; 120]);
        reads.push(Rd { pos: junction - 500, seq, cigar: cig, rev: i % 2 == 1 });
    }
    // 10 reads starting at the junction (left clip at junction).
    for i in 0..10 {
        let cig = vec![Cigar::SoftClip(120), Cigar::Match(500)];
        let mut seq = vec![b'A'; 120];
        seq.extend_from_slice(&r[junction..junction + 500]);
        reads.push(Rd { pos: junction, seq, cigar: cig, rev: i % 2 == 1 });
    }
    // Same-side (left) clips at 990 and 1010, within `--clip-merge` (50). With
    // `--chunk-size 1000` the records start in different chunks; the clips must
    // still merge into one region after the chunk merge.
    for start in [990usize, 1_010] {
        for i in 0..10 {
            let cig = vec![Cigar::SoftClip(120), Cigar::Match(400)];
            let mut seq = vec![b'A'; 120];
            seq.extend_from_slice(&r[start..start + 400]);
            reads.push(Rd { pos: start, seq, cigar: cig, rev: i % 2 == 1 });
        }
    }
    // A separate locus whose 10 reads carry only a 30 bp clip (below min len).
    for i in 0..10 {
        let cig = vec![Cigar::Match(200), Cigar::SoftClip(30)];
        let mut seq = r[2_800..3_000].to_vec();
        seq.extend_from_slice(&[b'A'; 30]);
        reads.push(Rd { pos: 2_800, seq, cigar: cig, rev: i % 2 == 1 });
    }
    let hot = |bed: &str| -> Vec<Vec<String>> {
        bed.lines()
            .filter(|l| l.contains("\tclip_hotspot\t"))
            .map(|l| l.split('\t').map(|s| s.to_string()).collect())
            .collect()
    };
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/qc-scan/test");

    // The default `--clip-high-s 100` is unreachable on this toy genome (no
    // background clips, so `s ~ log10(k)`); lower it so the red path is still
    // exercised here.
    let qp = QcParams { clip_high_s: 5.0, ..params(1, None) };
    let (_vcf, _cands, bed, json) = run_case_p("clip_junction", &r, reads.clone(), qp.clone());
    let h = hot(&bed);
    // The 30 bp locus gives nothing; the boundary and junction loci do.
    assert_eq!(h.len(), 2, "short clips must not form a hotspot:\n{bed}");
    let junction_h = h.iter().find(|x| x[1] == "1999").unwrap_or_else(|| panic!("{bed}"));
    assert_eq!((junction_h[5].as_str(), junction_h[6].as_str()), ("red", "clip_both"), "{junction_h:?}");
    assert_eq!((junction_h[1].as_str(), junction_h[2].as_str()), ("1999", "2001"), "{junction_h:?}");
    assert!(junction_h[9].contains("left=10;right=10"), "{junction_h:?}");
    // The two same-side clips across the chunk boundary are ONE merged region.
    let boundary_h: Vec<&Vec<String>> = h.iter().filter(|x| x[1] == "990").collect();
    assert_eq!(boundary_h.len(), 1, "boundary clips not merged:\n{bed}");
    assert_eq!((boundary_h[0][2].as_str(), boundary_h[0][6].as_str()), ("1011", "clip_left"), "{bed}");
    assert!(boundary_h[0][9].contains("left=20"), "{:?}", boundary_h[0]);
    assert_eq!(json["alerts"]["by_type"]["clip_hotspot"]["red"], 2, "{json}");

    // With the length filter off the 30 bp locus also forms a hotspot.
    let q0 = QcParams { clip_min_len: 0, ..qp.clone() };
    let (_v, _c, bed0, _j) = run_case_p("clip_junction_nolen", &r, reads.clone(), q0);
    let h0 = hot(&bed0);
    assert_eq!(h0.len(), 3, "{bed0}");
    assert!(h0.iter().any(|x| x[6] == "clip_right" && x[1] == "2999"), "{bed0}");

    // Threads 1 vs 8 (chunk boundary at 1000, between the two same-side clips)
    // produce identical BED and identical normalised JSON.
    let q8 = QcParams { threads: 8, chunk_size: Some(1_000), ..qp.clone() };
    let (_v8, _c8, bed8, _j8) = run_case_p("clip_junction_t8", &r, reads, q8);
    assert_eq!(bed, bed8, "BED differs between --threads 1 and 8");
    assert_eq!(
        json_without_volatile(&base.join("clip_junction/out.json")),
        json_without_volatile(&base.join("clip_junction_t8/out.json")),
        "JSON differs between --threads 1 and 8"
    );
}

/// Depth bands on a 60 kb reference at 40×: a 3 kb block with 40× extra
/// high-MAPQ reads (`depth_high`, orange), a 3 kb block whose reads are split
/// off as MAPQ 0 (`depth_half`, orange) and a 2 kb no-read block (`depth_low`,
/// red). Threads 1 vs 8 produce identical BED and normalised JSON.
#[test]
fn qc_depth_bands() {
    let l = 60_000usize;
    let r = random_ref(l, 123);
    // (pos, cigar, seq, rev, mapq)
    let mut reads: Vec<(usize, Vec<Cigar>, Vec<u8>, bool, u8)> = Vec::new();
    let mut i = 0usize;
    // Base 40× everywhere except the reads that would cover the two blocks.
    let omitted = |s: usize| (39_000..42_000).contains(&s) || (49_000..53_000).contains(&s);
    for s in (0..=l - 1_000).step_by(25) {
        if omitted(s) {
            continue;
        }
        reads.push((s, vec![Cigar::Match(1_000)], r[s..s + 1_000].to_vec(), i % 2 == 1, 60));
        i += 1;
    }
    // 80× block [20000, 23000): 40× of extra high-MAPQ reads.
    for s in (19_000..=22_975).step_by(25) {
        reads.push((s, vec![Cigar::Match(1_000)], r[s..s + 1_000].to_vec(), i % 2 == 1, 60));
        i += 1;
    }
    // 0.5× block [50000, 53000): 20× of MAPQ-0 reads (a split tandem copy).
    for s in (49_000..=52_950).step_by(50) {
        reads.push((s, vec![Cigar::Match(1_000)], r[s..s + 1_000].to_vec(), i % 2 == 1, 0));
        i += 1;
    }

    let dir = workdir("depth_bands");
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr\n{}\n", std::str::from_utf8(&r).unwrap())).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), format!("chr\t{l}\t5\t{l}\t{}\n", l + 1)).unwrap();
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", l));
    let bam_path = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
    reads.sort_by_key(|r| r.0);
    for (n, (pos, cigar, seq, rev, mapq)) in reads.iter().enumerate() {
        let mut rec = Record::new();
        let qual = vec![30u8; seq.len()];
        rec.set(format!("r{n}").as_bytes(), Some(&CigarString(cigar.clone())), seq, &qual);
        rec.set_tid(0);
        rec.set_pos(*pos as i64);
        rec.set_mapq(*mapq);
        rec.set_flags(if *rev { 0x10 } else { 0 });
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        w.write(&rec).unwrap();
    }
    drop(w);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

    let p1 = dir.join("t1/out");
    let p8 = dir.join("t8/out");
    let s1 = run_qc(&bam_path, &fa, &p1, &params(1, None)).unwrap();
    let s8 = run_qc(&bam_path, &fa, &p8, &QcParams { threads: 8, chunk_size: Some(1_500), ..params(1, None) }).unwrap();
    assert!(s8.chunks > 10);
    assert_eq!(s1.contigs.0[0].1.depth_skipped, None);
    assert_eq!(s1.contigs.0[0].1.depth_mode, 40);
    assert_eq!(s1.reads.mean_aligned_span, 1_000.0);

    let bed = std::fs::read_to_string(p1.with_extension("bed")).unwrap();
    let regs: Vec<Vec<&str>> =
        bed.lines().filter(|l| l.contains("\tdepth_")).map(|l| l.split('\t').collect()).collect();
    let covering = |label: &str, s: u64, e: u64| {
        regs.iter().any(|x| {
            x[3] == label && x[1].parse::<u64>().unwrap() <= s && x[2].parse::<u64>().unwrap() >= e
        })
    };
    // 80× high-MAPQ block: depth_high, orange (no cluster/clip to upgrade it).
    let hi = regs.iter().find(|x| x[3] == "depth_high").unwrap_or_else(|| panic!("no depth_high:\n{bed}"));
    assert!(covering("depth_high", 20_000, 23_000), "{bed}");
    assert_eq!((hi[5], hi[4], hi[6]), ("orange", "0", "high"), "{hi:?}");
    assert!(hi[9].contains("windows=") && hi[9].contains("baseline="), "{hi:?}");
    // 0.5× MAPQ-0 block: depth_half, orange.
    assert!(covering("depth_half", 50_000, 53_000), "{bed}");
    assert!(regs.iter().any(|x| x[3] == "depth_half" && x[5] == "orange"), "{regs:?}");
    // No-read block: depth_low, red.
    assert!(covering("depth_low", 40_000, 42_000), "{bed}");
    let lo = regs.iter().find(|x| x[3] == "depth_low").unwrap();
    assert_eq!(lo[5], "red", "{lo:?}");

    // Findings: one per depth region, with the expected `likely` hint. (The new
    // `depth_step` signal also adds a finding at the circular-origin ramp.)
    let fbed = std::fs::read_to_string(p1.with_extension("findings.bed")).unwrap();
    let fl: Vec<&str> = fbed.lines().filter(|l| !l.starts_with('#')).collect();
    let band: Vec<&str> = fl
        .iter()
        .copied()
        .filter(|l| l.contains("depth_high:") || l.contains("depth_low:") || l.contains("depth_half:"))
        .collect();
    assert_eq!(band.len(), 3, "{fbed}");
    let likelies: Vec<&str> = band.iter().map(|l| l.split('\t').nth(6).unwrap()).collect();
    assert!(
        likelies.contains(&"collapsed_repeat")
            && likelies.contains(&"extra_sequence")
            && likelies.contains(&"false_duplication"),
        "{fbed}"
    );
    // The new `depth_step` signal adds its own finding at the circular-origin
    // coverage ramp (lowcov at the contig end).
    let origin = fl
        .iter()
        .find(|l| l.contains("depth_step:") && l.split('\t').nth(1).unwrap().parse::<u64>().unwrap() < 500);
    assert!(origin.is_some(), "no origin depth_step finding:\n{fbed}");

    // Threads 1 vs 8 agree.
    assert_eq!(bed, std::fs::read_to_string(p8.with_extension("bed")).unwrap());
    assert_eq!(fbed, std::fs::read_to_string(p8.with_extension("findings.bed")).unwrap());
    assert_eq!(json_without_volatile(&p1.with_extension("json")), json_without_volatile(&p8.with_extension("json")));
}

/// Misjoin coverage (formerly the `bridge_break` test): a 40 kb reference with
/// ~5 kb reads at 40x, modelled as reads tiling a larger window and clipped at
/// the contig ends (as in real data). A clean reference gives no `clip_hotspot`;
/// a misjoin at 20 kb (reads spanning it clip at the junction) gives one red
/// `clip_both` hotspot whose finding is `likely = misjoin`. Threads 1 vs 8 agree.
#[test]
fn qc_misjoin_clip_hotspot_and_finding() {
    const L: usize = 40_000;
    const JL: usize = 20_000;
    const RL: i64 = 5_000;
    let r = random_ref(L, 555);

    // One read from a window start `s` (may be negative / past the end): the
    // parts outside the contig are soft-clipped. A read whose window spans
    // `JL` is split into a left and a right record (as a real split alignment).
    let mut i = 0usize;
    let mut reads = |break_at: Option<usize>| -> Vec<Rd> {
        let mut v: Vec<Rd> = Vec::new();
        let push = |v: &mut Vec<Rd>, a: i64, b: i64, lclip: i64, rclip: i64, i: &mut usize| {
            if b <= a {
                return;
            }
            let (a, b) = (a as usize, b as usize);
            let mut cigar = Vec::new();
            if lclip > 0 {
                cigar.push(Cigar::SoftClip(lclip as u32));
            }
            cigar.push(Cigar::Match((b - a) as u32));
            if rclip > 0 {
                cigar.push(Cigar::SoftClip(rclip as u32));
            }
            let mut seq = vec![b'A'; lclip.max(0) as usize];
            seq.extend_from_slice(&r[a..b]);
            seq.extend(std::iter::repeat(b'A').take(rclip.max(0) as usize));
            v.push(Rd { pos: a, cigar, seq, rev: *i % 2 == 1 });
            *i += 1;
        };
        let mut s = -RL;
        while s < L as i64 {
            let e = s + RL;
            let a = s.max(0);
            let b = e.min(L as i64);
            let br = break_at.filter(|&br| s < br as i64 && e > br as i64);
            match br {
                Some(br) => {
                    let br = br as i64;
                    // A split alignment: the left record is soft-clipped on the
                    // right at the junction, the right record on the left.
                    push(&mut v, a, br, a - s, e - br, &mut i);
                    push(&mut v, br, b, br - s, e - b, &mut i);
                }
                None => push(&mut v, a, b, a - s, e - b, &mut i),
            }
            s += 125;
        }
        v
    };

    // The default `--clip-high-s 100` is unreachable on this toy genome (the
    // end clips raise the background, so `s` stays small); lower it so the red
    // path is still exercised.
    let qp = QcParams { clip_high_s: 5.0, ..params(1, None) };
    let (_v, _c, bed, json) = run_case_p("misjoin_clean", &r, reads(None), qp.clone());
    assert!(!bed.contains("\tclip_hotspot\t"), "clean reference must give no clip_hotspot:\n{bed}");
    assert_eq!(json["alerts"]["by_type"]["clip_hotspot"]["red"], 0, "{json}");

    let mis = reads(Some(JL));
    let (_vcf, _cands, bed, json) = run_case_p("misjoin", &r, mis.clone(), qp.clone());
    let hits: Vec<Vec<&str>> =
        bed.lines().filter(|l| l.contains("\tclip_hotspot\t")).map(|l| l.split('\t').collect()).collect();
    assert_eq!(hits.len(), 1, "{bed}");
    let h = &hits[0];
    assert_eq!(h[6], "clip_both", "{h:?}");
    let start: u64 = h[1].parse().unwrap();
    let end: u64 = h[2].parse().unwrap();
    assert!(start <= JL as u64 && end >= JL as u64, "clip_hotspot {start}-{end} not at {JL}: {h:?}");
    assert!(h[9].contains("s=") && h[9].contains("f="), "{h:?}");
    assert_eq!(json["alerts"]["by_type"]["clip_hotspot"]["red"], 1, "{json}");

    // Findings: the junction is one finding with `likely = misjoin`.
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/qc-scan/test");
    let fbed = std::fs::read_to_string(base.join("misjoin/out.findings.bed")).unwrap();
    let fl: Vec<Vec<&str>> =
        fbed.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').collect()).collect();
    assert_eq!(fl.len(), 1, "{fbed}");
    assert_eq!((fl[0][6], fl[0][3]), ("misjoin", "finding_1"), "{fbed}");
    let (fs, fe): (u64, u64) = (fl[0][1].parse().unwrap(), fl[0][2].parse().unwrap());
    assert!(fs <= JL as u64 && fe >= JL as u64, "finding {fl:?} not at {JL}");

    // Threads 1 vs 8 agree (different workdir, so compare BED bodies).
    let q8 = QcParams { threads: 8, chunk_size: Some(1_500), ..qp.clone() };
    let (_v8, _c8, bed8, _j8) = run_case_p("misjoin_t8", &r, mis, q8);
    assert_eq!(bed, bed8, "BED differs between --threads 1 and 8");
    assert_eq!(fbed, std::fs::read_to_string(base.join("misjoin_t8/out.findings.bed")).unwrap());
    assert_eq!(
        json_without_volatile(&base.join("misjoin/out.json")),
        json_without_volatile(&base.join("misjoin_t8/out.json"))
    );
}

/// Depth steps: a 60 kb contig at ~40× on `[0, 30 kb)` and ~10× on
/// `[30 kb, 60 kb)` gives a `depth_step` at the boundary plus a paired
/// `depth_segment`; `depth_segment` is output only (BED, no alert type), and
/// threads 1 vs 8 produce identical BED, findings and normalised JSON.
#[test]
fn qc_depth_step_region_and_thread_identity() {
    const L: usize = 60_000;
    const STEP: usize = 30_000;
    let r = random_ref(L, 4242);
    let mut reads: Vec<Rd> = Vec::new();
    let mut i = 0usize;
    // ~40× left of STEP.
    for s in (0..STEP - 1_000).step_by(25) {
        reads.push(Rd { pos: s, cigar: vec![Cigar::Match(1_000)], seq: r[s..s + 1_000].to_vec(), rev: i % 2 == 1 });
        i += 1;
    }
    // ~10× right of STEP (the last read is clipped at the contig end).
    for s in (STEP..=L - 1_000).step_by(100) {
        reads.push(Rd { pos: s, cigar: vec![Cigar::Match(1_000)], seq: r[s..s + 1_000].to_vec(), rev: i % 2 == 1 });
        i += 1;
    }
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/qc-scan/test");
    let (_v, _c, bed, json) = run_case_p("depth_step", &r, reads.clone(), params(1, None));
    let steps: Vec<Vec<&str>> =
        bed.lines().filter(|l| l.contains("\tdepth_step\t")).map(|l| l.split('\t').collect()).collect();
    assert!(!steps.is_empty(), "no depth_step in:\n{bed}");
    assert!(
        steps.iter().any(|x| x[6] == "step_down"
            && (x[1].parse::<i64>().unwrap() - STEP as i64).abs() < 2_000
            && matches!(x[5], "red" | "orange" | "yellow")),
        "no step_down near {STEP}:\n{bed}"
    );
    // depth_step is counted under alerts (PASS or lowcov); depth_segment is not.
    let ds = &json["alerts"]["by_type"]["depth_step"];
    let total = |v: &serde_json::Value| {
        v["red"].as_u64().unwrap_or(0) + v["orange"].as_u64().unwrap_or(0) + v["yellow"].as_u64().unwrap_or(0)
    };
    assert!(total(ds) + total(&ds["lowcov"]) >= 1, "{ds}");
    assert!(json["alerts"]["by_type"].get("depth_segment").is_none(), "{json}");
    if bed.contains("\tdepth_segment\t") {
        let fbed = std::fs::read_to_string(base.join("depth_step/out.findings.bed")).unwrap();
        assert!(!fbed.contains("depth_segment:"), "depth_segment must not be a finding member:\n{fbed}");
    }

    // Threads 1 vs 8 agree (different workdir, so compare BED bodies).
    let q8 = QcParams { threads: 8, chunk_size: Some(1_500), ..params(1, None) };
    let (_v8, _c8, bed8, _j8) = run_case_p("depth_step_t8", &r, reads, q8);
    assert_eq!(bed, bed8, "BED differs between --threads 1 and 8");
    assert_eq!(
        std::fs::read_to_string(base.join("depth_step/out.findings.bed")).unwrap(),
        std::fs::read_to_string(base.join("depth_step_t8/out.findings.bed")).unwrap(),
        "findings differ between --threads 1 and 8"
    );
    assert_eq!(
        json_without_volatile(&base.join("depth_step/out.json")),
        json_without_volatile(&base.join("depth_step_t8/out.json")),
        "JSON differs between --threads 1 and 8"
    );
}

/// The `.depth.bedgraph` track is *aligned-base* depth only (reads with a base,
/// deletions not counted, like IGV / `samtools depth`), while qc's DP keeps
/// counting reads whose only contribution is a deletion. All 20 reads carry the
/// same `100M 50D 100M` alignment, so the bedGraph dips to 0 over `200..250`
/// while the DP-based JSON `depth_mean` (20× over 250 reference bp of 600) does
/// not.
#[test]
fn qc_bedgraph_aligned_depth_excludes_deletions() {
    let l = 600usize;
    let r = random_ref(l, 123);
    let (del_start, del_len) = (200usize, 50usize);
    let cig = vec![Cigar::Match(100), Cigar::Del(del_len as u32), Cigar::Match(100)];
    let reads = (0..20)
        .map(|i| Rd { pos: 100, seq: read_seq(&r, 100, &cig, &[]), cigar: cig.clone(), rev: i % 2 == 1 })
        .collect();
    let (vcf, cands, _bed, json) = run_case_p("bedgraph_aligned", &r, reads, params(1, None));

    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/qc-scan/test");
    let bg = std::fs::read_to_string(base.join("bedgraph_aligned/out.depth.bedgraph")).unwrap();
    let depth_at = |p: u64| -> u64 {
        bg.lines()
            .find_map(|ln| {
                let f: Vec<&str> = ln.split('\t').collect();
                let (s, e) = (f[1].parse::<u64>().unwrap(), f[2].parse::<u64>().unwrap());
                (p >= s && p < e).then(|| f[3].parse().unwrap())
            })
            .unwrap_or_else(|| panic!("no bedGraph row covering {p}:\n{bg}"))
    };
    // Aligned bases outside the deletion: 20 reads; inside it: 0.
    assert_eq!(depth_at(150), 20, "{bg}");
    assert_eq!(depth_at(del_start as u64 + del_len as u64 / 2), 0, "{bg}");
    assert_eq!(depth_at(275), 20, "{bg}");

    // DP (base or deletion) does not dip: each read's deletion still contributes,
    // so the contig mean is 20 × 250 / 600 (the deletion is covered), not the
    // aligned-base mean 20 × 200 / 600.
    assert_eq!(json["contigs"]["chr"]["depth_mean"], 8.333, "{json}");
    assert_eq!(json["contigs"]["chr"]["depth_mode"], 20, "{json}");
    // The deletion is still called (deletions count toward DP / support).
    assert!(cands.iter().any(|c| c["type"] == "DEL" && c["dp"] == "20"), "{cands:?}");
    assert!(vcf.iter().any(|v| v[7].contains("TYPE=DEL")), "{vcf:?}");
}

/// `qc.gc.bedgraph`: one row per fixed window, 0-based half-open, contigs in
/// the BAM header order (the depth track's order), the last window short, and
/// non-ACGT bases excluded from both counts (a window that is all `N` is
/// 0.0000). `--gc-window 0` writes no file. Values match a direct computation
/// from the FASTA.
#[test]
fn qc_gc_bedgraph() {
    let dir = workdir("gc_bedgraph");
    let a = random_ref(600, 404);
    let mut b = random_ref(250, 505);
    for x in &mut b[100..200] {
        *x = b'N';
    }
    // FASTA order is chrB, chrA; the BAM header (and both bedGraphs) is
    // chrA, chrB.
    let fa = dir.join("ref.fa");
    let mut text = String::new();
    for (name, seq) in [("chrB", &b), ("chrA", &a)] {
        text.push('>');
        text.push_str(name);
        text.push('\n');
        text.push_str(std::str::from_utf8(seq).unwrap());
        text.push('\n');
    }
    std::fs::write(&fa, text).unwrap();
    // ">chrB\n" (6) + 250 bp + newline = chrA header at 257, its sequence at 263.
    std::fs::write(dir.join("ref.fa.fai"), "chrB\t250\t6\t250\t251\nchrA\t600\t263\t600\t601\n").unwrap();

    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chrA").push_tag(b"LN", 600));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chrB").push_tag(b"LN", 250));
    let bam = dir.join("reads.bam");
    let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();
    let cig = CigarString(vec![Cigar::Match(100)]);
    for (tid, seq) in [(0i32, &a), (1i32, &b)] {
        for i in 0..6 {
            let mut rec = Record::new();
            rec.set(format!("r{tid}_{i}").as_bytes(), Some(&cig), &seq[100..200], &vec![30u8; 100]);
            rec.set_tid(tid);
            rec.set_pos(100);
            rec.set_mapq(60);
            rec.set_flags(if i % 2 == 1 { 0x10 } else { 0 });
            rec.set_mtid(-1);
            rec.set_mpos(-1);
            w.write(&rec).unwrap();
        }
    }
    drop(w);
    bam::index::build(&bam, None, bam::index::Type::Bai, 1).unwrap();

    let prefix = dir.join("out");
    run_qc(&bam, &fa, &prefix, &params(1, None)).unwrap();
    let gc_text = std::fs::read_to_string(prefix.with_extension("gc.bedgraph")).unwrap();
    let rows: Vec<(String, u64, u64, f64)> = gc_text
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 4, "row is not `contig start end gc`: {l:?}");
            (f[0].to_string(), f[1].parse().unwrap(), f[2].parse().unwrap(), f[3].parse().unwrap())
        })
        .collect();
    // No header line; same contig order as the depth bedGraph.
    let order = |t: &str| {
        let mut v: Vec<String> = Vec::new();
        for l in t.lines() {
            let c = l.split('\t').next().unwrap().to_string();
            if v.last() != Some(&c) {
                v.push(c);
            }
        }
        v
    };
    let depth_text = std::fs::read_to_string(prefix.with_extension("depth.bedgraph")).unwrap();
    assert_eq!(order(&gc_text), order(&depth_text));
    assert_eq!(order(&gc_text), ["chrA".to_string(), "chrB".to_string()]);

    // Per contig: start at 0, no gaps, end == contig length, values in [0, 1].
    for (name, len) in [("chrA", 600u64), ("chrB", 250u64)] {
        let rs: Vec<_> = rows.iter().filter(|r| r.0 == name).collect();
        assert!(!rs.is_empty(), "{gc_text}");
        assert_eq!(rs[0].1, 0, "{name}: {gc_text}");
        for pair in rs.windows(2) {
            assert_eq!(pair[0].2, pair[1].1, "{name}: gap/overlap in {gc_text}");
        }
        assert_eq!(rs.last().unwrap().2, len, "{name}: {gc_text}");
        assert!(rs.iter().all(|r| r.2 - r.1 == 100 || r.2 == len), "{name}: {gc_text}");
        assert!(rs.iter().all(|r| (0.0..=1.0).contains(&r.3)), "{name}: {gc_text}");
    }

    // A known region matches a direct FASTA computation, formatted to 4 decimals.
    let direct_gc = |s: &[u8]| {
        let (gc, acgt) = s.iter().fold((0u64, 0u64), |(g, n), &b| match b.to_ascii_uppercase() {
            b'G' | b'C' => (g + 1, n + 1),
            b'A' | b'T' => (g, n + 1),
            _ => (g, n),
        });
        if acgt == 0 { 0.0 } else { gc as f64 / acgt as f64 }
    };
    let want = direct_gc(&a[0..100]);
    let line = gc_text.lines().find(|l| l.starts_with("chrA\t0\t100\t")).unwrap_or_else(|| panic!("{gc_text}"));
    assert_eq!(line.split('\t').nth(3).unwrap(), format!("{want:.4}"), "{gc_text}");
    let want = direct_gc(&b[200..250]);
    assert_eq!(
        gc_text.lines().find(|l| l.starts_with("chrB\t200\t250\t")).unwrap_or_else(|| panic!("{gc_text}")),
        format!("chrB\t200\t250\t{want:.4}")
    );
    // All-N window: no A/C/G/T -> 0.0000.
    assert!(gc_text.lines().any(|l| l == "chrB\t100\t200\t0.0000"), "{gc_text}");

    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(prefix.with_extension("json")).unwrap()).unwrap();
    assert_eq!(json["gc_window"], 100, "{json}");
    assert_eq!(json["params"]["gc_window"], 100, "{json}");

    // --gc-window 0: no file.
    let p0 = dir.join("nogc/out");
    run_qc(&bam, &fa, &p0, &QcParams { gc_window: 0, ..params(1, None) }).unwrap();
    assert!(!p0.with_extension("gc.bedgraph").exists());
    let j0: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(p0.with_extension("json")).unwrap()).unwrap();
    assert_eq!(j0["gc_window"], 0, "{j0}");
}

/// Gray tier (homopolymer indels): a *minority / strand-biased* yellow INS/DEL
/// whose INFO `HP` is at least `--gray-hp-min` (default 5) is `TIER=GRAY` with
/// `HPGRAY`, still written to the VCF but never a finding member and not counted
/// in the red/orange/yellow alerts. `--no-gray-homopolymer` restores yellow
/// (which does form a finding). A majority HP indel stays RED and a minority SNP
/// stays YELLOW either way.
#[test]
fn qc_gray_homopolymer_minority_indel() {
    const L: usize = 4_000;
    const RL: usize = 1_000;
    const STEP: usize = 25;
    const HP_RED: usize = 1_200; // run of 8 `A`, deleted by every spanning read
    const HP_GRAY: usize = 3_000; // run of 8 `A`, deleted by a minority only
    const HP_LEN: usize = 8;
    const SNP_MINOR: usize = 2_000;

    let mut r = random_ref(L, 321);
    for hp in [HP_RED, HP_GRAY] {
        r[hp - 1] = b'C';
        r[hp + HP_LEN] = b'G';
        for x in &mut r[hp..hp + HP_LEN] {
            *x = b'A';
        }
    }

    let mut reads: Vec<Rd> = Vec::new();
    let mut rng = Lcg(777);
    let mut i = 0usize;
    let mut s = 0usize;
    while s + RL <= L {
        let e = s + RL;
        let rev = i % 2 == 1;
        let pair = i / 2;
        let mut seq = Vec::new();
        let mut ops: Vec<Cigar> = Vec::new();
        let mut m = 0u32;
        let mut p = s;
        while p < e {
            // One `A` from a run, at a per-read offset so the reads delete
            // different `A`s (all normalise to one allele). The majority run is
            // deleted by every spanning read, the gray run by a minority.
            let spans = |hp: usize| s < hp && e > hp + HP_LEN;
            let del = (spans(HP_RED) && p == HP_RED + pair % HP_LEN)
                || (pair % 4 == 0 && spans(HP_GRAY) && p == HP_GRAY + pair % HP_LEN);
            if del {
                if m > 0 {
                    ops.push(Cigar::Match(m));
                    m = 0;
                }
                ops.push(Cigar::Del(1));
                p += 1;
                continue;
            }
            // 0.5 % random substitution noise (the planted test's background),
            // so the beta-binomial fit is not dominated by the planted sites.
            let mut b = if rng.next() % 200 == 0 { b"ACGT"[(rng.next() % 4) as usize] } else { r[p] };
            if p == SNP_MINOR && pair % 4 == 0 {
                b = alt_of(b);
            }
            seq.push(b);
            m += 1;
            p += 1;
        }
        if m > 0 {
            ops.push(Cigar::Match(m));
        }
        reads.push(Rd { pos: s, cigar: ops, seq, rev });
        i += 1;
        s += STEP;
    }
    // `--unlinked-action flag`: the lone minority SNP is not dropped if noise
    // partners are tested (mirrors the planted-events test).
    let p = |threads, chunk| QcParams {
        unlinked_action: cramhook::qc::UnlinkedAction::Flag,
        ..params(threads, chunk)
    };

    let find = |vcf: &[Vec<String>], pos: usize| -> Option<Vec<String>> {
        vcf.iter().find(|v| v[1] == pos.to_string()).cloned()
    };
    let covers = |fbed: &str, pos: usize| -> bool {
        fbed.lines().filter(|l| !l.starts_with('#')).any(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let (s, e) = (f[1].parse::<u64>().unwrap(), f[2].parse::<u64>().unwrap());
            s <= pos as u64 && (pos as u64) < e
        })
    };

    // Default: gray on.
    let (vcf, cands, _bed, json) = run_case_p("gray_on", &r, reads.clone(), p(1, None));
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/qc-scan/test");
    let fbed = std::fs::read_to_string(base.join("gray_on/out.findings.bed")).unwrap();

    // Majority HP deletion: RED, never gray.
    let red = find(&vcf, HP_RED).unwrap_or_else(|| panic!("no DEL at {HP_RED}: {vcf:?}"));
    assert!(red[7].contains("TIER=RED") && red[7].contains("TYPE=DEL"), "{red:?}");
    assert!(!red[7].contains("HPGRAY"), "{red:?}");
    // Minority HP deletion: GRAY + HPGRAY, FILTER PASS, reason kept.
    let gray = find(&vcf, HP_GRAY).unwrap_or_else(|| panic!("no DEL at {HP_GRAY}: {vcf:?}"));
    assert_eq!(gray[6], "PASS", "{gray:?}");
    assert!(gray[7].contains("TIER=GRAY") && gray[7].contains(";HPGRAY") && gray[7].contains("REASON=minority"), "{gray:?}");
    // Minority SNP: YELLOW, never gray.
    let snp = find(&vcf, SNP_MINOR + 1).unwrap_or_else(|| panic!("no SNP at {SNP_MINOR}: {vcf:?}"));
    assert!(snp[7].contains("TIER=YELLOW") && snp[7].contains("TYPE=SNP") && !snp[7].contains("HPGRAY"), "{snp:?}");
    // candidates.tsv shows tier gray and the HPGRAY flag.
    let gc = cands
        .iter()
        .find(|c| c["type"] == "DEL" && c["pos0"].parse::<usize>().map(|p| p.abs_diff(HP_GRAY) <= 2).unwrap_or(false))
        .unwrap_or_else(|| panic!("no gray DEL candidate: {cands:?}"));
    assert_eq!((gc["tier"].as_str(), gc["reason"].as_str()), ("gray", "minority"), "{gc:?}");
    assert!(gc["flags"].contains("HPGRAY"), "{gc:?}");
    // No finding holds the gray site; a finding does hold the red site.
    assert!(!covers(&fbed, HP_GRAY), "gray site became a finding member:\n{fbed}");
    assert!(covers(&fbed, HP_RED), "red site lost its finding:\n{fbed}");
    // Alerts: gray counted apart; yellow is the SNP only.
    assert_eq!(json["alerts"]["gray"], 1, "{json}");
    assert_eq!(json["alerts"]["by_type"]["DEL"]["gray"], 1, "{json}");
    assert_eq!(json["alerts"]["yellow"], 1, "{json}");

    // --no-gray-homopolymer: the HP indel is yellow again and forms a finding.
    let off = QcParams { gray_homopolymer: false, ..p(1, None) };
    let (vcf_off, _c_off, _bed_off, json_off) = run_case_p("gray_off", &r, reads.clone(), off);
    let gray_off = find(&vcf_off, HP_GRAY).unwrap_or_else(|| panic!("no DEL at {HP_GRAY}: {vcf_off:?}"));
    assert!(gray_off[7].contains("TIER=YELLOW") && !gray_off[7].contains("HPGRAY"), "{gray_off:?}");
    let fbed_off = std::fs::read_to_string(base.join("gray_off/out.findings.bed")).unwrap();
    assert!(covers(&fbed_off, HP_GRAY), "yellow HP indel lost its finding:\n{fbed_off}");
    assert_eq!(json_off["alerts"]["gray"], 0, "{json_off}");
    assert_eq!(json_off["alerts"]["yellow"], 2, "{json_off}");

    // Threads 1 vs 8 agree.
    let q8 = QcParams { threads: 8, chunk_size: Some(400), ..p(1, None) };
    let (vcf8, _c8, _b8, _j8) = run_case_p("gray_on_t8", &r, reads, q8);
    assert_eq!(vcf, vcf8, "VCF differs between --threads 1 and 8");
    assert_eq!(fbed, std::fs::read_to_string(base.join("gray_on_t8/out.findings.bed")).unwrap());
    assert_eq!(
        json_without_volatile(&base.join("gray_on/out.json")),
        json_without_volatile(&base.join("gray_on_t8/out.json"))
    );
}

/// A `weak_majority` HP indel (strand-balanced AF > 0.5 that misses the red bar)
/// stays YELLOW even with gray on: a majority homopolymer indel may mean the
/// consensus itself is wrong, so it is not downgraded to gray (`--gray-hp-min`
/// does not apply to it). It also still forms a finding.
#[test]
fn qc_gray_weak_majority_stays_yellow() {
    const L: usize = 3_000;
    const HP: usize = 1_600; // run of 8 `A`
    const HP_LEN: usize = 8;

    let mut r = random_ref(L, 909);
    r[HP - 1] = b'C';
    r[HP + HP_LEN] = b'G';
    for x in &mut r[HP..HP + HP_LEN] {
        *x = b'A';
    }

    // 6 reads at 6x over 1000..2000; 4 carry one homopolymer deletion each
    // (different `A`, all normalise to one allele) -> AF 4/6 = 0.667, strand
    // balanced (2 fwd / 2 rev): the majority gate fires but Wilson LB < 0.5, so
    // the tier is yellow `weak_majority`.
    let mut reads = Vec::new();
    for i in 0..6usize {
        let cigar = if i < 4 {
            let off = i % HP_LEN;
            let m1 = (HP + off - 1_000) as u32;
            vec![Cigar::Match(m1), Cigar::Del(1), Cigar::Match(1_000 - m1 - 1)]
        } else {
            vec![Cigar::Match(1_000)]
        };
        reads.push(Rd { pos: 1_000, seq: read_seq(&r, 1_000, &cigar, &[]), cigar, rev: i % 2 == 1 });
    }
    let (vcf, cands, _bed, json) = run_case_p("gray_weakmajority", &r, reads, params(1, None));

    let v = vcf.iter().find(|v| v[1] == HP.to_string()).unwrap_or_else(|| panic!("no DEL at {HP}: {vcf:?}"));
    assert!(v[7].contains("TYPE=DEL") && v[7].contains("HP=8"), "{v:?}");
    assert!(v[7].contains("TIER=YELLOW") && v[7].contains("REASON=weak_majority") && !v[7].contains("HPGRAY"), "{v:?}");
    let c = cands
        .iter()
        .find(|c| c["type"] == "DEL" && c["pos0"] == (HP - 1).to_string())
        .unwrap_or_else(|| panic!("no weak_majority DEL candidate: {cands:?}"));
    assert_eq!((c["tier"].as_str(), c["reason"].as_str()), ("yellow", "weak_majority"), "{c:?}");
    assert_eq!(json["alerts"]["gray"], 0, "{json}");
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/qc-scan/test");
    let fbed = std::fs::read_to_string(base.join("gray_weakmajority/out.findings.bed")).unwrap();
    assert!(fbed.lines().filter(|l| !l.starts_with('#')).any(|l| {
        let f: Vec<&str> = l.split('\t').collect();
        f[1].parse::<u64>().unwrap() <= HP as u64 && (HP as u64) < f[2].parse::<u64>().unwrap()
    }), "weak_majority HP indel lost its finding:\n{fbed}");
}

/// A minority **insertion** in a homopolymer is gray too (not just deletions):
/// AF ~0.25, HP >= 5.
#[test]
fn qc_gray_minority_homopolymer_ins() {
    const HP: usize = 2_001;
    const HP_LEN: usize = 8;
    let mut r = random_ref(4_000, 111);
    r[HP - 1] = b'C';
    r[HP + HP_LEN] = b'G';
    for x in &mut r[HP..HP + HP_LEN] {
        *x = b'A';
    }
    // 1/4 of reads insert one `A` in the run (normalises to anchor HP-1).
    let reads = hp_reads(&r, |pair, _rev, p| {
        if pair % 4 == 0 && p == HP + pair % HP_LEN { Edit::InsA } else { Edit::None }
    });
    let qp = QcParams { unlinked_action: cramhook::qc::UnlinkedAction::Flag, ..params(1, None) };
    let (vcf, cands, _bed, json) = run_case_p("gray_minority_ins", &r, reads, qp);

    let v = vcf.iter().find(|v| v[1] == HP.to_string()).unwrap_or_else(|| panic!("no INS at {HP}: {vcf:?}"));
    assert!(v[7].contains("TYPE=INS") && v[7].contains("HP=8"), "{v:?}");
    assert!(v[7].contains("TIER=GRAY") && v[7].contains(";HPGRAY") && v[7].contains("REASON=minority"), "{v:?}");
    let c = cands
        .iter()
        .find(|c| c["type"] == "INS" && c["pos0"] == (HP - 1).to_string())
        .unwrap_or_else(|| panic!("no INS candidate: {cands:?}"));
    assert_eq!((c["tier"].as_str(), c["reason"].as_str()), ("gray", "minority"), "{c:?}");
    assert_eq!(json["alerts"]["by_type"]["INS"]["gray"], 1, "{json}");
}

/// The `--gray-hp-min` boundary: with the default 5 a run of 4 stays yellow and
/// a run of 5 is gray; `--gray-hp-min 4` makes both gray, `--gray-hp-min 6`
/// makes neither gray.
#[test]
fn qc_gray_hp_min_boundary_and_custom() {
    const RUN4: usize = 1_201; // homopolymer length 4
    const LEN4: usize = 4;
    const RUN5: usize = 3_001; // homopolymer length 5
    const LEN5: usize = 5;
    let mut r = random_ref(4_000, 222);
    for (hp, len) in [(RUN4, LEN4), (RUN5, LEN5)] {
        r[hp - 1] = b'C';
        r[hp + len] = b'G';
        for x in &mut r[hp..hp + len] {
            *x = b'A';
        }
    }
    let reads = hp_reads(&r, |pair, _rev, p| {
        if pair % 4 != 0 {
            return Edit::None;
        }
        if p == RUN4 + pair % LEN4 {
            Edit::Del
        } else if p == RUN5 + pair % LEN5 {
            Edit::Del
        } else {
            Edit::None
        }
    });
    let base = |t, c| QcParams { unlinked_action: cramhook::qc::UnlinkedAction::Flag, ..params(t, c) };
    let info = |vcf: &[Vec<String>], pos: usize| -> String {
        vcf.iter().find(|v| v[1] == pos.to_string()).map(|v| v[7].clone()).unwrap_or_default()
    };

    // Default `--gray-hp-min 5`: HP4 yellow, HP5 gray.
    let (v, _c, _b, _j) = run_case_p("gray_hp_def", &r, reads.clone(), base(1, None));
    assert!(info(&v, RUN4).contains("TIER=YELLOW") && info(&v, RUN4).contains("HP=4") && !info(&v, RUN4).contains("HPGRAY"), "{v:?}");
    assert!(info(&v, RUN5).contains("TIER=GRAY") && info(&v, RUN5).contains("HP=5"), "{v:?}");

    // `--gray-hp-min 4`: both gray.
    let (v, _c, _b, _j) = run_case_p("gray_hp_4", &r, reads.clone(), QcParams { gray_hp_min: 4, ..base(1, None) });
    assert!(info(&v, RUN4).contains("TIER=GRAY") && info(&v, RUN4).contains("HPGRAY"), "{v:?}");
    assert!(info(&v, RUN5).contains("TIER=GRAY"), "{v:?}");

    // `--gray-hp-min 6`: neither gray.
    let (v, _c, _b, _j) = run_case_p("gray_hp_6", &r, reads.clone(), QcParams { gray_hp_min: 6, ..base(1, None) });
    assert!(info(&v, RUN4).contains("TIER=YELLOW") && !info(&v, RUN4).contains("HPGRAY"), "{v:?}");
    assert!(info(&v, RUN5).contains("TIER=YELLOW") && !info(&v, RUN5).contains("HPGRAY"), "{v:?}");
}

/// A strand-biased yellow HP indel (AF >= 0.40, alt on one strand only) is gray.
#[test]
fn qc_gray_strand_bias_indel() {
    const HP: usize = 1_601;
    const HP_LEN: usize = 8;
    let mut r = random_ref(3_000, 333);
    r[HP - 1] = b'C';
    r[HP + HP_LEN] = b'G';
    for x in &mut r[HP..HP + HP_LEN] {
        *x = b'A';
    }
    // Every forward read carries the deletion (AF ~0.5) and no reverse read
    // does: the strand gate fails -> yellow `strand_bias` (not `weak_majority`).
    let reads = hp_reads(&r, |pair, rev, p| {
        if !rev && p == HP + pair % HP_LEN { Edit::Del } else { Edit::None }
    });
    let (vcf, cands, _bed, json) = run_case_p("gray_strand_bias", &r, reads, params(1, None));

    let v = vcf.iter().find(|v| v[1] == HP.to_string()).unwrap_or_else(|| panic!("no DEL at {HP}: {vcf:?}"));
    assert!(v[7].contains("TYPE=DEL") && v[7].contains("TIER=GRAY") && v[7].contains(";STRAND"), "{v:?}");
    assert!(v[7].contains("REASON=strand_bias") && v[7].contains(";HPGRAY"), "{v:?}");
    let c = cands
        .iter()
        .find(|c| c["type"] == "DEL" && c["pos0"] == (HP - 1).to_string())
        .unwrap_or_else(|| panic!("no strand_bias DEL candidate: {cands:?}"));
    assert_eq!((c["tier"].as_str(), c["reason"].as_str()), ("gray", "strand_bias"), "{c:?}");
    assert_eq!(json["alerts"]["gray"], 1, "{json}");
}

/// A gray call below the absolute low-coverage depth is FILTER `LowCov` and is
/// counted under `alerts.lowcov.gray`, not `alerts.gray`.
#[test]
fn qc_gray_lowcov_counts() {
    const HP: usize = 2_001;
    const HP_LEN: usize = 8;
    let mut r = random_ref(4_000, 555);
    r[HP - 1] = b'C';
    r[HP + HP_LEN] = b'G';
    for x in &mut r[HP..HP + HP_LEN] {
        *x = b'A';
    }
    let reads = hp_reads(&r, |pair, _rev, p| {
        if pair % 4 == 0 && p == HP + pair % HP_LEN { Edit::Del } else { Edit::None }
    });
    // A high absolute threshold makes the ~40x gray call low-coverage by `abs`.
    let qp = QcParams { lowcov_abs: Some(100), unlinked_action: cramhook::qc::UnlinkedAction::Flag, ..params(1, None) };
    let (vcf, _c, _bed, json) = run_case_p("gray_lowcov", &r, reads, qp);

    let v = vcf.iter().find(|v| v[1] == HP.to_string()).unwrap_or_else(|| panic!("no DEL at {HP}: {vcf:?}"));
    assert_eq!(v[6], "LowCov", "{v:?}");
    assert!(v[7].contains(";LOWCOV") && v[7].contains("TIER=GRAY") && v[7].contains(";HPGRAY"), "{v:?}");
    assert_eq!(json["alerts"]["gray"], 0, "{json}");
    assert_eq!(json["alerts"]["lowcov"]["gray"], 1, "{json}");
    assert_eq!(json["alerts"]["by_type"]["DEL"]["lowcov"]["gray"], 1, "{json}");
}

/// A gray call inside another finding's span is not a member: the finding's
/// counts are unchanged. With gray off the same site is yellow and joins,
/// extending the finding by one `DEL` member.
#[test]
fn qc_gray_not_absorbed_by_finding() {
    const RED: usize = 1_211; // majority HP deletion -> red DEL finding
    const GRAY: usize = 1_261; // minority HP deletion, 50 bp away
    const HP_LEN: usize = 8;
    let mut r = random_ref(4_000, 444);
    for hp in [RED, GRAY] {
        r[hp - 1] = b'C';
        r[hp + HP_LEN] = b'G';
        for x in &mut r[hp..hp + HP_LEN] {
            *x = b'A';
        }
    }
    let reads = hp_reads(&r, |pair, _rev, p| {
        if p == RED + pair % HP_LEN {
            Edit::Del
        } else if pair % 4 == 0 && p == GRAY + pair % HP_LEN {
            Edit::Del
        } else {
            Edit::None
        }
    });
    let base = |t, c| QcParams { unlinked_action: cramhook::qc::UnlinkedAction::Flag, ..params(t, c) };
    let (_v, _c, _b, _j) = run_case_p("gray_find_on", &r, reads.clone(), base(1, None));
    let (_v, _c, _b, _j) = run_case_p("gray_find_off", &r, reads.clone(), QcParams { gray_homopolymer: false, ..base(1, None) });
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scratch/qc-scan/test");
    let f_on = std::fs::read_to_string(dir.join("gray_find_on/out.findings.bed")).unwrap();
    let f_off = std::fs::read_to_string(dir.join("gray_find_off/out.findings.bed")).unwrap();
    // The finding that carries a `DEL` member, as (start, end, extra).
    let del_finding = |fbed: &str| -> Option<(u64, u64, String)> {
        fbed.lines().filter(|l| !l.starts_with('#')).find_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let extra = f.get(7).copied().unwrap_or("");
            (extra.contains("members=") && extra.contains("DEL:"))
                .then(|| (f[1].parse().unwrap(), f[2].parse().unwrap(), extra.to_string()))
        })
    };
    // Gray on: the gray site is not a member and the finding stays DEL:1.
    assert!(!f_on.contains("DEL:2"), "gray site joined a finding:\n{f_on}");
    let (s_on, e_on, x_on) = del_finding(&f_on).unwrap_or_else(|| panic!("no DEL finding:\n{f_on}"));
    assert!(x_on.contains("members=DEL:1") && e_on < GRAY as u64, "{s_on} {e_on} {x_on}");
    // Gray off: the yellow site is inside the finding's span and joins it.
    let (s_off, e_off, x_off) = del_finding(&f_off).unwrap_or_else(|| panic!("no DEL finding:\n{f_off}"));
    assert!(x_off.contains("members=DEL:2"), "{x_off}");
    assert!(s_off <= GRAY as u64 && (GRAY as u64) < e_off, "gray site not inside the finding span {s_off} {e_off}");
    assert_eq!(s_on, s_off, "finding start moved: {s_on} vs {s_off}");
}
