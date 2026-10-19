//! Round-trip tests for a lifted input (`A→B→A`, `A→B→A→B`).
//!
//! For every deterministic synthetic scenario and both read models we lift
//! `old.cram` (A) to the new assembly (B), back to A, and to B again, all through
//! the real CLI (so the chained `@PG`s are genuine). Each step is compared to a
//! full remap of the reads onto that step's target using the same equivalence
//! criteria as `translate_synthetic` (`kind`: identical / equivalent / mismatch,
//! `AS` required on both, `NM`/`de` reported as metrics). The `→A` oracle is the
//! original `old.cram`: minimap2 is deterministic (verified: rerunning
//! `minimap2 -ax map-hifi --secondary=no -Y old.fa reads` reproduces `old.cram`
//! 100 % identical), so a full remap onto A is exactly the original mapping.
//!
//! A mismatch fails the test with a diagnosis; the test is not loosened to hide
//! one. The random `SR1..SR5` scenarios are in the ignored
//! `roundtrip_random_scenarios`.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::process::Command;

use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{self, Read};

/// Deterministic scenarios (the same set `translate_synthetic` runs).
const SIDS: &[&str] = &[
    "S1", "S2", "S3", "S4", "S5", "S6", "S6b", "S7", "S8", "S9", "S10a", "S10b", "S11", "S12",
    "S13", "S13b", "S14", "S14b", "S15", "S16", "S17", "S17b", "S18", "S19", "S20", "S21",
];
const MODELS: &[&str] = &["hifi", "ont"];
/// Random edit combinations; kept out of the default run.
const SIDS_RANDOM: &[&str] = &["SR1", "SR2", "SR3", "SR4", "SR5"];

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn setup() {
    let base = root().join("data/synthetic");
    assert!(
        base.exists(),
        "data/synthetic is missing: run `pixi run sim` and `pixi run crams-synthetic`"
    );
}

// ---------------------------------------------------------------------------
// record view + equivalence criteria (mirrors tests/translate_synthetic.rs)

#[derive(Clone, Debug, PartialEq, Eq)]
struct RV {
    qname: Vec<u8>,
    flag: u16,
    contig: String,
    pos: i64,
    mapq: u8,
    cigar: String,
    aux: BTreeMap<String, String>,
}

impl RV {
    fn primary(&self) -> bool {
        self.flag & (0x100 | 0x800) == 0
    }
    fn supplementary(&self) -> bool {
        self.flag & 0x800 != 0
    }
    fn reverse(&self) -> bool {
        self.flag & 0x10 != 0
    }
    fn alignment_key(&self) -> (String, i64, bool, String, u8) {
        (self.contig.clone(), self.pos, self.reverse(), self.cigar.clone(), self.mapq)
    }
}

fn aux_to_string(a: &Aux<'_>) -> String {
    match a {
        Aux::Char(v) => format!("A:{}", *v as char),
        Aux::I8(v) => format!("i:{v}"),
        Aux::U8(v) => format!("i:{v}"),
        Aux::I16(v) => format!("i:{v}"),
        Aux::U16(v) => format!("i:{v}"),
        Aux::I32(v) => format!("i:{v}"),
        Aux::U32(v) => format!("i:{v}"),
        Aux::Float(v) => format!("f:{v}"),
        Aux::Double(v) => format!("f:{v}"),
        Aux::String(s) => format!("Z:{s}"),
        other => format!("{other:?}"),
    }
}

fn distill(rec: &bam::Record, names: &[String]) -> RV {
    let tid = rec.tid();
    let contig = if tid >= 0 {
        names.get(tid as usize).cloned().unwrap_or_default()
    } else {
        String::new()
    };
    let mut aux = BTreeMap::new();
    for item in rec.aux_iter() {
        let (tag, val) = item.expect("aux");
        let tag = String::from_utf8_lossy(tag).to_string();
        // MD is a decode artifact, XL the translation marker, s2 stripped from
        // translated records (M10b).
        if tag == "MD" || tag == "XL" || tag == "s2" {
            continue;
        }
        aux.insert(tag, aux_to_string(&val));
    }
    RV {
        qname: rec.qname().to_vec(),
        flag: rec.flags(),
        contig,
        pos: rec.pos(),
        mapq: rec.mapq(),
        cigar: rec.cigar().to_string(),
        aux,
    }
}

fn read_cram(path: &Path, reference: &Path) -> Vec<RV> {
    let mut reader =
        bam::Reader::from_path(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    reader.set_reference(reference).unwrap();
    let names: Vec<String> =
        reader.header().target_names().iter().map(|n| String::from_utf8_lossy(n).to_string()).collect();
    reader.records().map(|r| distill(&r.unwrap(), &names)).collect()
}

fn by_name(recs: Vec<RV>) -> HashMap<Vec<u8>, Vec<RV>> {
    let mut m: HashMap<Vec<u8>, Vec<RV>> = HashMap::new();
    for r in recs {
        m.entry(r.qname.clone()).or_default().push(r);
    }
    m
}

fn parse_ops(s: &str) -> Option<Vec<Cigar>> {
    CigarString::try_from(s).ok().map(|c| c.0)
}

fn cq_len(ops: &[Cigar]) -> u64 {
    use Cigar::*;
    ops.iter()
        .map(|c| match c {
            Match(l) | Ins(l) | SoftClip(l) | Equal(l) | Diff(l) => *l as u64,
            _ => 0,
        })
        .sum()
}

fn cr_len(ops: &[Cigar]) -> u64 {
    use Cigar::*;
    ops.iter()
        .map(|c| match c {
            Match(l) | Del(l) | RefSkip(l) | Equal(l) | Diff(l) => *l as u64,
            _ => 0,
        })
        .sum()
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Identical,
    Equivalent,
    Mismatch,
}

/// Identical, equivalent, or a real mismatch (same definition as
/// `translate_synthetic::kind`: `AS` required on both; `NM`/`de` not part of it).
fn kind(a: &RV, b: &RV) -> Kind {
    if a.alignment_key() == b.alignment_key() {
        return Kind::Identical;
    }
    if a.contig != b.contig || a.pos != b.pos || a.reverse() != b.reverse() || a.mapq != b.mapq {
        return Kind::Mismatch;
    }
    let (Some(oa), Some(ob)) = (parse_ops(&a.cigar), parse_ops(&b.cigar)) else {
        return Kind::Mismatch;
    };
    if a.pos + cr_len(&oa) as i64 != b.pos + cr_len(&ob) as i64 || cq_len(&oa) != cq_len(&ob) {
        return Kind::Mismatch;
    }
    let (Some(aas), Some(bas)) = (a.aux.get("AS"), b.aux.get("AS")) else {
        return Kind::Mismatch;
    };
    if aas != bas {
        return Kind::Mismatch;
    }
    Kind::Equivalent
}

#[derive(Default, Clone, Copy)]
struct Agree {
    identical: usize,
    equivalent: usize,
    mismatch: usize,
    nm_differ: usize,
    de_differ: usize,
}

impl Agree {
    fn total(&self) -> usize {
        self.identical + self.equivalent + self.mismatch
    }
    fn agreeing(&self) -> usize {
        self.identical + self.equivalent
    }
    fn fraction(&self) -> f64 {
        if self.total() == 0 {
            1.0
        } else {
            self.agreeing() as f64 / self.total() as f64
        }
    }
}

fn tally(k: Kind, a: &RV, b: &RV, out: &mut Agree, name: &[u8], tag: &str, diag: &mut Vec<String>) {
    match k {
        Kind::Identical => out.identical += 1,
        Kind::Equivalent => {
            out.equivalent += 1;
            if a.aux.get("NM") != b.aux.get("NM") {
                out.nm_differ += 1;
            }
            if a.aux.get("de") != b.aux.get("de") {
                out.de_differ += 1;
            }
        }
        Kind::Mismatch => {
            out.mismatch += 1;
            diag.push(format!(
                "  {}{tag} | ours {} {}:{} {} {} | oracle {} {}:{} {} {}",
                String::from_utf8_lossy(name),
                a.contig,
                a.pos,
                if a.reverse() { "-" } else { "+" },
                a.cigar,
                a.mapq,
                b.contig,
                b.pos,
                if b.reverse() { "-" } else { "+" },
                b.cigar,
                b.mapq
            ));
        }
    }
}

/// Compare one lifted CRAM to its full-remap oracle (primary + supplementary).
fn compare(ours_path: &Path, ours_ref: &Path, oracle_path: &Path, oracle_ref: &Path) -> (Agree, Vec<String>) {
    let ours = by_name(read_cram(ours_path, ours_ref));
    let oracle = by_name(read_cram(oracle_path, oracle_ref));
    let mut agree = Agree::default();
    let mut diag = Vec::new();
    for (name, recs) in &ours {
        let empty: Vec<RV> = Vec::new();
        let orecs = oracle.get(name).unwrap_or(&empty);
        // Primary.
        for r in recs.iter().filter(|r| r.primary()) {
            match orecs.iter().find(|o| o.primary()) {
                Some(o) => tally(kind(r, o), r, o, &mut agree, name, "", &mut diag),
                None => {
                    agree.mismatch += 1;
                    diag.push(format!(
                        "  {} | ours {} {}:{} {} | oracle: no primary",
                        String::from_utf8_lossy(name),
                        r.contig,
                        r.pos,
                        if r.reverse() { "-" } else { "+" },
                        r.cigar
                    ));
                }
            }
        }
        // Supplementary: each of ours must match a distinct oracle record.
        let mut used = vec![false; orecs.len()];
        for r in recs.iter().filter(|r| r.supplementary()) {
            let mut matched = None;
            for (i, o) in orecs.iter().enumerate() {
                if !used[i] && o.supplementary() && kind(r, o) != Kind::Mismatch {
                    used[i] = true;
                    matched = Some(o);
                    break;
                }
            }
            match matched {
                Some(o) => tally(kind(r, o), r, o, &mut agree, name, " (supp)", &mut diag),
                None => {
                    agree.mismatch += 1;
                    diag.push(format!(
                        "  {} (supp) | ours {} {}:{} {} | oracle: no match",
                        String::from_utf8_lossy(name),
                        r.contig,
                        r.pos,
                        if r.reverse() { "-" } else { "+" },
                        r.cigar
                    ));
                }
            }
        }
        for (i, o) in orecs.iter().enumerate() {
            if o.supplementary() && !used[i] {
                agree.mismatch += 1;
                diag.push(format!(
                    "  {} (supp) | ours: missing | oracle {} {}:{} {}",
                    String::from_utf8_lossy(name),
                    o.contig,
                    o.pos,
                    if o.reverse() { "-" } else { "+" },
                    o.cigar
                ));
            }
        }
    }
    // Reads only in the oracle.
    for (name, orecs) in &oracle {
        if !ours.contains_key(name) {
            agree.mismatch += orecs.len();
            diag.push(format!("  {} | ours: missing", String::from_utf8_lossy(name)));
        }
    }
    (agree, diag)
}

// ---------------------------------------------------------------------------
// CLI driver

struct LiftStats {
    translated: u64,
    primary: u64,
}

/// Run `cramhook lift` through the real CLI (genuine chained `@PG`s) and read
/// the stats. Asserts the contract held.
fn cli_lift(dir: &Path, tag: &str, input: &Path, old_fa: &Path, new_fa: &Path) -> (PathBuf, LiftStats) {
    let out = dir.join(format!("{tag}.cram"));
    let stats = dir.join(format!("{tag}.json"));
    let o = Command::new(env!("CARGO_BIN_EXE_cramhook"))
        .args(["lift", "--in"])
        .arg(input)
        .args(["--old-ref"])
        .arg(old_fa)
        .args(["--new-ref"])
        .arg(new_fa)
        .args(["--out"])
        .arg(&out)
        .args(["--threads", "4", "--min-input-mb", "0", "--stats"])
        .arg(&stats)
        .output()
        .expect("running cramhook lift");
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(o.status.success(), "{tag}: lift failed: {err}");
    let j: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&stats).unwrap()).expect("stats json");
    assert_eq!(j["contract_ok"], serde_json::Value::Bool(true), "{tag}: contract failed: {err}");
    (
        out,
        LiftStats {
            translated: j["translated_reads"].as_u64().unwrap(),
            primary: j["records_in"]["primary"].as_u64().unwrap(),
        },
    )
}

struct Row {
    label: String,
    s1: LiftStats,
    a1: Agree,
    s2: LiftStats,
    a2: Agree,
    s3: LiftStats,
    a3: Agree,
}

fn run_scenario(sid: &str, model: &str, diag: &mut Vec<String>) -> Row {
    let d = root().join("data/synthetic").join(sid);
    let old_fa = d.join("old.fa");
    let new_fa = d.join("new.fa");
    let old_cram = d.join(model).join("old.cram");
    let oracle_b = d.join(model).join("oracle.cram");
    assert!(
        old_cram.exists() && oracle_b.exists() && old_fa.exists() && new_fa.exists(),
        "{sid}/{model}: missing fixtures: run `pixi run sim` and `pixi run crams-synthetic`"
    );
    let out_dir = root().join("target/test-lift-roundtrip").join(sid).join(model);
    std::fs::create_dir_all(&out_dir).unwrap();

    let (ab, s1) = cli_lift(&out_dir, "ab", &old_cram, &old_fa, &new_fa);
    let (a1, d1) = compare(&ab, &new_fa, &oracle_b, &new_fa);
    push_diag(diag, sid, model, "A→B", d1);

    let (aba, s2) = cli_lift(&out_dir, "aba", &ab, &new_fa, &old_fa);
    // →A: the original mapping is the full-remap oracle (minimap2 is deterministic).
    let (a2, d2) = compare(&aba, &old_fa, &old_cram, &old_fa);
    push_diag(diag, sid, model, "B→A", d2);

    let (abab, s3) = cli_lift(&out_dir, "abab", &aba, &old_fa, &new_fa);
    let (a3, d3) = compare(&abab, &new_fa, &oracle_b, &new_fa);
    push_diag(diag, sid, model, "A→B→A→B", d3);

    Row {
        label: format!("{sid:5} {model:5}"),
        s1,
        a1,
        s2,
        a2,
        s3,
        a3,
    }
}

fn push_diag(diag: &mut Vec<String>, sid: &str, model: &str, step: &str, ds: Vec<String>) {
    if !ds.is_empty() {
        diag.push(format!("{sid}/{model} {step}: {} mismatches", ds.len()));
        diag.extend(ds.into_iter().take(20));
    }
}

fn frac(s: &LiftStats) -> f64 {
    if s.primary == 0 {
        1.0
    } else {
        s.translated as f64 / s.primary as f64
    }
}

fn run_all(sids: &[&str]) -> (Vec<Row>, Vec<String>) {
    // One thread per read model, so the 2 (model) x N (scenario) passes overlap.
    std::thread::scope(|scope| {
        let handles: Vec<_> = MODELS
            .iter()
            .map(|model| {
                scope.spawn(move || {
                    let mut rows = Vec::new();
                    let mut diag = Vec::new();
                    for sid in sids {
                        rows.push(run_scenario(sid, model, &mut diag));
                    }
                    (rows, diag)
                })
            })
            .collect();
        let mut rows = Vec::new();
        let mut diag = Vec::new();
        for h in handles {
            let (r, d) = h.join().expect("round-trip worker panicked");
            rows.extend(r);
            diag.extend(d);
        }
        rows.sort_by(|a, b| a.label.cmp(&b.label));
        (rows, diag)
    })
}

fn report_and_assert(rows: &[Row], diag: &[String], what: &str) {
    println!("\n=== {what}: round-trip table ===");
    println!("scenario model  xl-frac A→B / B→A / A→B→A→B | agreement A→B / B→A / A→B→A→B (ident/equiv/mism)");
    let mut bad = 0;
    for r in rows {
        if r.a1.mismatch + r.a2.mismatch + r.a3.mismatch > 0 {
            bad += 1;
        }
        println!(
            "{}  {:.3} {:.3} {:.3} | {:.5}({}/{}/{}) {:.5}({}/{}/{}) {:.5}({}/{}/{})",
            r.label,
            frac(&r.s1),
            frac(&r.s2),
            frac(&r.s3),
            r.a1.fraction(),
            r.a1.identical,
            r.a1.equivalent,
            r.a1.mismatch,
            r.a2.fraction(),
            r.a2.identical,
            r.a2.equivalent,
            r.a2.mismatch,
            r.a3.fraction(),
            r.a3.identical,
            r.a3.equivalent,
            r.a3.mismatch,
        );
    }
    let nm: usize = rows.iter().map(|r| r.a1.nm_differ + r.a2.nm_differ + r.a3.nm_differ).sum();
    let de: usize = rows.iter().map(|r| r.a1.de_differ + r.a2.de_differ + r.a3.de_differ).sum();
    println!("equivalent-pair tag drift: NM {} / de {}", nm, de);
    if !diag.is_empty() {
        println!("\n=== {what}: mismatches ===");
        for d in diag.iter().take(60) {
            println!("{d}");
        }
    }
    assert!(diag.is_empty(), "{what}: {} mismatch diagnoses (see above)", diag.len());
    assert_eq!(bad, 0);
}

#[test]
fn roundtrip_deterministic_scenarios() {
    setup();
    let (rows, diag) = run_all(SIDS);
    report_and_assert(&rows, &diag, "A→B→A / A→B→A→B");
}

#[test]
#[ignore = "random combinations; run with --ignored"]
fn roundtrip_random_scenarios() {
    setup();
    let (rows, diag) = run_all(SIDS_RANDOM);
    report_and_assert(&rows, &diag, "A→B→A / A→B→A→B (random)");
}
