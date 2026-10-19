//! Integration tests for `cramhook translate` (M3).
//!
//! For every synthetic scenario and both read models we run the translator and
//! check the read-count invariant, the S1 identity property, the S6/S6b repeat
//! safety property and the agreement with `oracle.cram`. Missing fixtures are a
//! hard error.
//!
//! Generate the data with `pixi run sim && pixi run crams-synthetic`.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{self, Format, Header, Read, Record, Writer};
use serde::Deserialize;

use cramhook::{TranslateConfig, translate};

const SIDS: &[&str] = &[
    "S1", "S2", "S3", "S4", "S5", "S6", "S6b", "S7", "S8", "S9", "S10a", "S10b", "S11", "S12",
    "S13", "S13b", "S14", "S14b", "S15", "S16", "S17", "S17b", "S18", "S19", "S20", "S21",
];
const MODELS: &[&str] = &["hifi", "ont"];

/// Random 3-5 edit combinations (seed indices SR1..SR5). Kept out of the default
/// run; execute with `cargo test --release -- --ignored` (see
/// `oracle_agreement_random_scenarios`).
const SIDS_RANDOM: &[&str] = &["SR1", "SR2", "SR3", "SR4", "SR5"];

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

#[derive(Deserialize)]
struct Truth {
    repeats_old: Vec<Repeat>,
    contigs_old: BTreeMap<String, u64>,
}

#[derive(Deserialize)]
struct Repeat {
    name: String,
    contig: String,
    start: u64,
    end: u64,
}

/// A record distilled for comparison.
#[derive(Clone, Debug, PartialEq, Eq)]
struct RV {
    qname: Vec<u8>,
    flag: u16,
    contig: String,
    pos: i64,
    mapq: u8,
    cigar: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
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
    /// The five fields compared against the oracle.
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
        // MD is a decode artifact; XL is the translation marker; s2 is the
        // old-assembly competitor score stripped from translated records (M10b).
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
        seq: rec.seq().as_bytes(),
        qual: rec.qual().to_vec(),
        aux,
    }
}

fn read_bam(path: &Path) -> Vec<RV> {
    let mut reader =
        bam::Reader::from_path(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let names: Vec<String> = reader
        .header()
        .target_names()
        .iter()
        .map(|n| String::from_utf8_lossy(n).to_string())
        .collect();
    let mut out = Vec::new();
    for r in reader.records() {
        let r = r.unwrap();
        out.push(distill(&r, &names));
    }
    out
}

fn read_cram(path: &Path, reference: &Path) -> Vec<RV> {
    let mut reader =
        bam::Reader::from_path(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    reader.set_reference(reference).unwrap();
    let names: Vec<String> = reader
        .header()
        .target_names()
        .iter()
        .map(|n| String::from_utf8_lossy(n).to_string())
        .collect();
    let mut out = Vec::new();
    for r in reader.records() {
        let r = r.unwrap();
        out.push(distill(&r, &names));
    }
    out
}

fn by_name(recs: Vec<RV>) -> HashMap<Vec<u8>, Vec<RV>> {
    let mut m: HashMap<Vec<u8>, Vec<RV>> = HashMap::new();
    for r in recs {
        m.entry(r.qname.clone()).or_default().push(r);
    }
    m
}

/// Raw BAM/CRAM records (no distillation), for tag-level assertions.
fn read_raw(path: &Path, reference: &Path) -> Vec<Record> {
    let mut reader =
        bam::Reader::from_path(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    reader.set_reference(reference).unwrap();
    reader.records().map(|r| r.unwrap()).collect()
}

fn default_cfg(xl: bool) -> TranslateConfig {
    TranslateConfig { xl_tag: xl, ..Default::default() }
}

fn out_prefix(tag: &str, sid: &str, model: &str) -> PathBuf {
    root().join("target/test-translate").join(tag).join(sid).join(model)
}

fn translated_bam(tag: &str, sid: &str, model: &str) -> PathBuf {
    let mut s = out_prefix(tag, sid, model).into_os_string();
    s.push(".translated.bam");
    PathBuf::from(s)
}

fn assets(sid: &str, model: &str) -> (PathBuf, PathBuf, PathBuf) {
    let dir = root().join("data/synthetic").join(sid);
    (dir.join(model).join("old.cram"), dir.join("old.fa"), dir.join("new.fa"))
}

fn setup() {
    let base = root().join("data/synthetic");
    assert!(
        base.exists(),
        "data/synthetic is missing: run `pixi run sim` and `pixi run crams-synthetic`"
    );
}

/// Run translate for one scenario/model and return (stats, translated records).
fn run_case(tag: &str, sid: &str, model: &str) -> (cramhook::Stats, Vec<RV>) {
    let (cram, old_fa, new_fa) = assets(sid, model);
    assert!(
        cram.exists() && old_fa.exists() && new_fa.exists(),
        "{sid}/{model}: missing fixtures: run `pixi run sim` and `pixi run crams-synthetic`"
    );
    let stats = translate(&cram, &old_fa, &new_fa, &out_prefix(tag, sid, model), &default_cfg(true), 4)
        .unwrap_or_else(|e| panic!("{sid}/{model}: translate failed: {e:#}"));
    let translated = read_bam(&translated_bam(tag, sid, model));
    (stats, translated)
}

fn load_truth(sid: &str) -> Truth {
    let p = root().join("data/synthetic").join(sid).join("truth.json");
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

/// Old-coordinate window of a translated record (Plan 4.2 rule 7).
fn old_window(rec: &RV, margin: u64, contig_len: u64) -> (u64, u64) {
    let (_q, r, lead, trail) = cigar_spans(&rec.cigar);
    let ref_end = rec.pos + r as i64;
    let s = (rec.pos - lead as i64 - margin as i64).max(0) as u64;
    let e = (ref_end + trail as i64 + margin as i64).min(contig_len as i64).max(s as i64) as u64;
    (s, e)
}

fn cigar_spans(cigar: &str) -> (u64, u64, u64, u64) {
    let mut q = 0u64;
    let mut r = 0u64;
    let mut ops: Vec<(char, u64)> = Vec::new();
    let mut num = String::new();
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            let n: u64 = num.parse().unwrap_or(0);
            num.clear();
            ops.push((c, n));
        }
    }
    for (op, n) in &ops {
        match op {
            'M' | '=' | 'X' => {
                q += n;
                r += n;
            }
            'I' => q += n,
            'D' | 'N' => r += n,
            _ => {}
        }
    }
    let lead = ops.first().filter(|(o, _)| *o == 'S').map(|(_, n)| *n).unwrap_or(0);
    let trail = ops.last().filter(|(o, _)| *o == 'S').map(|(_, n)| *n).unwrap_or(0);
    (q, r, lead, trail)
}

// ---------------------------------------------------------------------------
// Test 1: correctness properties (must pass)

#[test]
fn invariants_s1_and_repeat_windows() {
    setup();
    for sid in SIDS {
        let truth = load_truth(sid);
        for model in MODELS {
            let (stats, translated) = run_case("invariants", sid, model);

            // Read-count invariant: translated reads + deferred reads == primary in.
            assert_eq!(
                stats.translated_reads + stats.deferred_reads,
                stats.records_in.primary,
                "{sid}/{model}: read-count invariant"
            );
            // Each read name appears once: the output has no duplicate qnames.
            let mut seen = HashSet::new();
            for r in &translated {
                if r.primary() {
                    assert!(seen.insert(r.qname.clone()), "{sid}/{model}: duplicate read");
                }
            }

            // S1: identical assemblies, no diff-level deferral, byte-identical
            // to the input apart from XL (and the stripped MD/s2 tags).
            if *sid == "S1" {
                assert!(
                    !stats.deferred_by_class.contains_key("diff"),
                    "S1/{model}: diff-level deferrals {:?}",
                    stats.deferred_by_class
                );
                assert!(
                    stats.origin_crossing_translated > 0,
                    "S1/{model}: no origin-crossing split read translated"
                );
                check_s1_identical(model, &translated);
            }

            // S2/S3/S4 must translate something (rotation / '-' / rename).
            if matches!(*sid, "S2" | "S3" | "S4") {
                assert!(stats.translated_reads > 0, "{sid}/{model}: translated nothing");
            }

            // S6/S6b: no translated window may overlap a chromosomal IS copy.
            if matches!(*sid, "S6" | "S6b") {
                let is_copies: Vec<&Repeat> =
                    truth.repeats_old.iter().filter(|r| r.name.starts_with("IS")).collect();
                for rec in translated.iter().filter(|r| r.primary()) {
                    let Some(&len) = truth.contigs_old.get(&rec.contig) else { continue };
                    let (ws, we) = old_window(rec, 500, len);
                    for rep in &is_copies {
                        if rep.contig == rec.contig {
                            assert!(
                                !(ws < rep.end && rep.start < we),
                                "{sid}/{model}: translated {:?} window {ws}..{we} overlaps {} {}..{}",
                                String::from_utf8_lossy(&rec.qname),
                                rep.name,
                                rep.start,
                                rep.end
                            );
                        }
                    }
                }
            }
        }
    }
}

/// On S1 every translated record must equal its input record, apart from XL
/// (and the stripped MD/s2 tags).
fn check_s1_identical(model: &str, translated: &[RV]) {
    let base = root().join("data/synthetic").join("S1");
    let input = read_cram(&base.join(model).join("old.cram"), &base.join("old.fa"));
    let input_by = by_name(input);
    for rec in translated {
        let candidates = input_by
            .get(&rec.qname)
            .unwrap_or_else(|| panic!("S1/{model}: translated read not in input"));
        let matched = candidates.iter().any(|c| c == rec);
        assert!(
            matched,
            "S1/{model}: translated record for {:?} is not identical to any input record\n  ours:  {rec:?}\n  input: {candidates:?}",
            String::from_utf8_lossy(&rec.qname)
        );
    }
}

/// M10b: a translated record carries `XL:A:T` and neither the decode `MD`
/// artifact nor the old-assembly competitor score `s2`. The input still has
/// `s2`, so the check is non-vacuous. Remapped (`XL:A:R`) records keep `s2`
/// (covered by the lift tests).
#[test]
fn translated_records_tagged_xl_and_drop_s2() {
    setup();
    let sid = "S1";
    let model = "hifi";
    let (cram, old_fa, _new_fa) = assets(sid, model);
    let input = read_raw(&cram, &old_fa);
    assert!(
        input.iter().any(|r| r.aux(b"s2").is_ok()),
        "{sid}/{model}: input has no s2 tag; the strip check would be vacuous"
    );
    assert!(
        input.iter().any(|r| r.aux(b"MD").is_ok()),
        "{sid}/{model}: input has no MD tag; the strip check would be vacuous"
    );

    run_case("tags", sid, model);
    let out = translated_bam("tags", sid, model);
    let mut reader = bam::Reader::from_path(&out)
        .unwrap_or_else(|e| panic!("open {}: {e}", out.display()));
    let mut n = 0u64;
    for rec in reader.records() {
        let rec = rec.unwrap();
        n += 1;
        assert!(rec.aux(b"s2").is_err(), "translated record keeps s2: {rec:?}");
        assert!(rec.aux(b"MD").is_err(), "translated record keeps MD: {rec:?}");
        match rec.aux(b"XL") {
            Ok(Aux::Char(b'T')) => {}
            other => panic!("translated record XL is not T: {other:?}"),
        }
    }
    assert!(n > 0, "{sid}/{model}: no translated records produced");
}

// ---------------------------------------------------------------------------
// Test 2: oracle agreement (Plan acceptance: >= 99.9%), with equivalence

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Identical,
    Equivalent,
    Mismatch,
}

fn parse_ops(s: &str) -> Option<Vec<Cigar>> {
    CigarString::try_from(s).ok().map(|c| c.0)
}

/// Query bases consumed by all ops (including soft clips).
fn cq_len(ops: &[Cigar]) -> u64 {
    use Cigar::*;
    ops.iter()
        .map(|c| match c {
            Match(l) | Ins(l) | SoftClip(l) | Equal(l) | Diff(l) => *l as u64,
            _ => 0,
        })
        .sum()
}

/// Reference bases consumed by all ops.
fn cr_len(ops: &[Cigar]) -> u64 {
    use Cigar::*;
    ops.iter()
        .map(|c| match c {
            Match(l) | Del(l) | RefSkip(l) | Equal(l) | Diff(l) => *l as u64,
            _ => 0,
        })
        .sum()
}

/// Identical, equivalent, or a real mismatch.
///
/// Equivalent = CIGAR differs, but contig, pos, END (reference end), strand,
/// MAPQ, read length (query length incl. soft clips) and `AS` are all equal.
/// `NM` and `de` are NOT part of the definition; they are reported as metrics
/// for equivalent pairs instead.
fn kind(a: &RV, b: &RV) -> Kind {
    if a.alignment_key() == b.alignment_key() {
        return Kind::Identical;
    }
    if a.contig != b.contig || a.pos != b.pos || a.reverse() != b.reverse() || a.mapq != b.mapq {
        return Kind::Mismatch;
    }
    // Same reference end (pos + ref span) and query length (incl. soft clips).
    let (Some(oa), Some(ob)) = (parse_ops(&a.cigar), parse_ops(&b.cigar)) else {
        return Kind::Mismatch;
    };
    let a_end = a.pos + cr_len(&oa) as i64;
    let b_end = b.pos + cr_len(&ob) as i64;
    if a_end != b_end || cq_len(&oa) != cq_len(&ob) {
        return Kind::Mismatch;
    }
    // Equal AS (must be present on both; it is minimap2's score for the same
    // read/ref content, carried over unchanged on our side).
    let (Some(aas), Some(bas)) = (a.aux.get("AS"), b.aux.get("AS")) else {
        return Kind::Mismatch;
    };
    if aas != bas {
        return Kind::Mismatch;
    }
    Kind::Equivalent
}

/// Whether the `tag` values of two records are equal (missing == missing).
fn aux_eq(a: &RV, b: &RV, tag: &str) -> bool {
    a.aux.get(tag) == b.aux.get(tag)
}

#[derive(Default)]
struct Counts {
    identical: usize,
    equivalent: usize,
    /// Equivalences found on records translated from a `+` block (a bug).
    equivalent_plus: usize,
    mismatch: usize,
    /// NM/de agreement among equivalent pairs (metrics, not part of the gate).
    nm_equal: usize,
    nm_differ: usize,
    de_equal: usize,
    de_differ: usize,
}

impl Counts {
    fn total(&self) -> usize {
        self.identical + self.equivalent + self.equivalent_plus + self.mismatch
    }
    fn agreeing(&self) -> usize {
        self.identical + self.equivalent + self.equivalent_plus
    }
    /// NM/de metrics for the equivalent pairs (flipped and '+' block).
    fn tally_equiv_tags(&mut self, a: &RV, b: &RV) {
        if aux_eq(a, b, "NM") {
            self.nm_equal += 1;
        } else {
            self.nm_differ += 1;
        }
        if aux_eq(a, b, "de") {
            self.de_equal += 1;
        } else {
            self.de_differ += 1;
        }
    }
}

#[test]
fn oracle_agreement_all_scenarios() {
    setup();
    let mut table = Vec::new();
    let mut mismatches = Vec::new();
    let mut failures = Vec::new();

    for sid in SIDS {
        for model in MODELS {
            let (stats, translated) = run_case("oracle", sid, model);
            let (cram, old_fa, new_fa) = assets(sid, model);
            let oracle_path = cram.parent().unwrap().join("oracle.cram");
            let oracle = read_cram(&oracle_path, &new_fa);
            let input = read_cram(&cram, &old_fa);
            let oracle_by = by_name(oracle);
            let ours_by = by_name(translated);
            let input_by = by_name(input);

            // Sanity: AS is carried over unchanged from old.cram.
            for (name, recs) in &ours_by {
                let irecs = input_by.get(name);
                for r in recs {
                    let i = irecs.and_then(|v| v.iter().find(|i| i.supplementary() == r.supplementary()));
                    match i {
                        Some(i) if r.aux.get("AS") == i.aux.get("AS") => {}
                        Some(i) => failures.push(format!(
                            "{sid}/{model} {}: AS changed by translation: ours {:?} input {:?}",
                            String::from_utf8_lossy(name),
                            r.aux.get("AS"),
                            i.aux.get("AS")
                        )),
                        None => failures.push(format!(
                            "{sid}/{model} {}: no input record for AS sanity",
                            String::from_utf8_lossy(name)
                        )),
                    }
                }
            }

            let mut pc = Counts::default();
            let mut sc = Counts::default();
            compare(
                &ours_by,
                &oracle_by,
                &input_by,
                &mut pc,
                &mut sc,
                &mut mismatches,
                sid,
                model,
            );

            let pfrac = if pc.total() == 0 { 1.0 } else { pc.agreeing() as f64 / pc.total() as f64 };
            let sfrac = if sc.total() == 0 { 1.0 } else { sc.agreeing() as f64 / sc.total() as f64 };
            table.push(format!(
                "{sid:5} {model:5} xl {:5} defer(read {:4} diff {:4})  primary ident {:4} equiv {:3} mism {:3} => {pfrac:.5}  \
                 (nm diff {} de diff {} of {} equiv; equiv on '+' blocks {})  supp {}/{} = {sfrac:.4}",
                stats.translated_reads,
                *stats.deferred_by_class.get("read").unwrap_or(&0),
                *stats.deferred_by_class.get("diff").unwrap_or(&0),
                pc.identical,
                pc.equivalent,
                pc.mismatch,
                pc.nm_differ,
                pc.de_differ,
                pc.equivalent + pc.equivalent_plus,
                pc.equivalent_plus,
                sc.agreeing(),
                sc.total(),
            ));

            if pc.total() > 0 && pfrac < 0.999 {
                failures.push(format!(
                    "{sid}/{model}: only {pfrac:.5} of {} translated primary records are identical+equivalent (need >= 0.999)",
                    pc.total()
                ));
            }
            if pc.equivalent_plus > 0 {
                failures.push(format!(
                    "{sid}/{model}: {} primary records are merely equivalent on a '+' block translation (bug)",
                    pc.equivalent_plus
                ));
            }
            if sc.total() > 0 && sfrac < 0.999 {
                failures.push(format!(
                    "{sid}/{model}: only {sfrac:.5} of {} translated supplementary records are identical+equivalent",
                    sc.total()
                ));
            }
        }
    }

    println!("\n=== translate per-scenario oracle table ===");
    for row in &table {
        println!("{row}");
    }
    if !mismatches.is_empty() {
        println!("\n=== mismatching reads (scenario/model name | ours | oracle) ===");
        for m in &mismatches {
            println!("{m}");
        }
    }
    println!("\n{} scenario/model runs below 99.9% or with '+' block equivalence", failures.len());
    for f in &failures {
        println!("  {f}");
    }
    assert!(
        failures.is_empty(),
        "oracle agreement failures in {} scenario/model runs (see table and mismatches above)",
        failures.len()
    );
}

/// Translate one scenario/model and return (agreement fraction, primary total,
/// mismatches, first mismatch lines). Used by the ignored random / expected-fail
/// runners below.
fn agreement(sid: &str, model: &str, tag: &str) -> (f64, usize, usize, Vec<String>) {
    let (_stats, translated) = run_case(tag, sid, model);
    let (cram, old_fa, new_fa) = assets(sid, model);
    let oracle = read_cram(&cram.parent().unwrap().join("oracle.cram"), &new_fa);
    let input = read_cram(&cram, &old_fa);
    let mut pc = Counts::default();
    let mut sc = Counts::default();
    let mut mismatches = Vec::new();
    compare(
        &by_name(translated),
        &by_name(oracle),
        &by_name(input),
        &mut pc,
        &mut sc,
        &mut mismatches,
        sid,
        model,
    );
    let frac = if pc.total() == 0 { 1.0 } else { pc.agreeing() as f64 / pc.total() as f64 };
    (frac, pc.total(), pc.mismatch, mismatches)
}

/// Random edit combinations (`SR1..SR5`). Ignored by default; run explicitly
/// with `cargo test --release -- --ignored oracle_agreement_random_scenarios`.
#[test]
#[ignore = "random combinations; run with --ignored"]
fn oracle_agreement_random_scenarios() {
    setup();
    let mut fails = Vec::new();
    for sid in SIDS_RANDOM {
        for model in MODELS {
            let (frac, total, mismatch, mismatches) = agreement(sid, model, "oracle_rand");
            println!("{sid}/{model}: agreement={frac:.5} n={total} mismatch={mismatch}");
            for m in mismatches.iter().take(5) {
                println!("    {m}");
            }
            if total > 0 && frac < 0.999 {
                fails.push(format!("{sid}/{model}: agreement={frac:.5} (< 0.999)"));
            }
        }
    }
    assert!(fails.is_empty(), "random-scenario oracle agreement failures: {fails:#?}");
}

fn compare(
    ours: &HashMap<Vec<u8>, Vec<RV>>,
    oracle: &HashMap<Vec<u8>, Vec<RV>>,
    input: &HashMap<Vec<u8>, Vec<RV>>,
    pc: &mut Counts,
    sc: &mut Counts,
    mismatches: &mut Vec<String>,
    sid: &str,
    model: &str,
) {
    let empty: Vec<RV> = Vec::new();
    for (name, recs) in ours {
        let orecs = oracle.get(name).unwrap_or(&empty);
        let irecs = input.get(name).unwrap_or(&empty);
        for r in recs.iter().filter(|r| r.primary()) {
            let op = orecs.iter().find(|o| o.primary());
            match op {
                Some(o) => match kind(r, o) {
                    Kind::Identical => pc.identical += 1,
                    Kind::Equivalent => {
                        pc.tally_equiv_tags(r, o);
                        let flipped = irecs
                            .iter()
                            .find(|i| i.primary())
                            .map(|i| i.reverse() != r.reverse())
                            .unwrap_or(false);
                        if flipped {
                            pc.equivalent += 1;
                        } else {
                            pc.equivalent_plus += 1;
                        }
                    }
                    Kind::Mismatch => {
                        pc.mismatch += 1;
                        mismatches.push(mismatch_line(sid, model, name, r, o, ""));
                    }
                },
                None => {
                    pc.mismatch += 1;
                    mismatches.push(format!(
                        "{sid}/{model} {} | ours {} {}:{}-{} | oracle: no primary",
                        String::from_utf8_lossy(name),
                        r.contig,
                        r.pos,
                        if r.reverse() { "-" } else { "+" },
                        r.cigar
                    ));
                }
            }
        }
        for r in recs.iter().filter(|r| r.supplementary()) {
            let op = orecs.iter().find(|o| o.supplementary());
            match op {
                Some(o) => match kind(r, o) {
                    Kind::Identical => sc.identical += 1,
                    Kind::Equivalent => sc.equivalent += 1,
                    Kind::Mismatch => {
                        sc.mismatch += 1;
                        mismatches.push(mismatch_line(sid, model, name, r, o, " (supp)"));
                    }
                },
                None => sc.mismatch += 1,
            }
        }
    }
}

fn mismatch_line(sid: &str, model: &str, name: &[u8], r: &RV, o: &RV, tag: &str) -> String {
    format!(
        "{sid}/{model} {}{tag} | ours {} {}:{} {} {} | oracle {} {}:{} {} {}",
        String::from_utf8_lossy(name),
        r.contig,
        r.pos,
        if r.reverse() { "-" } else { "+" },
        r.cigar,
        r.mapq,
        o.contig,
        o.pos,
        if o.reverse() { "-" } else { "+" },
        o.cigar,
        o.mapq
    )
}

// ---------------------------------------------------------------------------
// Unit-ish tests for the equivalence helper and the split ordering fixes

fn rv(cigar: &str, aux: &[(&str, &str)]) -> RV {
    RV {
        qname: b"r".to_vec(),
        flag: 0,
        contig: "chr".into(),
        pos: 100,
        mapq: 60,
        cigar: cigar.into(),
        seq: Vec::new(),
        qual: Vec::new(),
        aux: aux.iter().map(|(k, v)| (k.to_string(), v.to_string())).collect(),
    }
}

#[test]
fn kind_requires_present_as() {
    let full = [("AS", "i:10"), ("NM", "i:2")];
    let a = rv("10M", &full);
    let b = rv("5=5X", &full);
    assert_eq!(kind(&a, &b), Kind::Equivalent);
    // AS must be present on both
    assert_eq!(kind(&a, &rv("5=5X", &[("NM", "i:2")])), Kind::Mismatch);
    assert_eq!(kind(&rv("10M", &[("NM", "i:2")]), &b), Kind::Mismatch);
    assert_eq!(kind(&a, &rv("5=5X", &[])), Kind::Mismatch);
    // different AS -> mismatch
    assert_eq!(kind(&a, &rv("5=5X", &[("AS", "i:11"), ("NM", "i:2")])), Kind::Mismatch);
    // NM and de are not part of the definition -> still equivalent
    assert_eq!(kind(&a, &rv("5=5X", &[("AS", "i:10"), ("NM", "i:3")])), Kind::Equivalent);
    assert_eq!(
        kind(&a, &rv("5=5X", &[("AS", "i:10"), ("NM", "i:2"), ("de", "i:1")])),
        Kind::Equivalent
    );
    // different reference end or read length -> mismatch
    assert_eq!(kind(&a, &rv("11M", &[("AS", "i:10")])), Kind::Mismatch);
    assert_eq!(kind(&a, &rv("9M1D", &[("AS", "i:10")])), Kind::Mismatch);
}

struct RecSpec {
    flags: u16,
    pos: i64,
    sa: Option<&'static str>,
}

fn pseudo_random_seq(n: usize) -> String {
    let mut x: u64 = 0x1234_5678_9abc_def0;
    (0..n)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[((x >> 33) & 3) as usize] as char
        })
        .collect()
}

fn write_spec_bam(path: &Path, specs: &[RecSpec]) {
    let mut header = Header::new();
    header.push_record(
        HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"),
    );
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", 5000));
    header.push_record(
        HeaderRecord::new(b"PG")
            .push_tag(b"ID", "minimap2")
            .push_tag(b"PN", "minimap2")
            .push_tag(b"VN", "2.31")
            .push_tag(b"CL", "minimap2 -ax map-hifi --secondary=no -Y old.fa reads.fq"),
    );
    let mut w = Writer::from_path(path, &header, Format::Bam).unwrap();
    let cigar = CigarString(vec![Cigar::Match(100)]);
    let seq = vec![b'A'; 100];
    let qual = vec![30u8; 100];
    for s in specs {
        let mut r = Record::new();
        r.set(b"r", Some(&cigar), &seq, &qual);
        r.set_tid(0);
        r.set_pos(s.pos);
        r.set_flags(s.flags);
        r.set_mapq(60);
        r.push_aux(b"NM", Aux::I32(0)).unwrap();
        r.push_aux(b"AS", Aux::I32(100)).unwrap();
        r.push_aux(b"de", Aux::Float(0.001)).unwrap();
        if let Some(sa) = s.sa {
            r.push_aux(b"SA", Aux::String(sa)).unwrap();
        }
        w.write(&r).unwrap();
    }
}

fn run_spec_case(tag: &str, specs: &[RecSpec]) -> cramhook::Stats {
    let dir = root().join("target/test-translate").join(tag);
    std::fs::create_dir_all(&dir).unwrap();
    let seq = pseudo_random_seq(5000);
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    std::fs::write(&old_fa, format!(">chr\n{seq}\n")).unwrap();
    std::fs::write(&new_fa, format!(">chr\n{seq}\n")).unwrap();
    let bam = dir.join("old.bam");
    write_spec_bam(&bam, specs);
    translate(&bam, &old_fa, &new_fa, &dir.join("out"), &default_cfg(false), 2)
        .unwrap_or_else(|e| panic!("{tag}: translate failed: {e:#}"))
}

/// Build a BAM with one primary record carrying a duplicated aux tag `dup`
/// (minimap2 never writes duplicates; invalid SAM).
fn write_dup_tag_bam(path: &Path, dup: &[u8]) {
    let mut header = Header::new();
    header.push_record(
        HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"),
    );
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", 5000));
    header.push_record(
        HeaderRecord::new(b"PG")
            .push_tag(b"ID", "minimap2")
            .push_tag(b"PN", "minimap2")
            .push_tag(b"VN", "2.31")
            .push_tag(b"CL", "minimap2 -ax map-hifi --secondary=no -Y old.fa reads.fq"),
    );
    let mut w = Writer::from_path(path, &header, Format::Bam).unwrap();
    let cigar = CigarString(vec![Cigar::Match(100)]);
    let seq = vec![b'A'; 100];
    let qual = vec![30u8; 100];
    let mut r = Record::new();
    r.set(b"r", Some(&cigar), &seq, &qual);
    r.set_flags(0);
    r.set_tid(0);
    r.set_pos(0);
    r.set_mapq(60);
    r.push_aux(b"NM", Aux::I32(0)).unwrap();
    r.push_aux(b"AS", Aux::I32(100)).unwrap();
    r.push_aux(b"de", Aux::Float(0.001)).unwrap();
    r.push_aux(b"s1", Aux::I32(100)).unwrap();
    r.push_aux(b"s2", Aux::I32(0)).unwrap();
    if dup == b"MD" {
        r.push_aux(b"MD", Aux::String("100")).unwrap();
        r.push_aux_unchecked(b"MD", Aux::String("100")).unwrap();
    } else {
        r.push_aux_unchecked(dup, Aux::I32(7)).unwrap();
    }
    w.write(&r).unwrap();
}

/// M10b: a record with a duplicated aux tag is invalid SAM and must be deferred
/// as `tags`, never translated (both a duplicated `s2` and a duplicated `MD`).
#[test]
fn duplicate_aux_tag_defers_as_tags() {
    let dups: [&[u8]; 2] = [b"s2", b"MD"];
    for dup in dups {
        let tag = String::from_utf8_lossy(dup).to_string();
        let dir = root().join("target/test-translate").join(format!("dup_tag_{tag}"));
        std::fs::create_dir_all(&dir).unwrap();
        let seq = pseudo_random_seq(5000);
        let old_fa = dir.join("old.fa");
        let new_fa = dir.join("new.fa");
        std::fs::write(&old_fa, format!(">chr\n{seq}\n")).unwrap();
        std::fs::write(&new_fa, format!(">chr\n{seq}\n")).unwrap();
        let bam = dir.join("old.bam");
        write_dup_tag_bam(&bam, dup);
        let stats = translate(&bam, &old_fa, &new_fa, &dir.join("out"), &default_cfg(false), 2)
            .unwrap_or_else(|e| panic!("duplicate {tag}: translate failed: {e:#}"));
        assert_eq!(stats.translated_reads, 0, "duplicate {tag}: must not translate");
        assert_eq!(
            *stats.deferred_by_reason.get("tags").unwrap_or(&0),
            1,
            "duplicate {tag}: expected tags deferral, got {:?}",
            stats.deferred_by_reason
        );
    }
}

/// A supplementary-first true origin candidate followed by a primary without
/// `SA` must be deferred as inconsistent, not translated.
#[test]
fn ordering_supplementary_candidate_then_primary_without_sa() {
    let stats = run_spec_case(
        "order1",
        &[
            // candidate: reaches the contig end, valid single SA entry
            RecSpec { flags: 0x800, pos: 4900, sa: Some("chr,1,+,100M,60,0;") },
            // primary with no SA
            RecSpec { flags: 0, pos: 0, sa: None },
        ],
    );
    assert_eq!(stats.translated_reads, 0, "the primary must not be translated");
    assert_eq!(
        *stats.deferred_by_reason.get("split_sa_inconsistent").unwrap_or(&0),
        1,
        "expected split_sa_inconsistent, got {:?}",
        stats.deferred_by_reason
    );
}

/// A supplementary-first NON-origin read is decided on the first record and not
/// buffered; the primary is deferred with the same reason.
#[test]
fn ordering_supplementary_non_origin_not_buffered() {
    let stats = run_spec_case(
        "order2",
        &[
            // reaches neither contig boundary -> not an origin candidate
            RecSpec { flags: 0x800, pos: 500, sa: Some("chr,1,+,100M,60,0;") },
            RecSpec { flags: 0, pos: 0, sa: Some("chr,501,+,100M,60,0;") },
        ],
    );
    assert_eq!(stats.translated_reads, 0);
    assert_eq!(
        *stats.deferred_by_reason.get("split_not_origin").unwrap_or(&0),
        1,
        "expected split_not_origin, got {:?}",
        stats.deferred_by_reason
    );
}
