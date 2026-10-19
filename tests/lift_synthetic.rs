//! Integration tests for `cramhook lift` (M4).
//!
//! Runs the full pipeline on a few synthetic scenarios and checks the output is
//! a valid sorted/indexed CRAM with a correct header, read accounting and `XL`
//! tags, plus the contract-violation fallback and failure paths (reader thread
//! error, injected child exit codes) using env-override wrappers.
//!
//! Regenerate fixtures with `pixi run sim && pixi run crams-synthetic`.

use std::collections::HashSet;
use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::sync::{Mutex, MutexGuard};
use std::thread;
use std::time::Duration;

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{self, Format, Header, Read, Record, Writer};

use cramhook::contract::read_fai;
use cramhook::{LiftOptions, TranslateConfig, lift};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Env overrides are process-global, so all lift tests run one at a time.
fn lock() -> MutexGuard<'static, ()> {
    static LOCK: Mutex<()> = Mutex::new(());
    LOCK.lock().unwrap_or_else(|e| e.into_inner())
}

fn setup() {
    let base = root().join("data/synthetic");
    assert!(
        base.exists(),
        "data/synthetic is missing: run `pixi run sim` and `pixi run crams-synthetic`"
    );
}

fn default_cfg(xl: bool) -> TranslateConfig {
    TranslateConfig { xl_tag: xl, ..Default::default() }
}

struct Case {
    out: PathBuf,
    stats: cramhook::Stats,
}

fn run_lift(sid: &str, model: &str, xl: bool) -> Case {
    let dir = root().join("data/synthetic").join(sid);
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let cram = dir.join(model).join("old.cram");
    assert!(
        cram.exists() && old_fa.exists() && new_fa.exists(),
        "{sid}/{model}: missing fixtures: run `pixi run crams-synthetic`"
    );
    let out_dir = root().join("target/test-lift").join(sid).join(model);
    std::fs::create_dir_all(&out_dir).unwrap();
    let out = out_dir.join("out.cram");
    let stats_path = out_dir.join("stats.json");
    let opts = LiftOptions { strict: false, keep_tmp: false, min_input_mb: Some(0.0) };
    let stats = lift(
        &cram,
        &old_fa,
        &new_fa,
        &out,
        Some(&stats_path),
        4,
        &default_cfg(xl),
        &opts,
        None,
    )
    .unwrap_or_else(|e| panic!("{sid}/{model}: lift failed: {e:#}"));
    Case { out, stats }
}

/// Read the output CRAM and check header/order/tags/counts.
fn check_output(case: &Case, sid: &str, model: &str) {
    let dir = root().join("data/synthetic").join(sid);
    let new_fa = dir.join("new.fa");
    let new_fai = read_fai(&new_fa.with_extension("fa.fai")).unwrap();
    let expected_names: Vec<String> = new_fai.iter().map(|(n, _)| n.clone()).collect();

    let mut idx = case.out.as_os_str().to_owned();
    idx.push(".crai");
    assert!(Path::new(&idx).exists(), "{sid}/{model}: missing .crai index");

    let mut reader = bam::Reader::from_path(&case.out).unwrap();
    reader.set_reference(&new_fa).unwrap();

    let names: Vec<String> = reader
        .header()
        .target_names()
        .iter()
        .map(|n| String::from_utf8_lossy(n).to_string())
        .collect();
    assert_eq!(names, expected_names, "{sid}/{model}: @SQ differs from new.fa .fai");

    let mut seen: HashSet<Vec<u8>> = HashSet::new();
    let mut primaries = 0u64;
    let mut xl_t = 0u64;
    let mut xl_r = 0u64;
    let mut records = 0u64;
    let mut last: Option<(i32, i64)> = None;
    for r in reader.records() {
        let r = r.unwrap();
        records += 1;
        let key = (r.tid(), r.pos());
        if r.tid() >= 0 {
            if let Some(prev) = last {
                assert!(prev <= key, "{sid}/{model}: output not coordinate-sorted at {r:?} ({prev:?} > {key:?})");
            }
            last = Some(key);
        }
        let xl = match r.aux(b"XL") {
            Ok(Aux::Char(c)) => c,
            other => panic!("{sid}/{model}: missing/invalid XL tag: {other:?}"),
        };
        match xl {
            b'T' => {
                xl_t += 1;
                // M10b: translated records drop the old-assembly competitor
                // score s2. (MD is also stripped on write, but htslib regenerates
                // it when the output CRAM is decoded, so it is not asserted here.)
                assert!(
                    r.aux(b"s2").is_err(),
                    "{sid}/{model}: translated record keeps s2: {r:?}"
                );
            }
            b'R' => {
                xl_r += 1;
                // Remapped records keep minimap2's fresh s2 (mapped records).
                if r.tid() >= 0 {
                    assert!(
                        r.aux(b"s2").is_ok(),
                        "{sid}/{model}: remapped record lost s2: {r:?}"
                    );
                }
            }
            other => panic!("{sid}/{model}: unexpected XL value {other}"),
        }
        if r.flags() & (0x100 | 0x800) == 0 {
            primaries += 1;
            assert!(seen.insert(r.qname().to_vec()), "{sid}/{model}: duplicate primary read");
        }
    }

    assert_eq!(primaries, case.stats.records_in.primary, "{sid}/{model}: primary count");
    assert_eq!(primaries, case.stats.out_primary_records, "{sid}/{model}: primary invariant");
    assert_eq!(xl_t, case.stats.translated_records, "{sid}/{model}: XL:A:T count");
    assert_eq!(xl_r, records - case.stats.translated_records, "{sid}/{model}: XL:A:R count");
}

#[test]
fn lift_s1_hifi() {
    let _g = lock();
    setup();
    let c = run_lift("S1", "hifi", true);
    check_output(&c, "S1", "hifi");
    assert_eq!(c.stats.translated_reads + c.stats.remapped_reads, c.stats.records_in.primary);
}

#[test]
fn lift_s3_ont() {
    let _g = lock();
    setup();
    let c = run_lift("S3", "ont", true);
    check_output(&c, "S3", "ont");
}

#[test]
fn lift_s12_ont() {
    let _g = lock();
    setup();
    let c = run_lift("S12", "ont", true);
    check_output(&c, "S12", "ont");
}

/// M10b: `--no-xl-tag` (`TranslateConfig::xl_tag = false`) writes no `XL` tag on
/// any translated or remapped record, and the effective value is in stats.
#[test]
fn lift_no_xl_tag_omits_xl() {
    let _g = lock();
    setup();
    let case = run_lift("S1", "hifi", false);
    assert!(!case.stats.params.xl_tag, "params.xl_tag must record the effective value");
    let new_fa = root().join("data/synthetic/S1/new.fa");
    let mut reader = bam::Reader::from_path(&case.out).unwrap();
    reader.set_reference(&new_fa).unwrap();
    let mut n = 0u64;
    for rec in reader.records() {
        let rec = rec.unwrap();
        n += 1;
        assert!(rec.aux(b"XL").is_err(), "XL present despite xl_tag=false: {rec:?}");
    }
    assert!(n > 0, "no output records produced");
}

/// M9: with default settings a small input takes the tiny-input fallback: the
/// diff is skipped, every read is remapped, the stats record the `tiny_input`
/// reason, and `contract_ok` stays true. The output must still be a valid
/// sorted/indexed CRAM with each read exactly once as primary.
#[test]
fn lift_tiny_input_fallback_default() {
    let _g = lock();
    setup();
    let sid = "S3";
    let model = "ont"; // lr:hq threshold 20 MB > this CRAM's size
    let dir = root().join("data/synthetic").join(sid);
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let cram = dir.join(model).join("old.cram");
    let out_dir = root().join("target/test-lift").join(sid).join("ont-tiny");
    std::fs::create_dir_all(&out_dir).unwrap();
    let out = out_dir.join("out.cram");
    let stats_path = out_dir.join("stats.json");
    // min_input_mb: None -> per-preset default (here lr:hq = 20 MB).
    // keep_tmp -> also check the deferred TSV records the `tiny_input` reason.
    let opts = LiftOptions { strict: false, keep_tmp: true, min_input_mb: None };
    let stats = lift(
        &cram,
        &old_fa,
        &new_fa,
        &out,
        Some(&stats_path),
        4,
        &default_cfg(true),
        &opts,
        None,
    )
    .unwrap_or_else(|e| panic!("{sid}/{model}: tiny-input lift failed: {e:#}"));

    assert!(stats.contract_ok, "tiny-input fallback must keep contract_ok=true");
    assert!(
        stats.contract_reason.is_none(),
        "tiny-input fallback must keep contract_reason null"
    );
    let fb = stats.fallback.clone().expect("tiny-input fallback must be recorded");
    assert_eq!(fb.reason, "tiny_input");
    assert_eq!(fb.input_bytes, std::fs::metadata(&cram).unwrap().len());
    assert!(fb.input_bytes < fb.threshold_bytes, "input not below threshold: {fb:?}");
    assert_eq!(fb.threshold_bytes, 20_000_000, "lr:hq default threshold");
    assert_eq!(stats.params.min_input_mb, 20.0, "effective threshold in params");
    assert_eq!(stats.translated_reads, 0, "tiny-input fallback must translate nothing");
    assert_eq!(stats.remapped_reads, stats.records_in.primary);

    // Reason histogram: every primary deferred as `tiny_input`, none as `contract`.
    assert_eq!(
        *stats.deferred_by_reason.get("tiny_input").unwrap_or(&0),
        stats.records_in.primary,
        "all primaries must be deferred as tiny_input: {:?}",
        stats.deferred_by_reason
    );
    assert!(
        !stats.deferred_by_reason.contains_key("contract"),
        "tiny-input fallback must not tag reads as contract: {:?}",
        stats.deferred_by_reason
    );
    assert_eq!(*stats.deferred_by_class.get("fallback").unwrap_or(&0), stats.records_in.primary);

    let text = std::fs::read_to_string(&stats_path).unwrap();
    assert!(text.contains("\"tiny_input\""), "stats.json missing fallback reason: {text}");

    // The deferred TSV written by --keep-tmp must carry the `tiny_input` reason.
    let tsv = out_dir.join("out.cram.tmp").join("deferred.deferred.tsv");
    let tsv_text = std::fs::read_to_string(&tsv)
        .unwrap_or_else(|e| panic!("reading deferred TSV {}: {e}", tsv.display()));
    assert!(tsv_text.contains("tiny_input"), "deferred TSV missing tiny_input: {tsv_text}");
    assert_eq!(
        tsv_text.lines().filter(|l| l.contains("\ttiny_input\t")).count() as u64,
        stats.records_in.primary,
        "deferred TSV row count"
    );

    // Output validity: sorted/indexed, header, each primary exactly once.
    check_output(&Case { out, stats }, sid, model);
}

/// M9: tiny-input fallback must keep the record-level contract checks. A
/// secondary record aborts in `--strict` mode and preserves the previous output;
/// without `--strict` the contract reason wins over `tiny_input`.
#[test]
fn lift_tiny_input_secondary_strict_errors() {
    let _g = lock();
    setup();
    let d = out_dir("tiny_secondary_strict");
    let (old_fa, new_fa, bam) = write_contract_input(&d, true);
    let out = d.join("out.cram");
    std::fs::write(&out, b"PREVIOUS").unwrap();
    let out_crai = out.with_extension("cram.crai");
    std::fs::write(&out_crai, b"OLDCRAI").unwrap();

    let mut cfg = default_cfg(false);
    cfg.strict = true;
    let opts = LiftOptions { strict: true, keep_tmp: false, min_input_mb: None };
    let e = format!(
        "{:#}",
        lift(&bam, &old_fa, &new_fa, &out, None, 4, &cfg, &opts, None)
            .expect_err("strict tiny-input with a secondary record must error")
    );
    assert!(e.contains("contract violation"), "error: {e}");
    assert!(e.contains("mid-stream") && e.contains("secondary"), "record-level check not hit: {e}");
    assert_eq!(std::fs::read(&out).unwrap(), b"PREVIOUS", "previous CRAM clobbered");
    assert_eq!(std::fs::read(&out_crai).unwrap(), b"OLDCRAI", "previous index clobbered");
    let (staged, staged_crai) = cramhook::lift::staged_paths(&out);
    assert!(!staged.exists(), "staged cram remains");
    assert!(!staged_crai.exists(), "staged crai remains");
}

/// M9: without `--strict`, a tiny-input input with a secondary record reports
/// `contract_ok=false` and the contract reason; the fallback metadata must not
/// claim `tiny_input` and reads defer as `contract`.
#[test]
fn lift_tiny_input_secondary_nonstrict_contract_wins() {
    let _g = lock();
    setup();
    let d = out_dir("tiny_secondary_nonstrict");
    let (old_fa, new_fa, bam) = write_contract_input(&d, true);
    let out = d.join("out.cram");
    let stats = lift(
        &bam,
        &old_fa,
        &new_fa,
        &out,
        None,
        4,
        &default_cfg(false),
        &LiftOptions { strict: false, keep_tmp: false, min_input_mb: None },
        None,
    )
    .unwrap_or_else(|e| panic!("tiny-input contract fallback failed: {e:#}"));
    assert!(!stats.contract_ok, "contract violation must set contract_ok=false");
    assert!(
        stats.contract_reason.as_deref().unwrap_or("").contains("secondary"),
        "contract reason: {:?}",
        stats.contract_reason
    );
    assert!(stats.fallback.is_none(), "contract reason must win over tiny_input");
    assert_eq!(stats.translated_reads, 0);
    assert_eq!(stats.remapped_reads, stats.records_in.primary);
    assert!(!stats.deferred_by_reason.contains_key("tiny_input"), "{:?}", stats.deferred_by_reason);
    assert_eq!(*stats.deferred_by_reason.get("contract").unwrap_or(&0), stats.records_in.primary);
}

/// M9: tiny-input fallback must catch a hard-clipped primary (a record-level
/// contract violation) too. Strict aborts; non-strict reports the contract.
#[test]
fn lift_tiny_input_hardclip_strict_errors() {
    let _g = lock();
    setup();
    let d = out_dir("tiny_hardclip_strict");
    let (old_fa, new_fa, bam) = write_contract_input(&d, false);
    let out = d.join("out.cram");
    std::fs::write(&out, b"PREVIOUS").unwrap();
    let mut cfg = default_cfg(false);
    cfg.strict = true;
    let e = format!(
        "{:#}",
        lift(
            &bam,
            &old_fa,
            &new_fa,
            &out,
            None,
            4,
            &cfg,
            &LiftOptions { strict: true, keep_tmp: false, min_input_mb: None },
            None,
        )
        .expect_err("strict tiny-input with a hard-clipped primary must error")
    );
    assert!(e.contains("contract violation"), "error: {e}");
    assert!(e.contains("hard-clipped"), "record-level check not hit: {e}");
    assert_eq!(std::fs::read(&out).unwrap(), b"PREVIOUS", "previous CRAM clobbered");
}

/// M9: non-strict hard-clipped primary under the tiny-input fallback.
#[test]
fn lift_tiny_input_hardclip_nonstrict_contract_wins() {
    let _g = lock();
    setup();
    let d = out_dir("tiny_hardclip_nonstrict");
    let (old_fa, new_fa, bam) = write_contract_input(&d, false);
    let out = d.join("out.cram");
    let stats = lift(
        &bam,
        &old_fa,
        &new_fa,
        &out,
        None,
        4,
        &default_cfg(false),
        &LiftOptions { strict: false, keep_tmp: false, min_input_mb: None },
        None,
    )
    .unwrap_or_else(|e| panic!("tiny-input hard-clip fallback failed: {e:#}"));
    assert!(!stats.contract_ok);
    assert!(
        stats.contract_reason.as_deref().unwrap_or("").contains("hard-clipped"),
        "contract reason: {:?}",
        stats.contract_reason
    );
    assert!(stats.fallback.is_none(), "contract reason must win over tiny_input");
    assert_eq!(stats.translated_reads, 0);
}

// ---------------------------------------------------------------------------
// failure paths

fn wrappers_dir() -> PathBuf {
    let d = root().join("target/test-lift/wrappers");
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write_script(name: &str, body: &str) -> PathBuf {
    let p = wrappers_dir().join(name);
    std::fs::write(&p, body).unwrap();
    let mut perm = std::fs::metadata(&p).unwrap().permissions();
    perm.set_mode(0o755);
    std::fs::set_permissions(&p, perm).unwrap();
    p
}

fn out_dir(tag: &str) -> PathBuf {
    let d = root().join("target/test-lift").join(tag);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Build a CRAM whose minimap2 `@PG` says `--secondary=yes`.
fn make_secondary_yes_cram(old_cram: &Path, old_fa: &Path, out: &Path) {
    let h = Command::new("samtools")
        .args(["view", "-h", "--no-PG"])
        .arg(old_cram)
        .output()
        .expect("samtools view -h");
    assert!(h.status.success(), "samtools view -h failed");
    let sam = String::from_utf8(h.stdout).unwrap().replace("--secondary=no", "--secondary=yes");
    assert!(sam.contains("--secondary=yes"), "CL did not contain --secondary=no");
    let mut child = Command::new("samtools")
        .args(["view", "-C", "--no-PG", "-T"])
        .arg(old_fa)
        .args(["-o"])
        .arg(out)
        .arg("-")
        .stdin(Stdio::piped())
        .spawn()
        .expect("samtools view -C");
    child.stdin.take().unwrap().write_all(sam.as_bytes()).unwrap();
    assert!(child.wait().unwrap().success(), "samtools view -C failed");
}

/// The raw header text of a CRAM/BAM.
fn header_text(path: &Path) -> String {
    let r = bam::Reader::from_path(path).unwrap();
    String::from_utf8_lossy(r.header().as_bytes()).to_string()
}

/// Rewrite `in_cram` so its `cramhook lift` `@PG` has no argv-encoding marker
/// (the legacy format an older build wrote).
fn make_legacy_lift_cram(in_cram: &Path, ref_fa: &Path, out: &Path) {
    let h = Command::new("samtools")
        .args(["view", "-h", "--no-PG"])
        .arg(in_cram)
        .output()
        .expect("samtools view -h");
    assert!(h.status.success(), "samtools view -h failed");
    let sam = String::from_utf8(h.stdout).unwrap().replace(" argv=sh1", "");
    assert!(!sam.contains("argv=sh1"), "marker not removed");
    let mut child = Command::new("samtools")
        .args(["view", "-C", "--no-PG", "-T"])
        .arg(ref_fa)
        .args(["-o"])
        .arg(out)
        .arg("-")
        .stdin(Stdio::piped())
        .spawn()
        .expect("samtools view -C");
    child.stdin.take().unwrap().write_all(sam.as_bytes()).unwrap();
    assert!(child.wait().unwrap().success(), "samtools view -C failed");
}

/// Write `in_cram` with an extra `@PG` line inserted at the end of the header.
fn make_old_name_cram(in_cram: &Path, ref_fa: &Path, out: &Path, pg: &str) {
    let h = Command::new("samtools")
        .args(["view", "-h", "--no-PG"])
        .arg(in_cram)
        .output()
        .expect("samtools view -h");
    assert!(h.status.success(), "samtools view -h failed");
    let text = String::from_utf8(h.stdout).unwrap();
    let mut lines: Vec<&str> = text.lines().collect();
    let pos = lines.iter().position(|l| !l.starts_with('@')).unwrap_or(lines.len());
    lines.insert(pos, pg);
    let sam = lines.join("\n") + "\n";
    let mut child = Command::new("samtools")
        .args(["view", "-C", "--no-PG", "-T"])
        .arg(ref_fa)
        .args(["-o"])
        .arg(out)
        .arg("-")
        .stdin(Stdio::piped())
        .spawn()
        .expect("samtools view -C");
    child.stdin.take().unwrap().write_all(sam.as_bytes()).unwrap();
    assert!(child.wait().unwrap().success(), "samtools view -C failed");
}

#[test]
fn lift_contract_fallback_all_remapped() {
    let _g = lock();
    setup();
    let dir = root().join("data/synthetic").join("S1");
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let old_cram = dir.join("hifi").join("old.cram");
    let out_dir = out_dir("contract");
    let bad = out_dir.join("secondary_yes.cram");
    make_secondary_yes_cram(&old_cram, &old_fa, &bad);
    let out = out_dir.join("out.cram");
    let opts = LiftOptions { strict: false, keep_tmp: false, min_input_mb: Some(0.0) };
    let stats = lift(&bad, &old_fa, &new_fa, &out, None, 4, &default_cfg(true), &opts, None)
        .unwrap_or_else(|e| panic!("lift fallback failed: {e:#}"));
    assert!(!stats.contract_ok, "expected contract violation");
    assert!(
        stats.contract_reason.as_deref().unwrap_or("").contains("secondary=yes"),
        "header fallback reason: {:?}",
        stats.contract_reason
    );
    assert_eq!(stats.preset, "map-hifi", "preset must be recovered from @PG");
    assert_eq!(stats.translated_reads, 0, "fallback must translate nothing");
    assert_eq!(stats.remapped_reads, stats.records_in.primary);

    let mut reader = bam::Reader::from_path(&out).unwrap();
    reader.set_reference(&new_fa).unwrap();
    for r in reader.records() {
        let r = r.unwrap();
        match r.aux(b"XL") {
            Ok(Aux::Char(b'R')) => {}
            other => panic!("fallback record not XL:A:R: {other:?}"),
        }
    }
}

/// The reader thread failing on a bad @SQ must not hang: it owns minimap2's
/// stdout (so the fd closes), the pass stops, and the children are killed.
#[test]
fn lift_reader_failure_does_not_hang() {
    let _g = lock();
    setup();
    let dir = root().join("data/synthetic").join("S1");
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let old_cram = dir.join("hifi").join("old.cram");
    // Emit a wrong @SQ header, then a lot to stdout; with the fd held open by
    // the main thread this would deadlock.
    let wrapper = write_script(
        "wrong_sq.sh",
        "#!/bin/sh\nprintf '@HD\\tVN:1.6\\tSO:unsorted\\n@SQ\\tSN:WRONG\\tLN:1\\n'\ndd if=/dev/zero bs=1048576 count=16 2>/dev/null\ncat > /dev/null\n",
    );
    let tag = out_dir("reader_fail");
    let out = tag.join("out.cram");
    let (tx, rx) = mpsc::channel();
    let mut cfg = default_cfg(true);
    cfg.xl_tag = false;
    let opts = LiftOptions { strict: false, keep_tmp: false, min_input_mb: Some(0.0) };
    let out2 = out.clone();
    let old_fa2 = old_fa.clone();
    let new_fa2 = new_fa.clone();
    let old_cram2 = old_cram.clone();
    let wrapper2 = wrapper.clone();
    thread::spawn(move || {
        unsafe { std::env::set_var("CRAMHOOK_MINIMAP2", &wrapper2) };
        let r = lift(&old_cram2, &old_fa2, &new_fa2, &out2, None, 4, &cfg, &opts, None);
        unsafe { std::env::remove_var("CRAMHOOK_MINIMAP2") };
        tx.send(r).ok();
    });
    let got = rx.recv_timeout(Duration::from_secs(90));
    let r = match got {
        Ok(r) => r,
        Err(_) => panic!("lift hung on reader-thread failure"),
    };
    let e = format!("{:#}", r.expect_err("reader failure must be an error"));
    assert!(
        e.contains("@SQ") || e.contains("reader"),
        "unexpected error: {e}"
    );
    assert!(!out.exists(), "staged output was not cleaned up");
    let (staged, staged_crai) = cramhook::lift::staged_paths(&out);
    assert!(!staged.exists(), "staged cram remains");
    assert!(!staged_crai.exists(), "staged crai remains");
}

fn run_injected(exe_env: &str, script_body: &str, tag: &str) -> (String, PathBuf) {
    let dir = root().join("data/synthetic").join("S1");
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let old_cram = dir.join("hifi").join("old.cram");
    let wrapper = write_script(&format!("{tag}.sh"), script_body);
    let d = out_dir(tag);
    let out = d.join("out.cram");
    // A previous good output (and its index) must survive a failed run.
    std::fs::write(&out, b"PREVIOUS").unwrap();
    let mut out_crai = out.as_os_str().to_owned();
    out_crai.push(".crai");
    let out_crai = PathBuf::from(out_crai);
    std::fs::write(&out_crai, b"OLDCRAI").unwrap();
    unsafe { std::env::set_var(exe_env, &wrapper) };
    let r = lift(&old_cram, &old_fa, &new_fa, &out, None, 4, &default_cfg(false), &LiftOptions { strict: false, keep_tmp: false, min_input_mb: Some(0.0) }, None);
    unsafe { std::env::remove_var(exe_env) };
    let e = format!("{:#}", r.expect_err("injected failure must be an error"));
    (e, out)
}

#[test]
fn lift_minimap2_exit43_reports_status() {
    let _g = lock();
    setup();
    let (e, out) = run_injected("CRAMHOOK_MINIMAP2", "#!/bin/sh\nexit 43\n", "mm43");
    assert!(e.contains("minimap2 exited with status 43"), "error: {e}");
    assert_eq!(std::fs::read(&out).unwrap(), b"PREVIOUS", "previous output was clobbered");
    assert_eq!(std::fs::read(out.with_extension("cram.crai")).unwrap(), b"OLDCRAI", "previous index was clobbered");
    let (staged, staged_crai) = cramhook::lift::staged_paths(&out);
    assert!(!staged.exists(), "staged cram remains");
    assert!(!staged_crai.exists(), "staged crai remains");
}

#[test]
fn lift_minimap2_stderr_surfaced_on_failure() {
    let _g = lock();
    setup();
    let script = "#!/bin/sh\necho 'minimap2 boom: index unreadable' >&2\nexit 43\n";
    let (e, _out) = run_injected("CRAMHOOK_MINIMAP2", script, "mm43err");
    assert!(e.contains("minimap2 exited with status 43"), "error: {e}");
    assert!(e.contains("minimap2 boom: index unreadable"), "captured stderr missing: {e}");
}

#[test]
fn lift_samtools_stderr_surfaced_on_failure() {
    let _g = lock();
    setup();
    let script = "#!/bin/sh\ncat > /dev/null\necho 'samtools boom: reference missing' >&2\nexit 42\n";
    let (e, _out) = run_injected("CRAMHOOK_SAMTOOLS", script, "sort42err");
    assert!(e.contains("samtools sort exited with status 42"), "error: {e}");
    assert!(e.contains("samtools boom: reference missing"), "captured stderr missing: {e}");
}

#[test]
fn lift_sort_exit42_reports_status() {
    let _g = lock();
    setup();
    let (e, out) = run_injected("CRAMHOOK_SAMTOOLS", "#!/bin/sh\ncat > /dev/null\nexit 42\n", "sort42");
    assert!(e.contains("samtools sort exited with status 42"), "error: {e}");
    assert_eq!(std::fs::read(&out).unwrap(), b"PREVIOUS", "previous output was clobbered");
    assert_eq!(std::fs::read(out.with_extension("cram.crai")).unwrap(), b"OLDCRAI", "previous index was clobbered");
    let (staged, staged_crai) = cramhook::lift::staged_paths(&out);
    assert!(!staged.exists(), "staged cram remains");
    assert!(!staged_crai.exists(), "staged crai remains");
}

/// If promotion fails, the previous OUT/OUT.crai pair must be restored intact.
#[test]
fn lift_promotion_failure_restores_previous_output() {
    let _g = lock();
    setup();
    let dir = root().join("data/synthetic").join("S1");
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let old_cram = dir.join("hifi").join("old.cram");
    let d = out_dir("promote_fail");
    let out = d.join("out.cram");
    std::fs::write(&out, b"PREVIOUS").unwrap();
    // A directory at the index target makes the index rename fail (EISDIR).
    let out_crai = out.with_extension("cram.crai");
    std::fs::create_dir_all(&out_crai).unwrap();

    let r = lift(&old_cram, &old_fa, &new_fa, &out, None, 4, &default_cfg(true), &LiftOptions { strict: false, keep_tmp: false, min_input_mb: Some(0.0) }, None);
    let e = format!("{:#}", r.expect_err("promotion must fail"));
    assert!(e.contains("renaming"), "error: {e}");
    assert_eq!(std::fs::read(&out).unwrap(), b"PREVIOUS", "previous CRAM was not restored");
    assert!(out_crai.is_dir(), "previous index target was disturbed");
    assert!(!out.with_extension("cram.bak").exists(), "left a .bak file");
    assert!(!out_crai.with_extension("crai.bak").exists(), "left a .crai.bak");
    let (staged, staged_crai) = cramhook::lift::staged_paths(&out);
    assert!(!staged.exists(), "staged cram remains");
    assert!(!staged_crai.exists(), "staged crai remains");
    std::fs::remove_dir_all(&out_crai).ok();
}

/// Build a tiny BAM with a contract-valid header plus a secondary record, to
/// trigger a mid-stream contract violation.
fn write_midstream_input(path: &Path) {
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", 5000));
    header.push_record(
        HeaderRecord::new(b"PG")
            .push_tag(b"ID", "minimap2")
            .push_tag(b"PN", "minimap2")
            .push_tag(b"VN", "2.31")
            .push_tag(b"CL", "minimap2 -ax map-hifi --secondary=no -Y old.fa reads.fq"),
    );
    let mut w = Writer::from_path(path, &header, Format::Bam).unwrap();
    let seq = vec![b'A'; 20];
    let qual = vec![30u8; 20];
    // primary, unmapped
    let mut p = Record::new();
    p.set(b"r1", None, &seq, &qual);
    p.set_tid(-1);
    p.set_pos(-1);
    p.set_flags(0x4);
    w.write(&p).unwrap();
    // secondary record -> mid-stream contract violation
    let mut s = Record::new();
    let c = CigarString(vec![Cigar::Match(20)]);
    s.set(b"r2", Some(&c), &seq, &qual);
    s.set_tid(0);
    s.set_pos(100);
    s.set_flags(0x100);
    w.write(&s).unwrap();
}

/// Build a tiny contract-valid BAM with a record-level contract violation but
/// no header violation: `secondary` selects one secondary record (0x100),
/// otherwise one hard-clipped primary. Also writes matching old.fa/new.fa
/// (chr, 5000 bp). Returns `(old_fa, new_fa, bam)`.
fn write_contract_input(dir: &Path, secondary: bool) -> (PathBuf, PathBuf, PathBuf) {
    std::fs::create_dir_all(dir).unwrap();
    let seq = "ACGT".repeat(1250);
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    std::fs::write(&old_fa, format!(">chr\n{seq}\n")).unwrap();
    std::fs::write(&new_fa, format!(">chr\n{seq}\n")).unwrap();
    let bam = dir.join("old.bam");
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr").push_tag(b"LN", 5000));
    header.push_record(
        HeaderRecord::new(b"PG")
            .push_tag(b"ID", "minimap2")
            .push_tag(b"PN", "minimap2")
            .push_tag(b"VN", "2.31")
            .push_tag(b"CL", "minimap2 -ax map-hifi --secondary=no -Y old.fa reads.fq"),
    );
    let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();
    let bases = vec![b'A'; 20];
    let qual = vec![30u8; 20];
    if secondary {
        // primary, unmapped (a legitimate read to remap)
        let mut p = Record::new();
        p.set(b"r1", None, &bases, &qual);
        p.set_tid(-1);
        p.set_pos(-1);
        p.set_flags(0x4);
        w.write(&p).unwrap();
        // secondary record -> record-level contract violation
        let mut s = Record::new();
        let c = CigarString(vec![Cigar::Match(20)]);
        s.set(b"r2", Some(&c), &bases, &qual);
        s.set_tid(0);
        s.set_pos(100);
        s.set_flags(0x100);
        w.write(&s).unwrap();
    } else {
        // hard-clipped primary -> record-level contract violation
        let mut p = Record::new();
        let c = CigarString(vec![Cigar::HardClip(5), Cigar::Match(20)]);
        p.set(b"r1", Some(&c), &bases, &qual);
        p.set_tid(0);
        p.set_pos(100);
        p.set_flags(0);
        w.write(&p).unwrap();
    }
    (old_fa, new_fa, bam)
}

#[test]
fn lift_midstream_fallback_stats() {
    let _g = lock();
    setup();
    let d = out_dir("midstream");
    let seq = "ACGT".repeat(1250);
    let old_fa = d.join("old.fa");
    let new_fa = d.join("new.fa");
    std::fs::write(&old_fa, format!(">chr\n{seq}\n")).unwrap();
    std::fs::write(&new_fa, format!(">chr\n{seq}\n")).unwrap();
    let bam = d.join("old.bam");
    write_midstream_input(&bam);
    let out = d.join("out.cram");
    let stats = lift(&bam, &old_fa, &new_fa, &out, None, 4, &default_cfg(false), &LiftOptions { strict: false, keep_tmp: false, min_input_mb: Some(0.0) }, None)
        .unwrap_or_else(|e| panic!("mid-stream fallback failed: {e:#}"));
    assert!(!stats.contract_ok, "mid-stream fallback must set contract_ok=false");
    assert!(
        stats.contract_reason.as_deref().unwrap_or("").contains("secondary"),
        "mid-stream reason: {:?}",
        stats.contract_reason
    );
    assert_eq!(stats.translated_reads, 0);
    assert_eq!(stats.remapped_reads, stats.records_in.primary);
    assert_eq!(stats.preset, "map-hifi");
}

/// Supplementary sets of translated (`XL:T`) reads, per query name.
fn supp_by_name(path: &Path, ref_fa: &Path) -> std::collections::HashMap<Vec<u8>, HashSet<(i32, i64)>> {
    let mut r = bam::Reader::from_path(path).unwrap();
    r.set_reference(ref_fa).unwrap();
    let mut m: std::collections::HashMap<Vec<u8>, HashSet<(i32, i64)>> = Default::default();
    for rec in r.records() {
        let rec = rec.unwrap();
        if rec.flags() & 0x800 != 0 {
            m.entry(rec.qname().to_vec()).or_default().insert((rec.tid(), rec.pos()));
        }
    }
    m
}

fn translated_xl_names(path: &Path, ref_fa: &Path) -> HashSet<Vec<u8>> {
    let mut r = bam::Reader::from_path(path).unwrap();
    r.set_reference(ref_fa).unwrap();
    let mut s = HashSet::new();
    for rec in r.records() {
        let rec = rec.unwrap();
        if rec.flags() & (0x100 | 0x800) == 0 {
            if let Ok(Aux::Char(b'T')) = rec.aux(b"XL") {
                s.insert(rec.qname().to_vec());
            }
        }
    }
    s
}

/// The clip screen must not drop a supplementary record: for every translated
/// (`XL:T`) read, the output's supplementary set equals the oracle's. This is
/// the S5c/S6c chimera regression (bad reads were previously translated without
/// their p3 tail supplementary).
#[test]
fn lift_s5c_s6c_supplementary_matches_oracle() {
    let _g = lock();
    setup();
    for sid in ["S5c", "S6c"] {
        for model in ["hifi", "ont"] {
            let dir = root().join("data/synthetic").join(sid);
            let new_fa = dir.join("new.fa");
            let oracle = dir.join(model).join("oracle.cram");
            assert!(oracle.exists(), "{sid}/{model}: missing oracle.cram");
            let case = run_lift(sid, model, true);
            let names = translated_xl_names(&case.out, &new_fa);
            assert!(!names.is_empty(), "{sid}/{model}: no translated XL:T reads to compare");
            let t_supp = supp_by_name(&case.out, &new_fa);
            let o_supp = supp_by_name(&oracle, &new_fa);
            let empty = HashSet::new();
            for name in &names {
                let a = t_supp.get(name).unwrap_or(&empty);
                let b = o_supp.get(name).unwrap_or(&empty);
                assert_eq!(
                    a,
                    b,
                    "{sid}/{model}: supplementary differs for {:?}",
                    String::from_utf8_lossy(name)
                );
            }
        }
    }
}

/// Reviewer's repro: a 4 kb deletion in old `B` creates a new adjacency in new
/// `B` while `diff.novel` stays empty. A chimeric read whose clip spans that
/// join must not be translated without its supplementary; the edit-site k-mers
/// make the screen defer it (`clip_hit`). Before the fix these reads were
/// translated and lost their supplementary alignment.
#[test]
fn lift_clip_join_edit_site_supplementary_matches_oracle() {
    let _g = lock();
    let dir = root().join("tests/fixtures/clip_join");
    let old_bam = dir.join("old.bam");
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let oracle = dir.join("oracle.cram");
    assert!(old_bam.exists() && oracle.exists(), "missing clip_join fixture");
    let out_dir = root().join("target/test-lift/clip_join");
    std::fs::create_dir_all(&out_dir).unwrap();
    let out = out_dir.join("out.cram");
    let cfg = TranslateConfig { xl_tag: true, ..Default::default() };
    let stats = lift(
        &old_bam,
        &old_fa,
        &new_fa,
        &out,
        None,
        4,
        &cfg,
        &LiftOptions { strict: false, keep_tmp: false, min_input_mb: Some(0.0) },
        None,
    )
    .unwrap_or_else(|e| panic!("clip_join lift failed: {e:#}"));

    // The screen must fire on the new adjacency (before the fix it translated 3
    // chimeric reads without their supplementary).
    let clip_hits = stats.deferred_by_reason.get("clip_hit").copied().unwrap_or(0)
        + stats.deferred_by_reason.get("split_clip_hit").copied().unwrap_or(0);
    assert!(clip_hits >= 3, "clip screen did not fire on the edit site: {clip_hits} hits");

    // No translated (XL:T) read may differ from the oracle's supplementary set.
    let t_supp = supp_by_name(&out, &new_fa);
    let o_supp = supp_by_name(&oracle, &new_fa);
    let empty = HashSet::new();
    for name in &translated_xl_names(&out, &new_fa) {
        assert_eq!(
            t_supp.get(name).unwrap_or(&empty),
            o_supp.get(name).unwrap_or(&empty),
            "supplementary differs for {:?}",
            String::from_utf8_lossy(name)
        );
    }
}

/// Primary `(tid, pos, reverse, mapq, cigar)` per read.
fn primary_sig(
    path: &Path,
    ref_fa: &Path,
) -> std::collections::BTreeMap<Vec<u8>, (i32, i64, bool, u8, String)> {
    let mut r = bam::Reader::from_path(path).unwrap();
    r.set_reference(ref_fa).unwrap();
    let mut m = std::collections::BTreeMap::new();
    for rec in r.records() {
        let rec = rec.unwrap();
        if rec.flags() & (0x100 | 0x800) == 0 {
            let cigar = rec.cigar().to_string();
            let prev = m.insert(
                rec.qname().to_vec(),
                (rec.tid(), rec.pos(), rec.flags() & 0x10 != 0, rec.mapq(), cigar),
            );
            assert!(prev.is_none(), "duplicate primary read");
        }
    }
    m
}

/// Supplementary `(tid, pos, reverse, mapq, cigar)` set per read, MAPQ >= 20
/// (the repo's `compare.py` `SUPP_MIN_MAPQ` filter).
type SuppSig = std::collections::BTreeMap<Vec<u8>, std::collections::BTreeSet<(i32, i64, bool, u8, String)>>;

fn supp_sig(path: &Path, ref_fa: &Path) -> SuppSig {
    let mut r = bam::Reader::from_path(path).unwrap();
    r.set_reference(ref_fa).unwrap();
    let mut m: SuppSig = Default::default();
    for rec in r.records() {
        let rec = rec.unwrap();
        if rec.flags() & 0x800 != 0 && rec.mapq() >= 20 {
            m.entry(rec.qname().to_vec()).or_default().insert((
                rec.tid(),
                rec.pos(),
                rec.flags() & 0x10 != 0,
                rec.mapq(),
                rec.cigar().to_string(),
            ));
        }
    }
    m
}

/// All `@PG` ids and their `PP`, in order.
fn pg_ids(path: &Path) -> Vec<(String, Option<String>)> {
    let r = bam::Reader::from_path(path).unwrap();
    let text = String::from_utf8_lossy(r.header().as_bytes()).to_string();
    let mut out = Vec::new();
    for line in text.split('\n') {
        if !line.starts_with("@PG\t") {
            continue;
        }
        let mut id = None;
        let mut pp = None;
        for f in line.split('\t').skip(1) {
            if let Some(v) = f.strip_prefix("ID:") {
                id = Some(v.to_string());
            }
            if let Some(v) = f.strip_prefix("PP:") {
                pp = Some(v.to_string());
            }
        }
        out.push((id.expect("@PG without ID"), pp));
    }
    out
}

/// Run `cramhook lift` through the real CLI (so its `@PG` CL is genuine).
fn cli_lift(in_cram: &Path, old_fa: &Path, new_fa: &Path, out: &Path, extra: &[&str]) -> String {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_cramhook"));
    cmd.args(["lift", "--in"])
        .arg(in_cram)
        .args(["--old-ref"])
        .arg(old_fa)
        .args(["--new-ref"])
        .arg(new_fa)
        .args(["--out"])
        .arg(out)
        .args(["--threads", "4", "--min-input-mb", "0"])
        .args(extra);
    let o = cmd.output().expect("running cramhook lift");
    assert!(o.status.success(), "cramhook lift failed: {}", String::from_utf8_lossy(&o.stderr));
    String::from_utf8_lossy(&o.stderr).to_string()
}

/// Number of primary records carrying `XL:A:T` (translated reads).
fn translated_primaries(path: &Path, ref_fa: &Path) -> u64 {
    let mut r = bam::Reader::from_path(path).unwrap();
    r.set_reference(ref_fa).unwrap();
    let mut n = 0u64;
    for rec in r.records() {
        let rec = rec.unwrap();
        if rec.flags() & (0x100 | 0x800) == 0 && matches!(rec.aux(b"XL"), Ok(Aux::Char(b'T'))) {
            n += 1;
        }
    }
    n
}

/// Read a FASTA into `(name, sequence)`.
fn read_fasta(path: &Path) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            out.push((name.split_whitespace().next().unwrap().to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.extend_from_slice(line.trim().as_bytes());
        }
    }
    out
}

fn write_fasta(path: &Path, recs: &[(String, Vec<u8>)]) {
    let mut s = String::new();
    for (name, seq) in recs {
        s.push('>');
        s.push_str(name);
        s.push('\n');
        for chunk in seq.chunks(60) {
            s.push_str(std::str::from_utf8(chunk).unwrap());
            s.push('\n');
        }
    }
    std::fs::write(path, s).unwrap();
}

/// Revision C for the chain test: `new.fa` (B) plus a duplicated 5 kb segment of
/// its longest contig, as a new contig. This creates real competitors for the
/// reads in that region, so the second lift's competitor handling is exercised.
fn write_revision_c(b: &Path, c: &Path) {
    let mut recs = read_fasta(b);
    let (longest, _) = recs.iter().max_by_key(|(_, s)| s.len()).unwrap().clone();
    let seq = recs.iter().find(|(n, _)| *n == longest).unwrap().1.clone();
    let start = seq.len() / 3;
    let seg = seq[start..start + 5000].to_vec();
    recs.push(("dup1".to_string(), seg));
    write_fasta(c, &recs);
}

/// Chained input through the CLI: lift A->B, B->C (revision C with an inserted
/// duplicate) and a third lift C->C. The chain must translate (not fall back to
/// a full remap), keep unique `@PG` ids with an unambiguous `PP` chain, and its
/// placements must match a direct full remap of the reads onto C by the repo's
/// equivalence criteria (contig/pos/strand/cigar/mapq, plus supplementary sets).
#[test]
fn lift_chain_three_lifts_matches_full_remap() {
    let _g = lock();
    setup();
    let dir = root().join("data/synthetic").join("S1");
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let old_cram = dir.join("hifi").join("old.cram");
    let reads = dir.join("hifi").join("reads.fastq.gz");
    assert!(reads.exists(), "S1/hifi: missing reads.fastq.gz");

    let d = out_dir("chain3");
    let new2_fa = d.join("new2.fa");
    write_revision_c(&new_fa, &new2_fa);
    // An index for the CRAM decode / lift.
    assert!(
        Command::new("samtools").args(["faidx"]).arg(&new2_fa).status().unwrap().success(),
        "samtools faidx failed"
    );

    let out1 = d.join("out1.cram");
    let out2 = d.join("out2.cram");
    let out3 = d.join("out3.cram");

    // A->B, B->C, then C->C (a third lift, to exercise the chained header).
    let r1 = cli_lift(&old_cram, &old_fa, &new_fa, &out1, &[]);
    assert!(r1.contains("contract ok: true"), "lift1 contract: {r1}");
    // A newly written lift @PG carries the argv-encoding marker.
    assert!(
        header_text(&out1).contains("argv=sh1"),
        "lift1 output has no argv encoding marker"
    );
    let r2 = cli_lift(&out1, &new_fa, &new2_fa, &out2, &[]);
    assert!(r2.contains("contract ok: true"), "lift2 contract: {r2}");
    let r3 = cli_lift(&out2, &new2_fa, &new2_fa, &out3, &[]);
    assert!(r3.contains("contract ok: true"), "lift3 contract: {r3}");

    // The chain must translate, not remap everything.
    let t2 = translated_primaries(&out2, &new2_fa);
    assert!(t2 > 0, "lift2 translated nothing");
    assert!(r2.contains("translated:"), "no lift2 summary");
    assert!(!r2.contains("contract reason:"), "lift2 reported a contract reason: {r2}");

    // Unique @PG ids and every PP resolves to an existing id (an unambiguous
    // chain), across three lifts.
    let pgs = pg_ids(&out3);
    let ids: std::collections::HashSet<&str> = pgs.iter().map(|(id, _)| id.as_str()).collect();
    assert_eq!(ids.len(), pgs.len(), "duplicate @PG id in {:?}", pgs);
    for (id, pp) in &pgs {
        if let Some(pp) = pp {
            assert!(ids.contains(pp.as_str()), "@PG {id} PP:{pp} does not resolve");
        }
    }
    let cramhook_ids: Vec<&str> =
        pgs.iter().filter(|(id, _)| id.starts_with("cramhook")).map(|(id, _)| id.as_str()).collect();
    assert_eq!(cramhook_ids.len(), 3, "expected three cramhook @PGs, got {:?}", pgs);

    // Direct full remap of the reads onto C.
    let oracle = d.join("oracle_C.cram");
    let pipeline = format!(
        "set -o pipefail; minimap2 -ax map-hifi -t 4 --secondary=no -Y {c} {r} \
         | samtools sort -@ 4 -O CRAM --reference {c} -o {o} -",
        c = new2_fa.display(),
        r = reads.display(),
        o = oracle.display()
    );
    let st = Command::new("bash").arg("-c").arg(pipeline).status().unwrap();
    assert!(st.success(), "building the C oracle failed");

    // Equivalence per the repo's criteria: primary (contig/pos/strand/cigar/mapq)
    // and supplementary sets (MAPQ >= 20). The second lift and the third (C->C
    // identity) lift must both match the direct full remap.
    let b = primary_sig(&oracle, &new2_fa);
    let sb = supp_sig(&oracle, &new2_fa);
    let a = primary_sig(&out2, &new2_fa);
    assert_eq!(a, b, "chain primary placements differ from the full remap");
    let sa = supp_sig(&out2, &new2_fa);
    assert_eq!(sa, sb, "chain supplementary sets differ from the full remap");

    let t3 = translated_primaries(&out3, &new2_fa);
    assert!(t3 > 0, "third lift translated nothing");
    assert_eq!(primary_sig(&out3, &new2_fa), b, "third lift primary placements differ");
    assert_eq!(supp_sig(&out3, &new2_fa), sb, "third lift supplementary sets differ");
}

/// Fail-safe chain rule: inherited translation is only allowed for an
/// all-default chain. A previous lift that used an eligibility/safety option
/// (here `--max-s2-ratio`) makes the contract fail, and a non-default current
/// lift does too; both fully remap. A default->default chain translates.
#[test]
fn lift_chain_requires_all_default_options() {
    let _g = lock();
    setup();
    let dir = root().join("data/synthetic").join("S1");
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let old_cram = dir.join("hifi").join("old.cram");

    let chain = |tag: &str, first: &[&str], second: &[&str]| -> (u64, String) {
        let d = out_dir(&format!("chain_default/{tag}"));
        let o1 = d.join("o1.cram");
        let o2 = d.join("o2.cram");
        cli_lift(&old_cram, &old_fa, &new_fa, &o1, first);
        let err = cli_lift(&o1, &new_fa, &new_fa, &o2, second);
        (translated_primaries(&o2, &new_fa), err)
    };

    let (dd, e) = chain("default_default", &[], &[]);
    assert!(dd > 0, "default -> default must translate");
    assert!(e.contains("contract ok: true"), "{e}");

    let (rd, e) = chain("relaxed_default", &["--max-s2-ratio", "1"], &[]);
    assert!(e.contains("contract ok: false"), "relaxed-lift-1 chain must fail the contract: {e}");
    assert_eq!(rd, 0, "relaxed -> default must remap");

    let (ds, e) = chain("default_stricter", &[], &["--max-s2-ratio", "0.4"]);
    assert!(e.contains("contract ok: false"), "non-default current lift must fail the contract: {e}");
    assert_eq!(ds, 0, "default -> non-default current must remap");
}

/// Supplementary placements `(tid, pos, reverse, mapq)` per read. The clip_join
/// oracle uses hard clips where `-Y` writes soft clips, so cigars are not
/// comparable there; placements are what "loses its supplementary" means.
type SuppPlacements =
    std::collections::BTreeMap<Vec<u8>, std::collections::BTreeSet<(i32, i64, bool, u8)>>;

fn supp_placements(path: &Path, ref_fa: &Path) -> SuppPlacements {
    let mut r = bam::Reader::from_path(path).unwrap();
    r.set_reference(ref_fa).unwrap();
    let mut m: SuppPlacements = Default::default();
    for rec in r.records() {
        let rec = rec.unwrap();
        if rec.flags() & 0x800 != 0 {
            m.entry(rec.qname().to_vec())
                .or_default()
                .insert((rec.tid(), rec.pos(), rec.flags() & 0x10 != 0, rec.mapq()));
        }
    }
    m
}

/// clip_join: the clip screen is a safety check. A default->default two-lift
/// chain must equal a direct full remap including the supplementaries of
/// junction45/50/55; a first lift that weakened the screen makes the second
/// lift's contract fail, so it fully remaps (the only safe behaviour).
#[test]
fn lift_chain_clip_join_default_equals_full_remap() {
    let _g = lock();
    let cj = root().join("tests/fixtures/clip_join");
    let old_bam = cj.join("old.bam");
    let old_fa = cj.join("old.fa");
    let new_fa = cj.join("new.fa");
    let oracle = cj.join("oracle.cram");
    assert!(oracle.exists(), "missing clip_join oracle");
    let oracle_supp = supp_placements(&oracle, &new_fa);

    let run_chain = |tag: &str, first: &[&str]| -> (String, PathBuf) {
        let d = out_dir(&format!("chain_clip/{tag}"));
        let o1 = d.join("o1.cram");
        let o2 = d.join("o2.cram");
        cli_lift(&old_bam, &old_fa, &new_fa, &o1, first);
        let err = cli_lift(&o1, &new_fa, &new_fa, &o2, &[]);
        (err, o2)
    };

    // Default -> default: contract ok and supplementaries preserved.
    let (err, o2) = run_chain("default", &[]);
    assert!(err.contains("contract ok: true"), "default chain contract: {err}");
    let a = supp_placements(&o2, &new_fa);
    assert_eq!(a, oracle_supp, "default chain lost supplementaries vs the full remap");
    for name in ["junction45", "junction50", "junction55"] {
        assert!(
            !a.get(name.as_bytes()).is_none_or(|s| s.is_empty()),
            "{name} lost its supplementary"
        );
    }

    // A first lift with a weakened clip screen: the chain is not all-default, so
    // the second lift fails the contract and fully remaps (matching the oracle).
    let (err, o2) = run_chain("relaxed", &["--clip-min-hits", "100000"]);
    assert!(err.contains("contract ok: false"), "relaxed-screen chain must fail the contract: {err}");
    assert!(err.contains("--clip-min-hits"), "reason must name the option: {err}");
    assert_eq!(supp_placements(&o2, &new_fa), oracle_supp, "relaxed-screen chain must remap");
}

/// clip_join, sharpened: `old -> old` (identity) has no edit site, so the clip
/// screen does not fire and lift 1 translates junction45/50/55. The `old -> new`
/// lift 2 must then defer those three as `clip_hit` (the join is a new edit site)
/// and remap them, restoring their supplementaries equal to a full remap.
#[test]
fn lift_chain_clip_join_old_to_old_then_new_defers_clip_hits() {
    let _g = lock();
    let cj = root().join("tests/fixtures/clip_join");
    let old_bam = cj.join("old.bam");
    let old_fa = cj.join("old.fa");
    let new_fa = cj.join("new.fa");
    let oracle = cj.join("oracle.cram");
    assert!(oracle.exists(), "missing clip_join oracle");
    let oracle_supp = supp_placements(&oracle, &new_fa);

    let d = out_dir("chain_clip/old_old");
    let o1 = d.join("o1.cram");
    let o2 = d.join("o2.cram");

    // Identity old->old: no edit site, so the three junction reads translate.
    let e1 = cli_lift(&old_bam, &old_fa, &old_fa, &o1, &[]);
    assert!(e1.contains("contract ok: true"), "lift1 contract: {e1}");
    assert_eq!(translated_primaries(&o1, &old_fa), 3, "old->old must translate the 3 junction reads");
    let p1 = primary_sig(&o1, &old_fa);
    for name in ["junction45", "junction50", "junction55"] {
        assert!(p1.contains_key(name.as_bytes()), "{name} not translated by lift1");
    }

    // old->new: the join is a new edit site, so the clip screen fires and the
    // three reads are deferred `clip_hit` and remapped (supplementaries restored).
    let e2 = cli_lift(&o1, &old_fa, &new_fa, &o2, &[]);
    assert!(e2.contains("contract ok: true"), "lift2 contract: {e2}");
    assert!(e2.contains("clip_hit"), "lift2 must defer clip hits: {e2}");
    let a = supp_placements(&o2, &new_fa);
    assert_eq!(a, oracle_supp, "lift2 must restore the supplementaries equal to the full remap");
    for name in ["junction45", "junction50", "junction55"] {
        assert!(
            a.get(name.as_bytes()).is_some_and(|s| !s.is_empty()),
            "{name} lost its supplementary"
        );
    }
}

/// A legacy lift `@PG` (an older build's `.join(" ")` CL, no argv marker) must
/// not be chained: the contract fails and the input is fully remapped. The
/// reviewer's legacy CL can hide a gated option behind an apostrophe, so the
/// marker is required rather than trusting the CL.
#[test]
fn lift_chain_rejects_legacy_lift_pg_without_marker() {
    let _g = lock();
    setup();
    let dir = root().join("data/synthetic").join("S1");
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let old_cram = dir.join("hifi").join("old.cram");

    let d = out_dir("chain_legacy");
    let out1 = d.join("out1.cram");
    let legacy = d.join("legacy.cram");
    let out2 = d.join("out2.cram");

    let r1 = cli_lift(&old_cram, &old_fa, &new_fa, &out1, &[]);
    assert!(r1.contains("contract ok: true"), "lift1 contract: {r1}");
    make_legacy_lift_cram(&out1, &new_fa, &legacy);
    assert!(!header_text(&legacy).contains("argv=sh1"), "marker not stripped");

    let err = cli_lift(&legacy, &new_fa, &new_fa, &out2, &[]);
    assert!(
        err.contains("contract ok: false"),
        "a legacy lift @PG must fail the contract: {err}"
    );
    assert!(err.contains("older cramhook"), "reason must name the legacy build: {err}");
    assert_eq!(
        translated_primaries(&out2, &new_fa),
        0,
        "a legacy chain must remap everything"
    );
}

/// A lift `@PG` under the tool's old name is a hard error in the entry points,
/// even without `--strict`: neither `lift` nor `translate` may fall back to a
/// full remap or write output. Checked for `PN:cramlift` and for an `ID` prefix.
#[test]
fn lift_and_translate_reject_old_name_lift_pg() {
    let _g = lock();
    setup();
    let dir = root().join("data/synthetic").join("S1");
    let old_fa = dir.join("old.fa");
    let new_fa = dir.join("new.fa");
    let old_cram = dir.join("hifi").join("old.cram");
    let d = out_dir("old_name");

    let cases = [
        (
            "pn",
            "@PG\tID:legacy\tPN:cramlift\tVN:0.1.0\tCL:/x/cramlift lift --in old.cram",
        ),
        (
            "id",
            "@PG\tID:cramlift.1\tPN:something-else\tVN:0.1.0\tCL:other tool",
        ),
    ];
    for (tag, pg) in cases {
        let input = d.join(format!("{tag}.cram"));
        make_old_name_cram(&old_cram, &old_fa, &input, pg);

        let lift_out = d.join(format!("{tag}.lift.cram"));
        let o = Command::new(env!("CARGO_BIN_EXE_cramhook"))
            .args(["lift", "--in"])
            .arg(&input)
            .args(["--old-ref"])
            .arg(&old_fa)
            .args(["--new-ref"])
            .arg(&new_fa)
            .args(["--out"])
            .arg(&lift_out)
            .args(["--threads", "2", "--min-input-mb", "0"])
            .output()
            .expect("running cramhook lift");
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(!o.status.success(), "{tag}: lift must fail; stderr: {err}");
        assert!(err.contains("old name"), "{tag}: lift stderr must name the old tool: {err}");
        assert!(!lift_out.exists(), "{tag}: lift wrote output before rejecting the input");

        let prefix = d.join(format!("{tag}.tr"));
        let o = Command::new(env!("CARGO_BIN_EXE_cramhook"))
            .args(["translate", "--in"])
            .arg(&input)
            .args(["--old-ref"])
            .arg(&old_fa)
            .args(["--new-ref"])
            .arg(&new_fa)
            .args(["--out-prefix"])
            .arg(&prefix)
            .args(["--threads", "2"])
            .output()
            .expect("running cramhook translate");
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(!o.status.success(), "{tag}: translate must fail; stderr: {err}");
        assert!(err.contains("old name"), "{tag}: translate stderr must name the old tool: {err}");
        for ext in ["translated.bam", "deferred.fq", "deferred.tsv", "stats.json"] {
            let p = PathBuf::from(format!("{}.{}", prefix.display(), ext));
            assert!(!p.exists(), "{tag}: translate wrote {} before rejecting", p.display());
        }
    }
}
