//! Integration tests for `--progress` (lift, qc, fish): the stderr stream is
//! well formed (plan first, every phase started and ended in plan order,
//! monotonic `frac` in [0, 1], `done` last) and the outputs do not change.
//!
//! Uses the synthetic S1 fixtures (`pixi run sim && pixi run crams-synthetic`).

use std::path::{Path, PathBuf};
use std::process::Command;

use rust_htslib::bam::{self, Read};
use serde_json::Value;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn s1() -> PathBuf {
    let d = root().join("data/synthetic/S1");
    assert!(
        d.join("hifi/old.cram").exists(),
        "data/synthetic is missing: run `pixi run sim` and `pixi run crams-synthetic`"
    );
    d
}

fn workdir(name: &str) -> PathBuf {
    let d = root().join("target/test-progress").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Run the CLI; returns stderr. Updates are not rate limited so the small
/// fixtures produce some.
fn cramhook(args: &[&str]) -> String {
    let o = Command::new(env!("CARGO_BIN_EXE_cramhook"))
        .args(args)
        .env("CRAMHOOK_PROGRESS_INTERVAL_MS", "0")
        .output()
        .expect("running cramhook");
    let err = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(o.status.success(), "cramhook {args:?} failed: {err}");
    err
}

fn p(path: &Path) -> &str {
    path.to_str().unwrap()
}

/// Check the progress lines of `stderr`; returns (phases, number of updates).
fn check_stream(stderr: &str, cmd: &str) -> (Vec<String>, usize) {
    let ev: Vec<Value> = stderr
        .lines()
        .filter(|l| l.starts_with('{'))
        .map(|l| serde_json::from_str(l).unwrap_or_else(|e| panic!("bad JSON {l:?}: {e}")))
        .collect();
    assert!(ev.len() >= 2, "no progress stream in:\n{stderr}");
    let mut last_t = 0.0;
    for v in &ev {
        assert_eq!(v["progress"], 1, "{v}");
        assert_eq!(v["cmd"], cmd, "{v}");
        let t = v["elapsed_s"].as_f64().unwrap_or_else(|| panic!("no elapsed_s: {v}"));
        assert!(t >= last_t, "elapsed_s decreased: {v}");
        last_t = t;
    }
    assert_eq!(ev[0]["event"], "plan", "plan must come first: {}", ev[0]);
    let plan: Vec<String> =
        ev[0]["phases"].as_array().unwrap().iter().map(|x| x.as_str().unwrap().to_string()).collect();
    assert!(!plan.is_empty());
    let mut i = 1;
    let mut updates = 0;
    for ph in &plan {
        assert_eq!(ev[i]["event"], "start", "expected start of {ph}: {}", ev[i]);
        assert_eq!(ev[i]["phase"], ph.as_str());
        i += 1;
        let mut last = 0.0;
        while ev[i]["event"] == "update" {
            assert_eq!(ev[i]["phase"], ph.as_str(), "{}", ev[i]);
            let f = ev[i]["frac"].as_f64().unwrap();
            assert!((0.0..=1.0).contains(&f), "frac out of range: {}", ev[i]);
            assert!(f >= last, "frac not monotonic: {}", ev[i]);
            last = f;
            updates += 1;
            i += 1;
        }
        assert_eq!(ev[i]["event"], "end", "expected end of {ph}: {}", ev[i]);
        assert_eq!(ev[i]["phase"], ph.as_str());
        assert!(ev[i].get("frac").is_none());
        i += 1;
    }
    assert_eq!(ev[i]["event"], "done", "{}", ev[i]);
    assert_eq!(i, ev.len() - 1, "lines after done");
    (plan, updates)
}

fn read_text(path: PathBuf) -> String {
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// JSON without the run-dependent keys.
fn json_stable(path: PathBuf) -> Value {
    fn strip(v: &mut Value) {
        if let Value::Object(m) = v {
            m.retain(|k, _| k != "timing_s" && k != "timings_s");
            m.values_mut().for_each(strip);
        }
    }
    let mut v: Value = serde_json::from_str(&read_text(path)).unwrap();
    strip(&mut v);
    v
}

fn with_prefix(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{suffix}", prefix.display()))
}

#[test]
fn progress_qc() {
    let s = s1();
    let d = workdir("qc");
    let cram = s.join("hifi/old.cram");
    let fa = s.join("old.fa");
    let base = ["qc", "--in", p(&cram), "--ref", p(&fa), "--threads", "4", "--chunk-size", "20000"];
    let a = d.join("off");
    let b = d.join("on");
    let err_off = cramhook(&[&base[..], &["--out-prefix", p(&a)]].concat());
    assert!(!err_off.lines().any(|l| l.starts_with('{')), "progress lines without --progress");
    let err = cramhook(&[&base[..], &["--out-prefix", p(&b), "--progress"]].concat());
    let (plan, updates) = check_stream(&err, "qc");
    assert_eq!(plan, ["scan", "circularity", "call", "linkage", "output"]);
    assert!(updates > 0, "no updates:\n{err}");
    for suf in [".vcf", ".bed", ".findings.bed", ".depth.bedgraph", ".haplotypes.tsv"] {
        assert_eq!(read_text(with_prefix(&a, suf)), read_text(with_prefix(&b, suf)), "{suf} differs");
    }
    assert_eq!(json_stable(with_prefix(&a, ".json")), json_stable(with_prefix(&b, ".json")));
}

#[test]
fn progress_fish_with_bait() {
    let s = s1();
    let d = workdir("fish");
    let cram = s.join("hifi/old.cram");
    let fa = s.join("old.fa");
    // Bait: a renamed 20 kb piece of the first contig.
    let text = read_text(fa.clone());
    let seq: String = text.lines().skip(1).take_while(|l| !l.starts_with('>')).collect();
    let bait = d.join("bait.fa");
    std::fs::write(&bait, format!(">bait1\n{}\n", &seq[10_000..30_000])).unwrap();
    let base = [
        "fish", "--in", p(&cram), "--ref", p(&fa), "--threads", "4", "--bait", p(&bait),
        "--bait-min-input-mb", "0",
    ];
    let a = d.join("off");
    let b = d.join("on");
    let err_off = cramhook(&[&base[..], &["--out-prefix", p(&a)]].concat());
    assert!(!err_off.lines().any(|l| l.starts_with('{')), "progress lines without --progress");
    let err = cramhook(&[&base[..], &["--out-prefix", p(&b), "--progress"]].concat());
    let (plan, _) = check_stream(&err, "fish");
    assert_eq!(plan, ["pass", "bait_lift", "payload", "output"]);
    // The internal bait lift reports as the `bait_lift` phase, not as lift lines.
    assert!(!err.contains("\"cmd\":\"lift\""));
    assert_eq!(read_text(with_prefix(&a, ".fish.tsv")), read_text(with_prefix(&b, ".fish.tsv")));
    let fq = |prefix: &Path| {
        let mut r = rust_htslib::bgzf::Reader::from_path(with_prefix(prefix, ".fastq.gz")).unwrap();
        let mut v = Vec::new();
        std::io::Read::read_to_end(&mut r, &mut v).unwrap();
        v
    };
    assert_eq!(fq(&a), fq(&b), "FASTQ differs");
    assert_eq!(json_stable(with_prefix(&a, ".fish.json")), json_stable(with_prefix(&b, ".fish.json")));
}

/// Decoded records, order-independent (records at one position may interleave
/// differently between runs).
fn records(cram: &Path, fa: &Path) -> Vec<bam::Record> {
    let mut r = bam::Reader::from_path(cram).unwrap();
    r.set_reference(fa).unwrap();
    let mut out: Vec<bam::Record> = r.records().map(|x| x.unwrap()).collect();
    out.sort_by(|x, y| {
        (x.qname(), x.flags(), x.tid(), x.pos()).cmp(&(y.qname(), y.flags(), y.tid(), y.pos()))
    });
    out
}

#[test]
fn progress_lift() {
    let s = s1();
    let d = workdir("lift");
    let cram = s.join("hifi/old.cram");
    let old_fa = s.join("old.fa");
    let new_fa = s.join("new.fa");
    let base = [
        "lift", "--in", p(&cram), "--old-ref", p(&old_fa), "--new-ref", p(&new_fa), "--threads", "4",
        "--min-input-mb", "0",
    ];
    let a = d.join("off.cram");
    let b = d.join("on.cram");
    let sa = d.join("off.json");
    let sb = d.join("on.json");
    let err_off = cramhook(&[&base[..], &["--out", p(&a), "--stats", p(&sa)]].concat());
    assert!(!err_off.lines().any(|l| l.starts_with('{')), "progress lines without --progress");
    let err = cramhook(&[&base[..], &["--out", p(&b), "--stats", p(&sb), "--progress"]].concat());
    let (plan, updates) = check_stream(&err, "lift");
    assert_eq!(plan, ["diff", "translate_pass", "writer_drain", "sort_index"]);
    assert!(updates > 0, "no updates:\n{err}");
    assert_eq!(json_stable(sa), json_stable(sb));
    assert_eq!(records(&a, &new_fa), records(&b, &new_fa), "lifted records differ");
}

#[test]
fn progress_error_line() {
    let s = s1();
    let d = workdir("error");
    let o = Command::new(env!("CARGO_BIN_EXE_cramhook"))
        .args(["qc", "--in", p(&d.join("missing.cram")), "--ref", p(&s.join("old.fa"))])
        .args(["--out-prefix", p(&d.join("x")), "--progress"])
        .output()
        .unwrap();
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    let ev: Vec<Value> =
        err.lines().filter(|l| l.starts_with('{')).map(|l| serde_json::from_str(l).unwrap()).collect();
    assert_eq!(ev.first().unwrap()["event"], "plan");
    let last = ev.last().unwrap();
    assert_eq!(last["event"], "error");
    assert!(last["message"].as_str().is_some_and(|m| !m.is_empty()));
    assert!(err.contains("Error:"), "the usual error message is still printed: {err}");
}

/// Run the CLI with stderr on a pipe whose read end is already closed (a GUI
/// that went away): every stderr write fails with EPIPE. Returns the exit status.
fn cramhook_closed_stderr(args: &[&str]) -> std::process::ExitStatus {
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    Command::new(env!("CARGO_BIN_EXE_cramhook"))
        .args(args)
        .stdout(std::process::Stdio::null())
        .stderr(writer)
        .status()
        .expect("running cramhook")
}

#[test]
fn progress_closed_stderr_does_not_fail() {
    let s = s1();
    let d = workdir("closed_stderr");
    let cram = s.join("hifi/old.cram");
    let old_fa = s.join("old.fa");
    let new_fa = s.join("new.fa");

    // lift: a large --min-input-mb makes the tiny-input announcement fire
    // before the stats are written.
    let lift = [
        "lift", "--in", p(&cram), "--old-ref", p(&old_fa), "--new-ref", p(&new_fa), "--threads", "4",
        "--min-input-mb", "100000",
    ];
    let (la, lb) = (d.join("lift_ok.cram"), d.join("lift_closed.cram"));
    let (sa, sb) = (d.join("lift_ok.json"), d.join("lift_closed.json"));
    let err = cramhook(&[&lift[..], &["--out", p(&la), "--stats", p(&sa)]].concat());
    assert!(err.contains("tiny input"), "the tiny-input announcement must fire: {err}");
    let st = cramhook_closed_stderr(&[&lift[..], &["--out", p(&lb), "--stats", p(&sb), "--progress"]].concat());
    assert!(st.success(), "lift with a closed stderr: {st:?}");
    let stats = json_stable(sb.clone());
    assert_eq!(stats["fallback"]["reason"], "tiny_input");
    assert_eq!(json_stable(sa), stats);
    assert_eq!(records(&la, &new_fa), records(&lb, &new_fa), "lifted records differ");

    // qc
    let qc = ["qc", "--in", p(&cram), "--ref", p(&old_fa), "--threads", "4"];
    let (qa, qb) = (d.join("qc_ok"), d.join("qc_closed"));
    cramhook(&[&qc[..], &["--out-prefix", p(&qa)]].concat());
    let st = cramhook_closed_stderr(&[&qc[..], &["--out-prefix", p(&qb), "--progress"]].concat());
    assert!(st.success(), "qc with a closed stderr: {st:?}");
    assert_eq!(read_text(with_prefix(&qa, ".vcf")), read_text(with_prefix(&qb, ".vcf")));
    assert_eq!(json_stable(with_prefix(&qa, ".json")), json_stable(with_prefix(&qb, ".json")));

    // fish
    let fish = ["fish", "--in", p(&cram), "--ref", p(&old_fa), "--threads", "4"];
    let (fa, fb) = (d.join("fish_ok"), d.join("fish_closed"));
    cramhook(&[&fish[..], &["--out-prefix", p(&fa)]].concat());
    let st = cramhook_closed_stderr(&[&fish[..], &["--out-prefix", p(&fb), "--progress"]].concat());
    assert!(st.success(), "fish with a closed stderr: {st:?}");
    assert_eq!(read_text(with_prefix(&fa, ".fish.tsv")), read_text(with_prefix(&fb, ".fish.tsv")));
    assert_eq!(json_stable(with_prefix(&fa, ".fish.json")), json_stable(with_prefix(&fb, ".fish.json")));

    // An error still exits 1 (not a panic) with a closed stderr.
    let st = cramhook_closed_stderr(&["qc", "--in", p(&d.join("missing.cram")), "--ref", p(&old_fa), "--out-prefix", p(&d.join("x")), "--progress"]);
    assert_eq!(st.code(), Some(1), "{st:?}");
}
