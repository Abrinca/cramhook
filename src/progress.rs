//! `--progress`: machine-readable progress on stderr, one JSON object per line
//! (see docs/reference/output-formats.md, "Progress stream").
//!
//! A [`Progress`] handle is cheap to clone and safe to share between threads.
//! When disabled ([`Progress::off`], the default) every method is a no-op
//! behind one `Option` check, so the outputs are unchanged.

use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use rust_htslib::bam::HeaderView;

/// Human-readable stderr line, best effort: like `eprintln!`, but a closed or
/// broken stderr (e.g. a GUI that closed its pipe) is ignored instead of
/// panicking.
#[macro_export]
macro_rules! note {
    ($($arg:tt)*) => {{
        use ::std::io::Write as _;
        let _ = ::std::writeln!(::std::io::stderr().lock(), $($arg)*);
    }};
}

/// Version of the progress stream (the `progress` key of every line).
pub const PROGRESS_VERSION: u32 = 1;

/// Progress reporter; a no-op when disabled.
#[derive(Clone, Default)]
pub struct Progress(Option<Arc<Inner>>);

struct Inner {
    cmd: &'static str,
    t0: Instant,
    state: Mutex<State>,
}

struct State {
    out: Box<dyn Write + Send>,
    phase: Option<&'static str>,
    /// Last emitted whole percent of the running phase.
    pct: u32,
    /// Exact record total of the running record pass, if known.
    total: Option<u64>,
    /// Last emitted record count (kept monotonic).
    records: u64,
}

impl std::fmt::Debug for Progress {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.0 {
            Some(i) => write!(f, "Progress({})", i.cmd),
            None => write!(f, "Progress(off)"),
        }
    }
}

impl Progress {
    /// Disabled reporter.
    pub fn off() -> Self {
        Progress(None)
    }

    /// Reporter for `cmd` writing to stderr.
    pub fn stderr(cmd: &'static str) -> Self {
        Self::to_writer(cmd, Box::new(std::io::stderr()))
    }

    /// Reporter for `cmd` writing to `out` (one `write_all` per line).
    pub fn to_writer(cmd: &'static str, out: Box<dyn Write + Send>) -> Self {
        Progress(Some(Arc::new(Inner {
            cmd,
            t0: Instant::now(),
            state: Mutex::new(State { out, phase: None, pct: 0, total: None, records: 0 }),
        })))
    }

    pub fn is_on(&self) -> bool {
        self.0.is_some()
    }

    /// Record total of an alignment file for a record pass (`None` when off
    /// or unknown; computed only when on).
    pub fn record_total(&self, path: &std::path::Path) -> Option<u64> {
        if self.is_on() { crate::cramcount::record_count(path) } else { None }
    }

    /// First line: the phases in the order they run.
    pub fn plan(&self, phases: &[&str]) {
        let Some(i) = &self.0 else { return };
        let list: Vec<String> = phases.iter().map(|p| format!("\"{p}\"")).collect();
        let mut st = i.state.lock().unwrap();
        i.emit(&mut st, &format!("\"event\":\"plan\",\"phases\":[{}]", list.join(",")));
    }

    /// Start `phase`. Starting the phase that is already running is a no-op
    /// (a restarted pass continues its phase; `frac` stays monotonic).
    pub fn start(&self, phase: &'static str) {
        self.start_with(phase, None);
    }

    /// Start a record pass whose total record count is `total` (`None`: the
    /// update lines carry an estimate).
    pub fn start_records(&self, phase: &'static str, total: Option<u64>) {
        self.start_with(phase, total);
    }

    fn start_with(&self, phase: &'static str, total: Option<u64>) {
        let Some(i) = &self.0 else { return };
        let mut st = i.state.lock().unwrap();
        if st.phase == Some(phase) {
            return;
        }
        st.phase = Some(phase);
        st.pct = 0;
        st.total = total;
        st.records = 0;
        i.emit(&mut st, &format!("\"phase\":\"{phase}\",\"event\":\"start\""));
    }

    /// Report the fraction done of the running phase (clamped to [0, 1]). A
    /// line is written only when the whole percent `floor(100 * frac)` exceeds
    /// the last one written for the phase.
    pub fn update(&self, frac: f64) {
        if self.0.is_some() {
            self.update_with(frac, None);
        }
    }

    /// [`update`](Self::update) for a record pass: `records` read so far.
    pub fn update_records(&self, frac: f64, records: u64) {
        if self.0.is_some() {
            self.update_with(frac, Some(records));
        }
    }

    fn update_with(&self, frac: f64, records: Option<u64>) {
        let Some(i) = &self.0 else { return };
        let mut frac = if frac.is_finite() { frac.clamp(0.0, 1.0) } else { return };
        let mut st = i.state.lock().unwrap();
        let Some(phase) = st.phase else { return };
        let records = records.map(|r| r.max(st.records));
        // An exact total that is exceeded was wrong: estimate for the rest of
        // the phase instead of writing an inconsistent line.
        if let (Some(r), Some(n)) = (records, st.total)
            && r > n
        {
            st.total = None;
        }
        // With an exact total, the records read also count as progress, so a
        // tail without positions (the unplaced unmapped reads) still advances;
        // and 100 % means every record was read.
        if let (Some(r), Some(n)) = (records, st.total.filter(|&n| n > 0)) {
            frac = frac.max((r as f64 / n as f64).min(1.0));
            if r < n {
                frac = frac.min(0.999);
            }
        }
        let pct = whole_pct(frac) as u32;
        if pct <= st.pct {
            return;
        }
        st.pct = pct;
        // Truncated, not rounded, so `1.000` only appears when done.
        let shown = (frac * 1000.0 + 1e-9).floor() / 1000.0;
        let mut body = format!("\"phase\":\"{phase}\",\"event\":\"update\",\"frac\":{shown:.3}");
        if let Some(r) = records {
            st.records = r;
            let (n, exact) = match st.total {
                Some(n) => (n, true),
                None => ((r as f64 / frac).round() as u64, false),
            };
            body.push_str(&format!(",\"records\":{r},\"records_total\":{n},\"total_exact\":{exact}"));
        }
        i.emit(&mut st, &body);
    }

    /// Exact record total of the running phase, if any.
    fn phase_total(&self) -> Option<u64> {
        let i = self.0.as_ref()?;
        i.state.lock().unwrap().total
    }

    /// End the running phase.
    pub fn end(&self) {
        let Some(i) = &self.0 else { return };
        let mut st = i.state.lock().unwrap();
        if let Some(phase) = st.phase.take() {
            i.emit(&mut st, &format!("\"phase\":\"{phase}\",\"event\":\"end\""));
        }
    }

    /// Last line of a successful run (ends a phase still running).
    pub fn done(&self) {
        self.end();
        let Some(i) = &self.0 else { return };
        let mut st = i.state.lock().unwrap();
        i.emit(&mut st, "\"event\":\"done\"");
    }

    /// Last line of a failed run; the error is still reported as usual.
    pub fn error(&self, message: &str) {
        let Some(i) = &self.0 else { return };
        let mut st = i.state.lock().unwrap();
        let msg = serde_json::to_string(message).unwrap_or_else(|_| "\"\"".to_string());
        i.emit(&mut st, &format!("\"event\":\"error\",\"message\":{msg}"));
    }

    /// Emit `done` or `error` for a finished run and pass the result through.
    pub fn finish<T>(&self, r: anyhow::Result<T>) -> anyhow::Result<T> {
        match &r {
            Ok(_) => self.done(),
            Err(e) => self.error(&format!("{e:#}")),
        }
        r
    }
}

/// `floor(100 * frac)`, robust to rounding (`0.29 * 100` is `28.999…`).
fn whole_pct(frac: f64) -> u64 {
    (frac * 100.0 + 1e-9).floor() as u64
}

impl Inner {
    /// Write one whole line (`body` is the JSON between the common prefix and
    /// `elapsed_s`). The time is taken under the lock, so `elapsed_s` never
    /// decreases. Write errors are ignored: progress must never fail a run.
    fn emit(&self, st: &mut State, body: &str) {
        let elapsed = self.t0.elapsed().as_secs_f64();
        let line = format!(
            "{{\"progress\":{PROGRESS_VERSION},\"cmd\":\"{}\",{body},\"elapsed_s\":{elapsed:.2}}}\n",
            self.cmd
        );
        let _ = st.out.write_all(line.as_bytes());
        let _ = st.out.flush();
    }
}

/// Fraction of a chunked record pass: bp of finished chunks / total bp, plus
/// the records read in them. Shared by the worker threads.
pub struct ChunkFrac<'a> {
    p: &'a Progress,
    total: u64,
    done: AtomicU64,
    records: AtomicU64,
}

impl<'a> ChunkFrac<'a> {
    pub fn new(p: &'a Progress, total: u64) -> Self {
        ChunkFrac { p, total, done: AtomicU64::new(0), records: AtomicU64::new(0) }
    }

    /// A chunk of `bp` reference bases holding `records` records is finished
    /// (the unplaced task: `bp` 0).
    pub fn add(&self, bp: u64, records: u64) {
        if self.p.is_on() && self.total > 0 {
            let r = self.records.fetch_add(records, Ordering::Relaxed) + records;
            let d = self.done.fetch_add(bp, Ordering::Relaxed) + bp;
            self.p.update_records(d as f64 / self.total as f64, r);
        }
    }
}

/// Fraction of a pass over a coordinate-sorted input: reference position
/// (earlier contigs + pos) / total reference length, plus the records read.
/// With an exact record total (set by `start_records` before the pass) the
/// fraction is the larger of that and records / total, so the unplaced tail
/// also advances it.
pub struct RefFrac<'a> {
    p: &'a Progress,
    offsets: Vec<u64>,
    total: u64,
    /// Reference bp at which the next whole percent is reached.
    next: u64,
    /// Last reference position seen.
    bp: u64,
    records: u64,
    /// Exact record total of the phase, if known.
    rec_total: Option<u64>,
    /// Record count at which the next whole percent is reached.
    next_rec: u64,
}

impl<'a> RefFrac<'a> {
    /// `None` when `p` is disabled or the reference is empty.
    pub fn new(p: &'a Progress, lens: &[u64]) -> Option<Self> {
        if !p.is_on() {
            return None;
        }
        let mut offsets = Vec::with_capacity(lens.len());
        let mut total = 0u64;
        for &l in lens {
            offsets.push(total);
            total += l;
        }
        if total == 0 {
            return None;
        }
        let rec_total = p.phase_total().filter(|&n| n > 0);
        let mut f = RefFrac { p, offsets, total, next: 0, bp: 0, records: 0, rec_total, next_rec: 0 };
        f.set_next(0);
        Some(f)
    }

    pub fn from_header(p: &'a Progress, h: &HeaderView) -> Option<Self> {
        if !p.is_on() {
            return None;
        }
        let lens: Vec<u64> = (0..h.target_count()).map(|t| h.target_len(t).unwrap_or(0)).collect();
        Self::new(p, &lens)
    }

    /// Thresholds (reference bp, records) of the percent after `pct`.
    fn set_next(&mut self, pct: u64) {
        self.next = ((pct + 1) * self.total).div_ceil(100);
        self.next_rec = match self.rec_total {
            Some(n) => ((pct + 1) * n).div_ceil(100),
            None => u64::MAX,
        };
    }

    /// A record at (`tid`, `pos`) was read.
    #[inline]
    pub fn record(&mut self, tid: i32, pos: i64) {
        self.records += 1;
        let mut hit = self.records >= self.next_rec;
        if tid >= 0
            && pos >= 0
            && let Some(&off) = self.offsets.get(tid as usize)
        {
            self.bp = off + pos as u64;
            hit |= self.bp >= self.next;
        }
        if hit {
            let mut frac = self.bp as f64 / self.total as f64;
            if let Some(n) = self.rec_total {
                frac = frac.max(self.records as f64 / n as f64);
            }
            self.p.update_records(frac, self.records);
            self.set_next(whole_pct(frac.min(1.0)));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Default)]
    struct Buf(Arc<Mutex<Vec<u8>>>);
    impl Write for Buf {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(b);
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    impl Buf {
        fn lines(&self) -> Vec<serde_json::Value> {
            let s = String::from_utf8(self.0.lock().unwrap().clone()).unwrap();
            s.lines().map(|l| serde_json::from_str(l).unwrap()).collect()
        }
        fn fracs(&self) -> Vec<f64> {
            self.lines().iter().filter_map(|v| v["frac"].as_f64()).collect()
        }
        fn events(&self) -> Vec<String> {
            self.lines().iter().map(|v| v["event"].as_str().unwrap().to_string()).collect()
        }
    }

    fn capture(cmd: &'static str) -> (Progress, Buf) {
        let b = Buf::default();
        (Progress::to_writer(cmd, Box::new(b.clone())), b)
    }

    #[test]
    fn off_is_noop() {
        let p = Progress::off();
        assert!(!p.is_on());
        p.plan(&["a"]);
        p.start("a");
        p.update(0.5);
        p.update_records(0.5, 3);
        p.end();
        p.done();
        assert!(RefFrac::new(&p, &[100]).is_none());
        ChunkFrac::new(&p, 10).add(5, 1);
        assert_eq!(p.record_total(std::path::Path::new("/nonexistent")), None);
    }

    #[test]
    fn line_format() {
        let (p, b) = capture("qc");
        p.plan(&["scan", "call"]);
        p.start("scan");
        p.update(0.4216);
        p.end();
        p.start_records("call", Some(1000));
        p.update_records(0.5, 480);
        p.end();
        p.error("bad \"x\"\nline");
        p.done();
        let raw = String::from_utf8(b.0.lock().unwrap().clone()).unwrap();
        let first = raw.lines().next().unwrap();
        assert!(first.starts_with("{\"progress\":1,\"cmd\":\"qc\",\"event\":\"plan\",\"phases\":[\"scan\",\"call\"],\"elapsed_s\":"), "{first}");
        let l = b.lines();
        assert_eq!(l.len(), 9);
        assert_eq!(l[1]["phase"], "scan");
        assert_eq!(l[1]["event"], "start");
        assert_eq!(l[2]["event"], "update");
        assert_eq!(l[2]["frac"].as_f64().unwrap(), 0.421); // truncated
        assert!(l[2].get("records").is_none());
        assert_eq!(l[3]["event"], "end");
        assert_eq!(l[5]["records"], 480);
        assert_eq!(l[5]["records_total"], 1000);
        assert_eq!(l[5]["total_exact"], true);
        assert_eq!(l[7]["message"], "bad \"x\"\nline");
        assert_eq!(l[8]["event"], "done");
        assert!(l.iter().all(|v| v["progress"] == 1 && v["cmd"] == "qc" && v["elapsed_s"].is_number()));
    }

    #[test]
    fn whole_percent_steps_and_monotonic() {
        let (p, b) = capture("lift");
        p.update(0.5); // no phase: dropped
        p.start("pass");
        p.update(0.009); // 0 %
        p.update(0.10); // 10 %
        p.update(0.105); // still 10 %
        p.update(0.05); // backwards
        p.update(0.119); // 11 %
        p.update(0.43); // jump 11 -> 43 %: one line
        p.update(0.99975); // 99 %, shown truncated
        p.update(7.0); // clamped to 100 %
        p.update(1.0); // 100 % again
        p.update(f64::NAN);
        p.start("pass"); // already running: no-op
        p.end();
        p.end(); // nothing running
        assert_eq!(b.fracs(), vec![0.1, 0.119, 0.43, 0.999, 1.0]);
        assert_eq!(b.events(), ["start", "update", "update", "update", "update", "update", "end"]);
        // A new phase starts again at 0 %.
        p.start("next");
        p.update(0.02);
        assert_eq!(b.fracs().last(), Some(&0.02));
    }

    #[test]
    fn estimated_total_without_exact() {
        let (p, b) = capture("qc");
        p.start_records("linkage", None);
        p.update_records(0.25, 100);
        p.update_records(0.5, 90); // records never go back
        let l = b.lines();
        assert_eq!(l[1]["records_total"], 400);
        assert_eq!(l[1]["total_exact"], false);
        assert_eq!(l[2]["records"], 100);
        assert_eq!(l[2]["records_total"], 200);
    }

    #[test]
    fn exact_total_records_advance_and_gate_100() {
        let (p, b) = capture("qc");
        p.start_records("scan", Some(100));
        p.update_records(0.2, 50); // records ahead of the position
        p.update_records(1.0, 80); // all positions done, records not: 99 %
        p.update_records(1.0, 100); // every record read: 100 %
        let l = b.lines();
        let fr = b.fracs();
        assert_eq!(fr, vec![0.5, 0.999, 1.0]);
        assert_eq!(l[3]["records"], 100);
    }

    #[test]
    fn ref_frac_unplaced_tail_reaches_100() {
        // 100 placed records over the reference, then 100 unplaced unmapped.
        let (p, b) = capture("lift");
        p.start_records("translate_pass", Some(200));
        let mut rf = RefFrac::new(&p, &[1000]).unwrap();
        for i in 0..100 {
            rf.record(0, i * 10);
        }
        let placed = b.fracs().len();
        for _ in 0..100 {
            rf.record(-1, -1);
        }
        p.end();
        let l = b.lines();
        let upd: Vec<&serde_json::Value> = l.iter().filter(|v| v["event"] == "update").collect();
        assert!(upd.len() > placed, "the tail did not advance");
        let pct: Vec<u64> = upd.iter().map(|v| whole_pct(v["frac"].as_f64().unwrap())).collect();
        assert!(pct.windows(2).all(|w| w[0] < w[1]), "{pct:?}");
        assert!(upd.len() <= 100);
        let last = upd.last().unwrap();
        assert_eq!(last["frac"].as_f64().unwrap(), 1.0);
        assert_eq!(last["records"], 200);
        assert_eq!(last["records_total"], 200);
    }

    #[test]
    fn exceeded_exact_total_switches_to_estimate() {
        let (p, b) = capture("qc");
        p.start_records("scan", Some(0)); // e.g. an index reporting 0
        p.update_records(0.5, 10);
        p.update_records(0.8, 30);
        p.end();
        p.start_records("pass", Some(20));
        p.update_records(0.3, 10);
        p.update_records(0.6, 25); // more than the "exact" total
        p.update_records(0.9, 28);
        let u: Vec<serde_json::Value> = b.lines().into_iter().filter(|v| v["event"] == "update").collect();
        let got: Vec<(u64, u64, bool)> = u
            .iter()
            .map(|v| (v["records"].as_u64().unwrap(), v["records_total"].as_u64().unwrap(), v["total_exact"].as_bool().unwrap()))
            .collect();
        assert_eq!(got, vec![(10, 20, false), (30, 38, false), (10, 20, true), (25, 42, false), (28, 31, false)]);
        assert!(got.iter().all(|&(r, n, exact)| !exact || r <= n));
    }

    #[test]
    fn done_closes_open_phase() {
        let (p, b) = capture("fish");
        p.start("pass");
        p.done();
        assert_eq!(b.events(), ["start", "end", "done"]);
    }

    #[test]
    fn ref_frac_positions() {
        let (p, b) = capture("lift");
        p.start_records("translate_pass", Some(2000));
        let mut rf = RefFrac::new(&p, &[600, 400]).unwrap();
        rf.record(-1, 5); // counted, no position
        rf.record(1, 200); // (600 + 200) / 1000
        rf.record(1, 205); // still 80 %
        rf.record(1, 211); // 81 %
        assert_eq!(b.fracs(), vec![0.8, 0.811]);
        let recs: Vec<u64> = b.lines().iter().filter_map(|v| v["records"].as_u64()).collect();
        assert_eq!(recs, vec![2, 4]);
        assert!(RefFrac::new(&p, &[0, 0]).is_none());

        // Every record of a dense pass: one line per whole percent.
        let (p, b) = capture("lift");
        p.start_records("translate_pass", None);
        let mut rf = RefFrac::new(&p, &[1000]).unwrap();
        for pos in 0..1000 {
            rf.record(0, pos);
        }
        let fr = b.fracs();
        assert_eq!(fr.len(), 99);
        assert_eq!(fr[0], 0.01);
        assert_eq!(fr[98], 0.99);
    }

    #[test]
    fn chunk_frac_threads() {
        let (p, b) = capture("qc");
        p.start_records("scan", Some(400));
        let cf = ChunkFrac::new(&p, 1000);
        std::thread::scope(|s| {
            for _ in 0..4 {
                s.spawn(|| {
                    for _ in 0..25 {
                        cf.add(10, 4);
                    }
                });
            }
        });
        p.end();
        let fr = b.fracs();
        let pct: Vec<u64> = fr.iter().map(|&f| whole_pct(f)).collect();
        assert!(fr.len() >= 10, "too few updates: {fr:?}");
        assert!(pct.windows(2).all(|w| w[0] < w[1]), "whole percents not strictly increasing: {pct:?}");
        assert!(fr.iter().all(|f| (0.0..=1.0).contains(f)));
        let recs: Vec<u64> = b.lines().iter().filter_map(|v| v["records"].as_u64()).collect();
        assert!(recs.windows(2).all(|w| w[0] <= w[1]) && *recs.last().unwrap() <= 400, "{recs:?}");
        assert_eq!(cf.done.load(Ordering::Relaxed), 1000);
        assert_eq!(cf.records.load(Ordering::Relaxed), 400);
    }
}
