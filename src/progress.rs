//! `--progress`: machine-readable progress on stderr, one JSON object per line
//! (see docs/output-format.md, "Progress stream").
//!
//! A [`Progress`] handle is cheap to clone and safe to share between threads.
//! When disabled ([`Progress::off`], the default) every method is a no-op
//! behind one `Option` check, so the outputs are unchanged.

use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

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
/// At most one `update` line per phase in this interval (~5 per second).
/// `CRAMHOOK_PROGRESS_INTERVAL_MS` overrides it (tests only).
const MIN_INTERVAL: Duration = Duration::from_millis(200);
/// Minimum `frac` increase between two `update` lines of a phase.
const MIN_STEP: f64 = 0.01;

/// Progress reporter; a no-op when disabled.
#[derive(Clone, Default)]
pub struct Progress(Option<Arc<Inner>>);

struct Inner {
    cmd: &'static str,
    t0: Instant,
    interval: Duration,
    state: Mutex<State>,
}

struct State {
    out: Box<dyn Write + Send>,
    phase: Option<&'static str>,
    frac: f64,
    last: Instant,
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
        let interval = std::env::var("CRAMHOOK_PROGRESS_INTERVAL_MS")
            .ok()
            .and_then(|v| v.parse().ok())
            .map(Duration::from_millis)
            .unwrap_or(MIN_INTERVAL);
        Self::with_interval(cmd, out, interval)
    }

    fn with_interval(cmd: &'static str, out: Box<dyn Write + Send>, interval: Duration) -> Self {
        let t0 = Instant::now();
        Progress(Some(Arc::new(Inner {
            cmd,
            t0,
            interval,
            state: Mutex::new(State { out, phase: None, frac: 0.0, last: t0 }),
        })))
    }

    pub fn is_on(&self) -> bool {
        self.0.is_some()
    }

    /// First line: the phases in the order they run.
    pub fn plan(&self, phases: &[&str]) {
        let Some(i) = &self.0 else { return };
        let list: Vec<String> = phases.iter().map(|p| format!("\"{p}\"")).collect();
        let mut st = i.state.lock().unwrap();
        i.emit(&mut st, &format!("\"event\":\"plan\",\"phases\":[{}]", list.join(",")), None);
    }

    /// Start `phase`. Starting the phase that is already running is a no-op
    /// (a restarted pass continues its phase; `frac` stays monotonic).
    pub fn start(&self, phase: &'static str) {
        let Some(i) = &self.0 else { return };
        let mut st = i.state.lock().unwrap();
        if st.phase == Some(phase) {
            return;
        }
        st.phase = Some(phase);
        st.frac = 0.0;
        st.last = Instant::now();
        i.emit(&mut st, &format!("\"phase\":\"{phase}\",\"event\":\"start\""), None);
    }

    /// Report the fraction done of the running phase (clamped to [0, 1]).
    /// Rate limited; lines that would not advance `frac` are dropped.
    pub fn update(&self, frac: f64) {
        if self.0.is_some() {
            self.update_at(frac, None);
        }
    }

    /// `now` defaults to the time the lock is taken, so `elapsed_s` never
    /// decreases between lines written by different threads.
    fn update_at(&self, frac: f64, now: Option<Instant>) {
        let Some(i) = &self.0 else { return };
        let frac = if frac.is_finite() { frac.clamp(0.0, 1.0) } else { return };
        let mut st = i.state.lock().unwrap();
        let Some(phase) = st.phase else { return };
        let now = now.unwrap_or_else(Instant::now);
        if frac < st.frac + MIN_STEP || now.saturating_duration_since(st.last) < i.interval {
            return;
        }
        st.frac = frac;
        st.last = now;
        i.emit(
            &mut st,
            &format!("\"phase\":\"{phase}\",\"event\":\"update\",\"frac\":{frac:.3}"),
            Some(now),
        );
    }

    /// End the running phase.
    pub fn end(&self) {
        let Some(i) = &self.0 else { return };
        let mut st = i.state.lock().unwrap();
        if let Some(phase) = st.phase.take() {
            i.emit(&mut st, &format!("\"phase\":\"{phase}\",\"event\":\"end\""), None);
        }
    }

    /// Last line of a successful run (ends a phase still running).
    pub fn done(&self) {
        self.end();
        let Some(i) = &self.0 else { return };
        let mut st = i.state.lock().unwrap();
        i.emit(&mut st, "\"event\":\"done\"", None);
    }

    /// Last line of a failed run; the error is still reported as usual.
    pub fn error(&self, message: &str) {
        let Some(i) = &self.0 else { return };
        let mut st = i.state.lock().unwrap();
        let msg = serde_json::to_string(message).unwrap_or_else(|_| "\"\"".to_string());
        i.emit(&mut st, &format!("\"event\":\"error\",\"message\":{msg}"), None);
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

impl Inner {
    /// Write one whole line (`body` is the JSON between the common prefix and
    /// `elapsed_s`). Write errors are ignored: progress must never fail a run.
    fn emit(&self, st: &mut State, body: &str, now: Option<Instant>) {
        let now = now.unwrap_or_else(Instant::now);
        let elapsed = now.saturating_duration_since(self.t0).as_secs_f64();
        let line = format!(
            "{{\"progress\":{PROGRESS_VERSION},\"cmd\":\"{}\",{body},\"elapsed_s\":{elapsed:.2}}}\n",
            self.cmd
        );
        let _ = st.out.write_all(line.as_bytes());
        let _ = st.out.flush();
    }
}

/// Fraction of a chunked pass: bp of finished chunks / total bp. Shared by the
/// worker threads.
pub struct ChunkFrac<'a> {
    p: &'a Progress,
    total: u64,
    done: AtomicU64,
}

impl<'a> ChunkFrac<'a> {
    pub fn new(p: &'a Progress, total: u64) -> Self {
        ChunkFrac { p, total, done: AtomicU64::new(0) }
    }

    /// A chunk of `bp` reference bases is finished.
    pub fn add(&self, bp: u64) {
        if self.p.is_on() && self.total > 0 {
            let d = self.done.fetch_add(bp, Ordering::Relaxed) + bp;
            self.p.update(d as f64 / self.total as f64);
        }
    }
}

/// Fraction of a pass over a coordinate-sorted input: reference position
/// (earlier contigs + pos) / total reference length. The unplaced tail does not
/// advance it; the phase `end` closes it.
pub struct RefFrac<'a> {
    p: &'a Progress,
    offsets: Vec<u64>,
    total: u64,
    step: u64,
    next: u64,
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
        Some(RefFrac { p, offsets, total, step: (total / 200).max(1), next: 0 })
    }

    pub fn from_header(p: &'a Progress, h: &HeaderView) -> Option<Self> {
        if !p.is_on() {
            return None;
        }
        let lens: Vec<u64> = (0..h.target_count()).map(|t| h.target_len(t).unwrap_or(0)).collect();
        Self::new(p, &lens)
    }

    /// A record at (`tid`, `pos`) was read.
    #[inline]
    pub fn record(&mut self, tid: i32, pos: i64) {
        if tid < 0 || pos < 0 {
            return;
        }
        let Some(&off) = self.offsets.get(tid as usize) else { return };
        let bp = off + pos as u64;
        if bp >= self.next {
            self.next = bp + self.step;
            self.p.update(bp as f64 / self.total as f64);
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
        p.end();
        p.done();
        assert!(RefFrac::new(&p, &[100]).is_none());
        ChunkFrac::new(&p, 10).add(5);
    }

    #[test]
    fn line_format() {
        let (p, b) = capture("qc");
        p.plan(&["scan", "call"]);
        p.start("scan");
        let t = p.0.as_ref().unwrap().state.lock().unwrap().last;
        p.update_at(0.4216, Some(t + Duration::from_secs(1)));
        p.end();
        p.error("bad \"x\"\nline");
        p.done();
        let raw = String::from_utf8(b.0.lock().unwrap().clone()).unwrap();
        let first = raw.lines().next().unwrap();
        assert!(first.starts_with("{\"progress\":1,\"cmd\":\"qc\",\"event\":\"plan\",\"phases\":[\"scan\",\"call\"],\"elapsed_s\":"), "{first}");
        let l = b.lines();
        assert_eq!(l.len(), 6);
        assert_eq!(l[1]["phase"], "scan");
        assert_eq!(l[1]["event"], "start");
        assert_eq!(l[2]["event"], "update");
        assert_eq!(l[2]["frac"].as_f64().unwrap(), 0.422);
        assert_eq!(l[3]["event"], "end");
        assert_eq!(l[4]["message"], "bad \"x\"\nline");
        assert_eq!(l[5]["event"], "done");
        assert!(l.iter().all(|v| v["progress"] == 1 && v["cmd"] == "qc" && v["elapsed_s"].is_number()));
    }

    #[test]
    fn rate_limit_and_monotonic() {
        let (p, b) = capture("lift");
        p.update(0.5); // no phase: dropped
        p.start("pass");
        let t = p.0.as_ref().unwrap().state.lock().unwrap().last;
        let ms = |m: u64| Some(t + Duration::from_millis(m));
        p.update_at(0.10, ms(100)); // too soon after start
        p.update_at(0.10, ms(250)); // emitted
        p.update_at(0.50, ms(300)); // too soon
        p.update_at(0.105, ms(600)); // step < 0.01
        p.update_at(0.05, ms(700)); // backwards
        p.update_at(0.12, ms(800)); // emitted
        p.update_at(7.0, ms(1100)); // clamped to 1
        p.update_at(f64::NAN, ms(1400));
        p.start("pass"); // already running: no-op
        p.end();
        p.end(); // nothing running
        let fr: Vec<f64> = b.lines().iter().filter_map(|v| v["frac"].as_f64()).collect();
        assert_eq!(fr, vec![0.1, 0.12, 1.0]);
        let ev: Vec<String> = b.lines().iter().map(|v| v["event"].as_str().unwrap().to_string()).collect();
        assert_eq!(ev, ["start", "update", "update", "update", "end"]);
    }

    #[test]
    fn done_closes_open_phase() {
        let (p, b) = capture("fish");
        p.start("pass");
        p.done();
        let ev: Vec<String> = b.lines().iter().map(|v| v["event"].as_str().unwrap().to_string()).collect();
        assert_eq!(ev, ["start", "end", "done"]);
    }

    #[test]
    fn ref_frac_positions() {
        let (p, b) = capture("qc");
        p.start("pass");
        p.0.as_ref().unwrap().state.lock().unwrap().last = Instant::now() - Duration::from_secs(10);
        let mut rf = RefFrac::new(&p, &[600, 400]).unwrap();
        rf.record(-1, 5);
        rf.record(1, 200); // (600 + 200) / 1000
        let fr: Vec<f64> = b.lines().iter().filter_map(|v| v["frac"].as_f64()).collect();
        assert_eq!(fr, vec![0.8]);
        assert!(RefFrac::new(&p, &[0, 0]).is_none());
    }

    #[test]
    fn chunk_frac_threads() {
        // No rate limit, so the fast loop below must emit updates.
        let b = Buf::default();
        let p = Progress::with_interval("qc", Box::new(b.clone()), Duration::ZERO);
        p.start("scan");
        let cf = ChunkFrac::new(&p, 1000);
        std::thread::scope(|s| {
            for _ in 0..4 {
                s.spawn(|| {
                    for _ in 0..25 {
                        cf.add(10);
                    }
                });
            }
        });
        p.end();
        let fr: Vec<f64> = b.lines().iter().filter_map(|v| v["frac"].as_f64()).collect();
        assert!(fr.len() >= 10, "too few updates: {fr:?}");
        assert!(fr.windows(2).all(|w| w[0] < w[1]), "{fr:?}");
        assert!(fr.iter().all(|f| (0.0..=1.0).contains(f)));
        assert_eq!(cf.done.load(Ordering::Relaxed), 1000);
    }
}
