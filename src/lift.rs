//! M4: the full `cramhook lift` pipeline.
//!
//! At t=0 the remap `minimap2` and the merging `samtools sort` are spawned (so
//! minimap2's index build overlaps the assembly diff). Then a single streaming
//! pass writes translated records (uncompressed BAM) into sort's stdin and
//! deferred reads (FASTQ) into minimap2's stdin, while a reader thread drains
//! minimap2's SAM stdout into the same BAM writer. Records are tagged `XL:A:T`
//! (translated) / `XL:A:R` (remapped) by default; `--no-xl-tag` disables it.
//!
//! The output is staged (`OUT.tmp.cram` + `.tmp.cram.crai`) and renamed into
//! place only after all checks pass, so a failed run leaves no partial output
//! and never destroys a previous good one. Every failure path kills and reaps
//! both children, joins the reader thread and reports their exit statuses.
//!
//! Speed: the CRAM reader uses an htslib thread pool, the intermediate BAM is
//! uncompressed, and the diff/translate/remap phases overlap.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Read as _, Seek as _, SeekFrom, Write};
use std::os::unix::fs::OpenOptionsExt;
use std::os::unix::io::AsRawFd;
use std::os::unix::process::ExitStatusExt;
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, ChildStdout, Command, ExitStatus, Stdio};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result, anyhow, bail};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Format, Read};
use rust_htslib::tpool::ThreadPool;

use crate::contract::{
    check_contract, check_sq_against_fai, parse_minimap2_cl, read_fai, reject_old_name_lift,
};
use crate::diff::{Diff, Preset, compute_diff_with, write_outputs};
use crate::fasta::{Fasta, ensure_fai};
use crate::version::versioned;
use crate::translate::{
    Ctx, Fallback, LiftSink, Params, PassOutcome, Reason, Stats, TranslateConfig, build_output_header,
    de_threshold, ensure_parent, fai_path, open_reader, run_all_defer, run_translate_pass, scoring,
    suffix,
};

/// Options for `cramhook lift` beyond the shared translate config.
#[derive(Clone, Copy, Debug)]
pub struct LiftOptions {
    pub strict: bool,
    pub keep_tmp: bool,
    /// M9: tiny-input fallback threshold in MB (10^6 bytes). `None` uses the
    /// per-read-preset default (lr:hq 20, map-hifi 3); `Some(0.0)` disables the
    /// fallback. Overrides via `--min-input-mb`.
    pub min_input_mb: Option<f64>,
}

/// Default tiny-input threshold in MB (10^6 bytes) for a read preset. A full
/// remap of a CRAM this small costs less than the fixed lift overhead
/// (diff + index), so lifting is not worth it.
fn default_min_input_mb(preset: Preset) -> f64 {
    match preset {
        Preset::LrHq => 20.0,
        Preset::MapHifi => 3.0,
    }
}

/// Per-phase thread split: (CRAM decode, remap minimap2, sort, diff minimap2).
/// Overridable for experiments with `CRAMHOOK_THREAD_SPLIT=decode,remap,sort,diff`.
fn split_threads(n: usize) -> (usize, usize, usize, usize) {
    if let Ok(v) = std::env::var("CRAMHOOK_THREAD_SPLIT") {
        let t: Vec<usize> = v.split(',').filter_map(|x| x.trim().parse().ok()).collect();
        if t.len() == 4 {
            return (t[0].max(1), t[1].max(1), t[2].max(1), t[3].max(1));
        }
    }
    let n = n.max(1);
    // The remap minimap2 dominates; give it all threads and split the rest.
    ((n / 2).max(1), n, (n / 2).max(1), (n / 2).max(1))
}

fn executable(env: &str, default: &str) -> String {
    std::env::var(env).unwrap_or_else(|_| default.to_string())
}

fn fd_path(fd: i32) -> PathBuf {
    PathBuf::from(format!("/proc/self/fd/{fd}"))
}

/// The actual staging paths used during a run: `OUT.tmp.cram` and its index.
/// Public so tests assert the real paths instead of hard-coding them.
pub fn staged_paths(out: &Path) -> (PathBuf, PathBuf) {
    let mut s = out.as_os_str().to_owned();
    s.push(".tmp.cram");
    let staged = PathBuf::from(s);
    let mut i = staged.as_os_str().to_owned();
    i.push(".crai");
    (staged, PathBuf::from(i))
}

fn status_desc(s: &ExitStatus) -> String {
    if let Some(code) = s.code() {
        format!("status {code}")
    } else if let Some(sig) = s.signal() {
        format!("signal {sig}")
    } else {
        format!("{s}")
    }
}

// ---------------------------------------------------------------------------
// reader thread

/// Shared slot for a reader-thread failure, checked by the sink so the main
/// pass stops promptly instead of hanging.
type SharedError = Arc<Mutex<Option<String>>>;

struct RemapResult {
    primary_names: HashSet<Vec<u8>>,
    primary_count: u64,
}

/// Drain minimap2's SAM stdout into the BAM writer. Owns the stdout pipe so it
/// is closed when the thread exits (even on error), which unblocks minimap2.
fn spawn_reader(
    mm_stdout: ChildStdout,
    writer: Arc<Mutex<bam::Writer>>,
    new_names: Vec<String>,
    xl_tag: bool,
    shared_err: SharedError,
) -> JoinHandle<Result<RemapResult>> {
    let fd = mm_stdout.as_raw_fd();
    thread::spawn(move || {
        let _own = mm_stdout; // closed when this closure returns
        let outcome = (|| -> Result<RemapResult> {
            let mut rdr = bam::Reader::from_path(fd_path(fd))
                .context("opening minimap2 SAM stdout")?;
            let names: Vec<String> = rdr
                .header()
                .target_names()
                .iter()
                .map(|n| String::from_utf8_lossy(n).to_string())
                .collect();
            if names != new_names {
                bail!("minimap2 remap @SQ order differs from new.fa .fai: {names:?} vs {new_names:?}");
            }
            let mut rec = bam::Record::new();
            let mut primary_names = HashSet::new();
            let mut primary_count = 0u64;
            while let Some(res) = rdr.read(&mut rec) {
                res.context("reading minimap2 SAM")?;
                if xl_tag {
                    let _ = rec.push_aux(b"XL", Aux::Char(b'R'));
                }
                if rec.flags() & (0x100 | 0x800) == 0 {
                    primary_count += 1;
                    primary_names.insert(rec.qname().to_vec());
                }
                writer.lock().unwrap().write(&rec)?;
            }
            Ok(RemapResult { primary_names, primary_count })
        })();
        if let Err(e) = &outcome {
            *shared_err.lock().unwrap() = Some(format!("{e:#}"));
        }
        outcome
    })
}

// ---------------------------------------------------------------------------
// sink

#[derive(Default)]
struct SinkState {
    deferred_names: HashSet<Vec<u8>>,
}

struct MergeSink {
    writer: Arc<Mutex<bam::Writer>>,
    fq: Option<BufWriter<ChildStdin>>,
    tee: Option<BufWriter<File>>,
    tsv: Option<BufWriter<File>>,
    state: Arc<Mutex<SinkState>>,
    shared_err: SharedError,
}

fn write_fastq_entry(w: &mut impl Write, qname: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    w.write_all(b"@")?;
    w.write_all(qname)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")?;
    Ok(())
}

impl MergeSink {
    fn check_err(&self) -> Result<()> {
        if let Some(e) = self.shared_err.lock().unwrap().as_ref() {
            bail!("minimap2 reader failed: {e}");
        }
        Ok(())
    }
}

impl LiftSink for MergeSink {
    fn translated(&mut self, rec: &bam::Record) -> Result<()> {
        self.check_err()?;
        let mut w = self.writer.lock().unwrap();
        w.write(rec)?;
        Ok(())
    }
    fn fastq(&mut self, qname: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
        self.check_err()?;
        if let Some(f) = self.fq.as_mut() {
            write_fastq_entry(f, qname, seq, qual)?;
        }
        if let Some(f) = self.tee.as_mut() {
            write_fastq_entry(f, qname, seq, qual)?;
        }
        Ok(())
    }
    fn deferred(&mut self, qname: &[u8], reason: Reason) -> Result<()> {
        self.state.lock().unwrap().deferred_names.insert(qname.to_vec());
        if let Some(t) = self.tsv.as_mut() {
            let mut row = Vec::new();
            row.extend_from_slice(qname);
            row.push(b'\t');
            row.extend_from_slice(reason.name().as_bytes());
            row.push(b'\t');
            row.extend_from_slice(reason.class().name().as_bytes());
            row.push(b'\n');
            t.write_all(&row)?;
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<()> {
        if let Some(mut f) = self.fq.take() {
            f.flush()?;
            drop(f); // closes minimap2's stdin
        }
        if let Some(mut f) = self.tee.take() {
            f.flush()?;
        }
        if let Some(mut t) = self.tsv.take() {
            t.flush()?;
        }
        Ok(())
    }
}

/// cutime + cstime from /proc/self/stat (USER_HZ = 100 on Linux).
fn children_cpu_s() -> f64 {
    let Ok(s) = std::fs::read_to_string("/proc/self/stat") else {
        return 0.0;
    };
    let Some(pos) = s.rfind(')') else { return 0.0 };
    let rest: Vec<&str> = s[pos + 1..].split_whitespace().collect();
    let get = |i: usize| rest.get(i).and_then(|x| x.parse::<u64>().ok()).unwrap_or(0);
    (get(13) + get(14)) as f64 / 100.0
}

// ---------------------------------------------------------------------------
// children

struct Children {
    mm: Option<Child>,
    sort: Option<Child>,
    mm_stdin: Option<ChildStdin>,
    mm_stdout: Option<ChildStdout>,
    sort_stdin: Option<ChildStdin>,
    t0_mm: Instant,
    finished: bool,
    mm_stderr_path: Option<PathBuf>,
    sort_stderr_path: Option<PathBuf>,
    /// Captured child stderr text, filled in by `reap`/`wait`.
    mm_stderr: Option<String>,
    sort_stderr: Option<String>,
}

impl Children {
    fn empty() -> Children {
        Children {
            mm: None,
            sort: None,
            mm_stdin: None,
            mm_stdout: None,
            sort_stdin: None,
            t0_mm: Instant::now(),
            finished: false,
            mm_stderr_path: None,
            sort_stderr_path: None,
            mm_stderr: None,
            sort_stderr: None,
        }
    }

    /// Kill and reap both children, returning their exit statuses.
    fn reap(&mut self) -> (Option<ExitStatus>, Option<ExitStatus>) {
        self.finished = true;
        drop(self.mm_stdin.take());
        drop(self.sort_stdin.take());
        let mm = self.mm.as_mut().and_then(|c| {
            let _ = c.kill();
            c.wait().ok()
        });
        let sort = self.sort.as_mut().and_then(|c| {
            let _ = c.kill();
            c.wait().ok()
        });
        self.mm_stderr = read_stderr(&self.mm_stderr_path);
        self.sort_stderr = read_stderr(&self.sort_stderr_path);
        (mm, sort)
    }

    /// Close stdin and wait for both children to exit on their own (no kill).
    fn wait(&mut self) -> (Option<ExitStatus>, Option<ExitStatus>) {
        self.finished = true;
        drop(self.mm_stdin.take());
        drop(self.sort_stdin.take());
        let mm = self.mm.as_mut().and_then(|c| c.wait().ok());
        let sort = self.sort.as_mut().and_then(|c| c.wait().ok());
        self.mm_stderr = read_stderr(&self.mm_stderr_path);
        self.sort_stderr = read_stderr(&self.sort_stderr_path);
        (mm, sort)
    }
}

impl Drop for Children {
    /// Safety net for any `?` path that never reached `fail`/`wait`: never leave
    /// a spawned child running or unreaped.
    fn drop(&mut self) {
        if !self.finished {
            let _ = self.reap();
        }
    }
}

/// Temp file for a child's stderr. A file (not a pipe) needs no reader thread and
/// can never block the child on a full pipe; the text is read back after exit.
///
/// `create_new` never follows a symlink at the target and never truncates an
/// existing file, so a planted symlink in the temp dir cannot redirect the write;
/// mode `0o600` keeps the child's output private. A name collision just retries.
fn stderr_capture_file(tag: &str) -> Result<(File, PathBuf)> {
    let dir = std::env::temp_dir();
    for attempt in 0..16u32 {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let path = dir.join(format!(
            "cramhook-{tag}-{}-{nanos}-{attempt}.stderr",
            std::process::id()
        ));
        match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .mode(0o600)
            .open(&path)
        {
            Ok(f) => return Ok((f, path)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(anyhow!("creating {}: {e}", path.display())),
        }
    }
    bail!("could not create a unique stderr capture file in {}", dir.display())
}

/// Best-effort: an `Stdio` the child writes its stderr to, plus the file to read
/// back. Falls back to `inherit()` when no temp file can be created, so a run is
/// never failed just because capture is unavailable. `None` path = not captured.
fn capture_stderr(tag: &str) -> (Stdio, Option<PathBuf>) {
    match stderr_capture_file(tag) {
        Ok((f, p)) => (Stdio::from(f), Some(p)),
        Err(_) => (Stdio::inherit(), None),
    }
}

/// Read and remove a captured stderr file. `None` when absent/unreadable/empty.
///
/// Bounded and lossy: a child that floods stderr cannot make this allocate
/// without limit, and one non-UTF-8 byte cannot discard the whole diagnostic.
/// The result is a bounded *head* followed by a bounded *tail* of the file. The
/// head keeps the index-ready line [`parse_index_ready_s`] looks for; the tail
/// keeps the real error text of a failing child. Real minimap2/samtools output is
/// orders of magnitude smaller than either cap, so it is returned verbatim.
fn read_stderr(path: &Option<PathBuf>) -> Option<String> {
    const HEAD_CAP: u64 = 1 << 20; // 1 MiB
    const TAIL_CAP: u64 = 1 << 16; // 64 KiB
    let p = path.as_ref()?;
    let mut f = match File::open(p) {
        Ok(f) => f,
        Err(_) => {
            let _ = std::fs::remove_file(p);
            return None;
        }
    };
    let total = f.metadata().map(|m| m.len()).unwrap_or(0);
    let mut head = Vec::new();
    if (&mut f).take(HEAD_CAP).read_to_end(&mut head).is_err() {
        let _ = std::fs::remove_file(p);
        return None;
    }
    let mut text = String::from_utf8_lossy(&head).into_owned();
    if total > HEAD_CAP {
        // Read the tail too (the failure snippet takes the end of the string).
        // Never start before the head ends, so head and tail cannot overlap.
        let start = total.saturating_sub(TAIL_CAP).max(head.len() as u64);
        let mut tail = Vec::new();
        let tail_ok = if f.seek(SeekFrom::Start(start)).is_ok() {
            f.take(TAIL_CAP).read_to_end(&mut tail).is_ok()
        } else {
            false
        };
        if tail_ok {
            if start > head.len() as u64 {
                text.push_str(&format!("\n... [{} bytes omitted] ...\n", start - head.len() as u64));
            }
            text.push_str(&String::from_utf8_lossy(&tail));
        }
    }
    let _ = std::fs::remove_file(p);
    Some(text).filter(|t| !t.trim().is_empty())
}

/// Char-boundary-safe tail of `s` (never panics on UTF-8).
fn tail_chars(s: &str, max: usize) -> &str {
    if s.len() <= max {
        return s;
    }
    let mut start = s.len() - max;
    while !s.is_char_boundary(start) {
        start += 1;
    }
    &s[start..]
}

/// Captured child stderr, formatted for an error message: `": <text>"`, or `""`
/// when there is nothing to show. Capped so a chatty child can't flood the error.
fn stderr_snippet(s: Option<&str>) -> String {
    match s.map(str::trim).filter(|t| !t.is_empty()) {
        Some(t) => format!(": {}", tail_chars(t, 2000)),
        None => String::new(),
    }
}

/// minimap2 logs index build completion to stderr as
/// `[M::main::<secs>*<rate>] loaded/built the index for N target sequence(s)`.
/// `<secs>` is seconds since minimap2 started, i.e. the spawn→index-ready time
/// that is otherwise unobservable. Takes the last such line (split indexes with
/// `-I` can emit more than one).
fn parse_index_ready_s(stderr: Option<&str>) -> Option<f64> {
    let s = stderr?;
    let mut out = None;
    for line in s.lines() {
        if !line.contains("loaded/built the index") {
            continue;
        }
        let Some((_, rest)) = line.split_once("::") else { continue };
        let Some((_, rest)) = rest.split_once("::") else { continue };
        let num: String =
            rest.chars().take_while(|c| c.is_ascii_digit() || *c == '.').collect();
        if let Ok(v) = num.parse::<f64>() {
            out = Some(v);
        }
    }
    out
}

fn spawn_children(
    preset: Preset,
    new_ref: &Path,
    exe_mm: &str,
    exe_sort: &str,
    decode_t: usize,
    remap_t: usize,
    sort_t: usize,
    staged: &Path,
    staged_crai: &Path,
) -> Result<Children> {
    let _ = decode_t;
    // Holder created first: any failure below drops it, killing/reaping whatever
    // was already spawned.
    let mut ch = Children::empty();
    ch.t0_mm = Instant::now();
    // Set the path before the spawn so a failed spawn still cleans up the file.
    let (mm_err, mm_err_path) = capture_stderr("minimap2");
    ch.mm_stderr_path = mm_err_path;
    ch.mm = Some(
        Command::new(exe_mm)
            .args([
                "-ax",
                preset.minimap2_name(),
                "-t",
                &remap_t.to_string(),
                "--secondary=no",
                "-Y",
                &new_ref.display().to_string(),
                "-",
            ])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(mm_err)
            .spawn()
            .with_context(|| format!("spawning remap minimap2 ({exe_mm})"))?,
    );
    {
        ch.mm_stdin = Some(ch.mm.as_mut().unwrap().stdin.take().context("minimap2 stdin")?);
        ch.mm_stdout = Some(ch.mm.as_mut().unwrap().stdout.take().context("minimap2 stdout")?);
    }

    let mut o = staged.as_os_str().to_owned();
    o.push("##idx##");
    o.push(staged_crai.as_os_str());
    let o = PathBuf::from(o);
    let (sort_err, sort_err_path) = capture_stderr("samtools");
    ch.sort_stderr_path = sort_err_path;
    ch.sort = Some(
        Command::new(exe_sort)
            .args([
                "sort",
                "-@",
                &sort_t.to_string(),
                "-m",
                "1G",
                // The tail after minimap2 exits is dominated by CRAM encoding,
                // not by the (in-memory) sort merge or the index. CRAM level 1
                // roughly halves that encode tail for ~1% larger output; the
                // decoded records are byte-identical (see tests).
                "-l",
                "1",
                "-O",
                "CRAM",
                "--reference",
                &new_ref.display().to_string(),
                "--write-index",
                "-o",
                &o.display().to_string(),
                "-",
            ])
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(sort_err)
            .spawn()
            .with_context(|| format!("spawning samtools sort ({exe_sort})"))?,
    );
    ch.sort_stdin = Some(
        ch.sort
            .as_mut()
            .unwrap()
            .stdin
            .take()
            .context("sort stdin")?,
    );
    Ok(ch)
}

// ---------------------------------------------------------------------------
// session

enum RunMode<'a> {
    Normal { preset: Preset, diff: &'a Diff },
    AllDefer { preset: Preset, reason: Option<String> },
}

enum LiftOutcome {
    Done(Stats),
    MidStream(String),
}

struct Common<'a> {
    out: &'a Path,
    staged: &'a Path,
    staged_crai: &'a Path,
    header_bytes: &'a [u8],
    new_fai: &'a [(String, u64)],
    new_tids: &'a HashMap<String, i32>,
    old: &'a Fasta,
    new: &'a Fasta,
    old_names: &'a [String],
    tcfg: &'a TranslateConfig,
    diff_s: f64,
    t0: Instant,
    deferred_prefix: Option<&'a Path>,
    cram_input: bool,
    /// Resolved divergence threshold (`--max-de` or the preset default).
    max_de: f64,
}

impl Common<'_> {
    fn remove_staged(&self) {
        let _ = std::fs::remove_file(self.staged);
        let _ = std::fs::remove_file(self.staged_crai);
    }
}

/// One attempt. On any failure: kill+reap the children, join the reader thread,
/// remove the staged output, and return an error naming the child statuses.
fn run_session(
    common: &Common<'_>,
    mut children: Children,
    reader: &mut bam::Reader,
    mode: RunMode<'_>,
    decode_t: usize,
) -> Result<LiftOutcome> {
    // htslib thread pool for CRAM decode. `open_reader` already set `set_threads`;
    // sharing an explicit pool is what speeds up CRAM decode, but `set_thread_pool`
    // leaves a plain BAM reader at EOF, so only apply it to CRAM input.
    let pool = if common.cram_input {
        match ThreadPool::new(decode_t.max(1) as u32) {
            Ok(p) => Some(p),
            Err(e) => return Err(fail(common, &mut children, None, anyhow!("creating thread pool: {e}"))),
        }
    } else {
        None
    };
    if let Some(p) = &pool {
        if let Err(e) = reader.set_thread_pool(p) {
            return Err(fail(common, &mut children, None, anyhow!("attaching thread pool: {e}")));
        }
    }

    // Output BAM writer (uncompressed) into sort's stdin.
    let sort_in_fd = children.sort_stdin.as_ref().unwrap().as_raw_fd();
    let header_view = bam::HeaderView::from_bytes(common.header_bytes);
    let header = bam::Header::from_template(&header_view);
    let mut writer = match bam::Writer::from_path(fd_path(sort_in_fd), &header, Format::Bam) {
        Ok(w) => w,
        Err(e) => {
            let err = fail(common, &mut children, None, anyhow!("creating BAM writer: {e}"));
            return Err(err);
        }
    };
    if let Err(e) = writer.set_compression_level(bam::CompressionLevel::Uncompressed) {
        let err = fail(common, &mut children, None, anyhow!("setting uncompressed level: {e}"));
        return Err(err);
    }
    let writer = Arc::new(Mutex::new(writer));

    // Reader thread (owns minimap2's stdout).
    let shared_err: SharedError = Arc::new(Mutex::new(None));
    let new_names: Vec<String> = common.new_fai.iter().map(|(n, _)| n.clone()).collect();
    let thread = match children.mm_stdout.take() {
        Some(out) => spawn_reader(out, writer.clone(), new_names, common.tcfg.xl_tag, shared_err.clone()),
        None => {
            let err = fail(common, &mut children, None, anyhow!("missing minimap2 stdout"));
            return Err(err);
        }
    };

    // Sink.
    let state = Arc::new(Mutex::new(SinkState::default()));
    let (tee, tsv) = match common.deferred_prefix {
        Some(prefix) => {
            let fq = File::create(suffix(prefix, ".deferred.fq")).ok().map(BufWriter::new);
            let tsv = File::create(suffix(prefix, ".deferred.tsv")).ok().map(BufWriter::new);
            (fq, tsv)
        }
        None => (None, None),
    };
    let sink = MergeSink {
        writer: writer.clone(),
        fq: children.mm_stdin.take().map(BufWriter::new),
        tee,
        tsv,
        state: state.clone(),
        shared_err: shared_err.clone(),
    };

    // Run the pass.
    let t_pass = Instant::now();
    let pass = match &mode {
        RunMode::Normal { preset, diff } => {
            let ctx = Ctx {
                old: common.old,
                new: common.new,
                diff,
                old_names: common.old_names,
                new_tids: common.new_tids,
                cfg: common.tcfg,
                thr: common.max_de,
            };
            run_translate_pass(reader, sink, ctx, *preset, t_pass)
        }
        RunMode::AllDefer { reason, .. } => {
            let ctx = Ctx {
                old: common.old,
                new: common.new,
                diff: &Diff::default(),
                old_names: common.old_names,
                new_tids: common.new_tids,
                cfg: common.tcfg,
                thr: 0.0,
            };
            // A tiny-input fallback (`reason` is `None`) still validates the
            // record-level contract so a real violation takes precedence.
            let validate = reason.is_none();
            let defer_reason = if validate { Reason::TinyInput } else { Reason::Contract };
            run_all_defer(reader, sink, ctx, reason.clone(), defer_reason, validate, t_pass)
        }
    };
    let pass = match pass {
        Ok(p) => p,
        Err(e) => return Err(fail(common, &mut children, Some(thread), e)),
    };
    let pass_s = t_pass.elapsed().as_secs_f64();

    let mid = match pass.outcome {
        PassOutcome::Completed => None,
        PassOutcome::ContractViolation(r) => Some(r),
    };
    if let Some(reason) = mid {
        // Abort: reap, join, drop staged; the caller may restart in all-defer.
        let _ = thread.join();
        let _ = children.reap();
        common.remove_staged();
        return Ok(LiftOutcome::MidStream(reason));
    }

    let mut stats = pass.stats;
    stats.timings_s.diff = common.diff_s;
    stats.timings_s.translate_pass = pass_s;
    stats.timings_s.stream = pass_s;
    if let RunMode::AllDefer { preset, .. } = &mode {
        stats.preset = preset.minimap2_name().to_string();
    }

    // Drain + close.
    let t_drain = Instant::now();
    let remap = match thread.join() {
        Ok(Ok(r)) => r,
        Ok(Err(e)) => return Err(fail(common, &mut children, None, e)),
        Err(_) => {
            return Err(fail(common, &mut children, None, anyhow!("reader thread panicked")));
        }
    };
    stats.timings_s.writer_drain = t_drain.elapsed().as_secs_f64();
    drop(writer);
    drop(children.sort_stdin.take());

    // Reap and check exit codes.
    let t_sort = Instant::now();
    let (mm_status, sort_status) = children.wait();
    stats.timings_s.sort_index = t_sort.elapsed().as_secs_f64();
    stats.timings_s.minimap2_total = children.t0_mm.elapsed().as_secs_f64();
    stats.timings_s.minimap2_index =
        parse_index_ready_s(children.mm_stderr.as_deref()).unwrap_or(0.0);
    if let Some(s) = &sort_status {
        if !s.success() {
            common.remove_staged();
            return Err(anyhow!(
                "samtools sort exited with {}{}",
                status_desc(s),
                stderr_snippet(children.sort_stderr.as_deref())
            ));
        }
    }
    if let Some(s) = &mm_status {
        if !s.success() {
            common.remove_staged();
            return Err(anyhow!(
                "minimap2 exited with {}{}",
                status_desc(s),
                stderr_snippet(children.mm_stderr.as_deref())
            ));
        }
    }

    // Final check (without re-reading the output).
    stats.remapped_reads = remap.primary_count;
    stats.out_primary_records = stats.translated_reads + remap.primary_count;
    if stats.out_primary_records != stats.records_in.primary {
        common.remove_staged();
        bail!(
            "read-count invariant violated: translated {} + remapped {} != primary in {}",
            stats.translated_reads,
            remap.primary_count,
            stats.records_in.primary
        );
    }
    let deferred = state.lock().unwrap().deferred_names.clone();
    if deferred != remap.primary_names {
        common.remove_staged();
        let only_out = remap.primary_names.difference(&deferred).count();
        let only_deferred = deferred.difference(&remap.primary_names).count();
        bail!(
            "remapped primary names != deferred names ({only_out} only-out, {only_deferred} only-deferred)"
        );
    }

    // Success: promote the staged output safely. The previous good pair is
    // moved aside first, the index is renamed before the CRAM, and any failure
    // restores the backups so `OUT`/`OUT.crai` are never left inconsistent.
    if let Err(e) = promote(common.staged, common.staged_crai, common.out) {
        common.remove_staged();
        return Err(e);
    }

    stats.timings_s.total = common.t0.elapsed().as_secs_f64();
    stats.timings_s.children_cpu = children_cpu_s();
    stats.params = Params::new(common.tcfg, common.max_de);
    Ok(LiftOutcome::Done(stats))
}

/// Cleanup on failure: reap the children, join the reader thread, remove the
/// staged output, and append child statuses to the error.
fn fail(
    common: &Common<'_>,
    children: &mut Children,
    thread: Option<JoinHandle<Result<RemapResult>>>,
    primary: anyhow::Error,
) -> anyhow::Error {
    let (mm, sort) = children.reap();
    if let Some(t) = thread {
        let _ = t.join();
    }
    common.remove_staged();
    let mut parts = Vec::new();
    if let Some(s) = &mm {
        if !s.success() {
            parts.push(format!(
                "minimap2 exited with {}{}",
                status_desc(s),
                stderr_snippet(children.mm_stderr.as_deref())
            ));
        }
    }
    if let Some(s) = &sort {
        if !s.success() {
            parts.push(format!(
                "samtools sort exited with {}{}",
                status_desc(s),
                stderr_snippet(children.sort_stderr.as_deref())
            ));
        }
    }
    if parts.is_empty() {
        primary
    } else {
        anyhow!("{primary:#} ({})", parts.join("; "))
    }
}

/// Promote `staged`/`staged_crai` to `out`/`out.crai` without ever leaving an
/// inconsistent pair. Existing regular files are moved to `.bak`, the index is
/// renamed before the CRAM, and on any failure the backups are restored.
fn promote(staged: &Path, staged_crai: &Path, out: &Path) -> Result<()> {
    let out_crai = suffix(out, ".crai");
    let out_bak = suffix(out, ".bak");
    let crai_bak = suffix(&out_crai, ".bak");

    let mut moved_out = false;
    let mut moved_crai = false;
    if out.is_file() {
        std::fs::rename(out, &out_bak)
            .with_context(|| format!("backing up {}", out.display()))?;
        moved_out = true;
    }
    if out_crai.is_file() {
        if let Err(e) = std::fs::rename(&out_crai, &crai_bak) {
            if moved_out {
                let _ = std::fs::rename(&out_bak, out);
            }
            return Err(anyhow!("backing up {}: {e}", out_crai.display()));
        }
        moved_crai = true;
    }

    // Index first, then the CRAM.
    let res = std::fs::rename(staged_crai, &out_crai)
        .map_err(|e| anyhow!("renaming {} to {}: {e}", staged_crai.display(), out_crai.display()))
        .and_then(|_| {
            std::fs::rename(staged, out)
                .map_err(|e| anyhow!("renaming {} to {}: {e}", staged.display(), out.display()))
        });

    match res {
        Ok(()) => {
            if moved_out {
                let _ = std::fs::remove_file(&out_bak);
            }
            if moved_crai {
                let _ = std::fs::remove_file(&crai_bak);
            }
            Ok(())
        }
        Err(e) => {
            let _ = std::fs::remove_file(&out_crai);
            if moved_crai {
                let _ = std::fs::rename(&crai_bak, &out_crai);
            }
            if moved_out {
                let _ = std::fs::rename(&out_bak, out);
            }
            Err(e)
        }
    }
}

// ---------------------------------------------------------------------------
// entry point

/// Best-effort preset extraction from a (possibly invalid) minimap2 `@PG` CL.
pub(crate) fn best_effort_preset(header_text: &str) -> Option<Preset> {
    let map = |p: &str| match p {
        "map-hifi" => Some(Preset::MapHifi),
        "lr:hq" => Some(Preset::LrHq),
        _ => None,
    };
    for l in crate::contract::parse_header(header_text) {
        if l.rec_type != "PG" {
            continue;
        }
        let is_mm = l.field("PN") == Some("minimap2")
            || l.field("ID").is_some_and(|id| id.starts_with("minimap2"));
        if !is_mm {
            continue;
        }
        let Some(cl) = l.field("CL") else { continue };
        if let Ok(p) = parse_minimap2_cl(cl) {
            return Some(p);
        }
        let toks: Vec<&str> = cl.split_ascii_whitespace().collect();
        for i in 0..toks.len() {
            let t = toks[i];
            if t == "-ax" || t == "-x" {
                if let Some(p) = toks.get(i + 1).and_then(|p| map(p)) {
                    return Some(p);
                }
            } else if let Some(rest) = t.strip_prefix("-x") {
                if let Some(p) = map(rest) {
                    return Some(p);
                }
            }
        }
    }
    None
}

/// Run `cramhook lift`. Returns the stats.
#[allow(clippy::too_many_arguments)]
pub fn lift(
    in_cram: &Path,
    old_ref: &Path,
    new_ref: &Path,
    out: &Path,
    stats_path: Option<&Path>,
    threads: usize,
    tcfg: &TranslateConfig,
    opts: &LiftOptions,
    preset_override: Option<Preset>,
) -> Result<Stats> {
    let t0 = Instant::now();
    tcfg.validate()?;
    crate::translate::warn_if_unsafe_clip_screen(tcfg);
    ensure_parent(out)?;
    ensure_fai(old_ref)?;
    ensure_fai(new_ref)?;
    let old_fai = read_fai(&fai_path(old_ref))?;
    let new_fai = read_fai(&fai_path(new_ref))?;
    let new_tids: HashMap<String, i32> =
        new_fai.iter().enumerate().map(|(i, (n, _))| (n.clone(), i as i32)).collect();

    let mut reader = open_reader(in_cram, old_ref, threads)?;
    let header_text = String::from_utf8_lossy(reader.header().as_bytes()).to_string();
    let old_names: Vec<String> = reader
        .header()
        .target_names()
        .iter()
        .map(|n| String::from_utf8_lossy(n).to_string())
        .collect();
    check_sq_against_fai(&header_text, &old_fai)?;
    // A lift `@PG` under the old name is a hard error, not a remap fallback.
    reject_old_name_lift(&header_text).map_err(anyhow::Error::msg)?;

    let argv: String = crate::translate::recorded_argv();

    // Contract: preset from the @PG when parseable; otherwise fall back to the
    // `--preset` override, else error. A contract violation becomes an all-defer
    // fallback with the actual reason kept in the stats.
    let (effective_preset, mut contract_reason, chained) = match check_contract(&header_text) {
        Ok(c) => (c.preset, None, c.chained),
        Err(r) => {
            if opts.strict {
                bail!("input contract violation: {r}");
            }
            let preset = best_effort_preset(&header_text).or(preset_override).ok_or_else(|| {
                anyhow!(
                    "input contract violation ({r}) and no preset could be recovered from the @PG; pass --preset"
                )
            })?;
            (preset, Some(r), false)
        }
    };
    // A chained input is only accepted when the current lift's own
    // eligibility/safety options are at their defaults too; otherwise the input
    // is fully remapped (fail-safe), exactly like a contract violation.
    if contract_reason.is_none() && chained && !tcfg.eligibility_is_default() {
        let r = "chained input: the current lift must use default eligibility/safety options"
            .to_string();
        if opts.strict {
            bail!("input contract violation: {r}");
        }
        contract_reason = Some(r);
    }

    // M9 tiny-input fallback: lifting a very small CRAM is not worth its fixed
    // overhead (diff + sort/index), so remap every read instead. A contract
    // violation already remaps everything and takes precedence. The output goes
    // through the same all-defer path (each read exactly once, staged rename).
    let effective_min_input_mb =
        opts.min_input_mb.unwrap_or_else(|| default_min_input_mb(effective_preset));
    let input_bytes = std::fs::metadata(in_cram)
        .with_context(|| format!("stat input {}", in_cram.display()))?
        .len();
    let threshold_bytes = (effective_min_input_mb * 1_000_000.0).round() as u64;
    let tiny_input =
        contract_reason.is_none() && effective_min_input_mb > 0.0 && input_bytes < threshold_bytes;

    // Rule 7d scoring comes from the resolved read preset (runtime state).
    let mut tcfg_local = *tcfg;
    tcfg_local.scoring = scoring(effective_preset);
    // Trust an inherited `XL:A:T` record's stripped `s2` only for an all-default
    // chain (the contract already verified the earlier lifts; the check above
    // verified this lift).
    tcfg_local.inherited_xl_ok = chained && contract_reason.is_none();
    let tcfg = &tcfg_local;

    let remap_cmd = format!(
        "minimap2 -ax {} -t {} --secondary=no -Y {} -",
        effective_preset.minimap2_name(),
        threads,
        new_ref.display()
    );
    let ds = format!(
        "old-ref={} new-ref={} remap={}",
        old_ref.display(),
        new_ref.display(),
        remap_cmd
    );
    let header_bytes = build_output_header(&header_text, &new_fai, &argv, Some(&ds));

    let tmp_dir: Option<PathBuf> = if opts.keep_tmp {
        let mut d = out.as_os_str().to_owned();
        d.push(".tmp");
        let d = PathBuf::from(d);
        std::fs::create_dir_all(&d).context("creating keep-tmp dir")?;
        Some(d)
    } else {
        None
    };
    let deferred_prefix = tmp_dir.as_ref().map(|d| d.join("deferred"));

    let old = Fasta::load(old_ref)?;
    let new = Fasta::load(new_ref)?;

    let (staged, staged_crai) = staged_paths(out);
    let (decode_t, remap_t, sort_t, diff_t) = split_threads(threads);
    let exe_mm = executable("CRAMHOOK_MINIMAP2", "minimap2");
    let exe_sort = executable("CRAMHOOK_SAMTOOLS", "samtools");

    let common = Common {
        out,
        staged: &staged,
        staged_crai: &staged_crai,
        header_bytes: &header_bytes,
        new_fai: &new_fai,
        new_tids: &new_tids,
        old: &old,
        new: &new,
        old_names: &old_names,
        tcfg,
        diff_s: 0.0,
        t0,
        deferred_prefix: deferred_prefix.as_deref(),
        cram_input: in_cram
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("cram")),
        max_de: tcfg.max_de.unwrap_or_else(|| de_threshold(effective_preset)),
    };

    let empty_diff = Diff::default();

    // Spawn the children BEFORE the diff so minimap2's index build overlaps it.
    let mut children = spawn_children(
        effective_preset,
        new_ref,
        &exe_mm,
        &exe_sort,
        decode_t,
        remap_t,
        sort_t,
        &staged,
        &staged_crai,
    )?;

    let (diff, diff_s) = if contract_reason.is_none() && !tiny_input {
        let t = Instant::now();
        match compute_diff_with(
            &old,
            &new,
            old_ref,
            new_ref,
            effective_preset,
            diff_t,
            tcfg.diff_params(),
        ) {
            Ok(o) => {
                let s = t.elapsed().as_secs_f64();
                if let Some(d) = &tmp_dir {
                    if let Err(e) = write_outputs(&d.join("diff"), &o) {
                        let common2 = Common { diff_s: 0.0, ..common };
                        return Err(fail(&common2, &mut children, None, e));
                    }
                }
                (o.diff, s)
            }
            Err(e) => {
                // Diff failed: cancel the children.
                let common2 = Common { diff_s: 0.0, ..common };
                return Err(fail(&common2, &mut children, None, e));
            }
        }
    } else {
        (empty_diff, 0.0)
    };
    let common = Common { diff_s, ..common };

    // Normal or all-defer first attempt. All-defer covers both a contract
    // violation (with its reason) and the tiny-input fallback (reason `None`,
    // so `contract_ok` stays true).
    let mode = if let Some(reason) = &contract_reason {
        RunMode::AllDefer { preset: effective_preset, reason: Some(reason.clone()) }
    } else if tiny_input {
        RunMode::AllDefer { preset: effective_preset, reason: None }
    } else {
        RunMode::Normal { preset: effective_preset, diff: &diff }
    };
    let mut outcome = run_session(&common, children, &mut reader, mode, decode_t)?;

    // Mid-stream contract violation: restart in all-defer mode.
    if let LiftOutcome::MidStream(reason) = &outcome {
        if opts.strict {
            bail!("input contract violation (mid-stream): {reason}");
        }
        let mut reader2 = open_reader(in_cram, old_ref, threads)?;
        let children2 = spawn_children(
            effective_preset,
            new_ref,
            &exe_mm,
            &exe_sort,
            decode_t,
            remap_t,
            sort_t,
            &staged,
            &staged_crai,
        )?;
        outcome = run_session(
            &common,
            children2,
            &mut reader2,
            RunMode::AllDefer { preset: effective_preset, reason: Some(reason.clone()) },
            decode_t,
        )?;
    }

    let mut stats = match outcome {
        LiftOutcome::Done(s) => s,
        LiftOutcome::MidStream(r) => bail!("mid-stream contract violation persisted: {r}"),
    };
    stats.params.min_input_mb = effective_min_input_mb;
    stats.params.threads = threads;
    // Contract reason wins: only record (and announce) the fallback when the run
    // did not turn into a contract fallback.
    if tiny_input && stats.contract_ok {
        eprintln!(
            "tiny input: {} bytes < {} bytes (preset {}); remapping all reads",
            input_bytes,
            threshold_bytes,
            effective_preset.minimap2_name()
        );
        stats.fallback = Some(Fallback {
            reason: "tiny_input".to_string(),
            input_bytes,
            threshold_bytes,
        });
    }
    if let Some(p) = stats_path {
        write_stats(p, &stats)?;
    }
    Ok(stats)
}

fn write_stats(path: &Path, stats: &Stats) -> Result<()> {
    let f = File::create(path).with_context(|| format!("creating stats file {}", path.display()))?;
    let mut w = BufWriter::new(f);
    serde_json::to_writer_pretty(&mut w, &versioned(stats)).context("writing stats")?;
    w.flush().context("flushing stats")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Real minimap2 2.31 stderr for a 3-target build.
    const MM_STDERR: &str = "\
[M::mm_idx_gen::0.001*1.74] collected minimizers
[M::mm_idx_gen::0.001*2.07] sorted minimizers
[M::main::0.250*1.78] loaded/built the index for 3 target sequence(s)
[M::mm_mapopt_update::0.260*1.73] mid_occ = 50
[M::main] Real time: 0.300 sec; CPU: 0.600 sec; Peak RSS: 0.036 GB
";

    #[test]
    fn parse_index_ready_from_loaded_built_line() {
        let v = parse_index_ready_s(Some(MM_STDERR)).expect("must parse the loaded/built line");
        assert!((v - 0.250).abs() < 1e-9, "got {v}, want 0.250 (not Real time 0.300)");
    }

    #[test]
    fn parse_index_ready_absent() {
        assert_eq!(parse_index_ready_s(None), None);
        assert_eq!(parse_index_ready_s(Some("")), None);
        assert_eq!(parse_index_ready_s(Some("[M::main] no index line here\n")), None);
        // A `loaded/built` line without a parseable `::<secs>` prefix yields None.
        assert_eq!(parse_index_ready_s(Some("[M::main] loaded/built the index for 1\n")), None);
    }

    #[test]
    fn tail_chars_short_returns_input() {
        assert_eq!(tail_chars("abc", 8), "abc");
        assert_eq!(tail_chars("abcdefghij", 4), "ghij");
    }

    #[test]
    fn tail_chars_multibyte_safe() {
        // 5 two-byte chars: `len() - max` (= 5) lands mid-codepoint.
        let s = "αβγδε";
        assert_eq!(s.len(), 10);
        let got = tail_chars(s, 5);
        assert!(s.ends_with(got), "not a suffix: {got:?}");
        assert!(got.ends_with("δε"), "got {got:?}");
    }

    #[test]
    fn stderr_snippet_empty_and_text() {
        assert_eq!(stderr_snippet(None), "");
        assert_eq!(stderr_snippet(Some("  \n")), "");
        let got = stderr_snippet(Some("  boom  \n"));
        assert!(got.starts_with(": "), "missing prefix: {got:?}");
        assert!(got.contains("boom"), "missing text: {got:?}");
    }

    #[test]
    fn stderr_capture_file_is_private() {
        use std::os::unix::fs::PermissionsExt;
        let (f, path) = stderr_capture_file("unittest-mode").unwrap();
        drop(f);
        let mode = std::fs::metadata(&path).unwrap().permissions().mode() & 0o777;
        assert_eq!(mode, 0o600, "capture file mode is {mode:o}, want 600");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn read_stderr_is_lossy_and_removes_the_file() {        let (f, path) = stderr_capture_file("unittest-lossy").unwrap();
        drop(f);
        std::fs::write(&path, [b'o', b'k', 0xff, b'\n']).unwrap();
        let got = read_stderr(&Some(path.clone())).expect("captured text");
        assert!(got.contains("ok"), "lossy decode lost the ASCII text: {got:?}");
        assert!(!path.exists(), "capture file was not removed");
    }

    #[test]
    fn read_stderr_keeps_bounded_head_and_tail() {
        const HEAD: usize = 1 << 20;
        const TAIL: usize = 1 << 16;
        let (f, path) = stderr_capture_file("unittest-cap").unwrap();
        drop(f);
        // A head marker + the index-ready line (inside the head cap), then far
        // more than the cap of noise, then the real failure at the very end.
        let mut body = String::with_capacity(HEAD + (2 << 20) + 128);
        body.push_str("HEAD-MARKER\n");
        body.push_str("[M::main::0.250*1.0] loaded/built the index for 3 target sequence(s)\n");
        body.push_str(&"h".repeat(HEAD));
        body.push_str(&"n".repeat(2 << 20));
        body.push_str("FINAL-ERROR-MARKER\n");
        assert!(body.len() > HEAD + TAIL, "test body must exceed the caps: {}", body.len());
        std::fs::write(&path, &body).unwrap();

        let got = read_stderr(&Some(path.clone())).expect("captured text");
        assert!(got.starts_with("HEAD-MARKER"), "head was not preserved: {:?}", &got[..20]);
        assert!(got.contains("loaded/built the index"), "head lost the index line");
        assert!(got.contains("bytes omitted"), "no truncation marker");
        assert!(got.ends_with("FINAL-ERROR-MARKER\n"), "tail lost the error: {:?}", &got[got.len() - 32..]);
        assert!(got.len() <= HEAD as usize + TAIL as usize + 64, "not bounded: {}", got.len());
        // The failure snippet surfaces the real tail, not the head.
        let snip = stderr_snippet(Some(&got));
        assert!(snip.contains("FINAL-ERROR-MARKER"), "{snip:?}");
        assert!(!path.exists(), "capture file was not removed");
    }

    #[test]
    fn read_stderr_none_for_missing_and_empty() {
        assert_eq!(read_stderr(&None), None);
        let p = std::env::temp_dir().join("cramhook-unittest-does-not-exist.stderr");
        let _ = std::fs::remove_file(&p);
        assert!(!p.exists(), "fixture path unexpectedly exists: {}", p.display());
        assert_eq!(read_stderr(&Some(p)), None);
        let (f, path) = stderr_capture_file("unittest-empty").unwrap();
        drop(f);
        assert_eq!(read_stderr(&Some(path.clone())), None);
        assert!(!path.exists(), "empty capture file was not removed");
    }

    #[test]
    fn stderr_snippet_caps_long_text_and_keeps_prefix() {
        let long = format!("start{}end-marker", "x".repeat(5000));
        let s = stderr_snippet(Some(&long));
        assert!(s.starts_with(": "), "{s:?}");
        assert_eq!(s.len(), 2 + 2000, "cap not applied");
        assert!(s.ends_with("end-marker"), "wrong end kept: {s:?}");
    }
}
