#!/usr/bin/env python3
"""Evaluate `cramhook lift` on the real cases (docs/plans/Plan-v1.md §7, M5).

Usage:
  python scripts/real_eval.py [--root data/real] [--only SUBSTR[,SUBSTR]]
      [--threads 16] [--tag NAME] [--extra "--max-clip 500 --margin 300"]
      [--bench N] [--out data/real/eval]

For each data/real/<case>/ the script runs `cramhook lift`, compares the result
to oracle.cram with scripts/compare.py, collects agreement/perf metrics, and
writes <out>/<tag>/summary.{json,md}.  With ``--bench N`` it benchmarks the
oracle full-remap command and the lift command with hyperfine.
"""
import argparse
import gzip
import json
import os
import re
import shlex
import shutil
import socket
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import compare, write_result, read_fai_lengths  # noqa: E402

BIN = ROOT / "target" / "release" / "cramhook"
MANIFEST = "manifest.tsv"
GROUP_PLATFORM = {"pacbio": "PacBio", "ont": "ONT"}
T_MISMATCH_GATE = 0.001  # 0.1 %


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def run(cmd, cwd=None, timeout=None):
    if isinstance(cmd, str):
        argv = ["bash", "-o", "pipefail", "-c", cmd]
    else:
        argv = cmd
    return subprocess.run(argv, cwd=str(cwd) if cwd else None,
                          capture_output=True, text=True, timeout=timeout)


def read_json(path):
    path = Path(path)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def read_manifest(root):
    """case -> manifest row (platform, preset, ...)."""
    rows = {}
    path = root / MANIFEST
    if not path.exists():
        return rows
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            rows[parts[0]] = dict(zip(header, parts))
    return rows


def case_dir_ok(d):
    return all((d / f).exists() for f in ("old.fa", "new.fa", "old.cram", "oracle.cram",
                                          "oracle.json"))


def discover(root, only):
    manifest = read_manifest(root)
    cases = []
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        if not case_dir_ok(d):
            continue
        if only and not any(s in d.name for s in only):
            continue
        row = manifest.get(d.name, {})
        cases.append({
            "case": d.name,
            "dir": d,
            "platform": row.get("platform", ""),
            "preset": row.get("preset", ""),
        })
    return cases


def lift_flags_supported():
    """Set of option names `cramhook lift --help` advertises (best effort)."""
    p = run([str(BIN), "lift", "--help"])
    flags = set()
    for m in re.finditer(r"(--[a-z0-9-]+)", p.stdout or ""):
        flags.add(m.group(1))
    return flags


def contig_count(fa):
    try:
        return len(read_fai_lengths(fa))
    except Exception:
        return None


def top_reasons(stats, k=3):
    reasons = (stats or {}).get("deferred_by_reason")
    if not isinstance(reasons, dict):
        return []
    return sorted(reasons.items(), key=lambda kv: (-kv[1], kv[0]))[:k]


def tool_version(name):
    """First line of `<tool> --version`, or None if unavailable."""
    try:
        p = run([name, "--version"])
    except Exception:
        return None
    if p.returncode != 0:
        return None
    out = (p.stdout or "").strip().splitlines()
    return out[0] if out else None


def file_sig(path):
    """Resolved size+mtime_ns of a file (follows symlinks)."""
    try:
        p = Path(path).resolve()
        st = p.stat()
        return {"path": str(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    except OSError:
        return None


FASTA_EXT = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")
FASTQ_EXT = (".fq", ".fastq", ".fq.gz", ".fastq.gz")


def _is_fasta(p):
    return p.lower().endswith(FASTA_EXT)


def _is_fastq(p):
    return p.lower().endswith(FASTQ_EXT)


def cmd_ref_reads(cmd):
    """(refs, reads) paths named by the oracle pipeline.

    Refs = the minimap2 positional FASTA plus every ``samtools --reference``
    value; reads = the minimap2 positional FASTQ.
    """
    refs = [m.group(1) for m in re.finditer(r"--reference\s+(\S+)", cmd)]
    reads = []
    head = cmd.split("|", 1)[0]
    try:
        toks = shlex.split(head)
    except ValueError:
        toks = head.split()
    for t in toks[1:]:
        if t.startswith("-"):
            continue
        if _is_fastq(t):
            reads.append(t)
        elif _is_fasta(t):
            refs.append(t)
    return refs, reads


def full_bench_input_sigs(cmd):
    """Resolved signatures of every input the full cmd reads (+ ref .fai)."""
    refs, reads = cmd_ref_reads(cmd)
    out = []
    seen = set()
    for p in reads + refs:
        if p in seen:
            continue
        seen.add(p)
        out.append(file_sig(p))
        if _is_fasta(p):
            out.append(file_sig(p + ".fai"))
    return out


def _under_root(path):
    """Lexical (symlink-preserving) containment in the repo root."""
    ap = Path(os.path.abspath(str(path)))
    return ap == ROOT or ROOT in ap.parents


def _parent_inside_root(p):
    """The directory that would receive ``<p>.fai`` must live inside the repo.

    The named path itself may be a symlink to elsewhere (reads/refs are symlinked
    into case dirs); only the resolved parent directory matters for writes.
    """
    try:
        return _under_root(Path(p).parent.resolve())
    except OSError:
        return False


def cmd_inputs_outside_root(cmd):
    """Input paths unsafe to benchmark: a reference/reads whose resolved parent
    directory is outside the repo (so a samtools ``.fai`` write would land in a
    read-only source tree), or a ``<ref>.fai`` symlink resolving outside it."""
    bad = []
    refs, reads = cmd_ref_reads(cmd)
    for p in refs + reads:
        if not _parent_inside_root(p):
            bad.append(p)
            continue
        if _is_fasta(p):
            fai = Path(str(p) + ".fai")
            if fai.is_symlink():
                try:
                    if not _under_root(fai.resolve()):
                        bad.append(str(fai))
                except OSError:
                    bad.append(str(fai))
    return bad


def localise_cmd(cmd, case):
    """Point the cmd's ref/reads at the case-local symlinks when they match.

    The oracle cmd records the resolved source paths; re-running it could make
    samtools write a ``.fai`` into the read-only source tree.  Case dirs carry
    ``old.fa``/``new.fa``/``reads.fastq.gz`` symlinks, so map resolved targets
    back onto them.  Unknown paths are left untouched (and later refused by
    ``cmd_inputs_outside_root``).
    """
    d = case["dir"]
    repl = {}
    for name in ("old.fa", "new.fa", "reads.fastq.gz"):
        lp = d / name
        try:
            repl[str(lp.resolve())] = str(lp)
        except OSError:
            continue
    for src, dst in repl.items():
        cmd = cmd.replace(src, dst)
    return cmd


def fastq_lengths(path, want_names=False):
    """Stream a FASTQ (gzip if suffixed .gz).

    Returns (n_records, total_bases, {qname: seq_len}) where the qname map is
    only built when ``want_names`` is true (first record per name wins).  Raises
    ``ValueError`` on a malformed record (bad @ header, bad + line, or
    seq/qual length mismatch).
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    lengths = {}
    n = 0
    total = 0
    with opener(path, "rt") as fh:
        while True:
            head = fh.readline()
            if not head:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            n += 1
            if not head.startswith("@"):
                raise ValueError(f"record {n} in {path}: header does not start with '@'")
            if not plus.startswith("+"):
                raise ValueError(f"record {n} in {path}: '+' line missing")
            if not seq or not plus or not qual or not qual.endswith("\n"):
                raise ValueError(f"record {n} in {path}: truncated record")
            s, q = seq.rstrip("\n"), qual.rstrip("\n")
            if len(s) != len(q):
                raise ValueError(f"record {n} in {path}: seq/qual length mismatch "
                                 f"({len(s)} != {len(q)})")
            total += len(s)
            if want_names:
                name = head[1:].strip().split()[0]
                lengths.setdefault(name, len(s))
    return n, total, lengths


def rewrite_oracle_output(cmd, out_cram):
    """Rewrite the `-o <path>` of the oracle pipeline to ``out_cram``."""
    new, n = re.subn(r"(-o\s+)(\S+)", lambda m: m.group(1) + str(out_cram), cmd, count=1)
    if n != 1:
        raise ValueError("could not find -o <path> in oracle cmd")
    return new


# ----------------------------------------------------------------------------
# benchmarking
# ----------------------------------------------------------------------------
def hyperfine_mean(cmd, runs, json_path):
    p = subprocess.run(
        ["hyperfine", "--warmup", "1", "--runs", str(runs),
         "--export-json", str(json_path), "--shell", "bash", cmd],
        capture_output=True, text=True)
    if p.returncode != 0:
        return None, None, p.stderr.strip()[-500:]
    data = read_json(json_path) or {}
    results = data.get("results") or []
    if not results:
        return None, None, "hyperfine produced no results"
    return results[0].get("mean"), results[0].get("stddev"), None


def full_bench_key(case, full_cmd, threads, runs):
    """Cache key for the expensive oracle full-remap benchmark.

    The full remap is a pure function of this cmd + the environment + the input
    files, so any change to those invalidates the cached timing.  Only the files
    the cmd actually reads are signed (plus the reference's .fai).
    """
    return {
        "full_cmd": full_cmd,
        "minimap2": tool_version("minimap2"),
        "samtools": tool_version("samtools"),
        "inputs": full_bench_input_sigs(full_cmd),
        "threads": threads,
        "runs": runs,
        "hostname": socket.gethostname(),
    }


def bench_case(case, out_dir, threads, runs, extra, rebench_full=False):
    """Benchmark full remap (oracle cmd) vs lift. Returns (dict, error).

    The full-remap timing is cached per case in ``<case_dir>/bench_full.json``
    (keyed by the exact cmd, tool versions, input size+mtime, threads, runs and
    host).  A matching cache is reused; ``rebench_full`` forces a re-measure.
    """
    d = case["dir"]
    bench_dir = out_dir / "bench"
    bench_dir.mkdir(parents=True, exist_ok=True)
    full_out = bench_dir / f"{case['case']}.full.cram"
    lift_out = bench_dir / f"{case['case']}.lift.cram"

    oracle = read_json(d / "oracle.json") or {}
    oracle_cmd = oracle.get("cmd")
    if not oracle_cmd:
        return None, "oracle.json has no cmd"
    # Point the recorded source paths at the case-local symlinks, then rewrite
    # the output; refuse unsafe (external) inputs so a rerun cannot write into
    # the read-only source tree.
    oracle_cmd = localise_cmd(oracle_cmd, case)
    try:
        remap = rewrite_oracle_output(oracle_cmd, full_out)
    except ValueError as exc:
        return None, str(exc)
    full_cmd = f"set -o pipefail; {remap} && samtools index {full_out}"
    unsafe = cmd_inputs_outside_root(full_cmd)
    if unsafe:
        return None, f"refusing benchmark: cmd input outside repo root: {', '.join(unsafe)}"
    lift_cmd = (f"{BIN} lift --in {d / 'old.cram'} --old-ref {d / 'old.fa'} "
                f"--new-ref {d / 'new.fa'} --out {lift_out} --threads {threads}"
                + (f" {extra}" if extra else ""))

    cache_path = d / "bench_full.json"
    key = full_bench_key(case, full_cmd, threads, runs)
    cached = read_json(cache_path) or {}
    full_cached = (not rebench_full and cached.get("key") == key
                   and cached.get("mean") is not None)
    full_mean = full_sd = None
    try:
        if full_cached:
            full_mean = cached.get("mean")
            full_sd = cached.get("sd")
            print(f"  full remap: cached ({full_mean:.3f}s)")
        else:
            full_mean, full_sd, err = hyperfine_mean(
                full_cmd, runs, bench_dir / f"{case['case']}.full.json")
            if err:
                return None, f"full remap bench: {err}"
            try:
                cache_path.write_text(json.dumps(
                    {"key": key, "mean": full_mean, "sd": full_sd,
                     "measured": time.strftime("%Y-%m-%dT%H:%M:%S")}, indent=2) + "\n")
            except OSError:
                pass

        lift_mean, lift_sd, err = hyperfine_mean(
            lift_cmd, runs, bench_dir / f"{case['case']}.lift.json")
        if err:
            return None, f"lift bench: {err}"
    finally:
        for f in (full_out, lift_out, Path(str(full_out) + ".crai"),
                  Path(str(lift_out) + ".crai"),
                  bench_dir / f"{case['case']}.full.json",
                  bench_dir / f"{case['case']}.lift.json"):
            try:
                f.unlink()
            except OSError:
                pass

    speedup = (full_mean / lift_mean) if (full_mean and lift_mean) else None
    return {
        "full_mean_s": full_mean,
        "full_sd_s": full_sd,
        "lift_mean_s": lift_mean,
        "lift_sd_s": lift_sd,
        "speedup": speedup,
        "full_cached": full_cached,
    }, None


# ----------------------------------------------------------------------------
# per-case evaluation
# ----------------------------------------------------------------------------
def deferred_breakdown(case, cram, stats):
    """Bases/reads deferred, split by reason and class, from the keep-tmp files.

    The lift tmp dir is ``<out>.tmp`` (see src/lift.rs) holding the deferred
    FASTQ (one record per deferred read) and TSV (qname, reason, class).  The
    FASTQ bases are authoritative; every FASTQ qname must have exactly one valid
    TSV row and the read count must reconcile with ``stats.json``.  On any
    mismatch the breakdown is unavailable (``{"breakdown_error": ...}``) rather
    than silently partial.  The caller owns tmp-dir cleanup.
    """
    d = case["dir"]
    tmp = Path(str(cram) + ".tmp")
    if not tmp.is_dir():
        return {"breakdown_error": "lift tmp dir missing"}

    # The lift writes <tmp>/deferred.deferred.{fq,tsv} (prefix "deferred" plus
    # the ".deferred.*" suffix); accept the undoubled spelling too.
    def tmp_file(ext):
        for name in (f"deferred{ext}", f"deferred.deferred{ext}"):
            p = tmp / name
            if p.exists():
                return p
        return tmp / f"deferred{ext}"

    fq_path, tsv_path = tmp_file(".fq"), tmp_file(".tsv")
    if not fq_path.exists():
        return {"breakdown_error": f"deferred FASTQ missing ({fq_path.name})"}
    if not tsv_path.exists():
        return {"breakdown_error": f"deferred TSV missing ({tsv_path.name})"}

    try:
        n_fq, fq_bases, lengths = fastq_lengths(fq_path, want_names=True)
    except ValueError as exc:
        return {"breakdown_error": f"deferred FASTQ malformed: {exc}"}
    if len(lengths) != n_fq:
        return {"breakdown_error":
                f"duplicate qnames in deferred FASTQ ({n_fq} records, "
                f"{len(lengths)} unique)"}

    rows = {}
    problems = []
    for line in tsv_path.read_text().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 3:
            problems.append(f"malformed row {line!r}")
            continue
        qname, reason, cls = parts
        if qname in rows:
            problems.append(f"duplicate TSV row for {qname}")
            continue
        if qname not in lengths:
            problems.append(f"TSV qname absent from FASTQ: {qname}")
            continue
        rows[qname] = (reason, cls)
    if problems:
        return {"breakdown_error": "; ".join(problems[:5])}
    missing = set(lengths) - set(rows)
    if missing:
        return {"breakdown_error":
                f"{len(missing)} FASTQ qname(s) without a TSV row "
                f"(e.g. {sorted(missing)[:3]})"}

    expected = (stats or {}).get("deferred_reads")
    if expected is None:
        return {"breakdown_error": "stats.json has no deferred_reads"}
    if n_fq != expected:
        return {"breakdown_error":
                f"deferred reads {n_fq} != stats deferred_reads {expected}"}

    by_reason, by_reason_reads = {}, {}
    by_class, by_class_reads = {}, {}
    for qname, (reason, cls) in rows.items():
        length = lengths[qname]
        by_reason[reason] = by_reason.get(reason, 0) + length
        by_reason_reads[reason] = by_reason_reads.get(reason, 0) + 1
        by_class[cls] = by_class.get(cls, 0) + length
        by_class_reads[cls] = by_class_reads.get(cls, 0) + 1

    # Total input reads/bases from the original FASTQ (same read set as the
    # primary records in old.cram, cheaper than decoding CRAM).
    total_reads = total_bases = None
    reads_fq = d / "reads.fastq.gz"
    if reads_fq.exists():
        try:
            total_reads, total_bases, _ = fastq_lengths(reads_fq)
        except ValueError as exc:
            return {"breakdown_error": f"input FASTQ malformed: {exc}"}

    return {
        "deferred_bases_by_reason": by_reason,
        "deferred_reads_by_reason": by_reason_reads,
        "deferred_bases_by_class": by_class,
        "deferred_reads_by_class": by_class_reads,
        "deferred_reads": n_fq,
        "deferred_fq_reads": n_fq,
        "deferred_bases": fq_bases,       # FASTQ is authoritative
        "deferred_fq_bases": fq_bases,
        "total_reads": total_reads,
        "total_bases": total_bases,
    }


def evaluate_case(case, out_dir, tag, threads, extra, flags):
    d = case["dir"]
    target = out_dir / tag
    target.mkdir(parents=True, exist_ok=True)
    cram = target / f"{case['case']}.cram"
    stats_path = target / f"{case['case']}.stats.json"

    result = {"case": case["case"], "preset": case["preset"],
              "platform": case["platform"], "lift_rc": None, "error": None}

    extra_unsupported = []
    if extra:
        toks = shlex.split(extra)
        for i, t in enumerate(toks):
            if t.startswith("--") and t not in flags:
                extra_unsupported.append(t)
    result["extra_unsupported"] = extra_unsupported

    lift_cmd = [str(BIN), "lift", "--in", str(d / "old.cram"),
                "--old-ref", str(d / "old.fa"), "--new-ref", str(d / "new.fa"),
                "--out", str(cram), "--threads", str(threads),
                "--stats", str(stats_path), "--xl-tag", "--keep-tmp"]
    if extra:
        lift_cmd += shlex.split(extra)

    tmp = Path(str(cram) + ".tmp")
    try:
        t0 = time.monotonic()
        p = run(lift_cmd)
        lift_wall = time.monotonic() - t0
        result["lift_rc"] = p.returncode
        result["lift_wall_s"] = lift_wall
        if p.returncode != 0:
            result["error"] = (f"lift failed rc={p.returncode}: "
                               f"{p.stderr.strip()[-400:]}")
            return result, None

        stats = read_json(stats_path) or {}
        result["params"] = stats.get("params") if isinstance(stats, dict) else None

        # deferral by bases (from the keep-tmp deferred FASTQ/TSV)
        breakdown = deferred_breakdown(case, cram, stats)
        if "breakdown_error" in breakdown:
            result["breakdown_error"] = breakdown["breakdown_error"]
            result["deferred_bases_by_reason"] = None
        else:
            result.update(breakdown)
            tb = breakdown.get("total_bases")
            db = breakdown.get("deferred_bases")
            result["deferred_bases_fraction"] = (db / tb) if (db is not None and tb) \
                else None
            result["translated_bases_fraction"] = (1.0 - result["deferred_bases_fraction"]) \
                if result["deferred_bases_fraction"] is not None else None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # compare lift vs oracle
    res = compare(cram, d / "new.fa", d / "oracle.cram", d / "new.fa")
    write_result(res, target / f"{case['case']}.cmp")

    # contigs
    c_old = contig_count(d / "old.fa")
    c_new = contig_count(d / "new.fa")
    result["contigs_old"] = c_old
    result["contigs_new"] = c_new
    result["contigs_differ"] = c_old != c_new

    # stats-derived
    rec_in = (stats.get("records_in") or {}) if isinstance(stats, dict) else {}
    prim_in = rec_in.get("primary") if isinstance(rec_in, dict) else None
    translated = stats.get("translated_reads") if isinstance(stats, dict) else None
    result["translated_reads"] = translated
    result["translated_fraction"] = (translated / prim_in) if (translated is not None
                                                               and prim_in) else None
    result["deferred_by_class"] = stats.get("deferred_by_class")
    result["deferred_by_reason"] = stats.get("deferred_by_reason")
    result["top_reasons"] = top_reasons(stats)

    # agreement
    def cat(name):
        return res["primary"].get(name, {"n": 0, "n_identical": 0,
                                         "n_equivalent": 0, "n_mismatch": 0})
    tcat, rcat = cat("T"), cat("R")
    result["t"] = {"n": tcat["n"], "identical": tcat["n_identical"],
                   "equivalent": tcat["n_equivalent"], "mismatch": tcat["n_mismatch"]}
    result["t_mismatch_rate"] = (tcat["n_mismatch"] / tcat["n"]) if tcat["n"] else 0.0
    result["r"] = {"n": rcat["n"], "identical": rcat["n_identical"],
                   "equivalent": rcat["n_equivalent"], "mismatch": rcat["n_mismatch"]}
    so = res["supplementary"]["overall"]
    result["supplementary"] = dict(so)
    stcat = res["supplementary"].get("T", {"n": 0, "equal": 0,
                                           "equivalent": 0, "differ": 0})
    result["t_supp"] = dict(stcat)
    sn = stcat["n"]
    result["t_supp_fraction"] = ((stcat["equal"] + stcat["equivalent"]) / sn) if sn else 1.0
    dp = res["depth"]["overall"]
    result["depth"] = {"fraction_identical": dp["fraction_identical"],
                       "max_abs_diff": dp["max_abs_diff"],
                       "n_positions": dp["n_positions"]}
    result["invariants_ok"] = bool(res["invariants_ok"])
    result["mismatch_total"] = res["mismatch_total"]

    # timing vs recorded oracle
    oracle = read_json(d / "oracle.json") or {}
    result["oracle_wall_s"] = oracle.get("wall_s")
    if result["oracle_wall_s"]:
        result["speedup_oracle"] = result["oracle_wall_s"] / lift_wall
    return result, res


# ----------------------------------------------------------------------------
# summary
# ----------------------------------------------------------------------------
def _dist(values):
    if not values:
        return None
    return {"n": len(values), "min": min(values), "median": statistics.median(values),
            "max": max(values), "sorted": sorted(values)}


def speedup_groups(results):
    groups = {}
    for grp, pred in (
        ("PacBio", lambda r: GROUP_PLATFORM.get(r["platform"]) == "PacBio"),
        ("ONT", lambda r: GROUP_PLATFORM.get(r["platform"]) == "ONT"),
        ("cons-to-*", lambda r: "__cons-to-" in r["case"]),
        ("*-to-cons", lambda r: "__cons-to-" not in r["case"]),
    ):
        vals = [r["speedup"] for r in results if r.get("speedup") and pred(r)]
        groups[grp] = _dist(vals)
    return groups


def best_speedup(r):
    return r.get("speedup_bench") or r.get("speedup_oracle")


def deferral_by_bases_lines(results):
    """Markdown section: deferral aggregated by read bases (remap cost driver).

    An empty reason dict is a valid breakdown (zero deferrals): such cases still
    count toward the case count and the input-base denominator.  Only genuinely
    missing/failed breakdowns are excluded (and listed).
    """
    agg = {}  # reason -> [reads, bases]
    total_deferred = 0
    total_input = 0
    n_cases = 0
    excluded = []
    for r in results:
        if r.get("error"):
            excluded.append((r["case"], r["error"]))
            continue
        if r.get("breakdown_error"):
            excluded.append((r["case"], f"breakdown: {r['breakdown_error']}"))
            continue
        if "deferred_bases_by_reason" not in r or r["deferred_bases_by_reason"] is None:
            excluded.append((r["case"], "no deferral breakdown"))
            continue
        n_cases += 1
        rb = r["deferred_bases_by_reason"]
        rr = r.get("deferred_reads_by_reason") or {}
        for reason, bases in rb.items():
            a = agg.setdefault(reason, [0, 0])
            a[0] += rr.get(reason, 0)
            a[1] += bases
        total_deferred += r.get("deferred_bases") or 0
        total_input += r.get("total_bases") or 0

    lines = ["", "## deferral by bases", ""]
    if n_cases == 0:
        lines.append("- no deferral breakdown recorded (lift tmp files missing?)")
        for case, why in excluded:
            lines.append(f"- {case}: {why}")
        lines.append("")
        return lines

    pct_def = f"{total_deferred / total_input:.3%}" if total_input else "n/a"
    lines.append(f"- cases: {n_cases}; input bases: {total_input:,}; "
                 f"deferred bases: {total_deferred:,} ({pct_def})")
    for case, why in excluded:
        lines.append(f"- excluded {case}: {why}")

    families = {}
    for reason, (reads, bases) in agg.items():
        fam = "split_*" if reason.startswith("split") else "non-split"
        f = families.setdefault(fam, [0, 0])
        f[0] += reads
        f[1] += bases
    for fam in ("split_*", "non-split"):
        if fam not in families:
            continue
        reads, bases = families[fam]
        p = f"{bases / total_deferred:.1%}" if total_deferred else "n/a"
        lines.append(f"- family {fam}: {reads:,} reads, {bases:,} bases ({p} of deferred bases)")

    lines += ["", "| reason | reads | bases | % deferred bases | % input bases |",
              "|---|---|---|---|---|"]
    for reason, (reads, bases) in sorted(agg.items(), key=lambda kv: (-kv[1][1], kv[0])):
        pd = f"{bases / total_deferred:.2%}" if total_deferred else "n/a"
        pi = f"{bases / total_input:.3%}" if total_input else "n/a"
        lines.append(f"| {reason} | {reads:,} | {bases:,} | {pd} | {pi} |")

    lines += ["", "per-case top-3 reasons by bases:", ""]
    for r in results:
        if r.get("error"):
            continue
        if r.get("breakdown_error"):
            lines.append(f"- {r['case']}: n/a ({r['breakdown_error']})")
            continue
        if not isinstance(r.get("deferred_bases_by_reason"), dict):
            continue
        top = sorted(r["deferred_bases_by_reason"].items(),
                     key=lambda kv: (-kv[1], kv[0]))[:3]
        s = ", ".join(f"{k}={v:,}" for k, v in top) or "—"
        lines.append(f"- {r['case']}: {s}")
    lines.append("")
    return lines


def write_summary(results, out_dir, tag, args, partial=False):
    head = f"# cramhook real-data eval (tag `{tag}`)"
    if partial:
        head += " — PARTIAL (no summary.json; results from *.result.json only)"
    lines = [head, ""]
    lines.append(f"- cases: {len(results)}; threads: {args.threads}; "
                 f"bench runs: {args.bench if args.bench else 'off'}")
    lines.append(f"- lift flags supported: {', '.join(sorted(args.flags))}")
    if args.extra:
        lines.append(f"- extra args: `{args.extra}`")
    lines.append("")
    lines.append("| case | preset | contigs | translated | deferred bases % | "
                 "top deferrals | T ident/equiv/mism | T supp | depth | "
                 "wall lift/full | speedup |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in results:
        if r.get("error"):
            lines.append(f"| {r['case']} | {r['preset']} | — | — | — | — | — | — | — | "
                         f"ERROR: {r['error']} | — |")
            continue
        contigs = f"{r['contigs_old']}→{r['contigs_new']}"
        if r.get("contigs_differ"):
            contigs += " !"
        tf = f"{r['translated_fraction']:.3%}" if r.get("translated_fraction") is not None \
            else "n/a"
        dbf = f"{r['deferred_bases_fraction']:.3%}" \
            if r.get("deferred_bases_fraction") is not None else "n/a"
        reasons = ", ".join(f"{k}={v}" for k, v in (r.get("top_reasons") or [])) or "—"
        t = r["t"]
        tm = f"{t['identical']}/{t['equivalent']}/{t['mismatch']}"
        tsupfrac = r.get("t_supp_fraction", 1.0)
        tsup = f"{tsupfrac:.4%}"
        t_sup = r.get("t_supp", {}) or {}
        differ = t_sup.get("differ", 0)
        differ_raw = t_sup.get("differ_raw")
        if tsupfrac < 1.0:
            tsup = "!" + tsup + f" ({differ} differ)"
        if differ_raw is not None and differ_raw != differ:
            ignored = (t_sup.get("low_mapq_ignored_test", 0)
                       + t_sup.get("low_mapq_ignored_oracle", 0))
            tsup += f" ({ignored} lowQ ignored)"
        dep = r["depth"]
        depth = f"{dep['fraction_identical']:.4%} max{dep['max_abs_diff']}"
        full = r.get("full_wall_s")
        fulls = f"{full:.3f}" if full else "n/a"
        sp = best_speedup(r)
        sps = f"{sp:.2f}x" if sp else "n/a"
        lines.append(f"| {r['case']} | {r['preset']} | {contigs} | {tf} | {dbf} | "
                     f"{reasons} | {tm} | {tsup} | {depth} | "
                     f"{r['lift_wall_s']:.3f}/{fulls} | {sps} |")

    lines += deferral_by_bases_lines(results)

    dist = speedup_groups([dict(r, speedup=best_speedup(r)) for r in results])
    lines += ["", "## speedup distribution", ""]
    allv = [best_speedup(r) for r in results if best_speedup(r)]
    overall = _dist(allv)
    if overall:
        lines.append(f"- overall: n={overall['n']} min={overall['min']:.2f}x "
                     f"median={overall['median']:.2f}x max={overall['max']:.2f}x")
    else:
        lines.append("- overall: no speedups recorded (run with --bench N or ensure "
                     "oracle.json has wall_s)")
    for grp in ("PacBio", "ONT", "cons-to-*", "*-to-cons"):
        d = dist.get(grp)
        if not d:
            continue
        vals = " ".join(f"{v:.2f}" for v in d["sorted"])
        lines.append(f"- {grp}: n={d['n']} min={d['min']:.2f}x "
                     f"median={d['median']:.2f}x max={d['max']:.2f}x  [{vals}]")
    lines.append("")

    md = "\n".join(lines)
    (out_dir / tag / "summary.md").write_text(md + "\n")

    summary = {
        "tag": tag,
        "threads": args.threads,
        "extra": args.extra,
        "bench_runs": args.bench,
        "lift_flags": sorted(args.flags),
        "results": results,
        "speedup_distribution": dist,
    }
    (out_dir / tag / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return md


def load_tag_results(tag_dir):
    """Assemble the case results for a tag for ``--summary-only``.

    If ``summary.json`` exists its case set is authoritative: results come from
    a per-case ``*.result.json`` only when that file is newer than summary.json,
    otherwise from the summary entry.  Extra result files (cases not in the
    summary) are never silently merged in — they are reported.  With no
    summary.json, all result files are used and the output is marked PARTIAL.

    Returns ``(results, partial, notes)`` or ``(None, partial, notes)``.
    """
    meta_path = tag_dir / "summary.json"
    meta = read_json(meta_path)
    res_files = {}
    for p in tag_dir.glob("*.result.json"):
        res_files[p.name[: -len(".result.json")]] = p
    notes = []
    if isinstance(meta, dict) and meta.get("results"):
        meta_mtime = meta_path.stat().st_mtime if meta_path.exists() else 0.0
        by_case = {r.get("case"): r for r in meta["results"]}
        out = []
        for case in [r.get("case") for r in meta["results"]]:
            f = res_files.get(case)
            if f is not None and f.stat().st_mtime > meta_mtime:
                d = read_json(f)
                if isinstance(d, dict):
                    out.append(d)
                    notes.append(f"{case}: *.result.json")
                    continue
            out.append(by_case[case])
            notes.append(f"{case}: summary.json")
        for case in sorted(c for c in res_files if c not in by_case):
            notes.append(f"{case}: ignored (not in summary.json)")
        return out, False, notes
    if res_files:
        out = []
        for case in sorted(res_files):
            d = read_json(res_files[case])
            if isinstance(d, dict):
                out.append(d)
                notes.append(f"{case}: *.result.json")
        if out:
            return out, True, notes
    return None, True, ["no results found"]


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    global BIN
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default=str(ROOT / "data" / "real"))
    ap.add_argument("--only", default=None,
                    help="comma-separated case substrings")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--tag", default="default")
    ap.add_argument("--extra", default=None)
    ap.add_argument("--bench", type=int, default=0)
    ap.add_argument("--rebench-full", action="store_true",
                    help="ignore the cached full-remap benchmark and re-measure it")
    ap.add_argument("--summary-only", action="store_true",
                    help="re-render summary.md/json for --tag from the per-case "
                         "result JSONs already in --out; runs no lift/bench")
    ap.add_argument("--binary", default=str(BIN),
                    help="cramhook binary to use (default: target/release/cramhook)")
    ap.add_argument("--out", default=str(ROOT / "data" / "real" / "eval"))
    args = ap.parse_args()

    out_dir = Path(args.out)
    if args.summary_only:
        tag_dir = out_dir / args.tag
        results, partial, notes = load_tag_results(tag_dir)
        if results is None:
            print(f"error: no per-case results or summary.json under {tag_dir}",
                  file=sys.stderr)
            return 2
        meta = read_json(tag_dir / "summary.json") or {}
        args.threads = meta.get("threads", args.threads)
        args.extra = meta.get("extra", args.extra)
        args.bench = meta.get("bench_runs", args.bench)
        args.flags = set(meta.get("lift_flags") or [])
        print("summary-only sources:")
        for n in notes:
            print(f"  {n}")
        md = write_summary(results, out_dir, args.tag, args, partial=partial)
        print(md)
        print(f"wrote {tag_dir / 'summary.md'} and {tag_dir / 'summary.json'}")
        return 0

    BIN = Path(args.binary)
    if not BIN.exists():
        print(f"error: {BIN} not found; build with `pixi run cargo build --release`",
              file=sys.stderr)
        return 2
    only = [s for s in args.only.split(",") if s] if args.only else None
    cases = discover(Path(args.root), only)
    if not cases:
        print(f"no real cases found under {args.root}", file=sys.stderr)
        return 2
    args.flags = lift_flags_supported()
    (out_dir / args.tag).mkdir(parents=True, exist_ok=True)

    results = []
    failures = []
    for case in cases:
        print(f"=== {case['case']} ===")
        try:
            r, _ = evaluate_case(case, out_dir, args.tag, args.threads,
                                 args.extra, args.flags)
        except Exception as exc:
            r = {"case": case["case"], "preset": case["preset"],
                 "platform": case["platform"], "error": repr(exc),
                 "lift_rc": None, "invariants_ok": False}
        if r.get("error"):
            failures.append(f"{case['case']}: {r['error']}")
        else:
            if not r.get("invariants_ok"):
                failures.append(f"{case['case']}: invariants_ok false")
            if r.get("t_mismatch_rate", 0) > T_MISMATCH_GATE:
                failures.append(f"{case['case']}: T mismatch rate "
                                f"{r['t_mismatch_rate']:.4%} > 0.1%")
            if r.get("t_supp_fraction", 1.0) < 1.0:
                failures.append(
                    f"{case['case']}: T supplementary agreement "
                    f"{r['t_supp_fraction']:.4%} < 100% "
                    f"({r.get('t_supp', {}).get('differ', 0)} differ)")
            if args.bench:
                bench, err = bench_case(case, out_dir, args.threads, args.bench,
                                        args.extra, rebench_full=args.rebench_full)
                if err:
                    r["bench_error"] = err
                    print(f"  bench failed: {err}", file=sys.stderr)
                else:
                    r.update({f"bench_{k}": v for k, v in bench.items()})
                    r["full_wall_s"] = bench["full_mean_s"]
                    r["speedup_bench"] = bench["speedup"]
        if not r.get("full_wall_s") and r.get("oracle_wall_s"):
            r["full_wall_s"] = r["oracle_wall_s"]
        results.append(r)
        (out_dir / args.tag / f"{case['case']}.result.json").write_text(
            json.dumps(r, indent=2) + "\n")

    md = write_summary(results, out_dir, args.tag, args)
    print()
    print(md)
    print(f"wrote {out_dir / args.tag / 'summary.md'} and "
          f"{out_dir / args.tag / 'summary.json'}")
    if failures:
        print("\nFAILURES:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
