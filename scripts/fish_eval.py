#!/usr/bin/env python3
"""Evaluate ``cramhook fish`` on the spiked assemblies from ``scripts/qc_spike.py``.

``cramhook fish`` selects the reads that disagree with the assembly and writes
``P.fastq.gz``, ``P.fish.tsv`` and ``P.fish.json`` (spec:
``scratch/manager/fish-step1-spec.md``).  ``P.fish.tsv`` columns::

    read  length  reason  finding  contig  pos  strand  record  detail  sa

with reasons ``unmapped clip large_ins large_del linked locus`` and ``finding``
the ``finding_<n>`` label of ``P.findings.bed`` (or ``.``).  This script measures
whether fish selects the *right* reads, using the planted events of
``scripts/qc_spike.py`` as truth.

A read attributed to one of an event's matched findings is scored **wrong** only
when its aligned footprint touches *neither* the planted event interval *nor* the
matched finding's own span(s) (both rows of a wrapped finding).  fish attributes
a clip/indel event to every finding whose span contains the event position
± ``--finding-pad`` (1000 bp, ``Findings::within`` in ``src/fish.rs``); testing a
read against that window would be tautological (it is the attribution criterion
itself), so the footprint is tested against the finding's span directly.  The
historic "wrong near the event only" count is kept as ``wrong_near_event``.

For every spiked replicate (and the per-sample control) it runs ``cramhook qc``
into ``scratch/fish-eval/<tag>/<case>/<rep>/qc.*`` and ``cramhook fish`` into
``.../<rep>/fish.fish.*`` (default categories), then scores the outputs.  A
second fish configuration can be run and scored under a subdir label with
``--fish-args`` / ``--fish-label``.  Without ``--run`` only existing outputs are
scored.  Event -> finding matching reuses ``scripts/qc_spike_eval.py``.

Usage::

  pixi run python scripts/fish_eval.py --run --tag f1 [--threads 16]
      [--fish-args '--scope genome --locus' --fish-label alt] [--self-test]

Output: ``scratch/fish-eval/<tag>/per_event.tsv``, ``genome_wide.tsv``,
``far_not_selected.tsv`` and ``summary.md``.
"""
import argparse
import gzip
import json
import shlex
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
import qc_spike as qs           # noqa: E402
import qc_spike_eval as qse     # noqa: E402  (load_events, parse_findings,
#                                                 finding_call, large_hit,
#                                                 genome_bp, overlap_len, run_qc)

BIN = ROOT / "target" / "release" / "cramhook"
SPIKE_ROOT = qs.SPIKE_ROOT
EVAL_ROOT = ROOT / "scratch" / "fish-eval"

# Truth: a must-find read covers >= min(TRUTH_MUST_CAP, L) of the event's
# original interval; a neutral read 1..min(cap, L)-1; wrong = 0 overlap.
TRUTH_MUST_CAP = 100

# ``--finding-pad`` fish uses to attribute a clip/indel event to a finding
# (``src/fish.rs`` ``Findings::within``); 1 kb default.  Not used for the wrong
# test: mirroring it would select exactly the attributed reads again, so the
# finding span itself is used (see the module docstring).
FINDING_PAD = 1000

# fish reasons (no ``no_bridge``: fish has no such category).
REASON_ORDER = ["unmapped", "clip", "large_ins", "large_del", "linked", "locus"]
FINDING_REASONS = frozenset(["clip", "large_ins", "large_del", "linked", "locus"])
KIND_ORDER = ["del", "ins_random", "dup", "collapse", "inv", "swap"]
RECALL_KINDS = ("del", "collapse")

FISH_COLS = ["read", "length", "reason", "finding", "contig", "pos", "strand",
             "record", "detail", "sa"]


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def pct(k, n):
    return "n/a" if not n else f"{100.0 * k / n:.1f}%"


def ratio(k, n):
    return "n/a" if not n else f"{100.0 * k / n:.1f}%"


def median_int(xs):
    return int(statistics.median(xs)) if xs else 0


def read_json(path):
    p = Path(path)
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text())
    except Exception:
        return {}


def mapped_path(d):
    """The mapped CRAM (or BAM) inside a rep/control dir, or None."""
    d = Path(d)
    for name in ("mapped.cram", "mapped.bam"):
        if (d / name).exists():
            return d / name
    return None


def resolve_original_fasta(sample, sdir):
    """Path of the unmodified assembly FASTA for a sample, or None."""
    for p in [sdir / "control" / "build.json"] + sorted(sdir.glob("rep*/build.json")):
        if not p.exists():
            continue
        info = read_json(p)
        for key in ("ref", "reference", "fasta", "new_fa", "original_fasta",
                    "orig_fa", "new"):
            v = info.get(key)
            if v and Path(v).exists():
                return Path(v)
        if info.get("oracle"):
            cand = Path(info["oracle"]).parent / "new.fa"
            if cand.exists():
                return cand
    d = Path(sample.get("dir", ""))
    cand = (d if d.is_absolute() else ROOT / d) / "new.fa"
    return cand if cand.exists() else None


def genome_bp_of(prefix, run_dir):
    bp = qse.genome_bp(prefix)
    if bp:
        return bp
    p = Path(run_dir) / "summary.json"
    if p.exists():
        cl = read_json(p).get("contig_lengths")
        if cl:
            return sum(cl.values())
    return 0


def fish_prefix(run_dir, label):
    """fish ``--out-prefix``: ``<run_dir>/fish`` or ``<run_dir>/<label>/fish``."""
    run_dir = Path(run_dir)
    return (run_dir / label / "fish") if label else (run_dir / "fish")


def fish_tsv_path(prefix):
    """The existing fish.tsv for a prefix (``.fish.tsv`` or ``.tsv``), or None."""
    for suffix in (".fish.tsv", ".tsv"):
        p = Path(str(prefix) + suffix)
        if p.exists():
            return p
    return None


def fish_tsv_exists(prefix):
    return fish_tsv_path(prefix) is not None


def fish_json_path(prefix):
    """The existing fish JSON for a prefix (``.fish.json`` or ``.json``), or None."""
    for suffix in (".fish.json", ".json"):
        p = Path(str(prefix) + suffix)
        if p.exists():
            return p
    return None


def fish_fastq_path(prefix):
    """The existing fish FASTQ for a prefix (``.fastq.gz`` or ``.fastq``), or None."""
    for suffix in (".fastq.gz", ".fastq"):
        p = Path(str(prefix) + suffix)
        if p.exists():
            return p
    return None


def fastq_lengths(path):
    """``{read name: sequence length}`` from a (possibly gzipped) FASTQ."""
    out = {}
    if path is None:
        return out
    p = Path(path)
    if not p.exists():
        return out
    try:
        fh = gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")
        with fh:
            while True:
                head = fh.readline()
                if not head:
                    break
                if not head.startswith("@"):
                    continue
                seq = fh.readline()
                fh.readline()
                fh.readline()
                name = head[1:].strip().split()[0] if len(head) > 1 else ""
                if name:
                    out.setdefault(name, len(seq.rstrip("\n")))
    except OSError:
        return out
    return out


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------
def parse_fish(path):
    """Parse ``P.fish.tsv`` into a list of dicts.  Missing file -> []."""
    out = []
    if path is None:
        return out
    p = Path(path)
    if not p.exists():
        return out
    lines = p.read_text().splitlines()
    if not lines:
        return out
    header = lines[0].split("\t")
    if header[:2] != ["read", "length"]:
        header = FISH_COLS
        body = lines
    else:
        body = lines[1:]
    idx = {h: i for i, h in enumerate(header)}

    def col(f, name, default):
        i = idx.get(name, default)
        return f[i] if 0 <= i < len(f) else ""

    for line in body:
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 3:
            continue
        try:
            length = int(float(col(f, "length", 1) or 0))
        except ValueError:
            length = 0
        out.append({
            "read": col(f, "read", 0), "length": length,
            "reason": col(f, "reason", 2), "finding": col(f, "finding", 3),
            "contig": col(f, "contig", 4), "pos": col(f, "pos", 5),
            "strand": col(f, "strand", 6), "record": col(f, "record", 7),
            "detail": col(f, "detail", 8), "sa": col(f, "sa", 9),
        })
    return out


def input_totals(qc_prefix):
    """Input read and base totals from the qc JSON's ``reads`` block."""
    d = read_json(str(qc_prefix) + ".json")
    r = d.get("reads", {}) if isinstance(d, dict) else {}
    reads = r.get("total") or r.get("records", {}).get("primary") or 0
    bp = r.get("bases_total") or 0
    return int(reads), int(bp)


# ---------------------------------------------------------------------------
# read-based truth (pysam)
# ---------------------------------------------------------------------------
def _open(path, ref):
    if path is None:
        return None
    p = Path(path)
    try:
        if str(p).endswith(".cram"):
            return pysam.AlignmentFile(str(p), "rc", reference_filename=str(ref))
        return pysam.AlignmentFile(str(p), "rb")
    except (ValueError, OSError):
        return None


def _skip_record(r):
    """Skip secondary / duplicate / qcfail / unmapped (keep primary + supp)."""
    return r.is_unmapped or r.is_secondary or r.is_duplicate or r.is_qcfail


def _blocks(read):
    try:
        return read.get_blocks()
    except Exception:
        return []


def _overlap_bp(blocks, lo, hi):
    return sum(qse.overlap_len(a, b, lo, hi) for a, b in blocks)


def control_overlaps(af, contig, lo, hi):
    """``{read name: union bp in [lo, hi)}`` (only overlap > 0).

    Primary + supplementary; secondary/dup/qcfail/unmapped skipped.  A read's
    overlap is the *union* of the aligned M/=/X blocks of ALL its counted
    records (primary + supplementary), so two disjoint 60 bp alignments of one
    read inside a 1 kb event count as 120 bp, not 60, and overlapping
    alignments are not double counted.
    """
    out = {}
    if af is None or hi <= lo:
        return out
    try:
        it = af.fetch(contig, lo, hi)
    except (ValueError, KeyError):
        return out
    blocks_by_name = defaultdict(list)
    for r in it:
        if _skip_record(r):
            continue
        for a, b in _blocks(r):
            if a < hi and lo < b:
                blocks_by_name[r.query_name].append((max(a, lo), min(b, hi)))
    for n, ivs in blocks_by_name.items():
        out[n] = sum(b - a for a, b in qse.merge_intervals(ivs))
    return out


def locate_names(cram_path, ref, names, cache):
    """``{name: "unmapped" | "contig:pos1" | "absent"}`` for the requested names.

    One full pass (read names have no index); cached per (cram, name-set).
    """
    if not names:
        return {}
    key = (str(cram_path), tuple(sorted(names)))
    if key not in cache:
        locs = {}
        af = _open(cram_path, ref)
        if af is not None:
            try:
                for r in af.fetch(until_eof=True):
                    n = r.query_name
                    if n not in names or n in locs:
                        continue
                    if r.is_unmapped:
                        locs[n] = "unmapped"
                    elif not (r.is_secondary or r.is_duplicate or r.is_qcfail):
                        locs[n] = f"{r.reference_name}:{r.reference_start + 1}"
            finally:
                af.close()
        cache[key] = locs
    return {n: cache[key].get(n, "absent") for n in names}


# ---------------------------------------------------------------------------
# running the tools
# ---------------------------------------------------------------------------
def run_fish(cram, ref, prefix, findings, haplotypes, threads, force,
             extra_args=None):
    """Run ``cramhook fish``; rerun when ``fish.tsv`` is missing."""
    if cram is None or ref is None:
        return False
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if not force and fish_tsv_exists(prefix):
        return True
    for o in prefix.parent.glob(prefix.name + ".*"):
        try:
            o.unlink()
        except OSError:
            pass
    cmd = [str(BIN), "fish", "--in", str(cram), "--ref", str(ref),
           "--out-prefix", str(prefix), "--threads", str(threads)]
    if findings and Path(findings).exists():
        cmd += ["--findings", str(findings)]
    if haplotypes and Path(haplotypes).exists():
        cmd += ["--haplotypes", str(haplotypes)]
    if extra_args:
        cmd += shlex.split(extra_args)
    p = subprocess.run(["pixi", "run"] + cmd, cwd=str(ROOT),
                       capture_output=True, text=True)
    if p.returncode != 0:
        print(f"  fish failed: {' '.join(cmd)}\n    {p.stderr.strip()[-400:]}",
              file=sys.stderr)
        return False
    return fish_tsv_exists(prefix)


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------
def _fish_sets(frows):
    sel_any, reason_reads = set(), defaultdict(set)
    for r in frows:
        n = r["read"]
        sel_any.add(n)
        reason_reads[r["reason"]].add(n)
    return sel_any, reason_reads


def fish_bp_totals(prefix, frows):
    """``(total_reads, total_bp, by_reason, source)`` for one fish run.

    ``fish.tsv``'s ``length`` column is the *event* length for clip/indel rows,
    not the read length, so bp is taken from ``P.fish.json``'s ``reads`` block
    (``reads``/``bp`` and per reason); if that is missing, read lengths are
    parsed from the FASTQ and combined with the reasons in ``fish.tsv``
    (per-reason bp = sum of read lengths of reads having that reason).

    ``by_reason`` values are ``(reads, bp)``.
    """
    if not fish_tsv_exists(prefix):
        return 0, 0, {rr: (0, 0) for rr in REASON_ORDER}, "none"
    fj = read_json(fish_json_path(prefix)) if fish_json_path(prefix) else {}
    rb = fj.get("reads") if isinstance(fj, dict) else None
    if isinstance(rb, dict) and "bp" in rb:
        by = {}
        for rr, v in (rb.get("by_reason") or {}).items():
            try:
                by[rr] = (int(v.get("reads", 0)), int(v.get("bp", 0)))
            except (TypeError, ValueError):
                by[rr] = (0, 0)
        for rr in REASON_ORDER:
            by.setdefault(rr, (0, 0))
        try:
            return int(rb.get("reads", 0)), int(rb.get("bp", 0)), by, "fish.json"
        except (TypeError, ValueError):
            pass
    lens = fastq_lengths(fish_fastq_path(prefix))
    sel_any, reason_reads = _fish_sets(frows)
    by = {rr: (len(reason_reads.get(rr, set())),
               sum(lens.get(n, 0) for n in reason_reads.get(rr, set())))
          for rr in REASON_ORDER}
    return (len(sel_any), sum(lens.get(n, 0) for n in sel_any), by, "fastq")


def score_rep(case, rep, sample, rep_dir, qc_prefix, control_dir, original,
              run_dir, label, cache, warns):
    findings = qse.parse_findings(str(qc_prefix) + ".findings.bed")
    fpref = fish_prefix(run_dir, label)
    frows = parse_fish(fish_tsv_path(fpref))
    if not fish_tsv_exists(fpref):
        warns.add(f"missing {Path(str(fpref)).name}.fish.tsv")
    elif not Path(str(fpref) + ".fastq.gz").exists():
        warns.add(f"missing {Path(str(fpref)).name}.fastq.gz")
    events = qse.load_events(Path(rep_dir) / "events.tsv")
    genome_bp = genome_bp_of(qc_prefix, run_dir)
    input_reads, input_bp = input_totals(qc_prefix)
    label_name = label or "default"

    sel_any, reason_reads = _fish_sets(frows)
    tot_reads, tot_bp, reason_bp, bp_source = fish_bp_totals(fpref, frows)

    spiked_af = _open(mapped_path(rep_dir), Path(rep_dir) / "spiked.fa")
    control_af = _open(mapped_path(control_dir), original) \
        if control_dir is not None and original is not None else None
    if control_af is None:
        warns.add(f"no readable control mapping/reference for {case}")

    event_rows = []
    must_union = set()
    far_names = set()
    for ev in events:
        fmatch = [f for f in findings if qse.large_hit(qse.finding_call(f), ev)]
        labels = {f["label"] for f in fmatch}
        sel_finding = {r["read"] for r in frows
                       if r["finding"] in labels and r["reason"] in FINDING_REASONS}
        lo = ev.get("orig_start0")
        hi = ev.get("orig_end0")
        if lo is None or hi is None:
            lo = hi = 0
        length = hi - lo
        rec = ev["kind"] in RECALL_KINDS

        # The matched finding's own span(s), in the spiked assembly, where fish
        # attributed the reads.  Collect every row carrying a matched id: a wrapped
        # finding is two rows sharing one label, and only one of them need match the
        # event (`fmatch`), so the other row would otherwise be missed.
        finding_spans = [(f["start"], f["end"]) for f in findings
                         if f["label"] in labels and f["contig"] == ev["contig"]]

        far_not = set()
        if rec:
            ov = control_overlaps(control_af, ev["contig"], lo, hi)
            cap = min(TRUTH_MUST_CAP, length) if length > 0 else TRUTH_MUST_CAP
            must = {n for n, o in ov.items() if o >= cap}
            neutral = {n for n, o in ov.items() if 1 <= o < cap}
            must_union |= must
            far_not = must - sel_any
            far_names |= far_not
            # Reads whose spiked-mapping footprint touches a matched finding span
            # (union of M/=/X blocks, primary+supplementary, as for `must`).
            near_finding = set()
            if spiked_af is not None:
                for fs, fend in finding_spans:
                    near_finding |= set(control_overlaps(
                        spiked_af, ev["contig"], fs, fend))
            sel_must = len(must & sel_finding)
            sel_neutral = len(neutral & sel_finding)
            wrong_near_event = len(sel_finding - must - neutral)
            sel_wrong = len(sel_finding - must - neutral - near_finding)
            recall_near = len(must & sel_finding)
            recall_total = len(must & sel_any)
            reason_must = {rr: len(must & reason_reads.get(rr, set()))
                           for rr in REASON_ORDER}
            n_must = len(must)
        else:
            sel_must = sel_neutral = sel_wrong = wrong_near_event = 0
            recall_near = recall_total = n_must = 0
            reason_must = {rr: 0 for rr in REASON_ORDER}

        event_rows.append({
            "case": case, "rep": rep, "platform": sample.get("platform", ""),
            "preset": sample.get("preset", ""), "label": label_name,
            "event_id": ev["event_id"], "kind": ev["kind"], "size": ev["size"],
            "contig": ev["contig"], "start0": ev["start0"], "end0": ev["end0"],
            "breakpoints": ",".join(str(x) for x in ev["breakpoints"]),
            "orig_start0": lo, "orig_end0": hi, "event_len": length,
            "matched": ";".join(sorted(labels)),
            "n_findings": len({f["label"] for f in fmatch}),
            "n_sel_finding": len(sel_finding),
            "n_must": n_must, "sel_must": sel_must, "sel_neutral": sel_neutral,
            "sel_wrong": sel_wrong, "wrong_near_event": wrong_near_event,
            "recall_near": recall_near,
            "recall_total": recall_total, "reason_must": reason_must,
            "far_not": far_not,
        })

    if far_names:
        locs = locate_names(mapped_path(rep_dir), Path(rep_dir) / "spiked.fa",
                            far_names, cache)
    else:
        locs = {}
    for row in event_rows:
        far = row["far_not"]
        row["far_not_selected"] = len(far)
        row["far_unmapped"] = sum(1 for n in far if locs.get(n) == "unmapped")
        elsewhere = Counter()
        detail = []
        for n in sorted(far):
            loc = locs.get(n, "absent")
            if loc == "unmapped":
                continue
            elsewhere[loc.split(":", 1)[0]] += 1
            detail.append((n, loc))
        row["far_elsewhere"] = sum(elsewhere.values())
        row["far_elsewhere_contigs"] = ",".join(
            f"{c}:{n}" for c, n in sorted(elsewhere.items()))
        row["far_detail"] = detail

    if spiked_af is not None:
        spiked_af.close()
    if control_af is not None:
        control_af.close()

    gw = {
        "case": case, "rep": rep, "platform": sample.get("platform", ""),
        "preset": sample.get("preset", ""), "label": label_name,
        "selected_reads": tot_reads,
        "selected_bp": tot_bp,
        "input_reads": input_reads, "input_bp": input_bp,
        "bp_source": bp_source,
        "by_reason": {rr: {"reads": reason_bp.get(rr, (0, 0))[0],
                           "bp": reason_bp.get(rr, (0, 0))[1]}
                      for rr in REASON_ORDER},
        "selected_must": len(sel_any & must_union),
        "selected_not_must": len(sel_any - must_union),
    }
    return {"case": case, "rep": rep, "platform": sample.get("platform", ""),
            "preset": sample.get("preset", ""), "label": label_name,
            "genome_bp": genome_bp, "events": event_rows, "gw": gw,
            "n_events": len(events), "n_fish_rows": len(frows)}


def score_control(case, sample, qc_prefix, run_dir, label, warns):
    fpref = fish_prefix(run_dir, label)
    frows = parse_fish(fish_tsv_path(fpref))
    if not fish_tsv_exists(fpref):
        warns.add(f"missing {Path(str(fpref)).name}.fish.tsv")
    elif not Path(str(fpref) + ".fastq.gz").exists():
        warns.add(f"missing {Path(str(fpref)).name}.fastq.gz")
    input_reads, input_bp = input_totals(qc_prefix)
    sel_any, reason_reads = _fish_sets(frows)
    tot_reads, tot_bp, reason_bp, bp_source = fish_bp_totals(fpref, frows)
    label_name = label or "default"
    gw = {
        "case": case, "rep": "control", "platform": sample.get("platform", ""),
        "preset": sample.get("preset", ""), "label": label_name,
        "selected_reads": tot_reads,
        "selected_bp": tot_bp,
        "input_reads": input_reads, "input_bp": input_bp,
        "bp_source": bp_source,
        "by_reason": {rr: {"reads": reason_bp.get(rr, (0, 0))[0],
                           "bp": reason_bp.get(rr, (0, 0))[1]}
                      for rr in REASON_ORDER},
        "selected_must": 0, "selected_not_must": len(sel_any),
    }
    return {"case": case, "platform": sample.get("platform", ""),
            "preset": sample.get("preset", ""), "label": label_name,
            "genome_bp": genome_bp_of(qc_prefix, run_dir), "gw": gw,
            "n_fish_rows": len(frows)}


# ---------------------------------------------------------------------------
# discovery / dataset scoring
# ---------------------------------------------------------------------------
def discover(samples, spike_root, only=None):
    for sample in samples:
        if only and not any(s in sample["case"] for s in only):
            continue
        sdir = Path(spike_root) / sample["case"]
        if not sdir.is_dir():
            continue
        reps = sorted((p.name for p in sdir.glob("rep*") if p.is_dir()),
                      key=lambda n: int(n.replace("rep", "") or 0))
        yield sample, sdir, reps


def base_outputs(prefix):
    """The qc outputs that must exist before fish can be run/scored."""
    return all(Path(str(prefix) + e).exists()
               for e in (".vcf", ".bed", ".findings.bed", ".json"))


def score_dataset(spike_root, eval_root, tag, samples, replicates=None, only=None,
                  labels=None, allow_missing_fish=False):
    qc_root = Path(eval_root) / tag
    labels = labels or [None]
    reps_res, ctrl_res, missing = [], [], []
    warns = set()
    cache = {}
    for sample, sdir, reps in discover(samples, spike_root, only):
        case = sample["case"]
        original = resolve_original_fasta(sample, sdir)
        if original is None:
            warns.add(f"no original FASTA for {case}")
        control_dir = sdir / "control"
        cprefix = qc_root / case / "control" / "qc"
        c_run_dir = qc_root / case / "control"
        for label in labels:
            if not base_outputs(cprefix):
                missing.append(f"{case}/control")
            elif not fish_tsv_exists(fish_prefix(c_run_dir, label)) \
                    and not allow_missing_fish:
                missing.append(f"{case}/control (missing fish.tsv"
                               + (f" [{label}]" if label else "") + ")")
            else:
                ctrl_res.append(score_control(case, sample, cprefix, c_run_dir,
                                              label, warns))
        for rep in reps:
            rn = rep.replace("rep", "")
            if replicates and int(rn or 0) not in replicates:
                continue
            rep_dir = sdir / rep
            prefix = qc_root / case / rep / "qc"
            run_dir = qc_root / case / rep
            for label in labels:
                if not base_outputs(prefix):
                    missing.append(f"{case}/{rep}")
                    break
                if not fish_tsv_exists(fish_prefix(run_dir, label)) \
                        and not allow_missing_fish:
                    missing.append(f"{case}/{rep} (missing fish.tsv"
                                   + (f" [{label}]" if label else "") + ")")
                    continue
                reps_res.append(score_rep(case, rep, sample, rep_dir, prefix,
                                          control_dir, original, run_dir, label,
                                          cache, warns))
    return reps_res, ctrl_res, missing, sorted(warns)


# ---------------------------------------------------------------------------
# aggregation / report
# ---------------------------------------------------------------------------
def agg_events(rows):
    g = {"n": 0, "sel_finding": [], "sel_finding_total": 0, "must_total": 0,
         "sel_must": 0, "sel_neutral": 0, "sel_wrong": 0,
         "wrong_near_event": 0, "recall_near": 0,
         "recall_total": 0, "far_not": 0, "far_unmapped": 0, "far_elsewhere": 0,
         "reason_must": defaultdict(int)}
    for r in rows:
        g["n"] += 1
        g["sel_finding"].append(r["n_sel_finding"])
        g["sel_finding_total"] += r["n_sel_finding"]
        g["must_total"] += r["n_must"]
        g["sel_must"] += r["sel_must"]
        g["sel_neutral"] += r["sel_neutral"]
        g["sel_wrong"] += r["sel_wrong"]
        g["wrong_near_event"] += r["wrong_near_event"]
        g["recall_near"] += r["recall_near"]
        g["recall_total"] += r["recall_total"]
        g["far_not"] += r["far_not_selected"]
        g["far_unmapped"] += r["far_unmapped"]
        g["far_elsewhere"] += r["far_elsewhere"]
        for k, v in r["reason_must"].items():
            g["reason_must"][k] += v
    return g


def agg_gw(gws):
    """Aggregate genome-wide per-rep dicts."""
    g = {"runs": len(gws), "selected_reads": 0, "selected_bp": 0,
         "input_reads": 0, "input_bp": 0, "selected_must": 0,
         "selected_not_must": 0,
         "by_reason": {rr: {"reads": 0, "bp": 0} for rr in REASON_ORDER}}
    for w in gws:
        g["selected_reads"] += w["selected_reads"]
        g["selected_bp"] += w["selected_bp"]
        g["input_reads"] += w["input_reads"]
        g["input_bp"] += w["input_bp"]
        g["selected_must"] += w["selected_must"]
        g["selected_not_must"] += w["selected_not_must"]
        for rr in REASON_ORDER:
            g["by_reason"][rr]["reads"] += w["by_reason"][rr]["reads"]
            g["by_reason"][rr]["bp"] += w["by_reason"][rr]["bp"]
    return g


def _ev_cols():
    return (["case", "rep", "platform", "preset", "label", "event_id", "kind",
             "size", "contig", "start0", "end0", "breakpoints", "orig_start0",
             "orig_end0", "event_len", "matched_findings", "n_findings",
             "n_sel_finding", "n_must", "sel_must", "sel_neutral", "sel_wrong",
             "wrong_near_event", "recall_near", "recall_total"]
            + [f"reason_must_{rr}" for rr in REASON_ORDER]
            + ["far_not_selected", "far_unmapped", "far_elsewhere",
               "far_elsewhere_contigs"])


def write_outputs(out_dir, tag, reps_res, ctrl_res, missing, warns):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    labels = []
    for res in reps_res + ctrl_res:
        if res["label"] not in labels:
            labels.append(res["label"])
    if not labels:
        labels = ["default"]

    # --- per-event TSV -----------------------------------------------------
    cols = _ev_cols()
    lines = ["\t".join(cols)]
    far_lines = ["label\tcase\trep\tevent_id\tread\tspiked_location"]
    for res in reps_res:
        for r in res["events"]:
            row = [res["case"], res["rep"], res["platform"], res["preset"],
                   res["label"], r["event_id"], r["kind"], r["size"], r["contig"],
                   r["start0"], r["end0"], r["breakpoints"], r["orig_start0"],
                   r["orig_end0"], r["event_len"], r["matched"], r["n_findings"],
                   r["n_sel_finding"], r["n_must"], r["sel_must"],
                   r["sel_neutral"], r["sel_wrong"], r["wrong_near_event"],
                   r["recall_near"], r["recall_total"]]
            row += [r["reason_must"].get(rr, 0) for rr in REASON_ORDER]
            row += [r["far_not_selected"], r["far_unmapped"], r["far_elsewhere"],
                    r["far_elsewhere_contigs"]]
            lines.append("\t".join(str(x) for x in row))
            for name, loc in r["far_detail"]:
                far_lines.append("\t".join([res["label"], res["case"], res["rep"],
                                            r["event_id"], name, loc]))
    (out_dir / "per_event.tsv").write_text("\n".join(lines) + "\n")
    (out_dir / "far_not_selected.tsv").write_text("\n".join(far_lines) + "\n")

    # --- genome-wide TSV ---------------------------------------------------
    gw_lines = ["label\tcase\trep\tplatform\trun\tselected_reads\tselected_bp\t"
                "input_reads\tinput_bp\tpct_reads\tpct_bp\tselected_must\t"
                "selected_not_must"
                + "".join(f"\t{rr}_reads\t{rr}_bp" for rr in REASON_ORDER)]
    for res in reps_res:
        gw_lines.append(_gw_tsv_row(res["gw"], "spiked"))
    for res in ctrl_res:
        gw_lines.append(_gw_tsv_row(res["gw"], "control"))
    (out_dir / "genome_wide.tsv").write_text("\n".join(gw_lines) + "\n")

    # --- summary.md --------------------------------------------------------
    md = [f"# cramhook fish eval (tag `{tag}`)", ""]
    md.append("Scores `P.fish.tsv` against the planted events of "
              "`scripts/qc_spike.py`.  Findings are matched to an event with the "
              "same rule as `scripts/qc_spike_eval.py` (`large_hit`, tolerance "
              "`max(50 bp, 1% of size)`).  Event->finding labels come from "
              "`P.findings.bed`; the input read/base totals come from "
              "`qc.json`'s `reads` block (`reads.total` / `bases_total`).")
    md.append("")
    md += ["## TRUTH DEFINITION (exactly as scored; per del/collapse event)", ""]
    md.append(f"- take a control read (primary+supplementary; secondary/dup/qcfail/"
              f"unmapped skipped) and its aligned `M/=/X` bp **inside** "
              f"`[orig_start0, orig_end0)`; `L = orig_end0 - orig_start0`, "
              f"`cap = min({TRUTH_MUST_CAP}, L)`.")
    md.append(f"- **must-find**: overlap `>= cap`; **neutral**: `1 <= overlap < cap` "
              f"(both measured against the planted interval).")
    md.append(f"- **overlap** is the *union* of the read's aligned `M/=/X` blocks "
              f"across all its counted records (primary+supplementary) inside the "
              f"interval (overlapping alignments are not double counted).")
    md.append("- **selected via a matched finding** = reads with a fish row whose "
              "`finding` is a matched finding of the event and whose reason is one "
              "of `clip`/`large_*`/`linked`/`locus` (an `unmapped` row has no finding "
              "id, so it is never counted as near).")
    md.append("- **wrong** (headline) = selected via a matched finding but its "
              "aligned footprint touches *neither* the planted event interval *nor* "
              "the matched finding's own span(s) (both rows of a wrapped finding). "
              "fish attributes a clip/indel event to every finding whose span "
              f"contains the event position ± `--finding-pad` "
              f"({FINDING_PAD} bp default) in `src/fish.rs`; scoring against that "
              "window would just re-select the attributed reads, so the footprint is "
              "tested against the finding's span directly (spiked coordinates).")
    md.append("- **wrong_near_event** (secondary) = the historic definition: "
              "selected via a matched finding with `0` bp overlap of the planted "
              "event interval (some of these sit on a linked far end of the same "
              "finding and belong to the problem region).")
    md.append("- **recall_near** = must-find reads selected via a matched finding "
              "(`clip`/`large_*`/`linked`/`locus` rows) / must-find reads.")
    md.append("- **recall_total** = must-find reads selected by ANY row (incl. "
              "`unmapped` and rows outside the event's findings) / must-find reads.")
    md.append("- **per-reason contribution** = must-find reads with at least one row "
              "of that reason (a read can be counted by several reasons).")
    md.append("- **far must-find reads not selected at all** = must-find reads with "
              "no fish row; we report their count and where they map in the spiked "
              "mapping (`unmapped` / elsewhere `contig:pos`).")
    md.append("- **precision** = among distinct reads selected via matched findings: "
              "must-find / neutral / wrong counts and the wrong share "
              "(`wrong` / all attributed reads).")
    md.append("")

    for label in labels:
        _emit_label(md, label, reps_res, ctrl_res)

    if missing:
        md += ["## missing outputs (not scored)", ""]
        for m in missing:
            md.append(f"- {m}")
        md.append("")
    if warns:
        md += ["## warnings", ""]
        for w in warns:
            md.append(f"- {w}")
        md.append("")

    (out_dir / "summary.md").write_text("\n".join(md) + "\n")

    # --- printed summary ---------------------------------------------------
    print(f"\n== cramhook fish eval (tag {tag}) ==")
    print(f"spiked reps scored: {len(reps_res)}; controls: {len(ctrl_res)}")
    if warns:
        print(f"warnings: {'; '.join(warns)}")
    for label in labels:
        tag_l = "" if label == "default" else f" [{label}]"
        for kind in KIND_ORDER:
            rows = [r for r in reps_res if r["label"] == label
                    for r in r["events"] if r["kind"] == kind]
            if not rows:
                continue
            g = agg_events(rows)
            rn = pct(g["recall_near"], g["must_total"])
            rt = pct(g["recall_total"], g["must_total"])
            print(f"  {label:8s}{tag_l} {kind:11s} events={g['n']:2d} "
                  f"sel_finding={g['sel_finding_total']:5d} must={g['must_total']:5d} "
                  f"recall_near={rn:>6s} recall_total={rt:>6s} "
                  f"far_not={g['far_not']:4d}")
    for label in labels:
        gw_sp = agg_gw([r["gw"] for r in reps_res if r["label"] == label])
        gw_ct = agg_gw([r["gw"] for r in ctrl_res if r["label"] == label])
        print(f"  [{label}] genome-wide spiked: {gw_sp['selected_reads']} reads "
              f"({ratio(gw_sp['selected_reads'], gw_sp['input_reads'])} of input), "
              f"must={gw_sp['selected_must']}; control: "
              f"{gw_ct['selected_reads']} reads "
              f"({ratio(gw_ct['selected_reads'], gw_ct['input_reads'])})")
    if missing:
        print(f"MISSING outputs for {len(missing)} run(s): "
              f"{', '.join(missing[:6])}"
              + (" ..." if len(missing) > 6 else ""), file=sys.stderr)
    print(f"wrote {out_dir}/per_event.tsv, genome_wide.tsv, "
          f"far_not_selected.tsv and summary.md")
    return out_dir


def _gw_tsv_row(w, run):
    row = [w["label"], w["case"], w["rep"], w["platform"], run,
           w["selected_reads"], w["selected_bp"], w["input_reads"], w["input_bp"],
           f"{100.0 * w['selected_reads'] / w['input_reads']:.3f}"
           if w["input_reads"] else "n/a",
           f"{100.0 * w['selected_bp'] / w['input_bp']:.3f}"
           if w["input_bp"] else "n/a",
           w["selected_must"], w["selected_not_must"]]
    for rr in REASON_ORDER:
        row += [w["by_reason"][rr]["reads"], w["by_reason"][rr]["bp"]]
    return "\t".join(str(x) for x in row)


def _ev_row_cells(kind, g, recall):
    # Wrong share is over every distinct read attributed to a matched finding.
    attributed = g["sel_finding_total"]
    cells = [
        kind, g["n"], median_int(g["sel_finding"]), g["sel_finding_total"],
        g["must_total"],
        g["recall_near"] if recall else "-",
        ratio(g["recall_near"], g["must_total"]) if recall else "-",
        g["recall_total"] if recall else "-",
        ratio(g["recall_total"], g["must_total"]) if recall else "-",
        g["sel_must"], g["sel_neutral"], g["sel_wrong"],
        ratio(g["sel_wrong"], attributed) if recall and attributed else "-",
        g["wrong_near_event"],
        g["far_not"], g["far_unmapped"], g["far_elsewhere"],
    ]
    return " | ".join(str(c) for c in cells)


EV_TABLE_HEADER = ("| kind | events | sel_finding med | sel_finding | must | "
                   "recall_near | recall_near % | recall_total | recall_total % | "
                   "sel_must | sel_neutral | sel_wrong | wrong_share | "
                   "wrong_near_event | far_not | "
                   "far_unmapped | far_elsewhere |")
EV_TABLE_SEP = "|---" * 17 + "|"


def _emit_label(md, label, reps_res, ctrl_res):
    reps = [r for r in reps_res if r["label"] == label]
    ctrls = [r for r in ctrl_res if r["label"] == label]
    all_events = [e for r in reps for e in r["events"]]
    head = "default" if label == "default" else label
    md += [f"## fish config `{head}`", ""]

    # per kind
    md += ["### per event kind", ""]
    md.append("`sel_finding` = distinct reads selected via a matched finding "
              "(median per event / total); `sel_any` = distinct selected reads; "
              "recall/precision columns are for `del`/`collapse` only.")
    md.append("")
    md.append(EV_TABLE_HEADER)
    md.append(EV_TABLE_SEP)
    for kind in KIND_ORDER:
        rows = [e for e in all_events if e["kind"] == kind]
        if not rows:
            continue
        md.append("| " + _ev_row_cells(kind, agg_events(rows),
                                       kind in RECALL_KINDS) + " |")
    md.append("")

    # per-reason recall contribution
    md += ["### per-reason recall contribution (must-find reads caught per reason)",
           ""]
    md.append("A must-find read counts once per reason that selected it; a read can "
              "appear under several reasons.  Only `del`/`collapse` events.")
    md.append("")
    md.append("| kind | " + " | ".join(REASON_ORDER) + " | must |")
    md.append("|---" * (len(REASON_ORDER) + 2) + "|")
    for kind in RECALL_KINDS:
        rows = [e for e in all_events if e["kind"] == kind]
        if not rows:
            continue
        g = agg_events(rows)
        md.append(f"| {kind} | "
                  + " | ".join(str(g["reason_must"].get(rr, 0))
                               for rr in REASON_ORDER)
                  + f" | {g['must_total']} |")
    md.append("")

    # per preset
    md += ["### per preset (platform)", ""]
    md.append("Same columns; `del`/`collapse` events only for recall/precision.")
    md.append("")
    md.append(EV_TABLE_HEADER)
    md.append(EV_TABLE_SEP)
    platforms = sorted({e["platform"] for e in all_events})
    for plat in platforms:
        rows = [e for e in all_events if e["platform"] == plat]
        md.append("| " + _ev_row_cells(plat, agg_events(rows), True) + " |")
    md.append("")

    # genome-wide
    md += ["### genome-wide extraction (fish selected reads / bp)", ""]
    md.append("Per-reason `reads/bp`; `% reads`/`% bp` are against the qc input "
              "totals.  Controls (unmodified assembly) are the noise baseline; for "
              "spiked runs `must-any` is the selected reads that are must-find for "
              "at least one `del`/`collapse` event.")
    md.append("")
    hdr = ("| run | platform | runs | sel reads | sel bp | % reads | % bp | "
           + " | ".join(f"{rr} R/B" for rr in REASON_ORDER) + " | must-any | "
           "not-must | must share |")
    md.append(hdr)
    md.append("|---" * (7 + len(REASON_ORDER) + 3) + "|")
    for plat in sorted({w["platform"] for w in
                        [r["gw"] for r in reps + ctrls]}):
        g = agg_gw([w["gw"] for w in reps if w["gw"]["platform"] == plat])
        if g["runs"]:
            md.append(_gw_md_row("spiked", plat, g))
        g = agg_gw([w["gw"] for w in ctrls if w["gw"]["platform"] == plat])
        if g["runs"]:
            md.append(_gw_md_row("control", plat, g))
    md.append("")


def _gw_md_row(run, plat, g):
    per = " | ".join(f"{g['by_reason'][rr]['reads']}/{g['by_reason'][rr]['bp']}"
                     for rr in REASON_ORDER)
    share = ratio(g["selected_must"],
                  g["selected_must"] + g["selected_not_must"])
    return (f"| {run} | {plat} | {g['runs']} | {g['selected_reads']} | "
            f"{g['selected_bp']} | {ratio(g['selected_reads'], g['input_reads'])} | "
            f"{ratio(g['selected_bp'], g['input_bp'])} | {per} | "
            f"{g['selected_must']} | {g['selected_not_must']} | {share} |")


# ---------------------------------------------------------------------------
# running the tools
# ---------------------------------------------------------------------------
def _second_label(fish_args, fish_label):
    """Label of the optional second fish config, or None when it is not run."""
    if fish_args is None and not fish_label:
        return None
    return fish_label or "alt"


def run(only, replicates, threads, force, tag, extra_args=None, fish_args=None,
        fish_label=None, allow_missing_fish=False):
    out_dir = EVAL_ROOT / tag
    out_dir.mkdir(parents=True, exist_ok=True)
    label2 = _second_label(fish_args, fish_label)
    labels = [None] + ([label2] if label2 else [])
    configs = [(None, None)] + ([(label2, fish_args)] if label2 else [])

    for sample, sdir, reps in discover(qs.SAMPLES, SPIKE_ROOT, only):
        case = sample["case"]
        preset = sample["preset"]
        original = resolve_original_fasta(sample, sdir)
        control_dir = sdir / "control"
        c_run_dir = out_dir / case / "control"
        cprefix = c_run_dir / "qc"
        c_run_dir.mkdir(parents=True, exist_ok=True)
        ccram = mapped_path(control_dir)
        if ccram and original is not None:
            qse.run_qc(ccram, original, preset, cprefix, threads, force,
                       extra_args)
            for label, fargs in configs:
                fpref = fish_prefix(c_run_dir, label)
                run_fish(ccram, original, fpref, str(cprefix) + ".findings.bed",
                         str(cprefix) + ".haplotypes.tsv", threads,
                         force or not fish_tsv_exists(fpref), fargs)
        else:
            print(f"{case}/control: no mapping or reference", file=sys.stderr)
        for rep in reps:
            rn = rep.replace("rep", "")
            if replicates and int(rn or 0) not in replicates:
                continue
            rep_dir = sdir / rep
            run_dir = out_dir / case / rep
            run_dir.mkdir(parents=True, exist_ok=True)
            prefix = run_dir / "qc"
            cram = mapped_path(rep_dir)
            if cram is None:
                print(f"{case}/{rep}: no mapping", file=sys.stderr)
                continue
            qse.run_qc(cram, rep_dir / "spiked.fa", preset, prefix, threads,
                       force, extra_args)
            for label, fargs in configs:
                suf = "" if label is None else f" [{label}]"
                fpref = fish_prefix(run_dir, label)
                ok = run_fish(cram, rep_dir / "spiked.fa", fpref,
                              str(prefix) + ".findings.bed",
                              str(prefix) + ".haplotypes.tsv", threads,
                              force or not fish_tsv_exists(fpref), fargs)
                if ok:
                    print(f"{case}/{rep}{suf}: qc+fish done")
                else:
                    print(f"{case}/{rep}{suf}: fish produced no fish.tsv",
                          file=sys.stderr)

    reps_res, ctrl_res, missing, warns = score_dataset(
        SPIKE_ROOT, EVAL_ROOT, tag, qs.SAMPLES, replicates, only, labels,
        allow_missing_fish)
    write_outputs(out_dir, tag, reps_res, ctrl_res, missing, warns)
    return 0


# ---------------------------------------------------------------------------
# self-test (synthetic data, no qc/fish binary)
# ---------------------------------------------------------------------------
def _write_bam(path, refs, recs):
    header = pysam.AlignmentHeader.from_references(
        [c for c, _ in refs], [n for _, n in refs])
    ordered = sorted(recs, key=lambda d: (0 if d.get("refid", 0) >= 0 else 1,
                                          d.get("refid", 0), d.get("pos0", 0)))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for d in ordered:
            a = pysam.AlignedSegment(header)
            a.query_name = d["name"]
            a.query_sequence = d["seq"]
            a.flag = d.get("flag", 0)
            a.reference_id = d.get("refid", 0)
            a.reference_start = d.get("pos0", 0)
            a.mapping_quality = 60
            a.cigarstring = d.get("cigar")
            out.write(a)
    pysam.index(str(path))


def _write_fa(path, name, seq):
    Path(path).write_text(f">{name}\n{seq}\n")


FISH_HEADER = "read\tlength\treason\tfinding\tcontig\tpos\tstrand\trecord\tdetail\tsa\n"


def self_test(_args):
    import shutil
    import tempfile
    problems = []

    # --- 1. fish.tsv parsing (header + bare rows) --------------------------
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        ft = tdp / "fish.fish.tsv"
        ft.write_text(
            FISH_HEADER
            + "r1\t500\tlarge_del\tfinding_1\tchr1\t5000\t+\tprimary\tlen=1000;at=5000\t.\n"
            + "r2\t400\tunmapped\t.\t.\t.\t.\tprimary\t.\t.\n")
        rows = parse_fish(ft)
        if len(rows) != 2 or rows[0]["length"] != 500 or rows[1]["reason"] != "unmapped" \
                or rows[1]["finding"] != ".":
            problems.append(f"parse_fish header: {rows}")
        ft2 = tdp / "bare.tsv"
        ft2.write_text("r3\t300\tclip\tfinding_1\tchr1\t10\t-\tsupplementary\t"
                       "side=left;len=60\tchr1:99:+\n")
        rows2 = parse_fish(ft2)
        if len(rows2) != 1 or rows2[0]["record"] != "supplementary" or rows2[0]["length"] != 300:
            problems.append(f"parse_fish bare: {rows2}")
        if parse_fish(tdp / "nope.tsv") != [] or parse_fish(None) != []:
            problems.append("parse_fish missing file not []")

    # --- 2. findings.bed parsing on a real file if one exists --------------
    real = sorted((ROOT / "scratch" / "qc").rglob("*.findings.bed"))
    if real:
        fs = qse.parse_findings(real[0])
        if not fs or not all("label" in f for f in fs):
            problems.append(f"parse_findings({real[0]}) bad")

    # --- 3. synthetic dataset: three classes + attribution -----------------
    td = tempfile.mkdtemp()
    try:
        tdp = Path(td)
        spike_root = tdp / "spike"
        eval_root = tdp / "eval"
        case = "syn"
        rep_dir = spike_root / case / "rep1"
        ctrl_dir = spike_root / case / "control"
        rep_dir.mkdir(parents=True)
        ctrl_dir.mkdir(parents=True)
        run_dir = eval_root / "t" / case / "rep1"
        crun_dir = eval_root / "t" / case / "control"
        run_dir.mkdir(parents=True)
        crun_dir.mkdir(parents=True)

        orig = "ACGT" * 5000                       # 20000 bp
        _write_fa(ctrl_dir / "new.fa", "chr1", orig)
        _write_fa(rep_dir / "spiked.fa", "chr1", orig[:5000] + orig[6000:])
        (ctrl_dir / "build.json").write_text(json.dumps({"ref": str(ctrl_dir / "new.fa")}))
        (rep_dir / "events.tsv").write_text(
            "event_id\tkind\tsize\tcontig\tstart0\tend0\tbreakpoints\torig_start0\torig_end0\n"
            "del_1000_0\tdel\t1000\tchr1\t5000\t5000\t5000\t5000\t6000\n")
        (rep_dir / "summary.json").write_text(json.dumps({"contig_lengths": {"chr1": 19000}}))

        def seq(n):
            return "A" * n

        # control (original assembly): overlap with [5000,6000)
        # r_union: two disjoint 60 bp alignments (primary + supplementary); the
        # union overlap is 120 bp, so it must be must-find (cap = 100), not neutral.
        _write_bam(ctrl_dir / "mapped.bam", [("chr1", 20000)], [
            {"name": "r_must", "pos0": 5100, "cigar": "400M", "seq": seq(400)},   # 400
            {"name": "r_must2", "pos0": 4800, "cigar": "500M", "seq": seq(500)},  # 300
            {"name": "r_union", "pos0": 5100, "cigar": "60M", "seq": seq(60)},    # 60
            {"name": "r_union", "pos0": 5300, "cigar": "60M", "seq": seq(60),
             "flag": 0x800},                                                     # +60
            {"name": "r_neutral", "pos0": 4980, "cigar": "50M", "seq": seq(50)},  # 30
            {"name": "r_zero", "pos0": 4000, "cigar": "500M", "seq": seq(500)},   # 0
            {"name": "r_far", "pos0": 5200, "cigar": "500M", "seq": seq(500)},    # 500
            {"name": "r_longclip", "pos0": 5000, "cigar": "5000M", "seq": seq(5000)},  # 1000
        ])
        # spiked: r_must unmapped; r_union/r_far map elsewhere and are not selected
        _write_bam(rep_dir / "mapped.bam", [("chr1", 19000)], [
            {"name": "r_must2", "pos0": 4700, "cigar": "500M", "seq": seq(500)},
            {"name": "r_neutral", "pos0": 4700, "cigar": "50M", "seq": seq(50)},
            {"name": "r_zero", "pos0": 3500, "cigar": "500M", "seq": seq(500)},
            {"name": "r_longclip", "pos0": 4700, "cigar": "5000M", "seq": seq(5000)},
            {"name": "r_union", "pos0": 16000, "cigar": "120M", "seq": seq(120)},
            {"name": "r_must", "flag": 0x4, "refid": -1, "pos0": -1, "cigar": None,
             "seq": seq(400)},
            {"name": "r_far", "pos0": 15000, "cigar": "500M", "seq": seq(500)},
            # footprint 4900-5000 touches finding_1's matching row (4950-5050) but not
            # the planted event (5000-6000): new rule -> not wrong, wrong_near_event.
            {"name": "r_gap", "pos0": 4900, "cigar": "100M", "seq": seq(100)},
            # footprint 100-200 touches only finding_1's wrapped far row (0-1000), the
            # row that does NOT match the event: also not wrong under the new rule.
            {"name": "r_wrap", "pos0": 100, "cigar": "100M", "seq": seq(100)},
        ])

        for base, reads_total, bases_total in ((run_dir, 100, 100000),
                                               (crun_dir, 200, 200000)):
            # finding_1 is a wrapped finding: two rows with one id, only the first
            # matches the planted event.
            (base / "qc.findings.bed").write_text(
                "#contig\tstart\tend\tlabel\tscore\ttier\tlikely\textra\n"
                "chr1\t4950\t5050\tfinding_1\t5\tred\tmissing_sequence\t"
                "evidence=large_del;members=large_del:1;lowcov=0\n"
                "chr1\t0\t1000\tfinding_1\t5\tred\tmissing_sequence\t"
                "evidence=large_del;members=large_del:1;lowcov=0\n")
            (base / "qc.json").write_text(json.dumps(
                {"assembly": {"length_bp": 19000},
                 "reads": {"total": reads_total, "bases_total": bases_total}}))
            (base / "qc.haplotypes.tsv").write_text("group\tread\n")
            for e in (".vcf", ".bed"):
                (base / ("qc" + e)).write_text("")
        # spiked fish.tsv: the `length` column is the EVENT length for clip/indel
        # rows (as the real tool writes it), so bp must come from fish.json.  r_must
        # is selected only via an unmapped row; r_longclip is a 5 kb read selected by
        # a 60 bp clip (tests that bp uses the read length, not the clip length).
        (run_dir / "fish.fish.tsv").write_text(
            FISH_HEADER
            + "r_must2\t1000\tlarge_del\tfinding_1\tchr1\t5000\t+\tprimary\tlen=1000;at=5000\t.\n"
            + "r_neutral\t1\tlocus\tfinding_1\tchr1\t5000\t+\tprimary\t.\t.\n"
            + "r_zero\t60\tclip\tfinding_1\tchr1\t5000\t+\tprimary\tside=left;len=60\t.\n"
            + "r_must\t400\tunmapped\t.\t.\t.\t.\tprimary\t.\t.\n"
            + "r_longclip\t60\tclip\tfinding_1\tchr1\t5000\t+\tprimary\tside=left;len=60\t.\n"
            + "r_gap\t60\tclip\tfinding_1\tchr1\t5000\t+\tprimary\tside=left;len=60\t.\n"
            + "r_wrap\t60\tclip\tfinding_1\tchr1\t150\t+\tprimary\tside=left;len=60\t.\n")
        (run_dir / "fish.fastq.gz").write_bytes(b"")
        (run_dir / "fish.fish.json").write_text(json.dumps({"reads": {
            "reads": 5, "bp": 6450, "by_reason": {
                "clip": {"reads": 2, "bp": 5500},
                "locus": {"reads": 1, "bp": 50},
                "large_del": {"reads": 1, "bp": 500},
                "unmapped": {"reads": 1, "bp": 400}}}}))
        # control fish.tsv: just some noise
        (crun_dir / "fish.fish.tsv").write_text(
            FISH_HEADER + "r_noise\t60\tclip\tfinding_1\tchr1\t10000\t+\tprimary\t"
            "side=left;len=60\t.\n")
        (crun_dir / "fish.fastq.gz").write_bytes(b"")
        (crun_dir / "fish.fish.json").write_text(json.dumps({"reads": {
            "reads": 1, "bp": 300, "by_reason": {"clip": {"reads": 1, "bp": 300}}}}))

        sample = {"case": case, "dir": str(ctrl_dir.parent), "platform": "ont",
                  "preset": "lr:hq"}
        reps_res, ctrl_res, missing, warns = score_dataset(
            spike_root, eval_root, "t", [sample])
        if missing:
            problems.append(f"synthetic missing: {missing}")
        if len(reps_res) != 1 or len(ctrl_res) != 1:
            problems.append(f"synthetic scored {len(reps_res)}/{len(ctrl_res)}")
        else:
            ev = reps_res[0]["events"][0]
            # r_zero (3500-4000) touches neither the event nor finding_1 -> wrong in
            # both rules; r_gap (4900-5000) touches finding_1's matching row and
            # r_wrap (100-200) touches only finding_1's *wrapped* far row -> both are
            # not wrong, but still wrong_near_event (the historic rule).
            expect = {"n_must": 5, "recall_near": 2, "recall_total": 3,
                      "sel_must": 2, "sel_neutral": 1, "sel_wrong": 1,
                      "wrong_near_event": 3,
                      "far_not_selected": 2, "far_unmapped": 0,
                      "far_elsewhere": 2}
            for k, v in expect.items():
                if ev.get(k) != v:
                    problems.append(f"synthetic event {k}={ev.get(k)} != {v}")
            if ev["reason_must"].get("unmapped") != 1:
                problems.append(f"reason_must unmapped {ev['reason_must']}")
            if ev["reason_must"].get("large_del") != 1:
                problems.append(f"reason_must large_del {ev['reason_must']}")
            if ev["reason_must"].get("clip") != 1:
                problems.append(f"reason_must clip {ev['reason_must']}")
            if ev["far_elsewhere_contigs"] != "chr1:2":
                problems.append(f"far_elsewhere_contigs {ev['far_elsewhere_contigs']}")
            gw = reps_res[0]["gw"]
            # bp comes from fish.json (6450), NOT the fish.tsv `length` sum (1521)
            # selected_reads/bp come from fish.json; must/not-must from the rows, so
            # r_gap and r_wrap (rows, not must-finds) move not-must 2 -> 4.
            if (gw["selected_reads"], gw["selected_bp"], gw["selected_must"],
                    gw["selected_not_must"]) != (5, 6450, 3, 4):
                problems.append(f"synthetic gw {gw}")
            if gw["by_reason"]["clip"] != {"reads": 2, "bp": 5500}:
                problems.append(f"synthetic gw clip {gw['by_reason']['clip']}")
            if gw["by_reason"]["unmapped"] != {"reads": 1, "bp": 400}:
                problems.append(f"synthetic gw unmapped {gw['by_reason']['unmapped']}")
            cw = ctrl_res[0]["gw"]
            if cw["selected_reads"] != 1 or cw["selected_must"] != 0:
                problems.append(f"synthetic control gw {cw}")

        # recall_near must NOT count the unmapped-caught must-find read
        rows = reps_res[0]["events"]
        if agg_events(rows)["recall_near"] >= agg_events(rows)["recall_total"]:
            problems.append("recall_near not < recall_total (unmapped attribution)")

        out_dir = tdp / "out"
        write_outputs(out_dir, "t", reps_res, ctrl_res, missing, warns)
        for f in ("per_event.tsv", "genome_wide.tsv", "far_not_selected.tsv",
                  "summary.md"):
            if not (out_dir / f).exists():
                problems.append(f"write_outputs missing {f}")
        sm = (out_dir / "summary.md").read_text()
        if "TRUTH DEFINITION" not in sm or "recall_total" not in sm:
            problems.append("summary.md missing truth definition")
        if "no_bridge" in sm:
            problems.append("summary.md still mentions no_bridge")
        pe = (out_dir / "per_event.tsv").read_text().splitlines()
        if len(pe) != 2:
            problems.append(f"per_event.tsv rows {len(pe)} != 2")
        far = (out_dir / "far_not_selected.tsv").read_text().splitlines()
        if not any("r_far" in l and "chr1:15001" in l for l in far):
            problems.append(f"far_not_selected.tsv missing r_far location: {far}")
        if not any("r_union" in l and "chr1:16001" in l for l in far):
            problems.append(f"far_not_selected.tsv missing r_union location: {far}")

        # --- 4. missing fish.tsv is an error unless --allow-missing-fish ---
        run2 = eval_root / "t2" / case / "rep1"
        crun2 = eval_root / "t2" / case / "control"
        run2.mkdir(parents=True)
        crun2.mkdir(parents=True)
        for src, dst in ((run_dir, run2), (crun_dir, crun2)):
            for f in src.iterdir():
                if f.name != "fish.fish.tsv":
                    (dst / f.name).write_bytes(f.read_bytes())
        rr, cc, miss2, _w = score_dataset(spike_root, eval_root, "t2", [sample])
        if not any("missing fish.tsv" in m for m in miss2) or rr or cc:
            problems.append(f"missing fish.tsv not an error: {miss2}")
        rr2, cc2, _m2, _w2 = score_dataset(spike_root, eval_root, "t2", [sample],
                                           allow_missing_fish=True)
        if not rr2 or not cc2:
            problems.append("allow_missing_fish did not score")

        # --- 5. bp falls back to FASTQ read lengths when fish.json is absent --
        run3 = eval_root / "t3" / case / "rep1"
        crun3 = eval_root / "t3" / case / "control"
        run3.mkdir(parents=True)
        crun3.mkdir(parents=True)
        for src, dst in ((run_dir, run3), (crun_dir, crun3)):
            for f in src.iterdir():
                if f.name != "fish.fish.json":
                    (dst / f.name).write_bytes(f.read_bytes())
        with gzip.open(run3 / "fish.fastq.gz", "wt") as fh:
            for name, n in (("r_must2", 500), ("r_neutral", 50), ("r_zero", 500),
                            ("r_must", 400), ("r_longclip", 5000)):
                fh.write(f"@{name}\n{'A' * n}\n+\n{'I' * n}\n")
        with gzip.open(crun3 / "fish.fastq.gz", "wt") as fh:
            fh.write("@r_noise\n" + "A" * 300 + "\n+\n" + "I" * 300 + "\n")
        rr3, _cc3, miss3, _w3 = score_dataset(spike_root, eval_root, "t3", [sample])
        if not rr3:
            problems.append(f"fastq fallback did not score: {miss3}")
        elif rr3[0]["gw"]["selected_bp"] != 6450:
            problems.append(
                f"fastq-fallback bp {rr3[0]['gw']['selected_bp']} != 6450")
        elif rr3[0]["gw"]["by_reason"]["clip"]["bp"] != 5500:
            problems.append(
                f"fastq-fallback clip bp {rr3[0]['gw']['by_reason']['clip']}")
    finally:
        shutil.rmtree(td, ignore_errors=True)

    if problems:
        print("SELF-TEST FAILED:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 1
    print("self-test OK: fish.tsv parsing, truth classes (must/neutral/wrong) with "
          "union overlap across records, wrong excludes reads on a matched "
          "finding's spans incl. a wrapped finding's other row "
          "(wrong_near_event keeps them), unmapped != finding "
          "attribution (recall_near < recall_total), far-not-selected locations, bp "
          "from fish.json (read lengths) and the FASTQ fallback, report writing")
    return 0


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="store_true",
                    help="run cramhook qc then cramhook fish; without it, score "
                         "existing outputs")
    ap.add_argument("--tag", default="f1")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--samples", default=None, help="comma-separated case substrings")
    ap.add_argument("--replicates", default=None,
                    help="comma-separated replicate numbers (default: all found)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--extra-args", default=None,
                    help="extra arguments appended to every `cramhook qc` call")
    ap.add_argument("--fish-args", default=None,
                    help="extra arguments appended to the optional second "
                         "`cramhook fish` call (e.g. '--scope genome --locus')")
    ap.add_argument("--fish-label", default=None,
                    help="subdir label for the second fish config (default: alt)")
    ap.add_argument("--allow-missing-fish", action="store_true",
                    help="score runs whose fish.tsv is missing as zero rows "
                         "(default: treat it as an error and skip the run)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test(args)
    only = [s for s in args.samples.split(",") if s] if args.samples else None
    reps = [int(x) for x in args.replicates.split(",") if x] if args.replicates else []
    (EVAL_ROOT / args.tag).mkdir(parents=True, exist_ok=True)
    if args.run:
        return run(only, reps, args.threads, args.force, args.tag, args.extra_args,
                   args.fish_args, args.fish_label, args.allow_missing_fish)
    label2 = _second_label(args.fish_args, args.fish_label)
    labels = [None] + ([label2] if label2 else [])
    reps_res, ctrl_res, missing, warns = score_dataset(
        SPIKE_ROOT, EVAL_ROOT, args.tag, qs.SAMPLES, reps, only, labels,
        args.allow_missing_fish)
    write_outputs(EVAL_ROOT / args.tag, args.tag, reps_res, ctrl_res, missing, warns)
    return 0


if __name__ == "__main__":
    sys.exit(main())
