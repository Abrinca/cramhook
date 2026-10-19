#!/usr/bin/env python3
"""Evaluate `cramhook qc` calls against the truth from scripts/qc_truth.py.

Usage:
  pixi run python scripts/qc_eval.py --truth scratch/qc/truth \
      --runs scratch/qc/runs --sets m6large,real,val [--cases SUB[,SUB]] \
      [--run] [--threads 16] [--extra-args "..."] [--out scratch/qc/eval/<tag>] \
      [--tag NAME] [--self-test]

With ``--run`` the tool is executed per case as

  pixi run /usr/bin/time -v ./target/release/cramhook qc \
      --in <case>/old.cram --ref <case>/old.fa --preset <manifest preset> \
      --out-prefix scratch/qc/runs/<tag>/<set>/<case>/qc --threads T \
      --candidates [<extra-args>]

and wall time / max RSS are recorded into ``<prefix>.run.json``.  Without
``--run`` the existing outputs under ``--runs`` are evaluated.

Matching (per Plan-v3 step 1) is positional; the event kind is not scored
(a type-concordant ``family_recall`` is reported separately):
  * truth SNP: any VCF call within +/-2 bp;
  * truth indel (INS/DEL/LINS/LDEL): any VCF call within +/-max(5, hp, str) bp;
  * a tiered BED region (RED/ORANGE/YELLOW in 0-based column 5) matches any
    truth event whose interval (start0/end0) it overlaps within +/-10 bp;
    untiered BED records (e.g. ``not_assessable``/``low_coverage``) are
    annotations, not calls.
The VCF tier is now the INFO field ``TIER=RED|ORANGE|YELLOW``; FILTER is
``PASS``/``LowCov`` (a confidence flag, ``LowCov`` = site DP below the absolute
low-coverage threshold).  ``TIER=GRAY`` (homopolymer indels, a tier below
yellow) is a reported call but **not** a scored call: it is skipped like an
untiered record (no recall, no FP), so scoring stays comparable.  A truth event
is recalled at tier T if a call of tier >= T matches it (RED > ORANGE > YELLOW),
regardless of FILTER; recall is reported over all calls and separately over PASS
calls only.  A call is a false positive if no truth event and no
``unaligned``/``multi`` region is within tolerance; calls inside such regions are
counted separately, and LowCov FPs are counted apart from PASS FPs.
"""
import argparse
import json
import re
import shlex
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SET_DIRS = {
    "m6large": ROOT / "data" / "m6large" / "real",
    "real": ROOT / "data" / "real",
    "val": ROOT / "data" / "val" / "real",
}
SET_ORDER = ["m6large", "real", "val"]
MANIFEST = "manifest.tsv"
BIN = ROOT / "target" / "release" / "cramhook"

TIER_RANK = {"RED": 3, "ORANGE": 2, "YELLOW": 1}
TIERS = set(TIER_RANK)
# GRAY = homopolymer indels, a tier below yellow: it is a reported VCF call but
# not a scored call, so it is treated like an untiered record (no recall, no FP)
# to keep scoring comparable.
GRAY_TIER = "GRAY"
TYPE_ORDER = ["SNP", "INS", "DEL", "LINS", "LDEL"]
FAMILIES = ["SNP", "INS", "DEL"]
BED_FAMILY = {"large_del": "DEL", "large_ins": "INS", "snp_cluster": "SNP"}
HP_BUCKETS = ["1", "2", "3", "4", "5", "6", "7", "8+"]

TOL_SNP = 2
TOL_INDEL_MIN = 5
TOL_BED = 10

# "truth by read support" view: one class per truth event, from the run's
# candidate AF and the unaligned/multi regions.  `unknown` when candidates.tsv
# is missing (the class depends on the run, so it is computed per run).
SUPPORT_ORDER = ["repeat_unaligned", "clear", "mixed", "minority",
                 "no_support", "unknown"]
SUPPORT_UNKNOWN = "unknown"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def read_manifest(root):
    rows = {}
    path = Path(root) / MANIFEST
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


def family(typ):
    if typ == "SNP":
        return "SNP"
    if typ in ("INS", "LINS"):
        return "INS"
    return "DEL"


def hp_bucket(hp):
    if hp <= 0:
        return "0"
    if hp >= 8:
        return "8+"
    return str(hp)


def safe_div(a, b):
    return (a / b) if b else None


def pct(x):
    return "n/a" if x is None else f"{100.0 * x:.1f}%"


# ---------------------------------------------------------------------------
# truth / calls
# ---------------------------------------------------------------------------
def load_truth(case_dir):
    events = []
    with open(case_dir / "truth.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            d = dict(zip(header, parts))
            ev = {
                "contig": d["contig"], "pos0": int(d["pos0"]), "type": d["type"],
                "len": int(d["len"]), "ref": d.get("ref", ""), "alt": d.get("alt", ""),
                "hp": int(d.get("hp", 0) or 0), "str": int(d.get("str", 0) or 0),
            }
            s0, e0 = d.get("start0", ""), d.get("end0", "")
            if s0 not in (None, "") and e0 not in (None, ""):
                ev["start0"] = int(s0)
                ev["end0"] = int(e0)
            events.append(ev)
    regions = []
    rpath = case_dir / "regions.bed"
    if rpath.exists():
        for line in rpath.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 4:
                continue
            label = f[3]
            if label in ("unaligned", "multi"):
                regions.append((f[0], int(f[1]), int(f[2]), label))
    summary = {}
    spath = case_dir / "summary.json"
    if spath.exists():
        try:
            summary = json.loads(spath.read_text())
        except Exception:
            summary = {}
    return events, regions, summary


def parse_vcf(path):
    calls = []
    if not Path(path).exists():
        return calls
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 8:
            continue
        # The tier lives in INFO (TIER=...); FILTER is PASS or LowCov.
        m = re.search(r"(?:^|;)TIER=([A-Za-z]+)", f[7])
        if not m:
            continue
        tier = m.group(1).upper()
        # GRAY is not a scored call (below yellow): skip it like an untiered
        # record so recall and FP stay comparable with runs that have no gray.
        if tier == GRAY_TIER:
            continue
        if tier not in TIERS:
            continue
        filt = f[6].split(";")[0].strip().upper()
        lowcov = filt == "LOWCOV"
        m = re.search(r"(?:^|;)TYPE=([A-Za-z]+)", f[7])
        if m:
            typ = m.group(1).upper()
        else:
            ref, alt = f[3], f[4].split(",")[0]
            typ = "INS" if len(alt) > len(ref) else ("DEL" if len(alt) < len(ref) else "SNP")
        calls.append({"source": "vcf", "contig": f[0], "pos0": int(f[1]) - 1,
                      "type": typ, "tier": tier, "lowcov": lowcov})
    return calls


def parse_bed(path):
    """Parse ``qc.bed`` into ``(calls, annotations)``.

    The tier is the value in column 6 (0-based index 5; spec order: contig,
    start, end, label, score, tier, ...).  A record whose tier is not
    RED/ORANGE/YELLOW is an annotation (e.g. ``not_assessable``), not a call.
    """
    calls, annotations = [], []
    if not Path(path).exists():
        return calls, annotations
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 4:
            f = line.split()
        if len(f) < 4:
            continue
        try:
            start, end = int(f[1]), int(f[2])
        except ValueError:
            continue
        label = f[3]
        if label == "depth_segment":
            # Output only (never a finding or an alert): not a scored call.
            continue
        tier = f[5].strip().upper() if len(f) > 5 else None
        extra = f[9] if len(f) > 9 else ""
        lowcov = "lowcov=1" in extra
        rec = {"contig": f[0], "start": start, "end": end, "label": label,
               "tier": tier, "lowcov": lowcov}
        if tier in TIERS:
            calls.append(dict(rec, source="bed"))
        else:
            annotations.append(rec)
    return calls, annotations


def parse_findings(path):
    """Parse ``qc.findings.bed``: one dict per finding (tier/likely/extra).

    The findings layer is optional: a missing file (older outputs) returns [].
    """
    out = []
    path = Path(path)
    if not path.exists():
        return out
    lines = path.read_text().splitlines()
    if not lines:
        return out
    header = lines[0].lstrip("#").split("\t")
    idx = {h: i for i, h in enumerate(header)}

    def col(f, name, default):
        i = idx.get(name, default)
        return f[i] if i < len(f) else ""

    for line in lines[1:]:
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 4:
            continue
        extra = col(f, "extra", 7)
        lowcov = "lowcov=1" in extra
        try:
            start, end = int(col(f, "start", 1)), int(col(f, "end", 2))
        except ValueError:
            continue
        out.append({"contig": col(f, "contig", 0), "start": start, "end": end,
                    "label": col(f, "label", 3), "score": int(col(f, "score", 4) or 0),
                    "tier": col(f, "tier", 5).upper(),
                    "likely": col(f, "likely", 6), "lowcov": lowcov})
    return out


def truth_interval(e):
    """Affected old interval (0-based, half-open).

    Uses the explicit ``start0``/``end0`` recorded by ``qc_truth.py`` so that a
    right-anchored origin deletion is ``[0, len)`` (not ``[1, len+1)``) and a
    whole-contig deletion ends at the contig end.  Only truth files lacking the
    columns (legacy) fall back to deriving the interval from ``pos0``/``ref``.
    """
    if "start0" in e and "end0" in e:
        return (e["start0"], e["end0"])
    if family(e["type"]) == "DEL" and e["len"] > 0:
        return (e["pos0"] + 1, e["pos0"] + 1 + e["len"])
    return (e["pos0"], e["pos0"] + 1)


def truth_tolerance(e):
    """Tolerance for a truth event (truth-dependent)."""
    if family(e["type"]) == "SNP":
        return TOL_SNP
    return max(TOL_INDEL_MIN, e["hp"], e.get("str", 0))


def call_matches_truth(call, e):
    """Positional match only (Plan-v3: the event kind is not scored)."""
    if call["contig"] != e["contig"]:
        return False
    if call["source"] == "vcf":
        return abs(call["pos0"] - e["pos0"]) <= truth_tolerance(e)
    a, b = truth_interval(e), (call["start"], call["end"])
    return a[0] - TOL_BED < b[1] and b[0] - TOL_BED < a[1]


def call_family_concordant(call, e):
    """Whether a call also agrees on the event family (optional column)."""
    fam = family(e["type"])
    if call["source"] == "vcf":
        return family(call["type"]) == fam
    return BED_FAMILY.get(call["label"]) == fam


def truth_in_regions(e, regions):
    """Whether a truth event's interval overlaps any (contig, start, end) region."""
    a = truth_interval(e)
    for r in regions:
        if r["contig"] != e["contig"]:
            continue
        if a[0] < r["end"] and r["start"] < a[1]:
            return True
    return False


def call_near_region(call, regions):
    a = (call["pos0"], call["pos0"] + 1) if call["source"] == "vcf" \
        else (call["start"], call["end"])
    for contig, s, e, label in regions:
        if contig != call["contig"]:
            continue
        if a[0] - TOL_BED < e and s - TOL_BED < a[1]:
            return label
    return None


def classify_truth(events, regions, cands):
    """Classify each truth event by read support (one class string per event).

    ``cands`` is ``{contig: [(pos0, af)]}`` from qc.candidates.tsv, or ``None``
    when the file is missing (then every event is ``unknown``).

    A class depends on the run, so it is derived from the run's candidates:
      * ``repeat_unaligned``: the event is within +/-10 bp of an
        ``unaligned``/``multi`` region of regions.bed (takes precedence);
      * otherwise by AF = max ``af`` over candidate rows within the event's
        tolerance (SNP +/-2, indel +/-max(5, hp, str)):
        ``clear`` AF > 0.6, ``mixed`` 0.4-0.6, ``minority`` 0 < AF < 0.4,
        ``no_support`` when there is no candidate within tolerance.
    """
    classes = []
    for e in events:
        if cands is None:
            classes.append(SUPPORT_UNKNOWN)
            continue
        # same half-open +/-10 bp overlap predicate as call_near_region / matching
        a = truth_interval(e)
        if any(r[0] == e["contig"] and a[0] - TOL_BED < r[2]
               and r[1] - TOL_BED < a[1] for r in regions):
            classes.append("repeat_unaligned")
            continue
        p = e["pos0"]
        tol = truth_tolerance(e)
        af = 0.0
        for q, a in cands.get(e["contig"], ()):
            if abs(q - p) <= tol and a > af:
                af = a
        if af > 0.6:
            classes.append("clear")
        elif af >= 0.4:
            classes.append("mixed")
        elif af > 0:
            classes.append("minority")
        else:
            classes.append("no_support")
    return classes


def evaluate_case(events, regions, calls_vcf, calls_bed, annotations, contig_lengths,
                  cands=None):
    """Return (res, missed, fps).

    ``cands`` is the per-contig candidate AF map (``{contig: [(pos0, af)]}``)
    used for the "truth by read support" view, or ``None`` when candidates.tsv
    is unavailable (classes then become ``unknown``).
    """
    calls = calls_vcf + calls_bed
    classes = classify_truth(events, regions, cands)
    truth_by_contig = defaultdict(list)
    for i, e in enumerate(events):
        truth_by_contig[e["contig"]].append(i)
    calls_by_contig = defaultdict(list)
    for i, c in enumerate(calls):
        calls_by_contig[c["contig"]].append(i)
    vcf_by_contig = defaultdict(list)
    for c in calls_vcf:
        vcf_by_contig[c["contig"]].append(c)

    best = [0] * len(events)
    best_pass = [0] * len(events)
    best_vcf = [0] * len(events)
    best_fam = [0] * len(events)
    exact = [False] * len(events)
    for i, e in enumerate(events):
        for ci in calls_by_contig.get(e["contig"], ()):
            c = calls[ci]
            if not call_matches_truth(c, e):
                continue
            r = TIER_RANK.get(c["tier"], 0)
            if r > best[i]:
                best[i] = r
            if not c.get("lowcov") and r > best_pass[i]:
                best_pass[i] = r
            if call_family_concordant(c, e) and r > best_fam[i]:
                best_fam[i] = r
            if c["source"] == "vcf" and c["pos0"] == e["pos0"]:
                exact[i] = True
        for c in vcf_by_contig.get(e["contig"], ()):
            if not call_matches_truth(c, e):
                continue
            r = TIER_RANK.get(c["tier"], 0)
            if r > best_vcf[i]:
                best_vcf[i] = r

    is_tp = [False] * len(calls)
    for ci, c in enumerate(calls):
        for i in truth_by_contig.get(c["contig"], ()):
            if call_matches_truth(c, events[i]):
                is_tp[ci] = True
                break

    # --- truth-side aggregation -------------------------------------------
    truth_counts = {t: 0 for t in TYPE_ORDER}
    truth_hp = {fam: defaultdict(int) for fam in FAMILIES}
    recall = {t: {"red": 0, "orange": 0, "all": 0} for t in TYPE_ORDER}
    recall_pass = {t: {"red": 0, "orange": 0, "all": 0} for t in TYPE_ORDER}
    fam_recall = {t: {"red": 0, "orange": 0, "all": 0} for t in TYPE_ORDER}
    exact_counts = {t: 0 for t in TYPE_ORDER}
    # Annotations are split by BED label: only ``not_assessable`` regions count
    # as "not assessable" (e.g. ``low_coverage`` is a separate annotation).
    ann_by_label = defaultdict(list)
    for a in annotations:
        ann_by_label[a["label"]].append(a)
    in_not_assessable = 0
    truth_in_annotation = defaultdict(int)
    missed = []
    for i, e in enumerate(events):
        t = e["type"]
        truth_counts[t] = truth_counts.get(t, 0) + 1
        truth_hp[family(t)][hp_bucket(e["hp"])] += 1
        for name, thresh in (("red", 3), ("orange", 2), ("all", 1)):
            if best[i] >= thresh:
                recall[t][name] += 1
            if best_pass[i] >= thresh:
                recall_pass[t][name] += 1
            if best_fam[i] >= thresh:
                fam_recall[t][name] += 1
        if exact[i]:
            exact_counts[t] = exact_counts.get(t, 0) + 1
        if best[i] == 0:
            missed.append(dict(e, support_class=classes[i]))
        for label, regs in ann_by_label.items():
            if truth_in_regions(e, regs):
                truth_in_annotation[label] += 1
        if truth_in_regions(e, ann_by_label.get("not_assessable", [])):
            in_not_assessable += 1

    # --- call-side aggregation --------------------------------------------
    fp = {"RED": 0, "ORANGE": 0, "YELLOW": 0}
    fp_lowcov = {"RED": 0, "ORANGE": 0, "YELLOW": 0}
    in_region = {"vcf": 0, "bed": 0}
    fps = []
    for ci, c in enumerate(calls):
        if is_tp[ci]:
            continue
        if call_near_region(c, regions) is not None:
            in_region[c["source"]] += 1
            continue
        tier = c["tier"] or "YELLOW"
        fp[tier] = fp.get(tier, 0) + 1
        if c.get("lowcov"):
            fp_lowcov[tier] = fp_lowcov.get(tier, 0) + 1
        fps.append(dict(c))

    genome_bp = sum(contig_lengths.values())
    fp_total = sum(fp.values())
    fp_lowcov_total = sum(fp_lowcov.values())
    total_truth = sum(truth_counts.values())

    def overall(rec):
        return {k: safe_div(sum(rec[t][k] for t in TYPE_ORDER), total_truth)
                for k in ("red", "orange", "all")}

    def support_table(best_arr):
        """Per-class truth counts and cumulative recall, by read support."""
        tbl = {c: {"n": 0, "red": 0, "red_orange": 0, "any": 0}
               for c in SUPPORT_ORDER}
        for i in range(len(events)):
            t = tbl[classes[i]]
            t["n"] += 1
            if best_arr[i] >= 3:
                t["red"] += 1
            if best_arr[i] >= 2:
                t["red_orange"] += 1
            if best_arr[i] >= 1:
                t["any"] += 1
        return tbl

    res = {
        "truth_counts": truth_counts,
        "truth_total": total_truth,
        "truth_hp": {fam: dict(truth_hp[fam]) for fam in FAMILIES},
        "recall": recall,
        "recall_overall": overall(recall),
        "recall_pass": recall_pass,
        "recall_pass_overall": overall(recall_pass),
        "family_recall": fam_recall,
        "family_recall_overall": overall(fam_recall),
        "exact": exact_counts,
        "exact_total": sum(exact_counts.values()),
        "exact_overall": safe_div(sum(exact_counts.values()), total_truth),
        "fp": fp,
        "fp_total": fp_total,
        "fp_lowcov": fp_lowcov,
        "fp_lowcov_total": fp_lowcov_total,
        "fp_per_mb": safe_div(fp_total, genome_bp / 1e6) if genome_bp else None,
        "in_region": in_region,
        "in_region_total": in_region["vcf"] + in_region["bed"],
        "calls": {"vcf": len(calls_vcf), "bed": len(calls_bed)},
        "annotations": len(annotations),
        "annotation_bp": sum(a["end"] - a["start"] for a in annotations),
        "annotation_bp_by_label": {
            label: sum(a["end"] - a["start"] for a in regs)
            for label, regs in ann_by_label.items()},
        "annotation_counts": {label: len(regs) for label, regs in ann_by_label.items()},
        "not_assessable_bp": sum(
            a["end"] - a["start"] for a in ann_by_label.get("not_assessable", [])),
        "truth_in_not_assessable": in_not_assessable,
        "truth_in_annotation": dict(truth_in_annotation),
        "support": support_table(best),
        "support_vcf": support_table(best_vcf),
        "support_pass": support_table(best_pass),
        "support_classes": classes,
        "genome_bp": genome_bp,
    }
    return res, missed, fps


# ---------------------------------------------------------------------------
# running the tool
# ---------------------------------------------------------------------------
def parse_hms(value):
    """Parse GNU time's ``[h:]mm:ss[.ff]`` elapsed value into seconds."""
    parts = value.strip().split(":")
    try:
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        if len(parts) == 2:
            return int(parts[0]) * 60 + float(parts[1])
        return float(parts[0])
    except ValueError:
        return None


def parse_time_verbose(stderr):
    wall = rss = None
    # GNU time: "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.23" — the
    # label itself contains colons, so match the parenthesised unit list first.
    m = re.search(r"Elapsed \(wall clock\) time \([^)]*\):\s*([0-9:.]+)", stderr)
    if m:
        wall = parse_hms(m.group(1))
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", stderr)
    if m:
        rss = int(m.group(1))
    return wall, rss


def run_case(set_name, case, preset, runs_root, tag, threads, extra_args, binary):
    case_dir = SET_DIRS[set_name] / case
    out_dir = runs_root / tag / set_name / case
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "qc"
    # delete stale outputs so a failed rerun cannot be scored from old files
    for f in out_dir.glob("qc.*"):
        try:
            f.unlink()
        except OSError:
            pass
    inner = ["/usr/bin/time", "-v", str(binary), "qc",
             "--in", str(case_dir / "old.cram"), "--ref", str(case_dir / "old.fa"),
             "--preset", preset, "--out-prefix", str(prefix),
             "--threads", str(threads), "--candidates"]
    if extra_args:
        inner += shlex.split(extra_args)
    t0 = time.monotonic()
    p = subprocess.run(["pixi", "run"] + inner, cwd=str(ROOT),
                       capture_output=True, text=True)
    wall_mono = time.monotonic() - t0
    wall, rss = parse_time_verbose(p.stderr)
    # /usr/bin/time -v writes its report as tab-indented lines on stderr; drop
    # them so the recorded tail is the tool's own output, and flatten newlines.
    clean = "\n".join(l for l in p.stderr.splitlines() if not l.startswith("\t")).strip()
    tail = " ".join(clean.split())[-400:]
    meta = {"rc": p.returncode, "wall_s": wall if wall is not None else round(wall_mono, 3),
            "max_rss_kb": rss, "stderr_tail": tail if p.returncode != 0 else ""}
    try:
        Path(str(prefix) + ".run.json").write_text(json.dumps(meta, indent=2) + "\n")
    except OSError:
        pass
    return meta


def read_run_meta(prefix):
    p = Path(str(prefix) + ".run.json")
    if p.exists():
        try:
            return json.loads(p.read_text())
        except Exception:
            return {}
    return {}


# ---------------------------------------------------------------------------
# nearest candidate (qc.candidates.tsv)
# ---------------------------------------------------------------------------
def load_candidates(path):
    path = Path(path)
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    if not lines:
        return []
    header = lines[0].rstrip("\n").split("\t")
    ci = pi = pi0 = None
    for idx, h in enumerate(header):
        hl = h.strip().lower()
        if ci is None and hl in ("contig", "chrom", "chromosome", "ref", "reference"):
            ci = idx
        if pi0 is None and hl in ("pos0", "start"):
            pi0 = idx
        if pi is None and hl in ("pos", "position"):
            pi = idx
    if pi0 is None and pi is None:
        return []
    out = []
    for line in lines[1:]:
        f = line.rstrip("\n").split("\t")
        try:
            if pi0 is not None:
                pos0 = int(float(f[pi0]))
            else:
                pos0 = int(float(f[pi])) - 1
        except (ValueError, IndexError):
            continue
        contig = f[ci] if ci is not None and ci < len(f) else None
        out.append((contig, pos0, line))
    return out


def nearest_candidate(cands, contig, pos0):
    best = None
    for c, p, line in cands:
        if c is not None and contig is not None and c != contig:
            continue
        d = abs(p - pos0)
        if best is None or d < best[0]:
            best = (d, line)
    return best


def load_candidates_af(path):
    """``{contig: [(pos0, af)]}`` from qc.candidates.tsv.

    Returns ``None`` (not an error) when the file is missing, empty/headerless,
    or lacks the required columns, so the by-support view falls back to
    ``unknown``.  A valid header-only file returns ``{}`` (no candidates).
    ``pos0``/``start`` are 0-based and preferred; ``pos``/``position`` are
    1-based and decremented, like ``load_candidates``.
    """
    path = Path(path)
    if not path.exists():
        return None
    lines = path.read_text().splitlines()
    if not lines:
        return None
    header = lines[0].rstrip("\n").split("\t")
    ci = pi0 = pi1 = ai = None
    for idx, h in enumerate(header):
        hl = h.strip().lower()
        if ci is None and hl in ("contig", "chrom", "chromosome", "ref", "reference"):
            ci = idx
        if pi0 is None and hl in ("pos0", "start"):
            pi0 = idx
        if pi1 is None and hl in ("pos", "position"):
            pi1 = idx
        if ai is None and hl == "af":
            ai = idx
    if ci is None or (pi0 is None and pi1 is None) or ai is None:
        return None
    out = defaultdict(list)
    for line in lines[1:]:
        f = line.rstrip("\n").split("\t")
        try:
            if pi0 is not None:
                pos0 = int(float(f[pi0]))
            else:
                pos0 = int(float(f[pi1])) - 1
            out[f[ci]].append((pos0, float(f[ai])))
        except (ValueError, IndexError):
            continue
    return out


# ---------------------------------------------------------------------------
# aggregation / reporting
# ---------------------------------------------------------------------------
def summarise_results(results):
    """Aggregate per (set, platform); invalid cases are excluded."""
    groups = {}
    for r in results:
        if not r.get("valid", True):
            continue
        key = (r["set"], r.get("platform", ""))
        g = groups.setdefault(key, {
            "cases": 0, "truth_counts": {t: 0 for t in TYPE_ORDER},
            "truth_hp": {fam: defaultdict(int) for fam in FAMILIES},
            "recall": {t: {"red": 0, "orange": 0, "all": 0} for t in TYPE_ORDER},
            "recall_pass": {t: {"red": 0, "orange": 0, "all": 0} for t in TYPE_ORDER},
            "fam_recall": {t: {"red": 0, "orange": 0, "all": 0} for t in TYPE_ORDER},
            "exact": 0, "exact_t": {t: 0 for t in TYPE_ORDER}, "truth_total": 0,
            "fp": {"RED": 0, "ORANGE": 0, "YELLOW": 0}, "fp_total": 0,
            "fp_lowcov": {"RED": 0, "ORANGE": 0, "YELLOW": 0}, "fp_lowcov_total": 0,
            "genome_bp": 0, "in_region": 0, "calls": 0,
            "annotations": 0, "annotation_bp": 0, "annotation_bp_by_label": {},
            "annotation_counts": {}, "not_assessable_bp": 0, "truth_in_na": 0,
            "truth_in_annotation": {},
            "support": {c: {"n": 0, "red": 0, "red_orange": 0, "any": 0}
                        for c in SUPPORT_ORDER},
            "support_vcf": {c: {"n": 0, "red": 0, "red_orange": 0, "any": 0}
                            for c in SUPPORT_ORDER},
            "support_pass": {c: {"n": 0, "red": 0, "red_orange": 0, "any": 0}
                             for c in SUPPORT_ORDER},
            "walls": [], "rss": [],
            "findings_total": 0,
            "findings_tier": {"RED": 0, "ORANGE": 0, "YELLOW": 0},
            "findings_likely": defaultdict(int),
        })
        g["cases"] += 1
        g["findings_total"] += r.get("findings_total", 0)
        for t, n in r.get("findings_tier", {}).items():
            g["findings_tier"][t] = g["findings_tier"].get(t, 0) + n
        for lk, n in r.get("findings_likely", {}).items():
            g["findings_likely"][lk] += n
        g["truth_total"] += r["truth_total"]
        g["exact"] += r["exact_total"]
        for t in TYPE_ORDER:
            g["truth_counts"][t] += r["truth_counts"].get(t, 0)
            g["exact_t"][t] += r["exact"].get(t, 0)
        for fam in FAMILIES:
            for b, n in r["truth_hp"].get(fam, {}).items():
                g["truth_hp"][fam][b] += n
        for t in TYPE_ORDER:
            for tier in ("red", "orange", "all"):
                g["recall"][t][tier] += r["recall"].get(t, {}).get(tier, 0)
                g["recall_pass"][t][tier] += r["recall_pass"].get(t, {}).get(tier, 0)
                g["fam_recall"][t][tier] += r["family_recall"].get(t, {}).get(tier, 0)
        for tier in g["fp"]:
            g["fp"][tier] += r["fp"].get(tier, 0)
            g["fp_lowcov"][tier] += r.get("fp_lowcov", {}).get(tier, 0)
        g["fp_total"] += r["fp_total"]
        g["fp_lowcov_total"] += r.get("fp_lowcov_total", 0)
        g["genome_bp"] += r["genome_bp"]
        g["in_region"] += r["in_region_total"]
        g["calls"] += r["calls"]["vcf"] + r["calls"]["bed"]
        g["annotations"] += r.get("annotations", 0)
        g["annotation_bp"] += r.get("annotation_bp", 0)
        for label, bp in r.get("annotation_bp_by_label", {}).items():
            g["annotation_bp_by_label"][label] = g["annotation_bp_by_label"].get(label, 0) + bp
        for label, n in r.get("annotation_counts", {}).items():
            g["annotation_counts"][label] = g["annotation_counts"].get(label, 0) + n
        g["not_assessable_bp"] += r.get("not_assessable_bp", 0)
        g["truth_in_na"] += r.get("truth_in_not_assessable", 0)
        for label, n in r.get("truth_in_annotation", {}).items():
            g["truth_in_annotation"][label] = g["truth_in_annotation"].get(label, 0) + n
        for cls in SUPPORT_ORDER:
            for key in ("n", "red", "red_orange", "any"):
                g["support"][cls][key] += r.get("support", {}).get(cls, {}).get(key, 0)
                g["support_vcf"][cls][key] += r.get("support_vcf", {}).get(cls, {}).get(key, 0)
                g["support_pass"][cls][key] += r.get("support_pass", {}).get(cls, {}).get(key, 0)
        if r.get("wall_s") is not None:
            g["walls"].append(r["wall_s"])
        if r.get("max_rss_kb") is not None:
            g["rss"].append(r["max_rss_kb"])
    return groups


def write_per_case(results, out_dir):
    cols = (["set", "case", "platform", "preset", "valid"]
            + [f"truth_{t}" for t in TYPE_ORDER] + ["truth_total"]
            + ["recall_red", "recall_orange", "recall_all",
               "recall_pass_red", "recall_pass_orange", "recall_pass_all",
               "family_recall_all", "exact_total", "exact_frac"]
            + ["fp_RED", "fp_ORANGE", "fp_YELLOW", "fp_total", "fp_lowcov", "fp_per_mb",
               "in_region_vcf", "in_region_bed", "calls_vcf", "calls_bed",
               "annotations", "annotation_bp", "not_assessable_bp",
               "truth_in_not_assessable"]
            + [f"sup_{c}" for c in SUPPORT_ORDER]
            + ["genome_bp", "wall_s", "max_rss_kb", "error"])
    lines = ["\t".join(cols)]
    for r in results:
        rc = r["recall_overall"]
        rcp = r["recall_pass_overall"]
        row = ([r["set"], r["case"], r.get("platform", ""), r.get("preset", ""),
                1 if r.get("valid", True) else 0]
               + [r["truth_counts"].get(t, 0) for t in TYPE_ORDER]
               + [r["truth_total"],
                  f"{rc['red']:.4f}" if rc["red"] is not None else "n/a",
                  f"{rc['orange']:.4f}" if rc["orange"] is not None else "n/a",
                  f"{rc['all']:.4f}" if rc["all"] is not None else "n/a",
                  f"{rcp['red']:.4f}" if rcp["red"] is not None else "n/a",
                  f"{rcp['orange']:.4f}" if rcp["orange"] is not None else "n/a",
                  f"{rcp['all']:.4f}" if rcp["all"] is not None else "n/a",
                  f"{r['family_recall_overall']['all']:.4f}"
                  if r["family_recall_overall"]["all"] is not None else "n/a",
                  r["exact_total"],
                  f"{r['exact_overall']:.4f}" if r["exact_overall"] is not None else "n/a",
                  r["fp"]["RED"], r["fp"]["ORANGE"], r["fp"]["YELLOW"], r["fp_total"],
                  r.get("fp_lowcov_total", 0),
                  f"{r['fp_per_mb']:.2f}" if r["fp_per_mb"] is not None else "n/a",
                  r["in_region"]["vcf"], r["in_region"]["bed"],
                  r["calls"]["vcf"], r["calls"]["bed"],
                  r.get("annotations", 0), r.get("annotation_bp", 0),
                  r.get("not_assessable_bp", 0),
                  r.get("truth_in_not_assessable", 0)]
               + [r.get("support", {}).get(c, {}).get("n", 0) for c in SUPPORT_ORDER]
               + [r["genome_bp"],
                  f"{r['wall_s']:.2f}" if r.get("wall_s") is not None else "",
                  r.get("max_rss_kb") if r.get("max_rss_kb") is not None else "",
                  r.get("error", "")])
        lines.append("\t".join(str(x) for x in row))
    (out_dir / "per_case.tsv").write_text("\n".join(lines) + "\n")


def write_summary_tsv(groups, out_dir):
    cols = (["set", "platform", "cases"]
            + [f"truth_{t}" for t in TYPE_ORDER] + ["truth_total"]
            + [f"hp_{b}" for b in HP_BUCKETS]
            + [f"recall_red_{t}" for t in TYPE_ORDER]
            + [f"recall_orange_{t}" for t in TYPE_ORDER]
            + [f"recall_all_{t}" for t in TYPE_ORDER]
            + [f"recall_pass_red_{t}" for t in TYPE_ORDER]
            + [f"recall_pass_orange_{t}" for t in TYPE_ORDER]
            + [f"recall_pass_all_{t}" for t in TYPE_ORDER]
            + [f"family_recall_all_{t}" for t in TYPE_ORDER]
            + ["exact_total", "exact_frac"]
            + ["fp_RED", "fp_ORANGE", "fp_YELLOW", "fp_total", "fp_lowcov", "fp_per_mb",
               "in_region", "calls", "annotations", "annotation_bp",
               "not_assessable_bp", "truth_in_not_assessable",
               "wall_median_s", "wall_max_s", "max_rss_kb"])
    lines = ["\t".join(cols)]
    for (set_name, platform), g in sorted(groups.items()):
        tt = sum(g["truth_counts"].values())
        hp_all = defaultdict(int)
        for fam in FAMILIES:
            for b, n in g["truth_hp"][fam].items():
                hp_all[b] += n
        fp_per_mb = safe_div(g["fp_total"], g["genome_bp"] / 1e6) if g["genome_bp"] else None
        walls = g["walls"]
        row = ([set_name, platform, g["cases"]]
               + [g["truth_counts"][t] for t in TYPE_ORDER] + [tt]
               + [hp_all.get(b, 0) for b in HP_BUCKETS]
               + [f"{safe_div(g['recall'][t]['red'], g['truth_counts'][t]):.4f}"
                  if g["truth_counts"][t] else "n/a" for t in TYPE_ORDER]
               + [f"{safe_div(g['recall'][t]['orange'], g['truth_counts'][t]):.4f}"
                  if g["truth_counts"][t] else "n/a" for t in TYPE_ORDER]
               + [f"{safe_div(g['recall'][t]['all'], g['truth_counts'][t]):.4f}"
                  if g["truth_counts"][t] else "n/a" for t in TYPE_ORDER]
               + [f"{safe_div(g['recall_pass'][t]['red'], g['truth_counts'][t]):.4f}"
                  if g["truth_counts"][t] else "n/a" for t in TYPE_ORDER]
               + [f"{safe_div(g['recall_pass'][t]['orange'], g['truth_counts'][t]):.4f}"
                  if g["truth_counts"][t] else "n/a" for t in TYPE_ORDER]
               + [f"{safe_div(g['recall_pass'][t]['all'], g['truth_counts'][t]):.4f}"
                  if g["truth_counts"][t] else "n/a" for t in TYPE_ORDER]
               + [f"{safe_div(g['fam_recall'][t]['all'], g['truth_counts'][t]):.4f}"
                  if g["truth_counts"][t] else "n/a" for t in TYPE_ORDER]
               + [g["exact"],
                  f"{safe_div(g['exact'], tt):.4f}" if tt else "n/a"]
               + [g["fp"]["RED"], g["fp"]["ORANGE"], g["fp"]["YELLOW"], g["fp_total"],
                  g["fp_lowcov_total"],
                  f"{fp_per_mb:.2f}" if fp_per_mb is not None else "n/a",
                  g["in_region"], g["calls"], g["annotations"], g["annotation_bp"],
                  g["not_assessable_bp"], g["truth_in_na"],
                  f"{statistics.median(walls):.2f}" if walls else "",
                  f"{max(walls):.2f}" if walls else "",
                  max(g["rss"]) if g["rss"] else ""])
        lines.append("\t".join(str(x) for x in row))
    (out_dir / "summary.tsv").write_text("\n".join(lines) + "\n")


def write_by_support(groups, out_dir):
    """Per set x platform x support-class recall (all calls and VCF-only)."""
    def frac(k, n):
        return f"{safe_div(k, n):.4f}" if n else "n/a"

    cols = ["set", "platform", "class", "n",
            "red", "red_orange", "any",
            "n_vcf", "red_vcf", "red_orange_vcf", "any_vcf",
            "red_pass", "red_orange_pass", "any_pass"]
    lines = ["\t".join(cols)]
    for (set_name, platform), g in sorted(groups.items()):
        for cls in SUPPORT_ORDER:
            s, v = g["support"][cls], g["support_vcf"][cls]
            p = g["support_pass"][cls]
            if s["n"] == 0 and v["n"] == 0:
                continue
            lines.append("\t".join(str(x) for x in [
                set_name, platform, cls, s["n"],
                frac(s["red"], s["n"]), frac(s["red_orange"], s["n"]),
                frac(s["any"], s["n"]),
                v["n"], frac(v["red"], v["n"]), frac(v["red_orange"], v["n"]),
                frac(v["any"], v["n"]),
                frac(p["red"], p["n"]), frac(p["red_orange"], p["n"]),
                frac(p["any"], p["n"])]))
    (out_dir / "by_support.tsv").write_text("\n".join(lines) + "\n")


def write_report(results, groups, out_dir, args):
    lines = [f"# cramhook qc eval (tag `{args.tag}`)", ""]
    lines.append(f"- truth: `{args.truth}`; runs: `{args.runs}`; "
                 f"sets: {args.sets}; threads: {args.threads}")
    if args.extra_args:
        lines.append(f"- extra args: `{args.extra_args}`")
    lines.append(f"- cases: {len(results)}")

    # Headline view: truth by read support.
    lines += ["", "## truth by read support (headline)", ""]
    lines.append("Each truth event is classified by read support: max candidate AF "
                 "within its tolerance (`clear` >0.6, `mixed` 0.4-0.6, "
                 "`minority` 0<AF<0.4, `no_support`); events within +/-10 bp of an "
                 "`unaligned`/`multi` region are `repeat_unaligned` (takes precedence); "
                 "`unknown` when qc.candidates.tsv is missing. "
                 "`red`/`red+orange`/`any` are cumulative recall; VCF-only columns "
                 "restrict matching to VCF site calls; `pass` columns restrict to "
                 "FILTER=PASS calls (`LowCov` excluded).")
    lines.append("")
    lines.append("| set | platform | class | n | red | red+orange | any | "
                 "vcf red | vcf red+orange | vcf any | n_vcf | "
                 "pass red | pass red+orange | pass any |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for (set_name, platform), g in sorted(groups.items()):
        for cls in SUPPORT_ORDER:
            s, v, p = g["support"][cls], g["support_vcf"][cls], g["support_pass"][cls]
            if s["n"] == 0 and v["n"] == 0:
                continue
            lines.append(f"| {set_name} | {platform} | {cls} | {s['n']} | "
                         f"{pct(safe_div(s['red'], s['n']))} | "
                         f"{pct(safe_div(s['red_orange'], s['n']))} | "
                         f"{pct(safe_div(s['any'], s['n']))} | "
                         f"{pct(safe_div(v['red'], v['n']))} | "
                         f"{pct(safe_div(v['red_orange'], v['n']))} | "
                         f"{pct(safe_div(v['any'], v['n']))} | {v['n']} | "
                         f"{pct(safe_div(p['red'], p['n']))} | "
                         f"{pct(safe_div(p['red_orange'], p['n']))} | "
                         f"{pct(safe_div(p['any'], p['n']))} |")

    lines += ["", "## per set x platform", ""]
    lines.append("Recall columns are cumulative over all calls (regardless of "
                 "FILTER); the `pass` columns restrict to FILTER=PASS calls "
                 "(`LowCov` calls excluded), so `pass` <= `all`.")
    lines.append("")
    lines.append("| set | platform | cases | truth SNP/INS/DEL/LINS/LDEL | "
                 "recall red/orange/all | pass red/orange/all | family recall all | "
                 "exact | FP R/O/Y | FP lowcov | FP/Mb | wall med/max s | RSS MB |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for (set_name, platform), g in sorted(groups.items()):
        tt = sum(g["truth_counts"].values())
        tc = "/".join(str(g["truth_counts"][t]) for t in TYPE_ORDER)
        rec = []
        for tier in ("red", "orange", "all"):
            a = sum(g["recall"][t][tier] for t in TYPE_ORDER)
            rec.append(f"{100.0 * a / tt:.1f}%" if tt else "n/a")
        rp = []
        for tier in ("red", "orange", "all"):
            a = sum(g["recall_pass"][t][tier] for t in TYPE_ORDER)
            rp.append(f"{100.0 * a / tt:.1f}%" if tt else "n/a")
        fam = sum(g["fam_recall"][t]["all"] for t in TYPE_ORDER)
        fam_s = f"{100.0 * fam / tt:.1f}%" if tt else "n/a"
        exact_s = f"{g['exact']} ({100.0 * g['exact'] / tt:.1f}%)" if tt else "n/a"
        fp = g["fp"]
        fp_per_mb = safe_div(g["fp_total"], g["genome_bp"] / 1e6) if g["genome_bp"] else None
        walls = g["walls"]
        wm = f"{statistics.median(walls):.2f}/{max(walls):.2f}" if walls else "n/a"
        rss = f"{max(g['rss']) / 1024:.0f}" if g["rss"] else "n/a"
        fpm = f"{fp_per_mb:.2f}" if fp_per_mb is not None else "n/a"
        lines.append(f"| {set_name} | {platform} | {g['cases']} | {tc} | "
                     f"{'/'.join(rec)} | {'/'.join(rp)} | {fam_s} | {exact_s} | "
                     f"{fp['RED']}/{fp['ORANGE']}/{fp['YELLOW']} | "
                     f"{g['fp_lowcov_total']} | "
                     f"{fpm} | {wm} | {rss} |")

    lines += ["", "## recall by type (red / red+orange / all) + family-concordant all", ""]
    lines.append("`pass all` is cumulative recall over FILTER=PASS calls only; "
                 "the difference from `all` is the share recalled only by LowCov calls.")
    lines.append("")
    lines.append("| set | platform | type | truth | red | red+orange | all | "
                 "pass all | family all | exact |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for (set_name, platform), g in sorted(groups.items()):
        for t in TYPE_ORDER:
            n = g["truth_counts"][t]
            if n == 0:
                continue
            lines.append(f"| {set_name} | {platform} | {t} | {n} | "
                         f"{pct(safe_div(g['recall'][t]['red'], n))} | "
                         f"{pct(safe_div(g['recall'][t]['orange'], n))} | "
                         f"{pct(safe_div(g['recall'][t]['all'], n))} | "
                         f"{pct(safe_div(g['recall_pass'][t]['all'], n))} | "
                         f"{pct(safe_div(g['fam_recall'][t]['all'], n))} | "
                         f"{pct(safe_div(g['exact_t'][t], n))} |")

    lines += ["", "## annotation regions (not-assessable vs other labels)", ""]
    lines.append("`not_assessable` counts only records whose BED label is "
                 "`not_assessable`; `annotation_bp` is the total over all untiered "
                 "records (incl. `low_coverage`). The `by label (bp)` column gives "
                 "the per-label bp breakdown.")
    lines.append("")
    lines.append("| set | platform | annotation regions | annotation bp | "
                 "by label (bp) | not-assessable regions | not-assessable bp | "
                 "truth in not-assessable |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for (set_name, platform), g in sorted(groups.items()):
        na_regions = g["annotation_counts"].get("not_assessable", 0)
        by_label = ", ".join(f"{lbl}={bp}" for lbl, bp
                             in sorted(g["annotation_bp_by_label"].items()))
        lines.append(f"| {set_name} | {platform} | {g['annotations']} | "
                     f"{g['annotation_bp']} | {by_label or '-'} | {na_regions} | "
                     f"{g['not_assessable_bp']} | {g['truth_in_na']} |")

    lines += ["", "## truth homopolymer buckets (by family)", ""]
    lines.append("| set | platform | family | " + " | ".join(HP_BUCKETS) + " |")
    lines.append("|---" * (3 + len(HP_BUCKETS)) + "|")
    for (set_name, platform), g in sorted(groups.items()):
        for fam in FAMILIES:
            hb = g["truth_hp"][fam]
            if not hb:
                continue
            lines.append(f"| {set_name} | {platform} | {fam} | "
                         + " | ".join(str(hb.get(b, 0)) for b in HP_BUCKETS) + " |")

    lines += ["", "## false positives per tier", ""]
    lines.append("`LowCov` counts the subset of FPs that carry FILTER=LowCov "
                 "(or BED `lowcov=1`); they are part of the tier totals.")
    lines.append("")
    lines.append("| set | platform | RED | ORANGE | YELLOW | total | LowCov | per Mb | "
                 "in unaligned/multi | calls |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for (set_name, platform), g in sorted(groups.items()):
        fp_per_mb = safe_div(g["fp_total"], g["genome_bp"] / 1e6) if g["genome_bp"] else None
        fpm = f"{fp_per_mb:.2f}" if fp_per_mb is not None else "n/a"
        lines.append(f"| {set_name} | {platform} | {g['fp']['RED']} | "
                     f"{g['fp']['ORANGE']} | {g['fp']['YELLOW']} | {g['fp_total']} | "
                     f"{g['fp_lowcov_total']} | {fpm} | {g['in_region']} | {g['calls']} |")

    lines += ["", "## findings likely distribution per set", ""]
    lines.append("Findings group the tiered calls into one problem location each; "
                 "this is the `likely` hint distribution (no truth for kinds here). "
                 "The `(none)`-free total and per-tier counts are also shown.")
    lines.append("")
    lines.append("| set | platform | findings | per Mb | red | orange | yellow | likely |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for (set_name, platform), g in sorted(groups.items()):
        per_mb = safe_div(g["findings_total"], g["genome_bp"] / 1e6) if g["genome_bp"] else None
        per_mb_s = f"{per_mb:.2f}" if per_mb is not None else "n/a"
        ft = g["findings_tier"]
        lk = ", ".join(f"{k} {v}" for k, v in sorted(g["findings_likely"].items()))
        lines.append(f"| {set_name} | {platform} | {g['findings_total']} | {per_mb_s} | "
                     f"{ft.get('RED', 0)} | {ft.get('ORANGE', 0)} | "
                     f"{ft.get('YELLOW', 0)} | {lk or '-'} |")

    invalid = [r for r in results if not r.get("valid", True)]
    if invalid:
        lines += ["", "## invalid cases (excluded from aggregates)", ""]
        for r in invalid:
            lines.append(f"- {r['set']}/{r['case']}: {r.get('error', 'invalid')}")
    lines.append("")
    (out_dir / "Report.md").write_text("\n".join(lines) + "\n")


def write_missed_fp(results, out_dir, cand_cache):
    mlines = ["set\tcase\tcontig\tpos0\ttype\tlen\tref\talt\thp\tstr\t"
              "support_class\tnearest_dist\tnearest_candidate"]
    flines = ["set\tcase\tsource\tcontig\tpos0\tstart\tend\ttype\tlabel\ttier\tlowcov"]
    for r in results:
        for m in r["missed"]:
            cands = cand_cache.get((r["set"], r["case"]), [])
            near = nearest_candidate(cands, m["contig"], m["pos0"]) if cands else None
            dist = near[0] if near else ""
            line = near[1] if near else ""
            mlines.append("\t".join(str(x) for x in [
                r["set"], r["case"], m["contig"], m["pos0"], m["type"], m["len"],
                m.get("ref", ""), m.get("alt", ""), m["hp"], m.get("str", 0),
                m.get("support_class", ""), dist, line]))
        for c in r["fp_list"]:
            pos0 = c.get("pos0", "")
            flines.append("\t".join(str(x) for x in [
                r["set"], r["case"], c["source"], c["contig"], pos0,
                c.get("start", ""), c.get("end", ""), c.get("type", ""),
                c.get("label", ""), c.get("tier", ""),
                1 if c.get("lowcov") else 0]))
    (out_dir / "missed.tsv").write_text("\n".join(mlines) + "\n")
    (out_dir / "fp.tsv").write_text("\n".join(flines) + "\n")


# ---------------------------------------------------------------------------
# self-test
# ---------------------------------------------------------------------------
def _unit_support(work):
    """Unit checks for the by-support classification (incl. missing candidates)."""
    problems = []
    regions = [("chr1", 100000, 101000, "unaligned")]

    def ev(contig, pos0, typ="SNP", hp=1, st=0, start0=None, end0=None):
        return {"contig": contig, "pos0": pos0, "type": typ, "len": 1,
                "ref": "A", "alt": "C", "hp": hp, "str": st,
                "start0": pos0 if start0 is None else start0,
                "end0": pos0 + 1 if end0 is None else end0}

    # repeat_unaligned takes precedence over a high-AF candidate
    got = classify_truth([ev("chr1", 100500)], regions, {"chr1": [(100500, 0.9)]})
    if got != ["repeat_unaligned"]:
        problems.append(f"repeat precedence: {got} != ['repeat_unaligned']")
    # region +/-10 bp boundaries (half-open, matching call_near_region)
    if classify_truth([ev("chr1", 101009)], regions, {"chr1": [(101009, 0.9)]}) \
            != ["repeat_unaligned"]:
        problems.append("region end -1 (inside +10) not repeat_unaligned")
    if classify_truth([ev("chr1", 101010)], regions, {"chr1": [(101010, 0.9)]}) \
            != ["clear"]:
        problems.append("region end +10 (half-open) not clear")
    if classify_truth([ev("chr1", 99990)], regions, {"chr1": [(99990, 0.9)]}) \
            != ["repeat_unaligned"]:
        problems.append("region start -10 edge not repeat_unaligned")
    if classify_truth([ev("chr1", 99989)], regions, {"chr1": [(99989, 0.9)]}) \
            != ["clear"]:
        problems.append("just past region start -10 not clear")
    # classification uses the truth interval, not just the anchor: a long
    # deletion whose anchor is left of the region but whose interval overlaps it
    big_del = ev("chr1", 90000, "LDEL", start0=90001, end0=101000)
    got = classify_truth([big_del], regions, {"chr1": [(90000, 0.9)]})
    if got != ["repeat_unaligned"]:
        problems.append(f"interval overlaps region: {got} != ['repeat_unaligned']")
    # ...and a SNP anchored inside the region but with the interval just past the
    # half-open boundary is not repeat_unaligned
    if classify_truth([ev("chr1", 101010)], regions, {"chr1": [(101010, 0.9)]}) \
            != ["clear"]:
        problems.append("anchor-only region test not replaced")
    # missing candidates file -> unknown (must not fail)
    if classify_truth([ev("chr1", 100)], regions, None) != ["unknown"]:
        problems.append("missing candidates did not yield unknown")
    # AF class boundaries (no region nearby; events spaced beyond tolerance)
    events = [ev("c", b) for b in (0, 100, 200, 300, 400)]
    cands = {"c": [(0, 0.4), (100, 0.6), (200, 0.601), (300, 0.399), (400, 0.0)]}
    got = classify_truth(events, [], cands)
    want = ["mixed", "mixed", "clear", "minority", "no_support"]
    if got != want:
        problems.append(f"AF boundaries: {got} != {want}")
    # candidate outside tolerance -> no_support
    got = classify_truth([ev("c", 1000)], [], {"c": [(2000, 0.9)]})
    if got != ["no_support"]:
        problems.append(f"out-of-tolerance candidate: {got} != ['no_support']")
    # indel tolerance uses max(5, hp, str)
    got = classify_truth([ev("c", 1000, "INS", hp=8)], [], {"c": [(1006, 0.9)]})
    if got != ["clear"]:
        problems.append(f"indel hp tolerance: {got} != ['clear']")

    # --- candidates file parsing (pos vs pos0; empty vs header-only) -------
    p1 = work / "cand_pos.tsv"
    p1.write_text("contig\tpos\taf\nchr1\t101\t0.9\n")
    c1 = load_candidates_af(p1)
    if c1 != {"chr1": [(100, 0.9)]}:
        problems.append(f"1-based `pos` column: {c1}")
    p2 = work / "cand_pos0.tsv"
    p2.write_text("contig\tpos0\tpos\taf\nchr1\t100\t101\t0.9\n")
    c2 = load_candidates_af(p2)
    if c2 != {"chr1": [(100, 0.9)]}:
        problems.append(f"`pos0` not preferred over `pos`: {c2}")
    # empty file -> None -> unknown
    p3 = work / "cand_empty.tsv"
    p3.write_text("")
    c3 = load_candidates_af(p3)
    if c3 is not None:
        problems.append(f"empty candidates file -> {c3!r}, expected None")
    if classify_truth([ev("c", 100)], [], c3) != ["unknown"]:
        problems.append("empty candidates file did not yield unknown")
    # valid header-only file -> {} -> no_support
    p4 = work / "cand_header.tsv"
    p4.write_text("contig\tpos0\taf\n")
    c4 = load_candidates_af(p4)
    if c4 != {}:
        problems.append(f"header-only candidates file -> {c4!r}, expected {{}}")
    if classify_truth([ev("c", 100)], [], c4) != ["no_support"]:
        problems.append("header-only candidates file did not yield no_support")
    return problems


def self_test(args):
    work = ROOT / "scratch" / "qc" / "selftest_eval"
    case_dir = work / "truth" / "selfset" / "case1"
    case_dir.mkdir(parents=True, exist_ok=True)

    events = [
        {"contig": "chr1", "pos0": 100, "type": "SNP", "len": 1, "ref": "A", "alt": "C",
         "hp": 1, "str": 0, "start0": 100, "end0": 101},
        {"contig": "chr1", "pos0": 200, "type": "DEL", "len": 2, "ref": "TA", "alt": "T",
         "hp": 1, "str": 0, "start0": 201, "end0": 203},
        {"contig": "chr1", "pos0": 300, "type": "INS", "len": 3, "ref": "G", "alt": "GTTT",
         "hp": 4, "str": 0, "start0": 301, "end0": 301},
        {"contig": "chr2", "pos0": 500, "type": "LDEL", "len": 100, "ref": "A" + "C" * 100,
         "alt": "A", "hp": 2, "str": 0, "start0": 501, "end0": 601},
        {"contig": "chr2", "pos0": 900, "type": "LINS", "len": 300, "ref": "T",
         "alt": "T" + "A" * 300, "hp": 1, "str": 0, "start0": 901, "end0": 901},
    ]
    with open(case_dir / "truth.tsv", "w") as fh:
        fh.write("contig\tpos0\ttype\tlen\tref\talt\thp\tstr\tstart0\tend0\n")
        for e in events:
            fh.write(f"{e['contig']}\t{e['pos0']}\t{e['type']}\t{e['len']}\t"
                     f"{e['ref']}\t{e['alt']}\t{e['hp']}\t{e['str']}\t"
                     f"{e['start0']}\t{e['end0']}\n")
    with open(case_dir / "regions.bed", "w") as fh:
        fh.write("chr1\t100000\t101000\tunaligned\n")
        fh.write("chr1\t200000\t201000\tmulti\n")
    contig_lengths = {"chr1": 1_000_000, "chr2": 1_000_000}
    (case_dir / "summary.json").write_text(json.dumps(
        {"contig_lengths": contig_lengths}) + "\n")

    # Fake calls: every truth event as TIER=RED (positional match), 3 fake FPs
    # (one of them FILTER=LowCov), one tiered BED call inside the unaligned
    # region (tier in col 6), plus a not_assessable and a low_coverage
    # annotation (no valid tier -> annotations, not calls).
    vcf = ["##fileformat=VCFv4.2"]
    for e in events:
        typ = family(e["type"]) if e["type"] != "SNP" else "SNP"
        vcf.append(f"{e['contig']}\t{e['pos0'] + 1}\t.\t{e['ref']}\t{e['alt']}\t.\t"
                   f"PASS\tTIER=RED;TYPE={typ}")
    vcf.append("chr1\t5001\t.\tA\tC\t.\tPASS\tTIER=RED;TYPE=SNP")
    vcf.append("chr1\t6001\t.\tA\tC\t.\tLowCov\tTIER=ORANGE;TYPE=SNP")
    vcf.append("chr2\t7001\t.\tA\tC\t.\tPASS\tTIER=YELLOW;TYPE=SNP")
    (work / "qc.vcf").write_text("\n".join(vcf) + "\n")
    (work / "qc.bed").write_text(
        "chr1\t100500\t100600\tsnp_cluster\t0\tYELLOW\tx\t0.9\t1.0\tsites=3\n"  # tier col6 -> call
        "chr1\t100\t160\tnot_assessable\t0\tnone\tlow_dp\t.\t.\t.\n"          # annotation
        "chr1\t200\t260\tlow_coverage\t0\tnone\tlow_dp\t.\t.\t.\n")           # annotation
    calls_vcf = parse_vcf(work / "qc.vcf")
    calls_bed, annotations = parse_bed(work / "qc.bed")
    cands = {"chr1": [(100, 0.7), (200, 0.5), (300, 0.2)], "chr2": [(900, 0.9)]}
    res, missed, fps = evaluate_case(
        events, [("chr1", 100000, 101000, "unaligned"),
                 ("chr1", 200000, 201000, "multi")],
        calls_vcf, calls_bed, annotations, contig_lengths, cands=cands)

    problems = _unit_support(work)
    # GNU time parsing: m:ss and h:mm:ss forms
    if parse_hms("0:01.23") != 1.23:
        problems.append(f"parse_hms m:ss -> {parse_hms('0:01.23')}")
    if parse_hms("1:02:03") != 3723.0:
        problems.append(f"parse_hms h:mm:ss -> {parse_hms('1:02:03')}")
    w = parse_time_verbose("Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.23\n")[0]
    if w != 1.23:
        problems.append(f"parse_time_verbose m:ss -> {w}")
    w = parse_time_verbose("\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02:03\n")[0]
    if w != 3723.0:
        problems.append(f"parse_time_verbose h:mm:ss -> {w}")

    # tier comes from INFO TIER=..., FILTER is PASS/LowCov
    if {c["tier"] for c in calls_vcf} != {"RED", "ORANGE", "YELLOW"}:
        problems.append(f"VCF tiers from INFO: {sorted({c['tier'] for c in calls_vcf})}")
    lowcov_vcf = [c for c in calls_vcf if c.get("lowcov")]
    if len(lowcov_vcf) != 1 or lowcov_vcf[0]["tier"] != "ORANGE":
        problems.append(f"FILTER=LowCov not parsed (got {lowcov_vcf})")
    if any(c["source"] != "vcf" or "lowcov" not in c for c in calls_vcf):
        problems.append("VCF calls missing lowcov flag")

    if res["recall_overall"]["red"] != 1.0:
        problems.append(f"recall red {res['recall_overall']['red']} != 1.0")
    if res["recall_overall"]["all"] != 1.0:
        problems.append(f"recall all {res['recall_overall']['all']} != 1.0")
    if res["recall_pass_overall"]["red"] != 1.0:
        problems.append(f"recall pass red {res['recall_pass_overall']['red']} != 1.0")
    if res["recall_pass_overall"]["all"] != 1.0:
        problems.append(f"recall pass all {res['recall_pass_overall']['all']} != 1.0")
    if res["family_recall_overall"]["all"] != 1.0:
        problems.append(f"family recall all {res['family_recall_overall']['all']} != 1.0")
    if res["exact_total"] != len(events):
        problems.append(f"exact {res['exact_total']} != {len(events)}")
    if missed:
        problems.append(f"{len(missed)} truth events missed")
    if res["fp_total"] != 3:
        problems.append(f"FP {res['fp_total']} != 3")
    if res["fp_lowcov_total"] != 1:
        problems.append(f"fp_lowcov_total {res['fp_lowcov_total']} != 1")
    if res["fp_lowcov"]["ORANGE"] != 1:
        problems.append(f"fp_lowcov ORANGE {res['fp_lowcov']['ORANGE']} != 1")
    if res["in_region_total"] != 1:
        problems.append(f"in-region {res['in_region_total']} != 1")
    if len(calls_bed) != 1:
        problems.append("tiered BED call (col 6) not parsed as a call")
    if len(annotations) != 2:
        problems.append("BED annotations not separated from calls")
    if res["not_assessable_bp"] != 60:
        problems.append(f"not_assessable_bp {res['not_assessable_bp']} != 60")
    if res["annotation_bp"] != 120:
        problems.append(f"annotation_bp {res['annotation_bp']} != 120")
    if res["annotation_bp_by_label"].get("low_coverage") != 60:
        problems.append(f"low_coverage bp {res['annotation_bp_by_label'].get('low_coverage')} != 60")
    if res["truth_in_not_assessable"] != 1:
        problems.append(f"truth_in_not_assessable {res['truth_in_not_assessable']} != 1")
    if res["truth_in_annotation"].get("low_coverage") != 1:
        problems.append(f"truth_in_annotation[low_coverage] "
                        f"{res['truth_in_annotation'].get('low_coverage')} != 1")
    if res["truth_in_annotation"].get("not_assessable") != 1:
        problems.append(f"truth_in_annotation[not_assessable] "
                        f"{res['truth_in_annotation'].get('not_assessable')} != 1")

    # start0/end0 must be parsed from truth.tsv and drive truth_interval
    loaded, _, _ = load_truth(case_dir)
    if not loaded or truth_interval(loaded[0]) != (100, 101):
        problems.append(f"load_truth start0/end0 not parsed: "
                        f"{truth_interval(loaded[0]) if loaded else None}")
    # right-anchored origin deletion: interval [0, len), not [1, len+1)
    origin = {"contig": "chr1", "pos0": 0, "type": "LDEL", "len": 5,
              "ref": "AAAAAA", "alt": "A", "hp": 1, "str": 0, "start0": 0, "end0": 5}
    if truth_interval(origin) != (0, 5):
        problems.append(f"origin deletion interval {truth_interval(origin)} != (0, 5)")
    # whole-contig deletion: interval must not extend past the contig
    whole = {"contig": "chr1", "pos0": 0, "type": "LDEL", "len": 10,
             "ref": "A" * 10, "alt": "", "hp": 1, "str": 0, "start0": 0, "end0": 10}
    if truth_interval(whole) != (0, 10):
        problems.append(f"whole-contig deletion interval {truth_interval(whole)} != (0, 10)")

    # by-support classes are computed from the run's candidates
    want_classes = ["clear", "mixed", "minority", "no_support", "clear"]
    if res["support_classes"] != want_classes:
        problems.append(f"support_classes {res['support_classes']} != {want_classes}")
    if res["support"]["clear"] != {"n": 2, "red": 2, "red_orange": 2, "any": 2}:
        problems.append(f"support[clear] = {res['support']['clear']}")
    if sum(v["n"] for v in res["support"].values()) != len(events):
        problems.append("support class counts do not sum to truth_total")
    if res["support_vcf"]["clear"]["red"] != 2:
        problems.append(f"support_vcf[clear].red = {res['support_vcf']['clear']['red']}")

    if problems:
        print("SELF-TEST FAILED:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 1
    print(f"self-test OK: recall_red={res['recall_overall']['red']} "
          f"recall_all={res['recall_overall']['all']} "
          f"recall_pass_all={res['recall_pass_overall']['all']} "
          f"family_recall_all={res['family_recall_overall']['all']} "
          f"exact={res['exact_total']}/{len(events)} fp={res['fp_total']} "
          f"fp_lowcov={res['fp_lowcov_total']} "
          f"in_region={res['in_region_total']} calls_vcf={len(calls_vcf)} "
          f"calls_bed={len(calls_bed)} annotations={len(annotations)} "
          f"annot_bp={res['annotation_bp']} na_bp={res['not_assessable_bp']} "
          f"time_parse=ok")
    return 0


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", default=str(ROOT / "scratch" / "qc" / "truth"))
    ap.add_argument("--runs", default=str(ROOT / "scratch" / "qc" / "runs"))
    ap.add_argument("--sets", default="m6large,real,val")
    ap.add_argument("--cases", default=None, help="comma-separated case substrings")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--extra-args", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--tag", default="default")
    ap.add_argument("--binary", default=str(BIN))
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test(args)

    truth_root = Path(args.truth)
    runs_root = Path(args.runs)
    out_dir = Path(args.out) if args.out else ROOT / "scratch" / "qc" / "eval" / args.tag
    out_dir.mkdir(parents=True, exist_ok=True)
    sets = [s.strip() for s in args.sets.split(",") if s.strip()]
    only = [s for s in args.cases.split(",") if s] if args.cases else None

    results = []
    cand_cache = {}
    for set_name in sets:
        set_truth = truth_root / set_name
        if not set_truth.is_dir():
            print(f"no truth for set {set_name} under {set_truth}", file=sys.stderr)
            continue
        manifest = read_manifest(SET_DIRS.get(set_name, set_truth))
        for case_dir in sorted(p for p in set_truth.iterdir() if p.is_dir()):
            case = case_dir.name
            if only and not any(s in case for s in only):
                continue
            if not (case_dir / "truth.tsv").exists():
                continue
            events, regions, summary = load_truth(case_dir)
            contig_lengths = summary.get("contig_lengths") or {}
            if not contig_lengths:
                fai = SET_DIRS.get(set_name, set_truth) / case / "old.fa.fai"
                if fai.exists():
                    for line in fai.read_text().splitlines():
                        f = line.split("\t")
                        if len(f) >= 2:
                            contig_lengths[f[0]] = int(f[1])

            preset = manifest.get(case, {}).get("preset", "")
            prefix = runs_root / args.tag / set_name / case / "qc"
            if args.run:
                meta = run_case(set_name, case, preset, runs_root, args.tag,
                                args.threads, args.extra_args, Path(args.binary))
                if meta.get("rc") != 0:
                    print(f"  run failed {set_name}/{case}: rc={meta['rc']} "
                          f"{meta.get('stderr_tail', '')}", file=sys.stderr)
            else:
                meta = read_run_meta(prefix)

            have_outputs = all(Path(str(prefix) + ext).exists()
                               for ext in (".vcf", ".bed", ".json"))
            rc_ok = meta.get("rc") == 0
            valid = bool(rc_ok and have_outputs)
            calls_vcf, calls_bed, annotations = [], [], []
            cands = None
            findings = []
            error = ""
            # Parse outputs only for a valid run; a parse failure invalidates the
            # case (recorded as a reason) instead of crashing the whole run.
            if valid:
                try:
                    calls_vcf = parse_vcf(str(prefix) + ".vcf")
                    calls_bed, annotations = parse_bed(str(prefix) + ".bed")
                    cand_cache[(set_name, case)] = load_candidates(
                        str(prefix) + ".candidates.tsv")
                    # missing candidates.tsv -> None -> classes become `unknown`
                    cands = load_candidates_af(str(prefix) + ".candidates.tsv")
                    findings = parse_findings(str(prefix) + ".findings.bed")
                except Exception as exc:
                    valid = False
                    error = f"qc output parse error: {exc!r}"
            if not valid and not error:
                if not meta:
                    error = "no run metadata (run with --run)"
                elif meta.get("rc") != 0:
                    error = f"qc rc={meta['rc']}: {meta.get('stderr_tail', '')}"
                else:
                    error = "missing qc.vcf/qc.bed/qc.json"

            res, missed, fps = evaluate_case(events, regions, calls_vcf, calls_bed,
                                             annotations, contig_lengths, cands=cands)
            flikely = defaultdict(int)
            ftier = {"RED": 0, "ORANGE": 0, "YELLOW": 0}
            # A wrapped finding is two BED rows with the same id: count findings,
            # not rows.
            seen = set()
            for f in findings:
                if f["label"] in seen:
                    continue
                seen.add(f["label"])
                ftier[f["tier"]] = ftier.get(f["tier"], 0) + 1
                flikely[f["likely"]] += 1
            res.update({"set": set_name, "case": case,
                        "platform": manifest.get(case, {}).get("platform", ""),
                        "preset": preset, "valid": valid,
                        "wall_s": meta.get("wall_s"), "max_rss_kb": meta.get("max_rss_kb"),
                        "missed": missed, "fp_list": fps,
                        "findings_total": len(seen), "findings_tier": ftier,
                        "findings_likely": dict(flikely)})
            if error:
                res["error"] = " ".join(error.split())
            results.append(res)
            print(f"{'ok ' if valid else 'INVALID'} {set_name}/{case}: "
                  f"truth={res['truth_total']} "
                  f"recall_all={pct(res['recall_overall']['all'])} "
                  f"fp={res['fp_total']} in_region={res['in_region_total']} "
                  f"calls={res['calls']['vcf'] + res['calls']['bed']}")

    groups = summarise_results(results)
    write_per_case(results, out_dir)
    write_summary_tsv(groups, out_dir)
    write_by_support(groups, out_dir)
    write_missed_fp(results, out_dir, cand_cache)
    write_report(results, groups, out_dir, args)
    print(f"\nwrote {out_dir}/"
          f"{{per_case.tsv,summary.tsv,by_support.tsv,Report.md,missed.tsv,fp.tsv}}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
