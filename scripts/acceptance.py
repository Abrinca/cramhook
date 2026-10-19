#!/usr/bin/env python3
"""Acceptance checks for `cramhook lift` on the synthetic scenarios (docs/plans/Plan-v1.md §7).

Usage:
  python scripts/acceptance.py [--only S1,S3] [--threads 16] [--out data/synthetic/acceptance]

For each data/synthetic/<S>/<model> it runs `cramhook lift`, compares the result
to the oracle with scripts/compare.py, and evaluates the `docs/plans/Plan-v1.md` section 7
criteria C1-C5 (C6 is a reported metric, not a gate).  It writes a markdown
summary and JSON into the output directory and exits non-zero if any gate fails.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare as cmp_mod  # noqa: E402
from compare import compare, write_result  # noqa: E402

SYNTH = ROOT / "data" / "synthetic"
# Allow a non-default binary (e.g. a scratch CARGO_TARGET_DIR build) to be used.
BIN = Path(os.environ.get("CRAMHOOK_BIN", ROOT / "target" / "release" / "cramhook"))
PRESETS = {"hifi": "map-hifi", "ont": "lr:hq"}
WINDOW_MARGIN = 500


# ----------------------------------------------------------------------------
# process helpers
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


def discover(only=None):
    cases = []
    if not SYNTH.is_dir():
        return cases
    for s in sorted(p for p in SYNTH.iterdir() if p.is_dir()):
        if only and s.name not in only:
            continue
        truth = s / "truth.json"
        old_ref, new_ref = s / "old.fa", s / "new.fa"
        if not truth.exists() or not old_ref.exists() or not new_ref.exists():
            continue
        for model in sorted(p for p in s.iterdir() if p.is_dir()):
            if model.name not in PRESETS:
                continue
            if not (model / "old.cram").exists() or not (model / "oracle.cram").exists():
                continue
            cases.append({"scenario": s.name, "model": model.name,
                          "dir": model, "truth": truth,
                          "old": old_ref, "new": new_ref})
    return cases


# ----------------------------------------------------------------------------
# C5 helpers: identical translated records
# ----------------------------------------------------------------------------
def _tag_value(r, tag):
    v = r.get_tag(tag)
    if isinstance(v, (list, tuple)):
        return tuple(v)
    return v


def record_sig(r, drop=("XL", "MD", "s2")):
    """Hashable signature of a record ignoring XL, regenerated MD and ``s2``.

    M10b: translated (``XL:T``) records drop ``s2`` (minimap2's competing chain
    score on the old assembly, which cannot be recomputed by translation), so
    ``s2`` is excluded from both sides of the C5 "identical to old" comparison.
    """
    try:
        tags = tuple(sorted((t, _tag_value(r, t)) for t, _ in r.get_tags()
                            if t not in drop))
    except Exception:
        tags = ()
    quals = r.query_qualities
    return (r.flag, r.reference_name, r.reference_start, r.mapping_quality,
            r.cigarstring, r.query_sequence,
            tuple(quals) if quals is not None else None, tags)


def translated_identical_to_old(lift, old, ref):
    """Check every XL:T record in ``lift`` has a matching record in ``old``.

    Returns (ok, n_primary_checked, n_all_checked, n_missing, missing_names).
    """
    def load(path):
        af = pysam.AlignmentFile(str(path), "rc", reference_filename=str(ref))
        out = {}
        try:
            for r in af.fetch(until_eof=True):
                sig = record_sig(r)
                out.setdefault(r.query_name, Counter())[sig] += 1
        finally:
            af.close()
        return out

    old_by_name = load(old)
    missing = []
    n_all = n_primary = 0
    af = pysam.AlignmentFile(str(lift), "rc", reference_filename=str(ref))
    try:
        for r in af.fetch(until_eof=True):
            if not (r.has_tag("XL") and str(r.get_tag("XL")) == "T"):
                continue
            n_all += 1
            if not r.flag & 0x100 and not r.flag & 0x800:
                n_primary += 1
            sig = record_sig(r)
            avail = old_by_name.get(r.query_name)
            if not avail or avail.get(sig, 0) <= 0:
                missing.append(r.query_name)
            else:
                avail[sig] -= 1
    finally:
        af.close()
    return (not missing), n_primary, n_all, len(missing), sorted(set(missing))


def xl_audit(path, ref):
    """Audit XL tags: every record must be T or R.

    Returns a dict with totals, T/R counts split primary/supplementary, and the
    number of records with a missing or invalid tag.
    """
    out = {"n_total": 0, "n_missing": 0, "n_bad": 0,
           "primary": Counter(), "supplementary": Counter(), "bad_values": Counter()}
    af = pysam.AlignmentFile(str(path), "rc", reference_filename=str(ref))
    try:
        for r in af.fetch(until_eof=True):
            out["n_total"] += 1
            if not r.has_tag("XL"):
                out["n_missing"] += 1
                continue
            v = str(r.get_tag("XL"))
            if v not in ("T", "R"):
                out["n_bad"] += 1
                out["bad_values"][v] += 1
                continue
            if r.flag & 0x100:
                out.setdefault("secondary", Counter())[v] += 1
            elif r.flag & 0x800:
                out["supplementary"][v] += 1
            else:
                out["primary"][v] += 1
    finally:
        af.close()
    out["complete"] = out["n_missing"] == 0 and out["n_bad"] == 0
    out["bad_values"] = dict(out["bad_values"])
    out["primary"] = dict(out["primary"])
    out["supplementary"] = dict(out["supplementary"])
    out["secondary"] = dict(out.get("secondary", {}))
    return out


def all_records_xl(path, ref, value):
    """True if every record carries XL:<value>; returns (ok, n, n_bad)."""
    n = bad = 0
    af = pysam.AlignmentFile(str(path), "rc", reference_filename=str(ref))
    try:
        for r in af.fetch(until_eof=True):
            n += 1
            if not (r.has_tag("XL") and str(r.get_tag("XL")) == value):
                bad += 1
    finally:
        af.close()
    return bad == 0, n, bad


def supp_differences(test, test_ref, oracle, ref, cat=None):
    """Read names whose supplementary sets differ beyond equivalence.

    With ``cat`` set, restrict to reads whose test XL category equals ``cat``.
    """
    tp, _tn, ts, _tsec, txl, tdet = cmp_mod.collect(test, test_ref)
    op, _on, os_, _osec, _oxl, odet = cmp_mod.collect(oracle, ref)
    diffs = []
    for name in set(tp) | set(op):
        if cat is not None and txl.get(name, "none") != cat:
            continue
        a = ts.get(name, set())
        b = os_.get(name, set())
        if a == b:
            continue
        status = cmp_mod._supp_status(
            a, b,
            (tdet.get(name) or {}).get("supp", []),
            (odet.get(name) or {}).get("supp", []))
        # Only true "differ" cases count: "equal" means the sets agree once
        # low-MAPQ supplementary records are ignored (M7b), and "equivalent"
        # means they agree up to AS/CIGAR placement.
        if status == "differ":
            diffs.append((name, sorted(a), sorted(b)))
    return diffs


# ----------------------------------------------------------------------------
# C6 helpers: unnecessary diff-level deferrals
# ----------------------------------------------------------------------------
def edit_intervals(truth):
    """Old-coordinate intervals of planted edits, grouped by old contig."""
    contigs = truth.get("contigs_old", {})
    out = {}
    for e in truth.get("edits", []):
        t = e.get("type")
        oc = e.get("old_contig")
        if oc is None:
            continue  # add_contig has no old location
        if t in ("remove_contig", "reverse_complement"):
            L = contigs.get(oc)
            if L is not None:
                out.setdefault(oc, []).append((0, L))
            continue
        if t == "rotation":
            L = contigs.get(oc)
            if L is not None:
                out.setdefault(oc, []).append((0, 0))
                out.setdefault(oc, []).append((L, L))
            off = e.get("offset")
            if off is not None:
                out.setdefault(oc, []).append((off, off))
            continue
        s = e.get("old_start")
        if s is None:
            continue
        en = e.get("old_end")
        if en is None or en <= s:
            en = s + 1
        out.setdefault(oc, []).append((s, en))
    return out


def _overlap(a, b):
    return a[0] < b[1] and b[0] < a[1]


# Edit types that carry no old-assembly location by design.
OLD_LESS_TYPES = {"reorder", "add_contig"}


def edit_recognised(e):
    """True when `edit_intervals` can locate this planted edit in old coords.

    Guards the truth schema: a point edit that uses `contig`/`pos` instead of
    `old_contig`/`old_start` is silently skipped by `edit_intervals` and must be
    reported, not ignored.
    """
    t = e.get("type")
    if t in OLD_LESS_TYPES:
        return True
    if e.get("old_contig") is None:
        return False
    if t in ("rename", "remove_contig", "reverse_complement", "rotation"):
        return True
    return e.get("old_start") is not None


def unrecognised_edits(truth):
    return [i for i, e in enumerate(truth.get("edits", [])) if not edit_recognised(e)]


def justified_t2(taint_bed, repeats_old):
    """T2 taint intervals that overlap a repeats_old copy (old coords)."""
    out = []
    taint_bed = Path(taint_bed)
    if not taint_bed.exists():
        return out
    for line in taint_bed.read_text().splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        c, s, e, reason = parts[0], int(parts[1]), int(parts[2]), parts[3]
        if reason != "T2":
            continue
        for rep in repeats_old:
            if rep["contig"] == c and _overlap((s, e), (rep["start"], rep["end"])):
                out.append((c, s, e))
                break
    return out


def read_old_windows(old_cram, old_ref):
    """name -> (contig, start, end, left_clip, right_clip) for old primaries."""
    out = {}
    af = pysam.AlignmentFile(str(old_cram), "rc", reference_filename=str(old_ref))
    try:
        for r in af.fetch(until_eof=True):
            if r.is_unmapped or r.flag & 0x100 or r.flag & 0x800:
                continue
            tup = r.cigartuples or []
            lc = tup[0][1] if tup and tup[0][0] in (4, 5) else 0
            rc = tup[-1][1] if tup and tup[-1][0] in (4, 5) else 0
            out[r.query_name] = (r.reference_name, r.reference_start,
                                 r.reference_end, lc, rc)
    finally:
        af.close()
    return out


def c6_metric(deferred_tsv, old_cram, old_ref, taint_bed, truth):
    """(n_unnecessary, n_diff, fraction) for diff-class deferrals."""
    edits = edit_intervals(truth)
    justified = justified_t2(taint_bed, truth.get("repeats_old", []))
    windows = read_old_windows(old_cram, old_ref)
    n_diff = n_unnecessary = 0
    examples = []
    if not Path(deferred_tsv).exists():
        return n_unnecessary, n_diff, None, examples
    for line in Path(deferred_tsv).read_text().splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        name, reason, cls = parts[0], parts[1], parts[2]
        if cls != "diff":
            continue
        n_diff += 1
        w = windows.get(name)
        if w is None:
            n_unnecessary += 1
            examples.append(name)
            continue
        contig, start, end, lc, rc = w
        ws, we = start - lc - WINDOW_MARGIN, end + rc + WINDOW_MARGIN
        unnecessary = True
        for iv in edits.get(contig, []):
            if _overlap((ws, we), iv):
                unnecessary = False
                break
        if unnecessary:
            for (c, ts, te) in justified:
                if c == contig and _overlap((ws, we), (ts, te)):
                    unnecessary = False
                    break
        if unnecessary:
            n_unnecessary += 1
            if len(examples) < 20:
                examples.append(name)
    frac = (n_unnecessary / n_diff) if n_diff else None
    return n_unnecessary, n_diff, frac, examples


# ----------------------------------------------------------------------------
# per-case evaluation
# ----------------------------------------------------------------------------
def _canonical_lift(model_dir, old_ref, new_ref, threads):
    return ["lift", "--in", str(model_dir / "old.cram"),
            "--old-ref", str(old_ref),
            "--new-ref", str(new_ref),
            "--out", str(model_dir / "lift.cram"),
            "--threads", str(threads),
            "--stats", str(model_dir / "lift.stats.json"),
            "--xl-tag",
            # synthetic inputs are tiny; test the translate path, not the M9 fallback
            "--min-input-mb", "0"]


def evaluate_case(case, out_dir, threads):
    scen, model = case["scenario"], case["model"]
    model_dir = case["dir"]
    old_ref, new_ref = case["old"], case["new"]
    truth = read_json(case["truth"]) or {}
    result = {"scenario": scen, "model": model, "c1": False, "c2": False,
              "c2b": False, "c3": False, "c4": None, "c5": None, "c6": None,
              "error": None}
    report = [f"## {scen}/{model}", ""]

    # 0. truth schema: every planted edit must be locatable in old coordinates.
    unrec = unrecognised_edits(truth)
    result["edits_planted"] = len(truth.get("edits", []))
    result["edits_unrecognised"] = unrec
    if unrec:
        report.append(f"- EDIT TRUTH: {len(unrec)} of {result['edits_planted']} "
                      f"planted edits are not recognised (indices {unrec})")
        print(f"  WARNING {scen}/{model}: {len(unrec)} planted edit(s) "
              f"unrecognised: {unrec}", file=sys.stderr)

    # 1. run lift (always rerun)
    lift_cmd = [str(BIN)] + _canonical_lift(model_dir, old_ref, new_ref, threads)
    t0 = time.monotonic()
    p = run(lift_cmd)
    wall = time.monotonic() - t0
    result["lift_wall_s"] = wall
    result["lift_rc"] = p.returncode
    if p.returncode != 0:
        result["error"] = f"lift failed rc={p.returncode}: {p.stderr.strip()[-500:]}"
        report.append(f"lift FAILED (rc={p.returncode})")
        report.append("```\n" + p.stderr.strip()[-2000:] + "\n```")
        result["report"] = "\n".join(report)
        return result
    lift = model_dir / "lift.cram"
    stats = read_json(model_dir / "lift.stats.json") or {}

    # 2. compare lift vs oracle
    res = compare(lift, new_ref, model_dir / "oracle.cram", new_ref)
    write_result(res, model_dir / "cmp_lift")

    flat = res["primary"]["overall"]
    n = flat["n"] or 0

    # XL audit: every record must be T or R; reconcile with lift stats.
    rec_in = (stats.get("records_in") or {}) if isinstance(stats, dict) else {}
    prim_in = rec_in.get("primary") if isinstance(rec_in, dict) else None
    translated = stats.get("translated_reads") if isinstance(stats, dict) else None
    deferred_reads = stats.get("deferred_reads") if isinstance(stats, dict) else None
    xl = xl_audit(lift, new_ref)
    result["xl"] = xl
    n_prim_t = xl["primary"].get("T", 0)
    n_prim_r = xl["primary"].get("R", 0)
    result["xl_primary_T"] = n_prim_t
    result["xl_primary_R"] = n_prim_r
    accounting_ok = True
    accounting_msg = ""
    if translated is None or deferred_reads is None:
        accounting_ok = False
        accounting_msg = "lift stats missing translated_reads/deferred_reads"
    elif n_prim_t != translated or n_prim_r != deferred_reads:
        accounting_ok = False
        accounting_msg = (f"XL=T primaries {n_prim_t} != translated_reads {translated} "
                          f"or XL=R primaries {n_prim_r} != deferred_reads {deferred_reads}")
    result["accounting_ok"] = accounting_ok
    result["accounting_msg"] = accounting_msg

    # C1: invariants + exit code + XL audit/accounting integrity
    result["invariants_ok"] = bool(res["invariants_ok"])
    result["c1"] = (bool(res["invariants_ok"]) and p.returncode == 0
                    and xl["complete"] and accounting_ok)

    # C2: XL=T agreement, non-vacuous, and no worse than baseline noise.
    tcat = res["primary"].get("T", {"n": 0, "n_identical": 0,
                                    "n_equivalent": 0, "n_mismatch": 0})
    tn = tcat["n"]
    agree = (tcat["n_identical"] + tcat["n_equivalent"]) / tn if tn else 1.0
    noise = read_json(model_dir / "cmp_noise.json")
    noise_n = ((noise or {}).get("primary", {}).get("overall", {}) or {}).get("n")
    noise_mm = (noise or {}).get("mismatch_total")
    noise_rate = (noise_mm / noise_n) if (noise and noise_n) else None
    result["noise_mismatch"] = noise_mm
    result["noise_rate"] = noise_rate
    result["c2_noise_ok"] = noise_rate is not None
    mm_rate = (tcat["n_mismatch"] / tn) if tn else 1.0
    result["c2_mismatch_rate"] = mm_rate
    result["c2"] = (tn > 0 and agree >= 0.999
                    and noise_rate is not None
                    and mm_rate <= noise_rate)  # Plan 7.2, literal: not worse than noise
    result["c2_fraction"] = agree
    result["c2_n"] = tn
    result["c2_mismatch"] = tcat["n_mismatch"]

    # C2b: supplementary agreement for translated (XL=T) reads must be 100 %.
    stcat = res["supplementary"].get("T", {"n": 0, "equal": 0,
                                           "equivalent": 0, "differ": 0})
    sn = stcat["n"]
    sdiff = stcat["differ"]
    result["c2b_n"] = sn
    result["c2b_equal"] = stcat["equal"]
    result["c2b_equivalent"] = stcat["equivalent"]
    result["c2b_differ"] = sdiff
    result["c2b_fraction"] = ((stcat["equal"] + stcat["equivalent"]) / sn) if sn else 1.0
    result["c2b"] = (sdiff == 0 and sn > 0)
    result["c2b_differing"] = []
    if sdiff > 0:
        result["c2b_differing"] = [
            x[0] for x in supp_differences(lift, new_ref, model_dir / "oracle.cram",
                                           new_ref, cat="T")][:50]

    # C3: depth
    dp = res["depth"]["overall"]
    result["c3"] = dp["fraction_identical"] >= 0.999 and dp["max_abs_diff"] <= 2
    result["depth_fraction"] = dp["fraction_identical"]
    result["depth_max_abs_diff"] = dp["max_abs_diff"]

    # R-category agreement (reported)
    rcat = res["primary"].get("R", {"n": 0, "n_identical": 0, "n_equivalent": 0})
    result["r_n"] = rcat["n"]
    result["r_identical"] = rcat["n_identical"]
    result["r_equivalent"] = rcat["n_equivalent"]

    # translated fraction and deferrals (defensive)
    result["translated_reads"] = translated
    result["deferred_reads"] = deferred_reads
    result["translated_fraction"] = (translated / prim_in) if (translated is not None
                                                               and prim_in) else None
    result["deferred_by_reason"] = stats.get("deferred_by_reason")
    result["deferred_by_class"] = stats.get("deferred_by_class")
    oracle_meta = read_json(model_dir / "oracle.json") or {}
    result["oracle_wall_s"] = oracle_meta.get("wall_s")
    if result["oracle_wall_s"]:
        result["wall_ratio"] = wall / result["oracle_wall_s"]

    # C5 (S1 only): diff, taint, diff-level deferrals, identical translation
    if scen == "S1":
        with tempfile.TemporaryDirectory(dir=str(out_dir)) as td:
            pfx = Path(td) / "diff"
            d = run([str(BIN), "diff", "--old-ref", str(old_ref),
                     "--new-ref", str(new_ref),
                     "--preset", PRESETS[model], "--out-prefix", str(pfx),
                     "--threads", str(threads)])
            c5 = {"diff_rc": d.returncode}
            if d.returncode == 0:
                c5["blocks_ok"], c5["blocks"] = _blocks_ok(
                    Path(str(pfx) + ".blocks.bed"), truth, new_ref)
                c5["taint"] = _taint_lines(Path(str(pfx) + ".taint.bed"))
                c5["diff_deferrals"] = (stats.get("deferred_by_class") or {}).get("diff", 0) \
                    if isinstance(stats, dict) else 0
                ok, prim_checked, all_checked, missing, names = \
                    translated_identical_to_old(lift, model_dir / "old.cram", old_ref)
                c5["translated_checked"] = prim_checked
                c5["translated_records_checked"] = all_checked
                c5["translated_missing"] = missing
                c5["translated_missing_names"] = names[:50]
                c5["identical_ok"] = ok
                # non-vacuous: primary translated reads must be accounted for
                c5["checked_matches_stats"] = (
                    translated is not None and prim_checked == translated
                    and prim_checked > 0)
                c5["ok"] = (c5["blocks_ok"] and not c5["taint"]
                            and c5["diff_deferrals"] == 0 and ok
                            and c5["checked_matches_stats"])
            else:
                c5["ok"] = False
            result["c5"] = c5["ok"]
            result["c5_detail"] = c5
        report.append(f"- C5: {'PASS' if result['c5'] else 'FAIL'} {result.get('c5_detail', {})}")

    # C4 (S1/hifi, S1/ont): contract violation -> full remap fallback
    if scen == "S1" and model in ("hifi", "ont"):
        result["c4"], result["c4_detail"] = _c4_case(case, out_dir, threads)
        report.append(f"- C4: {'PASS' if result['c4'] else 'FAIL'} "
                      f"{result.get('c4_detail', {})}")

    # C6 (reported, not a gate)
    with tempfile.TemporaryDirectory(dir=str(out_dir)) as td:
        td = Path(td)
        # diff for taint regions
        pfx = td / "diff"
        run([str(BIN), "diff", "--old-ref", str(old_ref),
             "--new-ref", str(new_ref), "--preset", PRESETS[model],
             "--out-prefix", str(pfx), "--threads", str(threads)])
        # translate for deferred.tsv
        tpfx = td / "tr"
        run([str(BIN), "translate", "--in", str(model_dir / "old.cram"),
             "--old-ref", str(old_ref),
             "--new-ref", str(new_ref),
             "--out-prefix", str(tpfx), "--threads", str(threads), "--xl-tag"])
        nu, nd, frac, ex = c6_metric(
            Path(str(tpfx) + ".deferred.tsv"), model_dir / "old.cram",
            old_ref, Path(str(pfx) + ".taint.bed"), truth)
        result["c6_unnecessary"] = nu
        result["c6_diff_deferrals"] = nd
        result["c6_fraction"] = frac
        result["c6_examples"] = ex
        report.append(f"- C6 (reported): {nu}/{nd} unnecessary"
                      + (f" ({frac:.3%})" if frac is not None else ""))

    # human-readable per-case report
    report.append("")
    report.append(f"- C1 invariants_ok={result['invariants_ok']} lift_rc={p.returncode} "
                  f"xl_complete={xl['complete']} xl_T={n_prim_t} xl_R={n_prim_r} "
                  f"accounting_ok={accounting_ok}"
                  + (f" ({accounting_msg})" if accounting_msg else ""))
    report.append(f"- C2 XL=T agreement={agree:.4%} (n={tn}, equivalent={tcat['n_equivalent']}, "
                  f"mismatch={tcat['n_mismatch']}); mismatch_rate={mm_rate:.4%} "
                  f"vs max(0.1%, noise={('%.4f%%' % (100 * noise_rate)) if noise_rate is not None else 'n/a'})")
    report.append(f"- C2b XL=T supplementary agreement={result['c2b_fraction']:.4%} "
                  f"(n={sn}, equal={result['c2b_equal']}, equivalent={result['c2b_equivalent']}, "
                  f"differ={sdiff})")
    if result["c2b_differing"]:
        report.append(f"  differing translated reads ({sdiff} total, first "
                      f"{len(result['c2b_differing'])} listed):")
        for nm in result["c2b_differing"]:
            report.append(f"  - {nm}")
    report.append(f"- C3 depth identical={dp['fraction_identical']:.4%} "
                  f"max_abs_diff={dp['max_abs_diff']}")
    report.append(f"- R agreement identical={result['r_identical']}/{result['r_n']} "
                  f"(+{result['r_equivalent']} equivalent)")
    if result["translated_fraction"] is not None:
        report.append(f"- translated={result['translated_reads']} "
                      f"({result['translated_fraction']:.3%})")
    else:
        report.append(f"- translated={result['translated_reads']}")
    report.append(f"- deferrals by class={result['deferred_by_class']} "
                  f"by reason={result['deferred_by_reason']}")
    report.append(f"- wall {wall:.3f}s vs oracle {result['oracle_wall_s']}s "
                  f"(ratio {result.get('wall_ratio')})")
    # mismatch listing for C2
    mism = res.get("mismatches_full", [])
    if mism:
        report.append("")
        report.append(f"### C2 mismatches ({len(mism)} total, first "
                      f"{min(len(mism), 50)} listed)")
        for m in mism[:50]:
            report.append(f"- {m['name']} [{m['category']}] test={m['test']} "
                          f"oracle={m['oracle']}")
    result["report"] = "\n".join(report)
    return result


def _blocks_ok(blocks_bed, truth, new_ref_path):
    """Exactly one identity block per old contig spanning it, mapped to itself."""
    contigs = truth.get("contigs_old", {})
    if not Path(blocks_bed).exists():
        return False, {"error": "no blocks.bed"}
    rows = []
    for line in Path(blocks_bed).read_text().splitlines():
        parts = line.split()
        if len(parts) < 6:
            continue
        rows.append(parts)
    seen = set()
    for r in rows:
        oc, s, e, nc, ns, strand = r[0], int(r[1]), int(r[2]), r[3], int(r[4]), r[5]
        L = contigs.get(oc)
        if L is None or oc != nc or s != 0 or e != L or ns != 0 or strand != "+":
            return False, {"bad_row": r}
        if oc in seen:
            return False, {"duplicate": oc}
        seen.add(oc)
    ok = seen == set(contigs)
    return ok, {"n_blocks": len(rows), "contigs": sorted(seen)}


def _taint_lines(taint_bed):
    if not Path(taint_bed).exists():
        return []
    return [ln for ln in Path(taint_bed).read_text().splitlines() if ln.split()[:1]]


def _c4_case(case, out_dir, threads):
    """Build a --secondary=yes CRAM, lift it, require full remap fallback."""
    model_dir = case["dir"]
    old_ref, new_ref = case["old"], case["new"]
    scen, model = case["scenario"], case["model"]
    d = out_dir / scen / model
    d.mkdir(parents=True, exist_ok=True)
    sec = d / "c4_secondary.cram"
    preset = PRESETS[model]
    cmd = (f"minimap2 -ax {preset} -t {threads} --secondary=yes -Y "
           f"{old_ref} {model_dir / 'reads.fastq.gz'} | "
           f"samtools sort -O CRAM --reference {old_ref} -o {sec}")
    p = run(cmd)
    detail = {"cram_rc": p.returncode}
    if p.returncode != 0:
        detail["error"] = p.stderr.strip()[-500:]
        return False, detail
    run(["samtools", "index", str(sec)])
    liftout = d / "c4_lift.cram"
    q = run([str(BIN), "lift", "--in", str(sec),
             "--old-ref", str(old_ref), "--new-ref", str(new_ref),
             "--out", str(liftout), "--threads", str(threads),
             "--stats", str(d / "c4_lift.stats.json"), "--xl-tag",
             "--min-input-mb", "0"])
    detail["lift_rc"] = q.returncode
    if q.returncode != 0:
        detail["error"] = q.stderr.strip()[-500:]
        return False, detail
    all_r, n, n_bad = all_records_xl(liftout, new_ref, "R")
    detail["all_xl_R"] = all_r
    detail["n_records"] = n
    detail["n_not_R"] = n_bad
    res = compare(liftout, new_ref, model_dir / "oracle.cram", new_ref)
    write_result(res, d / "cmp_c4")
    flat = res["primary"]["overall"]
    frac = (flat["n_identical"] / flat["n"] if flat["n"] else 1.0)
    detail["primary_identical_fraction"] = frac
    detail["primary_mismatch"] = flat["n_mismatch"]
    detail["invariants_ok"] = bool(res["invariants_ok"])
    mism = res.get("mismatches_full", [])
    detail["mismatches_listed"] = len(mism)
    # supplementary agreement (equal or equivalent)
    so = res["supplementary"]["overall"]
    sn = so["n"]
    supp_frac = ((so["equal"] + so["equivalent"]) / sn) if sn else 1.0
    detail["supplementary_n"] = sn
    detail["supplementary_equal"] = so["equal"]
    detail["supplementary_equivalent"] = so["equivalent"]
    detail["supplementary_differ"] = so["differ"]
    detail["supplementary_fraction"] = supp_frac
    supp_diffs = []
    if supp_frac < 1.0:
        supp_diffs = supp_differences(liftout, new_ref, model_dir / "oracle.cram", new_ref)
        detail["supplementary_difference_names"] = [x[0] for x in supp_diffs[:50]]
    detail["supplementary_differences_listed"] = len(supp_diffs)
    # Plan 7.4: output equals a full remap (equivalent alignments allowed for supplementaries)
    ok = (all_r and bool(res["invariants_ok"]) and frac == 1.0
          and supp_frac == 1.0)
    return ok, detail


# ----------------------------------------------------------------------------
# summary
# ----------------------------------------------------------------------------
def _mark(v):
    if v is None:
        return "—"
    return "PASS" if v else "FAIL"


def write_summary(results, out_dir):
    header = ("| scenario/model | C1 | C2 | C2b | C3 | C4 | C5 | C6 (reported) |\n"
              "|---|---|---|---|---|---|---|---|")
    lines = ["# cramhook acceptance", "", header]
    for r in results:
        c6 = "—"
        if r.get("c6_diff_deferrals") is not None:
            frac = r.get("c6_fraction")
            c6 = (f"{r['c6_unnecessary']}/{r['c6_diff_deferrals']}"
                  + (f" ({frac:.3%})" if frac is not None else ""))
        lines.append(f"| {r['scenario']}/{r['model']} | {_mark(r['c1'])} | "
                     f"{_mark(r['c2'])} | {_mark(r.get('c2b'))} | {_mark(r['c3'])} | "
                     f"{_mark(r['c4'])} | {_mark(r['c5'])} | {c6} |")
    lines.append("")
    for r in results:
        if r.get("report"):
            lines.append(r["report"])
            lines.append("")
    md = "\n".join(lines)
    (out_dir / "acceptance.md").write_text(md + "\n")

    gates = [r for r in results if r.get("c1") is False or r.get("c2") is False
             or r.get("c2b") is False or r.get("c3") is False or r.get("c4") is False
             or r.get("c5") is False]
    summary = {
        "n_cases": len(results),
        "n_gate_failures": len(gates),
        "results": results,
    }
    (out_dir / "acceptance.json").write_text(json.dumps(summary, indent=2) + "\n")
    return md, summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", default=None,
                    help="comma-separated scenario names, e.g. S1,S3")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--out", default=str(SYNTH / "acceptance"))
    args = ap.parse_args()

    if not BIN.exists():
        print(f"error: {BIN} not found; build it with `pixi run cargo build --release`",
              file=sys.stderr)
        return 2
    only = set(args.only.split(",")) if args.only else None
    cases = discover(only)
    if not cases:
        print("no synthetic cases found", file=sys.stderr)
        return 2
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for case in cases:
        print(f"=== {case['scenario']}/{case['model']} ===")
        try:
            results.append(evaluate_case(case, out_dir, args.threads))
        except Exception as exc:
            results.append({"scenario": case["scenario"], "model": case["model"],
                            "c1": False, "c2": False, "c3": False, "c4": None,
                            "c5": None, "c6": None, "error": repr(exc),
                            "report": None})
    md, summary = write_summary(results, out_dir)
    print()
    print(md.split("\n\n")[0])
    print(f"\nwrote {out_dir / 'acceptance.md'} and {out_dir / 'acceptance.json'}")
    n_unrec = sum(1 for r in results if r.get("edits_unrecognised"))
    if n_unrec:
        print(f"{n_unrec} case(s) have planted edits that are not recognised",
              file=sys.stderr)
        return 1
    if summary["n_gate_failures"]:
        print(f"{summary['n_gate_failures']} case(s) failed a gate", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
