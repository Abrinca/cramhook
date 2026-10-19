#!/usr/bin/env python3
"""Local release check on the public FAM29080 test data (Plan-v4 step 3).

1. ``cargo build --release`` + ``cargo test --release`` (unit + synthetic suites).
2. On the test data (``pixi run fetch-testdata`` or ``build-testdata`` first):

   * ``ont_lift`` (case A, ONT): lift cons.cram -> polished, qc on the lifted CRAM,
     fish with findings + haplotypes; plus a full remap onto polished.fa.
   * ``ont_bait`` (case B, ONT): fish --bait phage.fa on the lifted CRAM (polished
     lacks the phage); truth = a full remap onto polished.fa + phage.
   * ``hifi_bait`` (case B, HiFi): qc + fish --bait on asm.cram (asm.fa lacks the
     phage); truth = a full remap onto asm.fa + phage.

3. Invariants (nonzero exit if one fails): the lift equals the full remap by the
   ``scripts/compare.py`` criteria, qc reports a red/orange finding covering
   >= 80 % of the prophage on both platforms, fish --bait recalls the truth phage
   reads, and every output kept in ``testdata/expected/`` exists.
4. Diff report against ``testdata/expected/`` (never fails).  ``--accept``
   overwrites ``testdata/expected/`` with the current outputs (only when all
   invariants pass).

Outputs go to ``testdata/data/run/`` (git-ignored).  Commands run with cwd
``testdata/data`` and relative paths, so paths inside the outputs are stable.

Usage::

  pixi run release-check [--accept] [--no-cargo-test] [--threads N]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import DATA, EXPECTED, ROOT, data_dir, is_unpacked, load_manifest  # noqa: E402
from compare import compare  # noqa: E402  (scripts/compare.py)

BIN = ROOT / "target" / "release" / "cramhook"
RUN = DATA / "run"

PHAGE = "FAM29080-p1-1.1_scf7"
# Prophage (integrated copy of the phage) in each assembly, 0-based half-open;
# the same intervals as scripts/fish_phage_eval.py.
PROPHAGE = {
    "ont": ("1", 543393, 563743),
    "hifi": ("FAM29080-p1-1.1_scf1", 273076, 293426),
}
PROPHAGE_TIERS = ("red", "orange")
MIN_PROPHAGE_COVER = 0.8
# fish --bait recall of the truth phage reads (any-record truth).  Measured on
# v1 of the data: ONT 22/22, HiFi 35/35 (100 %).  Earlier full-depth ONT runs
# saw ~70 %, so ONT keeps a wider margin.
MIN_BAIT_RECALL = {"ont": 0.80, "hifi": 0.95}

# Files kept in testdata/expected/<case>/ (small text; JSON is normalised).
KEEP = {
    "ont_lift": ["lift.stats.json", "qc.vcf", "qc.bed", "qc.findings.bed",
                 "qc.json", "qc.haplotypes.tsv", "fish.fish.tsv",
                 "fish.fish.json", "summary.json"],
    "ont_bait": ["fish.fish.tsv", "fish.fish.json", "summary.json"],
    "hifi_bait": ["qc.vcf", "qc.bed", "qc.findings.bed", "qc.json",
                  "qc.haplotypes.tsv", "fish.fish.tsv", "fish.fish.json",
                  "summary.json"],
}
# JSON keys dropped before storing/comparing (timings, machine-dependent).
VOLATILE_KEYS = {"timing_s", "timings_s", "threads"}
TIERS = ["red", "orange", "yellow", "gray"]


# ---------------------------------------------------------------------------
# running
# ---------------------------------------------------------------------------
class Runner:
    def __init__(self, log_path):
        self.log = open(log_path, "w")
        self.steps = []

    def run(self, label, cmd, shell=False, cwd=DATA):
        t0 = time.perf_counter()
        text = cmd if shell else " ".join(str(c) for c in cmd)
        self.log.write(f"$ (cd {cwd}) {text}\n")
        self.log.flush()
        p = subprocess.run(
            ["bash", "-o", "pipefail", "-c", cmd] if shell else [str(c) for c in cmd],
            cwd=cwd, stdout=self.log, stderr=subprocess.STDOUT)
        dt = time.perf_counter() - t0
        self.log.write(f"# rc={p.returncode} {dt:.1f}s\n\n")
        self.log.flush()
        self.steps.append((label, dt, p.returncode))
        print(f"  {label:<34} {dt:6.1f}s{'' if p.returncode == 0 else f'  FAILED rc={p.returncode}'}",
              flush=True)
        if p.returncode != 0:
            raise StepFailed(f"{label} failed (rc={p.returncode}); see {rel(self.log.name)}")
        return dt


class StepFailed(Exception):
    pass


def rel(p):
    try:
        return str(Path(p).resolve().relative_to(ROOT))
    except ValueError:
        return str(p)


def write_fasta_concat(out, *fastas):
    """``out`` = the records of ``fastas`` (re-wrapped; inputs may lack a final newline)."""
    with open(out, "w") as fh:
        for fa in fastas:
            with pysam.FastxFile(str(fa)) as fx:
                for r in fx:
                    fh.write(f">{r.name}{(' ' + r.comment) if r.comment else ''}\n")
                    s = r.sequence
                    for i in range(0, len(s), 80):
                        fh.write(s[i:i + 80] + "\n")


def remap(runner, label, preset, ref, reads, out, threads):
    """The input-contract mapping command (a full remap)."""
    runner.run(label, f"minimap2 -ax {preset} -t {threads} --secondary=no -Y {ref} {reads} "
                      f"| samtools sort -@ 4 -O CRAM --reference {ref} -o {out} "
                      f"&& samtools index {out}", shell=True)


# ---------------------------------------------------------------------------
# evaluation helpers
# ---------------------------------------------------------------------------
def read_findings(path):
    out = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        out.append({"contig": f[0], "start": int(f[1]), "end": int(f[2]),
                    "label": f[3], "tier": f[5], "likely": f[6] if len(f) > 6 else ""})
    return out


def prophage_cover(f, plat):
    """Fraction of the prophage interval covered by finding ``f`` (0 on another contig)."""
    c, s, e = PROPHAGE[plat]
    if f["contig"] != c:
        return 0.0
    return max(0, min(e, f["end"]) - max(s, f["start"])) / (e - s)


def prophage_finding(findings_path, plat):
    """The finding that best covers the prophage (red/orange first, then the
    largest cover), with its cover fraction; None if none overlaps."""
    hits = [dict(f, prophage_cover=round(prophage_cover(f, plat), 4))
            for f in read_findings(findings_path)]
    hits = [f for f in hits if f["prophage_cover"] > 0]
    hits.sort(key=lambda f: (f["tier"] not in PROPHAGE_TIERS, -f["prophage_cover"]))
    return hits[0] if hits else None


def prophage_ok(f, plat):
    """qc must report the prophage as one red/orange finding spanning >= 80 % of it
    (the whole element is one finding: junction, depth steps and clips are linked)."""
    return (f is not None and f["tier"] in PROPHAGE_TIERS
            and prophage_cover(f, plat) >= MIN_PROPHAGE_COVER)


def truth_reads(cram, ref):
    """Read names with a counted record on the phage contig (primary, any)."""
    prim, anyr = set(), set()
    with pysam.AlignmentFile(str(cram), "rc", reference_filename=str(ref)) as af:
        for r in af.fetch(PHAGE):
            if r.is_unmapped or r.is_secondary or r.is_duplicate or r.is_qcfail:
                continue
            anyr.add(r.query_name)
            if not r.is_supplementary:
                prim.add(r.query_name)
    return prim, anyr


def fish_rows(tsv):
    rows = []
    lines = Path(tsv).read_text().splitlines()
    hdr = lines[0].split("\t")
    for line in lines[1:]:
        if line and not line.startswith("#"):
            rows.append(dict(zip(hdr, line.split("\t"))))
    return rows


def bait_summary(case_dir, plat, truth_cram, truth_ref):
    prim, anyr = truth_reads(truth_cram, truth_ref)
    rows = fish_rows(case_dir / "fish.fish.tsv")
    sel = {r["read"] for r in rows}
    bait = {r["read"] for r in rows if r["reason"] == "bait"}
    fj = json.loads((case_dir / "fish.fish.json").read_text())
    return {
        "truth_primary": len(prim), "truth_any": len(anyr),
        "selected_reads": len(sel), "bait_reads": len(bait),
        "truth_selected": len(sel & anyr),
        "recall_any": round(len(sel & anyr) / len(anyr), 4) if anyr else None,
        "recall_primary": round(len(sel & prim) / len(prim), 4) if prim else None,
        "bait_not_truth": len(bait - anyr),
        "bait_lift": {k: fj.get("bait", {}).get(k)
                      for k in ("translated_reads", "remapped_reads", "fallback")},
    }


def lift_summary(case_dir, ref):
    res = compare(case_dir / "lift.cram", ref, case_dir / "oracle.cram", ref)
    pr = res["primary"]

    def cat(c):
        d = pr.get(c) or {}
        return {k: d.get(k, 0) for k in ("n", "n_identical", "n_equivalent", "n_mismatch")}
    so = res["supplementary"]["overall"]
    dp = res["depth"]["overall"]
    st = json.loads((case_dir / "lift.stats.json").read_text())
    return {
        "translated_reads": st["translated_reads"],
        "remapped_reads": st["remapped_reads"],
        "invariants_ok": res["invariants_ok"],
        "mismatch_total": res["mismatch_total"],
        "primary_T": cat("T"), "primary_R": cat("R"),
        "supplementary": {k: so[k] for k in ("n", "equal", "equivalent", "differ")},
        "depth_max_abs_diff": dp["max_abs_diff"],
        "depth_fraction_identical": dp["fraction_identical"],
    }


# ---------------------------------------------------------------------------
# the data checks
# ---------------------------------------------------------------------------
def run_data(runner, threads):
    D = data_dir().name           # inputs, relative to DATA
    if RUN.exists():
        shutil.rmtree(RUN)
    for c in KEEP:
        (RUN / c).mkdir(parents=True)
    b = BIN
    t = str(threads)
    ont, hifi, phage = f"{D}/ont", f"{D}/hifi", f"{D}/phage.fa"
    summ = {}

    # case A: ONT cons -> polished
    a = "run/ont_lift"
    runner.run("ont_lift: lift", [b, "lift", "--in", f"{ont}/cons.cram", "--old-ref", f"{ont}/cons.fa",
                                  "--new-ref", f"{ont}/polished.fa", "--out", f"{a}/lift.cram",
                                  "--threads", t, "--stats", f"{a}/lift.stats.json",
                                  "--min-input-mb", "0"])
    runner.run("ont_lift: qc", [b, "qc", "--in", f"{a}/lift.cram", "--ref", f"{ont}/polished.fa",
                                "--out-prefix", f"{a}/qc", "--threads", t])
    runner.run("ont_lift: fish", [b, "fish", "--in", f"{a}/lift.cram", "--ref", f"{ont}/polished.fa",
                                  "--out-prefix", f"{a}/fish", "--threads", t,
                                  "--findings", f"{a}/qc.findings.bed",
                                  "--haplotypes", f"{a}/qc.haplotypes.tsv"])
    remap(runner, "ont_lift: full remap (oracle)", "lr:hq", f"{ont}/polished.fa",
          f"{ont}/reads.fastq.gz", f"{a}/oracle.cram", threads)
    t0 = time.perf_counter()
    summ["ont_lift"] = lift_summary(DATA / a, DATA / ont / "polished.fa")
    summ["ont_lift"]["prophage_finding"] = prophage_finding(DATA / a / "qc.findings.bed", "ont")
    print(f"  {'ont_lift: compare.py':<34} {time.perf_counter() - t0:6.1f}s", flush=True)

    # case B, ONT: bait on the lifted (polished, phage-less) CRAM
    o = "run/ont_bait"
    runner.run("ont_bait: fish --bait", [b, "fish", "--in", f"{a}/lift.cram", "--ref", f"{ont}/polished.fa",
                                         "--out-prefix", f"{o}/fish", "--threads", t,
                                         "--findings", f"{a}/qc.findings.bed",
                                         "--haplotypes", f"{a}/qc.haplotypes.tsv", "--bait", phage])
    write_fasta_concat(DATA / o / "truth_ref.fa", DATA / ont / "polished.fa", DATA / phage)
    remap(runner, "ont_bait: full remap (truth)", "lr:hq", f"{o}/truth_ref.fa",
          f"{ont}/reads.fastq.gz", f"{o}/truth.cram", threads)
    summ["ont_bait"] = bait_summary(DATA / o, "ont", DATA / o / "truth.cram", DATA / o / "truth_ref.fa")

    # case B, HiFi
    h = "run/hifi_bait"
    runner.run("hifi_bait: qc", [b, "qc", "--in", f"{hifi}/asm.cram", "--ref", f"{hifi}/asm.fa",
                                 "--out-prefix", f"{h}/qc", "--threads", t])
    runner.run("hifi_bait: fish --bait", [b, "fish", "--in", f"{hifi}/asm.cram", "--ref", f"{hifi}/asm.fa",
                                          "--out-prefix", f"{h}/fish", "--threads", t,
                                          "--findings", f"{h}/qc.findings.bed",
                                          "--haplotypes", f"{h}/qc.haplotypes.tsv", "--bait", phage])
    write_fasta_concat(DATA / h / "truth_ref.fa", DATA / hifi / "asm.fa", DATA / phage)
    remap(runner, "hifi_bait: full remap (truth)", "map-hifi", f"{h}/truth_ref.fa",
          f"{hifi}/reads.fastq.gz", f"{h}/truth.cram", threads)
    summ["hifi_bait"] = bait_summary(DATA / h, "hifi", DATA / h / "truth.cram", DATA / h / "truth_ref.fa")
    summ["hifi_bait"]["prophage_finding"] = prophage_finding(DATA / h / "qc.findings.bed", "hifi")

    for case, s in summ.items():
        (RUN / case / "summary.json").write_text(json.dumps(s, indent=2) + "\n")
    return summ


def invariants(summ):
    """``[(name, ok, measured, required)]``."""
    out = []
    lf = summ["ont_lift"]
    T, R = lf["primary_T"], lf["primary_R"]
    out.append(("lift: compare.py invariants (names, uniqueness, no secondary/orphans)",
                lf["invariants_ok"], lf["invariants_ok"], True))
    out.append(("lift: translated reads (non-vacuous)", T["n"] > 0, T["n"], "> 0"))
    out.append(("lift: primary mismatches vs full remap", lf["mismatch_total"] == 0,
                f"{lf['mismatch_total']} (T {T['n_identical']}+{T['n_equivalent']}/{T['n']}, "
                f"R {R['n_identical']}+{R['n_equivalent']}/{R['n']} identical+equivalent)", "0"))
    out.append(("lift: supplementary differ vs full remap", lf["supplementary"]["differ"] == 0,
                lf["supplementary"]["differ"], "0"))
    out.append(("lift: depth max |diff| vs full remap", lf["depth_max_abs_diff"] <= 1,
                lf["depth_max_abs_diff"], "<= 1"))
    for case, plat in (("ont_lift", "ont"), ("hifi_bait", "hifi")):
        f = summ[case]["prophage_finding"]
        c, s, e = PROPHAGE[plat]
        ok = prophage_ok(f, plat)
        got = (f"{f['label']} {f['tier']} {f['likely']} {f['contig']}:{f['start']}-{f['end']}, "
               f"covers {f['prophage_cover']:.1%}" if f else "none")
        out.append((f"qc {plat}: finding over prophage {c}:{s}-{e}", ok, got,
                    f"{'/'.join(PROPHAGE_TIERS)}, covers >= {MIN_PROPHAGE_COVER:.0%}"))
    for case, plat in (("ont_bait", "ont"), ("hifi_bait", "hifi")):
        b = summ[case]
        r = b["recall_any"]
        out.append((f"fish --bait {plat}: truth recall", r is not None and r >= MIN_BAIT_RECALL[plat],
                    f"{r} ({b['truth_selected']}/{b['truth_any']})", f">= {MIN_BAIT_RECALL[plat]}"))
    return out


# ---------------------------------------------------------------------------
# normalised outputs, diff report, accept
# ---------------------------------------------------------------------------
def strip_volatile(x):
    if isinstance(x, dict):
        return {k: strip_volatile(v) for k, v in x.items() if k not in VOLATILE_KEYS}
    if isinstance(x, list):
        return [strip_volatile(v) for v in x]
    return x


def normalised(case, name):
    p = RUN / case / name
    if not p.exists():
        return None
    text = p.read_text()
    if name.endswith(".json"):
        text = json.dumps(strip_volatile(json.loads(text)), indent=2) + "\n"
    return text


def vcf_sites(text):
    out = {}
    for line in (text or "").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        out[(f[0], int(f[1]), f[3], f[4])] = (info.get("TIER", "?").lower(), f[6])
    return out


def bed_tiers(text):
    c = Counter()
    for line in (text or "").splitlines():
        if line and not line.startswith("#"):
            f = line.split("\t")
            c[f"{f[3]}:{f[5]}"] += 1
    return c


def findings_of(text):
    out = []
    for line in (text or "").splitlines():
        if line and not line.startswith("#"):
            f = line.split("\t")
            out.append((f[0], int(f[1]), int(f[2]), f[5], f[6] if len(f) > 6 else ""))
    return out


def fmt_counter(old, new):
    keys = sorted(set(old) | set(new))
    return ", ".join(f"{k} {old.get(k, 0)}->{new.get(k, 0)}" for k in keys
                     if old.get(k, 0) != new.get(k, 0))


def diff_json(a, b, path=""):
    """Leaf differences between two JSON values: ``[(path, old, new)]``."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b), key=str):
            out += diff_json(a.get(k), b.get(k), f"{path}.{k}" if path else str(k))
        return out
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += diff_json(x, y, f"{path}[{i}]")
        return out
    return [] if a == b else [(path, a, b)]


def diff_report():
    print("\n== diff vs testdata/expected (informational) ==")
    if not EXPECTED.exists():
        print("  no testdata/expected/ yet (run with --accept to create it)")
        return
    n_same = n_total = 0
    for case, names in KEEP.items():
        lines = []
        for name in names:
            n_total += 1
            new = normalised(case, name)
            ep = EXPECTED / case / name
            old = ep.read_text() if ep.exists() else None
            if old == new:
                n_same += 1
                continue
            if old is None:
                lines.append(f"{name}: new file (not in expected)")
                continue
            if new is None:
                lines.append(f"{name}: missing in this run")
                continue
            lines += [f"{name}: {d}" for d in describe(name, old, new)]
        if lines:
            print(f"  [{case}]")
            for line in lines:
                print(f"    {line}")
    if n_same == n_total:
        print(f"  no changes ({n_total} files identical)")
    else:
        print(f"  {n_same}/{n_total} files identical")


def describe(name, old, new):
    out = []
    if name.endswith(".vcf"):
        a, b = vcf_sites(old), vcf_sites(new)
        ca = Counter(f"{t}/{flt}" for t, flt in a.values())
        cb = Counter(f"{t}/{flt}" for t, flt in b.values())
        if ca != cb:
            out.append(f"sites per tier/FILTER: {fmt_counter(ca, cb)}")
        added = [k for k in b if k not in a]
        lost = [k for k in a if k not in b]
        changed = [k for k in a if k in b and a[k] != b[k]]
        if added or lost or changed:
            out.append(f"sites added {len(added)}, lost {len(lost)}, tier/FILTER changed {len(changed)}"
                       + "".join(f"\n      {k[0]}:{k[1]} {a[k][0]}->{b[k][0]}" for k in changed[:5]))
        if not out:
            out.append("text changed (same sites and tiers)")
    elif name.endswith("findings.bed"):
        a, b = findings_of(old), findings_of(new)

        def match(f, pool):
            return next((g for g in pool if g[0] == f[0] and g[1] < f[2] and f[1] < g[2]), None)
        for f in b:
            g = match(f, a)
            if g is None:
                out.append(f"finding added {f[0]}:{f[1]}-{f[2]} {f[3]} {f[4]}")
            elif (g[3], g[4]) != (f[3], f[4]):
                out.append(f"finding {f[0]}:{f[1]}-{f[2]} {g[3]}/{g[4]} -> {f[3]}/{f[4]}")
            elif g[1:3] != f[1:3]:
                out.append(f"finding {f[0]}:{g[1]}-{g[2]} -> {f[1]}-{f[2]} ({f[3]})")
        for g in a:
            if match(g, b) is None:
                out.append(f"finding lost {g[0]}:{g[1]}-{g[2]} {g[3]} {g[4]}")
        if not out:
            out.append("text changed (same findings)")
    elif name.endswith(".bed"):
        ca, cb = bed_tiers(old), bed_tiers(new)
        out.append(f"regions per label:tier: {fmt_counter(ca, cb)}" if ca != cb
                   else "regions changed (same counts per label:tier)")
    elif name == "fish.fish.tsv":
        def per_reason(text):
            sets = {}
            for line in text.splitlines()[1:]:
                f = line.split("\t")
                sets.setdefault(f[2], set()).add(f[0])
            return sets
        a, b = per_reason(old), per_reason(new)
        for r in sorted(set(a) | set(b)):
            x, y = a.get(r, set()), b.get(r, set())
            if x != y:
                out.append(f"reason {r}: {len(x)} -> {len(y)} reads (+{len(y - x)} -{len(x - y)})")
        if not out:
            out.append("rows changed (same reads per reason)")
    elif name.endswith(".json"):
        d = diff_json(json.loads(old), json.loads(new))
        if name == "qc.json":
            # keep the report compact: alert counts and the headline numbers first
            d.sort(key=lambda x: (not x[0].startswith(("alerts.", "findings", "qc_qv")), x[0]))
        for p, a, b in d[:12]:
            out.append(f"{p}: {a} -> {b}")
        if len(d) > 12:
            out.append(f"... {len(d) - 12} more changed keys")
    elif name.endswith(".tsv"):
        a, b = set(old.splitlines()), set(new.splitlines())
        out.append(f"rows +{len(b - a)} -{len(a - b)}")
    return out


def missing_outputs():
    return [f"{case}/{name}" for case, names in KEEP.items() for name in names
            if not (RUN / case / name).exists()]


def accept():
    # read everything first: the old baseline is only replaced by a complete one
    texts = {(case, name): normalised(case, name)
             for case, names in KEEP.items() for name in names}
    missing = [f"{c}/{n}" for (c, n), t in texts.items() if t is None]
    if missing:
        raise SystemExit(f"not accepting, outputs missing: {', '.join(missing)}")
    if EXPECTED.exists():
        shutil.rmtree(EXPECTED)
    for (case, name), text in texts.items():
        (EXPECTED / case).mkdir(parents=True, exist_ok=True)
        (EXPECTED / case / name).write_text(text)
    print(f"\naccepted: wrote {rel(EXPECTED)}/ from this run")


def cargo_test_totals(log_path, pos):
    """Sum the ``test result:`` lines cargo printed since log offset ``pos``."""
    tot = Counter()
    with open(log_path) as fh:
        fh.seek(pos)
        for line in fh:
            if line.startswith("test result:"):
                for part in line.split(";"):
                    w = part.split()
                    for i, tok in enumerate(w[:-1]):
                        if tok.isdigit() and w[i + 1] in ("passed", "failed", "ignored"):
                            tot[w[i + 1]] += int(tok)
    return f"{tot['passed']} passed, {tot['failed']} failed, {tot['ignored']} ignored"


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accept", action="store_true",
                    help="overwrite testdata/expected/ with this run's outputs")
    ap.add_argument("--no-cargo-test", action="store_true", help="skip cargo test")
    ap.add_argument("--threads", type=int, default=min(16, os.cpu_count() or 1))
    args = ap.parse_args(argv)
    t_start = time.perf_counter()

    man = load_manifest()
    if not is_unpacked(man):
        sys.exit(f"test data {man['name']} (sha256 {man['sha256'][:12]}...) not in "
                 f"{rel(DATA)}/: run `pixi run fetch-testdata` (or `pixi run build-testdata`)")
    DATA.mkdir(parents=True, exist_ok=True)
    RUN.mkdir(parents=True, exist_ok=True)
    runner = Runner(DATA / "release-check.log")
    failures = []

    print("== build and tests ==")
    try:
        runner.run("cargo build --release", ["cargo", "build", "--release"], cwd=ROOT)
        if not args.no_cargo_test:
            pos = runner.log.tell()
            try:
                runner.run("cargo test --release", ["cargo", "test", "--release"], cwd=ROOT)
            finally:
                print("    " + cargo_test_totals(Path(runner.log.name), pos))
    except StepFailed as e:
        failures.append(str(e))
    if failures and not BIN.exists():
        print("\n".join(f"FAIL: {f}" for f in failures))
        return 1

    print(f"== test data ({man['name']}, {args.threads} threads) ==")
    summ = None
    try:
        summ = run_data(runner, args.threads)
    except StepFailed as e:
        failures.append(str(e))

    if summ is not None:
        missing = missing_outputs()
        if missing:
            failures.append(f"outputs missing: {', '.join(missing)}")
        print("\n== invariants ==")
        for name, ok, got, need in invariants(summ):
            print(f"  {'ok  ' if ok else 'FAIL'} {name}: {got} (required {need})")
            if not ok:
                failures.append(name)
        diff_report()

    print(f"\nrelease-check: {time.perf_counter() - t_start:.0f}s total; log {rel(runner.log.name)}")
    if failures:
        print(f"FAILED ({len(failures)}):")
        for f in failures:
            print(f"  - {f}")
        if args.accept:
            print("not accepting: fix the failures first")
        return 1
    if args.accept:
        accept()
    print("release-check: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
