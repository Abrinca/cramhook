#!/usr/bin/env python3
"""One-command real-data regression smoke check (M9).

Runs a small, fixed set of real cases through ``scripts/real_eval.py`` and
prints a compact pass/fail table.  The cases are chosen to be fast while
covering the risky lift paths.  Datasets are local-only, so cases whose
directory/files are absent are reported as SKIP and skipped.

Usage:
  python scripts/real_smoke.py [--threads N] [--binary PATH] [--no-lock]

Each root is evaluated with a single ``real_eval.py --bench 0`` call (no
benchmarking), tag ``smoke``, output under ``scratch/real-smoke/<root-name>/``.
By default every call is wrapped in ``flock scratch/m6.lock`` so it serialises
against other heavy repo tasks; ``--no-lock`` disables that.

Exit codes: 0 all gates pass, 1 a real_eval call failed or a gate failed,
2 every case was missing.
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REAL_EVAL = ROOT / "scripts" / "real_eval.py"
TAG = "smoke"
BENCH = 0
T_MISMATCH_GATE = 0.001  # must match real_eval.py T_MISMATCH_GATE
REQUIRED = ("old.fa", "new.fa", "old.cram", "oracle.cram", "oracle.json")
LOCK = ROOT / "scratch" / "m6.lock"

# (root, case, one-line reason).  Fast, but each covers a risky lift path.
CASES = [
    ("data/v2/real", "ont_FAM29080__cons-to-polished",
     "ONT use case; many plasmids / origin-crossing reads"),
    ("data/v2/real", "ont_FAM20862__polished-to-rotated",
     "rotation + reverse-complemented chromosome, z-drop rule"),
    ("data/v2/real", "pacbio_FAM20848-p1-1.1__cons-to-finalized",
     "PacBio use case"),
    ("data/m6large/real", "ont_FAM22000__minipolish_04-to-polished",
     "split/junction-heavy raw->polished"),
    ("data/real", "pacbio_FAM14197-p1-1.1__flye_02-to-cons",
     "both contigs reversed; PacBio raw->cons"),
]


def root_key(root):
    """Filesystem-safe, collision-free name for a root path."""
    return str(root).replace("/", "_").replace("\\", "_")


def missing_reason(case_dir):
    """None if the case is usable, else a short SKIP reason."""
    if not case_dir.is_dir():
        return "case directory not present"
    missing = [f for f in REQUIRED if not (case_dir / f).exists()]
    if missing:
        return "missing " + ", ".join(missing)
    return None


def load_results(out_dir):
    """case -> result dict, from summary.json (preferred) then *.result.json."""
    tag_dir = out_dir / TAG
    results = {}
    summary = tag_dir / "summary.json"
    if summary.exists():
        try:
            data = json.loads(summary.read_text())
            for r in data.get("results", []):
                if isinstance(r, dict) and r.get("case"):
                    results[r["case"]] = r
        except (OSError, ValueError):
            pass
    for p in sorted(tag_dir.glob("*.result.json")):
        case = p.name[: -len(".result.json")]
        if case in results:
            continue
        try:
            results[case] = json.loads(p.read_text())
        except (OSError, ValueError):
            pass
    return results


def case_passes(r):
    """Mirror real_eval.py's per-case correctness gate."""
    if not r or r.get("error"):
        return False
    if not r.get("invariants_ok"):
        return False
    if (r.get("t_mismatch_rate") or 0.0) > T_MISMATCH_GATE:
        return False
    if r.get("t_supp_fraction", 1.0) < 1.0:
        return False
    return True


def fmt_pct(v, nd=2):
    return f"{v * 100:.{nd}f}" if isinstance(v, (int, float)) else "n/a"


def fmt_num(v, nd=2):
    return f"{v:.{nd}f}" if isinstance(v, (int, float)) else "n/a"


def make_command(root, cases, args, out_dir):
    cmd = [sys.executable, str(REAL_EVAL), "--root", root,
           "--only", ",".join(cases), "--tag", TAG,
           "--out", str(out_dir), "--bench", str(BENCH)]
    if args.threads is not None:
        cmd += ["--threads", str(args.threads)]
    if args.binary:
        cmd += ["--binary", args.binary]
    if args.lock:
        LOCK.parent.mkdir(parents=True, exist_ok=True)
        LOCK.touch(exist_ok=True)
        cmd = ["flock", str(LOCK), *cmd]
    return cmd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threads", type=int, default=None,
                    help="pass --threads to real_eval (default: real_eval's own)")
    ap.add_argument("--binary", default=None,
                    help="pass --binary PATH to real_eval "
                         "(default: target/release/cramhook)")
    ap.add_argument("--no-lock", dest="lock", action="store_false",
                    help="do not wrap each real_eval call in flock")
    ap.set_defaults(lock=True)
    args = ap.parse_args()

    # Group cases by root, keeping the declared order.
    roots = []
    for root, case, _reason in CASES:
        if root not in roots:
            roots.append(root)

    skipped = {}         # case -> reason (dataset dir/files not present locally)
    present_cases = []   # cases whose dataset exists and were handed to real_eval
    case_rows = {}       # case -> result dict (only present cases with a result)
    any_call_failed = False

    t_start = time.monotonic()
    for root in roots:
        group = [(c, reason) for (r, c, reason) in CASES if r == root]
        present = []
        for case, reason in group:
            why = missing_reason(ROOT / root / case)
            if why is None:
                present.append(case)
                present_cases.append(case)
            else:
                print(f"SKIP {case}: {why}")
                skipped[case] = why
        if not present:
            continue

        out_dir = ROOT / "scratch" / "real-smoke" / root_key(root)
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = make_command(root, present, args, out_dir)
        print(f"=== {root} ({len(present)} case{'s' if len(present) != 1 else ''}) ===")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            any_call_failed = True
            if proc.stdout:
                print(proc.stdout.rstrip())
            if proc.stderr:
                print(proc.stderr.rstrip(), file=sys.stderr)
            print(f"  real_eval rc={proc.returncode}")
        results = load_results(out_dir)
        for case in present:
            r = results.get(case)
            if isinstance(r, dict):
                case_rows[case] = r
            else:
                print(f"  no readable result for {case}", file=sys.stderr)

    total_wall = time.monotonic() - t_start

    # Compact table, one row per declared case (SKIP -> n/a).
    header = ["case", "translated %", "T mism", "T supp %",
              "depth ident %", "max ddiff", "lift wall s", "result"]
    rows = []
    for _root, case, _reason in CASES:
        if case in skipped:
            rows.append([case, *(["n/a"] * 6), "SKIP"])
            continue
        r = case_rows.get(case)
        if r is None:
            # Present dataset but no result: not a skip, it must fail.
            rows.append([case, *(["n/a"] * 6), "FAIL"])
            continue
        if r.get("error"):
            rows.append([case, *(["n/a"] * 6), "FAIL"])
            continue
        t = r.get("t") or {}
        dep = r.get("depth") or {}
        rows.append([
            case,
            fmt_pct(r.get("translated_fraction")),
            str(t.get("mismatch", "n/a")),
            fmt_pct(r.get("t_supp_fraction"), 4),
            fmt_pct(dep.get("fraction_identical"), 4),
            str(dep.get("max_abs_diff", "n/a")),
            fmt_num(r.get("lift_wall_s")),
            "PASS" if case_passes(r) else "FAIL",
        ])

    widths = [len(h) for h in header]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))
    line = "  ".join(h.ljust(widths[i]) for i, h in enumerate(header))
    sep = "  ".join("-" * widths[i] for i in range(len(header)))
    print()
    print(line)
    print(sep)
    for row in rows:
        print("  ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)))

    print()
    if not present_cases:
        print("all cases missing (local-only datasets); nothing to run")
        return 2

    no_result = [c for c in present_cases if c not in case_rows]
    gate_failures = [c for c, r in case_rows.items() if not case_passes(r)]
    print(f"total wall time: {total_wall:.1f}s")
    if any_call_failed or no_result or gate_failures:
        if no_result:
            print(f"no result: {', '.join(sorted(no_result))}", file=sys.stderr)
        if gate_failures:
            print(f"gate failures: {', '.join(sorted(gate_failures))}", file=sys.stderr)
        if any_call_failed:
            print("one or more real_eval calls failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
