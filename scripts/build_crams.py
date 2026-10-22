#!/usr/bin/env python3
"""Build old.cram / oracle.cram / noise.cram for synthetic and real cases (docs/project/plans/Plan-v1.md §6.1/6.2).

Usage:
  python scripts/build_crams.py {synthetic|real} [--threads 16] [--force] [--no-noise] [--only SUBSTR]

Standard mapping command shape (run with bash -o pipefail):
  minimap2 -ax $PRESET -t N --secondary=no -Y REF reads.fastq.gz \
    | samtools sort -@ 4 -O CRAM --reference REF -o OUT.cram
  samtools index OUT.cram
"""
import argparse
import json
import os
import shutil
import socket
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import compare, write_result  # noqa: E402

PRESETS = {"hifi": "map-hifi", "ont": "lr:hq"}
SORT_THREADS = 4
NOISE_CONTIG = "dummy_noise"
NOISE_LEN = 50_000
NOISE_SEED = 12345


def minimap2_version():
    return subprocess.run(["minimap2", "--version"], capture_output=True, text=True,
                          check=True).stdout.strip()


def samtools_version():
    out = subprocess.run(["samtools", "--version"], capture_output=True, text=True,
                         check=True).stdout
    return out.splitlines()[0].strip() if out else ""


def ensure_ref(fa):
    """Ensure a .fai exists next to the given FASTA path.

    ``samtools faidx`` writes ``<path>.fai`` next to the *given* path, so indexing a
    symlink keeps the .fai inside the case dir even when the target lives in a read-only
    source tree.  Never index the resolved target directly.
    """
    fa = Path(fa)
    if Path(str(fa) + ".fai").exists():
        return fa
    if os.access(fa.parent, os.W_OK):
        subprocess.run(["samtools", "faidx", str(fa)], check=True)
        return fa
    print(f"    cannot write .fai next to {fa}; copying into case dir")
    tmp = fa.with_name(fa.name + ".copy")
    shutil.copyfile(fa.resolve(), tmp)
    tmp.replace(fa)  # replaces the symlink entry with a real file
    subprocess.run(["samtools", "faidx", str(fa)], check=True)
    return fa


def make_noise_ref(new_ref, out):
    """new.fa + a 50 kb random dummy contig (written atomically)."""
    rng = np.random.default_rng(NOISE_SEED)
    bases = "".join(np.array(list("ACGT"))[rng.integers(0, 4, NOISE_LEN)])
    text = Path(new_ref).read_text()
    if not text.endswith("\n"):
        text += "\n"
    tmp = Path(str(out) + ".tmp")
    with open(tmp, "w") as fh:
        fh.write(text)
        fh.write(f">{NOISE_CONTIG}\n")
        for i in range(0, len(bases), 60):
            fh.write(bases[i:i + 60] + "\n")
    tmp.replace(out)


def run_mapping(name, ref, reads, preset, outdir, threads, force):
    cram = outdir / f"{name}.cram"
    crai = outdir / f"{name}.cram.crai"
    meta = outdir / f"{name}.json"
    if not force and cram.exists() and crai.exists() and meta.exists():
        print(f"  {name}: exists, skip")
        return
    meta.unlink(missing_ok=True)  # completion marker; written last
    tmp_cram = outdir / f"{name}.tmp.cram"
    tmp_crai = outdir / f"{name}.tmp.cram.crai"
    cmd = (f"minimap2 -ax {preset} -t {threads} --secondary=no -Y {ref} {reads} | "
           f"samtools sort -@ {SORT_THREADS} -O CRAM --reference {ref} -o {tmp_cram}")
    print(f"  {name}: {cmd}")
    started = datetime.now().astimezone().isoformat()
    t0 = time.monotonic()
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
    subprocess.run(["samtools", "index", str(tmp_cram)], check=True)
    wall = time.monotonic() - t0
    tmp_cram.replace(cram)
    tmp_crai.replace(crai)
    tmp_meta = meta.with_suffix(".json.tmp")
    tmp_meta.write_text(json.dumps({
        "cmd": cmd,
        "started": started,
        "wall_s": wall,
        "minimap2_version": minimap2_version(),
        "samtools_version": samtools_version(),
        "host": socket.gethostname(),
        "threads": threads,
    }, indent=2) + "\n")
    tmp_meta.replace(meta)


def do_compare(prefix, test, test_ref, oracle, ref, force):
    jpath = Path(str(prefix) + ".json")
    others = [Path(str(prefix) + s) for s in (".txt", ".mismatches.tsv")]
    if not force and jpath.exists() and all(p.exists() for p in others):
        jm = jpath.stat().st_mtime
        inputs = [Path(test), Path(oracle)]
        if all(p.exists() and p.stat().st_mtime <= jm for p in inputs):
            print(f"  {Path(prefix).name}: exists, skip")
            return
    print(f"  {Path(prefix).name}: comparing")
    write_result(compare(test, test_ref, oracle, ref), prefix)


def build_case(case, args):
    outdir = case["outdir"]
    outdir.mkdir(parents=True, exist_ok=True)
    old_ref = ensure_ref(case["old"])
    new_ref = ensure_ref(case["new"])
    threads = args.threads

    run_mapping("old", old_ref, case["reads"], case["preset"], outdir, threads, args.force)
    run_mapping("oracle", new_ref, case["reads"], case["preset"], outdir, threads, args.force)

    if not args.no_noise:
        noise_ref = outdir / "noise_ref.fa"
        if args.force or not noise_ref.exists():
            make_noise_ref(new_ref, noise_ref)
        noise_ref = ensure_ref(noise_ref)
        run_mapping("noise", noise_ref, case["reads"], case["preset"], outdir, threads, args.force)

    do_compare(outdir / "cmp_self", outdir / "oracle.cram", new_ref,
               outdir / "oracle.cram", new_ref, args.force)
    if not args.no_noise:
        do_compare(outdir / "cmp_noise", outdir / "noise.cram", noise_ref,
                   outdir / "oracle.cram", new_ref, args.force)


def discover_synthetic(root):
    syn = root / "synthetic"
    cases = []
    if not syn.is_dir():
        return cases
    for s in sorted(p for p in syn.iterdir() if p.is_dir()):
        old, new = s / "old.fa", s / "new.fa"
        if not old.exists() or not new.exists():
            continue
        for model in sorted(p for p in s.iterdir() if p.is_dir()):
            reads = model / "reads.fastq.gz"
            if not reads.exists():
                continue
            preset = PRESETS.get(model.name)
            if preset is None:
                print(f"skip unknown model: {s.name}/{model.name}")
                continue
            cases.append({
                "name": f"{s.name}/{model.name}",
                "old": old, "new": new, "reads": reads,
                "preset": preset, "outdir": model,
            })
    return cases


def discover_real(root):
    real = root / "real"
    cases = []
    manifest = real / "manifest.tsv"
    if not manifest.exists():
        return cases
    with open(manifest) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if not row or not row[0]:
                continue
            d = dict(zip(header, row))
            case_dir = real / d["case"]
            cases.append({
                "name": d["case"],
                "old": case_dir / "old.fa", "new": case_dir / "new.fa",
                "reads": case_dir / "reads.fastq.gz",
                "preset": d["preset"], "outdir": case_dir,
            })
    return cases


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["synthetic", "real"])
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-noise", action="store_true",
                    help="skip the noise-reference CRAM (not needed for M6 real cases)")
    ap.add_argument("--only", default=None, help="only cases whose name contains this substring")
    ap.add_argument("--root", default="data", help=argparse.SUPPRESS)
    args = ap.parse_args()

    root = Path(args.root).resolve()
    cases = discover_synthetic(root) if args.mode == "synthetic" else discover_real(root)
    if args.only:
        cases = [c for c in cases if args.only in c["name"]]
    if not cases:
        print(f"no cases found under {root / args.mode}")
        return

    failures = []
    for case in cases:
        print(f"=== {case['name']} ===")
        try:
            build_case(case, args)
        except Exception as exc:
            print(f"FAILED {case['name']}: {exc}", file=sys.stderr)
            failures.append(case["name"])
    if failures:
        print(f"{len(failures)} case(s) failed: {', '.join(failures)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
