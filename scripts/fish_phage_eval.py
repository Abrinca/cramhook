#!/usr/bin/env python3
"""Evaluate ``cramhook fish`` (and ``--bait``) on the phage demo cases.

Two built-in case sets:

* ``--dataset ont`` (default): ``data/v2/real/ont_<S>__polished-to-pbfinal``
  contain an ONT polished assembly (``old.fa``/``old.cram``) that lacks the phage
  contig, a PacBio-finalized assembly (``new.fa``) whose last contig *is* the
  phage, and an oracle remap of the same reads against ``new.fa``.
* ``--dataset pb``: PacBio HiFi reads mapped to the *raw* autocycler assembly
  (``assembly.fasta``/``assembly.cram``, which lacks the phage); the oracle is
  the same reads remapped to the finalized assembly and is built by
  ``--build-oracle`` (into ``scratch/bait/pb/``).  The finalized assemblies in
  the external tree were regenerated without the phage contig, so the oracle
  reference is built from ``scratch/pacbio-phage-ds/<S>-p1-1.1/with-phage.fasta``
  (or the external finalized) and ``bait.fa`` is appended when that base still
  lacks the phage contig.

In both, the phage is also integrated in the chromosome as a (near-identical)
prophage, so most reads from the free phage map cleanly to the prophage region
in the old assembly and ``fish`` never sees them as "disagreeing".

Per case this script:

1. derives the **truth** reads from the oracle CRAM (primary record, and
   separately *any* record, on the phage contig) and reports how many of them in
   the old CRAM are unmapped / clipped >= 50 bp / inside the prophage region
   without a clip;
2. runs several ``cramhook fish`` configurations and scores them against the
   truth (selected reads/bp, per reason, recall, and selected-not-truth);
3. reassembles each configuration's FASTQ (and the truth reads alone, as the
   upper bound) with ``myloasm`` (``--hifi`` for HiFi cases), maps the contigs
   back to the phage with ``minimap2 -cx asm20``, and reports coverage /
   identity / circularity;
4. writes ``summary.md`` plus machine-readable TSVs under
   ``scratch/bait/eval/<tag>/``.

Outputs are written only under ``scratch/bait/eval/<tag>/`` and
``scratch/bait/pb/``.  ``data/``, ``ref/`` and the external PacBio dataset root
(``$CRAMHOOK_DATA_PACBIO``) are never written to.

Usage::

  # ONT cases
  pixi run python scripts/fish_phage_eval.py --tag nb1 \\
      --cramhook /tmp/opencode/cramhook-nobait --configs abc

  # PacBio HiFi cases, all configs including --bait
  pixi run python scripts/fish_phage_eval.py --tag pb1 --dataset pb \\
      --cramhook /tmp/opencode/cramhook-bait --configs abcd --build-oracle

  # arbitrary cases from a JSON spec (@file or inline)
  pixi run python scripts/fish_phage_eval.py --tag x --cases @cases.json
"""
import argparse
import glob
import gzip
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BIN = ROOT / "target" / "release" / "cramhook"
EVAL_ROOT = ROOT / "scratch" / "bait" / "eval"

# fish reason order.  ``bait`` (added by the planned ``--bait`` option, config d)
# is not listed here; ``reason_cols`` appends any extra reason actually seen, so
# configs a-c keep exactly these columns and config d adds ``bait``.
REASONS = ["unmapped", "clip", "large_ins", "large_del", "linked", "locus"]

# Built-in demo cases.  A case dict has these fields (paths are relative to the
# repo root unless absolute):
#   name, old_cram, old_fa        old assembly + reads mapped to it
#   oracle_cram, oracle_fa        the reads mapped to the finalized assembly
#   bait                          phage/plasmid FASTA (its header names the phage)
#   prophage                      [contig, start0, end0] in the old assembly
#   preset                        "ont" (lr:hq) or "hifi" (map-hifi; myloasm --hifi)
#   qc_findings, qc_haplotypes    fish inputs from a previous qc run
# Optional: finalized_fa (source for --build-oracle), reads_fastq, phage.
#
# ONT cases: the old assembly is a polished assembly lacking the phage contig.
# The prophage intervals are in the old-assembly frame (0-based, half-open):
# FAM20853's is 40,884 bp (same length as the free phage); FAM29080's is a
# 20,350 bp region of a 20,262 bp phage (an inverted segment swap).
BUILTIN_CASES = [
    {
        "name": "FAM20853",
        "old_cram": "data/v2/real/ont_FAM20853__polished-to-pbfinal/old.cram",
        "old_fa": "data/v2/real/ont_FAM20853__polished-to-pbfinal/old.fa",
        "oracle_cram": "data/v2/real/ont_FAM20853__polished-to-pbfinal/oracle.cram",
        "oracle_fa": "data/v2/real/ont_FAM20853__polished-to-pbfinal/new.fa",
        "bait": "scratch/bait/FAM20853/bait.fa",
        "prophage": ["1", 201976, 242860],
        "preset": "ont",
        "qc_findings": "scratch/bait/FAM20853/qc.findings.bed",
        "qc_haplotypes": "scratch/bait/FAM20853/qc.haplotypes.tsv",
    },
    {
        "name": "FAM29080",
        "old_cram": "data/v2/real/ont_FAM29080__polished-to-pbfinal/old.cram",
        "old_fa": "data/v2/real/ont_FAM29080__polished-to-pbfinal/old.fa",
        "oracle_cram": "data/v2/real/ont_FAM29080__polished-to-pbfinal/oracle.cram",
        "oracle_fa": "data/v2/real/ont_FAM29080__polished-to-pbfinal/new.fa",
        "bait": "scratch/bait/FAM29080/bait.fa",
        "prophage": ["1", 543393, 563743],
        "preset": "ont",
        "qc_findings": "scratch/bait/FAM29080/qc.findings.bed",
        "qc_haplotypes": "scratch/bait/FAM29080/qc.haplotypes.tsv",
    },
]

# PacBio HiFi demo cases (--dataset pb): the old assembly is the *raw*
# autocycler assembly (assembly.fasta / assembly.cram, which lacks the phage);
# the oracle is the same reads remapped to the finalized assembly and is built
# by --build-oracle.  Everything is written under scratch/bait/pb/; the external
# PacBio dataset root (``$CRAMHOOK_DATA_PACBIO``) is read-only.
PB_ROOT = os.environ.get("CRAMHOOK_DATA_PACBIO", "")
PB_QC_GLOBS = {
    "FAM20853": "FAM20853-p1-1.1/3_autocycler/qc/qc-eced56fa*",
    "FAM29080": "FAM29080-p1-1.1/3_autocycler/qc/qc-471702c9*",
}
PB_PROPHAGE = {
    "FAM20853": ["FAM20853-p1-1.1_scf1", 1862513, 1903397],
    "FAM29080": ["FAM29080-p1-1.1_scf1", 273076, 293426],
}
# Finalized assembly to build the oracle against, and the reads FASTQ.  The
# finalized assemblies in the external PacBio dataset tree were regenerated on
# 2026-10-14 and lost the phage contig; build_oracle appends ``bait`` when the
# base lacks the phage contig (detected by the bait header name).  For FAM29080
# the user provided a with-phage copy (scratch/pacbio-phage-ds/with-phage.fasta);
# FAM20853 falls back to the external finalized + bait.  ``None`` reads means
# "read the path from the CRAM @PG CL line".
PB_FINALIZED = {
    "FAM20853": "scratch/pacbio-phage-ds/FAM20853-p1-1.1/with-phage.fasta",
    "FAM29080": "scratch/pacbio-phage-ds/FAM29080-p1-1.1/with-phage.fasta",
}
PB_FINALIZED_FALLBACK = {
    "FAM20853": f"{PB_ROOT}/FAM20853-p1-1.1/3_autocycler/consensus_assembly.finalized.fasta",
    "FAM29080": f"{PB_ROOT}/FAM29080-p1-1.1/3_autocycler/consensus_assembly.finalized.fasta",
}
PB_READS = {
    "FAM20853": None,
    "FAM29080": "scratch/pacbio-phage-ds/FAM29080-p1-1.1/m84151_260508_204155_s4.hifi_reads.bc2048--bc2048.hifi_reads.fastq.gz",
}
PB_READS_FALLBACK = {
    "FAM20853": None,
    "FAM29080": f"{PB_ROOT}/FAM29080-p1-1.1/0_reads/m84151_260508_204155_s4.hifi_reads.bc2048--bc2048.hifi_reads.fastq.gz",
}


def resolve_path(p):
    p = Path(p)
    return p if p.is_absolute() else (ROOT / p)


def display_path(p):
    """Path relative to the repo root when possible, else absolute."""
    p = Path(p)
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def normalize_case(c):
    """Resolve path fields and fill defaults for a case dict."""
    out = dict(c)
    for k in ("old_cram", "old_fa", "oracle_cram", "oracle_fa", "finalized_fa",
              "reads_fastq", "bait", "qc_findings", "qc_haplotypes"):
        if out.get(k):
            out[k] = resolve_path(out[k])
    if out.get("prophage") is not None:
        out["prophage"] = list(out["prophage"])
    out.setdefault("preset", "ont")
    return out


def pb_cases():
    """The built-in PacBio HiFi cases (qc dir resolved by glob)."""
    if not PB_ROOT:
        sys.exit("CRAMHOOK_DATA_PACBIO is not set; point it at the external PacBio dataset root")
    out = []
    for name, rel in PB_QC_GLOBS.items():
        matches = sorted(glob.glob(f"{PB_ROOT}/{rel}"))
        if not matches:
            continue
        qc = Path(matches[-1])
        fin = PB_FINALIZED[name]
        if not resolve_path(fin).exists():
            fin = PB_FINALIZED_FALLBACK[name]
        reads = PB_READS.get(name)
        if reads and not resolve_path(reads).exists():
            reads = PB_READS_FALLBACK.get(name)
        out.append(normalize_case({
            "name": name,
            "old_cram": str(qc / "assembly.cram"),
            "old_fa": str(qc / "assembly.fasta"),
            "oracle_cram": f"scratch/bait/pb/{name}/oracle.cram",
            "oracle_fa": f"scratch/bait/pb/{name}/finalized.fa",
            "finalized_fa": fin,
            "reads_fastq": reads,
            "bait": f"scratch/bait/{name}/bait.fa",
            "prophage": PB_PROPHAGE[name],
            "preset": "hifi",
            "qc_findings": f"scratch/bait/pb/{name}/qc.findings.bed",
            "qc_haplotypes": f"scratch/bait/pb/{name}/qc.haplotypes.tsv",
        }))
    return out


def load_cases(spec, builtin):
    """Resolve the ``--cases`` argument.

    ``spec`` is either a JSON case spec (``@file.json`` or an inline JSON list)
    or a comma-separated list of built-in case names (empty = all built-ins).
    """
    spec = (spec or "").strip()
    if spec.startswith("@") or spec.startswith("[") or spec.startswith("{"):
        data = json.loads(Path(spec[1:]).read_text()) if spec.startswith("@") \
            else json.loads(spec)
        if isinstance(data, dict):
            data = [data]
        return [normalize_case(c) for c in data]
    want = {x.strip() for x in spec.split(",") if x.strip()}
    return [normalize_case(c) for c in builtin
            if not want or c["name"] in want]


def reads_from_cram_pg(cram):
    """Reads FASTQ path from the minimap2 @PG CL line of a CRAM, or None."""
    try:
        hdr = pysam.AlignmentFile(str(cram), "rc",
                                  reference_filename=None).header.to_dict()
    except Exception:
        return None
    for pg in hdr.get("PG", []):
        cl = pg.get("CL", "")
        if "minimap2" in cl and ("-ax" in cl or "-a" in cl):
            toks = shlex.split(cl)
            return toks[-1] if toks else None
    return None


def fasta_has_contig(path, name):
    """True if a FASTA contains a contig whose first token is ``name``."""
    try:
        with open(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if line[1:].strip().split()[0] == name:
                        return True
    except OSError:
        return False
    return False


def build_oracle(case, threads, force, log):
    """Build ``case['oracle_cram']`` from the reads and the finalized assembly.

    Copies the finalized FASTA into scratch (so the ``.fai``/``.crai`` land
    there).  If the base lacks the phage contig (the finalized assemblies were
    regenerated without it), appends ``bait`` so the reference is the intended
    "finalized incl. phage".  Then runs exactly::

      minimap2 -ax map-hifi -t N --secondary=no -Y <reference> <reads>
        | samtools sort -@4 -O CRAM --reference <reference> -o <oracle.cram>
    """
    oc = Path(case["oracle_cram"])
    crai = Path(str(oc) + ".crai")
    if oc.exists() and crai.exists() and not force:
        return oc
    fin = Path(case["finalized_fa"])
    dst = Path(case["oracle_fa"])
    bait = Path(case["bait"])
    phage = case.get("phage") or first_header(bait)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.exists() or force:
        shutil.copyfile(fin, dst)
        if phage and not fasta_has_contig(dst, phage):
            with open(dst, "a") as out:
                out.write("\n")
                with open(bait) as b:
                    shutil.copyfileobj(b, out)
            print(f"[{case['name']}] oracle reference: appended bait "
                  f"`{phage}` to {display_path(dst)} (base lacked it)",
                  file=sys.stderr)
            if log is not None:
                log.write(f"# appended bait {phage} from {bait} to {dst}\n")
    if not Path(str(dst) + ".fai").exists() or force:
        run(["pixi", "run", "samtools", "faidx", str(dst)], log=log)
    reads = case.get("reads_fastq") or reads_from_cram_pg(case["old_cram"])
    if not reads or not Path(reads).exists():
        raise SystemExit(f"{case['name']}: cannot find reads FASTQ "
                         f"(set reads_fastq in the case spec)")
    oc.parent.mkdir(parents=True, exist_ok=True)
    cmd = (f"pixi run minimap2 -ax map-hifi -t {threads} --secondary=no -Y "
           f"{shlex.quote(str(dst))} {shlex.quote(str(reads))} | "
           f"pixi run samtools sort -@4 -O CRAM --reference "
           f"{shlex.quote(str(dst))} -o {shlex.quote(str(oc))}")
    t0 = time.perf_counter()
    p = subprocess.run(cmd, shell=True, cwd=str(ROOT), capture_output=True,
                       text=True)
    dt = time.perf_counter() - t0
    if log is not None:
        log.write(f"$ {cmd}\n# rc={p.returncode} elapsed={dt:.2f}s\n")
        if p.stderr:
            log.write(p.stderr[-4000:])
        log.write("\n")
    if p.returncode != 0:
        raise SystemExit(f"{case['name']}: oracle build failed")
    run(["pixi", "run", "samtools", "index", str(oc)], log=log)
    print(f"[{case['name']}] built oracle {oc} ({dt:.1f}s)", file=sys.stderr)
    return oc


def oracle_reference(case, out_dir, log):
    """A FASTA containing the phage contig, for decoding ``oracle_cram``.

    The prebuilt ONT oracle CRAMs still reference the phage contig, but the
    finalized assembly they point at was regenerated without it, so a reference
    lacking the phage is repaired by appending ``bait`` (into ``out_dir``).
    """
    fa = Path(case["oracle_fa"])
    bait = Path(case["bait"])
    phage = case.get("phage") or first_header(bait)
    if not phage or fasta_has_contig(fa, phage):
        return fa
    ref = out_dir / case["name"] / "oracle_ref.fa"
    ref.parent.mkdir(parents=True, exist_ok=True)
    if not ref.exists() or not Path(str(ref) + ".fai").exists():
        shutil.copyfile(fa, ref)
        with open(ref, "a") as out:
            out.write("\n")
            with open(bait) as b:
                shutil.copyfileobj(b, out)
        run(["pixi", "run", "samtools", "faidx", str(ref)], log=log)
        print(f"[{case['name']}] oracle reference repaired: appended bait "
              f"`{phage}` to {display_path(ref)}", file=sys.stderr)
    return ref


# Configurations.  ``extra`` are the fish arguments added to the base command
# (which always carries --findings/--haplotypes from qc, i.e. the default
# findings scope).  ``d`` additionally needs --bait and a bait-capable binary.
CONFIGS = {
    "a": {"desc": "default (findings scope, qc findings+haplotypes)", "extra": []},
    "b": {"desc": "a + --locus (default flank 5000)", "extra": ["--locus"]},
    "c": {"desc": "a + --scope genome", "extra": ["--scope", "genome"]},
    "d": {"desc": "a + --bait <bait.fa>", "extra": [], "bait": True},
}


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def pct(k, n):
    return "n/a" if not n else f"{100.0 * k / n:.1f}%"


def run(cmd, cwd=ROOT, timeout=None, log=None):
    """Run a command, return ``(CompletedProcess, elapsed_s)``.

    ``cmd`` may mix str and Path.  If ``log`` is a writable file object the
    command line, return code and output are appended to it.
    """
    argv = [str(c) for c in cmd]
    t0 = time.perf_counter()
    p = subprocess.run(argv, cwd=str(cwd), capture_output=True, text=True,
                       timeout=timeout)
    dt = time.perf_counter() - t0
    if log is not None:
        log.write(f"$ {' '.join(argv)}\n# rc={p.returncode} elapsed={dt:.2f}s\n")
        if p.stdout:
            log.write(p.stdout)
        if p.stderr:
            log.write(p.stderr)
        log.write("\n")
    return p, dt


def read_json(path):
    p = Path(path)
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text())
    except Exception:
        return {}


def first_header(path):
    """First FASTA header name (without ``>``) of a file, or None."""
    p = Path(path)
    if not p.exists():
        return None
    try:
        with open(p, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    return line[1:].strip().split()[0]
    except OSError:
        return None
    return None


def fasta_headers(path):
    out = []
    with open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                out.append(line[1:].strip())
    return out


def fasta_lengths(path):
    """``{contig name: length}`` for a FASTA."""
    lens, name, n = {}, None, 0
    with open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name = line[1:].strip().split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens


def open_cram(path, ref):
    path = Path(path)
    if str(path).endswith(".cram"):
        return pysam.AlignmentFile(str(path), "rc",
                                   reference_filename=str(ref))
    return pysam.AlignmentFile(str(path), "rb")


def counted(r):
    """Primary + supplementary (skip secondary/duplicate/qcfail/unmapped)."""
    return not (r.is_unmapped or r.is_secondary or r.is_duplicate
                or r.is_qcfail)


def clip_len(r):
    c = r.cigar
    if not c:
        return 0
    return max([ln for op, ln in c if op == 4] + [0])


def intervals_union(ivs):
    if not ivs:
        return 0
    ivs = sorted(ivs)
    total = 0
    cs, ce = ivs[0]
    for a, b in ivs[1:]:
        if a <= ce:
            ce = max(ce, b)
        else:
            total += ce - cs
            cs, ce = a, b
    return total + (ce - cs)


# ---------------------------------------------------------------------------
# step 1: truth
# ---------------------------------------------------------------------------
def truth_reads(oracle_cram, oracle_fa, phage):
    """``(primary_names, any_names)`` on the phage contig in oracle.cram."""
    prim, anyr = set(), set()
    with open_cram(oracle_cram, oracle_fa) as af:
        for r in af.fetch(phage):
            anyr.add(r.query_name)
            if not (r.is_unmapped or r.is_secondary or r.is_supplementary
                    or r.is_duplicate or r.is_qcfail):
                prim.add(r.query_name)
    return prim, anyr


def old_status(old_cram, old_fa, truth_any, pro_contig, pro_start, pro_end):
    """Per-truth-read status in the old CRAM.

    Returns ``(unmapped, clipped50, in_prophage, prophage_no_clip)`` name sets
    restricted to ``truth_any``.  A read is *clipped* when any counted record
    has a terminal soft clip >= 50 bp, and *prophage without clip* when a
    counted mapped record overlaps the prophage region and the read has no such
    clip.
    """
    unmapped, clipped, inpro, inpro_no = set(), set(), set(), set()
    with open_cram(old_cram, old_fa) as af:
        for r in af.fetch(until_eof=True):
            n = r.query_name
            if n not in truth_any or r.is_secondary or r.is_duplicate \
                    or r.is_qcfail:
                continue
            if r.is_unmapped:
                unmapped.add(n)
                continue
            cl = clip_len(r)
            if cl >= 50:
                clipped.add(n)
            if (r.reference_name == pro_contig
                    and r.reference_start < pro_end
                    and pro_start < r.reference_end):
                inpro.add(n)
                if cl < 50:
                    inpro_no.add(n)
    return unmapped, clipped, inpro, inpro_no


def extract_truth_fastq(old_cram, old_fa, names, out_fq, log):
    """FASTQ of the truth reads from the old CRAM (primary records, original
    orientation): ``samtools view -N names -F 0x904 | samtools fastq``."""
    out_fq = Path(out_fq)
    out_fq.parent.mkdir(parents=True, exist_ok=True)
    names_file = out_fq.parent / "truth_names.txt"
    names_file.write_text("\n".join(sorted(names)) + "\n")
    view = ["pixi", "run", "samtools", "view", "-h", "-N", str(names_file),
            "-F", "0x904", str(old_cram)]
    fastq = ["pixi", "run", "samtools", "fastq", "-"]
    t0 = time.perf_counter()
    with open(out_fq, "wb") as raw:
        p1 = subprocess.Popen([str(c) for c in view], cwd=str(ROOT),
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.run([str(c) for c in fastq], cwd=str(ROOT),
                            stdin=p1.stdout, stdout=raw,
                            stderr=subprocess.PIPE)
        p1.stdout.close()
        p1_err = p1.stderr.read().decode(errors="replace")
        p1.wait()
    dt = time.perf_counter() - t0
    if log is not None:
        log.write(f"$ {' '.join(str(c) for c in view)} | "
                  f"{' '.join(str(c) for c in fastq)}\n")
        log.write(f"# rc={p1.returncode}/{p2.returncode} elapsed={dt:.2f}s\n")
        if p1_err:
            log.write(p1_err)
        if p2.stderr:
            log.write(p2.stderr.decode(errors="replace"))
        log.write("\n")
    # BGZF-compress the plain FASTQ (samtools fastq wrote plain FASTQ).
    gz = out_fq.with_suffix(out_fq.suffix + ".gz")
    with open(out_fq, "rb") as fin, gzip.open(gz, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    out_fq.unlink()
    return gz, dt, p1.returncode, p2.returncode


# ---------------------------------------------------------------------------
# step 2: fish
# ---------------------------------------------------------------------------
def fish_paths(prefix):
    prefix = Path(prefix)
    return {
        "tsv": Path(str(prefix) + ".fish.tsv"),
        "json": Path(str(prefix) + ".fish.json"),
        "fastq": Path(str(prefix) + ".fastq.gz"),
    }


def parse_fish(tsv):
    rows = []
    p = Path(tsv)
    if not p.exists():
        return rows
    lines = p.read_text().splitlines()
    if not lines:
        return rows
    hdr = lines[0].split("\t")
    if hdr[:2] != ["read", "length"]:
        return rows
    idx = {h: i for i, h in enumerate(hdr)}
    for line in lines[1:]:
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        rows.append({k: (f[i] if i < len(f) else "") for k, i in idx.items()})
    return rows


def fish_selected(tsv):
    sel, by_reason = set(), defaultdict(set)
    for r in parse_fish(tsv):
        n = r.get("read", "")
        if not n:
            continue
        sel.add(n)
        by_reason[r.get("reason", "")].add(n)
    return sel, by_reason


def fish_totals(prefix):
    """``(reads, bp, {reason: (reads, bp)})`` from ``P.fish.json``."""
    j = read_json(fish_paths(prefix)["json"])
    rb = j.get("reads") if isinstance(j, dict) else None
    by = {rr: (0, 0) for rr in REASONS}
    if isinstance(rb, dict):
        for rr, v in (rb.get("by_reason") or {}).items():
            try:
                by[rr] = (int(v.get("reads", 0)), int(v.get("bp", 0)))
            except (TypeError, ValueError):
                by[rr] = (0, 0)
        try:
            return int(rb.get("reads", 0)), int(rb.get("bp", 0)), by
        except (TypeError, ValueError):
            pass
    sel, by_reason = fish_selected(str(prefix) + ".fish.tsv")
    lens = fastq_lengths(fish_paths(prefix)["fastq"])
    for rr in by_reason:
        by.setdefault(rr, (0, 0))
    by = {rr: (len(by_reason.get(rr, set())),
               sum(lens.get(n, 0) for n in by_reason.get(rr, set())))
          for rr in by}
    return len(sel), sum(lens.get(n, 0) for n in sel), by


def fastq_lengths(path):
    out = {}
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
                name = head[1:].strip().split()[0]
                if name:
                    out.setdefault(name, len(seq.rstrip("\n")))
    except OSError:
        return out
    return out


def has_bait(binpath):
    try:
        p = subprocess.run([str(binpath), "fish", "--help"], cwd=str(ROOT),
                           capture_output=True, text=True)
    except OSError:
        return False
    return "--bait" in (p.stdout + p.stderr)


def run_fish(binpath, old_cram, old_fa, extra, prefix, threads, force, log):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = fish_paths(prefix)
    if not force and paths["tsv"].exists():
        return 0.0, True
    for p in prefix.parent.glob(prefix.name + ".*"):
        try:
            p.unlink()
        except OSError:
            pass
    cmd = [binpath, "fish", "--in", old_cram, "--ref", old_fa,
           "--out-prefix", prefix, "--threads", str(threads)] + extra
    p, dt = run(cmd, log=log)
    return dt, paths["tsv"].exists()


# ---------------------------------------------------------------------------
# step 3: reassembly
# ---------------------------------------------------------------------------
def find_assembly(outdir):
    outdir = Path(outdir)
    if not outdir.is_dir():
        return None
    for name in ("assembly_primary.fa", "assembly_primary.fasta"):
        if (outdir / name).exists():
            return outdir / name
    cands = sorted(outdir.glob("*primary*.fa")) + sorted(outdir.glob("*.fa"))
    for p in cands:
        if p.parent == outdir:
            return p
    return None


def run_myloasm(fastq, outdir, threads, force, log, hifi=False):
    outdir = Path(outdir)
    asm = find_assembly(outdir)
    if not force and asm is not None:
        return asm, 0.0
    if outdir.exists() and force:
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = ["pixi", "run", "-e", "testing-myloasm", "myloasm", fastq,
           "-o", outdir, "-t", str(threads)]
    if hifi:
        cmd.append("--hifi")
    p, dt = run(cmd, log=log)
    return find_assembly(outdir), dt


def parse_paf(text):
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 12:
            continue
        de = None
        nm = None
        for tag in f[12:]:
            if tag.startswith("de:f:"):
                try:
                    de = float(tag[5:])
                except ValueError:
                    de = None
            elif tag.startswith("NM:i:"):
                try:
                    nm = int(tag[5:])
                except ValueError:
                    nm = None
        try:
            d = {
                "qname": f[0], "qlen": int(f[1]), "qstart": int(f[2]),
                "qend": int(f[3]), "strand": f[4], "tname": f[5],
                "tlen": int(f[6]), "tstart": int(f[7]), "tend": int(f[8]),
                "nmatch": int(f[9]), "alnlen": int(f[10]), "mapq": int(f[11]),
                "de": de, "nm": nm,
            }
        except ValueError:
            continue
        if d["de"] is not None:
            d["identity"] = 1.0 - d["de"]
        elif d["alnlen"]:
            d["identity"] = d["nmatch"] / d["alnlen"]
        else:
            d["identity"] = float("nan")
        rows.append(d)
    return rows


def map_to_phage(asm_fa, phage_fa, log):
    """Map assembly contigs to the phage; return per-contig coverage/identity."""
    if asm_fa is None:
        return {}, 0.0, None
    p, dt = run(["pixi", "run", "minimap2", "-cx", "asm20", "--secondary=no",
                 phage_fa, asm_fa], log=log)
    paf = parse_paf(p.stdout)
    phage_len = 0
    try:
        with open(phage_fa) as fh:
            seq = []
            for line in fh:
                if line.startswith(">"):
                    if seq:
                        break
                    continue
                seq.append(line.strip())
            phage_len = len("".join(seq))
    except OSError:
        pass
    per = {}
    for a in paf:
        d = per.setdefault(a["qname"], {"ivs": [], "id_wsum": 0.0,
                                        "aln_sum": 0})
        d["ivs"].append((a["tstart"], a["tend"]))
        if a["alnlen"]:
            d["id_wsum"] += a["identity"] * a["alnlen"]
            d["aln_sum"] += a["alnlen"]
    res = {}
    for name, d in per.items():
        cov_bp = intervals_union(d["ivs"])
        res[name] = {
            "cov_bp": cov_bp,
            "cov": (cov_bp / phage_len) if phage_len else 0.0,
            "identity": (d["id_wsum"] / d["aln_sum"]) if d["aln_sum"] else 0.0,
            "n_align": len(d["ivs"]),
        }
    return res, dt, phage_len


def reassemble(fastq, outdir, phage_fa, threads, force, log, hifi=False):
    """Run myloasm and score the contigs against the phage."""
    out = {"fastq": str(fastq), "n_contigs": 0, "best": None, "best_len": 0,
           "best_cov": 0.0, "best_identity": 0.0, "covers95": False,
           "circular": None, "small_best": None, "small_best_len": 0,
           "small_best_cov": 0.0, "small_best_identity": 0.0,
           "small_best_circular": None, "small_covers95": False,
           "myloasm_s": 0.0, "minimap2_s": 0.0, "total_cov": 0.0,
           "error": None}
    if fastq is None or not Path(fastq).exists():
        out["error"] = "no fastq"
        return out
    asm, dt_asm = run_myloasm(fastq, outdir, threads, force, log, hifi=hifi)
    out["myloasm_s"] = dt_asm
    if asm is None:
        out["error"] = "no myloasm assembly"
        return out
    headers = fasta_headers(asm)
    lengths = fasta_lengths(asm)
    out["n_contigs"] = len(headers)

    def circ_of(name):
        hdr = next((h for h in headers if h.split()[0] == name), "")
        if "circular-yes" in hdr:
            return "yes"
        if "circular-no" in hdr:
            return "no"
        if "circular-possibly" in hdr:
            return "possibly"
        return None

    per, dt_map, phage_len = map_to_phage(asm, phage_fa, log)
    out["minimap2_s"] = dt_map
    out["phage_len"] = phage_len
    if per:
        out["total_cov"] = min(1.0, sum(d["cov_bp"] for d in per.values())
                               / phage_len) if phage_len else 0.0
        best_name = max(per, key=lambda n: per[n]["cov_bp"])
        b = per[best_name]
        out["best"] = best_name
        out["best_len"] = lengths.get(best_name, 0)
        out["best_cov"] = b["cov"]
        out["best_identity"] = b["identity"]
        out["covers95"] = b["cov"] >= 0.95
        out["circular"] = circ_of(best_name)
        # "phage-scale" contigs: exclude whole-chromosome-scale contigs (which
        # trivially contain the prophage) when judging phage recovery.
        small = {n: d for n, d in per.items()
                 if lengths.get(n, 10 ** 18) <= 3 * phage_len}
        if small:
            sn = max(small, key=lambda n: small[n]["cov_bp"])
            s = small[sn]
            out["small_best"] = sn
            out["small_best_len"] = lengths.get(sn, 0)
            out["small_best_cov"] = s["cov"]
            out["small_best_identity"] = s["identity"]
            out["small_best_circular"] = circ_of(sn)
            out["small_covers95"] = s["cov"] >= 0.95
    return out


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def reason_cols(rows):
    """Reason column order: the known reasons plus any extra ones seen."""
    cols = list(REASONS)
    for r in rows:
        if r.get("skipped"):
            continue
        for k in (r.get("by_reason") or {}):
            if k not in cols:
                cols.append(k)
    return cols


def write_case_outputs(case_dir, tag, case, truth, fish_rows, reasm_rows,
                       binaries, threads):
    case_dir.mkdir(parents=True, exist_ok=True)
    md = []

    # ---- truth ------------------------------------------------------------
    md += [f"## {case['name']}", ""]
    md.append(f"- preset: `{truth.get('preset', 'ont')}`")
    md.append(f"- old assembly: `{truth['old_fa']}` (+ `{truth['old_cram']}`)")
    md.append(f"- oracle: `{truth['oracle_cram']}` mapped to "
              f"`{truth['oracle_fa']}`")
    if truth.get("finalized_fa"):
        md.append(f"  - oracle reference built from `{truth['finalized_fa']}`; "
                  f"`bait.fa` is appended when that base lacks the phage contig")
    md.append(f"- phage contig: `{truth['phage']}` "
              f"({truth['phage_len']} bp), bait `{truth['bait_fa']}`")
    md.append(f"- prophage in old assembly: `{truth['pro_contig']}:"
              f"{truth['pro_start']}-{truth['pro_end']}`")
    md.append("")
    md += ["### truth reads (oracle on the phage contig)", ""]
    md.append("| definition | reads | in old.cram unmapped | clipped >= 50 bp | "
              "in prophage | in prophage, no clip |")
    md.append("|---|---|---|---|---|---|")
    for label, key in (("primary", "primary"), ("any record", "any")):
        md.append(f"| {label} | {len(truth[key])} | {truth['old'][key]['unmapped']} "
                  f"| {truth['old'][key]['clipped']} "
                  f"| {truth['old'][key]['inpro']} "
                  f"| {truth['old'][key]['inpro_no_clip']} |")
    md.append("")
    hard = truth["old"]["any"]["inpro_no_clip"]
    md.append(f"The **{hard}** truth reads that map inside the prophage without a "
              f"clip are the hard case: the free phage is absent from the old "
              f"assembly but the integrated copy is (nearly) identical, so these "
              f"reads look like ordinary chromosome reads to `fish`.")
    missing = len(truth["any"]) - truth.get("n_extracted", 0)
    if missing:
        md.append("")
        md.append(f"Note: {missing} truth reads have no primary mapped record in "
                  f"`old.cram`, so `samtools view -F 0x904` omits them from the "
                  f"truth FASTQ (upper bound).")
    md.append("")

    # ---- fish -------------------------------------------------------------
    rcols = reason_cols(fish_rows)
    md += ["### fish configurations", ""]
    md.append("Base command for every config is `fish --in old.cram --ref old.fa "
              "--out-prefix <cfg>/fish --threads N --findings qc.findings.bed "
              "--haplotypes qc.haplotypes.tsv`; the extra flags below are added.")
    md.append("")
    md.append("| cfg | extra flags | selected reads | selected bp | "
              + " | ".join(f"{r} R/B" for r in rcols)
              + " | recall (primary) | recall (any) | not truth | time s |")
    md.append("|---" * (8 + len(rcols)) + "|")
    for r in fish_rows:
        if r.get("skipped"):
            cells = [r["cfg"], r["desc"], f"skipped: {r['skipped']}"] \
                + [""] * len(rcols) + ["", "", "", ""]
            md.append("| " + " | ".join(cells) + " |")
            continue
        per = " | ".join(
            f"{r['by_reason'].get(rr, (0, 0))[0]}/"
            f"{r['by_reason'].get(rr, (0, 0))[1]}" for rr in rcols)
        md.append(f"| {r['cfg']} | `{' '.join(r['extra']) or '-'}` "
                  f"| {r['selected_reads']} | {r['selected_bp']} | {per} "
                  f"| {r['recall_primary']} ({pct(r['tp'], r['truth_primary'])}) "
                  f"| {r['recall_any']} ({pct(r['ta'], r['truth_any'])}) "
                  f"| {r['not_truth']} | {r['time_s']:.1f} |")
    md.append("")
    md.append("`recall (primary)` / `recall (any)` = truth reads selected / truth "
              "reads, for the primary and any-record truth definitions.  "
              "`not truth` = selected reads that are not truth under the "
              "any-record definition.")
    md.append("")

    # ---- reassembly -------------------------------------------------------
    md += ["### reassembly (myloasm + minimap2 asm20 to the phage)", ""]
    md.append("| cfg | input | contigs | best contig | len | phage cov | identity "
              "| covers >= 95% | circular | phage-scale contig | len | cov | "
              "identity | circular | myloasm s | minimap2 s |")
    md.append("|---" * 16 + "|")
    for r in reasm_rows:
        cov = f"{100.0 * r['best_cov']:.1f}%" if r["best"] else "n/a"
        ident = f"{100.0 * r['best_identity']:.2f}%" if r["best"] else "n/a"
        circ = r["circular"] or "n/a"
        best = r["best"] or r.get("error") or "-"
        scov = (f"{100.0 * r['small_best_cov']:.1f}%"
                if r.get("small_best") else "n/a")
        sident = (f"{100.0 * r['small_best_identity']:.2f}%"
                  if r.get("small_best") else "n/a")
        md.append(
            f"| {r['cfg']} | {r['label']} | {r['n_contigs']} | `{best}` "
            f"| {r.get('best_len', 0)} | {cov} | {ident} "
            f"| {'yes' if r['covers95'] else 'no'} | {circ} "
            f"| `{r.get('small_best') or '-'}` | {r.get('small_best_len', 0)} "
            f"| {scov} | {sident} "
            f"| {r.get('small_best_circular') or 'n/a'} "
            f"| {r['myloasm_s']:.1f} | {r['minimap2_s']:.1f} |")
    md.append("")
    md.append("`best contig` = the assembled contig covering the most of the "
              "phage (this can be the whole chromosome, which contains the "
              "prophage); `phage cov` is the fraction of the phage covered by "
              "that single contig.  The **phage-scale contig** columns restrict "
              "the search to contigs no longer than 3x the phage, i.e. a contig "
              "that could actually *be* the phage.  `circular` is myloasm's "
              "`circular-yes`/`-no`/`-possibly` header flag.")
    md.append("")

    (case_dir / "summary.md").write_text("\n".join(md) + "\n")
    return md


def write_overall(out_dir, tag, cases, truth_rows, fish_rows, reasm_rows,
                  binaries, threads, configs):
    md = [f"# fish phage eval (tag `{tag}`)", ""]
    md.append(f"Binary (configs a-c): `{binaries['cramhook']}`")
    if binaries.get("bait_binary"):
        md.append(f"Binary (config d, --bait): `{binaries['bait_binary']}`")
    md.append(f"Threads: {threads}; configs: {', '.join(configs)}")
    md.append("")
    md.append("Truth = reads with a primary record (and separately *any* record) "
              "on the phage contig in `oracle.cram`.  Configs are the default "
              "findings-scope fish (`a`), `a`+`--locus` (`b`), `a`+`--scope "
              "genome` (`c`), and `a`+`--bait` (`d`).")
    md.append("")

    md += ["## Truth summary", ""]
    md.append("| case | preset | phage | phage bp | truth primary | truth any | "
              "prophage bp | truth in prophage no clip |")
    md.append("|---|---|---|---|---|---|---|---|")
    for t in truth_rows:
        md.append(f"| {t['name']} | {t.get('preset', 'ont')} "
                  f"| `{t['phage']}` | {t['phage_len']} "
                  f"| {len(t['primary'])} | {len(t['any'])} "
                  f"| {t['pro_end'] - t['pro_start']} "
                  f"| {t['old']['any']['inpro_no_clip']} |")
    md.append("")

    md += ["## fish recall (any-record truth)", ""]
    md.append("| case | cfg | extra | selected | recall any | recall primary "
              "| not truth |")
    md.append("|---|---|---|---|---|---|---|")
    for r in fish_rows:
        if r.get("skipped"):
            md.append(f"| {r['case']} | {r['cfg']} | {r['desc']} | skipped | | | |")
            continue
        md.append(f"| {r['case']} | {r['cfg']} | `{' '.join(r['extra']) or '-'}` "
                  f"| {r['selected_reads']} | {r['recall_any']} "
                  f"| {r['recall_primary']} | {r['not_truth']} |")
    md.append("")

    md += ["## reassembly (single contig covering the phage)", ""]
    md.append("| case | cfg | contigs | best cov | identity | covers >= 95% | "
              "circular | phage-scale cov | phage-scale covers >= 95% | "
              "phage-scale circular |")
    md.append("|---" * 10 + "|")
    for r in reasm_rows:
        cov = f"{100.0 * r['best_cov']:.1f}%" if r["best"] else "n/a"
        ident = f"{100.0 * r['best_identity']:.2f}%" if r["best"] else "n/a"
        circ = r["circular"] or "n/a"
        scov = (f"{100.0 * r['small_best_cov']:.1f}%"
                if r.get("small_best") else "n/a")
        md.append(f"| {r['case']} | {r['cfg']} | {r['n_contigs']} | {cov} "
                  f"| {ident} | {'yes' if r['covers95'] else 'no'} | {circ} "
                  f"| {scov} "
                  f"| {'yes' if r.get('small_covers95') else 'no'} "
                  f"| {r.get('small_best_circular') or 'n/a'} |")
    md.append("")
    md.append("`phage-scale` restricts to contigs no longer than 3x the phage, "
              "so the chromosome (which contains the prophage) does not count.")
    md.append("")

    # ---- conclusion -------------------------------------------------------
    md += ["## Conclusion", ""]
    for t in truth_rows:
        case_fish = [r for r in fish_rows if r["case"] == t["name"]
                     and not r.get("skipped")]
        case_re = [r for r in reasm_rows if r["case"] == t["name"]]
        truth_re = next((r for r in case_re if r["cfg"] == "truth"), None)
        re_by_cfg = {r["cfg"]: r for r in case_re}
        md.append(f"**{t['name']}.** Truth: {len(t['any'])} reads on the phage, "
                  f"of which {t['old']['any']['inpro_no_clip']} map inside the "
                  f"prophage without a clip in `old.cram` (the reads `fish` "
                  f"cannot see from the old assembly alone).")
        for r in case_fish:
            re = re_by_cfg.get(r["cfg"])
            re_txt = ""
            if re:
                if re.get("small_best"):
                    re_txt = (f"; reassembly recovers {100.0 * re['small_best_cov']:.0f}% "
                              f"of the phage on a phage-scale contig "
                              f"({re['small_best_len']} bp, "
                              f"circular={re['small_best_circular']})")
                else:
                    re_txt = "; reassembly produces no phage-scale contig"
            md.append(f"  - config `{r['cfg']}` ({' '.join(r['extra']) or 'default'}): "
                      f"{r['selected_reads']} reads selected, "
                      f"recall {r['recall_any']} (any) / "
                      f"{r['recall_primary']} (primary){re_txt}.")
        if truth_re:
            md.append(f"  - truth upper bound: myloasm from the truth reads gives "
                      f"{truth_re['n_contigs']} contigs; the phage-scale contig "
                      f"covers {100.0 * truth_re['small_best_cov']:.1f}% of the "
                      f"phage at {100.0 * truth_re['small_best_identity']:.2f}% "
                      f"identity (myloasm circular="
                      f"{truth_re['small_best_circular']}).")
        md.append("")
    md.append("The truth upper bound shows whether a single phage-scale contig can "
              "recover the phage at all; the configs show whether `fish` selects "
              "the reads needed to get there.  A config whose reassembly recovers "
              "a phage-scale contig covering >= 95% of the phage at high identity "
              "is the practical winner.")
    md.append("")

    (out_dir / "summary.md").write_text("\n".join(md) + "\n")
    return out_dir / "summary.md"


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def case_truth(case, tag, out_dir, threads, force, log):
    bait_fa = Path(case["bait"])
    oracle_cram = Path(case["oracle_cram"])
    oracle_fa = oracle_reference(case, out_dir, log)
    old_cram = Path(case["old_cram"])
    old_fa = Path(case["old_fa"])
    phage = case.get("phage") or first_header(bait_fa) or first_header(oracle_fa)
    prim, anyr = truth_reads(oracle_cram, oracle_fa, phage)
    pro_contig, pro_start, pro_end = case["prophage"]
    unm, clip, inpro, inpro_no = old_status(old_cram, old_fa, anyr, pro_contig,
                                            pro_start, pro_end)
    old = {}
    for label, names in (("primary", prim), ("any", anyr)):
        old[label] = {
            "unmapped": len(names & unm),
            "clipped": len(names & clip),
            "inpro": len(names & inpro),
            "inpro_no_clip": len(names & inpro_no),
        }
    phage_len = 0
    try:
        with open(bait_fa) as fh:
            seq = []
            for line in fh:
                if line.startswith(">"):
                    if seq:
                        break
                    continue
                seq.append(line.strip())
        phage_len = len("".join(seq))
    except OSError:
        pass
    truth = {
        "name": case["name"], "phage": phage, "phage_len": phage_len,
        "bait_fa": display_path(bait_fa),
        "primary": prim, "any": anyr, "old": old,
        "pro_contig": pro_contig, "pro_start": pro_start, "pro_end": pro_end,
        "preset": case.get("preset", "ont"),
        "old_cram": display_path(old_cram), "old_fa": display_path(old_fa),
        "oracle_cram": display_path(oracle_cram),
        "oracle_fa": display_path(oracle_fa),
        "finalized_fa": (display_path(case["finalized_fa"])
                         if case.get("finalized_fa") else None),
    }
    # truth FASTQ (upper bound).
    truth_fq = out_dir / case["name"] / "truth" / "truth.fastq"
    gz, dt, rc1, rc2 = extract_truth_fastq(old_cram, old_fa, anyr, truth_fq, log)
    truth["fastq"] = str(gz)
    truth["extract_s"] = dt
    truth["extract_rc"] = (rc1, rc2)
    truth["n_extracted"] = len(fastq_lengths(gz))
    return truth


def run_case(case, tag, out_dir, args, truth, configs, log):
    old_cram = Path(case["old_cram"])
    old_fa = Path(case["old_fa"])
    bait_fa = Path(case["bait"])
    qc_findings = Path(case["qc_findings"]) if case.get("qc_findings") else None
    qc_hap = Path(case["qc_haplotypes"]) if case.get("qc_haplotypes") else None
    hifi = str(case.get("preset", "")).lower() in ("hifi", "map-hifi")
    case_dir = out_dir / case["name"]

    base = []
    if qc_findings and qc_findings.exists():
        base += ["--findings", str(qc_findings)]
    if qc_hap and qc_hap.exists():
        base += ["--haplotypes", str(qc_hap)]

    fish_rows, reasm_rows = [], []

    for cfg in configs:
        spec = CONFIGS[cfg]
        extra = list(base) + list(spec["extra"])
        binpath = args.cramhook
        if spec.get("bait"):
            extra = list(base) + ["--bait", str(bait_fa)]
            binpath = args.bait_binary
        row = {"case": case["name"], "cfg": cfg, "desc": spec["desc"],
               "extra": [a for a in spec["extra"] if a != "<bait.fa>"]
               or (["--bait", display_path(bait_fa)]
                   if spec.get("bait") else [])}
        if spec.get("bait") and not (args.bait_binary and has_bait(args.bait_binary)):
            row["skipped"] = "no --bait in binary"
            fish_rows.append(row)
            continue
        if cfg == "b" and not (qc_findings and qc_findings.exists()):
            row["skipped"] = "no findings (--locus requires --findings)"
            fish_rows.append(row)
            continue
        prefix = case_dir / cfg / "fish"
        dt, ok = run_fish(binpath, old_cram, old_fa, extra, prefix,
                          args.threads, args.force, log)
        row["time_s"] = dt
        if not ok:
            row["skipped"] = "fish produced no tsv"
            fish_rows.append(row)
            continue
        sel, _ = fish_selected(fish_paths(prefix)["tsv"])
        reads, bp, by = fish_totals(prefix)
        row.update({
            "selected_reads": reads, "selected_bp": bp, "by_reason": by,
            "truth_primary": len(truth["primary"]),
            "truth_any": len(truth["any"]),
            "tp": len(sel & truth["primary"]),
            "ta": len(sel & truth["any"]),
            "not_truth": len(sel - truth["any"]),
            "recall_primary": pct(len(sel & truth["primary"]),
                                  len(truth["primary"])),
            "recall_any": pct(len(sel & truth["any"]), len(truth["any"])),
            "selected_set": sorted(sel),
        })
        fish_rows.append(row)

        # reassembly of this config's FASTQ
        if not args.skip_reassembly:
            r = reassemble(fish_paths(prefix)["fastq"], case_dir / cfg / "myloasm",
                           bait_fa, args.threads, args.force, log, hifi=hifi)
            r.update({"case": case["name"], "cfg": cfg, "label": "fish fastq"})
            reasm_rows.append(r)

    # truth upper bound reassembly
    if truth.get("fastq") and not args.skip_reassembly:
        r = reassemble(truth["fastq"], case_dir / "truth" / "myloasm", bait_fa,
                       args.threads, args.force, log, hifi=hifi)
        r.update({"case": case["name"], "cfg": "truth", "label": "truth reads"})
        reasm_rows.append(r)

    return fish_rows, reasm_rows


def write_tsvs(out_dir, truth_rows, fish_rows, reasm_rows):
    out_dir.mkdir(parents=True, exist_ok=True)
    # truth
    lines = ["case\tphage\tphage_bp\tpro_contig\tpro_start\tpro_end\t"
             "truth_primary\ttruth_any\t"
             "old_unmapped_primary\told_clipped_primary\told_inpro_primary\t"
             "old_inpro_noclip_primary\t"
             "old_unmapped_any\told_clipped_any\told_inpro_any\t"
             "old_inpro_noclip_any\ttruth_reads_extracted\textract_s"]
    for t in truth_rows:
        op, oa = t["old"]["primary"], t["old"]["any"]
        lines.append("\t".join(str(x) for x in [
            t["name"], t["phage"], t["phage_len"], t["pro_contig"],
            t["pro_start"], t["pro_end"], len(t["primary"]), len(t["any"]),
            op["unmapped"], op["clipped"], op["inpro"], op["inpro_no_clip"],
            oa["unmapped"], oa["clipped"], oa["inpro"], oa["inpro_no_clip"],
            t.get("n_extracted", ""), f"{t.get('extract_s', 0):.2f}"]))
    (out_dir / "truth.tsv").write_text("\n".join(lines) + "\n")

    # fish
    rcols = reason_cols(fish_rows)
    lines = ["case\tcfg\textra\tselected_reads\tselected_bp\t"
             + "\t".join(f"{r}_reads\t{r}_bp" for r in rcols)
             + "\ttruth_primary\ttruth_any\ttruth_primary_selected\t"
             "truth_any_selected\tnot_truth\trecall_primary\trecall_any\ttime_s"]
    for r in fish_rows:
        if r.get("skipped"):
            lines.append("\t".join([r["case"], r["cfg"],
                                    " ".join(r["extra"]), "skipped"]))
            continue
        row = [r["case"], r["cfg"], " ".join(r["extra"]), r["selected_reads"],
               r["selected_bp"]]
        for rr in rcols:
            row += [r["by_reason"].get(rr, (0, 0))[0],
                    r["by_reason"].get(rr, (0, 0))[1]]
        row += [r["truth_primary"], r["truth_any"], r["tp"], r["ta"],
                r["not_truth"], r["recall_primary"], r["recall_any"],
                f"{r['time_s']:.2f}"]
        lines.append("\t".join(str(x) for x in row))
    (out_dir / "fish.tsv").write_text("\n".join(lines) + "\n")

    # reassembly
    lines = ["case\tcfg\tlabel\tn_contigs\tbest_contig\tbest_len\tbest_cov\t"
             "best_identity\tcovers95\tcircular\tsmall_best_contig\t"
             "small_best_len\tsmall_best_cov\tsmall_best_identity\t"
             "small_covers95\tsmall_best_circular\ttotal_cov\tphage_bp\t"
             "myloasm_s\tminimap2_s\terror"]
    for r in reasm_rows:
        lines.append("\t".join(str(x) for x in [
            r["case"], r["cfg"], r["label"], r["n_contigs"], r["best"] or "",
            r.get("best_len", 0),
            f"{r['best_cov']:.4f}", f"{r['best_identity']:.4f}",
            int(bool(r["covers95"])),
            r["circular"] or "",
            r.get("small_best") or "", r.get("small_best_len", 0),
            f"{r.get('small_best_cov', 0):.4f}",
            f"{r.get('small_best_identity', 0):.4f}",
            int(bool(r.get("small_covers95"))),
            r.get("small_best_circular") or "",
            f"{r.get('total_cov', 0):.4f}", r.get("phage_len", ""),
            f"{r['myloasm_s']:.2f}", f"{r['minimap2_s']:.2f}",
            r.get("error") or ""]))
    (out_dir / "reassembly.tsv").write_text("\n".join(lines) + "\n")


def parse_configs(s):
    out = []
    for ch in s.replace(",", "").replace(" ", ""):
        if ch in CONFIGS and ch not in out:
            out.append(ch)
    return out or ["a", "b", "c"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default=None,
                    help="output tag: scratch/bait/eval/<tag>/ (required unless "
                         "--self-test)")
    ap.add_argument("--cramhook", default=str(DEFAULT_BIN),
                    help="cramhook binary for configs a-c (e.g. the no-bait copy)")
    ap.add_argument("--bait-binary", default=None,
                    help="cramhook binary with --bait for config d (skipped if "
                         "the binary has no --bait)")
    ap.add_argument("--configs", default="abc",
                    help="configs to run, e.g. 'abc' or 'a,b,c,d'")
    ap.add_argument("--dataset", choices=["ont", "pb"], default="ont",
                    help="built-in case set: 'ont' (polished assemblies) or "
                         "'pb' (PacBio HiFi, raw autocycler assembly); ignored "
                         "for a JSON case spec")
    ap.add_argument("--cases", default=None,
                    help="comma-separated built-in case names (default: all of "
                         "the dataset), or a JSON case spec: @file.json or an "
                         "inline JSON list of case objects")
    ap.add_argument("--build-oracle", action="store_true",
                    help="build the oracle CRAM (reads remapped to the finalized "
                         "assembly) for cases that define finalized_fa")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--force", action="store_true",
                    help="re-run fish and myloasm even if outputs exist")
    ap.add_argument("--skip-reassembly", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()

    if not args.tag:
        ap.error("--tag is required (unless --self-test)")

    configs = parse_configs(args.configs)
    builtin = pb_cases() if args.dataset == "pb" else BUILTIN_CASES
    cases = load_cases(args.cases, builtin)
    if not cases:
        ap.error(f"no cases selected (dataset={args.dataset}, "
                 f"--cases={args.cases!r})")
    # Config d uses --bait; default the bait binary to the main binary.
    args.bait_binary = args.bait_binary or args.cramhook
    out_dir = EVAL_ROOT / args.tag
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = out_dir / "run.log"
    with open(log_path, "a") as log:
        log.write(f"\n===== {time.strftime('%Y-%m-%d %H:%M:%S')} tag={args.tag} "
                  f"configs={''.join(configs)} dataset={args.dataset} "
                  f"cases={[c['name'] for c in cases]} =====\n")
        binaries = {"cramhook": args.cramhook,
                    "bait_binary": args.bait_binary}
        if args.build_oracle:
            for case in cases:
                if case.get("finalized_fa"):
                    build_oracle(case, args.threads, args.force, log)
        all_truth, all_fish, all_reasm = [], [], []
        for case in cases:
            print(f"[{case['name']}] truth ...", file=sys.stderr)
            truth = case_truth(case, args.tag, out_dir, args.threads,
                               args.force, log)
            print(f"[{case['name']}] truth: primary={len(truth['primary'])} "
                  f"any={len(truth['any'])} extracted={truth.get('n_extracted')}",
                  file=sys.stderr)
            all_truth.append(truth)
            fr, rr = run_case(case, args.tag, out_dir, args, truth, configs, log)
            all_fish += fr
            all_reasm += rr
            for r in fr:
                if r.get("skipped"):
                    print(f"[{case['name']}] cfg {r['cfg']}: skipped "
                          f"({r['skipped']})", file=sys.stderr)
                else:
                    print(f"[{case['name']}] cfg {r['cfg']}: "
                          f"{r['selected_reads']} selected, recall any "
                          f"{r['recall_any']}", file=sys.stderr)

        write_tsvs(out_dir, all_truth, all_fish, all_reasm)
        for t in all_truth:
            write_case_outputs(out_dir / t["name"], args.tag,
                               next(c for c in cases if c["name"] == t["name"]),
                               t,
                               [r for r in all_fish if r["case"] == t["name"]],
                               [r for r in all_reasm if r["case"] == t["name"]],
                               binaries, args.threads)
        summary = write_overall(out_dir, args.tag, cases, all_truth, all_fish,
                                all_reasm, binaries, args.threads, configs)
    print(f"wrote {summary}")
    print(f"wrote {out_dir}/truth.tsv, fish.tsv, reassembly.tsv, summary.md")
    return 0


# ---------------------------------------------------------------------------
# self-test (synthetic, no data or tools)
# ---------------------------------------------------------------------------
def self_test():
    import tempfile
    problems = []

    # PAF parsing + per-contig coverage/identity
    paf = (
        "ctg1\t5000\t0\t5000\t+\tphage\t1000\t0\t1000\t1000\t1000\t60\tde:f:0.01\n"
        "ctg2\t2000\t0\t400\t+\tphage\t1000\t100\t500\t400\t400\t60\tde:f:0.05\n"
        "ctg2\t2000\t0\t400\t+\tphage\t1000\t600\t1000\t400\t400\t60\tde:f:0.05\n"
    )
    rows = parse_paf(paf)
    if len(rows) != 3:
        problems.append(f"parse_paf rows {len(rows)}")
    per = {}
    for a in rows:
        d = per.setdefault(a["qname"], {"ivs": [], "id_wsum": 0.0, "aln_sum": 0})
        d["ivs"].append((a["tstart"], a["tend"]))
        d["id_wsum"] += a["identity"] * a["alnlen"]
        d["aln_sum"] += a["alnlen"]
    cov1 = intervals_union(per["ctg1"]["ivs"])
    cov2 = intervals_union(per["ctg2"]["ivs"])
    if cov1 != 1000 or cov2 != 800:
        problems.append(f"coverage union {cov1} {cov2}")
    id1 = per["ctg1"]["id_wsum"] / per["ctg1"]["aln_sum"]
    if abs(id1 - 0.99) > 1e-9:
        problems.append(f"identity {id1}")

    # fish.tsv parsing
    with tempfile.TemporaryDirectory() as td:
        tsv = Path(td) / "fish.fish.tsv"
        tsv.write_text(
            "read\tlength\treason\tfinding\tcontig\tpos\tstrand\trecord\tdetail\tsa\n"
            "r1\t500\tclip\tfinding_2\t1\t543393\t+\tprimary\tside=left;len=60\t.\n"
            "r1\t500\tlinked\tfinding_2\t1\t543393\t+\tprimary\tgroup=1\t.\n"
            "r2\t400\tunmapped\t.\t.\t.\t.\tprimary\t.\t.\n")
        sel, by = fish_selected(tsv)
        if sel != {"r1", "r2"} or by["clip"] != {"r1"} or by["linked"] != {"r1"}:
            problems.append(f"fish_selected {sel} {dict(by)}")

    # interval union edge cases
    if intervals_union([(0, 10), (5, 20), (30, 40)]) != 30:
        problems.append("intervals_union overlap")
    if intervals_union([]) != 0:
        problems.append("intervals_union empty")

    # config parsing
    if parse_configs("a,b,c,d") != ["a", "b", "c", "d"]:
        problems.append("parse_configs comma")
    if parse_configs("abcd") != ["a", "b", "c", "d"]:
        problems.append("parse_configs string")

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        return 1
    print("self-test OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
