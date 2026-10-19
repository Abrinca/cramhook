#!/usr/bin/env python3
"""Populate an ONT test dataset (symlinks + manifest) under ``--root``. Run from repo root.

Default (``--root data/m6 --kind all``) reproduces the M6 set: the 11 ONT samples with
4 cases each (``cons-to-polished`` + 3 raw ``<primary>-to-cons``), reads > 300 MB
downsampled to ~100 MB.

``--kind raw-to-polished`` builds the "raw primary -> polished assembly" cases
(``<primary>-to-polished``) instead, e.g. the full-size M6-large set under
``data/m6large``.

``--kind m7`` builds the M7 "real use case" set (e.g. under ``data/v2``): 6 ONT
samples × 3 cases plus 6 PacBio cases, all with full-size reads.

* The source tree ($CRAMHOOK_DATA_ONT, the ONT dataset root) is READ-ONLY.  Primaries are
  ungzipped into ``<root>/assemblies/`` (or reused from an existing copy under
  ``data/m6/assemblies/``); downsamples are written once to ``<root>/reads/``.
* ``--no-downsample`` always symlinks the original full-size reads.
* Existing minimap2 outputs are reused *per mapping* (``old`` / ``oracle``) whenever the
  resolved reference and reads are byte-identical, the preset and minimap2 version match,
  and the source lies in data/real, data/val/real or another case of the same set.  Reused
  CRAMs are hardlinked and their JSON ``cmd`` rewritten case-local via ``rewrite_cmd``.

Idempotent: rerunning rebuilds the manifest/symlinks and refreshes reused CRAMs.
"""
import argparse
import gzip
import hashlib
import json
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _dataset_root(var):
    """A local source dataset root from the environment (never hard-coded)."""
    value = os.environ.get(var)
    if not value:
        sys.exit(f"{var} is not set; point it at a local source dataset root")
    return Path(value)


SRC = None          # resolved in main(); only modes that read ONT data need it
DATA = ROOT / "data"

PLATFORM, PRESET, FINAL = "ont", "lr:hq", "polished"
N_PRIMARIES, LEN_TOL, SEED = 3, 0.05, 6
DOWNSAMPLE_MIN_BYTES = 300 * 1024 * 1024   # > 300 MB (measured with du -L)
TARGET_BYTES = 100 * 1024 * 1024            # ~100 MB gz

KIND_ALL, KIND_RAW = "all", "raw-to-polished"
KIND_M7 = "m7"

# M7 "real use case" set sources (resolved in main(); only --kind m7 needs it).
PACBIO_SRC = None
ONT_M7_SAMPLES = ["FAM20853", "FAM20856", "FAM20862", "FAM20866", "FAM22000",
                  "FAM29080"]
# PacBio samples for M7.  Selection rule: the 6 largest read files in PACBIO_SRC
# (<S>-p1-1.1/0_reads/*.fastq.gz, by gzipped size) whose sample is not already used by
# data/real or data/val/real (which supply FAM14197/FAM18987/FAM20843/FAM28538/FAM29030/
# FAM29046 -p1-1.1 to the M6 sets).  data-2025-12-03 reads are not available locally, so
# this is derived from PACBIO_SRC alone.
PACBIO_M7_SAMPLES = ["FAM20855", "FAM29028", "FAM20848", "FAM20840", "FAM20866",
                     "FAM20839"]
ROTATION_SEED = 7

SAMPLES = ["FAM14197", "FAM20838", "FAM20843", "FAM20847", "FAM20851", "FAM20853",
           "FAM20856", "FAM20862", "FAM20866", "FAM22000", "FAM29080"]

# Source manifests used to recover the previously-chosen primaries.  The m6 manifest
# records the seed-6 picks for the 5 new samples; data/real and data/val/real cover the
# other 6.
PRIMARIES_MANIFESTS = [DATA / "real" / "manifest.tsv",
                       DATA / "val" / "real" / "manifest.tsv",
                       DATA / "m6" / "real" / "manifest.tsv"]

# Target dataset paths (mutated from --root in main()).
DATASET = DATA / "m6"
OUT = DATASET / "real"
READS_OUT = DATASET / "reads"
ASM_OUT = DATASET / "assemblies"
ASM_SEARCH = [ASM_OUT, DATA / "m6" / "assemblies"]

_CONTENT_CACHE = {}


def total_len(fa):
    with (gzip.open(fa, "rt") if fa.suffix == ".gz" else open(fa)) as fh:
        return sum(len(line.strip()) for line in fh if not line.startswith(">"))


def link(target, path):
    """Symlink path -> resolved target (replace any existing entry)."""
    path.unlink(missing_ok=True)
    path.symlink_to(Path(target).resolve())


def existing_primaries():
    """sample -> sorted primary names used by the pre-existing ONT cases."""
    chosen = {}
    for manifest in PRIMARIES_MANIFESTS:
        if not manifest.exists():
            continue
        for line in manifest.read_text().splitlines()[1:]:
            row = line.split("\t")
            if len(row) < 7:
                continue
            case, platform, sample = row[0], row[1], row[3]
            if platform != "ont" or "__cons-to-" in case:
                continue
            chosen.setdefault(sample, set()).add(case.split("__", 1)[1].split("-to-")[0])
    return {s: sorted(v) for s, v in chosen.items()}


def chosen_primaries(sample, cons_len, previously_used):
    """The 3 primaries for a sample: reuse known picks, else filter + seed-6 sample."""
    if sample in previously_used:
        return previously_used[sample]
    primaries = sorted((SRC / sample / "2_assemblies").glob("*.fasta.gz"))
    ok = [p for p in primaries if abs(total_len(p) - cons_len) <= LEN_TOL * cons_len]
    rng = random.Random(SEED)
    return sorted(p.name.removesuffix(".fasta.gz")
                  for p in rng.sample(ok, min(N_PRIMARIES, len(ok))))


def ensure_local_assembly(sample, primary):
    """Path to an ungzipped primary: reuse an existing copy, else ungzip into ASM_OUT."""
    for d in ASM_SEARCH:
        local = d / sample / f"{primary}.fasta"
        if local.exists():
            return local
    src = SRC / sample / "2_assemblies" / f"{primary}.fasta.gz"
    assert src.exists(), f"missing primary {src}"
    local = ASM_OUT / sample / f"{primary}.fasta"
    local.parent.mkdir(parents=True, exist_ok=True)
    tmp = local.with_name(local.name + ".tmp")
    with gzip.open(src, "rb") as fin, open(tmp, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    tmp.replace(local)
    return local


def gz_bytes(path):
    return Path(path).resolve().stat().st_size


def read_stats(fastq_gz):
    """(n_reads, n_bases) for a (gzipped) FASTQ, via seqkit stats."""
    out = subprocess.run(["seqkit", "stats", "-T", "-j", "8", str(fastq_gz)],
                         capture_output=True, text=True, check=True).stdout
    fields = out.splitlines()[1].split("\t")
    return int(fields[3]), int(fields[4])


def downsample(sample, reads):
    """Write READS_OUT/<sample>.fastq.gz (~100 MB) once; return (path, fraction)."""
    target = READS_OUT / f"{sample}.fastq.gz"
    fraction = TARGET_BYTES / gz_bytes(reads)
    if target.exists():
        return target, fraction
    tmp = READS_OUT / f"{sample}.fastq.gz.tmp"
    cmd = f"seqtk sample -s {SEED} {reads} {fraction:.10f} | gzip > {tmp}"
    print(f"  downsampling {sample}: {cmd}")
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
    tmp.replace(target)
    return target, fraction


def _content_hash(path):
    rp = Path(path).resolve()
    st = rp.stat()
    key = (str(rp), st.st_size, st.st_mtime_ns)
    h = _CONTENT_CACHE.get(key)
    if h is None:
        h = hashlib.sha256()
        with open(rp, "rb") as fh:
            for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
                h.update(chunk)
        h = h.hexdigest()
        _CONTENT_CACHE[key] = h
    return st.st_size, h


def same_content(a, b):
    ra, rb = Path(a).resolve(), Path(b).resolve()
    if ra == rb:
        return True
    if ra.stat().st_size != rb.stat().st_size:
        return False
    return _content_hash(ra) == _content_hash(rb)


FASTA_EXT = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")
FASTQ_EXT = (".fq", ".fastq", ".fq.gz", ".fastq.gz")


def _is_fasta(p):
    return p.lower().endswith(FASTA_EXT)


def _is_fastq(p):
    return p.lower().endswith(FASTQ_EXT)


def parse_cmd_inputs(cmd):
    """(refs, reads) sequence paths named in a recorded minimap2/samtools cmd."""
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


def rewrite_cmd(cmd, ref, reads, out):
    """Repoint every reference arg, the reads arg and ``-o`` at case-local paths.

    Idempotent: applying it to an already-local cmd is a no-op.
    """
    refs, reads_toks = parse_cmd_inputs(cmd)
    new = cmd
    for t in dict.fromkeys(reads_toks):
        new = new.replace(t, str(reads))
    for t in dict.fromkeys(refs):
        new = new.replace(t, str(ref))
    new = re.sub(r"(-o\s+)(\S+)", lambda m: m.group(1) + str(out), new, count=1)
    return new


# ---------------------------------------------------------------------------
# mapping command / provenance validation
# ---------------------------------------------------------------------------
_FASTA_MD5_CACHE = {}
_CRAM_HDR_CACHE = {}


def _split_pipeline(cmd):
    parts = cmd.split("|", 1)
    return parts[0], (parts[1] if len(parts) > 1 else "")


def parse_minimap2_cmd(cmd):
    """(preset, ref, reads) for a strict minimap2 alignment command, else None.

    Accepts only the binary ``minimap2``, ``-t <int>``, ``--secondary=no``, ``-Y``,
    an optional ``-a`` and exactly one ``-x <preset>`` (the ``-ax``/``-xa`` short
    forms are tolerated), followed by exactly the reference and reads paths.  Any
    other option is rejected as uncertain.
    """
    try:
        toks = shlex.split(cmd)
    except ValueError:
        return None
    if not toks or Path(toks[0]).name != "minimap2":
        return None
    preset = None
    positionals = []
    i = 1
    while i < len(toks):
        t = toks[i]
        if t == "-t":
            if i + 1 >= len(toks) or not toks[i + 1].lstrip("-").isdigit():
                return None
            i += 2
            continue
        if t in ("--secondary=no", "-Y"):
            i += 1
            continue
        if t.startswith("-") and not t.startswith("--"):
            flags = t[1:]
            if not flags or set(flags) - {"a", "x"}:
                return None
            if "x" in flags:
                if i + 1 >= len(toks):
                    return None
                preset = toks[i + 1]
                i += 2
            else:
                i += 1
            continue
        if t.startswith("--"):
            return None
        positionals.append(t)
        i += 1
    if preset is None or len(positionals) != 2:
        return None
    return preset, positionals[0], positionals[1]


def fasta_md5s(fa):
    """Ordered ``[(name, length, md5)]`` for a FASTA; md5 of the uppercased sequence."""
    rp = Path(fa).resolve()
    st = rp.stat()
    key = (str(rp), st.st_size, st.st_mtime_ns)
    if key in _FASTA_MD5_CACHE:
        return _FASTA_MD5_CACHE[key]
    out = []
    name, h, ln = None, None, 0
    opener = gzip.open if rp.suffix == ".gz" else open
    with opener(rp, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, ln, h.hexdigest()))
                name = line[1:].split()[0]
                h, ln = hashlib.md5(), 0
            else:
                s = line.strip()
                h.update(s.upper().encode())
                ln += len(s)
        if name is not None:
            out.append((name, ln, h.hexdigest()))
    _FASTA_MD5_CACHE[key] = out
    return out


def cram_header(cram):
    """``(sq, pg_lines)`` parsed from a CRAM/BAM header (cached)."""
    rp = Path(cram).resolve()
    st = rp.stat()
    key = (str(rp), st.st_size, st.st_mtime_ns)
    if key in _CRAM_HDR_CACHE:
        return _CRAM_HDR_CACHE[key]
    text = subprocess.run(["samtools", "view", "-H", str(rp)],
                          capture_output=True, text=True, check=True).stdout
    sq, pg = [], []
    for line in text.splitlines():
        if line.startswith("@SQ"):
            fields = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
            sq.append((fields.get("SN"), int(fields.get("LN", -1)), fields.get("M5")))
        elif line.startswith("@PG"):
            pg.append(line)
    _CRAM_HDR_CACHE[key] = (sq, pg)
    return sq, pg


def sq_matches_fasta(cram, ref):
    """True when every @SQ name/length/M5 equals the FASTA's sequences."""
    sq, _ = cram_header(cram)
    want = fasta_md5s(ref)
    return len(sq) == len(want) and all(
        sn == wn and ln == wl and (m5 or "") == wm
        for (sn, ln, m5), (wn, wl, wm) in zip(sq, want))


def pg_minimap2_reads(cram):
    """Reads path from the @PG minimap2 CL, or None."""
    _, pg = cram_header(cram)
    for line in pg:
        fields = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
        if fields.get("ID") == "minimap2" or fields.get("PN") == "minimap2":
            parsed = parse_minimap2_cmd(fields.get("CL", ""))
            if parsed:
                return parsed[2]
    return None


def mapping_contract_ok(case_dir, mapping, ref, reads, mm2, preset):
    """``(ok, reason)`` for an existing ``<mapping>`` in ``case_dir`` against ref/reads.

    Requires the full command shape/preset (``preset`` is the case's manifest preset,
    e.g. ``lr:hq`` for ONT or ``map-hifi`` for PacBio), matching minimap2 ref+reads and
    samtools ``--reference``, the recorded minimap2 version, CRAM @SQ M5/name/length
    provenance, and current inputs (or, for legacy cases without ``inputs.json``, a
    @PG minimap2 reads path that resolves to, or has the content of, the reads).
    """
    cram = case_dir / f"{mapping}.cram"
    crai = case_dir / f"{mapping}.cram.crai"
    jp = case_dir / f"{mapping}.json"
    if not (cram.exists() and crai.exists() and jp.exists()):
        return False, "missing artifacts"
    meta = read_json_file(jp)
    if not meta or not meta.get("cmd"):
        return False, "unreadable json"
    if meta.get("minimap2_version") != mm2:
        return False, "minimap2 version mismatch"
    head, tail = _split_pipeline(meta["cmd"])
    parsed = parse_minimap2_cmd(head)
    if parsed is None:
        return False, "unexpected minimap2 command shape"
    cmd_preset, cmd_ref, cmd_reads = parsed
    if cmd_preset != preset:
        return False, "preset mismatch"
    if not re.search(r"(^|\s)samtools(\s|$)", tail) or not re.search(r"(^|\s)sort(\s|$)", tail):
        return False, "missing samtools sort"
    m = re.search(r"--reference\s+(\S+)", tail)
    if not m:
        return False, "missing samtools --reference"
    if not same_content(ref, cmd_ref):
        return False, "minimap2 ref != requested ref"
    if not same_content(ref, m.group(1)):
        return False, "samtools --reference != requested ref"
    if not same_content(reads, cmd_reads):
        return False, "minimap2 reads != requested reads"
    if not sq_matches_fasta(cram, ref):
        return False, "@SQ M5/name/length mismatch"
    if (case_dir / "inputs.json").exists():
        stored = read_inputs_json(case_dir)
        if stored != input_signatures(case_dir / "old.fa", case_dir / "new.fa",
                                      case_dir / "reads.fastq.gz"):
            return False, "inputs.json does not match current inputs"
    else:
        pg_reads = pg_minimap2_reads(cram)
        if pg_reads is None:
            return False, "no @PG minimap2 reads path"
        if not (Path(pg_reads).resolve() == Path(reads).resolve()
                or same_content(pg_reads, reads)):
            return False, "@PG minimap2 reads mismatch"
    return True, "ok"


def invalidate_mapping(case_dir, mapping):
    """Delete one mapping's completion artifacts so build_crams rebuilds it."""
    removed = []
    for name in (f"{mapping}.cram", f"{mapping}.cram.crai", f"{mapping}.json"):
        p = case_dir / name
        if p.exists():
            p.unlink()
            removed.append(name)
    return removed


def invalidate_outputs(case_dir):
    """Delete completed CRAMs/JSONs (and stale reference .fai/copied .fa) so
    build_crams rebuilds from current inputs."""
    removed = []
    for name in ("old.cram", "old.cram.crai", "oracle.cram", "oracle.cram.crai",
                 "old.json", "oracle.json"):
        p = case_dir / name
        if p.exists():
            p.unlink()
            removed.append(name)
    # Reference indexes derived from a changed reference must not survive; drop
    # them plus any case-local *copied* reference (symlinks are refreshed below).
    for name in ("old.fa.fai", "new.fa.fai"):
        p = case_dir / name
        if p.exists() or p.is_symlink():
            p.unlink()
            removed.append(name)
    for name in ("old.fa", "new.fa"):
        p = case_dir / name
        if p.is_symlink() or not p.exists():
            continue
        p.unlink()
        removed.append(name)
    return removed


def _has_outputs(case_dir):
    return any((case_dir / f).exists() for f in
               ("old.cram", "old.cram.crai", "oracle.cram", "oracle.cram.crai",
                "old.json", "oracle.json"))


def input_signatures(old, new, reads):
    """Resolved target path + size + mtime_ns for each intended input."""
    def sig(p):
        rp = Path(p).resolve()
        st = rp.stat()
        return {"path": str(rp), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    return {"old": sig(old), "new": sig(new), "reads": sig(reads)}


def read_inputs_json(case_dir):
    return read_json_file(case_dir / "inputs.json")


def write_inputs_json(case_dir, sigs):
    (case_dir / "inputs.json").write_text(json.dumps(sigs, indent=2) + "\n")


def read_json_file(path):
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return None


def links_match_intended(case_dir, old, new, reads):
    """True when the case-local symlinks already resolve to the intended inputs."""
    want = {"old.fa": Path(old).resolve(), "new.fa": Path(new).resolve(),
            "reads.fastq.gz": Path(reads).resolve()}
    for name, target in want.items():
        lp = case_dir / name
        if not lp.exists():
            return False
        try:
            if lp.resolve() != target:
                return False
        except OSError:
            return False
    return True


def candidate_dirs():
    """Existing case dirs to reuse mappings from.

    Searched newest-source-first: this M7 set also reuses the full-size M6-large /
    M6 case CRAMs, then the M6/data real sets, then the current OUT.
    """
    dirs = []
    seen = set()
    for r in (DATA / "m6large" / "real", DATA / "m6" / "real", DATA / "real",
              DATA / "val" / "real", OUT):
        if not r.is_dir():
            continue
        for d in sorted(p for p in r.iterdir() if p.is_dir()):
            rp = d.resolve()
            if rp in seen:
                continue
            if not ((d / "old.fa").exists() or (d / "new.fa").exists()):
                continue
            seen.add(rp)
            dirs.append(d)
    return dirs


# Mappings are reusable by content, independent of role: a target ``old`` may be
# satisfied by a source ``oracle`` (and vice versa) as long as the reference and reads
# match and the contract holds.
MAPPING_ROLES = (("old", "old.fa"), ("oracle", "new.fa"))


def find_source(target_ref, target_reads, case_dir, candidates, mm2, preset):
    """First ``(source case, source role)`` with a fully-valid, byte-identical mapping."""
    for src in candidates:
        if src.resolve() == case_dir.resolve():
            continue
        sreads = src / "reads.fastq.gz"
        if not sreads.exists():
            continue
        for role, ref_name in MAPPING_ROLES:
            sref = src / ref_name
            if not sref.exists():
                continue
            if not same_content(target_ref, sref) or not same_content(target_reads, sreads):
                continue
            ok, _ = mapping_contract_ok(src, role, sref, sreads, mm2, preset)
            if ok:
                return src, role
    return None, None


def reuse_mapping(mapping, case_dir, ref, reads, mm2, candidates, preset):
    """Hardlink a matching mapping into case_dir and rewrite its JSON.

    ``mapping`` is the *target* name (``old``/``oracle``); the source role may differ.
    Returns True if a source was found and linked.
    """
    src, role = find_source(ref, reads, case_dir, candidates, mm2, preset)
    if src is None:
        return False
    for suffix in (".cram", ".cram.crai"):
        dst = case_dir / f"{mapping}{suffix}"
        dst.unlink(missing_ok=True)
        Path(dst).hardlink_to(src / f"{role}{suffix}")
    data = json.loads((src / f"{role}.json").read_text())
    data["cmd"] = rewrite_cmd(data["cmd"], ref, reads, case_dir / f"{mapping}.cram")
    data["reused_from"] = f"{src}:{role}"
    (case_dir / f"{mapping}.json").write_text(json.dumps(data, indent=2) + "\n")
    print(f"  {case_dir.name}: reused {mapping} from {src}:{role}")
    return True


def ensure_mapping(case_dir, mapping, ref, reads, mm2, candidates, preset):
    """Keep, reuse or invalidate one mapping; returns ``(action, reason)``.

    A compatible existing mapping is kept.  Otherwise a fully-valid reuse source is
    hardlinked; if none exists the mapping's completion artifacts are deleted so
    build_crams rebuilds it.
    """
    ok, reason = mapping_contract_ok(case_dir, mapping, ref, reads, mm2, preset)
    if ok:
        return "kept", reason
    if reuse_mapping(mapping, case_dir, ref, reads, mm2, candidates, preset):
        return "reused", reason
    removed = invalidate_mapping(case_dir, mapping)
    return ("invalidated" if removed else "absent"), reason


def fix_json(root):
    """Repoint existing case JSONs at case-local ref/reads/-o (idempotent)."""
    root = Path(root).resolve()
    if not root.is_dir():
        print(f"error: {root} is not a directory", file=sys.stderr)
        return 2
    fixed = 0
    for case_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        changed = False
        specs = [("old.json", case_dir / "old.fa", case_dir / "old.cram"),
                 ("oracle.json", case_dir / "new.fa", case_dir / "oracle.cram")]
        for jname, ref, out in specs:
            jp = case_dir / jname
            if not jp.exists():
                continue
            try:
                data = json.loads(jp.read_text())
            except Exception as exc:
                print(f"  {case_dir.name}/{jname}: unreadable ({exc})")
                continue
            before = data.get("cmd")
            if not before:
                continue
            after = rewrite_cmd(before, ref, case_dir / "reads.fastq.gz", out)
            if after != before:
                data["cmd"] = after
                jp.write_text(json.dumps(data, indent=2) + "\n")
                changed = True
                print(f"{case_dir.name}/{jname}:")
                print(f"  before: {before}")
                print(f"  after:  {after}")
        fixed += 1 if changed else 0
    print(f"fixed {fixed} case(s) under {root}")
    return 0


def build_cases(sample, cons, final, sample_reads, kind, previously_used):
    """(kind, old, new) triples for a sample under the selected case kind."""
    if kind == KIND_ALL:
        cases = [("cons-to-polished", cons, final)]
        for p in chosen_primaries(sample, total_len(cons), previously_used):
            cases.append((f"{p}-to-cons", ensure_local_assembly(sample, p), cons))
        return cases
    # KIND_RAW: raw primary -> polished assembly.
    return [(f"{p}-to-polished", ensure_local_assembly(sample, p), final)
            for p in chosen_primaries(sample, total_len(cons), previously_used)]


# ---------------------------------------------------------------------------
# M7 rotated assembly (deterministic synthetic dnaA-style reorientation)
# ---------------------------------------------------------------------------
_COMPLEMENT = {}
for _a, _b in zip("ACGTN", "TGCAN"):
    _COMPLEMENT[_a] = _b
    _COMPLEMENT[_a.lower()] = _b.lower()


def _revcomp(seq):
    """Reverse complement; ACGTN (any case) mapped, everything else -> N."""
    return "".join(_COMPLEMENT.get(c, "N") for c in reversed(seq))


def _unrotate(seq, off):
    """Inverse of ``seq[off:] + seq[:off]``."""
    if off == 0:
        return seq
    return seq[len(seq) - off:] + seq[:len(seq) - off]


def _read_fasta(path):
    """``[(name, uppercased sequence)]`` preserving order; name = first header token."""
    contigs = []
    name, parts = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    contigs.append((name, "".join(parts).upper()))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        contigs.append((name, "".join(parts).upper()))
    return contigs


def _circular_contigs(cons):
    """{contig name: bool} from the autocycler consensus headers (``circular=true``)."""
    circ = {}
    with open(cons) as fh:
        for line in fh:
            if line.startswith(">"):
                circ[line[1:].split()[0]] = "circular=true" in line
    return circ


def rotation_plan(polished, cons):
    """Deterministic rotate/reorient plan for the polished assembly.

    Contig order/names follow ``polished``; circularity follows ``cons`` headers.
    ``random.Random(ROTATION_SEED)`` draws one offset per circular contig in file
    order (``off = rng.randint(int(0.1*L), int(0.9*L))``, ``seq = seq[off:] + seq[:off]``);
    the single largest contig (ties -> first) is then reverse-complemented.  Includes a
    unit-style self-check that lengths are preserved and the transform is invertible.
    """
    circ = _circular_contigs(cons)
    rng = random.Random(ROTATION_SEED)
    plan = []
    for name, seq in _read_fasta(polished):
        length = len(seq)
        if circ.get(name, False):
            off = rng.randint(int(0.1 * length), int(0.9 * length))
            out = seq[off:] + seq[:off]
        else:
            off, out = 0, seq
        plan.append({"name": name, "orig": seq, "length": length,
                     "offset": off, "revcomp": False, "out": out})
    largest = max(range(len(plan)), key=lambda i: (plan[i]["length"], -i))
    plan[largest]["revcomp"] = True
    plan[largest]["out"] = _revcomp(plan[largest]["out"])
    for e in plan:  # self-check: invertible and length-preserving
        assert len(e["out"]) == e["length"], e["name"]
        back = _revcomp(e["out"]) if e["revcomp"] else e["out"]
        assert _unrotate(back, e["offset"]) == e["orig"], e["name"]
    return plan


def _sha256_streamed(path):
    """sha256 of a file's bytes, read in bounded chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _rotated_text(plan):
    """The exact FASTA text a plan produces (headers + 80-col wrapped sequences)."""
    parts = []
    for e in plan:
        parts.append(f">{e['name']}\n")
        s = e["out"]
        for i in range(0, len(s), 80):
            parts.append(s[i:i + 80] + "\n")
    return "".join(parts)


def _write_text_atomic(path, text):
    """Atomically write ``text`` via an exclusively-created temp file + ``os.replace``."""
    path = Path(path)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _assert_under_dataset(path, what):
    """Require ``path.resolve()`` to stay inside the resolved --root (DATASET)."""
    root = Path(DATASET).resolve()
    rp = Path(path).resolve()
    try:
        rp.relative_to(root)
    except ValueError:
        raise SystemExit(f"error: {what} {path} resolves to {rp}, outside dataset {root}")
    return rp


def _assert_safe_dest(path):
    """Reject an m7 write destination that is a symlink resolving outside DATASET."""
    if Path(path).is_symlink():
        _assert_under_dataset(path, "destination")


def ensure_polished_rotated(sample, polished, cons):
    """Write ``assemblies/<sample>/polished_rotated.fasta`` (+ JSON); return (path, meta).

    Atomic (.tmp + rename) and idempotent.  The sidecar records a streamed sha256 of the
    polished source; the FASTA is only rewritten when its existing bytes differ from the
    recomputed output, so a byte-identical regeneration (e.g. adding the sha256 to an old
    sidecar) leaves the FASTA mtime - and thus ``inputs.json`` - untouched and only
    refreshes the sidecar.
    """
    out_dir = ASM_OUT / sample
    if out_dir.is_symlink():
        raise SystemExit(f"error: assemblies dir {out_dir} is a symlink; refusing to write")
    _assert_under_dataset(out_dir, "assemblies dir")
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "polished_rotated.fasta"
    jp = out_dir / "polished_rotated.json"
    _assert_safe_dest(out)
    _assert_safe_dest(jp)
    plan = rotation_plan(polished, cons)
    st = polished.resolve().stat()
    meta = {
        "seed": ROTATION_SEED,
        "source": {"path": str(polished.resolve()), "size": st.st_size,
                   "mtime_ns": st.st_mtime_ns, "sha256": _sha256_streamed(polished)},
        "contigs": [{"name": e["name"], "length": e["length"],
                     "offset": e["offset"], "revcomp": e["revcomp"]} for e in plan],
    }
    if out.exists() and jp.exists() and read_json_file(jp) == meta:
        return out, meta
    text = _rotated_text(plan)
    if out.exists() and out.read_bytes() == text.encode():
        _write_text_atomic(jp, json.dumps(meta, indent=2) + "\n")
        print(f"  refreshed {jp}")
        return out, meta
    _write_text_atomic(out, text)
    _write_text_atomic(jp, json.dumps(meta, indent=2) + "\n")
    print(f"  wrote {out}")
    return out, meta


def _materialize(case, platform, preset, sample, old, new, reads, mm2, candidates, rows):
    """Symlink one case's inputs, build the manifest row and reconcile its mappings."""
    case_dir = OUT / case
    if case_dir.is_symlink():
        raise SystemExit(f"error: case dir {case_dir} is a symlink; refusing to populate")
    _assert_under_dataset(case_dir, "case dir")
    case_dir.mkdir(parents=True, exist_ok=True)
    for name in ("inputs.json", "old.json", "oracle.json"):
        _assert_safe_dest(case_dir / name)

    sigs = input_signatures(old, new, reads)
    if _has_outputs(case_dir):
        prior = read_inputs_json(case_dir)
        if prior is not None:
            if prior != sigs:
                removed = invalidate_outputs(case_dir)
                (case_dir / "inputs.json").unlink(missing_ok=True)
                print(f"  {case}: inputs changed; invalidated {removed or 'nothing'}")
        elif links_match_intended(case_dir, old, new, reads):
            write_inputs_json(case_dir, sigs)
            print(f"  {case}: adopted existing outputs (wrote inputs.json)")
        else:
            removed = invalidate_outputs(case_dir)
            print(f"  {case}: inputs mismatch; invalidated {removed or 'nothing'}")

    link(old, case_dir / "old.fa")
    link(new, case_dir / "new.fa")
    link(reads, case_dir / "reads.fastq.gz")
    rows.append([case, platform, preset, sample, old, new, reads])

    actions = {}
    for mapping, ref in (("old", case_dir / "old.fa"), ("oracle", case_dir / "new.fa")):
        action, reason = ensure_mapping(case_dir, mapping, ref,
                                        case_dir / "reads.fastq.gz",
                                        mm2, candidates, preset)
        actions[mapping] = action
        if action == "invalidated":
            print(f"  {case}: {mapping} {reason}; invalidated")

    if _has_outputs(case_dir):
        write_inputs_json(case_dir, sigs)
    return actions


def build_m7(mm2, candidates):
    """Build the M7 "real use case" set under OUT (full-size reads)."""
    OUT.mkdir(parents=True, exist_ok=True)
    ASM_OUT.mkdir(parents=True, exist_ok=True)
    rows, notes, summary, rotations = [], [], [], {}

    def add(case, platform, preset, sample, old, new, reads):
        actions = _materialize(case, platform, preset, sample, old, new, reads,
                               mm2, candidates, rows)
        summary.append((case, actions["old"], actions["oracle"]))

    def note(sample, reads):
        n, bases = read_stats(reads)
        notes.append((sample, "n", None, gz_bytes(reads), n, bases))
        print(f"{sample}: reads={n} bases={bases}")

    for s in ONT_M7_SAMPLES:
        d = SRC / s
        cons = d / "3_autocycler" / "consensus_assembly.fasta"
        polished = d / "3_autocycler" / "consensus_assembly_polished.fasta"
        pbfinal = (PACBIO_SRC / f"{s}-p1-1.1" / "3_autocycler"
                   / "consensus_assembly.finalized.fasta")
        reads = d / "0_reads" / f"{s}.fastq.gz"
        assert polished.exists(), f"missing polished {polished}"
        assert pbfinal.exists(), f"missing pbfinal {pbfinal}"
        rotated, meta = ensure_polished_rotated(s, polished, cons)
        rotations[s] = meta
        add(f"ont_{s}__cons-to-polished", "ont", "lr:hq", s, cons, polished, reads)
        add(f"ont_{s}__polished-to-rotated", "ont", "lr:hq", s, polished, rotated, reads)
        add(f"ont_{s}__polished-to-pbfinal", "ont", "lr:hq", s, polished, pbfinal, reads)
        note(s, reads)

    for s in PACBIO_M7_SAMPLES:
        sample = f"{s}-p1-1.1"
        d = PACBIO_SRC / sample
        reads_files = sorted((d / "0_reads").glob("*.fastq.gz"))
        assert len(reads_files) == 1, \
            f"{sample}: expected exactly 1 read file, got {len(reads_files)}"
        old = d / "3_autocycler" / "consensus_assembly.fasta"
        new = d / "3_autocycler" / "consensus_assembly.finalized.fasta"
        assert old.exists() and new.exists(), f"{sample}: missing assembly"
        add(f"pacbio_{sample}__cons-to-finalized", "pacbio", "map-hifi", sample,
            old, new, reads_files[0])
        note(sample, reads_files[0])

    _assert_safe_dest(OUT / "manifest.tsv")
    with open(OUT / "manifest.tsv", "w") as fh:
        fh.write("case\tplatform\tpreset\tsample\told\tnew\treads\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    print(f"{len(rows)} cases -> {OUT / 'manifest.tsv'}")

    width = max(len(c) for c, _, _ in summary)
    print(f"{'case':{width}}  {'old':11}  {'oracle':11}")
    for case, oa, oo in summary:
        print(f"{case:{width}}  {oa:11}  {oo:11}")

    write_readme_m7(notes, rotations)
    return 0


def write_readme_m7(notes, rotations):
    lines = [
        "# M7 real use-case dataset", "",
        "Full-size reads (no downsampling) for 4 case kinds:", "",
        "- `ont_<S>__cons-to-polished`: old = unpolished `consensus_assembly.fasta`, "
        "new = `consensus_assembly_polished.fasta` (ONT `lr:hq`).",
        "- `ont_<S>__polished-to-rotated`: old = polished, new = a deterministic "
        "`polished_rotated.fasta` (synthetic dnaA-style reorientation).",
        "- `ont_<S>__polished-to-pbfinal`: old = ONT polished, new = the PacBio "
        "`<S>-p1-1.1` `consensus_assembly.finalized.fasta` (same ONT reads).",
        "- `pacbio_<S>-p1-1.1__cons-to-finalized`: old = PacBio `consensus_assembly.fasta`,"
        " new = `consensus_assembly.finalized.fasta` (PacBio `map-hifi`).",
        "",
        "## Rotation rule (polished-to-rotated)", "",
        f"Seed `{ROTATION_SEED}`. Contig order/names follow the polished FASTA; circularity "
        "comes from the `consensus_assembly.fasta` headers. For each circular contig in "
        "file order, `off = rng.randint(int(0.1*L), int(0.9*L))` and the sequence becomes "
        "`seq[off:] + seq[:off]`; the single largest contig (ties -> first) is then "
        "reverse-complemented. Sequences are uppercased and wrapped at 80 columns. "
        "Per-contig parameters live in `assemblies/<S>/polished_rotated.json`.", "",
        "| sample | contig | length | offset | revcomp |",
        "|---|---|---|---|---|",
    ]
    for s in ONT_M7_SAMPLES:
        for c in rotations[s]["contigs"]:
            lines.append(f"| {s} | {c['name']} | {c['length']} | {c['offset']} | "
                         f"{'yes' if c['revcomp'] else 'no'} |")
    lines += [
        "",
        "## PacBio selection rule", "",
        "The 6 largest PacBio read files in `$CRAMHOOK_DATA_PACBIO` (by gzipped "
        "size) whose sample is not already used by `data/real` or `data/val/real`; "
        "`data-2025-12-03` reads are not available locally. Chosen: "
        + ", ".join(PACBIO_M7_SAMPLES) + ".",
        "",
        "## Sample reads", "",
        "| sample | downsampled | gz fraction | reads gz (MB) | reads | bases |",
        "|---|---|---|---|---|---|",
    ]
    for sample, down, frac, size, n, bases in notes:
        fs = f"{frac:.6f}" if frac is not None else "—"
        lines.append(f"| {sample} | {down} | {fs} | {size / 1e6:.1f} | {n} | {bases} |")
    lines.append("")
    _assert_safe_dest(DATASET / "README.md")
    (DATASET / "README.md").write_text("\n".join(lines))


def _resolve_inside_repo(raw, what):
    """Resolve a --root/--fix-root path and require it to stay inside the repo."""
    p = Path(raw)
    if not p.is_absolute():
        p = ROOT / p
    rp = p.resolve()
    try:
        rp.relative_to(ROOT)
    except ValueError:
        raise SystemExit(f"error: {what} {rp} is outside the repository {ROOT}")
    return rp


def _validate_samples(raw):
    """Known ONT sample IDs only; no path separators."""
    samples = (SAMPLES if not raw
               else [s.strip() for s in raw.split(",") if s.strip()])
    for s in samples:
        if s not in SAMPLES or "/" in s or "\\" in s:
            raise SystemExit(f"error: unknown or unsafe sample id {s!r}")
    return samples


def main():
    global DATASET, OUT, READS_OUT, ASM_OUT, ASM_SEARCH, SRC, PACBIO_SRC
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fix-json", action="store_true",
                    help="rewrite existing case JSON cmds in place (no population)")
    ap.add_argument("--fix-root", default=None,
                    help="case-dir root for --fix-json (default: <root>/real)")
    ap.add_argument("--root", default="data/m6",
                    help="dataset dir holding real/, reads/, assemblies/ (default: data/m6)")
    ap.add_argument("--kind", choices=[KIND_ALL, KIND_RAW, KIND_M7], default=KIND_ALL,
                    help="case kind to generate (default: all, i.e. the M6 set)")
    ap.add_argument("--samples", default=None,
                    help="comma-separated sample subset (default: all ONT samples)")
    ap.add_argument("--no-downsample", action="store_true",
                    help="always symlink the original full-size reads")
    args = ap.parse_args()

    if args.kind == KIND_M7 and args.samples is not None:
        raise SystemExit("error: --samples is not supported with --kind m7 "
                         "(M7 is always the full 24-case set)")

    DATASET = _resolve_inside_repo(args.root, "--root")
    OUT = DATASET / "real"
    READS_OUT = DATASET / "reads"
    ASM_OUT = DATASET / "assemblies"
    ASM_SEARCH = [ASM_OUT, DATA / "m6" / "assemblies"]

    if args.fix_json:
        fix_root = _resolve_inside_repo(args.fix_root if args.fix_root else OUT,
                                        "--fix-root")
        return fix_json(fix_root)

    SRC = _dataset_root("CRAMHOOK_DATA_ONT")
    if args.kind == KIND_M7:
        PACBIO_SRC = _dataset_root("CRAMHOOK_DATA_PACBIO")

    samples = _validate_samples(args.samples)
    mm2 = subprocess.run(["minimap2", "--version"], capture_output=True, text=True,
                         check=True).stdout.strip()
    previously_used = existing_primaries()
    OUT.mkdir(parents=True, exist_ok=True)
    READS_OUT.mkdir(parents=True, exist_ok=True)
    ASM_OUT.mkdir(parents=True, exist_ok=True)
    candidates = candidate_dirs()

    if args.kind == KIND_M7:
        return build_m7(mm2, candidates)

    rows = []
    notes = []
    for sample in samples:
        d = SRC / sample
        cons = d / "3_autocycler" / "consensus_assembly.fasta"
        final = d / "3_autocycler" / f"consensus_assembly_{FINAL}.fasta"
        reads = d / "0_reads" / f"{sample}.fastq.gz"

        if not args.no_downsample and gz_bytes(reads) > DOWNSAMPLE_MIN_BYTES:
            sample_reads, fraction = downsample(sample, reads)
            down = "y"
        else:
            sample_reads, fraction = reads, None
            down = "n"

        for kind, old, new in build_cases(sample, cons, final, sample_reads,
                                          args.kind, previously_used):
            case = f"{PLATFORM}_{sample}__{kind}"
            case_dir = OUT / case
            case_dir.mkdir(parents=True, exist_ok=True)

            # Compare intended input signatures against the persisted ones before
            # replacing links; any difference invalidates stale outputs.
            sigs = input_signatures(old, new, sample_reads)
            if _has_outputs(case_dir):
                prior = read_inputs_json(case_dir)
                if prior is not None:
                    if prior != sigs:
                        removed = invalidate_outputs(case_dir)
                        (case_dir / "inputs.json").unlink(missing_ok=True)
                        print(f"  {case}: inputs changed; invalidated {removed or 'nothing'}")
                elif links_match_intended(case_dir, old, new, sample_reads):
                    # Existing completed case predating inputs.json: adopt it.
                    write_inputs_json(case_dir, sigs)
                    print(f"  {case}: adopted existing outputs (wrote inputs.json)")
                else:
                    removed = invalidate_outputs(case_dir)
                    print(f"  {case}: inputs mismatch; invalidated {removed or 'nothing'}")

            link(old, case_dir / "old.fa")
            link(new, case_dir / "new.fa")
            link(sample_reads, case_dir / "reads.fastq.gz")
            rows.append([case, PLATFORM, PRESET, sample, old, new, sample_reads])

            # Generic per-mapping contract: keep a valid existing mapping, else reuse a
            # fully-valid source, else invalidate so build_crams rebuilds it.
            for mapping, ref in (("old", case_dir / "old.fa"),
                                 ("oracle", case_dir / "new.fa")):
                action, reason = ensure_mapping(case_dir, mapping, ref,
                                                case_dir / "reads.fastq.gz",
                                                mm2, candidates, PRESET)
                if action == "invalidated":
                    print(f"  {case}: {mapping} {reason}; invalidated")

            # Persist the input signatures once outputs are (re)established.
            if _has_outputs(case_dir):
                write_inputs_json(case_dir, sigs)

        # Per-sample read stats (after any downsampling).
        n, bases = read_stats(sample_reads)
        notes.append((sample, down, fraction, gz_bytes(sample_reads), n, bases))
        print(f"{sample}: {'downsampled' if down == 'y' else 'original'} "
              f"reads={n} bases={bases}")

    with open(OUT / "manifest.tsv", "w") as fh:
        fh.write("case\tplatform\tpreset\tsample\told\tnew\treads\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    print(f"{len(rows)} cases -> {OUT / 'manifest.tsv'}")

    write_readme(notes, args.kind)
    return 0


def write_readme(notes, kind):
    if kind == KIND_RAW:
        title = "raw primary -> polished ONT dataset"
        intro = [
            f"ONT `--kind {kind}` set: for each sample, `<primary>-to-polished` cases",
            "with old = an ungzipped raw primary, new = "
            "`3_autocycler/consensus_assembly_polished.fasta` and full-size reads.",
        ]
    else:
        title = "M6 ONT test dataset"
        intro = [
            "11 ONT samples from `$CRAMHOOK_DATA_ONT/<S>/`, 4 cases each:",
            "`cons-to-polished` and 3 `<primary>-to-cons` cases.  Preset `lr:hq`.",
        ]
    lines = [f"# {title}", "", *intro, "",
             "Sample reads (after any downsampling):", "",
             "| sample | downsampled | gz fraction | reads gz (MB) | reads | bases |",
             "|---|---|---|---|---|---|"]
    for sample, down, frac, size, n, bases in notes:
        fs = f"{frac:.6f}" if frac is not None else "—"
        lines.append(f"| {sample} | {down} | {fs} | {size / 1e6:.1f} | {n} | {bases} |")
    lines.append("")
    (DATASET / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
