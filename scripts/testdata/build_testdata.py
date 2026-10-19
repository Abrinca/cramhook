#!/usr/bin/env python3
"""Build the public test-data tarball from the local FAM29080 sources.

Deterministic: fixed subsampling seed and fractions, a fixed minimap2 thread
count (it appears in the @PG line), relative paths in every command line and a
tarball with sorted members, zero mtimes/owners and a zero gzip mtime.  Same
sources + same tool versions give the same sha256.

Writes ``dist/<name>.tar.gz`` (+ ``.sha256``), stores the sha256 in
``testdata/manifest.json`` and unpacks the tarball into ``testdata/data/`` (as
``pixi run fetch-testdata`` would).  Staging goes to ``scratch/testdata/build``.

Usage::

  pixi run build-testdata [--threads N] [--keep-staging]
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (ROOT, load_manifest, save_manifest, sha256_file,  # noqa: E402
                    unpack)

# Sources (local, not in the repo).  Sample FAM29080, approved to be public.
ONT_DIR = ROOT / "data/v2/real/ont_FAM29080__cons-to-polished"
PB_DIR = ROOT / "scratch/pacbio-phage-ds/FAM29080-p1-1.1"
SOURCES = {
    "ont_reads": ONT_DIR / "reads.fastq.gz",
    "ont_cons": ONT_DIR / "old.fa",          # autocycler consensus
    "ont_polished": ONT_DIR / "new.fa",      # polished; lacks the phage contig
    "hifi_reads": PB_DIR / "m84151_260508_204155_s4.hifi_reads.bc2048--bc2048.hifi_reads.fastq.gz",
    "hifi_asm": PB_DIR / "without-phage.fasta",
    # the 20,262 bp phage (contig 7 of the PacBio-finalized assembly)
    "phage": ROOT / "scratch/bait/pbds/FAM29080/bait.fa",
}

# Subsampling: ~25x of the 2.08 Mb genome (52 Mb of reads).  Full sets:
# ONT 668.3 Mb (321x), HiFi 430.0 Mb (207x).
SEED = 29080
FRACTION = {"ont": 0.078, "hifi": 0.121}
PRESET = {"ont": "lr:hq", "hifi": "map-hifi"}
# minimap2 -t appears in the @PG CL line; keep it fixed for a stable tarball.
MAP_THREADS = 8


def sh(cmd, cwd):
    shown = cwd.relative_to(ROOT) if cwd.is_relative_to(ROOT) else cwd
    print(f"$ (cd {shown} && {cmd})", flush=True)
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], cwd=cwd, check=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=8,
                    help="threads for bgzip and samtools (not minimap2)")
    ap.add_argument("--staging", default=str(ROOT / "scratch/testdata/build"))
    ap.add_argument("--keep-staging", action="store_true")
    args = ap.parse_args(argv)

    man = load_manifest()
    name = man["name"]
    for k, p in SOURCES.items():
        if not p.exists():
            sys.exit(f"missing source {k}: {p}")

    staging = Path(args.staging).resolve()
    top = staging / name
    if top.exists():
        shutil.rmtree(top)
    (top / "ont").mkdir(parents=True)
    (top / "hifi").mkdir(parents=True)

    def copy_fa(src, dst):
        shutil.copyfile(src, dst)
        sh(f"samtools faidx {dst.name}", dst.parent)

    copy_fa(SOURCES["ont_cons"], top / "ont/cons.fa")
    copy_fa(SOURCES["ont_polished"], top / "ont/polished.fa")
    copy_fa(SOURCES["hifi_asm"], top / "hifi/asm.fa")
    copy_fa(SOURCES["phage"], top / "phage.fa")

    for plat, ref in (("ont", "cons.fa"), ("hifi", "asm.fa")):
        d = top / plat
        src = SOURCES[f"{plat}_reads"].resolve()
        sh(f"seqtk sample -s {SEED} {src} {FRACTION[plat]} "
           f"| bgzip -l 9 -@ {args.threads} > reads.fastq.gz", d)
        cram = ref.replace(".fa", ".cram")
        # exactly the input contract command shape (docs/lift.md)
        sh(f"minimap2 -ax {PRESET[plat]} -t {MAP_THREADS} --secondary=no -Y "
           f"{ref} reads.fastq.gz | samtools sort -O CRAM --reference {ref} "
           f"-o {cram}", d)
        # htslib writes @SQ UR:<absolute reference path>; drop it (it leaks the
        # local path and makes the sha depend on the checkout).  -P: no extra
        # @PG, so the header stays the input-contract minimap2 + samtools sort.
        sh(f"samtools reheader -P -c \"sed 's/\\tUR:[^\\t]*//'\" {cram} > noUR.cram "
           f"&& mv noUR.cram {cram}", d)
        sh(f"samtools index {cram}", d)

    check_no_local_paths(top)

    tarball = ROOT / "dist" / f"{name}.tar.gz"
    tarball.parent.mkdir(parents=True, exist_ok=True)
    make_tarball(top, tarball)
    digest = sha256_file(tarball)
    (tarball.parent / (tarball.name + ".sha256")).write_text(
        f"{digest}  {tarball.name}\n")
    man["sha256"] = digest
    man["size"] = tarball.stat().st_size
    save_manifest(man)
    print(f"wrote {tarball.relative_to(ROOT)} "
          f"({tarball.stat().st_size / 1e6:.1f} MB) sha256 {digest}")
    print("updated testdata/manifest.json")

    unpack(tarball, digest)
    if not args.keep_staging:
        shutil.rmtree(top)
    return 0


def check_no_local_paths(top):
    """Abort if a header line of a packed file names a local path or the user."""
    import getpass
    import gzip

    needles = ["/home", "/tmp", str(ROOT), getpass.getuser()]

    def header_lines(p):
        if p.suffix == ".cram":
            out = subprocess.run(["samtools", "view", "-H", "--no-PG", str(p)],
                                 check=True, capture_output=True, text=True).stdout
            yield from out.splitlines()
        elif p.suffix == ".fa":
            with open(p) as fh:
                yield from (line for line in fh if line.startswith(">"))
        elif p.name.endswith(".fastq.gz"):
            with gzip.open(p, "rt") as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 0:
                        yield line

    bad = []
    for p in sorted(top.rglob("*")):
        if p.is_file():
            bad += [f"{p.relative_to(top)}: {line.strip()[:160]}"
                    for line in header_lines(p) if any(n in line for n in needles)]
    if bad:
        sys.exit("local paths in packed headers:\n  " + "\n  ".join(bad[:10]))
    print("checked: no local paths in packed headers")


def make_tarball(top, out):
    """A reproducible tar.gz of ``top`` (members sorted, mtime/owner zeroed)."""
    import gzip
    import tarfile

    def norm(ti):
        ti.mtime = 0
        ti.uid = ti.gid = 0
        ti.uname = ti.gname = ""
        ti.mode = 0o755 if ti.isdir() else 0o644
        return ti

    paths = [top] + sorted(p for p in top.rglob("*"))
    tmp = out.with_name(out.name + ".tmp")
    with open(tmp, "wb") as raw:
        # the members are already compressed: level 1 only wraps them
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0,
                           compresslevel=1) as gz:
            with tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tf:
                for p in paths:
                    tf.add(p, arcname=str(p.relative_to(top.parent)),
                           recursive=False, filter=norm)
    tmp.replace(out)


if __name__ == "__main__":
    sys.exit(main())
