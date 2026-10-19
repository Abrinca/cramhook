#!/usr/bin/env python3
"""Feasibility probe: optimistic upper bound on reads translatable old->new.

A read is translatable only if its padded window lies inside one identity block
(maximal exact-match run between old and new, in old coordinates). Taint and the
divergence/tag rules are ignored, so the fits fraction is an upper bound.
Run: pixi run python scripts/probe.py CASE_DIR ... --out probe.tsv
"""
import argparse
import bisect
import re
import subprocess
import sys
from pathlib import Path

import pysam

SC = 4  # soft-clip CIGAR op
CS_TOK = re.compile(r"[:\*+\-~][0-9A-Za-z]*")  # cs ops: :N, *xy, +seq, -seq, ~seq


def run_diff(case, threads):
    """minimap2 -cx asm5 --cs new.fa old.fa -> PAF, cached as CASE_DIR/probe_diff.paf.

    The PAF is written to probe_diff.paf.tmp and renamed only after minimap2 exits 0.
    The cache is reused only when it is newer than both old.fa and new.fa (mtime of
    resolved symlink targets).
    """
    paf = case / "probe_diff.paf"
    old = case / "old.fa"
    new = case / "new.fa"
    if paf.exists():
        pm = paf.stat().st_mtime
        if (Path(old).resolve().stat().st_mtime < pm
                and Path(new).resolve().stat().st_mtime < pm):
            return paf
    tmp = case / "probe_diff.paf.tmp"
    with open(tmp, "w") as out:
        rc = subprocess.run(["minimap2", "-cx", "asm5", "--cs", "--secondary=no",
                             "-t", str(threads), str(new), str(old)],
                            stdout=out).returncode
    if rc != 0:
        raise RuntimeError(f"minimap2 failed (exit {rc}); not writing {paf}")
    tmp.replace(paf)
    return paf


def derive_blocks(paf):
    """Parse primary PAF, walk cs -> (old_contig, old_s, old_e, new_contig, new_s, strand)."""
    blocks, n_aln = [], 0
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        tags = {t[:5]: t[5:] for t in f[12:]}
        if tags.get("tp:A:") != "P" or not tags.get("cs:Z:"):
            continue
        n_aln += 1
        qname, qstart, qend, strand = f[0], int(f[2]), int(f[3]), f[4]
        tname, tcur, qoff = f[5], int(f[7]), 0
        for m in CS_TOK.finditer(tags["cs:Z:"]):
            op, arg = m.group()[0], m.group()[1:]
            if op == ":":                      # match run -> one block
                n = int(arg)
                s = qstart + qoff if strand == "+" else qend - qoff - n
                blocks.append((qname, s, s + n, tname, tcur, strand))
                qoff += n
                tcur += n
            elif op == "*":                    # substitution: advance both by 1
                qoff += 1
                tcur += 1
            elif op == "+":                    # insertion in query (old) only
                qoff += len(arg)
            elif op == "-":                    # deletion: bases in target (new) only
                tcur += len(arg)
            # "~": intron, never emitted for asm5 -> ignored
    return blocks, n_aln


def build_index(blocks):
    """contig -> (sorted starts, prefix-max ends); handles overlapping blocks."""
    by = {}
    for c, s, e, *_ in blocks:
        by.setdefault(c, []).append((s, e))
    index = {}
    for c, iv in by.items():
        iv.sort()
        starts = [s for s, _ in iv]
        ends = [e for _, e in iv]
        for k in range(1, len(ends)):
            ends[k] = max(ends[k], ends[k - 1])
        index[c] = (starts, ends)
    return index


def fits(index, contig, s, e):
    """Some block on contig contains [s, e) entirely (optimistic over overlaps)."""
    if e <= s or contig not in index:
        return False
    starts, maxend = index[contig]
    k = bisect.bisect_right(starts, s) - 1
    return k >= 0 and maxend[k] >= e


def block_stats(blocks, old_lengths):
    """(n_blocks, block N50, fraction of old bases covered by the block union)."""
    lens = sorted((e - s for _, s, e, *_ in blocks), reverse=True)
    total = sum(lens)
    n50 = acc = 0
    for L in lens:
        acc += L
        if acc * 2 >= total:
            n50 = L
            break
    by = {}
    for c, s, e, *_ in blocks:
        by.setdefault(c, []).append((s, e))
    covered = 0
    for iv in by.values():
        iv.sort()
        cs = ce = -1
        for s, e in iv:
            if s > ce:
                covered += ce - cs if ce > cs else 0
                cs, ce = s, e
            else:
                ce = max(ce, e)
        covered += ce - cs if ce > cs else 0
    total_len = sum(old_lengths.values())
    return len(blocks), n50, covered / total_len if total_len else 0.0


def stream_case(case, index, old_lengths, args):
    """Stream old.cram, classify primary records, count window fits."""
    c = dict(n=0, unmapped=0, has_sa=0, low_mapq=0, long_clip=0,
             candidates=0, fits=0, not_fits=0)
    with pysam.AlignmentFile(str(case / "old.cram"), "rc",
                             reference_filename=str(case / "old.fa"),
                             threads=args.threads) as cram:
        for r in cram:
            if r.is_secondary or r.is_supplementary:
                continue
            c["n"] += 1
            if r.is_unmapped:
                c["unmapped"] += 1
            elif r.has_tag("SA"):
                c["has_sa"] += 1
            elif r.mapping_quality < args.min_mapq:
                c["low_mapq"] += 1
            else:
                cig = r.cigartuples or ()
                left = cig[0][1] if cig and cig[0][0] == SC else 0
                right = cig[-1][1] if cig and cig[-1][0] == SC else 0
                if left + right > args.max_clip:
                    c["long_clip"] += 1
                    continue
                c["candidates"] += 1
                ws = max(0, r.reference_start - left - args.margin)
                we = min(old_lengths.get(r.reference_name, 0),
                         r.reference_end + right + args.margin)
                c["fits" if fits(index, r.reference_name, ws, we) else "not_fits"] += 1
    return c


def process(case, args):
    blocks, n_aln = derive_blocks(run_diff(case, args.threads))
    with pysam.FastaFile(str(case / "old.fa")) as old:
        old_lengths = dict(zip(old.references, old.lengths))
    n_blocks, n50, frac = block_stats(blocks, old_lengths)
    c = stream_case(case, build_index(blocks), old_lengths, args)
    n = c["n"] or 1
    return {"case": case.name, "n_primary_reads": c["n"], "unmapped": c["unmapped"],
            "has_sa": c["has_sa"], "low_mapq": c["low_mapq"], "long_clip": c["long_clip"],
            "candidates": c["candidates"], "fits": c["fits"],
            "frac_fits_of_all": c["fits"] / n,
            "frac_readlevel_deferred": (c["n"] - c["candidates"]) / n,
            "n_aln": n_aln, "n_blocks": n_blocks, "block_n50": n50,
            "frac_old_in_blocks": frac}


COLS = ["case", "n_primary_reads", "unmapped", "has_sa", "low_mapq", "long_clip",
        "candidates", "fits", "frac_fits_of_all", "frac_readlevel_deferred",
        "n_aln", "n_blocks", "block_n50", "frac_old_in_blocks"]


def fmt(v):
    return f"{v:.4f}" if isinstance(v, float) else str(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cases", nargs="+", type=Path)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--margin", type=int, default=500)
    ap.add_argument("--max-clip", type=int, default=200)
    ap.add_argument("--min-mapq", type=int, default=60)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--debug-blocks", action="store_true")
    args = ap.parse_args()

    if args.debug_blocks:
        for case in args.cases:
            for b in sorted(derive_blocks(run_diff(case, args.threads))[0]):
                print(f"{case.name}\t{b[0]}\t{b[1]}\t{b[2]}\t{b[3]}\t{b[4]}\t{b[5]}")
        return

    rows = []
    for case in args.cases:
        missing = [n for n in ("old.fa", "new.fa", "old.cram") if not (case / n).exists()]
        if missing:
            sys.exit(f"{case}: missing {', '.join(missing)}")
        rows.append(process(case, args))

    with open(args.out, "w") as out:
        out.write("\t".join(COLS) + "\n")
        for row in rows:
            out.write("\t".join(fmt(row[c]) for c in COLS) + "\n")
    widths = [max(len(c), *(len(fmt(r[c])) for r in rows)) for c in COLS]
    print("  ".join(c.rjust(w) for c, w in zip(COLS, widths)))
    for row in rows:
        print("  ".join(fmt(row[c]).rjust(w) for c, w in zip(COLS, widths)))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
