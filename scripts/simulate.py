#!/usr/bin/env python3
"""Generate deterministic synthetic old->new assembly scenarios with reads.

Base genome: a 300 kb circular chromosome (chr) and two plasmids (p1 45 kb,
p2 25 kb) with planted repeats: an IS element (1.2 kb) with 4 copies on chr
(one reverse-complemented) and 1 copy on p1, plus an rRNA-like operon (5 kb)
with 3 copies on chr. All repeats are identical (except the revcomp copy).

Scenarios: S1-S12 and the variants (S5c/S6b/S6c/S6w/S10a/S10b) predate the
2026-10-17 addition of the remaining edit classes.  S13..S21 add inversion,
contig reorder, split, join, translocation, a deletion-created competitor and
hand-picked combinations; SR1..SR5 are random 3-5 edit combinations.  New
scenarios are appended to SCENARIOS so the seed indices of the earlier ones
(and therefore their bytes) are unchanged.

Each scenario writes <out>/<ID>/{old.fa,new.fa,truth.json} and, per read model
(hifi, ont), <out>/<ID>/<model>/reads.fastq.gz simulated from new.fa.

Fully deterministic: every scenario/model uses a fixed integer seed with
numpy.random.default_rng. Run with `pixi run python scripts/simulate.py`.
"""

import argparse
import gzip
import json
import os

import numpy as np

# ----------------------------------------------------------------------------
# sequence helpers

_IDX2BASE = np.frombuffer(b"ACGT", dtype=np.uint8)
_B2I = np.zeros(256, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _B2I[_b] = _i
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def random_seq(rng, n):
    return _IDX2BASE[rng.integers(0, 4, size=n)].tobytes().decode()


def write_fasta(path, contigs):
    with open(path, "w") as f:
        for name, seq in contigs.items():
            f.write(">%s circular=true\n" % name)
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


# ----------------------------------------------------------------------------
# coordinate mapping helpers (all intervals 0-based half-open)


def map_repeat(f, rep, L, flip=False):
    """Map one repeat record through position function f -> list of records."""
    s, e = rep["start"], rep["end"]
    s2, e2 = f(s), f(e)
    strand = rep["strand"]
    if flip:
        strand = "-" if strand == "+" else "+"
        pieces = [(min(s2, e2), max(s2, e2))]
    elif s2 <= e2:
        pieces = [(s2, e2)]
    else:  # interval crosses the origin after rotation
        pieces = [(s2, L), (0, e2)]
    out = []
    for i, (a, b) in enumerate(pieces):
        r = dict(rep)
        r["start"], r["end"], r["strand"] = int(a), int(b), strand
        if len(pieces) > 1:
            r["name"] = rep["name"] + ".%d" % (i + 1)
        out.append(r)
    return out


def apply_edits(seq, edits):
    """Apply point edits (snp/ins/del) to seq; return (new_seq, old->new map)."""
    cur = seq
    off = 0
    events = []  # (threshold, delta, inclusive)
    for ed in sorted(edits, key=lambda e: e["pos"]):
        p = ed["pos"]
        if ed["kind"] == "snp":
            i = p + off
            cur = cur[:i] + ed["base"] + cur[i + 1:]
        elif ed["kind"] == "ins":
            L = len(ed["seq"])
            i = p + off
            cur = cur[:i] + ed["seq"] + cur[i:]
            events.append((p, L, False))
            off += L
        elif ed["kind"] == "del":
            L = ed["length"]
            i = p + off
            cur = cur[:i] + cur[i + L:]
            events.append((p + L, -L, True))
            off -= L

    def map_pos(x):
        d = 0
        for t, delta, inc in events:
            if (x >= t) if inc else (x > t):
                d += delta
        return x + d

    return cur, map_pos


def copy_repeats(reps):
    return [dict(r) for r in reps]


def _insert_shift_repeats(reps, at, delta, contig="chr"):
    """Reproject repeat coordinates through inserting ``delta`` bases at ``at``
    on ``contig`` (a repeat straddling ``at`` is split)."""
    out = []
    for r in reps:
        if r["contig"] != contig or r["end"] <= at:
            out.append(dict(r))
        elif r["start"] >= at:
            out.append({**r, "start": r["start"] + delta, "end": r["end"] + delta})
        else:
            out.append({**r, "end": at})
            out.append({**r, "start": at + delta, "end": r["end"] + delta})
    return out


# ----------------------------------------------------------------------------
# base genome

CHR_LEN, P1_LEN, P2_LEN = 300_000, 45_000, 25_000
IS_LEN, RRNA_LEN = 1_200, 5_000
# chr repeat layout: IS at 25k/65k/130k/230k (one revcomp), rRNA at 100k/170k/270k
CHR_IS_POS = [(25_000, "+"), (65_000, "+"), (130_000, "+"), (230_000, "-")]
CHR_RRNA_POS = [100_000, 170_000, 270_000]
P1_IS_POS = 25_000


def build_base(rng):
    seqs = {
        "chr": random_seq(rng, CHR_LEN),
        "p1": random_seq(rng, P1_LEN),
        "p2": random_seq(rng, P2_LEN),
    }
    is_seq = random_seq(rng, IS_LEN)
    rrna_seq = random_seq(rng, RRNA_LEN)
    reps = []

    chr_ = seqs["chr"]
    for i, (pos, strand) in enumerate(CHR_IS_POS):
        seg = is_seq if strand == "+" else revcomp(is_seq)
        chr_ = chr_[:pos] + seg + chr_[pos + IS_LEN:]
        reps.append({"name": "IS%d" % (i + 1), "contig": "chr",
                     "start": pos, "end": pos + IS_LEN, "strand": strand})
    p1 = seqs["p1"]
    p1 = p1[:P1_IS_POS] + is_seq + p1[P1_IS_POS + IS_LEN:]
    reps.append({"name": "IS5", "contig": "p1", "start": P1_IS_POS,
                 "end": P1_IS_POS + IS_LEN, "strand": "+"})
    for i, pos in enumerate(CHR_RRNA_POS):
        chr_ = chr_[:pos] + rrna_seq + chr_[pos + RRNA_LEN:]
        reps.append({"name": "rRNA%d" % (i + 1), "contig": "chr",
                     "start": pos, "end": pos + RRNA_LEN, "strand": "+"})
    seqs["chr"] = chr_
    seqs["p1"] = p1
    return seqs, reps, is_seq


# ----------------------------------------------------------------------------
# scenario builders


def sc_s1(base, reps, is_seq, rng):
    return {"old": dict(base), "new": dict(base),
            "reps_old": copy_repeats(reps), "reps_new": copy_repeats(reps),
            "edits": [], "description": "identical assemblies"}


def sc_s2(base, reps, is_seq, rng):
    k = 123_457
    L = CHR_LEN
    new = dict(base)
    new["chr"] = base["chr"][k:] + base["chr"][:k]
    f = lambda x: (x - k) % L
    rn = [map_repeat(f, r, L)[0] for r in reps if r["contig"] == "chr"]
    rn += [dict(r) for r in reps if r["contig"] != "chr"]
    edits = [{"type": "rotation", "old_contig": "chr", "old_start": 0,
              "old_end": L, "new_contig": "chr", "new_start": 0,
              "new_end": L, "offset": k, "strand": "+"}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "chromosome rotated by %d bp" % k}


def sc_s3(base, reps, is_seq, rng):
    L = P1_LEN
    new = dict(base)
    new["p1"] = revcomp(base["p1"])
    rn = []
    for r in reps:
        if r["contig"] == "p1":
            rn += map_repeat(lambda x: L - x, r, L, flip=True)
        else:
            rn.append(dict(r))
    edits = [{"type": "reverse_complement", "old_contig": "p1",
              "old_start": 0, "old_end": L, "new_contig": "p1",
              "new_start": 0, "new_end": L, "strand": "-"}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "p1 reverse-complemented"}


def sc_s4(base, reps, is_seq, rng):
    new = {"chr": base["chr"], "plasmid_A": base["p1"], "p2": base["p2"]}
    rn = []
    for r in reps:
        r = dict(r)
        if r["contig"] == "p1":
            r["contig"] = "plasmid_A"
        rn.append(r)
    edits = [{"type": "rename", "old_contig": "p1", "old_start": None,
              "old_end": None, "new_contig": "plasmid_A", "new_start": None,
              "new_end": None}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "p1 renamed to plasmid_A"}


def _p3_plain(rng):
    return random_seq(rng, 30_000)


def sc_s5(base, reps, is_seq, rng):
    p3 = _p3_plain(rng)
    new = dict(base)
    new["p3"] = p3
    edits = [{"type": "add_contig", "old_contig": None, "old_start": None,
              "old_end": None, "new_contig": "p3", "new_start": 0,
              "new_end": len(p3)}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": copy_repeats(reps), "edits": edits,
            "description": "unique 30 kb plasmid p3 added"}


def _p3_with_is(rng, is_seq, diverged):
    p3 = _p3_plain(rng)
    if diverged:
        seg = diverge(is_seq, rng, sub=0.08, indel=0.02)
    else:
        seg = is_seq
    off = 10_000
    p3 = p3[:off] + seg + p3[off + IS_LEN:]
    rep = {"name": "IS6", "contig": "p3", "start": off,
           "end": off + len(seg), "strand": "+"}
    if diverged:
        rep["diverged"] = True
        rep["approx_identity"] = 0.9
    return p3, rep


def sc_s6(base, reps, is_seq, rng):
    p3, rep = _p3_with_is(rng, is_seq, diverged=False)
    new = dict(base)
    new["p3"] = p3
    rn = copy_repeats(reps) + [rep]
    edits = [
        {"type": "add_contig", "old_contig": None, "old_start": None,
         "old_end": None, "new_contig": "p3", "new_start": 0,
         "new_end": len(p3)},
        {"type": "repeat_copy_added", "old_contig": "chr", "old_start": 25_000,
         "old_end": 26_200, "new_contig": "p3", "new_start": 10_000,
         "new_end": 11_200, "note": "exact copy of IS1"},
    ]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "p3 added carrying an exact IS element copy"}


def sc_s6b(base, reps, is_seq, rng):
    p3, rep = _p3_with_is(rng, is_seq, diverged=True)
    new = dict(base)
    new["p3"] = p3
    rn = copy_repeats(reps) + [rep]
    edits = [
        {"type": "add_contig", "old_contig": None, "old_start": None,
         "old_end": None, "new_contig": "p3", "new_start": 0,
         "new_end": len(p3)},
        {"type": "repeat_copy_added", "old_contig": "chr", "old_start": 25_000,
         "old_end": 26_200, "new_contig": "p3", "new_start": 10_000,
         "new_end": 10_000 + rep["end"] - rep["start"],
         "note": "~90% identity diverged copy of IS1"},
    ]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "p3 added carrying a ~90% identity IS element copy"}


def sc_s5c(base, reps, is_seq, rng):
    """S5 genomes (identical to S5) plus a set of targeted chr-p3 chimeras."""
    res = sc_s5(base, reps, is_seq, rng)
    res["description"] = res["description"] + " (+ 10% targeted chr-p3 chimeras)"
    res["note"] = ("targeted chimeras: 10% extra reads join a 2-15 kb chr "
                   "fragment to a 30-3000 bp p3 fragment (log-uniform); p3 "
                   "length is recoverable from the read name suffix :len{L}")
    res["targeted_chimeras"] = True
    return res


def sc_s6c(base, reps, is_seq, rng):
    """S6 genomes (identical to S6) plus targeted chr-p3 chimeras; half of the
    p3 fragments deliberately overlap the IS copy carried on p3."""
    res = sc_s6(base, reps, is_seq, rng)
    is6 = next(r for r in res["reps_new"] if r["name"] == "IS6")
    res["p3_is"] = (is6["start"], is6["end"])
    res["description"] = res["description"] + " (+ 10% targeted chr-p3 chimeras)"
    res["note"] = ("targeted chimeras: 10% extra reads join a 2-15 kb chr "
                   "fragment to a 30-3000 bp p3 fragment (log-uniform); half "
                   "of the p3 fragments overlap the IS copy in p3; p3 length "
                   "is recoverable from the read name suffix :len{L}")
    res["targeted_chimeras"] = True
    return res


def sc_s6w(base, reps, is_seq, rng):
    """S6 (exact IS copy on added p3), but the IS copy spans p3's origin.

    p3 is linearised inside the IS: the second half of the IS sits at the very
    start of p3 and the first half at the very end, so the circular arc
    ``[L3-half, half)`` is an exact copy of the chromosomal IS. Reads are
    sampled circularly (``fragment``/``extract`` already wrap), so the exact
    copy is only visible across the origin -- a T2 windowing scheme without
    wrap-around misses it and leaves the chromosomal IS copy untainted.
    """
    L3 = 30_000
    half = IS_LEN // 2  # 600
    p3 = _p3_plain(rng)
    # arc [L3-half, half) == is_seq exactly: is_seq[:half] at the end of p3,
    # is_seq[half:] at the start.
    p3 = is_seq[half:] + p3[half:L3 - half] + is_seq[:half]
    rep_new = [
        {"name": "IS6.1", "contig": "p3", "start": L3 - half, "end": L3, "strand": "+"},
        {"name": "IS6.2", "contig": "p3", "start": 0, "end": half, "strand": "+"},
    ]
    new = dict(base)
    new["p3"] = p3
    rn = copy_repeats(reps) + rep_new
    edits = [
        {"type": "add_contig", "old_contig": None, "old_start": None,
         "old_end": None, "new_contig": "p3", "new_start": 0,
         "new_end": L3},
        {"type": "repeat_copy_added", "old_contig": "chr", "old_start": 25_000,
         "old_end": 26_200, "new_contig": "p3", "new_start": L3 - half,
         "new_end": half, "note": "exact copy of IS1 spanning p3's origin"},
    ]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "p3 added carrying an exact IS element copy spanning p3's origin",
            "note": "the IS copy wraps p3's origin: is_seq[half:] at the start, "
                    "is_seq[:half] at the end"}


def sc_s7(base, reps, is_seq, rng):
    new = {"chr": base["chr"], "p1": base["p1"]}
    rn = [dict(r) for r in reps if r["contig"] != "p2"]
    edits = [{"type": "remove_contig", "old_contig": "p2", "old_start": 0,
              "old_end": P2_LEN, "new_contig": None, "new_start": None,
              "new_end": None}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits, "description": "p2 removed"}


def sc_s8(base, reps, is_seq, rng):
    pos, ln = 200_000, 12_000
    ins = random_seq(rng, ln)
    new = dict(base)
    new["chr"] = base["chr"][:pos] + ins + base["chr"][pos:]
    f = lambda x: x if x < pos else x + ln
    rn = [map_repeat(f, r, len(new["chr"]))[0] for r in reps
          if r["contig"] == "chr"]
    rn += [dict(r) for r in reps if r["contig"] != "chr"]
    edits = [{"type": "insertion", "old_contig": "chr", "old_start": pos,
              "old_end": pos, "new_contig": "chr", "new_start": pos,
              "new_end": pos + ln, "length": ln}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "12 kb unique sequence inserted into chr at %d" % pos}


def sc_s9(base, reps, is_seq, rng):
    a, b = 150_000, 158_000
    ln = b - a
    new = dict(base)
    new["chr"] = base["chr"][:a] + base["chr"][b:]
    f = lambda x: x if x < b else x - ln
    rn = [map_repeat(f, r, len(new["chr"]))[0] for r in reps
          if r["contig"] == "chr"]
    rn += [dict(r) for r in reps if r["contig"] != "chr"]
    edits = [{"type": "deletion", "old_contig": "chr", "old_start": a,
              "old_end": b, "new_contig": "chr", "new_start": a,
              "new_end": a, "length": ln}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "chr[%d:%d] deleted" % (a, b)}


def _collapse(base, reps, rng):
    X, d = 195_000, 4_000
    seg = base["chr"][X:X + d]
    old_chr = base["chr"][:X] + seg + base["chr"][X:]
    old = dict(base)
    old["chr"] = old_chr
    f1 = lambda x: x if x < X else x + d  # base -> old
    reps_old = []
    for r in reps:
        if r["contig"] == "chr":
            reps_old += map_repeat(f1, r, len(old_chr))
        else:
            reps_old.append(dict(r))
    f2 = lambda x: x if x < X + d else x - d  # old -> new
    reps_new = []
    for r in reps_old:
        if r["contig"] == "chr":
            reps_new += map_repeat(f2, r, len(base["chr"]))
        else:
            reps_new.append(dict(r))
    return X, d, old, reps_old, reps_new


def sc_s10a(base, reps, is_seq, rng):
    X, d, old, reps_old, reps_new = _collapse(base, reps, rng)
    edits = [{"type": "collapsed_duplication", "old_contig": "chr",
              "old_start": X, "old_end": X + 2 * d, "new_contig": "chr",
              "new_start": X, "new_end": X + d,
              "copies_old": 2, "copies_new": 1}]
    return {"old": old, "new": dict(base), "reps_old": reps_old,
            "reps_new": reps_new, "edits": edits,
            "description": "4 kb tandem duplication collapsed to one copy"}


def sc_s10b(base, reps, is_seq, rng):
    X, d, old, reps_old, reps_new = _collapse(base, reps, rng)
    # reverse of S10a: old is single copy, new is duplicated
    edits = [{"type": "expanded_duplication", "old_contig": "chr",
              "old_start": X, "old_end": X + d, "new_contig": "chr",
              "new_start": X, "new_end": X + 2 * d,
              "copies_old": 1, "copies_new": 2}]
    return {"old": dict(base), "new": old, "reps_old": copy_repeats(reps),
            "reps_new": reps_old, "edits": edits,
            "description": "4 kb segment tandemly duplicated (one -> two copies)"}


def _pick_positions(rng, L, blocked, n, margin=200, sep=20):
    allowed = np.ones(L, dtype=bool)
    for s, e in blocked:
        allowed[max(0, s - margin):min(L, e + margin)] = False
    idx = np.nonzero(allowed)[0]
    perm = rng.permutation(len(idx))
    chosen = []
    for i in perm:
        p = int(idx[i])
        if all(abs(p - c) >= sep for c in chosen):
            chosen.append(p)
            if len(chosen) == n:
                break
    chosen.sort()
    return chosen


def sc_s11(base, reps, is_seq, rng):
    counts = {"chr": (20, 10), "p1": (6, 3), "p2": (4, 2)}
    new = {}
    reps_new = []
    edits = []
    for contig, (n_snp, n_ind) in counts.items():
        seq = base[contig]
        L = len(seq)
        blocked = [(r["start"], r["end"]) for r in reps
                   if r["contig"] == contig]
        pos = _pick_positions(rng, L, blocked, n_snp + n_ind)
        kinds = ["snp"] * n_snp + ["indel"] * n_ind
        rng.shuffle(kinds)
        ed = []
        for p, kind in zip(pos, kinds):
            if kind == "snp":
                old = seq[p]
                newb = "ACGT"[(("ACGT".index(old)) + int(rng.integers(1, 4))) % 4]
                ed.append({"pos": p, "kind": "snp", "base": newb})
            elif rng.random() < 0.5:
                ln = int(rng.integers(1, 6))
                ed.append({"pos": p, "kind": "ins", "seq": random_seq(rng, ln)})
            else:
                ln = int(rng.integers(1, 6))
                if p + ln > L:
                    ln = 1
                ed.append({"pos": p, "kind": "del", "length": ln})
        cur, map_pos = apply_edits(seq, ed)
        new[contig] = cur
        for r in reps:
            if r["contig"] != contig:
                continue
            rr = dict(r)
            rr["start"], rr["end"] = int(map_pos(r["start"])), int(map_pos(r["end"]))
            reps_new.append(rr)
        for e in ed:
            p = e["pos"]
            if e["kind"] == "snp":
                edits.append({"type": "snp", "old_contig": contig,
                              "old_start": p, "old_end": p + 1,
                              "new_contig": contig, "new_start": int(map_pos(p)),
                              "new_end": int(map_pos(p)) + 1})
            elif e["kind"] == "ins":
                edits.append({"type": "insertion", "old_contig": contig,
                              "old_start": p, "old_end": p,
                              "new_contig": contig, "new_start": int(map_pos(p)),
                              "new_end": int(map_pos(p)) + len(e["seq"])})
            else:
                edits.append({"type": "deletion", "old_contig": contig,
                              "old_start": p, "old_end": p + e["length"],
                              "new_contig": contig, "new_start": int(map_pos(p)),
                              "new_end": int(map_pos(p))})
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": reps_new, "edits": edits,
            "description": "30 SNPs and 15 small indels (polishing-like)"}


def sc_s12(base, reps, is_seq, rng):
    k = 77_777
    ins1, ins3 = random_seq(rng, 800), random_seq(rng, 2_500)
    ed_old = [{"pos": 45_000, "kind": "ins", "seq": ins1},
              {"pos": 145_000, "kind": "del", "length": 1_500},
              {"pos": 250_000, "kind": "ins", "seq": ins3}]
    mid, map_mid = apply_edits(base["chr"], ed_old)
    Lmid = len(mid)
    new_chr = mid[k:] + mid[:k]

    def f(x):
        return (map_mid(x) - k) % Lmid

    reps_new = []
    for r in reps:
        if r["contig"] == "chr":
            reps_new += map_repeat(f, r, Lmid)
        else:
            reps_new.append(dict(r))
    p3 = _p3_plain(rng)
    new = dict(base)
    new["chr"] = new_chr
    new["p3"] = p3
    edits = [
        {"type": "rotation", "old_contig": "chr", "old_start": 0,
         "old_end": CHR_LEN, "new_contig": "chr", "new_start": 0,
         "new_end": Lmid, "offset": k, "strand": "+",
         "note": "new = rotate(old_with_indels, %d)" % k},
        {"type": "insertion", "old_contig": "chr", "old_start": 45_000,
         "old_end": 45_000, "new_contig": "chr", "new_start": int(f(45_000)),
         "new_end": int(f(45_000)) + 800, "length": 800},
        {"type": "deletion", "old_contig": "chr", "old_start": 145_000,
         "old_end": 146_500, "new_contig": "chr", "new_start": int(f(145_000)),
         "new_end": int(f(145_000)), "length": 1_500},
        {"type": "insertion", "old_contig": "chr", "old_start": 250_000,
         "old_end": 250_000, "new_contig": "chr", "new_start": int(f(250_000)),
         "new_end": int(f(250_000)) + 2_500, "length": 2_500},
        {"type": "add_contig", "old_contig": None, "old_start": None,
         "old_end": None, "new_contig": "p3", "new_start": 0,
         "new_end": len(p3)},
    ]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": reps_new, "edits": edits,
            "description": "rotation + unique 30 kb p3 added + 3 medium indels"}


# ----------------------------------------------------------------------------
# added scenarios (user, 2026-10-17): more edit classes and combinations
#
# These builders express old->new as a list of "placed blocks": every block is
# (new_contig, old_contig, old_start, old_end, strand) in new-assembly order.
# A block describes where a slice of an old contig ends up (revcomp for '-').
# ``assemble`` materialises the new contigs and ``map_repeats_blocks``
# reprojects the planted repeats through the same blocks.  Point edits and
# whole-contig rotation/revcomp reuse the existing helpers.


def _invert_rep(r, a, b):
    """Reproject one repeat through an intra-contig inversion of [a, b)."""
    if r["end"] <= a or r["start"] >= b:
        return [dict(r)]
    out = []
    if r["start"] < a:  # part left of the inversion stays put
        out.append({**r, "end": min(r["end"], a)})
    s2, e2 = max(r["start"], a), min(r["end"], b)
    if s2 < e2:  # inverted part: [a + (b-e2), a + (b-s2)), strand flips
        out.append({**r, "start": a + (b - e2), "end": a + (b - s2),
                    "strand": "-" if r["strand"] == "+" else "+"})
    if r["end"] > b:  # part right of the inversion stays put
        out.append({**r, "start": max(r["start"], b)})
    return out


def _placed(blocks):
    """Turn (new_contig, old_contig, old_start, old_end, strand) tuples into
    dicts carrying the assembled ``new_start`` (blocks per new_contig must be
    listed in order and contiguous)."""
    out = []
    off = {}
    for nc, oc, s, e, strand in blocks:
        d = {"new_contig": nc, "old_contig": oc, "old_start": int(s),
             "old_end": int(e), "strand": strand}
        d["new_start"] = off.get(nc, 0)
        off[nc] = d["new_start"] + (int(e) - int(s))
        out.append(d)
    return out


def assemble(contigs, placed):
    """Materialise the new contigs from placed blocks (dict, ordered)."""
    new = {}
    for b in placed:
        seg = contigs[b["old_contig"]][b["old_start"]:b["old_end"]]
        if b["strand"] == "-":
            seg = revcomp(seg)
        new[b["new_contig"]] = new.get(b["new_contig"], "") + seg
    return new


def _back_project_blocks(pos, placed):
    """Map a coordinate in the assembly produced by ``placed`` back to the
    block's ``(old_contig, old_coord, strand)`` source coordinate (``None`` if
    no block covers it)."""
    for b in placed:
        blen = b["old_end"] - b["old_start"]
        if b["new_start"] <= pos < b["new_start"] + blen:
            off = pos - b["new_start"]
            if b["strand"] == "+":
                return b["old_contig"], b["old_start"] + off, "+"
            return b["old_contig"], b["old_end"] - 1 - off, "-"
    return None


def map_repeats_blocks(reps, placed):
    """Reproject repeats through placed blocks; a repeat may be split."""
    out = []
    for r in reps:
        for b in placed:
            if b["old_contig"] != r["contig"]:
                continue
            s = max(r["start"], b["old_start"])
            e = min(r["end"], b["old_end"])
            if s >= e:
                continue
            o1, o2 = s - b["old_start"], e - b["old_start"]
            blen = b["old_end"] - b["old_start"]
            if b["strand"] == "+":
                ns, ne, st = b["new_start"] + o1, b["new_start"] + o2, r["strand"]
            else:
                ns = b["new_start"] + (blen - o2)
                ne = b["new_start"] + (blen - o1)
                st = "-" if r["strand"] == "+" else "+"
            out.append({**r, "contig": b["new_contig"],
                        "start": int(ns), "end": int(ne), "strand": st})
    return out


def _full_blocks(contigs, order=None):
    names = order if order is not None else list(contigs.keys())
    return [(n, n, 0, len(contigs[n]), "+") for n in names]


def _all_contig_blocks(contigs, special):
    """Blocks over all contigs, replacing ``special`` (a callable returning the
    block list for its contig) where given."""
    blocks = []
    for name in contigs:
        if name in special:
            blocks += special[name](contigs)
        else:
            blocks.append((name, name, 0, len(contigs[name]), "+"))
    return blocks


def sc_s13(base, reps, is_seq, rng):
    a, b = 120_000, 150_000  # contains the IS3 copy (130k..131.2k)
    new = dict(base)
    new["chr"] = base["chr"][:a] + revcomp(base["chr"][a:b]) + base["chr"][b:]
    rn = []
    for r in reps:
        if r["contig"] == "chr":
            rn += _invert_rep(r, a, b)
        else:
            rn.append(dict(r))
    edits = [{"type": "inversion", "old_contig": "chr", "old_start": a,
              "old_end": b, "new_contig": "chr", "new_start": a,
              "new_end": b, "strand": "-"}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "chr[%d:%d] inverted (contains IS3)" % (a, b)}


def sc_s13b(base, reps, is_seq, rng):
    # one inversion end falls inside the IS3 copy -> the copy is split.
    a, b = 131_000, 162_000
    new = dict(base)
    new["chr"] = base["chr"][:a] + revcomp(base["chr"][a:b]) + base["chr"][b:]
    rn = []
    for r in reps:
        if r["contig"] == "chr":
            rn += _invert_rep(r, a, b)
        else:
            rn.append(dict(r))
    edits = [{"type": "inversion", "old_contig": "chr", "old_start": a,
              "old_end": b, "new_contig": "chr", "new_start": a,
              "new_end": b, "strand": "-",
              "note": "left end splits the IS3 copy"}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "chr[%d:%d] inverted with an end inside IS3" % (a, b)}


def sc_s14(base, reps, is_seq, rng):
    order = ["p1", "chr", "p2"]
    new = {k: base[k] for k in order}
    edits = [{"type": "reorder", "old_contig": None, "old_start": None,
              "old_end": None, "new_contig": None, "new_start": None,
              "new_end": None, "order": order}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": copy_repeats(reps), "edits": edits,
            "description": "contig order changed in new.fa (p1, chr, p2)"}


def sc_s14b(base, reps, is_seq, rng):
    order = ["plasmid_B", "chr", "p1"]
    new = {"plasmid_B": base["p2"], "chr": base["chr"], "p1": base["p1"]}
    rn = [({**r, "contig": "plasmid_B"} if r["contig"] == "p2" else dict(r))
          for r in reps]
    edits = [
        {"type": "rename", "old_contig": "p2", "old_start": None, "old_end": None,
         "new_contig": "plasmid_B", "new_start": None, "new_end": None},
        {"type": "reorder", "old_contig": None, "old_start": None, "old_end": None,
         "new_contig": None, "new_start": None, "new_end": None, "order": order},
    ]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "p2 renamed plasmid_B and contigs reordered"}


def sc_s15(base, reps, is_seq, rng):
    pos = 180_000  # between rRNA2 (170k) and IS4 (230k)
    blocks = _placed([
        ("chr_A", "chr", 0, pos, "+"),
        ("chr_B", "chr", pos, CHR_LEN, "+"),
        ("p1", "p1", 0, P1_LEN, "+"),
        ("p2", "p2", 0, P2_LEN, "+"),
    ])
    new = assemble(base, blocks)
    rn = map_repeats_blocks(reps, blocks)
    edits = [
        {"type": "split", "old_contig": "chr", "old_start": 0, "old_end": pos,
         "new_contig": "chr_A", "new_start": 0, "new_end": pos},
        {"type": "split", "old_contig": "chr", "old_start": pos,
         "old_end": CHR_LEN, "new_contig": "chr_B", "new_start": 0,
         "new_end": CHR_LEN - pos},
    ]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "chr split into chr_A/chr_B at %d" % pos}


def sc_s16(base, reps, is_seq, rng):
    # p1 and p2 each carry half of a 20 kb segment R; joining them reconstructs
    # R, which already exists intact on chr (the contig-join junction class).
    p1, p2, half = _R_pieces(base)
    old = {"chr": base["chr"], "p1": p1, "p2": p2}
    blocks = _placed([
        ("chr", "chr", 0, CHR_LEN, "+"),
        ("p12", "p1", 0, len(p1), "+"),
        ("p12", "p2", 0, len(p2), "+"),
    ])
    new = assemble(old, blocks)
    rn = map_repeats_blocks(reps, blocks)
    dup_old, dup_new = _R_pieces_records(half, P1_LEN)
    edits = [
        {"type": "join", "old_contig": "p1", "old_start": 0, "old_end": len(p1),
         "new_contig": "p12", "new_start": 0, "new_end": len(p1)},
        {"type": "join", "old_contig": "p2", "old_start": 0, "old_end": len(p2),
         "new_contig": "p12", "new_start": len(p1),
         "new_end": len(p1) + len(p2),
         "note": "the join reconstructs R, duplicating the intact chr copy"},
    ]
    return {"old": old, "new": new, "reps_old": copy_repeats(reps) + dup_old,
            "reps_new": rn + dup_new, "edits": edits,
            "description": "p1 + p2 joined into p12; the join reconstructs a "
                           "20 kb copy that duplicates an intact chr segment "
                           "(contig-join competitor)"}


def sc_s17(base, reps, is_seq, rng):
    a, b, d = 60_000, 70_000, 200_000  # block move within chr (contains IS2)
    blocks = _placed([
        ("chr", "chr", 0, a, "+"),
        ("chr", "chr", b, d, "+"),
        ("chr", "chr", a, b, "+"),
        ("chr", "chr", d, CHR_LEN, "+"),
        ("p1", "p1", 0, P1_LEN, "+"),
        ("p2", "p2", 0, P2_LEN, "+"),
    ])
    new = assemble(base, blocks)
    rn = map_repeats_blocks(reps, blocks)
    edits = [{"type": "translocation", "old_contig": "chr", "old_start": a,
              "old_end": b, "new_contig": "chr", "new_start": d - (b - a),
              "new_end": d, "strand": "+",
              "note": "chr[%d:%d] moved after %d" % (a, b, d)}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "chr[%d:%d] translocated within chr" % (a, b)}


def sc_s17b(base, reps, is_seq, rng):
    a, b = 150_000, 155_000  # unique chr block moved onto p1
    blocks = _placed([
        ("chr", "chr", 0, a, "+"),
        ("chr", "chr", b, CHR_LEN, "+"),
        ("p1", "p1", 0, P1_LEN, "+"),
        ("p1", "chr", a, b, "+"),
        ("p2", "p2", 0, P2_LEN, "+"),
    ])
    new = assemble(base, blocks)
    rn = map_repeats_blocks(reps, blocks)
    edits = [{"type": "translocation", "old_contig": "chr", "old_start": a,
              "old_end": b, "new_contig": "p1", "new_start": P1_LEN,
              "new_end": P1_LEN + (b - a), "strand": "+",
              "note": "chr block appended to p1"}]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "chr[%d:%d] translocated onto p1" % (a, b)}


# --- deletion-created competitor / contig join (the lift_junction class) ----
#
# The junction bug (scratch/manager/lift-junction-taint-spec.md,
# tests/lift_junction.rs) needs an *unchanged* copy that carries the read at a
# high old MAPQ plus an edit that makes a second copy contiguous without adding
# new bases.  A segment long enough to contain whole reads (20 kb, longer than
# a typical HiFi read) is therefore duplicated: one copy (RcopyA) stays intact,
# the other is either fragmented by two separators that the new assembly
# deletes (S18/S21) or split across two old contigs that the new assembly joins
# (S16/S20).  Reads wholly inside the segment are ambiguous (MAPQ 0) in the new
# assembly while lift, missing the junction taint, keeps their old MAPQ 60.

R_SEG_LEN = 20_000
R_SEG_POS = 180_000  # unique in the base chr (between rRNA2 at 170k and IS4)
R_SEP_LEN = 20_000   # long enough that minimap2 cannot chain across the gaps
                     # (s2/s1 < max_s2_ratio), so the old copy is not a
                     # `competitor` deferral and the read is translated.


def _R_seg(base):
    return base["chr"][R_SEG_POS:R_SEG_POS + R_SEG_LEN]


def _fragmented_copy(base, rng, R0=R_SEG_POS, Rlen=R_SEG_LEN, P=250_000,
                     sep_len=R_SEP_LEN):
    """Old chr with a 3-piece fragmented copy of R at P (two separators).

    The pieces are R/3 (~6.7 kb) so minimap2 finds no near-equal secondary for
    a typical read on the old fragmented copy (s2 < 0.5*s1), i.e. the read is
    *not* deferred as `competitor` in old; deleting the separators in new then
    creates a full competitor and the read must be remapped (the junction bug).
    """
    R = base["chr"][R0:R0 + Rlen]
    a, b = Rlen // 3, 2 * Rlen // 3
    sep1, sep2 = random_seq(rng, sep_len), random_seq(rng, sep_len)
    B = R[:a] + sep1 + R[a:b] + sep2 + R[b:Rlen]
    return base["chr"][:P] + B + base["chr"][P:], R


def _fragmented_blocks(P, sep_len, Lold, Rlen=R_SEG_LEN):
    """Blocks that delete the two separators of the fragmented copy."""
    a, b = Rlen // 3, 2 * Rlen // 3
    return _placed([
        ("chr", "chr", 0, P, "+"),
        ("chr", "chr", P, P + a, "+"),
        ("chr", "chr", P + a + sep_len, P + b + sep_len, "+"),
        ("chr", "chr", P + b + 2 * sep_len, P + Rlen + 2 * sep_len, "+"),
        ("chr", "chr", P + Rlen + 2 * sep_len, Lold, "+"),
        ("p1", "p1", 0, P1_LEN, "+"),
        ("p2", "p2", 0, P2_LEN, "+"),
    ])


def _R_records(R0=R_SEG_POS, Rlen=R_SEG_LEN, P=250_000, sep_len=R_SEP_LEN):
    """Repeat records for the two R copies (old and new)."""
    a, b = Rlen // 3, 2 * Rlen // 3
    old = [
        {"name": "RcopyA", "contig": "chr", "start": R0, "end": R0 + Rlen,
         "strand": "+", "note": "intact copy that carried the read"},
        {"name": "RcopyB1", "contig": "chr", "start": P, "end": P + a,
         "strand": "+", "note": "fragmented copy, left third"},
        {"name": "RcopyB2", "contig": "chr", "start": P + a + sep_len,
         "end": P + b + sep_len, "strand": "+",
         "note": "fragmented copy, middle third"},
        {"name": "RcopyB3", "contig": "chr", "start": P + b + 2 * sep_len,
         "end": P + Rlen + 2 * sep_len, "strand": "+",
         "note": "fragmented copy, right third"},
    ]
    new = [
        {"name": "RcopyA", "contig": "chr", "start": R0, "end": R0 + Rlen,
         "strand": "+", "note": "intact copy that carried the read"},
        {"name": "RcopyB", "contig": "chr", "start": P, "end": P + Rlen,
         "strand": "+", "note": "copy joined by deleting two separators"},
    ]
    return old, new


def _R_pieces(base):
    """(p1, p2) with a copy of R split in half across the two contigs."""
    R = _R_seg(base)
    half = R_SEG_LEN // 2
    return base["p1"] + R[:half], R[half:] + base["p2"], half


def _R_pieces_records(half, p1_split, old_chr="chr", new_chr="chr", new_off=0):
    old = [
        {"name": "RcopyA", "contig": old_chr, "start": R_SEG_POS,
         "end": R_SEG_POS + R_SEG_LEN, "strand": "+",
         "note": "intact copy"},
        {"name": "RcopyB1", "contig": "p1", "start": p1_split,
         "end": p1_split + half, "strand": "+", "note": "left half of the split copy"},
        {"name": "RcopyB2", "contig": "p2", "start": 0, "end": half,
         "strand": "+", "note": "right half of the split copy"},
    ]
    new = [
        {"name": "RcopyA", "contig": new_chr, "start": R_SEG_POS - new_off,
         "end": R_SEG_POS - new_off + R_SEG_LEN, "strand": "+",
         "note": "intact copy"},
        {"name": "RcopyB", "contig": "p12", "start": p1_split,
         "end": p1_split + R_SEG_LEN, "strand": "+",
         "note": "copy reconstructed by the contig join"},
    ]
    return old, new


def sc_s18(base, reps, is_seq, rng):
    R0, Rlen, P, sep_len = R_SEG_POS, R_SEG_LEN, 250_000, R_SEP_LEN
    a, b = Rlen // 3, 2 * Rlen // 3
    old_chr, R = _fragmented_copy(base, rng, R0, Rlen, P, sep_len)
    Lold = len(old_chr)
    old = {"chr": old_chr, "p1": base["p1"], "p2": base["p2"]}
    blocks = _fragmented_blocks(P, sep_len, Lold, Rlen)
    new = assemble(old, blocks)
    # Inserting the fragmented copy (Rlen + 2*sep_len bases) at P shifts every
    # chr position >= P; the truth must use the same offset as the sequence.
    reps_ins = _insert_shift_repeats(copy_repeats(reps), P, Rlen + 2 * sep_len)
    rn = map_repeats_blocks(reps_ins, blocks)
    dup_old, dup_new = _R_records(R0, Rlen, P, sep_len)
    reps_old = reps_ins + dup_old
    reps_new = rn + dup_new
    edits = [
        {"type": "deletion", "old_contig": "chr", "old_start": P + a,
         "old_end": P + a + sep_len, "new_contig": "chr",
         "new_start": P + a, "new_end": P + a, "length": sep_len,
         "note": "deletes separator 1 of the fragmented copy"},
        {"type": "deletion", "old_contig": "chr",
         "old_start": P + b + sep_len, "old_end": P + b + 2 * sep_len,
         "new_contig": "chr", "new_start": P + b, "new_end": P + b,
         "length": sep_len, "note": "deletes separator 2 of the fragmented copy"},
    ]
    return {"old": old, "new": new, "reps_old": reps_old, "reps_new": reps_new,
            "edits": edits,
            "description": "deletion joins a fragmented 20 kb copy that "
                           "duplicates an intact chr segment "
                           "(deletion-created competitor)"}


def sc_s19(base, reps, is_seq, rng):
    a, b, k = 120_000, 150_000, 77_777
    L = CHR_LEN
    seq = base["chr"]
    seq = seq[:a] + revcomp(seq[a:b]) + seq[b:]
    seq = seq[k:] + seq[:k]
    new_chr = revcomp(seq)
    cr = [r for r in reps if r["contig"] == "chr"]
    rn = [dict(r) for r in reps if r["contig"] != "chr"]
    cr = [x for r in cr for x in _invert_rep(r, a, b)]
    cr = [x for r in cr for x in map_repeat(lambda p: (p - k) % L, r, L)]
    cr = [x for r in cr for x in map_repeat(lambda p: L - p, r, L, flip=True)]
    rn += cr
    new = dict(base)
    new["chr"] = new_chr
    edits = [
        {"type": "inversion", "old_contig": "chr", "old_start": a, "old_end": b,
         "new_contig": "chr", "new_start": a, "new_end": b, "strand": "-"},
        {"type": "rotation", "old_contig": "chr", "old_start": 0, "old_end": L,
         "new_contig": "chr", "new_start": 0, "new_end": L, "offset": k,
         "strand": "+", "note": "applied after the inversion"},
        {"type": "reverse_complement", "old_contig": "chr", "old_start": 0,
         "old_end": L, "new_contig": "chr", "new_start": 0, "new_end": L,
         "strand": "-", "note": "applied last (whole contig)"},
    ]
    return {"old": dict(base), "new": new, "reps_old": copy_repeats(reps),
            "reps_new": rn, "edits": edits,
            "description": "chr inversion + rotation + whole-contig revcomp"}


def sc_s20(base, reps, is_seq, rng):
    # split chr (outside R) + join p1/p2 (reconstructing R, a duplicate of the
    # intact chr copy) + reorder.
    split_pos = 100_000
    p1, p2, half = _R_pieces(base)
    old = {"chr": base["chr"], "p1": p1, "p2": p2}
    blocks = _placed([
        ("chr_B", "chr", split_pos, CHR_LEN, "+"),
        ("p12", "p1", 0, len(p1), "+"),
        ("p12", "p2", 0, len(p2), "+"),
        ("chr_A", "chr", 0, split_pos, "+"),
    ])
    new = assemble(old, blocks)
    rn = map_repeats_blocks(reps, blocks)
    dup_old, dup_new = _R_pieces_records(half, P1_LEN, "chr", "chr_B", split_pos)
    edits = [
        {"type": "split", "old_contig": "chr", "old_start": 0,
         "old_end": split_pos, "new_contig": "chr_A", "new_start": 0,
         "new_end": split_pos},
        {"type": "split", "old_contig": "chr", "old_start": split_pos,
         "old_end": CHR_LEN, "new_contig": "chr_B", "new_start": 0,
         "new_end": CHR_LEN - split_pos},
        {"type": "join", "old_contig": "p1", "old_start": 0, "old_end": len(p1),
         "new_contig": "p12", "new_start": 0, "new_end": len(p1)},
        {"type": "join", "old_contig": "p2", "old_start": 0, "old_end": len(p2),
         "new_contig": "p12", "new_start": len(p1),
         "new_end": len(p1) + len(p2),
         "note": "the join reconstructs R, duplicating the intact chr copy"},
        {"type": "reorder", "old_contig": None, "old_start": None,
         "old_end": None, "new_contig": None, "new_start": None, "new_end": None,
         "order": ["p12", "chr_B", "chr_A"]},
    ]
    return {"old": old, "new": new, "reps_old": copy_repeats(reps) + dup_old,
            "reps_new": rn + dup_new, "edits": edits,
            "description": "split + join + reorder; the join reconstructs a "
                           "20 kb copy that duplicates an intact chr segment"}


def sc_s21(base, reps, is_seq, rng):
    R0, Rlen, P, sep_len = R_SEG_POS, R_SEG_LEN, 250_000, R_SEP_LEN
    a, b = Rlen // 3, 2 * Rlen // 3
    old_chr, R = _fragmented_copy(base, rng, R0, Rlen, P, sep_len)
    Lold = len(old_chr)
    old = {"chr": old_chr, "p1": base["p1"], "p2": base["p2"]}
    # 1. delete the two separators (deletion-created competitor)
    del_blocks = _fragmented_blocks(P, sep_len, Lold, Rlen)
    nm = assemble(old, del_blocks)
    # The sequence path (kept identical to the pre-fix generator, so the reads
    # and crams are byte-identical): the polishing positions are chosen from the
    # *old* repeat chain.
    reps_nm_seq = map_repeats_blocks(reps, del_blocks)
    dup_old, _ = _R_records(R0, Rlen, P, sep_len)
    # The truth path: inserting the fragmented copy shifts chr repeats >= P, so
    # derive the truth from the same offset as the sequence.
    reps_old = _insert_shift_repeats(copy_repeats(reps), P, Rlen + 2 * sep_len) \
        + dup_old
    # 2. translocation inside chr: move a unique tail block to the front
    Lnm = len(nm["chr"])
    a2, b2 = P + Rlen + 15_000, P + Rlen + 21_000
    tblocks = _placed([
        ("chr", "chr", a2, b2, "+"),
        ("chr", "chr", 0, a2, "+"),
        ("chr", "chr", b2, Lnm, "+"),
        ("p1", "p1", 0, P1_LEN, "+"),
        ("p2", "p2", 0, P2_LEN, "+"),
    ])
    new = assemble(nm, tblocks)
    reps_t_seq = map_repeats_blocks(reps_nm_seq, tblocks)
    # 3. SNPs / small indels (polishing-like) on chr
    seq = new["chr"]
    blocked = [(r["start"], r["end"]) for r in reps_t_seq if r["contig"] == "chr"]
    pos = _pick_positions(rng, len(seq), blocked, 8)
    ed = []
    for p in pos:
        if rng.random() < 0.5:
            oldb = seq[p]
            newb = "ACGT"[(("ACGT".index(oldb)) + int(rng.integers(1, 4))) % 4]
            ed.append({"pos": p, "kind": "snp", "base": newb})
        elif rng.random() < 0.5:
            ed.append({"pos": p, "kind": "ins",
                       "seq": random_seq(rng, int(rng.integers(1, 6)))})
        else:
            ed.append({"pos": p, "kind": "del", "length": int(rng.integers(1, 6))})
    new["chr"], map_pos = apply_edits(seq, ed)
    reps_t = map_repeats_blocks(
        map_repeats_blocks(reps_old, del_blocks), tblocks)
    reps_new = []
    for r in reps_t:
        if r["contig"] == "chr":
            reps_new.append({**r, "start": int(map_pos(r["start"])),
                             "end": int(map_pos(r["end"]))})
        else:
            reps_new.append(r)
    # Back-project the post-translocation coordinates to the original old chr:
    # nm = del_blocks(old), the pre-polish new chr = tblocks(nm).
    _, ta, _ = _back_project_blocks(a2, del_blocks)
    _, tb, _ = _back_project_blocks(b2 - 1, del_blocks)
    edits = [
        {"type": "deletion", "old_contig": "chr", "old_start": P + a,
         "old_end": P + a + sep_len, "new_contig": "chr",
         "new_start": P + a, "new_end": P + a, "length": sep_len},
        {"type": "deletion", "old_contig": "chr",
         "old_start": P + b + sep_len, "old_end": P + b + 2 * sep_len,
         "new_contig": "chr", "new_start": P + b, "new_end": P + b,
         "length": sep_len},
        {"type": "translocation", "old_contig": "chr", "old_start": ta,
         "old_end": tb + 1, "new_contig": "chr", "new_start": 0,
         "new_end": b2 - a2, "strand": "+",
         "note": "unique tail block moved to the front of chr"},
    ]
    # Polishing edits, in the established old/new interval schema and
    # back-projected to original old chr coordinates.
    for e in ed:
        p = e["pos"]  # pre-polish new chr coordinate
        _, nm_pos, _ = _back_project_blocks(p, tblocks)
        oc, op, _ = _back_project_blocks(nm_pos, del_blocks)
        if e["kind"] == "snp":
            old_end, new_end, extra = op + 1, p + 1, {"base": e["base"]}
        elif e["kind"] == "ins":
            old_end, new_end, extra = op, p + len(e["seq"]), {"seq": e["seq"]}
        else:
            old_end, new_end, extra = op + e["length"], p, {"length": e["length"]}
        edits.append({
            "type": {"snp": "snp", "ins": "insertion", "del": "deletion"}[e["kind"]],
            "old_contig": oc, "old_start": op, "old_end": old_end,
            "new_contig": "chr", "new_start": p, "new_end": new_end, **extra,
        })
    return {"old": old, "new": new, "reps_old": reps_old, "reps_new": reps_new,
            "edits": edits,
            "description": "translocation + deletion-created competitor + "
                           "SNPs/small indels"}


# --- random combinations ----------------------------------------------------


class _Asm:
    """Mutable assembly state used to compose random edit sequences."""

    def __init__(self, contigs, reps):
        self.c = dict(contigs)
        self.reps = [dict(r) for r in reps]
        self.edits = []
        self.used = {}

    def _used_ok(self, contig, s, e, margin):
        for (us, ue) in self.used.get(contig, []):
            if s < ue + margin and us - margin < e:
                return False
        return True

    def free_region(self, contig, length, rng, margin=500):
        L = len(self.c[contig])
        if L < length + 2 * margin:
            return None
        blocked = [(r["start"], r["end"]) for r in self.reps
                   if r["contig"] == contig]
        for _ in range(50):
            s = int(rng.integers(0, L - length + 1))
            e = s + length
            if any(s < be + margin and bs - margin < e for bs, be in blocked):
                continue
            if not self._used_ok(contig, s, e, margin):
                continue
            return s, e
        return None

    def mark(self, contig, s, e):
        self.used.setdefault(contig, []).append((s, e))

    def rotate(self, contig, k):
        L = len(self.c[contig])
        self.c[contig] = self.c[contig][k:] + self.c[contig][:k]
        self.reps = [y for r in self.reps
                     for y in (map_repeat(lambda x: (x - k) % L, r, L)
                               if r["contig"] == contig else [r])]
        self.edits.append({"type": "rotation", "old_contig": contig,
                           "old_start": 0, "old_end": L, "new_contig": contig,
                           "new_start": 0, "new_end": L, "offset": k,
                           "strand": "+"})

    def revcomp(self, contig):
        L = len(self.c[contig])
        self.c[contig] = revcomp(self.c[contig])
        self.reps = [y for r in self.reps
                     for y in (map_repeat(lambda x: L - x, r, L, flip=True)
                               if r["contig"] == contig else [r])]
        self.edits.append({"type": "reverse_complement", "old_contig": contig,
                           "old_start": 0, "old_end": L, "new_contig": contig,
                           "new_start": 0, "new_end": L, "strand": "-"})

    def rename(self, old, new):
        self.c = {new if k == old else k: v for k, v in self.c.items()}
        for r in self.reps:
            if r["contig"] == old:
                r["contig"] = new
        self.edits.append({"type": "rename", "old_contig": old, "old_start": None,
                           "old_end": None, "new_contig": new, "new_start": None,
                           "new_end": None})

    def reorder(self, order):
        self.c = {k: self.c[k] for k in order}
        self.edits.append({"type": "reorder", "old_contig": None, "old_start": None,
                           "old_end": None, "new_contig": None, "new_start": None,
                           "new_end": None, "order": order})

    def invert(self, contig, a, b):
        seq = self.c[contig]
        self.c[contig] = seq[:a] + revcomp(seq[a:b]) + seq[b:]
        self.reps = [y for r in self.reps
                     for y in (_invert_rep(r, a, b) if r["contig"] == contig
                               else [r])]
        self.edits.append({"type": "inversion", "old_contig": contig,
                           "old_start": a, "old_end": b, "new_contig": contig,
                           "new_start": a, "new_end": b, "strand": "-"})

    def translocate(self, contig, a, b, d):
        L = len(self.c[contig])
        if a < b <= d:
            pieces = [(0, a), (b, d), (a, b), (d, L)]
        elif d <= a < b:
            pieces = [(0, d), (a, b), (d, a), (b, L)]
        else:
            return False
        special = {contig: (lambda c, pieces=pieces:
                            [(contig, contig, s, e, "+")
                             for (s, e) in pieces if s < e])}
        blocks = _all_contig_blocks(self.c, special)
        placed = _placed(blocks)
        self.c = assemble(self.c, placed)
        self.reps = map_repeats_blocks(self.reps, placed)
        self.edits.append({"type": "translocation", "old_contig": contig,
                           "old_start": a, "old_end": b, "new_contig": contig,
                           "new_start": d - (b - a) if d > b else d,
                           "new_end": d if d > b else d + (b - a),
                           "strand": "+"})
        return True

    def insert(self, contig, pos, length, rng):
        ins = random_seq(rng, length)
        self.c[contig] = self.c[contig][:pos] + ins + self.c[contig][pos:]
        for r in self.reps:
            if r["contig"] == contig and r["start"] >= pos:
                r["start"] += length
                r["end"] += length
        self.edits.append({"type": "insertion", "old_contig": contig,
                           "old_start": pos, "old_end": pos,
                           "new_contig": contig, "new_start": pos,
                           "new_end": pos + length, "length": length})

    def delete(self, contig, a, b):
        length = b - a
        self.c[contig] = self.c[contig][:a] + self.c[contig][b:]
        out = []
        for r in self.reps:
            if r["contig"] != contig:
                out.append(r)
            elif r["end"] <= a:
                out.append(r)
            elif r["start"] >= b:
                out.append({**r, "start": r["start"] - length,
                            "end": r["end"] - length})
        self.reps = out
        self.edits.append({"type": "deletion", "old_contig": contig,
                           "old_start": a, "old_end": b, "new_contig": contig,
                           "new_start": a, "new_end": a, "length": length})

    def expand(self, contig, pos, d):
        seq = self.c[contig]
        self.c[contig] = seq[:pos] + seq[pos:pos + d] + seq[pos:]
        for r in self.reps:
            if r["contig"] == contig and r["start"] >= pos:
                r["start"] += d
                r["end"] += d
        self.edits.append({"type": "expanded_duplication", "old_contig": contig,
                           "old_start": pos, "old_end": pos + d,
                           "new_contig": contig, "new_start": pos,
                           "new_end": pos + 2 * d, "copies_old": 1, "copies_new": 2})

    def snps(self, contig, n, rng):
        seq = self.c[contig]
        blocked = [(r["start"], r["end"]) for r in self.reps
                   if r["contig"] == contig]
        pos = _pick_positions(rng, len(seq), blocked, n)
        ed = []
        for p in pos:
            x = rng.random()
            if x < 0.5:
                oldb = seq[p]
                newb = "ACGT"[(("ACGT".index(oldb)) + int(rng.integers(1, 4))) % 4]
                ed.append({"pos": p, "kind": "snp", "base": newb})
            elif x < 0.75:
                ed.append({"pos": p, "kind": "ins",
                           "seq": random_seq(rng, int(rng.integers(1, 6)))})
            else:
                ed.append({"pos": p, "kind": "del",
                           "length": int(rng.integers(1, 6))})
        cur, map_pos = apply_edits(seq, ed)
        self.c[contig] = cur
        for r in self.reps:
            if r["contig"] == contig:
                r["start"], r["end"] = int(map_pos(r["start"])), int(map_pos(r["end"]))
        for e in ed:
            p = e["pos"]
            t = {"snp": "snp", "ins": "insertion", "del": "deletion"}[e["kind"]]
            self.edits.append({"type": t, "old_contig": contig, "old_start": p,
                               "old_end": p + (1 if e["kind"] == "snp" else 0),
                               "new_contig": contig,
                               "new_start": int(map_pos(p)),
                               "new_end": int(map_pos(p))})

    def add_contig(self, name, seq):
        self.c[name] = seq
        self.edits.append({"type": "add_contig", "old_contig": None,
                           "old_start": None, "old_end": None, "new_contig": name,
                           "new_start": 0, "new_end": len(seq)})

    def remove_contig(self, name):
        del self.c[name]
        self.reps = [r for r in self.reps if r["contig"] != name]
        self.edits.append({"type": "remove_contig", "old_contig": name,
                           "old_start": 0, "old_end": None, "new_contig": None,
                           "new_start": None, "new_end": None})

    def split(self, contig, pos, name_a, name_b):
        seq = self.c.pop(contig)
        self.c[name_a] = seq[:pos]
        self.c[name_b] = seq[pos:]
        out = []
        for r in self.reps:
            if r["contig"] != contig:
                out.append(r)
            elif r["end"] <= pos:
                out.append({**r, "contig": name_a})
            elif r["start"] >= pos:
                out.append({**r, "contig": name_b, "start": r["start"] - pos,
                            "end": r["end"] - pos})
        self.reps = out
        self.edits.append({"type": "split", "old_contig": contig, "old_start": 0,
                           "old_end": pos, "new_contig": name_a, "new_start": 0,
                           "new_end": pos})

    def join(self, c1, c2, name):
        s1 = self.c.pop(c1)
        s2 = self.c.pop(c2)
        self.c[name] = s1 + s2
        out = []
        for r in self.reps:
            if r["contig"] == c1:
                out.append({**r, "contig": name})
            elif r["contig"] == c2:
                out.append({**r, "contig": name, "start": r["start"] + len(s1),
                            "end": r["end"] + len(s1)})
            else:
                out.append(r)
        self.reps = out
        self.edits.append({"type": "join", "old_contig": c1, "old_start": 0,
                           "old_end": len(s1), "new_contig": name, "new_start": 0,
                           "new_end": len(s1)})


# Random op classes; each returns True when it applied.
_RANDOM_OPS = ["rotation", "revcomp", "rename", "reorder", "inversion",
               "split", "join", "translocation", "insertion", "deletion",
               "expansion", "snps", "add_plasmid", "remove_plasmid"]


def _random_op(asm, op, rng, counter):
    names = list(asm.c.keys())
    contigs = [n for n in names if len(asm.c[n]) >= 20_000]
    if op == "rotation" and contigs:
        n = contigs[int(rng.integers(0, len(contigs)))]
        asm.rotate(n, int(rng.integers(1000, len(asm.c[n]) - 1000)))
        return True
    if op == "revcomp" and contigs:
        asm.revcomp(contigs[int(rng.integers(0, len(contigs)))])
        return True
    if op == "rename":
        n = names[int(rng.integers(0, len(names)))]
        asm.rename(n, "%s_R%d" % (n, counter))
        return True
    if op == "reorder" and len(names) >= 2:
        perm = list(rng.permutation(len(names)))
        asm.reorder([names[int(i)] for i in perm])
        return True
    if op == "inversion" and contigs:
        n = contigs[int(rng.integers(0, len(contigs)))]
        ln = int(rng.integers(8_000, 20_001))
        reg = asm.free_region(n, ln, rng)
        if reg is None:
            return False
        asm.invert(n, reg[0], reg[1])
        asm.mark(n, reg[0], reg[1])
        return True
    if op == "translocation" and contigs:
        n = contigs[int(rng.integers(0, len(contigs)))]
        ln = int(rng.integers(5_000, 15_001))
        reg = asm.free_region(n, ln, rng)
        if reg is None:
            return False
        dreg = asm.free_region(n, 1, rng)
        if dreg is None:
            return False
        d = dreg[0]
        if reg[0] <= d < reg[1]:
            return False
        ok = asm.translocate(n, reg[0], reg[1], d if d > reg[1] else d)
        if ok:
            asm.mark(n, reg[0], reg[1])
        return ok
    if op == "insertion" and contigs:
        n = contigs[int(rng.integers(0, len(contigs)))]
        pt = asm.free_region(n, 1, rng)
        if pt is None:
            return False
        asm.insert(n, pt[0], int(rng.integers(500, 5_001)), rng)
        return True
    if op == "deletion" and contigs:
        n = contigs[int(rng.integers(0, len(contigs)))]
        ln = int(rng.integers(500, 5_001))
        reg = asm.free_region(n, ln, rng)
        if reg is None:
            return False
        asm.delete(n, reg[0], reg[1])
        asm.mark(n, reg[0], reg[1])
        return True
    if op == "expansion" and contigs:
        n = contigs[int(rng.integers(0, len(contigs)))]
        ln = int(rng.integers(1_000, 4_001))
        reg = asm.free_region(n, ln, rng)
        if reg is None:
            return False
        asm.expand(n, reg[0], ln)
        asm.mark(n, reg[0], reg[1])
        return True
    if op == "snps" and contigs:
        n = contigs[int(rng.integers(0, len(contigs)))]
        asm.snps(n, int(rng.integers(5, 16)), rng)
        return True
    if op == "add_plasmid":
        asm.add_contig("p_add%d" % counter, random_seq(rng, int(rng.integers(15_000, 40_001))))
        return True
    if op == "remove_plasmid" and len(names) > 1:
        n = names[int(rng.integers(0, len(names)))]
        asm.remove_contig(n)
        return True
    return False


def sc_sr(base, reps, is_seq, rng):
    asm = _Asm(base, reps)
    n_ops = int(rng.integers(3, 6))
    done = 0
    for attempt in range(n_ops * 15):
        if done >= n_ops:
            break
        op = _RANDOM_OPS[int(rng.integers(0, len(_RANDOM_OPS)))]
        if _random_op(asm, op, rng, attempt):
            done += 1
    return {"old": dict(base), "new": asm.c, "reps_old": copy_repeats(reps),
            "reps_new": asm.reps, "edits": asm.edits,
            "description": "random combination of %d edits" % done,
            "note": "ops: " + ", ".join(e["type"] for e in asm.edits)}


SCENARIOS = [
    ("S1", sc_s1), ("S2", sc_s2), ("S3", sc_s3), ("S4", sc_s4),
    ("S5", sc_s5), ("S6", sc_s6), ("S6b", sc_s6b), ("S7", sc_s7),
    ("S8", sc_s8), ("S9", sc_s9), ("S10a", sc_s10a), ("S10b", sc_s10b),
    ("S11", sc_s11), ("S12", sc_s12),
    ("S5c", sc_s5c), ("S6c", sc_s6c), ("S6w", sc_s6w),
    # added 2026-10-17 (appended so existing seed indices are unchanged)
    ("S13", sc_s13), ("S13b", sc_s13b), ("S14", sc_s14), ("S14b", sc_s14b),
    ("S15", sc_s15), ("S16", sc_s16), ("S17", sc_s17), ("S17b", sc_s17b),
    ("S18", sc_s18), ("S19", sc_s19), ("S20", sc_s20), ("S21", sc_s21),
    ("SR1", sc_sr), ("SR2", sc_sr), ("SR3", sc_sr), ("SR4", sc_sr),
    ("SR5", sc_sr),
]

# S5c/S6c reuse the exact S5/S6 genome-generation seed so their old.fa/new.fa
# (and truth edits/repeats) are identical to those scenarios.
GEN_SEED_ALIAS = {"S5c": "S5", "S6c": "S6"}

# ----------------------------------------------------------------------------
# read simulation


def diverge(seq, rng, sub=0.08, indel=0.02):
    code = _B2I[np.frombuffer(seq.encode(), dtype=np.uint8)].copy()
    n = len(code)
    m = rng.random(n) < sub
    if m.any():
        code[m] = (code[m] + rng.integers(1, 4, size=int(m.sum()))) % 4
    mi = rng.random(n) < indel
    if mi.any():
        pos = np.nonzero(mi)[0]
        lens = rng.integers(1, 4, size=len(pos))
        idxs, vals = [], []
        for p, l in zip(pos, lens):
            l = int(l)
            idxs.extend([int(p) + 1] * l)
            vals.extend(rng.integers(0, 4, size=l))
        code = np.insert(code, idxs, vals)
    n2 = len(code)
    md = rng.random(n2) < indel
    code = code[~md]
    return _IDX2BASE[code].tobytes().decode()


def hp_runlen(code):
    n = len(code)
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    change = np.empty(n, dtype=bool)
    change[0] = True
    change[1:] = code[1:] != code[:-1]
    run_id = np.cumsum(change) - 1
    counts = np.bincount(run_id)
    return counts[run_id]


def inject_errors(seq, rng, model):
    code = _B2I[np.frombuffer(seq.encode(), dtype=np.uint8)].copy()
    n = len(code)
    if model == "hifi":
        sub_r = ins_r = del_r = 0.001
        hp_bias = False
    else:
        sub_r, ins_r, del_r = 0.0075, 0.00375, 0.00375
        hp_bias = True

    m = rng.random(n) < sub_r
    if m.any():
        code[m] = (code[m] + rng.integers(1, 4, size=int(m.sum()))) % 4

    hp = hp_runlen(code) >= 3 if hp_bias else np.zeros(n, dtype=bool)
    ins_p = ins_r * np.where(hp, 3.0, 1.0)
    mi = rng.random(n) < ins_p
    if mi.any():
        pos = np.nonzero(mi)[0]
        code = np.insert(code, pos + 1, rng.integers(0, 4, size=len(pos)))

    n2 = len(code)
    hp2 = hp_runlen(code) >= 3 if hp_bias else np.zeros(n2, dtype=bool)
    del_p = del_r * np.where(hp2, 3.0, 1.0)
    md = rng.random(n2) < del_p
    if md.any():
        code = code[~md]

    if model == "hifi":
        q = rng.integers(30, 41, size=len(code))
    else:
        q = rng.integers(10, 26, size=len(code))
    qual = (q + 33).astype(np.uint8).tobytes().decode()
    return _IDX2BASE[code].tobytes().decode(), qual


def sample_len(rng, model):
    if model == "hifi":
        m, s, mn = 15_000.0, 3_000.0, 2_000
    else:
        m, s, mn = 10_000.0, 8_000.0, 1_000
    sig2 = np.log(1.0 + (s / m) ** 2)
    mu = np.log(m) - sig2 / 2.0
    x = float(rng.lognormal(mu, np.sqrt(sig2)))
    return max(mn, int(round(x)))


def fragment(rng, contigs, names, cum, total, model):
    u = rng.random() * total
    ci = int(np.searchsorted(cum, u, side="right"))
    if ci >= len(names):
        ci = len(names) - 1
    name = names[ci]
    seq = contigs[name]
    L = len(seq)
    start = int(rng.integers(0, L))
    ln = sample_len(rng, model)
    if ln >= L:
        ln = L
        template = seq[start:] + seq[:start]
    else:
        end = start + ln
        template = seq[start:end] if end <= L else seq[start:] + seq[:end - L]
    strand = "+" if rng.random() < 0.5 else "-"
    eseq, equal = inject_errors(template, rng, model)
    if strand == "-":
        eseq, equal = revcomp(eseq), equal[::-1]
    return name, start, start + ln, strand, eseq, equal


def extract(seq, start, length):
    """Circular extraction of `length` bases starting at `start`."""
    L = len(seq)
    if length >= L:
        return seq[start:] + seq[:start]
    end = start + length
    if end <= L:
        return seq[start:end]
    return seq[start:] + seq[:end - L]


def chimp3_read(rng, new_contigs, model, counter, is_region):
    """Targeted chimera: chr fragment joined to a short p3 fragment."""
    chr_seq, p3_seq = new_contigs["chr"], new_contigs["p3"]
    Lc, L3 = len(chr_seq), len(p3_seq)

    clen = int(rng.integers(2_000, 15_001))
    cstart = int(rng.integers(0, Lc))
    cstrand = "+" if rng.random() < 0.5 else "-"
    cseq = extract(chr_seq, cstart, clen)

    plen = int(round(float(np.exp(rng.uniform(np.log(30.0), np.log(3000.0))))))
    plen = max(30, min(3000, plen))
    if is_region is not None and rng.random() < 0.5:
        s0, e0 = is_region
        lo = max(0, s0 - plen + 1)
        pstart = int(rng.integers(lo, e0))  # guarantees overlap with the IS copy
    else:
        pstart = int(rng.integers(0, L3))
    pstrand = "+" if rng.random() < 0.5 else "-"
    pseq = extract(p3_seq, pstart, plen)
    pend = pstart + plen

    if cstrand == "-":
        cseq = revcomp(cseq)
    if pstrand == "-":
        pseq = revcomp(pseq)
    joined = cseq + pseq if rng.random() < 0.5 else pseq + cseq

    eseq, equal = inject_errors(joined, rng, model)
    name = "r%06d:chimp3:chr:%d-%d:%s|p3:%d-%d:%s:len%d" % (
        counter, cstart, cstart + clen, cstrand, pstart, pend, pstrand, plen)
    return name, eseq, equal


def simulate_reads(new_contigs, model, rng, targeted=None):
    names = list(new_contigs.keys())
    lens = np.array([len(new_contigs[n]) for n in names], dtype=np.float64)
    cum = np.cumsum(lens)
    total = float(cum[-1])
    target = 50.0 * total

    reads = []
    counter = 0
    acc = 0.0
    while acc < target:
        c, s, e, st, seq, qual = fragment(rng, new_contigs, names, cum, total, model)
        reads.append(("r%06d:%s:%d-%d:%s" % (counter, c, s, e, st), seq, qual))
        counter += 1
        acc += len(seq)

    n = len(reads)
    for _ in range(int(round(0.005 * n))):
        f1 = fragment(rng, new_contigs, names, cum, total, model)
        f2 = fragment(rng, new_contigs, names, cum, total, model)
        nm = "r%06d:chim:%s:%d-%d:%s|%s:%d-%d:%s" % (
            counter, f1[0], f1[1], f1[2], f1[3], f2[0], f2[1], f2[2], f2[3])
        reads.append((nm, f1[4] + f2[4], f1[5] + f2[5]))
        counter += 1

    for _ in range(int(round(0.002 * n))):
        ln = int(rng.integers(2_000, 10_001))
        seq = random_seq(rng, ln)
        if model == "hifi":
            q = rng.integers(30, 41, size=ln)
        else:
            q = rng.integers(10, 26, size=ln)
        qual = (q + 33).astype(np.uint8).tobytes().decode()
        reads.append(("r%06d:junk" % counter, seq, qual))
        counter += 1

    if targeted is not None:
        for _ in range(int(round(0.10 * n))):
            nm, seq, qual = chimp3_read(rng, new_contigs, model, counter,
                                        targeted.get("is_region"))
            reads.append((nm, seq, qual))
            counter += 1

    order = rng.permutation(len(reads))
    reads = [reads[int(i)] for i in order]
    return reads


def write_fastq(path, reads):
    with open(path, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0,
                           compresslevel=6) as gz:
            for nm, seq, qual in reads:
                gz.write(("@%s\n%s\n+\n%s\n" % (nm, seq, qual)).encode())


# ----------------------------------------------------------------------------
# main


def seed_for(scen_index, model_index):
    # stable integer seeds; never hash() of strings
    return 20_240_001 + 1_000 * scen_index + 10 * model_index


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="data/synthetic")
    ap.add_argument("--only", default=None,
                    help="comma-separated scenario IDs, e.g. S1,S6b")
    args = ap.parse_args()

    wanted = None
    if args.only:
        wanted = set(x.strip() for x in args.only.split(",") if x.strip())

    base_rng = np.random.default_rng(123_456_789)
    base, reps, is_seq = build_base(base_rng)

    idx_of = {sid: i for i, (sid, _) in enumerate(SCENARIOS)}
    for si, (sid, fn) in enumerate(SCENARIOS):
        if wanted is not None and sid not in wanted:
            continue
        alias = GEN_SEED_ALIAS.get(sid, sid)
        rng = np.random.default_rng(seed_for(idx_of[alias], 0) + 7)
        res = fn(base, reps, is_seq, rng)
        outdir = os.path.join(args.out, sid)
        os.makedirs(outdir, exist_ok=True)

        write_fasta(os.path.join(outdir, "old.fa"), res["old"])
        write_fasta(os.path.join(outdir, "new.fa"), res["new"])

        truth = {
            "scenario": sid,
            "description": res["description"],
            "edits": res["edits"],
            "repeats_old": res["reps_old"],
            "repeats_new": res["reps_new"],
            "contigs_old": {k: len(v) for k, v in res["old"].items()},
            "contigs_new": {k: len(v) for k, v in res["new"].items()},
        }
        if res.get("note"):
            truth["note"] = res["note"]
        with open(os.path.join(outdir, "truth.json"), "w") as f:
            json.dump(truth, f, indent=2)
            f.write("\n")

        targeted = None
        if res.get("targeted_chimeras"):
            targeted = {"is_region": res.get("p3_is")}
        for mi, model in enumerate(("hifi", "ont")):
            mdir = os.path.join(outdir, model)
            os.makedirs(mdir, exist_ok=True)
            rrng = np.random.default_rng(seed_for(si, mi))
            reads = simulate_reads(res["new"], model, rrng, targeted)
            write_fastq(os.path.join(mdir, "reads.fastq.gz"), reads)
            total_bp = sum(len(s) for _, s, _ in reads)
            print("%-5s %-5s %6d reads %12d bp" % (sid, model, len(reads), total_bp))


if __name__ == "__main__":
    main()
