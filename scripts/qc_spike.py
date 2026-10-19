#!/usr/bin/env python3
"""Spike known structural/dosage errors into finalized bacterial assemblies (Plan-v3 step 2a).

For each sample x replicate this builds a *spiked* assembly by editing the
longest contig (the chromosome) of the finalized assembly ``new.fa`` with a
fixed list of events (deletions, random insertions, duplications, a collapsed
natural repeat, inversions, a segment swap, SNPs and small indels), writes the
truth (``truth.tsv`` for the small events, ``events.tsv`` for the large ones)
and remaps the **real reads** with minimap2 to the spiked assembly.

No read simulation: the reads are the ones the sample was assembled from.  The
unmodified assembly is mapped once per sample as a false-positive control
(reusing the existing ``oracle.cram`` when its reference is the same ``new.fa``).

Outputs (all under ``scratch/qc/spike/``):
  <case>/rep<n>/{spiked.fa, events.tsv, truth.tsv, regions.bed, summary.json,
                 mapped.cram, mapped.cram.crai, build.json}
  <case>/control/mapped.cram          (symlink to oracle.cram when it matches)
  <case>/repeats.tsv                  (self-alignment repeats, for the scorer)
  build_summary.json

Usage:
  pixi run python scripts/qc_spike.py [--samples SUB[,SUB]] [--replicates 1,2]
      [--threads 16] [--force] [--self-test]

Truth conventions match ``scripts/qc_truth.py`` (imported for ``read_fasta``,
``hp_run_containing``, ``str_span`` and ``event_to_row``): a row's ``ref`` is the
**spiked** (qc reference) base and ``alt`` the base the reads carry, in the
spiked assembly's 0-based coordinates.
"""
import argparse
import json
import random
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
import qc_truth as qt  # noqa: E402  (read_fasta, hp/str, event_to_row)

SPIKE_ROOT = ROOT / "scratch" / "qc" / "spike"

# The six finalized assemblies from the spec (case dir, platform, preset).
SAMPLES = [
    {"case": "ont_FAM22000__cons-to-polished",
     "dir": "data/real/ont_FAM22000__cons-to-polished", "platform": "ont", "preset": "lr:hq"},
    {"case": "ont_FAM14197__cons-to-polished",
     "dir": "data/val/real/ont_FAM14197__cons-to-polished", "platform": "ont", "preset": "lr:hq"},
    {"case": "ont_FAM29080__cons-to-polished",
     "dir": "data/val/real/ont_FAM29080__cons-to-polished", "platform": "ont", "preset": "lr:hq"},
    {"case": "pacbio_FAM14197-p1-1.1__cons-to-finalized",
     "dir": "data/real/pacbio_FAM14197-p1-1.1__cons-to-finalized", "platform": "pacbio", "preset": "map-hifi"},
    {"case": "pacbio_FAM20843-p1-1.1__cons-to-finalized",
     "dir": "data/real/pacbio_FAM20843-p1-1.1__cons-to-finalized", "platform": "pacbio", "preset": "map-hifi"},
    {"case": "pacbio_FAM29030-p1-1.1__cons-to-finalized",
     "dir": "data/val/real/pacbio_FAM29030-p1-1.1__cons-to-finalized", "platform": "pacbio", "preset": "map-hifi"},
]

# Structural events, ordered by descending footprint so the big ones place first.
# (kind, size).  ``swap`` occupies two 10 kb blocks 100-300 kb apart.
STRUCT_SPECS = [
    ("swap", 10000),
    ("del", 20000), ("inv", 20000),
    ("dup", 10000),
    ("del", 5000), ("ins_random", 5000), ("dup", 5000), ("inv", 5000),
    ("del", 1000), ("ins_random", 1000), ("dup", 1000),
    ("del", 200), ("ins_random", 200), ("del", 50), ("ins_random", 50),
]
LARGE_GAP = 60000      # min spacing between structural events (+ event size)
SMALL_GAP = 2000       # min spacing for SNP / small-indel sites (see notes)
MARGIN = 20000         # min distance of any event from the contig ends
SWAP_MIN_GAP = 100000  # swap: 100-300 kb between the two exchanged blocks
SWAP_MAX_GAP = 300000
N_SNP = 30
N_INDEL = 30
INDEL_MAX = 10
HP_MIN = 4
COLLAPSE_MIN_LEN = 1000
COLLAPSE_MIN_IDENT = 0.95
GENOME_SEED_BASE = 20261012

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def revcomp(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def rand_seq(n, gc, rng):
    """Seeded random sequence with the requested GC fraction."""
    out = []
    for _ in range(n):
        if rng.random() < gc:
            out.append(rng.choice("CG"))
        else:
            out.append(rng.choice("AT"))
    return "".join(out)


# ---------------------------------------------------------------------------
# repeat finding (self-alignment) -- the only minimap2 use during planning
# ---------------------------------------------------------------------------
def find_repeats(contig_fa, threads):
    """Return distinct repeat-copy pairs ``[(contig, s0, e0, s1, e1, ident)]``.

    ``minimap2 -X -c -x asm5`` skips the exact diagonal but still reports
    near-diagonal hits with indels, so hits whose two intervals overlap are
    dropped.  Only copies of >= 1 kb and >= 95 % identity are kept.
    """
    cmd = ["minimap2", "-X", "-c", "-x", "asm5", "-t", str(threads),
           str(contig_fa), str(contig_fa)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"minimap2 self failed rc={p.returncode}: {p.stderr[-400:]}")
    pairs = []
    seen = set()
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qn, qs, qe, tn, ts, te = f[0], int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8])
        if qn != tn:
            continue
        tags = {t[:2]: t[5:] for t in f[12:] if len(t) > 4 and t[2] == ":"}
        de = float(tags.get("de", "1"))
        ident = 1.0 - de
        if ident < COLLAPSE_MIN_IDENT:
            continue
        if qs >= te or ts >= qe:          # distinct, non-overlapping copies
            a, b = sorted([(qs, qe), (ts, te)])
        else:
            continue
        if min(a[1] - a[0], b[1] - b[0]) < COLLAPSE_MIN_LEN:
            continue
        key = (a[0], a[1], b[0], b[1])
        if key in seen:
            continue
        seen.add(key)
        pairs.append((qn, a[0], a[1], b[0], b[1], round(ident, 4)))
    pairs.sort(key=lambda r: r[3] - r[1])  # closest copies first
    return pairs


def choose_collapse(pairs, lo, hi):
    """Pick the closest pair whose deletion region fits inside [lo, hi]."""
    for (contig, s0, e0, s1, e1, ident) in pairs:  # pairs sorted by span
        a, b = s0, s1
        L = min(e0 - s0, e1 - s1)
        if a >= lo and b + L <= hi:
            return {"contig": contig, "a": a, "b": b, "L": L, "ident": ident}
    return None


# ---------------------------------------------------------------------------
# planning: choose positions and build edits
# ---------------------------------------------------------------------------
class Placer:
    """Greedy non-overlapping placement on the original contig coordinates."""

    def __init__(self, L, rng, gap):
        self.L = L
        self.rng = rng
        self.gap = gap
        self.lo, self.hi = MARGIN, L - MARGIN
        self.reserved = []               # (start, end) influence intervals
        self.cands = list(range(self.lo, self.hi, 500))
        self.rng.shuffle(self.cands)

    def reserve(self, s, e):
        self.reserved.append((s, e))

    def fits(self, s, e):
        if s < self.lo or e > self.hi:
            return False
        for rs, re in self.reserved:
            if s - self.gap < re and rs - self.gap < e:
                return False
        return True

    def place(self, size, footprint=None):
        """First legal start for an event of ``size`` (footprint = occupied len)."""
        fp = size if footprint is None else footprint
        for p in self.cands:
            if self.fits(p, p + fp):
                self.reserve(p, p + fp)
                return p
        return None

    def place_specific(self, p, footprint):
        if self.fits(p, p + footprint):
            self.reserve(p, p + footprint)
            return True
        return False


def gc_content(seq):
    seq = seq.upper()
    acgt = sum(1 for c in seq if c in "ACGT")
    gc = sum(1 for c in seq if c in "GC")
    return (gc / acgt) if acgt else 0.5


def plan_structural(seq, rng, collapse):
    """Return the structural event list (edits on original coordinates)."""
    L = len(seq)
    pl = Placer(L, rng, LARGE_GAP)
    events = []
    if collapse is not None:
        pl.reserve(collapse["a"], collapse["b"] + collapse["L"])
        events.append({"kind": "collapse", "size": collapse["L"],
                       "a": collapse["a"], "b": collapse["b"],
                       "edit": (collapse["b"], collapse["b"] + collapse["L"], ""),
                       "ident": collapse["ident"]})
    for kind, size in STRUCT_SPECS:
        if kind == "swap":
            S = size
            placed = None
            for a in pl.cands:
                g = rng.randrange(SWAP_MIN_GAP, SWAP_MAX_GAP + 1)
                if pl.fits(a, a + 2 * S + g):
                    placed = (a, a + S + g)
                    break
            if placed is None:
                raise RuntimeError("could not place swap")
            a, b = placed
            inf = pl.place_specific(a, 2 * S + (b - a - S))
            assert inf
            repl = seq[b:b + S] + seq[a + S:b] + seq[a:a + S]
            events.append({"kind": "swap", "size": S, "a": a, "b": b,
                           "S": S, "edit": (a, b + S, repl)})
            continue
        if kind == "del":
            p = pl.place(size)
            if p is None:
                raise RuntimeError(f"could not place del {size}")
            events.append({"kind": "del", "size": size, "p": p,
                           "edit": (p, p + size, "")})
        elif kind == "ins_random":
            p = pl.place(size)
            if p is None:
                raise RuntimeError(f"could not place ins_random {size}")
            ins = rand_seq(size, gc_content(seq), rng)
            events.append({"kind": "ins_random", "size": size, "p": p,
                           "ins": ins, "edit": (p, p, ins)})
        elif kind == "dup":
            p = pl.place(size, footprint=2 * size)
            if p is None:
                raise RuntimeError(f"could not place dup {size}")
            block = seq[p:p + size]
            events.append({"kind": "dup", "size": size, "p": p,
                           "edit": (p, p + size, block + block)})
        elif kind == "inv":
            p = pl.place(size)
            if p is None:
                raise RuntimeError(f"could not place inv {size}")
            events.append({"kind": "inv", "size": size, "p": p,
                           "edit": (p, p + size, revcomp(seq[p:p + size]))})
        else:
            raise ValueError(kind)
    return events, pl


def plan_small(seq, rng, structural_reserved):
    """Place 30 SNPs and 30 small indels (half inside homopolymers >= 4)."""
    events = []
    placer = Placer(len(seq), rng, SMALL_GAP)
    for (s, e) in structural_reserved:          # keep clear of structural events
        placer.reserve(s - 3000, e + 3000)
    # --- SNPs: isolated, non-homopolymer positions -------------------------
    snps = 0
    for p in list(placer.cands):
        if snps >= N_SNP:
            break
        if seq[p] not in "ACGT":
            continue
        if seq[p] == seq[p - 1] or seq[p] == seq[p + 1]:
            continue
        if not placer.place_specific(p, 1):
            continue
        nb = rng.choice([b for b in "ACGT" if b != seq[p]])
        events.append({"kind": "snp", "type": "SNP", "pos": p,
                       "ref": seq[p], "alt": nb, "edit": (p, p + 1, nb)})
        snps += 1
    if snps != N_SNP:
        raise RuntimeError(f"placed only {snps}/{N_SNP} SNPs")

    # --- small indels ------------------------------------------------------
    # Homopolymer subset: pick runs of ACGT, length >= HP_MIN.
    runs = []
    i = MARGIN
    while i < len(seq) - MARGIN:
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if seq[i] in "ACGT" and j - i >= HP_MIN:
            runs.append((i, j))
        i = j
    rng.shuffle(runs)
    n_hp = N_INDEL // 2
    n_free = N_INDEL - n_hp
    placed = 0
    for (a, b) in runs:
        if placed >= n_hp:
            break
        runlen = b - a
        base = seq[a]
        if placed % 2 == 0:      # reads insertion: delete k from the assembly
            if runlen < HP_MIN + 1:          # leave a run of >= HP_MIN behind
                continue
            k = rng.randint(1, min(INDEL_MAX, runlen - HP_MIN))
            q = a
            if not placer.place_specific(q - 1, k + 1):
                continue
            events.append({"kind": "indel", "type": "INS", "pos": q,
                           "seq": seq[q:q + k], "hp": True,
                           "edit": (q, q + k, "")})
        else:                    # reads deletion: insert k copies of the base
            k = rng.randint(1, INDEL_MAX)
            q = a
            ins = base * k
            if not placer.place_specific(q - 1, k + 1):
                continue
            events.append({"kind": "indel", "type": "DEL", "pos": q,
                           "seq": ins, "hp": True,
                           "edit": (q, q, ins)})
        placed += 1
    if placed < n_hp:
        raise RuntimeError(f"placed only {placed}/{n_hp} homopolymer indels")

    n_free_done = 0
    for p in placer.cands:
        if n_free_done >= n_free:
            break
        if seq[p] not in "ACGT":
            continue
        if qt.hp_run_containing(seq, p) > 1:
            continue
        k = rng.randint(1, INDEL_MAX)
        if n_free_done % 2 == 0:  # reads insertion: delete k bases
            if not placer.place_specific(p, k + 1):
                continue
            events.append({"kind": "indel", "type": "INS", "pos": p,
                           "seq": seq[p:p + k], "hp": False,
                           "edit": (p, p + k, "")})
        else:                     # reads deletion: insert k random bases
            ins = rand_seq(k, gc_content(seq), rng)
            if not placer.place_specific(p, k + 1):
                continue
            events.append({"kind": "indel", "type": "DEL", "pos": p,
                           "seq": ins, "hp": False,
                           "edit": (p, p, ins)})
        n_free_done += 1
    if n_free_done != n_free:
        raise RuntimeError(f"placed only {n_free_done}/{n_free} free indels")
    return events


def apply_edits(seq, events):
    """Apply all event edits; return (spiked_seq, sorted_edits).

    Events do not overlap.  ``final_point(x)`` maps an original coordinate to
    the spiked assembly's coordinate (offset accumulated from edits to its left).
    """
    edits = []
    for i, e in enumerate(events):
        s, en, repl = e["edit"]
        edits.append((s, en, repl, i))
    edits.sort(key=lambda t: (t[0], t[1]))
    # sanity: no overlap
    prev_end = -1
    for s, en, repl, i in edits:
        if s < prev_end:
            raise RuntimeError(f"overlapping edits at {s} (< {prev_end})")
        prev_end = en

    def final_point(x):
        return x + sum(len(repl) - (en - s) for (s, en, repl, i) in edits if s < x)

    parts = []
    cur = 0
    for s, en, repl, i in edits:
        parts.append(seq[cur:s])
        parts.append(repl)
        cur = en
    parts.append(seq[cur:])
    return "".join(parts), edits, final_point


def placements(events, spiked, final_point):
    """Attach the spiked-coordinate fields to every event."""
    for e in events:
        k = e["kind"]
        if k == "del":
            j = final_point(e["p"])
            e.update(start0=j, end0=j, breakpoints=[j])
        elif k == "ins_random":
            fp = final_point(e["p"])
            e.update(start0=fp, end0=fp + e["size"], breakpoints=[fp, fp + e["size"]])
        elif k == "dup":
            fp = final_point(e["p"])
            e.update(start0=fp, end0=fp + 2 * e["size"],
                     breakpoints=[fp + e["size"], fp + 2 * e["size"]])
        elif k == "inv":
            fi = final_point(e["p"])
            e.update(start0=fi, end0=fi + e["size"], breakpoints=[fi, fi + e["size"]])
        elif k == "swap":
            af, bf = final_point(e["a"]), final_point(e["b"])
            e.update(start0=af, end0=bf + e["size"],
                     breakpoints=[af, af + e["size"], bf, bf + e["size"]])
        elif k == "collapse":
            cf = final_point(e["a"])
            e.update(start0=cf, end0=cf + e["size"], breakpoints=[final_point(e["b"])])
        elif k == "snp":
            e["pos0"] = final_point(e["pos"])
        elif k == "indel":
            e["pos0"] = final_point(e["pos"])
    return events


# Columns of ``events.tsv`` (large events). `orig_start0`/`orig_end0` are the
# edited interval in original coordinates (the deleted block for del/collapse).
EVENTS_HEADER = ["event_id", "kind", "size", "contig", "start0", "end0",
                 "breakpoints", "orig_start0", "orig_end0"]


def build_truth_rows(contig, spiked, events):
    """Small-event truth rows (qc_truth format) and large-event rows."""
    truth = []
    large = []
    for idx, e in enumerate(events):
        k = e["kind"]
        if k == "snp":
            p = e["pos0"]
            # REF = spiked (qc reference) base, ALT = original base (what the reads show).
            truth.append((contig, p, "SNP", 1, spiked[p], e["ref"],
                          qt.hp_run_containing(spiked, p), qt.str_span(spiked, p),
                          p, p + 1))
        elif k == "indel":
            p = e["pos0"]
            if e["type"] == "INS":       # reads have extra bases -> truth INS
                row = qt.event_to_row(contig, spiked,
                                      {"kind": "INS", "tpos": p, "seq": e["seq"]})
            else:                        # reads are missing bases -> truth DEL
                row = qt.event_to_row(contig, spiked,
                                      {"kind": "DEL", "tpos": p, "seq": e["seq"]})
            if row is None:
                raise RuntimeError(f"event_to_row dropped {e}")
            truth.append(row)
        else:
            eid = f"{k}_{e['size']}_{idx}"
            # `orig_start0`/`orig_end0`: the edited interval in the *original*
            # (unspiked) coordinates. For `del`/`collapse` this is exactly the
            # deleted block, which the scorer needs to find cross-mapping reads.
            os_, oe_ = e["edit"][0], e["edit"][1]
            large.append((eid, k, e["size"], contig, e["start0"], e["end0"],
                          ",".join(str(x) for x in e["breakpoints"]), os_, oe_))
    return truth, large


# ---------------------------------------------------------------------------
# IO / pipeline
# ---------------------------------------------------------------------------
def wrap(s, w=60):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


def write_fasta(path, seqs, order):
    with open(path, "w") as fh:
        for name in order:
            fh.write(f">{name}\n{wrap(seqs[name])}\n")


def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def read_sq(cram):
    """{(SN): (LN, M5)} for a BAM/CRAM header."""
    p = subprocess.run(["samtools", "view", "-H", str(cram)],
                       capture_output=True, text=True)
    if p.returncode != 0:
        return None
    out = {}
    for line in p.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        d = {}
        for f in line.split("\t")[1:]:
            if ":" in f:
                d[f.split(":", 1)[0]] = f.split(":", 1)[1]
        if "SN" in d and "LN" in d:
            out[d["SN"]] = (int(d["LN"]), d.get("M5", ""))
    return out


def dict_md5(fa):
    p = subprocess.run(["samtools", "dict", str(fa)], capture_output=True, text=True)
    if p.returncode != 0:
        return None
    out = {}
    for line in p.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        d = {}
        for f in line.split("\t")[1:]:
            if ":" in f:
                d[f.split(":", 1)[0]] = f.split(":", 1)[1]
        if "SN" in d and "LN" in d:
            out[d["SN"]] = (int(d["LN"]), d.get("M5", ""))
    return out


DATA_ROOT = (ROOT / "data").resolve()


def _remove_output(path):
    """Unlink a regular file or a symlink (never follow it)."""
    p = Path(path)
    if p.is_symlink():
        p.unlink()
    elif p.exists():
        p.unlink()


def _guard_output(path):
    """Refuse to write through a symlink or under the read-only ``data/`` tree."""
    p = Path(path)
    if p.is_symlink():
        raise RuntimeError(f"refusing to write through symlink {p} -> {p.resolve()}")
    rp = p.resolve()
    if rp == DATA_ROOT or DATA_ROOT in rp.parents:
        raise RuntimeError(f"refusing to write under data/: {rp}")


def map_reads(spiked_fa, reads, preset, threads, out_cram):
    """minimap2 | samtools sort -> CRAM (+ index).  Returns wall seconds."""
    _guard_output(out_cram)
    _guard_output(str(out_cram) + ".crai")
    _remove_output(out_cram)
    # Index the reference only when no index exists, and never under data/ (read-only):
    # the control maps to data/.../new.fa, whose .fai is already present.
    fai = Path(str(spiked_fa) + ".fai")
    if not fai.exists():
        _guard_output(fai)
        subprocess.run(["samtools", "faidx", str(spiked_fa)], check=True,
                       capture_output=True, text=True)
    cmd = (f"minimap2 -ax {preset} -Y --secondary=no -t {int(threads)} "
           f"{spiked_fa} {reads} | "
           f"samtools sort -@ 8 -o {out_cram} --reference {spiked_fa}")
    t0 = time.monotonic()
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"mapping failed rc={p.returncode}: {p.stderr[-600:]}")
    subprocess.run(["samtools", "index", str(out_cram)], check=True,
                   capture_output=True, text=True)
    return time.monotonic() - t0


def plan_rep(seqs, contig, collapse_spec, seed):
    """Deterministically plan and apply one replicate's events."""
    rng = random.Random(GENOME_SEED_BASE + seed)
    struct_events, placer = plan_structural(seqs[contig], rng, collapse_spec)
    small_events = plan_small(seqs[contig], rng, placer.reserved)
    all_events = struct_events + small_events
    spiked_chr, edits, final_point = apply_edits(seqs[contig], all_events)
    placements(all_events, spiked_chr, final_point)
    truth, large = build_truth_rows(contig, spiked_chr, all_events)
    return all_events, spiked_chr, truth, large


def build_rep(sample, case_dir, rep, seed, threads, force):
    out = SPIKE_ROOT / sample["case"] / f"rep{rep}"
    out.mkdir(parents=True, exist_ok=True)
    spiked_fa = out / "spiked.fa"
    done = all((out / f).exists() for f in
               ("spiked.fa", "events.tsv", "truth.tsv", "mapped.cram"))
    if done and not force:
        try:
            info = json.loads((out / "build.json").read_text())
        except Exception:
            info = {}
        if info.get("mapping_time_s") is not None:
            return info, "cached"

    seqs = qt.read_fasta(case_dir / "new.fa")
    if not seqs:
        raise RuntimeError("empty new.fa")
    order = list(seqs)
    contig = max(order, key=lambda c: len(seqs[c]))
    L = len(seqs[contig])
    lo, hi = MARGIN, L - MARGIN

    # self-alignment repeats of the whole assembly (per sample, reused by scorer)
    rep_fa = out / "_chr.fa"
    write_fasta(rep_fa, {contig: seqs[contig]}, [contig])
    pairs = find_repeats(rep_fa, threads)
    rep_rows = [(c, s0, e0, s1, e1, ident) for (c, s0, e0, s1, e1, ident) in pairs]
    write_tsv(SPIKE_ROOT / sample["case"] / "repeats.tsv",
              ["contig", "start0", "end0", "mate_start", "mate_end", "identity"],
              rep_rows)
    collapse_spec = choose_collapse(pairs, lo, hi)

    all_events, spiked_chr, truth, large = plan_rep(seqs, contig, collapse_spec, seed)

    # write spiked assembly (only the chromosome changed)
    spiked = dict(seqs)
    spiked[contig] = spiked_chr
    write_fasta(spiked_fa, spiked, order)
    write_tsv(out / "truth.tsv",
              ["contig", "pos0", "type", "len", "ref", "alt", "hp", "str",
               "start0", "end0"], truth)
    write_tsv(out / "events.tsv", EVENTS_HEADER, large)
    (out / "regions.bed").write_text("")
    (out / "summary.json").write_text(json.dumps({
        "contig_lengths": {c: len(s) for c, s in spiked.items()},
        "contig_order": order, "chromosome": contig,
        "platform": sample["platform"], "preset": sample["preset"]}) + "\n")

    t_map = map_reads(spiked_fa, case_dir / "reads.fastq.gz", sample["preset"],
                      threads, out / "mapped.cram")
    rep_fa.unlink(missing_ok=True)

    per_kind = {}
    for e in all_events:
        if e["kind"] in ("snp", "indel"):
            key = "snp" if e["kind"] == "snp" else f"indel_{e['type']}"
        else:
            key = f"{e['kind']}_{e['size']}"
        per_kind[key] = per_kind.get(key, 0) + 1
    info = {
        "case": sample["case"], "rep": rep, "seed": seed,
        "platform": sample["platform"], "preset": sample["preset"],
        "contig": contig, "contig_len": L,
        "collapse": collapse_spec,
        "events_per_kind": per_kind,
        "n_events": len(all_events),
        "n_truth_rows": len(truth), "n_large_events": len(large),
        "spiked_bp": sum(len(s) for s in spiked.values()),
        "mapping_time_s": round(t_map, 2),
        "spiked_fa_bytes": spiked_fa.stat().st_size,
        "cram_bytes": (out / "mapped.cram").stat().st_size,
    }
    (out / "build.json").write_text(json.dumps(info, indent=2) + "\n")
    return info, "ok"


def regenerate_truth(sample, case_dir, rep, seed):
    """Rewrite truth.tsv and events.tsv for a cached rep (no remapping)."""
    out = SPIKE_ROOT / sample["case"] / f"rep{rep}"
    info = json.loads((out / "build.json").read_text())
    collapse = info.get("collapse")
    seqs = qt.read_fasta(case_dir / "new.fa")
    contig = info["contig"]
    all_events, spiked_chr, truth, large = plan_rep(seqs, contig, collapse, seed)
    cached = qt.read_fasta(out / "spiked.fa")
    if cached.get(contig) != spiked_chr:
        raise RuntimeError(f"{sample['case']} rep{rep}: regenerated chromosome "
                           f"differs from cached spiked.fa -- refusing to rewrite truth")
    write_tsv(out / "truth.tsv",
              ["contig", "pos0", "type", "len", "ref", "alt", "hp", "str",
               "start0", "end0"], truth)
    write_tsv(out / "events.tsv", EVENTS_HEADER, large)
    return len(truth)


def build_control(sample, case_dir, threads, force):
    out = SPIKE_ROOT / sample["case"] / "control"
    out.mkdir(parents=True, exist_ok=True)
    dst = out / "mapped.cram"
    if (dst.is_symlink() or dst.exists()) and not force:
        reused = dst.is_symlink() and dst.resolve().name.startswith("oracle")
        return {"reused_oracle": reused, "cram_bytes": dst.stat().st_size}, "cached"
    oracle = case_dir / "oracle.cram"
    same = False
    if oracle.exists():
        sq_o = read_sq(oracle)
        sq_n = dict_md5(case_dir / "new.fa")
        same = bool(sq_o and sq_n and sq_o == sq_n)
    if same:
        _remove_output(dst)
        dst.symlink_to(oracle)
        idx = out / "mapped.cram.crai"
        _remove_output(idx)
        idx.symlink_to(oracle.with_suffix(".cram.crai"))
        info = {"reused_oracle": True, "oracle": str(oracle),
                "cram_bytes": oracle.stat().st_size}
        (out / "build.json").write_text(json.dumps(info, indent=2) + "\n")
        return info, "oracle"
    # remap to the unmodified assembly: never write through the control symlinks
    _remove_output(dst)
    _remove_output(out / "mapped.cram.crai")
    t_map = map_reads(case_dir / "new.fa", case_dir / "reads.fastq.gz",
                      sample["preset"], threads, dst)
    info = {"reused_oracle": False, "mapping_time_s": round(t_map, 2),
            "cram_bytes": dst.stat().st_size}
    (out / "build.json").write_text(json.dumps(info, indent=2) + "\n")
    return info, "mapped"


def run(only, replicates, threads, force):
    summary = {}
    for sample in SAMPLES:
        if only and not any(s in sample["case"] for s in only):
            continue
        case_dir = ROOT / sample["dir"]
        if not (case_dir / "new.fa").exists():
            summary[sample["case"]] = {"error": "missing new.fa"}
            continue
        print(f"=== {sample['case']} ({sample['platform']}, {sample['preset']})")
        cinfo, cstatus = build_control(sample, case_dir, threads, force)
        print(f"  control: {cstatus} {cinfo}")
        reps = {}
        for rep in replicates:
            info, status = build_rep(sample, case_dir, rep, rep, threads, force)
            reps[rep] = info
            pk = info.get("events_per_kind", {})
            col = info.get("collapse")
            cdesc = "no" if not col else f"a={col['a']} b={col['b']} L={col['L']} ident={col['ident']}"
            print(f"  rep{rep}: {status} events={info.get('n_events')} "
                  f"map={info.get('mapping_time_s')}s "
                  f"cram={info.get('cram_bytes', 0) / 1e6:.0f}MB collapse={cdesc}")
            print(f"    per kind: {pk}")
        summary[sample["case"]] = {"control": cinfo, "reps": reps,
                                   "platform": sample["platform"],
                                   "preset": sample["preset"]}
    SPIKE_ROOT.mkdir(parents=True, exist_ok=True)
    (SPIKE_ROOT / "build_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nwrote {SPIKE_ROOT}/build_summary.json")
    return 0


# ---------------------------------------------------------------------------
# self-test (no minimap2, no reads)
# ---------------------------------------------------------------------------
def _plant_repeat(seq, rng, L_block=2000):
    """Plant an identical tandem copy; return (seq, a, b, L)."""
    a = 60000
    block = seq[a:a + L_block]
    b = a + L_block
    return seq[:b] + block + seq[b:], a, b, L_block


def self_test(_args):
    problems = []
    rng = random.Random(1234)
    base = rand_seq(3_500_000, 0.5, rng)
    # a long homopolymer for the indel tests
    base = base[:150000] + "A" * 12 + base[150012:]
    seq, ca, cb, cL = _plant_repeat(base, rng)
    collapse = {"contig": "chr", "a": ca, "b": cb, "L": cL, "ident": 1.0}

    rng2 = random.Random(1234)
    struct_events, placer = plan_structural(seq, rng2, collapse)
    small_events = plan_small(seq, rng2, placer.reserved)
    events = struct_events + small_events
    spiked, edits, final_point = apply_edits(seq, events)
    placements(events, spiked, final_point)
    truth, large = build_truth_rows("chr", spiked, events)

    def ev(kind, size=None):
        for e in events:
            if e["kind"] == kind and (size is None or e.get("size") == size):
                return e
        return None

    # every structural kind present
    for k, s in STRUCT_SPECS:
        if ev(k, s) is None:
            problems.append(f"missing structural event {k} {s}")
    if ev("collapse") is None:
        problems.append("missing collapse")

    # --- structural coordinate checks -------------------------------------
    for e in struct_events:
        k = e["kind"]
        if k == "del":
            p, s, j = e["p"], e["size"], e["start0"]
            if not (spiked[j - 1] == seq[p - 1] and spiked[j] == seq[p + s]):
                problems.append(f"del {s}: junction mismatch at {j}")
        elif k == "ins_random":
            fp, s = e["start0"], e["size"]
            if not (spiked[fp - 1] == seq[e["p"] - 1] and spiked[fp + s] == seq[e["p"]]
                    and len(spiked[fp:fp + s]) == s):
                problems.append(f"ins_random {s}: mismatch at {fp}")
        elif k == "dup":
            fp, s = e["start0"], e["size"]
            if not (spiked[fp:fp + s] == seq[e["p"]:e["p"] + s]
                    and spiked[fp + s:fp + 2 * s] == seq[e["p"]:e["p"] + s]):
                problems.append(f"dup {s}: copies mismatch at {fp}")
        elif k == "inv":
            fi, s = e["start0"], e["size"]
            if spiked[fi:fi + s] != revcomp(seq[e["p"]:e["p"] + s]):
                problems.append(f"inv {s}: mismatch at {fi}")
        elif k == "swap":
            a, b, S = e["a"], e["b"], e["S"]
            af, bf = final_point(a), final_point(b)
            if not (spiked[af:af + S] == seq[b:b + S]
                    and spiked[af + S:bf] == seq[a + S:b]
                    and spiked[bf:bf + S] == seq[a:a + S]):
                problems.append(f"swap {S}: blocks mismatch")
        elif k == "collapse":
            cf, j, L = e["start0"], e["breakpoints"][0], e["size"]
            if spiked[cf:cf + L] != seq[e["a"]:e["a"] + L]:
                problems.append("collapse: remaining copy mismatch")
            if not (spiked[j - 1] == seq[e["b"] - 1] and spiked[j] == seq[e["b"] + L]):
                problems.append("collapse: junction mismatch")

    # --- small events ------------------------------------------------------
    n_snp = sum(1 for e in events if e["kind"] == "snp")
    n_ins = sum(1 for e in events if e["kind"] == "indel" and e["type"] == "INS")
    n_del = sum(1 for e in events if e["kind"] == "indel" and e["type"] == "DEL")
    if n_snp != N_SNP:
        problems.append(f"planned {n_snp} SNPs != {N_SNP}")
    if n_ins + n_del != N_INDEL:
        problems.append(f"planned {n_ins + n_del} indels != {N_INDEL}")
    typ_counts = {}
    for r in truth:
        typ_counts[r[2]] = typ_counts.get(r[2], 0) + 1
    if typ_counts.get("SNP", 0) != N_SNP:
        problems.append(f"truth SNPs {typ_counts.get('SNP')} != {N_SNP}")
    if typ_counts.get("INS", 0) != n_ins or typ_counts.get("DEL", 0) != n_del:
        problems.append(f"truth indels {typ_counts} != {n_ins}/{n_del}")
    hp_indels = sum(1 for r in truth if r[2] in ("INS", "DEL") and r[6] >= HP_MIN)
    if hp_indels < N_INDEL // 2:
        problems.append(f"homopolymer indels {hp_indels} < {N_INDEL // 2}")
    # each SNP row: REF = spiked base, ALT = original base (reads), and REF != ALT
    for e in events:
        if e["kind"] != "snp":
            continue
        p = e["pos0"]
        row = next((r for r in truth if r[2] == "SNP" and r[1] == p), None)
        if row is None or row[4] != spiked[p] or row[5] != seq[e["pos"]]:
            problems.append(f"SNP truth row at {p} wrong: {row}")
        elif row[4] == row[5]:
            problems.append(f"SNP truth row at {p} has REF == ALT ({row[4]})")

    # --- events.tsv has every structural event -----------------------------
    if len(large) != len(struct_events):
        problems.append(f"events.tsv {len(large)} != {len(struct_events)}")

    if problems:
        print("SELF-TEST FAILED:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 1
    print(f"self-test OK: contigs=1 coll@{ca}+{cb} L={cL} "
          f"struct={len(struct_events)} snp={n_snp} indel={n_ins}INS/{n_del}DEL "
          f"hp_indels={hp_indels} truth_rows={len(truth)} events={len(large)} "
          f"spiked_len={len(spiked)}")
    return 0


def run_truth_only(only, replicates):
    """Rewrite truth.tsv for cached reps without remapping (spiked.fa unchanged)."""
    for sample in SAMPLES:
        if only and not any(s in sample["case"] for s in only):
            continue
        case_dir = ROOT / sample["dir"]
        for rep in replicates:
            out = SPIKE_ROOT / sample["case"] / f"rep{rep}"
            if not (out / "spiked.fa").exists():
                print(f"  skip {sample['case']} rep{rep}: no spiked.fa")
                continue
            n = regenerate_truth(sample, case_dir, rep, rep)
            print(f"  {sample['case']} rep{rep}: rewrote truth.tsv + events.tsv ({n} rows)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", default=None, help="comma-separated case substrings")
    ap.add_argument("--replicates", default="1,2")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--truth-only", action="store_true",
                    help="regenerate truth.tsv from cached spiked.fa (no remapping)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test(args)
    only = [s for s in args.samples.split(",") if s] if args.samples else None
    reps = [int(x) for x in args.replicates.split(",") if x]
    if args.truth_only:
        return run_truth_only(only, reps)
    return run(only, reps, args.threads, args.force)


if __name__ == "__main__":
    sys.exit(main())
