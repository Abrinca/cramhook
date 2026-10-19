#!/usr/bin/env python3
"""Build truth for the `cramhook qc` command (Plan-v3 step 1).

Truth = errors of the old (raw) assembly, i.e. the differences old -> new in
**old (target)** coordinates, taken from

    minimap2 -c -x asm5 --cs=long -t T old.fa new.fa

Only primary alignments (``tp:A:P``) are used.  The ``cs`` string is parsed into
SNP / INS / DEL (small, < 50 bp) / LINS / LDEL (>= 50 bp) events at old
coordinates, using the VCF anchor convention (0-based anchor base precedes the
indel) and left-normalising indels against the old sequence.  ``hp`` is the
homopolymer run length in old: starting at ``pos0 + 1`` for indels, the run
containing ``pos0`` for SNPs.

The ``cs`` tag is always in the target (old) forward orientation, even when the
query strand is ``-`` (verified with a reverse-complemented query); so the
strand column is not needed for parsing.

Also writes the old intervals covered by no primary alignment (``unaligned``)
and by more than one primary alignment (``multi``) into ``regions.bed``.

Usage:
  pixi run python scripts/qc_truth.py [--sets m6large,real,val]
      [--cases SUBSTR[,SUBSTR]] [--out scratch/qc/truth] [--threads 16]
      [--force] [--self-test]

Per case outputs (under ``--out``): ``<set>/<case>/{truth.tsv,regions.bed,summary.json}``.
Cases with missing ``old.fa``/``new.fa`` are skipped and reported.  Existing
outputs are reused unless ``--force`` is given.
"""
import argparse
import json
import random
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SET_DIRS = {
    "m6large": ROOT / "data" / "m6large" / "real",
    "real": ROOT / "data" / "real",
    "val": ROOT / "data" / "val" / "real",
}
SET_ORDER = ["m6large", "real", "val"]
MANIFEST = "manifest.tsv"
LARGE = 50                    # >= 50 bp -> LINS / LDEL
CS_BASES = set("ACGTNacgtn")
TYPES = ["SNP", "INS", "DEL", "LINS", "LDEL"]


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def read_manifest(root):
    """case -> manifest row (platform, preset, ...)."""
    rows = {}
    path = Path(root) / MANIFEST
    if not path.exists():
        return rows
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            rows[parts[0]] = dict(zip(header, parts))
    return rows


def read_fasta(path):
    """Stream a FASTA into an ordered ``{name: uppercase sequence}`` dict."""
    seqs = {}
    name = None
    chunks = []
    with open(path, "r") as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def hp_run_starting_at(seq, s):
    """Length of the homopolymer run starting at 0-based ``s`` (0 if N/out)."""
    if s < 0 or s >= len(seq):
        return 0
    b = seq[s]
    if b not in "ACGT":
        return 0
    e = s
    while e < len(seq) and seq[e] == b:
        e += 1
    return e - s


def hp_run_containing(seq, p):
    """Length of the homopolymer run containing 0-based ``p`` (0 if N/out)."""
    if p < 0 or p >= len(seq):
        return 0
    b = seq[p]
    if b not in "ACGT":
        return 0
    a = p
    while a > 0 and seq[a - 1] == b:
        a -= 1
    e = p
    while e + 1 < len(seq) and seq[e + 1] == b:
        e += 1
    return e - a + 1


def str_span(seq, p):
    """STR span (bases covered by >= 3 tandem copies, period 2-4) at ``p``.

    Returns the span of the largest period-2..4 tandem repeat (>= 3 copies)
    whose region contains 0-based position ``p``, or 0 if there is none.

    A region ``[a, b)`` is period-``k`` when ``seq[i] == seq[i-k]`` for every
    ``i`` in ``[a+k, b)``; such a region corresponds to a maximal run of true
    values in ``ok[i] = (seq[i] == seq[i+k])`` as ``[L, R]`` -> ``[L, R+k+1)``.
    Only blocks intersecting the ``k+1``-wide window ``[p-k, p]`` can contain
    ``p``, so this stays local.
    """
    n = len(seq)
    if p < 0 or p >= n:
        return 0
    best = 0
    for k in (2, 3, 4):
        if n < 2 * k:
            continue
        last = n - k - 1          # ok[] is defined for 0 <= i <= last

        def okval(i):
            # non-ACGT (e.g. N) is a boundary: both compared bases must be ACGT
            return 0 <= i <= last and seq[i] in "ACGT" and seq[i] == seq[i + k]

        lo, hi = max(0, p - k), min(last, p)
        i = lo
        while i <= hi:
            if not okval(i):
                i += 1
                continue
            L = i
            while L - 1 >= 0 and okval(L - 1):
                L -= 1
            R = i
            while R + 1 <= last and okval(R + 1):
                R += 1
            if L <= p <= R + k:                     # region [L, R+k+1) contains p
                span = R + k + 1 - L
                if span >= 3 * k:
                    best = max(best, span)
            i = R + 1
    return best


def hp_bucket(hp):
    if hp <= 0:
        return "0"
    if hp >= 8:
        return "8+"
    return str(hp)


# ---------------------------------------------------------------------------
# cs parsing -> old-coordinate events
# ---------------------------------------------------------------------------
def parse_cs_events(cs, tstart, old_seq):
    """Parse a minimap2 ``cs`` string into raw events in old coordinates.

    Returns a list of dicts: ``{"kind": SNP|INS|DEL, "tpos": int, ...}`` where
    ``tpos`` is the old (target) coordinate of the event (for an indel, the
    position of the first affected target base; the anchor is ``tpos - 1``).
    """
    out = []
    i = 0
    n = len(cs)
    tpos = tstart
    while i < n:
        c = cs[i]
        if c == ":":                                   # short match run
            j = i + 1
            while j < n and cs[j].isdigit():
                j += 1
            tpos += int(cs[i + 1:j] or "0")
            i = j
        elif c == "=":                                 # long match run
            j = i + 1
            while j < n and cs[j] in CS_BASES:
                j += 1
            tpos += j - (i + 1)
            i = j
        elif c == "*":                                 # substitution (SNP)
            if i + 2 >= n:
                break
            out.append({"kind": "SNP", "tpos": tpos,
                        "ref": cs[i + 1].upper(), "alt": cs[i + 2].upper()})
            tpos += 1
            i += 3
        elif c == "+":                                 # insertion in query
            j = i + 1
            while j < n and cs[j] in CS_BASES:
                j += 1
            out.append({"kind": "INS", "tpos": tpos, "seq": cs[i + 1:j].upper()})
            i = j
        elif c == "-":                                 # deletion from query
            j = i + 1
            while j < n and cs[j] in CS_BASES:
                j += 1
            out.append({"kind": "DEL", "tpos": tpos, "seq": cs[i + 1:j].upper()})
            tpos += j - (i + 1)
            i = j
        elif c == "~":                                 # intron (not for asm5)
            j = i + 1
            while j < n and cs[j] not in ":=*+-~":
                j += 1
            i = j
        else:
            i += 1
    return out


def event_to_row(contig, old_seq, ev):
    """Convert a raw cs event to a truth row (VCF anchor, left-normalised).

    Returns ``(contig, pos0, type, len, ref, alt, hp, str, start0, end0)`` or
    ``None``.  ``start0``/``end0`` are the affected old interval (0-based,
    half-open): for a SNP the single base, for a deletion the deleted bases, and
    for an insertion a degenerate point at the insertion coordinate.  The
    interval is computed from the *actual* affected old coordinate, so a
    right-anchored origin deletion is ``[0, size)`` (never ``[1, size+1)``) and a
    whole-contig deletion ends at the contig end.

    Indels use the VCF anchor convention: ``pos0`` is the 0-based anchor base
    preceding the indel.  An indel that left-normalises to position 0 has no
    preceding base, so it uses the *following* base as a right anchor (the VCF
    rule for indels at the start of a contig): an insertion is
    ``REF=old[0], ALT=ins+old[0]`` and a deletion ``REF=old[0:L+1], ALT=old[L]``.
    A whole-contig deletion is ``REF=old, ALT=''``.  ``hp``/``str`` are taken at
    the actual affected coordinate.
    """
    n = len(old_seq)
    if ev["kind"] == "SNP":
        p = ev["tpos"]
        if p < 0 or p >= n:
            return None
        return (contig, p, "SNP", 1, ev["ref"], ev["alt"],
                hp_run_containing(old_seq, p), str_span(old_seq, p), p, p + 1)

    seq = ev["seq"]
    if not seq:
        return None
    p = ev["tpos"]
    size = len(seq)
    if ev["kind"] == "INS":
        ins = seq
        # left-normalise: shift the insertion left while equivalent
        while p > 0 and old_seq[p - 1] == ins[-1]:
            p -= 1
            ins = ins[-1] + ins[:-1]
        if p > 0:
            anchor = old_seq[p - 1]
            pos0, ref, alt = p - 1, anchor, anchor + ins
        else:                                          # contig start: right anchor
            anchor = old_seq[0]
            pos0, ref, alt = 0, anchor, ins + anchor
        typ = "LINS" if size >= LARGE else "INS"
        return (contig, pos0, typ, size, ref, alt,
                hp_run_starting_at(old_seq, p), str_span(old_seq, p), p, p)

    # DEL: affected interval is the deleted bases [p, p + size)
    d = seq
    while p > 0 and old_seq[p - 1] == d[-1]:
        p -= 1
        d = d[-1] + d[:-1]
    if p > 0:
        anchor = old_seq[p - 1]
        pos0, ref, alt = p - 1, anchor + d, anchor
    elif p + size < n:                                 # contig start: right anchor
        right = old_seq[p + size]
        pos0, ref, alt = 0, d + right, right
    else:                                              # whole contig deleted
        pos0, ref, alt = 0, d, ""
    typ = "LDEL" if size >= LARGE else "DEL"
    return (contig, pos0, typ, size, ref, alt,
            hp_run_starting_at(old_seq, p), str_span(old_seq, p),
            p, p + size)


# ---------------------------------------------------------------------------
# coverage regions (unaligned / multi)
# ---------------------------------------------------------------------------
def compute_regions(contig_lengths, contig_order, intervals):
    """Old intervals with depth 0 (unaligned) or depth >= 2 (multi), merged."""
    out = []
    for contig in contig_order:
        length = contig_lengths.get(contig, 0)
        diff = {}
        for s, e in intervals.get(contig, []):
            if e <= s:
                continue
            diff[s] = diff.get(s, 0) + 1
            diff[e] = diff.get(e, 0) - 1
        if not diff:
            if length > 0:
                out.append((contig, 0, length, "unaligned"))
            continue
        segs = []
        cur, depth = 0, 0
        for pos in sorted(diff):
            if pos > cur:
                segs.append((cur, pos, depth))
            depth += diff[pos]
            cur = pos
        if cur < length:
            segs.append((cur, length, depth))
        run_label = run_start = run_end = None
        for s, e, dep in segs:
            label = "unaligned" if dep == 0 else ("multi" if dep >= 2 else None)
            if label is None:
                if run_label is not None:
                    out.append((contig, run_start, run_end, run_label))
                run_label = None
                continue
            if run_label == label and s == run_end:
                run_end = e
            else:
                if run_label is not None:
                    out.append((contig, run_start, run_end, run_label))
                run_label, run_start, run_end = label, s, e
        if run_label is not None:
            out.append((contig, run_start, run_end, run_label))
    return out


# ---------------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------------
def run_minimap2(old_fa, new_fa, threads):
    cmd = ["minimap2", "-c", "-x", "asm5", "--cs=long",
           "-t", str(threads), str(old_fa), str(new_fa)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"minimap2 failed rc={p.returncode}: {p.stderr.strip()[-500:]}")
    return cmd, p.stdout


def parse_paf(paf_text, old, contig_order, contig_lengths):
    """Primary-only PAF -> (rows, intervals, stats)."""
    rows = []
    intervals = {}
    total = used = 0
    strands = {"+": 0, "-": 0}
    for line in paf_text.splitlines():
        if not line or line[0] == "\t":
            continue
        f = line.split("\t")
        if len(f) < 12:
            continue
        total += 1
        tags = {}
        for t in f[12:]:
            if len(t) >= 5 and t[2] == ":" and t[4] == ":":
                tags[t[:2]] = t[5:]
        if tags.get("tp") != "P":
            continue
        tname = f[5]
        old_seq = old.get(tname)
        if old_seq is None:
            continue
        tstart, tend = int(f[7]), int(f[8])
        strand = f[4]
        cs = tags.get("cs")
        if not cs:
            continue
        used += 1
        strands[strand] = strands.get(strand, 0) + 1
        intervals.setdefault(tname, []).append((tstart, tend))
        for ev in parse_cs_events(cs, tstart, old_seq):
            row = event_to_row(tname, old_seq, ev)
            if row is not None:
                rows.append(row)

    # dedupe exact duplicates (e.g. two queries covering the same old region)
    seen = set()
    uniq = []
    for r in rows:
        if r in seen:
            continue
        seen.add(r)
        uniq.append(r)
    idx = {c: i for i, c in enumerate(contig_order)}
    uniq.sort(key=lambda r: (idx.get(r[0], 1 << 30), r[1], r[2]))
    regions = compute_regions(contig_lengths, contig_order, intervals)
    stats = {"alignments_total": total, "alignments_used": used, "strands": strands}
    return uniq, regions, stats


def run_pipeline(old_fa, new_fa, threads):
    """Run minimap2 + parse; returns a dict with rows/regions/sequence info."""
    t0 = time.monotonic()
    old = read_fasta(old_fa)
    contig_order = list(old)
    contig_lengths = {c: len(s) for c, s in old.items()}
    cmd, paf = run_minimap2(old_fa, new_fa, threads)
    runtime = time.monotonic() - t0
    rows, regions, stats = parse_paf(paf, old, contig_order, contig_lengths)
    stats.update({"runtime_s": runtime, "minimap2_cmd": cmd,
                  "contig_lengths": contig_lengths, "contig_order": contig_order,
                  "old": old, "rows": rows, "regions": regions})
    return stats


def summarise(rows, regions, contig_lengths):
    counts = {t: 0 for t in TYPES}
    hp_by_family = {"SNP": {}, "INS": {}, "DEL": {}}
    str_by_family = {"SNP": {}, "INS": {}, "DEL": {}}
    for row in rows:
        typ, hp, st = row[2], row[6], row[7]
        counts[typ] = counts.get(typ, 0) + 1
        fam = "SNP" if typ == "SNP" else ("INS" if typ in ("INS", "LINS") else "DEL")
        b = hp_bucket(hp)
        hp_by_family[fam][b] = hp_by_family[fam].get(b, 0) + 1
        sb = hp_bucket(st)
        str_by_family[fam][sb] = str_by_family[fam].get(sb, 0) + 1
    unaligned_bp = sum(e - s for _c, s, e, lab in regions if lab == "unaligned")
    multi_bp = sum(e - s for _c, s, e, lab in regions if lab == "multi")
    return {
        "counts": counts,
        "hp_by_family": hp_by_family,
        "str_by_family": str_by_family,
        "unaligned_bp": unaligned_bp,
        "multi_bp": multi_bp,
        "total_bp": sum(contig_lengths.values()),
    }


def write_case(out_dir, res, rows, regions):
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "truth.tsv", "w") as fh:
        fh.write("contig\tpos0\ttype\tlen\tref\talt\thp\tstr\tstart0\tend0\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(out_dir / "regions.bed", "w") as fh:
        for contig, s, e, label in regions:
            fh.write(f"{contig}\t{s}\t{e}\t{label}\n")
    (out_dir / "summary.json").write_text(json.dumps(res, indent=2) + "\n")


def build_case(set_name, case, case_dir, manifest_row, out_root, threads, force):
    old_fa, new_fa = case_dir / "old.fa", case_dir / "new.fa"
    out_dir = out_root / set_name / case
    cached = all((out_dir / f).exists() for f in ("truth.tsv", "regions.bed", "summary.json"))
    if cached and not force:
        try:
            res = json.loads((out_dir / "summary.json").read_text())
        except Exception:
            res = None
        return "cached", res
    if not old_fa.exists() or not new_fa.exists():
        missing = [n for n, f in (("old.fa", old_fa), ("new.fa", new_fa)) if not f.exists()]
        return "missing", {"error": f"missing {','.join(missing)}"}
    pipe = run_pipeline(old_fa, new_fa, threads)
    s = summarise(pipe["rows"], pipe["regions"], pipe["contig_lengths"])
    res = {
        "set": set_name,
        "case": case,
        "platform": manifest_row.get("platform", ""),
        "preset": manifest_row.get("preset", ""),
        "counts": s["counts"],
        "hp_by_family": s["hp_by_family"],
        "str_by_family": s["str_by_family"],
        "unaligned_bp": s["unaligned_bp"],
        "multi_bp": s["multi_bp"],
        "total_bp": s["total_bp"],
        "contig_lengths": pipe["contig_lengths"],
        "alignments_total": pipe["alignments_total"],
        "alignments_used": pipe["alignments_used"],
        "strands": pipe["strands"],
        "runtime_s": round(pipe["runtime_s"], 3),
        "minimap2_cmd": pipe["minimap2_cmd"],
    }
    write_case(out_dir, res, pipe["rows"], pipe["regions"])
    return "ok", res


# ---------------------------------------------------------------------------
# self-test
# ---------------------------------------------------------------------------
def _find_real_old_fa():
    for s in SET_ORDER:
        d = SET_DIRS[s]
        if not d.exists():
            continue
        for c in sorted(p for p in d.iterdir() if p.is_dir()):
            f = c / "old.fa"
            if f.exists():
                return f
    return None


def _aligned_span(old_fa, probe_fa, threads):
    """Primary alignment span of ``old_fa`` against an identical query."""
    _, paf = run_minimap2(old_fa, probe_fa, threads)
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) >= 12 and any(t == "tp:A:P" for t in f[12:]):
            return int(f[7]), int(f[8])
    return None


def _plant_and_check(slice_seq, contig, work, threads):
    """Plant known events in ``slice_seq``, run the pipeline, assert recovery."""
    W = len(slice_seq)
    rng = random.Random(20261011)

    def wrap(s, w=60):
        return "\n".join(s[i:i + w] for i in range(0, len(s), w))

    old_fa = work / "st_old.fa"
    old_fa.write_text(f">{contig}\n{wrap(slice_seq)}\n")
    # find the span minimap2 actually aligns (its ends can be clipped), and only
    # plant events comfortably inside it
    probe_fa = work / "st_probe.fa"
    probe_fa.write_text(f">{contig}\n{wrap(slice_seq)}\n")
    span = _aligned_span(old_fa, probe_fa, threads)
    if span is None:
        return ["could not align the slice to itself"], None, {}, [], [], 0, 0
    lo, hi = span[0] + 1000, span[1] - 1000
    if hi - lo < 50000:
        return [f"aligned span too small to plant events: {span}"], None, {}, [], [], 0, 0

    q = list(slice_seq)
    big_del_at = lo + (hi - lo) // 3
    big_ins_at = lo + 2 * (hi - lo) // 3

    def far(p, points, d):
        return all(abs(p - x) >= d for x in points)

    # --- homopolymer runs -> 3 insertions + 3 deletions --------------------
    runs = []
    i = lo
    while i < hi:
        j = i
        while j < W and q[j] == q[i]:
            j += 1
        if q[i] in "ACGT" and j - i >= 5:
            runs.append((i, j))
        i = j
    hp_ops = []
    picked = []
    for r in runs:
        if len(picked) >= 6:
            break
        if abs(r[0] - big_del_at) < 3000 or abs(r[0] - big_ins_at) < 3000:
            continue
        if not far(r[0], [p[0] for p in picked], 3000):
            continue
        picked.append(r)
    for n, (a, b) in enumerate(picked):
        hp_ops.append(("INS" if n % 2 == 0 else "DEL", a, b))

    # --- SNPs at isolated, non-homopolymer positions -----------------------
    isolated = [i for i in range(lo, hi)
                if q[i] in "ACGT" and q[i] != q[i - 1] and q[i] != q[i + 1]]
    avoid = [big_del_at, big_ins_at] + [r[0] for r in picked]
    planted_snp = []
    chosen = []
    step = (hi - lo) // 12
    for k in range(1, 12):
        target = lo + k * step
        for p in isolated:
            if p < target:
                continue
            if abs(p - big_del_at) < 2100 or abs(p - big_ins_at) < 2100:
                continue
            if not far(p, [r[0] for r in picked], 2000):
                continue
            if not far(p, chosen, 500):
                continue
            chosen.append(p)
            break
    for p in chosen[:8]:
        nb = next(b for b in "ACGT"
                  if b != q[p] and b != q[p - 1] and b != q[p + 1])
        planted_snp.append((p, q[p], nb))
        q[p] = nb

    # --- apply length-changing ops from the highest coordinate down --------
    ops = []
    for kind, a, b in hp_ops:
        ops.append((b if kind == "INS" else a, "ins" if kind == "INS" else "del", 1, q[a]))
    ops.append((big_del_at, "del", 100, None))
    ins_seq = "".join(rng.choice("ACGT") for _ in range(300))
    ops.append((big_ins_at, "ins", 300, ins_seq))
    for pos, kind, length, payload in sorted(ops, key=lambda o: -o[0]):
        if kind == "ins":
            q[pos:pos] = [payload] if length == 1 else list(payload)
        else:
            del q[pos:pos + (1 if length == 1 else length)]
    mutated = "".join(q)
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    rc = "".join(comp.get(c, "N") for c in reversed(mutated))

    problems = []

    def check_run(name, query, expect_strand):
        """Run the pipeline on one query contig and assert every planted event."""
        new_fa = work / f"st_new_{name}.fa"
        new_fa.write_text(f">{contig}_{name}\n{wrap(query)}\n")
        pipe = run_pipeline(old_fa, new_fa, threads)
        rows = pipe["rows"]
        if pipe["strands"].get(expect_strand, 0) < 1:
            problems.append(f"{name}: query not aligned on '{expect_strand}' strand")

        def found(typ, pred):
            return any(r[0] == contig and r[2] == typ and pred(r) for r in rows)

        for p, _old, _new in planted_snp:
            if not found("SNP", lambda r, p=p: r[1] == p):
                problems.append(f"{name}: SNP at {p} not recovered")
        for kind, a, b in hp_ops:
            typ = "INS" if kind == "INS" else "DEL"
            if not found(typ, lambda r, a=a, b=b: a - 1 <= r[1] <= b):
                problems.append(f"{name}: {typ} homopolymer [{a},{b}) not recovered")
        if not found("LDEL", lambda r, d=big_del_at: r[1] + 1 < d + 100 and d < r[1] + 1 + r[3]):
            problems.append(f"{name}: LDEL at {big_del_at} not recovered")
        if not found("LINS", lambda r, i=big_ins_at: abs(r[1] - i) <= 10):
            problems.append(f"{name}: LINS at {big_ins_at} not recovered")
        return pipe

    pipe_fwd = check_run("fwd", mutated, "+")
    pipe_rc = check_run("rc", rc, "-")
    info = {"fwd_strands": pipe_fwd["strands"], "rc_strands": pipe_rc["strands"]}
    return problems, pipe_fwd, info, planted_snp, hp_ops, big_del_at, big_ins_at


def _unit_contig_start():
    """Direct checks of contig-start and whole-contig indel handling.

    minimap2 clips/offsets indels at the very start of a contig, so these are
    checked on ``event_to_row`` directly (deterministic).
    """
    problems = []
    seq = "AAAACCCCGGGGTTTTACGTACGT"

    # insertion at the very start: right anchor is the following base (seq[0])
    row = event_to_row("c", seq, {"kind": "INS", "tpos": 0, "seq": "A"})
    if row is None or row[1] != 0 or row[4] != seq[0] or row[5] != "A" + seq[0]:
        problems.append(f"contig-start insertion mishandled: {row}")
    elif row[6] != hp_run_starting_at(seq, 0) or row[7] != str_span(seq, 0):
        problems.append("contig-start insertion hp/str not from affected coord")
    elif (row[8], row[9]) != (0, 0):
        problems.append(f"contig-start insertion interval {row[8:]} != (0, 0)")

    # deletion of the first two bases: right anchor is seq[2]
    row = event_to_row("c", seq, {"kind": "DEL", "tpos": 0, "seq": "AA"})
    if row is None or row[1] != 0 or row[4] != "AAA" or row[5] != seq[2]:
        problems.append(f"contig-start deletion mishandled: {row}")
    elif (row[8], row[9]) != (0, 2):
        problems.append(f"contig-start deletion interval {row[8:]} != (0, 2)")

    # a deletion that left-normalises to the start from a homopolymer
    row = event_to_row("c", seq, {"kind": "DEL", "tpos": 3, "seq": "A"})
    if row is None or row[1] != 0 or row[4] != "AA" or row[5] != seq[1]:
        problems.append(f"left-shifted-to-start deletion mishandled: {row}")
    elif (row[8], row[9]) != (0, 1):
        problems.append(f"left-shifted-to-start deletion interval {row[8:]} != (0, 1)")

    # whole-contig deletion: no anchor at all, interval ends at the contig end
    row = event_to_row("c", seq, {"kind": "DEL", "tpos": 0, "seq": seq})
    if row is None or row[1] != 0 or row[3] != len(seq) or row[4] != seq or row[5] != "":
        problems.append(f"whole-contig deletion mishandled: {row}")
    elif (row[8], row[9]) != (0, len(seq)):
        problems.append(f"whole-contig deletion interval {row[8:]} != (0, {len(seq)})")

    # a normal (non-origin) deletion keeps its affected interval
    row = event_to_row("c", seq, {"kind": "DEL", "tpos": 16, "seq": "AC"})
    if row is None or (row[8], row[9]) != (16, 18):
        problems.append(f"normal deletion interval {row[8:] if row else None} != (16, 18)")
    return problems


def _unit_str_span():
    """N (non-ACGT) must be a repeat boundary, not repeat evidence."""
    problems = []
    # a pure ACGT period-2 repeat of >= 3 copies is still detected
    if str_span("ATATATAT", 2) != 8:
        problems.append(f"ACGT STR span {str_span('ATATATAT', 2)} != 8")
    # a run of N must not register as a repeat
    if str_span("NNNNNN", 2) != 0:
        problems.append(f"N-only STR span {str_span('NNNNNN', 2)} != 0")
    # N flanking an otherwise periodic region must truncate the span
    if str_span("ATATATNNNN", 2) != 6:
        problems.append(f"N-truncated STR span {str_span('ATATATNNNN', 2)} != 6")
    # a single SNP inside an N-containing region is not repeat evidence
    if str_span("ANANAN", 1) != 0:
        problems.append(f"N-interleaved STR span {str_span('ANANAN', 1)} != 0")
    return problems


def self_test(args):
    src = _find_real_old_fa()
    if src is None:
        print("self-test: no real old.fa found under data/", file=sys.stderr)
        return 2
    seqs = read_fasta(src)
    # pick the longest contig and slice 200 kb from its middle
    contig = max(seqs, key=lambda c: len(seqs[c]))
    full = seqs[contig]
    W = 200000
    if len(full) < W + 2000:
        print(f"self-test: contig {contig} too short ({len(full)})", file=sys.stderr)
        return 2
    start = (len(full) - W) // 2
    slice_seq = full[start:start + W]

    work = ROOT / "scratch" / "qc" / "selftest"
    work.mkdir(parents=True, exist_ok=True)
    problems, pipe, info, snps, hp_ops, bd, bi = _plant_and_check(slice_seq, contig, work, args.threads)
    problems = problems + _unit_contig_start() + _unit_str_span()
    if pipe is None:
        print("SELF-TEST FAILED:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 1

    # write the truth of the synthetic case for inspection
    res = {"self_test": True, "source_old_fa": str(src), "contig": contig,
           "slice_start": start, "slice_len": W,
           "alignments_used": pipe["alignments_used"], "strands": info,
           "planted_snps": [p for p, _o, _n in snps],
           "planted_hp_ops": hp_ops, "big_del_at": bd, "big_ins_at": bi,
           "problems": problems}
    s = summarise(pipe["rows"], pipe["regions"], pipe["contig_lengths"])
    res.update({k: s[k] for k in ("counts", "unaligned_bp", "multi_bp", "total_bp")})
    write_case(work / "truth", res, pipe["rows"], pipe["regions"])

    if problems:
        print("SELF-TEST FAILED:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 1
    print(f"self-test OK: source={src} contig={contig} slice={start}..{start + W} "
          f"snps={len(snps)} hp_ops={len(hp_ops)} large=2 contig_start=ok "
          f"fwd_strands={info['fwd_strands']} rc_strands={info['rc_strands']} "
          f"counts={s['counts']}")
    return 0


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sets", default="m6large,real,val")
    ap.add_argument("--cases", default=None, help="comma-separated case substrings")
    ap.add_argument("--out", default=str(ROOT / "scratch" / "qc" / "truth"))
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test(args)

    out_root = Path(args.out)
    sets = [s.strip() for s in args.sets.split(",") if s.strip()]
    only = [s for s in args.cases.split(",") if s] if args.cases else None

    agg = {}
    skipped = []
    for set_name in sets:
        set_dir = SET_DIRS.get(set_name)
        if set_dir is None:
            print(f"unknown set: {set_name}", file=sys.stderr)
            continue
        manifest = read_manifest(set_dir)
        cases = list(manifest) if manifest else \
            [p.name for p in sorted(set_dir.iterdir()) if p.is_dir()] if set_dir.exists() else []
        for case in cases:
            if only and not any(s in case for s in only):
                continue
            row = manifest.get(case, {})
            status, res = build_case(set_name, case, set_dir / case, row,
                                     out_root, args.threads, args.force)
            if status == "missing":
                skipped.append(f"{set_name}/{case}: {res['error']}")
                print(f"skip {set_name}/{case}: {res['error']}")
                continue
            a = agg.setdefault(set_name, {"cases": 0, "counts": {t: 0 for t in TYPES},
                                          "unaligned_bp": 0, "multi_bp": 0,
                                          "total_bp": 0, "runtime_s": 0.0})
            a["cases"] += 1
            for t in TYPES:
                a["counts"][t] += res.get("counts", {}).get(t, 0)
            a["unaligned_bp"] += res.get("unaligned_bp", 0)
            a["multi_bp"] += res.get("multi_bp", 0)
            a["total_bp"] += res.get("total_bp", 0)
            a["runtime_s"] += res.get("runtime_s", 0.0)
            print(f"{status:6s} {set_name}/{case}: "
                  f"SNP={res.get('counts', {}).get('SNP', 0)} "
                  f"INS={res.get('counts', {}).get('INS', 0)} "
                  f"DEL={res.get('counts', {}).get('DEL', 0)} "
                  f"LINS={res.get('counts', {}).get('LINS', 0)} "
                  f"LDEL={res.get('counts', {}).get('LDEL', 0)} "
                  f"unaligned={res.get('unaligned_bp', 0)} "
                  f"multi={res.get('multi_bp', 0)} "
                  f"t={res.get('runtime_s', 0)}s")

    print("\n=== truth summary ===")
    print(f"{'set':8s} {'cases':>5s} {'SNP':>8s} {'INS':>7s} {'DEL':>7s} "
          f"{'LINS':>6s} {'LDEL':>6s} {'unaligned_bp':>13s} {'multi_bp':>10s} "
          f"{'total_bp':>12s} {'runtime_s':>9s}")
    for set_name in sets:
        a = agg.get(set_name)
        if not a:
            continue
        print(f"{set_name:8s} {a['cases']:5d} {a['counts']['SNP']:8d} {a['counts']['INS']:7d} "
              f"{a['counts']['DEL']:7d} {a['counts']['LINS']:6d} {a['counts']['LDEL']:6d} "
              f"{a['unaligned_bp']:13d} {a['multi_bp']:10d} {a['total_bp']:12d} "
              f"{a['runtime_s']:9.1f}")
    if skipped:
        print(f"\nskipped {len(skipped)} case(s):")
        for s in skipped:
            print("  " + s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
