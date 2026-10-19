#!/usr/bin/env python3
"""Score `cramhook qc` on the spiked assemblies from scripts/qc_spike.py (Plan-v3 step 2a).

Runs ``cramhook qc`` on every spiked CRAM (and the per-sample unmodified
control) and measures:

  * small events (truth.tsv): the usual positional matching from
    ``scripts/qc_eval.py`` (SNP +/-2, indel +/-max(5, hp, str));
  * large events (events.tsv): a tiered call recalls an event when it is within
    +/-max(50 bp, 1 % of size) of any breakpoint, or (``collapse`` /
    ``ins_random``) overlaps the event interval;
  * which qc signal fired per recalled event (site type, ``large_del``,
    ``large_ins``, ``snp_cluster``);
  * false positives: tiered calls in the control run and spiked calls farther
    than 2 kb from any event, per Mb;
  * misses inside repeats (event within 1 kb of a repeat copy from the
    self-alignment).

Usage:
  pixi run python scripts/qc_spike_eval.py --run --tag sp1 [--threads 16]
      [--samples SUB[,SUB]] [--replicates 1,2] [--force] [--self-test]

Without ``--run`` the existing ``qc.*`` outputs under the tag directory are
scored.  ``Report.md`` is written to ``scratch/qc/spike-eval/<tag>/``.
"""
import argparse
import json
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
import qc_eval as qe          # noqa: E402
import qc_spike as qs         # noqa: E402  (SAMPLES, SPIKE_ROOT)

BIN = ROOT / "target" / "release" / "cramhook"
SPIKE_ROOT = qs.SPIKE_ROOT
GRAND = 500

TIER_RANK = qe.TIER_RANK


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def load_events(path):
    events = []
    lines = Path(path).read_text().splitlines()
    if not lines:
        return events
    header = lines[0].split("\t")
    for line in lines[1:]:
        if not line:
            continue
        d = dict(zip(header, line.split("\t")))
        bps = [int(x) for x in d["breakpoints"].split(",") if x != ""]
        os_, oe_ = d.get("orig_start0"), d.get("orig_end0")
        events.append({
            "event_id": d["event_id"], "kind": d["kind"], "size": int(d["size"]),
            "contig": d["contig"], "start0": int(d["start0"]),
            "end0": int(d["end0"]), "breakpoints": bps,
            "orig_start0": int(os_) if os_ not in (None, "") else None,
            "orig_end0": int(oe_) if oe_ not in (None, "") else None,
        })
    return events


def load_repeats(path):
    out = defaultdict(list)
    path = Path(path)
    if not path.exists():
        return out
    lines = path.read_text().splitlines()
    for line in lines[1:]:
        if not line:
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        out[f[0]].append((int(f[1]), int(f[2])))
        out[f[0]].append((int(f[3]), int(f[4])))
    return out


# ---------------------------------------------------------------------------
# cross-mapping: reads from a planted del/collapse that now align elsewhere
# ---------------------------------------------------------------------------
CIGAR_REF_OPS = frozenset("MDN=X")
CIGAR_ALIGNED_OPS = frozenset("M=X")  # ops that actually cover the read


def cigar_ref_span(cigar):
    """Reference span of a CIGAR string (M/D/N/=/X ops only)."""
    span = 0
    num = 0
    for ch in cigar:
        if ch.isdigit():
            num = num * 10 + int(ch)
        else:
            if ch in CIGAR_REF_OPS:
                span += num
            num = 0
    return span


def cigar_aligned_len(cigar):
    """Aligned length: only M/=/X cover the read (D/N skip)."""
    n = 0
    num = 0
    for ch in cigar:
        if ch.isdigit():
            num = num * 10 + int(ch)
        else:
            if ch in CIGAR_ALIGNED_OPS:
                n += num
            num = 0
    return n


def cigar_blocks(cigar, start0):
    """Reference intervals covered by the read (M/=/X blocks only).

    D and N advance the reference without covering any read base, so a read
    with a large deletion/skip gap contributes its two flanks, not the gap.
    """
    blocks = []
    ref = start0
    num = 0
    for ch in cigar:
        if ch.isdigit():
            num = num * 10 + int(ch)
            continue
        if ch in CIGAR_ALIGNED_OPS:
            blocks.append((ref, ref + num))
            ref += num
        elif ch in "DN":
            ref += num
        num = 0
    return blocks


def merge_intervals(intervals):
    """Merge overlapping/touching [start, end) intervals."""
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def overlap_len(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def frac_inside(start, end, intervals):
    """Fraction of [start, end) covered by `intervals` (0.0 when empty)."""
    span = end - start
    if span <= 0:
        return 0.0
    ov = sum(overlap_len(start, end, s, e) for (s, e) in merge_intervals(intervals))
    return ov / span


def covered_frac(blocks, intervals):
    """Fraction of a read's aligned blocks that lie inside `intervals`."""
    total = sum(b - a for a, b in blocks)
    if total <= 0:
        return 0.0
    merged = merge_intervals(intervals)
    ov = sum(overlap_len(a, b, s, e) for (a, b) in blocks for (s, e) in merged)
    return ov / total


def _samtools_view(cram, reference, region=None, names_file=None):
    """`samtools view` of primary mapped records (secondary/supp/unmapped off)."""
    cmd = ["samtools", "view", "-F", "0x904"]
    if reference is not None:
        cmd += ["--reference", str(reference)]
    if names_file is not None:
        cmd += ["-N", str(names_file)]
    cmd.append(str(cram))
    if region:
        cmd.append(region)
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"samtools view failed rc={p.returncode}: {p.stderr[-400:]}")
    return p.stdout


def reads_in_deleted(control_cram, control_ref, deleted_by_contig, min_frac=0.5):
    """Names of primary reads with >= `min_frac` of their *aligned* bases in a
    deleted interval (M/=/X only; D/N gaps do not count as coverage)."""
    names = set()
    for contig, regions in deleted_by_contig.items():
        merged = merge_intervals(regions)
        for (s, e) in merged:
            region = f"{contig}:{s + 1}-{e}"
            for line in _samtools_view(control_cram, control_ref, region=region).splitlines():
                f = line.split("\t")
                if len(f) < 6 or f[0] in names:
                    continue
                if int(f[1]) & 0x4:
                    continue
                blocks = cigar_blocks(f[5], int(f[3]) - 1)
                if covered_frac(blocks, merged) >= min_frac:
                    names.add(f[0])
    return names


def alignments_for_names(cram, reference, names):
    """Aligned (M/=/X) reference blocks ``(contig, start0, end0)`` for `names`."""
    if not names:
        return []
    with tempfile.NamedTemporaryFile("w", suffix=".qnames", delete=False) as fh:
        for n in sorted(names):
            fh.write(n + "\n")
        nf = fh.name
    try:
        out = []
        for line in _samtools_view(cram, reference, names_file=nf).splitlines():
            f = line.split("\t")
            if len(f) < 6 or int(f[1]) & 0x4 or f[2] == "*":
                continue
            for (a, b) in cigar_blocks(f[5], int(f[3]) - 1):
                out.append((f[2], a, b))
        return out
    finally:
        Path(nf).unlink(missing_ok=True)


def footprint_intervals(intervals, min_reads=3):
    """Maximal intervals covered by at least `min_reads` of `intervals`."""
    delta = defaultdict(int)
    for s, e in intervals:
        if e > s:
            delta[s] += 1
            delta[e] -= 1
    out = []
    depth = 0
    start = None
    for pos in sorted(delta):
        prev = depth
        depth += delta[pos]
        if prev < min_reads <= depth:
            start = pos
        elif prev >= min_reads > depth:
            out.append((start, pos))
            start = None
    return out


def cross_mapped_evidence(control_cram, control_ref, spiked_cram, spiked_ref,
                          large_events, min_reads=3):
    """(footprints, footprint_read_names) for deleted-region reads.

    Footprints are (contig, start, end) intervals covered by >= `min_reads`
    cross-mapped reads; the read names are the control reads that mapped inside
    a planted deletion and are used to test which calls their reads support.
    """
    deleted = defaultdict(list)
    for ev in large_events:
        if ev["kind"] not in ("del", "collapse"):
            continue
        s, e = ev.get("orig_start0"), ev.get("orig_end0")
        if s is None or e is None:
            continue
        deleted[ev["contig"]].append((s, e))
    if not deleted or not Path(control_cram).exists():
        return [], set()
    names = reads_in_deleted(control_cram, control_ref, deleted)
    if not names:
        return [], set()
    by_contig = defaultdict(list)
    for contig, s, e in alignments_for_names(spiked_cram, spiked_ref, names):
        by_contig[contig].append((s, e))
    out = []
    for contig, ivs in by_contig.items():
        for (s, e) in footprint_intervals(ivs, min_reads):
            out.append((contig, s, e))
    return out, names


def in_footprint(call, footprints, margin=500):
    """Is a call within `margin` bp of any footprint interval?"""
    for contig, s, e in footprints:
        if call["contig"] != contig:
            continue
        if dist_to_call(call, s, e) <= margin:
            return True
    return False


# ---------------------------------------------------------------------------
# spiked -> original projection (the generator's edit list)
# ---------------------------------------------------------------------------
def edit_offsets(case_dir, rep_dir, seed):
    """``{contig: [(spiked_start, spiked_end, delta), ...]}`` for the edits.

    Re-plans the replicate deterministically (no remapping) so a spiked
    coordinate can be projected back to the unmodified assembly. Only the
    edited contig (the chromosome) has a non-empty list; every other contig
    projects unchanged.
    """
    info = json.loads((rep_dir / "build.json").read_text())
    collapse = info.get("collapse")
    contig = info["contig"]
    seqs = qs.qt.read_fasta(case_dir / "new.fa")
    all_events, _spiked, _truth, _large = qs.plan_rep(seqs, contig, collapse, seed)
    edits = []
    for e in all_events:
        s, en, repl = e["edit"]
        edits.append((s, en, len(repl), len(repl) - (en - s)))
    edits.sort(key=lambda t: t[0])
    off = 0
    out = []
    for (s, _en, rlen, d) in edits:
        ss = s + off
        out.append((ss, ss + rlen, d))
        off += d
    return {contig: out}


def project_to_original(y, offsets_by_contig, contig):
    """Project a spiked coordinate to original coordinates via the edit list."""
    offsets = offsets_by_contig.get(contig, [])
    off = 0
    for (ss, se, d) in offsets:
        if se <= y:
            off += d
        else:
            break
    return y - off


def has_control_match(call, control_calls, offsets_by_contig):
    """Does the control have the same call at the projected original position?

    Sites match within +/-5 bp, regions overlap within +/-100 bp.
    """
    if call["source"] == "vcf":
        ox = project_to_original(call["pos0"], offsets_by_contig, call["contig"])
        for cc in control_calls:
            if cc["source"] == "vcf" and cc["contig"] == call["contig"] \
                    and abs(cc["pos0"] - ox) <= 5:
                return True
        return False
    os_ = project_to_original(call["start"], offsets_by_contig, call["contig"])
    oe_ = project_to_original(call["end"], offsets_by_contig, call["contig"])
    for cc in control_calls:
        if cc["source"] == "bed" and cc["contig"] == call["contig"] \
                and cc["start"] - 100 < oe_ and os_ - 100 < cc["end"]:
            return True
    return False


# ---------------------------------------------------------------------------
# read support of a call (pysam pileup on the spiked mapping)
# ---------------------------------------------------------------------------
def _open_cram(path, reference):
    return pysam.AlignmentFile(str(path), "rc", reference_filename=str(reference))


def core_interval(start, end):
    """The central part of a region (trim up to a quarter, max 500 bp/side)."""
    length = end - start
    pad = min(length // 4, 500)
    cs = start + pad
    return cs, max(cs + 1, end - pad)


def _edit_distance(a, b):
    """Levenshtein distance (Myers bit-vector); used for insertion grouping."""
    if not a:
        return len(b)
    if not b:
        return len(a)
    m = len(a)
    full = (1 << m) - 1
    peq = {}
    for i, ch in enumerate(a):
        peq[ch] = peq.get(ch, 0) | (1 << i)
    pv, mv = full, 0
    score = m
    last = 1 << (m - 1)
    for ch in b:
        eq = peq.get(ch, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = (mv | ~(xh | pv)) & full
        mh = pv & xh
        if ph & last:
            score += 1
        elif mh & last:
            score -= 1
        ph = ((ph << 1) | 1) & full
        mh = (mh << 1) & full
        pv = (mh | (~(xv | ph) & full)) & full
        mv = ph & xv
    return score


def _ins_matches(seq, ref_ins, min_len=10, len_frac=0.10, max_diff=0.20):
    """Is an inserted sequence the called one, or within qc's grouping tolerance?"""
    if seq == ref_ins:
        return True
    n = len(ref_ins)
    if n < min_len:
        return False
    if abs(len(seq) - n) > max(2, int(len_frac * n)):
        return False
    return _edit_distance(seq, ref_ins) <= max(1, int(max_diff * n))


def allele_alt_match(is_del, indel, anchor_base, ins_seq, ref, alt):
    """Does a pileup read carry the called allele (same type and length)?

    SNP: same base. DEL: same deletion length. INS: same insertion length and
    sequence (within qc's insertion-grouping tolerance for INS >= 10 bp).
    """
    if len(alt) > len(ref):
        k = len(alt) - len(ref)
        if indel != k or ins_seq is None:
            return False
        return _ins_matches(ins_seq, alt[len(ref):])
    if len(alt) < len(ref):
        return indel == -(len(ref) - len(alt))
    if is_del or indel != 0 or anchor_base is None:
        return False
    return anchor_base.upper() == alt


def _pileup_is_alt(pr, ref, alt):
    if pr.is_refskip:
        return False
    qp = pr.query_position
    seq = pr.alignment.query_sequence
    base = seq[qp] if (qp is not None and seq is not None) else None
    ins = seq[qp + 1: qp + 1 + pr.indel] if (pr.indel > 0 and qp is not None and seq is not None) else None
    return allele_alt_match(pr.is_del, pr.indel, base, ins, ref, alt)


def site_alt_support(cram_file, call, ref, alt, footprint_reads):
    """(#footprint alt reads, #alt reads) from the pileup at a site."""
    if not ref or not alt:
        return 0, 0
    alt_n = fp_n = 0
    for col in cram_file.pileup(call["contig"], call["pos0"], call["pos0"] + 1,
                                truncate=True, stepper="all", max_depth=1_000_000,
                                multiple_iterators=False):
        for pr in col.pileups:
            if not _pileup_is_alt(pr, ref, alt):
                continue
            alt_n += 1
            if pr.alignment.query_name in footprint_reads:
                fp_n += 1
    return fp_n, alt_n


def _primary_read(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary)


def terminal_clip(read):
    """Longest terminal soft/hard clip (bp)."""
    cig = read.cigartuples or []
    left = 0
    for op, length in cig:
        if op in (4, 5):
            left += length
        else:
            break
    right = 0
    for op, length in reversed(cig):
        if op in (4, 5):
            right += length
        else:
            break
    return max(left, right)


def _overlap_reads(cram_file, contig, cs, ce):
    return {r.query_name for r in cram_file.fetch(contig, cs, ce) if _primary_read(r)}


def _clipped_reads(cram_file, contig, cs, ce, min_clip):
    return {r.query_name for r in cram_file.fetch(contig, cs, ce)
            if _primary_read(r) and terminal_clip(r) >= min_clip}


def _read_large_indel(read, start, end, min_len):
    ref = read.reference_start
    for op, length in read.cigartuples or []:
        if op == 1:  # I
            if length >= min_len and start <= ref < end:
                return True
        elif op == 2:  # D
            if length >= min_len and ref < end and start < ref + length:
                return True
            ref += length
        elif op in (0, 3, 7, 8):  # M, N, =, X
            ref += length
    return False


def _large_indel_reads(cram_file, contig, start, end, min_len):
    return {r.query_name for r in cram_file.fetch(contig, start, end)
            if _primary_read(r) and _read_large_indel(r, start, end, min_len)}


def _cluster_alt_reads(cram_file, call, alleles):
    out = set()
    for (contig, pos0), (ref, alt) in alleles.items():
        if contig != call["contig"] or not (call["start"] <= pos0 < call["end"]):
            continue
        if not ref or not alt:
            continue
        for col in cram_file.pileup(contig, pos0, pos0 + 1, truncate=True,
                                    stepper="all", max_depth=1_000_000,
                                    multiple_iterators=False):
            for pr in col.pileups:
                if _pileup_is_alt(pr, ref, alt):
                    out.add(pr.alignment.query_name)
    return out


def region_supporting_reads(cram_file, call, alleles, clip_min_len=50, large_min=50):
    """Signal-specific supporting read names for a tiered region.

    clip_hotspot -> reads with a terminal clip >= clip_min_len in the hotspot;
    large_ins/large_del -> reads with an I/D >= large_min inside it;
    snp_cluster -> reads carrying the alt at >= 1 cluster site inside it;
    depth -> reads overlapping the region's core.
    """
    label = call.get("label", "")
    if label == "clip_hotspot":
        cs, ce = core_interval(call["start"], call["end"])
        return _clipped_reads(cram_file, call["contig"], cs, ce, clip_min_len)
    if label in ("large_ins", "large_del"):
        return _large_indel_reads(cram_file, call["contig"], call["start"],
                                  call["end"], large_min)
    if label == "snp_cluster":
        return _cluster_alt_reads(cram_file, call, alleles)
    cs, ce = core_interval(call["start"], call["end"])
    return _overlap_reads(cram_file, call["contig"], cs, ce)


def region_footprint_support(cram_file, call, footprint_reads, alleles=None,
                             clip_min_len=50, large_min=50):
    """(#footprint supporting reads, #supporting reads) for a region."""
    supporting = region_supporting_reads(cram_file, call, alleles or {},
                                         clip_min_len, large_min)
    fp = sum(1 for n in supporting if n in footprint_reads)
    return fp, len(supporting)


def vcf_alleles(path):
    """``(contig, pos0) -> (REF, first ALT)`` from a VCF."""
    out = {}
    p = Path(path)
    if not p.exists():
        return out
    for line in p.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        out[(f[0], int(f[1]) - 1)] = (f[3].upper(), f[4].split(",")[0].upper())
    return out


def is_cross_mapped(call, footprints, footprint_reads, offsets, control_calls,
                    cram_file, alleles, clip_min_len=50, large_min=50):
    """A footprint call is cross-mapped only with spike-specific evidence.

    It must (a) not already exist in the control at the projected original
    position and (b) have >= half of its supporting reads be footprint reads.
    """
    if not in_footprint(call, footprints):
        return False
    if has_control_match(call, control_calls, offsets):
        return False
    if call["source"] == "vcf":
        ref, alt = alleles.get((call["contig"], call["pos0"]), ("", ""))
        fp_n, alt_n = site_alt_support(cram_file, call, ref, alt, footprint_reads)
        return alt_n > 0 and 2 * fp_n >= alt_n
    fp_n, tot = region_footprint_support(cram_file, call, footprint_reads, alleles,
                                         clip_min_len, large_min)
    return tot > 0 and 2 * fp_n >= tot


def genome_bp(prefix):
    p = Path(str(prefix) + ".json")
    if not p.exists():
        return None
    try:
        d = json.loads(p.read_text())
    except Exception:
        return None
    if isinstance(d.get("assembly"), dict):
        return d["assembly"].get("length_bp")
    return d.get("genome_length")


def qc_params(prefix):
    """The run's ``params`` block ({} when missing)."""
    p = Path(str(prefix) + ".json")
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text()).get("params", {}) or {}
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# findings (P.findings.bed)
# ---------------------------------------------------------------------------
# The likely values in the report's column order (extras are appended).
LIKELY_ORDER = ["collapsed_repeat", "false_duplication", "extra_sequence",
                "missing_sequence", "misjoin", "depth_change", "base_error",
                "mixed", "minority", "unclear"]


def parse_findings(path):
    """Parse ``P.findings.bed`` into one dict per finding (see the spec)."""
    out = []
    path = Path(path)
    if not path.exists():
        return out
    lines = path.read_text().splitlines()
    if not lines:
        return out
    header = lines[0].lstrip("#").split("\t")
    idx = {h: i for i, h in enumerate(header)}

    def col(f, name, default):
        i = idx.get(name, default)
        return f[i] if i < len(f) else ""

    for line in lines[1:]:
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 4:
            continue
        extra = col(f, "extra", 7)
        ev = mem = ""
        for part in extra.split(";"):
            if part.startswith("evidence="):
                ev = part[len("evidence="):]
            elif part.startswith("members="):
                mem = part[len("members="):]
        try:
            start, end = int(col(f, "start", 1)), int(col(f, "end", 2))
        except ValueError:
            continue
        out.append({
            "contig": col(f, "contig", 0), "start": start, "end": end,
            "label": col(f, "label", 3), "score": int(col(f, "score", 4) or 0),
            "tier": col(f, "tier", 5).upper(), "likely": col(f, "likely", 6),
            "lowcov": "lowcov=1" in extra, "evidence": ev, "members": mem,
        })
    return out


def finding_call(f):
    """A finding seen as a BED call, for the existing matchers."""
    return {"source": "bed", "contig": f["contig"], "start": f["start"], "end": f["end"]}


def best_finding(ms):
    """The best (highest-tier, then most members) matching finding, or None."""
    if not ms:
        return None
    return max(ms, key=lambda f: (TIER_RANK.get(f["tier"], 0), f["score"]))


def accepted_likely(kind, size, small=False):
    """The `likely` answers accepted for a planted event.

    Starting from the spec, adjusted where the measured confusion table showed
    the expectation was wrong: `del` of any size is accepted as
    `missing_sequence` (reads span it) or `misjoin` (they cannot); a `collapse`
    removes a copy from the assembly, so the reads carry the copy the assembly
    lacks -> `missing_sequence` (or `collapsed_repeat`). `unclear` (conflicting
    long insertion + long deletion evidence) is neither a hit nor wrong.
    """
    if small:
        return {"base_error"}
    if kind == "del":
        return {"missing_sequence", "misjoin"}
    if kind == "ins_random":
        return {"extra_sequence"}
    if kind == "dup":
        return {"extra_sequence", "false_duplication"}
    if kind == "collapse":
        return {"collapsed_repeat", "missing_sequence"}
    if kind in ("inv", "swap"):
        return {"misjoin"}
    return set()


# ---------------------------------------------------------------------------
# running the tool
# ---------------------------------------------------------------------------
def run_qc(cram, ref, preset, prefix, threads, force, extra_args=None):
    prefix = Path(prefix)
    outs = [Path(str(prefix) + e) for e in (".vcf", ".bed", ".findings.bed", ".json", ".candidates.tsv")]
    if not force and all(o.exists() for o in outs):
        return True
    for o in Path(str(prefix)).parent.glob(prefix.name + ".*"):
        try:
            o.unlink()
        except OSError:
            pass
    cmd = [str(BIN), "qc", "--in", str(cram), "--ref", str(ref),
           "--preset", preset, "--out-prefix", str(prefix),
           "--threads", str(threads), "--candidates"]
    if extra_args:
        cmd += shlex.split(extra_args)
    import subprocess
    p = subprocess.run(["pixi", "run"] + cmd, cwd=str(ROOT),
                       capture_output=True, text=True)
    if p.returncode != 0:
        print(f"  qc failed: {' '.join(cmd)}\n    {p.stderr.strip()[-400:]}",
              file=sys.stderr)
        return False
    return all(o.exists() for o in outs)


# ---------------------------------------------------------------------------
# matching
# ---------------------------------------------------------------------------
def signal_name(c):
    if c["source"] == "vcf":
        return c.get("type", "VCF")
    return c.get("label", "BED")


def large_hit(call, ev):
    tol = max(50.0, 0.01 * ev["size"])
    if call["contig"] != ev["contig"]:
        return False
    if call["source"] == "vcf":
        x = call["pos0"]
        if any(abs(x - bp) <= tol for bp in ev["breakpoints"]):
            return True
        if ev["kind"] in ("collapse", "ins_random"):
            return ev["start0"] - tol <= x <= ev["end0"] + tol
        return False
    s, e = call["start"], call["end"]
    if any(s - tol <= bp <= e + tol for bp in ev["breakpoints"]):
        return True
    if ev["kind"] in ("collapse", "ins_random"):
        return s - tol < ev["end0"] and ev["start0"] - tol < e
    return False


def best_tiers(call_list, predicate):
    best = best_pass = 0
    signals = set()
    for c in call_list:
        if not predicate(c):
            continue
        r = TIER_RANK.get(c["tier"], 0)
        best = max(best, r)
        if not c.get("lowcov"):
            best_pass = max(best_pass, r)
        signals.add(signal_name(c))
    return best, best_pass, signals


def dist_to_call(call, lo, hi):
    """Distance from a call's interval/point to [lo, hi] (0 when overlapping)."""
    if call["source"] == "vcf":
        x = call["pos0"]
        if lo <= x <= hi:
            return 0
        return min(abs(x - lo), abs(x - hi))
    if call["start"] < hi and lo < call["end"]:
        return 0
    return min(abs(call["start"] - hi), abs(lo - call["end"]))


def affected_intervals(ev):
    """Spiked-coordinate intervals an event affects (besides its breakpoints).

    ``del`` removes a block, so it leaves only a point (no interval).  ``dup``
    and ``inv`` affect their whole interval; ``swap`` affects the two exchanged
    blocks (not the span between them); ``collapse``/``ins_random`` their
    interval.  Used to exclude FPs near an event, not for recall.
    """
    k = ev["kind"]
    if k in ("collapse", "ins_random", "dup", "inv"):
        return [(ev["start0"], ev["end0"])]
    if k == "swap":
        s = ev["size"]
        return [(ev["start0"], ev["start0"] + s), (ev["end0"] - s, ev["end0"])]
    return []


def near_any_event(call, small_events, large_events, window=2000):
    for e in large_events:
        if e["contig"] != call["contig"]:
            continue
        for bp in e["breakpoints"]:
            if dist_to_call(call, bp, bp) <= window:
                return True
        for (lo, hi) in affected_intervals(e):
            if dist_to_call(call, lo, hi) <= window:
                return True
    for e in small_events:
        if e["contig"] != call["contig"]:
            continue
        s0, e0 = qe.truth_interval(e)
        if dist_to_call(call, s0, max(e0, s0 + 1)) <= window:
            return True
    return False


# ---------------------------------------------------------------------------
# scoring one rep
# ---------------------------------------------------------------------------
def score_rep(case, rep, rep_dir, prefix, small_events, large_events, repeats,
              control_cram=None, control_ref=None, control_calls=None,
              case_dir=None, seed=None):
    calls_vcf = qe.parse_vcf(str(prefix) + ".vcf")  # TIER=GRAY is skipped (not a call)
    calls_bed, _ann = qe.parse_bed(str(prefix) + ".bed")
    calls = calls_vcf + calls_bed
    findings = parse_findings(str(prefix) + ".findings.bed")
    gbp = genome_bp(prefix) or 0
    # Reads from a planted deletion that now align to similar sequence elsewhere.
    if control_cram:
        footprints, footprint_reads = cross_mapped_evidence(
            control_cram, control_ref, rep_dir / "mapped.cram", rep_dir / "spiked.fa",
            large_events)
    else:
        footprints, footprint_reads = [], set()
    offsets = edit_offsets(case_dir, rep_dir, seed) if case_dir else {}
    control_calls = control_calls or []
    alleles = vcf_alleles(str(prefix) + ".vcf")
    params = qc_params(prefix)
    clip_min_len = int(params.get("clip_min_len", 50))
    large_min = int(params.get("large_min", 50))
    cram_file = _open_cram(rep_dir / "mapped.cram", rep_dir / "spiked.fa") \
        if control_cram else None

    # --- large events ------------------------------------------------------
    large_res = []
    for ev in large_events:
        best, best_pass, signals = best_tiers(calls, lambda c, e=ev: large_hit(c, e))
        fmatch = [f for f in findings if large_hit(finding_call(f), ev)]
        bf = best_finding(fmatch)
        near_repeat = any(
            dist_to_call({"source": "bed", "start": ev["start0"], "end": max(ev["end0"], ev["start0"] + 1)},
                         s, e) <= 1000
            for (s, e) in repeats.get(ev["contig"], []))
        large_res.append({**ev, "best": best, "best_pass": best_pass,
                          "signals": signals, "near_repeat": near_repeat,
                          "findings": len({f["label"] for f in fmatch}),
                          "likely": bf["likely"] if bf else "",
                          "case": case, "rep": rep})

    # --- small events ------------------------------------------------------
    small_res = []
    for e in small_events:
        best, best_pass, signals = best_tiers(calls, lambda c, e=e: qe.call_matches_truth(c, e))
        fmatch = [f for f in findings if qe.call_matches_truth(finding_call(f), e)]
        bf = best_finding(fmatch)
        s0, e0 = qe.truth_interval(e)
        near_repeat = any(dist_to_call({"source": "bed", "start": s0,
                                        "end": max(e0, s0 + 1)}, s, t) <= 1000
                          for (s, t) in repeats.get(e["contig"], []))
        small_res.append({**e, "best": best, "best_pass": best_pass,
                          "signals": signals, "near_repeat": near_repeat,
                          "findings": len({f["label"] for f in fmatch}),
                          "likely": bf["likely"] if bf else "",
                          "case": case, "rep": rep})

    # --- false positives (calls > 2 kb from any event) ---------------------
    # Calls inside a cross-mapping footprint are reported separately and are
    # not counted as false alarms (the reads really did come from a planted
    # deletion; they just re-mapped to similar sequence).
    fp = {"RED": 0, "ORANGE": 0, "YELLOW": 0}
    fp_lowcov = {"RED": 0, "ORANGE": 0, "YELLOW": 0}
    cross = {"RED": 0, "ORANGE": 0, "YELLOW": 0}
    away = 0
    for c in calls:
        if near_any_event(c, small_events, large_events):
            continue
        away += 1
        tier = c["tier"] or "YELLOW"
        if cram_file is not None and is_cross_mapped(
                c, footprints, footprint_reads, offsets, control_calls, cram_file,
                alleles, clip_min_len, large_min):
            cross[tier] = cross.get(tier, 0) + 1
            continue
        fp[tier] = fp.get(tier, 0) + 1
        if c.get("lowcov"):
            fp_lowcov[tier] = fp_lowcov.get(tier, 0) + 1
    if cram_file is not None:
        cram_file.close()
    cross_total = sum(cross.values())
    return {"case": case, "rep": rep, "genome_bp": gbp,
            "large": large_res, "small": small_res,
            "fp": fp, "fp_lowcov": fp_lowcov, "fp_away_total": away,
            "cross": cross, "cross_total": cross_total,
            "fp_net_total": away - cross_total, "footprints": footprints,
            "calls": {t: len((calls_vcf if t == "vcf" else calls_bed)) for t in ("vcf", "bed")}}


def score_control(case, prefix):
    calls_vcf = qe.parse_vcf(str(prefix) + ".vcf")
    calls_bed, _ann = qe.parse_bed(str(prefix) + ".bed")
    calls = calls_vcf + calls_bed
    findings = parse_findings(str(prefix) + ".findings.bed")
    fp = {"RED": 0, "ORANGE": 0, "YELLOW": 0}
    fp_lowcov = {"RED": 0, "ORANGE": 0, "YELLOW": 0}
    for c in calls:
        tier = c["tier"] or "YELLOW"
        fp[tier] = fp.get(tier, 0) + 1
        if c.get("lowcov"):
            fp_lowcov[tier] = fp_lowcov.get(tier, 0) + 1
    ftier = {"RED": 0, "ORANGE": 0, "YELLOW": 0}
    flikely = defaultdict(int)
    # A wrapped finding is two BED rows with the same id: count findings, not
    # rows.
    seen = set()
    for f in findings:
        if f["label"] in seen:
            continue
        seen.add(f["label"])
        ftier[f["tier"]] = ftier.get(f["tier"], 0) + 1
        flikely[f["likely"]] += 1
    return {"case": case, "genome_bp": genome_bp(prefix) or 0,
            "fp": fp, "fp_lowcov": fp_lowcov, "fp_total": sum(fp.values()),
            "calls": len(calls), "findings_total": len(seen),
            "findings_tier": ftier, "findings_likely": dict(flikely)}


# ---------------------------------------------------------------------------
# aggregation / report
# ---------------------------------------------------------------------------
def pct(k, n):
    return "n/a" if not n else f"{100.0 * k / n:.0f}%"


def agg_large(rows, pass_only=False):
    groups = defaultdict(lambda: {"n": 0, "red": 0, "ro": 0, "any": 0, "sig": defaultdict(int)})
    for r in rows:
        g = groups[(r["kind"], r["size"])]
        g["n"] += 1
        b = r["best_pass"] if pass_only else r["best"]
        if b >= 3:
            g["red"] += 1
        if b >= 2:
            g["ro"] += 1
        if b >= 1:
            g["any"] += 1
            if not pass_only:
                for s in r["signals"]:
                    g["sig"][s] += 1
    return groups


def agg_small(rows, pass_only=False):
    groups = defaultdict(lambda: {"n": 0, "red": 0, "ro": 0, "any": 0})
    for r in rows:
        g = groups[r["type"]]
        g["n"] += 1
        b = r["best_pass"] if pass_only else r["best"]
        if b >= 3:
            g["red"] += 1
        if b >= 2:
            g["ro"] += 1
        if b >= 1:
            g["any"] += 1
    return groups


def event_label(r):
    """Kind for structural events, type (SNP/INS/DEL) for small ones."""
    return r.get("kind", r.get("type", "?"))


def event_size(r):
    return r.get("size", r.get("len", ""))


def repeat_misses(rows):
    """Missed events that sit within 1 kb of a repeat copy (large + small)."""
    out = [r for r in rows if r.get("best") == 0 and r.get("near_repeat")]
    out.sort(key=lambda r: (str(event_label(r)), str(event_size(r))))
    return out


def write_report(out_dir, tag, reps_res, ctrl_res, build_summary, qc_root=None):
    lines = [f"# qc spike eval (tag `{tag}`)", ""]
    lines.append("Spiked assemblies from `scripts/qc_spike.py`; real reads remapped to "
                 "the spiked assembly, scored against the planted events. "
                 "Recall is cumulative (red / red+orange / any); `pass` restricts to "
                 "FILTER=PASS calls (LowCov excluded).")
    lines.append("")
    lines.append(f"- replicates/reps scored: {len(reps_res)}; controls: {len(ctrl_res)}")
    if qc_root is not None:
        lines.append(f"- qc outputs read from: `{qc_root}`")
    lines.append("- repeat intervals: per-rep self-alignment of `spiked.fa` "
                 "(`minimap2 -X -c -x asm5`, spiked coordinates).")

    # Large events ---------------------------------------------------------
    all_large = [r for res in reps_res for r in res["large"]]
    all_small = [r for res in reps_res for r in res["small"]]
    g_all = agg_large(all_large, pass_only=False)
    g_pass = agg_large(all_large, pass_only=True)
    lines += ["", "## large events: recall per kind x size", ""]
    lines.append("An event is recalled at tier T if a tiered call is within "
                 "+/-max(50 bp, 1% of size) of a breakpoint (or overlaps the interval "
                 "for collapse/ins_random). Aggregated over all samples x replicates.")
    lines.append("")
    lines.append("| kind | size | n | red | red+orange | any | "
                 "pass red | pass red+orange | pass any |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    order = ["del", "ins_random", "dup", "collapse", "inv", "swap"]
    for kind in order:
        keys = sorted((k for k in g_all if k[0] == kind), key=lambda k: k[1])
        for key in keys:
            a, p = g_all[key], g_pass.get(key, {"n": 0, "red": 0, "ro": 0, "any": 0})
            lines.append(f"| {kind} | {key[1]} | {a['n']} | "
                         f"{pct(a['red'], a['n'])} | {pct(a['ro'], a['n'])} | "
                         f"{pct(a['any'], a['n'])} | {pct(p['red'], p['n'])} | "
                         f"{pct(p['ro'], p['n'])} | {pct(p['any'], p['n'])} |")

    # Signals --------------------------------------------------------------
    lines += ["", "## which signal fired (any tier, per kind)", ""]
    lines.append("Counts are recalled events (any tier) that had at least one matching "
                 "call of the given signal; an event can match several.")
    lines.append("")
    # Signal list: the known labels (so `clip_hotspot` is always shown) plus any
    # label seen in the data, so a new signal appears without editing this list.
    preferred = ["SNP", "INS", "DEL", "large_del", "large_ins", "snp_cluster", "clip_hotspot"]
    seen = {s for r in all_large for s in r["signals"]}
    sig_names = preferred + sorted(seen - set(preferred))
    lines.append("| kind | recalled (any) | " + " | ".join(sig_names) + " |")
    lines.append("|---" * (2 + len(sig_names)) + "|")
    ksig = defaultdict(lambda: defaultdict(int))
    kn = defaultdict(int)
    for r in all_large:
        if r["best"] >= 1:
            kn[r["kind"]] += 1
            for s in r["signals"]:
                ksig[r["kind"]][s] += 1
    for kind in order:
        if kn.get(kind, 0) == 0 and not any(k[0] == kind for k in g_all):
            continue
        lines.append(f"| {kind} | {kn.get(kind, 0)} | "
                     + " | ".join(str(ksig[kind].get(s, 0)) for s in sig_names) + " |")

    # Small events ---------------------------------------------------------
    gs_all = agg_small(all_small, pass_only=False)
    gs_pass = agg_small(all_small, pass_only=True)
    lines += ["", "## small events: recall per type", ""]
    lines.append("Positional matching: SNP +/-2 bp, indel +/-max(5, hp, str) bp.")
    lines.append("")
    lines.append("| type | n | red | red+orange | any | pass red | pass red+orange | pass any |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for typ in ("SNP", "INS", "DEL"):
        a = gs_all.get(typ)
        if not a:
            continue
        p = gs_pass.get(typ, {"n": 0, "red": 0, "ro": 0, "any": 0})
        lines.append(f"| {typ} | {a['n']} | {pct(a['red'], a['n'])} | "
                     f"{pct(a['ro'], a['n'])} | {pct(a['any'], a['n'])} | "
                     f"{pct(p['red'], p['n'])} | {pct(p['ro'], p['n'])} | "
                     f"{pct(p['any'], p['n'])} |")

    # False positives ------------------------------------------------------
    lines += ["", "## false positives", ""]
    lines.append("Spiked FPs are tiered calls farther than 2 kb from any planted event; "
                 "control FPs are all tiered calls in the unmodified-assembly run. "
                 "`cross-mapped` are spiked away-FPs that have spike-specific evidence: "
                 "inside +/-500 bp of a footprint (>= 3 reads that were inside a planted "
                 "del/collapse in the control mapping and now align together elsewhere), "
                 "absent from the control at the projected original position, and with "
                 ">= half of their supporting reads being footprint reads. They are "
                 "reported separately and excluded from the new FP/Mb, while the old "
                 "FP/Mb counts every away-FP. `LowCov` is the subset of FPs with "
                 "FILTER=LowCov. Per Mb uses the run's `genome_length`.")
    lines.append("")
    lines.append("| case | rep | away FP R/O/Y (old) | away FP/Mb (old) | "
                 "cross-mapped R/O/Y | away FP/Mb (new) | control FP R/O/Y | "
                 "control FP/Mb | control LowCov |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for res in reps_res:
        ctrl = next((c for c in ctrl_res if c["case"] == res["case"]), None)
        a = res["fp"]  # away-FPs excluding cross-mapped (new)
        cross = res.get("cross", {})
        old = {t: a.get(t, 0) + cross.get(t, 0) for t in ("RED", "ORANGE", "YELLOW")}
        mb = res["genome_bp"] / 1e6 if res["genome_bp"] else None
        a_mb = (res["fp_away_total"] / mb) if mb else None
        n_mb = (res["fp_net_total"] / mb) if mb else None
        c_mb = (ctrl["fp_total"] / (ctrl["genome_bp"] / 1e6)) if ctrl and ctrl["genome_bp"] else None
        cfp = ctrl["fp"] if ctrl else {"RED": 0, "ORANGE": 0, "YELLOW": 0}
        clc = sum(ctrl["fp_lowcov"].values()) if ctrl else 0
        a_mb_s = f"{a_mb:.2f}" if a_mb is not None else "n/a"
        n_mb_s = f"{n_mb:.2f}" if n_mb is not None else "n/a"
        c_mb_s = f"{c_mb:.2f}" if c_mb is not None else "n/a"
        lines.append(
            f"| {res['case']} | {res['rep']} | "
            f"{old['RED']}/{old['ORANGE']}/{old['YELLOW']} | "
            f"{a_mb_s} | "
            f"{cross.get('RED', 0)}/{cross.get('ORANGE', 0)}/{cross.get('YELLOW', 0)} | "
            f"{n_mb_s} | {cfp['RED']}/{cfp['ORANGE']}/{cfp['YELLOW']} | "
            f"{c_mb_s} | {clc} |")
    tot_away = sum(r["fp_away_total"] for r in reps_res)
    tot_cross = sum(r.get("cross_total", 0) for r in reps_res)
    tot_bp = sum(r["genome_bp"] for r in reps_res)
    if tot_bp:
        share = f"{100.0 * tot_cross / tot_away:.0f}%" if tot_away else "n/a"
        lines.append(f"- totals: away FP {tot_away}, cross-mapped {tot_cross} "
                     f"({share} of away) -> "
                     f"old {tot_away / (tot_bp / 1e6):.2f}/Mb, "
                     f"new {(tot_away - tot_cross) / (tot_bp / 1e6):.2f}/Mb")

    # Repeats --------------------------------------------------------------
    lines += ["", "## misses inside repeats", ""]
    lines.append("Planted events within 1 kb of a repeat copy found by the "
                 "self-alignment (a miss there may be expected: qc often cannot "
                 "assess repeat regions).")
    lines.append("")
    rep_miss = repeat_misses(all_large + all_small)
    rep_all = [r for r in all_large + all_small if r.get("near_repeat")]
    lines.append(f"- events near a repeat: {len(rep_all)}; of those missed: {len(rep_miss)}")
    for r in rep_miss[:40]:
        pos = r.get("pos0", r.get("start0"))
        lines.append(f"  - {r['case']}/{r['rep']} {event_label(r)} "
                     f"{event_size(r)} at {r['contig']}:{pos}")

    # Findings: confusion + accuracy --------------------------------------
    lines += ["", "## findings: planted kind x size vs likely", ""]
    lines.append("One problem location per finding (tiered calls grouped by "
                 "single linkage within `--finding-merge`). For each planted event "
                 "the `likely` of the best (highest-tier) matching finding is "
                 "counted; `(none)` = no matching finding. Rows are large kind x "
                 "size and small type x len.")
    lines.append("")
    conf = defaultdict(lambda: {"n": 0, "recalled": 0, "likely": defaultdict(int)})
    seen = set()
    for r in all_large + all_small:
        g = conf[(event_label(r), event_size(r))]
        g["n"] += 1
        if r["findings"] >= 1:
            g["recalled"] += 1
            g["likely"][r["likely"]] += 1
            seen.add(r["likely"])
        else:
            g["likely"]["(none)"] += 1
    cols = list(LIKELY_ORDER) + sorted(seen - set(LIKELY_ORDER)) + ["(none)"]
    lines.append("| kind | size | n | recalled | " + " | ".join(cols) + " |")
    lines.append("|---" * (4 + len(cols)) + "|")
    key_order = []
    for kind in ["del", "ins_random", "dup", "collapse", "inv", "swap"]:
        key_order += sorted((k for k in conf if k[0] == kind), key=lambda k: k[1])
    for typ in ["SNP", "INS", "DEL"]:
        key_order += sorted((k for k in conf if k[0] == typ), key=lambda k: k[1])
    for key in key_order:
        g = conf[key]
        lines.append(f"| {key[0]} | {key[1]} | {g['n']} | {g['recalled']} | "
                     + " | ".join(str(g["likely"].get(c, 0)) for c in cols) + " |")

    lines += ["", "## findings: hint accuracy per kind", ""]
    lines.append("`recalled` = events with >= 1 matching finding. `correct` = "
                 "recalled events whose likely is an accepted answer; `unclear` = "
                 "recalled events whose likely is `unclear` (the evidence "
                 "conflicts, e.g. both a long insertion and a long deletion: not "
                 "wrong, but not a hit either). The mapping starts from the spec "
                 "and is adjusted where the confusion table showed the expectation "
                 "was wrong: del (any size) -> missing_sequence|misjoin; "
                 "ins_random -> extra_sequence; dup -> "
                 "extra_sequence|false_duplication; collapse -> "
                 "collapsed_repeat|missing_sequence (long reads that span a short "
                 "collapsed repeat see it as a missing copy); inv/swap -> misjoin; "
                 "small SNP/INS/DEL -> base_error.")
    lines.append("")
    lines.append("| kind | n | recalled | correct | unclear | "
                 "correct/recalled | unclear/recalled | correct/all |")
    lines.append("|---|---|---|---|---|---|---|---|")
    acc = defaultdict(lambda: {"n": 0, "recalled": 0, "accepted": 0, "unclear": 0})
    for r in all_large:
        a = acc[r["kind"]]
        a["n"] += 1
        if r["findings"] >= 1:
            a["recalled"] += 1
            if r["likely"] in accepted_likely(r["kind"], r["size"]):
                a["accepted"] += 1
            elif r["likely"] == "unclear":
                a["unclear"] += 1
    for r in all_small:
        a = acc[r["type"]]
        a["n"] += 1
        if r["findings"] >= 1:
            a["recalled"] += 1
            if r["likely"] in accepted_likely(r["type"], r.get("len", 0), small=True):
                a["accepted"] += 1
            elif r["likely"] == "unclear":
                a["unclear"] += 1
    kind_order = ["del", "ins_random", "dup", "collapse", "inv", "swap",
                  "SNP", "INS", "DEL"]
    for kind in kind_order:
        if kind not in acc:
            continue
        a = acc[kind]
        lines.append(f"| {kind} | {a['n']} | {a['recalled']} | {a['accepted']} | "
                     f"{a['unclear']} | "
                     f"{pct(a['accepted'], a['recalled'])} | "
                     f"{pct(a['unclear'], a['recalled'])} | "
                     f"{pct(a['accepted'], a['n'])} |")

    lines += ["", "## findings: separate findings per planted event", ""]
    lines.append("How many findings one planted event produces (ideally 1), over "
                 "recalled events. `mean` = mean findings per recalled event.")
    lines.append("")
    lines.append("| kind | recalled | 1 | 2 | 3+ | mean |")
    lines.append("|---|---|---|---|---|---|")
    per = defaultdict(lambda: {"n": 0, "one": 0, "two": 0, "more": 0, "sum": 0})
    for r in all_large + all_small:
        if r["findings"] < 1:
            continue
        g = per[event_label(r)]
        g["n"] += 1
        g["sum"] += r["findings"]
        if r["findings"] == 1:
            g["one"] += 1
        elif r["findings"] == 2:
            g["two"] += 1
        else:
            g["more"] += 1
    for kind in kind_order:
        if kind not in per:
            continue
        g = per[kind]
        lines.append(f"| {kind} | {g['n']} | {g['one']} | {g['two']} | "
                     f"{g['more']} | {g['sum'] / g['n']:.2f} |")

    # Controls: findings per Mb by tier + likely ---------------------------
    lines += ["", "## controls: findings per Mb by tier and likely", ""]
    lines.append("All findings in the unmodified-assembly control runs are false "
                 "alarms. `per Mb` uses the summed control assembly length.")
    lines.append("")
    tot_f = sum(c["findings_total"] for c in ctrl_res)
    tot_bp = sum(c["genome_bp"] for c in ctrl_res)
    tier_tot = {"RED": 0, "ORANGE": 0, "YELLOW": 0}
    likely_tot = defaultdict(int)
    for c in ctrl_res:
        for t, n in c.get("findings_tier", {}).items():
            tier_tot[t] = tier_tot.get(t, 0) + n
        for lk, n in c.get("findings_likely", {}).items():
            likely_tot[lk] += n
    fpm = tot_f / (tot_bp / 1e6) if tot_bp else None
    fpm_s = f"{fpm:.2f}" if fpm is not None else "n/a"
    lines.append(f"- controls: {len(ctrl_res)}; findings {tot_f}; per Mb {fpm_s}")
    lines.append(f"- by tier: red {tier_tot['RED']}, orange {tier_tot['ORANGE']}, "
                 f"yellow {tier_tot['YELLOW']}")
    lines.append("- by likely: "
                 + ", ".join(f"{lk} {n}" for lk, n in sorted(likely_tot.items())))
    lines.append("")
    lines.append("| case | findings | per Mb | red | orange | yellow | likely |")
    lines.append("|---|---|---|---|---|---|---|")
    for c in ctrl_res:
        mb = (c["findings_total"] / (c["genome_bp"] / 1e6)) if c["genome_bp"] else None
        mb_s = f"{mb:.2f}" if mb is not None else "n/a"
        ft = c.get("findings_tier", {})
        lk = ", ".join(f"{k} {v}" for k, v in sorted(c.get("findings_likely", {}).items()))
        lines.append(f"| {c['case']} | {c['findings_total']} | {mb_s} | "
                     f"{ft.get('RED', 0)} | {ft.get('ORANGE', 0)} | "
                     f"{ft.get('YELLOW', 0)} | {lk} |")

    lines.append("")
    (out_dir / "Report.md").write_text("\n".join(lines) + "\n")

    # per-event TSV for inspection
    rows = ["case\trep\tclass\tkind_id\tkind\tsize\tcontig\tstart0\tend0\t"
            "breakpoints\tbest_tier\tbest_pass\tsignals\tnear_repeat\tfindings\tlikely"]
    for res in reps_res:
        for r in res["large"]:
            rows.append("\t".join(str(x) for x in [
                r["case"], r["rep"], "large", r["event_id"], r["kind"], r["size"],
                r["contig"], r["start0"], r["end0"], ",".join(map(str, r["breakpoints"])),
                r["best"], r["best_pass"], ";".join(sorted(r["signals"])),
                1 if r["near_repeat"] else 0, r["findings"], r["likely"]]))
        for r in res["small"]:
            rows.append("\t".join(str(x) for x in [
                r["case"], r["rep"], "small", r["type"], r["type"], r.get("len", ""),
                r["contig"], r["pos0"], "", "", r["best"], r["best_pass"],
                ";".join(sorted(r["signals"])), 1 if r["near_repeat"] else 0,
                r["findings"], r["likely"]]))
    (out_dir / "per_event.tsv").write_text("\n".join(rows) + "\n")

    # build summary
    if build_summary:
        lines2 = [f"# build summary (tag `{tag}`)", ""]
        for case, s in build_summary.items():
            if "error" in s:
                lines2.append(f"- {case}: {s['error']}")
                continue
            lines2.append(f"- **{case}** ({s.get('platform')}, {s.get('preset')})")
            cinfo = s.get("control", {})
            lines2.append(f"  - control: {'oracle reuse' if cinfo.get('reused_oracle') else 'mapped'}"
                          f" ({cinfo.get('cram_bytes', 0) / 1e6:.0f} MB)")
            for rep, info in sorted(s.get("reps", {}).items()):
                col = info.get("collapse")
                cdesc = ("none" if not col else
                         f"a={col['a']} b={col['b']} L={col['L']} id={col['ident']}")
                lines2.append(f"  - rep{rep}: contig {info.get('contig')} "
                              f"({info.get('contig_len')} bp), events={info.get('n_events')}, "
                              f"map {info.get('mapping_time_s')}s, "
                              f"CRAM {info.get('cram_bytes', 0) / 1e6:.0f} MB, "
                              f"collapse {cdesc}")
                lines2.append(f"    - per kind: {info.get('events_per_kind')}")
        (out_dir / "BuildSummary.md").write_text("\n".join(lines2) + "\n")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def discover(only):
    for sample in qs.SAMPLES:
        if only and not any(s in sample["case"] for s in only):
            continue
        sdir = SPIKE_ROOT / sample["case"]
        if not sdir.is_dir():
            continue
        reps = sorted(p.name for p in sdir.glob("rep*") if p.is_dir())
        yield sample, sdir, reps


def has_outputs(prefix):
    # `.findings.bed` is required: scoring a run without it would silently count
    # zero findings, so a missing file must invalidate the run instead.
    return all(Path(str(prefix) + e).exists()
               for e in (".vcf", ".bed", ".findings.bed", ".json"))


def any_qc_outputs(root):
    root = Path(root)
    return root.is_dir() and any(root.rglob("qc.json"))


def resolve_qc_root(eval_root, tag, qc_tag, run_flag):
    """Where the qc.* outputs live: the tag dir when running, else reuse."""
    if qc_tag:
        return eval_root / qc_tag
    primary = eval_root / tag
    if run_flag or any_qc_outputs(primary):
        return primary
    sibs = [p for p in eval_root.iterdir()
            if p.is_dir() and p.name != tag and any_qc_outputs(p)]
    if sibs:
        sibs.sort(key=lambda p: max(f.stat().st_mtime for f in p.rglob("qc.json")),
                  reverse=True)
        print(f"[reuse] --tag {tag}: no qc outputs; reusing {sibs[0].name}")
        return sibs[0]
    return primary


def rep_repeats(rep_dir, cache_dir, threads, force):
    """Repeat copies of this rep's spiked.fa, in spiked coordinates (cached)."""
    cache = Path(cache_dir) / "repeats.spiked.tsv"
    if not force and cache.exists():
        return load_repeats(cache)
    pairs = qs.find_repeats(rep_dir / "spiked.fa", threads)
    cache.parent.mkdir(parents=True, exist_ok=True)
    with open(cache, "w") as fh:
        fh.write("contig\tstart0\tend0\tmate_start\tmate_end\tidentity\n")
        for (c, s0, e0, s1, e1, ident) in pairs:
            fh.write(f"{c}\t{s0}\t{e0}\t{s1}\t{e1}\t{ident}\n")
    return load_repeats(cache)


def run(only, replicates, threads, force, tag, run_flag, qc_tag, extra_args=None):
    eval_root = ROOT / "scratch" / "qc" / "spike-eval"
    out_dir = eval_root / tag
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_root = resolve_qc_root(eval_root, tag, qc_tag, run_flag)
    missing = []
    reps_res, ctrl_res = [], []
    for sample, sdir, reps in discover(only):
        case = sample["case"]
        case_dir = ROOT / sample["dir"]
        preset = sample["preset"]
        # control
        cprefix = qc_root / case / "control" / "qc"
        ctrl_ok = (run_qc(sdir / "control" / "mapped.cram", case_dir / "new.fa",
                          preset, cprefix, threads, force, extra_args) if run_flag
                   else has_outputs(cprefix))
        ctrl_calls = []
        if ctrl_ok:
            ctrl_res.append(score_control(case, cprefix))
            ctrl_calls = (qe.parse_vcf(str(cprefix) + ".vcf")
                          + qe.parse_bed(str(cprefix) + ".bed")[0])
            print(f"control {case}: fp={ctrl_res[-1]['fp_total']}")
        else:
            missing.append(f"{case}/control")
            print(f"  missing control outputs: {cprefix}.*", file=sys.stderr)
        for rep in reps:
            rep_num = rep.replace("rep", "")
            if replicates and int(rep_num) not in replicates:
                continue
            rep_dir = sdir / rep
            small, regions, _sum = qe.load_truth(rep_dir)
            large = load_events(rep_dir / "events.tsv")
            prefix = qc_root / case / rep / "qc"
            if run_flag:
                ok = run_qc(rep_dir / "mapped.cram", rep_dir / "spiked.fa", preset,
                            prefix, threads, force, extra_args)
            else:
                ok = has_outputs(prefix)
            if not ok:
                missing.append(f"{case}/{rep}")
                print(f"  missing qc outputs: {prefix}.*", file=sys.stderr)
                continue
            repeats = rep_repeats(rep_dir, out_dir / case / rep, threads, force)
            res = score_rep(case, rep, rep_dir, prefix, small, large, repeats,
                            sdir / "control" / "mapped.cram", case_dir / "new.fa",
                            control_calls=ctrl_calls, case_dir=case_dir,
                            seed=int(rep_num))
            reps_res.append(res)
            n_rec = sum(1 for r in res["large"] if r["best"] >= 1)
            print(f"{case}/{rep}: large recalled(any)={n_rec}/{len(large)} "
                  f"away_fp={res['fp_away_total']}")

    build_summary = {}
    bp = SPIKE_ROOT / "build_summary.json"
    if bp.exists():
        try:
            build_summary = json.loads(bp.read_text())
        except Exception:
            build_summary = {}
    write_report(out_dir, tag, reps_res, ctrl_res, build_summary, qc_root)
    if missing:
        print(f"\nMISSING qc outputs for {len(missing)} run(s): {', '.join(missing)}",
              file=sys.stderr)
    print(f"\nwrote {out_dir}/Report.md (+ per_event.tsv, BuildSummary.md)")
    return 0


# ---------------------------------------------------------------------------
# self-test
# ---------------------------------------------------------------------------
def self_test(_args):
    problems = []
    # large-event matching, both sources and interval kinds
    ev = {"event_id": "del_1000_0", "kind": "del", "size": 1000, "contig": "c",
          "start0": 5000, "end0": 5000, "breakpoints": [5000]}
    vcf_near = {"source": "vcf", "contig": "c", "pos0": 5040, "type": "SNP",
                "tier": "RED", "lowcov": False}
    vcf_far = {"source": "vcf", "contig": "c", "pos0": 5051, "type": "SNP",
               "tier": "RED", "lowcov": False}
    # tol = max(50, 10) = 50 -> 5050 in, 5051 out
    if not large_hit(vcf_near, ev):
        problems.append("vcf within tol not matched")
    if large_hit(vcf_far, ev):
        problems.append("vcf beyond tol matched")
    # size sets the tolerance: 20 kb -> 1% = 200 bp
    big = {"event_id": "del_20000_0", "kind": "del", "size": 20000, "contig": "c",
           "start0": 5000, "end0": 5000, "breakpoints": [5000]}
    if not large_hit({"source": "vcf", "contig": "c", "pos0": 5200, "type": "DEL",
                      "tier": "ORANGE", "lowcov": False}, big):
        problems.append("1% size tolerance not applied")
    # BED_call near a breakpoint
    if not large_hit({"source": "bed", "contig": "c", "start": 4900, "end": 5100,
                      "label": "large_del", "tier": "RED", "lowcov": False}, big):
        problems.append("bed region over breakpoint not matched")
    # interval kind (ins_random): a call inside the interval counts
    ins = {"event_id": "ins_random_1000_0", "kind": "ins_random", "size": 1000,
           "contig": "c", "start0": 7000, "end0": 8000, "breakpoints": [7000, 8000]}
    if not large_hit({"source": "vcf", "contig": "c", "pos0": 7500, "type": "DEL",
                      "tier": "YELLOW", "lowcov": False}, ins):
        problems.append("interval overlap not matched for ins_random")
    # a plain del must NOT match a call strictly inside (no breakpoint nearby)
    if large_hit({"source": "vcf", "contig": "c", "pos0": 7500, "type": "DEL",
                  "tier": "YELLOW", "lowcov": False},
                 {"event_id": "del2", "kind": "del", "size": 1000, "contig": "c",
                  "start0": 5000, "end0": 5000, "breakpoints": [5000, 5100]}):
        problems.append("del matched a call far from every breakpoint")

    # tier aggregation + PASS split
    calls = [{"source": "vcf", "contig": "c", "pos0": 5000, "type": "DEL",
              "tier": "ORANGE", "lowcov": True},
             {"source": "bed", "contig": "c", "start": 4990, "end": 5060,
              "label": "large_del", "tier": "RED", "lowcov": False}]
    best, best_pass, sig = best_tiers(calls, lambda c: large_hit(c, ev))
    if best != 3 or best_pass != 3:
        problems.append(f"best tiers {best}/{best_pass} != 3/3")
    if sig != {"DEL", "large_del"}:
        problems.append(f"signals {sig} != {{DEL, large_del}}")

    # false-positive window: 2 kb
    small = [{"contig": "c", "pos0": 100000, "type": "SNP", "len": 1, "ref": "A",
              "alt": "C", "hp": 1, "str": 0, "start0": 100000, "end0": 100001}]
    if not near_any_event({"source": "vcf", "contig": "c", "pos0": 101999, "tier": "RED"},
                          small, []):
        problems.append("call at +1999 not counted near")
    if not near_any_event({"source": "vcf", "contig": "c", "pos0": 102000, "tier": "RED"},
                          small, []):
        problems.append("call at +2000 not counted near")
    if near_any_event({"source": "vcf", "contig": "c", "pos0": 102002, "tier": "RED"},
                      small, []):
        problems.append("call at +2002 (past window) counted near")

    # cross-mapping footprint logic (synthetic)
    if footprint_intervals([(100, 200), (150, 250)]) != []:
        problems.append("footprint with 2 reads should be empty")
    if footprint_intervals([(100, 200), (150, 250), (120, 180)]) != [(150, 180)]:
        problems.append("footprint with 3 overlapping reads")
    fp2 = footprint_intervals([(0, 100), (10, 90), (20, 80),
                               (500, 600), (510, 590), (520, 580)])
    if fp2 != [(20, 80), (520, 580)]:
        problems.append(f"footprint clusters: {fp2}")
    if abs(frac_inside(120, 220, [(100, 200)]) - 0.8) > 1e-9:
        problems.append("frac_inside 0.8")
    if abs(frac_inside(150, 250, [(100, 200)]) - 0.5) > 1e-9:
        problems.append("frac_inside 0.5")
    if frac_inside(200, 300, [(100, 200)]) != 0.0:
        problems.append("frac_inside disjoint != 0")
    if cigar_ref_span("10S20M5I5D30M") != 55:
        problems.append(f"cigar_ref_span {cigar_ref_span('10S20M5I5D30M')}")
    fp_calls = [{"source": "vcf", "contig": "c", "pos0": 650}]
    if not in_footprint(fp_calls[0], [("c", 500, 600)], margin=500):
        problems.append("in_footprint +50 (within margin) not matched")
    if in_footprint({"source": "vcf", "contig": "c", "pos0": 1101},
                    [("c", 500, 600)], margin=500):
        problems.append("in_footprint far call matched")

    # large-gap CIGAR: D/N advance the reference without covering the read.
    if cigar_blocks("100M20000D100M", 0) != [(0, 100), (20100, 20200)]:
        problems.append(f"cigar_blocks large D: {cigar_blocks('100M20000D100M', 0)}")
    if cigar_aligned_len("100M20000D100M") != 200:
        problems.append(f"cigar_aligned_len large D: {cigar_aligned_len('100M20000D100M')}")
    if cigar_aligned_len("10S20M5I5D30M") != 50:
        problems.append(f"cigar_aligned_len clips/ins: {cigar_aligned_len('10S20M5I5D30M')}")
    # a read whose only overlap with the deleted interval is the skipped gap is
    # not "inside" it (the old ref-span rule wrongly selected it).
    if covered_frac(cigar_blocks("100M20000D100M", 0), [(100, 20100)]) != 0.0:
        problems.append("covered_frac counted the skipped gap")
    if abs(covered_frac(cigar_blocks("100M20000D100M", 0), [(0, 100)]) - 0.5) > 1e-9:
        problems.append("covered_frac exact half")

    # pre-existing control call inside a footprint: the projection finds it, so
    # the call is not cross-mapping. Offsets are keyed by edited contig; a
    # contig without edits (a plasmid) projects unchanged.
    offsets = {"chr": [(5000, 17051, 12051)]}  # one 12051 bp insertion
    ccalls = [{"source": "vcf", "contig": "chr", "pos0": 631673, "tier": "YELLOW"},
              {"source": "bed", "contig": "chr", "start": 787900, "end": 788200,
               "label": "clip_hotspot", "tier": "YELLOW"},
              {"source": "vcf", "contig": "plasmid", "pos0": 1000, "tier": "YELLOW"}]
    pre = {"source": "vcf", "contig": "chr", "pos0": 643724}
    if project_to_original(643724, offsets, "chr") != 631673:
        problems.append(f"project_to_original: {project_to_original(643724, offsets, 'chr')}")
    if project_to_original(643724, offsets, "plasmid") != 643724:
        problems.append("unedited contig was shifted")
    if not has_control_match(pre, ccalls, offsets):
        problems.append("pre-existing control call inside a footprint not detected")
    if not has_control_match({"source": "vcf", "contig": "plasmid", "pos0": 1000},
                             ccalls, offsets):
        problems.append("identical plasmid control call not matched")
    if has_control_match({"source": "vcf", "contig": "chr", "pos0": 700000}, ccalls, offsets):
        problems.append("absent call wrongly matched to the control")
    if not has_control_match({"source": "bed", "contig": "chr", "start": 800050,
                              "end": 800150}, ccalls, offsets):
        problems.append("overlapping control region not matched")
    if core_interval(1000, 1100) != (1025, 1075):
        problems.append(f"core_interval small: {core_interval(1000, 1100)}")
    if core_interval(0, 10000) != (500, 9500):
        problems.append(f"core_interval large: {core_interval(0, 10000)}")

    # indel support matches the called allele (type + length + inserted sequence)
    if not allele_alt_match(False, 1, "C", "A", "C", "CA"):
        problems.append("called INS not matched")
    if allele_alt_match(False, 2, "C", "TT", "C", "CA"):
        problems.append("different INS length wrongly matched")
    if allele_alt_match(False, 1, "C", "G", "C", "CA"):
        problems.append("different INS sequence wrongly matched")
    if not allele_alt_match(False, -2, None, None, "ACT", "A"):
        problems.append("called DEL not matched")
    if allele_alt_match(False, -1, None, None, "ACT", "A"):
        problems.append("different DEL length wrongly matched")
    if not allele_alt_match(False, 0, "G", None, "A", "G"):
        problems.append("called SNP not matched")
    if allele_alt_match(False, 0, "A", None, "A", "G"):
        problems.append("REF base wrongly matched as alt")
    # insertion-grouping tolerance for INS >= 10 bp
    if not _ins_matches("A" * 9 + "C", "A" * 10):
        problems.append("1-edit 10 bp INS not grouped")
    if _ins_matches("C" * 10, "A" * 10):
        problems.append("unrelated 10 bp INS wrongly grouped")
    if _ins_matches("A", "C"):
        problems.append("short INS grouped (must match exactly)")
    # C->CA support: two footprint reads insert TT, one non-footprint inserts A
    alt_reads = [allele_alt_match(False, 2, "C", "TT", "C", "CA"),
                 allele_alt_match(False, 2, "C", "TT", "C", "CA"),
                 allele_alt_match(False, 1, "C", "A", "C", "CA")]
    fp_reads = {0, 1}
    alt_n = sum(alt_reads)
    fp_n = sum(1 for i, v in enumerate(alt_reads) if v and i in fp_reads)
    if (fp_n, alt_n) != (0, 1):
        problems.append(f"C->CA support {fp_n}/{alt_n} != 0/1")

    # region support is signal-specific
    class _Read:
        def __init__(self, name, start, cig):
            self.query_name, self.reference_start, self.cigartuples = name, start, cig
            self.is_unmapped = self.is_secondary = self.is_supplementary = False

    class _Aln:
        def __init__(self, name, seq):
            self.query_name, self.query_sequence = name, seq

    class _PR:
        def __init__(self, name, seq, qp, indel=0, is_del=False):
            self.alignment, self.query_position = _Aln(name, seq), qp
            self.indel, self.is_del, self.is_refskip = indel, is_del, False

    class _Col:
        def __init__(self, prs):
            self.pileups = prs

    class _Cram:
        def __init__(self, reads=(), pileups=None):
            self.reads, self.pileups = list(reads), pileups or {}

        def fetch(self, contig, cs, ce):
            for r in self.reads:
                end = r.reference_start + sum(l for op, l in (r.cigartuples or [])
                                              if op in (0, 2, 3, 7, 8))
                if r.reference_start < ce and cs < end:
                    yield r

        def pileup(self, contig, pos0, end, **k):
            return [_Col(self.pileups.get(pos0, []))]

    clip = {"source": "bed", "contig": "c", "label": "clip_hotspot",
            "start": 100, "end": 200}
    cram = _Cram([_Read("clip", 100, [(4, 60), (0, 50)]), _Read("noclip", 120, [(0, 50)])])
    if region_supporting_reads(cram, clip, {}) != {"clip"}:
        problems.append(f"clip_hotspot supporting: {region_supporting_reads(cram, clip, {})}")
    if region_footprint_support(cram, clip, {"clip"}, {}) != (1, 1):
        problems.append("clip_hotspot footprint fraction wrong")

    ldel = {"source": "bed", "contig": "c", "label": "large_del", "start": 100, "end": 250}
    cram = _Cram([_Read("bigdel", 100, [(0, 50), (2, 100), (0, 50)]),
                  _Read("smalldel", 100, [(0, 50), (2, 20), (0, 50)])])
    if region_supporting_reads(cram, ldel, {}) != {"bigdel"}:
        problems.append(f"large_del supporting: {region_supporting_reads(cram, ldel, {})}")

    depth = {"source": "bed", "contig": "c", "label": "depth_low", "start": 100, "end": 200}
    cram = _Cram([_Read("a", 100, [(0, 50)]), _Read("b", 120, [(0, 50)])])
    if region_supporting_reads(cram, depth, {}) != {"a", "b"}:
        problems.append("depth region must use the overlap rule")

    cluster = {"source": "bed", "contig": "c", "label": "snp_cluster",
               "start": 50, "end": 200}
    cram = _Cram(pileups={100: [_PR("alt", "G", 0), _PR("ref", "A", 0)]})
    if _cluster_alt_reads(cram, cluster, {("c", 100): ("A", "G")}) != {"alt"}:
        problems.append("snp_cluster supporting reads wrong")

    # affected intervals: a call inside an inversion / duplication interval (far
    # from its breakpoints) is not an FP; a call inside a swap block is not an FP
    # but a call in the span *between* the two swap blocks is.
    inv = {"event_id": "inv_20000_0", "kind": "inv", "size": 20000, "contig": "c",
           "start0": 300000, "end0": 320000, "breakpoints": [300000, 320000]}
    if affected_intervals(inv) != [(300000, 320000)]:
        problems.append(f"inv affected_intervals {affected_intervals(inv)}")
    if not near_any_event({"source": "vcf", "contig": "c", "pos0": 310000, "tier": "RED"},
                          [], [inv]):
        problems.append("call inside an inversion interval counted as away")
    dup = {"event_id": "dup_10000_0", "kind": "dup", "size": 10000, "contig": "c",
           "start0": 400000, "end0": 420000, "breakpoints": [410000, 420000]}
    if not near_any_event({"source": "vcf", "contig": "c", "pos0": 405000, "tier": "RED"},
                          [], [dup]):
        problems.append("call inside a duplication interval counted as away")
    sw = {"event_id": "swap_10000_0", "kind": "swap", "size": 10000, "contig": "c",
          "start0": 500000, "end0": 750000, "breakpoints": [500000, 510000, 740000, 750000]}
    if affected_intervals(sw) != [(500000, 510000), (740000, 750000)]:
        problems.append(f"swap affected_intervals {affected_intervals(sw)}")
    if not near_any_event({"source": "vcf", "contig": "c", "pos0": 745000, "tier": "RED"},
                          [], [sw]):
        problems.append("call inside a swap block counted as away")
    if near_any_event({"source": "vcf", "contig": "c", "pos0": 620000, "tier": "RED"},
                      [], [sw]):
        problems.append("call in the span between swap blocks counted as near")

    # repeat-miss list must accept small (type/len) and large (kind/size) rows
    mixed = [
        {"case": "c", "rep": "rep1", "kind": "inv", "size": 5000, "contig": "x",
         "start0": 100, "end0": 5100, "best": 0, "near_repeat": True, "pos0": 100},
        {"case": "c", "rep": "rep1", "type": "SNP", "len": 1, "contig": "x",
         "pos0": 200, "start0": 200, "end0": 201, "best": 0, "near_repeat": True},
        {"case": "c", "rep": "rep1", "kind": "del", "size": 50, "contig": "x",
         "start0": 300, "end0": 300, "best": 3, "near_repeat": True, "pos0": 300},
    ]
    rm = repeat_misses(mixed)
    if len(rm) != 2:
        problems.append(f"repeat_misses picked {len(rm)} != 2")
    if {event_label(r) for r in rm} != {"inv", "SNP"}:
        problems.append(f"repeat_misses labels {[event_label(r) for r in rm]}")
    if event_size(rm[0]) not in (5000, 1):
        problems.append("repeat_misses size fallback broken")

    # findings: parsing, best selection, accepted answers, matching.  Fixtures
    # live in a fresh temp dir so consecutive runs cannot see stale files.
    fdir = Path(tempfile.mkdtemp())
    try:
        fpath = fdir / "qc.findings.bed"
        fpath.write_text(
            "#contig\tstart\tend\tlabel\tscore\ttier\tlikely\textra\n"
            "c\t1000\t1100\tfinding_1\t2\tred\tmisjoin\t"
            "evidence=clips,clips_both;members=DEL:1,clip_hotspot:1;lowcov=0\n"
            "c\t5000\t5100\tfinding_2\t1\tyellow\tminority\t"
            "evidence=site_minor;members=SNP:1;lowcov=1\n")
        fs = parse_findings(fpath)
        if len(fs) != 2 or fs[0]["likely"] != "misjoin" or fs[1]["tier"] != "YELLOW":
            problems.append(f"parse_findings: {fs}")
        if not fs[1]["lowcov"] or fs[0]["members"] != "DEL:1,clip_hotspot:1":
            problems.append(f"parse_findings extra: {fs}")
        if fs[0]["score"] != 2 or fs[0]["start"] != 1000 or fs[0]["end"] != 1100:
            problems.append(f"parse_findings columns: {fs[0]}")
        bf = best_finding(fs)
        if bf is None or bf["likely"] != "misjoin":
            problems.append(f"best_finding: {bf}")
        if best_finding([]) is not None:
            problems.append("best_finding([]) should be None")
        if accepted_likely("del", 1000) != {"missing_sequence", "misjoin"}:
            problems.append("accepted_likely del <=1kb")
        if accepted_likely("del", 20000) != {"missing_sequence", "misjoin"}:
            problems.append("accepted_likely del >=20kb")
        if accepted_likely("collapse", 1000) != {"collapsed_repeat", "missing_sequence"}:
            problems.append("accepted_likely collapse")
        if accepted_likely("dup", 5000) != {"extra_sequence", "false_duplication"}:
            problems.append("accepted_likely dup")
        if accepted_likely("SNP", 1, small=True) != {"base_error"}:
            problems.append("accepted_likely small")
        if not large_hit(finding_call(fs[1]), ev):
            problems.append("finding_call/large_hit not matched")

        # has_outputs must require qc.findings.bed (a missing file invalidates the run)
        hp = fdir / "hasout" / "qc"
        hp.parent.mkdir(parents=True, exist_ok=True)
        for e in (".vcf", ".bed", ".json"):
            Path(str(hp) + e).write_text("")
        if has_outputs(hp):
            problems.append("has_outputs accepted a run without qc.findings.bed")
        Path(str(hp) + ".findings.bed").write_text("")
        if not has_outputs(hp):
            problems.append("has_outputs rejected a complete run")

        # a report with zero away calls must not divide by zero
        zrow = {"case": "z", "rep": "rep1", "genome_bp": 1_000_000,
                "large": [], "small": [],
                "fp": {"RED": 0, "ORANGE": 0, "YELLOW": 0},
                "fp_lowcov": {"RED": 0, "ORANGE": 0, "YELLOW": 0},
                "fp_away_total": 0, "cross": {"RED": 0, "ORANGE": 0, "YELLOW": 0},
                "cross_total": 0, "fp_net_total": 0, "footprints": [],
                "calls": {"vcf": 0, "bed": 0}}
        try:
            (fdir / "zero").mkdir(parents=True, exist_ok=True)
            write_report(fdir / "zero", "t", [zrow], [], {})
            if "away FP 0" not in (fdir / "zero" / "Report.md").read_text():
                problems.append("zero-away report missing totals")
        except ZeroDivisionError:
            problems.append("write_report divided by zero away calls")
    finally:
        shutil.rmtree(fdir, ignore_errors=True)

    if problems:
        print("SELF-TEST FAILED:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 1
    print("self-test OK: large matching, tier/PASS split, signals, FP window, "
          "inv/dup/swap intervals, repeat-miss fallbacks")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="store_true",
                    help="run cramhook qc; without it, score existing outputs")
    ap.add_argument("--tag", default="sp1")
    ap.add_argument("--qc-tag", default=None,
                    help="tag dir holding the qc.* outputs (default: --tag; "
                         "without --run, falls back to an existing tag dir)")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--samples", default=None)
    ap.add_argument("--replicates", default=None,
                    help="comma-separated replicate numbers (default: all found)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--extra-args", default=None,
                    help="extra arguments appended to every `cramhook qc` call "
                         "(e.g. '--unlinked-action drop')")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test(args)
    only = [s for s in args.samples.split(",") if s] if args.samples else None
    reps = [int(x) for x in args.replicates.split(",") if x] if args.replicates else []
    return run(only, reps, args.threads, args.force, args.tag, args.run, args.qc_tag, args.extra_args)


if __name__ == "__main__":
    sys.exit(main())
