#!/usr/bin/env python3
"""Compare a lifted CRAM against the oracle mapping (docs/project/plans/Plan-v1.md §6.3).

Usage:
  python scripts/compare.py --test A.cram --test-ref A_REF.fa \
      --oracle O.cram --ref NEW.fa --out PREFIX

Writes PREFIX.json and PREFIX.txt. Exposes compare(test, test_ref, oracle, ref).
"""
import argparse
import json
import os
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
import pysam

FIELDS = ["contig", "pos", "strand", "cigar", "mapq"]
CATS = ["T", "R", "none"]
MAX_MISMATCHES = 5000

# M7b (user decision): supplementary records with MAPQ below this are ignored on
# both sides before comparing supplementary sets.  Such records carry no
# placement information, and minimap2 adds or drops them depending on the
# orientation of the read, so their presence/absence is not a real disagreement.
SUPP_MIN_MAPQ = 20

# pysam.cigartuples operation codes.
REF_CONSUME_CODES = frozenset({0, 2, 3, 7, 8})    # M D N = X
QUERY_CONSUME_CODES = frozenset({0, 1, 4, 7, 8})  # M I S = X


def _fai_path(ref):
    return Path(str(ref) + ".fai")


def ensure_fai(ref):
    if not _fai_path(ref).exists():
        subprocess.run(["samtools", "faidx", str(ref)], check=True)
    return ref


def read_fai_lengths(ref):
    ensure_fai(ref)
    lengths = {}
    with open(str(ref) + ".fai") as fh:
        for line in fh:
            parts = line.split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


def rec_info(r):
    """(contig, pos, strand, cigar, mapq); unmapped -> (None, -1, False, None, 0)."""
    if r.is_unmapped:
        return (None, -1, False, None, 0)
    return (r.reference_name, r.reference_start, r.is_reverse, r.cigarstring, r.mapping_quality)


def xl_cat(r):
    try:
        if r.has_tag("XL"):
            return str(r.get_tag("XL"))
    except Exception:
        pass
    return "none"


def _tag(r, name):
    """Value of a record tag, or None if absent/unparseable."""
    try:
        return r.get_tag(name)
    except Exception:
        return None


def _record_key(detail):
    """Hashable equivalence key for a record.

    ``(contig, pos, END, strand, mapq, read_length, AS)``.  ``END`` is the
    reference end (pos + reference span); ``read_length`` includes soft clips.
    Returns ``None`` for unmapped records or when the CIGAR or AS tag is
    missing (a missing AS means the record is not equivalent to anything).
    ``NM`` and ``de`` are deliberately not part of the key; they are reported
    as metrics for equivalent pairs instead.
    """
    if not detail:
        return None
    contig, pos, strand, _cigar, mapq = detail["info"]
    tuples = detail.get("tuples")
    as_tag = detail.get("AS")
    if contig is None or not tuples or as_tag is None:
        return None
    span = sum(l for c, l in tuples if c in REF_CONSUME_CODES)
    qlen = sum(l for c, l in tuples if c in QUERY_CONSUME_CODES)
    return (contig, pos, pos + span, strand, mapq, qlen, as_tag)


def _rec_detail(r):
    if r.is_unmapped:
        return {"info": rec_info(r), "tuples": None, "AS": None, "NM": None, "de": None}
    return {"info": rec_info(r), "tuples": list(r.cigartuples or []),
            "AS": _tag(r, "AS"), "NM": _tag(r, "NM"), "de": _tag(r, "de")}


def _equiv_metrics(td, od):
    """(nm_equal, de_equal, test_NM, oracle_NM, test_de, oracle_de)."""
    tn = td.get("NM") if td else None
    on = od.get("NM") if od else None
    tde = td.get("de") if td else None
    ode = od.get("de") if od else None
    return (tn == on), (tde == ode), tn, on, tde, ode


def _open(path, ref):
    return pysam.AlignmentFile(str(path), "rc", reference_filename=str(ref))


def collect(path, ref):
    """primary dict name->info, primary name Counter, supp dict name->set, n_secondary,
    xl dict, detail dict name->{primary, supp list}."""
    primary = {}
    names = Counter()
    supp = {}
    xl = {}
    detail = {}
    n_secondary = 0
    af = _open(path, ref)
    try:
        for r in af.fetch(until_eof=True):
            if r.flag & 0x100:
                n_secondary += 1
                continue
            name = r.query_name
            d = detail.setdefault(name, {"primary": None, "supp": []})
            if r.flag & 0x800:
                supp.setdefault(name, set()).add(rec_info(r))
                d["supp"].append(_rec_detail(r))
                continue
            primary[name] = rec_info(r)
            names[name] += 1
            xl[name] = xl_cat(r)
            d["primary"] = _rec_detail(r)
    finally:
        af.close()
    return primary, names, supp, n_secondary, xl, detail


def depth_arrays(path, ref, lengths, min_mapq=20):
    """Per-contig depth (numpy) from primary mapped records with MAPQ >= min_mapq, M/=/X/D counted (samtools depth -J)."""
    diff = {c: np.zeros(n + 1, dtype=np.int64) for c, n in lengths.items()}
    af = _open(path, ref)
    try:
        for r in af.fetch(until_eof=True):
            if r.flag & 0x100 or r.flag & 0x800:
                continue
            if r.is_unmapped or r.mapping_quality < min_mapq:
                continue
            arr = diff.get(r.reference_name)
            if arr is None:
                continue
            n = len(arr) - 1  # test may be on a same-named contig of different length
            # count M/=/X and D (like `samtools depth -J`); N and clips/insertions don't cover ref
            p = r.reference_start
            for op, ln in r.cigartuples:
                if op in (0, 2, 7, 8):
                    arr[min(p, n)] += 1
                    arr[min(p + ln, n)] -= 1
                if op in (0, 2, 3, 7, 8):
                    p += ln
    finally:
        af.close()
    return {c: np.cumsum(a)[:-1] for c, a in diff.items()}


def _empty_stats():
    d = {"n": 0}
    for f in FIELDS:
        d[f] = 0
    d["all"] = 0
    d["n_identical"] = 0
    d["n_equivalent"] = 0
    d["n_mismatch"] = 0
    d["equiv_nm_equal"] = 0
    d["equiv_nm_differ"] = 0
    d["equiv_de_equal"] = 0
    d["equiv_de_differ"] = 0
    return d


def _empty_supp():
    return {"n": 0, "equal": 0, "equivalent": 0, "differ": 0,
            "differ_raw": 0, "low_mapq_ignored_test": 0,
            "low_mapq_ignored_oracle": 0}


def _supp_status(a_set, b_set, a_list, b_list, min_mapq=SUPP_MIN_MAPQ):
    """equal / equivalent / differ for one read's supplementary record sets.

    Supp records with MAPQ < ``min_mapq`` are dropped from both sides first
    (M7b); pass ``min_mapq=0`` for the raw, unfiltered comparison.
    """
    a_set = {i for i in a_set if i[4] >= min_mapq}
    b_set = {i for i in b_set if i[4] >= min_mapq}
    a_list = [d for d in a_list if d["info"][4] >= min_mapq]
    b_list = [d for d in b_list if d["info"][4] >= min_mapq]
    if a_set == b_set:
        return "equal"
    ka = [_record_key(d) for d in a_list]
    kb = [_record_key(d) for d in b_list]
    if any(k is None for k in ka) or any(k is None for k in kb):
        return "differ"  # a missing AS tag is never equivalent
    return "equivalent" if Counter(ka) == Counter(kb) else "differ"


def compare(test, test_ref, oracle, ref):
    lengths = read_fai_lengths(ref)
    t_primary, t_names, t_supp, t_sec, t_xl, t_detail = collect(test, test_ref)
    o_primary, o_names, o_supp, o_sec, _, o_detail = collect(oracle, ref)

    def inv(names, n_sec):
        return {
            "n_primary": int(sum(names.values())),
            "n_unique_names": len(names),
            "all_unique": all(v == 1 for v in names.values()),
            "n_secondary": n_sec,
        }

    orphan_supp_test = len(set(t_supp) - set(t_names))
    orphan_supp_oracle = len(set(o_supp) - set(o_names))
    invariants = {
        "test": inv(t_names, t_sec),
        "oracle": inv(o_names, o_sec),
        "names_only_in_test": len(set(t_names) - set(o_names)),
        "names_only_in_oracle": len(set(o_names) - set(t_names)),
        "same_names": set(t_names) == set(o_names),
        "orphan_supplementary_test": orphan_supp_test,
        "orphan_supplementary_oracle": orphan_supp_oracle,
    }
    invariants_ok = bool(
        invariants["test"]["all_unique"]
        and invariants["oracle"]["all_unique"]
        and invariants["same_names"]
        and t_sec == 0
        and o_sec == 0
        and orphan_supp_test == 0
        and orphan_supp_oracle == 0
    )

    per = {c: _empty_stats() for c in CATS}
    overall = _empty_stats()
    mismatches = []
    equivalents = []
    total_mismatch = 0

    # union of primary read names: a name in only one file is a mismatch.
    all_names = list(t_primary)
    all_names += [n for n in o_primary if n not in t_primary]

    for name in all_names:
        info = t_primary.get(name)
        o = o_primary.get(name)
        if info is not None and o is not None:
            cat = t_xl.get(name, "none")
            if cat not in CATS:
                cat = "none"
            eq = {
                "contig": info[0] == o[0],
                "pos": info[1] == o[1],
                "strand": info[2] == o[2],
                "cigar": info[3] == o[3],
                "mapq": info[4] == o[4],
            }
            category = ",".join(f for f in FIELDS if not eq[f]) or "equal"
        elif info is not None:
            cat = t_xl.get(name, "none")
            if cat not in CATS:
                cat = "none"
            eq = {f: False for f in FIELDS}
            category = "missing_in_oracle"
        else:
            cat = "missing_in_test"
            eq = {f: False for f in FIELDS}
            category = "missing_in_test"
        per.setdefault(cat, _empty_stats())
        all_eq = all(eq.values())
        # equivalent: contig/pos/END/strand/mapq, read length and AS all agree
        # (CIGAR differs); a missing AS is never equivalent.  NM and de are not
        # part of the definition, only reported below.
        is_equiv = False
        nm_eq = de_eq = None
        pair_extra = {}
        if info is not None and o is not None and not all_eq:
            core_eq = eq["contig"] and eq["pos"] and eq["strand"] and eq["mapq"]
            if core_eq and not eq["cigar"]:
                td = (t_detail.get(name) or {}).get("primary")
                od = (o_detail.get(name) or {}).get("primary")
                kt = _record_key(td)
                is_equiv = kt is not None and kt == _record_key(od)
                if is_equiv:
                    nm_eq, de_eq, tnm, onm, tde, ode = _equiv_metrics(td, od)
                    pair_extra = {"test_NM": tnm, "oracle_NM": onm,
                                  "test_de": tde, "oracle_de": ode}
        for d in (per[cat], overall):
            d["n"] += 1
            for f in FIELDS:
                if eq[f]:
                    d[f] += 1
            if all_eq:
                d["all"] += 1
                d["n_identical"] += 1
            elif is_equiv:
                d["n_equivalent"] += 1
                d["equiv_nm_equal" if nm_eq else "equiv_nm_differ"] += 1
                d["equiv_de_equal" if de_eq else "equiv_de_differ"] += 1
            else:
                d["n_mismatch"] += 1
        pair = {
            "name": name,
            "category": category,
            "test": list(info) if info is not None else None,
            "oracle": list(o) if o is not None else None,
        }
        if all_eq:
            continue
        if is_equiv:
            pair.update(pair_extra)
            equivalents.append(pair)
        else:
            total_mismatch += 1
            mismatches.append(pair)

    supp_per = {c: _empty_supp() for c in CATS}
    supp_overall = _empty_supp()
    for name in all_names:
        cat = t_xl.get(name, "none")
        if cat not in supp_per:
            cat = "none"
        a = t_supp.get(name, set())
        b = o_supp.get(name, set())
        al = (t_detail.get(name) or {}).get("supp", [])
        bl = (o_detail.get(name) or {}).get("supp", [])
        status = _supp_status(a, b, al, bl)
        raw = _supp_status(a, b, al, bl, min_mapq=0)
        ignored_t = sum(1 for d in al if d["info"][4] < SUPP_MIN_MAPQ)
        ignored_o = sum(1 for d in bl if d["info"][4] < SUPP_MIN_MAPQ)
        for d in (supp_per[cat], supp_overall):
            d["n"] += 1
            d[status] += 1
            d["low_mapq_ignored_test"] += ignored_t
            d["low_mapq_ignored_oracle"] += ignored_o
            if raw == "differ":
                d["differ_raw"] += 1

    t_depth = depth_arrays(test, test_ref, lengths)
    o_depth = depth_arrays(oracle, ref, lengths)
    contigs = {}
    tot_n = tot_id = tot_gt2 = maxd = 0
    for c, n in lengths.items():
        d = np.abs(t_depth[c] - o_depth[c])
        ident = int((d == 0).sum())
        mx = int(d.max()) if n else 0
        gt2 = int((d > 2).sum())
        contigs[c] = {
            "n": n,
            "n_identical": ident,
            "fraction_identical": (ident / n) if n else 1.0,
            "max_abs_diff": mx,
            "n_gt2": gt2,
        }
        tot_n += n
        tot_id += ident
        tot_gt2 += gt2
        maxd = max(maxd, mx)
    depth = {
        "overall": {
            "n_positions": tot_n,
            "n_identical": tot_id,
            "fraction_identical": (tot_id / tot_n) if tot_n else 1.0,
            "max_abs_diff": maxd,
            "n_gt2": tot_gt2,
        },
        "contigs": contigs,
    }

    return {
        "test": str(test),
        "test_ref": str(test_ref),
        "oracle": str(oracle),
        "ref": str(ref),
        "invariants": invariants,
        "invariants_ok": invariants_ok,
        "primary": {**per, "overall": overall},
        "mismatch_total": total_mismatch,
        "mismatches": mismatches[:MAX_MISMATCHES],
        "mismatches_full": mismatches,
        "mismatches_truncated_in_json": len(mismatches) > MAX_MISMATCHES,
        "equivalent_total": len(equivalents),
        "equivalents": equivalents[:MAX_MISMATCHES],
        "equivalents_full": equivalents,
        "equivalents_truncated_in_json": len(equivalents) > MAX_MISMATCHES,
        "supplementary": {**supp_per, "overall": supp_overall},
        "depth": depth,
    }


def _pct(x, n):
    return 100.0 * x / n if n else 100.0


def summary_text(res):
    inv = res["invariants"]
    pr = res["primary"]
    sp = res["supplementary"]
    dp = res["depth"]["overall"]
    lines = [
        f"test={res['test']} oracle={res['oracle']} ref={res['ref']}",
        "reads: test primary={} secondary={} unique={} | oracle primary={} secondary={} unique={}".format(
            inv["test"]["n_primary"], inv["test"]["n_secondary"], inv["test"]["all_unique"],
            inv["oracle"]["n_primary"], inv["oracle"]["n_secondary"], inv["oracle"]["all_unique"]),
        "names: same={} only_in_test={} only_in_oracle={}".format(
            inv["same_names"], inv["names_only_in_test"], inv["names_only_in_oracle"]),
        "invariants_ok={} orphan_supplementary: test={} oracle={}".format(
            res["invariants_ok"], inv["orphan_supplementary_test"],
            inv["orphan_supplementary_oracle"]),
    ]
    for c in CATS:
        d = pr.get(c)
        if not d or d["n"] == 0:
            continue
        lines.append(
            "primary[{}]: n={} identical={:.4f}% equivalent={:.4f}% "
            "identical+equivalent={:.4f}% mismatch={} "
            "(equiv NM eq/diff {}/{}, de eq/diff {}/{}) "
            "(contig {:.3f} pos {:.3f} strand {:.3f} cigar {:.3f} mapq {:.3f})".format(
                c, d["n"], _pct(d["n_identical"], d["n"]),
                _pct(d["n_equivalent"], d["n"]),
                _pct(d["n_identical"] + d["n_equivalent"], d["n"]), d["n_mismatch"],
                d["equiv_nm_equal"], d["equiv_nm_differ"],
                d["equiv_de_equal"], d["equiv_de_differ"],
                _pct(d["contig"], d["n"]), _pct(d["pos"], d["n"]), _pct(d["strand"], d["n"]),
                _pct(d["cigar"], d["n"]), _pct(d["mapq"], d["n"])))
    o = pr.get("overall")
    if o and o["n"]:
        lines.append(
            "primary[overall]: n={} identical={:.4f}% equivalent={:.4f}% "
            "identical+equivalent={:.4f}% mismatch={} "
            "(equiv NM eq/diff {}/{}, de eq/diff {}/{})".format(
                o["n"], _pct(o["n_identical"], o["n"]), _pct(o["n_equivalent"], o["n"]),
                _pct(o["n_identical"] + o["n_equivalent"], o["n"]), o["n_mismatch"],
                o["equiv_nm_equal"], o["equiv_nm_differ"],
                o["equiv_de_equal"], o["equiv_de_differ"]))
    so = sp["overall"]
    if so["n"]:
        lines.append(
            "supplementary[overall]: n={} equal={} equivalent={} differ={} "
            "differ_raw={} low_mapq_ignored(test={} oracle={}) "
            "({:.4f}% equal, {:.4f}% equal+equivalent)".format(
                so["n"], so["equal"], so["equivalent"], so["differ"],
                so["differ_raw"], so["low_mapq_ignored_test"],
                so["low_mapq_ignored_oracle"],
                _pct(so["equal"], so["n"]),
                _pct(so["equal"] + so["equivalent"], so["n"])))
    lines.append("depth: positions={} identical={:.4f}% max_abs_diff={} n_gt2={}".format(
        dp["n_positions"], 100.0 * dp["fraction_identical"], dp["max_abs_diff"], dp["n_gt2"]))
    return "\n".join(lines)


def _atomic_write(path, text):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


TSV_COLS = ["name", "category",
            "test_contig", "test_pos", "test_strand", "test_cigar", "test_mapq",
            "oracle_contig", "oracle_pos", "oracle_strand", "oracle_cigar", "oracle_mapq"]
TSV_EQ_EXTRA = ["test_NM", "test_de", "oracle_NM", "oracle_de"]


def _pair_fields(m):
    fields = [str(m.get("name", "")), str(m.get("category", ""))]
    for rec in (m.get("test"), m.get("oracle")):
        if rec is None:
            fields += ["", "", "", "", ""]
        else:
            fields += ["" if v is None else str(v) for v in rec]
    return fields


def _pairs_tsv(rows, path):
    out = ["\t".join(TSV_COLS)]
    for m in rows:
        out.append("\t".join(_pair_fields(m)))
    _atomic_write(path, "\n".join(out) + "\n")


def _equiv_tsv(rows, path):
    out = ["\t".join(TSV_COLS + TSV_EQ_EXTRA)]
    for m in rows:
        fields = _pair_fields(m)
        for k in TSV_EQ_EXTRA:
            v = m.get(k)
            fields.append("" if v is None else str(v))
        out.append("\t".join(fields))
    _atomic_write(path, "\n".join(out) + "\n")


def write_result(res, prefix):
    prefix = Path(prefix)
    full = res.get("mismatches_full", res.get("mismatches", []))
    eq_full = res.get("equivalents_full", [])
    out = {k: v for k, v in res.items()
           if k not in ("mismatches_full", "equivalents_full")}
    out["mismatches"] = full[:MAX_MISMATCHES]
    out["mismatches_truncated_in_json"] = len(full) > MAX_MISMATCHES
    out["equivalents"] = eq_full[:MAX_MISMATCHES]
    out["equivalents_truncated_in_json"] = len(eq_full) > MAX_MISMATCHES
    _pairs_tsv(full, Path(str(prefix) + ".mismatches.tsv"))
    _equiv_tsv(eq_full, Path(str(prefix) + ".equivalent.tsv"))
    text = summary_text(res)
    _atomic_write(Path(str(prefix) + ".txt"), text + "\n")
    _atomic_write(Path(str(prefix) + ".json"), json.dumps(out, indent=2) + "\n")  # completion marker, last
    print(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--test", required=True)
    ap.add_argument("--test-ref", required=True)
    ap.add_argument("--oracle", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    write_result(compare(args.test, args.test_ref, args.oracle, args.ref), args.out)


if __name__ == "__main__":
    main()
