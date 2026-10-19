//! Integration tests for `cramhook diff` on the synthetic S1–S12 scenarios.
//!
//! Uses the library directly (`run_diff`) and writes outputs under
//! `target/test-diff/<S>/<label>/`. The truth in `truth.json` drives the
//! expectations; the invariant (every block is an exact match) is checked for
//! every scenario and is the strongest correctness signal.
//!
//! Missing fixtures are a hard error: generate them with `pixi run sim`.

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::path::PathBuf;

use serde::Deserialize;

use cramhook::{Block, Diff, Fasta, Preset, run_diff};

const SIDS: &[&str] = &[
    "S1", "S2", "S3", "S4", "S5", "S6", "S6b", "S6w", "S7", "S8", "S9", "S10a", "S10b", "S11",
    "S12",
];

/// Both read presets. The T2 probe and the invariant must hold for each.
const PRESETS: [(Preset, &str); 2] = [(Preset::MapHifi, "map-hifi"), (Preset::LrHq, "lr-hq")];

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

#[derive(Deserialize)]
struct Truth {
    #[allow(dead_code)]
    scenario: String,
    #[allow(dead_code)]
    description: String,
    #[serde(default)]
    edits: Vec<Edit>,
    repeats_old: Vec<Repeat>,
    contigs_old: BTreeMap<String, u64>,
    #[allow(dead_code)]
    contigs_new: BTreeMap<String, u64>,
}

#[derive(Deserialize)]
struct Edit {
    #[serde(rename = "type")]
    kind: String,
    #[serde(default)]
    old_contig: Option<String>,
    #[serde(default)]
    old_start: Option<u64>,
    #[serde(default)]
    old_end: Option<u64>,
    #[serde(default)]
    new_contig: Option<String>,
    #[serde(default)]
    new_start: Option<u64>,
    #[serde(default)]
    new_end: Option<u64>,
    #[serde(default)]
    offset: Option<u64>,
    #[allow(dead_code)]
    #[serde(default)]
    length: Option<u64>,
    #[allow(dead_code)]
    #[serde(default)]
    note: Option<String>,
}

#[derive(Deserialize)]
struct Repeat {
    name: String,
    contig: String,
    start: u64,
    end: u64,
    #[allow(dead_code)]
    #[serde(default)]
    strand: String,
}

struct Case {
    old: Fasta,
    new: Fasta,
    truth: Truth,
    diff: Diff,
}

/// Run one scenario with one preset. `tag` keeps parallel outputs separate.
/// Panics (rather than skipping) when the fixture is absent.
fn case(sid: &str, preset: Preset, tag: &str) -> Case {
    let base = root().join("data/synthetic");
    assert!(
        base.exists(),
        "data/synthetic is missing: run `pixi run sim` to generate the fixtures"
    );
    let dir = base.join(sid);
    let old_path = dir.join("old.fa");
    let new_path = dir.join("new.fa");
    assert!(
        old_path.exists() && new_path.exists(),
        "{sid}: missing old.fa/new.fa under {}: run `pixi run sim`",
        dir.display()
    );
    let out = root().join("target/test-diff").join(sid).join(tag);
    let run = run_diff(&old_path, &new_path, preset, 8, &out)
        .unwrap_or_else(|e| panic!("{sid} ({tag}): diff failed: {e:#}"));
    let old = Fasta::load(&old_path).unwrap();
    let new = Fasta::load(&new_path).unwrap();
    let truth: Truth =
        serde_json::from_str(&std::fs::read_to_string(dir.join("truth.json")).unwrap()).unwrap();
    Case { old, new, truth, diff: run.diff }
}

fn blocks_of<'a>(diff: &'a Diff, contig: &str) -> &'a [Block] {
    diff.blocks.get(contig).map(|v| v.as_slice()).unwrap_or(&[])
}

fn taint_positions(diff: &Diff, reason_t2: bool, contig: &str) -> BTreeSet<u64> {
    let map = if reason_t2 { &diff.taint_t2 } else { &diff.taint_t1 };
    let mut s = BTreeSet::new();
    if let Some(v) = map.get(contig) {
        for iv in v {
            for p in iv.start.0..iv.end.0 {
                s.insert(p);
            }
        }
    }
    s
}

fn novel_positions(diff: &Diff, contig: &str) -> BTreeSet<u64> {
    let mut s = BTreeSet::new();
    if let Some(v) = diff.novel.get(contig) {
        for iv in v {
            for p in iv.start.0..iv.end.0 {
                s.insert(p);
            }
        }
    }
    s
}

/// The invariant: the block set must be non-empty, cover a positive number of
/// old bases, and every block must be an exact (reverse-complemented for `-`)
/// match between the assemblies. Returns the number of bad blocks.
fn check_invariant(case: &Case) -> usize {
    let total: u64 = case.diff.blocks.values().flatten().map(|b| b.len()).sum();
    assert!(!case.diff.blocks.is_empty(), "no identity blocks");
    assert!(total > 0, "zero total block coverage");
    let mut bad = 0;
    for (contig, v) in &case.diff.blocks {
        let old = case.old.seq(contig).expect("old contig");
        for b in v {
            let o = &old[b.old.start.0 as usize..b.old.end.0 as usize];
            let n = case.new.seq(&b.new_contig).expect("new contig");
            let nv = &n[b.new_start.0 as usize..(b.new_start.0 + b.old.len()) as usize];
            let ok = match b.strand {
                cramhook::Strand::Plus => o == nv,
                cramhook::Strand::Minus => o == &revcomp(nv),
            };
            if !ok {
                bad += 1;
                eprintln!(
                    "invariant fail on {contig} {}..{} -> {} {} {}",
                    b.old.start.0, b.old.end.0, b.new_contig, b.new_start.0, b.strand.as_char()
                );
            }
        }
    }
    bad
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// normalisation of indel edges (shared with the test expectations)

/// Positions an insertion of `x` originally at `p` may occupy after left/right
/// normalisation on the old sequence (bcftools-style rotation).
fn ins_positions(old: &[u8], p: u64, x: &[u8]) -> BTreeSet<u64> {
    let l = x.len() as u64;
    let mut poss = BTreeSet::new();
    if l == 0 {
        return poss;
    }
    let mut seen: HashSet<(u64, u64)> = HashSet::new();
    let mut stack = vec![(p, 0u64)];
    seen.insert((p, 0));
    while let Some((pos, r)) = stack.pop() {
        poss.insert(pos);
        let last = x[((r + l - 1) % l) as usize];
        if pos > 0 && old[(pos - 1) as usize] == last {
            let st = (pos - 1, (r + l - 1) % l);
            if seen.insert(st) {
                stack.push(st);
            }
        }
        let first = x[r as usize];
        if pos < old.len() as u64 && old[pos as usize] == first {
            let st = (pos + 1, (r + 1) % l);
            if seen.insert(st) {
                stack.push(st);
            }
        }
    }
    poss
}

/// Start positions a deletion of `old[a..b)` may occupy after normalisation.
fn del_starts(old: &[u8], a: u64, b: u64) -> BTreeSet<u64> {
    let l = b - a;
    let mut starts = BTreeSet::new();
    let mut seen: HashSet<u64> = HashSet::new();
    let mut stack = vec![a];
    seen.insert(a);
    while let Some(aa) = stack.pop() {
        starts.insert(aa);
        let bb = aa + l;
        if aa > 0 && old[(aa - 1) as usize] == old[(bb - 1) as usize] && seen.insert(aa - 1) {
            stack.push(aa - 1);
        }
        if bb < old.len() as u64 && old[bb as usize] == old[aa as usize] && seen.insert(aa + 1) {
            stack.push(aa + 1);
        }
    }
    starts
}

/// The rotation edge in old coordinates: the first old position whose
/// (indel-adjusted) coordinate reaches the rotation offset.
fn rotation_edge(old_len: u64, k: u64, edits: &[Edit], contig: &str) -> u64 {
    let events: Vec<(u64, i64, bool)> = edits
        .iter()
        .filter(|e| e.old_contig.as_deref() == Some(contig))
        .filter_map(|e| match e.kind.as_str() {
            "insertion" => {
                let l = (e.new_end.unwrap() - e.new_start.unwrap()) as i64;
                Some((e.old_start.unwrap(), l, false))
            }
            "deletion" => {
                let l = (e.old_end.unwrap() - e.old_start.unwrap()) as i64;
                Some((e.old_end.unwrap(), -l, true))
            }
            _ => None,
        })
        .collect();
    let f = |x: u64| -> i64 {
        let mut d = 0i64;
        for &(t, delta, inc) in &events {
            if if inc { x >= t } else { x > t } {
                d += delta;
            }
        }
        x as i64 + d
    };
    for o in 0..=old_len {
        if f(o) >= k as i64 {
            return o;
        }
    }
    old_len
}

/// Internal block starts and ends (contig-end boundaries excluded).
fn block_boundaries(diff: &Diff, contig: &str, len: u64) -> (BTreeSet<u64>, BTreeSet<u64>) {
    let mut starts = BTreeSet::new();
    let mut ends = BTreeSet::new();
    for b in blocks_of(diff, contig) {
        if b.old.start.0 > 0 && b.old.start.0 < len {
            starts.insert(b.old.start.0);
        }
        if b.old.end.0 > 0 && b.old.end.0 < len {
            ends.insert(b.old.end.0);
        }
    }
    (starts, ends)
}

/// All valid edge positions for every planted edit on `contig`.
fn valid_edges(case: &Case, contig: &str) -> BTreeSet<u64> {
    let old = case.old.seq(contig).unwrap();
    let mut valid = BTreeSet::new();
    for e in &case.truth.edits {
        match e.kind.as_str() {
            "snp" if e.old_contig.as_deref() == Some(contig) => {
                valid.insert(e.old_start.unwrap());
                valid.insert(e.old_start.unwrap() + 1);
            }
            "insertion" if e.old_contig.as_deref() == Some(contig) => {
                let nc = e.new_contig.as_deref().unwrap();
                let x = &case.new.seq(nc).unwrap()
                    [e.new_start.unwrap() as usize..e.new_end.unwrap() as usize];
                valid.extend(ins_positions(old, e.old_start.unwrap(), x));
            }
            "deletion" if e.old_contig.as_deref() == Some(contig) => {
                let a = e.old_start.unwrap();
                let b = e.old_end.unwrap();
                for s in del_starts(old, a, b) {
                    valid.insert(s);
                    valid.insert(s + (b - a));
                }
            }
            "rotation" if e.old_contig.as_deref() == Some(contig) => {
                valid.insert(rotation_edge(
                    old.len() as u64,
                    e.offset.unwrap(),
                    &case.truth.edits,
                    contig,
                ));
            }
            _ => {}
        }
    }
    valid
}

fn edits_on<'a>(case: &'a Case, contig: &str) -> Vec<&'a Edit> {
    case.truth
        .edits
        .iter()
        .filter(|e| e.old_contig.as_deref() == Some(contig))
        .collect()
}

/// Insertion sequence `x` (from the new assembly) for an insertion edit.
fn insertion_seq<'a>(case: &'a Case, e: &Edit) -> &'a [u8] {
    let nc = e.new_contig.as_deref().unwrap();
    &case.new.seq(nc).unwrap()[e.new_start.unwrap() as usize..e.new_end.unwrap() as usize]
}

// ---------------------------------------------------------------------------
// tests

#[test]
fn all_scenarios_invariant() {
    for (preset, pname) in PRESETS {
        for sid in SIDS {
            let c = case(sid, preset, &format!("inv-{pname}"));
            assert_eq!(check_invariant(&c), 0, "{sid} ({pname}): invariant violations");
        }
    }
}

#[test]
fn s1_identical() {
    let c = case("S1", Preset::MapHifi, "s1");
    assert_eq!(check_invariant(&c), 0);
    for (contig, len) in &c.truth.contigs_old {
        let bl = blocks_of(&c.diff, contig);
        assert_eq!(bl.len(), 1, "{contig}: expected one block, got {}", bl.len());
        let b = &bl[0];
        assert_eq!((b.old.start.0, b.old.end.0), (0, *len), "{contig}");
        assert_eq!(b.new_contig, *contig, "{contig}");
        assert_eq!(b.new_start.0, 0, "{contig}");
        assert_eq!(b.strand, cramhook::Strand::Plus, "{contig}");
    }
    assert!(c.diff.taint_t1.values().all(|v| v.is_empty()), "S1 T1 taint");
    assert!(c.diff.taint_t2.values().all(|v| v.is_empty()), "S1 T2 taint");
    assert!(c.diff.novel.values().all(|v| v.is_empty()), "S1 novel");
}

#[test]
fn s2_rotation() {
    let c = case("S2", Preset::MapHifi, "s2");
    assert_eq!(check_invariant(&c), 0);
    let k = c.truth.edits[0].offset.unwrap();
    let len = c.truth.contigs_old["chr"];
    let bl = blocks_of(&c.diff, "chr");
    assert_eq!(bl.len(), 2, "S2 chr should have 2 blocks: {bl:?}");
    for b in bl {
        assert_eq!(b.strand, cramhook::Strand::Plus);
        assert_eq!(b.new_start.0, (b.old.start.0 + len - k) % len, "block {b:?}");
        for p in [b.old.start.0, b.old.end.0 - 1] {
            let mapped = b.map_pos(cramhook::OldPos(p)).unwrap().0;
            assert_eq!(mapped, (p + len - k) % len, "p={p}");
        }
    }
}

#[test]
fn s3_reverse_complement() {
    let c = case("S3", Preset::MapHifi, "s3");
    assert_eq!(check_invariant(&c), 0);
    let len = c.truth.contigs_old["p1"];
    let bl = blocks_of(&c.diff, "p1");
    assert_eq!(bl.len(), 1, "S3 p1 blocks: {bl:?}");
    let b = &bl[0];
    assert_eq!(b.strand, cramhook::Strand::Minus);
    assert_eq!(b.new_contig, "p1");
    for p in [0u64, 1, len / 2, len - 1] {
        assert_eq!(b.map_pos(cramhook::OldPos(p)).unwrap().0, len - 1 - p, "p={p}");
    }
}

#[test]
fn s4_rename() {
    let c = case("S4", Preset::MapHifi, "s4");
    assert_eq!(check_invariant(&c), 0);
    assert!(!blocks_of(&c.diff, "p1").is_empty());
    for b in blocks_of(&c.diff, "p1") {
        assert_eq!(b.new_contig, "plasmid_A");
    }
}

#[test]
fn s5_unique_plasmid_added() {
    for (preset, pname) in PRESETS {
        let c = case("S5", preset, &format!("s5-{pname}"));
        assert_eq!(check_invariant(&c), 0, "S5 ({pname})");
        let p3 = c.truth.contigs_new["p3"];
        let novel = novel_positions(&c.diff, "p3");
        assert_eq!(novel, (0..p3).collect::<BTreeSet<_>>(), "S5 ({pname}): novel != whole p3");
        assert!(
            c.diff.taint_t2.values().all(|v| v.is_empty()),
            "S5 ({pname}): T2 must be empty"
        );
    }
}

#[test]
fn s6_s6b_repeat_copy_taint() {
    for (preset, pname) in PRESETS {
        for sid in ["S6", "S6b", "S6w"] {
            let c = case(sid, preset, &format!("s6-{pname}"));
            assert_eq!(check_invariant(&c), 0, "{sid} ({pname})");
            for rep in &c.truth.repeats_old {
                if !rep.name.starts_with("IS") {
                    continue;
                }
                let t2 = taint_positions(&c.diff, true, &rep.contig);
                let covered = (rep.start..rep.end).filter(|p| t2.contains(p)).count() as u64;
                let frac = covered as f64 / (rep.end - rep.start) as f64;
                println!("{sid} ({pname}) {} {}: T2 coverage {:.4}", rep.name, rep.contig, frac);
                assert!(
                    frac >= 0.9,
                    "{sid} ({pname}) {} on {}: T2 only covers {:.4} of the copy",
                    rep.name,
                    rep.contig,
                    frac
                );
            }
        }
    }
}

#[test]
fn s7_removed_plasmid_taint() {
    let c = case("S7", Preset::MapHifi, "s7");
    assert_eq!(check_invariant(&c), 0);
    let len = c.truth.contigs_old["p2"];
    let t1 = taint_positions(&c.diff, false, "p2");
    assert_eq!(t1, (0..len).collect::<BTreeSet<_>>(), "S7 T1 must cover all of p2");
}

#[test]
fn s10_duplication_taint() {
    for (preset, pname) in PRESETS {
        for sid in ["S10a", "S10b"] {
            let c = case(sid, preset, &format!("s10-{pname}"));
            assert_eq!(check_invariant(&c), 0, "{sid} ({pname})");
            let ed = c
                .truth
                .edits
                .iter()
                .find(|e| e.kind == "collapsed_duplication" || e.kind == "expanded_duplication")
                .expect("duplication edit");
            let (s, e) = (ed.old_start.unwrap(), ed.old_end.unwrap());
            let mut taint = taint_positions(&c.diff, false, "chr");
            taint.extend(taint_positions(&c.diff, true, "chr"));
            let overlap = (s..e).filter(|p| taint.contains(p)).count() as u64;
            println!("{sid} ({pname}) duplicated segment {s}..{e}: taint overlap {overlap}/{}", e - s);
            assert!(overlap > 0, "{sid} ({pname}): no taint overlaps the duplicated segment");
        }
    }
}

/// Blocks must cover every old base except the edit bases (SNP bases, deleted
/// bases); insertions add new sequence only.
fn check_block_coverage(case: &Case, contig: &str, len: u64) {
    let mut ivs: Vec<(u64, u64)> =
        blocks_of(&case.diff, contig).iter().map(|b| (b.old.start.0, b.old.end.0)).collect();
    ivs.sort_unstable();
    let total: u64 = ivs.iter().map(|(a, b)| b - a).sum();
    let mut union = 0u64;
    let mut cur: Option<(u64, u64)> = None;
    for (a, b) in &ivs {
        match cur {
            None => cur = Some((*a, *b)),
            Some((ca, cb)) if *a <= cb => cur = Some((ca, cb.max(*b))),
            Some((ca, cb)) => {
                union += cb - ca;
                cur = Some((*a, *b));
            }
        }
    }
    if let Some((ca, cb)) = cur {
        union += cb - ca;
    }
    let removed: u64 = edits_on(case, contig)
        .iter()
        .map(|e| match e.kind.as_str() {
            "snp" => 1,
            "deletion" => e.old_end.unwrap() - e.old_start.unwrap(),
            _ => 0,
        })
        .sum();
    let expected = len - removed;
    assert_eq!(union, expected, "{contig}: blocks do not cover all non-edit old bases");
    assert_eq!(total, expected, "{contig}: blocks overlap or leave gaps");
}

#[test]
fn s8_s9_s11_s12_edges_match_edits() {
    for sid in ["S8", "S9", "S11", "S12"] {
        let c = case(sid, Preset::MapHifi, "edges");
        assert_eq!(check_invariant(&c), 0, "{sid}");
        for (contig, len) in &c.truth.contigs_old {
            check_block_coverage(&c, contig, *len);
            let (starts, ends) = block_boundaries(&c.diff, contig, *len);
            let obs: BTreeSet<u64> = starts.union(&ends).copied().collect();
            let valid = valid_edges(&c, contig);
            let unexplained: Vec<u64> = obs.difference(&valid).copied().collect();
            assert!(
                unexplained.is_empty(),
                "{sid} {contig}: block edges not explained by a planted edit: {unexplained:?}\nvalid={valid:?}"
            );

            for e in edits_on(&c, contig) {
                match e.kind.as_str() {
                    "snp" => {
                        let p = e.old_start.unwrap();
                        assert!(
                            ends.contains(&p) && starts.contains(&(p + 1)),
                            "{sid} {contig}: SNP at {p} must end a block at {p} and start one at {} (starts={starts:?}, ends={ends:?})",
                            p + 1
                        );
                    }
                    "insertion" => {
                        let p = e.old_start.unwrap();
                        let x = insertion_seq(&c, e);
                        let positions = ins_positions(c.old.seq(contig).unwrap(), p, x);
                        assert!(
                            positions.iter().any(|q| starts.contains(q) && ends.contains(q)),
                            "{sid} {contig}: insertion at {p} has no position where a block both ends and starts (reachable={positions:?}, starts={starts:?}, ends={ends:?})"
                        );
                    }
                    "deletion" => {
                        let (a, b) = (e.old_start.unwrap(), e.old_end.unwrap());
                        let l = b - a;
                        let sts = del_starts(c.old.seq(contig).unwrap(), a, b);
                        assert!(
                            sts.iter().any(|s| ends.contains(s) && starts.contains(&(s + l))),
                            "{sid} {contig}: deletion {a}..{b} has no consistent normalised edge pair (reachable={sts:?}, starts={starts:?}, ends={ends:?})"
                        );
                    }
                    "rotation" => {
                        let r = rotation_edge(
                            c.old.seq(contig).unwrap().len() as u64,
                            e.offset.unwrap(),
                            &c.truth.edits,
                            contig,
                        );
                        assert!(
                            ends.contains(&r) && starts.contains(&r),
                            "{sid} {contig}: rotation edge {r} does not end and start a block (starts={starts:?}, ends={ends:?})"
                        );
                    }
                    _ => {}
                }
            }
        }
    }
}

#[test]
fn s6b_report_coverage() {
    // Always prints the S6b coverage numbers for both presets; the assertion
    // lives in s6_s6b_repeat_copy_taint.
    for (preset, pname) in PRESETS {
        let c = case("S6b", preset, &format!("s6b-report-{pname}"));
        for rep in &c.truth.repeats_old {
            if !rep.name.starts_with("IS") {
                continue;
            }
            let t2 = taint_positions(&c.diff, true, &rep.contig);
            let covered = (rep.start..rep.end).filter(|p| t2.contains(p)).count() as u64;
            println!(
                "S6b ({pname}) {} {}: {}/{} = {:.4}",
                rep.name,
                rep.contig,
                covered,
                rep.end - rep.start,
                covered as f64 / (rep.end - rep.start) as f64
            );
        }
    }
}
