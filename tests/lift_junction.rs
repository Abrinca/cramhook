//! Phase-1 regressions for the lift junction-taint bug (manager spec
//! `scratch/manager/lift-junction-taint-spec.md`).
//!
//! T2 taint is projected from *novel* new positions: those not covered exactly
//! once by identity blocks (`src/diff.rs`, `compute_diff_with` step 5). An edit
//! that makes existing sequence more similar to a read **without adding new
//! bases** (a deletion that joins a fragmented copy, a contig join) adds no
//! novel interval, so T2 never probes the junction. A read on the unchanged
//! duplicate copy is then translated with its stale old-assembly MAPQ instead of
//! being remapped, while a full remap gives MAPQ 0 (the joined copy is now a
//! competent competitor).
//!
//! Each fixture is a tiny assembly pair plus one perfect read. `old.bam` is a
//! real minimap2 mapping to the old assembly; `oracle.cram` is a real full remap
//! to the new assembly. The test lifts `old.bam` and requires the result to equal
//! the oracle (contig/pos/strand/MAPQ/CIGAR). Regenerate fixtures with
//! `pixi run python scratch/lift-junction/make_fixtures.py`.
//!
//! The `deletion_join` and `contig_join` tests reproduced the bug before the T2
//! junction taint; they now pass.

use std::path::{Path, PathBuf};

use rust_htslib::bam::{self, Read};

use cramhook::diff::{DiffParams, Preset, run_diff_with};
use cramhook::{LiftOptions, TranslateConfig, lift};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn case_dir(name: &str) -> PathBuf {
    let d = root().join("tests/fixtures/lift-junction").join(name);
    assert!(d.join("old.bam").exists(), "missing fixture {name}: run make_fixtures.py");
    d
}

/// All primary records as `(contig, pos, reverse, mapq, cigar)`, sorted.
fn primary_records(path: &Path, ref_fa: &Path) -> Vec<(String, i64, bool, u8, String)> {
    let mut reader = bam::Reader::from_path(path).unwrap();
    reader.set_reference(ref_fa).unwrap();
    let names: Vec<String> = reader
        .header()
        .target_names()
        .iter()
        .map(|n| String::from_utf8_lossy(n).to_string())
        .collect();
    let mut out = Vec::new();
    for rec in reader.records() {
        let rec = rec.unwrap();
        if rec.flags() & (0x100 | 0x800) != 0 {
            continue;
        }
        let contig = names.get(rec.tid() as usize).cloned().unwrap_or_default();
        out.push((
            contig,
            rec.pos(),
            rec.is_reverse(),
            rec.mapq(),
            format!("{}", rec.cigar()),
        ));
    }
    out.sort();
    out
}

/// Lift the fixture and return the primary records of the output plus the stats.
fn lift_fixture(name: &str) -> (Vec<(String, i64, bool, u8, String)>, cramhook::Stats) {
    let d = case_dir(name);
    let out_dir = root().join("target/test-lift-junction").join(name);
    std::fs::create_dir_all(&out_dir).unwrap();
    let out = out_dir.join("out.cram");
    let cfg = TranslateConfig { xl_tag: true, ..Default::default() };
    let stats = lift(
        &d.join("old.bam"),
        &d.join("old.fa"),
        &d.join("new.fa"),
        &out,
        None,
        2,
        &cfg,
        &LiftOptions { strict: false, keep_tmp: false, min_input_mb: Some(0.0) },
        None,
    )
    .unwrap_or_else(|e| panic!("{name}: lift failed: {e:#}"));
    (primary_records(&out, &d.join("new.fa")), stats)
}

/// Assert the lifted record set equals the full-remap oracle (the bug makes
/// these differ in MAPQ, and for a contig join in contig too).
fn assert_matches_oracle(name: &str) {
    let d = case_dir(name);
    let (got, _) = lift_fixture(name);
    let want = primary_records(&d.join("oracle.cram"), &d.join("new.fa"));
    assert!(!want.is_empty(), "{name}: oracle has no records");
    assert_eq!(got, want, "{name}: lifted != full remap");
}

/// A deletion of two separators joins B's fragmented copy of the read, so new B
/// becomes as good as the unchanged A copy. Full remap: MAPQ 0. Before the T2
/// junction taint, lift translated the A read at MAPQ 60; now the junction taint
/// defers it and the output equals the oracle.
#[test]
fn lift_deletion_join_matches_oracle() {
    assert_matches_oracle("deletion_join");
}

/// Two old contigs X and Y are joined in the new assembly, reconstructing the
/// read's sequence on X. Full remap moves the read to X at MAPQ 0; before the
/// junction taint lift kept it on the unchanged A copy at MAPQ 60.
#[test]
fn lift_contig_join_matches_oracle() {
    assert_matches_oracle("contig_join");
}

/// P1a: the competitor copy lives on a **circular** old contig and spans its
/// origin, so in the linear old index the read's two halves are split; rotating
/// the contig makes the copy contiguous. minimap2 indexes contigs linearly, so
/// a wrap is a junction like any other (no circular-header skip). A full remap
/// gives MAPQ 0; the junction taint must defer the read.
#[test]
fn lift_origin_competitor_matches_oracle() {
    let (got, stats) = lift_fixture("origin_competitor");
    let d = case_dir("origin_competitor");
    assert_eq!(got, primary_records(&d.join("oracle.cram"), &d.join("new.fa")));
    assert!(
        stats.deferred_by_reason.get("taint_junction").copied().unwrap_or(0) >= 1,
        "origin-spanning competitor not caught: {:?}",
        stats.deferred_by_reason
    );
}

/// The reverse-complemented variant (the read's competitor is on the new `-`
/// strand) must be caught too; before the fix it also gave a contig/strand
/// disagreement.
#[test]
fn lift_origin_competitor_minus_matches_oracle() {
    let (got, stats) = lift_fixture("origin_competitor_minus");
    let d = case_dir("origin_competitor_minus");
    assert_eq!(got, primary_records(&d.join("oracle.cram"), &d.join("new.fa")));
    assert!(
        stats.deferred_by_reason.get("taint_junction").copied().unwrap_or(0) >= 1,
        "reverse-complemented competitor not caught: {:?}",
        stats.deferred_by_reason
    );
}

/// P1b: a block nested inside another (`[0,20000)`, `[1000,2000)`,
/// `[20000,40000)`) hid the `[0,20000)->[20000,40000)` adjacency from a
/// windows(2) scan over start-sorted blocks. The endpoint-event search finds it.
/// End-to-end the read on the unchanged copy was translated (MAPQ 60) where the
/// oracle gives 0.
#[test]
fn lift_nested_junction_matches_oracle() {
    let (got, stats) = lift_fixture("nested_junction");
    let d = case_dir("nested_junction");
    assert_eq!(got, primary_records(&d.join("oracle.cram"), &d.join("new.fa")));
    assert!(
        stats.deferred_by_reason.get("taint_junction").copied().unwrap_or(0) >= 1,
        "nested-block junction not caught: {:?}",
        stats.deferred_by_reason
    );
}

/// P2a: all new contigs are treated as circular, so the last-block→first-block
/// adjacency at the new origin is a junction too. The fixture's new contig is a
/// 30 kb arc of the 40 kb old contig: the seam joins old[10000)->old[20000)
/// (the deletion across the new origin), alongside the wrap at new 20000.
#[test]
fn diff_seam_junction_found() {
    let fx = root().join("tests/fixtures/lift-junction/seam_junction");
    let tmp = root().join("target/test-lift-junction/seam");
    std::fs::create_dir_all(&tmp).unwrap();
    let run = run_diff_with(
        &fx.join("old.fa"),
        &fx.join("new.fa"),
        Preset::LrHq,
        2,
        &tmp.join("d"),
        DiffParams::default(),
    )
    .unwrap();
    assert_eq!(run.summary.junctions_detected, 2, "expected the wrap and the seam junction");
    assert_eq!(run.summary.junction_rotation_wraps, 1);
    let ivs = &run.diff.junctions["C"];
    assert!(
        ivs.iter().any(|iv| iv.start.0 == 0 && iv.end.0 == 50),
        "seam junction at the new origin missing: {ivs:?}"
    );
    assert!(
        ivs.iter().any(|iv| iv.start.0 == 19950 && iv.end.0 == 20050),
        "wrap junction missing: {ivs:?}"
    );
}

/// Safety of the self-hit exclusion: excluding the window that aligns back to
/// the junction's own source must not leave the source side unhandled. A read
/// that spans the deleted sequence on the source copy is deferred by the
/// deleted-sequence rules (`taint_t1`/`no_block`), and the read on the unchanged
/// competitor copy is caught by the junction taint; both match the full remap.
#[test]
fn lift_deletion_source_safety() {
    let (got, stats) = lift_fixture("deletion_source_safety");
    let d = case_dir("deletion_source_safety");
    assert_eq!(got, primary_records(&d.join("oracle.cram"), &d.join("new.fa")));
    assert!(
        stats.deferred_by_reason.get("taint_junction").copied().unwrap_or(0) >= 1,
        "competitor read not caught: {:?}",
        stats.deferred_by_reason
    );
    let source_handled = stats.deferred_by_reason.get("taint_t1").copied().unwrap_or(0)
        + stats.deferred_by_reason.get("no_block").copied().unwrap_or(0);
    assert!(
        source_handled >= 1,
        "source-spanning read not deferred by the deleted-sequence rules: {:?}",
        stats.deferred_by_reason
    );
}

/// A block *move* breaks asm5 colinearity, so the moved block cannot be matched
/// at its new place and the aligner leaves a novel gap there. T2 therefore taints
/// the junction and the read is remapped: this class does not escape. Kept as a
/// documented non-escape and a guard that the fix does not over-taint it.
#[test]
fn lift_rearrange_join_matches_oracle() {
    assert_matches_oracle("rearrange_join");
}

/// A pure circular rotation is a representation change: no new competitor, so
/// lift must equal the full remap. The junction taint (phase 2) must not create
/// a junction at a rotation origin.
#[test]
fn lift_rotation_control_matches_oracle() {
    assert_matches_oracle("rotation_control");
}

/// A circular wrap is a junction like any other, independent of the old FASTA
/// header (minimap2 indexes the old contig linearly, so a copy split by the old
/// origin can become contiguous). For the pure rotation fixture the junction
/// window maps back only to its own source, so the self-hit exclusion leaves no
/// junction taint.
#[test]
fn diff_rotation_wrap_is_a_junction() {
    let fx = root().join("tests/fixtures/lift-junction/rotation_control");
    let tmp = root().join("target/test-lift-junction/rotation_wrap");
    std::fs::create_dir_all(&tmp).unwrap();
    let run = run_diff_with(
        &fx.join("old.fa"),
        &fx.join("new.fa"),
        Preset::LrHq,
        2,
        &tmp.join("d"),
        DiffParams::default(),
    )
    .unwrap();
    assert!(run.summary.junctions_detected >= 1, "rotation wrap not found");
    assert_eq!(
        run.summary.junctions_detected, run.summary.junction_rotation_wraps,
        "all junctions here should be wraps"
    );
    assert!(run.diff.junctions.values().any(|v| !v.is_empty()), "wrap interval missing");
    assert!(
        !run.diff.taint_junction.values().any(|v| !v.is_empty()),
        "pure rotation: only self-hits, no junction taint expected"
    );
}
