# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1] - 2026-10-21

### Added

- `cramhook fish --clip-cap N` / `--clip-bin BP` (default off / 2000): cap the clip-only selected
  reads per reference bin, for high-depth samples where a repeat or misjoin contributes thousands of
  clipped reads. A read selected by any other reason is never dropped; within a bin the `N` longest
  clips are kept (ties by read name), so the choice is identical for every `--threads` and input
  order. On hotspot samples this removes most of the clip volume (measured: 36 % fewer ONT selected
  reads at cap 10 on the spiked eval) while every planted event kept at least `N` reads on it. The
  JSON reports `cap` (`reads`, `bp`) and `params.clip_cap`/`clip_bin`.

### Fixed

- `docs/output-format.md` said the current format version was 2 while every output carries 1.

## [0.1.0] - 2026-10-19

### Added

- `cramhook lift`: translate long-read alignments from an old bacterial assembly to a new one
  (identity blocks, T1/T2/junction taint, clip screen) and remap the rest, writing one
  coordinate-sorted, indexed CRAM.
- `cramhook qc`: tiered site and region calls (`red`/`orange`/`yellow`/`gray`) where the reads
  disagree with the assembly, as VCF, BED, findings BED, bedGraph, candidates TSV, haplotypes TSV
  and JSON.
- `cramhook qc`: `q.gc.bedgraph`, a GC-fraction track in fixed `--gc-window` windows (default
  100 bp; `0` disables the file), computed from the `--ref` FASTA qc already loads. Non-ACGT
  bases are excluded from both counts; the JSON gains `gc_window`.
- `cramhook qc --crispr`: annotate CRISPR arrays detected in the reference with `pilercr-rs`
  (MIT, a Rust port of PILER-CR): untiered `crispr` BED rows and INFO `CRISPR` flags on sites
  inside an array; the top-level `crispr` JSON block (`arrays`, `bp`) is written either way
  (zeros without the flag). Annotation only: arrays never change a TIER, are not finding members
  and are not counted in `alerts`.
- `cramhook fish`: extract reads behind the calls (clips, large indels, unmapped, linked, locus),
  including `--bait` (lift against assembly + bait sequence) and `--bait-records any|primary`.
- Debug commands `cramhook diff` (identity blocks, taint, novel regions, junctions) and
  `cramhook translate` (classification only).
- `--progress` on `lift`, `qc` and `fish`: machine-readable progress on stderr, one JSON object per
  line (`plan`, per-phase `start`/`update`/`end`, `done` or `error`; stream version 1, see
  [docs/output-format.md](docs/output-format.md#progress-stream---progress)). `update` lines carry
  `frac`, one line per whole percent at most, and in record passes `records`, `records_total` and
  `total_exact` (alignment records; exact totals from the CRAM container headers or the BAM index;
  with an exact total, `frac` also follows the records read, so the unmapped tail reaches 100 %).
  Off by default; without it the outputs are unchanged.

### Changed

- Renamed from cramlift to cramhook.

### Output format

- Output format version 1 (`format_version` in JSON, `##format_version` in VCF). This release adds
  the qc `P.gc.bedgraph` file and the qc JSON `gc_window` key, plus the `crispr` BED label and VCF
  INFO `CRISPR` flag with `--crispr`. The qc JSON `crispr` block and `params.crispr` key are not
  conditional on the flag: they are always present (zeros / `false` without `--crispr`), so the
  schema is stable.
