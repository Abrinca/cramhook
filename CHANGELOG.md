# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-10-19

### Added

- `cramhook lift`: translate long-read alignments from an old bacterial assembly to a new one
  (identity blocks, T1/T2/junction taint, clip screen) and remap the rest, writing one
  coordinate-sorted, indexed CRAM.
- `cramhook qc`: tiered site and region calls (`red`/`orange`/`yellow`/`gray`) where the reads
  disagree with the assembly, as VCF, BED, findings BED, bedGraph, candidates TSV, haplotypes TSV
  and JSON.
- `cramhook fish`: extract reads behind the calls (clips, large indels, unmapped, linked, locus),
  including `--bait` (lift against assembly + bait sequence) and `--bait-records any|primary`.
- Debug commands `cramhook diff` (identity blocks, taint, novel regions, junctions) and
  `cramhook translate` (classification only).

### Changed

- Renamed from cramlift to cramhook.

### Output format

- Output format version 1 (format_version in JSON, ##format_version in VCF).
