# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.2] - 2026-11-18

### Added

- `cramhook fish` hairpin (fold-back) labelling: when a record's terminal soft clip matches an
  opposite-strand `SA` segment on the same contig with >= 50 % reference overlap of the shorter
  arm, adjacent original-orientation query intervals (`|gap| <= max(200 bp, 2% of read length)`),
  and each arm >= 20 % of the read and >= `--clip-min-len`, the clip becomes a `hairpin` event
  with arm intervals in `detail` (`arm1=q<s>-<e>,r<s>-<e>,<strand>;arm2=…`) unless >=
  `--clip-min-len` of it stays unexplained. Hairpin rows are selected by default and are not thinned
  by `--clip-cap`; `--no-hairpins` drops them without falling back to `clip`. Detection thresholds
  are exposed as hidden flags (`--hairpin-min-overlap`, `--hairpin-max-gap`,
  `--hairpin-max-gap-frac`, `--hairpin-min-arm-frac`). `P.fish.json` reports
  `reads.by_reason.hairpin` and `hairpin_detected` (`reads`, `records`). Measured: ONT FAM22000
  flagged 127/140 mapped fold-back reads (plus 3 extra hits confirmed by self-alignment); HiFi
  pb-FAM22000 flagged 7/7 known fold-backs.
- `cramhook fish` `<prefix>.lengths.tsv`: always-written read-length histogram (`bin_lo`, `reads`,
  `bp`) with `[0, 100)` followed by log-spaced bins at 20 per decade (`round(100 * 10^(k/20))`) up to
  the bin holding the longest read. Each read counts once from its primary or unmapped record
  (skipping secondary, supplementary, duplicate and QC-fail), so `sum(reads)` and `sum(bp)` match
  `P.fish.json` `records.primary + records.unmapped` and `mapped_bp + unmapped_bp`.
- `P.fish.json` `origin_crossing_skipped` (`reads`, `records`, `clips`, `clips_exact`,
  `clips_near_end`): counts origin-crossing clips skipped on circular contigs (before `--scope`
  filtering).

### Changed

- `cramhook fish --circ-end-tol BP` (default `100`, capped at half the contig length; `0` = exact
  ends only): on circular contigs, terminal clips stopping within `BP` of a contig end are now also
  skipped as origin crossings when `SA` has a same-contig, same-strand segment within `BP` of the
  opposite end (linear contigs are unchanged; reported in `P.fish.json` `params.circ_end_tol`).
  This changes default `fish` clip selection: on ONT FAM22000, the uncapped clip pool drops from
  2,030 to 872 reads, with all 1,158 removed reads being near-end wrap-arounds.

### Fixed

- `cramhook fish` and `cramhook qc` on an indexed BAM with zero records: htslib fails to seek to
  coordinate-less records (`*`) on an empty BAI; both commands now detect an empty index, emit
  `warning: <in>: input has no records` once on stderr, and write their normal empty/zero outputs.

### Output format

- Output format version remains `1` (`format_version` in JSON, `##format_version` in VCF). This
  release makes additive `fish` changes: the always-written `P.lengths.tsv` histogram, the `hairpin`
  reason in `P.fish.tsv`, and new `P.fish.json` fields (`hairpin_detected`,
  `origin_crossing_skipped`, `reads.by_reason.hairpin`, and `params.{hairpins, circ_end_tol,
  hairpin_min_overlap, hairpin_max_gap, hairpin_max_gap_frac, hairpin_min_arm_frac}`).
- `docs/reference/output-formats.md`: the format version is bumped only by breaking changes; additive
  changes (new keys, columns, reasons or files) keep it. The doc said every change bumps it, which
  0.1.1 and this release do not.

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

- `docs/reference/output-formats.md` said the current format version was 2 while every output carries 1.

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
  [docs/reference/output-formats.md](https://github.com/Abrinca/cramhook/blob/main/docs/reference/output-formats.md#progress-stream---progress)). `update` lines carry
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
