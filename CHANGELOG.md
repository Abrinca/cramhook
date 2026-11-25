# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `cramhook fish` structure detection from a read's own segments (the record plus its `SA`
  entries; same-contig, same-strand segments of >= 500 bp that are back to back in the read, query
  gap or overlap <= max(200 bp, 2 % of the read)):
  - **`concatemer`**: a read of a circular contig that goes round it, `units` (summed aligned query
    length / contig length) >= 1.2; `period` = the contig length. A concatemer is explained by the
    assembly: it is counted but **not selected** by default; **`--keep-concatemers`** selects it
    with reason `concatemer`. A clip covering the next unit becomes a `concatemer` event (also at
    the contig end, before the origin-crossing skip) unless >= `--clip-min-len` of it stays
    unexplained; a remaining real clip keeps its `clip` row.
  - **`tandem`**: consecutive segments that jump back >= 1 kb on the reference (modulo the contig
    length on circular contigs, so a plain origin crossing is never a tandem); `period` = the
    jump, `units` = aligned query length / period >= 1.2. The assembly does not explain it (a
    collapsed tandem repeat, an excised circle), so it is **selected by default** (reason `tandem`,
    detail `period=P;units=U`); **`--no-tandems`** drops it without falling back to `clip`.
  - Precedence: `hairpin` > `concatemer` > `tandem` > origin-crossing skip; one structure per read.
    The label does not depend on `--no-clips`, `--scope` or `--no-hairpins` (those only gate
    selection), so `P.structure.tsv` can list more hairpin reads than `hairpin_detected` (which
    keeps its 0.1.2 meaning). A hairpin read's replaced clips revert to `clip` events.
    Neither reason counts against or is thinned by `--clip-cap`.
  - New `<prefix>.structure.tsv` (always written, header-only when empty): one row per hairpin,
    concatemer or tandem read, selected or not (`read, length, structure, contig, strand, period,
    units, arms, selected`).
  - `P.fish.json`: `concatemer_detected` / `tandem_detected` (`reads`), `reads.by_reason.concatemer`
    / `.tandem` when selected, the thresholds and `keep_concatemers` / `tandems` in `params`, and
    per contig `origin_crossing_reads` (reads with an origin-crossing clip on that contig that is
    not part of a concatemer / tandem chain: evidence it is circular; independent of `--no-clips`
    and `--scope`), `concatemer_reads` (reads with a circular-concatemer call on that contig) and
    `multimer` (circular contigs <= 50 kb: `reads`, `reads_monomer`,
    `reads_dimer`, `reads_trimer`, `reads_4plus_units`, `multimer_bp_frac`,
    `frac_multimer_given_long`; else `null`).
  - Hidden thresholds: `--concatemer-min-segment` (500), `--concatemer-min-units` (1.2),
    `--tandem-min-period` (1000), `--chain-max-gap` (200), `--chain-max-gap-frac` (0.02),
    `--concatemer-max-contig` (50000).

### Changed

- `cramhook fish`: the clips of a circular concatemer are no longer counted as skipped origin
  crossings (`origin_crossing_skipped`); they are `concatemer` events (selected only with
  `--keep-concatemers`). A tandem read that was selected as `clip` keeps that row for any clip the
  tandem does not explain, and its other clips are `tandem` rows.

### Output format

- Output format version remains `1`; every change is additive (`P.structure.tsv`, new `P.fish.json`
  keys, `concatemer` and `tandem` reasons).

## [0.1.3] - 2026-11-20

### Changed

- `cramhook qc`: a linear FASTA header (`circular=false` / `topology=linear`) now counts. Such a
  contig wraps only when reads span its end-to-start junction (verdict `circular`, `circular_by:
  verdict`); with an `unclear` or `linear` verdict it stays linear (`circular_by: header`).
  Previously the header was ignored and an `unclear` verdict wrapped the contig. Circular headers,
  `--linear` and contigs without a header flag are unchanged.
- `cramhook fish`: circularity follows the `--ref` FASTA header (`circular=` / `topology=`) like
  `qc`, so a linear header turns off the origin-crossing clip skip for that contig. `--linear`
  still overrides; contigs without a (consistent) header flag stay circular. `P.fish.json` gains
  `contigs` (`name`, `length_bp`, `circular`, `circular_by`: `flag`, `header` or `default`).

### Output format

- Output format version remains `1` (`format_version` in JSON, `##format_version` in VCF). The only
  schema change is additive: `P.fish.json` gains `contigs`. The `qc` linear-header change alters
  which contigs wrap, not the meaning of `circular` / `circular_by`.

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
