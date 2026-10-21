# Output formats

Reference for every public output file of `qc`, `fish` and `lift`. The VCF and BED files are read
by igv.js (and any genome browser); the JSON files are meant for programs.

## Format version

All outputs share one format version number, currently **1**:

- every JSON file carries `"format_version": 1` as its first key;
- the qc VCF carries a `##format_version=1` header line;
- BED, bedGraph and TSV files share the same version number but carry no version line of their own
  (BED and TSV files start with one `#`-prefixed or plain column-name line; the bedGraph has none).

Any change to fields, keys, columns or the set of output files bumps the version and is recorded in
[CHANGELOG.md](../CHANGELOG.md). Coordinates in BED, bedGraph and TSV files are 0-based, half-open;
VCF `POS` is 1-based. One exception: the `sa` column of `P.fish.tsv` copies positions from the SAM
`SA` tag unchanged, so they are **1-based** (the `pos` column of the same file is 0-based).

| command | file | format |
|---|---|---|
| `lift` | `new.cram` + `.crai` | [CRAM](#lift-cram) |
| `lift` | `--stats` file | [JSON](#lift-stats-json) |
| `qc` | `P.vcf` | [VCF](#pvcf) |
| `qc` | `P.bed` | [BED, regions](#pbed) |
| `qc` | `P.findings.bed` | [BED, findings](#pfindingsbed) |
| `qc` | `P.json` | [JSON](#pjson) |
| `qc` | `P.depth.bedgraph` | [bedGraph](#pdepthbedgraph) |
| `qc` | `P.gc.bedgraph` | [bedGraph](#pgcbedgraph) |
| `qc` | `P.candidates.tsv` (with `--candidates`) | [TSV](#pcandidatestsv) |
| `qc` | `P.haplotypes.tsv` | [TSV](#phaplotypestsv) |
| `fish` | `P.fastq.gz` | [FASTQ](#pfastqgz) |
| `fish` | `P.fish.tsv` | [TSV](#pfishtsv) |
| `fish` | `P.fish.json` | [JSON](#pfishjson) |
| `fish` | `P.bait_ref.fa`, `P.bait.cram` (with `--bait`) | [FASTA, CRAM](#pbait_reffa-pbaitcram) |

With `--progress`, `lift`, `qc` and `fish` also write a [progress stream](#progress-stream---progress)
to stderr.

The debug subcommands write JSON with the same `format_version` first key (`cramhook diff`:
`P.summary.json`; `cramhook translate`: `P.stats.json`, the same keys as the lift stats); their
other files are not part of the public format.

## lift

### lift CRAM

- A coordinate-sorted CRAM with index, written at compression level 1. Files are ~1 % larger than
  at the default level, but encoding is faster.
- The header has an `@SQ` line for every contig of `new.fa`. It keeps the original `@PG` lines and
  adds a `cramhook` `@PG` line with the command line and the remap command, plus an `argv=sh1`
  marker in its `DS` (a unique id per lift: `cramhook`, `cramhook.1`, …).
- Every record carries `XL:A:T` (translated) or `XL:A:R` (remapped) unless `--no-xl-tag`; `s2` is
  removed from translated records. See [lift.md](lift.md#tags-on-translated-reads).

### lift stats JSON

Written with `--stats`: counts, deferrals by reason, parameters, per-phase timings. Top-level keys:

| key | contents |
|---|---|
| `format_version` | 1 |
| `contract_ok`, `contract_reason` | input-contract result; the reason is `null` when it passed |
| `preset` | `map-hifi` or `lr:hq` |
| `records_in` | `primary`, `supplementary`, `secondary`, `unmapped` records read |
| `translated_reads`, `translated_records` | reads / records translated |
| `deferred_reads` | reads sent to the remap |
| `remapped_reads` | primary records produced by the remap minimap2 |
| `out_primary_records` | primary records written (translated + remapped) |
| `deferred_by_reason`, `deferred_by_class` | deferral counts by reason (see [eligibility rules](lift.md#eligibility-rules)) and by class |
| `origin_crossing_translated`, `origin_crossing_deferred` | origin-crossing split reads |
| `qual_missing` | records without QUAL |
| `fallback` | `null`, or the tiny-input fallback: `reason` (`tiny_input`), `input_bytes`, `threshold_bytes`. A tiny-input fallback leaves `contract_ok` true and `contract_reason` null |
| `timings_s` | `diff`, `stream`, `minimap2_index`, `translate_pass`, `writer_drain`, `sort_index`, `minimap2_total`, `children_cpu`, `clip_screen`, `total` (seconds) |
| `params` | effective parameters: `min_mapq`, `xl_tag`, `max_clip`, `margin`, `max_de`, `t2_win`, `t2_step`, `t2_flank`, `t2_junction_flank`, `clip_k`, `clip_min_hits`, `clip_min_len`, `clip_flank`, `clip_subtract_old`, `unsafe_no_clip_screen`, `max_s2_ratio`, `max_minus_indel`, `max_zdrop`, `min_input_mb`, `threads` |

## qc

With `--out-prefix P`, qc writes the files below. Example lines from a PacBio run (console
summary, `P.vcf`, `P.bed`, `P.findings.bed`):

```
alerts: 0 red, 1 orange, 5 yellow, 0 gray (PASS); lowcov: 0 red, 0 orange, 0 yellow, 0 gray
contig_1  2180910  .  G  T  440.3  PASS  TIER=YELLOW;TYPE=SNP;AF=0.1762;DP=227;AD=187,40;ADF=104,30;ADR=83,10;DP4=104,83,30,10;HP=1;STR=0;CLASS=SNP;PV=9.239e-45;QV=2.450e-38;MQ0F=0.176;REASON=minority  GT:AD:DP  1:187,40:227
contig_2  8691  9978  clip_hotspot  0  orange  clip_both  0.756  1.000  left=64;right=64;ratio_left=0.365;ratio_right=0.365;s=38.61;f=0.3650
contig_2  8691  9978  finding_5  1  orange  misjoin  evidence=clips,clips_both;members=clip_hotspot:1;lowcov=0
```

### P.vcf

VCFv4.2, one sample column (`sample`), tiered site calls only. Header lines: `##fileformat`,
`##format_version`, `##source=cramhook qc <version>`, `##reference`, `##cramhook_preset`,
`##contig` per contig. INFO `TIER` (`RED`/`ORANGE`/`YELLOW`/`GRAY`); FILTER `PASS`/`LowCov`; INFO
`TYPE`, `AF`, `DP`, `AD`/`ADF`/`ADR`/`DP4`, `HP`, `STR`, `CLASS`, `PV`, `QV`, `MQ0F`, `REASON`,
plus the flags `MIXED`, `STRAND`, `ORIGIN`, `LOWCOV`, `LINKED`/`UNLINKED`, `GROUPED`, `HPGRAY`
(gray) and `CRISPR` (with `--crispr`); FORMAT `GT:AD:DP`. QUAL is the phred-scaled p-value (capped
at 999).

| field | Number, Type | meaning |
|---|---|---|
| FILTER `LowCov` | | DP below `--lowcov-abs`: too few reads to trust the call |
| INFO `TIER` | 1, String | what the reads show: `RED`, `ORANGE`, `YELLOW` or `GRAY` (homopolymer indel, below yellow) |
| INFO `TYPE` | 1, String | `SNP`, `INS` or `DEL` |
| INFO `AF` | A, Float | alt reads / DP |
| INFO `DP` | 1, Integer | reads with a base or a deletion at the position (unlike bcftools, deletions are counted); secondary, duplicate and QC-fail records excluded, supplementary included |
| INFO `AD` | R, Integer | ref and alt read counts (indel ref = DP minus all indel alleles at the anchor) |
| INFO `ADF`, `ADR` | R, Integer | ref and alt read counts, forward / reverse strand |
| INFO `DP4` | 4, Integer | ref-forward, ref-reverse, alt-forward, alt-reverse reads |
| INFO `HP` | 1, Integer | homopolymer run length (SNP: run containing the position; indel: run starting after the anchor) |
| INFO `STR` | 1, Integer | tandem-repeat span (period 2-4, ≥ 3 copies) at the position / after the anchor; 0 if none |
| INFO `CLASS` | 1, String | background class used for the test |
| INFO `PV` | 1, Float | beta-binomial upper-tail p-value against the class background |
| INFO `QV` | 1, Float | Benjamini-Hochberg q-value (m = assessable positions of the class) |
| INFO `MQ0F` | 1, Float | share of MAPQ 0 reads in DP |
| INFO `MIXED` | 0, Flag | alt fraction in the mixed band |
| INFO `HPGRAY` | 0, Flag | gray: a minority/strand-biased yellow INS/DEL in a homopolymer ≥ `--gray-hp-min` (`TIER=GRAY`; not counted in red/orange/yellow) |
| INFO `STRAND` | 0, Flag | strand-balance gate failed |
| INFO `LOWCOV` | 0, Flag | low coverage: DP below `--lowcov-abs` (abs) or below `--lowcov-rel` × the contig depth mode (rel); only abs sets `FILTER=LowCov` |
| INFO `LINKED` | 1, Integer | haplotype group id from the linkage second pass |
| INFO `UNLINKED` | 0, Flag | tested against ≥ 2 partner sites within the linkage window and linked to none |
| INFO `GROUPED` | 1, Integer | exact insertion alleles merged into this candidate (only when > 1; header line only when insertion grouping is on) |
| INFO `ORIGIN` | 0, Flag | indel before the first reference base (right-padded: anchor base is the last REF/ALT base) |
| INFO `CRISPR` | 0, Flag | the site's 0-based position (POS − 1; for indels the anchor) lies inside a CRISPR array detected by pilercr-rs; annotation only, TIER and FILTER are unchanged (header line only with `--crispr`) |
| INFO `REASON` | 1, String | tier rule that fired |
| FORMAT `GT` | 1, String | genotype (haploid, always `1`) |
| FORMAT `AD` | R, Integer | ref and alt read counts |
| FORMAT `DP` | 1, Integer | depth as INFO `DP` |

### P.bed

Tiered regions **plus** untiered annotations (`not_assessable`, `low_coverage`, and `crispr` with
`--crispr`, tier `none`). Header
line `#contig start end label score tier reason depth_ratio highmapq_share extra`. Columns contig,
start, end, label, score, **tier (col 6, 1-based)**, reason, depth_ratio, highmapq_share, extra
(`NA` for a missing ratio, `.` for an empty extra). Labels: `snp_cluster`, `large_ins`,
`large_del`, `clip_hotspot`, `depth_low`, `depth_half`, `depth_high`, `depth_step`,
`linked_haplotype`, `split_junction`, plus the untiered `not_assessable`, `low_coverage` and
`crispr`. A `crispr` row (with `--crispr`) has score 0, tier `none`, reason `.`, `.` for both
depth columns and `copies=<n>;repeat=<mean bp>;spacer=<mean bp>;consensus=<letters>` in extra; it
is **not** a finding member and is not counted in `alerts`. `depth_segment` (the paired interval
between two `depth_step` ends) is written to the BED but is **not** a finding member and is not
counted in `alerts`. A `clip_hotspot` whose clipped reads
continue on another contig gets `;clip_to=<contig>:<start>-<end>;clip_to_reads=<n>` (one pair per
side, `|`-joined).

### P.findings.bed

One finding per group of tiered calls within 1 kb. The two ends of a `split_junction` cluster and of
a `depth_segment` within `--finding-link-max` (default 100 kb; the shorter way round on a circular
contig) are **linked**: the two ends pool, together with every tiered call inside the span, into one
finding (tier = the highest member tier; `--finding-link-max 0` disables linking and gives the plain
grouping). Header line `#contig start end label score tier likely extra`. Columns contig, start,
end, label (`finding_<n>`), score, tier, `likely` (see [qc.md](qc.md#likely-hint)), extra. The
extra column holds `evidence=<list>;members=<label>:<count>,…;lowcov=<0|1>`. A finding that wraps
the origin of a circular contig is written as two rows with the same `finding_<n>` id. A finding
that holds one end of a cross-contig clip link also carries
`;link=finding_<k>@<contig>:<s>-<e>,reads=<n>,f=<f>,partner_depth=<d>` (k = the finding at the other
end, `.` when none; `d` = depth at that other end, no decimals for a whole count else one decimal;
several links are `|`-joined).

### P.json

`reads`, `assembly`, `alerts` (including a `gray` count for homopolymer indels, kept separate from
red/orange/yellow; `alerts.lowcov` has it too), `findings`, `warnings`, `contigs` (per contig
`circular` and `circular_by` — `verdict`, `header` or `flag` — say whether the contig wraps and which
source decided it, precedence `--linear` → FASTA header `circular=true`/`topology=circular` →
verdict, with a `circular`/`unclear` verdict wrapping and `linear` not; also `circularity`: verdict
circular/linear/unclear from reads spanning the end→start junction, the FASTA header flag, and the
junction as `clean`, `overlap_<n>` or `gap_<n>`), `linkage`, and `qc_qv` =
`max(0, −10·log10((trusted red+orange site calls + 1) / (confident bases + 1)))`, finite with zero
calls. `xlinks` counts the cross-contig clip links; `crispr` counts the arrays detected with
`--crispr` (always present, zeros without the flag). `qc_qv` is a summary of qc's own calls, not
comparable to Merqury or Inspector QV.

All top-level keys, in order:

| key | contents |
|---|---|
| `format_version` | 1 |
| `analysis_type` | `cramhook qc` |
| `version` | tool version |
| `preset` | `map-hifi` or `lr:hq` |
| `inputs` | `in`, `ref` paths |
| `indel_min_af`, `lowcov_abs`, `gc_window` | effective indel AF prefilter, absolute low-coverage depth and the GC track window in bp (`0` = no `P.gc.bedgraph`) |
| `reads` | `total`, `mapped`, `unmapped`, `unmapped_fraction`, `bases_total`, `bases_mapped`, `bases_unmapped`, `mean_aligned_span`, `records` (`primary`, `supplementary`, `secondary`, `duplicate`, `qcfail`, `unmapped`, `hard_clipped`) |
| `assembly` | `length_bp`, `contigs`, `assessable_bp`, `assessable_fraction`, `low_coverage_bp`, `low_coverage_fraction`, `confident_bp` |
| `crispr` | `arrays` and `bp` (total array span) detected in the reference with `--crispr`; zeros without the flag, so the schema is stable |
| `alerts` | `red`, `orange`, `yellow`, `gray`, `lowcov` (same tiers) and `by_type` (the same counts per site type and region label) |
| `findings` | `red`, `orange`, `yellow`, `by_likely`, and `lowcov` (the same for low-coverage findings) |
| `xlinks` | number of cross-contig clip links |
| `warnings` | list of `code`, `contigs`, `bp`, `fraction`, `message` |
| `qc_qv` | see above |
| `contigs` | per contig name, in assembly order: `length_bp`, `circular`, `circular_by`, `aligned_bases`, `depth_mean`, `depth_mode`, `depth_mode_excludes_zero`, `covered_fraction`, `copy_number`, `assessable_bp`, `low_coverage_bp`, `low_coverage_fraction`, `alerts`, `findings`, `warnings` (`code`, `reason`, `bp`, `fraction`, `message`), `edge_clip_records`, `depth_step_origin_dropped`, `depth_skipped`, `circularity` (`verdict`, `header_circular`, `spanning_reads`, `end_reads`, `start_reads`, `clipped_at_end`, `clipped_at_start`, `junction`, `junction_offset_median`, `junction_offset_consistent`) |
| `candidates` | pre-filtered candidates per type: `SNP`, `INS`, `DEL`, `total` |
| `background` | per class: `class`, the beta-binomial fit, `fit_sites`, `m`, `fallback`, `candidates` |
| `params` | effective qc parameters |
| `threads`, `chunks` | parallelism used |
| `timing_s` | `scan`, `qc`, `circularity`, `total` (seconds) |
| `linkage` | `regions`, `pairs_tested`, `groups`, `linked_sites`, `unlinked_sites`, `clusters_coincident`, `timing_s` |
| `circularity` | contigs per verdict: `circular`, `linear`, `unclear` |

### P.depth.bedgraph

Per-position **aligned-base** depth over the assembly, run-length encoded and 0-based (BED): reads
with a base at the position, deletions not counted — the same coverage IGV shows and the
`samtools depth -a` default. qc's own signals and INFO `DP` stay on the wider depth that also counts
reads whose only contribution is a deletion. Columns contig, start, end, depth; no header line.

### P.gc.bedgraph

GC fraction over the assembly in fixed, non-overlapping windows of `--gc-window` bp (default
100), computed from the `--ref` FASTA qc already loads. One row per window, 0-based half-open:
`contig`, start, end, `gc` = (G+C)/(A+C+G+T) over the window as a fraction with 4 decimals
(`0.5120`). Contigs are in the reference/header order of the depth track; the last window of a
contig is shorter (end = contig length), and a contig shorter than the window gets one short
window. Bases that are not A/C/G/T (`N`, ambiguity codes) are excluded from both numerator and
denominator, and a window with no A/C/G/T is `0.0000`. No header line. `--gc-window 0` writes no
file.

### P.candidates.tsv

Written with `--candidates`: every pre-filtered candidate with all statistics, one row each, with a
header line. Columns `contig`, `pos0`, `type`, `len`, `ref`, `alt`, `dp`, `dpf`, `dpr`,
`alt_reads`, `alt_fwd`, `alt_rev`, `ref_fwd`, `ref_rev`, `af`, `hp`, `str`, `mq0f`, `class`, `pv`,
`qv`, `wilson_lb`, `strand_ok`, `significant`, `tier`, `reason`, `flags` (`MIXED`, `HPGRAY`,
`STRAND`, `ORIGIN`, comma-joined, or `.`), `lowcov` (reason or `.`), plus `grouped` (exact alleles
merged, 1 when not grouped) only when insertion grouping is on (`--ins-group-min-len` > 0).

### P.haplotypes.tsv

One row per (haplotype group, read), with a header line. Columns `group` (group id), `contig`,
`span`, `read` (read name), `strand`, `alleles` (alleles carried), `sa` (SA target). This is what
`fish` selects from.

## fish

### P.fastq.gz

BGZF; each selected read once, in its original orientation (reverse-complemented / reverse-QUAL for
reverse records; missing QUAL becomes `!`), sorted by read name.

### P.fish.tsv

One row per (read, reason, finding), with a header line, sorted by read name, then reason, then
finding (numeric `finding_<n>`, `.` last). Columns `read`, `length`, `reason`, `finding`
(`finding_<n>` or `.`), `contig`, `pos`, `strand`, `record`, `detail`, `sa` (the read's `SA`
targets as `contig:pos:strand`, or `.`). Per read and reason, it keeps the largest event per
distinct finding (including `.` for events outside every finding; ties: smallest position, primary
over supplementary). A read with evidence in two findings, or inside one and outside all, gets one
row each. The `detail` values per reason are listed in [fish.md](fish.md#selection).

### P.fish.json

The selected reads and bp per reason and in total, the input record counts, the categories and
thresholds used, and timing. Top-level keys:

| key | contents |
|---|---|
| `format_version` | 1 |
| `reads` | `reads`, `bp` (distinct selected reads, summed length) and `by_reason` (per reason: `reads`, `bp`; a read can count under more than one) |
| `cap` | reads dropped by `--clip-cap`: `reads`, `bp` (summed read length); zeros when the cap is off. Dropped reads are not in `reads` or the outputs; see [fish.md](fish.md#capping-clip-reads) |
| `records` | input record counts: `total`, `primary`, `supplementary`, `secondary`, `duplicate`, `qcfail`, `unmapped`, `mapped_bp`, `unmapped_bp` |
| `params` | `scope`, `min_tier`, `finding_pad`, `unmapped`, `clips`, `indels`, `linked`, `bait`, `locus`, `locus_flank`, `clip_min_len`, `indel_min_len`, `indel_merge_gap`, `indel_piece_min`, `clip_cap` (0 = off), `clip_bin` |
| `bait` | `null` without `--bait`; see below |
| `notes` | categories turned off or skipped (e.g. `linked` without `--haplotypes`) |
| `timing_s` | `pass`, `total` (seconds) |

With `--bait`, the `bait` block holds the bait contigs (`contigs`: `name`, `length`, `reads`, `bp`,
`mapq0_reads` each), the overall selected `reads`/`bp` and `mapq0_reads` count, the lift stats
(`translated_reads`, `remapped_reads`, `contract_ok`, `contract_reason`, `fallback`) and
`timing_s`, plus `mode` (`any`/`primary`) and `supplementary_only_reads` (the reads selected only
through a supplementary bait record).

The TSV and the decompressed FASTQ are byte-identical across thread counts (only the JSON's timing
differs).

### P.bait_ref.fa, P.bait.cram

With `--bait`: `P.bait_ref.fa` (assembly + bait contigs, indexed) and `P.bait.cram` (+ index), the
lift of the input onto that reference (a [lift CRAM](#lift-cram)), kept for inspection in IGV.

## Progress stream (`--progress`)

With `--progress`, `lift`, `qc` and `fish` write machine-readable progress to **stderr**, one JSON
object per line, so a GUI can drive a progress bar. Without the flag nothing is written and every
output is unchanged; with it the output files are the same too (only the `lift` `@PG` `CL`, which
records the command line, shows the flag).

```
{"progress":1,"cmd":"qc","event":"plan","phases":["scan","circularity","call","linkage","output"],"elapsed_s":0.00}
{"progress":1,"cmd":"qc","phase":"scan","event":"start","elapsed_s":0.00}
{"progress":1,"cmd":"qc","phase":"scan","event":"update","frac":0.420,"records":36940,"records_total":87951,"total_exact":true,"elapsed_s":1.31}
{"progress":1,"cmd":"qc","phase":"scan","event":"end","elapsed_s":2.05}
…
{"progress":1,"cmd":"qc","event":"done","elapsed_s":2.40}
```

| key | meaning |
|---|---|
| `progress` | version of the stream, currently **1** (independent of `format_version`) |
| `cmd` | `lift`, `qc` or `fish` |
| `event` | `plan`, `start`, `update`, `end`, `done` or `error` |
| `phases` | (`plan`) the phase names in the order they run |
| `phase` | (`start`, `update`, `end`) the running phase |
| `frac` | (`update`) fraction of the phase done, in [0, 1], truncated to 3 decimals (so `1.000` only when the phase's work is complete) |
| `records` | (`update` of a record pass) alignment records read so far in this phase; never decreases |
| `records_total` | (with `records`) the records the phase will read: exact when `total_exact` is true, else the estimate `round(records / frac)` |
| `total_exact` | (with `records`) `true` when `records_total` is counted from the file, `false` for the estimate |
| `message` | (`error`) the error message |
| `elapsed_s` | seconds since the command started, 2 decimals, never decreasing |

Order: the `plan` line comes first; then every phase of the plan, in plan order, as one `start`,
zero or more `update` and one `end` line; the last line is `done`. On a failure the stream stops
with an `error` line (possibly in the middle of a phase) and the command exits non-zero with its
usual error message.

**Whole-percent rule:** an `update` line is written only when the whole percent `floor(100 · frac)`
exceeds the last one written for that phase, so a phase has at most 100 `update` lines (1 % …
100 %) and a jump (41 → 43 %) gives one line with the current `frac`. The number of lines therefore
depends only on the progress itself, not on timing. An update at 100 % is written when the work
measured by `frac` completes before the phase's own wrap-up; it is not forced, so the last
`update` of a phase may be below 100 % (the `end` line closes it).

**Records, not reads:** `records` counts alignment records (primary, supplementary, secondary,
unmapped), so a read with a supplementary alignment counts twice — the same unit as `samtools view
-c`. The exact total of an input comes from the CRAM container headers (CRAM 3.x, read without
decoding, in well under 1 ms) or the BAM index statistics, and is computed only with `--progress`;
when it is not available (another CRAM version, no BAM index or an index without its optional
statistics, an unreadable file) or a pass reads only part of the file, `records_total` is the
estimate. Should a pass read more records than an exact total, the rest of that phase switches to
the estimate (`total_exact:false`), so `records` never exceeds an exact `records_total`. Parallel chunked passes add the records of
each finished chunk, so `records` advances in chunk steps.

In a pass with an exact total (lift `translate_pass`, qc `scan`, fish `pass` and `payload`), `frac`
is the larger of the position-based fraction (table below) and `records / records_total`, so the
unplaced unmapped reads at the end of a sorted file still advance the bar; 100 % is reached exactly
when the last record is read (`records == records_total`), and below that `frac` stays under 1.
Phases with an estimated total keep the position-based `frac`.

| command | phases | `frac` |
|---|---|---|
| `lift` | `diff`, `translate_pass`, `writer_drain`, `sort_index` (the `timings_s` names of the stats) | `translate_pass`: reference position of the input record (length of the earlier contigs + position) / total reference length, from the header, or records / total when larger. Records: every input record, total exact. `writer_drain`: remapped primary records received from minimap2 / reads sent to the remap (one primary record per deferred read, so `records` / `records_total` are the same unit as `frac`; total exact = `deferred_reads` of the stats). `diff` (minimap2) and `sort_index` (samtools): no `frac`. With the tiny-input or contract fallback `diff` is skipped but still reported as a `start` + `end` pair |
| `qc` | `scan`, `circularity`, `call`, `linkage`, `output` | `scan`, `circularity`, `linkage`: reference bp of the finished chunks / total, with records. `scan` reads every record (total exact; records / total when larger); `circularity` and `linkage` read only the placed records of all / some contigs (total estimated). `call` (site calling, regions) and `output` (findings, writing the files): no `frac` |
| `fish` | `pass`, `bait_lift` (only with `--bait`), `payload`, `output` | `pass` (selection), `payload` (FASTQ payload of the selected reads): reference bp of the finished chunks / total, or records / total when larger (every record, total exact). `bait_lift` (the internal lift onto assembly + bait, reported as this one phase) and `output`: no `frac` |

Each line is written whole, under a lock, so it never interleaves with other stderr output of
cramhook. stderr also carries the usual human-readable messages (warnings, the final summary): **lines
not starting with `{` are not progress lines** and should be skipped. minimap2 and samtools run with
their stderr captured, not inherited. All stderr writes are best effort: if the reader closes the
pipe, the command still runs to completion and writes its outputs (with or without `--progress`).
