# Results

> Historical measurement log (lift, milestones M1-M11, 2026-10-07..11). Current qc/lift numbers
> are in `docs/qc-validation.md` and `docs/decisions.md`.

## M1 — test data, oracle, baseline noise, feasibility probe (2026-10-07)

Host: 16 cores, minimap2/samtools from the pixi env, `-t 16`, sort `-@ 4`.

### Baseline noise (oracle vs. mapping to `new.fa` + 50 kb random `dummy_noise` contig)

All 28 synthetic sample sets (14 scenarios × hifi/ont) and all 24 real cases: **100.0000 %** primary
agreement (contig, pos, strand, CIGAR, MAPQ), 100 % supplementary sets equal, depth identical at
100 % of positions, `invariants_ok=True`. Self-comparisons (oracle vs. oracle) are also 100 %.
Negative controls (`old.cram` vs. oracle): S2 hifi 20.7 % agreement, S11 ont 14.2 %. So compare.py
discriminates.

Consequence: a dummy contig does not perturb minimap2 on genomes this size. Criterion 7.2's "not
worse than baseline noise" therefore means **zero tolerance** under this definition. If M4 shows
small inherent disagreements (e.g. from minimizer frequency changes caused by real edits), a
stronger noise definition is needed. Do not loosen the criterion silently.

### Full-remap cost (the thing to beat)

Wall time of `minimap2 | samtools sort -O CRAM` (oracle.json): 2–27 s per real case.
Breakdown on the two largest cases:

| case | minimap2 → SAM | sort → CRAM | asm5 diff | `samtools view -c` old.cram |
|---|---|---|---|---|
| pacbio_FAM20843 cons-to-finalized (148k reads) | 22.4 s | 2.1 s | 0.85 s | 0.12 s |
| ont_FAM22000 cons-to-polished (112k reads) | 23.0 s | 3.7 s | 0.49 s | 0.59 s |

minimap2 is ~85–90 % of the cost. Fixed costs cramhook cannot avoid are diff + T2 mapping +
CRAM decode + sort/encode, ≈ 3–6 s.

### Feasibility probe (`scripts/probe.py`, optimistic upper bound)

Translatable = primary, mapped, no SA, MAPQ ≥ 60, clip ≤ 200, and window (±500 bp) inside one
identity block. Taint, `de` and tag rules are ignored. Full table: `data/real/probe_real.tsv`.

| group | fits (fraction of all reads) | read-level deferred | notes |
|---|---|---|---|
| PacBio, all 12 cases | 0.979–0.993 | 0.4–1.3 % | 2–28 blocks |
| ONT FAM20843 (4 cases) | 0.61–0.94 | 5–7 % | cons-to-polished: 94 blocks (polishing) → 0.80 |
| ONT FAM22000 (4 cases) | 0.57–0.63 | 32–36 % | **34 % of reads have SA** |
| ONT FAM20838 (4 cases) | 0.22–0.43 | 52–67 % | **45–62 % of reads have SA** |

Rough time estimate (fixed ≈ 4 s + deferred fraction × minimap2 time):

- PacBio cons-to-finalized: ≈ 4 + 0.01·22 + index ≈ 5–7 s vs 24 s → **~70–75 % saved**.
- ONT FAM22000: ≈ 4 + 0.4·23 ≈ 13–15 s vs 27 s → **~45–50 % saved**, borderline.
- ONT FAM20838: ≈ 60–75 % of reads remapped → **< 50 % saved**.

Caveats: per-read cost is not uniform (split/long reads likely cost more than average), and the
probe ignores taint, so real fractions will be lower.

Verdict: no abort. PacBio clearly clears 7.1. ONT is limited by **split reads (read-level), not by
the diff**. In ont_FAM20838 cons-to-polished, 5,689 primaries have SA. 99.7 % of them have all SA
records on the same contig, and 70 % have the primary within 2 kb of a contig end. So these are mostly
**reads crossing the origin of circular contigs** (four plasmids of 20–65 kb plus long ONT reads).
Translating origin-crossing split reads as a unit (all records inside identity blocks that map
consistently) would recover most of this. That is out of scope for v1 (section 9: "split reads are
always remapped"), so it needs an explicit decision.

### Other observations

- Synthetic S12: the 1.5 kb deletion block ends at 144999, not 145000. minimap2 places ambiguous
  indels one base off the planted coordinate. M2 tests must allow indel-placement ambiguity
  (equivalent positions) rather than exact planted coordinates.
- Synthetic reads crossing contig ends produce split alignments: ~12 % has_SA in synthetic samples,
  mostly on the plasmids.
- Real FASTA dirs were writable, so `.fai` files were written next to the source assemblies in
  `$DATA_ROOT/ont` and `$DATA_ROOT/pacbio`.

## M4 — full pipeline (2026-10-08)

Design: remap minimap2 (`-` stdin) and `samtools sort --write-index` spawned at t=0 (index build
overlaps the diff). Translate pass: htslib thread pool for decode; translated records as
**uncompressed BAM** into sort's stdin; deferred reads as FASTQ into minimap2's stdin; a reader
thread merges minimap2's SAM into the same writer. No compressed intermediates, no temp files.
Thread split for `--threads N`: decode N/2, minimap2 N, sort N/2, diff N/2.

I/O floors (`--threads 16`): decode `view -@16 -u` 0.45 s (PacBio FAM20843) / 0.88 s (ONT FAM22000);
sort → CRAM 1.7 / 2.6 s; index 0.2 / 0.4 s. `decode_md=0` drops NM too (needed for SA rewrite and the
output tags), so it is not used; with the default, the output tags decode identically to the oracle's.

Synthetic acceptance (`pixi run accept`, 28 runs): C1–C5 PASS everywhere (C2 literal: ≥ 99.9 % and
not worse than the baseline noise, i.e. 0 mismatches; C4 literal: equal to a full remap). C6
unnecessary diff-level deferrals: 0 % except S11 (0.8–1.4 %) and S12 (0.8–2.1 %).

Real cases (`--threads 16`, wall time of the `lift` process vs. oracle full remap):

| case | lift | full remap | speedup | translated | T agreement (ident / +equiv) | mismatches | depth identical (max diff) |
|---|---|---|---|---|---|---|---|
| pacbio FAM20843 cons-to-finalized | 3.2 s | 23.8 s | 7.4× | 99.1 % | 1.00000 / 1.00000 | 0 | 100 % (0) |
| pacbio FAM14197 flye_02-to-cons (both contigs reversed) | 2.6 s | 12.3 s | 4.7× | 97.8 % | 0.98282 / 0.99981 | 13 | 99.915 % (2) |
| ont FAM22000 cons-to-polished | 6.3 s | 26.9 s | 4.3× | 89.5 % | 1.00000 / 1.00000 | 0 | 100 % (0) |
| ont FAM20838 cons-to-polished | 3.3 s | 7.3 s | 2.2× | 55.9 % | 1.00000 / 1.00000 | 0 | 100 % (0) |

FAM14197's 13 mismatches are all on reversed contigs: same contig, pos, strand, MAPQ and AS, NM ±1
(co-optimal alignments that minimap2 finds differently on the reverse strand). Proper `hyperfine`
benchmarks of all 24 cases follow in M5.

## M5 — real data, tuning, validation (2026-10-08)

Final settings = defaults (margin 500, de 0.01/0.03, clip 200 + clip screen k=17/H=2, s2/s1 < 0.5,
'-' block indel < 50 bp, T2 win/step/flank 10k/2k/10k with circular windows). Benchmarks:
`real_eval.py --bench 3` (hyperfine, warmup 1; full remap = oracle command + `samtools index`;
lift = `cramhook lift --threads 16`). Summaries: `data/real/eval/final/summary.md` (tuning, 24 cases)
and `data/val/real/eval/final/summary.md` (hold-out, 24 cases from 6 other samples, seed 2).

**Correctness, all 48 cases:** invariants OK; translated-read supplementary agreement 100 %; depth
(deletions counted) identical at 100.0000 % of positions, max diff 0. Non-equivalent translated
primaries: 9 (tuning) / 3 (hold-out), all CIGAR-only on '-' blocks. The orientation study
(`scratch/orient/Report.md`) shows minimap2 itself produces the same differences when mapping to the
reverse-complemented assembly (accepted by the user as minimap2's direction effect).

**Speedup vs full remap (distribution):**

| set | group | n | min | median | max |
|---|---|---|---|---|---|
| tuning | PacBio | 12 | 3.56× | 3.95× | 5.94× |
| tuning | ONT | 12 | 0.49× | 1.18× | 2.52× |
| hold-out | PacBio | 12 | 1.13× | 3.26× | 4.62× |
| hold-out | ONT | 12 | 1.83× | 2.71× | 5.17× |
| both | cons-to-* | 12 | 0.96× | 4.0× | 5.94× |
| both | *-to-cons | 36 | 0.49× | 2.5× | 5.07× |

31 of 48 cases reach ≥ 2× (≥ 50 % compute saved, criterion 7.1). Translated fraction: tuning 83.7 %,
hold-out 82.2 %.

**Cases slow for biological/data reasons (not rule limitations):**
- ont FAM20838 (0.99–1.57×): a 7.8 kb circle between two IS copies is amplified (43 % of reads),
  and one polishing SNP inside the IS makes the copies more alike (T2), plus many chimeric/split reads.
- ont FAM20843 (0.49–1.17×): tiny dataset (4,914 reads, full remap 1.9 s); fixed costs dominate.
  ilesta_02 additionally has a large diff (9 → 1 contigs, 2.6 s of T2).
- ont FAM22000 *-to-cons (1.6–1.7×): primary assemblies differ structurally on the plasmids;
  ~50 % of reads are split/junction reads there. cons-to-polished is 2.5×.
- pacbio FAM29030 canu/pbipa (1.1–1.3×), FAM29046 minipolish (1.35×), ont FAM20866 raven (1.8×):
  the primary assemblies carry differing repeat copies, so T2 taints 17–48 k reads (correctly).

**Tuning:** see `docs/decisions.md` (M5). The clip limit was not raised (unsafe on real data even with
the screen). Candidates that were neutral on tuning (margin 250, de 0.02/0.05) were worse on hold-out
ONT, so the defaults stay. New rules found via validation: competitor (s2/s1) for MAPQ instability,
and the '-' block large-indel rule for structural direction effects.

## M6 step 1 — ONT baseline on the M6 set (2026-10-09)

Set: `data/m6/real` (44 cases = 11 ONT samples × cons→polished + 3 raw→cons; reads > 300 MB
downsampled to ~100 MB, seed 6; see `data/m6/README.md`). Binary = M5 final (no rule changes).
`real_eval.py --bench 3 --threads 16`; summary `data/m6/real/eval/baseline/summary.md`.

**Correctness:** 44/44 invariants OK, translated supplementary agreement 100 %, depth identical
100 % (max diff 0). One non-equivalent translated primary (FAM14197 minipolish_05, '-' block).

**Speed:** median 1.99× (0.43–3.98×); cons→polished median 2.45×, raw→cons median 1.86×.
Cases with a full remap > 3.5 s (n = 36): median 2.10×.

**Deferral by bases** (all 44 cases; 24.6 % of input bases deferred):
split_* family 10.8 % of input bases (split_junction 2.9, split_no_block 2.8, split_not_origin 2.7,
split_diff_contigs 0.9, other 1.5); no_block 4.5; taint_t2 3.4; de 2.0; competitor 1.5; mapq 0.8;
clip 0.8; minus_indel 0.5; unmapped 0.2.

**Fixed cost:** lift ≈ 1.1 s + 0.8 × (deferred fraction × full remap) + 0.6 × diff time
(least squares). At ~100 MB the fixed part (diff ≈ 0.5–0.7 s, sort/index ≈ 0.5 s) caps the
speedup at roughly 3× even with zero deferrals (median ceiling 3.3× for cases with a full
remap > 3.5 s). The cap rises with dataset size; deferrals dominate at full size.
Analysis: `scratch/manager/m6/fit.py`.

## M6 — full-size ONT raw → polished baseline (2026-10-09)

Set: `data/m6large/real` (18 cases: 6 samples with reads ≥ ~500 MB, 0.67–2.14 Gb, full size;
3 raw primaries → polished each). M5 binary, `real_eval.py --bench 3 --threads 16`. Run in two parts
(the first process died after 15 cases); merged summary `data/m6large/real/eval/merged/summary.md`.

**Speed:** median 2.43× (1.57–4.45×); full remap 18–54 s, lift 7–31 s.

**Correctness:** invariants OK and supplementary agreement 100 % in all 18 cases. **New:** in
FAM20856 flye_03/pflye_04 the depth is not identical everywhere (max diff 1, still rounds to
100.0000 %); T mismatches 4/3 include 1 pos (1 bp shift, 133S4M1I vs 134S4M1D) and 1 MAPQ (60 vs 59).
These reads are probably on '-' blocks (both cases have ~116 k reads translated through '-' blocks);
not analysed yet. M5 had 0 pos/MAPQ mismatches.

**Phases (s, e.g. FAM20853 myloasm_03: lift 7.2 s of 24.1 s):** diff 0.5–0.9 (1.9 raven);
translate pass 1.5–6 (15.5 raven); writer drain 2–9; sort + index 2.1–4.5. The translate pass is
2–3 s longer on cases with many '-' block reads (FAM20856 flye 4.7 vs myloasm 2.1 s; hypothesis:
the CIGAR reversal/left-align is CPU-bound in the single pass thread).

**Deferral by bases** (22.5 % of input bases deferred): taint_t2 5.9 %, no_block 4.8 %, split_* 5.8 %
(junction 1.8, not_origin 1.6, no_block 1.5), de 1.9 %, competitor 1.6 %, mapq 0.9 %,
minus_indel 0.8 %, clip 0.6 %.

## M6 — speed changes on the full-size set (2026-10-09)

Merged: the '-' block reverse-complement now works on the 4-bit encoding (`m6-translate`), and
CRAM output is written at `-l 1` (`m6-sort`). Records are identical to the baseline (body md5
checked on 5 cases); `cargo test --release` and `pixi run accept` pass.
`real_eval --bench 3` on `data/m6large/real` (tag `m6opt`; the full-remap timings are cached from
the baseline):

- **Speedup median 2.43× → 2.76×**, min 1.57 → 1.72×, max 4.45 → 5.13×. Total lift time over the
  18 cases 216.5 → 189.3 s (−12.6 %; median per-case ratio 0.88).
- Correctness unchanged (same residuals as the baseline: the FAM20856 '-' block direction effects
  plus one MAPQ near-tie at s2/s1 = 0.498; see `scratch/m6-mismatch/Report.md`).
- Speedups are against the full remap at samtools' default CRAM level; lift writes at level 1
  (+~1.3 % size).

## M7 — real-use-case validation (2026-10-10, measure only)

Build `74327e5` (= M6). `real_eval.py --bench 3 --threads 16`, tag `m7`; full remaps at samtools'
default CRAM level, lift at `-l 1`.

**Use-case set** `data/v2/real` (24 cases, full-size reads; `populate_m6.py --kind m7`, see
`data/v2/README.md`):

| kind | n | speedup median (range) | deferred bases, median per case |
|---|---|---|---|
| ONT cons → polished | 6 | 4.1× (2.8–5.8×) | 8.3 % |
| ONT polished → rotated (seed 7; chromosome reverse-complemented) | 6 | 3.5× (2.1–4.8×) | 9.9 % |
| ONT polished → PacBio-finalized (same ONT reads) | 6 | 3.3× (2.2–5.6×) | 13.0 % |
| PacBio cons → finalized (6 largest unused samples) | 6 | 4.8× (3.4–5.9×) | 3.1 % |
| all | 24 | 3.9× (2.1–5.9×) | 10.4 % of all input bases (M6 raw → polished: 22.5 %) |

Deferral by bases: `split_*` 44.8 % of deferred bases, `no_block` only 1.8 % of input.

**Correctness:** invariants OK in all 24 cases. Two cases fail gates for the first time:
`ont_FAM20862__polished-to-rotated` and `__polished-to-pbfinal` (same reads, chromosome on a `-`
block): translated supplementary agreement 99.999 % (2 reads differ), depth identical at 99.81 %
of positions (max diff 1; 4,858 positions, 4,856 from one read). `scratch/m7-mismatch/Report.md`:
all four affected reads reproduce exactly when minimap2 maps them to the reverse-complemented old
assembly, i.e. minimap2's own orientation effect (one read's z-drop split fires only on one
strand in a ~450 bp noisy stretch). All four passed every rule with margin. Not fixed (measure
only). Candidate signal: the old alignment's max lr:hq score drop ≥ 400 (minimap2's own z-drop),
2 reads per case (0.001 %). FAM20856 rotated/pbfinal show the known M6 residuals (max depth diff 1).

**M8 gate (small-edit patching), `scratch/m7-deferral/Report.md`:** recoverable `no_block` bases
(edits ≤ 5 bp) are 0.76 % of input (gate ≥ 10 %), best kind cons → polished 1.5 %; estimated lift
time saving 2.4 % (gate ≥ 20 %), best kind 5.4 %. Rotation and PacBio cases recover 0: their
`no_block` reads sit at structural breakpoints. **Gate not met.**

**Other sets, current build:**
- `data/m6large/real` (18 ONT raw → polished): median 2.74× (1.72–5.13×), m6opt 2.76×; total lift
  200.0 s vs 204.3 s. Same residuals as M6.
- M5 PacBio (24 cases of `data/real` + `data/val/real`): median 3.88× → **4.80×** (1.20–6.45×;
  tuning 3.95 → 4.92×, hold-out 3.26 → 3.60×). All gates pass.

## M7b — z-drop rule and supplementary gate (2026-10-10)

`real_eval --threads 16`, tag `m7b2` (`data/v2` with `--bench 3`; other sets correctness only,
single lift run). `cargo test --release` (120) and `pixi run accept` pass.

- **Correctness, 90 real cases** (`data/v2` 24, `data/m6large` 18, `data/real` 24, `data/val/real`
  24): all gates pass. FAM20862 rotated/pbfinal: supplementary 100 % (1 MAPQ-1 supplementary
  ignored), depth 99.9999 % (2 positions, the accepted 2 bp shift; before 99.81 %). Max depth diff
  is 1 everywhere.
- **Cost:** `zdrop` defers 16 / 11 / 0 / 11 reads (≤ 0.002 % of bases) per set. Total lift time
  on `data/v2` 178.4 s (m7) → 174.6 s (within noise).
- **`data/v2` speedup:** median 3.96× (2.10–6.04×); cons → polished 4.0×, polished → rotated 3.4×,
  polished → PacBio 3.3×, PacBio cons → finalized 4.9×.

## M9 — hardening (2026-10-10)

`cargo test` (release and debug: 125 tests), `pixi run accept`, and `pixi run real-smoke`
(5/5 PASS, 187 s) pass. Tiny-input fallback on ONT FAM20847 ilesta_02 (16.8 MB): 1.08 s vs
1.42 s full remap (normal lift 3.26 s); output identical to the full remap.

## M10 — remap overlap via minimap2 `-K` (2026-10-11, measure only)

`scratch/m10/Report.md`. 5 full-size cases, `--threads 16`, 3 runs, `-K` injected via
`CRAMHOOK_MINIMAP2`. Best `-K 20M`: median −4.6 % lift time (−3.8 … −10.5 %); with 12 remap
threads −7.5 %. **Gate (≥ 10 %) not met.** Overlap works (post-pass minimap2 wait 1–10 s →
~0.4 s), but minimap2's threads then preempt the single-threaded pass, which slows by about as
much. Records identical for every K (order-insensitive); with `-K ≤ 50M` the order of
equal-coordinate ties varies run to run (already at the default `-K` when deferred input
> 500 Mb). Peak RSS equal or lower.

Pass profile (FAM20866 rotated, default `-K`): the pass thread is 99 % CPU-bound (decode floor
1.1 s vs 4.1 s pass): `translate_one` 29 %, record rebuild (`Record::set` + `bam_set1`) 34 %,
memmove 10 %, `zdrop_statistic` 3 %.

## M11 — faster streaming pass (2026-10-11)

In-place record rewrite. Median of 3, `--threads 16`, base (M10) → new: translate pass −11 %
median (FAM20866 rotated 4.21 → 3.69 s, FAM22000 −4 %, PacBio FAM20855 −9 %, FAM14197 −13 %),
total lift −3.7 % median (−0.3 … −4.9 %). Records identical (order-insensitive md5) and stats
identical. `-K` re-test on the new pass (`scratch/m10/Report.md`, "M11 re-test"): best
`-K 20M` + 12 remap threads −7.3 % median, gate not met. `cargo test` (127, release + debug),
`pixi run accept`, `pixi run real-smoke` (5/5) pass.

## M10b — tag reliability (2026-10-11)

9,170,646 translated primaries + 79,571 supplementaries in 90 real cases vs the full remap
(`scratch/m10b/Report.md`). % equal on records with an identical alignment (81.7 %; the rest
differ only in CIGAR, all on `-` blocks): `NM` 100, `ms` 99.99997, `AS` 99.99999, `nn`/`tp` 100,
`cm` 99.9991, `s1` 99.9993, `rl` 99.9994, `de` 99.99999, `SA` 100, MAPQ 100 (2 of 9.17 M differ
overall). **`s2` 99.150 %** (99.857 % on `+` blocks; below 99.9 % in 63 of 90 cases; median
|Δ| 69, max 3533) → removed from translated records. `XL` is now on by default.
