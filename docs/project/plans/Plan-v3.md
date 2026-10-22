# cramlift — Plan v3: `cramlift qc` and `cramlift fish` (draft, 2026-10-11)

> The tool was renamed from cramlift to cramhook on 2026-10-19.

> Historical design document (Plan v3, 2026-10-11). Implemented; the decisions and measured
> results are in `docs/project/decisions.md` and `docs/validation/qc-validation.md`.

Status: draft, refined with the user after reading CRAQ, Inspector, Flagger, Clair3, NucFreq and
bcftools (reports: `scratch/manager/v3refs/*.md`; paper + code reading only, nothing run).
Numbers below are **starting values** to be measured, not final thresholds. Plan-v2 is complete.

## Goals

- **`cramlift qc`**: a fast assembly sanity check. From a read-to-assembly alignment
  (BAM/CRAM, mapping done upstream with minimap2), flag suspicious sites and regions as
  **yellow / orange / red**. A triage tool: it raises a flag that something is off at a position;
  it does **not** classify the structural variant, call variants, realign, or recalibrate.
- **`cramlift fish`**: extract reads that do not agree with the assembly, for follow-up work
  (taxonomy, plasmid rescue, aligning elsewhere).

## Architecture

Two subcommands on one shared scanning engine; `fish` consumes `qc`'s findings instead of
recomputing them.

```
src/scan.rs   one pass over the alignment → genome-wide per-position counters + per-read facts
src/qc.rs     counters → signals → tiers → VCF / BED / JSON / depth
src/fish.rs   read selection: unmapped, clips/indels, bait remap, regions from qc
```

- **`scan`** walks each record's CIGAR against the reference; no htslib pileup API (too slow),
  no local realignment (Clair3 disables it for long reads: lower precision, 3× runtime; BAQ is off
  in bcftools' long-read presets). Decoding and record parsing reuse the lift pass.
- **`qc` is strictly single-pass.** The scan is built so a later region pass over flagged regions
  (≈ 1 s of decode, shared with `fish`) can be added for per-read checks; added only if
  measurement shows clusters + depth miss or misflag collapses.
- **Memory:** plain genome-wide arrays (bacterial genomes, ~5 Mb + plasmids); no windowing.
- **Parallelism:** split by contig/region via the index, one reader per chunk, merge counters.
- **`fish` from `qc`:** `fish --regions qc.bed` extracts reads in a second pass. Default stays
  two commands.

### Input

- Any coordinate-sorted minimap2 BAM/CRAM; lift's strict contract does not apply to `qc`.
  `--preset` (HiFi or accurate ONT) selects thresholds.
- **Record filter:** skip secondary (0x100; the input may lack `--secondary=no`), duplicate
  (0x400) and QC-fail (0x200). **Count supplementary (0x800)** in depth, base counts and clips:
  each covers a different part of the read, and dropping them creates false depth dips and clip
  hotspots at every origin crossing and misjoin. Origin-crossing splits are recognised (lift's
  logic). MAPQ is taken per record (a supplementary can be MAPQ0 while its primary is unique).
- **Hard clips:** `qc` accepts them and counts `H` as clip length, and notes that the input is
  not usable for `fish`. **`fish` requires `-Y`** (soft-clipped supplementaries, full `SEQ` in
  every record), like lift, and stops with a clear error on a hard-clipped record.

### `scan` counters

Per position (all cheap):
- A/C/G/T **per strand** (needed for ONT: strand-biased errors at Dam GATC / Dcm CCWGG motifs);
- deletions covering the position (counted in depth), deletions starting after it, insertions
  after it (indels anchored on the preceding base, Clair3 layout); dominant indel length;
- left-clip and right-clip starts, separately;
- supplementary depth, MAPQ0 depth, high-MAPQ depth;
- large indels (≥ 50 bp) after merging nearby indels within one read (Inspector).

Per read: unmapped, clip lengths, clipped at both ends, longest indel, `SA` present.
From the assembly: homopolymer / tandem-repeat run length per position.

## `qc`: signals and tiers (starting values)

| Signal | Rule | Source |
| --- | --- | --- |
| Site disagreement | by alt fraction, see tier table below | Inspector, CRAQ, Clair3 |
| Minority allele | candidate pre-filter: ≥ 2 alt reads and SNP ≥ 8 %, indel ≥ 8–10 % (HiFi) / 15–20 % (ONT); then beta-binomial test + strand balance | Clair3, bcftools |
| Minority-SNP cluster | baseline ≥ 2 reads, ≥ 10 %, ~5 sites within 500 bp, as a true sliding window that wraps across the origin; cluster + ~2× depth → red; cluster at normal depth → orange | NucFreq, Flagger |
| Clip hotspot | clipped reads within ±25–50 bp ÷ depth on the aligned side, left/right separate, positions merged within 50 bp; > 0.75 → red, ~0.15–0.75 graded; share of reads clipped at both ends as supporting evidence | CRAQ |
| Depth | baseline = mode of depth per contig/plasmid; bands ≤ 0.25× / ~0.5× / 1× / ≥ 2×; ~0.5× only flagged where high-MAPQ share ≤ 0.25, ≥ 2× only where ≥ 0.75; ~1 kb windows (HiFi) + 5–10 kb view; minimum region length instead of an HMM; lower expected depth at linear contig ends | Flagger |
| Large indels | ≥ 50 bp, supporting reads ÷ local depth, with a minimum read count | Inspector |
| Not assessable | high MAPQ0 share, contig and coverage edges: reported as such, never as clean | all |

**Site tiers by alt fraction** (haploid: a majority disagreeing with the assembly means the
assembly base is wrong; 40–60 % is a mixed population or no clear consensus):

| Tier | Rule |
| --- | --- |
| yellow | significant minority allele, alt < 40 % |
| orange | alt 40–60 %, or a red candidate with too little evidence |
| red | alt > 60 % **and** lower confidence bound of the alt fraction > 50 %, minimum read count |

**Statistics:**
- Error background from the data, per preset and class (substitution; indel per homopolymer run
  length). bcftools' homopolymer indel error formula can seed the indel rates.
- Errors are correlated (same base, same strand), so a plain binomial is far too optimistic at
  100–200×. Use **both**: a beta-binomial (spread fitted once per class by moments, genome-wide)
  for significance, ranking and FDR, and **strand balance** as a gate for orange and above and
  a flag on yellow. If the beta-binomial turns out slow, fall back to the strand requirement.
- Window densities against the genome-wide rate; robust depth statistics.

**Severity:** an explicit, testable rule table. One weak signal → yellow; a strong or repeated
signal → orange; a majority-disagreeing call or several independent signals in one region → red.
Thresholds per preset.

**Later (after measurement):** `SA`-based evidence attached to clip hotspots, only as a "split
here" flag (no SV typing); region pass with per-read local-mismatch counts.

## `qc`: outputs

- **VCF** for site findings: AD, ADF, ADR, DP4; our own DP definition in the header (bcftools'
  DP excludes deletions); one haploid sample column (works with `bcftools isec`); tier in
  `FILTER`, reasons in an INFO field; `MIXED` for alt 40–60 %.
- **BED** (or GFF) for region findings: label, reason, depth ratio, high-MAPQ share.
- **JSON summary:** mapped/unmapped reads (count, bp, %), per-contig and overall depth, plasmid
  copy numbers from depth ratios, assessable fraction, counts per tier, and **`qc_qv`**
  (−10·log10 of flagged errors over assessable bases, with a pseudo-count). `qc_qv` is a summary
  of our own calls, not comparable with Merqury or Inspector QV; README defines it.
- **Depth track** for IGV (run-length bedGraph), replacing a separate `samtools depth`.

## `fish`: read selection (candidates)

Each selectable and tunable on its own:
- unmapped reads;
- reads with long clips or indels (inside the read or at the ends);
- reads in `qc`-flagged regions (consensus disagreement);
- **bait** (optional extra sequences, e.g. new contigs): probe the bait against the assembly
  (T2-style) to find similar regions; candidates = reads there + unmapped/clipped/split reads +
  reads whose clips share k-mers with the bait (clip screen); remap only the candidates to
  assembly + bait with the same minimap2 command; extract a read if its primary lands on the
  bait or scores better there;
- (later, region pass) reads carrying ≥ 2 of a cluster's minority alleles: separates paralog or
  contaminant reads without a remap (NucFreq report).

Output: FASTQ plus a table with the reason(s) per read.

## Validation (measure first)

Success means a flag at (within a tolerance of) the right position; the kind of event is not
scored.

- **Known consensus errors:** reads mapped to a raw assembly; the diff to the polished/finalized
  assembly lists positions where the raw assembly is wrong. Sets: `data/m6large` raw → polished
  (18 cases), `data/real`/`data/val` raw → cons.
- **Error-spiking generator (in scope):** edit a finalized assembly (delete one copy of a repeat,
  duplicate a block, swap segments, insert/delete 50 bp–50 kb) and remap the **real reads** with
  minimap2; no read simulation. Built after the SNP/indel signals are measured, before depth and
  clip thresholds are tuned. Recall by event size; misses inside repeats reported separately.
- **Real-world repeat case:** Hall et al. (eLife 2024): two near-identical 7.5 kb paralogous
  repeats where Clair3 missed 45 of 47 SNPs; tests cluster and low-MAPQ signals.
- **False positives:** finalized assemblies (expected nearly clean).
- **Depth stability:** down-sample; calls should be stable above ~30×.
- **Spot checks (concordance, not precision; same reads):** bcftools
  `mpileup -X pacbio-ccs|ont-sup -d 0 …| call -m --ploidy 1` for red only (it never reports
  minority alleles; `-d 0` because the default depth cap of 250 silently disables indel calls
  on high-copy plasmids); Clair3 `--haploid_sensitive` for yellow/orange.
- **fish:** pairs where the extra sequence is known (e.g. `data/v2` cases whose new assembly has
  extra contigs): does `fish` extract the reads a full remap places there, at what false-positive
  rate?
- Speed target: seconds per bacterial genome at full depth.

## References

CRAQ (clip-based misjoins), NucFreq (minor-allele clusters), Flagger (depth-based classes),
Inspector (structural and small errors), Clair3 (candidate gate, count layout), bcftools
(annotations, VCF conventions). Cloned read-only in `ref/`; ideas only, no copied code.
Paper supplements were not available; several paper and code defaults disagree (the reports
take the code as what runs).

## Order

1. `scan` + `qc` site signals (SNP/indel, minority, clusters) on the raw → polished truth sets.
2. Error-spiking generator; depth, clip and large-indel signals.
3. `fish` (unmapped/clip/indel selection, bait remap, `--regions`).

Each step as before: measure → stop and report → the user decides scope.

## Open questions (later)

- Exact thresholds and window sizes per preset (all numbers above are starting values).
- Tier combination rules for several weak signals in one region.
- Depth noise at 1 kb windows (≈ 10–15 % at 50–100×): measure on `data/m6large`.
- BED vs GFF.
- Bait candidate selection and score comparison details.
- `fish` output layout (one FASTQ per reason or a single file).
