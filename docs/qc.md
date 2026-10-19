# qc: assembly sanity check

`cramhook qc` scans a read-to-assembly alignment and flags **sites** (VCF) and **regions** (BED)
where the reads disagree with the assembly, as `yellow`, `orange` or `red` (with `gray` below
yellow for homopolymer indels). It is a triage tool: it does not call genotypes, polish, or
classify structural variants beyond the `likely=` hint (a label that never changes the tier).

## Input

A coordinate-sorted, indexed BAM or CRAM mapped to the assembly with minimap2, using
`--preset lr:hq` (accurate ONT) or `map-hifi` (PacBio HiFi); the preset is read from the minimap2
`@PG` line unless `--preset` is given. Supplementary alignments count toward depth, base counts
and clips; secondary, duplicate and QC-fail records are skipped. Hard clips are accepted.

## Quick start

```bash
cramhook qc --in old.cram --ref old.fa --preset lr:hq --out-prefix q --candidates
```

| Flag | Default | Purpose |
|------|---------|---------|
| `--threads` | all cores | worker threads for the scan, circularity and linkage passes |
| `--candidates` | off | also write `q.candidates.tsv` |
| `--progress` | off | machine-readable progress on stderr, one JSON object per line (phases `scan`, `circularity`, `call`, `linkage`, `output`); see [output-format.md](output-format.md#progress-stream---progress) |

The thresholds (`--min-depth`, `--fdr`, `--clip-*`, `--depth-*`, `--link-*`, …) are listed by
`cramhook qc --help`.

## Outputs

With `--out-prefix q`, qc writes:

| file | contents |
|---|---|
| `q.vcf` | tiered site calls (`TIER`, FILTER `PASS`/`LowCov`) |
| `q.bed` | tiered regions plus untiered annotations (`not_assessable`, `low_coverage`) |
| `q.findings.bed` | one finding per group of tiered calls, with the `likely` hint |
| `q.json` | run summary: reads, assembly, alerts, findings, warnings, contigs, linkage, `qc_qv` |
| `q.depth.bedgraph` | per-position aligned-base depth |
| `q.candidates.tsv` (with `--candidates`) | every pre-filtered candidate with all statistics |
| `q.haplotypes.tsv` | one row per (haplotype group, read); `fish` selects from it |

Every column, INFO field and JSON key, with example lines, is in
[output-format.md](output-format.md#qc).

## Signals

Every signal reduces to "k reads disagree out of n": a **significance** score `s` (each signal has
a low bar — the candidate gate, i.e. it is reported at all — and a high bar) and a disagreement
**fraction** `f`. The tier is one rule for all signals:

| | `f <= 0.5` | `f > 0.5` |
|---|---|---|
| `s >= high bar` | **orange** | **red** |
| reported, `s < high bar` | yellow | yellow |

Red always needs `f > 0.5` (for sites, the Wilson lower bound > 0.5). A finding keeps the highest
member tier.

| signal | detects | tier rule (defaults) |
|---|---|---|
| site calls | a position where the reads disagree with the assembly base or carry an indel | candidates need only the prefilters (≥ 2 alt reads; AF ≥ 0.08 SNP, ≥ 0.15 `lr:hq` / 0.08 `map-hifi` indel). Red: strand-balanced and AF > 0.5 with Wilson lower bound > 0.50 and ≥ 5 alt reads (the old AF > 0.60 bound is gone). A strand-balanced AF > 0.5 that misses the red bar → yellow (`weak_majority`). Significant AF 0.40–0.50 → orange (`mixed`); significant AF < 0.40 → yellow (`minority`); significant AF ≥ 0.40 failing strand balance → yellow + `STRAND`; non-significant → no tier. A yellow INS/DEL with reason `minority` or `strand_bias` whose INFO `HP` is ≥ `--gray-hp-min` (default 5) is **gray** (`HPGRAY`), not yellow |
| SNP clusters | ≥ 5 qualifying SNP sites in a 500 bp window (each AF ≥ 0.10, ≥ 2 alt reads, p ≤ 1e-3, strand-balanced); sites below the absolute low-coverage depth and not-assessable positions are excluded | every window is at the high bar; `f = (r − 1) / r` with `r` = region mean depth ÷ contig depth mode: orange, red when `f > 0.5` (`r > 2`) |
| large indels | a merged indel ≥ 50 bp with ≥ 3 supporting reads and support ≥ 0.30 of depth | every call is at the high bar; `f = support / depth`: orange, red when `f > 0.5` |
| clip hotspots | clipped read ends (≥ 50 bp, ≥ 3 reads, chained within 50 bp) on the aligned side | reported when the ratio ≥ `--clip-yellow` (0.15); `s = −log10 p` (Poisson tail vs the local clip rate) ≥ `--clip-high-s` (100) → orange, red when `f > 0.5`, else yellow. `f =` cluster reads ÷ aligned-side depth, or `(left + right clips) / depth` for a two-sided `clip_both` hotspot |
| split-read junctions | a read whose consecutive alignments join two loci ≥ `--junction-min-dist` (1 kb) apart on one contig (read-coordinate gap ≤ `--junction-max-gap` 200 bp; on a circular contig either breakpoint within `--junction-edge-excl` 2 kb of an end is skipped, the origin wrap included). Clusters: both breakpoints within `--junction-tol` (200 bp; ONT breakpoints jitter), forward and reverse reads merged (a reverse-complement traversal of the same junction is one event, keeping both breakpoint sides); support = distinct reads; `orient` is `tandem` (end→start), `deletion` (start→end) or `inverted_l`/`inverted_r`; `fwd=`/`rev=` count the reads by traversal direction relative to the canonical `(pos,side)` order (`fwd` = visits the lower breakpoint first), not by read strand | two `split_junction` regions per cluster (one at each breakpoint, ± tol), fraction = support ÷ mean breakpoint depth. Reported with ≥ `--junction-min-reads` (5) reads and fraction ≥ `--junction-yellow` (0.02); every reported call is at the high bar: orange, red when `f > 0.5` |
| depth bands | 500 bp windows whose depth leaves the local ±50 kb baseline | every run is at the high bar; `f = 1 − r` for `depth_low` (≤ 0.25×, always red) and `depth_half` (0.35–0.65×, high-MAPQ share ≤ 0.25), `f = (r − 1) / r` for `depth_high` (≥ 1.75×, share ≥ 0.75): orange, red when `f > 0.5` (half when `r < 0.5`; high when `r > 2`) |
| depth steps | the per-position depth change `Δ = mean[i, i+w) − mean[i−w, i)` (`--step-window` 100 bp), scored over ±`--step-span` (25 kb) on a 1 kb grid as `z = (Δ − median) / (1.4826 · max(MAD, sqrt(local_depth/w)))` (the floor is inside the `1.4826 ×`); peaks are the local maxima of `abs(z)` suppressed ±`--step-nms` (500 bp) and refined at 1 bp, and circular contigs wrap. On a **wrapping** contig a peak whose 100 bp interval touches the origin (`start == 0` or `end == len`) is an aligner/representation artefact — a read crossing the arbitrarily chosen origin is split into a primary and a supplementary alignment, soft-clipped at both ends, so the read counts just inside each end differ — and is **dropped**; it is removed *after* suppression so it cannot unmask a neighbour it suppressed, and `depth_segment`s pair from the remaining steps (the per-contig count is `depth_step_origin_dropped` in the JSON). Plus an output-only `depth_segment`: the peaks are grouped into maximal same-sign runs, each represented by its strongest step, and a segment links an up run and a down run (gain up→down, dip/loss down→up). Runs alternate in sign, so there are two alternating matchings (start at run 0 or run 1; on a circular contig the last and first runs are neighbours, so a matching may wrap the origin; on a linear contig they are not). The matching that pairs more runs wins, then the smaller total forward segment length (`n` = contig length, so a wrapping arc is the long way); ties go to the run-0 matching | `f = abs(Δ) / max(depth left, depth right)`: red if `abs(z) ≥ --step-red` (20) and `f > 0.5`; orange if `abs(z) ≥ 20` and `f ≤ 0.5`; yellow if `abs(z) ≥ --step-yellow` (8). A `depth_segment` whose inside/flank ratio contradicts its direction (gain `≤ 1`, loss `≥ 1`) is dropped; its tier is the lower of the two end tiers |
| linked haplotypes | a second pass records, per read, the allele at every candidate SNP site; pairs of sites within 5 kb (≥ 3 alt/alt reads, concordance ≥ 0.8 both ways, phasing ≤ 0.1/0.2, pooled Fisher p ≤ 1e-6; per strand, Fisher p ≤ `--link-strand-p` 0.01 with ≥ 2 alt/alt reads) are **linked**; connected components are haplotype groups | `linked_haplotype` region over the group span; every group is at the high bar, tier by `f = max(group mean AF, (r − 1) / r)`: orange, red when `f > 0.5`. Member sites get INFO `LINKED`; tested minority sites linked to none get `UNLINKED` (or are dropped, `--unlinked-action`) |

## Tiers and FILTER

The tier says what the reads show; FILTER says how far to trust it. A site below `--lowcov-abs`
(20× `lr:hq`, 10× `map-hifi`) or below `--lowcov-rel` (0.25× the contig depth mode) is marked INFO
`LOWCOV`, and the absolute case sets `FILTER=LowCov`. The tier is never capped by depth: a real
error seen in few reads stays `TIER=RED` with `FILTER=LowCov`, and `bcftools view -f PASS` gives
the confident subset. A strand-balanced AF > 0.5 that misses the red bar (Wilson lower bound
≤ 0.5 or fewer than 5 alt reads) becomes yellow (`weak_majority`). A significant site with AF
≥ 0.40 that fails strand balance becomes yellow with a `STRAND` flag, even one that would
otherwise be red or orange; a non-significant strand-failing site gets no tier.
`--no-yellow-strand-gate` restores the flag-only rule for minority sites (AF < 0.40, which then
stay yellow with `STRAND`).

### Gray

**Gray** sits below yellow and covers homopolymer indels: a yellow call with reason `minority` or
`strand_bias` (a minority / strand-biased signal), of type INS or DEL, whose INFO `HP` (the indel's
forward homopolymer run after the anchor) is ≥ `--gray-hp-min` (default 5) is written as `TIER=GRAY`
with `HPGRAY` (the `REASON` is kept). It still gets its normal FILTER (`PASS`/`LowCov`) and stays in
the VCF, but it is **not** a finding member (never creates or joins a finding, never absorbed) and is
**not** counted in the red/orange/yellow `alerts` or tier counts (the JSON carries a separate `gray`
count per contig and top level; `qc_qv` uses red/orange only). SNPs, red/orange calls and
`weak_majority` (strand-balanced AF > 0.5 that missed the red bar) are never gray: a majority
homopolymer indel may mean the consensus itself is wrong, so it stays yellow. The point is to keep
low-quality homopolymer indels from noise visible without letting them dominate triage or the FP
score. `--no-gray-homopolymer` disables it (such calls stay yellow and do form findings); it is on
for both presets. Linkage and SNP clusters are unaffected: both work on SNP sites only, so a gray
INS/DEL can never be a linked group member or a cluster site.

## `likely` hint

The `likely` column of `findings.bed` comes from a fixed rule table — `collapsed_repeat`,
`false_duplication`, `extra_sequence`, `missing_sequence`, `misjoin`, `depth_change`, `base_error`,
`mixed`, `minority`, `subpopulation`, `unclear` — and never changes a tier. On the planted errors it
is right for 100 % of duplications and inversions, 92–98 % of deletions, random insertions, swaps,
SNPs and small indels, and 10 of 12 collapses ([qc-validation.md](qc-validation.md), §4). A
`linked_haplotype` finding is `subpopulation` when no depth or structural evidence explains it; a
high-depth / copy-number signal (a trusted `depth_high`, or the group's own depth ratio ≥ 1.75 with
high-MAPQ share ≥ 0.75) gives `collapsed_repeat`. `depth_half` and `depth_low` keep their own hints
(`false_duplication`, `extra_sequence`), and a `depth_step` that is the finding's only evidence is
`depth_change`.

## Linkage second pass

`qc` re-reads the file once more to separate a coherent second haplotype from independent
errors. It collects, per read (primary and supplementary records consolidated to one read), the
allele at every candidate SNP site, and counts site pairs genome-wide per contig. A pair of sites
within `--link-window` (5,000 bp) is linked when the pooled table has `a >= --link-min-reads` (3),
concordance `a/(a+b)` and `a/(a+c) >= --link-min-conc` (0.8), phasing `c/(c+d) <= 0.1` and
`b/(a+b) <= 0.2` and a one-sided Fisher exact p `<= --link-p` (1e-6), and the same holds within the
forward reads and within the reverse reads (each with `a >= 2` and its own Fisher p <=
`--link-strand-p`, 0.01: the strand check only guards against strand-specific errors, significance
comes from the pooled `--link-p` 1e-6).

Connected components of linked pairs are haplotype groups, but a merge is rejected when a
**decisive incompatibility** joins the two groups: a tested pair with
`a+b+c >= --link-conflict-min-reads` (6, enough reads carry the alt at X or Y) whose alt/alt share
`a/(a+b+c) < --link-conflict-max-share` (0.3, the two minority alleles mostly sit on different
reads). This stops a merge only when a pair within the window shows such a decisive conflict; when
that evidence is absent (an untested pair, or an underpowered but compatible one) the guard cannot
prevent two haplotypes from being bridged. An ambiguous middle site with no conflicting pair
therefore does not split a group. A tested minority site (AF < 0.40, ≥ 2 partners within the window
sharing ≥ 10 reads) in no group is either flagged (INFO `UNLINKED`, tier unchanged;
`--unlinked-action flag`) or dropped (default `drop`: tier none, reason `unlinked`, kept in
`candidates.tsv`). Group regions carry the reads' SA targets / clip sides as `elsewhere=` /
`clipped=` when ≥ 30 % of the group reads fall inside one 1 kb SA window or share a clip side.

## Performance and validation

About 1 s per bacterial genome (0.36–4.99 s over the 66 raw → polished cases, median 1.18 s,
16 threads; `s33`); the earlier 15-assembly comparison measured ≈ 12 s for Clair3 and ≈ 4 s for
Sniffles2. Planted errors: every large event kind × size is recalled at red+orange (100 %), and at
red for all but the two larger duplications and one collapse; SNPs/INS/DEL at 99–100 % red.
Raw → polished truth: clear-support events are recalled at red for 94–100 % (57 % on real ONT).
Most extra calls are yellow. Details and numbers: [qc-validation.md](qc-validation.md).

## Known limitations

- Most ONT false alarms are YELLOW; treat it as "look here", not as an error.
- By default (`--unlinked-action drop`) a minority site that is tested against nearby candidate
  sites and linked to none leaves the VCF (it stays in `candidates.tsv`). This removes about a
  quarter of the yellow false alarms on the raw → polished cases (1,142 of 4,346 in the run that
  measured it, before the gray tier) at the cost of 4 of their truth events, but it also drops a
  real minority variant that stands alone. Use `--unlinked-action flag` to keep such sites (yellow,
  INFO `UNLINKED`).
- The origin of circular contigs is not checked: the use case runs dnaapler, which rotates each
  assembly to start at dnaA, so breaks there are not covered.
- Reads of a large deletion can cross-map and cause off-target false alarms.
- Short, partial collapses in raw assemblies only reach orange.
- Low-coverage calls are reported separately (`FILTER=LowCov`) and never promoted or capped by depth.
- Not a variant caller: no SV classification, realignment or recalibration; `qc_qv` is not
  comparable to Merqury or Inspector QV.
