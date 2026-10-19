# `cramhook qc` — validation and quality controls

How `cramhook qc` is validated, and how far its output can be trusted. `qc` is a
triage tool: it scans a read-to-assembly alignment and flags suspicious **sites**
(VCF) and **regions** (BED) as `yellow`, `orange` or `red` (`gray` sits below
yellow and marks homopolymer indels). It does **not** classify structural
variants, call genotypes, realign or recalibrate (`docs/plans/Plan-v3.md`). A flag means
"the reads disagree with the assembly here"; the tier says how strongly, and
`FILTER` (`PASS`/`LowCov`) says how much the local depth supports it.

Sources of every number below: `docs/decisions.md` (from "qc split-read junction
signal" onward), `docs/plans/Plan-v3.md`, and the scratch reports
`scratch/qc/eval/s33/Report.md` (raw → polished, current baseline tag `s33`),
`scratch/qc/spike-eval/sp31/Report.md` (planted errors, current baseline tag
`sp31`), `scratch/qc-gray/Report.md` (gray tier, measured at `s28` vs `s26`),
`scratch/bait/eval/pbds2/summary.md` (phage `fish` demo, tag `pbds2`) and
`scratch/cmp-tools/Report.md` (Clair3 / Sniffles2 comparison, an **earlier** run,
commit `901e2cf5`; see §5–§6). Older tags (`s26`–`s32`, `sp24`–`sp30`) are kept
only where a comparison is historical. Numbers are reproducible with the
commands in §9. Each table names the tag it comes from.

## 1. Summary

**What qc is validated for.** That a flag lands at (within tolerance of) the
right position, and that the tier and `FILTER` are meaningful. The *kind* of
event is not scored (that is the `likely=` hint's job, §4).

**Headline results.**

- **Planted errors** (`scripts/qc_spike.py`, 6 samples × 2 replicates = 12 runs,
  6 unmodified controls; `sp31`): every large event kind × size is recalled at
  **red+orange** (100 %), and at **red** for all of deletions, random
  insertions, inversions and swaps; duplications reach red on 100 % / 50 % /
  67 % of the 1 kb / 5 kb / 10 kb events, and the one 1,447 bp collapse only
  reaches orange (the other collapses are red). Small events are recalled at red
  for 100 % of SNPs, 99 % of INS and 100 % of DEL (§4).
- **Raw → polished truth** (66 real cases: raw assembly vs its polished/consensus
  version; `s33`): truth events with clear read support (`clear` class) are
  recalled at red for 94.4 % (m6large ONT, 90 events), 100 % (val ONT, 85),
  100 % (val PacBio, 83), 100 % (real PacBio, 14) and 57.4 % (real ONT, 68)
  (§3). Events with only minority or no read support are largely not recalled by
  design.
- **Gray tier** (homopolymer indels; `qc-gray/Report.md`): on the 66 raw →
  polished cases the change moves 421 ONT and 103 PacBio VCF calls from yellow to
  gray and drops the yellow findings from 878 → 528 (ONT) and 142 → 51 (PacBio);
  red/orange recall and false alarms are unchanged. Gray calls stay in the VCF
  but are not findings and are not counted in the red/orange/yellow totals (§3).
- **Cost.** On the 66 raw → polished cases (16 threads) qc takes 0.36–4.99 s
  (median 1.18 s, RSS 244–919 MB; `s33`). The earlier comparison on 15 polished
  assemblies gave ≈1 s for qc against ≈12 s for Clair3 and ≈4 s for Sniffles2
  (§6).

**Main limits.** Truth that the reads do not support cannot be recalled; ONT
yellow false alarms are dominated by repeat/dosage ambiguity and strand-biased
errors; the origin of circular contigs is not checked; reads of a large planted
deletion cross-map and cause off-target false alarms; short, partial collapses in
raw assemblies reach only orange; low-coverage calls are reported separately and
are never promoted or capped by depth (§8).

## 2. Test suites

Run with `pixi run cargo test --release`. All suites pass (0 failed, 3 ignored;
429 passed).

| Suite | File | Tests | Covers |
|---|---|---|---|
| Unit | `src/lib.rs` | 298 | scan counters, statistics (incl. the Myers edit distance), tiers, `likely`, contract, diff/translate, VCF/BED/JSON writers |
| Synthetic integration | `tests/diff_synthetic.rs` | 11 | `diff` identity blocks / taint / novel (upstream of qc) |
| Synthetic integration | `tests/lift_synthetic.rs` | 24 | `lift` end-to-end on synthetic scenarios, chained lifted input |
| Synthetic integration | `tests/lift_junction.rs` | 10 | T2 junction taint (deletion/join, origin competitor, nested junction, seam) |
| Synthetic integration | `tests/lift_roundtrip.rs` | 1 (+1 ignored) | A→B→A and A→B→A→B over S1–S21 × hifi/ont; random SR1–SR5 behind `--ignored` |
| Synthetic integration | `tests/translate_synthetic.rs` | 7 (+1 ignored) | translate pass, split reads, contract; random SR1–SR5 behind `--ignored` |
| Synthetic integration | `tests/qc_synthetic.rs` | 26 | `qc` signals and outputs on in-process synthetic BAMs |
| Synthetic integration | `tests/qc_circularity.rs` | 17 | circularity verdicts, wrapping and cross-contig clip links |
| Synthetic integration | `tests/qc_linkage.rs` | 4 | the linkage second pass |
| Synthetic integration | `tests/fish.rs` | 31 | `fish` selection, outputs, `--bait`, `--bait-records` |

`tests/qc_synthetic.rs` builds reference + BAMs in process and asserts on the
VCF, BED, `candidates.tsv` and JSON. It contains the qc-specific regression
scenarios (small-input edge cases, indels, clusters, clip hotspots, depth bands,
depth steps, junctions, gray indels, linked findings) and the
**thread-determinism** checks: cases are run with `--threads 1` and
`--threads 8` (different chunk sizes) and must agree, but each test compares a
different subset. `qc_ins_grouping` compares the parsed VCF and `candidates.tsv`
plus the BED; `qc_clip_hotspot_junction_and_length_filter`,
`qc_depth_bands`, `qc_depth_step_region_and_thread_identity`,
`qc_bedgraph_aligned_depth_excludes_deletions` and
`qc_misjoin_clip_hotspot_and_finding` compare the BED (and `findings.bed`) plus
the normalised JSON; the main `qc_planted_events_and_thread_identity` compares
the VCF, BED, `candidates.tsv`, `depth.bedgraph` and `gc.bedgraph` byte for byte
plus the normalised JSON.

**Script self-tests** (`--self-test` on each script; see §9):

| Script | Self-test checks |
|---|---|
| `scripts/qc_truth.py` | truth construction on a real `old.fa` slice: SNP/INS/DEL/LINS/LDEL counts, homopolymer handling, contig-start, forward and reverse-complement strand |
| `scripts/qc_eval.py` | recall/FP/family/exact matching and PASS/low-coverage split on synthetic truth + calls |
| `scripts/qc_spike.py` | the generator without minimap2/reads: every structural kind present, structural coordinates, SNP/indel counts, and the SNP truth rows (REF = spiked base, ALT = original base, REF ≠ ALT) |
| `scripts/qc_spike_eval.py` | large-event matching, tolerance, tier/PASS split, signal names, FP window, inversion/duplication/swap intervals, finding parsing and `likely` acceptance |

All four self-tests pass and are repeatable.

## 3. Truth set 1 — raw → polished (66 cases)

**What it is.** For each case, truth = the differences between the raw assembly
and its polished/consensus version, in **raw (old) coordinates**, from
`minimap2 -c -x asm5 --cs=long old.fa new.fa` (primary alignments only),
parsed by `scripts/qc_truth.py` into SNP / INS / DEL (< 50 bp) / LINS / LDEL
(≥ 50 bp). Truth is read against `qc` calls by position; the event kind is not
scored.

**Sets and platforms (66 cases):**

| Set | Platform | Cases | Inputs |
|---|---|---|---|
| `m6large` | ONT | 18 | raw assembler output → polished |
| `real` | ONT | 12 | raw → consensus |
| `real` | PacBio | 12 | raw → consensus |
| `val` | ONT | 12 | raw → consensus |
| `val` | PacBio | 12 | raw → consensus |

**Caveats.** A raw → polished difference is not necessarily a majority error in
the raw assembly: it can be a mixed plasmid population, a repeat copy paired
differently by the two assemblies, or a low-depth polishing change. `qc_eval.py`
therefore classifies each truth event by the read support it actually sees
(`scripts/qc_eval.py`, `s33`):

- `clear` — max candidate AF > 0.6 (majority disagreement);
- `mixed` — AF 0.4–0.6;
- `minority` — 0 < AF < 0.4;
- `no_support` — no candidate at all (the reads do not show the difference);
- `repeat_unaligned` — event within ±10 bp of an `unaligned`/`multi` region
  (takes precedence over the AF classes);
- `unknown` — `qc.candidates.tsv` missing (none in these runs).

**Recall per read-support class and tier** (`s33`; `red` / `red+orange` / `any`
are cumulative; `any` = any tier):

| set | platform | class | n events | red | red+orange | any |
|---|---|---|---|---|---|---|
| m6large | ONT | repeat_unaligned | 137 | 85.4 % | 99.3 % | 99.3 % |
| m6large | ONT | clear | 90 | 94.4 % | 94.4 % | 100.0 % |
| m6large | ONT | mixed | 52 | 9.6 % | 15.4 % | 100.0 % |
| m6large | ONT | minority | 58 | 3.4 % | 5.2 % | 10.3 % |
| m6large | ONT | no_support | 3 | 0.0 % | 0.0 % | 0.0 % |
| real | ONT | repeat_unaligned | 75 | 40.0 % | 97.3 % | 97.3 % |
| real | ONT | clear | 68 | 57.4 % | 57.4 % | 98.5 % |
| real | ONT | mixed | 203 | 1.0 % | 12.8 % | 55.2 % |
| real | ONT | minority | 119 | 0.0 % | 0.8 % | 5.0 % |
| real | ONT | no_support | 42 | 0.0 % | 0.0 % | 4.8 % |
| real | PacBio | clear | 14 | 100.0 % | 100.0 % | 100.0 % |
| val | ONT | repeat_unaligned | 95 | 91.6 % | 98.9 % | 98.9 % |
| val | ONT | clear | 85 | 100.0 % | 100.0 % | 100.0 % |
| val | ONT | mixed | 8 | 25.0 % | 25.0 % | 62.5 % |
| val | ONT | minority | 3 | 0.0 % | 0.0 % | 0.0 % |
| val | ONT | no_support | 4 | 0.0 % | 0.0 % | 0.0 % |
| val | PacBio | repeat_unaligned | 127 | 63.0 % | 100.0 % | 100.0 % |
| val | PacBio | clear | 83 | 100.0 % | 100.0 % | 100.0 % |
| val | PacBio | mixed | 22 | 18.2 % | 100.0 % | 100.0 % |
| val | PacBio | minority | 6 | 83.3 % | 83.3 % | 100.0 % |
| val | PacBio | no_support | 2 | 0.0 % | 0.0 % | 0.0 % |

The whole-set (all classes) recall and false-alarm rate per set × platform
(`s33`; `pass` restricts to `FILTER=PASS` calls, so LowCov is excluded):

| set | platform | cases | truth SNP/INS/DEL/LINS/LDEL | recall red / red+orange / any | pass red / red+orange / any | FP R/O/Y | FP lowcov | FP/Mb |
|---|---|---|---|---|---|---|---|---|
| m6large | ONT | 18 | 247/59/21/11/2 | 61.5 % / 68.2 % / 83.5 % | 61.5 % / 68.2 % / 83.5 % | 35/431/370 | 0 | 17.94 |
| real | ONT | 12 | 253/154/92/5/3 | 14.0 % / 27.4 % / 51.3 % | 10.8 % / 20.5 % / 22.9 % | 34/498/1847 | 260 | 62.96 |
| real | PacBio | 12 | 0/13/0/1/0 | 100.0 % / 100.0 % / 100.0 % | 100.0 % / 100.0 % / 100.0 % | 6/11/22 | 0 | 1.07 |
| val | ONT | 12 | 130/37/19/8/1 | 89.2 % / 92.8 % / 94.4 % | 89.2 % / 92.8 % / 94.4 % | 27/357/120 | 0 | 16.57 |
| val | PacBio | 12 | 194/15/28/3/0 | 71.7 % / 98.8 % / 99.2 % | 71.7 % / 98.8 % / 99.2 % | 157/169/810 | 0 | 43.65 |

**False alarms per tier** (`s33`; calls with no truth event and no
`unaligned`/`multi` region within tolerance; `LowCov` is a subset of the tier
totals):

| set | platform | RED | ORANGE | YELLOW | total | LowCov | per Mb | calls |
|---|---|---|---|---|---|---|---|---|
| m6large | ONT | 35 | 431 | 370 | 836 | 0 | 17.94 | 1568 |
| real | ONT | 34 | 498 | 1847 | 2379 | 260 | 62.96 | 2851 |
| real | PacBio | 6 | 11 | 22 | 39 | 0 | 1.07 | 63 |
| val | ONT | 27 | 357 | 120 | 504 | 0 | 16.57 | 1013 |
| val | PacBio | 157 | 169 | 810 | 1136 | 0 | 43.65 | 1692 |

The totals above already exclude the depth-step origin artefacts (below). Two
output changes moved the false-alarm counts from the earlier `s28` numbers
without changing recall: `s32` drops a `depth_step` at the origin of a wrapping
contig (114 calls = 57 peaks removed; FP red/orange/yellow 277/1504/3199 →
259/1466/3169 over the 66 cases; the origin-drop decision in
`docs/decisions.md`), and `s33` adds the cross-contig clip links (`clip_to=`/
`link=` notes, JSON `xlinks`) with the VCF byte-identical and `.bed`/
`.findings.bed` identical once those note fields are removed: 56 links in 17
cases, no tier or `likely` change (the cross-contig clip links decision). The
finding note carries `reads=`, `f=` and `partner_depth=`; no BEDPE is written
(`P.links.bedpe` was tried and dropped).

The false alarms are dominated by YELLOW. On ONT, the ORANGE false alarms were
analysed separately (`scratch/worker-orangefp/Report.md`, s6b, N = 1349 PASS
ONT ORANGE site calls, an earlier run): 99.9 % (1348/1349) sit in a locus where
the assembly is collapsed or duplicated (local depth ≈ 1.2–2.4× the contig depth
mode) and the disagreeing reads carry a second, internally consistent haplotype
— genuine repeat/dosage ambiguity, not base-call errors (methylation-motif
artefacts 0, homopolymer-indel artefacts 0, unexplained 0). Recommendation in
that report: keep but relabel as repeat/ambiguous.

**Gray tier.** A minority or strand-biased yellow INS/DEL whose reference
homopolymer run (INFO `HP`) is ≥ `--gray-hp-min` (5) is written as `TIER=GRAY`
(INFO flag `HPGRAY`), stays in the VCF with its normal FILTER, is
never a finding member and is not counted in the red/orange/yellow totals or in
`qc_qv` (`qc-gray/Report.md`). On the 66 cases (`s28` vs `s26`):

| platform | RED | ORANGE | YELLOW | GRAY |
|---|---|---|---|---|
| ONT | 788 | 1118 | 3109 → 2688 | 0 → 421 |
| PacBio | 464 | 156 | 940 → 837 | 0 → 103 |

Findings per tier moved from 878 → 528 (ONT) and 142 → 51 (PacBio) yellow; red
and orange findings were unchanged. 29 truth events that `s26` recalled with a
yellow call are now missed (all INS/DEL, `hp` 5–13, mostly `minority`); red/orange
recall and false alarms are unchanged. On the planted runs (`sp26` vs `sp24`)
gray moves 78 ONT and 23 PacBio yellow calls out and planted recall is unchanged
(any-tier 456/456, 100 %). `--no-gray-homopolymer` disables it.

## 4. Truth set 2 — planted errors (tag `sp31`)

**Generator.** `scripts/qc_spike.py` edits the longest contig (the chromosome)
of a *finalized* assembly with a fixed event list, writes truth in spiked
coordinates, and remaps the **real reads** (no read simulation) to the spiked
assembly with minimap2. The unmodified assembly is mapped once per sample as a
false-positive control.

**Events.** Large (`events.tsv`): deletions 50 / 200 / 1000 / 5000 / 20000 bp;
random insertions 50 / 200 / 1000 / 5000 bp; tandem duplications 1000 / 5000 /
10000 bp; one collapsed natural repeat per sample (1002–1447 bp); inversions
5000 / 20000 bp; one segment swap (10000 bp). Small (`truth.tsv`): SNPs and
small insertions/deletions, spread over homopolymer lengths. Truth rows use the
`qc_truth.py` convention: `ref` = spiked (qc reference) base, `alt` = the base
the reads carry.

**Runs.** 6 samples × 2 replicates = 12 scored spiked runs, plus 6 unmodified
controls. Recall is cumulative; an event is recalled at tier T if a tiered call
is within ±max(50 bp, 1 % of size) of a breakpoint (or overlaps the interval for
`collapse`/`ins_random`).

**Large events, recall per kind × size.** Every row is 100 % at red+orange and
any tier; the table gives the **red** column (`sp31`):

| kind | sizes (bp) | events per size | red |
|---|---|---|---|
| del | 50, 200, 1000, 5000, 20000 | 12 each | 100 % each |
| ins_random | 50, 200, 1000, 5000 | 12 each | 100 % each |
| dup | 1000 / 5000 / 10000 | 12 each | 100 % / 50 % / 67 % |
| collapse | 1002, 1030, 1215, 1350, 1447 | 2, 2, 2, 4, 2 | 100 % each except 1447 bp: 0 % |
| inv | 5000, 20000 | 12 each | 100 % each |
| swap | 10000 | 12 | 100 % |

The `pass` columns equal the unrestricted columns (no LowCov calls in these
runs). `bridge_break` no longer exists: large events are caught by `clip_hotspot`,
`large_ins`/`large_del`, `depth_step` and `split_junction` (`sp31`, "which signal
fired"). The duplications and the 1,447 bp collapse that stay below red are
recalled at orange (red+orange 100 %).

**Small events, recall per type** (`sp31`; positional matching):

| type | n events | red | red+orange | any |
|---|---|---|---|---|
| SNP | 360 | 100 % | 100 % | 100 % |
| INS | 192 | 99 % | 99 % | 100 % |
| DEL | 168 | 100 % | 100 % | 100 % |

**False alarms in the controls** (unmodified assemblies; every tiered call is a
false alarm; `sp31`):

- 6 controls, 16 findings, 0.94 findings per Mb; by tier: red 2, orange 6,
  yellow 8; by `likely`: `depth_change` 5, `misjoin` 4, `minority` 2 and one each
  of `base_error`, `collapsed_repeat`, `extra_sequence`, `mixed`,
  `subpopulation`. (The 5 `depth_change` findings were `unclear` before the
  `likely=depth_change` rule; `sp26` had 22 control findings, `sp31` has 16
  because the depth-step origin drop removed 6.)
- Away from planted events (calls farther than 2 kb from any event): 323
  away-false-alarm calls over the 12 spiked replicates → **9.52 /Mb**; 110 of
  them (34 %) are recognised as **cross-mapped** (spike-specific evidence: near
  a deletion/collapse footprint and absent from the control), leaving
  **6.28 /Mb**. The spiked runs with the most away false alarms are
  `pacbio_FAM29030` and `ont_FAM22000` (`sp31`, false-positives table).
- Misses inside repeats: 97 planted events lie within 1 kb of a repeat copy;
  0 of them were missed.

**Accuracy of the `likely` hint** (a rule table that never changes a tier;
`correct` = an accepted answer, `unclear` = conflicting evidence, e.g. a long
insertion *and* a long deletion; `sp31`):

| kind | n | correct/recalled | unclear/recalled |
|---|---|---|---|
| del | 60 | 97 % | 0 % |
| ins_random | 48 | 96 % | 0 % |
| dup | 36 | 100 % | 0 % |
| collapse | 12 | 83 % | 17 % |
| inv | 24 | 100 % | 0 % |
| swap | 12 | 92 % | 0 % |
| SNP | 360 | 96 % | 0 % |
| INS | 192 | 97 % | 0 % |
| DEL | 168 | 98 % | 0 % |

The `likely` hint is looser than before because the linked-findings rule pools a
small error with the region that links it, and the small error then takes the
region's hint (`docs/decisions.md`, "linked findings"). Since `sp27`/`s29`, a
finding whose only evidence is a `depth_step` is labelled `depth_change` instead
of `unclear` (the evidence flag is `depth_step`; `depth_segment` is output only):
130 findings in `s29` and 33 in `sp27` move, with `.bed`/`.vcf` tiers, recall and
false alarms unchanged (`docs/decisions.md`, "`likely=depth_change`"). The hint
table above is unchanged between `sp26` and `sp31`.

**Findings per planted event** (`sp31`): one, except swaps (mean 3.33 findings —
4 breakpoints). Deletions, insertions, duplications, collapses, inversions and
all small events average 1.00.

## 5. Polished assemblies — red alarms (earlier comparison)

**This section and §6 describe an earlier run** (`scratch/cmp-tools/Report.md`,
qc at commit `901e2cf5`, 2026-10-13), *before* `bridge_break` was removed and
before the one-tier rule, linked findings and the gray tier. The qc call names
below are therefore historical; the read-level case reviews remain valid as
descriptions of the loci, and the tool costs are still indicative. The current
validation is §3–§4.

On the 15 polished assemblies of §6, qc raised **13 RED calls** (3 red sites +
10 red regions), reported as **7 RED locations**. The two higher-confidence
groups are positions where an independent caller agrees at high AF; the rest are
qc-only regions.

| case | contig:pos | qc call(s) | other tools | judgement |
|---|---|---|---|---|
| ont_FAM20862__flye_04 | 1:682877 | RED SNP, AF 0.67, DP 389 | Clair3 sensitive SNP AF 0.68 | real remaining issue |
| ont_FAM20862__flye_04 | 1:835158 | RED SNP, AF 0.64, DP 19987 | Sniffles2 DUP VAF 1.00 | qc site inside a flagged region |
| ont_FAM20862__flye_04 | 1:834501–837500 | RED `depth_high` + `clip_hotspot` + `bridge_break` | Sniffles2 DUP VAF 1.00 | qc region + Sniffles SV agree |
| ont_FAM22000__minipolish_04 | 4:4917 | RED SNP, AF 0.79, DP 3238 | Clair3 precise + sensitive SNP AF 0.79 | real remaining issue |
| ont_FAM20838__cons-to-polished | 4:13943–25171 | RED `depth_high` + `bridge_break` | none | real: collapsed / mis-joined multi-copy repeat (inspected) |
| pacbio_FAM29030-p1-1.1 | scf3:1001–2940 | RED `depth_low` + `bridge_break` | none | no assembly error; unusual coverage, flagging intended (inspected) |
| pacbio_FAM29046-p1-1.1 | scf4:17014–18682 | RED `bridge_break` + `clip_hotspot` | none | real: unsupported ~669 bp segment (inspected) |

The three qc-only regions were inspected read by read
(`scratch/worker-red3/Report.md`):

- **FAM20838 contig 4 (real problem).** A 35 kb multi-copy element (~60 copies
  per genome). A 2.9 kb unit appears twice in the assembly, 7.8 kb apart; 787
  reads carry deletions of 7,791–7,798 bp there and none the insertion. The raw
  metaMDBG assembly keeps this part as a separate 7.8 kb contig; the finalized
  assembly joined it into contig 4. A collapsed / mis-joined repeat.
- **FAM29046 scf4 (real problem, medium confidence).** 110 reads (76 % of the
  reads there) are clipped at 17,513 and continue at 18,182: a 669 bp segment
  that only 60 of 1,616 reads carry (0.17× depth), flanked by an exact 670 bp
  inverted repeat. Either a spurious insertion or a minor plasmid variant kept as
  the consensus.
- **FAM29030 scf3 (no assembly error; unusual coverage, flagging intended).** A
  circular 45 kb contig whose sequence is confirmed by reads everywhere. Read 5′
  starts concentrate near position 2,436 (about 340 reads in that 1 kb bin), and
  depth before it drops to 0.27×; the biological mechanism is unresolved.
  `bridge_break` counted reads that end within 500 bp as not bridging, and
  `depth_low` reacts to the dip, so both turned red although no read contradicts
  the assembly. This is intended: an abrupt coverage change like this is unusual
  and worth a look (user decision, `docs/decisions.md`).

Of the 7 red locations, 4 are supported by Clair3 or Sniffles2 (one of them, the
FAM20862 SNP at 835,158, lies inside the FAM20862 region), 2 more were judged
real by read-level inspection, and 1 is not an assembly error but an unusual
coverage pattern that qc is meant to flag.

## 6. Comparison with Clair3 and Sniffles2 (earlier comparison)

Source: `scratch/cmp-tools/Report.md` (33 inputs = 15 polished assemblies,
12 spiked replicates, 6 raw assemblies; each processed by all five runs; qc at
committed HEAD `901e2cf5`, Clair3 2.0.3, Sniffles2 2.8.1; 16 threads). **The
qc-side counts predate the bridge_break removal, the one-tier rule and the gray
tier, so they are historical; re-run the comparison to refresh them. The tool
costs are still indicative.**

**Scale / cost** (per genome):

| group | tool | n runs | wall median / max (s) | RSS median / max (MB) |
|---|---|---|---|---|
| polished | qc | 15 | 0.7 / 3.7 | 583 / 847 |
| polished | Sniffles2 | 15 | 3.9 / 25.9 | 69 / 102 |
| polished | Clair3 precise | 15 | 12.1 / 73.3 | 1059 / 1692 |
| polished | Clair3 sensitive | 15 | 11.6 / 75.0 | 1059 / 1691 |
| spike | qc | 12 | 0.9 / 1.3 | 580 / 698 |
| raw | qc | 6 | 1.8 / 2.1 | 716 / 827 |
| raw | Clair3 sensitive | 6 | 31.2 / 37.0 | 1198 / 1244 |

**Polished assemblies** (near-clean). Over the 15 assemblies qc flagged 528
sites and 50 regions (13 RED, 7 RED clusters); Clair3 `--haploid_precise` called
18 variants, `--haploid_sensitive` 1021; Sniffles2 17 SVs (mosaic 2). qc and
Clair3-sensitive agree on 120/528 qc sites (±5 bp); qc regions and Sniffles
overlap on 6/50. Most extra Clair3-sensitive calls are low-AF 0/1 sites below
qc's thresholds, and the Sniffles calls are dominated by whole-contig `DUP`
artefacts (VAF ≈ 1.0).

**Planted errors.** Small events are recalled at ≈100 % by both tools for every
type; for large events qc recalled 100 % of every kind × size via
`clip_hotspot` / `large_*` (with the historical `bridge_break`), and Sniffles2
also recalls ≈100 % of the same events (`--mosaic` is much worse, usually
0–25 %, but 100 % for 5 kb duplications and the 1,447 bp collapse).

**False alarms away from planted events** (12 spiked replicates; all rates
pooled over the 12 replicates, total calls ÷ 33.92 Mb): qc **8 red** calls →
**0.24 /Mb**, **70 red+orange** → **2.1 /Mb**, **316 all-tier** → **9.3 /Mb**;
Clair3 precise **0.12 /Mb** (4 calls), Clair3 sensitive **3.18 /Mb** (108),
Sniffles2 **0.44 /Mb** (15), Sniffles2 `--mosaic` **0.06 /Mb** (2). (The
current planted FP numbers are in §4.)

Comparing **tier for tier**: qc red-only (0.24 /Mb) is about **2×** Clair3
precise (0.12 /Mb); qc red+orange (2.1 /Mb) is **below** Clair3 sensitive
(3.18 /Mb); qc all tiers (9.3 /Mb) is ~**2.9×** Clair3 sensitive, ~**21×**
Sniffles2, ~**79×** Clair3 precise and ~**158×** Sniffles2 `--mosaic`. So qc buys
its 100 % large-event recall with a higher flag rate, and its conservative tiers
are comparable to the variant callers, while its YELLOW tier carries most of the
extra noise.

**Raw assemblies.** At the six ONT raw → polished cases with the most
orange+yellow false alarms (1155 calls), Clair3 precise calls essentially none
(0) and Clair3 sensitive calls 218/1155 ≈ 19 %, always at AF ≈ 0.4 (median
0.40). About 19 % of the qc false alarms are low-AF sites Clair3-sensitive also
sees (minority/mixed, often in dense clusters); about 81 % are qc-only
(systematic ONT error / repeat-region flags that a variant caller filters out).
qc's role is to flag, not to genotype.

## 7. Review process

Policy: an independent `reviewer` agent checks the diff before commit
(`scratch/manager/Handoff.md`; `docs/decisions.md`). Bugs found this way, with
the regression tests they led to:

- **`AF > 1` from indel normalisation** (commit `5feb6ff`). Grouping similar
  insertions in a ±5 bp window pooled alleles from different anchors, so support
  and depth no longer came from the same reads (`AF > 1`, overlapping REF/ALT).
  Fixed by grouping only same-anchor alleles (`--ins-group-min-len 0` restores
  exact-sequence counting). `qc_synthetic.rs::assert_af_le_one` is an assertion
  helper (AF ≤ 1 and alt ≤ DP for every candidate) reused by several synthetic
  tests; the regression test for this bug is `qc_ins_grouping`, which asserts the
  three noisy copies merge into one INS candidate with `alt_reads=9`/`dp=10` and
  `GROUPED=3` (TIER=RED), that threads 1 and 8 agree, and that
  `--ins-group-min-len 0` splits them into three low-support alleles with no
  `grouped` column or INFO field.
- **Wrong SNP truth ALT in the generator** (commit `2b598c6`). The spiking
  generator's SNP truth rows initially carried the wrong ALT; the generator
  self-test now asserts, for every SNP row, `REF = spiked base`,
  `ALT = original base (what the reads carry)` and `REF ≠ ALT`
  (`scripts/qc_spike.py::self_test`).
- **Myers edit-distance underflow** (commit `4b2bae1`). The banded bit-parallel
  edit distance (Myers with Ukkonen banding) underflowed an unsigned block
  cursor (`bl -= 1`) when pruning the only active block; debug builds panicked
  (minimal repro `[b'A'; 64]` vs `[b'C'; 64]`, band 0). Found by an independent
  reviewer with a randomised property test over 2,045 sequence pairs
  (`scratch/reviewer-myers/README.md`); fixed by returning `None` before the
  unsigned decrement, with a regression test in `src/qc/stats.rs`.
- **Indel support counted from the wrong reads** (step 1, commit `dc41563`).
  Left-normalising an indel could move it before the supporting read's first
  aligned base, or into a CIGAR `N` skip, so its support was divided by the depth
  of unrelated reads (`AF = 2.0`); two indels of one read could also collapse into
  one allele and count twice. Fixed by counting an indel only from reads with a
  depth-contributing segment at the normalised anchor, and once per read;
  regression tests `qc_normalised_before_read_start`, `qc_normalised_into_refskip`
  and `qc_multiple_indels_one_read_count_once`.
- **Panic on a large insertion at a contig end** (step 1, commit `dc41563`): the
  region interval fell outside the contig. Regression test
  `qc_terminal_large_insertion`.

Newer features carry their own regression tests: `qc_bedgraph_aligned_depth_excludes_deletions`
(the depth track counts aligned bases only, matching `samtools depth -a`),
`qc_depth_step_region_and_thread_identity`, `qc_misjoin_clip_hotspot_and_finding`
(the `misjoin` hint now comes from clip evidence) and the `qc_gray_*` tests
(gray calls stay out of findings, tiers and `qc_qv`).

## 8. Known limitations

- **ONT yellow noise.** Most ONT false alarms are YELLOW (370 in m6large ONT,
  1847 in real ONT, 120 in val ONT; `s33`). An earlier run measured that the
  strand gate removed ~85 % of ~29,000 ONT site false positives, and the ONT
  ORANGE false alarms are 99.9 % repeat/dosage ambiguity
  (`scratch/worker-orangefp/Report.md`). Treat YELLOW as "look here", not as an
  error.
- **Origin of circular contigs is not checked.** By decision, the use case runs
  `dnaapler`, which rotates each assembly to start at dnaA; an assembler may get
  the region around the original start wrong, but there is no reason to expect
  errors at the start of dnaA (`docs/decisions.md`, "qc scope decisions").
  Breaks at the origin of circular contigs are therefore not covered by the
  bridging-read check (and the `bridge_break` signal no longer exists). The
  depth-step signal is different: on a wrapping contig an origin `depth_step` is
  an aligner representation artefact and is dropped (`depth_step_origin_dropped`
  in the JSON), not reported.
- **Cross-mapping after large deletions.** Reads of a large planted deletion can
  cross-map and produce off-target false alarms. `qc_spike_eval.py` identifies
  these as a `cross-mapped` subset: 110 of the 323 away false alarms (34 %) on
  the 12 spiked replicates (`sp31`), concentrated in the `pacbio_FAM29030`
  runs.
- **Short, partial collapses in raw assemblies only reach orange.** In the raw
  assemblies, collapses are short (0.5–2.5 kb) and partial (~1.5× depth), and
  their SNP clusters already flag them orange, not red (`docs/decisions.md`,
  "qc depth bands"). This is not about collapses longer than the reads, and the
  planted short collapses are unaffected: they all reached red or orange (§4).
- **Duplications can stay below red.** The planted 5 kb and 10 kb tandem
  duplications reach red for only 50 % / 67 % of events (the rest are orange;
  red+orange is 100 %; §4).
- **Low coverage.** A site below `--lowcov-abs` (lr:hq 20×, map-hifi 10×) or
  `--lowcov-rel` (0.25× the contig depth mode) is marked `INFO LOWCOV`, and the
  absolute case sets `FILTER=LowCov`. The tier is never capped: a real error
  seen in few reads stays `TIER=RED` with `FILTER=LowCov`. On `real/ont` the
  all-call recall is 51.3 % but only 22.9 % over `FILTER=PASS` calls
  (`s33`), so about half of the recall comes from low-coverage calls; every other
  set is unaffected. Use `bcftools view -f PASS` for the confident subset.
- **Split junctions near the ends of circular contigs.** `split_junction` drops every
  junction with a breakpoint within `--junction-edge-excl` (2 kb) of an end of a wrapping
  contig, to keep the origin join out. On circular contigs below 200 kb this blanks a median
  18–20 % of the length, and on raw → polished it suppresses 21 clusters that pass the gate;
  3 of them are not origin-related (e.g. `utg000001c:944–7372`, 796 reads). An origin-only
  rule plus a read-continuation test was tried and not adopted (`docs/decisions.md`,
  "split_junction near contig ends: origin-only rule tried, not adopted").
- **Not a variant caller.** `qc` does not classify structural variants, realign
  or recalibrate; `qc_qv` in the JSON is a summary of qc's own calls and is not
  comparable with Merqury or Inspector QV.

## 9. How to reproduce

Build and run the tests:

```sh
pixi run cargo build --release
pixi run cargo test --release
```

Script self-tests (all pass):

```sh
pixi run python scripts/qc_truth.py      --self-test
pixi run python scripts/qc_eval.py       --self-test
pixi run python scripts/qc_spike.py      --self-test
pixi run python scripts/qc_spike_eval.py --self-test
```

Truth set 1 (raw → polished), build truth then run qc and evaluate:

```sh
pixi run python scripts/qc_truth.py --sets m6large,real,val
pixi run python scripts/qc_eval.py \
    --truth scratch/qc/truth --runs scratch/qc/runs \
    --sets m6large,real,val --run --tag s33 \
    --out scratch/qc/eval/s33 --threads 16
```

Truth set 2 (planted errors), generate the spiked assemblies then score
(`scripts/qc_spike.py` builds them; `qc_spike_eval.py --run --tag sp31` scores
them):

```sh
pixi run python scripts/qc_spike.py
pixi run python scripts/qc_spike_eval.py --run --tag sp31
```

`fish` selection and the phage demo:

```sh
pixi run python scripts/fish_eval.py --run --tag f3
pixi run python scripts/fish_phage_eval.py --tag pbds2 --dataset pb --configs abcd
```

Comparison with Clair3 and Sniffles2 (needs the `testing-clair3` and
`testing-sniffles2` pixi environments; the stored report is the earlier run of
§5–§6):

```sh
pixi run -e default python scratch/cmp-tools/analyze.py
```

Run `qc` directly (`--preset map-hifi` is for PacBio HiFi):

```sh
pixi run ./target/release/cramhook qc --in X.cram --ref X.fa --preset lr:hq \
    --out-prefix P --candidates --threads 16
```
