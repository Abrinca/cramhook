# fish: extract reads that disagree with the assembly

`fish` is the extraction step. Two passes over the alignment: the first applies its own per-read
checks and selects reads into a reasons table, the second decodes the FASTQ payload of the selected
reads only. `fish` selects reads only; the reassembly is up to the caller.

```bash
cramhook fish --in old.cram --ref asm.fa --out-prefix P [--threads N] [--linear c1,c2] \
  [--findings P.findings.bed] [--haplotypes P.haplotypes.tsv] \
  [--scope findings|genome] [--min-tier yellow|orange|red] [--finding-pad 1000] \
  [--no-unmapped] [--no-clips] [--no-hairpins] [--keep-concatemers] [--no-tandems] \
  [--no-indels] [--no-linked] \
  [--locus [--locus-flank 5000]] [--bait bait.fa] [--bait-min-input-mb MB] \
  [--bait-records any|primary] \
  [--clip-min-len 50] [--circ-end-tol 100] [--indel-min-len 50] [--indel-merge-gap 50] [--indel-piece-min 10] \
  [--clip-cap N [--clip-bin 2000]] [--progress]
```

| Flag | Default | Purpose |
|------|---------|---------|
| `--linear c1,c2` | none | contigs treated as linear, overriding a FASTA `circular=`/`topology=` header; otherwise the FASTA header decides (`circular=true`/`topology=circular` or `circular=false`/`topology=linear`), and absent or conflicting headers default to circular |
| `--no-hairpins` | off | do not emit `hairpin` events (and do not fall back to `clip` for them, so a read selected only because of a hairpin is not selected; `hairpin_detected` in `P.fish.json` is still counted whenever clips are enabled and in scope) |
| `--keep-concatemers` | off | also select circular concatemers (reason `concatemer`); by default the assembly explains them, so they are counted and listed in `P.structure.tsv` but not selected; see [Concatemers and tandems](#concatemers-and-tandems) |
| `--no-tandems` | off | do not select `tandem` reads (no fallback to `clip`; `tandem_detected` still counts them) |
| `--circ-end-tol BP` | `100` | circular contigs: a clip stopping within `BP` of a contig end is an origin crossing (skipped) when `SA` has a same-contig, same-strand segment within `BP` of the opposite end (`0` = exact ends only; capped at half the contig length; the partner is read from `SA` and not checked to be a counted record or to cover the clipped bases) |
| `--clip-cap N` | `0` (off) | keep at most `N` clip-selected reads per bin of `--clip-bin` reference bp; see [Capping clip reads](#capping-clip-reads) |
| `--clip-bin BP` | `2000` | bin size for `--clip-cap`; `0` is an error (`--clip-bin must be positive`) |
| `--progress` | off | machine-readable progress on stderr, one JSON object per line (phases `pass`, `bait_lift` with `--bait`, `payload`, `output`); see [output-format.md](../reference/output-formats.md#progress-stream-progress) |

## Selection

Counted records skip secondary (0x100), duplicate (0x400) and QC-fail (0x200); primary and
supplementary are both checked. **`fish` requires a `minimap2 -Y` alignment** (supplementary
records must carry SEQ): any hard clip in a counted record is an error. A read is selected when any
of its counted records qualifies for any enabled category:

| reason | a read is selected when | detail |
|---|---|---|
| `unmapped` | it has an unmapped record (0x4), placed or unplaced | `.` |
| `clip` | a terminal soft clip ≥ `--clip-min-len` (not converted to `hairpin`); on a **circular** contig a clip at an exact contig end — or within `--circ-end-tol` bp with a same-contig, same-strand `SA` segment within `--circ-end-tol` bp of the opposite end — is an origin crossing (skipped), on linear contigs (`--linear` or a `circular=false`/`topology=linear` FASTA header; absent/conflicting headers default to circular) it counts | `side=left;len=N` / `side=right;len=N` |
| `hairpin` | a qualifying terminal clip (≥ `--clip-min-len`, in scope, not an origin crossing) is paired with an opposite-strand `SA` entry on the same contig whose reference overlap is ≥ 50 % of the shorter span, whose original-orientation query interval is adjacent to the record's (\|gap or overlap\| ≤ max(200, 0.02 × read length)) with each arm ≥ 20 % of the read length and ≥ `--clip-min-len`, and whose arm covers the clipped read region leaving residual < `--clip-min-len`. Any pair of counted records qualifies (primary/supplementary or supplementary/supplementary); the partner is taken from the record's `SA` tag (not checked to be a counted record) | `arm1=q<s>-<e>,r<s>-<e>,<strand>;arm2=q<s>-<e>,r<s>-<e>,<strand>` (`arm1` earlier in the read; `q` = 0-based half-open original-orientation query coords, `r` = 0-based half-open reference coords) |
| `concatemer` | with `--keep-concatemers`: a clip is covered by the read's circular-concatemer chain (see below) | `period=P;units=U` |
| `tandem` | a clip is covered by the read's tandem chain, unless `--no-tandems` (see below) | `period=P;units=U` |
| `large_ins` / `large_del` | a merged large indel (`--indel-merge-gap`, `--indel-piece-min`) ≥ `--indel-min-len` | `len=N;at=P` |
| `linked` | the read is listed in `--haplotypes` for a group whose span overlaps a finding of tier ≥ `--min-tier` (findings scope) or any group (genome scope) | `group=N` |
| `locus` | any counted record overlaps a finding (tier ≥ `--min-tier`) ± `--locus-flank` | `.` |
| `bait` | with `--bait`: a full remap to assembly + bait places one of the read's counted records on a bait contig (see below) | `contig=<bait>;mapq=<MAPQ>;s1=<s1>;s2=<s2>;record=primary\|supplementary` |

Scope `findings` (the default when `--findings` is given) counts `clip`/`hairpin`/`concatemer`/`tandem`/`large_*` only when
the event position (clip position, indel start) is within a finding of tier ≥ `--min-tier` ±
`--finding-pad`; scope `genome` counts them anywhere and `unmapped` is genome-wide in both.
When a hairpin partner leaves a residual ≥ `--clip-min-len` in the clipped region, or on the read's
other clip side, that clip stays `clip` (so a hairpin read with extra unaligned sequence keeps its
`clip` reason). With `--no-hairpins`, hairpin-converted clips are not emitted as either `hairpin` or
`clip` (though `hairpin_detected` in `P.fish.json` still counts in-scope, clip-enabled converted
events; both `hairpin_detected` and `origin_crossing_skipped` are 0 when `--no-clips` is set, and
`origin_crossing_skipped` is counted before `--scope` filtering). `--locus` and `--scope findings`
require `--findings`. There is no positive `--linked` flag: linked
selection is on by default and requires `--haplotypes` (if no file is given it is silently skipped
and recorded in the JSON; `--no-linked` turns it off). Reads that merely end at a position, or agree
with the assembly, are not selected; site-only and depth-only findings contribute only through
`--locus`.

## Concatemers and tandems

A read whose own segments (the record and its `SA` entries) show the same stretch of reference
repeated, back to back, is not a chimera: it is a read of a circular molecule, or of a tandem
repeat. Segments of >= 500 bp are ordered by their original-orientation query interval, and two
consecutive segments on the same contig and strand are *back to back* when their query gap or
overlap is <= max(200 bp, 2 % of the read length). In the read direction, let `d` be the next
segment's reference start minus the previous segment's reference end.

- **`concatemer`** (circular contigs only): `d` is 0 modulo the contig length `L` (within the
  same tolerance, capped at `L/4`): the read goes round the circle, e.g. `[x, L)` then `[0, y)`, or
  a full circle `[0, L)` followed by `[0, L)`. `period` = `L`, `units` = the summed aligned query
  length of the chain / `L`; called at **units >= 1.2**. Below that it is a plain origin crossing
  (skipped, as before). It is not selected by default because the assembly explains it; use
  `--keep-concatemers` to select it.
- **`tandem`**: `d <= -1000` (a back jump `J` >= 1 kb; on a circular contig `d` is first reduced
  modulo `L` into (-L/2, L/2], so an origin crossing is never a tandem and the jump must be
  `J < L/2` after the reduction). `period` = the median `J` of the chain,
  `units` = the chain's aligned query length / `period`; called at **units >= 1.2**. The assembly
  does not explain it (a collapsed tandem repeat, an excised circle, a duplicated prophage), so it
  is **selected by default**; `--no-tandems` drops it. Segments that are far apart in the read
  (an IS pair or an rRNA operon: same strand, back jump, unaligned query between them) are not
  back to back and give no call.
- **Precedence**: `hairpin` (opposite strands) > `concatemer` > `tandem` > the origin-crossing
  skip; a read has one structure, listed in `P.structure.tsv`. The label is a property of the
  read: hairpin geometry decides it whatever `--no-clips`, `--scope` and `--no-hairpins` say.
  Those flags only gate *selection*, so `P.structure.tsv` also lists hairpin reads that are out of
  scope or run under `--no-clips` / `--no-hairpins` with `selected = no`, whereas `hairpin_detected`
  keeps its 0.1.2 meaning (clips enabled, hairpin in scope) and can therefore be smaller than the
  number of `hairpin` rows. A hairpin read has no `concatemer` / `tandem` events: the clips they
  replaced revert to what they would have been without structure detection (`clip` events if
  clips are on and the clip is in scope, origin-crossing skips).
- **Only the best call converts clips.** A read can hold several chains (e.g. on two contigs);
  the best call (concatemer before tandem, then most units) is the read's structure and the only
  one that converts clips. Clips of a second chain stay `clip`.
- **Clips.** A clip becomes a `concatemer` / `tandem` event when the read's chain (which must
  contain the clipped record) covers the clipped read region, leaving a residual < `--clip-min-len`;
  this is checked after the hairpin rule and before the origin-crossing skip (a circular
  concatemer's clips all sit at the contig ends). A clip with a real residual stays `clip` (or an
  origin crossing), so such a read keeps its `clip` row. Like hairpins, a suppressed `concatemer` /
  `tandem` (default / `--no-tandems`) does not fall back to `clip`. Both are selected through
  clips, so `--no-clips` selects neither (detection, the counts and `P.structure.tsv` stay).
- **`--clip-cap`.** `concatemer` and `tandem` rows are never thinned and never consume cap slots.
- **Per-contig statistics.** `P.fish.json` `contigs[]` reports `multimer` (circular contigs up to
  50 kb: how many reads carry 1, 2, 3, 4+ units of the contig), `concatemer_reads` (reads with a
  circular-concatemer call on the contig) and `origin_crossing_reads` (reads with an
  origin-crossing clip on the contig that is not part of a concatemer / tandem chain: evidence that
  it is circular). The last two do not depend on `--no-clips`, `--scope` or the selection flags; see
  [output-format.md](../reference/output-formats.md#pfishjson).

Hidden thresholds (recorded in `params`): `--concatemer-min-segment` (500), `--concatemer-min-units`
(1.2), `--tandem-min-period` (1000), `--chain-max-gap` (200), `--chain-max-gap-frac` (0.02),
`--concatemer-max-contig` (50000).

## Capping clip reads

On a high-depth sample a few hotspots (repeats, rDNA, IS elements, misjoins) can contribute hundreds
to thousands of clipped reads per bin and flood the reassembly input: FAM20862 (~230x ONT) has
9,710 clip reads in a single 2 kb bin. **`--clip-cap N`** keeps at most `N` clip-selected reads per
`--clip-bin` bp bin of the reference (default bin 2000 bp). It is off by default and then the outputs
are unchanged; most bins are below the cap and unaffected.

This is a downsampling knob, not a filter: where a bin exceeds the cap, the kept clip coverage drops
to about `N` reads. Measured on the plant-error eval (genome scope, 72 deletions/collapses), cap 10
cut ONT selected reads by 36 % (bp 29 %) and cap 20 by 31 % (bp 24 %), and every planted event still
had at least 10 (cap 10) or 20 (cap 20) selected reads on it (median 132 / 140). Aggressive cuts are
therefore safe for reassembly but are not "free": raise `--clip-cap` (20-50) when the downstream
step wants more depth, and use `--clip-min-len` (instead) when the noise is uniform rather than
hotspot-driven.

- **Eligibility.** Only a read whose *sole* reason is `clip` can be dropped. A read selected by
  `hairpin`, `concatemer`, `tandem`, `unmapped`, `large_ins`, `large_del`, `linked`, `locus` or `bait` is never dropped, even
  when it also has a `clip` event (and `hairpin` rows do not consume cap slots). Every clip-only read
  is eligible, whatever the length of its clip.
- **Binning.** A read is placed by its best clip (the largest, then the smallest position): bin =
  (contig, position / `--clip-bin`).
- **Choice.** In a bin with more than `N` eligible reads, the `N` reads with the largest best-clip
  length are kept (ties by read name, ascending). A long clip spans furthest into sequence the
  assembly does not hold, so this favours exactly the reads worth re-examining, and a hotspot's
  short clips cannot crowd out its long ones.
- **Determinism.** The ranking uses only the clip lengths and the read names, so the kept set is
  identical for every `--threads` and input order.
- **Effect.** Dropped reads are not selected at all: they are absent from `P.fastq.gz`,
  `P.fish.tsv` and the `reads` / `by_reason` counts, and the payload pass skips them.
  `P.fish.json` reports the totals in `cap` (`reads`, `bp`; zeros when off) and `params` records
  `clip_cap` and `clip_bin`.

## `--bait`

**`--bait`** answers "which reads would a full remap to assembly + bait place on the bait?" for extra
sequences such as a phage or plasmid contig. It reuses the `lift` session: it writes `P.bait_ref.fa` =
the `--ref` assembly plus the bait contigs (bait contig names must not clash with assembly names — a
clash is an error, names are never renamed), then runs the same lift as `cramhook lift` (`--threads`
passed through) to write `P.bait.cram` (+ index). `lift` translates reads in regions that are
identical and unique in both references and remaps the rest with the same minimap2 command, so the
placement equals a full remap onto assembly + bait (T2 projection marks assembly regions similar to
the bait; clipped/unmapped/split reads are deferred). Every read with a counted record on a bait
contig is selected with reason `bait`; the detail names the best record (primary first, then highest
MAPQ) and its `MAPQ`, `s1`/`s2` and record kind. A read whose original primary record is duplicate
(0x400) or QC-fail (0x200) is **not** selected: `lift` remaps such a primary, so its remapped record
would otherwise look counted while the FASTQ payload pass ignores it. `--bait` is repeatable and
accepts a multi-record FASTA; `--bait-min-input-mb` overrides the bait lift's tiny-input fallback
threshold (0 disables the fallback, so reads are translated as well as remapped). All other
categories are still evaluated on the original `--in` alignment. `P.bait.cram` is kept so the
placement can be inspected in IGV. `MAPQ 0` (or `s2` ≈ `s1`) means the read fits the assembly
equally well (e.g. the integrated prophage copy): it is still selected, because for reassembly both
copies are the same sequence.

## `--bait-records`

**`--bait-records any|primary`** (default `any`) chooses which bait records select a read. Under
`any`, any counted record on a bait contig selects the read, primary or supplementary: a chimeric
read that reaches the bait only through a supplementary alignment (its primary record lies
elsewhere) is kept, and so are junction-spanning molecules. Under `primary`, a read is selected
only when its **primary** record is on a bait contig; reads that reach the bait solely through a
supplementary record are dropped. This matters when the supplementary-only molecules pull a
reassembly off a clean circle: on the FAM20853 demo the 9 such reads are the difference between B'
and the standard workflow's single circular contig (see
[case-study-phage-FAM20853.md](../case-studies/phage-FAM20853.md)). The `bait` block of `P.fish.json`
always records `mode` and the number of supplementary-only reads; an absent flag and an explicit
`--bait-records any` are equivalent.

## Outputs

| file | contents |
|---|---|
| `P.fastq.gz` | the selected reads (BGZF), each once, in original orientation, sorted by read name |
| `P.fish.tsv` | one row per (read, reason, finding) |
| `P.structure.tsv` | one row per hairpin, concatemer or tandem read, selected or not (`read, length, structure, contig, strand, period, units, arms, selected`); always written (header only when empty) |
| `P.lengths.tsv` | log-spaced read-length histogram (`bin_lo`, `reads`, `bp`) across every counted read (primary + unmapped); always written (header only when empty) |
| `P.fish.json` | selected reads and bp per reason and in total, reads dropped by `--clip-cap`, `hairpin_detected`, `concatemer_detected`, `tandem_detected` and `origin_crossing_skipped` counts, input record counts, per-contig circularity, `origin_crossing_reads` and `multimer` statistics (`contigs`), categories and thresholds, timing; a `bait` block with `--bait` |
| `P.bait_ref.fa` (with `--bait`) | assembly + bait contigs, indexed |
| `P.bait.cram` (with `--bait`) | the bait lift, with index |

The TSVs and the decompressed FASTQ are byte-identical across thread counts (only the JSON's timing
differs). Columns and keys: [output-format.md](../reference/output-formats.md#fish). An empty input
(0 records) emits `warning: <in>: input has no records` to stderr and writes empty `P.fastq.gz`,
header-only `P.fish.tsv`, `P.structure.tsv` and `P.lengths.tsv`, and zero counts in `P.fish.json`.

## Demo: recovering a phage missing from an assembly

A worked example of `qc` + `fish --bait` + reassembly. PacBio HiFi reads of two bacterial isolates
(FAM20853, FAM29080) were mapped to the automated `autocycler` assembly, which lacks a phage: the
phage is integrated in the chromosome as a prophage and is also present as free circular copies.
The curated assembly has the phage as an extra contig. The data are local, not in the repo.

A separate case study, [phage-FAM20853.md](../case-studies/phage-FAM20853.md), runs the complete
workflow (`lift` / `qc` / `fish --bait` / reassembly) on a different experiment: **ONT** reads of
FAM20853 only, not the HiFi demo above. It compares it with a standard minimap2 + Clair3 +
Sniffles2 script (wall/CPU time, read selection, reassembly outcome).

```bash
# 1. map the reads (-Y is required by fish: supplementary records must carry SEQ)
minimap2 -ax map-hifi --secondary=no -Y asm.fa reads.fq.gz \
  | samtools sort -O CRAM --reference asm.fa -o asm.cram
samtools index asm.cram

# 2. flag where the reads disagree with the assembly
cramhook qc --in asm.cram --ref asm.fa --out-prefix qc

# 3. select the reads, with a candidate sequence (here the phage) as bait
cramhook fish --in asm.cram --ref asm.fa --findings qc.findings.bed \
  --haplotypes qc.haplotypes.tsv --bait phage.fa --out-prefix fish

# 4. reassemble the selected reads
myloasm fish.fastq.gz --hifi -o reasm
```

`qc` reports each prophage as one finding spanning the whole element: the reads joining the
prophage end to its start (`split_junction`), the depth steps at both ends and the clips there are
linked, and everything in between is pooled. FAM20853: one orange finding at 1,862,313-1,903,595
(clips on both sides, `depth_high`, depth steps, junction); FAM29080: one red finding at
272,876-293,622 (also linked minority SNPs). The `likely` hint names the read pattern, not the
biology (`collapsed_repeat` here): an excised phage looks like a collapsed repeat with a junction.

`--bait` is needed because most phage reads map cleanly inside the prophage copy and show no clip, so
the default `fish` (findings scope) finds only 15 % / 31 % of them. The bait is a candidate sequence
(the numbers below use the phage contig itself): `--bait` asks which reads a full remap to
assembly + bait would place on the bait.

| | FAM20853 (40,865 bp phage) | FAM29080 (20,262 bp) |
|---|---|---|
| default fish: selected / truth found | 281 / 15.3 % | 101 / 31.2 % |
| `--locus` | 12,291 / 100 % | 6,414 / 100 % |
| `--bait` | 1,280 / 100 % | 251 / 100 % |
| fish `--bait` time (16 threads) | 3.5 s | 2.2 s |
| myloasm on `--bait` reads | one circular contig, 99.6 % of the phage at 100 % identity | one circular contig, 100 % at 100 % identity |
| myloasm on default-fish reads | a 40.9 kb contig covering the phage (depth 4, "possibly" circular) | one circular contig, 100 % (depth 7) |

Truth = reads that a full remap to the curated assembly places on the phage contig (1,180 and 218
reads). `--bait` recovers every truth read with few extra reads, and myloasm rebuilds each phage as a
circular contig at full read depth (the truth reads alone give a 96.5 % contig for FAM20853).
`--locus` (all reads within 5 kb of a finding) also catches every truth read because the finding now
spans the prophage, at the cost of ~50x more reads.

On ONT reads of the same isolates `--bait` recovered 70-74 % of the truth reads (interior reads tie
between the prophage and the phage copy, and minimap2 picks one at random), and the reassembly still
recovered each phage completely. On ONT `qc` flags both: FAM29080 as above, FAM20853 (free phage only
~1.1x deep, no depth signal) as one orange finding from the junction reads alone.

The evaluation script is `scripts/fish_phage_eval.py`.
