# fish: extract reads that disagree with the assembly

`fish` is the extraction step. Two passes over the alignment: the first applies its own per-read
checks and selects reads into a reasons table, the second decodes the FASTQ payload of the selected
reads only. `fish` selects reads only; the reassembly is up to the caller.

```bash
cramhook fish --in old.cram --ref asm.fa --out-prefix P [--threads N] [--linear c1,c2] \
  [--findings P.findings.bed] [--haplotypes P.haplotypes.tsv] \
  [--scope findings|genome] [--min-tier yellow|orange|red] [--finding-pad 1000] \
  [--no-unmapped] [--no-clips] [--no-indels] [--no-linked] \
  [--locus [--locus-flank 5000]] [--bait bait.fa] [--bait-min-input-mb MB] \
  [--bait-records any|primary] \
  [--clip-min-len 50] [--indel-min-len 50] [--indel-merge-gap 50] [--indel-piece-min 10] \
  [--progress]
```

| Flag | Default | Purpose |
|------|---------|---------|
| `--progress` | off | machine-readable progress on stderr, one JSON object per line (phases `pass`, `bait_lift` with `--bait`, `payload`, `output`); see [output-format.md](output-format.md#progress-stream---progress) |

## Selection

Counted records skip secondary (0x100), duplicate (0x400) and QC-fail (0x200); primary and
supplementary are both checked. **`fish` requires a `minimap2 -Y` alignment** (supplementary
records must carry SEQ): any hard clip in a counted record is an error. A read is selected when any
of its counted records qualifies for any enabled category:

| reason | a read is selected when | detail |
|---|---|---|
| `unmapped` | it has an unmapped record (0x4), placed or unplaced | `.` |
| `clip` | a terminal soft clip ≥ `--clip-min-len`; a clip at a contig end of a **circular** contig is an origin crossing (skipped), on `--linear` contigs it counts | `side=left;len=N` / `side=right;len=N` |
| `large_ins` / `large_del` | a merged large indel (`--indel-merge-gap`, `--indel-piece-min`) ≥ `--indel-min-len` | `len=N;at=P` |
| `linked` | the read is listed in `--haplotypes` for a group whose span overlaps a finding of tier ≥ `--min-tier` (findings scope) or any group (genome scope) | `group=N` |
| `locus` | any counted record overlaps a finding (tier ≥ `--min-tier`) ± `--locus-flank` | `.` |
| `bait` | with `--bait`: a full remap to assembly + bait places one of the read's counted records on a bait contig (see below) | `contig=<bait>;mapq=<MAPQ>;s1=<s1>;s2=<s2>;record=primary\|supplementary` |

Scope `findings` (the default when `--findings` is given) counts `clip`/`large_*` only when the
event position (clip position, indel start) is within a finding of tier ≥ `--min-tier` ±
`--finding-pad`; scope `genome` counts them anywhere and `unmapped` is genome-wide in both.
`--locus` and `--scope findings` require `--findings`. There is no positive `--linked` flag: linked
selection is on by default and requires `--haplotypes` (if no file is given it is silently skipped
and recorded in the JSON; `--no-linked` turns it off). Reads that merely end at a position, or agree
with the assembly, are not selected; site-only and depth-only findings contribute only through
`--locus`.

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
[case-study-phage-FAM20853.md](case-study-phage-FAM20853.md)). The `bait` block of `P.fish.json`
always records `mode` and the number of supplementary-only reads; an absent flag and an explicit
`--bait-records any` are equivalent.

## Outputs

| file | contents |
|---|---|
| `P.fastq.gz` | the selected reads (BGZF), each once, in original orientation, sorted by read name |
| `P.fish.tsv` | one row per (read, reason, finding) |
| `P.fish.json` | selected reads and bp per reason and in total, input record counts, categories and thresholds, timing; a `bait` block with `--bait` |
| `P.bait_ref.fa` (with `--bait`) | assembly + bait contigs, indexed |
| `P.bait.cram` (with `--bait`) | the bait lift, with index |

The TSV and the decompressed FASTQ are byte-identical across thread counts (only the JSON's timing
differs). Columns and keys: [output-format.md](output-format.md#fish).

## Demo: recovering a phage missing from an assembly

A worked example of `qc` + `fish --bait` + reassembly. PacBio HiFi reads of two bacterial isolates
(FAM20853, FAM29080) were mapped to the automated `autocycler` assembly, which lacks a phage: the
phage is integrated in the chromosome as a prophage and is also present as free circular copies.
The curated assembly has the phage as an extra contig. The data are local, not in the repo.

The full write-up — standard workflow versus `cramhook`, wall/CPU time, read selection and the
reassembly outcome — is [case-study-phage-FAM20853.md](case-study-phage-FAM20853.md).

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
