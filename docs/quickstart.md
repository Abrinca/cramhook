# Quickstart

This page walks through the whole workflow: map once, lift to a revised assembly, flag where the
reads disagree, extract the reads behind the flags, and reassemble them. You need `cramhook`,
`minimap2` and `samtools` ([install](install.md)), your reads, and the old and new assemblies as
plain FASTA.

## 1. Map once

Map the reads to the old assembly with exactly this command shape. This is the
[input contract](commands/lift.md#input-contract); `-Y` is required by `fish`.

=== "ONT"

    ```bash
    minimap2 -ax lr:hq -t 16 --secondary=no -Y old.fa reads.fastq.gz \
      | samtools sort -O CRAM --reference old.fa -o old.cram
    ```

=== "HiFi"

    ```bash
    minimap2 -ax map-hifi -t 16 --secondary=no -Y old.fa reads.fastq.gz \
      | samtools sort -O CRAM --reference old.fa -o old.cram
    ```

`lift` reads the preset from the minimap2 `@PG` line, and so does `qc`. The only options allowed
there are the preset, `-t`, `--secondary=no` and `-Y`; do not post-process the CRAM with other
tools.

## 2. Lift to the revised assembly

```bash
cramhook lift --in old.cram --old-ref old.fa --new-ref new.fa --out new.cram \
  --threads 16 --stats lift.stats.json
```

Reads in unchanged, unambiguous regions are translated, the rest are remapped with the same
minimap2 command. `new.cram` and its `.crai` are written together, and every record is tagged
`XL:A:T` (translated) or `XL:A:R` (remapped). `lift.stats.json` has the counts and the deferral
reasons.

!!! note
    If the input does not satisfy the [input contract](commands/lift.md#input-contract), `lift`
    remaps all reads and records the reason in the stats. Add `--strict` to make it an error
    instead. If the preset cannot be recovered from the `@PG` line either, `lift` errors unless you
    supply `--preset lr:hq|map-hifi`.

## 3. Flag where the reads disagree

```bash
cramhook qc --in new.cram --ref new.fa --out-prefix qc --threads 16
```

Start with `qc.findings.bed` (one finding per problem location, with tier and `likely` hint) and
`qc.json` (the run summary). Load `qc.vcf` and `qc.bed` in IGV to see the sites and regions.

## 4. Extract the reads

```bash
cramhook fish --in new.cram --ref new.fa --out-prefix fish --threads 16 \
  --findings qc.findings.bed --haplotypes qc.haplotypes.tsv
```

Add `--bait phage.fa` to also select the reads that a full remap to the assembly plus the bait
would place on it, for a candidate sequence that is missing from the assembly (a phage, a plasmid):

```bash
cramhook fish --in new.cram --ref new.fa --out-prefix fish --threads 16 \
  --findings qc.findings.bed --haplotypes qc.haplotypes.tsv --bait phage.fa
```

`fish.fastq.gz` holds the selected reads and `fish.fish.tsv` says why each was selected. Bait
contig names must not clash with the names in `new.fa`.

## 5. Reassemble

`fish` only selects reads; use the assembler of your choice on `fish.fastq.gz`, for example:

```bash
myloasm fish.fastq.gz -o reasm            # ONT
myloasm fish.fastq.gz --hifi -o reasm     # HiFi
```

## If the assembly did not change

Skip `lift`. Index the mapping and run `qc` and `fish` on it:

```bash
samtools index old.cram
cramhook qc   --in old.cram --ref old.fa --out-prefix qc --threads 16
cramhook fish --in old.cram --ref old.fa --out-prefix fish --threads 16 \
  --findings qc.findings.bed --haplotypes qc.haplotypes.tsv --bait phage.fa
```

`qc` and plain `fish` (without `--bait`) work on any coordinate-sorted, indexed BAM or CRAM mapped
with `minimap2 -Y`; `qc` takes the preset from the `@PG` line unless you pass
`--preset lr:hq|map-hifi`.

`fish --bait` is stricter: it runs `lift` internally, so it needs a supported preset (`lr:hq` or
`map-hifi`) that can be recovered from the mapping's minimap2 `@PG` header (`fish` has no
`--preset` option). If the mapping does not follow the [input contract](commands/lift.md#input-contract)
but a preset is still recoverable, the bait lift remaps all reads instead of translating them.

## Outputs at a glance

| command | file | what it is |
|---|---|---|
| `lift` | `new.cram` + `.crai` | the alignments on the new assembly; each record tagged `XL:A:T` (translated) or `XL:A:R` (remapped) |
| `lift` | `lift.stats.json` (`--stats`) | counts, deferrals by reason, parameters, timings |
| `qc` | `qc.vcf` | tiered site calls (`TIER`, FILTER `PASS`/`LowCov`) |
| `qc` | `qc.bed` | tiered regions (clip hotspots, large indels, depth changes, split-read junctions, linked haplotypes), plus untiered `not_assessable`/`low_coverage` and, with `--crispr`, `crispr` array annotations |
| `qc` | `qc.findings.bed` | one finding per problem location, with tier and `likely` hint |
| `qc` | `qc.json` | run summary: alert and finding counts, per-contig depth and circularity, `qc_qv` |
| `qc` | `qc.depth.bedgraph` | per-position aligned-base depth |
| `qc` | `qc.gc.bedgraph` | GC fraction per `--gc-window` bp window (default 100; `0` disables) |
| `qc` | `qc.haplotypes.tsv` | reads of each linked haplotype group (input to `fish`) |
| `qc` | `qc.candidates.tsv` (`--candidates`) | every pre-filtered candidate site with all statistics |
| `fish` | `fish.fastq.gz` | the selected reads |
| `fish` | `fish.fish.tsv` | why each read was selected |
| `fish` | `fish.fish.json` | selection summary per reason, bait block |
| `fish` | `fish.bait_ref.fa`, `fish.bait.cram` (`--bait`) | assembly + bait, and the reads lifted onto it |

The VCF and BED files load in IGV / igv.js; the JSON files are for programs. Every column and key
is documented in [Output formats](reference/output-formats.md).

## Next steps

- [lift](commands/lift.md): eligibility rules, lifted-input chains, all flags.
- [qc](commands/qc.md): signals, tiers, the `likely` hint, limitations.
- [fish](commands/fish.md): selection reasons, `--bait`, `--bait-records`, `--clip-cap`.
- [Case study](case-studies/phage-FAM20853.md): this workflow on a real ONT dataset, with timings.
- [Test data](project/testdata.md): a small public data set to try it on.
