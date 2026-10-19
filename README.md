<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/img/logo/logo-cramhook-dark-full.svg">
    <img src="docs/img/logo/logo-cramhook-light-full.svg" alt="cramhook" width="360">
  </picture>
</h1>

**A fast, light workflow for long-read bacterial assemblies** (ONT R10 sup reads mapped with
`lr:hq`, PacBio HiFi): check where the reads disagree with an assembly and pull out the reads that
explain it, without remapping everything after each assembly revision.

| command | what it does |
|---|---|
| `lift` | moves read alignments from an old assembly to a revised one (rotated, reoriented, polished, a contig added or removed). Reads in unchanged, unambiguous regions are translated, the rest are remapped with the same minimap2 command; typically **2–6× faster** than a full remap |
| `qc` | flags sites (VCF) and regions (BED) where the reads disagree with the assembly, tiered **red** / **orange** / **yellow** (plus **gray** for homopolymer indels), grouped into findings with a `likely` hint |
| `fish` | extracts the reads behind the flags (clips, large indels, unmapped, linked haplotypes) as FASTQ for reassembly; `--bait` adds the reads a full remap would place on a candidate sequence such as a phage |

Together they replace a full remap + Clair3 + Sniffles2 (+ a second remap to select reads). On one
ONT isolate with a phage missing from the new assembly, the whole workflow including reassembly
took **17.0 s instead of 98.7 s wall** (5.8×, 127.0 s instead of 644.6 s CPU) and recovered the
full phage sequence: see the [case study](docs/case-study-phage-FAM20853.md).

## Install

- **Binary:** download the static Linux binary from the
  [GitHub releases page](https://github.com/Abrinca/cramhook/releases) and put it on `PATH`.
- **From source:** `pixi run cargo build --release` (binary in `target/release/cramhook`).

Runtime dependencies: `minimap2` and `samtools` on `PATH` (used by `lift` and `fish --bait`;
tested with minimap2 2.31-r1302 and samtools/htslib 1.24).
`CRAMHOOK_MINIMAP2` and `CRAMHOOK_SAMTOOLS` override only the remap minimap2 and the output sort.

## Quick start

```bash
# map once, with exactly this command shape (the input contract; -Y is required by fish)
minimap2 -ax lr:hq -t 16 --secondary=no -Y old.fa reads.fastq.gz \
  | samtools sort -O CRAM --reference old.fa -o old.cram

# 1. the assembly was revised: move the alignments instead of remapping
cramhook lift --in old.cram --old-ref old.fa --new-ref new.fa --out new.cram \
  --threads 16 --stats lift.stats.json

# 2. flag where the reads disagree with the new assembly
cramhook qc --in new.cram --ref new.fa --out-prefix qc --threads 16

# 3. extract the reads behind the flags (optionally with a candidate sequence as bait)
cramhook fish --in new.cram --ref new.fa --out-prefix fish --threads 16 \
  --findings qc.findings.bed --haplotypes qc.haplotypes.tsv --bait phage.fa

# reassemble the selected reads with the assembler of your choice, e.g.
myloasm fish.fastq.gz -o reasm
```

Use `map-hifi` instead of `lr:hq` for PacBio HiFi. Without an assembly revision, skip `lift`: index
the mapping (`samtools index old.cram`) and run `qc` and `fish` on it directly. The input contract
(which minimap2 options, which `@PG` lines) is described in [docs/lift.md](docs/lift.md#input-contract).

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

The VCF and BED files load in IGV / igv.js; the JSON files are for programs. All of them share one
format version (`format_version` in JSON, `##format_version` in the VCF). Every column and key is
documented in [docs/output-format.md](docs/output-format.md).

## Documentation

- [docs/lift.md](docs/lift.md): how `lift` works, input contract, eligibility rules, lifted-input
  chains, CLI, performance and validation, debug subcommands
- [docs/qc.md](docs/qc.md): `qc` signals, tiers and FILTER, `likely` hint, linkage, performance,
  limitations
- [docs/fish.md](docs/fish.md): `fish` selection, `--bait`, `--bait-records`, the phage demo
- [docs/output-format.md](docs/output-format.md): every output file, the format version
- [docs/case-study-phage-FAM20853.md](docs/case-study-phage-FAM20853.md): the workflow versus
  remap + Clair3 + Sniffles2 on a real ONT dataset
- [docs/qc-validation.md](docs/qc-validation.md): `qc` validation on planted errors and
  raw → polished truth sets
- [docs/development.md](docs/development.md): tests and source layout
- [docs/results.md](docs/results.md), [docs/decisions.md](docs/decisions.md),
  [docs/plans/](docs/plans/): measured results, decision log, plans
- [CHANGELOG.md](CHANGELOG.md)

## Development

```bash
pixi run cargo build --release
pixi run cargo test               # unit + integration tests (debug or --release)
pixi run cargo test --release -- --ignored   # the ignored random SR1–SR5 scenarios
pixi run sim && pixi run crams-synthetic && pixi run accept   # synthetic acceptance
```

Real-data runs, the test targets and the source layout are in
[docs/development.md](docs/development.md).

## License

GPL-3.0-or-later; see [LICENSE](LICENSE).
