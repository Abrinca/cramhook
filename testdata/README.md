# Public test data

A slim, real data set for `pixi run release-check`: sample **FAM29080** (a bacterial isolate,
2.08 Mb: a 1.93 Mb chromosome and five plasmids, plus a 20,262 bp phage that is also integrated as
a prophage), ONT R10 and PacBio HiFi reads, made public with the author's approval.

```bash
pixi run fetch-testdata     # download + verify sha256 + unpack into testdata/data/ (git-ignored)
pixi run release-check      # cargo test + lift/qc/fish on the data + invariants + diff report
pixi run release-check --accept   # also overwrite testdata/expected/ (only if the invariants pass)
pixi run build-testdata     # rebuild the tarball from the local sources (author only)
```

`testdata/manifest.json` holds the download URL (the `testdata-v1` release asset on GitHub;
to use a local copy instead, run `pixi run fetch-testdata --from dist/cramhook-testdata-v1.tar.gz`) and the
tarball's sha256.

## Contents of `cramhook-testdata-v1.tar.gz` (78.9 MB)

| file | what it is |
|---|---|
| `ont/reads.fastq.gz` | ONT R10 sup reads, 6,090 reads, 53.4 Mb (~25x of the genome) |
| `ont/cons.fa` | autocycler consensus assembly (6 contigs, no phage) |
| `ont/polished.fa` | the same assembly after polishing (no phage) |
| `ont/cons.cram` | `ont/reads.fastq.gz` mapped to `cons.fa` |
| `hifi/reads.fastq.gz` | PacBio HiFi reads, 7,534 reads, 52.9 Mb (~25x) |
| `hifi/asm.fa` | autocycler assembly of the HiFi reads (6 contigs, no phage) |
| `hifi/asm.cram` | `hifi/reads.fastq.gz` mapped to `asm.fa` |
| `phage.fa` | the phage (`FAM29080-p1-1.1_scf7`, contig 7 of the PacBio-finalized assembly) |

Mean depth on the chromosome is 18x (ONT) and 22x (HiFi); the plasmids are deeper (ONT 66-242x,
HiFi 47-74x). At 18x, `qc` marks about half of the ONT chromosome as low coverage (`--lowcov-abs`
20x for `lr:hq`).

## Provenance and build

`scripts/testdata/build_testdata.py` (`pixi run build-testdata`) builds the tarball from the full
local read sets (ONT 668 Mb = 321x, HiFi 430 Mb = 207x):

```bash
seqtk sample -s 29080 <full ONT reads> 0.078  | bgzip -l 9 > ont/reads.fastq.gz
seqtk sample -s 29080 <full HiFi reads> 0.121 | bgzip -l 9 > hifi/reads.fastq.gz
# the input contract command (docs/commands/lift.md), run inside ont/ and hifi/
minimap2 -ax lr:hq   -t 8 --secondary=no -Y cons.fa reads.fastq.gz | samtools sort -O CRAM --reference cons.fa -o cons.cram
minimap2 -ax map-hifi -t 8 --secondary=no -Y asm.fa  reads.fastq.gz | samtools sort -O CRAM --reference asm.fa  -o asm.cram
# drop the @SQ UR:<absolute path> that htslib adds (no extra @PG), then index
samtools reheader -P -c "sed 's/\tUR:[^\t]*//'" cons.cram > noUR.cram && mv noUR.cram cons.cram   # same for asm.cram
samtools index cons.cram; samtools index asm.cram
```

Whole reads are subsampled (not regions). The build is deterministic for fixed tool versions
(minimap2 2.31-r1302, samtools/htslib 1.24, seqtk 1.5): sorted tar members, zero mtimes and owners,
fixed `-t 8` in the `@PG` line, no local paths in any packed header (checked by the builder). It writes `dist/cramhook-testdata-v1.tar.gz`, stores the sha256 in
`testdata/manifest.json` and unpacks into `testdata/data/`.

## Cases (`scripts/testdata/release_check.py`)

| case | runs | invariant |
|---|---|---|
| `ont_lift` (A) | `lift` cons.cram → polished.fa (`--min-input-mb 0`, so reads are translated), `qc` on the lifted CRAM, `fish --findings --haplotypes` | lift equals a full remap onto polished.fa (`scripts/compare.py`: 0 primary mismatches, 0 supplementary differences, depth max diff ≤ 1, ≥ 1 translated read); `qc` has a red/orange finding covering ≥ 80 % of the prophage `1:543393-563743` |
| `ont_bait` (B, ONT) | `fish --bait phage.fa` on the lifted CRAM (polished.fa lacks the phage) | recall of the truth phage reads ≥ 0.80 |
| `hifi_bait` (B, HiFi) | `qc` and `fish --bait phage.fa` on asm.cram (asm.fa lacks the phage) | `qc` has a red/orange finding covering ≥ 80 % of the prophage `FAM29080-p1-1.1_scf1:273076-293426`; recall ≥ 0.95 |

Truth = reads with a counted record on the phage contig in a full remap onto assembly + phage
(computed during the check). Measured on v1: lift 5,493 translated + 597 remapped reads, all
identical to the full remap; prophage findings orange on both platforms, each covering 100 % of the prophage (red at full HiFi depth);
`--bait` recall 22/22 (ONT) and 35/35 (HiFi).

## Expected outputs

`testdata/expected/<case>/` holds the small text outputs (qc VCF/BED/findings/JSON/haplotypes,
fish TSV/JSON, lift stats) plus a derived `summary.json`; JSON timing keys and `threads` are
removed. The diff report compares a run against them (per-tier site counts, findings added/lost/
changed tier, lift counts, fish reads per reason) and never fails. Outputs are identical across
thread counts. The run's full outputs (CRAMs, bedgraphs, FASTQs) are in `testdata/data/run/`.
