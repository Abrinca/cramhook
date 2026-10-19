# Development

```bash
pixi run cargo build --release
pixi run cargo test               # unit + integration tests (debug or --release)
pixi run cargo test --release -- --ignored   # the ignored random SR1–SR5 scenarios
pixi run sim && pixi run crams-synthetic && pixi run accept   # synthetic acceptance
pixi run real-smoke                                            # 5 real cases, ~3 min (local datasets)
pixi run real-eval -- --root data/m6large/real --bench 3      # real data (local datasets)
```

## Release check

```bash
pixi run fetch-testdata           # public FAM29080 test data into testdata/data/ (sha256-checked)
pixi run release-check            # cargo test --release + lift/qc/fish on the test data, ~2 min
pixi run release-check --accept   # also update testdata/expected/ (only if the invariants pass)
pixi run build-testdata           # rebuild the tarball from the local sources (author only)
```

`release-check` exits nonzero when a test or an invariant fails (lift equals a full remap,
qc flags the prophage on ONT and HiFi, `fish --bait` recall); the diff against
`testdata/expected/` is informational only. `--no-cargo-test` skips the cargo tests (~12 s then).
Cases, provenance and thresholds: [testdata/README.md](../testdata/README.md).

## Tests

`pixi run cargo test --release` runs 452 tests (3 ignored): 315 unit, 11 `diff_synthetic`,
25 `lift_synthetic`, 10 `lift_junction`, 1 `lift_roundtrip` (+1 ignored), 7 `translate_synthetic`
(+1 ignored), 26 `qc_synthetic`, 17 `qc_circularity`, 4 `qc_linkage`, 31 `fish` and 5 `progress`.

| Test target | Covers |
|---|---|
| `src/**` (unit) | diff, translate, contract, qc statistics/tiers/`likely`, output writers |
| `tests/diff_synthetic.rs` | `diff` identity blocks, taint (T1/T2) and novel regions |
| `tests/translate_synthetic.rs` | translate pass, split reads and contract over S1–S21 × hifi/ont |
| `tests/lift_synthetic.rs` | `lift` end-to-end, contract fallbacks, chained lifted input |
| `tests/lift_junction.rs` | junction-taint fixtures (deletion/join, origin competitor, nested junction, seam) |
| `tests/lift_roundtrip.rs` | A→B→A and A→B→A→B over S1–S21 × hifi/ont, each step vs its full remap |
| `tests/qc_synthetic.rs`, `qc_circularity.rs`, `qc_linkage.rs` | qc signals/outputs, circularity verdicts, cross-contig clip links, linkage pass |
| `tests/fish.rs` | fish selection, outputs, `--bait`, `--bait-records` |
| `tests/progress.rs` | `--progress` stream of `lift`, `qc`, `fish` (order, whole-percent `frac`, record counts vs `samtools view -c`, error line), unchanged outputs, and a run with a closed stderr pipe |

## Source layout

`ref/rust-htslib` holds the pinned rust-htslib 1.0.2 source used via `[patch]`.

| Path | Contents |
|------|----------|
| `src/diff.rs` | assembly diff: identity blocks, T1/T2 taint, T2-junction taint, junctions |
| `src/translate.rs` | eligibility rules, translation, clip screen, streaming pass |
| `src/lift.rs` | pipeline: child processes, merge/sort, checks, staging |
| `src/bamio.rs` | shared BAM/CRAM helpers: `is_cram`, field-mask `open_indexed`, `keep`, CIGAR ops, terminal clips, `aux_int` |
| `src/scan.rs` | shared alignment scan used by qc and fish (clips, indels, depth) |
| `src/qc/mod.rs` | `run_qc`: the qc pass that ties the modules together |
| `src/qc/{params,tiers,sites,sig,alerts,likely,findings,calling}.rs` | thresholds, tier rule, site types, significance, alert counts, `likely` table, findings, site calling and candidate collection |
| `src/qc/regions/{cluster,clip,junction,indel,depth,step,xlink}.rs` | region signals: SNP clusters, clip hotspots + cross-contig links, split junctions, large indels, depth bands and steps, xlinks |
| `src/qc/{circularity,junctions,link,output,stats,test_util}.rs` | circularity verdicts, split-read/cross-contig event extraction, linkage pass, VCF/BED/JSON/BEDGRAPH writers, statistics, test helpers |
| `src/fish.rs` | fish selection, reasons table, `--bait` lift |
| `src/contract.rs` | header and `@PG` contract |
| `src/version.rs` | the shared output `format_version` (see [output-format.md](output-format.md)) |
| `src/progress.rs` | the `--progress` stream: `Progress` handle, chunk and reference-position fractions |
| `src/cramcount.rs` | record count of a CRAM (container headers) or indexed BAM (index stats), for the `--progress` totals |
| `src/paf.rs`, `src/coords.rs`, `src/fasta.rs`, `src/main.rs` | PAF/`cs` parsing, typed coordinates, FASTA, CLI |
| `scripts/` | test data, oracle construction, comparison, acceptance, evaluation |
| `scripts/qc_truth.py` | qc truth set: raw → polished differences (SNP/INS/DEL/LINS/LDEL) |
| `scripts/qc_eval.py` | qc recall / false-alarm evaluation against truth, by read-support class |
| `scripts/qc_spike.py` | qc planted-error generator: spike events into a finalized assembly |
| `scripts/qc_spike_eval.py` | qc scoring of planted errors: recall, tiers, `likely` |
| `scripts/fish_eval.py` | fish selection evaluation against planted-event reads |
| `scripts/fish_phage_eval.py` | fish phage demo: recall vs the oracle, myloasm reassembly |
| `scripts/testdata/` | public test data: build, fetch, `release-check` |
| `testdata/` | test-data manifest (URL, sha256), expected outputs, README |
| `docs/plans/` (Plan-v1..v4), `docs/decisions.md`, `docs/results.md` | design history, decision log, measured results |
