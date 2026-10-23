# Plan v4: release-ready repo (2026-10-19)

## Vision

A standalone tool for long-read bacterial assemblies. The product is the workflow
(`lift` → `qc` → `fish`): a faster, lighter and more sensitive replacement for a full remap plus
Clair3 and Sniffles2. The user's assembly GUI is one consumer (already in use), not the reason the
repo exists. `lift` began as an experiment in saving remap time; most of the value is now in `qc`
and `fish`.

- **Scope:** ONT R10 sup (`lr:hq`) and PacBio HiFi; bacterial isolates. Simple mixtures and
  haploid eukaryotes are later experiments.
- **Callers** run it as a command line. The public outputs are JSON, VCF and BED (the VCF and BED
  feed igv.js). Plain crashes are enough for errors. The caller controls all mapping, so the input
  contract holds.
- **`fish` selects reads only.** The caller runs the reassembly. Baits come from the user (e.g.
  copied from igv.js); unmapped reads can be assembled directly.
- **Speed target** for lift + qc: < 10 s ideal, < 30 s fine; huge datasets may take longer. The
  ONT case study already takes ~8.6 s on 16 threads.
- **qc keeps being tuned** (not frozen).
- **Used only by the author:** no backwards compatibility, no CI, no testing on other distros.

## Release decisions

GPL-3.0; a static Linux binary attached to GitHub releases; the tool gets a new name; tests run
locally.

## Steps

| # | Work |
|---|---|
| 0 | Decisions: the new name is **cramhook** (free on crates.io, GitHub, bioconda and PyPI as of 2026-10-19); public samples listed below |
| 1 | Repo basics: LICENSE (GPL-3.0), `--version`, CHANGELOG, semver from 0.1.0 |
| 2 | Slim public test data from the public samples: reads cut to ~20–30×, both platforms, at least one phage case (FAM29080, 20 kb phage), plus expected outputs |
| 3 | `pixi run release-check` (local): unit + synthetic tests, lift → qc → fish smoke test on the public data, qc evaluations |
| 4 | Static musl binary, built locally and attached to a GitHub release; minimap2 and samtools stay runtime dependencies, with their tested versions documented |
| 5 | Short user-facing README; reference material moves to `docs/`; JSON, VCF and BED documented under one shared format version, with format changes recorded in the CHANGELOG |
| 6 | Optional: `--progress`, one JSON line per phase on stderr |
| 7 | **Done 2026-10-19:** renamed to **cramhook** (binary, crate, `@PG` ID, docs). No compatibility: a lift `@PG` under the old name is an error |
| 8 | Guard for continued qc tuning: `release-check` compares the evaluation numbers against the previous release; any change goes in the CHANGELOG |

## Dependency note

`ref/rust-htslib` is an unmodified clone of rust-htslib v1.0.2, used via `[patch.crates-io]`
because 1.0.2 is not on crates.io yet. When it is published, drop the `[patch]` and `ref/rust-htslib`
and use the normal crates.io dependency. This does not block the release (binaries on GitHub
releases only; crates.io is not a target). rust-htslib is MIT, compatible with GPL-3.0.

**Done 2026-10-23:** 1.0.2 is on crates.io; the `[patch]` is dropped and the crates.io build is
used. The clones under `ref/` (including `ref/htslib`) are kept as read-only reference material.

## Public samples

FAM14197, FAM20838, FAM20843, FAM20847, FAM20851, FAM20853, FAM20856, FAM20862, FAM20866, FAM22000,
FAM29080 (ONT and PacBio).

## Parked

From `scratch/mrtomrod/POSSIBLE_NEXT_STEPS.md`: reducing junction-taint cost (3), inter-contig junction signal (5),
fish adapter-clip filter (6), origin `split_junction` (7). From `scratch/mrtomrod/ideas-lift.md`: small-edit
patching, the generalised split-read unit, the parallel pass. Items 2 (reproducible case study) and
4 (data drift) fit into steps 2 and 3 when useful.

## Later

Experiments with simple mixtures and haploid eukaryotes (larger genomes, more contigs).
