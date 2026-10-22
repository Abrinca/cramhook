# Install

`cramhook` runs on Linux (x86_64).

## Binary

Download the static Linux binary from the
[GitHub releases page](https://github.com/Abrinca/cramhook/releases) and put it on your `PATH`.
The release archive `cramhook-<version>-x86_64-linux.tar.gz` contains the `cramhook` executable,
`README.md`, `CHANGELOG.md` and `LICENSE`; a `.sha256` file is published next to it.

```bash
sha256sum -c cramhook-<version>-x86_64-linux.tar.gz.sha256
tar xzf cramhook-<version>-x86_64-linux.tar.gz
install -m 755 cramhook-<version>-x86_64-linux/cramhook ~/.local/bin/
```

## Build from source

The repository ships a [pixi](https://pixi.sh) workspace that provides the Rust toolchain, the C
build dependencies, `minimap2` and `samtools`:

```bash
git clone https://github.com/Abrinca/cramhook
cd cramhook
pixi run cargo build --release
```

The binary is `target/release/cramhook`. To run it with the pixi-provided `minimap2` and
`samtools`, either call it through `pixi run` or put the environment on `PATH`:

```bash
export PATH="$PWD/.pixi/envs/default/bin:$PWD/target/release:$PATH"
```

## Runtime dependencies

`minimap2` and `samtools` must be on `PATH`. `lift` and `fish --bait` call them (`lift` also uses
`minimap2` for the assembly diff); `qc` and plain `fish` read the alignment directly.

| tool | tested with | pixi requirement |
|---|---|---|
| `minimap2` | 2.31-r1302 | `>=2.31` |
| `samtools` / htslib | 1.24 | `>=1.24` |

Two environment variables select other executables:

| variable | overrides |
|---|---|
| `CRAMHOOK_MINIMAP2` | only the minimap2 that remaps deferred reads |
| `CRAMHOOK_SAMTOOLS` | only the `samtools` that sorts the output |

The assembly diff and the indexing still use `PATH`.

## Verify the install

```bash
cramhook --version
minimap2 --version
samtools --version | head -1
cramhook --help
```

`cramhook --help` lists `lift`, `qc` and `fish`, plus the debug commands `diff` and `translate`.

For a full check on real data, use the public test data set:

```bash
pixi run fetch-testdata     # download, verify sha256, unpack into testdata/data/
pixi run release-check      # cargo test + lift/qc/fish on the data + invariants + diff report
```

## Public test data

A slim data set (sample FAM29080, ONT R10 and PacBio HiFi reads, assemblies and ready-made CRAMs,
78.9 MB) is published as the `testdata-v1` release asset on GitHub. You can map the reads yourself
and follow the [quickstart](quickstart.md) with it. Contents, provenance and the checks that
`release-check` runs are in [Test data](project/testdata.md).
