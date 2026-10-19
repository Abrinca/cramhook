#!/usr/bin/env bash
# Package the static binary: dist/cramhook-<version>-x86_64-linux.tar.gz, unpacking to
# cramhook-<version>-x86_64-linux/{cramhook,LICENSE,README.md,CHANGELOG.md}
# Run via `pixi run -e static dist` (builds first).
set -euo pipefail
cd "$(dirname "$0")/.."

version="$(sed -n 's/^version *= *"\(.*\)"/\1/p' Cargo.toml | head -1)"
name="cramhook-$version-x86_64-linux"
mkdir -p dist
rm -rf "dist/$name" "dist/$name.tar.gz"
mkdir -p "dist/$name"
cp target/x86_64-unknown-linux-musl/release/cramhook "dist/$name/cramhook"
strip "dist/$name/cramhook"
cp LICENSE README.md CHANGELOG.md "dist/$name/"
tar -C dist --sort=name --owner=0 --group=0 --numeric-owner -czf "dist/$name.tar.gz" "$name"
(cd dist && sha256sum "$name.tar.gz" > "$name.tar.gz.sha256")
ls -l "dist/$name" "dist/$name.tar.gz"
