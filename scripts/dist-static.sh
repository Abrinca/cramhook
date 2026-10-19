#!/usr/bin/env bash
# Package the static binary: dist/cramhook-<version>-x86_64-linux{,.tar.gz}
# Run via `pixi run -e static dist` (builds first).
set -euo pipefail
cd "$(dirname "$0")/.."

version="$(sed -n 's/^version *= *"\(.*\)"/\1/p' Cargo.toml | head -1)"
name="cramhook-$version-x86_64-linux"
mkdir -p dist
cp target/x86_64-unknown-linux-musl/release/cramhook "dist/$name"
strip "dist/$name"
tar -C dist -czf "dist/$name.tar.gz" "$name"
(cd dist && sha256sum "$name.tar.gz" > "$name.tar.gz.sha256")
ls -l "dist/$name" "dist/$name.tar.gz"
