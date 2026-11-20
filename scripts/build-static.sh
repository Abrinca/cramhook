#!/usr/bin/env bash
# Fully static x86_64 Linux (musl) build of cramhook. Run via
#   pixi run -e static build-static
# Output: target/x86_64-unknown-linux-musl/release/cramhook
set -euo pipefail
cd "$(dirname "$0")/.."

TARGET=x86_64-unknown-linux-musl
# Must match the rust pinned in pixi.toml [feature.static.dependencies].
RUST_VERSION=1.100.0
RUST_STD_SHA256=9d5211bc7b40195f9700a6f9133831832a32f50d0e2a2c0ae8af801e7656de5c

sysroot="$(rustc --print sysroot)"
have="$(rustc --version | cut -d' ' -f2)"
[[ "$have" == "$RUST_VERSION" ]] || { echo "rustc $have != $RUST_VERSION" >&2; exit 1; }

# conda-forge ships no musl rust-std; install the upstream one into this env.
# The stamp catches a stale copy left behind when pixi upgrades rust in place.
stamp="$sysroot/lib/rustlib/$TARGET/rust-std-version"
if [[ "$(cat "$stamp" 2>/dev/null)" != "$RUST_VERSION" ]]; then
  rm -rf "$sysroot/lib/rustlib/$TARGET"
  dl=scratch/static-build
  tarball="rust-std-$RUST_VERSION-$TARGET.tar.xz"
  mkdir -p "$dl"
  if [[ ! -f "$dl/$tarball" ]]; then
    curl -sSfL -o "$dl/$tarball.part" "https://static.rust-lang.org/dist/$tarball"
    echo "$RUST_STD_SHA256  $dl/$tarball.part" | sha256sum -c - \
      || { rm -f "$dl/$tarball.part"; exit 1; }
    mv "$dl/$tarball.part" "$dl/$tarball"
  fi
  echo "$RUST_STD_SHA256  $dl/$tarball" | sha256sum -c -
  tar -xJf "$dl/$tarball" -C "$dl"
  cp -r "$dl/rust-std-$RUST_VERSION-$TARGET/rust-std-$TARGET/lib/rustlib/$TARGET" "$sysroot/lib/rustlib/"
  echo "$RUST_VERSION" > "$stamp"
fi

zig_inc="$CONDA_PREFIX/lib/zig/libc/include"
export CC_x86_64_unknown_linux_musl="$PWD/scripts/zig-cc-musl.sh"
export AR_x86_64_unknown_linux_musl="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-ar"
# bindgen (hts-sys) must parse the htslib headers against musl, not the glibc sysroot.
export BINDGEN_EXTRA_CLANG_ARGS_x86_64_unknown_linux_musl="-nostdlibinc \
-isystem $zig_inc/x86_64-linux-musl -isystem $zig_inc/generic-musl \
-isystem $zig_inc/x86-linux-any -isystem $zig_inc/any-linux-any"
export ZIG_GLOBAL_CACHE_DIR="$PWD/target/zig-cache"

cargo build --release --locked --target "$TARGET" "$@"
