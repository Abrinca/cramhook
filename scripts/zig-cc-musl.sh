#!/usr/bin/env bash
# C compiler for x86_64-unknown-linux-musl: zig cc with its bundled musl.
# cc-rs passes a rust-style --target, which zig does not understand; drop it.
args=()
for a in "$@"; do
  case "$a" in
    --target=*) ;;
    *) args+=("$a") ;;
  esac
done
exec zig cc -target x86_64-linux-musl "${args[@]}"
