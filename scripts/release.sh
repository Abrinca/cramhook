#!/usr/bin/env bash
# Publish a cramhook release from this machine. Run via
#   pixi run release            # full release-check first
#   pixi run release --no-cargo-test   # passed through to release-check
#
# Preconditions (checked): on main, clean tree, Cargo.toml version X.Y.Z has a
# `## [X.Y.Z]` CHANGELOG.md section, no GitHub release vX.Y.Z yet, and tag
# vX.Y.Z (if it exists) points at HEAD. Then: release-check, static tarball,
# a confirmation prompt, tag (annotated, if missing), push main + tag, and
# `gh release create` with the changelog section as notes.
set -euo pipefail
cd "$(dirname "$0")/.."

die() { echo "release: $*" >&2; exit 1; }

version="$(sed -n 's/^version *= *"\(.*\)"/\1/p' Cargo.toml | head -1)"
[[ -n "$version" ]] || die "no version in Cargo.toml"
tag="v$version"
name="cramhook-$version-x86_64-linux"

# --- preconditions --------------------------------------------------------
[[ "$(git rev-parse --abbrev-ref HEAD)" == main ]] || die "not on main"
[[ -z "$(git status --porcelain)" ]] || die "working tree not clean"

notes="$(mktemp)"
trap 'rm -f "$notes"' EXIT
awk -v h="## [$version]" 'index($0, h) == 1 {f = 1; next} /^## \[/ {f = 0} f' \
  CHANGELOG.md > "$notes"
grep -q '[^[:space:]]' "$notes" || die "CHANGELOG.md has no '## [$version]' section"

if git rev-parse -q --verify "refs/tags/$tag" >/dev/null; then
  [[ "$(git rev-parse "$tag^{commit}")" == "$(git rev-parse HEAD)" ]] \
    || die "tag $tag exists but does not point at HEAD"
  have_tag=1
else
  have_tag=0
fi
if gh release view "$tag" >/dev/null 2>&1; then
  die "GitHub release $tag already exists"
fi

# --- check and build ------------------------------------------------------
pixi run release-check "$@"
pixi run -e static dist

bin="dist/$name/cramhook"
[[ "$("$bin" --version)" == "cramhook $version" ]] || die "$bin reports '$("$bin" --version)'"
file "$bin" | grep -q 'static' || die "$bin is not statically linked"
[[ -z "$(git status --porcelain)" ]] || die "release-check/dist changed tracked files; commit them first"

# --- publish --------------------------------------------------------------
echo
echo "== $tag: notes =="
cat "$notes"
echo "== assets =="
ls -l "dist/$name.tar.gz" "dist/$name.tar.gz.sha256"
echo
read -r -p "Publish $tag (push main + tag, create GitHub release)? [y/N] " ok
[[ "$ok" == [yY] ]] || die "aborted (nothing pushed)"

(( have_tag )) || git tag -a "$tag" -m "cramhook $version"
git push origin main
git push origin "$tag"
gh release create "$tag" --verify-tag \
  --title "cramhook $version" \
  --notes-file "$notes" \
  "dist/$name.tar.gz" "dist/$name.tar.gz.sha256"
