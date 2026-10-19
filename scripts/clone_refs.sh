#!/bin/bash
set -e

REFS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/ref"
mkdir -p "$REFS_DIR"

# Define repositories
declare -A REPOS=(
    ["minimap2"]="https://github.com/lh3/minimap2"
    ["leviosam2"]="https://github.com/milkschen/leviosam2"
    ["hts-specs"]="https://github.com/samtools/hts-specs"
    ["htslib"]="https://github.com/samtools/htslib"
    ["rust-htslib"]="https://github.com/rust-bio/rust-htslib"
    ["noodles"]="https://github.com/zaeleus/noodles"
)

# Clone or skip existing
for name in "${!REPOS[@]}"; do
    url="${REPOS[$name]}"
    dir="$REFS_DIR/$name"

    if [ -d "$dir" ]; then
        echo "Skipping $name (already exists)"
    else
        echo "Cloning $name..."
        git clone --depth 1 "$url" "$dir"
    fi
done

echo "Reference repositories ready."
