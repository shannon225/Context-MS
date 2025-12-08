#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Directories
indir="feature"
outdir="1_percolator"
mkdir -p "$outdir"

# Collect matching inputs from indir (safe even if none match)
inputs=( "$indir"/*_features_nontarget*.tsv )

# Exit early (or just warn) if no files found
if ((${#inputs[@]} == 0)); then
  echo "No files matched: $indir/*_features_nontarget*.tsv"
  exit 0
fi

for infile in "${inputs[@]}"; do
  base="$(basename "$infile" .tsv)"

  # Loop over seeds 1 through 5
  for seed in {1..5}; do
    out_weights="$outdir/${base}_seed${seed}_weights.txt"
    console="$outdir/${base}_seed${seed}_console.txt"

    /mnt/c/users/m334793/percolator/build/src/percolator \
      --seed "$seed" \
      --weights "$out_weights" \
      --protein-report-duplicates \
      --post-processing-tdc \
      --override \
      "$infile" \
      > "$console" 2>&1

    echo "Processed $(basename "$infile") with seed $seed → $out_weights"
  done
done
