#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

script_dir="$(cd -- "$(dirname "$0")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

percolator_bin="${repo_root}/percolator/build/src/percolator"

if [[ ! -x "$percolator_bin" ]]; then
  echo "Error: percolator binary not found at: $percolator_bin" >&2
  echo "Hint: cd third_party/percolator && mkdir -p build && cd build && cmake .. && make -j4" >&2
  exit 1
fi

indir="${repo_root}/features"
outdir="${repo_root}/results/1_percolator"
mkdir -p "$outdir"

inputs=( "$indir"/*_features_nontarget*.tsv )

if ((${#inputs[@]} == 0)); then
  echo "No files matched: $indir/*_features_nontarget*.tsv"
  exit 0
fi

for infile in "${inputs[@]}"; do
  base="$(basename "$infile" .tsv)"

  for seed in {1..2}; do
    out_weights="${outdir}/${base}_seed${seed}_weights.txt"
    console="${outdir}/${base}_seed${seed}_console.txt"

    "${percolator_bin}" \
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
