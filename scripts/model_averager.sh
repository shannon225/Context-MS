#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

script_dir="$(cd -- "$(dirname "$0")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

indir="${repo_root}/results/1_percolator"

# Loop over each weights file for each seed
for file in "$indir"/*features_nontarget_seed*_weights.txt; do
  bn=$(basename "$file")                                       # e.g., tcell100p_features_nontarget_seed1_weights.txt

  # Extract seed (robustly): take text after "_seed", then stop at next "_"
  seed_part=${bn#*_seed}                                       # "1_weights.txt"
  seed=${seed_part%%_*}                                        # "1"

  # Build a base WITHOUT the "_seed..." tail
  base_no_seed=${bn%_seed*}                                    # "tcell100p_features_nontarget"

  # Final output path (no literal '*')
  out="${indir}/${base_no_seed}_average_seed${seed}_weights.txt"

  # Compute the average across lines 6, 9, and 12 (within this single file)
  avg_line=$(awk '
    NR==6||NR==9||NR==12 {
      for (i=1; i<=NF; i++) { sum[i]+=$i; nf=(NF>nf?NF:nf); cnt++ }
    }
    END {
      for (i=1; i<=nf; i++) {
        val = (cnt ? sum[i]/cnt : 0)
        printf i<nf ? "%.6f\t" : "%.6f\n", val
      }
    }' "$file")

  # Replace lines 6, 9, and 12 with the averaged line
  awk -v avg="$avg_line" '
    NR==6 || NR==9 || NR==12 { print avg; next }
    { print }
  ' "$file" > "$out"

  echo "Wrote averaged file: $out"
done
