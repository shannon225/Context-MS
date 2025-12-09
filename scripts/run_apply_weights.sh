#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

script_dir="$(cd -- "$(dirname "$0")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

features_dir="${repo_root}/features"
weights_dir="${repo_root}/1_percolator"
outdir="${repo_root}/2_linearcombo"
mkdir -p "$outdir"

python_script="${repo_root}/python/apply_weights_to_features.py"

# Loop over all average weight files like:
#   tcell100p_features_nontarget_average_12_weights.txt
for wfile in "${weights_dir}"/*_average_*_weights.txt; do
  # If none match, the loop won’t run due to nullglob
  [ -e "$wfile" ] || { echo "No weight files found under ${weights_dir}."; break; }

  bn="$(basename "$wfile")"

  # Seed is the part after "_average_" up to next "_"
  # e.g., ..._average_12_weights.txt -> "12"
  seed_part="${bn#*_average_}"        # "12_weights.txt"
  seed="${seed_part%%_*}"             # "12"

  # Base prefix shared by target/nontarget:
  # strip everything from first "_features_nontarget" onward
  # e.g., "tcell100p_features_nontarget_average_12_weights.txt" -> "tcell100p"
  base_prefix="${bn%%_features_nontarget*}"

  # Corresponding features_target file
  feature_file="${features_dir}/${base_prefix}_features_target.tsv"
  if [[ ! -f "$feature_file" ]]; then
    echo "⚠️  Missing target features for base '${base_prefix}': ${feature_file} not found. Skipping seed ${seed}."
    continue
  fi

  out_prefix="${outdir}/${base_prefix}_scored_target_${seed}"

  echo "🔹 Processing base='${base_prefix}' seed=${seed}"
  echo "   weights:  $wfile"
  echo "   features: $feature_file"
  echo "   out:      ${out_prefix}*"

  python3 "$python_script" \
    --weights "$wfile" \
    --features "$feature_file" \
    --out "$out_prefix" \
    --write-label0 \
    --hist
done

echo "✅ Done. Results in ${outdir}/"