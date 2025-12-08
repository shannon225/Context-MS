#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

script_dir="$(cd -- "$(dirname "$0")" && pwd)"
cd "$script_dir"

# Input: label0 files from apply_weights
label0_dir="2_linearcombo"
# Output: pyIsoPEP results
pep_dir="3_pep"
mkdir -p "$pep_dir"

# Find all scored target label0 tables
label0_files=( "$label0_dir"/*_scored_target_*_label0.tsv )
if ((${#label0_files[@]} == 0)); then
  echo "❌ No *_scored_target_*_label0.tsv files found in ${label0_dir}/" >&2
  exit 1
fi

echo "🔹 Found ${#label0_files[@]} label0 files in ${label0_dir}/"

for f in "${label0_files[@]}"; do
  bn="$(basename "$f")"   # e.g. tcell100p_scored_target_seed1_label0.tsv

  # Extract seed number from the segment after '_scored_target_'
  seed_token="${bn#*_scored_target_}"   # "seed1_label0.tsv"
  seed_token="${seed_token%%_*}"        # "seed1"

  # If pattern is seed<NUM>, extract just the number, else keep token as is
  if [[ "$seed_token" =~ seed([0-9]+) ]]; then
    seed="${BASH_REMATCH[1]}"
  else
    seed="$seed_token"
  fi

  context_in="${pep_dir}/context_pep_seed${seed}.tsv"
  context_out="${pep_dir}/context_pep_seed${seed}_ipspline.txt"

  # Copy the label0 table into 3_pep as the context_pep file expected by pyIsoPEP
  echo "🔹 Preparing seed ${seed}:"
  echo "   label0: $f"
  echo "   ctx in: $context_in"
  echo "   ctx out: $context_out"

  cp "$f" "$context_in"

  # Run pyIsoPEP q2pep
  podman run -it --rm \
    -v "$pep_dir":/data \
    ghcr.io/statisticalbiotechnology/pyisotonicpep:main \
    q2pep \
      --cat-file "/data/$(basename "$context_in")" \
      --qcol qvalue \
      --score-col score \
      --label-col Label \
      --target-label 1 \
      --decoy-label 0 \
      --calc-q-from-fdr \
      --output "/data/$(basename "$context_out")"
done

echo "✅ pyIsoPEP q2pep completed for all seeds → ${pep_dir}/context_pep_seed*_ipspline.txt"
