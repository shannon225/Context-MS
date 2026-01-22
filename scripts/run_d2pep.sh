#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

# Input: label0 files from apply_weights
label0_dir="results/2_linearcombo"
# Output: pyIsoPEP results
pep_dir="results/3_pep"
mkdir -p "$pep_dir"

# Find all scored target label0 tables
label0_files=( "$label0_dir"/*_scored_target_*_label0.tsv )
if ((${#label0_files[@]} == 0)); then
  echo "❌ No *_scored_target_*_label0.tsv files found in ${label0_dir}/" >&2
  exit 1
fi

echo "🔹 Found ${#label0_files[@]} label0 files in ${label0_dir}/"

for f in "${label0_files[@]}"; do
  bn="$(basename "$f" .tsv)"   # e.g. tcell100p_scored_target_seed1_label0

  # Extract seed number from the segment after '_scored_target_'
  seed_token="${bn#*_scored_target_}"   # "seed1_label0.tsv"
  seed_token="${seed_token%%_*}"        # "seed1"

  # If pattern is seed<NUM>, extract just the number, else keep token as is
  if [[ "$seed_token" =~ seed([0-9]+) ]]; then
    seed="${BASH_REMATCH[1]}"
  else
    seed="$seed_token"
  fi

  context_in="${pep_dir}/${bn}_context_pep.tsv"
  context_out="${pep_dir}/${bn}_pep_ipspline.txt"

  echo "🔹 Preparing seed ${seed}:"
  echo "   label0: $f"
  echo "   ctx in: $context_in"
  echo "   ctx out: $context_out"

  # Copy the label0 table into 3_pep as the context_pep file expected by pyIsoPEP
  cp "$f" "$context_in"

  # Absolute path for Podman volume (important!)
  host_pep_dir="$(cd "$pep_dir" && pwd)"

  # Work around WSLg XDG_RUNTIME_DIR ownership issue
  if [[ -n "${XDG_RUNTIME_DIR:-}" && ! -O "$XDG_RUNTIME_DIR" ]]; then
    unset XDG_RUNTIME_DIR
  fi

  # Run pyIsoPEP q2pep
  podman run --rm \
    -v "$host_pep_dir":/data \
    ghcr.io/statisticalbiotechnology/pyisotonicpep:main \
    d2pep \
      --cat-file "/data/$(basename "$context_in")" \
      --score-col score \
      --label-col Label \
      --target-label 1 \
      --decoy-label 0 \
      --calc-q-from-fdr \
      --output "/data/$(basename "$context_out")"
done

# Plotting the scores and PEP values from the ipspline outputs
repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

indir="$repo_root/results/3_pep"
plot_script="${repo_root}/python/plot_pyispep_vs_score.py"

# Troubleshooting
# echo "DEBUG: repo_root=$repo_root"
# echo "DEBUG: pwd=$(pwd)"
# echo "DEBUG: indir=$indir"
# ls -la "$indir" || true

# Collect all ipspline outputs (with basename + seed)
files=( "${indir}"/*_pep_ipspline.txt )

if ((${#files[@]} == 0)); then
    echo "⚠️  No *_pep_ipspline.txt files found in ${indir}/"
    exit 0
fi

for file in "${files[@]}"; do
    echo "🔹 Plotting $file ..."
    python3 "$plot_script" "$file"
done

echo "✅ All seed plots completed from ${indir}/."
