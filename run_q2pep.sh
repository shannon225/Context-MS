#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

folder="/mnt/c/Users/m334793/git/automating_context/target_output"

# Only match target files; derive decoy from each
targets=( "$folder"/scored_nontarget_seed*_weights.{txt,tsv} )
if ((${#targets[@]} == 0)); then
  echo "No files found: $folder/scored_nontarget_seed*_weights.{txt,tsv}" >&2
  exit 1
fi

# Optional: ensure folder exists
[[ -d "$folder" ]] || { echo "Folder not found: $folder" >&2; exit 1; }

for target in "${targets[@]}"; do
  base_t="$(basename "$target")"

  # Extract digits after 'seed'
  if [[ "$base_t" =~ seed([0-9]+) ]]; then
    seed="${BASH_REMATCH[1]}"
  else
    echo "Could not extract seed from: $base_t" >&2
    continue
  fi

  # Derive the decoy path by swapping _target -> _decoy (extension preserved)
  decoy="${target/_target./_decoy.}"
  base_d="$(basename "$decoy")"
  if [[ ! -f "$decoy" ]]; then
    echo "Missing decoy for seed $seed: $decoy" >&2
    continue
  fi

  output="$folder/scored_target_pep_seed${seed}_ipspline.txt"
  [[ -f "$output" ]] && { echo "Exists, skipping: $output"; continue; }

  echo "Processing seed $seed → $output"

         # Run the container
    podman run -it --rm \
      -v "$folder":/data \
      ghcr.io/statisticalbiotechnology/pyisotonicpep:main \
      d2pep \
        --cat-file "/data/context_pep_seed${seed}.tsv" \
        --score-col LDA_score \
        --label-col Label \
        --target-label 1 \
        --decoy-label 0 \
        --regression-algo ispline \
        --calc-q-from-fdr \
        --output "/data/context_pep_seed${seed}_ipspline.txt"
done
