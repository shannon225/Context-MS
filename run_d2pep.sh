#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

# Directories
IN_DIR="$(pwd)/2_linearcombo"
OUT_DIR="$(pwd)/3_pep"
mkdir -p "$OUT_DIR"


# Helper scripts (kept next to this bash file)
SCRIPT_DIR="$(cd -- "$(dirname "$0")" && pwd)"
PLOT_SCRIPT="${SCRIPT_DIR}/plot_pyispep_vs_score.py"
COUNT_SCRIPT="${SCRIPT_DIR}/count_pyispep.py"

# Column names in your scored input files
SCORE_COL="score"    # change to "LDA_score" if that's what your table uses
LABEL_COL="Label"    # expected 0/1

# 1) RUN d2pep on each scored target file from 2_linearcombo
#    Example file: tcell30p_scored_target_seed1_label0.tsv
targets=( "$IN_DIR"/*_scored_target_seed*_label0.tsv )
if ((${#targets[@]} == 0)); then
  echo "No target files found: $IN_DIR/*_scored_target_seed*_label0.tsv" >&2
  exit 1
fi

for target in "${targets[@]}"; do
  tb="$(basename "$target")"  # e.g., tcell30p_scored_target_seed1_label0.tsv

  # Extract seed from "...seedNNN..."
  if [[ "$tb" =~ seed([0-9]+) ]]; then
    seed="${BASH_REMATCH[1]}"
  else
    echo "Could not extract seed from: $tb" >&2
    continue
  fi

  # Base prefix before "_scored_target_seed"
  base_prefix="${tb%%_scored_target_seed*}"  # "tcell30p"

  # Output PEP file in 3_pep
  pep_out="${OUT_DIR}/${base_prefix}_pep_seed${seed}_ipspline.txt"

  if [[ -f "$pep_out" ]]; then
    echo "PEP exists, skipping d2pep: $pep_out"
  else
    echo "🔹 Seed ${seed}"
    echo "   input  : $target"
    echo "   pep_out: $pep_out"

    # d2pep (pyisotonic i-spline, q from FDR); mount IN/OUT so d2pep sees paths
    podman run --rm \
      -v "$IN_DIR":/in \
      -v "$OUT_DIR":/out \
      ghcr.io/statisticalbiotechnology/pyisotonicpep:main \
      d2pep \
        --cat-file "/in/$(basename "$target")" \
        --score-col "$SCORE_COL" \
        --label-col "$LABEL_COL" \
        --target-label 1 \
        --decoy-label 0 \
        --regression-algo ispline \
        --calc-q-from-fdr \
        --output "/out/$(basename "$pep_out")"
  fi

  # 2) PLOT PEP vs score, reading the PEP file from 3_pep
  if [[ -f "$PLOT_SCRIPT" ]]; then
    echo "   📈 Plotting: ${pep_out}.score_vs_pep.png"
    python3 "$PLOT_SCRIPT" "$pep_out"
  else
    echo "   ⚠️ Plot script not found: $PLOT_SCRIPT (skipping plot)"
  fi
done

# 3) COUNT PEP/q < 0.01 across all seeds in 3_pep
if [[ -f "$COUNT_SCRIPT" ]]; then
  echo "🧮 Summarizing pyIsoPEP counts (0.01) in $OUT_DIR ..."
  python3 "$COUNT_SCRIPT" \
    --indir "$OUT_DIR" \
    --pattern "*_pep_seed*_ipspline.txt" \
    --thresh 0.01 0.05 \
    --out "$OUT_DIR/pyispep_counts_0p01_dedup.tsv"
else
  echo "⚠️ Count script not found: $COUNT_SCRIPT (skipping summary)"
fi

echo "✅ Done. Outputs in: $OUT_DIR/"
