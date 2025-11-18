#!/usr/bin/env bash
# run_context.sh
# Orchestrates:
#  1) run_percolator_on_background.sh  -> 1_percolator/
#  2) model_averager.sh                -> 1_percolator/ (in/out)
#  3) run_apply_weights.sh             -> 2_linearcombo/
#  4) run_d2pep.sh                     -> 3_pep/  (also plots & counts)
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd -- "$(dirname "$0")" && pwd)"

# Component scripts (as provided)
PERCO="${SCRIPT_DIR}/run_percolator_on_background.sh"   # 1_percolator/
AVER="${SCRIPT_DIR}/model_averager.sh"                 # in/out: 1_percolator/
APPLY="${SCRIPT_DIR}/run_apply_weights.sh"              # in: 1_percolator/  out: 2_linearcombo/
D2PEP="${SCRIPT_DIR}/run_d2pep.sh"                      # in: 2_linearcombo/ out: 3_pep/

# Make sure output dirs exist (the step scripts also ensure these)
mkdir -p "${SCRIPT_DIR}/1_percolator" "${SCRIPT_DIR}/2_linearcombo" "${SCRIPT_DIR}/3_pep"

# Ensure executability (safe even if already +x)
chmod +x "$PERCO" "$AVER" "$APPLY" "$D2PEP" || true

echo "────────────────────────────────────────────────────────"
echo "🟦 [1/4] Percolator on background → 1_percolator/"
echo "────────────────────────────────────────────────────────"
# NOTE: Seeds are currently controlled inside run_percolator_on_background.sh (defaults to 1..5).
#       Edit that script if you want a different seed range.
"$PERCO"

echo
echo "────────────────────────────────────────────────────────"
echo "🟨 [2/4] Model averager (within-file lines 6/9/12) → 1_percolator/"
echo "────────────────────────────────────────────────────────"
"$AVER"

echo
echo "────────────────────────────────────────────────────────"
echo "🟧 [3/4] Apply weights to target features → 2_linearcombo/"
echo "────────────────────────────────────────────────────────"
"$APPLY"

echo
echo "────────────────────────────────────────────────────────"
echo "🟥 [4/4] Run d2pep + plots + counts → 3_pep/"
echo "────────────────────────────────────────────────────────"
"$D2PEP"

echo
echo "✅ All done. Outputs:"
echo "   • 1_percolator/   (weights, averages)"
echo "   • 2_linearcombo/  (scored target *_label0.tsv)"
echo "   • 3_pep/          (pyIsoPEP tables, plots, multi-threshold counts)"
