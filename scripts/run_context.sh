#!/usr/bin/env bash
# run_context.sh
#
# Orchestrates the full Context pipeline:
#   1) run_percolator_on_background.sh  -> 1_percolator/
#   2) model_averager.sh                -> 1_percolator/ (in/out)
#   3) run_apply_weights.sh             -> 2_linearcombo/
#   4) run_fdr.sh                       -> 2_linearcombo/ (add BH FDR/q)
#   5) run_q2pep.sh                     -> 3_pep/ (pyIsoPEP d2pep)
#   6) run_plot_pyispep.sh              -> plots from 3_pep/
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd -- "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Component scripts
PERCO="${SCRIPT_DIR}/run_percolator_on_background.sh"   # uses: feature/ -> 1_percolator/
AVER="${SCRIPT_DIR}/model_averager.sh"                  # in/out: 1_percolator/
APPLY="${SCRIPT_DIR}/run_apply_weights.sh"              # in: feature/,1_percolator/ -> out: 2_linearcombo/
FDR="${SCRIPT_DIR}/run_fdr.sh"                          # in/out: 2_linearcombo/  (calls bh_fdr.py)
Q2PEP="${SCRIPT_DIR}/run_q2pep.sh"                      # in: 2_linearcombo/ -> out: 3_pep/
PLOT="${SCRIPT_DIR}/run_plot_pyispep.sh"                # in: 3_pep/

# Make sure output dirs exist (most scripts also mkdir as needed)
mkdir -p "${SCRIPT_DIR}/1_percolator" "${SCRIPT_DIR}/2_linearcombo" "${SCRIPT_DIR}/3_pep"

# Ensure executability
chmod +x "$PERCO" "$AVER" "$APPLY" "$FDR" "$Q2PEP" "$PLOT" || true

echo "────────────────────────────────────────────────────────"
echo "🟦 [1/6] Percolator on background → 1_percolator/"
echo "────────────────────────────────────────────────────────"
"$PERCO"

echo
echo "────────────────────────────────────────────────────────"
echo "🟨 [2/6] Model averager (within-file lines 6/9/12) → 1_percolator/"
echo "────────────────────────────────────────────────────────"
"$AVER"

echo
echo "────────────────────────────────────────────────────────"
echo "🟧 [3/6] Apply weights to target features → 2_linearcombo/"
echo "────────────────────────────────────────────────────────"
"$APPLY"

echo
echo "────────────────────────────────────────────────────────"
echo "🟫 [4/6] BH FDR + q-values (run_fdr.sh) → 2_linearcombo/"
echo "────────────────────────────────────────────────────────"
"$FDR"

echo
echo "────────────────────────────────────────────────────────"
echo "🟥 [5/6] Run q2pep (pyIsoPEP d2pep) → 3_pep/"
echo "────────────────────────────────────────────────────────"
"$Q2PEP"

echo
echo "🟪 [6/6] Plot pyIsoPEP outputs → 3_pep/"
echo "────────────────────────────────────────────────────────"
"$PLOT"

echo
echo "✅ All done. Outputs:"
echo "   • 1_p_
