#!/usr/bin/env bash
#SBATCH --job-name=context-aggregate
#SBATCH --account=naiss2025-22-1406
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/proj/proteoforma_nsc/Context/results/logs/context_aggregate_%x_%j.out

# for p in tcell0pt3p tcell1p tcell3p tcell10p tcell30p; do bash scripts/run_aggregate.sh "$p"; done
set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$REPO_ROOT"

PREFIX="${1:?usage: run_aggregate.sh <prefix> [engine]}"
ENGINE="${2:-percolator}"
BASE_OUTDIR="${OUTDIR:-results}"
ENGINE_OUTDIR="$BASE_OUTDIR/$ENGINE"
PYTHON="${AGGREGATE_PYTHON:-/proj/proteoforma_nsc/mamba_env/pymc_py12/bin/python3}"

mkdir -p "$BASE_OUTDIR/logs"

"$PYTHON" scripts/aggregate.py \
    --target "features/${PREFIX}_features_target.tsv" \
    --prefix "$PREFIX" \
    --engine "$ENGINE" \
    --outdir "$ENGINE_OUTDIR"
