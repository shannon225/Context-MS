#!/usr/bin/env bash
#SBATCH --job-name=context-reproducibility
#SBATCH --account=naiss2025-22-1406
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/proj/proteoforma_nsc/Context/results/logs/context_%x_%A_%a.out

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$REPO_ROOT"

PREFIX="${1:?usage: sbatch --array=1-100 run_reproducibility_hpc.sh <prefix> [engine]}"
ENGINE="${2:-percolator}"
BASE_OUTDIR="${OUTDIR:-results}"
ENGINE_OUTDIR="$BASE_OUTDIR/$ENGINE"
IMAGE="${CONTEXT_SIF:-$REPO_ROOT/context.sif}"

mkdir -p "$BASE_OUTDIR/logs"

if [[ ! -f "$IMAGE" ]]; then
    echo "Container image not found at $IMAGE" >&2
    echo "Build it once: apptainer build $IMAGE docker://ghcr.io/shannon225/context:containerize" >&2
    exit 1
fi

SEED="$SLURM_ARRAY_TASK_ID"
TAG="${PREFIX}.seed${SEED}"

apptainer run --bind "$REPO_ROOT:/work" --pwd /work "$IMAGE" \
    run \
        --nontarget    "features/${PREFIX}_features_nontarget.tsv" \
        --target       "features/${PREFIX}_features_target.tsv" \
        --prefix       "$PREFIX" \
        --outdir       "$ENGINE_OUTDIR" \
        --engine       "$ENGINE" \
        --seed         "$SEED" \
        --weights-out  "${TAG}.weights.txt" \
        --rescored-out "${TAG}.rescored_features.tsv" \
        --psm-out      "${TAG}.psm.target.txt" \
        --peptide-out  "${TAG}.peptide.target.txt"
