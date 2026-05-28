#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

N_SEEDS="${N_SEEDS:-100}"
PREFIXES=(tcell0pt3p tcell1p tcell3p tcell10p tcell30p)
OUTDIR="${OUTDIR:-results}"

mkdir -p "$OUTDIR/logs"

for PREFIX in "${PREFIXES[@]}"; do
    ARRAY_JID=$(sbatch --parsable \
        --job-name="ctx-${PREFIX}" \
        --array="1-${N_SEEDS}" \
        scripts/run_reproducibility_hpc.sh "$PREFIX")
    echo "[submit] ${PREFIX}: array job ${ARRAY_JID} (1-${N_SEEDS})"

    AGG_JID=$(sbatch --parsable \
        --job-name="ctx-agg-${PREFIX}" \
        --dependency="afterok:${ARRAY_JID}" \
        scripts/run_aggregate.sh "$PREFIX")
    echo "[submit] ${PREFIX}: aggregate job ${AGG_JID} (depends on ${ARRAY_JID})"
done
