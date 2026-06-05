#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

N_SEEDS="${N_SEEDS:-100}"
PREFIXES=(tcell0pt3p tcell1p tcell3p tcell10p tcell30p)
ENGINES="${ENGINES:-percolator mprophet}"
OUTDIR="${OUTDIR:-results}"
PYTHON="${COMBINE_PYTHON:-/proj/proteoforma_nsc/mamba_env/pymc_py12/bin/python3}"

mkdir -p "$OUTDIR/logs"

ALL_AGG_JIDS=()
for ENGINE in $ENGINES; do
    for PREFIX in "${PREFIXES[@]}"; do
        ARRAY_JID=$(sbatch --parsable \
            --job-name="ctx-${ENGINE}-${PREFIX}" \
            --array="1-${N_SEEDS}" \
            scripts/run_reproducibility_hpc.sh "$PREFIX" "$ENGINE")
        echo "[submit] ${ENGINE}/${PREFIX}: array job ${ARRAY_JID} (1-${N_SEEDS})"

        AGG_JID=$(sbatch --parsable \
            --job-name="ctx-agg-${ENGINE}-${PREFIX}" \
            --dependency="afterok:${ARRAY_JID}" \
            scripts/run_aggregate.sh "$PREFIX" "$ENGINE")
        echo "[submit] ${ENGINE}/${PREFIX}: aggregate job ${AGG_JID} (depends on ${ARRAY_JID})"
        ALL_AGG_JIDS+=("$AGG_JID")
    done
done

if ((${#ALL_AGG_JIDS[@]})); then
    DEP="afterok:$(IFS=:; echo "${ALL_AGG_JIDS[*]}")"
    COMB_JID=$(sbatch --parsable \
        --job-name="ctx-combine" \
        --dependency="$DEP" \
        --cpus-per-task=2 --mem=4G --time=00:05:00 \
        --output="$OUTDIR/logs/context_combine_%j.out" \
        --wrap="$PYTHON scripts/combine_boxplots.py --base-dir $OUTDIR")
    echo "[submit] combine job ${COMB_JID} (depends on all aggregates)"
fi
