#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

script_dir="$(cd -- "$(dirname "$0")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
cd "$repo_root"

indir="results/3_pep"
plot_script="${script_dir}/plot_pyispep_vs_score.py"

# Collect all ipspline outputs (with basename + seed)
files=( "${indir}"/*_context_pep_seed*_ipspline.txt )

if ((${#files[@]} == 0)); then
    echo "⚠️  No *_context_pep_seed*_ipspline.txt files found in ${indir}/"
    exit 0
fi

for file in "${files[@]}"; do
    echo "🔹 Plotting $file ..."
    python3 "$plot_script" "$file"
done

echo "✅ All seed plots completed from ${indir}/."
