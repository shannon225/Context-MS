#!/usr/bin/env bash
set -euo pipefail

indir="3_pep"
plot_script="plot_pyispep_vs_score.py"

# Loop over all files matching pattern
for file in "${indir}"/context_pep_seed*_ipspline.txt; do
    # Skip if no files match
    [ -e "$file" ] || { echo "No context_pep_seed*_ipspline.txt files found in ${indir}/"; break; }

    echo "🔹 Plotting $file ..."
    python3 "$plot_script" "$file"
done

echo "✅ All seed plots completed."
