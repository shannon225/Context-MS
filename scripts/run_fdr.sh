#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

script_dir="$(cd -- "$(dirname "$0")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
cd "$repo_root"

indir="2_linearcombo"
python_script="${repo_root}/python/bh_fdr.py"

column_name="pvalue"   # change this if your p-value column has a different name

if [[ ! -f "$python_script" ]]; then
  echo "❌ bh_fdr.py not found at: $python_script" >&2
  exit 1
fi

if [[ ! -d "$indir" ]]; then
  echo "❌ Directory not found: $indir" >&2
  exit 1
fi

echo "🔹 Running BH FDR/q-values on files in ${indir}/ using column '${column_name}' ..."
python3 "$python_script" --input-dir "$indir" --column "$column_name"

echo "✅ FDR/q-values added in-place in ${indir}/"
