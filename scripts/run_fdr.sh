#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

script_dir="$(cd -- "$(dirname "$0")" && pwd)"
cd "$script_dir"

indir="2_linearcombo"
python_script="${script_dir}/bh_fdr.py"
column_name="column"   # change this if your p-value column has a different name

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
