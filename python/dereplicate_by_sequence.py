#!/usr/bin/env python3
"""
dereplicate_by_sequence.py

Dereplicate a TSV so there is only 1 row per peptide sequence.
If multiple rows share the same sequence, keep the row with the highest "primary" score.

- Input: tab-delimited file with a sequence column and a primary score column
- Output: tab-delimited file with 1 row per sequence

Typical usage before run_context.sh:
  python3 scripts/dereplicate_by_sequence.py --in results/1_features/foo.tsv --out results/1_features_derep/foo.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import pandas as pd


def norm_col(s: str) -> str:
    return "".join(ch.lower() for ch in s.strip() if ch.isalnum())


def resolve_column(df: pd.DataFrame, candidates: list[str], what: str) -> str:
    """
    Find a column in df that matches any candidate, case/whitespace/punct-insensitive.
    """
    nmap = {norm_col(c): c for c in df.columns}
    for cand in candidates:
        key = norm_col(cand)
        if key in nmap:
            return nmap[key]
    # fallback: try contains match (useful if columns like "primary_score" exist)
    for cand in candidates:
        key = norm_col(cand)
        for nkey, orig in nmap.items():
            if key and key in nkey:
                return orig
    raise SystemExit(
        f"ERROR: Could not find {what} column. Looked for: {candidates}. "
        f"Available columns: {list(df.columns)[:30]}{' ...' if len(df.columns) > 30 else ''}"
    )


def read_tsv_robust(path: Path) -> pd.DataFrame:
    # robust reader consistent with your Context-style TSVs
    return pd.read_csv(
        path,
        sep="\t",
        engine="python",
        dtype=str,
        na_filter=False,
        on_bad_lines="error",
    )


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Dereplicate TSV to 1 row per sequence using max primary score.")
    p.add_argument("--in", dest="inp", required=True, type=Path, help="Input TSV (features file).")
    p.add_argument("--out", dest="out", required=True, type=Path, help="Output TSV.")
    p.add_argument("--sequence-col", default="sequence",
                   help='Sequence column name (default: "sequence"). Case-insensitive matching is applied.')
    p.add_argument("--primary-col", default="primary",
                   help='Primary score column name (default: "primary"). Case-insensitive matching is applied.')
    p.add_argument("--keep-ties", choices=["first", "last"], default="first",
                   help="If multiple rows have identical max primary for a sequence, keep first or last (default: first).")
    p.add_argument("--drop-empty-seq", action="store_true",
                   help="Drop rows where the sequence is empty/blank (default: keep them).")
    args = p.parse_args(argv)

    if not args.inp.exists():
        raise SystemExit(f"ERROR: input file not found: {args.inp}")

    df = read_tsv_robust(args.inp)

    seq_col = resolve_column(df, [args.sequence_col, "sequence","Sequence", "peptide", "Peptide", "ModifiedSequence", "modifiedsequence"], "sequence")
    prim_col = resolve_column(df, [args.primary_col, "Primary", "primary_score", "PrimaryScore", "primaryscore"], "primary")

    # Clean sequence values
    df[seq_col] = df[seq_col].astype(str).str.strip()
    if args.drop_empty_seq:
        df = df[df[seq_col] != ""].copy()

    # Coerce primary to numeric for ranking
    primary_num = pd.to_numeric(df[prim_col].astype(str).str.strip(), errors="coerce")

    # If primary is missing/unparseable, treat as -inf so it won't win over real values
    primary_num = primary_num.fillna(float("-inf"))

    # Stable tie-breaking: preserve input order and select first/last among equal maxima
    df["_primary_num__"] = primary_num
    df["_row_order__"] = range(len(df))

    # Sort so "best" row per sequence is first (or last) within group
    # Highest primary first; for ties choose by row order.
    ascending = [True, False, True]  # seq asc, primary desc, row_order asc
    if args.keep_ties == "last":
        ascending = [True, False, False]  # choose last row among ties

    df_sorted = df.sort_values(by=[seq_col, "_primary_num__", "_row_order__"], ascending=ascending, kind="mergesort")

    # Keep one per sequence (first after sort)
    out_df = df_sorted.drop_duplicates(subset=[seq_col], keep="first").copy()

    # Drop helper columns
    out_df.drop(columns=["_primary_num__", "_row_order__"], inplace=True)

    # Write output
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)

    # Diagnostics
    sys.stderr.write(
        f"[derep] in={args.inp.name} rows={len(df)} unique_sequences={out_df[seq_col].nunique()} out={args.out}\n"
    )
    sys.stderr.write(
        f"[derep] sequence_col={seq_col} primary_col={prim_col} tie_break={args.keep_ties}\n"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
