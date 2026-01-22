#!/usr/bin/env python3
"""
deduplicate_features_by_primary.py

Removes duplicate rows that share the same 'id' in Skyline/Percolator-style
feature files, keeping only the row with the highest 'primary' value.

If multiple rows have the same 'id' and 'primary', one is kept arbitrarily.

Usage:
    python deduplicate_features_by_primary.py \
        --indir . \
        --pattern "tcell100p_features*.tsv" \
        --outdir dedup_features

Outputs:
    e.g. tcell100p_features_dedup.tsv
          tcell100p_features_nontarget_dedup.tsv
"""

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import sys


def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only the row with the highest 'primary' value per id."""
    if "id" not in df.columns or "primary" not in df.columns:
        missing = [c for c in ("id", "primary") if c not in df.columns]
        raise SystemExit(f"ERROR: Missing required column(s): {', '.join(missing)}")

    # Convert primary to numeric, replace non-numeric with -inf
    df["primary"] = pd.to_numeric(df["primary"], errors="coerce").fillna(float("-inf"))

    # Sort descending by primary, then drop duplicates by id
    df = df.sort_values(["id", "primary"], ascending=[True, False])
    df = df.drop_duplicates(subset="id", keep="first")

    return df


def main():
    p = argparse.ArgumentParser(description="Deduplicate feature files by highest primary per id.")
    p.add_argument("--indir", type=Path, default=Path("."), help="Input directory (default: current)")
    p.add_argument("--pattern", type=str, default="*.tsv", help="Glob pattern for input files")
    p.add_argument("--outdir", type=Path, default=None, help="Output directory (default: <indir>/dedup)")
    args = p.parse_args()

    indir = args.indir
    outdir = args.outdir or (indir / "dedup")
    outdir.mkdir(parents=True, exist_ok=True)

    files = sorted(indir.glob(args.pattern))
    if not files:
        sys.exit(f"No files matched: {indir}/{args.pattern}")

    for f in files:
        print(f"🔹 Processing {f.name} ...", file=sys.stderr)
        try:
            df = pd.read_csv(f, sep="\t", dtype=str)
            df.columns = df.columns.str.strip()
        except Exception as e:
            print(f"❌ Failed to read {f.name}: {e}", file=sys.stderr)
            continue

        before = len(df)
        df_dedup = deduplicate(df)
        after = len(df_dedup)

        out_path = outdir / f"{f.stem}_dedup.tsv"
        df_dedup.to_csv(out_path, sep="\t", index=False)
        print(f"✅ Wrote {out_path.name} ({before} → {after} rows)", file=sys.stderr)

    print(f"\n🎉 Done. Deduplicated files written to: {outdir}")


if __name__ == "__main__":
    main()
