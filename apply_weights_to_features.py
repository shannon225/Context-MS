#!/usr/bin/env python3
"""
apply_weights_to_features.py

Compute LDA-style scores by applying a row of weights (with intercept) to a peptide-by-feature matrix.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import pandas as pd

INTERCEPT_KEYS = {"m0", "intercept", "bias", "constant", "c", "(intercept)"}


def read_weights_with_intercept(weights_path: Path, header_row: int = 4, weights_row: int = 6):
    """Read tab-delimited weights where header_row (1-based) has names and weights_row has coefficients."""
    # Read header row (1-based -> skip header_row-1)
    w_header = pd.read_csv(weights_path, sep="\t", skiprows=header_row - 1, nrows=1, header=0)
    w_header.columns = w_header.columns.str.strip()

    # Read weights row
    w_row = pd.read_csv(weights_path, sep="\t", skiprows=weights_row - 1, nrows=1, header=None)
    w_row.columns = w_header.columns
    w_row.columns = w_row.columns.str.strip()

    # Intercept
    intercept_cols = [c for c in w_row.columns if c.strip().lower() in INTERCEPT_KEYS]
    if not intercept_cols:
        raise ValueError(
            f"No intercept column found. Expected one of: {', '.join(sorted(INTERCEPT_KEYS))}"
        )
    intercept_col = intercept_cols[-1]
    intercept = float(w_row.iloc[0][intercept_col])

    # Feature weights (exclude intercept)
    weights = {c: float(w_row.iloc[0][c]) for c in w_row.columns if c != intercept_col}
    return weights, intercept, list(w_row.columns)


def read_features_robust(features_path: Path, metadata_cols: int, strict: bool = False) -> pd.DataFrame:
    """
    Robust TSV reader:
    - Gets the first non-empty, non-comment line as header
    - Enforces that header's column count
    - Skips malformed rows unless strict=True
    """
    # Grab header line
    header_line = None
    with open(features_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                header_line = line.rstrip("\n")
                break
    if header_line is None:
        raise SystemExit("Feature file appears empty or only comments.")
    cols = [c for c in header_line.split("\t") if c != ""]
    # Read with python engine; optionally skip bad lines
    on_bad = "error" if strict else "skip"
    fdf = pd.read_csv(
        features_path,
        sep="\t",
        engine="python",
        header=0,         # use file's first data header
        names=cols,       # enforce the header we just parsed
        usecols=cols,     # ignore any trailing blank columns
        dtype=str,
        na_filter=False,
        on_bad_lines=on_bad,
    )
    if fdf.shape[1] <= metadata_cols:
        raise SystemExit("Feature file has too few columns relative to --metadata-cols.")
    return fdf


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Apply weights (with intercept) to feature matrix.")
    p.add_argument("--weights", required=True, type=Path, help="Weights file (TSV).")
    p.add_argument("--features", required=True, type=Path, help="Feature file (TSV).")
    p.add_argument("--out", required=True, type=Path, help="Output path prefix (no extension).")
    p.add_argument("--metadata-cols", type=int, default=3,
                   help="Number of metadata columns to ignore from the start of features (default: 3).")
    p.add_argument("--header-row", type=int, default=4,
                   help="1-based row number containing the weights header (default: 4).")
    p.add_argument("--weights-row", type=int, default=6,
                   help="1-based row number containing the weights values (default: 6).")
    p.add_argument("--write-label0", action="store_true",
                   help="Also write a version with Label -1 changed to 0.")
    p.add_argument("--hist", action="store_true",
                   help="Also save a dual histogram image of LDA_score split by Label.")
    p.add_argument("--strict", action="store_true",
                   help="Fail on malformed feature lines instead of skipping them.")
    args = p.parse_args(argv)

    # Weights + intercept
    weights, intercept, weight_cols = read_weights_with_intercept(
        args.weights, header_row=args.header_row, weights_row=args.weights_row
    )

    # Features (robust read)
    fdf = read_features_robust(args.features, metadata_cols=args.metadata_cols, strict=args.strict)

    # --- Normalize Label column to numeric ints if present ---
    if "Label" in fdf.columns:
        fdf["Label"] = fdf["Label"].astype(str).str.strip()
        fdf["Label"] = pd.to_numeric(fdf["Label"], errors="coerce")
        fdf["Label"] = fdf["Label"].fillna(0).astype(int)

    # Split metadata vs features (ignore the first N metadata columns)
    meta_cols = args.metadata_cols
    feat_df = fdf.iloc[:, meta_cols:].copy()
    feat_df.columns = feat_df.columns.str.strip()

    # Align columns
    present = [c for c in weights if c in feat_df.columns]
    missing = [c for c in weights if c not in feat_df.columns]
    if not present:
        raise SystemExit("No overlapping feature columns between weights and features.")

    # Numeric coercion
    X = feat_df[present].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Score = dot + intercept
    coef = [weights[c] for c in present]
    fdf["score"] = (X * coef).sum(axis=1) + intercept

    # Write scored file
    out_tsv = args.out.with_suffix(".tsv")
    fdf.to_csv(out_tsv, sep="\t", index=False)

    # Optionally write label0 version
    if args.write_label0:
        fdf_label0 = fdf.copy()
        if "Label" in fdf_label0.columns:
            fdf_label0["Label"] = fdf_label0["Label"].where(fdf_label0["Label"] != -1, 0)
        out_label0 = args.out.with_name(args.out.name + "_label0").with_suffix(".tsv")
        fdf_label0.to_csv(out_label0, sep="\t", index=False)

    # Optional histogram (lazy import)
    if args.hist:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("WARN: matplotlib not installed; skipping histogram.", file=sys.stderr)
        else:
            if "Label" not in fdf.columns:
                print("WARN: 'Label' column not found; skipping histogram.", file=sys.stderr)
            else:
                targets = fdf.loc[fdf["Label"] == 1, "score"]
                decoys  = fdf.loc[fdf["Label"] == -1, "score"]

                plt.figure(figsize=(8, 5))
                plt.hist(decoys, bins=80, color="#2166ac", alpha=0.6, label="Decoy (Label = -1)")
                plt.hist(targets, bins=80, color="#b2182b", alpha=0.6, label="Target (Label = 1)")
                plt.xlabel("score")
                plt.ylabel("Peptide count")
                plt.title("Target vs Decoy Peptide Score Distribution")
                plt.legend()
                plt.grid(alpha=0.3)
                plt.tight_layout()
                out_png = args.out.with_name(args.out.name + "_hist").with_suffix(".png")
                plt.savefig(out_png, dpi=150)
                plt.close()

    # Diagnostics
    sys.stderr.write(
        f"[apply] weights={args.weights.name} features={args.features.name} -> out={out_tsv.name}\n"
    )
    sys.stderr.write(
        f"[apply] intercept={intercept} columns: provided={len(weights)} used={len(present)} missing={len(missing)}\n"
    )
    if missing:
        preview = ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else "")
        sys.stderr.write(f"[apply] missing (treated as 0): {preview}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
