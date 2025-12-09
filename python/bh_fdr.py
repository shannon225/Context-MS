#!/usr/bin/env python3
"""
bh_fdr.py

Apply Benjamini–Hochberg FDR and monotone q-values to a column of p-values
in all *_label0.tsv files in a given directory.

Usage:
    python3 bh_fdr.py --input-dir /path/to/2_linearcombo --column score
"""

import argparse
import glob
import os

import numpy as np
import pandas as pd

DEFAULT_EPSILON = 1e-20


def fdr_to_qvalue(fdr: np.ndarray) -> np.ndarray:
    """
    Take an array of (possibly non-monotone) FDR estimates and return
    monotone non-decreasing q-values, using the cumulative-min trick.
    """
    fdr = np.asarray(fdr, dtype=float)

    # Sort FDRs in ascending order
    order = np.argsort(fdr, kind="mergesort")
    fdr_sorted = fdr[order]

    # Enforce monotonicity from right to left, then flip back
    qvalues_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]

    # Revert to original order
    qvalues = np.empty_like(qvalues_sorted)
    qvalues[order] = qvalues_sorted

    # Avoid exactly 0 or 1
    qvalues = np.clip(qvalues, DEFAULT_EPSILON, 1.0 - DEFAULT_EPSILON)
    return qvalues


def bh_fdr_from_pvalues(pvals: np.ndarray):
    """
    Classic Benjamini–Hochberg procedure.

    Parameters
    ----------
    pvals : array-like
        P-values for each test.

    Returns
    -------
    fdr : np.ndarray
        Raw BH FDR estimates (before monotone correction), in original order.
    qvalues : np.ndarray
        Monotone BH q-values (after fdr_to_qvalue), in original order.
    """
    pvals = np.asarray(pvals, dtype=float)
    m = pvals.size
    if m == 0:
        return np.array([]), np.array([])

    # Sort p-values ascending
    order = np.argsort(pvals)
    p_sorted = pvals[order]

    # Raw BH FDR estimates: p(i) * m / i
    ranks = np.arange(1, m + 1, dtype=float)
    fdr_sorted = p_sorted * m / ranks

    # Monotone q-values
    q_sorted = fdr_to_qvalue(fdr_sorted)

    # Map both back to original order
    fdr = np.empty_like(fdr_sorted)
    qvalues = np.empty_like(q_sorted)
    fdr[order] = fdr_sorted
    qvalues[order] = q_sorted

    return fdr, qvalues


def process_file(path: str, column: str = "score"):
    """
    Read a TSV, compute BH FDR/q-values on the given column,
    and write back in place.
    """
    print(f"  → Processing {os.path.basename(path)}")
    df = pd.read_csv(path, sep="\t")

    if column not in df.columns:
        raise ValueError(
            f"Column '{column}' not found in {path}. "
            f"Available columns: {list(df.columns)}"
        )

    pvals = df[column].values
    fdr, qvals = bh_fdr_from_pvalues(pvals)

    df["fdr"] = fdr
    df["qvalue"] = qvals

    # Overwrite in-place (you can change this to write a new file if you prefer)
    df.to_csv(path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Apply BH FDR + q-values to *_label0.tsv files in a directory."
    )
    parser.add_argument(
        "--input-dir",
        "-i",
        required=True,
        help="Directory containing *_label0.tsv files (e.g. 2_linearcombo)",
    )
    parser.add_argument(
        "--column",
        "-c",
        default="score",
        help="Name of the column to use (default: 'score')",
    )
    args = parser.parse_args()

    pattern = os.path.join(args.input_dir, "*_label0.tsv")
    files = sorted(glob.glob(pattern))

    if not files:
        raise SystemExit(f"No files matching *_label0.tsv in {args.input_dir}")

    print(f"Found {len(files)} files in {args.input_dir} matching *_label0.tsv")
    print(f"Using column '{args.column}' for p-values")

    for path in files:
        process_file(path)


if __name__ == "__main__":
    main()
