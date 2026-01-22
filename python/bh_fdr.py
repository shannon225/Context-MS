#!/usr/bin/env python3
"""
bh_fdr.py

1) Compute empirical two-sided p-values from a score column, using
   Label == 0 as the null distribution.
2) Apply Benjamini–Hochberg FDR and monotone q-values to the resulting
   p-values for the target rows (Label == 1).
3) Write pvalue, fdr, and qvalue columns back into each *_label0.tsv file.

Usage:
    python3 bh_fdr.py --input-dir /path/to/2_linearcombo --column score --label-col Label
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


def empirical_two_sided_p(scores: np.ndarray, null_scores: np.ndarray) -> np.ndarray:
    """
    Empirical two-sided p-values using the null distribution defined by null_scores.

    For each score s:
        p = 2 * min( P_null(S >= s), P_null(S <= s) )

    where probabilities are estimated empirically from null_scores.
    """
    scores = np.asarray(scores, dtype=float)
    null_scores = np.asarray(null_scores, dtype=float)

    if null_scores.size == 0:
        raise ValueError("Null distribution is empty (no rows with label == 0).")

    n_null = null_scores.size
    null_sorted = np.sort(null_scores)

    # counts of null <= s
    le_counts = np.searchsorted(null_sorted, scores, side="right")
    # counts of null >= s
    ge_counts = n_null - np.searchsorted(null_sorted, scores, side="left")

    le_prop = le_counts / n_null
    ge_prop = ge_counts / n_null

    pvals = 2.0 * np.minimum(le_prop, ge_prop)
    pvals = np.clip(pvals, DEFAULT_EPSILON, 1.0)
    return pvals


def process_file(path: str, score_col: str = "score", label_col: str = "Label"):
    """
    Read a TSV, compute empirical p-values and BH FDR/q-values, and write back in place.

    - Null distribution: rows with label_col == 0.
    - P-values: empirical two-sided, computed for target rows (label_col == 1).
    - BH FDR/q-values: applied to target p-values only.
    - Decoys (label_col == 0): pvalue, fdr, qvalue set to 1.0.
    """
    print(f"  → Processing {os.path.basename(path)}")
    df = pd.read_csv(path, sep="\t")

    # Basic column checks
    missing = [c for c in (score_col, label_col) if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing columns {missing} in {path}. "
            f"Available columns: {list(df.columns)}"
        )

    scores = df[score_col].values
    labels = df[label_col].values

    mask_decoy = labels == 0
    mask_target = labels == 1

    null_scores = scores[mask_decoy]
    target_scores = scores[mask_target]

    if target_scores.size == 0:
        print("    ⚠️ No target (label == 1) rows found; skipping BH FDR for this file.")
        # Still write out pvalue/fdr/qvalue=1.0 so downstream doesn't break
        df["pvalue"] = 1.0
        df["fdr"] = 1.0
        df["qvalue"] = 1.0
        df.to_csv(path, sep="\t", index=False)
        return

    # 1) Empirical two-sided p-values for target scores
    pvals_target = empirical_two_sided_p(target_scores, null_scores)

    # Initialize full-length arrays with 1.0 for decoys
    pvals_full = np.ones_like(scores, dtype=float)
    fdr_full = np.ones_like(scores, dtype=float)
    qvals_full = np.ones_like(scores, dtype=float)

    # 2) BH FDR/q-values on target p-values only
    fdr_target, qvals_target = bh_fdr_from_pvalues(pvals_target)

    # Map target-only results back into full arrays
    pvals_full[mask_target] = pvals_target
    fdr_full[mask_target] = fdr_target
    qvals_full[mask_target] = qvals_target

    # Add columns to the DataFrame
    df["pvalue"] = pvals_full
    df["fdr"] = fdr_full
    df["qvalue"] = qvals_full

    # Overwrite in-place (you can change this to write a new file if you prefer)
    df.to_csv(path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute empirical p-values from score/label, then apply BH FDR/q-values "
            "to *_label0.tsv files in a directory."
        )
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
        help="Name of the score column to use (default: 'score')",
    )
    parser.add_argument(
        "--label-col",
        "-l",
        default="Label",
        help="Name of the label column (default: 'Label', with 0 = null, 1 = target)",
    )
    args = parser.parse_args()

    pattern = os.path.join(args.input_dir, "*_label0.tsv")
    files = sorted(glob.glob(pattern))

    if not files:
        raise SystemExit(f"No files matching *_label0.tsv in {args.input_dir}")

    print(f"Found {len(files)} files in {args.input_dir} matching *_label0.tsv")
    print(f"Using score column '{args.column}' and label column '{args.label_col}'")

    for path in files:
        process_file(path, score_col=args.column, label_col=args.label_col)


if __name__ == "__main__":
    main()
