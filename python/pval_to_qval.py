#!/usr/bin/env python3
"""
Usage:
    chmod +x pval_to_qval.py (run once to make executable)
    python3 pval_to_qval.py --input-dir /path/to/2_linearcombo --column column
References:
    [1] https://kaell.se/dsbook/statistics/multiple.html
    [2] https://github.com/statisticalbiotechnology/triqler/blob/master/triqler/qvality.py
"""

import argparse
import glob
import os
import numpy as np
import pandas as pd

DEFAULT_EPSILON = 1e-20
numLambda = 100
maxLambda = 0.95
def estimatePi0(pvalues, numBoot=100):
    pvalues = np.array(pvalues)
    lambdas, pi0s = list(), list()
    numPvals = len(pvalues)

    for lambdaIdx in range(numLambda + 1):
        l = ((lambdaIdx + 1.0) / numLambda) * maxLambda
        startIdx = np.searchsorted(pvalues, l, side='right')
        Wl = numPvals - startIdx
        pi0 = Wl / (1.0 - l) / numPvals
        if pi0 > 0.0:
            lambdas.append(l)
            pi0s.append(pi0)

    if len(pi0s) == 0:
        print(
            "Error in the input data: too good separation between targets and decoys.\n" \
            "Impossible to estimate pi0, setting pi0 = 1"
        )
        return 1.0

    minPi0 = min(pi0s)

    mse = [0.0] * len(pi0s)
    # Examine which lambda level is most stable under bootstrap
    for boot in range(numBoot):
        pBoot = bootstrap(pvalues, maxSize=numPvals)
        n = len(pBoot)
        for idx, l in enumerate(lambdas):
            startIdxBoot = np.searchsorted(pBoot, l, side='right')
            WlBoot = numPvals - startIdxBoot
            pi0Boot = WlBoot / (1.0 - l) / numPvals
            # Estimated mean-squared error.
            mse[idx] += (pi0Boot - minPi0) * (pi0Boot - minPi0)
    return max([min([pi0s[np.argmin(mse)], 1.0]), 0.0])


def bootstrap(allVals, maxSize=1000):
    return sorted(
        np.random.choice(allVals, min([len(allVals), maxSize]), replace=True)
    )


def pvalue_to_qvalue(pvalues, pi0):
    pvalues = np.asarray(pvalues)
    m = pvalues.size
    order = np.argsort(pvalues, kind="mergesort")
    pvalues_sorted = pvalues[order]
    ranks = np.arange(1, m + 1)
    fdr_sorted = (pi0 * m * pvalues_sorted) / ranks
    fdr_sorted = np.minimum(fdr_sorted, 1.0)
    qvalues_sorted = fdr_sorted.copy()
    qvalues_sorted = np.minimum.accumulate(qvalues_sorted[::-1])[::-1]

    fdr = np.empty_like(fdr_sorted)
    qvalues = np.empty_like(fdr_sorted)
    fdr[order] = fdr_sorted
    qvalues[order] = qvalues_sorted

    fdr = np.clip(fdr, DEFAULT_EPSILON, 1.0 - DEFAULT_EPSILON)
    qvalues = np.clip(qvalues, DEFAULT_EPSILON, 1.0 - DEFAULT_EPSILON)

    return fdr, qvalues


def fdr_to_qvalue(fdr):
    fdr = np.asarray(fdr)
    order = np.argsort(fdr, kind="mergesort")
    fdr_sorted = fdr[order]
    qvalues_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
    qvalues = np.empty_like(fdr_sorted)
    qvalues[order] = qvalues_sorted
    qvalues = np.clip(qvalues, DEFAULT_EPSILON, 1.0 - DEFAULT_EPSILON)

    return qvalues


def process_file(path: str, column: str):
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
    df = df.sort_values(column, kind="mergesort").reset_index(drop=True)
    pvals = df[column].values
    pi0 = estimatePi0(pvals)
    df["FDR"], df["q-value"] = pvalue_to_qvalue(df[column], pi0)

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
        default="p-value",
        help="Name of the p-value column to use (default: 'p-value')",
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