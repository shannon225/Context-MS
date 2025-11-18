#!/usr/bin/env python3
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def pick_column(df, preferred_exact, fallback_predicate):
    # Try exact (case-sensitive), then case-insensitive
    for name in preferred_exact:
        if name in df.columns: return name
    for name in preferred_exact:
        for c in df.columns:
            if c.lower() == name.lower():
                return c
    # Predicate fallback
    matches = [c for c in df.columns if fallback_predicate(c)]
    return matches[0] if matches else None

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_pyispep_vs_score.py <*_pep_seed*_ipspline.txt>")
        sys.exit(1)

    infile = Path(sys.argv[1])
    if not infile.exists():
        sys.exit(f"❌ File not found: {infile}")

    df = pd.read_csv(infile, sep="\t", comment="#", dtype=str)
    df.columns = df.columns.str.strip()

    # Score col: try LDA_score, Score, then a sensible fallback
    score_col = pick_column(
        df,
        preferred_exact=["LDA_score", "score", "Score"],
        fallback_predicate=lambda c: c.lower().endswith("score")
    )

    # PEP col: pyIsoPEP PEP (case/space tolerant)
    pep_col = pick_column(
        df,
        preferred_exact=["pyIsoPEP PEP"],
        fallback_predicate=lambda c: ("pyispep" in c.lower() and "pep" in c.lower() and "q" not in c.lower())
    )

    if not score_col or not pep_col:
        print(f"❌ Required columns not found.\nColumns available: {list(df.columns)}")
        sys.exit(1)

    # Convert numeric & drop NAs
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df[pep_col]   = pd.to_numeric(df[pep_col], errors="coerce")
    df = df.dropna(subset=[score_col, pep_col])

    # Scatter (log Y)
    plt.figure(figsize=(7, 5))
    plt.scatter(df[score_col], df[pep_col], s=10, alpha=0.6)
    plt.yscale("log")
    plt.xlabel(score_col)
    plt.ylabel(pep_col)
    plt.title(f"{infile.name}: {score_col} vs {pep_col}")
    plt.grid(alpha=0.3, which="both")
    plt.tight_layout()

    out_png = infile.with_suffix(".score_vs_pep.png")
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"📊 Saved plot to: {out_png}")

if __name__ == "__main__":
    main()
