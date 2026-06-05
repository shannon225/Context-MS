#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# /proj/proteoforma_nsc/mamba_env/pymc_py12/bin/python /proj/proteoforma_nsc/Context/scripts/combine_boxplots.py

ENGINES = ("percolator", "mprophet")
ENGINE_COLORS = {"percolator": "#1f77b4", "mprophet": "#d62728"}
METRICS = ("q-value", "PEP")
LEVELS = ("PSMs", "peptides")
LAYOUT = [
    [None,       "tcell30p", "tcell10p"],
    ["tcell3p",  "tcell1p",  "tcell0pt3p"],
]


def per_engine_grid(base_dir, engine, level, metric):
    figs_dir = base_dir / engine / "figs"
    fig, axes = plt.subplots(2, 3, figsize=(21, 12))
    for r, row in enumerate(LAYOUT):
        for c, prefix in enumerate(row):
            ax = axes[r, c]
            ax.axis("off")
            if prefix is None:
                continue
            img_path = figs_dir / f"{prefix}.boxplot.{level}.{metric}.png"
            if img_path.exists():
                ax.imshow(mpimg.imread(img_path))
            else:
                ax.text(0.5, 0.5, f"missing\n{img_path.name}",
                        ha="center", va="center")
    fig.tight_layout()
    out_path = figs_dir / f"combined.boxplot.{level}.{metric}.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {out_path}")


def _load_counts(base_dir, engine, prefix):
    p = base_dir / engine / f"{prefix}.counts.tsv"
    if not p.exists():
        return None
    return pd.read_csv(p, sep="\t")


def _comparison_subplot(ax, counts_by_engine, level, metric, prefix):
    metric_key = {"q-value": "q", "PEP": "pep"}[metric]

    thresholds = None
    series = {}
    for eng in ENGINES:
        cdf = counts_by_engine.get(eng)
        if cdf is None:
            continue
        sub = cdf[(cdf["level"] == level) & (cdf["metric"] == metric_key)]
        if sub.empty:
            continue
        thrs = sorted(sub["threshold"].unique())
        if thresholds is None:
            thresholds = thrs
        series[eng] = [
            sub.loc[sub["threshold"] == t, "n_accepted"].to_numpy()
            for t in thrs
        ]

    if not series or thresholds is None:
        ax.text(0.5, 0.5, f"no data\n{prefix}", ha="center", va="center")
        ax.axis("off")
        return

    n_thr = len(thresholds)
    n_eng = len(series)
    width = 0.8 / n_eng
    centers = np.arange(1, n_thr + 1)
    handles, labels = [], []
    for i, (eng, data) in enumerate(series.items()):
        offset = (i - (n_eng - 1) / 2) * width
        positions = centers + offset
        bp = ax.boxplot(
            data, positions=positions, widths=width * 0.9,
            patch_artist=True, manage_ticks=False,
        )
        color = ENGINE_COLORS[eng]
        for box in bp["boxes"]:
            box.set(facecolor=color, alpha=0.5, edgecolor=color)
        for median in bp["medians"]:
            median.set(color="black", linewidth=1.2)
        handles.append(bp["boxes"][0])
        labels.append(eng)

    ax.set_xticks(centers)
    ax.set_xticklabels([f"{metric}<{t}" for t in thresholds])
    ax.set_ylim(-2, 150)
    ax.set_ylabel(f"# {level} accepted")
    ax.set_title(prefix, fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(handles, labels, loc="upper right", fontsize=9)


def comparison_grid(base_dir, level, metric):
    cross_figs_dir = base_dir / "figs"
    cross_figs_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for r, row in enumerate(LAYOUT):
        for c, prefix in enumerate(row):
            ax = axes[r, c]
            if prefix is None:
                ax.axis("off")
                continue
            counts = {eng: _load_counts(base_dir, eng, prefix) for eng in ENGINES}
            _comparison_subplot(ax, counts, level, metric, prefix)
    fig.suptitle(f"Per-seed acceptance: percolator vs mprophet "
                 f"({level}, by {metric})", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_path = cross_figs_dir / f"combined.compare.boxplot.{level}.{metric}.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {out_path}")


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--base-dir", type=Path, default=Path("results"),
                   help="Base directory holding <base>/<engine>/ subtrees.")
    p.add_argument("--engines", nargs="+", default=list(ENGINES),
                   choices=list(ENGINES))
    args = p.parse_args(argv)

    if not args.base_dir.is_dir():
        sys.exit(f"base dir not found: {args.base_dir}")

    for engine in args.engines:
        for level in LEVELS:
            for metric in METRICS:
                per_engine_grid(args.base_dir, engine, level, metric)

    if len(args.engines) >= 2:
        for level in LEVELS:
            for metric in METRICS:
                comparison_grid(args.base_dir, level, metric)


if __name__ == "__main__":
    main()
