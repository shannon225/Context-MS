#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


_SEED_RE = re.compile(r"\.seed(\d+)\.")
_METRIC_LABELS = {"q": "q-value", "pep": "PEP"}
Q_COL = "q-value"
PEP_COL = "posterior_error_prob"
LABEL_COL = "Label"


def parse_percolator_weights(path):
    lines = [ln.rstrip("\n") for ln in Path(path).read_text().splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]
    if len(lines) % 3 != 0:
        raise ValueError(f"{path}: expected 3 lines per CV bin, got {len(lines)}")
    k = len(lines) // 3
    feature_names = None
    w_raw, b_raw = [], []
    for i in range(k):
        header = [c.strip() for c in lines[3 * i].split("\t")][:-1]
        if feature_names is None:
            feature_names = header
        elif header != feature_names:
            raise ValueError(f"{path}: feature-name drift between CV bins")
        raw = np.array([float(v) for v in lines[3 * i + 2].split("\t")])
        w_raw.append(raw[:-1])
        b_raw.append(float(raw[-1]))
    return feature_names, np.stack(w_raw, axis=0), np.array(b_raw)


def parse_pyprophet_weights(path):
    df = pd.read_csv(path, sep="\t")
    feat = [str(x) for x in df["feature"].tolist()]
    w = df["weight"].to_numpy(dtype=float)
    bias_mask = [name == "__bias__" for name in feat]
    if any(bias_mask):
        bias_idx = bias_mask.index(True)
        b = float(w[bias_idx])
        feat = [f for f, m in zip(feat, bias_mask) if not m]
        w = w[[not m for m in bias_mask]]
    else:
        b = 0.0
    return feat, w[np.newaxis, :], np.array([b])


def parse_weights(path, engine):
    if engine == "percolator":
        return parse_percolator_weights(path)
    if engine == "pyprophet":
        return parse_pyprophet_weights(path)
    raise ValueError(f"unknown engine {engine!r}")


def boxplot_counts(prefix, counts_df, level, metric, png_path, *,
                   n_targets, n_decoys, engine):
    sub = counts_df[(counts_df["level"] == level) & (counts_df["metric"] == metric)]
    thresholds = sorted(sub["threshold"].unique())
    data = [sub.loc[sub["threshold"] == t, "n_accepted"].to_numpy()
            for t in thresholds]
    label = _METRIC_LABELS[metric]
    fig, ax = plt.subplots(figsize=(max(6.0, 2.4 + 1.6 * len(thresholds)), 4.8))
    ax.boxplot(data, tick_labels=[f"{label}<{t}" for t in thresholds])
    ax.set_ylabel(f"# {level} accepted")
    ax.set_ylim(-2, 150)
    ax.set_title(f"{prefix} [{engine}] ({n_targets} targets, {n_decoys} decoys)",
                 fontsize=11)
    ax.grid(alpha=0.3)
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _seed_files(directory, pattern):
    out = []
    for p in sorted(directory.glob(pattern)):
        m = _SEED_RE.search(p.name)
        if m:
            out.append((int(m.group(1)), p))
    return out


def aggregate(target_path, *, outdir, prefix, engine,
              q_thresholds, pep_thresholds):
    outdir = Path(outdir)
    weights_dir = outdir / "weights"
    figs_dir = outdir / "figs"
    figs_dir.mkdir(parents=True, exist_ok=True)

    weights_files = _seed_files(
        weights_dir, f"{prefix}.seed*.weights.txt",
    )
    if not weights_files:
        sys.exit(
            f"No weights files matching "
            f"{prefix}.seed*.weights.txt under {weights_dir}"
        )

    feature_names = None
    w_list, b_list = [], []
    for _, wpath in weights_files:
        feats, w_raw, b_raw = parse_weights(wpath, engine)
        if feature_names is None:
            feature_names = feats
        w_list.append(w_raw)
        b_list.append(b_raw)
    w_mean = np.concatenate(w_list, axis=0).mean(axis=0)
    b_mean = float(np.concatenate(b_list).mean())

    avg_w_path = outdir / f"{prefix}.avg_weights.tsv"
    pd.DataFrame({"feature": feature_names + ["__bias__"],
                  "weight": np.concatenate([w_mean, [b_mean]])}
                ).to_csv(avg_w_path, sep="\t", index=False)
    print(f"[write] {avg_w_path} (over {len(weights_files)} seeds)")

    target_df = pd.read_csv(target_path, sep="\t")
    labels = target_df[LABEL_COL].astype(int).to_numpy()
    n_targets = int((labels == 1).sum())
    n_decoys = int((labels == -1).sum())

    psm_seeds = _seed_files(outdir, f"{prefix}.seed*.psm.target.txt")
    pep_seeds = _seed_files(outdir, f"{prefix}.seed*.peptide.target.txt")
    if not (psm_seeds and pep_seeds):
        sys.exit("per-seed target tables not found; run the per-seed jobs first")

    rows = []
    metric_specs = (("q", Q_COL, q_thresholds), ("pep", PEP_COL, pep_thresholds))
    for level, seeds in (("PSMs", psm_seeds), ("peptides", pep_seeds)):
        for seed, p in seeds:
            df_s = pd.read_csv(p, sep="\t")
            for metric, col, thrs in metric_specs:
                for thr in thrs:
                    rows.append({"seed": seed, "level": level,
                                 "metric": metric, "threshold": thr,
                                 "n_accepted": int((df_s[col] < thr).sum())})
    counts = pd.DataFrame(rows)
    counts_out = outdir / f"{prefix}.counts.tsv"
    counts.to_csv(counts_out, sep="\t", index=False)
    print(f"[write] {counts_out}")
    for level in ("PSMs", "peptides"):
        for metric, _, _ in metric_specs:
            metric_label = _METRIC_LABELS[metric]
            boxplot_counts(
                prefix, counts, level, metric,
                figs_dir / f"{prefix}.boxplot.{level}.{metric_label}.png",
                n_targets=n_targets, n_decoys=n_decoys, engine=engine,
            )


def build_parser():
    p = argparse.ArgumentParser(
        prog="aggregate.py",
        description="Aggregate per-seed context outputs into averaged weights, "
                    "counts, and boxplots, for one engine.",
    )
    p.add_argument("--target", type=Path, required=True,
                   help="Target panel feature TSV (used to count targets/decoys).")
    p.add_argument("--prefix", required=True,
                   help="Output prefix shared with the per-seed runs.")
    p.add_argument("--outdir", type=Path, required=True,
                   help="Directory containing the per-seed results (typically "
                        "results/<engine>).")
    p.add_argument("--engine", choices=("percolator", "pyprophet"),
                   default="percolator",
                   help="Selects the weight-file format parser.")
    p.add_argument("--q-thresholds", type=float, nargs="+", default=[0.01, 0.05],
                   help="q-value cutoffs for the per-seed acceptance boxplots.")
    p.add_argument("--pep-thresholds", type=float, nargs="+", default=[0.01, 0.05],
                   help="PEP cutoffs for the per-seed acceptance boxplots.")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if not args.target.is_file():
        sys.exit(f"Not found: {args.target}")
    aggregate(args.target, outdir=args.outdir, prefix=args.prefix,
              engine=args.engine,
              q_thresholds=args.q_thresholds, pep_thresholds=args.pep_thresholds)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
