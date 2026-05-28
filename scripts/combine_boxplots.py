from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# /proj/proteoforma_nsc/mamba_env/pymc_py12/bin/python /proj/proteoforma_nsc/Context/scripts/combine_boxplots.py
FIGS_DIR = Path("/proj/proteoforma_nsc/Context/results/figs")

LAYOUT = [
    [None, "tcell30p", "tcell10p"],
    ["tcell3p", "tcell1p", "tcell0pt3p"],
]


def combine(metric, out_name):
    fig, axes = plt.subplots(2, 3, figsize=(21, 12))
    for r, row in enumerate(LAYOUT):
        for c, prefix in enumerate(row):
            ax = axes[r, c]
            ax.axis("off")
            if prefix is None:
                continue
            img_path = FIGS_DIR / f"{prefix}.boxplot.peptides.{metric}.png"
            ax.imshow(mpimg.imread(img_path))
    fig.tight_layout()
    out_path = FIGS_DIR / out_name
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


if __name__ == "__main__":
    combine("PEP", "combined.boxplot.peptides.PEP.png")
    combine("q-value", "combined.boxplot.peptides.q-value.png")
