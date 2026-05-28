# `Context`: Confidence Estimation for Targeted Proteomics

## Overview
`Context` produces PSM- and peptide-level q-values and PEPs for targeted
mass-spec assays. Because the target panel is usually too small to train a
discriminant directly, `Context` trains the discriminant on the **non-target**
PSMs from the same run (a large, statistically well-behaved background) and
transfers the learned linear discriminant to score the target panel.
[`pyIsoPEP`](https://pypi.org/project/pyIsoPEP/) (`q2pep`) then computes
FDR + q-values from TDC scores and fits an I-spline q -> PEP, separately at
the PSM and peptide levels.

**`percolator`** (default) and **`pyprophet`** are supported for the discriminant step.

`Context` is distributed as both a [**container image**](https://github.com/shannon225/Context/pkgs/container/context)
(Podman-, Docker- and Apptainer-compatible) and a [**Python package**](https://pypi.org/project/context-ms)
exposing the same CLI.

---

## Installation

### PyPI

```bash
pip install context-ms
```

`Context` is a Python package with two interchangeable external engines:
Percolator and pyProphet. For PyPI installations, `Context` prefers a local
executable on `PATH` (`percolator` or `pyprophet`); when neither is present
it falls back to a container image via Podman or Docker (set with
`--container-cmd`, default `podman`):

* Percolator container: `ghcr.io/percolator/percolator:master`
* pyProphet container: `ghcr.io/pyprophet/pyprophet:latest`

### Container image
```bash
podman pull ghcr.io/shannon225/context:main
# or
apptainer build context.sif docker://ghcr.io/shannon225/context:main
```

---

## Input & output

**Input:** two Percolator-compatible tsv files for non-target and target PSMs/peptides.
They must share an **identical header** (same columns, same order).

**Output:** for `context run --nontarget NT --target TG --prefix P --outdir results/ --engine E`:

```
results/
  weights/P.weights.txt           # trained weights
  P.rescored_features.tsv         # target features rescored with the trained weights
  P.psm.target.txt                # PSM-level target output
  P.peptide.target.txt            # peptide-level target output
```

Each output path can be overridden individually via `--weights-out`,
`--rescored-out`, `--psm-out`, `--peptide-out`. Values may be plain file
names (written inside `--outdir`, or inside `--outdir/weights` for the
weights file) or absolute paths.

---

## Command-line reference

```bash
context run -h
```

| flag | default | description |
|------|---------|-------------|
| `--nontarget FILE` | *required* | nontarget (background) feature TSV used to train the engine |
| `--target FILE` | *required* | target panel feature TSV |
| `--prefix STR` | *required* | output prefix for results files |
| `--outdir DIR` | `results` | output directory |
| `--engine NAME` | `percolator` | `percolator` or `pyprophet` |
| `--seed INT` | `1` | seed; passed to Percolator. |
| `--container-cmd CMD` | `podman` | container runtime fallback (`podman` or `docker`) |
| `--weights-out FILE` | `<prefix>.weights.txt` | weights output file name (or path); relative paths land under `<outdir>/weights` |
| `--rescored-out FILE` | `<prefix>.rescored_features.tsv` | rescored-features output file name (or path); relative paths land under `<outdir>` |
| `--psm-out FILE` | `<prefix>.psm.target.txt` | PSM-level output file name (or path); relative paths land under `<outdir>` |
| `--peptide-out FILE` | `<prefix>.peptide.target.txt` | peptide-level output file name (or path); relative paths land under `<outdir>` |

---

## Examples
```bash
cd example
```
### PyPI

```bash
# Percolator
context run \
  --nontarget nontarget.tsv \
  --target    target.tsv \
  --prefix    run01 \
  --outdir    results_run01

# pyProphet
context run \
  --nontarget nontarget.tsv \
  --target    target.tsv \
  --prefix    run01 \
  --outdir    results_run01 \
  --engine    pyprophet
```

### Container image

```bash
# Podman
podman run --rm -v "$PWD:/work" -w /work \
  ghcr.io/shannon225/context:main \
  run --nontarget nontarget.tsv --target target.tsv \
      --prefix run01 --outdir results_run01 --engine pyprophet

# Apptainer
apptainer run --bind "$PWD:/work" --pwd /work context.sif \
  run --nontarget nontarget.tsv --target target.tsv \
      --prefix run01 --outdir results_run01
```

---

## Links
* **PyPI package:** <https://pypi.org/project/context-ms>
* **Container image:** <https://github.com/shannon225/Context/pkgs/container/context>
* **GitHub repository:** <https://github.com/shannon225/Context>
* **Percolator:** <http://percolator.ms>
* **pyProphet:** <https://github.com/PyProphet/pyprophet>
* **pyIsoPEP:** <https://github.com/statisticalbiotechnology/smooth_q_to_pep>
