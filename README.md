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

**`percolator`** (default) and **`mprophet`** are supported for the discriminant step.

`Context` is distributed as both a [**container image**](https://github.com/shannon225/Context/pkgs/container/context)
(Podman-, Docker- and Apptainer-compatible) and a [**Python package**](https://pypi.org/project/context-ms)
exposing the same CLI.

---

## Installation

### PyPI

```bash
pip install context-ms
```

`Context` is a Python package. The `mprophet` engine has no external
dependency beyond NumPy/pandas. The `percolator` engine prefers a local
`percolator` executable on `PATH` and falls back to a container image via
Podman or Docker (set with `--container-cmd`, default `podman`):

* Percolator container: `ghcr.io/percolator/percolator:master`

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

Weights file format:

* `percolator`: Percolator's native `--weights` output (3 lines per CV bin).
* `mprophet`: a two-column TSV with `feature` and `weight`, plus a final
  `__bias__` row carrying the LDA constant.

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
| `--engine NAME` | `percolator` | `percolator` or `mprophet` |
| `--seed INT` | `1` | seed; passed to Percolator or to mprophet's RNG |
| `--container-cmd CMD` | `podman` | container runtime fallback (`podman` or `docker`) |
| `--input-profile NAME` | `encyclopedia` | mprophet-only; feature-column selection profile (`auto`, `pin`, `encyclopedia`) |
| `--seed-coefficients NAME_OR_PATH` | `encyclopedia` | mprophet-only; built-in name (`encyclopedia`, `none`) or path to a JSON file mapping feature names to seed-model coefficients |
| `--weights-out FILE` | `<prefix>.weights.txt` | weights output file name (or path); relative paths land under `<outdir>/weights` |
| `--rescored-out FILE` | `<prefix>.rescored_features.tsv` | rescored-features output file name (or path); relative paths land under `<outdir>` |
| `--psm-out FILE` | `<prefix>.psm.target.txt` | PSM-level output file name (or path); relative paths land under `<outdir>` |
| `--peptide-out FILE` | `<prefix>.peptide.target.txt` | peptide-level output file name (or path); relative paths land under `<outdir>` |

### Input profiles (mprophet only)

* `auto` — if any column starts with `var_` or `main_var_`, keep only
  those (OpenSWATH/pyprophet convention); otherwise fall back to
  `encyclopedia`.
* `pin` — keep every column the pin convention exposes as a feature.
* `encyclopedia` — drop the metadata columns Encyclopedia's
  `MProphetFeatureReader` excludes by name (`pepLength`, `charge1..4`,
  `precursorMass`, `RTinMin`, `midTime`, `numberOfMatchingPeaksAboveThreshold`,
  `primary`, `TD`).

`--input-profile` has no effect on the percolator engine.

### Seed coefficients (mprophet only)

A JSON dictionary mapping feature column names to starting linear
coefficients for the seed LDA. Names that don't appear in the input
contribute 0; if all entries drop out, the seed model is disabled
and inner iter 0 falls back to ranking by the single best feature.
Pass `none` to disable the seed model unconditionally.

---

## Examples
```bash
cd example
```
### PyPI

```bash
# Percolator (default)
context run \
  --nontarget nontarget.tsv \
  --target    target.tsv \
  --prefix    run01 \
  --outdir    results_run01

# mProphet
context run \
  --nontarget nontarget.tsv \
  --target    target.tsv \
  --prefix    run01 \
  --outdir    results_run01 \
  --engine    mprophet
```

### Container image

```bash
# Podman
podman run --rm -v "$PWD:/work" -w /work \
  ghcr.io/shannon225/context:main \
  run --nontarget nontarget.tsv --target target.tsv \
      --prefix run01 --outdir results_run01 --engine mprophet

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
* **Encyclopedia (mProphet reference):** <https://bitbucket.org/searleb/encyclopedia>
* **pyIsoPEP:** <https://github.com/statisticalbiotechnology/smooth_q_to_pep>
