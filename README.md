# `Context`: Confidence Estimation for Targeted Proteomics

## Overview
`Context` produces PSM- and peptide-level q-values and PEPs for targeted
mass-spec assays. Because the target panel is usually too small to train
Percolator directly, `Context` trains Percolator on the **non-target** PSMs
from the same run (a large, statistically well-behaved background) and
transfers the learned linear discriminant to score the target panel.
[`pyIsoPEP`](https://pypi.org/project/pyIsoPEP/) (`q2pep`) then computes
FDR + q-values from TDC scores and fits an isotonic spline q -> PEP,
separately at the PSM and peptide levels.

`Context` is distributed as both a [**container image**](https://github.com/shannon225/Context/pkgs/container/context)
(Podman-, Docker- and Apptainer-compatible) and a [**Python package**](https://pypi.org/project/context/)
exposing the same CLI.

---

## Installation

### PyPI

```bash
pip install context
```

`Context` is a Python package with one external dependency: Percolator. For
PyPI installations, either a local `percolator` executable on `PATH` or Podman
must be available. `Context` uses the local executable when present; otherwise,
it runs `ghcr.io/percolator/percolator:master` through Podman.

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

**Output:** for `context run --nontarget NT --target TG --prefix P --outdir results/ --seed K`:

```
results/
  weights/P.seedK.weights.txt           # raw Percolator weights (3 CV-bin rows)
  P.seedK.rescored_features.tsv         # target features rescored with the seed weights
  P.seedK.psm.target.txt                # PSM-level target output
  P.seedK.peptide.target.txt            # peptide-level target output
```

---

## Command-line reference

```bash
context run -h
```

| flag | default | description |
|------|---------|-------------|
| `--nontarget FILE` | *required* | nontarget (background) feature TSV used to train Percolator |
| `--target FILE` | *required* | target panel feature TSV |
| `--prefix STR` | *required* | output prefix for results files |
| `--outdir DIR` | `results` | where outputs land |
| `--seed INT` | `1` | Percolator seed |
| `--container-cmd CMD` | `podman` | container runtime for the Percolator / `pyIsoPEP` fallback paths (`podman` or `docker`) |

---

## Examples

### PyPI

```bash
context run \
  --nontarget example/nontarget.tsv \
  --target    example/target.tsv \
  --prefix    run01 \
  --outdir    results
```

### Container image

```bash
# Podman
podman run --rm -v "$PWD:/work" -w /work \
  ghcr.io/shannon225/context:main \
  run --nontarget nontarget.tsv --target target.tsv \
      --prefix run01 --outdir results

# Apptainer
apptainer run --bind "$PWD:/work" --pwd /work context.sif \
  run --nontarget nontarget.tsv --target target.tsv \
      --prefix run01 --outdir results
```

---

## Links
* **PyPI package:** <https://pypi.org/project/context>
* **Container image:** <https://github.com/shannon225/Context/pkgs/container/context>
* **GitHub repository:** <https://github.com/shannon225/Context>
* **Percolator:** <http://percolator.ms>
* **pyIsoPEP:** <https://github.com/statisticalbiotechnology/smooth_q_to_pep>
