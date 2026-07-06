# `Context`: Confidence Estimation for Targeted Proteomics

## Overview
`Context` produces PSM- and peptide-level q-values and PEPs for targeted
mass-spec assays. Because the reference panel is usually too small to train a
discriminant directly, `Context` trains the discriminant on the **background**
PSMs from the same run (a large, statistically well-behaved background) and
transfers the learned linear discriminant to score the reference panel.
Two engines are supported, each paired with its own confidence-estimation
back-end:

* **`percolator`** (default) — Percolator trains the discriminant, and
  [`pyIsoPEP`](https://pypi.org/project/pyIsoPEP/) (`q2pep`) then computes
  FDR + q-values from TDC scores and fits an I-spline q -> PEP, separately
  at the PSM and peptide levels.
* **`mprophet`** — a Python re-implementation of the mProphet-based rescoring
  and statistical validation components built into EncyclopeDIA for PRM.
  q-values are Benjamini–Hochberg-adjusted p-values from a Gaussian decoy null,
  and PEPs are Storey's local FDR (`qvalue::lfdr`, KDE on logit-p-values, π₀
  floored at 0.05), exactly as EncyclopeDIA reports them. pyIsoPEP is **not**
  used on this pathway.

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

**Input:** two Percolator-compatible tsv files for background and reference PSMs/peptides.
They must share an **identical header** (same columns, same order).

**Output:** for `context run --background BG --reference REF --prefix P --outdir results/ --engine E`:

```
results/
  weights/P.weights.txt           # trained weights
  P.rescored_features.tsv         # reference features rescored with the trained weights
  P.psm.reference.txt             # PSM-level reference targets
  P.peptide.reference.txt         # peptide-level reference targets
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
| `--background FILE` | *required* | background feature TSV used to train the engine |
| `--reference FILE` | *required* | reference panel feature TSV |
| `--prefix STR` | *required* | output prefix for results files |
| `--outdir DIR` | `results` | output directory |
| `--engine NAME` | `percolator` | `percolator` or `mprophet` |
| `--seed INT` | `1` | seed; passed to Percolator or to mprophet's RNG |
| `--container-cmd CMD` | `podman` | container runtime fallback (`podman` or `docker`) |
| `--input-profile NAME` | `encyclopedia` | mprophet-only; feature-column selection profile (`auto`, `pin`, `encyclopedia`) |
| `--seed-coefficients NAME_OR_PATH` | `encyclopedia` | mprophet-only; built-in name (`encyclopedia`, `none`) or path to a JSON file mapping feature names to seed-model coefficients |
| `--weights-out FILE` | `<prefix>.weights.txt` | weights output file name (or path); relative paths land under `<outdir>/weights` |
| `--rescored-out FILE` | `<prefix>.rescored_features.tsv` | rescored-features output file name (or path); relative paths land under `<outdir>` |
| `--psm-out FILE` | `<prefix>.psm.reference.txt` | PSM-level output file name (or path); relative paths land under `<outdir>` |
| `--peptide-out FILE` | `<prefix>.peptide.reference.txt` | peptide-level output file name (or path); relative paths land under `<outdir>` |

### Input profiles (mprophet only)

* `auto` — if any column starts with `var_` or `main_var_`, keep only
  those (OpenSWATH/pyprophet convention); otherwise fall back to
  `encyclopedia`.
* `pin` — keep every column the pin convention exposes as a feature.
* `encyclopedia` — drop the metadata columns EncyclopeDIA's
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

### Confidence estimation

The two engines use different back-ends for q-values and PEPs:

* **`percolator`** — Percolator trains the discriminant, then pyIsoPEP
  computes q-values (from TDC counts) and PEPs (from an I-spline fit to
  q(score)). Called separately for the PSM-level and the peptide-level
  tables.
* **`mprophet`** — EncyclopeDIA's method throughout, at both
  the PSM level and (after per-sequence deduplication, highest score
  wins) the peptide level:
    * fit a Gaussian N(µ_d, σ_d) to the decoy scores;
    * `p = 1 - Φ((score - µ_d) / σ_d)` for both targets and decoys;
    * `q` = Benjamini–Hochberg-adjusted p-values;
    * `PEP` = Storey `qvalue::lfdr` (KDE on the logit-transformed
      p-values, Silverman bandwidth, monotone non-decreasing in p),
      with π₀ floored at 0.05.

  The same routine is used inside the mProphet training loop for the
  "passing targets" selection at each inner iteration and for the
  held-out evaluation that picks between the trained LDA and the seed
  model.

---

## Examples
```bash
cd example
```
### PyPI

```bash
# Percolator (default)
context run \
  --background background.tsv \
  --reference    reference.tsv \
  --prefix    run01 \
  --outdir    results_run01

# mProphet
context run \
  --background background.tsv \
  --reference    reference.tsv \
  --prefix    run01 \
  --outdir    results_run01 \
  --engine    mprophet
```

### Container image

```bash
# Podman
podman run --rm -v "$PWD:/work" -w /work \
  ghcr.io/shannon225/context:main \
  run --background background.tsv --reference reference.tsv \
      --prefix run01 --outdir results_run01 --engine mprophet

# Apptainer
apptainer run --bind "$PWD:/work" --pwd /work context.sif \
  run --background background.tsv --reference reference.tsv \
      --prefix run01 --outdir results_run01
```

---

## Links
* **PyPI package:** <https://pypi.org/project/context-ms>
* **Container image:** <https://github.com/shannon225/Context/pkgs/container/context>
* **GitHub repository:** <https://github.com/shannon225/Context>
* **Percolator:** <http://percolator.ms>
* **EncyclopeDIA:** <https://bitbucket.org/searleb/encyclopedia>
* **pyIsoPEP:** <https://github.com/statisticalbiotechnology/smooth_q_to_pep>
