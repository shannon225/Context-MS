from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
import os
from .io import numeric_features

PYPROPHET_IMAGE = "ghcr.io/pyprophet/pyprophet:latest"
SIGMA_DROP_THRESHOLD = 1e-6

def train(nontarget_path, weights_out, seed, *, feature_cols, label_col, id_col,
          container_cmd="podman", image=PYPROPHET_IMAGE, classifier="LDA",
          ss_num_iter=10):
    if classifier not in ("LDA", "SVM"):
        raise ValueError(
            f"pyprophet engine only supports LDA/SVM (linear); got {classifier!r}"
        )

    nontarget_path = Path(nontarget_path).resolve()
    weights_out = Path(weights_out).resolve()
    weights_out.parent.mkdir(parents=True, exist_ok=True)

    nt_df = pd.read_csv(nontarget_path, sep="\t")
    Xnt = numeric_features(nt_df, feature_cols).to_numpy(dtype=float)
    mu = Xnt.mean(axis=0)
    sigma = Xnt.std(axis=0, ddof=0)

    bad = np.where(sigma < SIGMA_DROP_THRESHOLD)[0]
    if bad.size:
        names = [feature_cols[i] for i in bad]
        raise RuntimeError(
            f"feature(s) {names} have ~zero variance in nontarget after "
            "pruning; this should have been dropped upstream."
        )

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        pp_input = tmp / "nontarget_pp.tsv"
        _write_pyprophet_tsv(
            nt_df, feature_cols, label_col, id_col,
            mu=mu, sigma=sigma, out_path=pp_input,
        )

        _invoke_pyprophet(
            pp_input, seed=seed,
            classifier=classifier, ss_num_iter=ss_num_iter,
            container_cmd=container_cmd, image=image,
        )

        weights_csv = pp_input.with_name(pp_input.stem + "_weights.csv")
        if not weights_csv.exists():
            sys.exit(f"pyprophet weights file not produced: {weights_csv}")

        w_std = _read_pyprophet_weights(weights_csv, feature_cols)

    w_raw = w_std / sigma
    _write_raw_weights(weights_out, feature_cols, w_raw, mu=mu, sigma=sigma)
    return feature_cols, w_raw, 0.0


def _write_pyprophet_tsv(df, feature_cols, label_col, id_col, *,
                          mu, sigma, out_path):
    """Convert pin to pyprophet-compatible tsv:
      id_col: transition_group_id
      label_col: decoy (Label == 1: decoy=0, Label == -1: decoy=1)
      first feature: main_var_<name>
      remaining feature: var_<name>
    """
    label_in = pd.to_numeric(df[label_col], errors="coerce").fillna(0).astype(int)
    decoy = (label_in != 1).astype(int)

    Xraw = numeric_features(df, feature_cols).to_numpy(dtype=float)
    Xstd = (Xraw - mu) / sigma

    out = pd.DataFrame({
        "transition_group_id": df[id_col].astype(str).values,
        "run_id": "run0",
        "group_id": df[id_col].astype(str).values,
        "decoy": decoy.values,
    })
    rename = {}
    for i, name in enumerate(feature_cols):
        col = f"main_var_{name}" if i == 0 else f"var_{name}"
        rename[name] = col
        out[col] = Xstd[:, i]
    out.to_csv(out_path, sep="\t", index=False)
    return rename


def _invoke_pyprophet(pp_input, *, seed, classifier, ss_num_iter,
                     container_cmd, image):
    env = {**os.environ, "PYTHONHASHSEED": str(seed)}

    if shutil.which("pyprophet"):
        cmd = [
            "pyprophet", "score",
            "--in", str(pp_input),
            "--classifier", classifier,
            "--level", "ms2",
            "--ss_num_iter", str(ss_num_iter),
        ]
        cwd = pp_input.parent
    else:
        in_dir = pp_input.parent
        cmd = [
            container_cmd, "run", "--rm",
            "-e", f"PYTHONHASHSEED={seed}",
            "-v", f"{in_dir}:/in",
            "-w", "/in",
            image,
            "pyprophet", "score",
            "--in", f"/in/{pp_input.name}",
            "--classifier", classifier,
            "--level", "ms2",
            "--ss_num_iter", str(ss_num_iter),
        ]
        cwd = None

    console = pp_input.with_suffix(".pyprophet.log")
    with console.open("w") as fh:
        rc = subprocess.run(
            cmd, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=cwd,
        ).returncode
    if rc != 0:
        persisted = Path(tempfile.gettempdir()) / f"{pp_input.stem}.pyprophet.log"
        shutil.copy(console, persisted)
        log_tail = console.read_text(errors="replace")
        sys.exit(
            f"pyprophet failed (rc={rc}). Log copied to {persisted}\n"
            f"--- pyprophet output ---\n{log_tail}"
        )


def _read_pyprophet_weights(weights_csv, feature_cols):
    df = pd.read_csv(weights_csv)
    if "level" in df.columns:
        df = df[df["level"] == "ms2"]
    lookup = {}
    for _, row in df.iterrows():
        name = str(row["score"])
        for prefix in ("main_var_", "var_"):
            if name.startswith(prefix):
                lookup[name[len(prefix):]] = float(row["weight"])
                break

    missing = [c for c in feature_cols if c not in lookup]
    if missing:
        raise ValueError(
            f"pyprophet weights missing entries for: {missing} "
            f"(weights file: {weights_csv})"
        )
    return np.array([lookup[c] for c in feature_cols], dtype=float)


def _write_raw_weights(weights_out, feature_cols, w_raw, *, mu, sigma):
    pd.DataFrame({
        "feature": list(feature_cols) + ["__bias__"],
        "weight":  list(w_raw)        + [0.0],
        "nontarget_mu":    list(mu)    + [np.nan],
        "nontarget_sigma": list(sigma) + [np.nan],
    }).to_csv(weights_out, sep="\t", index=False)


def parse_weights(path):
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
