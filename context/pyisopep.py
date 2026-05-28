from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd

PYISOPEP_IMAGE = "ghcr.io/statisticalbiotechnology/pyisopep:main"
Q_COL = "q-value"
PEP_COL = "posterior_error_prob"


def annotate_q_and_pep(df, *, score_col="score", label_col="Label",
                       container_cmd="podman", image=PYISOPEP_IMAGE):
    try:
        return _via_python_api(df, score_col=score_col, label_col=label_col)
    except ImportError:
        return _via_container(df, score_col=score_col, label_col=label_col,
                              container_cmd=container_cmd, image=image)


def _via_python_api(df, *, score_col, label_col):
    from pyIsoPEP.IsotonicPEP import IsotonicPEP
    lab = df[label_col].astype(int).to_numpy()
    lab_iso = np.where(lab == 1, 0, 1)
    obs = np.column_stack([df[score_col].to_numpy(dtype=float), lab_iso])
    _, q_arr, pep_arr, _ = IsotonicPEP().pep_regression(
        obs=obs, method="q2pep", calc_q_from_fdr=True, calc_q_from_pep=False,
    )
    tgt = df[lab == 1].copy().reset_index(drop=True)
    tgt[Q_COL] = q_arr
    tgt[PEP_COL] = pep_arr
    return tgt.sort_values(Q_COL, kind="mergesort").reset_index(drop=True)


def _via_container(df, *, score_col, label_col, container_cmd, image):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        in_tsv = tmp / "in.tsv"
        out_tsv = tmp / "out.tsv"
        df.to_csv(in_tsv, sep="\t", index=False)
        cmd = [
            container_cmd, "run", "--rm",
            "--entrypoint", "pyisopep",
            "-v", f"{tmp}:/{tmp}",
            image,
            "q2pep",
            "--cat-file", f"/{tmp}/in.tsv",
            "--score-col", score_col,
            "--label-col", label_col,
            "--target-label", "1",
            "--decoy-label", "-1",
            "--calc-q-from-fdr",
            "--output", f"/{tmp}/out.tsv",
        ]
        rc = subprocess.run(cmd).returncode
        if rc != 0:
            sys.exit(f"pyIsoPEP container failed (rc={rc})")
        out = pd.read_csv(out_tsv, sep="\t")
        out = out.rename(columns={
            "pyIsoPEP q-value from FDR": Q_COL,
            "pyIsoPEP PEP": PEP_COL,
        })
        if "pyIsoPEP FDR" in out.columns:
            out = out.drop(columns=["pyIsoPEP FDR"])
        return out.sort_values(Q_COL, kind="mergesort").reset_index(drop=True)
