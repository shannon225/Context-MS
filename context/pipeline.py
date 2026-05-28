from pathlib import Path
import numpy as np
from .io import (
    FeatureFile, OUT_SCORE, check_headers_match, coerce_label,
    numeric_features, psm_to_peptide, read_features_raw, to_final_output,
)
from .percolator import parse_weights, run_percolator
from .pyisopep import annotate_q_and_pep


def score_table(df, feature_cols, w, b):
    X = numeric_features(df, feature_cols).to_numpy(dtype=float)
    return X @ np.asarray(w, dtype=float) + float(b)


def score_and_calibrate(ff, w, b, *, container_cmd):
    rescored = ff.df.copy()
    rescored[OUT_SCORE] = score_table(rescored, ff.feature_cols, w, b)

    psm_tgt = annotate_q_and_pep(
        rescored, score_col=OUT_SCORE, label_col=ff.label_col,
        container_cmd=container_cmd,
    )
    psm_out = to_final_output(psm_tgt, ff)

    pep_input = psm_to_peptide(rescored, peptide_col=ff.peptide_col,
                                score_col=OUT_SCORE)
    pep_tgt = annotate_q_and_pep(
        pep_input, score_col=OUT_SCORE, label_col=ff.label_col,
        container_cmd=container_cmd,
    )
    pep_out = to_final_output(pep_tgt, ff)
    return rescored, psm_out, pep_out


def run(nontarget_path, target_path, *, outdir, prefix, seed,
            container_cmd):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    weights_dir = outdir / "weights"
    weights_dir.mkdir(parents=True, exist_ok=True)

    check_headers_match(nontarget_path, target_path)

    wpath = weights_dir / f"{prefix}.seed{seed}.weights.txt"
    print(f"[percolator] {prefix} seed={seed}", flush=True)
    run_percolator(nontarget_path, wpath, seed=seed, container_cmd=container_cmd)

    feature_names, w_raw, b_raw = parse_weights(wpath)
    w = w_raw.mean(axis=0)
    b = float(b_raw.mean())

    df = read_features_raw(target_path)
    ff = FeatureFile(df, feature_cols=feature_names)
    coerce_label(ff.df, ff.label_col)

    rescored, psm_out, pep_out = score_and_calibrate(
        ff, w, b, container_cmd=container_cmd,
    )

    rescored_path = outdir / f"{prefix}.seed{seed}.rescored_features.tsv"
    psm_path = outdir / f"{prefix}.seed{seed}.psm.target.txt"
    pep_path = outdir / f"{prefix}.seed{seed}.peptide.target.txt"
    rescored.to_csv(rescored_path, sep="\t", index=False)
    psm_out.to_csv(psm_path, sep="\t", index=False)
    pep_out.to_csv(pep_path, sep="\t", index=False)

    print(f"[write] {rescored_path.name}  ({len(rescored)} rows including decoys)")
    print(f"[write] {psm_path.name}        ({len(psm_out)} target PSMs)")
    print(f"[write] {pep_path.name}    ({len(pep_out)} target peptides)")

    return {
        "weights": wpath,
        "rescored": rescored_path,
        "psm": psm_path,
        "peptide": pep_path,
    }
