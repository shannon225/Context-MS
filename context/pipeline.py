from pathlib import Path
import tempfile
import numpy as np
import pandas as pd
from .io import (
    FeatureFile, OUT_SCORE, OUT_Q, OUT_PEP, check_headers_match, coerce_label,
    detect_feature_cols, near_constant_cols, numeric_features,
    psm_to_peptide, read_features_raw, select_training_features,
    to_final_output, write_pruned_tsv,
)
from . import percolator
from . import mprophet as mprophet_engine
from . import fdr as fdr_mod
from .pyisopep import annotate_q_and_pep


ENGINES = ("percolator", "mprophet")
DEFAULT_INPUT_PROFILE = "encyclopedia"
DEFAULT_SEED_COEFFICIENTS = "encyclopedia"


def score_table(df, feature_cols, w, b):
    X = numeric_features(df, feature_cols).to_numpy(dtype=float)
    return X @ np.asarray(w, dtype=float) + float(b)


def _score_all(ff, w, b):
    rescored = ff.df.copy()
    rescored[OUT_SCORE] = score_table(rescored, ff.feature_cols, w, b)
    return rescored


def _split_scores(df, label_col, score_col=OUT_SCORE):
    lab = pd.to_numeric(df[label_col], errors="coerce").fillna(0).astype(int)
    t_mask = (lab == 1).to_numpy()
    d_mask = (lab == -1).to_numpy()
    scores = df[score_col].to_numpy(dtype=float)
    return t_mask, d_mask, scores


def _percolator_q_and_pep(ff, rescored, container_cmd):
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
    return psm_out, pep_out


def _mprophet_q_and_pep_frame(df, label_col):
    t_mask, d_mask, scores = _split_scores(df, label_col)
    q, pep = fdr_mod.q_and_pep(scores[t_mask], scores[d_mask])
    tgt = df.loc[t_mask].copy().reset_index(drop=True)
    tgt[OUT_Q] = q
    tgt[OUT_PEP] = pep
    return tgt.sort_values(OUT_Q, kind="mergesort").reset_index(drop=True)


def _mprophet_q_and_pep(ff, rescored):
    psm_tgt = _mprophet_q_and_pep_frame(rescored, ff.label_col)
    psm_out = to_final_output(psm_tgt, ff)

    pep_input = psm_to_peptide(rescored, peptide_col=ff.peptide_col,
                                score_col=OUT_SCORE)
    pep_tgt = _mprophet_q_and_pep_frame(pep_input, ff.label_col)
    pep_out = to_final_output(pep_tgt, ff)
    return psm_out, pep_out


def _train_percolator(pruned_bg, weights_out, seed, container_cmd):
    percolator.run_percolator(
        pruned_bg, weights_out, seed=seed, container_cmd=container_cmd,
    )
    feature_names, w_raw, b_raw = percolator.parse_weights(weights_out)
    w = w_raw.mean(axis=0)
    b = float(b_raw.mean())
    return feature_names, w, b


def _train_mprophet(pruned_bg_path, weights_out, seed, *,
                     pin_feature_cols, input_profile, seed_coefficients_name):
    bg_df = read_features_raw(pruned_bg_path)
    label_col = pin_feature_cols_label(pruned_bg_path)
    training_cols = select_training_features(
        pin_feature_cols, profile=input_profile,
    )
    if not training_cols:
        raise ValueError(
            f"no training features remain after applying input profile "
            f"{input_profile!r}; check your feature column names"
        )
    seed_coeffs = mprophet_engine.load_seed_coefficients(seed_coefficients_name)
    feature_names, w, b = mprophet_engine.train(
        bg_df, training_cols, label_col, seed=seed,
        seed_coefficients=seed_coeffs,
    )
    mprophet_engine.write_weights(weights_out, feature_names, w, b)
    return feature_names, w, b


def pin_feature_cols_label(path):
    with open(path, encoding="utf-8") as f:
        header = [c.strip() for c in f.readline().rstrip("\n").split("\t")]
    return header[1]


def _resolve_out(name, *, base, default_name):
    if name is None:
        return base / default_name
    p = Path(name)
    if p.is_absolute():
        return p
    return base / p


def run(background_path, reference_path, *, outdir, prefix, seed, engine,
        container_cmd, psm_out=None, peptide_out=None,
        rescored_out=None, weights_out=None,
        input_profile=DEFAULT_INPUT_PROFILE,
        seed_coefficients=DEFAULT_SEED_COEFFICIENTS):
    if engine not in ENGINES:
        raise ValueError(f"unknown engine {engine!r}; choose from {ENGINES}")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    weights_dir = outdir / "weights"
    weights_dir.mkdir(parents=True, exist_ok=True)

    header = check_headers_match(background_path, reference_path)
    feature_cols = detect_feature_cols(header)

    bg_df = read_features_raw(background_path)
    ref_df = read_features_raw(reference_path)
    drop = near_constant_cols(bg_df, feature_cols)
    if drop:
        print(f"[prune] dropping {len(drop)} near-constant feature(s): {drop}",
              flush=True)
    kept = [c for c in feature_cols if c not in set(drop)]

    wpath = _resolve_out(
        weights_out, base=weights_dir,
        default_name=f"{prefix}.weights.txt",
    )
    wpath.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        pruned_bg = write_pruned_tsv(bg_df, drop, tmp / "background.tsv")
        pruned_ref = write_pruned_tsv(ref_df, drop, tmp / "reference.tsv")

        print(f"[{engine}] {prefix} seed={seed}", flush=True)
        if engine == "percolator":
            feature_names, w, b = _train_percolator(
                pruned_bg, wpath, seed, container_cmd,
            )
        else:
            feature_names, w, b = _train_mprophet(
                pruned_bg, wpath, seed,
                pin_feature_cols=kept,
                input_profile=input_profile,
                seed_coefficients_name=seed_coefficients,
            )

        df = read_features_raw(pruned_ref)

    ff = FeatureFile(df, feature_cols=feature_names)
    coerce_label(ff.df, ff.label_col)

    rescored = _score_all(ff, w, b)
    if engine == "percolator":
        psm_df, pep_df = _percolator_q_and_pep(ff, rescored, container_cmd)
    else:
        psm_df, pep_df = _mprophet_q_and_pep(ff, rescored)

    rescored_path = _resolve_out(
        rescored_out, base=outdir,
        default_name=f"{prefix}.rescored_features.tsv",
    )
    psm_path = _resolve_out(
        psm_out, base=outdir,
        default_name=f"{prefix}.psm.reference.txt",
    )
    pep_path = _resolve_out(
        peptide_out, base=outdir,
        default_name=f"{prefix}.peptide.reference.txt",
    )
    for p in (rescored_path, psm_path, pep_path):
        p.parent.mkdir(parents=True, exist_ok=True)
    rescored.to_csv(rescored_path, sep="\t", index=False)
    psm_df.to_csv(psm_path, sep="\t", index=False)
    pep_df.to_csv(pep_path, sep="\t", index=False)

    print(f"[write] {rescored_path}  ({len(rescored)} rows including decoys)")
    print(f"[write] {psm_path}        ({len(psm_df)} reference target PSMs)")
    print(f"[write] {pep_path}    ({len(pep_df)} reference target peptides)")

    return {
        "engine": engine,
        "weights": wpath,
        "rescored": rescored_path,
        "psm": psm_path,
        "peptide": pep_path,
    }
