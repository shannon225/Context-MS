"""
Context identifies the non-feature columns by position, following the Percolator convention:
    col[0]    id
    col[1]    label    (target = +1, decoy = -1)
    col[-2]   peptide
    col[-1]   proteins (multiple proteins are joined by comma)
ScanNr is looked up by name.

The mprophet engine accepts an additional --input-profile that further
filters the feature set by name. Percolator is unaffected.
"""
import pandas as pd
import numpy as np

OUT_SCORE = "score"
OUT_Q = "q-value"
OUT_PEP = "posterior_error_prob"

SIGMA_DROP_THRESHOLD = 1e-6

ENCYCLOPEDIA_METADATA_COLS = frozenset({
    "TD", "pepLength",
    "charge1", "charge2", "charge3", "charge4",
    "precursorMass", "RTinMin", "midTime",
    "numberOfMatchingPeaksAboveThreshold", "primary",
})

INPUT_PROFILES = ("auto", "pin", "encyclopedia")


def read_features_raw(path, *, multi_protein_join=","):
    with open(path, encoding="utf-8") as f:
        header = [c.strip() for c in f.readline().rstrip("\n").split("\t")]
        ncols = len(header)
        rows = []
        for raw in f:
            if not raw.strip():
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) > ncols:
                merged = multi_protein_join.join(parts[ncols - 1:])
                parts = parts[:ncols - 1] + [merged]
            if len(parts) != ncols:
                raise ValueError(
                    f"{path}: row has {len(parts)} fields, expected {ncols}"
                )
            rows.append(parts)
    return pd.DataFrame(rows, columns=header)


def read_header(path):
    with open(path, encoding="utf-8") as f:
        return [c.strip() for c in f.readline().rstrip("\n").split("\t")]


def check_headers_match(nontarget_path, target_path):
    nt = read_header(nontarget_path)
    tg = read_header(target_path)
    if nt != tg:
        raise ValueError(
            f"nontarget and target headers do not match.\n"
            f"  nontarget: {nontarget_path}\n  target:    {target_path}\n"
        )
    return nt


class FeatureFile:
    def __init__(self, df, feature_cols):
        self.df = df
        cols = list(df.columns)
        missing = [c for c in feature_cols if c not in cols]
        if missing:
            raise ValueError(
                f"Nontarget features not present in target header: {missing}"
            )
        self.feature_cols = list(feature_cols)
        self.id_col = cols[0]
        self.label_col = cols[1]
        self.peptide_col = cols[-2]
        self.proteins_col = cols[-1]
        lc = {c.lower(): c for c in cols}
        self.scan_col = lc.get("scannr") or lc.get("scan")


def coerce_label(df, label_col):
    df[label_col] = (
        pd.to_numeric(df[label_col].astype(str).str.strip(), errors="coerce")
          .fillna(0).astype(int)
    )
    return df


def numeric_features(df, feature_cols):
    return df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)


def detect_feature_cols(header):
    id_col, label_col = header[0], header[1]
    peptide_col, proteins_col = header[-2], header[-1]
    lc = {c.lower(): c for c in header}
    scan_col = lc.get("scannr") or lc.get("scan")
    non_feat = {id_col, label_col, peptide_col, proteins_col}
    if scan_col is not None:
        non_feat.add(scan_col)
    return [c for c in header if c not in non_feat]


def near_constant_cols(df, feature_cols, *, threshold=SIGMA_DROP_THRESHOLD):
    num = numeric_features(df, feature_cols)
    sigma = num.std(ddof=0)
    return [c for c in feature_cols if sigma[c] < threshold]


def select_training_features(feature_cols, *, profile="auto"):
    if profile not in INPUT_PROFILES:
        raise ValueError(
            f"unknown input profile {profile!r}; choose from {INPUT_PROFILES}"
        )
    if profile == "auto":
        if any(c.startswith(("var_", "main_var_")) for c in feature_cols):
            return [c for c in feature_cols
                    if c.startswith(("var_", "main_var_"))]
        profile = "encyclopedia"
    if profile == "pin":
        return list(feature_cols)
    return [c for c in feature_cols if c not in ENCYCLOPEDIA_METADATA_COLS]


def psm_to_peptide(df, peptide_col, score_col=OUT_SCORE):
    return (
        df.sort_values(score_col, ascending=False, kind="mergesort")
          .drop_duplicates(subset=peptide_col, keep="first")
          .reset_index(drop=True)
    )


def to_final_output(target_df, ff):
    scan_name = ff.scan_col
    scan_vals = (target_df[scan_name].values
                 if scan_name is not None else [""] * len(target_df))
    scan_out_name = scan_name if scan_name is not None else "ScanNr"
    out = pd.DataFrame({
        ff.id_col: target_df[ff.id_col].values,
        scan_out_name: scan_vals,
        OUT_SCORE: target_df[OUT_SCORE].values,
        OUT_Q: target_df[OUT_Q].values,
        OUT_PEP: target_df[OUT_PEP].values,
        ff.peptide_col: target_df[ff.peptide_col].values,
        ff.proteins_col: target_df[ff.proteins_col].values,
    })
    return out


def write_pruned_tsv(df, drop_cols, out_path):
    keep = [c for c in df.columns if c not in set(drop_cols)]
    df[keep].to_csv(out_path, sep="\t", index=False)
    return out_path
