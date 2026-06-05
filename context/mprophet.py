import json
from pathlib import Path
import numpy as np
import pandas as pd

from . import fdr as fdr_mod
from .io import numeric_features

SEED_DIR = Path(__file__).parent / "seeds"
BUILTIN_SEEDS = ("encyclopedia", "none")

OUTER_ITERS = 50
INNER_ITERS = 10
MAX_KEPT = 25
DECOY_CAP_RATIO = 10
RANGE_THRESHOLD = 1e-4
INITIAL_TOP_FRACTION = 0.15
INNER_TARGET_FDR = (None, 0.02, 0.01)  # iter 0 uses percentile; rest 1%


def load_seed_coefficients(name_or_path):
    if name_or_path is None or name_or_path == "none":
        return {}
    candidate = SEED_DIR / f"{name_or_path}.json"
    if candidate.is_file():
        path = candidate
    else:
        path = Path(name_or_path)
        if not path.is_file():
            raise FileNotFoundError(
                f"seed coefficients not found: {name_or_path!r} "
                f"(tried built-in {candidate} and {path})"
            )
    with path.open() as f:
        raw = json.load(f)
    return {k: float(v) for k, v in raw.items() if not k.startswith("_")}


def _seed_lda(feature_cols, seed_dict):
    coeffs = np.array([float(seed_dict.get(c, 0.0)) for c in feature_cols])
    if not np.any(coeffs):
        return None
    return (coeffs, 0.0)


def _is_usable(lda):
    if lda is None:
        return False
    w, b = lda
    if not (np.all(np.isfinite(w)) and np.isfinite(b)):
        return False
    return bool(np.any(w))


def _score(X, lda):
    w, b = lda
    return X @ w + b


def build_lda(positive, negative):
    pos = np.asarray(positive, dtype=float)
    neg = np.asarray(negative, dtype=float)
    if pos.ndim != 2 or neg.ndim != 2 or pos.shape[1] != neg.shape[1]:
        return None
    n_pos, n_neg = pos.shape[0], neg.shape[0]
    total = n_pos + n_neg
    if n_pos == 0 or n_neg == 0:
        return None
    pos_prior = n_pos / total
    neg_prior = n_neg / total

    p_range = pos.max(axis=0) - pos.min(axis=0)
    n_range = neg.max(axis=0) - neg.min(axis=0)
    use = (p_range > RANGE_THRESHOLD) & (n_range > RANGE_THRESHOLD)
    if not np.any(use):
        return None

    pos_u = pos[:, use]
    neg_u = neg[:, use]
    mean_pos = pos_u.mean(axis=0)
    mean_neg = neg_u.mean(axis=0)
    mean_all = mean_pos * pos_prior + mean_neg * neg_prior
    mean_sum = mean_pos + mean_neg
    mean_diff = mean_pos - mean_neg

    pos_c = pos_u - mean_all
    neg_c = neg_u - mean_all
    cov_pos = (pos_c.T @ pos_c) / pos_c.shape[0]
    cov_neg = (neg_c.T @ neg_c) / neg_c.shape[0]
    pooled = cov_pos * pos_prior + cov_neg * neg_prior

    try:
        inv = np.linalg.inv(pooled)
    except np.linalg.LinAlgError:
        return None

    coefs_used = inv @ mean_diff
    if not np.all(np.isfinite(coefs_used)):
        return None
    zero_pt = float(coefs_used @ mean_sum)
    constant = -np.log(neg_prior / pos_prior) - 0.5 * zero_pt

    coeffs = np.zeros(pos.shape[1])
    coeffs[use] = coefs_used
    return (coeffs, float(constant))


def average_ldas(ldas):
    if not ldas:
        return None
    w = np.mean([lda[0] for lda in ldas], axis=0)
    b = float(np.mean([lda[1] for lda in ldas]))
    return (w, b)


def _split_kfold(n, k, rng, max_fold_size):
    perm = rng.permutation(n)
    folds = [[] for _ in range(k)]
    for i, idx in enumerate(perm):
        slot = i % k
        if len(folds[slot]) < max_fold_size:
            folds[slot].append(int(idx))
    return [np.array(f, dtype=int) for f in folds]


def _passing_by_percentile(scores, frac):
    if len(scores) == 0 or frac <= 0:
        return np.array([], dtype=int)
    thr = np.quantile(scores, 1.0 - frac)
    return np.where(scores >= thr)[0]


def _passing_by_fdr(target_scores, decoy_scores, target_fdr, rng):
    if len(target_scores) == 0 or len(decoy_scores) == 0:
        return np.array([], dtype=int)
    q, _ = fdr_mod.qvalues_from_scores(target_scores, decoy_scores, rng=rng)
    return np.where(q < target_fdr)[0]


def _auto_main_score(target_X, decoy_X):
    t_range = target_X.max(axis=0) - target_X.min(axis=0)
    d_range = decoy_X.max(axis=0) - decoy_X.min(axis=0)
    valid = (t_range > 1e-6) & (d_range > 1e-6)
    if not np.any(valid):
        return 0
    diffs = np.abs(target_X.mean(axis=0) - decoy_X.mean(axis=0))
    pooled = (target_X.std(axis=0) + decoy_X.std(axis=0)) / 2.0
    score = np.where(valid, diffs / np.maximum(pooled, 1e-12), -1.0)
    return int(np.argmax(score))


def _initial_score(train_t, train_d, lda):
    if lda is not None:
        return _score(train_t, lda)
    col = _auto_main_score(train_t, train_d)
    return train_t[:, col]


def _inner_loop(train_t, train_d, seed_lda, *, inner_iters, decoy_cap_ratio,
                 rng):
    lda = seed_lda if _is_usable(seed_lda) else None
    prev = None
    best = 0
    for i in range(inner_iters):
        if i == 0:
            scores_t = _initial_score(train_t, train_d, lda)
            passing_idx = _passing_by_percentile(scores_t, INITIAL_TOP_FRACTION)
        else:
            if lda is None:
                break
            scores_t = _score(train_t, lda)
            scores_d = _score(train_d, lda)
            target_fdr = INNER_TARGET_FDR[i] if i < len(INNER_TARGET_FDR) else 0.01
            passing_idx = _passing_by_fdr(scores_t, scores_d, target_fdr, rng)

        if len(passing_idx) == 0 or (i > 2 and len(passing_idx) <= best):
            if prev is not None:
                lda = prev
            break
        best = len(passing_idx)

        pos = train_t[passing_idx]
        n_neg_keep = min(len(train_d), max(1, len(pos) * decoy_cap_ratio))
        neg = train_d[:n_neg_keep]
        new_lda = build_lda(pos, neg)
        if not _is_usable(new_lda):
            break
        prev = lda
        lda = new_lda
    return lda


def train(df, feature_cols, label_col, *, seed,
          seed_coefficients=None,
          outer_iters=OUTER_ITERS,
          inner_iters=INNER_ITERS,
          max_kept=MAX_KEPT,
          max_fold_size=None,
          decoy_cap_ratio=DECOY_CAP_RATIO):
    feature_cols = list(feature_cols)
    X = numeric_features(df, feature_cols).to_numpy(dtype=float)
    labels = pd.to_numeric(df[label_col], errors="coerce").fillna(0).astype(int).to_numpy()
    target_mask = labels == 1
    decoy_mask = labels == -1
    X_t = X[target_mask]
    X_d = X[decoy_mask]
    n_t, n_d = len(X_t), len(X_d)
    if n_t == 0 or n_d == 0:
        raise ValueError(
            f"mProphet needs both targets and decoys; got {n_t} targets, {n_d} decoys"
        )

    seed_dict = seed_coefficients or {}
    seed_lda = _seed_lda(feature_cols, seed_dict)
    fold_cap = max_fold_size if max_fold_size is not None else max(n_t, n_d)

    rng = np.random.default_rng(seed)
    fdr_rng = np.random.default_rng(seed + 1)
    kept = []  # (passing_count, lda)

    for it in range(outer_iters):
        t_folds = _split_kfold(n_t, 2, rng, fold_cap)
        d_folds = _split_kfold(n_d, 2, rng, fold_cap)
        train_t, test_t = X_t[t_folds[0]], X_t[t_folds[1]]
        train_d, test_d = X_d[d_folds[0]], X_d[d_folds[1]]
        if len(train_t) == 0 or len(train_d) == 0:
            continue

        lda = _inner_loop(train_t, train_d, seed_lda,
                          inner_iters=inner_iters,
                          decoy_cap_ratio=decoy_cap_ratio,
                          rng=fdr_rng)
        if lda is None:
            continue

        lda_pass = _passing_by_fdr(_score(test_t, lda),
                                    _score(test_d, lda), 0.01, fdr_rng)
        if _is_usable(seed_lda):
            seed_pass = _passing_by_fdr(_score(test_t, seed_lda),
                                         _score(test_d, seed_lda), 0.01, fdr_rng)
            if len(seed_pass) > len(lda_pass):
                kept.append((len(seed_pass), seed_lda))
                continue
        kept.append((len(lda_pass), lda))

    if not kept:
        if seed_lda is None:
            raise RuntimeError(
                "mProphet training failed: no usable models produced and no "
                "seed model available. Try a different seed-coefficients "
                "profile or check that targets/decoys are sufficient."
            )
        print("[mprophet] no trained model beat seed; falling back to seed",
              flush=True)
        avg = seed_lda
    else:
        kept.sort(key=lambda x: x[0], reverse=True)
        top = [lda for _, lda in kept[:max_kept]]
        avg = average_ldas(top)
        print(f"[mprophet] kept top {len(top)}/{len(kept)} models for averaging",
              flush=True)

    w, b = avg
    return feature_cols, w, b


def write_weights(out_path, feature_cols, w, b):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "feature": list(feature_cols) + ["__bias__"],
        "weight": list(map(float, w)) + [float(b)],
    }).to_csv(out_path, sep="\t", index=False)


def parse_weights(path):
    df = pd.read_csv(path, sep="\t")
    feat = [str(x) for x in df["feature"].tolist()]
    w = df["weight"].to_numpy(dtype=float)
    bias_mask = [name == "__bias__" for name in feat]
    if any(bias_mask):
        bi = bias_mask.index(True)
        b = float(w[bi])
        feat = [f for f, m in zip(feat, bias_mask) if not m]
        w = w[[not m for m in bias_mask]]
    else:
        b = 0.0
    return feat, w[np.newaxis, :], np.array([b])
