import numpy as np
from scipy.special import erfc

PI0_FLOOR = 0.05
RNG_SEED = 42
NUM_LAMBDA = 100
MAX_LAMBDA = 0.95
DEFAULT_NULL = "gaussian"


def gaussian_pvalues(target_scores, decoy_scores):
    target_scores = np.asarray(target_scores, dtype=float)
    decoy_scores = np.asarray(decoy_scores, dtype=float)
    if decoy_scores.size == 0:
        return np.ones_like(target_scores)
    mu = float(decoy_scores.mean())
    sigma = float(decoy_scores.std(ddof=0))
    if not np.isfinite(sigma) or sigma <= 0:
        return np.where(target_scores > mu, 0.0, 1.0)
    z = (target_scores - mu) / sigma
    p = 0.5 * erfc(z / np.sqrt(2.0))
    return np.clip(p, 0.0, 1.0)


def empirical_pvalues(target_scores, decoy_scores):
    target_scores = np.asarray(target_scores, dtype=float)
    decoy_scores = np.asarray(decoy_scores, dtype=float)
    n_d = decoy_scores.size
    if n_d == 0:
        return np.ones_like(target_scores)
    decoy_sorted = np.sort(decoy_scores)
    idx = np.searchsorted(decoy_sorted, target_scores, side="left")
    p = (n_d - idx + 1.0) / (n_d + 1.0)
    return np.clip(p, 0.0, 1.0)


def pvalues(target_scores, decoy_scores, *, null=DEFAULT_NULL):
    if null == "gaussian":
        return gaussian_pvalues(target_scores, decoy_scores)
    if null == "empirical":
        return empirical_pvalues(target_scores, decoy_scores)
    raise ValueError(f"unknown null {null!r}; choose 'gaussian' or 'empirical'")


def estimate_pi0(p, n_boot=100, rng=None):
    rng = np.random.default_rng(RNG_SEED) if rng is None else rng
    p = np.sort(np.asarray(p, dtype=float))
    m = len(p)
    if m == 0:
        return 1.0
    lams, pi0s = [], []
    for k in range(1, NUM_LAMBDA + 1):
        l = k / NUM_LAMBDA * MAX_LAMBDA
        Wl = m - np.searchsorted(p, l, side="right")
        v = Wl / ((1.0 - l) * m)
        if v > 0:
            lams.append(l); pi0s.append(v)
    if not pi0s:
        return 1.0
    pi0s = np.asarray(pi0s)
    min_pi0 = pi0s.min()
    mse = np.zeros_like(pi0s)
    for _ in range(n_boot):
        b = np.sort(rng.choice(p, m, replace=True))
        nb = len(b)
        for i, l in enumerate(lams):
            Wb = nb - np.searchsorted(b, l, side="right")
            mse[i] += (Wb / ((1.0 - l) * nb) - min_pi0) ** 2
    return float(np.clip(pi0s[int(np.argmin(mse))], 0.0, 1.0))


def storey_qvalue(p, pi0):
    p = np.asarray(p, dtype=float)
    m = p.size
    if m == 0:
        return np.array([])
    order = np.argsort(p, kind="mergesort")
    ps = p[order]
    ranks = np.arange(1, m + 1)
    fdr_s = pi0 * m * ps / ranks
    q_s = np.minimum.accumulate(fdr_s[::-1])[::-1]
    q = np.empty_like(q_s)
    q[order] = q_s
    return q


def qvalues_from_scores(target_scores, decoy_scores, *, rng=None,
                          null=DEFAULT_NULL):
    p = pvalues(target_scores, decoy_scores, null=null)
    pi0 = max(PI0_FLOOR, estimate_pi0(p, rng=rng))
    q = storey_qvalue(p, pi0)
    return q, pi0
