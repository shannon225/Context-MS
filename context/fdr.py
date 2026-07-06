"""FDR / q-value / PEP calculations for the mprophet engine,
matching Encyclopedia's method exactly.
"""
import numpy as np
from scipy.special import erfc
from scipy.stats import gaussian_kde

PI0_FLOOR = 0.05
SMALLEST_MEANINGFUL_P = 1e-8


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


def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float)
    m = p.size
    if m == 0:
        return np.array([])
    order = np.argsort(p, kind="mergesort")
    ps = p[order]
    ranks = np.arange(1, m + 1, dtype=float)
    adj = ps * m / ranks
    q_s = np.minimum.accumulate(adj[::-1])[::-1]
    q_s = np.clip(q_s, 0.0, 1.0)
    q = np.empty_like(q_s)
    q[order] = q_s
    return q


def estimate_pi0(p):
    p = np.asarray(p, dtype=float)
    m = p.size
    if m == 0:
        return 1.0
    lambdas = np.arange(1, 10) / 10.0
    pi0s = np.array([np.mean(p > lam) / (1.0 - lam) for lam in lambdas])
    return float(np.median(pi0s))


def local_fdr(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return np.array([])
    pi0 = max(PI0_FLOOR, estimate_pi0(p))

    eps = SMALLEST_MEANINGFUL_P
    x = np.log((p + eps) / (1.0 - p + eps))
    ex = np.exp(x)
    dx = ex / (1.0 + ex) ** 2

    try:
        density = gaussian_kde(x, bw_method="silverman")
        f = density(x)
    except np.linalg.LinAlgError:
        f = np.full(n, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        lfdr = np.where(f > 0, pi0 * dx / f, 1.0)
    lfdr = np.clip(lfdr, 0.0, 1.0)

    order = np.argsort(p, kind="mergesort")
    lfdr_sorted = lfdr[order]
    lfdr_sorted = np.maximum.accumulate(lfdr_sorted)
    out = np.empty(n)
    out[order] = lfdr_sorted
    return out


def qvalues_bh(target_scores, decoy_scores):
    p = gaussian_pvalues(target_scores, decoy_scores)
    return benjamini_hochberg(p)


def q_and_pep(target_scores, decoy_scores):
    p = gaussian_pvalues(target_scores, decoy_scores)
    q = benjamini_hochberg(p)
    pep = local_fdr(p)
    return q, pep
