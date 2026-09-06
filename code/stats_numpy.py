"""
NumPy-only statistics for the Cyno proteome half-life pipeline.

The original pipeline deliberately avoids SciPy so that results are reproducible
from NumPy alone (SciPy is not installed in the Cowork sandbox). This module
re-implements the handful of tests the pipeline needs:

  - rankdata (average ranks, tie-aware)
  - spearman  (rho + p via t-approximation)
  - mannwhitneyu (U + normal-approx p with tie correction)
  - cliffs_delta
  - ols  (multiple linear regression: beta, SE, t, p, R^2)
  - regularized incomplete beta -> Student-t two-sided p-value

All permutation-based callers must seed NumPy with default_rng(42).
"""

import numpy as np


# --------------------------------------------------------------------------
# Student-t two-sided p-value via the regularized incomplete beta function.
# (No SciPy: continued-fraction evaluation, Numerical-Recipes style.)
# --------------------------------------------------------------------------
def _betacf(a, b, x, itmax=300, eps=1e-12):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betainc(a, b, x):
    """Regularized incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    from math import lgamma, log, exp
    bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b)
             + a * log(x) + b * log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_sf_two_sided(t, df):
    """Two-sided Student-t p-value for statistic t with df degrees of freedom."""
    t = abs(float(t))
    if df <= 0:
        return float("nan")
    x = df / (df + t * t)
    return betainc(df / 2.0, 0.5, x)


def normal_sf_two_sided(z):
    from math import erfc, sqrt
    return erfc(abs(float(z)) / sqrt(2.0))


# --------------------------------------------------------------------------
# Ranks / Spearman / Mann-Whitney / Cliff's delta
# --------------------------------------------------------------------------
def rankdata(a):
    """Average ranks (ties shared), 1-based."""
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)
    sa = a[order]
    i = 0
    n = len(a)
    while i < n:
        j = i
        while j + 1 < n and sa[j + 1] == sa[i]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        ranks[order[i:j + 1]] = avg
        i = j + 1
    return ranks


def spearman(x, y):
    """Spearman rho and two-sided p (t-approximation, df = n-2)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    rx, ry = rankdata(x), rankdata(y)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    rho = float((rx * ry).sum() / denom) if denom > 0 else float("nan")
    if n > 2 and abs(rho) < 1.0:
        t = rho * np.sqrt((n - 2) / (1.0 - rho * rho))
        p = t_sf_two_sided(t, n - 2)
    else:
        p = 0.0 if abs(rho) == 1.0 else float("nan")
    return rho, p


def mannwhitneyu(x, y):
    """
    Two-sided Mann-Whitney U with normal approximation and tie correction.
    Returns (U for sample x, p). U here is U1 (rank-sum based on x).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n1, n2 = len(x), len(y)
    allv = np.concatenate([x, y])
    ranks = rankdata(allv)
    r1 = ranks[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    # tie correction
    _, counts = np.unique(allv, return_counts=True)
    tie = (counts ** 3 - counts).sum()
    N = n1 + n2
    sigma2 = (n1 * n2 / 12.0) * ((N + 1) - tie / (N * (N - 1)))
    sigma = np.sqrt(sigma2) if sigma2 > 0 else 0.0
    if sigma == 0:
        return float(u1), float("nan")
    z = (u1 - mu) / sigma  # no continuity correction (matches pipeline)
    p = normal_sf_two_sided(z)
    return float(u1), float(p)


def cliffs_delta(x, y):
    """Cliff's delta = (#x>y - #x<y) / (nx*ny), via rank-sum identity."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n1, n2 = len(x), len(y)
    ranks = rankdata(np.concatenate([x, y]))
    r1 = ranks[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2.0
    return float(2.0 * u1 / (n1 * n2) - 1.0)


# --------------------------------------------------------------------------
# Ordinary least squares
# --------------------------------------------------------------------------
def ols(X, y):
    """
    X: (n, k) design WITHOUT intercept column (added here).
    Returns dict with beta (incl intercept first), se, t, p, r2, adj_r2.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n = X.shape[0]
    Xd = np.column_stack([np.ones(n), X])
    k = Xd.shape[1]
    XtX = Xd.T @ Xd
    XtX_inv = np.linalg.inv(XtX)
    beta = XtX_inv @ (Xd.T @ y)
    resid = y - Xd @ beta
    ss_res = float(resid @ resid)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot
    df = n - k
    sigma2 = ss_res / df
    se = np.sqrt(np.diag(XtX_inv) * sigma2)
    tvals = beta / se
    pvals = np.array([t_sf_two_sided(t, df) for t in tvals])
    adj_r2 = 1.0 - (1.0 - r2) * (n - 1) / df
    return dict(beta=beta, se=se, t=tvals, p=pvals,
                r2=r2, adj_r2=adj_r2, df=df, n=n)
