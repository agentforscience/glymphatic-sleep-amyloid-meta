#!/usr/bin/env python3
"""
metalib.py - Core meta-analysis primitives for the slow-wave-sleep / amyloid-PET
meta-analysis.

Everything here operates on the Fisher-z scale (y) with known sampling variances (v).
Implemented from primary formulations rather than relying on PyMARE, because we need
(a) Knapp-Hartung-Sidik-Jonkman small-sample standard errors, (b) cluster-robust variance
estimation, and (c) exact permutation inference - none of which PyMARE 0.0.12 provides.
PyMARE is used in `validate_metalib.py` as an independent cross-check of the REML point
estimates.

References for the formulations used:
  - Fisher z transform:            z = artanh(r),  Var(z) = 1/(n-3)
  - Kendall tau -> Pearson r:      r = sin(pi*tau/2)             (Greiner's relation)
  - REML tau^2:                    Viechtbauer (2005), J Educ Behav Stat 30:261
  - Knapp-Hartung (HKSJ):          Knapp & Hartung (2003), Stat Med 22:2693
  - Prediction interval:           Higgins, Thompson & Spiegelhalter (2009) JRSS-A 172:137
  - Correlated-effects aggregation: Borenstein et al. (2009), Intro to Meta-Analysis, ch. 24
  - Cluster-robust variance (RVE): Hedges, Tipton & Johnson (2010), Res Synth Methods 1:39
  - Egger regression:              Egger et al. (1997), BMJ 315:629
  - Trim-and-fill (L0):            Duval & Tweedie (2000), Biometrics 56:455
  - TOST equivalence:              Lakens (2017), Soc Psychol Personal Sci 8:355
"""
from __future__ import annotations

import itertools
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Scale conversions
# ---------------------------------------------------------------------------

R_CLIP = 0.9999  # guard against artanh(+-1)


def r_to_z(r):
    """Fisher z transform of a correlation."""
    return np.arctanh(np.clip(np.asarray(r, dtype=float), -R_CLIP, R_CLIP))


def z_to_r(z):
    """Inverse Fisher transform."""
    return np.tanh(np.asarray(z, dtype=float))


def var_z(n):
    """Sampling variance of Fisher z for a correlation from n observations."""
    n = np.asarray(n, dtype=float)
    return 1.0 / (n - 3.0)


def tau_to_r(tau):
    """Kendall's tau -> Pearson r under bivariate normality (Greiner's relation)."""
    return np.sin(np.pi * np.asarray(tau, dtype=float) / 2.0)


def r_from_p_and_df(p_two_sided, df, sign):
    """Recover a (partial) correlation from a two-sided p-value and residual df.

    A regression coefficient's two-sided p-value is a monotone function of its partial
    correlation given df, via the t statistic:  t = r*sqrt(df) / sqrt(1-r^2)  =>
    r = t / sqrt(t^2 + df).  This is the standard "test-statistic recovery" conversion
    (Rosenthal & Rubin 1979; Borenstein et al. 2009 sec. 6) and is the only assumption-light
    way to bring unstandardised regression coefficients onto the correlation scale.

    Args:
        p_two_sided: reported two-sided p-value (0 < p < 1).
        df:          residual degrees of freedom of the coefficient's t test.
        sign:        +1 if the effect supports "worse sleep <-> more amyloid", else -1.

    Returns:
        Signed partial correlation r.
    """
    p = float(np.clip(p_two_sided, 1e-12, 0.999999))
    df = float(df)
    if df <= 0:
        return np.nan
    t = stats.t.isf(p / 2.0, df)          # |t| implied by the two-sided p
    r = t / np.sqrt(t * t + df)
    return float(np.sign(sign) * r)


# ---------------------------------------------------------------------------
# Random-effects model
# ---------------------------------------------------------------------------

def _tau2_dl(y, v):
    """DerSimonian-Laird moment estimator of tau^2 (used as the REML starting value)."""
    y, v = np.asarray(y, float), np.asarray(v, float)
    k = len(y)
    if k < 2:
        return 0.0
    w = 1.0 / v
    mu = np.sum(w * y) / np.sum(w)
    Q = float(np.sum(w * (y - mu) ** 2))
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    return max(0.0, (Q - (k - 1)) / c) if c > 0 else 0.0


def _tau2_reml(y, v, max_iter=200, tol=1e-10):
    """REML estimator of tau^2 by Fisher scoring (Viechtbauer 2005, eq. 12-13).

    Iterates
        tau2_{s+1} = tau2_s + [ sum w_i^2 ((y_i - mu)^2 - (v_i + tau2)) + tr-correction ]
                              / sum w_i^2
    with w_i = 1/(v_i + tau2), truncated at 0.
    """
    y, v = np.asarray(y, float), np.asarray(v, float)
    k = len(y)
    if k < 2:
        return 0.0
    t2 = _tau2_dl(y, v)
    for _ in range(max_iter):
        w = 1.0 / (v + t2)
        sw = np.sum(w)
        mu = np.sum(w * y) / sw
        # REML estimating equation numerator/denominator
        num = np.sum(w ** 2 * ((y - mu) ** 2 - v)) + np.sum(w ** 2) / sw
        den = np.sum(w ** 2)
        new = num / den
        new = max(0.0, float(new))
        if abs(new - t2) < tol:
            t2 = new
            break
        t2 = new
    return float(t2)


def re_meta(y, v, method="REML", hksj=True, level=0.95):
    """Fit a random-effects meta-analysis on the Fisher-z scale.

    Args:
        y:      array of effect estimates (Fisher z).
        v:      array of sampling variances.
        method: "REML" (default), "DL", or "FE" (fixed effect, tau^2 = 0).
        hksj:   apply the Knapp-Hartung-Sidik-Jonkman variance adjustment and use a
                t_{k-1} reference distribution. Strongly recommended when k is small,
                which it always is here.
        level:  confidence level.

    Returns:
        dict with the pooled estimate on both z and r scales, CI, p, tau^2, Q, I^2,
        H^2, and a 95% prediction interval.
    """
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    ok = np.isfinite(y) & np.isfinite(v) & (v > 0)
    y, v = y[ok], v[ok]
    k = len(y)
    out = {"k": int(k)}
    if k == 0:
        return {**out, "note": "no effects"}
    if k == 1:
        se = float(np.sqrt(v[0]))
        crit = stats.norm.isf((1 - level) / 2)
        out.update(est=float(y[0]), se=se, ci_low=float(y[0] - crit * se),
                   ci_high=float(y[0] + crit * se), tau2=0.0, Q=0.0, I2=0.0, H2=1.0,
                   p=float(2 * stats.norm.sf(abs(y[0]) / se)), method="single-effect",
                   pi_low=np.nan, pi_high=np.nan)
    else:
        if method == "FE":
            t2 = 0.0
        elif method == "DL":
            t2 = _tau2_dl(y, v)
        else:
            t2 = _tau2_reml(y, v)

        w = 1.0 / (v + t2)
        sw = np.sum(w)
        mu = float(np.sum(w * y) / sw)
        se = float(np.sqrt(1.0 / sw))

        # Knapp-Hartung-Sidik-Jonkman adjustment: rescale by the weighted residual MS
        if hksj and k >= 2:
            q_hksj = float(np.sum(w * (y - mu) ** 2) / (k - 1))
            se_adj = se * np.sqrt(max(q_hksj, 1e-12))
            # Standard safeguard: never let HKSJ shrink the SE below the model-based one
            se = max(se, se_adj)
            crit = stats.t.isf((1 - level) / 2, k - 1)
            pval = float(2 * stats.t.sf(abs(mu) / se, k - 1))
            meth = f"{method}+HKSJ"
        else:
            crit = stats.norm.isf((1 - level) / 2)
            pval = float(2 * stats.norm.sf(abs(mu) / se))
            meth = method

        # Heterogeneity (computed at the fixed-effect weights, standard practice)
        wf = 1.0 / v
        muf = np.sum(wf * y) / np.sum(wf)
        Q = float(np.sum(wf * (y - muf) ** 2))
        dfQ = k - 1
        I2 = float(max(0.0, 100.0 * (Q - dfQ) / Q)) if Q > 0 else 0.0
        H2 = float(Q / dfQ) if dfQ > 0 else np.nan
        Qp = float(stats.chi2.sf(Q, dfQ)) if dfQ > 0 else np.nan

        # 95% prediction interval (Higgins 2009): mu +- t_{k-2} * sqrt(tau2 + se^2)
        if k >= 3:
            tcrit_pi = stats.t.isf((1 - level) / 2, k - 2)
            half = tcrit_pi * np.sqrt(t2 + se ** 2)
            pi_low, pi_high = mu - half, mu + half
        else:
            pi_low = pi_high = np.nan

        out.update(est=mu, se=se, ci_low=float(mu - crit * se),
                   ci_high=float(mu + crit * se), tau2=float(t2), Q=Q, Q_p=Qp,
                   I2=I2, H2=H2, p=pval, method=meth,
                   pi_low=float(pi_low), pi_high=float(pi_high))

    out["r"] = float(z_to_r(out["est"]))
    out["r_ci_low"] = float(z_to_r(out["ci_low"]))
    out["r_ci_high"] = float(z_to_r(out["ci_high"]))
    out["r_pi_low"] = float(z_to_r(out["pi_low"])) if np.isfinite(out.get("pi_low", np.nan)) else None
    out["r_pi_high"] = float(z_to_r(out["pi_high"])) if np.isfinite(out.get("pi_high", np.nan)) else None
    out["n_total"] = None
    return out


# ---------------------------------------------------------------------------
# Dependency handling
# ---------------------------------------------------------------------------

def aggregate_correlated(y, v, rho=0.6):
    """Combine m dependent effects from one cluster into a single synthetic effect.

    Borenstein et al. (2009) ch. 24: the mean of m correlated estimates has variance
        v_bar = (1/m^2) * [ sum_i v_i + sum_{i != j} rho * sqrt(v_i v_j) ]
    where rho is the assumed correlation between the estimates. rho = 1 reproduces the
    single-effect variance (no precision gain); rho = 0 reproduces independence.

    Returns (y_bar, v_bar).
    """
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    m = len(y)
    if m == 0:
        return np.nan, np.nan
    if m == 1:
        return float(y[0]), float(v[0])
    ybar = float(np.mean(y))
    s = np.sqrt(v)
    cross = rho * (np.outer(s, s).sum() - np.sum(v))  # sum_{i != j} rho sqrt(vi vj)
    vbar = float((np.sum(v) + cross) / (m ** 2))
    return ybar, vbar


def rve_meta(y, v, cluster, rho=0.6, small_sample=True, level=0.95):
    """Cluster-robust (sandwich) variance estimation for an intercept-only meta-regression.

    Hedges, Tipton & Johnson (2010). Weights use the approximate inverse-variance form
    w_j = 1 / (m_j * (v_bar_j + tau2)) within cluster j. The small-sample correction uses
    df = m - 1 with an inflation factor of m/(m-1), where m is the number of clusters.
    Only meaningful with >= 4 clusters; we return a warning below that.
    """
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    cluster = np.asarray(cluster)
    clusters = list(dict.fromkeys(cluster.tolist()))
    m = len(clusters)

    # tau^2 from the cluster-aggregated effects (stable with tiny k)
    agg_y, agg_v = [], []
    for c in clusters:
        sel = cluster == c
        a, b = aggregate_correlated(y[sel], v[sel], rho=rho)
        agg_y.append(a)
        agg_v.append(b)
    t2 = _tau2_reml(np.array(agg_y), np.array(agg_v))

    num = den = 0.0
    for c in clusters:
        sel = cluster == c
        mj = int(sel.sum())
        wj = 1.0 / (mj * (np.mean(v[sel]) + t2))
        num += wj * np.sum(y[sel])
        den += wj * mj
    mu = num / den

    meat = 0.0
    for c in clusters:
        sel = cluster == c
        mj = int(sel.sum())
        wj = 1.0 / (mj * (np.mean(v[sel]) + t2))
        ej = np.sum(y[sel] - mu)
        meat += (wj * ej) ** 2
    var = meat / (den ** 2)
    if small_sample and m > 1:
        var *= m / (m - 1.0)
    se = float(np.sqrt(var))
    df = m - 1
    if df >= 1:
        crit = stats.t.isf((1 - level) / 2, df)
        p = float(2 * stats.t.sf(abs(mu) / se, df))
    else:
        crit, p = np.nan, np.nan
    return {"est": float(mu), "se": se, "df": int(df), "n_clusters": int(m),
            "tau2": float(t2), "ci_low": float(mu - crit * se) if df >= 1 else np.nan,
            "ci_high": float(mu + crit * se) if df >= 1 else np.nan, "p": p,
            "r": float(z_to_r(mu)),
            "r_ci_low": float(z_to_r(mu - crit * se)) if df >= 1 else None,
            "r_ci_high": float(z_to_r(mu + crit * se)) if df >= 1 else None,
            "warning": None if m >= 4 else
            f"RVE with only {m} clusters is unreliable (>=4 recommended); reported for completeness"}


# ---------------------------------------------------------------------------
# Contrasts
# ---------------------------------------------------------------------------

def paired_contrast(z_a, v_a, z_b, v_b, rho_within=0.5):
    """Within-study paired difference of two Fisher-z effects from the same sample.

    Two correlations measured on the same participants are themselves correlated. Their
    difference has variance
        Var(z_a - z_b) = v_a + v_b - 2*rho*sqrt(v_a*v_b)
    where rho is the correlation between the two estimates. rho is not identified from
    published summaries, so it is assumed and varied in sensitivity analysis. Larger rho
    means a *smaller* variance, i.e. rho = 0 is the conservative choice for detecting a
    difference; we use rho = 0.5 as the primary value.
    """
    d = float(z_a - z_b)
    vd = float(v_a + v_b - 2.0 * rho_within * np.sqrt(v_a * v_b))
    vd = max(vd, 1e-12)
    return d, vd


def between_stratum_contrast(res_a, res_b):
    """Difference between two independently-pooled strata, on the Fisher-z scale.

    Uses a Wald test with the two pooled standard errors. Valid only when the strata
    contain disjoint samples; when they share cohorts this is anti-conservative and the
    paired contrast should be preferred.
    """
    d = res_a["est"] - res_b["est"]
    se = np.sqrt(res_a["se"] ** 2 + res_b["se"] ** 2)
    zst = d / se
    return {"diff_z": float(d), "se": float(se), "z_stat": float(zst),
            "p": float(2 * stats.norm.sf(abs(zst))),
            "ci_low": float(d - 1.96 * se), "ci_high": float(d + 1.96 * se),
            "delta_r": float(z_to_r(res_a["est"]) - z_to_r(res_b["est"]))}


# ---------------------------------------------------------------------------
# Small-study bias
# ---------------------------------------------------------------------------

def egger_test(y, v):
    """Egger's regression test for funnel-plot asymmetry.

    Regresses the standard normal deviate (y/se) on precision (1/se); the intercept is the
    asymmetry test. Equivalent to weighted regression of y on se.
    """
    y, v = np.asarray(y, float), np.asarray(v, float)
    k = len(y)
    if k < 3:
        return {"k": k, "note": "k < 3, not computable"}
    se = np.sqrt(v)
    snd = y / se
    prec = 1.0 / se
    res = stats.linregress(prec, snd)
    return {"k": int(k), "intercept": float(res.intercept),
            "intercept_se": float(res.intercept_stderr), "slope": float(res.slope),
            "p": float(_intercept_p(res, k)),
            "interpretable": bool(k >= 5),
            "note": None if k >= 5 else "k < 5: Egger's test has essentially no power here"}


def _intercept_p(res, k):
    t = res.intercept / res.intercept_stderr if res.intercept_stderr > 0 else np.nan
    return 2 * stats.t.sf(abs(t), k - 2) if np.isfinite(t) else np.nan


def _l0_missing_left(y, v, max_iter=100):
    """Duval & Tweedie L0: estimated number of studies missing from the LEFT of the funnel.

    Iterates: trim the k0 largest effects, re-estimate the centre from the trimmed set,
    then recompute k0 from the signed ranks of the FULL set's residuals about that centre.
        Tn = sum of the ranks of the positive-signed residuals
        L0 = (4*Tn - n(n+1)) / (2n - 1)
    A positive L0 means positive effects dominate the funnel more than symmetry allows,
    i.e. small studies with SMALL effects appear to be missing.
    """
    y, v = np.asarray(y, float), np.asarray(v, float)
    n = len(y)
    k0 = 0
    order_desc = np.argsort(-y)          # largest effect first
    for _ in range(max_iter):
        keep = np.ones(n, bool)
        if k0 > 0:
            keep[order_desc[:k0]] = False
        if keep.sum() < 2:
            break
        mu = M_pool(y[keep], v[keep])
        d = y - mu
        nz = d != 0
        if nz.sum() == 0:
            new_k0 = 0
        else:
            ranks = np.zeros(n)
            ranks[nz] = stats.rankdata(np.abs(d[nz]))
            Tn = float(np.sum(ranks[d > 0]))
            m = int(nz.sum())
            new_k0 = int(max(0, round((4 * Tn - m * (m + 1)) / (2 * m - 1))))
        new_k0 = min(new_k0, n - 2)
        if new_k0 == k0:
            break
        k0 = new_k0
    keep = np.ones(n, bool)
    if k0 > 0:
        keep[order_desc[:k0]] = False
    mu = M_pool(y[keep], v[keep]) if keep.sum() >= 2 else M_pool(y, v)
    return k0, float(mu), order_desc[:k0]


def M_pool(y, v):
    """Random-effects centre used inside trim-and-fill (no HKSJ - a point estimate only)."""
    return re_meta(y, v, hksj=False)["est"]


def trim_and_fill(y, v, side="auto", max_iter=100):
    """Duval & Tweedie (2000) trim-and-fill with the L0 estimator.

    Estimates how many studies are missing from one side of the funnel plot, imputes their
    mirror images about the trimmed centre, and refits.

    side: "left"  - suspect small-effect studies are missing (funnel skewed positive)
          "right" - suspect large-effect studies are missing (funnel skewed negative)
          "auto"  - run both orientations and take whichever imputes more studies

    NOTE on power: trim-and-fill is generally regarded as uninformative below about ten
    studies. Every stratum here has 3-6 cohort families, so the result is reported for
    completeness and explicitly flagged as uninterpretable.
    """
    y, v = np.asarray(y, float), np.asarray(v, float)
    k = len(y)
    if k < 3:
        return {"k": int(k), "note": "k < 3, not computable"}

    kl, mul, idxl = _l0_missing_left(y, v)
    kr, mur, idxr = _l0_missing_left(-y, v)

    if side == "auto":
        side = "left" if kl >= kr else "right"
    if side == "left":
        k0, mu, idx = kl, mul, idxl
    else:
        k0, mu, idx = kr, -mur, idxr

    base = re_meta(y, v, hksj=False)
    out = {"k": int(k), "side": side, "n_imputed": int(k0),
           "unadjusted_z": float(base["est"]), "unadjusted_r": float(z_to_r(base["est"])),
           "interpretable": bool(k >= 10),
           "note": None if k >= 10 else
           "k < 10: trim-and-fill has essentially no power and should not be interpreted"}
    if k0 <= 0:
        out.update(adjusted_z=float(base["est"]), adjusted_r=float(z_to_r(base["est"])),
                   adjusted_ci_r=[float(z_to_r(base["ci_low"])),
                                  float(z_to_r(base["ci_high"]))])
        return out
    # mirror the trimmed studies about the trimmed centre and refit on the augmented set
    centre = mu
    mirrored = 2 * centre - y[idx]
    y_new = np.concatenate([y, mirrored])
    v_new = np.concatenate([v, v[idx]])
    adj = re_meta(y_new, v_new, hksj=False)
    out.update(adjusted_z=float(adj["est"]), adjusted_r=float(z_to_r(adj["est"])),
               adjusted_ci_r=[float(z_to_r(adj["ci_low"])), float(z_to_r(adj["ci_high"]))],
               trimmed_centre_z=float(centre))
    return out


# ---------------------------------------------------------------------------
# Exact / resampling inference for tiny k
# ---------------------------------------------------------------------------

def permutation_p(y, v, n_max_exact=2 ** 16):
    """Exact sign-flip permutation p-value for H0: pooled effect = 0.

    With k effects there are 2^k sign assignments; when 2^k is small we enumerate all of
    them (exact), otherwise we sample. The test statistic is the inverse-variance-weighted
    mean. This is valid under the (weak) assumption that the effects are symmetrically
    distributed about 0 under H0, and is far more trustworthy than a Wald p at k <= 6.
    """
    y, v = np.asarray(y, float), np.asarray(v, float)
    k = len(y)
    if k < 2:
        return {"k": k, "note": "k < 2"}
    w = 1.0 / v
    obs = abs(np.sum(w * y) / np.sum(w))
    if 2 ** k <= n_max_exact:
        signs = np.array(list(itertools.product([-1, 1], repeat=k)), float)
        exact = True
    else:
        rng = np.random.default_rng(42)
        signs = rng.choice([-1.0, 1.0], size=(n_max_exact, k))
        exact = False
    stat = np.abs((signs * y) @ w / np.sum(w))
    p = float((np.sum(stat >= obs - 1e-12)) / len(stat))
    return {"k": int(k), "observed_weighted_mean_z": float(np.sum(w * y) / np.sum(w)),
            "p_permutation": p, "exact": bool(exact), "n_permutations": int(len(stat))}


def tost_equivalence(est, se, df, bound_z):
    """Two one-sided tests of equivalence against +-bound_z on the Fisher-z scale.

    H0: |true effect| >= bound. Rejecting both one-sided tests concludes equivalence, i.e.
    that the true effect is smaller in magnitude than the bound.
    """
    t_lo = (est - (-bound_z)) / se     # test against lower bound
    t_hi = (est - bound_z) / se        # test against upper bound
    p_lo = float(stats.t.sf(t_lo, df))       # H0: effect <= -bound
    p_hi = float(stats.t.cdf(t_hi, df))      # H0: effect >= +bound
    p = max(p_lo, p_hi)
    return {"bound_z": float(bound_z), "bound_r": float(z_to_r(bound_z)),
            "p_lower": p_lo, "p_upper": p_hi, "p_tost": p,
            "equivalent_at_05": bool(p < 0.05)}


def detectable_effect(v_list, power=0.80, alpha=0.05, tau2=0.0):
    """Smallest true Fisher-z detectable with the given power by a random-effects model.

    Uses the fixed-weight approximation se = sqrt(1/sum(1/(v_i+tau2))); the minimum
    detectable effect is (z_{1-a/2} + z_{power}) * se.
    """
    v = np.asarray(v_list, float)
    if len(v) == 0:
        return {"note": "no effects"}
    se = float(np.sqrt(1.0 / np.sum(1.0 / (v + tau2))))
    mde_z = (stats.norm.isf(alpha / 2) + stats.norm.ppf(power)) * se
    return {"k": int(len(v)), "pooled_se_z": se, "mde_z": float(mde_z),
            "mde_r": float(z_to_r(mde_z)), "power": power, "alpha": alpha,
            "tau2_assumed": float(tau2)}


def power_for_effect(v_list, true_z, alpha=0.05, tau2=0.0):
    """Power of the pooled test to detect a given true Fisher-z effect."""
    v = np.asarray(v_list, float)
    se = float(np.sqrt(1.0 / np.sum(1.0 / (v + tau2))))
    lam = abs(true_z) / se
    crit = stats.norm.isf(alpha / 2)
    return float(stats.norm.sf(crit - lam) + stats.norm.cdf(-crit - lam))
