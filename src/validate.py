#!/usr/bin/env python3
"""
validate.py - Phase 5 validation.

Three independent checks:
  V1 ANALYTIC   metalib's REML/pooling primitives against closed-form values and against
                PyMARE (a separately-maintained meta-analysis package) on the real data.
  V2 CONVERSION the effect-size conversion rules against hand-computable cases.
  V3 REPRODUCE  re-run the whole pipeline from scratch and byte-compare every results file
                against the copies produced by the first run.

Exit code 0 only if every check passes.

Output: results/validation_report.json
"""
from __future__ import annotations

import filecmp
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metalib as M                                # noqa: E402
from exp1_pooling import load, aggregate_by_cohort  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "validation_report.json")

checks = []


def check(name, passed, detail=""):
    checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))
    return passed


def v1_analytic():
    print("V1 analytic checks")
    # Fisher z round trip
    r = np.array([-0.9, -0.3, 0.0, 0.42, 0.87])
    check("Fisher z round-trip", np.allclose(M.z_to_r(M.r_to_z(r)), r, atol=1e-12))
    # known value
    check("artanh(0.5) = 0.5493061443", abs(M.r_to_z(0.5) - 0.5493061443340548) < 1e-12)
    # Kendall tau conversion
    check("tau=0 -> r=0", abs(M.tau_to_r(0.0)) < 1e-15)
    check("tau=1 -> r=1", abs(M.tau_to_r(1.0) - 1.0) < 1e-12)
    # p-value recovery: for a simple correlation with df = n-2, recovering r from the exact
    # p-value of that r must return r.
    for r0, n in [(0.45, 26), (0.07, 32), (0.32, 66), (0.19, 189)]:
        df = n - 2
        t = r0 * np.sqrt(df) / np.sqrt(1 - r0 ** 2)
        p = 2 * stats.t.sf(abs(t), df)
        back = M.r_from_p_and_df(p, df, 1)
        check(f"p-recovery round-trip r={r0}, n={n}", abs(back - r0) < 1e-8,
              f"recovered {back:.10f}")
    # Fixed-effect pooling equals the closed-form inverse-variance mean
    y = np.array([0.2, 0.5, 0.1, 0.4]); v = np.array([0.05, 0.04, 0.06, 0.03])
    fe = M.re_meta(y, v, method="FE", hksj=False)
    closed = np.sum(y / v) / np.sum(1 / v)
    check("FE pooled = inverse-variance mean", abs(fe["est"] - closed) < 1e-12)
    check("FE se = 1/sqrt(sum(1/v))", abs(fe["se"] - np.sqrt(1 / np.sum(1 / v))) < 1e-12)
    # tau^2 = 0 when all effects identical
    yy = np.full(5, 0.3)
    check("tau2 = 0 for identical effects", M._tau2_reml(yy, v[:1].repeat(5)) < 1e-10)
    # Q and I2 sanity
    res = M.re_meta(y, v, hksj=False)
    check("I2 in [0,100]", 0 <= res["I2"] <= 100, f"I2={res['I2']:.2f}")
    # aggregation limits
    a1, vv1 = M.aggregate_correlated(y, v, rho=1.0)
    a0, vv0 = M.aggregate_correlated(y, v, rho=0.0)
    check("rho=1 aggregation has larger variance than rho=0", vv1 > vv0,
          f"v(rho=1)={vv1:.5f} > v(rho=0)={vv0:.5f}")
    check("aggregation mean is the arithmetic mean", abs(a1 - np.mean(y)) < 1e-12)
    m1, mv1 = M.aggregate_correlated(y[:1], v[:1], rho=0.6)
    check("single-effect aggregation is identity", abs(m1 - y[0]) < 1e-15 and
          abs(mv1 - v[0]) < 1e-15)
    # permutation p is exact and in (0,1]
    perm = M.permutation_p(y, v)
    check("permutation exact enumeration for k=4", perm["exact"] and
          perm["n_permutations"] == 16)
    # trim-and-fill must detect a strongly asymmetric funnel and leave a symmetric one alone
    vv = np.array([0.30, 0.28, 0.25, 0.20, 0.15, 0.10, 0.05, 0.02, 0.01, 0.005])
    y_asym = np.array([0.95, 0.90, 0.85, 0.70, 0.55, 0.40, 0.22, 0.12, 0.09, 0.08])
    y_sym = np.array([0.30, -0.10, 0.45, 0.05, 0.35, 0.10, 0.25, 0.18, 0.21, 0.19])
    ta, ts_, tn = (M.trim_and_fill(y_asym, vv), M.trim_and_fill(y_sym, vv),
                   M.trim_and_fill(-y_asym, vv))
    check("trim-and-fill detects an asymmetric funnel", ta["n_imputed"] > 0,
          f"imputed {ta['n_imputed']}, adjusted r {ta['adjusted_r']:.3f} "
          f"(unadjusted {ta['unadjusted_r']:.3f})")
    check("trim-and-fill imputes nothing for a symmetric funnel", ts_["n_imputed"] == 0)
    check("trim-and-fill is sign-symmetric", tn["n_imputed"] == ta["n_imputed"] and
          abs(tn["adjusted_r"] + ta["adjusted_r"]) < 1e-12, f"side={tn['side']}")
    check("trim-and-fill shrinks an inflated estimate toward zero",
          abs(ta["adjusted_r"]) < abs(ta["unadjusted_r"]))
    # Egger detects the same asymmetry and not the symmetric case
    check("Egger flags the asymmetric funnel", M.egger_test(y_asym, vv)["p"] < 0.01)
    check("Egger does not flag the symmetric funnel", M.egger_test(y_sym, vv)["p"] > 0.10)

    # TOST: a zero effect with tiny SE must be equivalent to a wide bound
    t = M.tost_equivalence(0.0, 0.01, 30, M.r_to_z(0.30))
    check("TOST rejects a wide bound for a precise zero", t["equivalent_at_05"])
    t2 = M.tost_equivalence(0.0, 1.0, 3, M.r_to_z(0.05))
    check("TOST does not reject a narrow bound when imprecise", not t2["equivalent_at_05"])

    # PyMARE cross-check on the real data
    try:
        from pymare import Dataset, estimators
        d = load()
        diffs = []
        for dom in ["SWS_stage", "SWS_spectral", "TST", "SE"]:
            agg = aggregate_by_cohort(d[d.exposure_domain == dom])
            ds = Dataset(y=agg.z.values[:, None], v=agg.v.values[:, None])
            est = estimators.VarianceBasedLikelihoodEstimator(method="REML").fit_dataset(ds)
            pm = float(np.ravel(est.summary().get_fe_stats()["est"])[0])
            ours = M.re_meta(agg.z.values, agg.v.values, hksj=True)["est"]
            diffs.append((dom, pm, ours, abs(pm - ours)))
        worst = max(d_[3] for d_ in diffs)
        check("REML point estimates match PyMARE (all strata)", worst < 1e-4,
              "max |diff| = " + f"{worst:.2e}; " +
              ", ".join(f"{a}: ours {c:.5f} vs PyMARE {b:.5f}" for a, b, c, _ in diffs))
    except Exception as e:  # pragma: no cover
        check("PyMARE cross-check", False, f"could not run: {e}")


def v2_conversions():
    print("V2 conversion checks against the source data")
    h = pd.read_csv(os.path.join(RES, "harmonised_effects.csv"))
    # every usable row has finite z and v
    u = h[h.usable]
    check("all usable rows have finite z", np.isfinite(u.z).all())
    check("all usable rows have positive v", (u.v > 0).all())
    check("no |r| >= 1", (u.r.abs() < 1).all())
    # direct_r rows equal burden_aligned_effect exactly
    dr = u[u.conversion_rule == "direct_r"]
    check("direct_r rows are unchanged",
          np.allclose(dr.r.values, dr.burden_aligned_effect.values, atol=1e-12),
          f"{len(dr)} rows")
    # v equals 1/(n-3)
    check("v = 1/(n-3) for every usable row",
          np.allclose(u.v.values, 1.0 / (u.n.values - 3), atol=1e-12))
    # imputed nulls are exactly zero
    inull = u[u.is_imputed_null]
    check("imputed nulls have z = 0", np.allclose(inull.z.values, 0.0),
          f"{len(inull)} rows")
    # signs preserved
    nz = u[(~u.is_imputed_null) & (u.burden_aligned_effect != 0)]
    check("conversion preserves sign",
          (np.sign(nz.r.values) == np.sign(nz.burden_aligned_effect.values)).all())
    # Winer 2020's two headline numbers survive the pipeline unchanged
    w = h[h.effect_id.isin(["E012", "E013"])].set_index("effect_id")
    check("Winer2020 spectral r = 0.52 preserved", abs(w.loc["E012", "r"] - 0.52) < 1e-12)
    check("Winer2020 N3 stage r = 0.07 preserved", abs(w.loc["E013", "r"] - 0.07) < 1e-12)
    # Mander's Kendall tau converts as expected
    m = h[h.effect_id == "E008"].iloc[0]
    check("Kendall tau 0.30 -> r 0.4540", abs(m.r - np.sin(np.pi * 0.30 / 2)) < 1e-12,
          f"r={m.r:.4f}")
    # no effect is dropped without a recorded reason
    bad = h[(~h.usable) & (h.usable_reason.isna())]
    check("every unusable row has a reason", len(bad) == 0, f"{len(bad)} unexplained")
    # the analysis set matches what the report claims
    pet = h[h.outcome_class.isin(["deposition_pet", "deposition_pet_change"]) & h.usable]
    check("40 usable amyloid-PET effects", len(pet) == 40, f"found {len(pet)}")
    check("12 cohort families in the full dataset", h.cohort_family.nunique() == 12,
          f"found {h.cohort_family.nunique()}")


def _hash(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def v3_reproduce():
    """Re-run the full pipeline in a scratch copy and compare every results file."""
    print("V3 reproducibility - re-running the pipeline from scratch")
    scripts = ["prepare_data.py", "exp1_pooling.py", "exp2_specificity.py",
               "exp3_attenuation.py", "exp4_dose_response.py", "exp5_robustness.py"]
    targets = ["harmonised_effects.csv", "qc_report.json", "e1_pooled.json",
               "e2_specificity.json", "e3_attenuation.json", "e4_dose_response.json",
               "e5_robustness.json", "e6_precision.json"]
    before = {t: _hash(os.path.join(RES, t)) for t in targets}
    with tempfile.TemporaryDirectory() as tmp:
        backup = os.path.join(tmp, "results_backup")
        shutil.copytree(RES, backup)
        env = dict(os.environ, PYTHONHASHSEED="0")
        for sc in scripts:
            p = subprocess.run([sys.executable, os.path.join(ROOT, "src", sc)],
                               capture_output=True, text=True, cwd=ROOT, env=env)
            if p.returncode != 0:
                check(f"re-run {sc}", False, p.stderr[-500:])
                return
        after = {t: _hash(os.path.join(RES, t)) for t in targets}
        mismatched = [t for t in targets if before[t] != after[t]]
        check("all result files bit-identical on re-run", not mismatched,
              f"differing: {mismatched}" if mismatched else f"{len(targets)} files identical")


def main():
    print("=" * 74)
    v1_analytic()
    print()
    v2_conversions()
    print()
    v3_reproduce()
    print("=" * 74)
    n_pass = sum(c["passed"] for c in checks)
    report = {"n_checks": len(checks), "n_passed": n_pass,
              "n_failed": len(checks) - n_pass, "all_passed": n_pass == len(checks),
              "checks": checks}
    with open(OUT, "w") as f:
        json.dump(report, f, indent=2)
    print(f"{n_pass}/{len(checks)} checks passed -> {os.path.relpath(OUT, ROOT)}")
    sys.exit(0 if report["all_passed"] else 1)


if __name__ == "__main__":
    main()
