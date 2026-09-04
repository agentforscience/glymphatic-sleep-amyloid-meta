#!/usr/bin/env python3
"""
exp5_robustness.py - E5 (robustness / sensitivity / bias) and E6 (precision).

E5 sensitivity analyses, all preregistered in planning.md section A.4/A.6:
  S1  exclude the imputed narrative nulls
  S2  exclude effects recovered from p-values (unstandardised betas)
  S3  vary the residual-df rule used in the p-value recovery
  S4  vary the within-cohort dependency rho in {0.3, 0.6, 0.8}
  S5  restrict to strictly cognitively normal samples
  S6  restrict to PSG-measured exposures (like-for-like measurement)
  S7  leave-one-cohort-family-out, per stratum
  S8  fixed-effect instead of random-effects
  S9  small-study bias: Egger regression and trim-and-fill
  S10 CSF secondary analysis, stratified by cognitive status (the sign-flip check)
  S11 benchmark against the published pooled estimates of Chen 2025

E6 precision:
  minimum detectable effect at 80% power per stratum, power to detect the r = 0.30-0.45
  effects reported by the Berkeley cohort, and TOST equivalence tests for the N3 stratum.

Output: results/e5_robustness.json, results/e6_precision.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metalib as M                                       # noqa: E402
from exp1_pooling import load, aggregate_by_cohort, RHO_PRIMARY, PET  # noqa: E402
from prepare_data import covariate_count                  # noqa: E402

np.random.seed(42)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARM = os.path.join(ROOT, "results", "harmonised_effects.csv")
OUT5 = os.path.join(ROOT, "results", "e5_robustness.json")
OUT6 = os.path.join(ROOT, "results", "e6_precision.json")

MAIN = ["SWS_stage", "SWS_spectral", "TST", "SE"]


def pool_stratum(sub, rho=RHO_PRIMARY):
    """Cohort-aggregate then pool; returns None for an empty stratum."""
    if len(sub) == 0:
        return None
    agg = aggregate_by_cohort(sub, rho=rho)
    out = M.re_meta(agg.z.values, agg.v.values, hksj=True)
    out["n_cohorts"] = int(len(agg))
    out["k_effects"] = int(len(sub))
    return out


def compact(res):
    """Trim a pooling result down to the fields worth tabulating."""
    if res is None:
        return None
    return {k: res.get(k) for k in
            ["k_effects", "n_cohorts", "r", "r_ci_low", "r_ci_high", "p", "tau2", "I2"]}


def run_all_strata(d, rho=RHO_PRIMARY):
    return {dom: compact(pool_stratum(d[d.exposure_domain == dom], rho)) for dom in MAIN}


def recompute_recovered(d, df_rule):
    """Re-derive r for p-recovered rows under an alternative residual-df rule."""
    d = d.copy()
    mask = d.conversion_rule == "recover_r_from_p_and_df"
    for i in d.index[mask]:
        n = float(d.at[i, "n"])
        ncov, _ = covariate_count(d.at[i, "adjustment"])
        if df_rule == "n_minus_2":
            dfree = n - 2
        elif df_rule == "n_minus_1_minus_2ncov":
            dfree = n - 1 - 2 * ncov
        else:                                   # primary: n - 1 - ncov
            dfree = n - 1 - ncov
        dfree = max(dfree, 2.0)
        r = M.r_from_p_and_df(d.at[i, "p_value"], dfree,
                              np.sign(d.at[i, "burden_aligned_effect"]) or 1.0)
        d.at[i, "r"] = r
        d.at[i, "z"] = float(M.r_to_z(r))
        d.at[i, "df_used"] = dfree
    return d


def main():
    d = load()
    e1 = json.load(open(os.path.join(ROOT, "results", "e1_pooled.json")))
    full = pd.read_csv(HARM)

    R = {"experiment": "E5 - robustness, sensitivity and bias diagnostics",
         "primary_for_reference": {dom: {
             "r": e1["strata"][dom]["primary"]["r"],
             "ci": [e1["strata"][dom]["primary"]["r_ci_low"],
                    e1["strata"][dom]["primary"]["r_ci_high"]],
             "p": e1["strata"][dom]["primary"]["p"]} for dom in MAIN}}

    # ---- S1 / S2 inclusion rules ---------------------------------------------------
    R["S1_exclude_imputed_nulls"] = {
        "rationale": ("Imputing z=0 for reported nulls is conservative in magnitude but adds "
                      "precision that the papers did not actually supply. Removing them tests "
                      "whether the conclusions depend on the imputation."),
        "results": run_all_strata(load(include_imputed_nulls=False))}
    R["S2_exclude_p_recovered_betas"] = {
        "rationale": ("Tests whether the t-statistic recovery of unstandardised regression "
                      "coefficients drives any conclusion. This drops the two largest cohorts "
                      "(A4, AIBL) from the TST stratum."),
        "results": run_all_strata(load(include_recovered=False))}
    R["S1_and_S2_both_excluded"] = {
        "rationale": "Most restrictive set: only effects reported natively as correlations.",
        "results": run_all_strata(load(include_imputed_nulls=False, include_recovered=False))}

    # ---- S3 df rule ----------------------------------------------------------------
    R["S3_df_rule_for_p_recovery"] = {"rationale":
        "The residual df used in the t-recovery is not always fully recoverable from the "
        "papers; this varies the rule across a plausible range.", "results": {}}
    for rule in ["n_minus_1_minus_ncov", "n_minus_2", "n_minus_1_minus_2ncov"]:
        dd = recompute_recovered(d, rule)
        R["S3_df_rule_for_p_recovery"]["results"][rule] = run_all_strata(dd)

    # ---- S4 rho ---------------------------------------------------------------------
    R["S4_dependency_rho"] = {"rationale":
        "rho is the assumed correlation between multiple effects from one cohort. rho=1 "
        "would mean they carry no extra information; rho=0 would treat them as independent.",
        "results": {f"rho={r}": run_all_strata(d, rho=r) for r in [0.3, 0.6, 0.8, 1.0]}}

    # ---- S5 cognitive status --------------------------------------------------------
    cn = d[d.cognitive_status.str.contains("cognitively_normal|cognitively_unimpaired",
                                           case=False, na=False)]
    R["S5_cognitively_normal_only"] = {
        "rationale": "The hypothesis is explicitly about cognitively normal adults.",
        "excluded_studies": sorted(set(d.study_id) - set(cn.study_id)),
        "results": run_all_strata(cn)}

    # ---- S6 PSG only ----------------------------------------------------------------
    R["S6_PSG_measured_exposures_only"] = {
        "rationale": ("N3 stage can only be measured by PSG, so comparing it against "
                      "self-reported TST/SE mixes measurement modalities."),
        "results": run_all_strata(d[d.exposure_modality == "PSG_lab"])}

    # ---- S7 leave-one-cohort-out ----------------------------------------------------
    R["S7_leave_one_cohort_family_out"] = {}
    for dom in MAIN:
        sub = d[d.exposure_domain == dom]
        rows = []
        for cf in sorted(sub.cohort_family.unique()):
            res = pool_stratum(sub[sub.cohort_family != cf])
            rows.append({"omitted_cohort_family": cf,
                         "r": res["r"] if res else None,
                         "ci": [res["r_ci_low"], res["r_ci_high"]] if res else None,
                         "p": res["p"] if res else None,
                         "n_cohorts_remaining": res["n_cohorts"] if res else 0})
        R["S7_leave_one_cohort_family_out"][dom] = rows

    # ---- S8 fixed effect -------------------------------------------------------------
    R["S8_fixed_effect"] = {}
    for dom in MAIN:
        agg = aggregate_by_cohort(d[d.exposure_domain == dom])
        fe = M.re_meta(agg.z.values, agg.v.values, method="FE", hksj=False)
        R["S8_fixed_effect"][dom] = {k: fe.get(k) for k in
                                     ["k", "r", "r_ci_low", "r_ci_high", "p"]}

    # ---- S9 small-study bias ---------------------------------------------------------
    R["S9_small_study_bias"] = {}
    for dom in MAIN:
        sub = d[d.exposure_domain == dom]
        agg = aggregate_by_cohort(sub)
        R["S9_small_study_bias"][dom] = {
            "level": "cohort-family aggregated",
            "k": int(len(agg)),
            "egger": M.egger_test(agg.z.values, agg.v.values),
            "trim_and_fill": M.trim_and_fill(agg.z.values, agg.v.values),
            "egger_all_effects": M.egger_test(sub.z.values, sub.v.values),
        }

    # ---- S10 CSF secondary ------------------------------------------------------------
    csf = full[(full.outcome_class == "soluble_csf") & full.usable].copy()
    R["S10_CSF_secondary"] = {
        "rationale": ("CSF Abeta42 inversely proxies deposition, and its sign relative to "
                      "burden flips between pre-deposition and post-deposition samples. "
                      "Pooled here only stratified by cognitive status, never with PET."),
        "n_effects": int(len(csf)),
        "by_cognitive_status": {},
        "note": None}
    for status, g in csf.groupby("cognitive_status"):
        R["S10_CSF_secondary"]["by_cognitive_status"][status] = {
            "studies": sorted(g.study_id.unique().tolist()),
            "k": int(len(g)),
            "by_domain": {dom: compact(pool_stratum(g[g.exposure_domain == dom]))
                          for dom in g.exposure_domain.unique()}}
    R["S10_CSF_secondary"]["note"] = (
        "Varga 2016 (cognitively normal, pre-deposition) and Liguori 2020 (mixed SCI/MCI/AD) "
        "give opposite-signed relationships for the same construct once aligned to burden; "
        "this is why CSF is excluded from the primary model.")

    # ---- S11 benchmark ---------------------------------------------------------------
    prior = pd.read_csv(os.path.join(ROOT, "datasets", "meta_analysis",
                                     "prior_meta_analyses.csv"))
    R["S11_benchmark_vs_published"] = {
        "chen2025_sleep_quality_vs_amyloid_PET": {
            "fisher_z": 0.153, "ci": [0.005, 0.302], "k": 9, "I2_pct": 99,
            "as_r": float(M.z_to_r(0.153))},
        "chen2025_sleep_duration_vs_amyloid_PET": {
            "reported_fisher_z": -0.002, "ci": [-0.003, -0.001], "k": 6,
            "caveat": ("these are raw regression-coefficient units mislabelled as Fisher z; "
                       "a Fisher z of -0.002 with that CI is not a credible pooled "
                       "correlation and is not comparable to our estimates")},
        "our_TST_vs_amyloid_PET": R["primary_for_reference"]["TST"],
        "our_SWS_stage_vs_amyloid_PET": R["primary_for_reference"]["SWS_stage"],
        "comparison": ("Our TST estimate (r = "
                       f"{R['primary_for_reference']['TST']['r']:.3f}) is of the same order as "
                       "Chen 2025's sleep-quality pooling (r = 0.152), which is reassuring for "
                       "the harmonisation; the N3-stage estimate is an order of magnitude "
                       "smaller and has no published counterpart."),
        "prior_rows_available": int(len(prior))}

    with open(OUT5, "w") as f:
        json.dump(R, f, indent=2, default=str)
    print(f"[exp5] wrote {OUT5}")

    # ================= E6 precision ==================================================
    P = {"experiment": "E6 - precision and detectable effects",
         "note": ("An honest null must state what it could have detected. MDE is the smallest "
                  "true Fisher-z the stratum's pooled precision could detect with 80% power at "
                  "alpha=0.05, computed at the cohort-aggregated variances.")}
    P["by_stratum"] = {}
    for dom in MAIN:
        agg = aggregate_by_cohort(d[d.exposure_domain == dom])
        tau2 = e1["strata"][dom]["primary"]["tau2"]
        mde0 = M.detectable_effect(agg.v.values, tau2=0.0)
        mdet = M.detectable_effect(agg.v.values, tau2=tau2)
        P["by_stratum"][dom] = {
            "n_cohorts": int(len(agg)),
            "mde_r_assuming_no_heterogeneity": mde0["mde_r"],
            "mde_r_at_estimated_tau2": mdet["mde_r"],
            "tau2_estimated": tau2,
            "power_to_detect_r_0.20": M.power_for_effect(agg.v.values, M.r_to_z(0.20), tau2=tau2),
            "power_to_detect_r_0.30": M.power_for_effect(agg.v.values, M.r_to_z(0.30), tau2=tau2),
            "power_to_detect_r_0.45_Berkeley_size":
                M.power_for_effect(agg.v.values, M.r_to_z(0.45), tau2=tau2),
        }
    # TOST equivalence on the N3 stratum
    n3 = e1["strata"]["SWS_stage"]["primary"]
    P["equivalence_tests_SWS_stage"] = {
        "estimate_r": n3["r"], "estimate_z": n3["est"], "se_z": n3["se"],
        "df": n3["k"] - 1,
        "tests": {f"bound_r_{b}": M.tost_equivalence(n3["est"], n3["se"], n3["k"] - 1,
                                                     float(M.r_to_z(b)))
                  for b in [0.10, 0.15, 0.20, 0.30, 0.45]},
        "interpretation": ("A significant TOST at bound b means the data are inconsistent with "
                           "a true association as large as r = b, i.e. positive evidence for a "
                           "small-or-null effect rather than mere absence of evidence.")}
    # same for the spectral stratum, which is the more interesting scientific claim
    sp = e1["strata"]["SWS_spectral"]["primary"]
    P["equivalence_tests_SWS_spectral"] = {
        "estimate_r": sp["r"], "se_z": sp["se"], "df": sp["k"] - 1,
        "tests": {f"bound_r_{b}": M.tost_equivalence(sp["est"], sp["se"], sp["k"] - 1,
                                                    float(M.r_to_z(b)))
                  for b in [0.20, 0.30, 0.45]}}

    with open(OUT6, "w") as f:
        json.dump(P, f, indent=2, default=str)
    print(f"[exp6] wrote {OUT6}\n")

    # -------- console summaries -------------------------------------------------------
    print("Sensitivity analyses - pooled r for the N3-stage stratum and comparators")
    hdr = f"{'analysis':<38}{'N3 stage':>12}{'spectral':>11}{'TST':>10}{'SE':>10}"
    print(hdr); print("-" * len(hdr))

    def line(name, res):
        def g(dom):
            v = res.get(dom)
            return f"{v['r']:+.3f}" if v else "   --"
        print(f"{name:<38}{g('SWS_stage'):>12}{g('SWS_spectral'):>11}{g('TST'):>10}{g('SE'):>10}")

    line("PRIMARY", {k: {"r": v["r"]} for k, v in R["primary_for_reference"].items()})
    line("S1 no imputed nulls", R["S1_exclude_imputed_nulls"]["results"])
    line("S2 no p-recovered betas", R["S2_exclude_p_recovered_betas"]["results"])
    line("S1+S2 correlations only", R["S1_and_S2_both_excluded"]["results"])
    for rule, res in R["S3_df_rule_for_p_recovery"]["results"].items():
        line(f"S3 df = {rule}", res)
    for k, res in R["S4_dependency_rho"]["results"].items():
        line(f"S4 {k}", res)
    line("S5 cognitively normal only", R["S5_cognitively_normal_only"]["results"])
    line("S6 PSG exposures only", R["S6_PSG_measured_exposures_only"]["results"])
    line("S8 fixed effect", R["S8_fixed_effect"])

    print("\nLeave-one-cohort-family-out, N3 stage stratum:")
    for row in R["S7_leave_one_cohort_family_out"]["SWS_stage"]:
        print(f"   omit {row['omitted_cohort_family']:<22} r = {row['r']:+.3f} "
              f"[{row['ci'][0]:+.3f},{row['ci'][1]:+.3f}]")
    print("Leave-one-cohort-family-out, spectral SWA stratum:")
    for row in R["S7_leave_one_cohort_family_out"]["SWS_spectral"]:
        print(f"   omit {row['omitted_cohort_family']:<22} r = {row['r']:+.3f} "
              f"[{row['ci'][0]:+.3f},{row['ci'][1]:+.3f}]")

    print("\nPrecision (E6):")
    for dom, v in P["by_stratum"].items():
        print(f"   {dom:<14} cohorts={v['n_cohorts']}  MDE r={v['mde_r_at_estimated_tau2']:.3f}"
              f"   power@r=.30: {v['power_to_detect_r_0.30']:.2f}"
              f"   power@r=.45: {v['power_to_detect_r_0.45_Berkeley_size']:.2f}")
    print("\nTOST equivalence, N3 stage stratum:")
    for b, t in P["equivalence_tests_SWS_stage"]["tests"].items():
        print(f"   {b:<14} p_TOST = {t['p_tost']:.4f}   equivalent: {t['equivalent_at_05']}")


if __name__ == "__main__":
    main()
