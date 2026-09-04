#!/usr/bin/env python3
"""
exp1_pooling.py - E1: stratified random-effects pooling of sleep-metric <-> amyloid-PET
associations, with cohort-family dependency handling.

Answers clause (i) of the hypothesis: is reduced N3 stage percentage associated with higher
amyloid-beta burden? And provides the four stratum estimates that E2 contrasts.

Three estimators are reported per stratum, in increasing order of defensibility:
  naive   - every extracted effect treated as independent (WRONG, shown for contrast only)
  agg     - one synthetic effect per cohort family (Borenstein correlated-effects), REML+HKSJ
  rve     - cluster-robust variance on cohort family (only interpretable with >= 4 clusters)

Output: results/e1_pooled.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metalib as M  # noqa: E402

np.random.seed(42)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARM = os.path.join(ROOT, "results", "harmonised_effects.csv")
OUT = os.path.join(ROOT, "results", "e1_pooled.json")

PET = ["deposition_pet", "deposition_pet_change"]
RHO_PRIMARY = 0.6

STRATA = {
    "SWS_stage":    "N3 / slow-wave-sleep stage amount (the hypothesis's exposure)",
    "SWS_spectral": "Low-frequency spectral slow-wave activity",
    "TST":          "Total sleep time (comparator 1)",
    "SE":           "Sleep efficiency (comparator 2)",
    "WASO":         "Wake after sleep onset (secondary)",
    "other_stage":  "Other sleep stages: N1/N2/REM amount (secondary)",
}


def load(include_imputed_nulls=True, include_recovered=True):
    """Load harmonised amyloid-PET effects under a given inclusion rule."""
    d = pd.read_csv(HARM)
    d = d[d.outcome_class.isin(PET) & d.usable].copy()
    if not include_imputed_nulls:
        d = d[~d.is_imputed_null]
    if not include_recovered:
        d = d[d.conversion_rule != "recover_r_from_p_and_df"]
    return d


def aggregate_by_cohort(sub: pd.DataFrame, rho=RHO_PRIMARY) -> pd.DataFrame:
    """Collapse a stratum to one synthetic effect per cohort family."""
    rows = []
    for cf, g in sub.groupby("cohort_family", sort=True):
        ybar, vbar = M.aggregate_correlated(g.z.values, g.v.values, rho=rho)
        rows.append({
            "cohort_family": cf,
            "z": ybar, "v": vbar, "m_effects": len(g),
            "n": int(g.n.max()),
            "studies": "|".join(sorted(g.study_id.unique())),
            "year_min": int(g.year.min()), "year_max": int(g.year.max()),
            "any_imputed_null": bool(g.is_imputed_null.any()),
            "any_recovered": bool((g.conversion_rule == "recover_r_from_p_and_df").any()),
            "designs": "|".join(sorted(g.design.unique())),
        })
    return pd.DataFrame(rows)


def analyse_stratum(sub: pd.DataFrame, label: str, rho=RHO_PRIMARY) -> dict:
    """Run naive, cohort-aggregated and cluster-robust pooling for one stratum."""
    res = {"exposure_domain": label, "description": STRATA.get(label, ""),
           "k_effects": int(len(sub)),
           "n_cohort_families": int(sub.cohort_family.nunique()),
           "cohort_families": sorted(sub.cohort_family.unique().tolist()),
           "studies": sorted(sub.study_id.unique().tolist()),
           "n_participants_sum_over_cohorts": int(
               sub.groupby("cohort_family").n.max().sum()),
           "k_imputed_null": int(sub.is_imputed_null.sum()),
           "k_recovered_from_p": int(
               (sub.conversion_rule == "recover_r_from_p_and_df").sum()),
           }
    if len(sub) == 0:
        res["note"] = "empty stratum"
        return res

    res["naive"] = M.re_meta(sub.z.values, sub.v.values, hksj=True)

    agg = aggregate_by_cohort(sub, rho=rho)
    res["cohort_aggregated_inputs"] = agg.round(5).to_dict("records")
    res["primary"] = M.re_meta(agg.z.values, agg.v.values, hksj=True)
    res["primary"]["rho_assumed"] = rho
    res["primary"]["estimator"] = (
        f"cohort-family aggregation (rho={rho}) + REML random effects + HKSJ")
    res["permutation"] = M.permutation_p(agg.z.values, agg.v.values)

    if sub.cohort_family.nunique() >= 2:
        res["rve"] = M.rve_meta(sub.z.values, sub.v.values,
                                sub.cohort_family.values, rho=rho)
    # design subgroups
    res["by_design"] = {}
    for des, g in sub.groupby("design"):
        a = aggregate_by_cohort(g, rho=rho)
        res["by_design"][des] = {
            "k_effects": int(len(g)), "n_cohorts": int(g.cohort_family.nunique()),
            "pooled": M.re_meta(a.z.values, a.v.values, hksj=True)}
    return res


def main():
    d = load()
    out = {
        "experiment": "E1 - stratified pooling of sleep metrics vs amyloid PET",
        "analysis_set": "outcome_class in {deposition_pet, deposition_pet_change}, usable effects",
        "n_effects": int(len(d)), "n_studies": int(d.study_id.nunique()),
        "n_cohort_families": int(d.cohort_family.nunique()),
        "rho_primary": RHO_PRIMARY,
        "estimator_note": (
            "Primary = one synthetic effect per cohort family (Borenstein correlated-effects "
            "variance, rho=0.6), then REML random-effects with Knapp-Hartung-Sidik-Jonkman "
            "small-sample standard errors and a t_{k-1} reference distribution. 'naive' pools "
            "every extracted effect as if independent and is reported only to quantify how "
            "much dependency handling matters."),
        "strata": {},
    }
    for dom in STRATA:
        out["strata"][dom] = analyse_stratum(d[d.exposure_domain == dom], dom)

    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"[exp1] wrote {OUT}")

    # console summary
    hdr = f"{'stratum':<14}{'k':>3}{'coh':>5}{'  primary r':>13}{'  95% CI':>20}{'  p':>8}{'  I2':>7}{'  naive r':>10}"
    print(hdr); print("-" * len(hdr))
    for dom, r in out["strata"].items():
        if "primary" not in r:
            print(f"{dom:<14}{r['k_effects']:>3}  (empty)"); continue
        p = r["primary"]
        print(f"{dom:<14}{r['k_effects']:>3}{r['n_cohort_families']:>5}"
              f"{p['r']:>13.3f}   [{p['r_ci_low']:>6.3f},{p['r_ci_high']:>6.3f}]"
              f"{p['p']:>8.3f}{p['I2']:>7.1f}{r['naive']['r']:>10.3f}")


if __name__ == "__main__":
    main()
