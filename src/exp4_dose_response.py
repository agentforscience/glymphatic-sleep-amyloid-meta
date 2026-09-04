#!/usr/bin/env python3
"""
exp4_dose_response.py - E4: convert the pooled standardised effects into the units the
hypothesis actually states - amyloid PET SUVR (and Centiloid) change per 10 percentage-point
reduction in N3 sleep percentage.

Chain of conversion (each step stated so it can be criticised):
    Delta_amyloid per 10pp N3  =  r_pooled  x  (10 / SD_N3%)  x  SD_amyloid

  r_pooled       from E1, cohort-aggregated random-effects estimate
  SD_N3%         between-person dispersion of N3 % of TST in the relevant age range,
                 derived from Sleep-EDF Expanded and age-weighted to match the amyloid
                 studies' mean ages
  SD_amyloid     between-person dispersion of global amyloid PET in cognitively unimpaired
                 older adults, reconstructed from Jack 2017's A-/A+ medians and IQRs as a
                 two-component mixture with A+ prevalence p

Uncertainty is propagated by Monte Carlo over all three inputs (10,000 draws, seed 42):
  - r_pooled via a t_{k-1} draw on the Fisher-z scale using the HKSJ standard error
  - SD_N3% via the chi-square sampling distribution of a standard deviation
  - SD_amyloid via a Beta prior on A+ prevalence and chi-square sampling of the component SDs

A separate longitudinal version expresses the effect as Centiloids/year, benchmarked against
the empirical accumulation rate in cognitively unimpaired adults (Carvalho 2024: 2.1 +- 5.3
CL/year; Hanseeuw 2021: +3 to +11 CL/year in amyloid-positive individuals).

FEASIBILITY CHECK: the analysis also asks whether a "10% reduction in N3 percentage" is even
attainable in the target population, using the empirical distribution of N3% by age.

Output: results/e4_dose_response.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metalib as M                       # noqa: E402
from exp1_pooling import load             # noqa: E402

SEED = 42
N_MC = 10000
rng = np.random.default_rng(SEED)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "e4_dose_response.json")
SUBJ = os.path.join(ROOT, "datasets", "sleep_edfx", "sleepedf_subject_level.csv")

IQR_TO_SD = 1.349  # IQR = 1.349 * SD for a normal distribution

# Jack 2017 anchors (datasets/meta_analysis/amyloid_pet_reference_scale.csv)
SUVR_ANEG_MED, SUVR_ANEG_IQR = 1.31, (1.26, 1.35)
SUVR_APOS_MED, SUVR_APOS_IQR = 1.57, (1.47, 1.77)
CL_ANEG_MED,  CL_ANEG_IQR = 9.0, (5.0, 12.0)
CL_APOS_MED,  CL_APOS_IQR = 31.0, (23.0, 48.0)
APOS_PREVALENCE = 0.30            # A+ prevalence among cognitively unimpaired older adults
APOS_PREV_BETA = (30.0, 70.0)     # Beta(30,70): mean 0.30, 95% ~ [0.215, 0.394]

# Longitudinal benchmark: Carvalho 2024 Table 1, cognitively unimpaired Mayo MCSA
DELTA_CL_PER_YEAR_MEAN, DELTA_CL_PER_YEAR_SD, DELTA_CL_N = 2.1, 5.3, 64
# Carvalho 2024's own age-equivalence anchor: 1.46 CL/year ~ 6 additional years of age
CL_PER_YEAR_PER_AGE_YEAR = 1.46 / 6.0


def mixture_sd(med_neg, sd_neg, med_pos, sd_pos, p):
    """SD of a two-component mixture with given component means/SDs and mixing weight p."""
    mu = (1 - p) * med_neg + p * med_pos
    var = ((1 - p) * (sd_neg ** 2 + med_neg ** 2) + p * (sd_pos ** 2 + med_pos ** 2)) - mu ** 2
    return float(np.sqrt(max(var, 1e-12))), float(mu)


def sd_draw(sd_hat, n, size):
    """Draw from the sampling distribution of an SD estimated from n observations."""
    df = max(n - 1, 1)
    return sd_hat * np.sqrt(rng.chisquare(df, size) / df)


def age_weighted_n3_sd(subj: pd.DataFrame, study_ages):
    """Between-person SD of N3 % of TST, weighted to the age distribution of the studies.

    Uses one record per subject (first night) to avoid double-counting, restricted to the
    Sleep-EDF age bands that overlap the amyloid studies' mean ages.
    """
    first = subj.sort_values("night_n").groupby(["cohort", "subj_num"], as_index=False).first()
    bands = [(40, 60), (60, 75), (75, 200)]
    out = {}
    for lo, hi in bands:
        g = first[(first.age >= lo) & (first.age < hi)]
        if len(g) >= 5:
            out[f"{lo}-{hi if hi < 200 else '+'}"] = {
                "n": int(len(g)), "mean": float(g.noct_n3_pct_tst.mean()),
                "sd": float(g.noct_n3_pct_tst.std(ddof=1)),
                "median": float(g.noct_n3_pct_tst.median()),
                "pct_at_or_above_10pp": float((g.noct_n3_pct_tst >= 10).mean() * 100)}
    # weight the bands by how many amyloid-study mean ages fall in each
    w = {k: 0 for k in out}
    for a in study_ages:
        for k in out:
            lo = float(k.split("-")[0])
            hi = 200.0 if k.endswith("+") else float(k.split("-")[1])
            if lo <= a < hi:
                w[k] += 1
    tot = sum(w.values())
    if tot == 0:
        w = {k: 1 for k in out}
        tot = len(out)
    sd_w = sum(out[k]["sd"] * w[k] for k in out) / tot
    n_w = int(sum(out[k]["n"] * w[k] for k in out) / tot)
    return {"bands": out, "weights": w, "weighted_sd": float(sd_w),
            "effective_n": n_w}


def mc_dose(pooled, sd_n3_hat, n_n3, scale="SUVR", delta_pp=10.0):
    """Monte-Carlo distribution of the amyloid change per `delta_pp` points of N3%."""
    k = pooled["k"]
    df = max(k - 1, 1)
    # r_pooled draws
    z = pooled["est"] + pooled["se"] * rng.standard_t(df, N_MC)
    r = np.tanh(z)
    # SD_N3 draws
    sd_n3 = sd_draw(sd_n3_hat, n_n3, N_MC)
    # SD_amyloid draws
    p = rng.beta(*APOS_PREV_BETA, N_MC)
    if scale == "SUVR":
        sd_neg = sd_draw((SUVR_ANEG_IQR[1] - SUVR_ANEG_IQR[0]) / IQR_TO_SD, 100, N_MC)
        sd_pos = sd_draw((SUVR_APOS_IQR[1] - SUVR_APOS_IQR[0]) / IQR_TO_SD, 100, N_MC)
        med_neg, med_pos = SUVR_ANEG_MED, SUVR_APOS_MED
    else:
        sd_neg = sd_draw((CL_ANEG_IQR[1] - CL_ANEG_IQR[0]) / IQR_TO_SD, 100, N_MC)
        sd_pos = sd_draw((CL_APOS_IQR[1] - CL_APOS_IQR[0]) / IQR_TO_SD, 100, N_MC)
        med_neg, med_pos = CL_ANEG_MED, CL_APOS_MED
    mu = (1 - p) * med_neg + p * med_pos
    var = ((1 - p) * (sd_neg ** 2 + med_neg ** 2) + p * (sd_pos ** 2 + med_pos ** 2)) - mu ** 2
    sd_amy = np.sqrt(np.maximum(var, 1e-12))

    delta = r * (delta_pp / sd_n3) * sd_amy
    return {
        "scale": scale, "delta_pp_N3": delta_pp,
        "point_estimate": float(np.median(delta)),
        "mean": float(np.mean(delta)),
        "ci_95": [float(np.percentile(delta, 2.5)), float(np.percentile(delta, 97.5))],
        "ci_80": [float(np.percentile(delta, 10)), float(np.percentile(delta, 90))],
        "p_direction_positive": float(np.mean(delta > 0)),
        "sd_amyloid_used_median": float(np.median(sd_amy)),
        "sd_n3_used_median": float(np.median(sd_n3)),
        "n_mc": N_MC,
    }


def main():
    e1 = json.load(open(os.path.join(ROOT, "results", "e1_pooled.json")))
    d = load()
    subj = pd.read_csv(SUBJ)

    study_ages = d.groupby("cohort_family").mean_age.mean().dropna().tolist()
    n3 = age_weighted_n3_sd(subj, study_ages)

    res = {
        "experiment": "E4 - dose-response conversion (hypothesis clause ii)",
        "seed": SEED, "n_monte_carlo": N_MC,
        "conversion_formula": "Delta_amyloid per 10pp N3 = r_pooled * (10 / SD_N3%) * SD_amyloid",
        "inputs": {
            "n3_dispersion": n3,
            "study_mean_ages_by_cohort": {k: float(v) for k, v in
                                          d.groupby("cohort_family").mean_age.mean()
                                          .dropna().items()},
            "amyloid_scale_anchors": {
                "source": "Jack 2017 (datasets/meta_analysis/amyloid_pet_reference_scale.csv)",
                "A_negative_SUVR_median_IQR": [SUVR_ANEG_MED, list(SUVR_ANEG_IQR)],
                "A_positive_SUVR_median_IQR": [SUVR_APOS_MED, list(SUVR_APOS_IQR)],
                "A_negative_CL_median_IQR": [CL_ANEG_MED, list(CL_ANEG_IQR)],
                "A_positive_CL_median_IQR": [CL_APOS_MED, list(CL_APOS_IQR)],
                "A_positive_prevalence_prior": f"Beta{APOS_PREV_BETA} (mean {APOS_PREVALENCE})"},
        },
    }

    sd_suvr, mu_suvr = mixture_sd(
        SUVR_ANEG_MED, (SUVR_ANEG_IQR[1] - SUVR_ANEG_IQR[0]) / IQR_TO_SD,
        SUVR_APOS_MED, (SUVR_APOS_IQR[1] - SUVR_APOS_IQR[0]) / IQR_TO_SD, APOS_PREVALENCE)
    sd_cl, mu_cl = mixture_sd(
        CL_ANEG_MED, (CL_ANEG_IQR[1] - CL_ANEG_IQR[0]) / IQR_TO_SD,
        CL_APOS_MED, (CL_APOS_IQR[1] - CL_APOS_IQR[0]) / IQR_TO_SD, APOS_PREVALENCE)
    res["derived_amyloid_dispersion"] = {
        "SUVR_mixture_mean": mu_suvr, "SUVR_mixture_sd": sd_suvr,
        "Centiloid_mixture_mean": mu_cl, "Centiloid_mixture_sd": sd_cl,
        "A_neg_to_A_pos_gap_SUVR": SUVR_APOS_MED - SUVR_ANEG_MED,
        "A_neg_to_A_pos_gap_CL": CL_APOS_MED - CL_ANEG_MED,
        "internal_consistency_check_CL_per_SUVR": float(
            (CL_APOS_MED - CL_ANEG_MED) / (SUVR_APOS_MED - SUVR_ANEG_MED))}

    # ---- dose estimates for each stratum ------------------------------------------
    res["dose_per_10pp_N3_equivalent"] = {}
    for dom in ["SWS_stage", "TST", "SE", "SWS_spectral"]:
        pooled = e1["strata"][dom]["primary"]
        suvr = mc_dose(pooled, n3["weighted_sd"], n3["effective_n"], "SUVR")
        cl = mc_dose(pooled, n3["weighted_sd"], n3["effective_n"], "Centiloid")
        entry = {
            "pooled_r": pooled["r"], "k_cohorts": pooled["k"],
            "delta_SUVR": suvr, "delta_Centiloid": cl,
            "pct_of_A_minus_to_A_plus_SUVR_gap": float(
                100 * suvr["point_estimate"] / (SUVR_APOS_MED - SUVR_ANEG_MED)),
            "pct_of_A_minus_to_A_plus_CL_gap": float(
                100 * cl["point_estimate"] / (CL_APOS_MED - CL_ANEG_MED)),
        }
        if dom != "SWS_stage":
            entry["note"] = ("expressed per 10 percentage points of N3 only to put all "
                             "strata on one axis; for this exposure the natural unit differs")
        res["dose_per_10pp_N3_equivalent"][dom] = entry

    # ---- longitudinal framing ------------------------------------------------------
    lon = e1["strata"]["SWS_stage"].get("by_design", {}).get("longitudinal")
    res["longitudinal_framing"] = {
        "benchmark_accumulation": {
            "mean_CL_per_year": DELTA_CL_PER_YEAR_MEAN, "sd_CL_per_year": DELTA_CL_PER_YEAR_SD,
            "n": DELTA_CL_N,
            "source": "Carvalho 2024, cognitively unimpaired Mayo MCSA, Table 1",
            "A_positive_range_CL_per_year": "+3 to +11 (Hanseeuw 2021)"},
        "age_equivalence_anchor": {
            "CL_per_year_per_extra_year_of_age": CL_PER_YEAR_PER_AGE_YEAR,
            "source": "Carvalho 2024: 1.46 CL/year effect described as equivalent to 6 "
                      "additional years of baseline age"},
    }
    if lon and lon["pooled"].get("k"):
        pl = lon["pooled"]
        delta_cl_yr = pl["r"] * (10.0 / n3["weighted_sd"]) * DELTA_CL_PER_YEAR_SD
        res["longitudinal_framing"]["SWS_stage_delta_CL_per_year_per_10pp"] = {
            "point_estimate": float(delta_cl_yr),
            "ci_from_r_ci": [float(pl["r_ci_low"] * (10 / n3["weighted_sd"]) * DELTA_CL_PER_YEAR_SD),
                             float(pl["r_ci_high"] * (10 / n3["weighted_sd"]) * DELTA_CL_PER_YEAR_SD)],
            "equivalent_extra_years_of_age": float(delta_cl_yr / CL_PER_YEAR_PER_AGE_YEAR),
            "k_cohorts": pl["k"],
            "caveat": "single cohort family contributes the longitudinal N3-stage effect"}

    # ---- feasibility of a 10-percentage-point reduction ----------------------------
    first = subj.sort_values("night_n").groupby(["cohort", "subj_num"], as_index=False).first()
    older = first[first.age >= 60]
    res["feasibility_of_the_stated_dose"] = {
        "question": ("Is a '10% reduction in N3 sleep percentage' attainable in the "
                     "population the hypothesis is about?"),
        "n_subjects_age_60plus": int(len(older)),
        "mean_N3_pct": float(older.noct_n3_pct_tst.mean()),
        "median_N3_pct": float(older.noct_n3_pct_tst.median()),
        "sd_N3_pct": float(older.noct_n3_pct_tst.std(ddof=1)),
        "pct_with_N3_at_least_10pp": float((older.noct_n3_pct_tst >= 10).mean() * 100),
        "pct_with_N3_below_5pp": float((older.noct_n3_pct_tst < 5).mean() * 100),
        "interpretation": None,  # filled below
    }
    f = res["feasibility_of_the_stated_dose"]
    f["interpretation"] = (
        f"Mean N3 is {f['mean_N3_pct']:.1f}% of TST in adults aged 60+ and "
        f"{100 - f['pct_with_N3_at_least_10pp']:.0f}% of them already have less than 10 "
        "percentage points of N3 in total. An absolute 10-percentage-point reduction is "
        "therefore not available as an exposure contrast for most of the target population; "
        "the stated dose spans roughly the entire observable range rather than a small "
        "perturbation within it.")
    # ratio of the stated dose to the population SD
    f["stated_dose_in_SD_units"] = float(10.0 / n3["weighted_sd"])

    with open(OUT, "w") as f_out:
        json.dump(res, f_out, indent=2, default=str)
    print(f"[exp4] wrote {OUT}\n")

    print(f"N3% dispersion (age-weighted to study ages): SD = {n3['weighted_sd']:.2f} pp "
          f"(effective n = {n3['effective_n']})")
    print(f"  => a 10-pp change is {10 / n3['weighted_sd']:.2f} SD of the exposure")
    print(f"Amyloid dispersion: SUVR SD = {sd_suvr:.3f}, Centiloid SD = {sd_cl:.1f} "
          f"(A-/A+ gap = {SUVR_APOS_MED - SUVR_ANEG_MED:.2f} SUVR / "
          f"{CL_APOS_MED - CL_ANEG_MED:.0f} CL)\n")
    print(f"{'stratum':<14}{'ΔSUVR/10pp':>14}{'95% CI':>26}{'ΔCL/10pp':>12}{'% of A-/A+ gap':>17}")
    print("-" * 84)
    for dom, v in res["dose_per_10pp_N3_equivalent"].items():
        s, c = v["delta_SUVR"], v["delta_Centiloid"]
        print(f"{dom:<14}{s['point_estimate']:>+14.4f}"
              f"   [{s['ci_95'][0]:+.4f}, {s['ci_95'][1]:+.4f}]"
              f"{c['point_estimate']:>12.2f}{v['pct_of_A_minus_to_A_plus_CL_gap']:>16.1f}%")
    print("\nFeasibility:", f["interpretation"])


if __name__ == "__main__":
    main()
