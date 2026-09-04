#!/usr/bin/env python3
"""
exp3_attenuation.py - E3: the measurement-reliability null model.

Competing explanation for any "metric A beats metric B" finding: metric A may simply be
measured more reliably. Classical attenuation says an observed correlation is
    r_obs = r_true * sqrt(ICC_exposure) * sqrt(ICC_outcome)
so if two exposures have identical TRUE correlations with amyloid, the ratio of their
OBSERVED correlations is sqrt(ICC_a / ICC_b) - the outcome term cancels because the outcome
is the same amyloid PET measure in both.

Empirical single-night test-retest reliabilities come from Sleep-EDF Expanded (75 subjects
with two consecutive nights), derived in the resource-gathering phase:
    N3 % of TST      ICC(2,1) = 0.775   -> attenuation factor 0.880
    sleep efficiency ICC(2,1) = 0.562   -> 0.750
    total sleep time ICC(2,1) = 0.225   -> 0.474

So a single night of PSG measures N3% far more precisely than it measures TST, and the
measurement-only prediction is that N3% should show a correlation about 1.86x larger than
TST even if the underlying biology is identical. The hypothesis's specificity claim therefore
requires the observed N3 advantage to EXCEED 1.86x, not merely to exist.

CRITICAL MODALITY CAVEAT, handled explicitly below: 6 of the 9 TST effects use *self-reported
habitual* sleep duration, not single-night PSG. Self-report is a retrospective trait estimate
with much higher test-retest reliability than one night of PSG, so the Sleep-EDF PSG ICC is
the wrong reliability for those rows. E3 therefore reports a PSG-only restriction as the
like-for-like comparison, and treats the all-modality version as an upper bound on how much
of the gap reliability could explain.

Output: results/e3_attenuation.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metalib as M                                        # noqa: E402
from exp1_pooling import load, aggregate_by_cohort, RHO_PRIMARY  # noqa: E402
from exp2_specificity import RHO_WITHIN_PRIMARY            # noqa: E402

np.random.seed(42)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REL_CSV = os.path.join(ROOT, "datasets", "derived", "night_to_night_reliability.csv")
OUT = os.path.join(ROOT, "results", "e3_attenuation.json")

# Map each exposure domain to the Sleep-EDF metric whose reliability applies to it.
# SWS_spectral has no hypnogram-derivable counterpart (spectral power needs the raw EEG,
# and we hold annotations only), so it is left unmapped and flagged.
DOMAIN_TO_METRIC = {
    "SWS_stage": "noct_n3_pct_tst",
    "TST":       "noct_tst_min",
    "SE":        "noct_sleep_eff_pct",
    "WASO":      "noct_waso_min",
}
# Literature value used only in a clearly-labelled sensitivity arm: test-retest reliability
# of self-reported habitual sleep duration over weeks-to-months is typically r ~ 0.75.
SELF_REPORT_TST_ICC = 0.75


def icc_with_ci(rel: pd.DataFrame, metric: str, level=0.95):
    """ICC point estimate plus a CI derived from the night1-night2 correlation.

    The ICC(2,1) and the between-night Pearson r are near-identical here, so we take the CI
    from the Fisher-z interval of the between-night r with n subjects, and carry its variance
    forward into the disattenuation.
    """
    row = rel.loc[metric]
    icc = float(row["icc_2_1"])
    n = float(row["n_subjects"])
    r_nn = float(row["pearson_r_night1_night2"])
    z = M.r_to_z(r_nn)
    se_z = 1.0 / np.sqrt(n - 3.0)
    crit = 1.959963985
    lo, hi = float(M.z_to_r(z - crit * se_z)), float(M.z_to_r(z + crit * se_z))
    # delta-method variance of r on the original scale
    var_r = ((1 - r_nn ** 2) ** 2) * se_z ** 2
    return {"metric": metric, "icc": icc, "icc_ci": [lo, hi], "var_icc": float(var_r),
            "attenuation_factor_sqrt_icc": float(np.sqrt(icc)), "n_subjects": int(n)}


def disattenuate(z, v, icc, var_icc):
    """Correct a Fisher-z effect for exposure measurement error; propagate ICC uncertainty.

    Works on the r scale via the delta method, then transforms back:
        r_dis = r / sqrt(icc)
        Var(r_dis) ~ Var(r)/icc + r^2 * Var(icc) / (4 * icc^3)
        z_dis = artanh(r_dis),  Var(z_dis) ~ Var(r_dis) / (1 - r_dis^2)^2
    """
    r = float(M.z_to_r(z))
    var_r = ((1 - r ** 2) ** 2) * v
    r_dis = r / np.sqrt(icc)
    r_dis = float(np.clip(r_dis, -0.999, 0.999))
    var_r_dis = var_r / icc + (r ** 2) * var_icc / (4 * icc ** 3)
    z_dis = float(M.r_to_z(r_dis))
    var_z_dis = float(var_r_dis / max((1 - r_dis ** 2) ** 2, 1e-9))
    return z_dis, var_z_dis, r_dis


def pooled_disattenuated(d, domain, iccs, modality_filter=None, rho=RHO_PRIMARY):
    """Pool a stratum after disattenuating each cohort-level effect."""
    sub = d[d.exposure_domain == domain]
    if modality_filter:
        sub = sub[sub.exposure_modality.isin(modality_filter)]
    if len(sub) == 0:
        return None
    agg = aggregate_by_cohort(sub, rho=rho)
    info = iccs[domain]
    zz, vv, rr = [], [], []
    for _, row in agg.iterrows():
        zd, vd, rd = disattenuate(row.z, row.v, info["icc"], info["var_icc"])
        zz.append(zd); vv.append(vd); rr.append(rd)
    pooled = M.re_meta(np.array(zz), np.array(vv), hksj=True)
    pooled["icc_used"] = info["icc"]
    pooled["n_cohorts"] = int(len(agg))
    pooled["per_cohort_disattenuated_r"] = [
        {"cohort_family": c, "r_obs": float(M.z_to_r(z)), "r_disattenuated": float(rd)}
        for c, z, rd in zip(agg.cohort_family, agg.z, rr)]
    return pooled


def paired_disattenuated(d, dom_a, dom_b, iccs, modality_filter=None,
                         rho=RHO_PRIMARY, rho_within=RHO_WITHIN_PRIMARY):
    """Within-cohort paired contrast after disattenuating both exposures."""
    def tab(dom):
        sub = d[d.exposure_domain == dom]
        if modality_filter:
            sub = sub[sub.exposure_modality.isin(modality_filter)]
        if len(sub) == 0:
            return None
        return aggregate_by_cohort(sub, rho=rho).set_index("cohort_family")
    A, B = tab(dom_a), tab(dom_b)
    if A is None or B is None:
        return {"note": "one side empty under this modality filter"}
    shared = sorted(set(A.index) & set(B.index))
    if not shared:
        return {"note": "no shared cohorts under this modality filter",
                "n_shared_cohorts": 0}
    dz, vz, detail = [], [], []
    for cf in shared:
        za, va, ra = disattenuate(A.loc[cf, "z"], A.loc[cf, "v"],
                                  iccs[dom_a]["icc"], iccs[dom_a]["var_icc"])
        zb, vb, rb = disattenuate(B.loc[cf, "z"], B.loc[cf, "v"],
                                  iccs[dom_b]["icc"], iccs[dom_b]["var_icc"])
        dd, vv = M.paired_contrast(za, va, zb, vb, rho_within)
        dz.append(dd); vz.append(vv)
        detail.append({"cohort_family": cf, "r_a_dis": ra, "r_b_dis": rb, "diff_z": dd})
    pooled = M.re_meta(np.array(dz), np.array(vz), hksj=True)
    return {"n_shared_cohorts": len(shared), "per_cohort": detail,
            "pooled_difference": pooled,
            "interpretation": "positive => N3 stage stronger than comparator AFTER "
                              "correcting both for measurement unreliability"}


def main():
    rel = pd.read_csv(REL_CSV).set_index("metric")
    d = load()
    e1 = json.load(open(os.path.join(ROOT, "results", "e1_pooled.json")))

    iccs = {dom: icc_with_ci(rel, met) for dom, met in DOMAIN_TO_METRIC.items()}

    res = {
        "experiment": "E3 - measurement-reliability null model (competing explanation)",
        "reliability_source": ("Sleep-EDF Expanded, 75 subjects x 2 consecutive nights; "
                               "single-night ICC(2,1) per metric"),
        "iccs": iccs,
        "unmapped_domains": {
            "SWS_spectral": ("no reliability estimate derivable from hypnogram annotations "
                             "(spectral power requires the raw EEG, which was not downloaded); "
                             "not disattenuated")},
        "modality_caveat": (
            f"{int((d[d.exposure_domain=='TST'].exposure_modality=='self_report').sum())} of "
            f"{int((d.exposure_domain=='TST').sum())} TST effects and "
            f"{int((d[d.exposure_domain=='SE'].exposure_modality=='self_report').sum())} of "
            f"{int((d.exposure_domain=='SE').sum())} SE effects use self-report, whose "
            "test-retest reliability is NOT the single-night PSG ICC. The PSG-only arm below "
            "is the like-for-like comparison."),
    }

    # ---- 1. The measurement-only prediction -------------------------------------------
    pred = {}
    for comp in ["TST", "SE"]:
        ratio = float(np.sqrt(iccs["SWS_stage"]["icc"] / iccs[comp]["icc"]))
        obs_a = e1["strata"]["SWS_stage"]["primary"]["r"]
        obs_b = e1["strata"][comp]["primary"]["r"]
        obs_ratio = float(obs_a / obs_b) if abs(obs_b) > 1e-6 else np.nan
        pred[f"SWS_stage_over_{comp}"] = {
            "predicted_observed_ratio_if_true_effects_equal": ratio,
            "observed_r_SWS_stage": obs_a, "observed_r_comparator": obs_b,
            "observed_ratio": obs_ratio,
            "excess_over_measurement_prediction": (
                float(obs_ratio / ratio) if np.isfinite(obs_ratio) else None),
            "verdict": ("hypothesis requires observed ratio > "
                        f"{ratio:.2f}; observed ratio is "
                        f"{obs_ratio:.2f}" if np.isfinite(obs_ratio) else "comparator ~0"),
        }
    res["measurement_only_prediction"] = pred

    # ---- 2. Disattenuated stratum estimates -------------------------------------------
    res["disattenuated_pooled"] = {}
    res["disattenuated_pooled_psg_only"] = {}
    for dom in DOMAIN_TO_METRIC:
        p = pooled_disattenuated(d, dom, iccs)
        if p:
            res["disattenuated_pooled"][dom] = {
                "observed_r": e1["strata"][dom]["primary"]["r"],
                "disattenuated_r": p["r"],
                "disattenuated_ci_r": [p["r_ci_low"], p["r_ci_high"]],
                "p": p["p"], "k": p["k"], "icc_used": p["icc_used"],
                "per_cohort": p["per_cohort_disattenuated_r"]}
        p2 = pooled_disattenuated(d, dom, iccs, modality_filter=["PSG_lab"])
        if p2:
            res["disattenuated_pooled_psg_only"][dom] = {
                "disattenuated_r": p2["r"],
                "disattenuated_ci_r": [p2["r_ci_low"], p2["r_ci_high"]],
                "p": p2["p"], "k": p2["k"]}

    # observed PSG-only pooled estimates (no disattenuation) for reference
    res["observed_pooled_psg_only"] = {}
    for dom in ["SWS_stage", "TST", "SE", "SWS_spectral"]:
        sub = d[(d.exposure_domain == dom) & (d.exposure_modality == "PSG_lab")]
        if len(sub):
            agg = aggregate_by_cohort(sub)
            r = M.re_meta(agg.z.values, agg.v.values, hksj=True)
            res["observed_pooled_psg_only"][dom] = {
                "k_effects": int(len(sub)), "n_cohorts": int(len(agg)),
                "r": r["r"], "ci_r": [r["r_ci_low"], r["r_ci_high"]], "p": r["p"],
                "I2": r["I2"]}

    # ---- 3. Disattenuated paired contrasts (the formal test) --------------------------
    res["disattenuated_paired_contrasts"] = {
        "SWS_stage_vs_TST_all_modalities": paired_disattenuated(d, "SWS_stage", "TST", iccs),
        "SWS_stage_vs_SE_all_modalities": paired_disattenuated(d, "SWS_stage", "SE", iccs),
        "SWS_stage_vs_TST_PSG_only": paired_disattenuated(
            d, "SWS_stage", "TST", iccs, modality_filter=["PSG_lab"]),
        "SWS_stage_vs_SE_PSG_only": paired_disattenuated(
            d, "SWS_stage", "SE", iccs, modality_filter=["PSG_lab"]),
    }

    # ---- 4. Self-report-reliability sensitivity ---------------------------------------
    iccs_sr = {k: dict(v) for k, v in iccs.items()}
    iccs_sr["TST"] = {**iccs["TST"], "icc": SELF_REPORT_TST_ICC,
                      "var_icc": iccs["TST"]["var_icc"]}
    res["sensitivity_self_report_reliability"] = {
        "assumed_icc_self_report_TST": SELF_REPORT_TST_ICC,
        "rationale": ("If self-reported habitual sleep duration is in fact a reliable trait "
                      "measure (ICC ~ 0.75) rather than a noisy single-night one (0.225), "
                      "TST's observed correlation needs far less upward correction, which "
                      "makes the N3 gap even harder to explain by measurement error."),
        "predicted_observed_ratio_if_true_effects_equal": float(
            np.sqrt(iccs["SWS_stage"]["icc"] / SELF_REPORT_TST_ICC)),
        "disattenuated_TST_r": pooled_disattenuated(d, "TST", iccs_sr)["r"],
    }

    # ---- 5. What would N3's true correlation have to be? ------------------------------
    tst_dis = res["disattenuated_pooled"]["TST"]["disattenuated_r"]
    needed_obs = tst_dis * np.sqrt(iccs["SWS_stage"]["icc"])
    res["required_effect_for_parity"] = {
        "disattenuated_TST_r": tst_dis,
        "observed_r_N3_needed_for_equal_true_effect": float(needed_obs),
        "observed_r_N3_actual": e1["strata"]["SWS_stage"]["primary"]["r"],
        "shortfall": float(needed_obs - e1["strata"]["SWS_stage"]["primary"]["r"]),
        "note": ("For N3 stage to have the same TRUE association with amyloid as TST, the "
                 "OBSERVED single-night N3 correlation would have to be this large. It is not.")}

    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)
    print(f"[exp3] wrote {OUT}\n")

    print("Measurement-only prediction (if true effects were equal):")
    for k, v in pred.items():
        print(f"  {k}: predicted observed ratio {v['predicted_observed_ratio_if_true_effects_equal']:.2f}x"
              f", actual {v['observed_ratio']:.2f}x")
    print("\nDisattenuated pooled r (all modalities):")
    for k, v in res["disattenuated_pooled"].items():
        print(f"  {k:<10} observed r={v['observed_r']:+.3f} -> true r={v['disattenuated_r']:+.3f} "
              f"[{v['disattenuated_ci_r'][0]:+.3f},{v['disattenuated_ci_r'][1]:+.3f}] (ICC={v['icc_used']:.3f})")
    print("\nObserved pooled r, PSG-measured exposures only (like-for-like):")
    for k, v in res["observed_pooled_psg_only"].items():
        print(f"  {k:<14} k={v['k_effects']} cohorts={v['n_cohorts']}  r={v['r']:+.3f} "
              f"[{v['ci_r'][0]:+.3f},{v['ci_r'][1]:+.3f}]  p={v['p']:.3f}")
    print("\nDisattenuated paired contrasts (positive favours the hypothesis):")
    for k, v in res["disattenuated_paired_contrasts"].items():
        pdp = v.get("pooled_difference")
        if pdp:
            print(f"  {k:<38} dz={pdp['est']:+.3f} [{pdp['ci_low']:+.3f},{pdp['ci_high']:+.3f}] p={pdp['p']:.3f}"
                  f"  (cohorts={v['n_shared_cohorts']})")
        else:
            print(f"  {k:<38} {v.get('note')}")


if __name__ == "__main__":
    main()
