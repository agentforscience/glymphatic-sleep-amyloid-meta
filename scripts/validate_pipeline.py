#!/usr/bin/env python3
"""End-to-end validation that the curated effect-size dataset supports the planned
meta-analysis. NOT the final analysis - a smoke test producing artifacts/pilot_meta.json.

Primary contrast (hypothesis): is the pooled slow-wave-sleep <-> amyloid-PET association
stronger than the total-sleep-time and sleep-efficiency associations?
"""
import json, numpy as np, pandas as pd, os
from pymare import Dataset, estimators

os.makedirs("artifacts", exist_ok=True)
d = pd.read_csv("datasets/meta_analysis/effect_sizes.csv")

R_TYPES = {"pearson_r", "partial_r", "kendall_tau", "r_from_partial_eta2"}
pet = d[d.outcome_class.isin(["deposition_pet", "deposition_pet_change"])].copy()
corr = pet[pet.effect_type.isin(R_TYPES)].copy()
# kendall tau -> pearson r (Greiner's relation)
corr["r"] = np.where(corr.effect_type == "kendall_tau",
                     np.sin(np.pi * corr.burden_aligned_effect / 2),
                     corr.burden_aligned_effect)
corr["z"] = np.arctanh(corr.r.clip(-0.999, 0.999))
corr["v"] = 1.0 / (corr.n - 3)

def pool(sub, label):
    if len(sub) < 2:
        return {"contrast": label, "k": int(len(sub)), "note": "too few effects to pool"}
    ds = Dataset(y=sub.z.values[:, None], v=sub.v.values[:, None])
    est = estimators.VarianceBasedLikelihoodEstimator(method="REML").fit_dataset(ds)
    fe = est.summary().get_fe_stats()
    b, lo, hi = (float(np.ravel(fe[k])[0]) for k in ("est", "ci_l", "ci_u"))
    tau2 = float(np.ravel(est.params_["tau2"])[0])
    w = 1 / sub.v.values
    Q = float(np.sum(w * (sub.z.values - np.average(sub.z.values, weights=w)) ** 2))
    dfree = len(sub) - 1
    I2 = max(0.0, 100 * (Q - dfree) / Q) if Q > 0 else 0.0
    return {"contrast": label, "k": int(len(sub)), "n_total": int(sub.n.sum()),
            "pooled_fisher_z": round(b, 4), "pooled_r": round(float(np.tanh(b)), 4),
            "ci_r": [round(float(np.tanh(lo)), 4), round(float(np.tanh(hi)), 4)],
            "tau2": round(tau2, 5), "Q": round(Q, 3), "I2_pct": round(I2, 1),
            "studies": sorted(sub.study_id.unique().tolist()),
            "cohort_families": sorted(sub.cohort_family.unique().tolist()),
            "n_independent_cohorts": int(sub.cohort_family.nunique()),
            "dependency_warning": ("effects share cohorts - naive pooling overstates precision"
                                   if sub.cohort_family.nunique() < len(sub) else None)}

results = {"generated_by": "scripts/validate_pipeline.py",
           "input": "datasets/meta_analysis/effect_sizes.csv",
           "n_effects_total": int(len(d)), "n_studies_total": int(d.study_id.nunique()),
           "n_amyloid_pet_effects": int(len(pet)),
           "n_correlation_effects_poolable": int(len(corr)), "pooled": []}

for dom, label in [("SWS_spectral", "spectral slow-wave activity vs amyloid PET"),
                   ("SWS_stage",    "N3/SWS stage duration vs amyloid PET"),
                   ("TST",          "total sleep time vs amyloid PET"),
                   ("SE",           "sleep efficiency vs amyloid PET")]:
    results["pooled"].append(pool(corr[corr.exposure_domain == dom], label))

# attenuation-corrected comparison using empirical single-night reliabilities
rel = pd.read_csv("datasets/derived/night_to_night_reliability.csv").set_index("metric")
ATT = {"SWS_stage": "noct_n3_pct_tst", "TST": "noct_tst_min", "SE": "noct_sleep_eff_pct"}
results["measurement_reliability_note"] = (
    "Single-night ICC(2,1) from Sleep-EDF Expanded (n=75 two-night subjects). "
    "Observed correlations are attenuated by sqrt(ICC); dividing by that factor gives a "
    "disattenuated estimate of the true association.")
results["attenuation"] = {
    dom: {"reliability_metric": m,
          "icc": float(rel.loc[m, "icc_2_1"]),
          "attenuation_factor": float(rel.loc[m, "attenuation_factor_sqrt_icc"])}
    for dom, m in ATT.items()}

with open("artifacts/pilot_meta.json", "w") as f:
    json.dump(results, f, indent=2)
print(json.dumps(results, indent=2))
