#!/usr/bin/env python3
"""
prepare_data.py - Harmonise the extracted effect sizes onto a single Fisher-z scale.

Implements the conversion rules preregistered in planning.md section A.2. Every row of the
input keeps its identity in the output, with explicit columns recording which rule was
applied and whether the row is usable in the primary model. Nothing is silently dropped.

Outputs:
    results/harmonised_effects.csv   one row per input effect, with z, v and provenance
    results/qc_report.json           data-quality checks and conversion audit
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metalib as M  # noqa: E402

SEED = 42
np.random.seed(SEED)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_CSV = os.path.join(ROOT, "datasets", "meta_analysis", "effect_sizes.csv")
OUT_CSV = os.path.join(ROOT, "results", "harmonised_effects.csv")
QC_JSON = os.path.join(ROOT, "results", "qc_report.json")

# Effect types already on the correlation scale
DIRECT_R = {"pearson_r", "partial_r", "r_from_partial_eta2"}
# Effect types requiring test-statistic recovery from (p, df)
RECOVER = {"unstd_beta", "unstd_beta_per_SD", "unstd_beta_per_hour",
           "unstd_beta_tertile", "unstd_beta_per_ordinal_category",
           "unstd_beta_group_contrast"}

# Number of covariates implied by each `adjustment` string. Explicit rather than parsed,
# so the count is auditable. `approx=True` marks strings ending in "etc" where the true
# covariate count is not fully recoverable from the paper's text.
COVARIATE_COUNTS = {
    "unadjusted":                                  (0, False),
    "age_sex":                                     (2, False),
    "age_sex_APOE":                                (3, False),
    "age_sex_education":                           (3, False),
    "age_sex_time":                                (3, False),
    "age_APOE4_baseline_amyloid":                  (3, False),
    "age_sex_apnea_interval":                      (4, False),
    "age_sex_APOE_education":                      (4, False),
    "age_sex_education_AHI":                       (4, False),
    "age_sex_APOE_time":                           (4, False),
    "age_sex_education_BMI_etc":                   (5, True),
    "age_sex_APOE4_DM_smoking_BMI_antidep_MMSE":   (8, False),
    "fully_adjusted":                              (8, True),
}
DEFAULT_COVARIATES = 4  # fallback if an adjustment string is unseen; flagged in output


def covariate_count(adjustment: str):
    """Return (n_covariates, is_approximate) for an adjustment description."""
    if not isinstance(adjustment, str):
        return DEFAULT_COVARIATES, True
    return COVARIATE_COUNTS.get(adjustment.strip(), (DEFAULT_COVARIATES, True))


def harmonise(df: pd.DataFrame) -> pd.DataFrame:
    """Apply the preregistered conversion rules; return the augmented dataframe."""
    rows = []
    for _, row in df.iterrows():
        etype = row["effect_type"]
        eff = row["burden_aligned_effect"]
        n = float(row["n"])
        ncov, approx = covariate_count(row.get("adjustment"))
        rec = {
            "conversion_rule": None, "r": np.nan, "z": np.nan, "v": np.nan,
            "df_used": np.nan, "n_covariates": ncov, "df_approximate": approx,
            "is_imputed_null": False, "usable": False, "usable_reason": None,
        }

        if etype in DIRECT_R:
            rec["conversion_rule"] = "direct_r"
            rec["r"] = float(eff)
            rec["df_used"] = n - 2
            rec["usable"] = True

        elif etype == "kendall_tau":
            rec["conversion_rule"] = "kendall_tau_to_r_sin"
            rec["r"] = float(M.tau_to_r(eff))
            rec["df_used"] = n - 2
            rec["usable"] = True

        elif etype in RECOVER:
            p = row.get("p_value")
            if pd.isna(p):
                rec["conversion_rule"] = "recover_from_p_FAILED_no_p"
                rec["usable_reason"] = "unstandardised beta with no reported p-value"
            else:
                dfree = n - 1.0 - ncov
                rec["conversion_rule"] = "recover_r_from_p_and_df"
                rec["df_used"] = dfree
                rec["r"] = M.r_from_p_and_df(p, dfree, np.sign(eff) if eff != 0 else 1.0)
                rec["usable"] = True

        elif etype == "narrative_null":
            rec["conversion_rule"] = "imputed_null_z0"
            rec["r"] = 0.0
            rec["df_used"] = n - 2
            rec["is_imputed_null"] = True
            rec["usable"] = True

        else:
            rec["conversion_rule"] = f"unhandled:{etype}"
            rec["usable_reason"] = f"no conversion rule for effect_type={etype}"

        if rec["usable"]:
            if n <= 4:
                rec["usable"] = False
                rec["usable_reason"] = "n <= 4, Fisher-z variance undefined"
            else:
                rec["z"] = float(M.r_to_z(rec["r"]))
                rec["v"] = float(M.var_z(n))

        rows.append(rec)

    aug = pd.DataFrame(rows, index=df.index)
    return pd.concat([df, aug], axis=1)


def qc_checks(df: pd.DataFrame, out: pd.DataFrame) -> dict:
    """Data-quality checks run before any modelling."""
    pet = out[out.outcome_class.isin(["deposition_pet", "deposition_pet_change"])]
    qc = {
        "seed": SEED,
        "input_file": os.path.relpath(IN_CSV, ROOT),
        "n_effects": int(len(df)),
        "n_studies": int(df.study_id.nunique()),
        "n_cohort_families": int(df.cohort_family.nunique()),
        "n_participants_unique_cohort_max": int(
            df.groupby("cohort_family").n.max().sum()),
        "duplicate_effect_ids": df.effect_id[df.effect_id.duplicated()].tolist(),
        "missing_n": int(df.n.isna().sum()),
        "missing_burden_aligned": int(df.burden_aligned_effect.isna().sum()),
        "missing_provenance": int(df.provenance.isna().sum()),
        "effects_with_r_out_of_range": out.effect_id[
            out.r.abs() > 1].tolist(),
        "effect_type_counts": df.effect_type.value_counts().to_dict(),
        "outcome_class_counts": df.outcome_class.value_counts().to_dict(),
        "conversion_rule_counts": out.conversion_rule.value_counts().to_dict(),
        "n_usable": int(out.usable.sum()),
        "n_unusable": int((~out.usable).sum()),
        "unusable_detail": out.loc[~out.usable, ["effect_id", "effect_type",
                                                 "usable_reason"]].to_dict("records"),
        "n_df_approximate": int(out.df_approximate.sum()),
        "amyloid_pet": {
            "n_effects": int(len(pet)),
            "n_usable": int(pet.usable.sum()),
            "n_imputed_null": int(pet.is_imputed_null.sum()),
            "n_recovered_from_p": int(
                (pet.conversion_rule == "recover_r_from_p_and_df").sum()),
            "by_exposure_domain": pet.groupby("exposure_domain").agg(
                k=("effect_id", "size"),
                k_usable=("usable", "sum"),
                n_cohorts=("cohort_family", "nunique"),
            ).reset_index().to_dict("records"),
        },
        "cohort_family_effect_counts": df.cohort_family.value_counts().to_dict(),
    }
    # Sanity assertion: signs of burden_aligned effects match the recorded z signs
    mism = out[(out.usable) & (~out.is_imputed_null) &
               (np.sign(out.burden_aligned_effect) != np.sign(out.r)) &
               (out.burden_aligned_effect != 0)]
    qc["sign_mismatches"] = mism.effect_id.tolist()
    return qc


def main():
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    df = pd.read_csv(IN_CSV)
    out = harmonise(df)
    qc = qc_checks(df, out)

    # Hard assertions - fail loudly rather than producing a silently wrong analysis
    assert not qc["duplicate_effect_ids"], f"duplicate effect ids: {qc['duplicate_effect_ids']}"
    assert not qc["effects_with_r_out_of_range"], "converted |r| > 1"
    assert not qc["sign_mismatches"], f"sign inversion in conversion: {qc['sign_mismatches']}"
    assert qc["missing_n"] == 0, "missing sample sizes"

    out.to_csv(OUT_CSV, index=False)
    with open(QC_JSON, "w") as f:
        json.dump(qc, f, indent=2, default=str)

    print(f"[prepare_data] wrote {OUT_CSV}  ({len(out)} rows, {int(out.usable.sum())} usable)")
    print(f"[prepare_data] wrote {QC_JSON}")
    print(json.dumps({k: qc[k] for k in
                      ["n_effects", "n_studies", "n_cohort_families", "n_usable",
                       "conversion_rule_counts", "amyloid_pet"]}, indent=2, default=str))


if __name__ == "__main__":
    main()
