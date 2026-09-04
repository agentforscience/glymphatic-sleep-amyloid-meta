#!/usr/bin/env python3
"""
exp2_specificity.py - E2: is the N3-stage association *stronger* than the total-sleep-time
and sleep-efficiency associations?

This is clause (iii) of the hypothesis, and it is a comparative claim, so it cannot be read
off any single pooled estimate. Two complementary tests are run:

  (a) BETWEEN-STRATUM: Wald contrast of the two pooled estimates from E1. Simple, but the
      strata share cohorts (Berkeley, AgeWell and MARS each contribute to several), so the
      independence assumption is violated and the test is anti-conservative.

  (b) WITHIN-STUDY PAIRED (primary): for every cohort family that reports BOTH exposures
      against amyloid PET in the SAME participants, form the difference
      d = z_stage - z_comparator with Var(d) = v_a + v_b - 2*rho_w*sqrt(v_a*v_b), then pool
      the differences with REML+HKSJ. This removes all between-cohort confounding from the
      exact claim the hypothesis makes and is, to our knowledge, not done in this literature.

Three confirmatory contrasts (Holm-corrected): stage vs TST, stage vs SE, stage vs spectral.

Output: results/e2_specificity.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metalib as M                        # noqa: E402
from exp1_pooling import load, aggregate_by_cohort, RHO_PRIMARY  # noqa: E402

np.random.seed(42)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "e2_specificity.json")

RHO_WITHIN_PRIMARY = 0.5          # correlation between two effects from the same sample
RHO_WITHIN_SENS = [0.0, 0.3, 0.7]

CONTRASTS = [
    ("SWS_stage", "TST",          "N3 stage vs total sleep time (hypothesis clause iii)"),
    ("SWS_stage", "SE",           "N3 stage vs sleep efficiency (hypothesis clause iii)"),
    ("SWS_stage", "SWS_spectral", "N3 stage vs spectral slow-wave activity (dissociation)"),
]


def cohort_table(d: pd.DataFrame, domain: str, rho=RHO_PRIMARY) -> pd.DataFrame:
    """Cohort-family-level synthetic effects for one exposure domain."""
    sub = d[d.exposure_domain == domain]
    if len(sub) == 0:
        return pd.DataFrame(columns=["cohort_family", "z", "v", "m_effects", "n", "studies"])
    return aggregate_by_cohort(sub, rho=rho)


def paired_test(d, dom_a, dom_b, rho_within=RHO_WITHIN_PRIMARY, rho_agg=RHO_PRIMARY):
    """Pool within-cohort differences z_a - z_b over cohorts reporting both exposures."""
    A = cohort_table(d, dom_a, rho_agg).set_index("cohort_family")
    B = cohort_table(d, dom_b, rho_agg).set_index("cohort_family")
    shared = sorted(set(A.index) & set(B.index))
    rows = []
    for cf in shared:
        dz, vz = M.paired_contrast(A.loc[cf, "z"], A.loc[cf, "v"],
                                   B.loc[cf, "z"], B.loc[cf, "v"], rho_within)
        rows.append({"cohort_family": cf, "n": int(max(A.loc[cf, "n"], B.loc[cf, "n"])),
                     "z_a": float(A.loc[cf, "z"]), "z_b": float(B.loc[cf, "z"]),
                     "r_a": float(M.z_to_r(A.loc[cf, "z"])),
                     "r_b": float(M.z_to_r(B.loc[cf, "z"])),
                     "diff_z": dz, "var_diff": vz,
                     "studies_a": A.loc[cf, "studies"], "studies_b": B.loc[cf, "studies"]})
    out = {"contrast": f"{dom_a} - {dom_b}", "rho_within_assumed": rho_within,
           "n_shared_cohorts": len(shared), "per_cohort": rows}
    if len(rows) == 0:
        out["note"] = "no cohort reports both exposures against amyloid PET"
        return out
    dz = np.array([r["diff_z"] for r in rows])
    vz = np.array([r["var_diff"] for r in rows])
    pooled = M.re_meta(dz, vz, hksj=True)
    out["pooled_difference"] = pooled
    out["pooled_difference"]["interpretation"] = (
        "positive => the N3-stage association is STRONGER than the comparator")
    out["permutation"] = M.permutation_p(dz, vz)
    out["n_cohorts_favouring_stage"] = int(np.sum(dz > 0))
    out["n_cohorts_favouring_comparator"] = int(np.sum(dz < 0))
    return out


def main():
    d = load()
    e1 = json.load(open(os.path.join(ROOT, "results", "e1_pooled.json")))

    res = {
        "experiment": "E2 - specificity contrasts (hypothesis clause iii)",
        "primary_test": "within-study paired difference of Fisher z, pooled REML+HKSJ",
        "rho_within_primary": RHO_WITHIN_PRIMARY,
        "contrasts": {},
    }

    raw_p = {}
    for a, b, label in CONTRASTS:
        key = f"{a}_vs_{b}"
        entry = {"label": label}
        # (a) between-stratum
        entry["between_stratum"] = M.between_stratum_contrast(
            e1["strata"][a]["primary"], e1["strata"][b]["primary"])
        entry["between_stratum"]["caveat"] = (
            "strata share cohort families; this test assumes independence and is "
            "anti-conservative. Prefer the paired result.")
        # (b) within-study paired (primary)
        entry["paired"] = paired_test(d, a, b)
        # rho_within sensitivity
        entry["paired_rho_sensitivity"] = {
            f"rho={rw}": (paired_test(d, a, b, rho_within=rw).get("pooled_difference") or {})
            for rw in RHO_WITHIN_SENS}
        res["contrasts"][key] = entry
        pd_ = entry["paired"].get("pooled_difference")
        raw_p[key] = pd_["p"] if pd_ else np.nan

    # Holm correction across the three confirmatory contrasts
    valid = {k: v for k, v in raw_p.items() if np.isfinite(v)}
    order = sorted(valid, key=lambda k: valid[k])
    m = len(order)
    holm, running = {}, 0.0
    for i, k in enumerate(order):
        adj = min(1.0, (m - i) * valid[k])
        running = max(running, adj)
        holm[k] = running
    res["holm_adjusted_p_paired"] = holm
    res["multiplicity_note"] = (
        "Holm step-down across the three preregistered confirmatory contrasts, "
        "applied to the within-study paired p-values (the primary test).")

    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=str)
    print(f"[exp2] wrote {OUT}\n")

    for k, e in res["contrasts"].items():
        p = e["paired"].get("pooled_difference")
        print(f"--- {k}: {e['label']}")
        print(f"    shared cohorts: {e['paired']['n_shared_cohorts']}"
              f"  ({e['paired'].get('n_cohorts_favouring_stage','-')} favour N3, "
              f"{e['paired'].get('n_cohorts_favouring_comparator','-')} favour comparator)")
        if p:
            print(f"    paired dz = {p['est']:+.3f}  95% CI [{p['ci_low']:+.3f}, "
                  f"{p['ci_high']:+.3f}]  p = {p['p']:.3f}  (Holm p = {holm.get(k, float('nan')):.3f})"
                  f"  perm p = {e['paired']['permutation']['p_permutation']:.3f}")
        bs = e["between_stratum"]
        print(f"    between-stratum dz = {bs['diff_z']:+.3f}  p = {bs['p']:.3f} "
              f"(anti-conservative)\n")


if __name__ == "__main__":
    main()
