#!/usr/bin/env python3
"""
make_figures.py - all figures for the report.

Design follows the project's data-viz guidance: form chosen before color; a single
categorical hue where there is one series (no legend box needed - the title names it);
blue as slot 1 and orange as slot 2 where two entities must be distinguished (that pair
clears the all-pairs CVD floors); recessive grid and axes; direct labels rather than a
number on every mark; text in ink tokens rather than series color.

Every figure also has a machine-readable table companion written to results/figure_data/,
so identity is never carried by color alone.

Outputs: figures/*.png (300 dpi) and results/figure_data/*.csv
"""
from __future__ import annotations

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metalib as M                                        # noqa: E402
from exp1_pooling import load, aggregate_by_cohort         # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
FDAT = os.path.join(ROOT, "results", "figure_data")
os.makedirs(FIG, exist_ok=True)
os.makedirs(FDAT, exist_ok=True)

# --- palette (reference instance from the data-viz guidance) -------------------------
SURFACE = "#fcfcfb"
INK1, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"   # blue, orange, aqua, yellow
GRID = "#e4e3df"
CRITICAL, GOOD = "#d03b3b", "#0ca30c"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK1,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlecolor": INK1,
    "grid.color": GRID, "grid.linewidth": 0.8, "figure.dpi": 110, "savefig.dpi": 300,
    "axes.spines.top": False, "axes.spines.right": False,
})

LABELS = {
    "SWS_stage": "N3 / SWS stage amount",
    "SWS_spectral": "Slow-wave activity (spectral)",
    "TST": "Total sleep time",
    "SE": "Sleep efficiency",
}
MAIN = ["SWS_stage", "SWS_spectral", "TST", "SE"]

# Short display names so forest-plot tick labels do not spill into the neighbouring panel.
SHORT = {"BerkeleyAgingCohort": "Berkeley", "AgeWell_Caen": "Age-Well", "MayoMCSA": "Mayo MCSA",
         "NYU_Barcelona": "NYU-Barc", "AU_community": "AU-community", "NYU_CBH": "NYU-CBH",
         "RomeTorVergata": "Rome-TV"}


def short(name):
    """Compact cohort-family label for axis ticks."""
    return SHORT.get(name, name)


def save(fig, name):
    path = os.path.join(FIG, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  [fig] {os.path.relpath(path, ROOT)}")


def diamond(ax, y, lo, hi, est, color, h=0.22):
    ax.add_patch(Polygon([[lo, y], [est, y + h], [hi, y], [est, y - h]],
                         closed=True, facecolor=color, edgecolor=SURFACE, lw=1.2, zorder=5))


# =====================================================================================
def fig1_forest(d, e1):
    """F1 - forest plot: cohort-family effects and pooled estimate, one panel per stratum."""
    fig, axes = plt.subplots(1, 4, figsize=(17, 5.4), sharex=True)
    fig.subplots_adjust(wspace=0.62)
    rows_out = []
    for ax, dom in zip(axes, MAIN):
        sub = d[d.exposure_domain == dom]
        agg = aggregate_by_cohort(sub).sort_values("z")
        pooled = e1["strata"][dom]["primary"]
        ys = np.arange(len(agg))[::-1]
        r = M.z_to_r(agg.z.values)
        lo = M.z_to_r(agg.z.values - 1.96 * np.sqrt(agg.v.values))
        hi = M.z_to_r(agg.z.values + 1.96 * np.sqrt(agg.v.values))
        # marker area proportional to precision (inverse variance)
        size = 40 + 700 * (1 / agg.v.values) / max((1 / agg.v.values).max(), 1e-9)
        accent = dom == "SWS_stage"
        col = S2 if accent else S1
        ax.hlines(ys, lo, hi, color=col, lw=2, alpha=0.85, zorder=3)
        ax.scatter(r, ys, s=size, color=col, edgecolor=SURFACE, lw=1.2, zorder=4)
        ax.axvline(0, color=INK3, lw=1, zorder=1)
        diamond(ax, -1.5, pooled["r_ci_low"], pooled["r_ci_high"], pooled["r"], INK1)
        if pooled.get("r_pi_low") is not None and abs(pooled["r_pi_low"]) < 0.99:
            ax.hlines(-1.5, pooled["r_pi_low"], pooled["r_pi_high"],
                      color=INK3, lw=1.4, ls=(0, (3, 2)), zorder=2)
        ax.set_yticks(list(ys) + [-1.5])
        ax.set_yticklabels([f"{short(c)}\nn={int(n)}, m={int(m)}"
                            for c, n, m in zip(agg.cohort_family, agg.n, agg.m_effects)]
                           + ["POOLED"], fontsize=8.5, linespacing=1.35)
        for t, lab in zip(ax.get_yticklabels(), list(agg.cohort_family) + ["POOLED"]):
            if lab == "POOLED":
                t.set_fontweight("bold"); t.set_color(INK1)
            else:
                t.set_color(INK2)
        ax.set_ylim(-2.3, len(agg) - 0.35)
        ax.set_title(f"{LABELS[dom]}\nr = {pooled['r']:+.3f} "
                     f"[{pooled['r_ci_low']:+.2f}, {pooled['r_ci_high']:+.2f}]   "
                     f"I2 = {pooled['I2']:.0f}%", fontsize=10.5,
                     color=S2 if accent else INK1)
        ax.grid(axis="x", alpha=0.6); ax.set_axisbelow(True)
        ax.set_xticks([-1, -0.5, 0, 0.5, 1])
        for c, z, v, m in zip(agg.cohort_family, agg.z, agg.v, agg.m_effects):
            rows_out.append({"stratum": dom, "cohort_family": c, "m_effects": int(m),
                             "r": float(M.z_to_r(z)),
                             "ci_low": float(M.z_to_r(z - 1.96 * np.sqrt(v))),
                             "ci_high": float(M.z_to_r(z + 1.96 * np.sqrt(v)))})
        rows_out.append({"stratum": dom, "cohort_family": "POOLED", "m_effects": pooled["k"],
                         "r": pooled["r"], "ci_low": pooled["r_ci_low"],
                         "ci_high": pooled["r_ci_high"]})
    axes[0].set_xlim(-1.05, 1.05)
    fig.suptitle("Sleep metrics and amyloid-PET burden: cohort-family effects and "
                 "random-effects pooling", fontsize=13.5, fontweight="bold", y=1.03)
    fig.text(0.5, 0.005, "correlation r   (positive = worse sleep is associated with more amyloid)",
             ha="center", fontsize=10, color=INK2)
    fig.text(0.5, -0.075,
             "Marker area is proportional to precision. Solid bar = 95% CI. Diamond = pooled "
             "estimate; dashed line = 95% prediction interval.\nm = number of extracted effects "
             "aggregated for that cohort family. The N3 stage panel (orange) is the "
             "hypothesis's own exposure.", ha="center", fontsize=8.5, color=INK2)
    pd.DataFrame(rows_out).to_csv(os.path.join(FDAT, "fig1_forest.csv"), index=False)
    save(fig, "fig1_forest_by_stratum.png")


# =====================================================================================
def fig2_specificity(e1, e2):
    """F2 - the specificity claim: pooled estimates side by side, plus paired contrasts."""
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 4.6),
                                   gridspec_kw={"width_ratios": [1.15, 1]})
    doms = ["SWS_stage", "TST", "SE", "SWS_spectral"]
    ys = np.arange(len(doms))[::-1]
    est = [e1["strata"][x]["primary"]["r"] for x in doms]
    lo = [e1["strata"][x]["primary"]["r_ci_low"] for x in doms]
    hi = [e1["strata"][x]["primary"]["r_ci_high"] for x in doms]
    cols = [S2 if x == "SWS_stage" else S1 for x in doms]
    axL.hlines(ys, lo, hi, color=cols, lw=3, alpha=0.9)
    axL.scatter(est, ys, s=140, color=cols, edgecolor=SURFACE, lw=1.4, zorder=4)
    axL.axvline(0, color=INK3, lw=1)
    for y, e in zip(ys, est):
        axL.annotate(f"{e:+.3f}", (e, y), textcoords="offset points", xytext=(0, 13),
                     ha="center", fontsize=9.5, color=INK1, fontweight="bold")
    axL.set_ylim(-0.55, len(doms) - 0.45)
    axL.set_yticks(ys)
    axL.set_yticklabels([f"{LABELS[x]}\n({e1['strata'][x]['n_cohort_families']} cohorts, "
                         f"{e1['strata'][x]['k_effects']} effects)" for x in doms], fontsize=9)
    axL.set_xlabel("pooled correlation r with amyloid PET")
    axL.set_title("Clause (i): is each metric associated with amyloid?")
    axL.grid(axis="x", alpha=0.6); axL.set_axisbelow(True)
    axL.set_xlim(-0.75, 0.75)

    keys = ["SWS_stage_vs_TST", "SWS_stage_vs_SE", "SWS_stage_vs_SWS_spectral"]
    names = ["N3 stage − total sleep time", "N3 stage − sleep efficiency",
             "N3 stage − spectral SWA"]
    ys2 = np.arange(len(keys))[::-1]
    d_, l_, h_, p_ = [], [], [], []
    for k in keys:
        pdiff = e2["contrasts"][k]["paired"]["pooled_difference"]
        d_.append(pdiff["est"]); l_.append(pdiff["ci_low"]); h_.append(pdiff["ci_high"])
        p_.append(pdiff["p"])
    axR.axvspan(-3, 0, color=CRITICAL, alpha=0.05, zorder=0)
    axR.hlines(ys2, l_, h_, color=S2, lw=3, alpha=0.9, zorder=3)
    axR.scatter(d_, ys2, s=140, color=S2, edgecolor=SURFACE, lw=1.4, zorder=4)
    axR.axvline(0, color=INK3, lw=1, zorder=2)
    for y, e, p in zip(ys2, d_, p_):
        axR.annotate(f"Δz = {e:+.3f},  p = {p:.2f}", (e, y), textcoords="offset points",
                     xytext=(0, 13), ha="center", fontsize=9.5, color=INK1)
    axR.set_ylim(-0.55, len(keys) - 0.45)
    axR.set_yticks(ys2); axR.set_yticklabels(names, fontsize=9)
    axR.set_xlabel("within-study paired difference in Fisher z")
    axR.set_title("Clause (iii): is N3 stage STRONGER than the comparators?")
    axR.grid(axis="x", alpha=0.6); axR.set_axisbelow(True)
    axR.set_xlim(-2.9, 2.9)
    fig.suptitle("The hypothesis's two testable comparative claims", fontsize=13,
                 fontweight="bold", y=1.04)
    fig.text(0.5, -0.10,
             "Left: pooled random-effects estimates with 95% CI, cohort-family aggregated. "
             "Right: within-study paired differences — only cohorts reporting BOTH exposures in "
             "the same participants contribute.\nThe red shading marks the region where N3 stage "
             "is WEAKER than its comparator, i.e. the direction opposite to the hypothesis. "
             "All three point estimates fall in it.", ha="center", fontsize=8.5, color=INK2)
    pd.DataFrame({"contrast": names, "paired_diff_z": d_, "ci_low": l_, "ci_high": h_,
                  "p": p_}).to_csv(os.path.join(FDAT, "fig2_specificity.csv"), index=False)
    save(fig, "fig2_specificity_contrasts.png")


# =====================================================================================
def fig3_attenuation(e3):
    """F3 - the measurement-reliability null model."""
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15.5, 4.4))
    fig.subplots_adjust(wspace=0.42)
    SHORTLAB = {"SWS_stage": "N3 stage", "SE": "Sleep efficiency", "TST": "Total sleep time",
                "WASO": "WASO"}

    # panel 1: reliabilities
    iccs = e3["iccs"]
    order = ["SWS_stage", "SE", "WASO", "TST"]
    vals = [iccs[k]["icc"] for k in order]
    errs = [[iccs[k]["icc"] - iccs[k]["icc_ci"][0] for k in order],
            [iccs[k]["icc_ci"][1] - iccs[k]["icc"] for k in order]]
    cols = [S2 if k == "SWS_stage" else S1 for k in order]
    xs = np.arange(len(order))
    a1.bar(xs, vals, color=cols, width=0.62, zorder=3)
    a1.errorbar(xs, vals, yerr=errs, fmt="none", ecolor=INK2, capsize=4, lw=1.2, zorder=4)
    for x, v in zip(xs, vals):
        a1.annotate(f"{v:.3f}", (x, v + 0.045), ha="center", fontsize=9, color=INK1)
    a1.set_xticks(xs); a1.set_xticklabels([SHORTLAB.get(k, k) for k in order],
                                          fontsize=9, rotation=18, ha="right")
    a1.set_ylabel("single-night ICC(2,1)")
    a1.set_ylim(0, 1.05)
    a1.set_title("How reliably is each metric\nmeasured in one night?", fontsize=10.5)
    a1.grid(axis="y", alpha=0.6); a1.set_axisbelow(True)

    # panel 2: predicted vs observed ratio
    pred = e3["measurement_only_prediction"]
    keys = list(pred)
    xs2 = np.arange(len(keys))
    predv = [pred[k]["predicted_observed_ratio_if_true_effects_equal"] for k in keys]
    obsv = [pred[k]["observed_ratio"] for k in keys]
    w = 0.34
    a2.bar(xs2 - w / 2 - 0.01, predv, w, color=S1, label="predicted if true effects equal",
           zorder=3)
    a2.bar(xs2 + w / 2 + 0.01, obsv, w, color=S2, label="observed", zorder=3)
    a2.axhline(1.0, color=INK3, lw=1, ls=(0, (3, 2)))
    for x, v in zip(xs2 - w / 2 - 0.01, predv):
        a2.annotate(f"{v:.2f}×", (x, v + 0.05), ha="center", fontsize=9, color=INK1)
    for x, v in zip(xs2 + w / 2 + 0.01, obsv):
        a2.annotate(f"{v:.2f}×", (x, v + 0.05), ha="center", fontsize=9, color=INK1)
    a2.set_xticks(xs2)
    a2.set_xticklabels(["N3 stage\nvs TST", "N3 stage\nvs sleep efficiency"], fontsize=9)
    a2.set_ylabel("ratio of observed correlations")
    a2.set_title("Measurement error alone predicts N3\nshould look STRONGER. It doesn't.",
                 fontsize=10.5)
    a2.legend(fontsize=8, frameon=False, loc="upper center", ncol=1)
    a2.grid(axis="y", alpha=0.6); a2.set_axisbelow(True)
    a2.set_ylim(0, max(predv) * 1.35)

    # panel 3: observed vs disattenuated
    dis = e3["disattenuated_pooled"]
    order3 = ["SWS_stage", "SE", "TST"]
    ys = np.arange(len(order3))[::-1]
    for y, k in zip(ys, order3):
        o, t = dis[k]["observed_r"], dis[k]["disattenuated_r"]
        a3.annotate("", xy=(t, y), xytext=(o, y),
                    arrowprops=dict(arrowstyle="-|>", color=INK3, lw=1.6,
                                    shrinkA=4, shrinkB=4))
        a3.scatter([o], [y], s=110, color=INK3, zorder=4)
        a3.scatter([t], [y], s=150, color=S2 if k == "SWS_stage" else S1,
                   edgecolor=SURFACE, lw=1.3, zorder=5)
        a3.annotate(f"{t:+.3f}", (t, y), textcoords="offset points", xytext=(0, 12),
                    ha="center", fontsize=9, color=INK1, fontweight="bold")
    a3.axvline(0, color=INK3, lw=1)
    a3.set_yticks(ys); a3.set_yticklabels([SHORTLAB[k] for k in order3], fontsize=9)
    a3.set_ylim(-0.6, len(order3) - 0.4)
    a3.set_title("Correcting for measurement error\nwidens the gap against N3", fontsize=10.5)
    a3.set_xlim(-0.05, 0.42)
    a3.grid(axis="x", alpha=0.6); a3.set_axisbelow(True)
    a3.set_xlabel("correlation r   (grey = observed → colour = disattenuated)", fontsize=9)

    fig.suptitle("E3: could differential measurement reliability explain the pattern? "
                 "No — it predicts the opposite.", fontsize=13, fontweight="bold", y=1.04)
    pd.DataFrame([{"metric": k, "icc": iccs[k]["icc"],
                   "observed_r": dis.get(k, {}).get("observed_r"),
                   "disattenuated_r": dis.get(k, {}).get("disattenuated_r")}
                  for k in order]).to_csv(os.path.join(FDAT, "fig3_attenuation.csv"),
                                          index=False)
    save(fig, "fig3_reliability_null_model.png")


# =====================================================================================
def fig4_dose(e4, subj):
    """F4 - the dose statement in physical units, and whether the dose is even attainable."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 4.8),
                                 gridspec_kw={"width_ratios": [1, 1.05]})
    dose = e4["dose_per_10pp_N3_equivalent"]
    order = ["SWS_stage", "TST", "SE", "SWS_spectral"]
    ys = np.arange(len(order))[::-1]
    est = [dose[k]["delta_SUVR"]["point_estimate"] for k in order]
    lo = [dose[k]["delta_SUVR"]["ci_95"][0] for k in order]
    hi = [dose[k]["delta_SUVR"]["ci_95"][1] for k in order]
    cols = [S2 if k == "SWS_stage" else S1 for k in order]
    gap = e4["derived_amyloid_dispersion"]["A_neg_to_A_pos_gap_SUVR"]
    a1.axvline(gap, color=CRITICAL, lw=1.6, ls=(0, (4, 2)), zorder=2)
    a1.annotate(f"amyloid-negative → amyloid-positive separation\n"
                f"({gap:.2f} SUVR = 22 Centiloid)", (gap, 0.985),
                xycoords=("data", "axes fraction"), textcoords="offset points",
                xytext=(-6, -4), ha="right", va="top", fontsize=8.5, color=CRITICAL)
    a1.hlines(ys, lo, hi, color=cols, lw=3, alpha=0.9, zorder=3)
    a1.scatter(est, ys, s=140, color=cols, edgecolor=SURFACE, lw=1.3, zorder=4)
    a1.axvline(0, color=INK3, lw=1)
    for y, e, k in zip(ys, est, order):
        a1.annotate(f"{e:+.4f} SUVR  ({dose[k]['pct_of_A_minus_to_A_plus_SUVR_gap']:.1f}% "
                    "of the gap)", (e, y), textcoords="offset points", xytext=(0, 12),
                    ha="center", fontsize=8.5, color=INK1)
    a1.set_ylim(-0.55, len(order) - 0.35)
    a1.set_yticks(ys); a1.set_yticklabels([LABELS[k] for k in order], fontsize=9)
    a1.set_xlabel("Δ amyloid PET SUVR per 10-percentage-point worse sleep metric")
    a1.set_title("Clause (ii): the effect in physical units")
    a1.grid(axis="x", alpha=0.6); a1.set_axisbelow(True)
    a1.set_xlim(-0.11, 0.30)

    first = subj.sort_values("night_n").groupby(["cohort", "subj_num"], as_index=False).first()
    older = first[first.age >= 60]
    younger = first[first.age < 60]
    bins = np.arange(0, 46, 2.5)
    a2.hist(younger.noct_n3_pct_tst, bins=bins, color=S1, alpha=0.55, label="age < 60",
            zorder=3)
    a2.hist(older.noct_n3_pct_tst, bins=bins, color=S2, alpha=0.85,
            label="age ≥ 60 (the target population)", zorder=4)
    a2.axvline(10, color=CRITICAL, lw=1.8, ls=(0, (4, 2)), zorder=5)
    f = e4["feasibility_of_the_stated_dose"]
    a2.annotate(f"the hypothesis's stated dose:\na 10-point reduction in N3%\n\n"
                f"{100 - f['pct_with_N3_at_least_10pp']:.0f}% of adults aged 60+ have\n"
                f"less than 10 points of N3 in total",
                (0.30, 0.62), xycoords="axes fraction", va="top", fontsize=9.5,
                color=CRITICAL)
    a2.set_xlabel("N3 as % of total sleep time (Sleep-EDF Expanded, one night per subject)")
    a2.set_ylabel("number of subjects")
    a2.set_title("Is a 10-point reduction in N3 even available?")
    a2.legend(fontsize=9, frameon=False, loc="upper right")
    a2.grid(axis="y", alpha=0.6); a2.set_axisbelow(True)

    fig.suptitle("E4: translating the pooled effect into the units the hypothesis states",
                 fontsize=13, fontweight="bold", y=1.03)
    pd.DataFrame([{"stratum": k, "delta_SUVR_per_10pp": dose[k]["delta_SUVR"]["point_estimate"],
                   "ci_low": dose[k]["delta_SUVR"]["ci_95"][0],
                   "ci_high": dose[k]["delta_SUVR"]["ci_95"][1],
                   "delta_Centiloid": dose[k]["delta_Centiloid"]["point_estimate"],
                   "pct_of_A_gap": dose[k]["pct_of_A_minus_to_A_plus_SUVR_gap"]}
                  for k in order]).to_csv(os.path.join(FDAT, "fig4_dose.csv"), index=False)
    save(fig, "fig4_dose_response_and_feasibility.png")


# =====================================================================================
def fig5_robustness(e5, e6):
    """F5 - specification curve for the N3 stratum, leave-one-out, and equivalence tests."""
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(16, 4.6),
                                     gridspec_kw={"width_ratios": [1.35, 1, 1]})
    specs, vals = [], []
    specs.append(("PRIMARY", e5["primary_for_reference"]["SWS_stage"]["r"]))
    specs.append(("S1 drop imputed nulls", e5["S1_exclude_imputed_nulls"]["results"]["SWS_stage"]["r"]))
    specs.append(("S2 drop p-recovered betas", e5["S2_exclude_p_recovered_betas"]["results"]["SWS_stage"]["r"]))
    specs.append(("S1+S2 correlations only", e5["S1_and_S2_both_excluded"]["results"]["SWS_stage"]["r"]))
    for rule, res in e5["S3_df_rule_for_p_recovery"]["results"].items():
        specs.append((f"S3 df: {rule}", res["SWS_stage"]["r"]))
    for rho, res in e5["S4_dependency_rho"]["results"].items():
        specs.append((f"S4 {rho}", res["SWS_stage"]["r"]))
    specs.append(("S5 cognitively normal only", e5["S5_cognitively_normal_only"]["results"]["SWS_stage"]["r"]))
    specs.append(("S6 PSG exposures only", e5["S6_PSG_measured_exposures_only"]["results"]["SWS_stage"]["r"]))
    specs.append(("S8 fixed effect", e5["S8_fixed_effect"]["SWS_stage"]["r"]))
    for row in e5["S7_leave_one_cohort_family_out"]["SWS_stage"]:
        specs.append((f"S7 omit {row['omitted_cohort_family']}", row["r"]))
    names = [s[0] for s in specs]; vals = [s[1] for s in specs]
    ys = np.arange(len(names))[::-1]
    a1.axvspan(-0.05, 0.05, color=GOOD, alpha=0.07, zorder=0)
    a1.scatter(vals, ys, s=70, color=[S2 if n == "PRIMARY" else S1 for n in names],
               edgecolor=SURFACE, lw=1, zorder=4)
    a1.axvline(0, color=INK3, lw=1)
    a1.axvline(e5["primary_for_reference"]["SWS_stage"]["r"], color=S2, lw=1,
               ls=(0, (3, 2)), zorder=2)
    a1.set_yticks(ys); a1.set_yticklabels(names, fontsize=8)
    a1.set_xlim(-0.25, 0.55)
    a1.set_title("Specification curve: N3 stage", fontsize=11)
    a1.grid(axis="x", alpha=0.6); a1.set_axisbelow(True)
    a1.set_xlabel("pooled r, N3 stage vs amyloid PET   (green band: |r| < 0.05)")

    # leave-one-cohort-out for the spectral stratum - where it DOES matter
    rows = e5["S7_leave_one_cohort_family_out"]["SWS_spectral"]
    labs = ["(none omitted)"] + [f"omit {short(r['omitted_cohort_family'])}" for r in rows]
    rr = [e5["primary_for_reference"]["SWS_spectral"]["r"]] + [r["r"] for r in rows]
    llo = [e5["primary_for_reference"]["SWS_spectral"]["ci"][0]] + [r["ci"][0] for r in rows]
    lhi = [e5["primary_for_reference"]["SWS_spectral"]["ci"][1]] + [r["ci"][1] for r in rows]
    ys2 = np.arange(len(labs))[::-1]
    cols2 = [INK1] + [CRITICAL if "Berkeley" in l else S1 for l in labs[1:]]
    a2.hlines(ys2, llo, lhi, color=cols2, lw=2.6, alpha=0.85, zorder=3)
    a2.scatter(rr, ys2, s=110, color=cols2, edgecolor=SURFACE, lw=1.2, zorder=4)
    a2.axvline(0, color=INK3, lw=1)
    for y, v in zip(ys2, rr):
        a2.annotate(f"{v:+.3f}", (v, y), textcoords="offset points", xytext=(0, 12),
                    ha="center", fontsize=9, color=INK1)
    a2.set_ylim(-0.55, len(labs) - 0.4)
    a2.set_yticks(ys2); a2.set_yticklabels(labs, fontsize=8.5)
    a2.set_xlim(-1.02, 1.02)
    a2.set_xlabel("pooled r, spectral SWA vs amyloid PET")
    a2.set_title("Spectral SWA depends entirely\non one cohort family", fontsize=11)
    a2.grid(axis="x", alpha=0.6); a2.set_axisbelow(True)

    # equivalence tests
    tests = e6["equivalence_tests_SWS_stage"]["tests"]
    bounds = [float(k.split("_")[-1]) for k in tests]
    ps = [tests[k]["p_tost"] for k in tests]
    cols3 = [GOOD if p < 0.05 else INK3 for p in ps]
    a3.bar(np.arange(len(bounds)), ps, color=cols3, width=0.6, zorder=3)
    a3.axhline(0.05, color=CRITICAL, lw=1.4, ls=(0, (4, 2)), zorder=4)
    a3.annotate("α = 0.05", (len(bounds) - 0.4, 0.058), ha="right", fontsize=8.5,
                color=CRITICAL)
    for x, p in enumerate(ps):
        a3.annotate(f"{p:.3f}", (x, p + 0.006), ha="center", fontsize=8.5, color=INK1)
    a3.set_xticks(np.arange(len(bounds)))
    a3.set_xticklabels([f"r = {b}" for b in bounds], fontsize=9)
    a3.set_ylabel("TOST p-value")
    a3.set_xlabel("equivalence bound")
    a3.set_title("Evidence FOR a small effect\n(green = bound rejected)", fontsize=11)
    a3.grid(axis="y", alpha=0.6); a3.set_axisbelow(True)

    fig.suptitle("E5/E6: does the N3 null survive every specification, and how small is it?",
                 fontsize=13, fontweight="bold", y=1.04)
    pd.DataFrame({"specification": names, "pooled_r_SWS_stage": vals}).to_csv(
        os.path.join(FDAT, "fig5_specification_curve.csv"), index=False)
    save(fig, "fig5_robustness_and_equivalence.png")


# =====================================================================================
def fig6_funnel(d, e1):
    """F6 - funnel plots at the cohort-family level."""
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.0), sharey=True)
    for ax, dom in zip(axes, MAIN):
        agg = aggregate_by_cohort(d[d.exposure_domain == dom])
        se = np.sqrt(agg.v.values)
        pooled = e1["strata"][dom]["primary"]["est"]
        smax = max(se.max() * 1.15, 0.05)
        ys = np.linspace(0.001, smax, 50)
        ax.plot(pooled - 1.96 * ys, ys, color=GRID, lw=1.2)
        ax.plot(pooled + 1.96 * ys, ys, color=GRID, lw=1.2)
        ax.fill_betweenx(ys, pooled - 1.96 * ys, pooled + 1.96 * ys, color=GRID, alpha=0.35)
        ax.axvline(pooled, color=INK3, lw=1, ls=(0, (3, 2)))
        col = S2 if dom == "SWS_stage" else S1
        ax.scatter(agg.z, se, s=110, color=col, edgecolor=SURFACE, lw=1.2, zorder=4)
        for z, s, c in zip(agg.z, se, agg.cohort_family):
            ax.annotate(short(c), (z, s), textcoords="offset points", xytext=(0, -13),
                        ha="center", fontsize=7, color=INK2)
        ax.margins(x=0.22)
        ax.invert_yaxis()
        ax.set_xlabel("Fisher z")
        ax.set_title(LABELS[dom], fontsize=10)
        ax.grid(alpha=0.4); ax.set_axisbelow(True)
    axes[0].set_ylabel("standard error (cohort-aggregated)")
    fig.suptitle("F6: funnel plots. With 3–6 cohort families per stratum, asymmetry "
                 "tests have essentially no power — shown for completeness.",
                 fontsize=12, fontweight="bold", y=1.06)
    save(fig, "fig6_funnel.png")


def main():
    d = load()
    e1 = json.load(open(os.path.join(ROOT, "results", "e1_pooled.json")))
    e2 = json.load(open(os.path.join(ROOT, "results", "e2_specificity.json")))
    e3 = json.load(open(os.path.join(ROOT, "results", "e3_attenuation.json")))
    e4 = json.load(open(os.path.join(ROOT, "results", "e4_dose_response.json")))
    e5 = json.load(open(os.path.join(ROOT, "results", "e5_robustness.json")))
    e6 = json.load(open(os.path.join(ROOT, "results", "e6_precision.json")))
    subj = pd.read_csv(os.path.join(ROOT, "datasets", "sleep_edfx",
                                    "sleepedf_subject_level.csv"))
    print("[figures] rendering")
    fig1_forest(d, e1)
    fig2_specificity(e1, e2)
    fig3_attenuation(e3)
    fig4_dose(e4, subj)
    fig5_robustness(e5, e6)
    fig6_funnel(d, e1)
    print("[figures] done")


if __name__ == "__main__":
    main()
