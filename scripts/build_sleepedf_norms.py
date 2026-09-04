#!/usr/bin/env python3
"""Merge Sleep-EDF hypnogram summaries with subject age/sex; derive
   (a) age-stratified normative distributions of N3%, TST, sleep efficiency
   (b) night-to-night (test-retest) reliability of each metric (ICC(2,1), Pearson r)
Nocturnal-window metrics are used throughout (comparable to in-lab PSG).
"""
import pandas as pd, numpy as np, os

os.makedirs("datasets/derived", exist_ok=True)
st = pd.read_csv("datasets/sleep_edfx/hypnogram_stages.csv")
d = st[st.cohort == "cassette"].copy()
meta = pd.read_excel("datasets/sleep_edfx/SC-subjects.xls").rename(
    columns={"subject": "subj_num", "night": "night_n", "sex (F=1)": "sex_f1"})
d = d.merge(meta[["subj_num", "night_n", "age", "sex_f1"]], on=["subj_num", "night_n"], how="left")
assert d.age.notna().all()
d["sex"] = d.sex_f1.map({1: "F", 2: "M"})
d.to_csv("datasets/sleep_edfx/sleepedf_subject_level.csv", index=False)
print(f"subject-level: {d.shape[0]} recordings, {d.subj_num.nunique()} subjects, age {d.age.min()}-{d.age.max()}")

METRICS = ["noct_n3_pct_tst", "noct_n3_min", "noct_n2_pct_tst", "noct_rem_pct_tst",
           "noct_tst_min", "noct_sleep_eff_pct", "noct_waso_min"]

rows = []
for lo, hi, lbl in [(20,40,"20-39"),(40,60,"40-59"),(60,75,"60-74"),(75,102,"75+"),(0,999,"all")]:
    s = d[(d.age >= lo) & (d.age < hi)]
    if len(s) < 3: continue
    for m in METRICS:
        rows.append({"age_band": lbl, "n_recordings": len(s), "n_subjects": s.subj_num.nunique(),
                     "metric": m, "mean": round(s[m].mean(), 3), "sd": round(s[m].std(ddof=1), 3),
                     "median": round(s[m].median(), 3), "p25": round(s[m].quantile(.25), 3),
                     "p75": round(s[m].quantile(.75), 3), "min": round(s[m].min(), 3),
                     "max": round(s[m].max(), 3)})
norms = pd.DataFrame(rows)
norms.to_csv("datasets/derived/n3_norms_by_age.csv", index=False)
print("\n--- normative N3% of TST by age band (Sleep-EDF Cassette) ---")
print(norms[norms.metric == "noct_n3_pct_tst"].to_string(index=False))

# night-to-night reliability
piv = d.pivot_table(index="subj_num", columns="night_n", values=METRICS)
rel = []
for m in METRICS:
    x, y = piv[(m, 1)], piv[(m, 2)]
    ok = x.notna() & y.notna(); x, y = x[ok].values, y[ok].values
    n, k = len(x), 2
    Y = np.column_stack([x, y]); gm = Y.mean()
    MSR = k * ((Y.mean(1) - gm) ** 2).sum() / (n - 1)
    MSC = n * ((Y.mean(0) - gm) ** 2).sum() / (k - 1)
    MSE = ((Y - Y.mean(1)[:, None] - Y.mean(0)[None, :] + gm) ** 2).sum() / ((n - 1) * (k - 1))
    icc = (MSR - MSE) / (MSR + (k - 1) * MSE + k * (MSC - MSE) / n)
    rel.append({"metric": m, "n_subjects": n,
                "pearson_r_night1_night2": round(float(np.corrcoef(x, y)[0, 1]), 4),
                "icc_2_1": round(float(icc), 4),
                "within_subject_sd": round(float(np.sqrt(MSE)), 4),
                "pooled_sd": round(float(Y.std(ddof=1)), 4),
                "attenuation_factor_sqrt_icc": round(float(np.sqrt(max(icc, 0))), 4)})
rel = pd.DataFrame(rel)
rel.to_csv("datasets/derived/night_to_night_reliability.csv", index=False)
print("\n--- single-night test-retest reliability (n=75 subjects with 2 nights) ---")
print(rel.to_string(index=False))
