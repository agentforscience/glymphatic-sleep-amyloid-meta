# Datasets

Four datasets support the analysis. Two are **curated by hand from the literature** and are the
primary analysis inputs; two are **derived from an openly downloadable PSG corpus** and supply
the calibration constants the hypothesis needs.

Bulk data is excluded from git (see `datasets/.gitignore`); the curated CSVs are small and are
tracked, because they are not reproducible by re-running a download.

---

## 1. `meta_analysis/effect_sizes.csv` — PRIMARY

Study-level effect sizes linking sleep measures to amyloid burden, extracted by hand from the
full texts in `papers/`.

- **49 effects** from **17 studies** across **12 independent cohort families**
- **Source**: manual extraction; every row carries a `provenance` field quoting the sentence or
  table cell it came from, plus a DOI and PMID verified against Europe PMC
- **License**: the numbers are facts extracted from published papers; the compilation is ours

### Columns that matter

| Column | Meaning |
|---|---|
| `exposure_domain` | `SWS_stage`, `SWS_spectral`, `TST`, `SE`, `WASO`, `other_stage`, `*_variability` — **the key stratifier** |
| `exposure_measure` | the exact metric as the paper defined it (units and band limits differ across studies) |
| `outcome_class` | `deposition_pet`, `deposition_pet_change`, `soluble_csf` |
| `effect_type` | `pearson_r`, `partial_r`, `kendall_tau`, `r_from_partial_eta2`, `unstd_beta_*`, `narrative_null` |
| `effect` | the effect as reported |
| `burden_aligned_effect` | sign-normalised so **positive = worse sleep goes with more amyloid** |
| `burden_aligned_note` | why the sign was flipped, or a warning about the row |
| `cohort_family` | **use this for dependency handling** — Berkeley, A4, AIBL, Age-Well etc. |
| `provenance` | verbatim quote supporting the number |

### Three traps this file is designed to make visible

1. **Cohort non-independence.** `Mander2015`, `Winer2019` and `Winer2020` are all
   `cohort_family = BerkeleyAgingCohort` with overlapping participants. Pooling them as three
   independent studies overstates precision badly.
2. **CSF sign instability.** For `outcome_class = soluble_csf`, low CSF Aβ42 indicates *more*
   deposited plaque in cohorts that include AD, but *less* soluble production in pre-deposition
   cohorts. Varga 2016 and Liguori 2020 therefore point opposite ways. Do not pool CSF with PET
   without stratifying on cognitive status.
3. **`narrative_null` rows are real data.** They encode nulls the papers reported in prose
   without a coefficient (e.g. Champetier 2024). Dropping them biases the pooled estimate
   upward.

### Load it
```python
import pandas as pd
d = pd.read_csv("datasets/meta_analysis/effect_sizes.csv")
pet = d[d.outcome_class.isin(["deposition_pet", "deposition_pet_change"])]
```

---

## 2. `meta_analysis/amyloid_pet_reference_scale.csv` and `meta_analysis/prior_meta_analyses.csv`

- **Reference scale** (12 rows): SUVR ↔ Centiloid anchors needed to express a standardized
  effect in the units the hypothesis states. Amyloid-positivity threshold 1.42 PiB SUVR =
  19 Centiloid; A− median 1.31 SUVR (Centiloid 9); A+ median 1.57 SUVR (Centiloid 31).
  Sources: Jack 2017 (*Lancet Neurol*), Hanseeuw 2021.
- **Prior meta-analyses** (7 rows): pooled benchmarks from Chen 2025 to contextualise our own
  estimates, including their I² = 99% for the sleep-quality/amyloid-PET pooling.

---

## 3. Sleep-EDF Expanded — `sleep_edfx/`

Openly downloadable PSG corpus from PhysioNet, used here **only for its hypnograms**. The full
corpus is ~8 GB of raw EDF; we download only the annotation files (~800 KB total), which is all
that is needed to compute sleep-stage percentages.

- **Source**: https://physionet.org/content/sleep-edfx/1.0.0/
- **Contents downloaded**: 153 Sleep-Cassette + 44 Sleep-Telemetry hypnograms, plus
  `SC-subjects.xls` / `ST-subjects.xls` (age, sex, lights-off time)
- **Sleep-Cassette**: 78 healthy subjects aged 25–101, ~2 nights each — this two-night design is
  what makes the reliability analysis possible
- **License**: Open Data Commons Attribution License v1.0
- **No amyloid data.** This corpus cannot test the hypothesis directly; it supplies the
  measurement-property constants the hypothesis's dose statement requires.

### Download
```bash
source .venv/bin/activate
python scripts/get_sleepedf.py         # hypnograms + subject metadata (~800 KB)
python scripts/parse_hypnograms.py     # -> sleep_edfx/hypnogram_stages.csv
python scripts/build_sleepedf_norms.py # -> derived/*.csv
```

### A caveat about the recording window
Sleep-Cassette records span roughly 20 hours of ambulatory recording, so metrics computed over
the whole file include daytime wake and naps and are not comparable to in-lab PSG. Every metric
is therefore computed twice: `full_*` over the entire recording and `noct_*` over a 12-hour
window anchored at the subject's lights-off time. **Use the `noct_*` columns.** The difference
is material — full-window sleep efficiency averages 78.5% (SD 18.2) versus 87.9% (SD 9.7)
for the nocturnal window.

---

## 4. Derived calibration tables — `derived/`

### `night_to_night_reliability.csv`
Single-night test-retest reliability from the 75 Sleep-Cassette subjects with two nights.

| metric | ICC(2,1) | attenuation factor √ICC |
|---|---|---|
| N3 % of TST | 0.775 | 0.880 |
| N3 minutes | 0.790 | 0.889 |
| sleep efficiency % | 0.562 | 0.750 |
| WASO min | 0.517 | 0.719 |
| N2 % of TST | 0.390 | 0.625 |
| REM % of TST | 0.225 | 0.474 |
| TST min | 0.225 | 0.474 |

This is the input to the attenuation null model: if the *true* correlations with amyloid were
identical across metrics, a single-night N3% measurement would still show an observed
correlation about **1.9x** that of a single-night TST measurement, purely because it is measured
more reliably. Any claim that N3 is *biologically* special has to clear this bar first.

### `n3_norms_by_age.csv`
Age-stratified distributions (mean, SD, median, IQR, range) for N3%, N3 min, N2%, REM%, TST,
sleep efficiency and WASO.

| age band | n recordings | N3 % of TST, mean (SD) |
|---|---|---|
| 20–39 | 39 | 17.30 (9.62) |
| 40–59 | 41 | 8.25 (7.24) |
| 60–74 | 41 | 9.29 (10.06) |
| 75+ | 32 | 5.25 (8.32) |
| all | 153 | 10.21 (9.85) |

These supply the denominator for the dose-response conversion. In adults over 60 the SD of N3%
is roughly 8–10 percentage points, so **a 10-percentage-point reduction in N3% is about 1.0–1.2
standard deviations** — a large shift, not a marginal one. That matters for interpreting what
"each 10% reduction" would have to mean in SUVR units.

### Sample records
`sleep_edfx/hypnogram_stages.csv`, first two Sleep-Cassette recordings (nocturnal window):

| record | subj | night | age | noct_tst_min | noct_sleep_eff_pct | noct_n3_pct_tst | noct_rem_pct_tst |
|---|---|---|---|---|---|---|---|
| SC4001EC | 0 | 1 | 33 | 326.5 | 90.57 | 33.69 | 19.14 |
| SC4002EC | 0 | 2 | 33 | 472.0 | 93.65 | 31.46 | 22.78 |

---

## Datasets we could not obtain

Every cohort in the primary literature — ADNI, A4, BLSA, AIBL, WRAP, Mayo Clinic Study of Aging,
Age-Well/IMAP, OASIS-3 — requires a data-use agreement with human review that cannot be
completed inside this pipeline. NSRR (`sleepdata.org`) was also unreachable from this
environment. **No open dataset pairs polysomnography with amyloid PET in the same
participants.** That is precisely why the study-level meta-analysis in `meta_analysis/` is the
primary analysis rather than an individual-participant-data analysis, and it is a real ceiling
on what can be concluded here.
