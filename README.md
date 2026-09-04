# Slow-wave sleep and amyloid-beta: a study-level meta-analysis

A meta-analysis testing whether reduced N3 (slow-wave) sleep is dose-dependently associated with
amyloid-beta burden on PET in cognitively normal adults, and whether that association is stronger
than the associations with total sleep time or sleep efficiency. 40 amyloid-PET effect sizes from
17 studies across 12 independent cohort families, with cohort-family dependency handling, a
measurement-reliability null model, and conversion of the pooled effect into SUVR/Centiloid units.

**Full write-up: [REPORT.md](REPORT.md)** · Preregistered plan: [planning.md](planning.md)

## Key findings

- **N3 stage amount is not associated with amyloid PET burden**: pooled r = **+0.026** (95% CI
  −0.153 to +0.203, p = 0.68, I² = 0% across 4 cohort families). Equivalence testing rejects a true
  correlation as large as r = 0.20 (TOST p = 0.026).
- **It is weaker, not stronger, than its comparators.** Total sleep time r = +0.163, sleep
  efficiency r = +0.090. All three preregistered within-study paired contrasts point the opposite
  way to the hypothesis (Δz = −0.147, −0.049, −0.177; none significant after Holm).
- **Measurement error predicts the opposite of what is seen.** One night of PSG measures N3%
  (ICC 0.775) far better than total sleep time (ICC 0.225), so N3 should show a correlation
  **1.86×** larger than TST even if the true effects were identical. The observed ratio is
  **0.16×** — about 11× below the measurement-only prediction. Disattenuating widens the gap
  against N3 (true r ≈ 0.03 vs 0.23).
- **The dose is negligible and the stated dose is unattainable.** A 10-percentage-point reduction
  in N3 maps to ΔSUVR = **+0.005** (95% CI −0.028 to +0.038) — 1.7% of the amyloid-negative /
  amyloid-positive separation. Meanwhile 78% of adults aged 60+ have less than 10 percentage
  points of N3 in total, so the stated contrast spans the whole observable range.
- **The positive slow-wave literature is one cohort.** Spectral slow-wave activity pools to
  r = +0.112 (I² = 77%); removing the Berkeley Aging Cohort Study drops it to **r = −0.025**.
  That stratum is too imprecise for any conclusion in either direction — its honest verdict is
  *insufficient evidence*, not *no effect*.

**Bottom line:** the hypothesis as literally worded is not supported. The stage/spectral
distinction, routinely collapsed in the secondary literature, is where the whole question turns.

## Reproducing

```bash
uv venv && source .venv/bin/activate && uv sync   # Python >= 3.10
python src/run_all.py                             # full pipeline + validation, ~29 s on CPU
```

No GPU required — this is a study-level meta-analysis. Seed 42; the only stochastic step is the
E4 Monte Carlo. `src/validate.py` re-runs the entire pipeline in a scratch copy and asserts that
all eight result files are **bit-identical**; it also cross-checks the hand-implemented REML
estimator against PyMARE (agreement to 4.4e-09). Current status: **39/39 checks pass**.

## Structure

```
planning.md              preregistered plan: motivation, direction budget, analysis spec
REPORT.md                the research report - read this
CODE_WALKTHROUGH.md      how the code works and how to extend it
literature_review.md     54-paper review (resource-gathering phase)
resources.md             catalogue of papers, datasets and code

src/
  metalib.py             meta-analysis primitives: Fisher z, REML tau^2, HKSJ, RVE,
                         correlated-effects aggregation, Egger, trim-and-fill, permutation, TOST
  prepare_data.py        effect harmonisation onto the Fisher-z scale + QC
  exp1_pooling.py        E1  stratified pooling with cohort-family dependency handling
  exp2_specificity.py    E2  between-stratum and within-study paired contrasts
  exp3_attenuation.py    E3  measurement-reliability null model
  exp4_dose_response.py  E4  Monte-Carlo conversion to SUVR/Centiloid + feasibility
  exp5_robustness.py     E5/E6  11 sensitivity analyses, bias diagnostics, power, equivalence
  make_figures.py        six figures + tabular companions
  validate.py            39 analytic, conversion and reproducibility checks
  run_all.py             runs everything in order

datasets/                effect sizes, reliability, N3 norms, PET reference scale, Sleep-EDF
results/                 all JSON/CSV outputs (see the appendix of REPORT.md)
figures/                 fig1-fig6 PNG at 300 dpi
papers/                  54 full texts + metadata
```

## Data sources

Effect sizes were hand-extracted from the 54 full-text papers in `papers/`, each with a verbatim
provenance quote and a Europe-PMC-verified DOI/PMID
(`datasets/meta_analysis/effect_sizes.csv`). Sleep-stage reliability and N3 norms were derived
from Sleep-EDF Expanded hypnograms (PhysioNet). SUVR/Centiloid anchors come from Jack et al. 2017.
No individual-participant data was available — ADNI, A4, BLSA, AIBL, WRAP, Mayo MCSA, Age-Well and
OASIS-3 all require data-use agreements with human review, which caps this at a study-level
analysis (see Limitations in REPORT.md).
