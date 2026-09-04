# Code Walkthrough

How the analysis is put together, what each piece does and why, and how to change it.

## Data flow

```
datasets/meta_analysis/effect_sizes.csv        (49 hand-extracted effects, 17 studies)
            │
            ▼  prepare_data.py  — apply the four conversion rules; QC; hard assertions
results/harmonised_effects.csv                 (49 rows + z, v, conversion_rule, usable, reason)
results/qc_report.json
            │
            ├─▶ exp1_pooling.py       ──▶ results/e1_pooled.json       (per-stratum pooled r)
            │        │
            │        ├─▶ exp2_specificity.py ──▶ results/e2_specificity.json
            │        ├─▶ exp3_attenuation.py ──▶ results/e3_attenuation.json
            │        │      ▲ datasets/derived/night_to_night_reliability.csv
            │        ├─▶ exp4_dose_response.py ──▶ results/e4_dose_response.json
            │        │      ▲ datasets/sleep_edfx/sleepedf_subject_level.csv
            │        │      ▲ amyloid_pet_reference_scale.csv (constants, transcribed in-file)
            │        └─▶ exp5_robustness.py ──▶ results/e5_robustness.json, e6_precision.json
            │
            └─▶ make_figures.py ──▶ figures/*.png + results/figure_data/*.csv
                     validate.py ──▶ results/validation_report.json
```

`exp1_pooling.py` exports `load()` and `aggregate_by_cohort()`, which every downstream experiment
imports — so the analysis set and the dependency-handling rule are defined in exactly one place.

Run everything with `python src/run_all.py`.

## `src/metalib.py` — the statistical primitives

Implemented from primary formulations rather than taken from a package, because no available
Python meta-analysis package provides HKSJ standard errors, cluster-robust variance and exact
permutation inference together. `validate.py` cross-checks the REML point estimates against
PyMARE on the real data (agreement to 4.4e-09).

### Scale conversions

| Function | Purpose |
|---|---|
| `r_to_z(r)` / `z_to_r(z)` | Fisher transform and inverse; clips at ±0.9999 to keep `artanh` finite |
| `var_z(n)` | `1/(n−3)`, the sampling variance of Fisher z |
| `tau_to_r(tau)` | `sin(πτ/2)` — Greiner's relation, exact under bivariate normality |
| `r_from_p_and_df(p, df, sign)` | test-statistic recovery, see below |

#### `r_from_p_and_df(p_two_sided, df, sign)`

**Purpose.** Recover a signed partial correlation from a reported two-sided p-value and residual
degrees of freedom, so that unstandardised regression coefficients can be pooled with correlations.

**Rationale.** A coefficient's p-value is a monotone function of its partial correlation given df:
`t = r√df / √(1−r²)`, so `r = t/√(t²+df)` with `t = F⁻¹(1−p/2; df)`. This is the standard
conversion (Rosenthal & Rubin 1979). It is the only assumption-light way to keep 12 of the 40
amyloid-PET effects — including A4 (n = 4,417) and AIBL (n = 189) — in a common metric. Raw betas
are never pooled in their native units.

**Verification.** `validate.py` round-trips it: for a simple correlation with df = n−2, recovering
r from the exact p-value of that r returns the original r to 1e-8, for four different (r, n) pairs.

```python
r = M.r_from_p_and_df(p_two_sided=0.042, df=64 - 1 - 3, sign=-1)   # Carvalho 2024, E020
```

### Random-effects model

#### `re_meta(y, v, method="REML", hksj=True, level=0.95)`

**Returns** a dict with `est`, `se`, `ci_low/high`, `p`, `tau2`, `Q`, `Q_p`, `I2`, `H2`,
`pi_low/high` (prediction interval) and the same quantities back-transformed to the r scale.

**Design decisions.**

- **REML τ² by Fisher scoring** (Viechtbauer 2005), started from the DerSimonian–Laird moment
  estimate, truncated at 0. REML rather than DL because DL is biased downward when k is small.
- **HKSJ on by default.** With k = 3–6, the standard Wald interval is far too narrow. HKSJ
  rescales the SE by the weighted residual mean square and uses a `t_{k−1}` reference. The
  implementation applies the standard safeguard of never letting HKSJ *shrink* the SE below the
  model-based one, so it is conservative in both directions.
- **Heterogeneity at fixed-effect weights**, which is the conventional definition of Q and I².
- **Prediction interval** (Higgins 2009) only when k ≥ 3, since it needs `t_{k−2}`.

#### `aggregate_correlated(y, v, rho)`

Collapses m dependent effects from one cohort into a single synthetic effect using the Borenstein
correlated-effects variance `v̄ = (1/m²)(Σvᵢ + Σ_{i≠j} ρ√(vᵢvⱼ))`. ρ = 1 reproduces the
single-effect variance (no precision gain from the extra effects); ρ = 0 reproduces independence.
Validated to satisfy both limits.

#### `rve_meta(y, v, cluster, rho)`

Hedges–Tipton–Johnson sandwich variance for an intercept-only meta-regression, with the `m/(m−1)`
small-sample correction and `df = m−1`. Emits an explicit warning below 4 clusters, where RVE is
known to be unreliable — every stratum here except TST triggers that warning, which is why RVE is
reported alongside rather than as the primary estimator.

### Contrasts

- `paired_contrast(z_a, v_a, z_b, v_b, rho_within)` — `Var(z_a − z_b) = v_a + v_b − 2ρ√(v_a v_b)`.
  ρ_w is the correlation between two estimates measured on the same participants; it is not
  identified from published summaries, so it is assumed (0.5) and varied (0, 0.3, 0.7). Note that
  larger ρ_w *shrinks* the variance, so ρ_w = 0 is the conservative choice for detecting a
  difference.
- `between_stratum_contrast(res_a, res_b)` — Wald test on two pooled estimates. Valid only for
  disjoint samples; here the strata share cohorts, so this is anti-conservative and is reported
  with that caveat while the paired version is primary.

### Bias and small-k inference

- `egger_test(y, v)` — regression of the standard normal deviate on precision; the intercept is
  the asymmetry test. Flags `interpretable=False` below k = 5.
- `trim_and_fill(y, v, side="auto")` — Duval & Tweedie L0. **Note:** the first implementation had
  a side-selection bug that made it silently return 0 imputed studies on a strongly asymmetric
  funnel. It was rewritten to run the L0 iteration in both orientations (`_l0_missing_left` on
  `y` and on `−y`) and take whichever imputes more. `validate.py` now asserts it detects a
  constructed asymmetric funnel (5 imputed, r 0.139 → 0.103), imputes nothing on a symmetric one,
  and is sign-symmetric. Flags `interpretable=False` below k = 10.
- `permutation_p(y, v)` — exact sign-flip test, enumerating all 2^k assignments when 2^k ≤ 65,536
  (always, here) and sampling otherwise. Far more trustworthy than a Wald p at k ≤ 6.
- `tost_equivalence(est, se, df, bound_z)` — two one-sided tests. Rejecting both concludes the true
  effect is smaller than the bound. This is what turns "no evidence of an effect" into "evidence
  of no meaningful effect".
- `detectable_effect(v, power, alpha, tau2)` / `power_for_effect(...)` — minimum detectable effect
  and power, so a null result can state what it could have found.

## `src/prepare_data.py` — harmonisation

Applies the four preregistered conversion rules, one branch per `effect_type`, and records for
every row which rule fired (`conversion_rule`), whether it is usable, and if not, why
(`usable_reason`). Nothing is dropped silently.

`COVARIATE_COUNTS` maps each `adjustment` string to an explicit covariate count so that
`df = n − 1 − n_covariates` is auditable rather than parsed heuristically; strings ending in
`etc` are marked `df_approximate` and sensitivity analysis S3 varies the rule.

Finishes with hard assertions — duplicate IDs, |r| > 1, sign inversion, missing n — so a corrupted
input fails loudly instead of producing a plausible-looking wrong answer.

## The experiment scripts

| Script | Answers | Key implementation note |
|---|---|---|
| `exp1_pooling.py` | clause (i); supplies the strata for everything else | reports `naive` / `primary` (cohort-aggregated) / `rve` side by side, plus design subgroups |
| `exp2_specificity.py` | clause (iii) | `paired_test()` intersects the cohort-family index of two strata, so only cohorts measuring both exposures contribute; Holm across the three confirmatory contrasts |
| `exp3_attenuation.py` | the competing explanation | `disattenuate()` works on the r scale by delta method and propagates ICC uncertainty (`Var(r)/icc + r²·Var(icc)/(4icc³)`), then transforms back to z; runs an all-modality and a PSG-only arm because 6 of 9 TST effects are self-report |
| `exp4_dose_response.py` | clause (ii) | 10,000-draw Monte Carlo over pooled r (t_{k−1} on the z scale), SD of N3% (χ² sampling of an SD) and SD of amyloid (Beta prior on A+ prevalence applied to a two-component mixture); includes the feasibility check |
| `exp5_robustness.py` | E5 + E6 | 11 sensitivity analyses; `recompute_recovered()` re-derives r under alternative df rules without touching the rest of the pipeline |

## Figures

`make_figures.py` follows the project's data-viz guidance: form chosen before color, one hue where
there is one series, blue/orange only where two entities must be distinguished (that pair clears
the all-pairs CVD floors), recessive grid, direct labels rather than a number on every mark, and
text in ink tokens rather than series color. Every figure writes a CSV companion to
`results/figure_data/` so identity is never carried by color alone.

## Reproducing and verifying

```bash
source .venv/bin/activate
python src/run_all.py          # prepare → E1..E6 → figures → validate   (~29 s, CPU only)
python src/validate.py         # validation alone; exits non-zero on any failure
```

`validate.py` runs three families of checks:

- **V1 analytic** — closed-form identities (Fisher round-trip, FE pooling = inverse-variance mean,
  aggregation limits at ρ = 0 and 1), behavioural properties (TOST rejects a wide bound for a
  precise zero but not a narrow bound when imprecise; trim-and-fill detects asymmetry and not
  symmetry), and a **PyMARE cross-check** of REML on the real data.
- **V2 conversion** — every usable row has finite z and positive v; `v = 1/(n−3)` exactly;
  direct-r rows unchanged; signs preserved; the two Winer 2020 headline numbers survive the
  pipeline unchanged; every unusable row carries a reason; the analysis set is the size the report
  claims.
- **V3 reproducibility** — re-runs the whole pipeline and SHA-256-compares all eight result files.

Current status: **39/39 pass**.

## Extending this

- **Add a study.** Append rows to `datasets/meta_analysis/effect_sizes.csv` with a
  `burden_aligned_effect`, a `cohort_family` and a provenance quote, then re-run `run_all.py`. If
  the study reports an effect type not yet handled, add a branch in `prepare_data.py::harmonise`
  and a covariate count in `COVARIATE_COUNTS`; the assertions will tell you if you missed one.
- **Change the dependency assumption.** `RHO_PRIMARY` in `exp1_pooling.py`; S4 already reports
  ρ ∈ {0.3, 0.6, 0.8, 1.0}.
- **Add a stratum.** Add the key to `STRATA` in `exp1_pooling.py` and to `MAIN`/`LABELS` in the
  downstream scripts.
- **Move to individual-participant data.** `metalib.re_meta` would be replaced by a mixed model;
  the harmonisation layer and the E3 reliability model transfer unchanged and would become
  considerably more informative.

## Known limitations of the code

- `_tau2_reml` uses fixed-point Fisher scoring with a 200-iteration cap; it converges in a handful
  of iterations on every stratum here but has no formal convergence guarantee for pathological
  inputs.
- `rve_meta` uses the approximate (rather than the correlated-effects-with-known-ρ) weighting from
  Hedges–Tipton–Johnson, which is standard practice but is one of several defensible choices.
- The amyloid PET reference-scale constants in `exp4_dose_response.py` are transcribed from
  `amyloid_pet_reference_scale.csv` as module-level constants rather than read from it, so that
  the Monte-Carlo section is self-documenting; if that CSV changes, update the constants.
