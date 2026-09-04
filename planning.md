# Research Plan and Direction Budget

> **Document status.** Sections 0 and A below were added by the `experiment_runner` phase
> (2026-09-04) and constitute the *preregistered* analysis specification. Everything from
> "What the literature already says" onward is the `resource_finder` phase's plan, retained
> verbatim; the direction budget (D1/D2/D3 kept, D4-D8 pruned) was **not** re-opened, because
> no new evidence overturned its ranking.

---

# PHASE 0 - Motivation & Novelty Assessment

## Why This Research Matters

Sleep is one of a very small number of *modifiable* candidate risk factors for Alzheimer's
disease, and the glymphatic hypothesis gives it a concrete mechanism: slow waves drive
CSF pulsation, CSF pulsation clears interstitial amyloid-β, so less slow-wave sleep should mean
more plaque. That story is now widely repeated in reviews, in press coverage, and in the design
of behavioural and acoustic-stimulation trials that aim to boost N3 sleep in order to slow
amyloid accumulation. If the human evidence linking *slow-wave sleep specifically* to amyloid
burden is weaker or more fragile than the narrative implies, then trials are being powered and
targeted on an effect that may not exist at the claimed magnitude - an expensive and
consequential error. The people who benefit from getting this right are trialists choosing an
intervention target, clinicians deciding what to tell cognitively normal older adults about
their sleep, and the field's own priors about which sleep metric to measure.

## Gap in Existing Work

From `literature_review.md` (54 full texts) four gaps are established, and every one of them is
load-bearing for the hypothesis:

1. **No meta-analysis has ever pooled slow-wave-specific exposures against amyloid PET.**
   Chen 2025 (the most recent and largest, 30 studies / 14,997 subjects) pooled only *sleep
   quality* and *sleep duration*; Cui 2022 and Bubu 2019/2020 covered obstructive sleep apnea.
   The exposure the hypothesis names has never been meta-analysed against the outcome it names.
2. **Stage and spectral "slow-wave sleep" are silently conflated.** Winer 2020 reports both in
   the *same* 32 participants: spectral <1 Hz SWA r = -0.52 (p = 0.002) versus N3 *stage*
   duration r = -0.07 (p = 0.72). A pooled estimate that merges these two exposures is
   uninterpretable, yet the distinction is routinely lost in secondary literature.
3. **Cohort non-independence is unhandled.** Mander 2015, Winer 2019 and Winer 2020 are all the
   Berkeley Aging Cohort Study with overlapping participants - 13 of our 49 extracted effects.
   Treating them as independent studies would inflate precision by roughly √3.
4. **Differential measurement reliability has never been tested as a competing explanation.**
   No paper in the corpus asks whether "metric A correlates more strongly than metric B" is
   simply "metric A is measured more reliably than metric B on a single night."

## Our Novel Contribution

Four things this project does that the literature has not:

- **C1. The first meta-analytic estimate of the N3-stage ↔ amyloid-PET association**, kept
  strictly separate from the spectral-SWA ↔ amyloid-PET association, with variance handled at
  the level of the *cohort family* rather than the publication.
- **C2. A direct, within-study paired test of the specificity claim.** Because several studies
  report N3 stage, TST and sleep efficiency against amyloid in the *same* participants, the
  contrast "is N3 stronger than TST?" can be tested as a paired difference within study rather
  than as a between-stratum comparison across different samples. This removes between-cohort
  confounding from the exact claim the hypothesis makes, and to our knowledge has not been done.
- **C3. A quantitative measurement-reliability null model.** Using single-night test-retest
  ICCs derived from Sleep-EDF Expanded (75 two-night subjects), we compute what correlation
  ratio between exposures would be expected *if the underlying true correlations were identical*
  and only measurement precision differed. The hypothesis's specificity claim only survives if
  the observed N3-vs-TST gap beats this null.
- **C4. The first explicit estimate of the hypothesis's own dose statement** - expected amyloid
  PET SUVR / Centiloid change per 10 percentage-point reduction in N3% - obtained by converting
  pooled standardised effects using external N3% dispersion and published SUVR↔Centiloid
  anchors, with the approximation's assumptions stated rather than hidden.

We also correct a defect in the benchmark meta-analysis: Chen 2025 reports "Fisher z = -0.002,
95% CI -0.003 to -0.001" for sleep duration vs amyloid PET, which are raw regression-coefficient
units mislabelled as Fisher z. Our pooling uses a genuinely common metric.

## Experiment Justification

| # | Experiment | Why it is necessary |
|---|---|---|
| **E1** | Stratified random-effects pooling of N3-stage, spectral-SWA, TST and SE against amyloid PET, with cohort-family aggregation | This *is* the hypothesis's first clause. Without it there is no meta-analytic estimate at all. Stratification is what prevents the stage/spectral conflation (gap 2); cohort aggregation is what prevents Berkeley counting three times (gap 3). |
| **E2** | Specificity contrasts: N3-stage vs TST, N3-stage vs SE, and (scientific interest) N3-stage vs spectral-SWA - run both between-stratum with cluster-robust variance and within-study paired | This is the hypothesis's third clause ("stronger than TST or SE alone"), which is a *comparative* claim and therefore cannot be answered by any single pooled estimate. The paired version (C2) is the only version free of between-cohort confounding. |
| **E3** | Reliability-attenuation null model | Without it, any observed difference between exposures is uninterpretable: it could be biology or it could be that N3% has ICC 0.775 while TST has ICC 0.225. E3 turns the difference into a testable excess over a measurement-only prediction (C3, gap 4). |
| **E4** | Dose-response conversion to ΔSUVR and ΔCentiloid per 10 pp N3, with uncertainty propagation | This is the hypothesis's second clause, stated in physical units. It has never been estimated (C4). It also converts a possibly-significant standardized effect into a magnitude that can be compared against the A-/A+ separation, which is what determines whether the effect is *clinically* meaningful rather than merely non-zero. |
| **E5** | Robustness: leave-one-cohort-out, ρ sensitivity for aggregation, inclusion/exclusion of imputed nulls and of p-value-recovered betas, Egger regression, trim-and-fill, permutation tests, and equivalence testing (TOST) of the N3 stratum | With k this small, a single influential cohort or a single analytic choice can drive the whole answer; robustness checks are what separate a defensible negative result from an underpowered one. TOST is specifically required because we anticipate a null: a non-significant pooled estimate is only evidence *for* absence if the CI excludes effects large enough to matter. |
| **E6** | Precision / detectable-effect analysis per stratum | An honest negative result must state what it could and could not have detected. This bounds the interpretation of E1-E2 and is the main guard against overclaiming "no effect" when the truthful claim is "no effect larger than r ≈ x". |

---

# SECTION A - Preregistered analysis specification (experiment_runner)

Fixed **before** running any model beyond the resource_finder smoke test.

## A.1 Primary analysis set

`datasets/meta_analysis/effect_sizes.csv` filtered to
`outcome_class ∈ {deposition_pet, deposition_pet_change}`. CSF outcomes (`soluble_csf`) are
**excluded from the primary model** and analysed only as a cognitive-status-stratified
secondary, because the sign of CSF Aβ42 relative to deposited burden reverses between
pre-deposition and post-deposition samples (Varga 2016 vs Liguori 2020).

All effects use `burden_aligned_effect`, signed so that **positive = worse sleep is associated
with more amyloid**, i.e. positive supports the hypothesis.

## A.2 Effect harmonisation rules (fixed in advance, applied uniformly)

Everything is pooled on the Fisher-z scale, `z = artanh(r)`, `v = 1/(n-3)`.

| `effect_type` | Rule | Justification |
|---|---|---|
| `pearson_r`, `partial_r`, `r_from_partial_eta2` | used directly as r | already the target metric; ηp² → r = √ηp² is exact for a 1-df contrast |
| `kendall_tau` | r = sin(πτ/2) | Greiner's relation, exact under bivariate normality |
| `unstd_beta_*` (per-SD, per-hour, per-tertile, per-ordinal-category, group-contrast) | **t-statistic recovery**: from the reported two-sided p and residual df, t = sign × F⁻¹(1 - p/2; df), then r_partial = t/√(t² + df) | A regression coefficient's p-value is a monotone function of its partial correlation given df. Recovering r from (p, df) is a standard and assumption-light conversion (Rosenthal & Rubin) and is the only way to keep 10 of 40 PET effects - including the two largest cohorts, A4 (n = 4,417) and AIBL (n = 189) - in a common metric. Raw betas are **never** pooled in their native units. |
| `narrative_null` | imputed z = 0, v = 1/(n-3) | These encode genuinely reported nulls ("no association with TST or SE"). Dropping them biases every pooled estimate away from zero. Imputing z = 0 is conservative in magnitude but *anti-conservative in precision*, so it is carried as the primary rule **and** removed in a preregistered sensitivity analysis. |

Residual df for recovery: `n - 1 - n_covariates`, with `n_covariates` taken from the
`adjustment` field; where the covariate count is not recoverable, df = n - 5 is used and the
row is flagged. Sensitivity to this choice is reported.

## A.3 Dependency handling

Primary: **aggregate to one effect per `cohort_family` per stratum** using the Borenstein
correlated-effects formula - mean of the z's, with
`v̄ = (1/m²)(Σvᵢ + Σ_{i≠j} ρ√(vᵢvⱼ))`, ρ = 0.6 (sensitivity ρ ∈ {0.3, 0.8}).
Then random-effects REML with **Knapp-Hartung-Sidik-Jonkman (HKSJ)** standard errors, which
is the recommended small-k adjustment. Naive per-effect pooling is reported alongside purely to
show how much dependency handling matters; it is *not* the headline number.
Where a stratum has ≥ 4 cohort families, cluster-robust variance estimation
(Hedges-Tipton-Johnson, with the small-sample df correction) is also reported.

## A.4 Statistical tests, thresholds, multiplicity

- Two-sided α = 0.05 throughout. Primary confirmatory tests: the three specificity contrasts
  (stage vs TST, stage vs SE, stage vs spectral). Holm correction across those three.
- Heterogeneity: τ² (REML), Q, I², and a 95% **prediction interval**.
- Small-study bias: Egger-type regression of z on its standard error, and trim-and-fill (L0),
  reported only where k ≥ 5 and labelled as uninterpretable below that.
- Given k ≤ 6 per stratum, every primary p-value is accompanied by an **exact permutation
  p-value** over sign flips of the aggregated effects.
- **Equivalence testing (TOST)** on the N3-stage stratum against bounds r = ±0.10 and ±0.20.

## A.5 Success criteria (what would support vs refute the hypothesis)

The hypothesis has three conjunctive clauses; it is preregistered as **supported only if all
three hold**:

| Clause | Supported if | Refuted if |
|---|---|---|
| (i) N3 stage % is associated with amyloid PET burden | pooled r > 0 with 95% CI excluding 0 | CI includes 0 and excludes r ≥ 0.20 (TOST) |
| (ii) the association is dose-dependent and measurable in SUVR | ΔSUVR per 10 pp N3 CI excludes 0 and is a non-trivial fraction of the A-/A+ separation | CI includes 0, or magnitude is < ~5% of the A-/A+ SUVR gap |
| (iii) it is *stronger* than TST and SE | paired Δz = z_N3 - z_TST > 0, CI excluding 0, **and** the excess survives the E3 reliability-attenuation null | Δz ≤ 0, or any positive gap is fully explained by differential ICC |

A negative result on any clause is reported as the finding, not treated as a failed experiment.

## A.6 What could go wrong, and the contingency

- *k too small for stable τ² estimation.* Mitigated by HKSJ + permutation + reporting τ² as
  unstable rather than interpreting it. Contingency: report fixed-effect estimates alongside.
- *t-recovery conversion is wrong for some beta.* Mitigated by running the entire primary
  analysis with those rows excluded (preregistered sensitivity S2); if the conclusion flips, the
  conclusion is reported as conversion-dependent.
- *Imputed nulls dominate a stratum.* Mitigated by sensitivity S1 (drop them) and by reporting
  the imputed fraction per stratum.
- *The answer is "we cannot tell".* This is a legitimate outcome and is what E6 is for; it would
  be reported as an evidence-insufficiency result with the detectable-effect bound stated.

## A.7 Milestones

| Milestone | Est. | Artifact |
|---|---|---|
| M1 harmonisation + QC | 25 min | `results/harmonised_effects.csv`, `results/qc_report.json` |
| M2 E1 pooling | 25 min | `results/e1_pooled.json` |
| M3 E2 specificity contrasts | 30 min | `results/e2_specificity.json` |
| M4 E3 attenuation null | 20 min | `results/e3_attenuation.json` |
| M5 E4 dose conversion | 20 min | `results/e4_dose_response.json` |
| M6 E5/E6 robustness + power | 35 min | `results/e5_robustness.json`, `results/e6_precision.json` |
| M7 figures | 20 min | `figures/*.png` |
| M8 REPORT.md + README.md | 30 min | `REPORT.md`, `README.md` |

Buffer for debugging: ~30% of the above.

---


**Hypothesis under test.** Reduced slow-wave sleep duration shows a dose-dependent association
with higher amyloid-beta burden in cognitively normal adults; each 10% reduction in N3 sleep
percentage is associated with a measurable increase in amyloid PET SUVR; and this relationship
is stronger than associations with total sleep time (TST) or sleep efficiency (SE) alone.

## What the literature already says (Phase 1 evidence, see `literature_review.md`)

Three facts constrain everything downstream:

1. **The stage/spectral dissociation.** In the single dataset that reports both, Winer et al.
   (2020, *Curr Biol*) found baseline proportion of <1 Hz slow-wave activity predicted the rate
   of cortical amyloid accumulation at r = -0.52 (p = 0.002), while **N3 stage duration in the
   very same participants was null (r = -0.07, p = 0.72)**. Montagne et al. (2026) similarly
   found N3% did not differ between amyloid-positive and amyloid-negative participants
   (ηp² = 0.003, p = 0.628) while TST did (ηp² = 0.065, p = 0.031). The hypothesis names the
   *stage* metric; the literature's positive signal sits on the *spectral* metric.

2. **The positive spectral literature is largely one cohort.** Mander 2015, Winer 2019 and
   Winer 2020 all draw on the Berkeley Aging Cohort Study with overlapping participants.
   Champetier et al. (2024, n = 127, Age-Well/IMAP) attempted the replication in an independent
   cohort with delta sub-band decomposition and found **no** association with amyloid burden at
   baseline or with accumulation over ~21 months.

3. **No prior meta-analysis has pooled slow-wave-specific effects against amyloid PET.**
   Chen et al. (2025, *Alzheimers Dement*, 30 studies / 14,997 subjects) pooled only *sleep
   quality* and *sleep duration*; Cui et al. (2022) and Bubu et al. (2020) covered obstructive
   sleep apnea. Chen et al. additionally reported I² = 99% for the sleep-quality/amyloid-PET
   pooling and appear to have pooled raw regression coefficients under a "Fisher z" label
   (e.g. z = -0.002, 95% CI -0.003 to -0.001), which is not a defensible common metric.

## Direction budget

Directions were scored on four criteria (1-5 each): **Evidence** (literature support that the
direction is answerable), **Relevance** (bears directly on the stated hypothesis),
**Information gain** (how much the answer moves belief), **Feasibility** (executable with the
resources actually staged in this workspace).

| # | Direction | Ev | Rel | IG | Feas | Total | Decision |
|---|-----------|----|-----|----|------|-------|----------|
| D1 | Stage-vs-spectral dissociation: pool N3/SWS stage %, low-frequency SWA, TST and SE against amyloid PET separately, with cohort-clustered variance | 5 | 5 | 5 | 5 | **20** | **KEEP** |
| D2 | Measurement-reliability confound: test whether exposure-metric differences are explained by differing single-night reliability rather than biology | 4 | 5 | 5 | 5 | **19** | **KEEP** |
| D3 | Dose-response conversion to SUVR/Centiloid units per 10% N3, plus heterogeneity, bias and small-study diagnostics | 4 | 5 | 4 | 4 | **17** | **KEEP** |
| D4 | Glymphatic MRI (DTI-ALPS) mediation meta-analysis | 4 | 2 | 3 | 3 | 12 | reject |
| D5 | CPAP / OSA-treatment intervention effects on amyloid | 3 | 2 | 3 | 3 | 11 | reject |
| D6 | CSF Aβ42 kinetics as the primary outcome | 4 | 2 | 2 | 4 | 12 | reject (retained as secondary stratum only) |
| D7 | Train an automated sleep stager on public PSG to re-derive N3% | 3 | 1 | 1 | 3 | 8 | reject |
| D8 | Individual-participant-data meta-analysis (ADNI / A4 / OASIS-3) | 5 | 5 | 5 | 1 | 16 | reject (blocked) |

### Kept directions

**D1 - Stage-vs-spectral dissociation (primary).**
Pool `SWS_stage`, `SWS_spectral`, `TST` and `SE` effects against amyloid PET as four separate
strata on the Fisher-z scale, then formally contrast the strata. Every effect carries a
`cohort_family` label so that Berkeley (Mander 2015 / Winer 2019 / Winer 2020) and Age-Well
(Champetier 2024 / Montagne 2026) are not treated as independent. The pilot run
(`artifacts/pilot_meta.json`) already shows the shape of the answer: spectral SWA r = 0.44 but
from a **single** cohort family; N3/SWS stage r = 0.01 across two cohort families; TST r = 0.30
across three independent cohort families. The headline deliverable is a direct test of the
hypothesis's specificity claim, which the preliminary evidence points *against* as literally
worded.

**D2 - Measurement-reliability confound.**
A stronger observed correlation for one exposure can arise purely because that exposure is
measured more reliably. We derived single-night test-retest reliabilities from Sleep-EDF
Expanded (75 subjects with two nights each): N3% ICC(2,1) = 0.775, sleep efficiency 0.562,
TST 0.225 (`datasets/derived/night_to_night_reliability.csv`). Attenuation factors √ICC are
0.88, 0.75 and 0.47. So if the *true* correlations were identical, single-night N3% would show
an observed correlation ~1.9x that of TST. This gives a quantitative null model: the hypothesis
survives only if the N3-vs-TST gap exceeds what differential reliability alone predicts.
The observed direction is the opposite (N3 ≈ 0, TST ≈ 0.30), which is a strong result in itself.

**D3 - Dose-response conversion and bias diagnostics.**
Convert pooled standardized effects into the units the hypothesis states. Using the
age-stratified N3% dispersion from Sleep-EDF (SD ≈ 8-10 percentage points in adults over 60,
`datasets/derived/n3_norms_by_age.csv`), a 10-percentage-point reduction in N3% is roughly
1.0-1.2 SD. Combined with the amyloid PET reference scale (amyloid-positivity threshold
1.42 PiB SUVR = 19 Centiloid; A- median 1.31 SUVR vs A+ median 1.57 SUVR,
`datasets/meta_analysis/amyloid_pet_reference_scale.csv`), a pooled r maps to an expected SUVR
shift per 10% N3. Accompany with τ², I², prediction intervals, leave-one-cohort-out sensitivity,
and small-study diagnostics (Egger-type regression, trim-and-fill; see `Lin2018`, `Moreno2009`).

### Why the rejected directions were pruned

- **D4 (DTI-ALPS)** measures a perivascular-diffusion proxy, not a sleep stage. It is one step
  removed from the hypothesis's exposure and the index itself is contested. Papers are retained
  as mechanistic background (Kamagata 2022, Huang 2024, Hsu 2023, Hong 2024) but not analysed.
- **D5 (CPAP/OSA)** changes the exposure to apnea severity and adds intermittent hypoxemia as a
  competing mechanism. Carvalho 2024 shows the sign of the SWA-amyloid relation can even
  *reverse* in OSA, so mixing these in would confound the primary contrast rather than inform it.
- **D6 (CSF Aβ42)** is not interchangeable with PET SUVR, and its sign relative to amyloid burden
  flips with disease stage: Varga 2016 (cognitively normal) found less SWS with *higher* CSF
  Aβ42, whereas Liguori 2020 (cohort including AD) found less N3 with *lower* CSF Aβ42. Kept in
  the dataset with an explicit `burden_aligned_note` and analysed only as a stage-stratified
  secondary sensitivity check.
- **D7** would consume most of the compute budget and produce no amyloid-linked inference: the
  public PSG corpora carry no amyloid labels.
- **D8** is the strongest design available in principle, but ADNI, A4 and OASIS-3 all require
  data-use agreements with human review that cannot be completed inside this pipeline. This is a
  genuine capability gap, not a modelling choice, and is recorded as such in `resources.md`.

## Concrete next steps for the experiment runner

1. Load `datasets/meta_analysis/effect_sizes.csv`. Filter to
   `outcome_class in {deposition_pet, deposition_pet_change}` for the primary analysis.
2. Convert every poolable effect to Fisher z. Conversions already supported:
   `pearson_r`, `partial_r`, `kendall_tau` (r = sin(πτ/2)), `r_from_partial_eta2` (r = √ηp²).
   Rows typed `unstd_beta_*` and `narrative_null` need a documented handling rule - do not
   silently drop them; the `narrative_null` rows encode genuine reported nulls and omitting
   them would bias the pooled estimate upward.
3. Fit random-effects models per `exposure_domain` with variance clustered on `cohort_family`
   (or run a one-effect-per-cohort-family sensitivity analysis - `artifacts/pilot_meta.json`
   flags each stratum where this matters).
4. Contrast the strata; report the specificity claim explicitly.
5. Apply the D2 attenuation null model and the D3 unit conversion.
6. Report heterogeneity, prediction intervals and small-study diagnostics.

Re-open the direction ranking only if new evidence overturns points 1-3 above; if that happens,
record the revision and its reason in `STATE.md`.
