# Literature Review: Slow-Wave Sleep, Glymphatic Clearance and Amyloid-Beta Burden

Scope: 54 full-text papers plus 18 abstract-only records, gathered 2026-09-04.
Search was run against Europe PMC (the paper-finder service at `localhost:8000` was not
running; see `resources.md` for the fallback strategy). Every numeric claim below is quoted
from the paper's own text and carries a provenance string in
`datasets/meta_analysis/effect_sizes.csv`.

---

## 1. Research area overview

The proposed causal chain has four links: (i) slow-wave sleep drives glymphatic/perivascular
convective flow; (ii) that flow clears interstitial amyloid-β; (iii) chronically reduced
slow-wave sleep therefore permits amyloid accumulation; (iv) accumulation is measurable as
increased amyloid PET SUVR. Links (i) and (ii) are well supported mechanistically. Link (iii)
in humans is where the evidence becomes contested, and the contest turns on *which* measure of
slow-wave sleep is used.

**Mechanism.** Xie et al. (2013, *Science*) showed convective interstitial-fluid exchange
increases during natural sleep and anaesthesia in mice, with a ~60% increase in interstitial
space and faster clearance of exogenous amyloid-β. Fultz et al. (2019, *Science*) provided the
human counterpart: during NREM sleep, coupled electrophysiological, haemodynamic and CSF
oscillations produce large macroscopic CSF flow pulsations locked to slow-wave activity.
Holth et al. (2019, *Science*) extended the same logic to tau. Rasmussen et al. (2018,
*Lancet Neurol*) is the standard review of the glymphatic pathway across disorders.

**Human causal manipulation.** Ju et al. (2017, *Brain*) is the closest thing to an experiment.
Seventeen healthy adults underwent acoustic slow-wave-activity disruption versus sham in a
crossover design. Specific SWA disruption correlated with an increase in CSF amyloid-β40
(r = 0.610, p = 0.009), and the authors state the effect "was specific for slow wave activity,
and not for sleep duration or efficiency." This is the single strongest piece of support for
the specificity part of the hypothesis - but the outcome is *soluble CSF Aβ over one night*,
not deposited plaque on PET.

---

## 2. Key papers

### 2.1 Mechanistic foundation

| Paper | Design | Key result |
|---|---|---|
| **Xie 2013**, *Science*, 10.1126/science.1241224 | mouse, 2-photon | Sleep/anaesthesia increases convective interstitial exchange; faster Aβ clearance |
| **Fultz 2019**, *Science*, 10.1126/science.aax5440 | human fMRI+EEG | NREM slow waves drive macroscopic CSF flow oscillations |
| **Holth 2019**, *Science*, 10.1126/science.aav2546 | mouse + human | Sleep-wake cycle regulates interstitial tau and CSF tau |
| **Rasmussen 2018**, *Lancet Neurol* | review | Glymphatic pathway across neurological disease |

### 2.2 Slow-wave-specific human studies (the hypothesis's core evidence)

**Mander et al. 2015**, *Nat Neurosci*, 10.1038/nn.4035. n = 26 cognitively normal older adults
(18 female, 75.1 ± 3.5 y), PSG + PiB PET. Proportion of medial-prefrontal NREM SWA in the
0.6-1 Hz band correlated with mPFC PiB DVR at r = -0.45 (p = 0.020; Kendall τ = -0.30,
p = 0.035). Critically the effect was **frequency-specific**: the 1-4 Hz band went the other way
(t = +1.85, p = 0.076), and REM delta 0.6-1 Hz was not significant (r = -0.31, p = 0.123).

**Winer et al. 2019**, *J Neurosci*, 10.1523/jneurosci.0503-19.2019. n = 31, Berkeley Aging
Cohort Study. Replicated Mander: proportion 0.6-1 Hz SWA vs cortical PiB r = -0.36 (p = 0.04),
rising to r = -0.43 after adjusting for age, sex, apnea risk and scan interval. Established a
double dissociation - SWA tracks amyloid, slow-oscillation/spindle coupling tracks tau
(SWA vs MTL tau r = 0.05, p = 0.80).

**Winer et al. 2020**, *Curr Biol*, 10.1016/j.cub.2020.08.017. n = 32, longitudinal. **The single
most important study for this hypothesis.** Baseline proportion <1 Hz SWA predicted the rate of
subsequent cortical Aβ accumulation (r = -0.52, p = 0.002). In the same participants:
sleep efficiency r = -0.35 (p = 0.05) prospectively and r = -0.43 (p = 0.01) cross-sectionally;
TST r = -0.36 (p = 0.04, did not survive correction); and **N3/SWS stage duration r = -0.07
(p = 0.72)**, N2 r = -0.32 (p = 0.08), REM r = -0.18 (p = 0.34). The authors conclude the signal
indicates "specificity to the efficiency of sleep, rather than any specific sleep stage."

**Varga et al. 2016**, *Sleep*, 10.5665/sleep.6240. n = 36 cognitively normal (66.8 ± 8.2 y),
PSG + lumbar puncture. %SWS vs CSF Aβ42 r = -0.36; SWS duration r = -0.35; total frontal SWA
r = -0.45 (p < 0.01). TST, N1, N2 and REM were all uncorrelated. Note the direction: in this
pre-deposition cohort *less* SWS went with *higher* CSF Aβ42, interpreted as elevated soluble
production rather than lower plaque load.

**Champetier et al. 2024**, *Sleep*, 10.1093/sleep/zsae012. n = 127 cognitively unimpaired,
Age-Well/IMAP. The most direct independent replication attempt of the Berkeley finding, with
delta decomposed into sub-bands: "Delta power was neither associated with amyloid burden at
baseline nor its accumulation over time, whatever the frequency sub-band." **Null.**

**Carvalho et al. 2024**, *Brain Commun*, 10.1093/braincomms/fcae354. n = 64, Mayo Clinic Study
of Aging, older adults with OSA, ≥2 serial PiB scans. Per 1 SD, slow-oscillation relative power
*increased* annualized PiB accumulation by 0.0033 (95% CI 0.0001-0.0064, p = 0.042) and
SO-slope by 0.0069 (0.0009-0.0129, p = 0.026) - **the opposite direction to the hypothesis** -
while delta-slope went the hypothesised way (-0.0082, -0.0143 to -0.0021, p = 0.009). In OSA,
elevated slow-wave measures may index unmet sleep pressure and hyper-excitability.

**Lucey et al. 2019**, *Sci Transl Med*, 10.1126/scitranslmed.aau6550. n = 38 with single-channel
EEG. NREM SWA inversely related to AD pathology, most strongly for **tau** and at the lowest
frequencies (1-2 Hz), with amyloid associations weaker.

**Stankeviciute et al. 2024**, *Alzheimers Dement (Amst)*, 10.1002/dad2.12616. n = 39. Decreased
SWS associated with higher temporal tau and lower cortical thickness; explicitly framed as
*amyloid-independent* - "adjustment for cortical amyloid burden did not alter the associations
between percentages of N1 and N3 sleep with tau."

**Montagne et al. 2026**, *Alzheimers Dement*, 10.1002/alz.71376. n = 76 with PSG (Aβ- n = 54,
Aβ+ n = 22), cognitively unimpaired, Age-Well. Amyloid-status ANCOVA adjusted for age, sex,
education and AHI: **N3% ηp² = 0.003, p = 0.628** (Aβ+ actually slightly higher: 17.4% vs
16.5%); N3 minutes ηp² = 0.003, p = 0.659; **TST ηp² = 0.065, p = 0.031** (369.0 vs 397.1 min);
SE ηp² ≈ 0, p = 1.00. The stage-specific claim fails while the duration claim holds.

**Jin et al. 2025**, *Alzheimers Dement*, 10.1002/alz.14495. n = 128. SWS% tertiles vs Aβ PET
SUVR were null in fully adjusted models (middle tertile β = -0.26, 95% CI -0.48 to 0.06,
p = 0.102; top tertile β = 0.15, -0.19 to 0.65, p = 0.427). REM latency, not SWS, was the metric
associated with Aβ PET SUVR (β = 0.08, 0.03-0.13, p = 0.002).

**Zavecz et al. 2023**, *BMC Medicine*, 10.1186/s12916-023-02811-z. n = 62. Reframes NREM SWA as
a *cognitive reserve* factor: SWA moderated the effect of Aβ status on memory
(std. β = 0.64, p = 0.042) rather than predicting amyloid level itself.

### 2.3 Sleep duration, efficiency and quality vs amyloid

| Paper | n | Exposure | Result |
|---|---|---|---|
| **Winer 2021**, JAMA Neurol | 4,417 (A4) | self-reported duration | shorter sleep → higher Aβ, β = -0.01 SUVR/h, p = .005 |
| **Spira 2013**, JAMA Neurol | 70 (BLSA) | self-reported duration | cortical DVR B = 0.08 (0.03-0.14), p = .005; precuneus B = 0.11, p = .007; partial r = 0.32 excluding MCI/dementia |
| **Insel 2021**, JAMA Netw Open | 4,425 (A4) | nighttime duration | null for global Aβ; regional effects for daytime sleep only |
| **Pivac 2024**, Alz Dement (Amst) | 189 (AIBL) | SE < 65% vs > 85% | faster accumulation, β = 2.91 ± 0.93, p = .002 |
| **Brown 2016**, Sleep | 184 (AIBL) | PSQI components | only sleep latency significant (B = 0.003, p = .02) |
| **Fenton 2023**, Brain Commun | 52 | actigraphic *variability* | SE variability b = 0.03, p < .01; duration variability b = 0.01, p = .04 |
| **Liguori 2020**, Alz Res Ther | 258 (mixed SCI→AD) | PSG | N3 vs CSF Aβ42 r = +0.362; TST +0.306; SE +0.403 |

### 2.4 Prior systematic reviews and meta-analyses

**Chen et al. 2025**, *Alzheimers Dement*, 10.1002/alz.70096. 30 studies, 14,997 subjects.
Poor sleep quality vs amyloid PET: Fisher z = 0.153 (k = 9, 95% CI 0.005-0.302, p = .042),
heterogeneity I² = 99%. Sleep duration vs amyloid PET: z = -0.002 (k = 6, -0.003 to -0.001,
p = .002). **Neither slow-wave sleep nor N3 percentage was pooled as an exposure.** The
reported "Fisher z" values for duration are implausible on that scale and appear to be pooled
raw regression coefficients - a metric-harmonisation weakness our analysis should avoid.

**Bubu et al. 2017** (*Sleep*, 476 citations) and **Bubu et al. 2020** (*Sleep Med Rev*) cover
sleep/cognition/AD risk and OSA respectively. **Cui et al. 2022** (*Front Aging Neurosci*)
meta-analyses AD biomarkers in OSA.

### 2.5 Measurement-methods papers

- **Levendowski 2017**, *J Clin Sleep Med* - accuracy, night-to-night variability and stability of frontopolar sleep EEG.
- **Nikkonen 2024**, *J Sleep Res* - multicentre sleep-stage scoring agreement (inter-scorer reliability ceiling).
- **Jack 2017**, *Lancet Neurol* - amyloid PET SUVR/Centiloid reference distributions and the 1.42 SUVR / 19 Centiloid positivity threshold.
- **Hanseeuw 2021**, *Eur J Nucl Med Mol Imaging* - Centiloid = 26 threshold for progression to dementia.
- **Lin 2018**, *Biometrics*; **Moreno 2009**, *BMC Med Res Methodol* - publication-bias quantification and regression-based adjustment.

---

## 3. Synthesis

### Common methodologies
- **Exposure**: overnight lab PSG (Mander, Varga, Winer, Champetier, Montagne, Jin, Liguori); single-channel/forehead EEG (Lucey 2019); actigraphy (Ju 2017 home nights, Fenton, Spira 2014); self-report PSQI/MOS/WHIIRS (Spira 2013, Sprecher, Brown, Winer 2021, Insel).
- **Outcome**: PiB or florbetapir PET expressed as DVR/SUVR/Centiloid, either global cortical composite or region-specific (mPFC, precuneus); CSF Aβ42 or Aβ40 by ELISA; a minority use plasma.
- **Spectral exposure definitions vary a lot** - "proportion 0.6-1 Hz SWA", "prop <1 Hz SWA", "SO relative power 0.5-0.9 Hz", "delta 1-3.9 Hz", "delta downslope". This is a major heterogeneity source and must be a moderator.

### Standard baselines and comparators
Sleep duration and sleep efficiency are the field's de-facto comparators - which is exactly what
the hypothesis proposes to beat. Prior pooled benchmarks to beat or contextualise:
Chen 2025 sleep-quality/amyloid-PET Fisher z = 0.153, and Winer 2021's β = -0.01 SUVR per hour
of sleep in n = 4,417.

### Evaluation metrics
Fisher-z-transformed correlation as the common metric; τ², Q and I² for heterogeneity;
prediction intervals; leave-one-cohort-out sensitivity; Egger-type regression and trim-and-fill
for small-study effects; and for the dose-response statement, SUVR or Centiloid change per
10 percentage points of N3.

### Datasets used in the literature
Berkeley Aging Cohort Study; A4 screening (n ≈ 4,400); BLSA neuroimaging substudy;
Wisconsin Registry for Alzheimer's Prevention; AIBL; Mayo Clinic Study of Aging;
Age-Well/IMAP (Caen); ADNI. **None of these are openly downloadable**; all require a
data-use agreement. This is why the analysis dataset here is a curated study-level effect-size
table rather than individual-participant data.

---

## 4. Gaps and opportunities

1. **No meta-analysis has ever pooled slow-wave-specific exposures against amyloid PET.**
   This is the gap the hypothesis occupies, and it is real.
2. **Stage vs spectral conflation.** Papers and reviews routinely write "slow-wave sleep" for
   both N3 stage percentage and low-frequency spectral power. Winer 2020 shows these dissociate
   *within one sample* (r = -0.52 vs r = -0.07). Any pooled estimate that merges them is
   uninterpretable.
3. **Cohort non-independence.** The positive spectral literature is concentrated in the Berkeley
   Aging Cohort Study across at least three publications with overlapping participants. Naive
   pooling would treat one cohort as three and badly overstate precision.
4. **Direction reversal by disease stage and by comorbidity.** CSF Aβ42's sign relative to
   amyloid burden flips between pre-deposition (Varga 2016) and post-deposition (Liguori 2020)
   samples; the SWA-amyloid sign flips in OSA (Carvalho 2024).
5. **Differential measurement reliability is an untested confound.** No paper in this corpus
   asks whether metric-to-metric differences in observed correlation strength are explained by
   metric-to-metric differences in single-night reliability. We can now test this directly.
6. **The "10% N3" dose statement has never been estimated.** No study reports amyloid PET SUVR
   change per unit of N3 percentage.

---

## 5. Recommendations for the experiment

**Primary dataset.** `datasets/meta_analysis/effect_sizes.csv` - 49 effects from 17 studies
across 12 independent cohort families, each with a verified DOI/PMID, a provenance quote, a
`burden_aligned_effect` normalised so that positive means "worse sleep ↔ more amyloid", and a
`cohort_family` label for dependency handling.

**Supporting datasets.**
`datasets/derived/night_to_night_reliability.csv` (empirical single-night ICCs for the
attenuation null model), `datasets/derived/n3_norms_by_age.csv` (N3% dispersion for the
dose-response conversion), `datasets/meta_analysis/amyloid_pet_reference_scale.csv` (SUVR ↔
Centiloid anchors), `datasets/meta_analysis/prior_meta_analyses.csv` (benchmarks).

**Recommended analysis.**
1. Random-effects (REML) pooling of Fisher-z per `exposure_domain`, restricted to amyloid PET.
2. Cluster or subsample on `cohort_family`; report both naive and clustered estimates.
3. Formally contrast `SWS_stage` vs `TST` vs `SE` - this is the hypothesis's specificity claim.
4. Contrast `SWS_spectral` vs `SWS_stage` - the dissociation is the real scientific story.
5. Apply the reliability-attenuation null model (D2) before attributing any gap to biology.
6. Convert to SUVR per 10% N3 using the norms and reference scale (D3).
7. Report τ², I², prediction intervals, leave-one-cohort-out, and small-study diagnostics.

**Methodological cautions.**
- Do not drop the `narrative_null` rows; they encode real reported nulls and omitting them
  biases the pooled estimate upward.
- Do not pool CSF Aβ42 with PET SUVR without stage stratification - the sign is not stable.
- Do not treat `unstd_beta_*` rows as if they were correlations; they need an explicit,
  documented conversion or a separate narrative synthesis.
- The likely honest conclusion is that the hypothesis **as literally worded** (N3 stage
  percentage, dose-dependent, stronger than TST/SE) is **not supported**, while a refined
  version restricted to low-frequency spectral SWA has suggestive but cohort-confounded support
  that an independent replication (Champetier 2024) failed to confirm. Report that outcome
  plainly rather than searching for a specification that rescues the original claim.
