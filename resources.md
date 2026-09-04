# Resources Catalog

Gathered 2026-09-04 for: *Glymphatic Clearance Efficiency During Sleep Stages and Its
Relationship to Amyloid-Beta Accumulation: A Meta-Analysis of Human Neuroimaging Studies.*

## Summary

| Resource | Count | Location |
|---|---|---|
| Papers with full text (PDF + extracted text) | 54 | `papers/` |
| Papers, abstract only | 20 | `paper_search_results/abstract_only.json` |
| Curated analysis datasets | 3 | `datasets/meta_analysis/` |
| Derived calibration datasets | 2 | `datasets/derived/` |
| Downloaded raw data | 197 hypnograms | `datasets/sleep_edfx/` |
| Cloned repositories | 4 | `code/` |
| Pipeline validation artifact | 1 | `artifacts/pilot_meta.json` |

---

## Papers

54 papers with full text, spanning glymphatic mechanism (Xie 2013, Fultz 2019, Holth 2019),
slow-wave-specific human studies (Mander 2015, Varga 2016, Ju 2017, Winer 2019/2020/2021,
Lucey 2019, Champetier 2024, Carvalho 2024, Montagne 2026, Jin 2025), sleep-duration and
sleep-quality studies (Spira 2013/2014/2021, Sprecher 2015/2017, Insel 2021, Brown 2016,
Pivac 2024, Fenton 2023, Naismith 2022), prior meta-analyses (Chen 2025, Cui 2022, Bubu
2019/2020), and measurement-methods papers (Jack 2017, Hanseeuw 2021, Levendowski 2017,
Nikkonen 2024, Lin 2018, Moreno 2009).

Full per-paper detail, including a "why relevant" line for each, is in `papers/README.md`.
Machine-readable bibliographic metadata is in **`papers/paper_metadata.csv`** (title, authors,
journal, year, DOI, PMID, PMCID, citations, OA status, file paths), sourced from Europe PMC
records rather than parsed from PDFs.

Extracted plain text for every paper is in `papers/fulltext/*.txt` — grep these rather than
re-reading PDFs.

---

## Datasets

| Name | Source | Size | Role | Location |
|---|---|---|---|---|
| **Effect sizes** | manual extraction from `papers/` | 49 effects, 17 studies, 12 cohort families | **primary analysis input** | `datasets/meta_analysis/effect_sizes.csv` |
| Amyloid PET reference scale | Jack 2017, Hanseeuw 2021 | 12 rows | SUVR ↔ Centiloid anchors for the dose-response conversion | `datasets/meta_analysis/amyloid_pet_reference_scale.csv` |
| Prior meta-analyses | Chen 2025 | 7 rows | pooled benchmarks to contextualise our estimates | `datasets/meta_analysis/prior_meta_analyses.csv` |
| Night-to-night reliability | derived from Sleep-EDF | 7 metrics, n = 75 subjects | attenuation null model (direction D2) | `datasets/derived/night_to_night_reliability.csv` |
| N3 norms by age | derived from Sleep-EDF | 5 age bands × 7 metrics | dose-response denominator (direction D3) | `datasets/derived/n3_norms_by_age.csv` |
| Sleep-EDF Expanded hypnograms | PhysioNet | 197 recordings, ~800 KB | raw input for the two derived tables | `datasets/sleep_edfx/` |

Distribution across the strata the hypothesis compares (amyloid-PET rows only):

| exposure_domain | cross-sectional PET | longitudinal PET change |
|---|---|---|
| `SWS_spectral` | 5 | 5 |
| `SWS_stage` | 5 | 1 |
| `TST` | 7 | 2 |
| `SE` | 4 | 2 |
| `WASO` | 1 | 1 |
| `other_stage` | 2 | 2 |

Full documentation, download instructions, loading code, caveats and sample records are in
**`datasets/README.md`**.

---

## Code Repositories

| Name | URL | Purpose | Status |
|---|---|---|---|
| PyMARE | github.com/neurostuff/PyMARE | meta-analysis engine | **installed + smoke-tested** |
| metafor | github.com/wviechtb/metafor | reference formulations (R) | doc only — R not installed |
| dosresmeta | github.com/alecri/dosresmeta | dose-response meta-analysis method (R) | doc only |
| YASA | github.com/raphaelvallat/yasa | sleep analysis toolbox | cloned, not installed (not needed) |

See `code/README.md`, which also records a PyMARE API gotcha (`get_fe_stats()` returns 2-D
arrays) and PyMARE's lack of cluster-robust variance — relevant because our effects are nested
within cohort families.

---

## Resource-gathering notes

### Search strategy
The paper-finder service at `localhost:8000` was **not running**, so the search fell back to
manual retrieval. Given the biomedical domain, Europe PMC was used as the primary index rather
than arXiv or Papers with Code:

1. 12 broad topic queries (`scripts/search_epmc.py`) → 345 unique records
2. 16 targeted queries by author and exact title for seminal work not surfaced by the broad
   sweep → 75 more
3. 12 supplemental queries for prior meta-analyses, PSG reliability and publication-bias methods
   → 61 more
4. 16 further queries for sleep-architecture studies named in Chen 2025's inclusion table → 25 more

Total pool 420 records; 74 selected against an explicit must-have list; 54 full texts obtained.

### Selection criteria
Priority went to (a) studies reporting a *slow-wave-specific* exposure against an amyloid
outcome, (b) studies reporting the comparator exposures the hypothesis must beat (TST, sleep
efficiency), (c) prior meta-analyses that set the benchmark, (d) mechanism papers establishing
the causal chain, and (e) measurement-methods papers needed for the reliability and
unit-conversion analyses. Negative and null findings were sought deliberately — Champetier 2024,
Montagne 2026 and Jin 2025 are all nulls, and excluding them would have produced a badly biased
corpus.

### Challenges encountered
- **paper-finder unavailable.** Fell back to Europe PMC + Unpaywall + PMC.
- **PMC blocks automated PDF retrieval.** Requests to `pmc.ncbi.nlm.nih.gov/articles/*/pdf/`
  return a JavaScript interstitial rather than the file. Worked around by resolving all
  Unpaywall OA locations first (recovering 33 publisher PDFs), then for the remainder extracting
  PMC full-text HTML and typesetting it locally (17 papers). Those 17 are flagged
  *(rendered from PMC full text)* in `papers/README.md`; their text is complete but figures and
  image-based tables are absent, so a figure-only number could have been missed.
- **Semantic Scholar rate-limited** (HTTP 429); used sparingly with backoff.
- **NSRR (`sleepdata.org`) unreachable** from this environment.
- **No R interpreter**, so metafor and dosresmeta are reference material only.

### Gaps and workarounds
- **No open dataset pairs PSG with amyloid PET.** ADNI, A4, BLSA, AIBL, WRAP, Mayo MCSA,
  Age-Well and OASIS-3 all require data-use agreements with human review. Consequence: the
  analysis must be a **study-level** meta-analysis, not individual-participant data. This caps
  what can be concluded — in particular, a genuine per-10%-N3 dose-response curve cannot be
  estimated directly and must be approximated by converting pooled standardized effects using
  external N3% dispersion and SUVR reference scales.
- **The `SWS_stage` vs amyloid-PET stratum is thin** (5 cross-sectional + 1 longitudinal effect,
  and only 3 are poolable as correlations). This is not an extraction failure — it reflects that
  few studies report N3 stage percentage against amyloid PET as a coefficient. The experiment
  runner should report this limitation rather than padding the stratum with CSF or
  sleep-quality proxies.
- **20 papers remain abstract-only.** None contribute an effect size; they were used for
  screening only.

---

## Recommendations for experiment design

1. **Primary dataset**: `datasets/meta_analysis/effect_sizes.csv`, filtered to
   `outcome_class in {deposition_pet, deposition_pet_change}`. Stratify by `exposure_domain`.
2. **Baselines / comparators**: `TST` and `SE` strata (the hypothesis's own comparators), plus
   the published benchmark of Chen 2025 (sleep-quality/amyloid-PET Fisher z = 0.153, I² = 99%).
3. **Metrics**: Fisher-z pooled correlation with REML τ²; Q and I²; prediction intervals;
   leave-one-cohort-out sensitivity; Egger-type regression and trim-and-fill; and for the dose
   statement, SUVR/Centiloid change per 10 percentage points of N3.
4. **Code to reuse**: PyMARE for pooling (already installed and validated);
   `scripts/validate_pipeline.py` as the working skeleton — it already loads the data, converts
   effect types, pools per stratum and flags cohort dependency.
5. **Non-negotiable analysis controls**:
   - cluster on `cohort_family` (Berkeley alone contributes 13 of 49 effects),
   - keep `narrative_null` rows in,
   - stratify CSF outcomes by cognitive status or exclude them from the primary model,
   - run the reliability-attenuation null model before attributing any N3-vs-TST gap to biology.

### What the preliminary run already shows

`artifacts/pilot_meta.json` (produced by `scripts/validate_pipeline.py`) contains a first pass:

| stratum | k | independent cohorts | pooled r | 95% CI | I² |
|---|---|---|---|---|---|
| spectral slow-wave activity vs amyloid PET | 5 | **1** | 0.444 | 0.367 – 0.515 | 0% |
| N3/SWS stage duration vs amyloid PET | 3 | 2 | **0.012** | −0.159 – 0.181 | 0% |
| total sleep time vs amyloid PET | 3 | 3 | 0.299 | 0.169 – 0.418 | 0% |
| sleep efficiency vs amyloid PET | 3 | 2 | 0.241 | −0.364 – 0.702 | 64.5% |

Read carefully, this points toward the hypothesis being **unsupported as literally worded**:
N3 stage percentage shows essentially no association with amyloid PET, and is *weaker* than both
TST and sleep efficiency rather than stronger. The one strong stratum is spectral slow-wave
activity — but all five of those effects come from a single cohort family, and the independent
replication attempt (Champetier 2024, n = 127) was null.

These are preliminary numbers from a smoke test on a small k, not the final analysis: the
cohort clustering, the null rows, the attenuation model and the bias diagnostics are all still
to be applied. But the experiment runner should expect to report a negative or heavily qualified
result and should design the analysis to characterise that finding rigorously, rather than
searching for a specification under which the original claim survives.
