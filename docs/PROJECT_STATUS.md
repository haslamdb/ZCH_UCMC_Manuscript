# ZCH-UCMC NICU Microbiome & Resistome Study - Project Status

**Project:** Multi-Site NICU Microbiome Analysis (Cincinnati & Hangzhou)
**Type:** Research Manuscript / Multi-Institutional Analysis
**Last Updated:** 2026-05-05

---

## Current Status

**Phase:** RADAR re-validation v2 — full clinical-covariate analysis layer in place
**Priority:** High (national meeting talk deadline ~1-2 weeks from 2026-05-05)

### Active sub-task: v2 analysis pipeline complete

The 642-sample biogpu/RADAR rerun has been aggregated, canonicalized, and run through a full clinical-covariate analysis pipeline. Ready for figure curation and manuscript update.

---

## Major work completed 2026-05-05

### 1. Sample-name canonicalization (data integrity fix)

Two collision bugs were silently dropping all 240 ZCH AMR samples from analysis:

- The metadata `sample_name` column was set from the unprefixed `Sample` column (`N01_1_2`) instead of `SampleID` (`ZJH_N01_1_2`), so ZCH metadata→AMR joins didn't match.
- ZCH single-digit subjects were zero-padded (`N1` → `N01`) at manifest stage, but the raw FASTQs use the single-digit form.

**Fixed by `scripts/canonicalize_zjh_names.py` and `scripts/canonicalize_radar_dirs.sh`** (renamed 48 radar dirs, rewrote sample names in AMR/summary/manifest/sample-key files; original copies preserved as `*.pre_canonicalize`). Companion fixes in `01_create_master_metadata.py` (drop the bad pad) and `02_quality_control.py` (group on SubjectID, not PatientID, since UCMC and ZCH share PatientID strings for N10–N59).

**Net effect:** v2 cohort went from a buggy UCMC-only 347 samples to the real 642 (359 UCMC + 283 ZCH). `data/sample_summary.tsv` was also augmented with zero-rows for the 54 radar-processed but strict-zero samples so they aren't dropped by selection bias.

### 2. Clinical-covariate analysis layer (the new manuscript scaffold)

Replaced the old "UCMC vs ZCH" framing with per-feature linear mixed-effects models on:
- BodySite (Axilla / Groin / Stool)
- SampleCollectionWeek (postnatal age proxy: Week 1 vs Week 3)
- PostNatalAbxCohort (No / Low / High infant antibiotic exposure)
- MaternalAntibiotics (Mat.Abx vs None)
- AnyMilk (Mother vs Donor)
- GestationTime (continuous)
- SubjectID random intercept

Same model run at four resolutions (helper at `scripts/python/_clinical_lmm.py`):
- **04d** gene_family (237 features)
- **04c** ARO id (117 features) — extends a descriptive panel with an inferential block
- **04b** MEGARes class (23 features)
- **12** clinical category (11 features) — uses the screening-category definitions from script 11

### 3. Clinical-screening view (where to swab)

`scripts/python/11_functional_class_screening.py` defines 11 clinical AMR categories (ESBL, carbapenemase, mecA/MRSA, aminoglycoside, fluoroquinolone, MLS, tetracycline, sulfonamide, trimethoprim, colistin/mcr, VRE) by `gene_symbol` regex + `subclass` fallback. For each category we compute per-subject × per-body-site presence (RPM ≥ 1 in the strict file, any sample of subject at site). Outputs include:

- Prevalence by site, week, location, and the full 3-way (Site × Week × Location) breakdown
- Venn diagrams (one per category, subject overlap across the 3 sites)
- Site-capture heatmap (% of category-positive subjects caught by each single/dual/triple site combination)

### 4. Diversity and volcano figures

- `10_clinical_diversity_figures.py` — boxplots of richness / total RPM / Shannon by Site × Week and Site × PostNatalAbxCohort; volcano panel (5 effects from gene-family LMM)

### 5. Intrinsic vs acquired split

- `09b_intrinsic_acquired_split.py` and a new section in `09_create_summary_report.py` — splits resistome by NCBI `resistance_type` (intrinsic vs acquired); reports median RPM and acquired-fraction trajectory by site × week.

---

## Headline biological findings (from the v2 layer)

1. **Body site dominates everything.** 220/237 gene families differ Stool vs Axilla; 181 differ Groin vs Axilla. Effect sizes are 3–10× larger than any antibiotic effect.
2. **Postnatal age (Week 3 vs Week 1) is the second-largest axis** — broad colonization expansion across 180/237 gene families, mostly increases.
3. **Infant antibiotics suppress the resistome.** 85/237 families significant; 77 of those decrease (consistent with broad flora depletion). MEGARes-class view: 13/13 significant classes decrease with high-abx exposure. Carbapenemase is the only category that *increases* with high abx.
4. **Maternal antibiotics, donor milk, gestational age** show essentially no class- or category-level signal. The resistome is shaped by exposure variables (where, when, current antibiotics), not by host/constitutional ones.
5. **Geographic asymmetry, opposite directions, in two specific markers:**
   - Carbapenemase (NDM family): 0% UCMC, **39% ZCH** subjects
   - Vancomycin (VRE): **29% UCMC**, 0% ZCH subjects
   ESBL is higher in UCMC (73% vs 41%), surprisingly.
6. **Body site determines what you screen.** Axilla catches 90% of mecA carriers (skin organism) but misses VRE and mcr entirely. Stool catches 100% of VRE and mcr carriers and most ESBL/carbapenemase. Groin is intermediate but never dominant. **Axilla + Stool together capture ≥95% of every category** — the optimal 2-swab combo.
7. **mecA reverses by week in stool** (38.8% → 21.6%) while every other category rises — gut commensals are displacing the early-life skin-derived flora.

---

## Pipeline state

| Script | Purpose | Resolution |
|---|---|---|
| `01_create_master_metadata.py` | Canonicalize + merge | — |
| `02_quality_control.py` | Sample/gene QC | — |
| `03_exploratory_analysis.py` | PCA + diversity | gene_family |
| `04_differential_abundance.py` | UCMC vs ZCH (legacy) | gene_family |
| `04b_megares_class_differential.py` | Clinical LMM | MEGARes class (23) |
| `04c_aro_panel.py` | Descriptive + clinical LMM | ARO (117) |
| `04d_gene_family_clinical.py` | Clinical LMM | gene_family (237) |
| `05_antibiotic_correlations.py` | Spearman abx-AMR | total RPM |
| `06_longitudinal_analysis.py` | Week 1 vs 3 paired | gene_family |
| `07_mixed_effects_models.py` | Total AMR LMM | total RPM |
| `08_publication_figures.py` | Pub figs | gene_family |
| `09_create_summary_report.py` | Markdown report | — |
| `09b_intrinsic_acquired_split.py` | Intrinsic vs acquired | resistance_type |
| `10_clinical_diversity_figures.py` | Diversity + volcano | gene_family |
| `11_functional_class_screening.py` | Per-subject × site presence | clinical category (11) |
| `12_category_clinical_lmm.py` | Clinical LMM | clinical category (11) |
| `_clinical_lmm.py` | Shared LMM helper | — |

Settings doc: `nicu_resistome_analysis_v2/ANALYSIS_PIPELINE_SETTINGS.md` (covers strict filter, canonicalization, ARO/MEGARes integration).

---

## Upcoming work

- [ ] Curate manuscript figures from the v2 outputs (1–2 weeks for talk)
- [ ] Reviewer response v2 incorporating the new clinical-covariate framing
- [ ] Decide whether 04 (legacy UCMC-vs-ZCH framing) should be deprecated or kept as a supplementary contrast

---

## Session Log

| Date | Work Completed |
|------|----------------|
| 2026-02-03 | Created project status file |
| 2026-05-05 | RADAR v2 batch (642 samples) launched, aggregator written |
| 2026-05-05 | Sample-name canonicalization fix (recovers full 642-sample cohort) |
| 2026-05-05 | Clinical-covariate LMM layer at 4 resolutions (gene_family / ARO / MEGARes class / clinical category); intrinsic-vs-acquired split; clinical screening view with body-site Venns and 6-panel Site × Week × Location prevalence |
