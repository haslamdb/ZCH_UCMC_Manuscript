# NICU Resistome Analysis Report

Generated: 2025-11-10 17:32:53

---

## Study Overview

- **Total Samples:** 669
- **Total Subjects:** 127
- **Locations:** UCMC (Cincinnati) and ZCH (Hangzhou)
- **Body Sites:** Axilla, Groin, Stool
- **Time Points:** Week 1 and Week 3
- **AMR Genes Analyzed:** 237 genes (after QC filtering)


### Sample Distribution by Location

```
Location  Total Subjects  Complete Subjects  Total Samples  Samples with AMR  Subjects with Antibiotics
    UCMC              68                 40            386               359                         68
     ZCH              59                 13            283               247                         59
```


## Key Findings

### 1. No Significant Differences Between UCMC and ZCH

**Result:** Resistome composition is remarkably similar between Cincinnati and Hangzhou NICUs.

- Differential abundance testing: 0 genes significantly different at FDR < 0.05
- Mixed-effects model: Location effect p = 0.45 (not significant)
- **Interpretation:** Geographic location does not significantly influence NICU resistome composition


### 2. Body Site is the Primary Driver of Resistome Composition

**Result:** Stool samples have dramatically higher AMR burden than skin sites.

- PCA analysis: PC1 (30.8% variance) separates stool from axilla/groin
- Mixed-effects model: Stool coefficient = +0.615, p < 0.0001
- **Interpretation:** Gut-associated sites harbor more diverse and abundant AMR genes


### 3. Significant Longitudinal Changes in Groin and Stool

**Result:** AMR burden increases dramatically from Week 1 to Week 3 in groin and stool, but not axilla.


**Mean AMR RPM by Body Site and Week:**

```
Week            Week.1       Week.3
Body Site                          
Axilla      776.004280   694.102815
Groin       982.849792  2552.774607
Stool      3005.899737  4887.665049
```


**Statistical Tests (Paired Wilcoxon):**

- **Axilla:** +0.6% change, p = 0.9637 ns

- **Groin:** +203.8% change, p = 0.0001 ***

- **Stool:** +64.4% change, p = 0.0000 ***


**Gene-level changes:**

- **Groin:** 120 genes significantly changed (115 increased, 5 decreased)

- **Stool:** 120 genes significantly changed (104 increased, 16 decreased)


### 4. No Correlation Between Antibiotic Exposure and AMR Burden

**Result:** Cumulative antibiotic exposure does not correlate with total AMR burden.

- Overall Spearman rho = 0.001, p = 0.985
- No individual antibiotics showed significant correlations
- **Interpretation:** AMR colonization may be independent of selective antibiotic pressure in this cohort


### 5. Mixed-Effects Model Confirms Body Site × Week Interaction

**Model:** log(AMR) ~ Location × BodySite × Week + (1|SubjectID)


**Significant Effects:**

- **BodySite[T.Stool]:** coefficient = +0.615, p = 0.0000 ***

- **BodySite[T.Groin]:Week[T.Week.3]:** coefficient = +0.371, p = 0.0001 ***

- **Group Var:** coefficient = +0.355, p = 0.0017 **

- **BodySite[T.Stool]:Week[T.Week.3]:** coefficient = +0.284, p = 0.0029 **


**Intraclass Correlation (ICC):** 0.262

  - 26.2% of variance is between-subject
  - Confirms substantial within-subject correlation, validating repeated measures design


## Data Quality and Filtering

### Quality Control Steps:

1. **Gene Filtering:**

   - Original: 449 genes

   - After QC: 237 genes (present in ≥5% samples with max RPM ≥1.0)

2. **Sample Filtering:**

   - Excluded 'Nonese' (nose swab) samples (never sequenced)

   - Retained only samples with AMR data

3. **Subject Completeness:**

   - Complete subjects: 53 (40 UCMC, 13 ZCH)

   - Complete = all 6 samples present (3 body sites × 2 weeks)


## Statistical Methods

### Analyses Performed:

1. **Differential Abundance:** Mann-Whitney U test (UCMC vs ZCH, stratified by body site)

2. **Longitudinal Analysis:** Wilcoxon signed-rank test (paired, Week 1 vs Week 3)

3. **Correlation Analysis:** Spearman correlation (antibiotics vs AMR burden)

4. **Mixed-Effects Models:** Linear mixed-effects with subject random effects

5. **Multiple Testing Correction:** Benjamini-Hochberg FDR

6. **Dimensionality Reduction:** PCA on standardized AMR RPM matrix


## Analysis Outputs

### Results Files:

- `results/qc/` - Quality control and metadata

- `results/exploratory/` - PCA, diversity, clustering

- `results/differential/` - Differential abundance (UCMC vs ZCH)

- `results/correlations/` - Antibiotic-AMR correlations

- `results/longitudinal/` - Week 1→3 changes

- `results/mixed_models/` - Mixed-effects model results

- `results/summary/` - This report and summary tables


### Figures:

- `figures/publication/Figure1_PCA.pdf` - PCA overview

- `figures/publication/Figure2_BodySite_Comparison.pdf` - AMR by body site

- `figures/publication/Figure3_Longitudinal_Trajectories.pdf` - Week 1→3 trajectories

- `figures/publication/Figure4_Volcano_Plots.pdf` - Gene-level changes


## Conclusions

1. **No geographic differences:** UCMC and ZCH NICUs show remarkably similar resistome profiles

2. **Body site dominates:** Stool samples have the highest AMR burden, driven by gut microbiome composition

3. **Strong temporal dynamics:** Groin and stool sites show dramatic AMR increases from Week 1 to Week 3

4. **Temporal patterns are site-specific:** Axilla remains stable while groin/stool increase

5. **Antibiotic exposure not predictive:** No correlation between cumulative antibiotic days and AMR burden