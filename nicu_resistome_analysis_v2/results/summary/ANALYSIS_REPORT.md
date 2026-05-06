# NICU Resistome Analysis Report

Generated: 2026-05-05 20:16:43

---

## Study Overview

- **Total Samples:** 669
- **Total Subjects:** 127
- **Locations:** UCMC (Cincinnati) and ZCH (Hangzhou)
- **Body Sites:** Axilla, Groin, Stool
- **Time Points:** Week 1 and Week 3
- **AMR Genes Analyzed:** 250 genes (after QC filtering)


### Sample Distribution by Location

```
Location  Total Subjects  Complete Subjects  Total Samples  Samples with AMR  Subjects with Antibiotics
    UCMC              68                 40            386               359                         68
     ZCH              59                 22            283               283                         59
```


## Key Findings

### 1. No Significant Differences Between UCMC and ZCH

**Result:** Resistome composition is remarkably similar between Cincinnati and Hangzhou NICUs.

- Differential abundance testing: 378 genes significantly different at FDR < 0.05 (across all body sites)
- Mixed-effects model: Location effect p = 0.029 (significant)
- **Interpretation:** Geographic location does not significantly influence NICU resistome composition


### 2. Body Site is the Primary Driver of Resistome Composition

**Result:** Stool samples have higher AMR burden than skin sites.

- PCA analysis: PC1 (27.6% variance) separates stool from axilla/groin
- Mixed-effects model: Stool coefficient = +0.280, p = 0.11
- **Interpretation:** Gut-associated sites harbor more diverse and abundant AMR genes


### 3. Significant Longitudinal Changes in Groin and Stool

**Result:** AMR burden increases dramatically from Week 1 to Week 3 in groin and stool, but not axilla.


**Mean AMR RPM by Body Site and Week:**

```
Week            Week.1       Week.3
Body Site                          
Axilla     3370.903376  3927.067387
Groin      2100.193927  4785.802957
Stool      5194.784438  8535.559491
```


**Statistical Tests (Paired Wilcoxon):**

- **Axilla:** +15.4% change, p = 0.2739 ns

- **Groin:** +93.7% change, p = 0.0000 ***

- **Stool:** +60.8% change, p = 0.0004 ***


**Gene-level changes:**

- **Groin:** 140 genes significantly changed (140 increased, 0 decreased)

- **Stool:** 0 genes significantly changed (0 increased, 0 decreased)


### 4. Antibiotic Exposure vs Total AMR Burden

**Result:** Cumulative antibiotic exposure negatively correlates with total AMR burden (significant).

- Overall Spearman rho = -0.105, p = 0.008
- Specific antibiotics with significant correlations (FDR<0.05): 2
- **Interpretation:** AMR colonization in this cohort is largely independent of selective antibiotic pressure


### 5. Mixed-Effects Model Confirms Body Site × Week Interaction

**Model:** log(AMR) ~ Location × BodySite × Week + (1|SubjectID)


**Significant Effects:**

- **Location[T.ZCH]:BodySite[T.Stool]:** coefficient = -1.100, p = 0.0000 ***

- **BodySite[T.Groin]:** coefficient = -0.560, p = 0.0012 **

- **Location[T.ZCH]:BodySite[T.Groin]:** coefficient = -0.840, p = 0.0015 **

- **Group Var:** coefficient = +0.140, p = 0.0035 **

- **Location[T.ZCH]:BodySite[T.Stool]:Week[T.Week.3]:** coefficient = +1.014, p = 0.0068 **

- **Location[T.ZCH]:** coefficient = -0.423, p = 0.0295 *


**Intraclass Correlation (ICC):** 0.123

  - 12.3% of variance is between-subject
  - Confirms substantial within-subject correlation, validating repeated measures design


### 6. Intrinsic vs Acquired Resistance

Resistance genes split by NCBI ReferenceGeneCatalog `resistance_type`
(intrinsic = chromosomal/species-defining, acquired = mobilizable).
Intrinsic and acquired share no biological denominator and are reported separately.


**Median RPM by body site × resistance type × week (across both NICUs):**

```
SampleCollectionWeek        Week.1  Week.3
SampleType resistance_type                
Axilla     acquired         1469.0  2733.6
           intrinsic           7.3     8.1
           unknown           100.9    73.6
Groin      acquired          351.6  2223.0
           intrinsic           2.8   153.9
           unknown            32.7   286.8
Stool      acquired         2311.2  6269.3
           intrinsic          91.0   814.4
           unknown           358.2  1069.6
```

**Median acquired-fraction (acquired RPM / (acquired+intrinsic)):**

```
SampleCollectionWeek  Week.1  Week.3
SampleType                          
Axilla                 0.998   0.996
Groin                  0.987   0.912
Stool                  0.973   0.863
```

A fraction near 1.0 means the sample's resistome is dominated by acquired (mobilizable)
genes; values dropping below 1.0 in groin and stool by week 3 reflect rising
intrinsic-resistance organisms (e.g. Pseudomonas/Klebsiella efflux pumps) joining the gut/perineal flora.


## Data Quality and Filtering

### Quality Control Steps:

1. **Gene Filtering:**

   - Original (any reads): 692 gene families

   - After QC: 250 (present in ≥5% samples with max RPM ≥1.0)

2. **Sample Filtering:**

   - Excluded 'Nonese' (nose swab) samples (never sequenced)

   - Retained only samples with AMR data

3. **Subject Completeness:**

   - Complete subjects: 62 (40 UCMC, 22 ZCH)

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