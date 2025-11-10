# Fluoroquinolone Resistance Mutation Analysis Plan

**Created**: 2025-11-10
**Status**: Ready to implement

---

## Overview

Analyze fluoroquinolone resistance mutations across NICU samples to determine if allele frequencies vary by:
1. **Location** (UCMC vs ZCH)
2. **Body site** (Axilla, Groin, Stool)
3. **Time** (Week 1 vs Week 3)
4. **Antibiotic exposure** (even though no FQ antibiotics were given)

---

## Data Files

### Input Data
1. `data/nicu_fq_allele_frequencies.csv` - 14,009 rows of allele-level data
   - Key columns: `sample_name`, `species`, `gene`, `position`, `has_resistance_mutation`, `total_resistant_frequency`, `dominant_mutant_frequency`
2. `data/fq_mutation_prevalence.tsv` - Summary of mutation prevalence (281 mutations)
3. `data/fq_resistance_by_species.tsv` - Species-level resistance classification (5,999 rows)
4. `results/qc/master_metadata_with_qc.tsv` - Sample metadata with antibiotic exposure

### Key Species with FQ Resistance Data

**Species with substantial resistance data (prioritize for analysis):**

1. **Staphylococcus aureus** - 339 samples with resistance
   - Key mutations: parE Y470N (293 samples, 99.8% freq), parE D/S432V, parC E/R80S, gyrA G/S84L/Y
   - High frequency mutations indicating established resistance

2. **Klebsiella pneumoniae** - 214 samples with resistance
   - Key mutations: gyrA K154R (177 samples, 99.5% freq - nearly universal!), gyrA D/S83A/L/F/N, gyrA D/S87Y/N/A
   - Widespread chromosomal mutations

3. **Enterococcus faecium** - 180 samples with resistance
   - Key mutation: gyrA E/S87Y (100% in all samples, 99.8% freq)
   - Fixed resistance mutation

4. **Enterococcus faecalis** - 132 samples with resistance
   - Key mutation: gyrA E/S87K (100% in all samples, 99.7% freq)
   - Fixed resistance mutation

5. **Escherichia coli** - 73 samples with resistance (but LOW frequency ~1-5%)
   - Key mutations: parE D/S458A, parC A/G78V (mostly wildtype with rare resistant alleles)
   - Limited clinical significance due to low frequencies

6. **Klebsiella oxytoca** - 62 samples with resistance
   - Key mutations: parC S80R/I/T (common QRDR mutation), gyrA D87A/G/Y, gyrA T83I/S
   - Moderate frequency mutations

7. **Pseudomonas aeruginosa** - 35 samples with resistance
   - Key mutations: gyrA D/S87Y/N, gyrA A/D83V
   - Moderate to high frequency

8. **Streptococcus mitis** - 33 samples with resistance
   - Key mutations: parC S79I (24 samples, 95.5% freq), gyrA S81F
   - Note: Different amino acid numbering than Gram-negatives

9. **Streptococcus oralis** - 17 samples with resistance
   - Key mutations: gyrA S81F/Y
   - Limited sample size

**Note**: Serratia marcescens and Streptococcus pyogenes do not have FQ resistance mutation data in this dataset, though they were significantly different in pathogen-antibiotic analysis. S. pneumoniae detected but no resistance mutations observed.

---

## Analysis Plan

### Phase 1: Data Preparation and QC

**Objective**: Merge FQ data with metadata, filter to resistance mutations only

**Script**: `scripts/python/10_fq_data_preparation.py`

#### Steps:
1. Load `nicu_fq_allele_frequencies.csv`
2. Merge with `master_metadata_with_qc.tsv` to add:
   - Location (UCMC/ZCH)
   - BodySite (Axilla/Groin/Stool)
   - Week (Week 1/Week 3)
   - SubjectID
   - Antibiotic exposure variables
3. **Filter to resistance mutations**:
   - `has_resistance_mutation == True`
   - `total_resistant_frequency > 0` (exclude wildtype-only positions)
4. Create derived variables:
   - `mutation_id` = `{species}_{gene}_{position}_{wildtype_aa}->{mutant_aa}`
   - `is_high_frequency` = `total_resistant_frequency >= 0.75` (≥75% mutant)
5. QC checks:
   - How many samples have FQ resistance data?
   - How many unique resistance mutations detected?
   - Distribution by species

#### Outputs:
- `data/fq_resistance_filtered.tsv` - Resistance mutations only
- `results/fq_resistance/qc_summary.tsv` - QC statistics
- `results/fq_resistance/mutation_catalog.tsv` - Unique mutations detected

---

### Phase 2: Species-Level Resistance Prevalence

**Objective**: Determine which species carry FQ resistance mutations and at what prevalence

**Script**: `scripts/python/11_fq_species_prevalence.py`

#### Analyses:
1. **Per-sample species resistance status**:
   - For each sample × species, calculate:
     - Number of resistance mutations
     - Mean resistance allele frequency
     - Presence of known high-impact mutations (e.g., E. coli gyrA S83L, parC S80I)

2. **Species prevalence by location**:
   - % samples with resistant E. coli (UCMC vs ZCH)
   - % samples with resistant K. pneumoniae
   - % samples with resistant S. aureus
   - etc.
   - **Test**: Fisher's exact test for prevalence differences

3. **Species prevalence by body site**:
   - Stratify by Axilla/Groin/Stool
   - **Test**: Chi-square test

#### Outputs:
- `results/fq_resistance/species_prevalence_by_location.tsv`
- `results/fq_resistance/species_prevalence_by_bodysite.tsv`
- `figures/fq_resistance/species_prevalence_barplot.pdf`

---

### Phase 3: Mutation-Level Analysis

**Objective**: Compare specific mutation frequencies across groups

**Script**: `scripts/python/12_fq_mutation_analysis.py`

#### Analyses:

##### 3.1 Top Mutations Overall
- Identify most prevalent resistance mutations (top 20-30)
- Focus on mutations detected in this dataset:
  - **S. aureus**: parE Y470N, parE D/S432V, parC E/R80S, gyrA G/S84L/Y
  - **K. pneumoniae**: gyrA K154R (nearly universal!), gyrA D/S83A/L/F/N, gyrA D/S87Y/N/A
  - **K. oxytoca**: parC S80R/I/T, gyrA D87A/G/Y, gyrA T83I
  - **E. faecium**: gyrA E/S87Y (fixed mutation)
  - **E. faecalis**: gyrA E/S87K (fixed mutation)
  - **E. coli**: parE D/S458A, parC A/G78V (low frequency, ~1-5%)
  - **P. aeruginosa**: gyrA D/S87Y/N, gyrA A/D83V
  - **S. mitis**: parC S79I, gyrA S81F
  - **S. oralis**: gyrA S81F/Y

##### 3.2 Location Comparison (UCMC vs ZCH)
For each mutation:
- Mean frequency at UCMC vs ZCH
- Prevalence (% samples with mutation)
- **Test**: Mann-Whitney U test (frequency), Fisher's exact (prevalence)
- **Multiple testing correction**: Benjamini-Hochberg FDR

**Stratify by body site** - critical!

##### 3.3 Body Site Comparison
- Compare mutation frequencies: Axilla vs Groin vs Stool
- **Hypothesis**: Stool has higher FQ resistance (gut microbiome reservoir)
- **Test**: Kruskal-Wallis test

##### 3.4 Longitudinal Changes (Week 1 → Week 3)
- For paired samples (complete subjects), test if mutation frequencies change
- **Test**: Wilcoxon signed-rank test (paired)
- Stratify by body site

#### Outputs:
- `results/fq_resistance/top_mutations.tsv`
- `results/fq_resistance/mutation_location_comparison.tsv`
- `results/fq_resistance/mutation_bodysite_comparison.tsv`
- `results/fq_resistance/mutation_longitudinal_changes.tsv`
- `figures/fq_resistance/mutation_heatmap_by_location.pdf`
- `figures/fq_resistance/mutation_frequency_boxplots.pdf`

---

### Phase 4: Gene-Level Resistance Burden

**Objective**: Summarize resistance across entire genes (gyrA, gyrB, parC, parE)

**Script**: `scripts/python/13_fq_gene_level_analysis.py`

#### Analyses:

##### 4.1 Composite Resistance Score
For each sample × species × gene:
- Sum all resistance mutation frequencies
- Example: If E. coli gyrA has S83L (freq=0.9) + D87N (freq=0.6) → score = 1.5

##### 4.2 Gene-Level Comparisons
- Compare composite scores: UCMC vs ZCH, by body site, Week 1 vs Week 3
- **Test**: Mann-Whitney U, Kruskal-Wallis, Wilcoxon signed-rank

##### 4.3 Multi-Gene Resistance Patterns
- Co-occurrence of mutations in gyrA + parC (double mutants)
- Triple/quadruple mutants (gyrA + gyrB + parC + parE)
- **Hypothesis**: Higher mutation burden → higher-level resistance

#### Outputs:
- `results/fq_resistance/gene_level_scores.tsv`
- `results/fq_resistance/multi_gene_patterns.tsv`
- `figures/fq_resistance/gene_scores_by_location.pdf`
- `figures/fq_resistance/mutation_burden_distribution.pdf`

---

### Phase 5: Antibiotic Exposure Correlations

**Objective**: Test if non-FQ antibiotics correlate with FQ resistance

**Script**: `scripts/python/14_fq_antibiotic_correlations.py`

#### Rationale:
Even though no fluoroquinolones were given, other antibiotics might:
- Select for multi-drug resistant strains (that also carry FQ resistance)
- Co-select via linked resistance genes

#### Analyses:

##### 5.1 Total Antibiotic Exposure vs FQ Resistance
- For each sample, correlate:
  - Total antibiotic days (any antibiotic)
  - Number of FQ resistance mutations
  - Mean FQ mutation frequency
- **Test**: Spearman correlation
- Stratify by body site

##### 5.2 Specific Antibiotic Classes
Test correlation with:
- Beta-lactams (ampicillin, piperacillin-tazobactam)
- Aminoglycosides (gentamicin)
- Glycopeptides (vancomycin)
- Carbapenems

##### 5.3 High vs Low Exposure Groups
- Stratify samples: No antibiotics vs Any antibiotics
- Compare FQ resistance burden
- **Test**: Mann-Whitney U test

#### Outputs:
- `results/fq_resistance/antibiotic_correlations.tsv`
- `results/fq_resistance/antibiotic_groups_comparison.tsv`
- `figures/fq_resistance/antibiotic_fq_scatterplots.pdf`

---

### Phase 6: Mixed-Effects Modeling

**Objective**: Comprehensive model accounting for repeated measures

**Script**: `scripts/python/15_fq_mixed_models.py`

#### Model Structure:
```python
# Per-sample FQ resistance burden
FQ_resistance_score ~ Location + BodySite + Week +
                       Location:BodySite + Location:Week + BodySite:Week +
                       (1|SubjectID)
```

#### Response Variables:
1. **Total FQ mutation count** - Number of resistance mutations per sample
2. **Mean FQ mutation frequency** - Average frequency across all mutations
3. **Species-specific scores** - E.g., E. coli gyrA resistance score

#### Fixed Effects to Test:
- **Location**: UCMC vs ZCH
- **BodySite**: Axilla vs Groin vs Stool
- **Week**: Week 1 vs Week 3
- **Interactions**: Location × BodySite, BodySite × Week

#### Random Effects:
- `(1|SubjectID)` - Account for repeated measures within subjects

#### Model Diagnostics:
- Residual plots
- Q-Q plots
- Intraclass correlation (ICC)

#### Outputs:
- `results/fq_resistance/mixed_model_coefficients.tsv`
- `results/fq_resistance/mixed_model_anova.tsv`
- `figures/fq_resistance/model_diagnostics.pdf`
- `figures/fq_resistance/predicted_values_by_factors.pdf`

---

### Phase 7: Publication Figures

**Objective**: Create publication-quality figures for manuscript

**Script**: `scripts/python/16_fq_publication_figures.py`

#### Figure 1: FQ Resistance Overview
**Panel A**: Heatmap of top 30 mutations × samples
- Rows: Mutations (grouped by species/gene)
- Columns: Samples (ordered by Location, BodySite, Week)
- Color: Mutation frequency (0-100%)
- Annotation: Location, BodySite, Week

**Panel B**: Species-level prevalence barplot
- X-axis: Species
- Y-axis: % samples with resistance mutations
- Fill: Location (UCMC vs ZCH)
- Facet: Body site

#### Figure 2: Location and Body Site Comparisons
**Panel A**: Boxplots of FQ mutation frequency by Location
- X-axis: Location (UCMC vs ZCH)
- Y-axis: Mean FQ mutation frequency per sample
- Facet: Body site
- Statistics: p-values from Mann-Whitney tests

**Panel B**: Boxplots by Body Site
- X-axis: Body site (Axilla, Groin, Stool)
- Y-axis: Mean FQ mutation frequency
- Fill: Location
- Statistics: p-values from Kruskal-Wallis

#### Figure 3: Longitudinal Trajectories
**Panel A**: Paired trajectories for complete subjects
- X-axis: Week 1 vs Week 3
- Y-axis: FQ resistance score
- Lines: Individual subjects
- Facet: Body site
- Color: Location

**Panel B**: Change in specific mutations (Week 1 → 3)
- Bar plot of log2 fold-change for top changing mutations
- Stratify by body site
- Highlight significant changes (FDR < 0.05)

#### Figure 4: Species-Specific Resistance Patterns
Focus on key pathogens with substantial resistance data: S. aureus, K. pneumoniae, K. oxytoca, Enterococcus spp., S. mitis:

**Panel A**: Mutation frequency heatmap (species-specific)
- Rows: Known resistance mutations for each species
- Columns: Samples with that species detected
- Color: Mutation frequency

**Panel B**: Multi-gene resistance patterns
- Venn diagram or UpSet plot showing co-occurrence of:
  - gyrA only
  - parC only
  - gyrA + parC
  - gyrA + parC + gyrB
  - etc.

#### Supplementary Figures:
- **S1**: Mutation catalog with clinical resistance classification
- **S2**: Antibiotic exposure vs FQ resistance scatterplots
- **S3**: Species-stratified longitudinal changes
- **S4**: Mixed-effects model predicted values

#### Outputs:
- `figures/fq_resistance/publication/Figure1_FQ_Overview.pdf`
- `figures/fq_resistance/publication/Figure2_Location_BodySite.pdf`
- `figures/fq_resistance/publication/Figure3_Longitudinal.pdf`
- `figures/fq_resistance/publication/Figure4_Species_Patterns.pdf`
- `figures/fq_resistance/supplementary/` - All supplementary figures

---

### Phase 8: Summary Report and Tables

**Objective**: Generate manuscript-ready tables and summary text

**Script**: `scripts/python/17_fq_summary_report.py`

#### Tables:

##### Table 1: FQ Resistance Mutation Catalog
Columns:
- Species
- Gene
- Position
- Mutation (e.g., S83L)
- Clinical significance (High/Medium/Low resistance)
- Prevalence (% samples)
- Mean frequency
- Location (UCMC vs ZCH prevalence)

##### Table 2: Species-Level Resistance Summary
Columns:
- Species
- N samples with species
- N samples with resistance mutations
- % with resistance
- Top mutations (most common)
- UCMC vs ZCH comparison (p-value)

##### Table 3: Differential Mutations by Location
For mutations with FDR < 0.05:
- Mutation
- UCMC mean frequency
- ZCH mean frequency
- Fold-change
- p-value
- FDR-adjusted p-value

##### Table 4: Longitudinal Changes
For mutations changing significantly Week 1 → 3:
- Mutation
- Body site
- Week 1 mean frequency
- Week 3 mean frequency
- % change
- p-value (paired Wilcoxon)
- FDR

##### Table 5: Mixed-Effects Model Results
- Factor
- Coefficient
- Standard error
- t-value
- p-value
- 95% CI

#### Summary Report:
- Markdown report with key findings
- Statistics for Results section
- Interpretation for Discussion

#### Outputs:
- `results/fq_resistance/summary/Table1_Mutation_Catalog.tsv`
- `results/fq_resistance/summary/Table2_Species_Summary.tsv`
- `results/fq_resistance/summary/Table3_Differential_Mutations.tsv`
- `results/fq_resistance/summary/Table4_Longitudinal_Changes.tsv`
- `results/fq_resistance/summary/Table5_Mixed_Model.tsv`
- `results/fq_resistance/summary/FQ_RESISTANCE_REPORT.md`

---

## Statistical Considerations

### Multiple Testing Correction
- **Problem**: Testing hundreds of mutations across multiple comparisons
- **Solution**: Benjamini-Hochberg FDR correction
- **Threshold**: FDR < 0.05 for significance

### Body Site Stratification
- **Critical**: Stool vs skin microbiomes differ dramatically
- **Approach**:
  - Stratify all analyses by body site OR
  - Include BodySite as covariate in models

### Paired vs Unpaired Tests
- **Complete subjects (n=53)**: Use paired tests (Wilcoxon signed-rank)
- **All samples**: Use unpaired tests (Mann-Whitney, Kruskal-Wallis)

### Normalization
- FQ mutation frequencies are already normalized (0-1 scale)
- For composite scores, may need log-transformation if highly skewed

### Missing Data
- Not all species detected in all samples
- Strategy: Analyze only samples where species is present (depth ≥ 10X)

---

## Key Hypotheses to Test

### H1: Location Effect
**Null**: FQ resistance mutation frequencies are equal between UCMC and ZCH
**Alternative**: Significant differences in resistance profiles
**Prediction**: Expect similarity (like AMR gene analysis)

### H2: Body Site Effect
**Null**: FQ resistance is equal across body sites
**Alternative**: Stool has higher FQ resistance than skin sites
**Prediction**: Stool > Groin/Axilla (gut reservoir hypothesis)

### H3: Temporal Dynamics
**Null**: FQ resistance does not change from Week 1 to Week 3
**Alternative**: Resistance increases over time
**Prediction**: Increase in groin/stool (matching AMR trend)

### H4: Antibiotic Selection
**Null**: Non-FQ antibiotic exposure does not correlate with FQ resistance
**Alternative**: Correlation exists (co-selection)
**Prediction**: Weak or no correlation (based on prior AMR findings)

---

## Expected Challenges and Solutions

### Challenge 1: Low Depth for Some Species
- **Problem**: Not all species detected in all samples
- **Solution**: Set minimum depth threshold (e.g., ≥10 reads per position)

### Challenge 2: Mutation Frequency vs Prevalence
- **Problem**: High frequency in few samples vs low frequency in many samples
- **Solution**: Report both metrics:
  - **Prevalence**: % samples with mutation
  - **Frequency**: Mean allele frequency when present

### Challenge 3: Clinical Interpretation
- **Problem**: Not all mutations have equal resistance impact
- **Solution**: Classify mutations:
  - **High impact**: Known to confer resistance (e.g., gyrA S83L)
  - **Low impact**: Polymorphisms with unclear significance
  - Use literature to categorize

### Challenge 4: Multiple Mutations per Gene
- **Problem**: Samples may have multiple mutations in same gene
- **Solution**: Create composite score (sum of frequencies)

---

## Timeline Estimate

Assuming ~2-3 hours per phase:

1. **Phase 1-2**: Data prep + species prevalence (~4 hours)
2. **Phase 3-4**: Mutation and gene-level analysis (~6 hours)
3. **Phase 5**: Antibiotic correlations (~2 hours)
4. **Phase 6**: Mixed models (~3 hours)
5. **Phase 7**: Publication figures (~4 hours)
6. **Phase 8**: Tables and report (~2 hours)

**Total**: ~20-25 hours of analysis time

---

## References for Clinical Mutations

Key papers to cite for FQ resistance mutations:

1. **Enterobacteriaceae (E. coli, Klebsiella spp.)**:
   - gyrA S83L, D87N; parC S80I (PMID: 8071223, 9322037)
   - K. pneumoniae gyrA K154R (PMID: 15583292)
   - K. oxytoca QRDR mutations (PMID: 15583292)

2. **S. aureus**:
   - gyrA S84L; parC S80F/Y, E80S (PMID: 8071224, 9322038)
   - parE mutations less well characterized

3. **Enterococcus spp.**:
   - gyrA S83I/Y, E87K/Y (PMID: 10103236)
   - Fixed resistance mutations common

4. **P. aeruginosa**:
   - gyrA T83I, D87N (PMID: 7721033)
   - parC S87L (PMID: 9758780)

5. **Streptococcus spp.**:
   - parC S79F/Y/I (PMID: 9620347)
   - gyrA S81F/Y (PMID: 10223938)
   - Different numbering than Gram-negatives

---

## Directory Structure

Create subdirectories for organization:

```
nicu_resistome_analysis/
├── results/
│   └── fq_resistance/
│       ├── qc/
│       ├── species_prevalence/
│       ├── mutation_analysis/
│       ├── gene_analysis/
│       ├── antibiotic_correlations/
│       ├── mixed_models/
│       └── summary/
└── figures/
    └── fq_resistance/
        ├── publication/
        └── supplementary/
```

---

## Next Steps

When you're ready to start:

1. Run Phase 1 (data preparation) first - this creates the filtered dataset
2. Review the filtered data to ensure it looks correct
3. Proceed through phases sequentially
4. Generate summary report after each phase to check results
5. Create figures last, once you know what the key findings are

**Good luck! Let me know if you have questions when you start implementing.**

---

**Created by**: Claude Code
**Date**: 2025-11-10
