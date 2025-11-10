# FQ Resistance Analysis - Quick Start Guide

## What I've Created for You

### 1. Comprehensive Analysis Plan
**File**: `FQ_RESISTANCE_ANALYSIS_PLAN.md`

This is your complete roadmap with:
- 8 analysis phases covering all your questions
- Statistical methods for each comparison
- Expected outputs and figures
- Estimated timeline (~20-25 hours total)
- Key hypotheses to test
- Solutions to anticipated challenges

### 2. Directory Structure
```
results/fq_resistance/
├── qc/                         # Quality control outputs
├── species_prevalence/         # Species-level resistance
├── mutation_analysis/          # Mutation-level comparisons
├── gene_analysis/              # Gene-level resistance burden
├── antibiotic_correlations/    # Antibiotic exposure correlations
├── mixed_models/               # Mixed-effects model results
└── summary/                    # Tables and final report

figures/fq_resistance/
├── publication/                # Main manuscript figures
└── supplementary/              # Supplementary figures
```

### 3. Starter Script
**File**: `scripts/python/10_fq_data_preparation.py`

Ready-to-run Phase 1 script that will:
- Load FQ allele frequency data (14,009 rows)
- Merge with sample metadata
- Filter to resistance mutations only
- Create derived variables (mutation_id, high frequency indicators)
- Generate QC summaries and mutation catalog

**To run**:
```bash
cd /home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis
python scripts/python/10_fq_data_preparation.py
```

## What the Analysis Will Answer

### Primary Questions:
1. **Does FQ resistance differ by location?** (UCMC vs ZCH)
   - Compare mutation frequencies and prevalence
   - Stratified by body site and species

2. **Does FQ resistance differ by body site?** (Axilla vs Groin vs Stool)
   - Hypothesis: Stool has higher resistance (gut reservoir)
   - Test with Kruskal-Wallis

3. **Does FQ resistance change over time?** (Week 1 → Week 3)
   - Paired analysis in complete subjects
   - Expect increases in groin/stool (matching AMR trends)

4. **Does antibiotic exposure correlate with FQ resistance?**
   - Even though no fluoroquinolones given
   - Test for co-selection with other antibiotics
   - Likely null result (based on AMR findings)

## Analysis Workflow

### Phase 1: Data Prep ✓ (Script ready!)
Run `10_fq_data_preparation.py` to create filtered dataset

### Phase 2: Species Prevalence
- Which species carry FQ resistance?
- Prevalence by location and body site

### Phase 3: Mutation-Level Analysis
- Compare specific mutations across groups
- Top clinically relevant mutations
- Differential mutations (UCMC vs ZCH)

### Phase 4: Gene-Level Resistance
- Composite resistance scores
- Multi-gene resistance patterns
- Double/triple/quadruple mutants

### Phase 5: Antibiotic Correlations
- Total antibiotic days vs FQ resistance
- Specific antibiotic classes
- High vs low exposure groups

### Phase 6: Mixed-Effects Models
- Comprehensive model with repeated measures
- Account for subject clustering
- Test all interactions

### Phase 7: Publication Figures
- Figure 1: FQ Resistance Overview (heatmap + prevalence)
- Figure 2: Location/Body Site Comparisons (boxplots)
- Figure 3: Longitudinal Trajectories (paired lines)
- Figure 4: Species-Specific Patterns (heatmap + co-occurrence)

### Phase 8: Summary Report
- Manuscript-ready tables
- Key statistics for Results section
- Summary report with interpretation

## Species with FQ Resistance Data

**Priority species (substantial resistance data):**

1. **Staphylococcus aureus** - 339 samples with resistance mutations
   - parE Y470N (99.8% freq), parE D/S432V, parC E/R80S, gyrA G/S84L/Y
   - High frequency = established resistance

2. **Klebsiella pneumoniae** - 214 samples with resistance mutations
   - gyrA K154R nearly universal (99.5% in 177 samples!)
   - gyrA D/S83A/L/F/N, gyrA D/S87Y/N/A

3. **Enterococcus faecium** - 180 samples with resistance
   - gyrA E/S87Y (100% prevalence, 99.8% freq) - FIXED MUTATION

4. **Enterococcus faecalis** - 132 samples with resistance
   - gyrA E/S87K (100% prevalence, 99.7% freq) - FIXED MUTATION

5. **Escherichia coli** - 73 samples (but LOW frequency ~1-5%)
   - parE D/S458A, parC A/G78V - mostly wildtype, rare resistant alleles
   - Limited clinical significance

6. **Klebsiella oxytoca** - 62 samples with resistance
   - parC S80R/I/T, gyrA D87A/G/Y, gyrA T83I

7. **Pseudomonas aeruginosa** - 35 samples with resistance
   - gyrA D/S87Y/N, gyrA A/D83V

8. **Streptococcus mitis** - 33 samples with resistance
   - parC S79I (95.5% freq), gyrA S81F

9. **Streptococcus oralis** - 17 samples with resistance
   - gyrA S81F/Y

**Note**: Serratia marcescens and Streptococcus pyogenes do NOT have FQ resistance data in this dataset (though they were significant in pathogen analysis). S. pneumoniae detected but no resistance mutations observed.

## Data Files Available

### Input Data:
- `data/nicu_fq_allele_frequencies.csv` - 14,009 allele observations
- `data/fq_mutation_prevalence.tsv` - 281 mutation summary
- `data/fq_resistance_by_species.tsv` - 5,999 species-level classifications
- `results/qc/master_metadata_with_qc.tsv` - Sample metadata with antibiotics

### Key Columns in Allele Data:
- `sample_name`, `species`, `gene`, `position`
- `has_resistance_mutation` - True/False (already classified!)
- `total_resistant_frequency` - Frequency of resistance allele (0-1)
- `dominant_mutant_frequency` - Frequency of dominant mutant
- `total_depth` - Sequencing depth
- `mutation_summary` - Human-readable mutation (e.g., "S83L(95%)")

## Statistical Approach

### Comparisons:
- **Location**: Mann-Whitney U test (UCMC vs ZCH)
- **Body site**: Kruskal-Wallis test (Axilla vs Groin vs Stool)
- **Longitudinal**: Wilcoxon signed-rank (paired, Week 1 vs 3)
- **Prevalence**: Fisher's exact test or Chi-square

### Multiple Testing:
- Benjamini-Hochberg FDR correction
- Threshold: FDR < 0.05

### Mixed Models:
```
FQ_resistance ~ Location + BodySite + Week +
                Location:BodySite + Location:Week + BodySite:Week +
                (1|SubjectID)
```

## Expected Findings (Predictions)

Based on your AMR gene analysis:

1. **No location effect** - UCMC ≈ ZCH (like AMR genes)
2. **Strong body site effect** - Stool >> Groin > Axilla
3. **Temporal increase** - Week 3 > Week 1 (groin/stool)
4. **No antibiotic correlation** - Null result (like AMR)

## Timeline Estimate

- **Phase 1-2**: 4 hours (data prep + species prevalence)
- **Phase 3-4**: 6 hours (mutation + gene analysis)
- **Phase 5**: 2 hours (antibiotic correlations)
- **Phase 6**: 3 hours (mixed models)
- **Phase 7**: 4 hours (publication figures)
- **Phase 8**: 2 hours (tables + report)

**Total**: ~20-25 hours

Can be broken into shorter sessions!

## Tips for Implementation

1. **Start with Phase 1** - Run the data prep script first to see if filtering works correctly
2. **Check outputs** - Review the mutation catalog and QC summary before proceeding
3. **One phase at a time** - Each phase builds on the previous
4. **Save checkpoints** - Keep intermediate files so you can go back if needed
5. **Generate summaries** - Create a small report after each phase
6. **Body site stratification** - This is critical! Always stratify or control for body site

## Questions to Consider

As you analyze the data:

1. Are there specific mutations that are UCMC-specific or ZCH-specific?
2. Do multi-gene resistance patterns differ by location?
3. Which species drive the body site differences?
4. Are there longitudinal trajectories specific to certain pathogens?
5. Do beta-lactams correlate with FQ resistance (co-selection hypothesis)?

## Next Steps

When you get home:

1. Read through `FQ_RESISTANCE_ANALYSIS_PLAN.md` in detail
2. Run `scripts/python/10_fq_data_preparation.py`
3. Review the QC outputs to understand the data
4. Decide which analyses are most important for your manuscript
5. Start implementing phases sequentially

Good luck! This should give you a comprehensive FQ resistance analysis to complement your AMR gene findings.

---

**Created**: 2025-11-10
**Status**: Ready to implement
