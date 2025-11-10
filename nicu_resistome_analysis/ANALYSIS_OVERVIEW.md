# NICU Resistome Analysis Overview

## Project Structure

```
nicu_resistome_analysis/
├── data/                      # Symlinks to merged AMR/stress data
├── scripts/
│   ├── python/               # Python analysis scripts
│   └── R/                    # R statistical scripts
├── results/
│   ├── qc/                   # Quality control outputs
│   ├── exploratory/          # EDA results
│   ├── differential/         # Differential abundance
│   ├── longitudinal/         # Longitudinal analysis
│   ├── correlations/         # Antibiotic-AMR correlations
│   └── mixed_models/         # Mixed-effects models
└── figures/
    ├── exploratory/
    ├── differential/
    ├── longitudinal/
    └── publication/
```

## Experimental Design

### Study Structure
- **Locations**: Cincinnati (UCMC) vs Hangzhou (ZCH)
- **Body Sites**: Axilla, Groin, Stool (3 sites)
- **Time Points**: Week 1 vs Week 3 (2 timepoints)
- **Expected samples per subject**: 6 (3 sites × 2 weeks)

### Current Data Summary
- **Total samples**: 771
- **Samples with AMR data**: 606 (78.6%)
- **Unique subjects**: 68
  - Cincinnati: 386 samples
  - Hangzhou: 385 samples
- **Complete subjects** (all 6 samples): 9
- **Incomplete subjects**: 59

### Sample Distribution
| Body Site | Count |
|-----------|-------|
| Axilla    | 221   |
| Groin     | 211   |
| Stool     | 237   |
| Nonese    | 102   | *(failed/missing samples)*

| Week   | Count |
|--------|-------|
| Week 1 | 400   |
| Week 3 | 371   |

### Antibiotic Exposure
- **Subjects with antibiotic data**: 68 (100%)
- **Any antibiotics at Week 1**: 11/68 (16%)
- **Any antibiotics at Week 3**: 4/68 (6%)
- **Mean exposure days**:
  - Week 1: 0.69 days
  - Week 3: 0.40 days

**Note**: Low overall antibiotic exposure in this cohort.

## Data Files

### Merged Data (symlinked from /fastpool/analysis/nicu_amr_stress/merged/)
1. **nicu_amr_stress_merged.tsv** - Combined AMR + stress genes (77,219 rows)
2. **nicu_amr_only.tsv** - AMR genes only (76,958 rows)
3. **nicu_stress_only.tsv** - Stress response genes (262 rows)
4. **sample_summary.tsv** - Per-sample statistics (642 samples)
5. **gene_prevalence.tsv** - Gene prevalence across samples
6. **fq_mutation_matrix.tsv** - Fluoroquinolone resistance mutations
7. **fq_resistance_summary.tsv** - FQ resistance summary

### Metadata Files
1. **master_metadata.tsv** - Comprehensive metadata with:
   - Sample information (PatientID, Location, BodySite, Week)
   - AMR statistics (total_amr_rpm, unique_amr_genes, num_drug_classes)
   - Antibiotic exposure (all antibiotics by week)
   - Quality flags (has_metadata, has_amr_data, subject_complete)

2. **sample_qc_report.tsv** - QC status for each sample

## Analysis Plan

### Phase 1: Quality Control ✅ SETUP COMPLETE
1. ✅ Create master metadata file
2. ⏳ Assess sequencing depth and AMR burden
3. ⏳ Filter low-abundance genes
4. ⏳ Identify complete subjects for paired analysis
5. ⏳ Generate QC figures

### Phase 2: Exploratory Data Analysis
1. PCA/PCoA on AMR gene matrix
   - Stratify by body site (critical!)
   - Color by Location, shape by Week
2. Hierarchical clustering with metadata annotation
3. Diversity metrics (Shannon, Simpson) by Location × Site × Week
4. Overall AMR burden comparisons

### Phase 3: Differential Abundance Analysis
**Body site-specific analyses** (essential due to different microbiomes):

For each body site (Axilla, Groin, Stool):
1. Cincinnati vs Hangzhou comparisons
   - Use DESeq2 or edgeR for robust testing
   - Account for paired samples where possible
2. Identify location-specific resistance signatures
3. Drug class enrichment analysis
4. Create volcano plots and heatmaps

### Phase 4: Longitudinal Analysis
**Week 1 → Week 3 changes within subjects**:
1. Paired comparisons within complete subjects (n=9)
2. Test Location × Time interaction
3. Identify genes that change over time
4. Trajectory plots for key resistance genes

### Phase 5: Antibiotic-AMR Correlations
1. Correlate antibiotic exposure with AMR gene abundance
2. Stratify by antibiotic class
3. Test if exposure predicts resistance gene load
4. Account for baseline differences by location

### Phase 6: Mixed-Effects Modeling
**Comprehensive model** accounting for repeated measures:
```
AMR_abundance ~ Location + BodySite + Week +
                Location:BodySite + Location:Week + BodySite:Week +
                (1|Subject)
```
- Test main effects and interactions
- Separate models for total AMR burden and individual genes
- Multiple testing correction (FDR)

### Phase 7: Stress Response Analysis
1. Compare stress gene categories (SOS, oxidative, heat shock)
2. Test stress-AMR co-occurrence
3. Location-specific stress responses

### Phase 8: Fluoroquinolone Resistance
1. Mutation frequency comparisons
2. Location-specific mutations
3. Association with FQ antibiotic exposure

## Key Statistical Considerations

### Normalization
- **Using RPM** (Reads Per Million) - already calculated
- RPM accounts for sequencing depth differences
- Appropriate for between-sample comparisons

### Important Caveats
1. **Body site heterogeneity**: Stool ≠ skin microbiomes
   - Must analyze separately or control for site
2. **Paired structure**: Multiple samples per subject
   - Requires mixed models or blocking
3. **Incomplete data**: Only 9/68 subjects have all 6 samples
   - Limits power for full factorial repeated-measures analysis
   - Can still do cross-sectional comparisons
4. **Low antibiotic exposure**: Limited power for exposure-resistance correlations
5. **Multiple testing burden**: Hundreds of genes × multiple comparisons
   - Strict FDR correction required

### Recommended Approaches
1. **For complete subjects (n=9)**: Full repeated-measures ANOVA or mixed models
2. **For all samples (n=606 with AMR data)**: Cross-sectional comparisons with site stratification
3. **Body site**: Always analyze separately or include as covariate
4. **Location**: Primary factor of interest
5. **Week**: Secondary factor, test for temporal changes

## Next Steps

Run quality control script to:
1. Generate QC figures
2. Filter low-abundance genes
3. Create analysis-ready datasets

Then proceed with exploratory analysis and hypothesis testing.

## Scripts Created

1. `scripts/python/01_create_master_metadata.py` - Merge all metadata ✅
2. `scripts/python/02_quality_control.py` - QC analysis (ready to run)
3. More scripts to be added for each analysis phase

---

**Created**: 2025-11-10
**Status**: Data integration complete, ready for QC and analysis
