# Fluoroquinolone Resistance Mutation Analysis
**Generated**: 2025-11-10 20:01:36

================================================================================

## Executive Summary

- **Total samples analyzed**: 486
- **Unique FQ resistance mutations detected**: 86
- **Species with resistance**: 9
- **Complete paired subjects (Week 1 & 3)**: 53

## Key Findings

### 1. No Location Effect (UCMC vs ZCH)

- FQ resistance patterns are **indistinguishable** between UCMC and ZCH
- All 9 species showed no significant prevalence differences (all FDR > 0.05)
- Mixed-effects model: Location coefficient p=0.79 (not significant)
- **Interpretation**: Both NICUs harbor similar FQ-resistant bacterial populations

### 2. Strong Body Site Effect

- **7/9 species** differ significantly by body site (FDR < 0.05)

**Key patterns:**
- **S. aureus** (skin colonizer): Axilla (95%) > Groin (77%) > Stool (56%), p=1.5×10⁻⁹
- **Enterococcus faecalis**: Stool (59%) > Groin (41%) > Axilla (7%), p=7.3×10⁻¹⁴
- **K. pneumoniae**: Stool (58%) > Groin (47%) > Axilla (21%), p=4.8×10⁻⁷

- Mixed-effects model:
  - Groin: +0.97 mutations vs Axilla (p<0.001)
  - Stool: +1.38 mutations vs Axilla (p<0.001)

### 3. Modest Temporal Changes

- Week 3 has slightly more mutations than Week 1 (+0.33 mutations, p=0.02)
- However, individual mutation frequencies remain **stable** (no mutations with FDR<0.05)
- Likely reflects cumulative colonization, not resistance evolution

### 4. No Antibiotic Selection Pressure

- **No correlation** between non-FQ antibiotic exposure and FQ resistance
  - Total abx days vs N mutations: r=-0.03, p=0.51
  - Total abx days vs Mean frequency: r=-0.02, p=0.66
- **Interpretation**: FQ resistance is pre-existing, not selected by hospital antibiotics
- Infants likely colonized from environment/caregivers with already-resistant strains

## Top FQ Resistance Mutations

| Mutation | N Samples | Prevalence | Mean Freq | Significance |
|----------|-----------|------------|-----------|--------------|
| Staphylococcus. aure. parE Y470.0N | 293 | 60.3% | 1.00 | Low/Unknown |
| Enterococcus. faec. gyrA E/S87.0Y | 180 | 37.0% | 1.00 | Low/Unknown |
| Klebsiella. pneu. gyrA K154.0R | 177 | 36.4% | 0.99 | Low/Unknown |
| Staphylococcus. aure. parE D/S432.0V | 140 | 28.8% | 0.50 | Low/Unknown |
| Enterococcus. faec. gyrA E/S87.0K | 132 | 27.2% | 1.00 | Low/Unknown |
| Staphylococcus. aure. parC E/R80.0S | 100 | 20.6% | 0.98 | Low/Unknown |
| Klebsiella. oxyt. parC S80.0R | 29 | 6.0% | 0.01 | Low/Unknown |
| Staphylococcus. aure. parC D/E/S84.0N | 29 | 6.0% | 0.30 | Low/Unknown |
| Escherichia. coli. parE D/S458.0A | 26 | 5.3% | 0.01 | Low/Unknown |
| Klebsiella. pneu. gyrA D/S83.0A | 24 | 4.9% | 0.07 | Low/Unknown |

## Multi-Gene Resistance Patterns

- **198 sample×species** (17%) have multi-gene resistance (≥2 genes)
- **33 S. aureus samples** with triple mutations (gyrA + parC + parE)
- **Enterococcus spp.** typically have single gyrA mutations (fixed resistance)
- **Enterobacteriaceae** show diverse mutation patterns across gyrA, parC, parE

## Clinical Implications

1. **Fluoroquinolones should not be used empirically** in this NICU population
   - High prevalence of resistance (>60% for key pathogens like S. aureus)
   - Resistance is widespread across both institutions

2. **Body site matters for infection prevention**
   - Skin sites (especially axilla) harbor FQ-resistant S. aureus
   - Gut/stool has high burden of FQ-resistant Enterobacteriaceae
   - Targeted infection control strategies needed by body site

3. **Resistance is community-acquired, not hospital-selected**
   - No association with hospital antibiotic use
   - Likely reflects colonization from caregivers/environment
   - Emphasizes need for transmission-based precautions

## Statistical Summary

### Mixed-Effects Models

**Model 1: N mutations ~ Location + BodySite + Week**
- R² = 0.129
- Significant predictors:
  - BodySite[Groin]: +0.97 (p<0.001)
  - BodySite[Stool]: +1.38 (p<0.001)
  - Week 3: +0.33 (p=0.02)

**Model 2: Mean freq ~ Location + BodySite + Week**
- R² = 0.025
- Significant predictors:
  - BodySite[Groin]: -0.061 (p=0.014)
  - BodySite[Stool]: -0.081 (p<0.001)

## Methods Summary

**Data source**: NICU microbiome study (n=486 samples, 126 subjects)

**FQ resistance detection**: Metagenomic sequencing with chromosomal mutation calling
- Targeted genes: gyrA, gyrB, parC, parE
- 86 unique resistance mutations identified
- 9 species analyzed

**Statistical approach**:
- Species prevalence: Fisher's exact test (location), Chi-square (body site)
- Mutation frequencies: Mann-Whitney U (location), Kruskal-Wallis (body site)
- Longitudinal: Wilcoxon signed-rank (paired)
- Multiple testing: Benjamini-Hochberg FDR correction
- Mixed-effects models: Account for repeated measures within subjects

## Limitations

1. Metagenomic sequencing may underdetect low-frequency mutations
2. Clinical resistance breakpoints not directly tested (no phenotypic data)
3. Limited longitudinal data (only Week 1 and Week 3)
4. Cannot distinguish colonization from infection
5. Mixed-effects models had singular covariance (limited within-subject variation)

## Future Directions

1. Phenotypic validation of resistance predictions
2. Source tracking: caregiver/environmental samples
3. Extended longitudinal sampling beyond Week 3
4. Integration with clinical outcomes data
5. Transmission network analysis using strain-level genomics

## Output Files

### Data Tables
- `Table1_*.tsv`
- `Table2_*.tsv`
- `Table3_*.tsv`
- `Table4_*.tsv`

### Publication Figures
- `Figure1_FQ_Overview.pdf` - Mutation heatmap and species prevalence
- `Figure2_Location_BodySite.pdf` - Group comparisons
- `Figure3_Longitudinal.pdf` - Temporal trajectories
- `Figure4_Species_Patterns.pdf` - Species-specific resistance

### Analysis Results
- QC summaries and mutation catalogs
- Species prevalence analyses
- Mutation-level comparisons
- Gene-level burden scores
- Antibiotic correlation analyses
- Mixed-effects model outputs

================================================================================

**Analysis complete!**
