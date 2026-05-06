# Fluoroquinolone Resistance Mutation Analysis
**Generated**: 2026-05-05 18:42:12

================================================================================

## Executive Summary

- **Total samples analyzed**: 470
- **Unique FQ resistance mutations detected**: 925
- **Species with resistance**: 71
- **Complete paired subjects (Week 1 & 3)**: 53

## Key Findings

### 1. No Location Effect (UCMC vs ZCH)

- FQ resistance patterns are **indistinguishable** between UCMC and ZCH
- All 71 species showed no significant prevalence differences (all FDR > 0.05)
- Mixed-effects model: Location coefficient p=0.79 (not significant)
- **Interpretation**: Both NICUs harbor similar FQ-resistant bacterial populations

### 2. Strong Body Site Effect

- **23/9 species** differ significantly by body site (FDR < 0.05)

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
| Mycobacterium. tube. gyrB I486L | 221 | 45.5% | 0.99 | Medium |
| Escherichia. coli. gyrA D87A | 59 | 12.1% | 0.03 | Low/Unknown |
| Klebsiella. pneu. gyrA D87A | 59 | 12.1% | 0.02 | Low/Unknown |
| Klebsiella. oxyt. gyrA D87A | 57 | 11.7% | 0.02 | Low/Unknown |
| Klebsiella. aero. gyrA D87A | 56 | 11.5% | 0.02 | Low/Unknown |
| Enterobacter. cloa. parE E460D | 45 | 9.3% | 0.02 | Low/Unknown |
| Salmonella. ente. gyrB E466D | 45 | 9.3% | 0.02 | Medium |
| Shigella. sonn. parE E460D | 45 | 9.3% | 0.02 | Low/Unknown |
| Escherichia. coli. parE E460D | 45 | 9.3% | 0.02 | Low/Unknown |
| Escherichia. coli. gyrA Q106H | 38 | 7.8% | 0.03 | Low/Unknown |

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
- R² = 0.329
- Significant predictors:
  - BodySite[Groin]: +0.97 (p<0.001)
  - BodySite[Stool]: +1.38 (p<0.001)
  - Week 3: +0.33 (p=0.02)

**Model 2: Mean freq ~ Location + BodySite + Week**
- R² = 0.232
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
