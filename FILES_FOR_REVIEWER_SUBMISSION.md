# Files to Include in Reviewer Response

## Summary
You performed comprehensive ARG analysis in response to the reviewer's suggestion. While it won't be in the main manuscript (to keep focus), you're offering to include 2-3 supplementary figures showing you did thorough analysis.

---

## Recommended Approach

**Strategy**: Acknowledge excellent suggestion → Show comprehensive analysis was done → Explain surprising findings (no geographic differences!) → Offer supplementary figures → Explain why not adding to main text (focus/scope)

---

## Three Figures to Include as Supplementary Material

### Option A: Two Core Figures (Minimal Response)

**Supplementary Figure S#**: Resistome composition by location and body site
- **File to submit**: `/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis/figures/publication/Figure2_BodySite_Comparison.pdf`
- **Figure legend**: "Antibiotic resistance gene burden by location, body site, and week. Total AMR reads per million (RPM) shown for UCMC (Cincinnati) and ZCH (Hangzhou) samples across three body sites (axilla, groin, stool) at Week 1 and Week 3. No significant differences were observed between locations (Mann-Whitney U, all p > 0.05 after FDR correction). Stool samples showed significantly higher AMR burden than skin sites (p < 0.0001). Sample sizes: UCMC n=359, ZCH n=247."

**Supplementary Figure S#**: Longitudinal AMR trajectories by body site
- **File to submit**: `/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis/figures/publication/Figure3_Longitudinal_Trajectories.pdf`
- **Figure legend**: "Longitudinal changes in antibiotic resistance gene burden from Week 1 to Week 3. Paired trajectories shown for subjects with complete data at both timepoints (n=53 subjects, 40 UCMC, 13 ZCH). Groin samples showed dramatic AMR increase (+204%, p=0.0001, paired Wilcoxon), stool samples increased moderately (+64%, p<0.0001), while axilla samples remained stable (+0.6%, p=0.96). Each line represents a single subject."

### Option B: Add Third Figure (Stronger Response - RECOMMENDED)

**Supplementary Figure S#**: Antibiotic exposure versus resistome burden
- **File to submit**: `/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis/figures/exploratory/antibiotic_amr_correlation.pdf`
- **Figure legend**: "Correlation between cumulative antibiotic exposure and antibiotic resistance gene burden. Scatter plots show total AMR RPM versus cumulative antibiotic days, stratified by body site. No significant correlations were observed (Spearman: Axilla rho=-0.093, p=0.19; Groin rho=0.073, p=0.30; Stool rho=-0.077, p=0.28), suggesting resistance colonization is independent of direct antibiotic selection pressure during NICU hospitalization."

**Why include the third figure?** It directly addresses the reviewer's concern about "significant differences in antibiotic use patterns" by showing that despite different prescribing, resistance patterns are similar AND resistance doesn't correlate with exposure.

---

## Key Statistics to Mention in Response

### No Geographic Differences
- Differential abundance: **0/237 genes** significant (FDR < 0.05)
- Mixed-effects model location effect: **p = 0.45**
- Interpretation: Remarkably similar resistomes despite different institutions

### No Antibiotic Correlation
- Overall Spearman correlation: **rho = 0.001, p = 0.985**
- Axilla: rho = -0.093, p = 0.19
- Groin: rho = 0.073, p = 0.30
- Stool: rho = -0.077, p = 0.28
- Interpretation: Resistance colonization independent of antibiotic pressure

### Body Site is Primary Driver
- PCA PC1 variance: **30.8%** (separates stool from skin)
- Mean AMR RPM: Stool = 4,888; Groin = 2,553; Axilla = 694
- Mixed model coefficient: **Stool +0.615, p < 0.0001**
- Interpretation: **7-fold** stronger effect than any other factor

### Temporal Dynamics
- Groin Week 1→3: **+204% increase, p = 0.0001**
- Stool Week 1→3: **+64% increase, p < 0.0001**
- Axilla Week 1→3: +0.6% increase, p = 0.96 (not significant)
- Gene-level: 120 genes increased in groin, 104 in stool

### K. pneumoniae FQ Resistance (addresses ESBL example)
- gyrA K154R mutation: **99.5%** prevalence in K. pneumoniae samples
- No location effect: **p = 0.79**
- Interpretation: Resistance is pre-existing and universal

---

## Text to Add to Results Section (if including supplementary figures)

**Option 1 - Brief mention (1 sentence):**
> "Comprehensive analysis of antibiotic resistance genes revealed no significant differences between UCMC and ZCH resistomes (0/237 genes significant, FDR < 0.05), with body site and temporal dynamics as the primary drivers (Supplementary Figures S#-S#)."

**Option 2 - More detailed (2-3 sentences):**
> "Comprehensive analysis of antibiotic resistance genes across 237 resistance genes and 669 samples revealed no significant differences between UCMC and ZCH resistomes (0/237 genes significant, FDR < 0.05), despite documented differences in institutional antibiotic prescribing patterns. Body site was the primary driver of resistome composition (stool >> groin > axilla, p < 0.0001), with dramatic temporal increases observed in groin (+204%, p = 0.0001) and stool (+64%, p < 0.0001) from Week 1 to Week 3. No correlation was observed between antibiotic exposure and resistance burden (Spearman rho = 0.001, p = 0.985), suggesting resistance colonization is largely independent of direct antibiotic selection pressure during NICU hospitalization (Supplementary Figures S#-S#)."

---

## Materials Available (if reviewer wants more detail)

### Additional Figures Available:
- PCA plot showing body site separation: `figures/publication/Figure1_PCA.pdf`
- Volcano plots of gene-level changes: `figures/publication/Figure4_Volcano_Plots.pdf`
- FQ resistance figures (4 publication figures): `figures/fq_resistance/publication/`
- Model diagnostics: `figures/mixed_models/`

### Summary Reports:
- Full analysis report: `results/summary/ANALYSIS_REPORT.md`
- Key findings table: `results/summary/key_findings.tsv`
- Top changing genes: `results/summary/top_changing_genes.tsv`
- Mixed model coefficients: `results/mixed_models/full_model_coefficients.tsv`
- FQ resistance report: `results/fq_resistance/summary/FQ_RESISTANCE_REPORT.md`

### Methods Summary:
- Sample size: 669 samples, 127 subjects
- Genes analyzed: 237 (after QC filtering from 449)
- Statistical tests: Mann-Whitney U, Wilcoxon signed-rank, Spearman correlation, linear mixed-effects models
- Multiple testing correction: Benjamini-Hochberg FDR
- ICC: 0.262 (26.2% between-subject variance)

---

## Study Design Limitation to Clarify

**Important**: We cannot directly link microbiome resistome to BSI isolate susceptibility because:
1. Microbiome sequencing = subset of infants
2. BSI data = all infants in NICU during study period
3. Only small number of microbiome-sequenced infants had BSI
4. Therefore, no robust statistical comparison possible

This is addressed in the reviewer response document under "Addressing the Specific Example: ESBL-Related ARGs and K. pneumoniae"

---

## Response Strategy Summary

**What you did**: Comprehensive ARG analysis (237 genes, 669 samples, 9 analysis scripts, 40+ figures)

**What you found**:
- Surprising! No geographic differences despite different antibiotic use
- No antibiotic correlation (resistance appears pre-existing)
- Body site dominates (stool >> skin)
- Dramatic temporal increases in groin/stool

**Why not in manuscript**: Gets off-topic from core convergence narrative

**What you're offering**: 2-3 supplementary figures + brief Results mention

**Key message**: We did thorough analysis → Found surprising results → Willing to share → Keeps manuscript focused

---

## Recommended Files to Attach to Reviewer Response

1. The reviewer response text (from `REVIEWER_RESPONSE_AMR_Analysis.md`)
2. The 2-3 supplementary figures (PDFs listed above)
3. Optional: One supplementary table with key statistics (could use `results/summary/key_findings.tsv` or create a cleaner version)

All files are publication-ready and formatted appropriately for submission.
