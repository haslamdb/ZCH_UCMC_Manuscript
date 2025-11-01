# Analysis Plan: Antibiotic Exposure and BSI Pathogen Abundance

## Objective
Correlate antibiotic exposures (week 1 and cumulative through week 3) with abundance of 8 significant BSI-associated pathogens.

---

## Data Preparation

### 1. Antibiotic Data
**Files:**
- `UC_antibiotics_parsed.csv`
- `ZCH_antibiotics_parsed.csv`

**Processing Steps:**
1. Load both antibiotic datasets
2. Combine UCMC and ZCH data
3. For each antibiotic:
   - Keep week 1 exposures as-is (`_w1` columns)
   - Create cumulative week 3: `cumulative_w3 = week1 + week3` (sum `_w1` + `_w2` columns)
   - Note: Our `_w2` columns actually represent week 3 data
4. Create final antibiotic matrix with:
   - Patient identifier (MRN)
   - Location (UCMC/ZCH)
   - All week 1 antibiotics (e.g., `Ampicillin_wk1`)
   - All cumulative week 3 antibiotics (e.g., `Ampicillin_cumulative_wk3`)

### 2. BSI Pathogen Data
**Need to identify:**
- What file contains the microbiome/pathogen abundance data?
- Likely from the R workspace: `/home/david/projects/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929/NICUData20250514`

**8 BSI Pathogens (from Figure 4):**
1. Staphylococcus aureus
2. Klebsiella pneumoniae
3. Klebsiella oxytoca
4. Enterococcus faecium
5. Enterococcus faecalis
6. Serratia marcescens
7. Escherichia coli
8. Streptococcus pyogenes

**Processing Steps:**
1. Extract abundance data for these 8 species
2. Filter to same patients in antibiotic dataset
3. Match by patient identifier (likely need MRN → Subject ID mapping)
4. Filter to appropriate time points (week 1 and week 3 samples)

---

## Linking Datasets

### Challenge: MRN vs Subject ID
- Antibiotic data uses **MRN** (e.g., 6402984, 65000272)
- Microbiome data likely uses **Subject ID** (e.g., N01, N02)

**Solutions:**
1. Check if Subject-MRN mapping exists in metadata files
2. Check `metadata/AllNICUSampleKey20250206.csv` for Subject column
3. May need to create mapping based on location and enrollment order

---

## Statistical Analysis Approaches

### Option 1: Spearman Correlation (Simple, Interpretable)
**For each antibiotic-pathogen pair:**
- Calculate Spearman's rank correlation coefficient
- Separate analysis for:
  - Week 1 antibiotics vs Week 1 pathogen abundance
  - Cumulative week 3 antibiotics vs Week 3 pathogen abundance
- Adjust for multiple comparisons (FDR correction)

**Pros:**
- Simple to interpret
- Non-parametric (handles zeros and skewed distributions)
- Direct measure of association

**Cons:**
- Doesn't account for confounders
- Doesn't separate location effects

---

### Option 2: Linear Mixed Models (More Rigorous)
**Model structure:**
```
Pathogen_Abundance ~ Antibiotic_Days + Location + GestationalAge + (1|Subject)
```

**For each pathogen:**
- Test association with each antibiotic separately
- Include random effect for subject (repeated measures)
- Control for location and gestational age
- Log-transform abundances (or use GLMM with appropriate family)

**Pros:**
- Controls for confounders
- Accounts for location differences
- Handles repeated measures

**Cons:**
- More complex
- May have convergence issues with sparse data

---

### Option 3: Elastic Net Regression (Exploratory)
**For each pathogen:**
- Use all antibiotics as predictors
- Elastic net penalty handles collinearity
- Identify which antibiotics are most predictive
- Separate models for week 1 and week 3

**Pros:**
- Handles multiple correlated predictors
- Built-in feature selection
- Can identify antibiotic combinations

**Cons:**
- Less interpretable
- Requires careful cross-validation

---

## Recommended Approach (Tiered)

### Tier 1: Exploratory Correlation Matrix
1. Create heatmap showing Spearman correlations
2. Rows: Antibiotics (week 1 and cumulative week 3)
3. Columns: 8 BSI pathogens (week 1 and week 3)
4. Color by correlation strength and significance
5. This gives quick overview of associations

### Tier 2: Individual Associations
For significant correlations (p < 0.05 after FDR):
1. Create scatter plots
2. Separate by location (UCMC vs ZCH)
3. Report correlation coefficients with 95% CI

### Tier 3: Adjusted Models (if needed)
For key findings from Tier 1-2:
1. Fit linear mixed models
2. Control for location, gestational age, sample type
3. Report adjusted effect sizes

---

## Visualization Plan

### Figure 1: Correlation Heatmap
- **Week 1 Panel**: Antibiotics (wk1) × Pathogens (wk1)
- **Week 3 Panel**: Antibiotics (cumulative wk3) × Pathogens (wk3)
- Color scale: Blue (negative) to Red (positive)
- Asterisks for significance (*, **, ***)

### Figure 2: Scatter Plots (Top Associations)
- Select top 6-8 significant associations
- Each panel: Antibiotic days (x) vs Pathogen abundance (y)
- Points colored by location (UCMC/ZCH)
- Add trend line with 95% CI
- Show Spearman ρ and p-value

### Figure 3: Network Diagram (Optional)
- Nodes: Antibiotics and Pathogens
- Edges: Significant correlations
- Edge thickness: Correlation strength
- Edge color: Positive (red) vs Negative (blue)

---

## Statistical Considerations

### Multiple Testing Correction
- Total tests: ~16 antibiotics × 8 pathogens × 2 time points = 256 tests
- Use FDR (Benjamini-Hochberg) correction
- Report both raw and adjusted p-values

### Handling Zeros
- Many patients have 0 antibiotic exposure
- Many samples have 0 pathogen abundance
- Options:
  1. Include zeros (Spearman handles this well)
  2. Analyze presence/absence separately
  3. Use zero-inflated models for adjusted analysis

### Location Stratification
- Consider separate analyses by location
- ZCH and UCMC have very different antibiotic practices
- May pool if effects are similar, stratify if heterogeneous

---

## Output Files

### Data Files
1. `antibiotic_pathogen_merged.csv` - Combined dataset
2. `antibiotic_pathogen_correlations.csv` - All correlation results
3. `antibiotic_pathogen_significant.csv` - FDR-significant associations

### Figures
1. `Antibiotic_Pathogen_Correlation_Heatmap.pdf`
2. `Antibiotic_Pathogen_Scatterplots.pdf`
3. `Antibiotic_Pathogen_Network.pdf` (optional)

### Statistics Tables
1. `Antibiotic_Pathogen_Stats_Week1.csv`
2. `Antibiotic_Pathogen_Stats_Week3.csv`

---

## Code Structure (To Implement)

### Script 1: Data Preparation
```r
# Load and merge antibiotic data
# Create cumulative week 3 exposures
# Load pathogen abundance data
# Link datasets by patient ID
# Output: merged dataset
```

### Script 2: Correlation Analysis
```r
# Calculate Spearman correlations
# Apply FDR correction
# Generate correlation matrices
# Output: correlation results
```

### Script 3: Visualization
```r
# Create correlation heatmap
# Generate scatter plots for top associations
# Optional: network diagram
# Output: publication-ready figures
```

---

## Questions to Resolve

1. **Data Location:**
   - Where is the pathogen abundance data?
   - Is it in the R workspace or separate CSV?

2. **Sample Matching:**
   - How do we link MRN to Subject ID?
   - Which samples to use (Stool only? All sample types?)

3. **Time Point Matching:**
   - Exact definition of "week 1" vs "week 3" samples
   - Should we use average abundance across sample types?

4. **Analysis Scope:**
   - Just Spearman correlations, or include adjusted models?
   - Analyze locations separately or combined?
   - Include sample type as variable?

---

## Expected Results

### Hypotheses:
1. **Negative associations**: Antibiotics that reduce specific pathogens
   - Example: Beta-lactams → reduced E. coli, Klebsiella

2. **Positive associations**: Antibiotics that may promote pathogen overgrowth
   - Example: Broad-spectrum → increased resistant organisms

3. **Cumulative effects**: Stronger associations at week 3 (cumulative) than week 1

4. **Location differences**: Different patterns at UCMC vs ZCH due to different antibiotic practices

---

## Timeline
1. **Data preparation** (30-60 min): Merge datasets, create cumulative exposures
2. **Correlation analysis** (15-30 min): Calculate correlations, FDR correction
3. **Visualization** (30-60 min): Generate heatmaps and scatter plots
4. **Review and interpretation** (30-60 min): Identify key findings

**Total: ~2-3 hours**

---

## Notes
- Start with simple Spearman correlations to get quick overview
- If interesting patterns emerge, can follow up with more sophisticated models
- Be prepared for sparse data (many zeros)
- Location differences may dominate signal
- Consider log-transforming abundances for visualization

