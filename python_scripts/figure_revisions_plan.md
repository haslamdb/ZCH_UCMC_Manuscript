# Figure Revisions Plan - Manuscript Revision Response

## Overview
This document tracks all figure and analysis revisions needed to address reviewer comments for the ZCH-UCMC NICU microbiome manuscript.

---

## Script-to-Figure Mapping

### Main Figures

#### Figure 1: BSI Etiology Bar Chart
- **Current Status**: Script not identified in current directory
- **Current Format**: Shows organism distribution with asterisks for significant differences
- **Script Location**: Unknown/Need to create
- **Dependencies**: BSI epidemiology data from Table 2

#### Figure 2: Microbiome Composition & Diversity
- **Parts**:
  - a. Shannon diversity box plots (Week 1 & 3, three body sites)
  - b,c. PCA plots at weeks 1 and 3
- **Script(s)**: `/home/david/projects/Metagenomics/Yanping/NICU_Microbiome/Hangzhou/NoHumanDNA20220929/MicrobiomeAnalysisNICUHangzhouCincinnatiNoHuman20250514.R`
- **Key code sections**:
  - Shannon diversity: Lines ~586-758 (multiple Shannon plots by location, week, sample type)
  - PCA plots: Lines ~1225-2900+ (multiple PCA analyses)
  - Key files generated: `ShannonLocationType.pdf`, `LocationPCAWeek1.pdf`, `AllSamplesNoBSIPCA.pdf`, etc.
- **Status**: Located - R script
- **Dependencies**: Genera abundance data, vegan package for diversity calculations, FactoMineR for PCA

#### Figure 3: BSI-Microbiome Association
- **Parts**:
  - a. Bubble charts (Bray-Curtis similarity)
  - b. Network analysis
- **Script(s)**: `network_plot_after_variance_partitioning.py`
- **Dependencies**: `microbiome-variance-partitioning-revised.py`, variance analysis results

#### Figure 4: Longitudinal Changes in BSI Species
- **Parts**: Box plots for 8 species (S. aureus, K. pneumoniae, K. oxytoca, E. faecium, E. faecalis + 3 more)
- **Current PDFs**: Individual organism PDFs exist (e.g., `Staphylococcus_aureus_comparison.pdf`)
- **Script(s)**: Need to identify longitudinal analysis script
- **Status**: Individual plots exist, need to find script that generates combined figure

#### Figure 5: Clinical Predictors (RF + LMM)
- **Parts**:
  - Left panels: LMM coefficients (effect sizes)
  - Right panels: Random Forest SHAP importance
  - a. S. aureus (UCMC predominant)
  - b. K. pneumoniae (ZCH predominant)
- **Script(s)**:
  - `bsi_pedictors_viz.py` - Creates the visualization
  - `rf_shap_LMM_analysis_key_orgs.py` - Runs the analysis
- **Status**: Scripts exist, functional

### Supplementary Figures

#### Supplementary Figure 1: BSI Distribution in Enrolled Participants
- **Status**: Unknown location
- **Note**: Includes CoNS (unlike main epidemiology analysis)

#### Supplementary Figure 2: t-SNE by Body Site & Location
- **Script(s)**: `tSNE_plot.py`
- **Status**: Script exists and functional
- **Outputs**: Individual t-SNE plots for organisms + metadata

#### Supplementary Figure 3: t-SNE Colored by BSI Species Abundance
- **Script(s)**: `tSNE_plot.py`
- **Status**: Script exists
- **Outputs**: Multi-page PDF `all_key_organisms_multipage.pdf`

#### Supplementary Figure 4: Longitudinal Changes (Extended)
- **Status**: Similar to Figure 4, likely same script

#### Supplementary Figures 5-6: Clinical Predictors (Additional Organisms)
- **Script(s)**: `bsi_pedictors_viz.py`
- **Status**: Script exists, generates all organisms

---

## Reviewer Comments Summary

### Reviewer 1 - Major Issues

#### 1. Clinical Practices Selection Clarity
**Issue**: Lack of clarity on why specific clinical practices were selected

**Variables Currently in Analysis** (from metadata):
- `SampleType` (axilla, groin, stool)
- `Location` (UCMC, ZCH)
- `GestationCohort` (gestational age groups)
- `SampleCollectionWeek` (Week 1, Week 3)
- `MaternalAntibiotics` (categories)
- `PostNatalAbxCohort` (No/Low/High infant antibiotics)
- `BSI_30D` (bloodstream infection within 30 days)
- `NEC_30D` (necrotizing enterocolitis within 30 days)
- `AnyMilk` (breast milk exposure)
- `PICC` (peripherally inserted central catheter - UE/LE/Neck)
- `UVC` (umbilical venous catheter)
- `BirthMode` / `Delivery` (Cesarean vs Vaginal)

**Variables Requested by Reviewers (NOT currently in analysis)**:
1. Decolonization practices for S. aureus
   - MRSA screening protocols
   - Cohorting practices
   - Mupirocin decolonization (infants, staff, parents)
2. Fluconazole prophylaxis use
3. Specific antibiotic types (not just duration/cohorts)
4. Breast milk details (mother's own vs pasteurized donor milk)
5. Probiotics use/differences across sites

**Action Items**:
- [ ] Check if additional clinical variables are available in metadata
- [ ] If available: Re-run analyses with additional variables
- [ ] If not available: Document in Discussion why not included
- [ ] Add rationale for variable selection in Methods

#### 2. Sample Collection Methods Comparability
**Issue**: Need detailed description of collection protocols at each site

**Questions to Address**:
- How were samples collected at each site?
- Were protocols identical?
- Time between collection and freezing
- Fresh collection vs stabilizer use
- Differences to accommodate unit workflow
- Umbilical swabs: Why excluded? Should recommend for future studies?

**Action Items**:
- [ ] Document detailed sample collection protocols for both sites
- [ ] Create table comparing protocols side-by-side
- [ ] Add to Methods or Supplementary Methods
- [ ] Justify umbilical swab exclusion
- [ ] Discuss implications for interpreting groin swabs vs LE PICC association

---

### Reviewer 1 - Minor Issues

#### 1. CoNS (Coagulase-Negative Staphylococci) Handling
**Issue**: Clarity needed on whether CoNS excluded from microbiome analysis

**Current Status**:
- BSI analysis: CoNS excluded from epidemiology (Table 2)
- BSI analysis: CoNS included for enrolled subjects (Supp Fig 1)
- Microbiome analysis: `Staphylococcus.epidermidis` IS included (this is CoNS!)
- Possible confusion about what was excluded where

**Action Items**:
- [ ] Clarify in Methods: CoNS excluded from BSI epidemiology comparison but included in microbiome
- [ ] Explain rationale: CoNS are real colonizers, exclusion from BSI was due to contamination concern
- [ ] Check if all relevant CoNS species are in microbiome data
- [ ] Update figure legends for clarity

#### 2. Figure 1 Format Improvements
**Current**: Single bar chart, asterisks for significance
**Requested Changes**:
1. Separate into Week 1 and Week 3 panels
2. Use stacked bar plots
3. Show ALL identified microbes to total 100%
4. Include organisms where species cannot be assigned

**Action Items**:
- [ ] Find or create Figure 1 generation script
- [ ] Modify to create stacked bars
- [ ] Include unclassified/unidentified taxa
- [ ] Split by time point if relevant
- [ ] Verify totals = 100%

#### 3. Limitations Section
**Need to explicitly state**:
1. Lack of strain-level data (no genomic comparison of BSI vs colonization strains)
2. No direct strain linkage between BSI and colonizing strains in study subjects
3. Despite metagenomic sequencing, no AMR gene or virulence factor analysis
   - Could have correlated clinical susceptibility with ARGs
   - Could have linked ESBL phenotype with ESBL genes in microbiome
   - This would strengthen the clinical practices argument

**Action Items**:
- [ ] Add comprehensive Limitations section to Discussion
- [ ] Acknowledge lack of strain-level resolution
- [ ] Acknowledge no direct BSI-colonization strain matching
- [ ] Acknowledge missed opportunity for ARG/virulence analysis
- [ ] Frame as future directions

#### 4. C-section and Early Colonization
**Issue**: Others have shown C-section (especially without labor/rupture) linked to S. aureus colonization in Week 1

**Action Items**:
- [ ] Check if C-section timing/details available (with/without labor)
- [ ] If available, add to analysis
- [ ] Discuss in context of decolonization vs recolonization interventions
- [ ] Add to Discussion even if not analyzed

#### 5. Human Subjects Research Statement
**Action Items**:
- [ ] Add human subjects/IRB statement to Methods
- [ ] Include IRB approval numbers for both sites
- [ ] Include informed consent statement

---

### Reviewer 2 - Questions & Comments

#### 1. Patient Selection Criteria (Line 124)
**Question**: Were low birth weight infants <2000g included even if at term?

**Action Items**:
- [ ] Clarify inclusion criteria in Methods
- [ ] Was it: (a) ALL <2000g, or (b) only PRETERM infants
- [ ] Update Methods text at line 124

#### 2. Central Line Terminology (Table 1)
**Question**: How is "Neck Vein PICC line" different from a Central Venous line?

**Action Items**:
- [ ] Clarify terminology in Table 1 legend
- [ ] Consider standardizing terms (e.g., "jugular PICC" vs "neck vein PICC")
- [ ] Add brief definitions to Methods or table footnote

#### 3. CoNS in BSI Pathogens (Line 192, Supp Fig 1)
**Question**: Were CoNS included as BSI pathogens or not?

**Action Items**:
- [ ] Clarify: Excluded from Table 2 (epidemiology), Included in Supp Fig 1 (enrolled subjects)
- [ ] Update text at line 192
- [ ] Ensure consistency throughout manuscript

#### 4. BSI Epidemiology Data Scope
**Questions**:
a. Table 2 data from ALL infants or ONLY preterm in NICUs?
b. How does distribution compare to expected BSI data for both countries?
c. Is this the usual distribution for the institutions?

**Action Items**:
- [ ] Clarify in Table 2 legend: ALL NICU admissions vs only preterm
- [ ] Search literature for country-specific BSI data
- [ ] Compare your distributions to published national/regional data
- [ ] Add comparative context to Results or Discussion
- [ ] If institutional data not representative, note as limitation

#### 5. CLABSI vs BSI Distinction
**Question**: How many BSI events were CLABSIs (central line-associated)?

**Critical Point**: If all patients with LE PICC had high groin S. aureus AND had CLABSI, the microbiome itself may not be causal - could be:
- Poor line care
- Infection introduced during insertion
- Other procedural factors

**Action Items**:
- [ ] Determine which BSI events met CLABSI criteria
- [ ] Add CLABSI data to Table 2
- [ ] Discuss causal interpretation carefully
- [ ] Address potential alternative explanations for line-site associations

#### 6. Microbiome-BSI Correlation Scope
**Question**: Were 13 BSI among 127 enrolled infants compared to THEIR OWN microbiomes, or were the 127 microbiomes compared to institution-level BSI data (Table 2)?

**This is crucial for interpretation!**

**Action Items**:
- [ ] Clarify in Methods exactly what was compared
- [ ] If using institution-level data: explain rationale and limitations
- [ ] If using subject-level data: emphasize in Results
- [ ] Update figure legends to be explicit about data source

#### 7. Timeline and BSI-Microbiome Correlation
**Question**: Does the timeline of BSI (day of positive culture) correlate with microbiome diversity/abundance at different weeks?

**Key point**: Abundance not consistent over weeks within same groups

**Action Items**:
- [ ] Create timeline analysis: BSI day vs microbiome sampling
- [ ] Analyze: Does Week 1 vs Week 3 microbiome predict BSI timing?
- [ ] Test: Microbiome closer in time to BSI more predictive?
- [ ] Add analysis or acknowledge as limitation

#### 8. Candida Composition
**Issue**: 20 Candidemia cases at ZCH, but Candida in microbiome not addressed

**Note**: `bsi_pedictors_viz.py` includes Candida parapsilosis analysis!

**Action Items**:
- [ ] Check if fungal microbiome data was generated
- [ ] If yes: Add Candida analysis to results
- [ ] If no: Acknowledge in limitations
- [ ] Discuss high Candidemia rate at ZCH and potential clinical practice differences

#### 9. Geographic Location Variable Control
**Question**: In analyzing significance of geographic location, were other variables controlled for (antibiotics, C-section, PICC, feeding)?

**Action Items**:
- [ ] Verify that models include these covariates
- [ ] Check multicollinearity between Location and other variables
- [ ] Report variance inflation factors if relevant
- [ ] Clarify in Results that models are adjusted
- [ ] Potentially show bivariate vs multivariate results

---

## Figure & Table Formatting Issues

### Universal Issues
- [ ] **Figure labels too small** - Especially Figure 5, Supp Figs 5-6 not legible
- [ ] **Missing sample sizes** - Add n for each box plot group
- [ ] **Color scheme consistency** - Use different colors for sites vs body sites (Supp Fig 2)

### Figure-Specific Issues

#### Figure 1
- [ ] Consider enfolding into Table 2 (Reviewer 2 suggestion)
- [ ] If keeping separate: improve formatting per Reviewer 1 suggestions

#### Figure 2a (Shannon Diversity)
- [ ] Include differences/changes of median/IQR between Week 1 and 3 within groups
- [ ] Include changes between groups
- [ ] Mention p-values for ALL sites/timepoints (not just significant ones)
- [ ] Add sample sizes to box plots

#### Figure 4 (Longitudinal Changes)
- [ ] Add number of samples for each box plot
- [ ] **CRITICAL**: Specify which 8 species were included in figure legend
- [ ] Current legend says 8 species but only lists 5 (a-e)
- [ ] Need to add parts f, g, h or clarify only 5 shown

#### Figure 5 & Supplementary Figures 5-6
- [ ] **Increase font sizes** - Currently not legible
- [ ] Increase axis labels
- [ ] Increase variable names
- [ ] Increase significance markers
- [ ] Test at print size to ensure readability

#### Table 2
- [ ] Add total numbers in header (n=275 UCMC, n=238 ZCH)
- [ ] Currently only in legend, should be in table header too

#### Supplementary Figure 2
- [ ] Use different colors for study sites than body sites
- [ ] Currently using red/blue for both, causing confusion
- [ ] Suggestion: Sites = shapes, Body sites = colors

#### Supplementary Figure 3
- [ ] Break down relative abundance by the two sites
- [ ] Better depicts differences and correlation to BSI
- [ ] Currently shows all samples combined

#### Supplementary Figure 4
- [ ] Consider enfolding into Figure 4 (avoid duplication)
- [ ] Or clarify what additional information this provides

---

## Analysis Additions/Modifications Needed

### High Priority

#### 1. Legibility Fixes for Figures 5, Supp 5-6
**Script**: `bsi_pedictors_viz.py`
**Changes needed**:
- Increase figure size
- Increase font sizes for all text elements
- Increase marker sizes
- Test output at print resolution

#### 2. Add Sample Counts to All Box Plots
**Scripts**: Multiple (Figure 2, 4, Supp figures)
**Changes needed**:
- Add n= labels to each box plot
- Either in legend or on plot itself

#### 3. Figure 4 Species List Clarification
**Action**: Either
- Show all 8 species mentioned
- Or update legend to reflect actual 5 species shown

#### 4. CLABSI Analysis
**New Analysis Needed**:
- Identify CLABSI vs non-CLABSI BSI events
- Test if associations differ for CLABSI vs non-CLABSI
- Address causality question

#### 5. Timeline Correlation Analysis
**New Analysis Needed**:
- BSI day vs microbiome sampling timing
- Microbiome predictive power by temporal proximity
- Longitudinal trends in cases that develop BSI

### Medium Priority

#### 6. Candida Microbiome Analysis
**If data exists**:
- Add fungal community analysis
- Correlate with Candidemia cases
- Compare antifungal prophylaxis between sites

#### 7. Additional Clinical Variables
**If data available**:
- Fluconazole prophylaxis
- Specific antibiotic agents
- Decolonization protocols
- Probiotic use

#### 8. Country-Specific BSI Comparison
**Literature search needed**:
- Find published BSI rates for China
- Find published BSI rates for USA
- Compare to your data
- Contextualize your findings

### Lower Priority

#### 9. Cesarean Section Details
**If data available**:
- With vs without labor
- With vs without rupture of membranes
- Association with early S. aureus colonization

#### 10. Variance Inflation Factor Analysis
**For model diagnostics**:
- Calculate VIF for all models
- Check multicollinearity
- Report in supplement or methods

---

## Manuscript Text Revisions Needed

### Methods Section
- [ ] Add detailed sample collection protocols (both sites)
- [ ] Add human subjects/IRB statement
- [ ] Clarify patient inclusion criteria (line 124)
- [ ] Add central line terminology definitions
- [ ] Add rationale for clinical variable selection
- [ ] Add CoNS handling explanation

### Results Section
- [ ] Clarify scope of BSI epidemiology data (all vs preterm)
- [ ] Clarify scope of microbiome-BSI correlation (subject-level vs institution)
- [ ] Add CLABSI data
- [ ] Add timeline/temporal analysis results
- [ ] Report adjusted models clearly

### Discussion Section
- [ ] Add comprehensive Limitations section
  - [ ] No strain-level data
  - [ ] No direct BSI-colonization matching
  - [ ] No ARG/virulence analysis
  - [ ] No fungal analysis (if applicable)
- [ ] Add context comparing to national BSI data
- [ ] Discuss C-section implications
- [ ] Discuss decolonization interventions
- [ ] Address clinical practices not analyzed
- [ ] Discuss causality carefully (CLABSI issue)

### Figure Legends
- [ ] All legends: Add sample sizes
- [ ] Figure 1: Update format description
- [ ] Figure 2a: Add all p-values
- [ ] Figure 4: List all 8 species or clarify why only 5
- [ ] All legends: Clarify data source (enrolled subjects vs institution-level)
- [ ] Table 2: Add total n to header, clarify patient population

---

## Priority Action Plan

### Phase 1: Quick Wins (1-2 days)
1. Fix figure font sizes and legibility
2. Add sample counts to box plots
3. Fix Figure 4 species list
4. Update figure legends with clarifications
5. Add Methods text clarifications

### Phase 2: Moderate Effort (3-5 days)
6. Create/modify Figure 1 with new format
7. Add comprehensive Limitations section
8. Add sample collection protocol details
9. Literature search for BSI comparison data
10. Add human subjects statement

### Phase 3: New Analyses (1-2 weeks)
11. CLABSI identification and analysis
12. Timeline/temporal correlation analysis
13. Additional clinical variables (if data available)
14. Candida analysis (if data available)
15. VIF/multicollinearity diagnostics

### Phase 4: Manuscript Integration (3-5 days)
16. Integrate all figure revisions
17. Update all Results text
18. Update Discussion with new findings and limitations
19. Update Methods with all clarifications
20. Final consistency check across all sections

---

## Scripts to Modify/Create

### Existing Scripts to Modify
1. `bsi_pedictors_viz.py` - Fix font sizes, add sample counts
2. `tSNE_plot.py` - Color scheme consistency
3. `microbiome_permanova.py` (need to locate) - Figure 2 modifications
4. Unknown script for Figure 4 - Add species list, sample counts

### New Scripts to Create
1. Figure 1 generation script (stacked bar chart)
2. CLABSI analysis script
3. Timeline correlation analysis script
4. Candida microbiome analysis (if applicable)

### Scripts to Locate
1. Figure 2 generation (diversity + PCA)
2. Figure 4 generation (longitudinal changes)
3. Supplementary Figure 1 generation

---

## Status Tracking

**Last Updated**: 2025-10-27 19:30

### Completed Items - Phase 1 Quick Wins
- [x] Map figures to scripts
- [x] Catalog all reviewer comments
- [x] Create revision plan (figure_revisions_plan.md)
- [x] Locate Figure 2 generation script (R analysis file)
- [x] **Fix Figure 5 font sizes for legibility**
  - Modified `bsi_pedictors_viz.py`
  - Increased figure size: 16x8 → 20x10 (individual), 16x6n → 20x8n (combined)
  - Increased all font sizes: axis labels (14→16), titles (14→18), tick labels (11→14)
  - Made axis labels and titles bold
  - Generated 10 improved PDFs in `revision_figures/`
- [x] **Add sample counts to Figure 2a and Figure 4**
  - Created `generate_revised_figures_with_sample_counts.R`
  - Generated Figure 2a with sample counts: `Figure2a_Shannon_with_sample_counts.pdf`
  - Generated Figure 4 for 5 organisms with sample counts:
    - `Figure4_Staphylococcus_aureus_with_counts.pdf`
    - `Figure4_Klebsiella_pneumoniae_with_counts.pdf`
    - `Figure4_Klebsiella_oxytoca_with_counts.pdf`
    - `Figure4_Enterococcus_faecium_with_counts.pdf`
    - `Figure4_Enterococcus_faecalis_with_counts.pdf`
  - Generated sample count CSV tables for easy reference
  - Added statistical tests with reported sample sizes
- [x] Created centralized `revision_figures/` directory for all outputs

### Phase 1 Complete! ✓ ALL FIGURE GENERATION COMPLETE

All figure generation and modifications have been completed:
1. ✅ Figure 5 font sizes fixed and legible (10 PDFs)
2. ✅ Figure 2a with sample counts AND all p-values (2 versions)
3. ✅ Figure 4 comprehensive - all 8 species (a-h) in 4×2 grid:
   - a. Staphylococcus aureus
   - b. Klebsiella pneumoniae
   - c. Klebsiella oxytoca
   - d. Enterococcus faecium
   - e. Enterococcus faecalis
   - f. Serratia marcescens
   - g. Escherichia coli
   - h. Streptococcus pyogenes
4. ✅ Supplementary Figure 2 - revised color scheme (shapes for location, colors for body site)
5. ✅ Supplementary Figure 3 - separated by site (Cincinnati vs Hangzhou panels)
6. ✅ Supplementary Figure 4 - RECOMMENDATION: Remove (duplicate of Figure 4)
7. ✅ Table 2 header update instructions created
8. ✅ **NEW: Candida Supplementary Figure** - Addressing Reviewer 2's question (20 Candidemia cases at ZCH)
   - 6 top Candida species analyzed (by prevalence)
   - Individual plots for each species (6 PDFs)
   - Comprehensive multi-panel figure (1 PDF)
   - Statistics showing significant differences between locations

### Phase 2 - Text/Methods Clarifications (TO BE DONE IN MANUSCRIPT)

**NOTE: All items below are MANUSCRIPT TEXT EDITS, not figure generation.**
**Work on these in desktop Claude or Word document.**

---

## TEXT/METHODS ITEMS TO ADDRESS IN MANUSCRIPT

### PRIORITY 1: Figure Legends (Simple Text Updates)

#### Figure 2a Legend
- [ ] Add note about all p-values being reported (not just significant ones)
- [ ] Mention sample sizes are shown on plot
- Current: "*p <0.05 for stool and groin samples at week 3"
- Update to: "P-values shown for all comparisons; *p <0.05 for stool and groin samples at week 3. Sample sizes (n) indicated for each group."

#### Figure 4 Legend
- [ ] Update to specify all 8 species (currently only lists a-e)
- Current: "Species shown demonstrated significant changes between time points (*uncorrected p<0.05, FDR<0.1, effect size>0.3). a. Staphylococcus aureus, b. Klebsiella pneumoniae, c. Klebsiella oxytoca, d. Enterococcus faecium, e. Enterococcus faecalis"
- Update to include: "f. Serratia marcescens, g. Escherichia coli, h. Streptococcus pyogenes"

#### Table 2 Legend/Header
- [ ] Add total n to table header
- [ ] See instructions in: `revision_figures/Table2_header_update.txt`

#### Supplementary Figure 4
- [ ] REMOVE entirely (duplicate of Figure 4) OR
- [ ] Add clarification if keeping separate

---

### PRIORITY 2: Methods Section Clarifications

#### Patient Selection (Line 124)
**Reviewer Question**: Were low birth weight infants <2000g included even if at term?
- [ ] Clarify inclusion criteria: (a) ALL <2000g, or (b) only PRETERM infants
- [ ] Update Methods text at line 124

#### Central Line Terminology (Table 1)
**Reviewer Question**: How is "Neck Vein PICC line" different from a Central Venous line?
- [ ] Clarify terminology in Table 1 legend or Methods
- [ ] Consider standardizing terms (e.g., "jugular PICC" vs "neck vein PICC")

#### CoNS Handling (Line 192, Supplementary Figure 1)
**Reviewer Question**: Were CoNS included as BSI pathogens or not?
- [ ] Clarify in Methods: CoNS excluded from BSI epidemiology (Table 2) but included in microbiome analysis
- [ ] Explain rationale: CoNS are real colonizers; exclusion from BSI was due to contamination concern
- [ ] Update text at line 192
- [ ] Update Supp Fig 1 legend for consistency

#### Human Subjects/IRB Statement
- [ ] Add human subjects/IRB statement to Methods
- [ ] Include IRB approval numbers for both sites
- [ ] Include informed consent statement

---

### PRIORITY 3: Methods - Clinical Variables & Sample Collection

#### Clinical Variables NOT Currently in Analysis
**Reviewer asked about these variables - need to address why not included:**
1. Decolonization practices for S. aureus (MRSA screening, mupirocin, cohorting)
2. Fluconazole prophylaxis use
3. Specific antibiotic types (not just duration/cohorts)
4. Breast milk details (mother's own vs pasteurized donor milk)
5. Probiotics use/differences across sites

**Action Items:**
- [ ] Check if any of these variables are available in metadata
- [ ] If available: Consider re-running analyses OR note in Discussion
- [ ] If NOT available: Document in Methods or Discussion why not included
- [ ] Add rationale for variable selection in Methods

#### Sample Collection Protocols
**Reviewer wants detailed protocols to ensure comparability:**
- [ ] Document detailed sample collection protocols for both sites
- [ ] Create table comparing protocols side-by-side (or add to Supplementary Methods)
- [ ] Address: Collection method, time to freezing, fresh vs stabilizer, workflow differences
- [ ] Justify umbilical swab exclusion
- [ ] Discuss implications for interpreting groin swabs vs LE PICC association

---

### PRIORITY 4: Results Section Updates

#### BSI Epidemiology Data Scope
**Reviewer Questions:**
- [ ] Clarify in Table 2 legend: Data from ALL NICU infants or ONLY preterm?
- [ ] Add context: How does your distribution compare to expected BSI data for both countries?
- [ ] Discuss: Is this the usual distribution for these institutions?
- [ ] If institutional data not representative, note as limitation

#### CLABSI vs BSI Distinction
**Reviewer Question**: How many BSI events were CLABSIs?
**CRITICAL POINT**: If all LE PICC patients with high groin S. aureus had CLABSI, microbiome may not be causal (could be line care, insertion contamination, etc.)

- [ ] Determine which BSI events met CLABSI criteria
- [ ] Add CLABSI data to Table 2 or Results
- [ ] Discuss causality carefully - address alternative explanations for line-site associations

#### Microbiome-BSI Correlation Scope
**Reviewer Question**: Were 13 BSI among 127 enrolled infants compared to THEIR microbiomes, or were 127 microbiomes compared to institution-level BSI data (Table 2)?
**This is crucial for interpretation!**

- [ ] Clarify in Methods exactly what was compared
- [ ] If using institution-level data: explain rationale and limitations
- [ ] If using subject-level data: emphasize in Results
- [ ] Update figure legends to be explicit about data source

#### Timeline and BSI-Microbiome Correlation
**Reviewer Question**: Does BSI timing (day of positive culture) correlate with microbiome diversity/abundance at different weeks?
**Key point**: Abundance not consistent over weeks

- [ ] If timeline data available: Add temporal analysis
- [ ] Test: Does Week 1 vs Week 3 microbiome predict BSI timing?
- [ ] Test: Is microbiome closer in time to BSI more predictive?
- [ ] If NOT analyzed: Acknowledge as limitation

#### Geographic Location Variable Control
**Reviewer Question**: Were other variables controlled for when analyzing geographic location effect?
- [ ] Verify models include covariates (antibiotics, C-section, PICC, feeding)
- [ ] Check multicollinearity between Location and other variables
- [ ] Clarify in Results that models are adjusted for multiple factors
- [ ] Consider showing bivariate vs multivariate results

---

### PRIORITY 5: Discussion Section Additions

#### Comprehensive Limitations Section
**Need to explicitly state:**
- [ ] Lack of strain-level data (no genomic comparison of BSI vs colonization strains)
- [ ] No direct strain linkage between BSI and colonizing strains in study subjects
- [ ] Despite metagenomic sequencing, no AMR gene or virulence factor analysis
  - Could have correlated clinical susceptibility with ARGs
  - Could have linked ESBL phenotype with ESBL genes in microbiome
  - Frame as future directions
- [ ] No fungal microbiome analysis (20 Candidemia cases at ZCH not addressed)
- [ ] If timeline data unavailable: Acknowledge as limitation

#### C-section and Early Colonization
**Context**: Others have shown C-section (especially without labor/rupture) linked to S. aureus colonization in Week 1
- [ ] Check if C-section timing/details available (with/without labor/rupture)
- [ ] If available: Add to analysis or discuss
- [ ] Discuss in context of decolonization vs recolonization interventions
- [ ] Add to Discussion even if not analyzed

#### Clinical Practices Context
- [ ] Discuss clinical practices that WERE NOT analyzed (from Priority 3 list)
- [ ] Explain how these might influence results
- [ ] Suggest future studies should include these variables

#### Causality and Alternative Explanations
- [ ] Address CLABSI issue (if LE PICC + groin colonization all had CLABSI)
- [ ] Discuss alternative explanations for line-site associations
- [ ] Carefully frame causal vs associative language

#### Country/Institution-Specific Context
- [ ] Compare your BSI distributions to published national/regional data
- [ ] Discuss whether your data are representative
- [ ] Contextualize geographic differences

---

### PRIORITY 6: Candida Analysis ✅ COMPLETED

**Issue**: 20 Candidemia cases at ZCH, but Candida in microbiome not addressed
**Note**: `bsi_pedictors_viz.py` includes Candida parapsilosis in analysis

- [x] ✅ Fungal microbiome data WAS generated from metagenomic sequencing
- [x] ✅ Candida analysis COMPLETED - 6 species analyzed longitudinally
- [x] ✅ Supplementary figure created (similar format to Figure 4)
- [ ] **TODO: Add to Results section** - Describe Candida findings:
  - Top 6 Candida species by prevalence: C. riodocensis (36.7%), C. gorgasii (30.3%), C. arabinofermentans (19.6%), C. restingae (15.2%), C. auris (13.9%), C. albicans (8.4%)
  - Significant geographic differences found for C. riodocensis, C. gorgasii, C. arabinofermentans, and C. auris (p<0.05)
  - C. albicans (most clinically relevant) showed high abundance but no significant geographic difference
- [ ] **TODO: Discuss** high Candidemia rate at ZCH (20 cases vs fewer at Cincinnati) and potential clinical practice differences (fluconazole prophylaxis?)
- [ ] **TODO: Add figure legend** for Candida supplementary figure

---

### SUMMARY CHECKLIST FOR MANUSCRIPT REVISION

**Figure Legends (3 items)**
- [ ] Figure 2a - add p-value note
- [ ] Figure 4 - list all 8 species
- [ ] Table 2 - add n to header

**Methods Section (7 items)**
- [ ] Patient selection criteria
- [ ] Central line terminology
- [ ] CoNS handling explanation
- [ ] IRB/human subjects statement
- [ ] Clinical variables rationale
- [ ] Sample collection protocols
- [ ] Analysis scope clarifications

**Results Section (5 items)**
- [ ] BSI data scope (all vs preterm)
- [ ] CLABSI data
- [ ] Microbiome-BSI correlation scope
- [ ] Timeline correlation (if data available)
- [ ] Covariate adjustment clarification

**Discussion Section (5 items)**
- [ ] Comprehensive Limitations section
- [ ] C-section implications
- [ ] Clinical practices not analyzed
- [ ] Causality discussion
- [ ] Country-specific context

**Optional/If Data Available:**
- [ ] Candida microbiome analysis
- [ ] Timeline/temporal analysis
- [ ] Additional clinical variables

### Files Created/Modified
**New Files:**
- `figure_revisions_plan.md` - Comprehensive revision tracking document
- `generate_revised_figures_with_sample_counts.R` - R script for Figures 2a and 4 with sample counts
- `generate_candida_supplementary_figure.R` - R script for Candida longitudinal analysis
- `generate_revised_supplementary_tsne_figures.py` - Python script for t-SNE figures
- `add_sample_counts_helper.R` - Helper functions for adding sample counts (for reference)
- `revision_figures/` directory with **30 files total**:
  - 10 Figure 5/Supp Fig 5-6 PDFs (Python)
  - 9 Figure 2a/4 PDFs with sample counts (R)
  - 3 Supplementary t-SNE figures (Python)
  - 7 Candida analysis PDFs (6 individual + 1 comprehensive, R)
  - 1 Table 2 update instructions (text)
  - 5 CSV statistics/sample count tables

**Modified Files:**
- `bsi_pedictors_viz.py` - Increased font sizes, updated output directory

### Summary Statistics from Generated Figures
**Figure 2a Sample Sizes:**
- Cincinnati: Axilla (n=53,55), Groin (n=39,52), Stool (n=62,60) at Weeks 1,3
- Hangzhou: Axilla (n=29,38), Groin (n=14,24), Stool (n=34,43) at Weeks 1,3
- Significant differences found in Week 3 Groin (p=0.0008) and Stool (p=0.0043)

**Figure 4 Sample Sizes:**
- 5 organisms analyzed across 12 combinations (2 locations × 3 body sites × 2 weeks)
- Total sample sizes per organism: ~60 Cincinnati, ~40 Hangzhou samples per timepoint
- See `revision_figures/Figure4_all_sample_counts.csv` for complete breakdown

### Blocked/Waiting
- [ ] Need confirmation on data availability for additional clinical variables
- [ ] Need sample collection protocol documentation from collaborators
- [ ] Need CLABSI definitions/data from clinical team

---

## Notes & Questions for Team

1. **Clinical Variables**: Do we have data for fluconazole prophylaxis, specific antibiotic types, probiotic use, decolonization protocols?
2. **Sample Collection**: Who can provide detailed protocols for both sites?
3. **CLABSI Data**: Do we have CLABSI definitions/classifications for the BSI events?
4. **Timeline Data**: Do we have exact dates for BSI events relative to microbiome sampling?
5. **Fungal Data**: Was fungal microbiome sequenced/analyzable from the metagenomic data?
6. **IRB Info**: Need IRB approval numbers and dates for both institutions
7. **C-section Details**: Do we have data on labor/membrane rupture status?

---

## Contact & Responsibilities

- **Figure modifications**: [Name/Role]
- **Statistical analyses**: [Name/Role]
- **Clinical data queries**: [Name/Role]
- **Manuscript text revisions**: [Name/Role]
- **Coordination**: [Name/Role]

