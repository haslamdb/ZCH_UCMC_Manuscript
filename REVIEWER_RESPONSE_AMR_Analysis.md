# Response to Reviewer: Antibiotic Resistance Gene Analysis

## Reviewer Comment:
> "Despite metagenomic sequencing, no analysis of Antibiotic resistance genes or virulence factors that could add more resolution to the colonization patterns was done. This could be especially insightful, given the significant differences in antibiotic use patterns. The information on BSI could be correlated with clinical susceptibility data. For example, ESBL phenotype-related ARGs being more abundant in the microbiome and linked with ESBL K. pneumoniae can be a powerful piece of evidence for the point the authors are making. If not done, this is a limitation that should be noted."

---

## Proposed Response:

We thank the reviewer for this excellent suggestion. In response to this comment, we performed a comprehensive analysis of antibiotic resistance genes (ARGs) across our metagenomic dataset. This analysis examined 237 resistance genes across 669 samples from both UCMC (Cincinnati) and ZCH (Hangzhou) NICUs, including longitudinal samples from 127 infants.

### Key Findings from Resistome Analysis:

**1. No Geographic Differences Despite Different Antibiotic Use Patterns**

Surprisingly, despite documented differences in antibiotic prescribing practices between the two institutions, we found no significant differences in resistome composition between UCMC and ZCH (differential abundance testing: 0/237 genes significant at FDR < 0.05; mixed-effects model location effect p = 0.45). This remarkable similarity in resistance gene profiles across geographically and clinically distinct NICUs suggests that NICU resistomes may be shaped more by universal factors (e.g., body site, temporal dynamics) than by institution-specific antibiotic pressure.

**2. No Correlation Between Antibiotic Exposure and Resistance Burden**

We tested for correlations between cumulative antibiotic exposure and AMR burden across all samples and found no significant associations (Spearman rho = 0.001, p = 0.985). This null finding was consistent across all body sites (Axilla: rho = -0.093, p = 0.19; Groin: rho = 0.073, p = 0.30; Stool: rho = -0.077, p = 0.28), suggesting that AMR colonization in this NICU cohort may be largely independent of direct antibiotic selection pressure during hospitalization.

**3. Body Site is the Primary Driver of Resistome Composition**

Principal component analysis revealed that body site accounts for the largest variance in resistome composition (PC1: 30.8%), with stool samples harboring dramatically higher AMR burden than skin sites (mean AMR RPM: Stool = 4,888; Groin = 2,553; Axilla = 694; mixed-effects model: Stool coefficient = +0.615, p < 0.0001). This body site effect was >7-fold stronger than any other factor in our analysis.

**4. Dramatic Temporal Increases in Groin and Stool Sites**

From Week 1 to Week 3, we observed dramatic increases in AMR burden in groin (+204%, p = 0.0001) and stool (+64%, p < 0.0001) samples, while axilla samples remained stable (+0.6%, p = 0.96). At the gene level, 120 resistance genes significantly increased in groin and 104 increased in stool, indicating widespread resistance gene expansion during early NICU colonization.

**5. Fluoroquinolone Resistance Analysis**

We also analyzed fluoroquinolone (FQ) resistance mutations across 9 bacterial species. Consistent with the ARG findings, we observed: (i) no location effect (p = 0.79), (ii) strong body site stratification, and (iii) high-frequency resistance mutations in key pathogens including *Klebsiella pneumoniae* gyrA K154R (present in 99.5% of samples with *K. pneumoniae*), *Staphylococcus aureus* parE Y470N (99.8%), and fixed mutations in *Enterococcus* species (99.7-99.8%).

---

### Why We Are Not Including This in the Main Manuscript:

While this resistome analysis provides valuable insights into NICU microbial ecology, we believe it extends beyond the core narrative of our manuscript, which focuses on the remarkable convergence of microbiome composition between geographically distinct NICUs. The ARG analysis—particularly the surprising finding that antibiotic resistance patterns are similar despite different prescribing practices—represents a distinct and complex story that would dilute the focus of our current work.

However, we recognize the scientific value of these findings and would be happy to include 2-3 key figures in the supplementary materials to demonstrate the comprehensiveness of our analysis. These would include:

1. **Supplementary Figure X**: Resistome composition by location and body site (showing no geographic differences)
   - **File**: `nicu_resistome_analysis/figures/publication/Figure2_BodySite_Comparison.pdf`
   - **Shows**: AMR burden (RPM) stratified by location (UCMC vs ZCH), body site, and week
   - **Key message**: Nearly identical resistome profiles between institutions despite different prescribing patterns

2. **Supplementary Figure Y**: Longitudinal trajectories of AMR burden by body site (showing dramatic increases in groin/stool)
   - **File**: `nicu_resistome_analysis/figures/publication/Figure3_Longitudinal_Trajectories.pdf`
   - **Shows**: Paired Week 1→Week 3 trajectories for each body site with statistical tests
   - **Key message**: Groin (+204%) and stool (+64%) show dramatic AMR expansion while axilla remains stable

3. **Supplementary Figure Z** (RECOMMENDED - addresses reviewer's antibiotic use concern directly): Antibiotic exposure vs. AMR burden correlations
   - **File**: `nicu_resistome_analysis/figures/exploratory/antibiotic_amr_correlation.pdf`
   - **Shows**: Scatter plots of cumulative antibiotic exposure vs. total AMR burden by body site
   - **Key message**: No correlation between antibiotic exposure and resistance (Spearman rho ≈ 0, all p > 0.19), suggesting resistance colonization is independent of direct selection pressure

A brief description of the resistome findings could be added to the Results section (1-2 sentences) with reference to the supplementary figures, such as:

> "Comprehensive analysis of antibiotic resistance genes revealed no significant differences between UCMC and ZCH resistomes (0/237 genes significant, FDR < 0.05), despite different institutional antibiotic prescribing patterns. Body site and temporal dynamics, rather than geographic location or antibiotic exposure, were the primary drivers of resistome composition (Supplementary Figures X-Y)."

---

### Addressing the Specific Example: ESBL-Related ARGs and K. pneumoniae

The reviewer mentions ESBL-related ARGs as a potential link to clinical K. pneumoniae BSI. While we recognize the value of this approach, direct correlation between metagenomic resistome data and clinical BSI isolate susceptibility was not possible in our study design. The microbiome sequencing was performed on a subset of infants enrolled in the study, while BSI data was collected from all infants in each NICU during the study period. Furthermore, only a small number of the infants with microbiome data developed BSI, precluding robust statistical comparisons between colonizing resistome patterns and bloodstream infection isolate phenotypes.

However, our fluoroquinolone resistance analysis at the population level demonstrates that resistance mutations in key NICU pathogens such as *K. pneumoniae* are nearly universal (gyrA K154R present in 99.5% of samples with *K. pneumoniae*) and do not differ by location (p = 0.79). This suggests that resistance phenotypes commonly observed in clinical NICU isolates likely reflect widespread pre-existing colonization patterns rather than institution-specific antibiotic selection. A future study design with prospective enrollment of infants who develop BSI would be needed to directly link metagenomic resistome signatures with clinical culture susceptibility profiles.

---

### Limitation Statement (if ARG analysis is not included):

If the editor/reviewer prefers that we do not include the ARG analysis even in supplementary materials, we will add the following limitation to the Discussion:

> "Our analysis focused on taxonomic composition and diversity metrics. While metagenomic sequencing also provides information on antibiotic resistance genes and virulence factors, comprehensive analysis of these functional elements was beyond the scope of the current study. Future work examining resistome dynamics in relation to clinical antibiotic susceptibility patterns and bloodstream infection etiology would provide valuable insights into NICU colonization-infection pathways."

---

## Summary for Reviewer:

- **Analysis performed**: ✓ Comprehensive resistome analysis completed
- **Key finding**: No geographic differences despite different antibiotic use (surprising!)
- **No antibiotic correlation**: Resistance appears independent of direct selection pressure
- **Body site dominates**: Stool >> Groin > Axilla (7-fold effect)
- **Temporal dynamics**: Dramatic increases in groin/stool Week 1→3
- **Recommendation**: Include 2-3 supplementary figures with brief Results mention, but keep detailed analysis out of main manuscript to maintain focus

---

## Materials Available for Supplementary Submission:

All analyses, figures, and summary tables are available in the `nicu_resistome_analysis/` directory and can be prepared for supplementary materials upon request.
