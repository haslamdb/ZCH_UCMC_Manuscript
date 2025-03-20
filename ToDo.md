# ToDo List for ZCH_UCMC_Manuscript

## 🐞 Bug Fixes
- [ ] Set up proper paths for databases in `bash_scripts/process_reads.sh`


## 🔧 Features to Improve

## High Priority Tasks

### Data Processing
- [ ] Create a consistent file naming scheme for processed data

### Analysis Improvements
- [ ] Refactor the code in `genus_analysis.R` to improve organization and readability
- [ ] Address redundant code in `genus_analysis.R` and `bsi_microbiome_comparison.R` 
- [ ] Consolidate duplicate effect size calculation functions across scripts
- [ ] Integrate the `nicu_utils.R` functions more consistently across all analysis scripts

### Statistical Analysis
- [ ] Perform false discovery rate correction more consistently across all comparisons
- [ ] Implement multivariate models to account for confounding factors - GLMM, Random Forest, SHaP
- [ ] Add power analysis to determine if sample sizes are adequate for conclusions
- [ ] Verify assumptions of statistical tests being used

### Visualization Enhancements
- [ ] Create consistent figure styles across all analysis scripts
- [ ] Generate publication-quality heatmaps for differentially abundant taxa
- [ ] Improve PCA plots to better distinguish between groups
- [ ] Develop interactive visualizations for key findings

## Medium Priority Tasks

### Documentation
- [ ] Add detailed docstring headers to all R functions
- [ ] Create a comprehensive documentation file explaining the analysis pipeline
- [ ] Document parameter choices and thresholds used in each analysis step
- [ ] Update README.md with clearer installation and usage instructions

### Code Organization
- [ ] Refactor repetitive code into functions in `nicu_utils.R`
- [ ] Standardize function parameter names across all scripts
- [ ] Add proper error handling to data import functions
- [ ] Create a unified configuration file for all analysis parameters

### Additional Analyses
- [ ] Compare microbial diversity results with published NICU studies
- [ ] Analyze correlations between gestational age and specific microbial taxa
- [ ] Implement longitudinal analysis to better track changes over time
- [ ] Add functional prediction analysis (e.g., PICRUSt2)

## Other Tasks

### Repository Structure
- [ ] Organize result outputs into subdirectories by analysis type
- [ ] Create a consistent naming convention for all output files
- [ ] Add version information to all scripts
- [ ] Set up continuous integration for testing code changes

### Feature Additions
- [ ] Add support for additional diversity metrics
- [ ] Implement machine learning approaches for predicting clinical outcomes
- [ ] Create a Shiny app for interactive exploration of results
- [ ] Add support for metabolomic data integration

### Optimization
- [ ] Optimize memory usage for large dataset processing
- [ ] Improve performance of taxonomic abundance calculations
- [ ] Parallelize computationally intensive operations
- [ ] Implement checkpointing for long-running analyses

## For Publication

### Manuscript Preparation
- [ ] Generate final versions of all figures for publication
- [ ] Create supplementary materials with detailed methods
- [ ] Update statistical analyses based on reviewer feedback
- [ ] Prepare data for public repository submission

### Data Sharing
- [ ] Anonymize patient metadata for public sharing
- [ ] Prepare processed data files for repository submission
- [ ] Document data processing steps for reproducibility
- [ ] Create R markdown documents that reproduce key analyses

## Future Directions
- [ ] Extend analysis to include additional NICU sites
- [ ] Add metagenomic analysis to complement 16S rRNA findings
- [ ] Investigate host-microbiome interactions with additional data types
- [ ] Develop predictive models for clinical outcomes based on microbiome profiles

## 📅 Future Plans
- [ ] Test SHAP force plots with new dataset
- [ ] Deploy microbiome analysis pipeline on cloud
