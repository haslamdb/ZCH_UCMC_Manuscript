# NICU Microbiome Project

A comprehensive project for analyzing microbiome composition in neonatal intensive care unit (NICU) patients from Cincinnati and Hangzhou, with a focus on factors influencing microbiome development and its relationship to bloodstream infection (BSI).

## Project Overview

This repository contains R and Python scripts for analyzing taxonomic profiles from Kraken2 and Bracken, with tools for statistical analysis and visualization of NICU patient microbiome data. The project compares microbiome composition between different body sites, gestational ages, antibiotic exposures, and locations (Cincinnati vs Hangzhou).

### Key Research Questions

1. How do factors like gestational age, postnatal age, and antibiotic exposure influence the early-life microbiome?
2. What are the differences in microbiome composition between Cincinnati and Hangzhou NICU patients?
3. What taxa are associated with bloodstream infections (BSI) in NICU patients?
4. How does microbiome development differ across various body sites (stool, axilla, groin)?

## Repository Structure

```
.
├── R_scripts/                      # R analysis scripts
│   ├── bsi_microbiome_comparison.R # Analysis comparing BSI and microbiome data
│   ├── dbRDA.R                     # Distance-based redundancy analysis
│   ├── genus_analysis.R            # Genus-level analysis
│   └── import_initial_species_analysis.r # Initial species data import and processing
│
├── python_scripts/                 # Python analysis scripts
│   ├── feature_selection_rf_microbiome_composition.py # Random Forest feature selection
│   ├── microbiome_transform.py     # Microbiome data transformation utilities
│   ├── microbiome-variance-partitioning-revised.py # Variance partitioning analysis
│   ├── network_plot_after_variance_partitioning.py # Network visualization
│   ├── rf_shap_LMM_analysis_key_orgs.py # Random Forest with SHAP and LMM analysis
│   ├── shap_feature_importance_revised.py # Enhanced SHAP analysis
│   ├── shap_feature_selection.py   # Feature selection with SHAP
│   ├── tSNE_plot.py                # t-SNE visualizations
│   └── zero-inflated-glmm.py       # Zero-inflated GLMM for microbiome data
│
├── bash_scripts/                   # Shell scripts for processing
│   └── process_reads.sh            # Script to process raw microbiome reads
│
├── utils/                          # Utility scripts and configuration
│   ├── config-file.txt             # Configuration parameters
│   ├── python-config-loader.py     # Python configuration utilities
│   ├── r-config-loader.r           # R configuration utilities
│   └── [other utility files]
│
└── README.md                       # This file
```

## Data Files

The analysis scripts reference the following key data files (not included in the repository):

- `AllNICUSampleKeyRevised20250206_for_HangzhouCincinnatiSamples.csv`: Sample metadata
- `BSIData20250206.csv`: Bloodstream infection data
- `HumanReactiveKraken2.csv`: List of human-reactive species to filter
- `NICUSpeciesReduced.csv`: Processed microbiome count data

## Getting Started

### Prerequisites

#### R Dependencies
```R
install.packages(c("tidyverse", "vegan", "labdsv", "RColorBrewer", "Heatplus", 
                   "ggplot2", "FactoMineR", "factoextra", "ade4", "scales", 
                   "pheatmap", "VennDiagram", "permute", "lme4"))
```

#### Python Dependencies
```bash
pip install pandas numpy scipy scikit-learn scikit-bio shap matplotlib seaborn
```

### Running the Analysis

#### R Analysis Pipeline
1. Process species data:
```R
source("R_scripts/import_initial_species_analysis.r")
```

2. Run BSI and microbiome comparison:
```R
source("R_scripts/bsi_microbiome_comparison.R")
```

3. Perform genus-level analysis:
```R
source("R_scripts/genus_analysis.R")
```

4. Run distance-based redundancy analysis:
```R
source("R_scripts/dbRDA.R")
```

#### Python Analysis Pipeline
1. Transform microbiome data:
```bash
python python_scripts/microbiome_transform.py
```

2. Run variance partitioning analysis:
```bash
python python_scripts/microbiome-variance-partitioning-revised.py
```

3. Perform feature selection and importance analysis:
```bash
python python_scripts/rf_shap_LMM_analysis_key_orgs.py
```

4. Generate visualizations:
```bash
python python_scripts/tSNE_plot.py
```

## Key Features

- **Taxonomic Analysis**: Process Kraken2 and Bracken outputs for taxonomic profiling
- **Statistical Methods**: 
  - PERMANOVA for community-level differences
  - Generalized Linear Mixed Models (GLMM) with subject-level random effects
  - Zero-inflated models for sparse microbiome data
- **Feature Selection**:
  - Random Forest importance metrics
  - SHAP (SHapley Additive exPlanations) values for interpretable ML
- **Visualization**:
  - t-SNE plots for microbiome composition
  - Network visualization of clinical factor interactions
  - Taxonomic heatmaps and PCA/PCoA plots

## Key Results

### Body Site Differences
The analysis reveals distinct microbial communities between stool, axilla, and groin samples, with each body site hosting specific signature taxa.

### Gestational Age Effects
Gestational age significantly impacts microbiome composition, with differences observed between different gestational age cohorts (23-27 weeks, 28-32 weeks, 33-36 weeks, 37-42 weeks).

### Antibiotic Effects
Postnatal antibiotic exposure significantly alters the microbiome, with reduced abundance of beneficial bacteria like Bifidobacterium and increased prevalence of potential pathogens.

### Location Differences
Significant differences in microbiome composition were observed between Cincinnati and Hangzhou NICU patients, potentially reflecting differences in clinical practices, environment, or genetics.

### Bloodstream Infection Associations
Several taxa were found to be differentially abundant between infants with and without bloodstream infections, suggesting potential biomarkers for infection risk.

## License

This project is proprietary research code. All rights reserved.

## Authors

- David Haslam (Primary Investigator)

## Acknowledgments

- Cincinnati Children's Hospital Medical Center
- [Collaborating institutions]
