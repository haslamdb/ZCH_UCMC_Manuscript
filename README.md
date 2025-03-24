# Analysis of ZCH and UCMC BSI-Microbiome Associations

## Introduction

This repository contains the code and data analysis pipeline for a comparative microbiome study of neonatal intensive care units (NICUs) in Cincinnati and Hangzhou. The project investigates the relationships between epidemiology of bloodstream infections (BSI) and infant microbiome composition and various clinical factors, including:

- Geographic location (Cincinnati vs. Hangzhou)
- Sample collection site (axilla, groin, stool)
- Sample collection time (Week 1, Week 3)
- Antibiotic exposure (postnatal and maternal)
- Gestational age
- Birth mode
- Breast milk intake
- Intravenous access method (PICC, UVC, peripheral IV)



The analysis examines taxonomic composition, diversity metrics, and associations with clinical variables through a combination of statistical approaches and machine learning techniques.

## Repository Structure

```
ZCH_UCMC_Manuscript/
├── README.md                  # This file
├── metadata/                  # Metadata directory
│   ├── AllNICUSampleKeyRevised*.csv   # Sample metadata
│   ├── BSIData*.csv                   # Bloodstream infection count data
│   └── HumanReactiveKraken2.csv       # Human reactive species list
├── data/                      # Data directory
│   └── Kraken2/               # Kraken2 taxonomic classifications
├── R_scripts/                 # R analysis scripts
├── python_scripts/            # Python analysis scripts
├── bash_scripts/              # Bash processing scripts
└── results/                   # Analysis results directory
    ├── figures/               # Generated figures
    └── tables/                # Generated data tables
```

## Analysis Pipeline

### 1. Raw Data Processing

#### `bash_scripts/process_reads.sh`

This script processes raw FASTQ sequencing files to prepare them for microbiome analysis:

```bash
bash bash_scripts/process_reads.sh <input_dir> <output_dir> <metadata_file>
```

Key steps:
- Removes host (human) DNA using kneaddata
- Performs taxonomic classification with Kraken2
- Estimates relative abundance with Bracken
- Organizes outputs into structured directories

The script leverages the kraken_tools pipeline (https://github.com/haslamdb/kraken_tools) for consistent processing. Output files include taxonomic classifications at different levels (e.g., species, genus) that serve as input for downstream analyses.

### 2. Core R Analysis Scripts

#### `R_scripts/nicu_utils.R`

This utility script provides common functions used throughout the R analyses:

- Custom color palettes and visualization themes
- Statistical functions for effect size calculations
- Data processing functions (normalization, filtering)
- Helper functions for fold change calculations

This script is sourced by other R scripts and serves as a foundation for consistent analysis approaches.

#### `R_scripts/import_initial_species_analysis.R`

This script performs the initial processing and exploratory analysis of microbiome data at the species level:

```R
source("R_scripts/import_initial_species_analysis.R")
```

Key analyses:
- Imports and merges sample metadata with Kraken2 taxonomic data
- Filters and normalizes species abundance data
- Calculates diversity metrics (Shannon, Simpson, Fisher's alpha)
- Performs statistical comparisons between groups
  - Location (Cincinnati vs. Hangzhou)
  - Antibiotic exposure
  - Gestational age cohorts
- Creates diversity visualization plots
- Performs ordination analyses (PCA) to visualize community structure

The output includes filtered species tables, diversity metrics for each sample, and visualization files that identify key differences between groups.

#### `R_scripts/bsi_microbiome_comparison.R`

This script focuses on the relationship between bloodstream infections (BSI) and the microbiome:

```R
source("R_scripts/bsi_microbiome_comparison.R")
```

Key analyses:
- Merges BSI data with microbiome species data
- Computes Bray-Curtis distance matrices to measure sample similarity
- Calculates observed-to-expected ratios for BSI organisms
- Performs PCA to visualize relationships between BSI organisms and microbiome composition
- Generates visualizations showing BSI-differential organisms
- Creates location-specific analyses of *Staphylococcus aureus*, *Klebsiella pneumoniae*, and other clinically important organisms

This analysis helps identify connections between bloodstream pathogens and microbiome composition, potentially revealing early warning markers or protective microbial signatures.

### 3. Machine Learning and Advanced Visualization

#### `python_scripts/rf_shap_LMM_analysis.py`

This script implements a machine learning approach to identify key clinical features associated with microbiome composition:

```bash
python python_scripts/rf_shap_LMM_analysis.py
```

Key analyses:
- Automatically selects appropriate data transformation based on microbiome properties
- Builds Random Forest regression models to predict organism abundance
- Calculates SHAP (SHapley Additive exPlanations) values to quantify feature importance
- Fits mixed-effects models with subject-level random effects to account for repeated measures
- Visualizes feature importance with SHAP plots
- Generates summary statistics and effect estimations for clinical variables

This analysis provides interpretable insights into which clinical factors most strongly influence the abundance of specific microbes.

#### `python_scripts/tSNE_plot.py`

This script creates dimensionality reduction visualizations to explore microbiome patterns:

```bash
python python_scripts/tSNE_plot.py
```

Key features:
- Generates t-SNE (t-Distributed Stochastic Neighbor Embedding) plots to visualize high-dimensional data
- Creates organism-specific visualizations colored by abundance
- Produces metadata-based visualizations (by sample type, location, etc.)
- Outputs both individual and combined visualizations in PDF and PNG formats
- Creates multi-page visualization documents for comprehensive exploration

These visualizations help identify clustering patterns in the data that may not be apparent in traditional statistical analyses.

### 4. Supplementary Analyses

#### Genus-Level Analysis

`R_scripts/genus_analysis.R` - Performs analysis at the genus taxonomic level:
- Filters and processes genus-level abundance data
- Runs generalized linear mixed-effects models (GLMM)
- Analyzes effects of antibiotics, gestational age, and maternal factors
- Calculates effect sizes between comparison groups

#### Cross-Validation and Feature Selection Scripts

- `python_scripts/microbiome_shap_analysis_Kfold.py`: Extends the SHAP analysis with K-fold cross-validation for more robust feature importance
- `python_scripts/feature_selection_rf.py`: Uses Random Forest for feature selection based on Bray-Curtis distances
- `python_scripts/shap_feature_selection.py`: Focuses on SHAP-based feature importance for microbiome differences

#### Advanced Statistical Modeling

`python_scripts/zero-inflated-glmm.py`: Implements zero-inflated generalized linear mixed models:
- Addresses the high proportion of zeros typical in microbiome data
- Implements a two-part model (presence/absence + abundance when present)
- Incorporates subject-level random effects
- Provides coefficient estimates and significance tests

#### Data Transformation Utilities

`python_scripts/microbiome_transform.py`: Module providing functions for microbiome data transformation:
- CLR (Centered Log-Ratio) transformation for compositional data
- TSS (Total Sum Scaling) normalization
- VST (Variance-Stabilizing Transformation)
- Log transformation with pseudocount handling
- Rarefaction (subsampling) to normalize sequencing depth

## Dependencies

### R Packages
- vegan: Ecological diversity analysis and ordination
- ggplot2: Data visualization
- FactoMineR and factoextra: PCA and visualization
- dplyr: Data manipulation
- lme4: Mixed-effects models
- NBZIMM: Zero-inflated negative binomial regression
- tidyverse: Data wrangling
- pheatmap: Heatmap visualization

### Python Packages
- pandas: Data manipulation
- numpy: Numerical operations
- scikit-learn: Machine learning models
- shap: Model interpretability and feature importance
- statsmodels: Statistical modeling
- matplotlib and seaborn: Visualization
- scipy: Scientific computing

### External Tools
- Kraken2: Taxonomic classification
- Bracken: Abundance estimation
- kneaddata: Host DNA removal

## Getting Started

### Setting Up the Repository

1. Clone this repository:
   ```bash
   git clone https://github.com/YOUR-USERNAME/ZCH_UCMC_Manuscript.git
   cd ZCH_UCMC_Manuscript
   ```

2. Create the required directory structure:
   ```bash
   mkdir -p data results/figures results/tables KrakenAlignments/Kraken2
   ```

3. Place your input data in the appropriate directories:
   - Sample metadata in `metadata/`
   - Kraken2 results in `KrakenAlignments/Kraken2/`

### Running the Full Analysis Pipeline

1. Process raw reads:
   ```bash
   bash bash_scripts/process_reads.sh <input_dir> <output_dir> <metadata_file>
   ```

2. Run species-level analysis:
   ```R
   source("R_scripts/import_initial_species_analysis.R")
   ```

3. Run BSI-microbiome comparison:
   ```R
   source("R_scripts/bsi_microbiome_comparison.R")
   ```

4. Run machine learning analysis:
   ```bash
   python python_scripts/rf_shap_LMM_analysis.py
   ```

5. Generate t-SNE visualizations:
   ```bash
   python python_scripts/tSNE_plot.py
   ```

6. Run additional analyses as needed:
   ```R
   source("R_scripts/genus_analysis.R")
   ```
   
   ```bash
   python python_scripts/feature_selection_rf.py
   python python_scripts/zero-inflated-glmm.py
   ```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

dbhaslam@interface-labs.com  
david.haslam@cchmc.org
