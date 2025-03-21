# ZCH_UCMC_Manuscript

## Repository for NICU Microbiome Analysis

This repository contains the code and documentation for microbiome data analysis comparing samples from NICUs in Cincinnati and Hangzhou. Bloodstream microbiology data is compared to microbiome composition. Microbiome analysis examines bacterial composition, diversity metrics, and associations with factors such as antibiotics exposure, gestational age, and sample collection site/timing.

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
│   ├── nicu_utils.R           # Utility functions for data analysis and visualization
│   ├── import_initial_species_analysis.R  # Primary species-level analysis script
│   ├── bsi_microbiome_comparison.R        # Correlation between BSI microbes and microbiome
│   └── genus_analysis.R       # Genus-level taxonomic analysis
├── python_scripts/            # Python analysis scripts
│   ├── microbiome_transform.py         # Transformation functions for microbiome data
│   ├── rf_shap_LMM_analysis.py         # Random Forest and SHAP analysis with mixed effects models
│   ├── microbiome_shap_analysis_Kfold.py # K-fold cross-validation version of SHAP analysis
│   ├── feature_selection_rf.py         # Random forest for feature selection
│   └── shap_feature_selection.py       # SHAP-based feature importance analysis
├── bash_scripts/              # Bash processing scripts
│   └── process_reads.sh       # Process raw reads with Kraken2 and Bracken
└── results/                   # Analysis results directory
    ├── figures/               # Generated figures
    └── tables/                # Generated data tables
```

## Data Processing Pipeline

The analysis pipeline consists of several stages:

1. **Raw Data Processing**: Process FASTQ files to filter host DNA and assign taxonomic classifications
2. **Taxonomic Classification**: Raw sequencing data is processed using Kraken2 to assign taxonomic classifications
3. **Data Import & Filtering**: Taxonomic data is imported, filtered, and normalized
4. **Diversity Analysis**: Alpha diversity metrics are calculated and compared across groups
5. **Differential Abundance Analysis**: Statistical testing to identify differentially abundant taxa
6. **Ordination Analysis**: PCA and other methods to visualize community similarities/differences
7. **BSI Correlation Analysis**: Comparison of bloodstream infection data with microbiome composition
8. **Machine Learning Analysis**: Random Forest and SHAP analysis to identify important clinical features
9. **Visualization**: Generation of figures for publication

## Script Descriptions

### R Scripts

#### `nicu_utils.R`

This script contains utility functions used throughout the analysis:

- Color palettes and visualization themes
- Statistical analysis functions for effect size calculations
- Data processing functions for normalization and filtering
- Custom plotting functions and parameters
- Helper functions for fold change calculations and transformations

#### `import_initial_species_analysis.R`

This script performs the core analysis of microbiome data from the NICU study at the species level:

- Imports sample metadata and taxonomic abundance data
- Filters and normalizes species counts
- Calculates diversity metrics (Shannon, Simpson, etc.)
- Performs statistical comparisons between groups (location, antibiotics, gestational age)
- Identifies differentially abundant species
- Creates publication-quality visualizations
- Runs ordination analysis to visualize community structure

#### `bsi_microbiome_comparison.R`

This script focuses on the relationship between bloodstream infections and microbiome composition:

- Merges BSI (bloodstream infection) data with microbiome data
- Calculates Bray-Curtis distances between samples
- Computes observed-to-expected ratios for BSI organisms
- Performs PCA analysis to visualize similarities/differences
- Creates heatmaps and boxplots of BSI-differential organisms
- Conducts statistical tests for BSI-microbiome associations

#### `genus_analysis.R`

This script performs genus-level analysis of the microbiome data:

- Filters and processes genus-level taxonomy data
- Runs generalized linear mixed-effects models (GLMM) for various comparisons
- Analyzes effects of antibiotics, gestational age, and maternal factors
- Creates genus-specific visualizations
- Performs PCA analysis for different sample types and conditions
- Calculates effect sizes and fold changes between groups

### Python Scripts

#### `microbiome_transform.py`

This module provides functions for transforming microbiome count data:

- CLR (Centered Log-Ratio) transformation for compositional data
- TSS (Total Sum Scaling) normalization
- VST (Variance-Stabilizing Transformation)
- Log transformation with pseudocount handling
- Rarefaction (subsampling) to normalize sequencing depth

#### `rf_shap_LMM_analysis.py`

This script performs advanced analysis combining machine learning and mixed-effects models:

- Automatically selects appropriate transformation for microbiome data
- Implements Random Forest regression for abundance prediction
- Calculates SHAP (SHapley Additive exPlanations) values to identify important features
- Fits mixed-effects models with subject-level random effects
- Creates visualizations of feature importance with SHAP plots
- Generates summary statistics and effect estimations

#### `microbiome_shap_analysis_Kfold.py`

An extension of the SHAP analysis with K-fold cross-validation:

- Implements K-fold cross-validation for more robust feature importance
- Aggregates SHAP values across folds to reduce overfitting
- Includes additional validation steps and model diagnostics

#### `feature_selection_rf.py`

This script uses Random Forest for feature selection in microbiome data:

- Computes Bray-Curtis distances between samples
- Uses feature importance from Random Forest to identify key predictors
- Creates visualizations of important features

#### `shap_feature_selection.py`

This script focuses on SHAP-based feature importance for microbiome differences:

- Implements SHAP analysis for feature importance
- Creates summary plots and dependence plots for key features
- Integrates with distance-based analyses

### Bash Scripts

#### `process_reads.sh`

This script processes the raw sequencing data:

- Takes raw FASTQ files as input
- Removes human DNA using kneaddata
- Performs taxonomic classification using Kraken2
- Estimates abundance with Bracken
- Organizes outputs into structured directories
- Uses the kraken_tools pipeline for consistent processing

## Dependencies

### R Packages
- vegan: For diversity analyses and ordination
- ggplot2: For creating visualizations
- reshape2: For data reshaping
- FactoMineR and factoextra: For PCA and visualization
- dplyr: For data manipulation
- sda: For effect size calculations
- lme4: For mixed-effects models
- NBZIMM: For zero-inflated negative binomial regression
- tidyverse: For data wrangling
- pheatmap: For heatmaps
- VennDiagram: For creating Venn diagrams
- EnhancedVolcano: For volcano plots

### Python Packages
- pandas: For data manipulation
- numpy: For numerical operations
- scikit-learn: For machine learning models
- shap: For model interpretability and feature importance
- statsmodels: For statistical modeling
- matplotlib and seaborn: For visualization
- scipy: For scientific computing and statistics

### External Tools
- Kraken2: For taxonomic classification
- Bracken: For abundance estimation
- kneaddata: For host DNA removal
- kraken_tools: For pipeline execution (https://github.com/haslamdb/kraken_tools)

## Usage Instructions

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

### Running the Analysis

1. Process raw reads (if needed):
   ```bash
   bash bash_scripts/process_reads.sh <input_dir> <output_dir> <metadata_file>
   ```

2. Run the species-level analysis:
   ```R
   source("R_scripts/import_initial_species_analysis.R")
   ```

3. Run the BSI-microbiome comparison:
   ```R
   source("R_scripts/bsi_microbiome_comparison.R")
   ```

4. Run the genus-level analysis:
   ```R
   source("R_scripts/genus_analysis.R")
   ```

5. Run the Python-based machine learning analysis:
   ```bash
   python python_scripts/rf_shap_LMM_analysis.py
   ```

6. Transform microbiome data using the transformation module:
   ```python
   import microbiome_transform as mt
   # Example usage
   transformed_data = mt.clr_transform(microbiome_df)
   ```

### Interpreting Results

The analysis generates several key results:

1. **Diversity Analysis**: Comparison of microbial diversity between different locations, antibiotic exposures, and gestational age groups
2. **Taxonomic Differences**: Identification of bacterial species that differ significantly between groups
3. **Community Structure**: PCA visualizations showing similarities and differences in microbial communities
4. **BSI Correlations**: Analysis of relationships between bloodstream infections and microbiome composition
5. **Effect Size Analysis**: Quantification of the magnitude of differences between comparison groups
6. **Feature Importance**: SHAP and Random Forest analyses highlighting key clinical predictors of microbiome composition
7. **Mixed Effects Models**: Statistical models accounting for subject-level random effects

## License

MIT License

Copyright (c) 2025 David Haslam

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Citation


## Contact

dbhaslam@interface-labs.com
david.haslam@cchmc.org
