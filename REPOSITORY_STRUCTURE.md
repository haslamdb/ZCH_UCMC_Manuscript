# Repository Structure

This repository contains code and analysis for the ZCH-UCMC NICU microbiome manuscript comparing neonatal intensive care units in Hangzhou (ZCH) and Cincinnati (UCMC).

## Directory Organization

```
ZCH_UCMC_Manuscript/
│
├── scripts/                      # All analysis scripts
│   ├── python/                   # Python scripts
│   │   ├── analysis/            # Statistical analysis and modeling
│   │   │   ├── feature_selection_rf.py
│   │   │   ├── microbiome_permanova.py
│   │   │   ├── microbiome_shap_analysis_Kfold.py
│   │   │   ├── rf_shap_LMM_analysis.py
│   │   │   ├── shap_feature_importance_revised.py
│   │   │   ├── zero-inflated-glmm.py
│   │   │   └── resistome/      # AMR-specific analyses
│   │   ├── visualization/      # Plotting and figure generation
│   │   │   ├── bsi_pedictors_viz.py
│   │   │   ├── network_plot_after_variance_partitioning.py
│   │   │   ├── plot_antibiotic_comparison.py
│   │   │   ├── tsne_multilayer_visualization.py
│   │   │   └── tSNE_plot.py
│   │   └── utils/              # Helper functions and utilities
│   │       ├── microbiome_transform.py
│   │       ├── parse_antibiotic_data.py
│   │       └── rename_antibiotic_columns.py
│   ├── R/                      # R analysis scripts
│   │   ├── shannon_diversity_figure.R
│   │   └── [other R scripts]
│   └── bash/                   # Shell scripts for automation
│
├── data/                        # Input data
│   ├── raw/                   # Original data files (gitignored)
│   ├── processed/             # Processed/transformed data
│   └── metadata/              # Sample metadata and clinical data
│
├── results/                     # Analysis outputs
│   ├── figures/               # Generated plots and visualizations
│   ├── tables/                # Summary statistics and CSVs
│   ├── models/                # Model outputs and coefficients
│   └── resistome/             # AMR analysis results
│
├── manuscript/                  # Paper-related files
│   ├── figures/               # Final figures for publication
│   └── supplementary/         # Supplementary materials
│
├── documentation/               # Project documentation
│   └── methods/               # Detailed method descriptions
│
└── notebooks/                   # Jupyter notebooks (if any)
```

## Key Analysis Scripts

### Python Analysis Scripts (`scripts/python/analysis/`)
- **rf_shap_LMM_analysis.py**: Combined Random Forest SHAP and Linear Mixed Model analysis
- **microbiome_shap_analysis_Kfold.py**: K-fold cross-validation with SHAP analysis
- **feature_selection_rf.py**: Random Forest feature selection
- **microbiome_permanova.py**: PERMANOVA analysis for microbiome data
- **zero-inflated-glmm.py**: Zero-inflated generalized linear mixed models

### Python Visualization Scripts (`scripts/python/visualization/`)
- **bsi_pedictors_viz.py**: BSI pathogen predictor visualizations
- **tsne_multilayer_visualization.py**: t-SNE plots with multiple layers
- **plot_antibiotic_comparison.py**: Antibiotic exposure comparisons
- **network_plot_after_variance_partitioning.py**: Network visualization

### R Scripts (`scripts/R/`)
- **shannon_diversity_figure.R**: Shannon diversity analysis and plotting
- Additional R scripts for statistical analysis

## Data Organization

### Raw Data (`data/raw/`)
- Original sequencing data and clinical metadata (gitignored)
- Should not be modified directly

### Processed Data (`data/processed/`)
- Transformed microbiome abundance data
- Feature matrices ready for analysis

### Metadata (`data/metadata/`)
- Sample information
- Clinical variables
- Antibiotic exposure data

## Results Organization

### Figures (`results/figures/`)
- All generated plots and visualizations
- Both exploratory and publication-ready figures

### Tables (`results/tables/`)
- CSV files with statistical results
- Summary statistics
- Model coefficients

### Models (`results/models/`)
- Saved model objects
- Cross-validation results
- Feature importance rankings

### Resistome Analysis (`results/resistome/`)
- AMR gene analysis results
- DIAMOND alignment outputs
- Resistance gene abundance tables

## Running Analyses

### Main Analysis Pipeline

1. **Data Processing**:
   ```bash
   python scripts/python/utils/microbiome_transform.py
   ```

2. **Statistical Analysis**:
   ```bash
   python scripts/python/analysis/rf_shap_LMM_analysis.py
   python scripts/python/analysis/microbiome_permanova.py
   ```

3. **Visualization**:
   ```bash
   python scripts/python/visualization/bsi_pedictors_viz.py
   Rscript scripts/R/shannon_diversity_figure.R
   ```

### Key Dependencies

- Python 3.8+
  - pandas, numpy, scikit-learn
  - shap, matplotlib, seaborn
  - statsmodels, scipy
- R 4.0+
  - ggplot2, dplyr, tidyr
  - vegan, lme4

## Git Workflow

The `.gitignore` file is configured to exclude:
- Raw data files (`data/raw/`)
- Generated figures (`*.pdf`, `*.png`)
- Large data files (`*.tsv`)
- Python cache files (`__pycache__/`, `*.pyc`)
- Results directories (tracked structure only)

## Legacy Directories

The following directories contain old analyses and will be removed in future cleanup:
- `python_scripts/` (contents moved to `scripts/python/`)
- `R_scripts/` (contents moved to `scripts/R/`)
- `results-LMM/` (contents moved to `results/models/`)
- `revision_figures/` (contents moved to `results/figures/`)

## Contact

For questions about this repository or the analyses, please contact the repository maintainers.