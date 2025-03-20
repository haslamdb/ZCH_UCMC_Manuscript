# ZCH_UCMC_Manuscript

## Repository for NICU Microbiome Analysis

This repository contains the code and documentation for microbiome data analysis comparing samples from NICUs in Cincinnati and Hangzhou. The analysis examines bacterial composition, diversity metrics, and associations with factors such as antibiotics exposure, gestational age, and sample collection site/timing. Bloodstream microbiology data is compared to microbiome composition.

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
8. **Visualization**: Generation of figures for publication

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

### Interpreting Results

The analysis generates several key results:

1. **Diversity Analysis**: Comparison of microbial diversity between different locations, antibiotic exposures, and gestational age groups
2. **Taxonomic Differences**: Identification of bacterial species that differ significantly between groups
3. **Community Structure**: PCA visualizations showing similarities and differences in microbial communities
4. **BSI Correlations**: Analysis of relationships between bloodstream infections and microbiome composition
5. **Effect Size Analysis**: Quantification of the magnitude of differences between comparison groups

## License

*Add your license information here*

## Citation

If you use this code in your research, please cite:

*Add citation information when available*

## Contact

*Add your contact information here*
