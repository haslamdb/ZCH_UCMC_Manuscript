# ZCH_UCMC_Manuscript

## Repository for NICU Microbiome Analysis

This repository contains the code and documentation for microbiome data analysis comparing samples from NICUs in Cincinnati and Hangzhou. Bloodstream microbiology is compared to microbiome composition. The analysis examines bacterial composition, diversity, and associations with factors such as antibiotics exposure, gestational age, and sample collection site/timing. 

## Repository Structure

```
ZCH_UCMC_Manuscript/
├── README.md                  # This file
├── metadata/                          # Metadata directory
│   ├── AllNICUSampleKeyRevised*.csv   # Sample metadata
│   ├── BSIData*.csv                   # BSI count data
│   └── HumanReactiveKraken2.csv       # Human reactive species list
├── data/                        
│   └── Kraken2/               # Kraken2 taxonomic classifications
├── R_scripts/                                  # R analysis scripts
│   └── nicu_utils.R                            # Utilities for data analysis and visualization
│   └── import_initial_species_analysis.R       # Main analysis script
│   └── bsi_microbiome_comparison.R             # Correlation between BSI microbes with skin and gut microbiome
├── bash_scripts/              # Bash processing scripts
└── results/                   # Analysis results 
    ├── figures/               # Generated figures
    └── tables/                # Generated data tables
```

## Data Processing Pipeline

The analysis pipeline consists of several stages:

1. **Taxonomic Classification**: Raw sequencing data is processed using Kraken2 to assign taxonomic classifications
2. **Data Import & Filtering**: Taxonomic data is imported, filtered, and normalized
3. **Diversity Analysis**: Alpha diversity metrics are calculated and compared across groups
4. **Differential Abundance Analysis**: Statistical testing to identify differentially abundant taxa
5. **Ordination Analysis**: PCA and other methods to visualize community similarities/differences
6. **Visualization**: Generation of figures for publication

## Script Descriptions

### R Scripts

#### `import_initial_analysis.R`

This script performs the core analysis of microbiome data from the NICU study. Key functions:

- Imports sample metadata and taxonomic abundance data
- Filters and normalizes species counts
- Calculates diversity metrics (Shannon, Simpson, etc.)
- Performs statistical comparisons between groups (location, antibiotics, gestational age)
- Identifies differentially abundant species
- Creates publication-quality visualizations
- Runs ordination analysis to visualize community structure

Usage:
```R
# Set the project directory in the script, then run:
source("R_scripts/import_initial_analysis.R")
```

Key outputs:
- Diversity statistics by group
- Differentially abundant species tables
- PCA plots showing sample relationships
- Boxplots of diversity metrics and species abundances

### Bash Scripts

*Note: Details about humann3_tools scripts will be added here.*

## Dependencies

### R Packages
- vegan: For diversity analyses and ordination
- ggplot2: For creating visualizations
- reshape2: For data reshaping
- FactoMineR and factoextra: For PCA and visualization
- dplyr: For data manipulation
- sda: For effect size calculations

### Other Dependencies
- Kraken2 (for taxonomic classification)
- *Additional dependencies from humann3_tools will be listed here*

## Usage Instructions

### Setting Up the Repository

1. Clone this repository:
   ```bash
   git clone https://github.com/YOUR-USERNAME/ZCH_UCMC_Manuscript.git
   cd ZCH_UCMC_Manuscript
   ```

2. Create the required directory structure:
   ```bash
   mkdir -p data results/figures results/tables KrakenAlignments/Kraken2 R_scripts bash_scripts
   ```

3. Place your input data in the appropriate directories:
   - Sample metadata in `data/`
   - Kraken2 results in `KrakenAlignments/Kraken2/`

### Running the Analysis

1. Update file paths in the R script to match your system
2. Run the main analysis script:
   ```R
   source("R_scripts/import_initial_analysis.R")
   ```

### Interpreting Results

The analysis generates several key results:

1. **Diversity Analysis**: Comparison of microbial diversity between different locations, antibiotic exposures, and gestational age groups
2. **Taxonomic Differences**: Identification of bacterial species that differ significantly between groups
3. **Community Structure**: PCA visualizations showing similarities and differences in microbial communities
4. **Species-Specific Analysis**: Detailed examination of key species of interest

## License

*Add your license information here*

## Citation

If you use this code in your research, please cite:

*Add citation information when available*

## Contact

*Add your contact information here*
