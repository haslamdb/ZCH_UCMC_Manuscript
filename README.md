# Kraken Tools

A comprehensive Python package for analyzing taxonomic profiles from Kraken2 and Bracken, with tools for data processing, statistical analysis, and visualization.

## Overview

Kraken Tools provides an end-to-end solution for microbiome analysis:

1. **Preprocessing**: Run quality control and host depletion with KneadData
2. **Taxonomic Classification**: Process sequences through Kraken2 and estimate abundances with Bracken
3. **Data Processing**: Merge, normalize, and filter taxonomic abundance data
4. **Exploratory Analysis**: Alpha/beta diversity, taxonomic visualization, PCA/PCoA analysis
5. **Differential Analysis**: Multiple methods for identifying significant differences between groups
   - PERMANOVA analysis for community-level differences
   - ALDEx2, ANCOM, ANCOM-BC for taxon-level differences
   - Feature selection with Random Forest
6. **Advanced Modeling**:
   - GLMM (Generalized Linear Mixed Models) for complex experimental designs
   - Random Forest with SHAP values for feature importance
   - Linear Mixed Models for specific taxa

## Installation

### Conda Environment Setup (Recommended)

```bash
# Clone the repository
git clone https://github.com/haslamdb/kraken_tools.git
cd kraken_tools

# Create a Conda environment
conda create -n kraken_tools python=3.12 -y
conda activate kraken_tools

# Install dependencies
conda install -c conda-forge -c bioconda pandas numpy scipy scikit-bio scikit-learn scikit-posthocs statsmodels matplotlib seaborn matplotlib-venn tqdm psutil shap
conda install -c bioconda kraken2 bracken kneaddata

# Install the package
pip install -e .
```

### Standard Installation

```bash
# Install directly using pip (once published to PyPI)
pip install kraken-tools

# Note: External tools (Kraken2, Bracken, KneadData) must be installed separately
```

## Workflow Overview

![Workflow Diagram](docs/workflow-diagram.png)

## Command-Line Usage

### 1. Full Pipeline (Raw Reads to Analysis)

Run the complete analysis pipeline in one command:

```bash
kraken-tools full-pipeline \
    --input-fastq reads_1.fastq.gz reads_2.fastq.gz \
    --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --sample-key metadata.csv \
    --output-dir results/ \
    --group-col "Group" \
    --min-abundance 0.01 \
    --min-prevalence 0.1 \
    --threads 8
```

### 2. Step by Step Analysis

#### Step 1: Preprocessing with KneadData

Quality control and host sequence removal:

```bash
kraken-tools preprocess \
    --input-fastq reads_1.fastq.gz reads_2.fastq.gz \
    --paired \
    --kneaddata-dbs /path/to/kneaddata_db \
    --output-dir results/preprocessed/ \
    --threads 8
```

#### Step 2: Taxonomic Classification (Kraken2 + Bracken)

Classify sequences and estimate abundances:

```bash
kraken-tools classify \
    --input-fastq clean_reads.fastq \
    --kraken-db /path/to/kraken_db \
    --bracken-db /path/to/kraken_db/database150mers.kmer_distrib \
    --output-dir results/taxonomy/ \
    --taxonomic-level S \
    --threads 8
```

#### Step 3: Data Processing

Merge, normalize, and filter abundance data:

```bash
kraken-tools process \
    --kreport-dir kraken_reports/ \
    --bracken-dir bracken_files/ \
    --sample-key metadata.csv \
    --output-dir results/processed/ \
    --min-abundance 0.01 \
    --min-prevalence 0.1
```

#### Step 4: Exploratory Analysis

Basic taxonomic and diversity analysis:

```bash
kraken-tools analyze \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/analysis/ \
    --group-col "Group"
```

#### Step 5: Differential Abundance Testing

Multiple methods for identifying significant taxa:

```bash
kraken-tools diff-abundance \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/diff_abundance/ \
    --group-col "Group" \
    --methods aldex2,ancom,ancom-bc
```

#### Step 6: PERMANOVA Analysis

Test for overall community differences:

```bash
kraken-tools permanova \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/permanova/ \
    --group-col "Group" \
    --categorical-vars "Treatment,TimePoint" \
    --distance-metric "bray"
```

#### Step 7: Feature Selection with Random Forest

Identify important variables driving microbiome differences:

```bash
kraken-tools feature-selection \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/feature_selection/ \
    --predictors "Treatment,TimePoint,Subject,Age" \
    --n-estimators 100
```

#### Step 8: GLMM Analysis

Complex modeling with mixed effects:

```bash
kraken-tools glmm \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/glmm/ \
    --formula "Count ~ Group + (1|Subject)" \
    --model negbin
```

#### Step 9: Random Forest with SHAP Analysis

Detailed feature importance analysis for specific taxa:

```bash
kraken-tools rf-shap \
    --abundance-file processed_abundance.tsv \
    --sample-key metadata.csv \
    --output-dir results/rf_shap/ \
    --target-taxa "Bacteroides.fragilis,Escherichia.coli" \
    --predictors "Treatment,TimePoint,Age" \
    --random-effects "Subject"
```

## Sample Key Format

The sample key CSV file should contain:

- A column with sample identifiers matching the file names in the input directories
- Additional columns for grouping and metadata

**Example**:
```csv
SampleName,Group,Treatment,TimePoint,Subject
sample1,Control,Placebo,Day0,Subject1
sample2,Treatment,Drug,Day0,Subject2
sample3,Control,Placebo,Day7,Subject1
sample4,Treatment,Drug,Day7,Subject2
```

## Output Structure

The output directory will have the following structure:

```
output_dir/
├── PreprocessedData/            # KneadData results
│   └── kneaddata_output/        # Clean reads after host removal
├── TaxonomyData/                # Taxonomic classification results
│   ├── kraken_reports/          # Kraken2 reports
│   └── bracken_output/          # Bracken abundance files
├── ProcessedData/               # Processed abundance data
│   ├── KrakenProcessed/         # Processed Kraken reports
│   └── BrackenProcessed/        # Processed Bracken files
├── ExploratoryAnalysis/         # Basic analysis results
│   ├── taxonomy_heatmap.svg     # Taxonomic heatmap
│   ├── taxonomy_pca.svg         # PCA/PCoA plots
│   └── diversity/               # Alpha/beta diversity
├── DifferentialAbundance/       # Differential abundance results
│   ├── aldex2_results.csv       # ALDEx2 results
│   ├── ancom_results.csv        # ANCOM results
│   └── ancom_bc_results.csv     # ANCOM-BC results
├── PERMANOVA/                   # PERMANOVA analysis
│   ├── permanova_results.csv    # Statistical results
│   └── pcoa_plots/              # PCoA visualization
├── FeatureSelection/            # Random Forest results
│   ├── feature_importance.pdf   # Feature importance plots
│   └── feature_importance.csv   # Feature ranking data
├── GLMM/                        # GLMM analysis results
│   ├── glmm_results.csv         # Model results
│   └── coefficient_plots/       # Coefficient visualizations
└── RF_SHAP/                     # Random Forest with SHAP
    ├── shap_summary/            # SHAP summary plots
    └── mixed_models/            # Linear mixed model results
```

## Python API Usage

For more flexibility, you can use the Python API:

```python
from kraken_tools import run_full_pipeline

# Run complete analysis
abundance_file, success = run_full_pipeline(
    sample_key="metadata.csv",
    kreport_dir="kraken_reports/",
    bracken_dir="bracken_files/",
    output_dir="results/",
    group_col="Group",
    min_abundance=0.01,
    min_prevalence=0.1,
    log_file="kraken_analysis.log"
)

# Run differential abundance analysis
from kraken_tools.analysis.differential import run_differential_abundance_analysis

results = run_differential_abundance_analysis(
    abundance_df=abundance_df,
    metadata_df=metadata_df,
    output_dir="diff_abundance/",
    group_col="Group",
    methods=["aldex2", "ancom-bc"],
    logger=logger
)
```

## Detailed Documentation

For more detailed information on each analysis step, see the individual documentation files:

- [Preprocessing with KneadData](docs/preprocessing.md)
- [Taxonomic Classification](docs/taxonomy_classification.md)
- [Data Processing](docs/data_processing.md)
- [Exploratory Analysis](docs/exploratory_analysis.md)
- [Differential Abundance Analysis](docs/differential_abundance.md)
- [PERMANOVA Analysis](docs/permanova.md)
- [Feature Selection](docs/feature_selection.md)
- [GLMM Analysis](docs/glmm_analysis.md)
- [Random Forest with SHAP](docs/rf_shap_analysis.md)

## Troubleshooting

See the [Troubleshooting Guide](docs/troubleshooting.md) for common issues and solutions.

## Citation

If you use Kraken Tools in your research, please cite:

- The original Kraken2 paper: Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol. 2019;20:257.
- Bracken: Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. PeerJ Comput Sci. 2017;3:e104.
- KneadData: McIver LJ, et al. bioBakery: a meta'omic analysis environment. Bioinformatics. 2018;34(7):1235-1237.
- This tool: Haslam D. (2025). Kraken Tools: A comprehensive framework for taxonomic analysis.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
