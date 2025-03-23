# Unified Configuration for Microbiome Analysis

This document provides an overview of the unified configuration system for your microbiome analysis pipeline.

## Overview

The unified configuration system consists of:

1. A central YAML configuration file (`config.yaml`)
2. Utility modules for loading the configuration in Python and R
3. Example scripts demonstrating how to use the configuration

This system ensures consistent parameters across all analyses, making your workflow more reproducible and easier to manage.

## Configuration File Structure

The `config.yaml` file contains standardized parameters for all analyses, organized into the following sections:

- **paths**: Project directory paths
- **input_files**: Paths to input data files
- **data_processing**: Parameters for data processing and filtering
- **statistics**: Statistical analysis parameters
- **mixed_effects**: Parameters for mixed-effects models
- **zero_inflated**: Parameters for zero-inflated models
- **visualization**: Visualization settings and styles
- **tsne**: t-SNE algorithm parameters
- **pca**: PCA parameters
- **database**: Reference database paths
- **key_organisms**: List of target microbes to analyze
- **metadata_features**: Metadata columns for modeling
- **analysis**: Analysis-specific parameters

## Using the Configuration

### In Python Scripts

```python
# Import the configuration loader
from utils.config import load_config

# Load the configuration
config = load_config()

# Access configuration values
input_file = config['input_files']['microbiome_counts']
results_dir = config['paths']['results_dir']
key_organisms = config['key_organisms']
```

### In R Scripts

```r
# Source the configuration loader
source("utils/config.R")

# Load the configuration
config <- load_config()

# Access configuration values
input_file <- config$input_files$microbiome_counts
results_dir <- config$paths$results_dir
key_organisms <- config$key_organisms
```

## Path Resolution

The configuration system supports variable references in paths using the `${variable_name}` syntax. For example:

```yaml
paths:
  project_dir: "."
  data_dir: "${project_dir}/data"
  results_dir: "${project_dir}/results"
```

This allows for flexible directory structures while maintaining relative paths.

## Customization

You can modify the `config.yaml` file to adjust parameters for your specific needs. Common customizations include:

- Updating file paths for your environment
- Adjusting statistical thresholds
- Modifying visualization parameters
- Adding or removing key organisms

## Directory Structure

```
project_root/
├── config.yaml                 # Unified configuration file
├── utils/
│   ├── config.py               # Python configuration loader
│   └── config.R                # R configuration loader
├── examples/
│   ├── python_example.py       # Example Python script
│   └── r_example.R             # Example R script
├── data/                       # Data directory
│   └── ...
└── results/                    # Results directory
    ├── figures/                # Generated figures
    ├── models/                 # Saved models
    └── reports/                # Analysis reports
```

## Benefits

Using the unified configuration system provides several benefits:

1. **Reproducibility**: Standardized parameters ensure consistent results
2. **Maintainability**: Central location for all parameters
3. **Flexibility**: Easy to adjust parameters without modifying scripts
4. **Documentation**: Configuration serves as documentation of analysis parameters
5. **Consistency**: Same parameters used across Python and R analyses

## Migration Guide

To migrate existing scripts to use the unified configuration:

1. Replace hardcoded parameters with references to the configuration
2. Ensure consistent parameter names across scripts
3. Update file paths to use the standardized directory structure
4. Add the configuration loading code at the beginning of each script
