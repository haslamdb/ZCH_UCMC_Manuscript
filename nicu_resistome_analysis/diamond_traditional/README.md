# DIAMOND Traditional Pipeline Analysis - NICU Resistome

**Date Moved:** November 13, 2025
**Original Location:** `/home/david/projects/benchmark_biogpu/`

This directory contains the complete DIAMOND-based resistome analysis for the NICU study (641 samples, 360 UCMC + 281 ZCH).

## Overview

DIAMOND blastx was used to align NICU infant microbiome samples against the AMR+Stress protein database (~14,280 genes). This analysis complements the BioGPU analysis and provides traditional bioinformatics pipeline results for comparison.

## Directory Structure

```
diamond_traditional/
├── README.md                           # This file
├── DIAMOND_ANALYSIS_SUMMARY.md         # Comprehensive analysis summary
├── scripts/                            # Analysis scripts (8 scripts)
│   ├── create_gene_annotation_mapping.py
│   ├── combine_diamond_results.py
│   ├── create_metadata_from_samples.py
│   ├── analyze_diamond_resistome.py
│   ├── differential_abundance_analysis.py
│   ├── antibiotic_resistance_correlations.py
│   ├── gene_antibiotic_correlations.py
│   └── plot_antibiotic_exposure_by_location.py
├── data/                               # Processed data (6 files)
│   ├── gene_annotations.tsv
│   ├── diamond_amr_combined.tsv        # 1.68M entries, 6,408 genes
│   ├── diamond_metadata.tsv
│   ├── diamond_metadata_updated.tsv    # Corrected metadata (641 samples)
│   ├── master_metadata_fixed.tsv
│   └── originally_missing_samples_antibiotic_status.tsv
├── results/                            # Analysis results (5 subdirectories)
│   ├── diamond_resistome/              # PCA, diversity, gene matrix
│   ├── differential_abundance/         # UCMC vs ZCH, body sites, weeks
│   ├── antibiotic_correlations/        # Total & specific antibiotic correlations
│   ├── gene_antibiotic_correlations/   # 19,106 significant gene×antibiotic pairs
│   └── antibiotic_exposure_comparison/ # Statistics for exposure differences
└── figures/                            # Publication-ready figures (3 subdirectories)
    ├── diamond_resistome/              # PCA, diversity, clustering plots
    ├── antibiotic_correlations/        # Correlation scatter plots
    └── antibiotic_exposure_comparison/ # 13 antibiotic comparison plots

Total size: ~40 MB
```

## Key Analyses

1. **Exploratory Analysis** (`analyze_diamond_resistome.py`)
   - PCA, diversity metrics, hierarchical clustering
   - 641 samples × 6,408 genes

2. **Differential Abundance** (`differential_abundance_analysis.py`)
   - UCMC vs ZCH: 2,666 significant genes
   - Body sites: 3,399 significant genes
   - Week 1 vs Week 3: 3,262 significant genes

3. **Antibiotic Correlations** (`antibiotic_resistance_correlations.py`)
   - Total antibiotic burden vs total AMR burden
   - Key finding: Negative correlation overall (rho = -0.150, p < 0.001)
   - Week 1: Strong negative correlation (rho = -0.322, p < 0.0001)

4. **Gene × Antibiotic Correlations** (`gene_antibiotic_correlations.py`)
   - 19,106 significant correlations (p < 0.05, uncorrected)
   - All 6,408 genes tested against 13 antibiotics

5. **Antibiotic Exposure Comparison** (`plot_antibiotic_exposure_by_location.py`)
   - 23 significant differences between UCMC and ZCH
   - Dramatic differences in prescribing practices

## Key Findings

- **ZCH uses more:** Penicillin, Moxalactam, Cefotaxime, Fluconazole, Meropenem
- **UCMC uses more:** Gentamicin, Nafcillin
- **Tetracycline resistance genes** highly enriched in UCMC vs ZCH
- **Quinolone resistance** increases 3-fold from Week 1 to Week 3
- **Negative antibiotic-AMR correlation** in early life (Week 1) suggests antibiotics reduce bacterial load

## Raw DIAMOND Results

Raw DIAMOND alignment results are stored in:
`/home/david/projects/benchmark_biogpu/results/traditional/`

- 642 sample directories (N01_1_2, ZJH_N01_1_2, etc.)
- Each contains: abundance tables, DIAMOND output, statistics, timing

## Related Analyses

- **BioGPU Analysis:** `../` (parent directory)
- **Original Metadata:** `/home/david/projects/ZCH_UCMC_Manuscript/metadata/`
- **Pipeline Comparison:** `/home/david/projects/benchmark_biogpu/` (original location)

## Running the Analysis

All scripts use absolute paths pointing to:
- Raw results: `/home/david/projects/benchmark_biogpu/results/traditional/`
- Metadata: This directory (`data/`)

To rerun any analysis, execute from the scripts directory:
```bash
cd scripts
python3 <script_name>.py
```

## Notes

- Sample naming: UCMC = N##_W_S, ZCH = ZJH_N##_W_S (zero-padded 01-09)
- Week encoding: 1 = Week.1, 2 = Week.3
- Site encoding: 2 = Axilla, 3 = Groin, 4 = Stool
- Antibiotic exposure: Cumulative (Week 1 = w1, Week 3 = w1 + w2)

## Documentation

See `DIAMOND_ANALYSIS_SUMMARY.md` for comprehensive documentation of all analyses, methods, and findings.

---

**Last Updated:** November 13, 2025
