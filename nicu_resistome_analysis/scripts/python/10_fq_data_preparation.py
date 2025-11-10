#!/usr/bin/env python3
"""
FQ Resistance Data Preparation - Phase 1
=========================================

Merges FQ allele frequency data with sample metadata and filters to resistance mutations only.

Input:
    - data/nicu_fq_allele_frequencies.csv
    - results/qc/master_metadata_with_qc.tsv

Output:
    - data/fq_resistance_filtered.tsv
    - results/fq_resistance/qc/qc_summary.tsv
    - results/fq_resistance/qc/mutation_catalog.tsv
"""

import pandas as pd
import numpy as np
import os

# Paths
DATA_DIR = "data"
RESULTS_DIR = "results/fq_resistance"
QC_DIR = os.path.join(RESULTS_DIR, "qc")

print("=" * 80)
print("FQ Resistance Data Preparation - Phase 1")
print("=" * 80)

# ============================================================================
# Step 1: Load data
# ============================================================================
print("\n[1/5] Loading data...")

# Load FQ allele frequencies
print("  - Loading FQ allele frequency data...")
fq_data = pd.read_csv(os.path.join(DATA_DIR, "nicu_fq_allele_frequencies.csv"))
print(f"    Loaded {len(fq_data):,} rows")

# Load metadata
print("  - Loading metadata...")
metadata = pd.read_csv(
    os.path.join("results/qc", "master_metadata_with_qc.tsv"),
    sep="\t"
)
print(f"    Loaded metadata for {len(metadata):,} samples")

# ============================================================================
# Step 2: Filter to resistance mutations only
# ============================================================================
print("\n[2/5] Filtering to resistance mutations...")

print(f"  - Original data: {len(fq_data):,} rows")

# Filter 1: has_resistance_mutation == True
fq_resistance = fq_data[fq_data['has_resistance_mutation'] == True].copy()
print(f"  - After filtering to resistance mutations: {len(fq_resistance):,} rows")

# Filter 2: total_resistant_frequency > 0 (exclude wildtype-only)
fq_resistance = fq_resistance[fq_resistance['total_resistant_frequency'] > 0].copy()
print(f"  - After filtering frequency > 0: {len(fq_resistance):,} rows")

# ============================================================================
# Step 3: Merge with metadata
# ============================================================================
print("\n[3/5] Merging with metadata...")

# Merge on sample_name
fq_merged = fq_resistance.merge(
    metadata[['sample_name', 'PatientID', 'Location', 'BodySite', 'Week',
              'subject_complete', 'total_abx_days']],
    on='sample_name',
    how='left'
)

print(f"  - Merged data: {len(fq_merged):,} rows")
print(f"  - Samples with metadata: {fq_merged['Location'].notna().sum():,}")

# ============================================================================
# Step 4: Create derived variables
# ============================================================================
print("\n[4/5] Creating derived variables...")

# Create mutation ID
fq_merged['mutation_id'] = (
    fq_merged['species'].str.replace(' ', '_') + "_" +
    fq_merged['gene'] + "_" +
    fq_merged['position'].astype(str) + "_" +
    fq_merged['wildtype_aa'] + "->" +
    fq_merged['dominant_mutant_aa']
)

# Create high frequency indicator (≥75% mutant)
fq_merged['is_high_frequency'] = fq_merged['total_resistant_frequency'] >= 0.75

# Create binary presence indicator
fq_merged['mutation_present'] = fq_merged['total_resistant_frequency'] > 0.1  # >10% frequency

print(f"  - Created mutation_id column")
print(f"  - High frequency mutations (≥75%): {fq_merged['is_high_frequency'].sum():,}")
print(f"  - Mutations present (>10%): {fq_merged['mutation_present'].sum():,}")

# ============================================================================
# Step 5: QC and summary statistics
# ============================================================================
print("\n[5/5] Generating QC summaries...")

# QC summary
qc_summary = pd.DataFrame({
    'metric': [
        'Total samples in FQ data',
        'Samples with metadata',
        'Unique resistance mutations',
        'Total resistance observations',
        'Mean mutations per sample',
        'Samples from UCMC',
        'Samples from ZCH',
        'Axilla samples',
        'Groin samples',
        'Stool samples',
        'Week 1 samples',
        'Week 3 samples',
        'Complete subjects'
    ],
    'value': [
        fq_merged['sample_name'].nunique(),
        fq_merged['Location'].notna().sum(),
        fq_merged['mutation_id'].nunique(),
        len(fq_merged),
        fq_merged.groupby('sample_name').size().mean(),
        (fq_merged['Location'] == 'UCMC').sum(),
        (fq_merged['Location'] == 'ZCH').sum(),
        (fq_merged['BodySite'] == 'Axilla').sum(),
        (fq_merged['BodySite'] == 'Groin').sum(),
        (fq_merged['BodySite'] == 'Stool').sum(),
        (fq_merged['Week'] == 'Week 1').sum(),
        (fq_merged['Week'] == 'Week 3').sum(),
        fq_merged[fq_merged['subject_complete'] == True]['sample_name'].nunique()
    ]
})

print("\nQC Summary:")
print(qc_summary.to_string(index=False))

# Mutation catalog
mutation_catalog = fq_merged.groupby(['species', 'gene', 'position', 'wildtype_aa',
                                      'dominant_mutant_aa', 'mutation_id']).agg({
    'sample_name': 'nunique',
    'total_resistant_frequency': ['mean', 'median', 'std', 'min', 'max'],
    'total_depth': ['mean', 'median']
}).reset_index()

# Flatten column names
mutation_catalog.columns = [
    'species', 'gene', 'position', 'wildtype_aa', 'mutant_aa', 'mutation_id',
    'n_samples', 'mean_frequency', 'median_frequency', 'std_frequency',
    'min_frequency', 'max_frequency', 'mean_depth', 'median_depth'
]

# Sort by prevalence
mutation_catalog = mutation_catalog.sort_values('n_samples', ascending=False)

print(f"\nTop 10 most prevalent mutations:")
print(mutation_catalog[['mutation_id', 'n_samples', 'mean_frequency']].head(10).to_string(index=False))

# Species summary
species_summary = fq_merged.groupby('species').agg({
    'sample_name': 'nunique',
    'mutation_id': 'nunique',
    'total_resistant_frequency': 'mean'
}).reset_index()
species_summary.columns = ['species', 'n_samples', 'n_unique_mutations', 'mean_resistance_freq']
species_summary = species_summary.sort_values('n_samples', ascending=False)

print(f"\nSpecies summary:")
print(species_summary.to_string(index=False))

# ============================================================================
# Save outputs
# ============================================================================
print("\n" + "=" * 80)
print("Saving outputs...")

# Save filtered data
output_file = os.path.join(DATA_DIR, "fq_resistance_filtered.tsv")
fq_merged.to_csv(output_file, sep='\t', index=False)
print(f"✓ Saved filtered data: {output_file}")
print(f"  ({len(fq_merged):,} rows × {len(fq_merged.columns)} columns)")

# Save QC summary
qc_file = os.path.join(QC_DIR, "qc_summary.tsv")
qc_summary.to_csv(qc_file, sep='\t', index=False)
print(f"✓ Saved QC summary: {qc_file}")

# Save mutation catalog
catalog_file = os.path.join(QC_DIR, "mutation_catalog.tsv")
mutation_catalog.to_csv(catalog_file, sep='\t', index=False)
print(f"✓ Saved mutation catalog: {catalog_file}")
print(f"  ({len(mutation_catalog):,} unique mutations)")

# Save species summary
species_file = os.path.join(QC_DIR, "species_summary.tsv")
species_summary.to_csv(species_file, sep='\t', index=False)
print(f"✓ Saved species summary: {species_file}")

print("\n" + "=" * 80)
print("Phase 1 complete!")
print("=" * 80)
print("\nNext step: Run scripts/python/11_fq_species_prevalence.py")
