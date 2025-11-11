#!/usr/bin/env python3
"""
Phase 2: Species-Level Resistance Prevalence
===========================================

Analyzes which species carry FQ resistance mutations and tests for prevalence differences
by location and body site.

Input:
    - data/fq_resistance_filtered.tsv

Outputs:
    - results/fq_resistance/species_prevalence/prevalence_by_location.tsv
    - results/fq_resistance/species_prevalence/prevalence_by_bodysite.tsv
    - results/fq_resistance/species_prevalence/prevalence_by_week.tsv
    - results/fq_resistance/species_prevalence/statistical_tests.tsv
    - figures/fq_resistance/species_prevalence_barplot.pdf
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Paths
DATA_DIR = "data"
RESULTS_DIR = "results/fq_resistance/species_prevalence"
FIGURES_DIR = "figures/fq_resistance"

# Create output directories
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

print("=" * 80)
print("Phase 2: Species-Level Resistance Prevalence")
print("=" * 80)

# ============================================================================
# Step 1: Load filtered FQ resistance data
# ============================================================================
print("\n[1/5] Loading filtered FQ resistance data...")

fq_data = pd.read_csv(os.path.join(DATA_DIR, "fq_resistance_filtered.tsv"), sep="\t")
print(f"  - Loaded {len(fq_data):,} resistance mutation observations")
print(f"  - {fq_data['sample_name'].nunique()} unique samples")
print(f"  - {fq_data['species'].nunique()} unique species")

# ============================================================================
# Step 2: Calculate per-sample species resistance status
# ============================================================================
print("\n[2/5] Calculating per-sample species resistance...")

# For each sample × species combination, summarize resistance
sample_species = fq_data.groupby(['sample_name', 'species', 'Location', 'BodySite', 'Week', 'SubjectID']).agg({
    'mutation_id': 'nunique',  # Number of unique mutations
    'total_resistant_frequency': 'mean',  # Mean mutation frequency
    'is_high_frequency': 'sum'  # Number of high-frequency mutations
}).reset_index()

sample_species.columns = [
    'sample_name', 'species', 'Location', 'BodySite', 'Week', 'SubjectID',
    'n_mutations', 'mean_mutation_freq', 'n_high_freq_mutations'
]

print(f"  - {len(sample_species):,} sample × species combinations")
print(f"  - Mean mutations per sample×species: {sample_species['n_mutations'].mean():.1f}")

# ============================================================================
# Step 3: Species prevalence by location
# ============================================================================
print("\n[3/5] Calculating species prevalence by location...")

# Get total number of samples per location
location_totals = fq_data.groupby('Location')['sample_name'].nunique().to_dict()
print(f"\n  Total samples per location:")
for loc, count in location_totals.items():
    print(f"    {loc}: {count}")

# Calculate prevalence by location
prevalence_by_location = sample_species.groupby(['species', 'Location']).agg({
    'sample_name': 'nunique',
    'n_mutations': 'mean',
    'mean_mutation_freq': 'mean'
}).reset_index()

prevalence_by_location.columns = ['species', 'Location', 'n_samples', 'mean_n_mutations', 'mean_mutation_freq']

# Add total samples per location and calculate percentage
prevalence_by_location['total_samples'] = prevalence_by_location['Location'].map(location_totals)
prevalence_by_location['prevalence_pct'] = (
    prevalence_by_location['n_samples'] / prevalence_by_location['total_samples'] * 100
)

# Pivot for easier comparison
prevalence_pivot_loc = prevalence_by_location.pivot(
    index='species',
    columns='Location',
    values='prevalence_pct'
).fillna(0)

print(f"\nSpecies prevalence by location (% of samples):")
print(prevalence_pivot_loc.to_string())

# Statistical tests (Fisher's exact test for each species)
print(f"\nTesting prevalence differences (Fisher's exact test)...")
test_results_location = []

for species in sample_species['species'].unique():
    # Create 2x2 contingency table
    ucmc_with = sample_species[(sample_species['species'] == species) &
                               (sample_species['Location'] == 'UCMC')]['sample_name'].nunique()
    ucmc_without = location_totals['UCMC'] - ucmc_with

    zch_with = sample_species[(sample_species['species'] == species) &
                              (sample_species['Location'] == 'ZCH')]['sample_name'].nunique()
    zch_without = location_totals['ZCH'] - zch_with

    contingency = [[ucmc_with, ucmc_without], [zch_with, zch_without]]

    # Fisher's exact test
    odds_ratio, p_value = stats.fisher_exact(contingency)

    test_results_location.append({
        'species': species,
        'UCMC_n': ucmc_with,
        'UCMC_total': location_totals['UCMC'],
        'UCMC_pct': (ucmc_with / location_totals['UCMC']) * 100,
        'ZCH_n': zch_with,
        'ZCH_total': location_totals['ZCH'],
        'ZCH_pct': (zch_with / location_totals['ZCH']) * 100,
        'odds_ratio': odds_ratio,
        'p_value': p_value
    })

test_results_location_df = pd.DataFrame(test_results_location)

# FDR correction
from statsmodels.stats.multitest import multipletests
test_results_location_df['fdr'] = multipletests(test_results_location_df['p_value'], method='fdr_bh')[1]

# Sort by p-value
test_results_location_df = test_results_location_df.sort_values('p_value')

print(f"\nTop species with location differences:")
print(test_results_location_df[['species', 'UCMC_pct', 'ZCH_pct', 'odds_ratio', 'p_value', 'fdr']].head(10).to_string(index=False))

# ============================================================================
# Step 4: Species prevalence by body site
# ============================================================================
print("\n[4/5] Calculating species prevalence by body site...")

# Get total number of samples per body site
bodysite_totals = fq_data.groupby('BodySite')['sample_name'].nunique().to_dict()
print(f"\n  Total samples per body site:")
for site, count in bodysite_totals.items():
    print(f"    {site}: {count}")

# Calculate prevalence by body site
prevalence_by_bodysite = sample_species.groupby(['species', 'BodySite']).agg({
    'sample_name': 'nunique',
    'n_mutations': 'mean',
    'mean_mutation_freq': 'mean'
}).reset_index()

prevalence_by_bodysite.columns = ['species', 'BodySite', 'n_samples', 'mean_n_mutations', 'mean_mutation_freq']

# Add total samples and calculate percentage
prevalence_by_bodysite['total_samples'] = prevalence_by_bodysite['BodySite'].map(bodysite_totals)
prevalence_by_bodysite['prevalence_pct'] = (
    prevalence_by_bodysite['n_samples'] / prevalence_by_bodysite['total_samples'] * 100
)

# Pivot for easier comparison
prevalence_pivot_site = prevalence_by_bodysite.pivot(
    index='species',
    columns='BodySite',
    values='prevalence_pct'
).fillna(0)

print(f"\nSpecies prevalence by body site (% of samples):")
print(prevalence_pivot_site.to_string())

# Chi-square test for each species
print(f"\nTesting body site differences (Chi-square test)...")
test_results_bodysite = []

for species in sample_species['species'].unique():
    # Count samples with this species per body site
    counts = []
    for site in ['Axilla', 'Groin', 'Stool']:
        with_species = sample_species[(sample_species['species'] == species) &
                                      (sample_species['BodySite'] == site)]['sample_name'].nunique()
        without_species = bodysite_totals[site] - with_species
        counts.append([with_species, without_species])

    # Chi-square test
    contingency = np.array(counts)
    chi2, p_value = stats.chi2_contingency(contingency)[:2]

    test_results_bodysite.append({
        'species': species,
        'Axilla_n': counts[0][0],
        'Axilla_pct': (counts[0][0] / bodysite_totals['Axilla']) * 100,
        'Groin_n': counts[1][0],
        'Groin_pct': (counts[1][0] / bodysite_totals['Groin']) * 100,
        'Stool_n': counts[2][0],
        'Stool_pct': (counts[2][0] / bodysite_totals['Stool']) * 100,
        'chi2': chi2,
        'p_value': p_value
    })

test_results_bodysite_df = pd.DataFrame(test_results_bodysite)

# FDR correction
test_results_bodysite_df['fdr'] = multipletests(test_results_bodysite_df['p_value'], method='fdr_bh')[1]

# Sort by p-value
test_results_bodysite_df = test_results_bodysite_df.sort_values('p_value')

print(f"\nTop species with body site differences:")
print(test_results_bodysite_df[['species', 'Axilla_pct', 'Groin_pct', 'Stool_pct', 'chi2', 'p_value', 'fdr']].head(10).to_string(index=False))

# ============================================================================
# Step 5: Create visualization
# ============================================================================
print("\n[5/5] Creating visualizations...")

# Figure: Species prevalence by location and body site
fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# Panel A: By location
ax = axes[0]
prevalence_plot_loc = prevalence_by_location.pivot_table(
    index='species',
    columns='Location',
    values='prevalence_pct',
    fill_value=0
)
prevalence_plot_loc = prevalence_plot_loc.sort_values('UCMC', ascending=False).head(10)

x = np.arange(len(prevalence_plot_loc))
width = 0.35

ax.bar(x - width/2, prevalence_plot_loc['UCMC'], width, label='UCMC', alpha=0.8)
ax.bar(x + width/2, prevalence_plot_loc['ZCH'], width, label='ZCH', alpha=0.8)

ax.set_xlabel('Species', fontsize=12, fontweight='bold')
ax.set_ylabel('Prevalence (% of samples)', fontsize=12, fontweight='bold')
ax.set_title('FQ Resistance Prevalence by Location', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(prevalence_plot_loc.index, rotation=45, ha='right')
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Panel B: By body site
ax = axes[1]
prevalence_plot_site = prevalence_by_bodysite.pivot_table(
    index='species',
    columns='BodySite',
    values='prevalence_pct',
    fill_value=0
)
prevalence_plot_site = prevalence_plot_site.loc[prevalence_plot_loc.index]  # Same order as panel A

x = np.arange(len(prevalence_plot_site))
width = 0.25

ax.bar(x - width, prevalence_plot_site['Axilla'], width, label='Axilla', alpha=0.8)
ax.bar(x, prevalence_plot_site['Groin'], width, label='Groin', alpha=0.8)
ax.bar(x + width, prevalence_plot_site['Stool'], width, label='Stool', alpha=0.8)

ax.set_xlabel('Species', fontsize=12, fontweight='bold')
ax.set_ylabel('Prevalence (% of samples)', fontsize=12, fontweight='bold')
ax.set_title('FQ Resistance Prevalence by Body Site', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(prevalence_plot_site.index, rotation=45, ha='right')
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "species_prevalence_barplot.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved figure: {FIGURES_DIR}/species_prevalence_barplot.pdf")

# ============================================================================
# Save outputs
# ============================================================================
print("\n" + "=" * 80)
print("Saving outputs...")

# Save prevalence tables
prevalence_by_location.to_csv(
    os.path.join(RESULTS_DIR, "prevalence_by_location.tsv"),
    sep='\t',
    index=False
)
print(f"✓ Saved: {RESULTS_DIR}/prevalence_by_location.tsv")

prevalence_by_bodysite.to_csv(
    os.path.join(RESULTS_DIR, "prevalence_by_bodysite.tsv"),
    sep='\t',
    index=False
)
print(f"✓ Saved: {RESULTS_DIR}/prevalence_by_bodysite.tsv")

# Save statistical test results
test_results_location_df.to_csv(
    os.path.join(RESULTS_DIR, "location_statistical_tests.tsv"),
    sep='\t',
    index=False
)
print(f"✓ Saved: {RESULTS_DIR}/location_statistical_tests.tsv")

test_results_bodysite_df.to_csv(
    os.path.join(RESULTS_DIR, "bodysite_statistical_tests.tsv"),
    sep='\t',
    index=False
)
print(f"✓ Saved: {RESULTS_DIR}/bodysite_statistical_tests.tsv")

# Save sample-level species resistance data
sample_species.to_csv(
    os.path.join(RESULTS_DIR, "sample_species_resistance.tsv"),
    sep='\t',
    index=False
)
print(f"✓ Saved: {RESULTS_DIR}/sample_species_resistance.tsv")

print("\n" + "=" * 80)
print("Phase 2 complete!")
print("=" * 80)

# Print key findings
print("\nKey Findings:")
print(f"  • {len(sample_species['species'].unique())} species with FQ resistance detected")
print(f"  • Most prevalent: {test_results_location_df.iloc[0]['species']} "
      f"(UCMC: {test_results_location_df.iloc[0]['UCMC_pct']:.1f}%, "
      f"ZCH: {test_results_location_df.iloc[0]['ZCH_pct']:.1f}%)")

n_sig_location = (test_results_location_df['fdr'] < 0.05).sum()
if n_sig_location > 0:
    print(f"  • {n_sig_location} species with significant location differences (FDR < 0.05)")
else:
    print(f"  • No species with significant location differences (FDR < 0.05)")

n_sig_bodysite = (test_results_bodysite_df['fdr'] < 0.05).sum()
if n_sig_bodysite > 0:
    print(f"  • {n_sig_bodysite} species with significant body site differences (FDR < 0.05)")
else:
    print(f"  • No species with significant body site differences (FDR < 0.05)")

print("\nNext step: Run scripts/python/12_fq_mutation_analysis.py")
