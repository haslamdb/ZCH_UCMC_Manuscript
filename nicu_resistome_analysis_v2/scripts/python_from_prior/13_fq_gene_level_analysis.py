#!/usr/bin/env python3
"""
Phase 4: Gene-Level Resistance Burden
=====================================

Analyzes resistance at the gene level by creating composite scores and identifying
multi-gene resistance patterns.

Input:
    - data/fq_resistance_filtered.tsv

Outputs:
    - results/fq_resistance/gene_analysis/gene_level_scores.tsv
    - results/fq_resistance/gene_analysis/multi_gene_patterns.tsv
    - results/fq_resistance/gene_analysis/gene_comparison_by_location.tsv
    - results/fq_resistance/gene_analysis/gene_comparison_by_bodysite.tsv
    - figures/fq_resistance/gene_scores_by_location.pdf
    - figures/fq_resistance/mutation_burden_distribution.pdf
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os

# Paths
DATA_DIR = "data"
RESULTS_DIR = "results/fq_resistance/gene_analysis"
FIGURES_DIR = "figures/fq_resistance"

# Create output directories
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

print("=" * 80)
print("Phase 4: Gene-Level Resistance Burden")
print("=" * 80)

# ============================================================================
# Step 1: Load data
# ============================================================================
print("\n[1/5] Loading filtered FQ resistance data...")

fq_data = pd.read_csv(os.path.join(DATA_DIR, "fq_resistance_filtered.tsv"), sep="\t")
print(f"  - Loaded {len(fq_data):,} observations")

# ============================================================================
# Step 2: Calculate composite resistance scores per gene
# ============================================================================
print("\n[2/5] Calculating gene-level composite resistance scores...")

# For each sample × species × gene, sum all mutation frequencies
gene_scores = fq_data.groupby(['sample_name', 'species', 'gene', 'Location', 'BodySite', 'Week', 'SubjectID']).agg({
    'total_resistant_frequency': 'sum',  # Sum of all mutation frequencies
    'mutation_id': 'nunique',  # Number of unique mutations
    'is_high_frequency': 'sum'  # Number of high-frequency mutations
}).reset_index()

gene_scores.columns = [
    'sample_name', 'species', 'gene', 'Location', 'BodySite', 'Week', 'SubjectID',
    'composite_score', 'n_mutations', 'n_high_freq'
]

print(f"  - {len(gene_scores):,} sample × species × gene combinations")
print(f"\nComposite score distribution:")
print(gene_scores['composite_score'].describe())

print(f"\nGene distribution:")
print(gene_scores['gene'].value_counts())

# Save gene scores
gene_scores.to_csv(
    os.path.join(RESULTS_DIR, "gene_level_scores.tsv"),
    sep='\t',
    index=False
)
print(f"\n✓ Saved: {RESULTS_DIR}/gene_level_scores.tsv")

# ============================================================================
# Step 3: Multi-gene resistance patterns
# ============================================================================
print("\n[3/5] Analyzing multi-gene resistance patterns...")

# For each sample × species, identify which genes have resistance
sample_species_genes = gene_scores.groupby(['sample_name', 'species', 'Location', 'BodySite']).agg({
    'gene': lambda x: ','.join(sorted(x)),  # List of genes
    'composite_score': 'sum',  # Total resistance burden
    'n_mutations': 'sum'  # Total number of mutations
}).reset_index()

sample_species_genes.columns = [
    'sample_name', 'species', 'Location', 'BodySite',
    'genes_with_resistance', 'total_resistance_score', 'total_mutations'
]

# Count number of genes with resistance
sample_species_genes['n_genes_with_resistance'] = sample_species_genes['genes_with_resistance'].apply(
    lambda x: len(x.split(','))
)

print(f"\nMulti-gene resistance patterns:")
print(sample_species_genes['n_genes_with_resistance'].value_counts().sort_index())

# Analyze co-occurrence patterns
print(f"\nMost common gene combinations:")
gene_combo_counts = sample_species_genes.groupby(['species', 'genes_with_resistance']).size().reset_index(name='n_samples')
gene_combo_counts = gene_combo_counts.sort_values('n_samples', ascending=False)
print(gene_combo_counts.head(20).to_string(index=False))

# Save multi-gene patterns
sample_species_genes.to_csv(
    os.path.join(RESULTS_DIR, "multi_gene_patterns.tsv"),
    sep='\t',
    index=False
)
print(f"\n✓ Saved: {RESULTS_DIR}/multi_gene_patterns.tsv")

# ============================================================================
# Step 4: Gene-level comparisons by location and body site
# ============================================================================
print("\n[4/5] Comparing gene scores by location and body site...")

# Location comparison
print("\n--- Location Comparison ---")
location_results = []

for gene in gene_scores['gene'].unique():
    gene_data = gene_scores[gene_scores['gene'] == gene]

    ucmc_scores = gene_data[gene_data['Location'] == 'UCMC']['composite_score']
    zch_scores = gene_data[gene_data['Location'] == 'ZCH']['composite_score']

    if len(ucmc_scores) >= 5 and len(zch_scores) >= 5:
        # Mann-Whitney U test
        statistic, p_value = stats.mannwhitneyu(ucmc_scores, zch_scores, alternative='two-sided')

        location_results.append({
            'gene': gene,
            'UCMC_n': len(ucmc_scores),
            'UCMC_mean': ucmc_scores.mean(),
            'UCMC_median': ucmc_scores.median(),
            'ZCH_n': len(zch_scores),
            'ZCH_mean': zch_scores.mean(),
            'ZCH_median': zch_scores.median(),
            'statistic': statistic,
            'p_value': p_value
        })

location_df = pd.DataFrame(location_results)
location_df['fdr'] = multipletests(location_df['p_value'], method='fdr_bh')[1]
location_df = location_df.sort_values('p_value')

print(location_df.to_string(index=False))

location_df.to_csv(
    os.path.join(RESULTS_DIR, "gene_comparison_by_location.tsv"),
    sep='\t',
    index=False
)
print(f"\n✓ Saved: {RESULTS_DIR}/gene_comparison_by_location.tsv")

# Body site comparison
print("\n--- Body Site Comparison ---")
bodysite_results = []

for gene in gene_scores['gene'].unique():
    gene_data = gene_scores[gene_scores['gene'] == gene]

    axilla_scores = gene_data[gene_data['BodySite'] == 'Axilla']['composite_score']
    groin_scores = gene_data[gene_data['BodySite'] == 'Groin']['composite_score']
    stool_scores = gene_data[gene_data['BodySite'] == 'Stool']['composite_score']

    if len(axilla_scores) >= 5 and len(groin_scores) >= 5 and len(stool_scores) >= 5:
        # Kruskal-Wallis test
        statistic, p_value = stats.kruskal(axilla_scores, groin_scores, stool_scores)

        bodysite_results.append({
            'gene': gene,
            'Axilla_n': len(axilla_scores),
            'Axilla_mean': axilla_scores.mean(),
            'Groin_n': len(groin_scores),
            'Groin_mean': groin_scores.mean(),
            'Stool_n': len(stool_scores),
            'Stool_mean': stool_scores.mean(),
            'statistic': statistic,
            'p_value': p_value
        })

bodysite_df = pd.DataFrame(bodysite_results)
bodysite_df['fdr'] = multipletests(bodysite_df['p_value'], method='fdr_bh')[1]
bodysite_df = bodysite_df.sort_values('p_value')

print(bodysite_df.to_string(index=False))

bodysite_df.to_csv(
    os.path.join(RESULTS_DIR, "gene_comparison_by_bodysite.tsv"),
    sep='\t',
    index=False
)
print(f"\n✓ Saved: {RESULTS_DIR}/gene_comparison_by_bodysite.tsv")

# ============================================================================
# Step 5: Create visualizations
# ============================================================================
print("\n[5/5] Creating visualizations...")

# Figure 1: Gene scores by location
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

genes = sorted(gene_scores['gene'].unique())[:4]  # Top 4 genes

for i, gene in enumerate(genes):
    ax = axes[i]
    gene_data = gene_scores[gene_scores['gene'] == gene]

    # Boxplot by location
    data_to_plot = [
        gene_data[gene_data['Location'] == 'UCMC']['composite_score'],
        gene_data[gene_data['Location'] == 'ZCH']['composite_score']
    ]

    bp = ax.boxplot(data_to_plot, tick_labels=['UCMC', 'ZCH'], patch_artist=True)

    # Color boxes
    for patch, color in zip(bp['boxes'], ['lightblue', 'lightcoral']):
        patch.set_facecolor(color)

    ax.set_ylabel('Composite Resistance Score', fontsize=10)
    ax.set_title(f"{gene}", fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "gene_scores_by_location.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/gene_scores_by_location.pdf")

# Figure 2: Mutation burden distribution
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel A: Distribution of number of genes with resistance
ax = axes[0]
burden_counts = sample_species_genes['n_genes_with_resistance'].value_counts().sort_index()
ax.bar(burden_counts.index, burden_counts.values, color='steelblue', alpha=0.7)
ax.set_xlabel('Number of Genes with Resistance', fontsize=12, fontweight='bold')
ax.set_ylabel('Number of Sample×Species', fontsize=12, fontweight='bold')
ax.set_title('Multi-Gene Resistance Burden', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Panel B: Total resistance score distribution
ax = axes[1]
ax.hist(sample_species_genes['total_resistance_score'], bins=30, color='coral', alpha=0.7, edgecolor='black')
ax.set_xlabel('Total Resistance Score', fontsize=12, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax.set_title('Distribution of Total Resistance Scores', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "mutation_burden_distribution.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/mutation_burden_distribution.pdf")

print("\n" + "=" * 80)
print("Phase 4 complete!")
print("=" * 80)

print("\nKey Findings:")
print(f"  • {len(gene_scores):,} sample×species×gene combinations analyzed")
print(f"  • Mean composite score: {gene_scores['composite_score'].mean():.2f}")
print(f"  • {(sample_species_genes['n_genes_with_resistance'] > 1).sum()} sample×species with multi-gene resistance")

n_sig_loc = (location_df['fdr'] < 0.05).sum()
n_sig_site = (bodysite_df['fdr'] < 0.05).sum()
print(f"  • {n_sig_loc} genes differ by location (FDR < 0.05)")
print(f"  • {n_sig_site} genes differ by body site (FDR < 0.05)")

print("\nNext step: Run scripts/python/14_fq_antibiotic_correlations.py")
