#!/usr/bin/env python3
"""
Phase 3: Mutation-Level Analysis
================================

Compares specific FQ resistance mutation frequencies across groups.

Input:
    - data/fq_resistance_filtered.tsv

Outputs:
    - results/fq_resistance/mutation_analysis/top_mutations.tsv
    - results/fq_resistance/mutation_analysis/mutation_location_comparison.tsv
    - results/fq_resistance/mutation_analysis/mutation_bodysite_comparison.tsv
    - results/fq_resistance/mutation_analysis/mutation_longitudinal_changes.tsv
    - figures/fq_resistance/mutation_heatmap_by_location.pdf
    - figures/fq_resistance/mutation_frequency_boxplots.pdf
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
RESULTS_DIR = "results/fq_resistance/mutation_analysis"
FIGURES_DIR = "figures/fq_resistance"

# Create output directories
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

print("=" * 80)
print("Phase 3: Mutation-Level Analysis")
print("=" * 80)

# ============================================================================
# Step 1: Load data
# ============================================================================
print("\n[1/6] Loading filtered FQ resistance data...")

fq_data = pd.read_csv(os.path.join(DATA_DIR, "fq_resistance_filtered.tsv"), sep="\t")
print(f"  - Loaded {len(fq_data):,} observations")
print(f"  - {fq_data['mutation_id'].nunique()} unique mutations")

# ============================================================================
# Step 2: Identify top mutations overall
# ============================================================================
print("\n[2/6] Identifying top mutations...")

# Summarize each mutation
mutation_summary = fq_data.groupby(['species', 'gene', 'position', 'mutation_id']).agg({
    'sample_name': 'nunique',
    'total_resistant_frequency': ['mean', 'median', 'std'],
    'is_high_frequency': 'sum'
}).reset_index()

mutation_summary.columns = [
    'species', 'gene', 'position', 'mutation_id',
    'n_samples', 'mean_freq', 'median_freq', 'std_freq', 'n_high_freq'
]

# Calculate prevalence
total_samples = fq_data['sample_name'].nunique()
mutation_summary['prevalence_pct'] = (mutation_summary['n_samples'] / total_samples) * 100

# Sort by prevalence
mutation_summary = mutation_summary.sort_values('n_samples', ascending=False)

print(f"\nTop 20 most prevalent mutations:")
print(mutation_summary[['mutation_id', 'n_samples', 'prevalence_pct', 'mean_freq']].head(20).to_string(index=False))

# Save top mutations
mutation_summary.to_csv(
    os.path.join(RESULTS_DIR, "top_mutations.tsv"),
    sep='\t',
    index=False
)
print(f"\n✓ Saved: {RESULTS_DIR}/top_mutations.tsv")

# ============================================================================
# Step 3: Location comparison (UCMC vs ZCH)
# ============================================================================
print("\n[3/6] Comparing mutations by location...")

# Get mutations present in at least 10 samples total
common_mutations = mutation_summary[mutation_summary['n_samples'] >= 10]['mutation_id'].tolist()
print(f"  - Testing {len(common_mutations)} mutations (≥10 samples)")

location_results = []

for mutation in common_mutations:
    mutation_data = fq_data[fq_data['mutation_id'] == mutation]

    # UCMC data
    ucmc_data = mutation_data[mutation_data['Location'] == 'UCMC']['total_resistant_frequency']
    # ZCH data
    zch_data = mutation_data[mutation_data['Location'] == 'ZCH']['total_resistant_frequency']

    if len(ucmc_data) >= 3 and len(zch_data) >= 3:  # Need at least 3 samples per group
        # Mann-Whitney U test
        statistic, p_value = stats.mannwhitneyu(ucmc_data, zch_data, alternative='two-sided')

        # Prevalence comparison (Fisher's exact)
        ucmc_total = fq_data[fq_data['Location'] == 'UCMC']['sample_name'].nunique()
        zch_total = fq_data[fq_data['Location'] == 'ZCH']['sample_name'].nunique()
        ucmc_with = len(ucmc_data)
        zch_with = len(zch_data)
        ucmc_without = ucmc_total - ucmc_with
        zch_without = zch_total - zch_with

        contingency = [[ucmc_with, ucmc_without], [zch_with, zch_without]]
        odds_ratio, p_prevalence = stats.fisher_exact(contingency)

        location_results.append({
            'mutation_id': mutation,
            'species': mutation_data['species'].iloc[0],
            'gene': mutation_data['gene'].iloc[0],
            'UCMC_n': len(ucmc_data),
            'UCMC_mean_freq': ucmc_data.mean(),
            'UCMC_median_freq': ucmc_data.median(),
            'ZCH_n': len(zch_data),
            'ZCH_mean_freq': zch_data.mean(),
            'ZCH_median_freq': zch_data.median(),
            'fold_change': ucmc_data.mean() / zch_data.mean() if zch_data.mean() > 0 else np.nan,
            'mw_statistic': statistic,
            'mw_pvalue': p_value,
            'prevalence_odds_ratio': odds_ratio,
            'prevalence_pvalue': p_prevalence
        })

location_df = pd.DataFrame(location_results)

# FDR correction
location_df['freq_fdr'] = multipletests(location_df['mw_pvalue'], method='fdr_bh')[1]
location_df['prevalence_fdr'] = multipletests(location_df['prevalence_pvalue'], method='fdr_bh')[1]

# Sort by p-value
location_df = location_df.sort_values('mw_pvalue')

print(f"\nMutations with location differences (top 10 by p-value):")
print(location_df[['mutation_id', 'UCMC_mean_freq', 'ZCH_mean_freq', 'fold_change', 'mw_pvalue', 'freq_fdr']].head(10).to_string(index=False))

n_sig = (location_df['freq_fdr'] < 0.05).sum()
print(f"\n  • {n_sig} mutations with significant frequency differences (FDR < 0.05)")

# Save location comparison
location_df.to_csv(
    os.path.join(RESULTS_DIR, "mutation_location_comparison.tsv"),
    sep='\t',
    index=False
)
print(f"✓ Saved: {RESULTS_DIR}/mutation_location_comparison.tsv")

# ============================================================================
# Step 4: Body site comparison
# ============================================================================
print("\n[4/6] Comparing mutations by body site...")

bodysite_results = []

for mutation in common_mutations:
    mutation_data = fq_data[fq_data['mutation_id'] == mutation]

    # Data per body site
    axilla_data = mutation_data[mutation_data['BodySite'] == 'Axilla']['total_resistant_frequency']
    groin_data = mutation_data[mutation_data['BodySite'] == 'Groin']['total_resistant_frequency']
    stool_data = mutation_data[mutation_data['BodySite'] == 'Stool']['total_resistant_frequency']

    # Need at least 3 samples per group
    if len(axilla_data) >= 3 and len(groin_data) >= 3 and len(stool_data) >= 3:
        # Kruskal-Wallis test
        statistic, p_value = stats.kruskal(axilla_data, groin_data, stool_data)

        bodysite_results.append({
            'mutation_id': mutation,
            'species': mutation_data['species'].iloc[0],
            'gene': mutation_data['gene'].iloc[0],
            'Axilla_n': len(axilla_data),
            'Axilla_mean_freq': axilla_data.mean(),
            'Groin_n': len(groin_data),
            'Groin_mean_freq': groin_data.mean(),
            'Stool_n': len(stool_data),
            'Stool_mean_freq': stool_data.mean(),
            'kw_statistic': statistic,
            'p_value': p_value
        })

bodysite_df = pd.DataFrame(bodysite_results)

# FDR correction
bodysite_df['fdr'] = multipletests(bodysite_df['p_value'], method='fdr_bh')[1]

# Sort by p-value
bodysite_df = bodysite_df.sort_values('p_value')

print(f"\nMutations with body site differences (top 10 by p-value):")
print(bodysite_df[['mutation_id', 'Axilla_mean_freq', 'Groin_mean_freq', 'Stool_mean_freq', 'p_value', 'fdr']].head(10).to_string(index=False))

n_sig = (bodysite_df['fdr'] < 0.05).sum()
print(f"\n  • {n_sig} mutations with significant body site differences (FDR < 0.05)")

# Save body site comparison
bodysite_df.to_csv(
    os.path.join(RESULTS_DIR, "mutation_bodysite_comparison.tsv"),
    sep='\t',
    index=False
)
print(f"✓ Saved: {RESULTS_DIR}/mutation_bodysite_comparison.tsv")

# ============================================================================
# Step 5: Longitudinal changes (Week 1 → Week 3)
# ============================================================================
print("\n[5/6] Analyzing longitudinal changes...")

# Get complete subjects (have both Week 1 and Week 3)
complete_subjects = fq_data[fq_data['subject_complete'] == True]['SubjectID'].unique()
print(f"  - {len(complete_subjects)} complete subjects for paired analysis")

longitudinal_results = []

for mutation in common_mutations:
    mutation_data = fq_data[(fq_data['mutation_id'] == mutation) &
                            (fq_data['SubjectID'].isin(complete_subjects))]

    # Need to have paired samples
    week1_data = mutation_data[mutation_data['Week'] == 'Week 1'].set_index('SubjectID')
    week3_data = mutation_data[mutation_data['Week'] == 'Week 3'].set_index('SubjectID')

    # Find subjects with both timepoints
    paired_subjects = week1_data.index.intersection(week3_data.index)

    if len(paired_subjects) >= 5:  # Need at least 5 paired samples
        week1_freqs = week1_data.loc[paired_subjects, 'total_resistant_frequency'].values
        week3_freqs = week3_data.loc[paired_subjects, 'total_resistant_frequency'].values

        # Wilcoxon signed-rank test (paired)
        try:
            statistic, p_value = stats.wilcoxon(week1_freqs, week3_freqs)
        except ValueError:
            # If all differences are zero
            statistic, p_value = 0, 1.0

        longitudinal_results.append({
            'mutation_id': mutation,
            'species': mutation_data['species'].iloc[0],
            'gene': mutation_data['gene'].iloc[0],
            'n_paired': len(paired_subjects),
            'Week1_mean_freq': week1_freqs.mean(),
            'Week3_mean_freq': week3_freqs.mean(),
            'mean_change': week3_freqs.mean() - week1_freqs.mean(),
            'pct_change': ((week3_freqs.mean() - week1_freqs.mean()) / week1_freqs.mean() * 100) if week1_freqs.mean() > 0 else 0,
            'statistic': statistic,
            'p_value': p_value
        })

longitudinal_df = pd.DataFrame(longitudinal_results)

# FDR correction
longitudinal_df['fdr'] = multipletests(longitudinal_df['p_value'], method='fdr_bh')[1]

# Sort by absolute change
longitudinal_df = longitudinal_df.sort_values('mean_change', key=abs, ascending=False)

print(f"\nTop mutations with longitudinal changes:")
print(longitudinal_df[['mutation_id', 'Week1_mean_freq', 'Week3_mean_freq', 'mean_change', 'p_value', 'fdr']].head(10).to_string(index=False))

n_sig = (longitudinal_df['fdr'] < 0.05).sum()
print(f"\n  • {n_sig} mutations with significant longitudinal changes (FDR < 0.05)")

# Save longitudinal analysis
longitudinal_df.to_csv(
    os.path.join(RESULTS_DIR, "mutation_longitudinal_changes.tsv"),
    sep='\t',
    index=False
)
print(f"✓ Saved: {RESULTS_DIR}/mutation_longitudinal_changes.tsv")

# ============================================================================
# Step 6: Create visualizations
# ============================================================================
print("\n[6/6] Creating visualizations...")

# Figure 1: Heatmap of top mutations by location and body site
print("  - Creating mutation heatmap...")

# Get top 30 mutations
top_mutations = mutation_summary.head(30)['mutation_id'].tolist()

# Create matrix: mutations × samples
# For visualization, sample randomly to avoid too many columns
sample_subset = fq_data['sample_name'].unique()[:100]  # First 100 samples
heatmap_data = fq_data[
    (fq_data['mutation_id'].isin(top_mutations)) &
    (fq_data['sample_name'].isin(sample_subset))
]

# Pivot to matrix
matrix = heatmap_data.pivot_table(
    index='mutation_id',
    columns='sample_name',
    values='total_resistant_frequency',
    fill_value=0
)

# Sort by total frequency
matrix = matrix.loc[matrix.sum(axis=1).sort_values(ascending=False).index]

# Create figure
fig, ax = plt.subplots(figsize=(16, 10))

sns.heatmap(
    matrix,
    cmap='YlOrRd',
    vmin=0,
    vmax=1,
    cbar_kws={'label': 'Mutation Frequency'},
    xticklabels=False,  # Too many samples to label
    yticklabels=True,
    ax=ax
)

ax.set_xlabel('Samples', fontweight='bold', fontsize=12)
ax.set_ylabel('Mutation', fontweight='bold', fontsize=12)
ax.set_title('Top 30 FQ Resistance Mutations Across Samples', fontweight='bold', fontsize=14)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "mutation_heatmap_by_location.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/mutation_heatmap_by_location.pdf")

# Figure 2: Boxplots of mutation frequencies by body site (for significant mutations)
print("  - Creating mutation frequency boxplots...")

sig_bodysite = bodysite_df[bodysite_df['fdr'] < 0.05].head(6)  # Top 6 significant

if len(sig_bodysite) > 0:
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, (idx, row) in enumerate(sig_bodysite.iterrows()):
        if i >= 6:
            break

        mutation = row['mutation_id']
        mutation_data = fq_data[fq_data['mutation_id'] == mutation]

        # Create boxplot
        ax = axes[i]
        data_to_plot = [
            mutation_data[mutation_data['BodySite'] == 'Axilla']['total_resistant_frequency'],
            mutation_data[mutation_data['BodySite'] == 'Groin']['total_resistant_frequency'],
            mutation_data[mutation_data['BodySite'] == 'Stool']['total_resistant_frequency']
        ]

        bp = ax.boxplot(data_to_plot, labels=['Axilla', 'Groin', 'Stool'], patch_artist=True)

        # Color boxes
        colors = ['lightblue', 'lightgreen', 'lightsalmon']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)

        ax.set_ylabel('Mutation Frequency', fontsize=10)
        ax.set_title(f"{mutation.replace('_', ' ')}\n(FDR={row['fdr']:.2e})", fontsize=9)
        ax.grid(axis='y', alpha=0.3)

    # Hide unused subplots
    for i in range(len(sig_bodysite), 6):
        axes[i].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "mutation_frequency_boxplots.pdf"), dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {FIGURES_DIR}/mutation_frequency_boxplots.pdf")
else:
    print("  (No significant body site differences to plot)")

print("\n" + "=" * 80)
print("Phase 3 complete!")
print("=" * 80)

print("\nKey Findings:")
print(f"  • {len(common_mutations)} mutations tested (≥10 samples)")
n_sig_loc = (location_df['freq_fdr'] < 0.05).sum()
n_sig_site = (bodysite_df['fdr'] < 0.05).sum()
n_sig_long = (longitudinal_df['fdr'] < 0.05).sum()
print(f"  • {n_sig_loc} mutations differ by location (FDR < 0.05)")
print(f"  • {n_sig_site} mutations differ by body site (FDR < 0.05)")
print(f"  • {n_sig_long} mutations change longitudinally (FDR < 0.05)")

print("\nNext step: Run scripts/python/13_fq_gene_level_analysis.py")
