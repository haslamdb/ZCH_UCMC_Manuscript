#!/usr/bin/env python3
"""
Phase 5: Antibiotic Exposure Correlations
=========================================

Tests whether non-FQ antibiotic exposure correlates with FQ resistance.

Input:
    - data/fq_resistance_filtered.tsv

Outputs:
    - results/fq_resistance/antibiotic_correlations/total_abx_correlation.tsv
    - results/fq_resistance/antibiotic_correlations/specific_abx_correlations.tsv
    - results/fq_resistance/antibiotic_correlations/abx_groups_comparison.tsv
    - figures/fq_resistance/antibiotic_fq_scatterplots.pdf
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Paths
DATA_DIR = "data"
RESULTS_DIR = "results/fq_resistance/antibiotic_correlations"
FIGURES_DIR = "figures/fq_resistance"

# Create output directories
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

print("=" * 80)
print("Phase 5: Antibiotic Exposure Correlations")
print("=" * 80)

# ============================================================================
# Step 1: Load data and create sample-level FQ resistance summary
# ============================================================================
print("\n[1/4] Loading data and calculating sample-level FQ resistance...")

fq_data = pd.read_csv(os.path.join(DATA_DIR, "fq_resistance_filtered.tsv"), sep="\t")
print(f"  - Loaded {len(fq_data):,} observations")

# Calculate per-sample FQ resistance burden
sample_fq = fq_data.groupby(['sample_name', 'Location', 'BodySite', 'Week', 'SubjectID',
                              'total_abx_days', 'any_abx']).agg({
    'mutation_id': 'nunique',  # Number of unique mutations
    'total_resistant_frequency': 'mean',  # Mean mutation frequency
    'is_high_frequency': 'sum'  # Number of high-frequency mutations
}).reset_index()

sample_fq.columns = [
    'sample_name', 'Location', 'BodySite', 'Week', 'SubjectID',
    'total_abx_days', 'any_abx',
    'n_fq_mutations', 'mean_fq_freq', 'n_high_freq_mutations'
]

print(f"  - {len(sample_fq)} samples with FQ resistance data")
print(f"\nAntibiotic exposure distribution:")
print(f"  Total abx days: {sample_fq['total_abx_days'].describe()}")
print(f"  Any antibiotic exposure: {sample_fq['any_abx'].sum()} / {len(sample_fq)} ({sample_fq['any_abx'].sum()/len(sample_fq)*100:.1f}%)")

# ============================================================================
# Step 2: Total antibiotic exposure vs FQ resistance
# ============================================================================
print("\n[2/4] Testing correlation: Total antibiotic days vs FQ resistance...")

# Overall correlation
overall_results = []

# Remove NaN values
valid_data = sample_fq[sample_fq['total_abx_days'].notna()]

# Correlation: total_abx_days vs n_fq_mutations
corr_mutations, p_mutations = stats.spearmanr(valid_data['total_abx_days'], valid_data['n_fq_mutations'])
overall_results.append({
    'comparison': 'Total abx days vs N mutations',
    'correlation': corr_mutations,
    'p_value': p_mutations,
    'n_samples': len(valid_data)
})

# Correlation: total_abx_days vs mean_fq_freq
corr_freq, p_freq = stats.spearmanr(valid_data['total_abx_days'], valid_data['mean_fq_freq'])
overall_results.append({
    'comparison': 'Total abx days vs Mean frequency',
    'correlation': corr_freq,
    'p_value': p_freq,
    'n_samples': len(valid_data)
})

print(f"\nOverall correlations:")
for result in overall_results:
    print(f"  {result['comparison']}: r={result['correlation']:.3f}, p={result['p_value']:.3e}")

# Stratified by body site
print(f"\nStratified by body site:")
bodysite_results = []

for site in ['Axilla', 'Groin', 'Stool']:
    site_data = valid_data[valid_data['BodySite'] == site]

    if len(site_data) >= 10:
        corr, p_val = stats.spearmanr(site_data['total_abx_days'], site_data['n_fq_mutations'])
        bodysite_results.append({
            'body_site': site,
            'comparison': 'Total abx days vs N mutations',
            'correlation': corr,
            'p_value': p_val,
            'n_samples': len(site_data)
        })
        print(f"  {site}: r={corr:.3f}, p={p_val:.3e} (n={len(site_data)})")

# Save total abx correlation results
total_abx_df = pd.DataFrame(overall_results + bodysite_results)
total_abx_df.to_csv(
    os.path.join(RESULTS_DIR, "total_abx_correlation.tsv"),
    sep='\t',
    index=False
)
print(f"\n✓ Saved: {RESULTS_DIR}/total_abx_correlation.tsv")

# ============================================================================
# Step 3: Specific antibiotic classes
# ============================================================================
print("\n[3/4] Testing specific antibiotic classes...")

# Load metadata with individual antibiotic columns
metadata = pd.read_csv("results/qc/master_metadata_with_qc.tsv", sep="\t")

# Merge with sample FQ data
sample_fq_full = sample_fq.merge(metadata, on='sample_name', how='left', suffixes=('', '_meta'))

# Define antibiotic groups
abx_groups = {
    'Beta-lactams': ['Ampicillin_w1', 'Penicillin_w1', 'Nafcillin_w1',
                     'Cefotaxime_w1', 'Ceftazidime_w1', 'Ceftriaxone_w1', 'Cefepime_w1',
                     'Piperacillin_Tazobactam_w1', 'Meropenem_w1'],
    'Aminoglycosides': ['Gentamicin_w1', 'Tobramycin_w1'],
    'Glycopeptides': ['Vancomycin_w1'],
    'Antifungals': ['Fluconazole_w1', 'Amphotericin_B_w1']
}

specific_results = []

for group_name, abx_list in abx_groups.items():
    # Check which antibiotics are in the data
    available_abx = [abx for abx in abx_list if abx in sample_fq_full.columns]

    if len(available_abx) > 0:
        # Create binary indicator: any antibiotic in this group
        sample_fq_full[f'{group_name}_any'] = (sample_fq_full[available_abx] > 0).any(axis=1)

        # Compare FQ resistance: exposed vs unexposed
        exposed = sample_fq_full[sample_fq_full[f'{group_name}_any'] == True]['n_fq_mutations']
        unexposed = sample_fq_full[sample_fq_full[f'{group_name}_any'] == False]['n_fq_mutations']

        if len(exposed) >= 5 and len(unexposed) >= 5:
            # Mann-Whitney U test
            statistic, p_value = stats.mannwhitneyu(exposed, unexposed, alternative='two-sided')

            specific_results.append({
                'antibiotic_group': group_name,
                'n_exposed': len(exposed),
                'mean_mutations_exposed': exposed.mean(),
                'n_unexposed': len(unexposed),
                'mean_mutations_unexposed': unexposed.mean(),
                'statistic': statistic,
                'p_value': p_value
            })

            print(f"\n{group_name}:")
            print(f"  Exposed (n={len(exposed)}): {exposed.mean():.2f} mutations")
            print(f"  Unexposed (n={len(unexposed)}): {unexposed.mean():.2f} mutations")
            print(f"  Mann-Whitney p={p_value:.3e}")

if len(specific_results) > 0:
    specific_df = pd.DataFrame(specific_results)
    from statsmodels.stats.multitest import multipletests
    specific_df['fdr'] = multipletests(specific_df['p_value'], method='fdr_bh')[1]

    specific_df.to_csv(
        os.path.join(RESULTS_DIR, "specific_abx_correlations.tsv"),
        sep='\t',
        index=False
    )
    print(f"\n✓ Saved: {RESULTS_DIR}/specific_abx_correlations.tsv")

# ============================================================================
# Step 4: High vs Low exposure groups
# ============================================================================
print("\n[4/4] Comparing high vs low antibiotic exposure...")

# Define groups
sample_fq['abx_group'] = 'No antibiotics'
sample_fq.loc[sample_fq['total_abx_days'] > 0, 'abx_group'] = 'Any antibiotics'

# Compare groups
no_abx = sample_fq[sample_fq['abx_group'] == 'No antibiotics']
any_abx_group = sample_fq[sample_fq['abx_group'] == 'Any antibiotics']

print(f"\nGroup sizes:")
print(f"  No antibiotics: n={len(no_abx)}")
print(f"  Any antibiotics: n={len(any_abx_group)}")

if len(no_abx) >= 5 and len(any_abx_group) >= 5:
    # Test for n_mutations
    stat_mut, p_mut = stats.mannwhitneyu(no_abx['n_fq_mutations'], any_abx_group['n_fq_mutations'])
    # Test for mean_freq
    stat_freq, p_freq = stats.mannwhitneyu(no_abx['mean_fq_freq'], any_abx_group['mean_fq_freq'])

    print(f"\nComparison (Mann-Whitney U):")
    print(f"  N mutations: No abx={no_abx['n_fq_mutations'].mean():.2f}, Any abx={any_abx_group['n_fq_mutations'].mean():.2f}, p={p_mut:.3e}")
    print(f"  Mean freq: No abx={no_abx['mean_fq_freq'].mean():.3f}, Any abx={any_abx_group['mean_fq_freq'].mean():.3f}, p={p_freq:.3e}")

    groups_df = pd.DataFrame([{
        'metric': 'N mutations',
        'no_abx_mean': no_abx['n_fq_mutations'].mean(),
        'any_abx_mean': any_abx_group['n_fq_mutations'].mean(),
        'statistic': stat_mut,
        'p_value': p_mut
    }, {
        'metric': 'Mean frequency',
        'no_abx_mean': no_abx['mean_fq_freq'].mean(),
        'any_abx_mean': any_abx_group['mean_fq_freq'].mean(),
        'statistic': stat_freq,
        'p_value': p_freq
    }])

    groups_df.to_csv(
        os.path.join(RESULTS_DIR, "abx_groups_comparison.tsv"),
        sep='\t',
        index=False
    )
    print(f"\n✓ Saved: {RESULTS_DIR}/abx_groups_comparison.tsv")

# ============================================================================
# Step 5: Create visualizations
# ============================================================================
print("\n[5/5] Creating visualizations...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: Total abx days vs N mutations (overall)
ax = axes[0, 0]
ax.scatter(valid_data['total_abx_days'], valid_data['n_fq_mutations'], alpha=0.5)
ax.set_xlabel('Total Antibiotic Days', fontsize=11, fontweight='bold')
ax.set_ylabel('Number of FQ Mutations', fontsize=11, fontweight='bold')
ax.set_title(f'Overall Correlation (r={corr_mutations:.3f}, p={p_mutations:.3e})', fontsize=12, fontweight='bold')
ax.grid(alpha=0.3)

# Panel B: Total abx days vs Mean frequency
ax = axes[0, 1]
ax.scatter(valid_data['total_abx_days'], valid_data['mean_fq_freq'], alpha=0.5, color='coral')
ax.set_xlabel('Total Antibiotic Days', fontsize=11, fontweight='bold')
ax.set_ylabel('Mean FQ Mutation Frequency', fontsize=11, fontweight='bold')
ax.set_title(f'Overall Correlation (r={corr_freq:.3f}, p={p_freq:.3e})', fontsize=12, fontweight='bold')
ax.grid(alpha=0.3)

# Panel C: Boxplot by body site
ax = axes[1, 0]
site_data_list = []
site_labels = []
for site in ['Axilla', 'Groin', 'Stool']:
    site_vals = valid_data[valid_data['BodySite'] == site]['n_fq_mutations']
    if len(site_vals) > 0:
        site_data_list.append(site_vals)
        site_labels.append(site)

if len(site_data_list) > 0:
    bp = ax.boxplot(site_data_list, tick_labels=site_labels, patch_artist=True)
    for patch, color in zip(bp['boxes'], ['lightblue', 'lightgreen', 'lightsalmon']):
        patch.set_facecolor(color)
    ax.set_ylabel('Number of FQ Mutations', fontsize=11, fontweight='bold')
    ax.set_title('FQ Resistance by Body Site', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

# Panel D: Abx exposure groups comparison
ax = axes[1, 1]
if len(no_abx) > 0 and len(any_abx_group) > 0:
    bp = ax.boxplot([no_abx['n_fq_mutations'], any_abx_group['n_fq_mutations']],
                     tick_labels=['No Antibiotics', 'Any Antibiotics'],
                     patch_artist=True)
    for patch, color in zip(bp['boxes'], ['lightgray', 'salmon']):
        patch.set_facecolor(color)
    ax.set_ylabel('Number of FQ Mutations', fontsize=11, fontweight='bold')
    ax.set_title(f'Antibiotic Exposure Groups (p={p_mut:.3e})', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "antibiotic_fq_scatterplots.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/antibiotic_fq_scatterplots.pdf")

print("\n" + "=" * 80)
print("Phase 5 complete!")
print("=" * 80)

print("\nKey Findings:")
print(f"  • Total abx days vs N mutations: r={corr_mutations:.3f}, p={p_mutations:.3e}")
print(f"  • Total abx days vs Mean freq: r={corr_freq:.3f}, p={p_freq:.3e}")
if len(specific_results) > 0:
    n_sig = (specific_df['fdr'] < 0.05).sum()
    print(f"  • {n_sig} antibiotic groups with significant associations (FDR < 0.05)")

print("\nNext step: Run scripts/python/15_fq_mixed_models.py")
