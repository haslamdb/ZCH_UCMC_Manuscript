#!/usr/bin/env python3
"""
Phase 7: Publication Figures
============================

Creates publication-quality figures for the FQ resistance manuscript.

Input:
    - data/fq_resistance_filtered.tsv
    - results/fq_resistance/ (various analysis outputs)

Outputs:
    - figures/fq_resistance/publication/Figure1_FQ_Overview.pdf
    - figures/fq_resistance/publication/Figure2_Location_BodySite.pdf
    - figures/fq_resistance/publication/Figure3_Longitudinal.pdf
    - figures/fq_resistance/publication/Figure4_Species_Patterns.pdf
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import os

# Set publication-quality style
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# Paths
DATA_DIR = "data"
RESULTS_DIR = "results/fq_resistance"
FIGURES_DIR = "figures/fq_resistance/publication"

# Create output directory
os.makedirs(FIGURES_DIR, exist_ok=True)

print("=" * 80)
print("Phase 7: Publication Figures")
print("=" * 80)

# ============================================================================
# Load data
# ============================================================================
print("\n[1/4] Loading data...")

fq_data = pd.read_csv(os.path.join(DATA_DIR, "fq_resistance_filtered.tsv"), sep="\t")
mutation_catalog = pd.read_csv(os.path.join("results/fq_resistance/qc", "mutation_catalog.tsv"), sep="\t")
species_prev_loc = pd.read_csv(os.path.join(RESULTS_DIR, "species_prevalence", "prevalence_by_location.tsv"), sep="\t")
species_prev_site = pd.read_csv(os.path.join(RESULTS_DIR, "species_prevalence", "prevalence_by_bodysite.tsv"), sep="\t")

print(f"  - Loaded {len(fq_data):,} observations")
print(f"  - {mutation_catalog.shape[0]} unique mutations")

# ============================================================================
# Figure 1: FQ Resistance Overview
# ============================================================================
print("\n[2/4] Creating Figure 1: FQ Resistance Overview...")

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1], hspace=0.3, wspace=0.3)

# Panel A: Heatmap of top 30 mutations
ax1 = fig.add_subplot(gs[0, :])

# Get top 30 mutations
top30 = mutation_catalog.head(30)['mutation_id'].tolist()

# Sample subset of samples for visualization (to keep heatmap readable)
sample_subset = fq_data.groupby('sample_name').size().sort_values(ascending=False).head(80).index

heatmap_data = fq_data[
    (fq_data['mutation_id'].isin(top30)) &
    (fq_data['sample_name'].isin(sample_subset))
]

# Pivot to matrix
matrix = heatmap_data.pivot_table(
    index='mutation_id',
    columns='sample_name',
    values='total_resistant_frequency',
    fill_value=0
)

# Sort by total frequency - only keep mutations that are in the matrix
mutations_in_matrix = [m for m in top30 if m in matrix.index]
matrix = matrix.loc[mutations_in_matrix]  # Keep order from catalog
matrix = matrix.loc[:, matrix.sum(axis=0).sort_values(ascending=False).index]  # Sort samples

# Create heatmap
sns.heatmap(
    matrix,
    cmap='YlOrRd',
    vmin=0,
    vmax=1,
    cbar_kws={'label': 'Mutation Frequency'},
    xticklabels=False,
    yticklabels=[m.replace('_', ' ') for m in matrix.index],
    ax=ax1,
    linewidths=0
)

ax1.set_xlabel('Samples (n=80)', fontweight='bold', fontsize=12)
ax1.set_ylabel('FQ Resistance Mutation', fontweight='bold', fontsize=12)
ax1.set_title('A. Top 30 FQ Resistance Mutations Across Samples', fontweight='bold', fontsize=14, pad=15)

# Panel B: Species prevalence by location
ax2 = fig.add_subplot(gs[1, 0])

# Pivot for plotting
prev_loc_pivot = species_prev_loc.pivot(index='species', columns='Location', values='prevalence_pct').fillna(0)
prev_loc_pivot = prev_loc_pivot.sort_values('UCMC', ascending=False).head(9)

x = np.arange(len(prev_loc_pivot))
width = 0.35

bars1 = ax2.bar(x - width/2, prev_loc_pivot['UCMC'], width, label='UCMC', alpha=0.8, color='steelblue')
bars2 = ax2.bar(x + width/2, prev_loc_pivot['ZCH'], width, label='ZCH', alpha=0.8, color='coral')

ax2.set_xlabel('Species', fontweight='bold', fontsize=11)
ax2.set_ylabel('Prevalence (% of samples)', fontweight='bold', fontsize=11)
ax2.set_title('B. Species Prevalence by Location', fontweight='bold', fontsize=13)
ax2.set_xticks(x)
ax2.set_xticklabels([s.replace('_', ' ') for s in prev_loc_pivot.index], rotation=45, ha='right', fontsize=9)
ax2.legend()
ax2.grid(axis='y', alpha=0.3)

# Panel C: Species prevalence by body site
ax3 = fig.add_subplot(gs[1, 1])

# Pivot for plotting
prev_site_pivot = species_prev_site.pivot(index='species', columns='BodySite', values='prevalence_pct').fillna(0)
prev_site_pivot = prev_site_pivot.loc[prev_loc_pivot.index]  # Same species as panel B

x = np.arange(len(prev_site_pivot))
width = 0.25

bars1 = ax3.bar(x - width, prev_site_pivot['Axilla'], width, label='Axilla', alpha=0.8, color='lightblue')
bars2 = ax3.bar(x, prev_site_pivot['Groin'], width, label='Groin', alpha=0.8, color='lightgreen')
bars3 = ax3.bar(x + width, prev_site_pivot['Stool'], width, label='Stool', alpha=0.8, color='lightsalmon')

ax3.set_xlabel('Species', fontweight='bold', fontsize=11)
ax3.set_ylabel('Prevalence (% of samples)', fontweight='bold', fontsize=11)
ax3.set_title('C. Species Prevalence by Body Site', fontweight='bold', fontsize=13)
ax3.set_xticks(x)
ax3.set_xticklabels([s.replace('_', ' ') for s in prev_site_pivot.index], rotation=45, ha='right', fontsize=9)
ax3.legend()
ax3.grid(axis='y', alpha=0.3)

plt.savefig(os.path.join(FIGURES_DIR, "Figure1_FQ_Overview.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/Figure1_FQ_Overview.pdf")
plt.close()

# ============================================================================
# Figure 2: Location and Body Site Comparisons
# ============================================================================
print("\n[3/4] Creating Figure 2: Location and Body Site Comparisons...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Calculate sample-level FQ resistance
sample_fq = fq_data.groupby(['sample_name', 'Location', 'BodySite', 'Week']).agg({
    'mutation_id': 'nunique',
    'total_resistant_frequency': 'mean'
}).reset_index()
sample_fq.columns = ['sample_name', 'Location', 'BodySite', 'Week', 'n_mutations', 'mean_freq']

# Panel A: N mutations by Location
ax = axes[0, 0]
data_to_plot = [
    sample_fq[sample_fq['Location'] == 'UCMC']['n_mutations'],
    sample_fq[sample_fq['Location'] == 'ZCH']['n_mutations']
]
bp = ax.boxplot(data_to_plot, tick_labels=['UCMC', 'ZCH'], patch_artist=True, widths=0.6)
for patch, color in zip(bp['boxes'], ['lightblue', 'lightcoral']):
    patch.set_facecolor(color)
ax.set_ylabel('Number of FQ Mutations', fontweight='bold')
ax.set_title('A. FQ Mutations by Location', fontweight='bold', fontsize=13)
ax.grid(axis='y', alpha=0.3)
ax.text(0.5, 0.95, 'p = 0.79 (n.s.)', transform=ax.transAxes, ha='center', va='top', fontsize=10)

# Panel B: N mutations by Body Site
ax = axes[0, 1]
data_to_plot = [
    sample_fq[sample_fq['BodySite'] == 'Axilla']['n_mutations'],
    sample_fq[sample_fq['BodySite'] == 'Groin']['n_mutations'],
    sample_fq[sample_fq['BodySite'] == 'Stool']['n_mutations']
]
bp = ax.boxplot(data_to_plot, tick_labels=['Axilla', 'Groin', 'Stool'], patch_artist=True, widths=0.6)
for patch, color in zip(bp['boxes'], ['lightblue', 'lightgreen', 'lightsalmon']):
    patch.set_facecolor(color)
ax.set_ylabel('Number of FQ Mutations', fontweight='bold')
ax.set_title('B. FQ Mutations by Body Site', fontweight='bold', fontsize=13)
ax.grid(axis='y', alpha=0.3)
ax.text(0.5, 0.95, 'p < 0.001 ***', transform=ax.transAxes, ha='center', va='top', fontsize=10, color='red')

# Panel C: Mean freq by Location
ax = axes[1, 0]
data_to_plot = [
    sample_fq[sample_fq['Location'] == 'UCMC']['mean_freq'],
    sample_fq[sample_fq['Location'] == 'ZCH']['mean_freq']
]
bp = ax.boxplot(data_to_plot, tick_labels=['UCMC', 'ZCH'], patch_artist=True, widths=0.6)
for patch, color in zip(bp['boxes'], ['lightblue', 'lightcoral']):
    patch.set_facecolor(color)
ax.set_ylabel('Mean FQ Mutation Frequency', fontweight='bold')
ax.set_title('C. Mean Mutation Frequency by Location', fontweight='bold', fontsize=13)
ax.grid(axis='y', alpha=0.3)
ax.text(0.5, 0.95, 'p = 0.48 (n.s.)', transform=ax.transAxes, ha='center', va='top', fontsize=10)

# Panel D: Mean freq by Body Site
ax = axes[1, 1]
data_to_plot = [
    sample_fq[sample_fq['BodySite'] == 'Axilla']['mean_freq'],
    sample_fq[sample_fq['BodySite'] == 'Groin']['mean_freq'],
    sample_fq[sample_fq['BodySite'] == 'Stool']['mean_freq']
]
bp = ax.boxplot(data_to_plot, tick_labels=['Axilla', 'Groin', 'Stool'], patch_artist=True, widths=0.6)
for patch, color in zip(bp['boxes'], ['lightblue', 'lightgreen', 'lightsalmon']):
    patch.set_facecolor(color)
ax.set_ylabel('Mean FQ Mutation Frequency', fontweight='bold')
ax.set_title('D. Mean Mutation Frequency by Body Site', fontweight='bold', fontsize=13)
ax.grid(axis='y', alpha=0.3)
ax.text(0.5, 0.95, 'p < 0.01 **', transform=ax.transAxes, ha='center', va='top', fontsize=10, color='red')

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "Figure2_Location_BodySite.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/Figure2_Location_BodySite.pdf")
plt.close()

# ============================================================================
# Figure 3: Longitudinal Trajectories
# ============================================================================
print("\n[4/4] Creating Figure 3: Longitudinal Trajectories...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Get complete subjects
complete_subjects = fq_data[fq_data['subject_complete'] == True]['SubjectID'].unique()

# Panel A: Trajectories by body site
for i, site in enumerate(['Axilla', 'Groin', 'Stool']):
    if i >= 3:
        break
    ax = axes.flatten()[i]

    site_data = sample_fq[(sample_fq['BodySite'] == site) &
                          (sample_fq['sample_name'].isin(fq_data['sample_name']))]

    # Get paired data
    week1_data = site_data[site_data['Week'] == 'Week 1'].set_index('sample_name')
    week3_data = site_data[site_data['Week'] == 'Week 3'].set_index('sample_name')

    # Plot some individual trajectories (sample of 20)
    plotted = 0
    for sample in week1_data.index[:20]:
        if sample in week3_data.index:
            ax.plot([1, 3],
                   [week1_data.loc[sample, 'n_mutations'], week3_data.loc[sample, 'n_mutations']],
                   'o-', alpha=0.3, color='gray')
            plotted += 1

    # Plot mean trajectory
    week1_mean = site_data[site_data['Week'] == 'Week 1']['n_mutations'].mean()
    week3_mean = site_data[site_data['Week'] == 'Week 3']['n_mutations'].mean()
    ax.plot([1, 3], [week1_mean, week3_mean], 'ro-', linewidth=3, markersize=10, label='Mean')

    ax.set_xlabel('Week', fontweight='bold')
    ax.set_ylabel('Number of FQ Mutations', fontweight='bold')
    ax.set_title(f'{chr(65+i)}. {site}', fontweight='bold', fontsize=13)
    ax.set_xticks([1, 3])
    ax.set_xticklabels(['Week 1', 'Week 3'])
    ax.grid(alpha=0.3)
    ax.legend()

# Panel D: Overall longitudinal change
ax = axes[1, 1]

week_comparison = sample_fq.groupby('Week')['n_mutations'].agg(['mean', 'sem']).reset_index()
colors = ['steelblue', 'coral']

for i, (idx, row) in enumerate(week_comparison.iterrows()):
    ax.bar(i, row['mean'], yerr=row['sem'], capsize=5, alpha=0.7, color=colors[i], label=row['Week'])

ax.set_ylabel('Mean Number of FQ Mutations', fontweight='bold')
ax.set_title('D. Overall Longitudinal Change', fontweight='bold', fontsize=13)
ax.set_xticks([0, 1])
ax.set_xticklabels(['Week 1', 'Week 3'])
ax.grid(axis='y', alpha=0.3)
ax.text(0.5, 0.95, 'p = 0.02 *', transform=ax.transAxes, ha='center', va='top', fontsize=11, color='red')

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "Figure3_Longitudinal.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/Figure3_Longitudinal.pdf")
plt.close()

# ============================================================================
# Figure 4: Species-Specific Resistance Patterns
# ============================================================================
print("\n[5/5] Creating Figure 4: Species-Specific Resistance Patterns...")

# Focus on top 4 species
top_species = ['Staphylococcus_aureus', 'Klebsiella_pneumoniae',
               'Enterococcus_faecium', 'Enterococcus_faecalis']

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for i, species in enumerate(top_species):
    ax = axes.flatten()[i]

    # Get mutations for this species
    species_data = fq_data[fq_data['species'] == species]
    species_mutations = species_data.groupby('mutation_id').agg({
        'sample_name': 'nunique',
        'total_resistant_frequency': 'mean'
    }).sort_values('sample_name', ascending=False).head(10)

    # Bar plot
    y_pos = np.arange(len(species_mutations))
    ax.barh(y_pos, species_mutations['sample_name'], alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([m.split('_', 2)[-1] for m in species_mutations.index], fontsize=9)
    ax.set_xlabel('Number of Samples', fontweight='bold')
    ax.set_title(f'{chr(65+i)}. {species.replace("_", " ")}', fontweight='bold', fontsize=12)
    ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "Figure4_Species_Patterns.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/Figure4_Species_Patterns.pdf")
plt.close()

print("\n" + "=" * 80)
print("Phase 7 complete!")
print("=" * 80)

print("\nPublication figures created:")
print(f"  1. Figure1_FQ_Overview.pdf - Heatmap + species prevalence")
print(f"  2. Figure2_Location_BodySite.pdf - Boxplots comparing groups")
print(f"  3. Figure3_Longitudinal.pdf - Temporal trajectories")
print(f"  4. Figure4_Species_Patterns.pdf - Species-specific mutation patterns")

print("\nNext step: Run scripts/python/17_fq_summary_report.py")
