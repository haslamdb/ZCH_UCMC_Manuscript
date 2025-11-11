#!/usr/bin/env python3
"""
Supplementary Analysis: FQ Allele Frequency by Clinical Factors
================================================================

Creates detailed visualizations of FQ resistance allele frequencies (not just prevalence)
broken down by clinical factors (location, body site, week).

Focus: Mean allele frequencies for species, genes, and specific mutations.

Input:
    - data/fq_resistance_filtered.tsv

Outputs:
    - figures/fq_resistance/supplementary/S1_Species_AllelFreq_BodySite.pdf
    - figures/fq_resistance/supplementary/S2_Gene_AlleleFreq_BodySite.pdf
    - figures/fq_resistance/supplementary/S3_TopMutations_AlleleFreq_BodySite.pdf
    - figures/fq_resistance/supplementary/S4_AlleleFreq_by_Location_Week.pdf
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Set style
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

# Paths
DATA_DIR = "data"
FIGURES_DIR = "figures/fq_resistance/supplementary"

# Create output directory
os.makedirs(FIGURES_DIR, exist_ok=True)

print("=" * 80)
print("Supplementary Analysis: FQ Allele Frequency by Clinical Factors")
print("=" * 80)

# ============================================================================
# Load data
# ============================================================================
print("\n[1/5] Loading data...")

fq_data = pd.read_csv(os.path.join(DATA_DIR, "fq_resistance_filtered.tsv"), sep="\t")
print(f"  - Loaded {len(fq_data):,} observations")

# ============================================================================
# Figure S1: Species-specific allele frequencies by body site
# ============================================================================
print("\n[2/5] Creating Figure S1: Species Allele Frequencies by Body Site...")

# Get top 6 species by prevalence
top_species = fq_data.groupby('species')['sample_name'].nunique().sort_values(ascending=False).head(6).index

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

for i, species in enumerate(top_species):
    ax = axes[i]

    # Get data for this species
    species_data = fq_data[fq_data['species'] == species]

    # Calculate mean allele frequency per sample × body site
    sample_means = species_data.groupby(['sample_name', 'BodySite'])['total_resistant_frequency'].mean().reset_index()

    # Create boxplot
    data_to_plot = []
    labels = []
    for site in ['Axilla', 'Groin', 'Stool']:
        site_vals = sample_means[sample_means['BodySite'] == site]['total_resistant_frequency']
        if len(site_vals) > 0:
            data_to_plot.append(site_vals)
            labels.append(f"{site}\n(n={len(site_vals)})")

    if len(data_to_plot) > 0:
        bp = ax.boxplot(data_to_plot, tick_labels=labels, patch_artist=True, widths=0.6)

        # Color boxes
        colors = ['lightblue', 'lightgreen', 'lightsalmon']
        for patch, color in zip(bp['boxes'], colors[:len(data_to_plot)]):
            patch.set_facecolor(color)

        ax.set_ylabel('Mean Allele Frequency', fontweight='bold', fontsize=10)
        ax.set_title(species.replace('_', ' '), fontweight='bold', fontsize=11)
        ax.set_ylim(0, 1.05)
        ax.grid(axis='y', alpha=0.3)

        # Add statistical test
        if len(data_to_plot) == 3:
            try:
                stat, p = stats.kruskal(*data_to_plot)
                if p < 0.001:
                    ax.text(0.5, 0.95, f'p < 0.001 ***', transform=ax.transAxes,
                           ha='center', va='top', fontsize=9, color='red')
                elif p < 0.05:
                    ax.text(0.5, 0.95, f'p = {p:.3f} *', transform=ax.transAxes,
                           ha='center', va='top', fontsize=9, color='red')
            except:
                pass

plt.suptitle('Species-Specific FQ Resistance Allele Frequencies by Body Site',
             fontweight='bold', fontsize=14, y=0.995)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "S1_Species_AlleleFreq_BodySite.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/S1_Species_AlleleFreq_BodySite.pdf")
plt.close()

# ============================================================================
# Figure S2: Gene-level allele frequencies by body site
# ============================================================================
print("\n[3/5] Creating Figure S2: Gene Allele Frequencies by Body Site...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

genes = ['gyrA', 'gyrB', 'parC', 'parE']

for i, gene in enumerate(genes):
    ax = axes[i]

    # Get data for this gene
    gene_data = fq_data[fq_data['gene'] == gene]

    if len(gene_data) > 0:
        # Calculate mean allele frequency per sample × body site
        sample_means = gene_data.groupby(['sample_name', 'BodySite'])['total_resistant_frequency'].mean().reset_index()

        # Create boxplot
        data_to_plot = []
        labels = []
        for site in ['Axilla', 'Groin', 'Stool']:
            site_vals = sample_means[sample_means['BodySite'] == site]['total_resistant_frequency']
            if len(site_vals) > 0:
                data_to_plot.append(site_vals)
                labels.append(f"{site}\n(n={len(site_vals)})")

        if len(data_to_plot) > 0:
            bp = ax.boxplot(data_to_plot, tick_labels=labels, patch_artist=True, widths=0.6)

            # Color boxes
            colors = ['lightblue', 'lightgreen', 'lightsalmon']
            for patch, color in zip(bp['boxes'], colors[:len(data_to_plot)]):
                patch.set_facecolor(color)

            ax.set_ylabel('Mean Allele Frequency', fontweight='bold', fontsize=10)
            ax.set_title(f'{gene} mutations', fontweight='bold', fontsize=12)
            ax.set_ylim(0, 1.05)
            ax.grid(axis='y', alpha=0.3)

            # Add statistical test
            if len(data_to_plot) == 3:
                try:
                    stat, p = stats.kruskal(*data_to_plot)
                    if p < 0.001:
                        ax.text(0.5, 0.95, f'p < 0.001 ***', transform=ax.transAxes,
                               ha='center', va='top', fontsize=9, color='red')
                    elif p < 0.05:
                        ax.text(0.5, 0.95, f'p = {p:.3f} *', transform=ax.transAxes,
                               ha='center', va='top', fontsize=9, color='red')
                except:
                    pass

plt.suptitle('Gene-Level FQ Resistance Allele Frequencies by Body Site',
             fontweight='bold', fontsize=14, y=0.995)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "S2_Gene_AlleleFreq_BodySite.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/S2_Gene_AlleleFreq_BodySite.pdf")
plt.close()

# ============================================================================
# Figure S3: Top mutation-specific allele frequencies by body site
# ============================================================================
print("\n[4/5] Creating Figure S3: Top Mutations Allele Frequencies by Body Site...")

# Get top 12 most prevalent mutations
top_mutations = fq_data.groupby('mutation_id')['sample_name'].nunique().sort_values(ascending=False).head(12).index

fig, axes = plt.subplots(3, 4, figsize=(16, 12))
axes = axes.flatten()

for i, mutation in enumerate(top_mutations):
    ax = axes[i]

    # Get data for this mutation
    mut_data = fq_data[fq_data['mutation_id'] == mutation]

    # Create boxplot by body site
    data_to_plot = []
    labels = []
    for site in ['Axilla', 'Groin', 'Stool']:
        site_vals = mut_data[mut_data['BodySite'] == site]['total_resistant_frequency']
        if len(site_vals) > 0:
            data_to_plot.append(site_vals)
            labels.append(f"{site}\n(n={len(site_vals)})")

    if len(data_to_plot) > 0:
        bp = ax.boxplot(data_to_plot, tick_labels=labels, patch_artist=True, widths=0.6)

        # Color boxes
        colors = ['lightblue', 'lightgreen', 'lightsalmon']
        for patch, color in zip(bp['boxes'], colors[:len(data_to_plot)]):
            patch.set_facecolor(color)

        # Format mutation name for title (shorten if needed)
        mut_name = mutation.replace('_', ' ')
        if len(mut_name) > 40:
            parts = mut_name.split()
            mut_name = ' '.join(parts[:3]) + '\n' + ' '.join(parts[3:])

        ax.set_ylabel('Allele Frequency', fontweight='bold', fontsize=9)
        ax.set_title(mut_name, fontsize=9, fontweight='bold')
        ax.set_ylim(0, 1.05)
        ax.grid(axis='y', alpha=0.3)

        # Add p-value if significant
        if len(data_to_plot) == 3:
            try:
                stat, p = stats.kruskal(*data_to_plot)
                if p < 0.001:
                    ax.text(0.5, 0.95, 'p<0.001***', transform=ax.transAxes,
                           ha='center', va='top', fontsize=8, color='red')
                elif p < 0.05:
                    ax.text(0.5, 0.95, f'p={p:.3f}*', transform=ax.transAxes,
                           ha='center', va='top', fontsize=8, color='red')
            except:
                pass

plt.suptitle('Top 12 FQ Resistance Mutations: Allele Frequencies by Body Site',
             fontweight='bold', fontsize=14, y=0.995)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "S3_TopMutations_AlleleFreq_BodySite.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/S3_TopMutations_AlleleFreq_BodySite.pdf")
plt.close()

# ============================================================================
# Figure S4: Allele frequencies by Location and Week
# ============================================================================
print("\n[5/5] Creating Figure S4: Allele Frequencies by Location and Week...")

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Panel row 1: Top 3 species by Location
for i, species in enumerate(top_species[:3]):
    ax = axes[0, i]

    species_data = fq_data[fq_data['species'] == species]
    sample_means = species_data.groupby(['sample_name', 'Location'])['total_resistant_frequency'].mean().reset_index()

    data_to_plot = []
    labels = []
    for loc in ['UCMC', 'ZCH']:
        loc_vals = sample_means[sample_means['Location'] == loc]['total_resistant_frequency']
        if len(loc_vals) > 0:
            data_to_plot.append(loc_vals)
            labels.append(f"{loc}\n(n={len(loc_vals)})")

    if len(data_to_plot) > 0:
        bp = ax.boxplot(data_to_plot, tick_labels=labels, patch_artist=True, widths=0.6)
        for patch, color in zip(bp['boxes'], ['steelblue', 'coral']):
            patch.set_facecolor(color)

        ax.set_ylabel('Mean Allele Frequency', fontweight='bold', fontsize=10)
        ax.set_title(f'{species.replace("_", " ")}\nby Location', fontweight='bold', fontsize=11)
        ax.set_ylim(0, 1.05)
        ax.grid(axis='y', alpha=0.3)

        # Add p-value
        if len(data_to_plot) == 2:
            try:
                stat, p = stats.mannwhitneyu(*data_to_plot)
                if p < 0.05:
                    ax.text(0.5, 0.95, f'p={p:.3f}', transform=ax.transAxes,
                           ha='center', va='top', fontsize=9, color='red' if p < 0.05 else 'black')
            except:
                pass

# Panel row 2: Top 3 species by Week
for i, species in enumerate(top_species[:3]):
    ax = axes[1, i]

    species_data = fq_data[fq_data['species'] == species]
    sample_means = species_data.groupby(['sample_name', 'Week'])['total_resistant_frequency'].mean().reset_index()

    data_to_plot = []
    labels = []
    for week in ['Week 1', 'Week 3']:
        week_vals = sample_means[sample_means['Week'] == week]['total_resistant_frequency']
        if len(week_vals) > 0:
            data_to_plot.append(week_vals)
            labels.append(f"{week}\n(n={len(week_vals)})")

    if len(data_to_plot) > 0:
        bp = ax.boxplot(data_to_plot, tick_labels=labels, patch_artist=True, widths=0.6)
        for patch, color in zip(bp['boxes'], ['lightsteelblue', 'salmon']):
            patch.set_facecolor(color)

        ax.set_ylabel('Mean Allele Frequency', fontweight='bold', fontsize=10)
        ax.set_title(f'{species.replace("_", " ")}\nby Week', fontweight='bold', fontsize=11)
        ax.set_ylim(0, 1.05)
        ax.grid(axis='y', alpha=0.3)

        # Add p-value
        if len(data_to_plot) == 2:
            try:
                stat, p = stats.mannwhitneyu(*data_to_plot)
                if p < 0.05:
                    ax.text(0.5, 0.95, f'p={p:.3f}', transform=ax.transAxes,
                           ha='center', va='top', fontsize=9, color='red' if p < 0.05 else 'black')
            except:
                pass

plt.suptitle('Species-Specific Allele Frequencies by Location and Week',
             fontweight='bold', fontsize=14, y=0.995)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "S4_AlleleFreq_by_Location_Week.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/S4_AlleleFreq_by_Location_Week.pdf")
plt.close()

# ============================================================================
# Create summary table of mean allele frequencies
# ============================================================================
print("\n[6/6] Creating summary table of mean allele frequencies...")

# Species × Body Site
species_site_summary = []
for species in top_species:
    species_data = fq_data[fq_data['species'] == species]
    sample_means = species_data.groupby(['sample_name', 'BodySite'])['total_resistant_frequency'].mean().reset_index()

    for site in ['Axilla', 'Groin', 'Stool']:
        site_vals = sample_means[sample_means['BodySite'] == site]['total_resistant_frequency']
        if len(site_vals) > 0:
            species_site_summary.append({
                'Species': species,
                'BodySite': site,
                'N_Samples': len(site_vals),
                'Mean_AlleleFreq': site_vals.mean(),
                'Median_AlleleFreq': site_vals.median(),
                'SD_AlleleFreq': site_vals.std()
            })

species_site_df = pd.DataFrame(species_site_summary)
species_site_df.to_csv(os.path.join("results/fq_resistance/summary", "Species_AlleleFreq_by_BodySite.tsv"),
                       sep='\t', index=False)
print(f"  ✓ Saved: results/fq_resistance/summary/Species_AlleleFreq_by_BodySite.tsv")

print("\n" + "=" * 80)
print("Supplementary Analysis Complete!")
print("=" * 80)

print("\nGenerated supplementary figures:")
print("  S1: Species-specific allele frequencies by body site")
print("  S2: Gene-level allele frequencies by body site")
print("  S3: Top 12 mutations allele frequencies by body site")
print("  S4: Allele frequencies by location and week")

print("\nKey insights:")
print("  • These show ALLELE FREQUENCIES (0-100% resistant alleles)")
print("  • Not prevalence (presence/absence)")
print("  • Species-specific and mutation-specific patterns by clinical factors")
print("  • Statistical tests included on each panel")
