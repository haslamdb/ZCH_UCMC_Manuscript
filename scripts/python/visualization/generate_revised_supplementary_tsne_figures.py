#!/usr/bin/env python3
"""
Generate revised Supplementary Figures 2 & 3 (t-SNE plots) with corrections:
- Supp Fig 2: Different colors for study sites vs body sites (use shapes for sites)
- Supp Fig 3: Break down by two sites separately
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import seaborn as sns
import microbiome_transform as mt

# Set style
plt.style.use('seaborn-v0_8-whitegrid')

# Load microbiome data
print("Loading data...")
microbiome_df = pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)
metadata_df = pd.read_csv("../metadata/AllNICUSampleKey20250206.csv", index_col=0)

# Transform microbiome data
microbiome_clr = mt.clr_transform(microbiome_df)

# Merge with metadata
data = microbiome_clr.merge(metadata_df, left_index=True, right_index=True, how='inner')
print(f"Data shape: {data.shape}")

# ============================================================================
# Supplementary Figure 2: t-SNE colored by Body Site AND Location
# FIXED: Use shapes for location, colors for body site
# ============================================================================

print("\n=== Generating Supplementary Figure 2 ===")

# Perform t-SNE
tsne = TSNE(n_components=2, perplexity=30, max_iter=1000, random_state=42)
tsne_results = tsne.fit_transform(microbiome_clr)

tsne_df = pd.DataFrame(tsne_results, columns=['tSNE1', 'tSNE2'], index=microbiome_clr.index)
tsne_with_meta = tsne_df.merge(metadata_df, left_index=True, right_index=True, how='inner')

# Filter for relevant samples
tsne_plot_data = tsne_with_meta[tsne_with_meta['SampleType'].isin(['Axilla', 'Groin', 'Stool'])]

# Create Supp Fig 2 with improved color scheme
# Body sites = colors, Locations = shapes
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

# Panel A: Colored by body site, shaped by location
body_site_colors = {'Axilla': '#e41a1c', 'Groin': '#377eb8', 'Stool': '#4daf4a'}
location_markers = {'Cincinnati': 'o', 'Hangzhou': '^'}

for site in ['Axilla', 'Groin', 'Stool']:
    for loc in ['Cincinnati', 'Hangzhou']:
        mask = (tsne_plot_data['SampleType'] == site) & (tsne_plot_data['Location'] == loc)
        ax1.scatter(tsne_plot_data.loc[mask, 'tSNE1'],
                   tsne_plot_data.loc[mask, 'tSNE2'],
                   c=body_site_colors[site],
                   marker=location_markers[loc],
                   s=100,
                   alpha=0.7,
                   label=f'{site}-{loc}',
                   edgecolors='black',
                   linewidths=0.5)

ax1.set_xlabel('t-SNE Component 1', fontsize=14, fontweight='bold')
ax1.set_ylabel('t-SNE Component 2', fontsize=14, fontweight='bold')
ax1.set_title('a. Samples colored by Body Site (shape = Location)', fontsize=16, fontweight='bold')
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, ncol=2)
ax1.grid(alpha=0.3)

# Panel B: Colored by location, shaped by body site
location_colors = {'Cincinnati': '#E69F00', 'Hangzhou': '#56B4E9'}
body_site_markers = {'Axilla': 'o', 'Groin': 's', 'Stool': '^'}

for loc in ['Cincinnati', 'Hangzhou']:
    for site in ['Axilla', 'Groin', 'Stool']:
        mask = (tsne_plot_data['SampleType'] == site) & (tsne_plot_data['Location'] == loc)
        ax2.scatter(tsne_plot_data.loc[mask, 'tSNE1'],
                   tsne_plot_data.loc[mask, 'tSNE2'],
                   c=location_colors[loc],
                   marker=body_site_markers[site],
                   s=100,
                   alpha=0.7,
                   label=f'{loc}-{site}',
                   edgecolors='black',
                   linewidths=0.5)

ax2.set_xlabel('t-SNE Component 1', fontsize=14, fontweight='bold')
ax2.set_ylabel('t-SNE Component 2', fontsize=14, fontweight='bold')
ax2.set_title('b. Samples colored by Location (shape = Body Site)', fontsize=16, fontweight='bold')
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('revision_figures/SuppFig2_tSNE_revised_color_scheme.pdf', dpi=300, bbox_inches='tight')
plt.savefig('revision_figures/SuppFig2_tSNE_revised_color_scheme.png', dpi=300, bbox_inches='tight')
plt.close()

print("✓ Saved Supplementary Figure 2 with revised color scheme")

# ============================================================================
# Supplementary Figure 3: t-SNE colored by organism abundance
# FIXED: Separate panels for Cincinnati and Hangzhou
# ============================================================================

print("\n=== Generating Supplementary Figure 3 ===")

# Key BSI-associated organisms
key_organisms = ["Staphylococcus.aureus",
                 "Klebsiella.pneumoniae",
                 "Klebsiella.oxytoca",
                 "Enterococcus.faecium",
                 "Enterococcus.faecalis",
                 "Serratia.marcescens",
                 "Escherichia.coli",
                 "Streptococcus.pyogenes"]

# Create multi-page PDF with side-by-side comparisons for each organism
from matplotlib.backends.backend_pdf import PdfPages

with PdfPages('revision_figures/SuppFig3_tSNE_by_organism_by_site.pdf') as pdf:
    for organism in key_organisms:
        if organism not in data.columns:
            print(f"Warning: {organism} not in data")
            continue

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

        # Get abundance data
        abundance_data = tsne_with_meta.copy()
        abundance_data['abundance'] = data.loc[abundance_data.index, organism]

        # Cincinnati panel
        cinci_data = abundance_data[abundance_data['Location'] == 'Cincinnati']
        scatter1 = ax1.scatter(cinci_data['tSNE1'],
                              cinci_data['tSNE2'],
                              c=cinci_data['abundance'],
                              cmap='YlOrRd',
                              s=80,
                              alpha=0.8,
                              edgecolors='black',
                              linewidths=0.5)
        ax1.set_xlabel('t-SNE Component 1', fontsize=12, fontweight='bold')
        ax1.set_ylabel('t-SNE Component 2', fontsize=12, fontweight='bold')
        ax1.set_title(f'Cincinnati (n={len(cinci_data)})', fontsize=14, fontweight='bold')
        ax1.grid(alpha=0.3)
        cbar1 = plt.colorbar(scatter1, ax=ax1)
        cbar1.set_label('Abundance (CLR)', fontsize=11)

        # Hangzhou panel
        hangz_data = abundance_data[abundance_data['Location'] == 'Hangzhou']
        scatter2 = ax2.scatter(hangz_data['tSNE1'],
                              hangz_data['tSNE2'],
                              c=hangz_data['abundance'],
                              cmap='YlOrRd',
                              s=80,
                              alpha=0.8,
                              edgecolors='black',
                              linewidths=0.5)
        ax2.set_xlabel('t-SNE Component 1', fontsize=12, fontweight='bold')
        ax2.set_ylabel('t-SNE Component 2', fontsize=12, fontweight='bold')
        ax2.set_title(f'Hangzhou (n={len(hangz_data)})', fontsize=14, fontweight='bold')
        ax2.grid(alpha=0.3)
        cbar2 = plt.colorbar(scatter2, ax=ax2)
        cbar2.set_label('Abundance (CLR)', fontsize=11)

        # Overall title
        fig.suptitle(f'{organism.replace(".", " ")}',
                    fontsize=16, fontweight='bold', y=0.98)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

        print(f"  ✓ Generated {organism}")

print("\n✓ Saved Supplementary Figure 3 with site-separated panels")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "="*70)
print("SUPPLEMENTARY FIGURES GENERATION COMPLETE")
print("="*70)
print("\nGenerated files in revision_figures/:")
print("  - SuppFig2_tSNE_revised_color_scheme.pdf")
print("  - SuppFig2_tSNE_revised_color_scheme.png")
print("  - SuppFig3_tSNE_by_organism_by_site.pdf (multi-page, 8 organisms)")
print("\nKey improvements:")
print("  ✓ Supp Fig 2: Different colors for sites vs body sites")
print("  ✓ Supp Fig 3: Separate panels for Cincinnati and Hangzhou")
print("="*70)
