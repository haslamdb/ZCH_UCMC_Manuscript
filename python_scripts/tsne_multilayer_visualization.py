#!/usr/bin/env python3
"""
Create tSNE visualizations that combine multiple metadata layers:
- Body site (shape)
- Timepoint/week (size or border)
- Organism abundance (color intensity)

This allows us to see if high-abundance clusters are associated with
specific body sites and collection timepoints.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from sklearn.manifold import TSNE
import seaborn as sns
import microbiome_transform as mt
from matplotlib.backends.backend_pdf import PdfPages

# Set style
plt.style.use('seaborn-v0_8-whitegrid')

# Load data
print("Loading data...")
microbiome_df = pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)
metadata_df = pd.read_csv("../metadata/AllNICUSampleKey20250206.csv", index_col=0)

# Transform microbiome data
microbiome_clr = mt.clr_transform(microbiome_df)

# Merge with metadata
data = microbiome_clr.merge(metadata_df, left_index=True, right_index=True, how='inner')
print(f"Data shape: {data.shape}")

# Key organisms to visualize
key_organisms = [
    "Staphylococcus.aureus",
    "Klebsiella.pneumoniae",
    "Klebsiella.oxytoca",
    "Enterococcus.faecium",
    "Enterococcus.faecalis",
    "Serratia.marcescens",
    "Escherichia.coli",
    "Streptococcus.pyogenes",
    "Candida.parapsilosis"
]

# Perform t-SNE once
print("Computing t-SNE...")
tsne = TSNE(n_components=2, perplexity=30, max_iter=1000, random_state=42)
tsne_results = tsne.fit_transform(microbiome_clr)

tsne_df = pd.DataFrame(tsne_results, columns=['tSNE1', 'tSNE2'], index=microbiome_clr.index)
tsne_with_meta = tsne_df.merge(metadata_df, left_index=True, right_index=True, how='inner')

# Filter for body sites and weeks we're interested in
plot_data = tsne_with_meta[tsne_with_meta['SampleType'].isin(['Axilla', 'Groin', 'Stool'])].copy()

# Create binary week groups (early vs late)
# Handle Week.1 and Week.3 format
plot_data['WeekGroup'] = plot_data['SampleCollectionWeek'].apply(
    lambda x: 'Week 1' if str(x) == 'Week.1' else 'Week 3'
)

print(f"Filtered data shape: {plot_data.shape}")
print(f"Week groups: {plot_data['WeekGroup'].value_counts()}")

# ============================================================================
# OPTION 1: Shape (body site) + Color intensity (abundance) + Border (week)
# ============================================================================

def create_multilayer_plot_option1(plot_data, data, organism, location=None):
    """
    Create a single plot combining:
    - Shape: Body site (circle, square, triangle)
    - Color intensity: Organism abundance (white to red gradient)
    - Border thickness: Week (thin=early, thick=late)
    """
    fig, ax = plt.subplots(figsize=(12, 10))

    # Filter by location if specified
    if location:
        subset = plot_data[plot_data['Location'] == location].copy()
        title_suffix = f" ({location})"
    else:
        subset = plot_data.copy()
        title_suffix = ""

    # Get abundance data
    subset['abundance'] = data.loc[subset.index, organism]

    # Define markers for body sites
    body_site_markers = {
        'Stool': 'o',      # circle
        'Groin': 's',      # square
        'Axilla': '^'      # triangle
    }

    # Define border widths for weeks
    week_linewidths = {
        'Week 1': 0.5,   # thin border
        'Week 3': 2.0     # thick border
    }

    # Normalize abundance for color mapping (0-1 range)
    vmin, vmax = subset['abundance'].min(), subset['abundance'].max()

    # Create colormap (white to red)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('white_red', ['#ffffff', '#d73027'])

    # Plot each combination of body site and week
    for body_site in ['Stool', 'Groin', 'Axilla']:
        for week_group in ['Week 1', 'Week 3']:
            mask = (subset['SampleType'] == body_site) & (subset['WeekGroup'] == week_group)
            subset_pts = subset[mask]

            if len(subset_pts) > 0:
                scatter = ax.scatter(
                    subset_pts['tSNE1'],
                    subset_pts['tSNE2'],
                    c=subset_pts['abundance'],
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    marker=body_site_markers[body_site],
                    s=120,
                    alpha=0.8,
                    edgecolors='black',
                    linewidths=week_linewidths[week_group],
                    label=f'{body_site}-{week_group}'
                )

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(f'{organism.replace(".", " ")} abundance (CLR)', fontsize=12, fontweight='bold')

    # Create custom legend
    # Body site markers
    body_site_legend = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=10, label='Stool', markeredgecolor='black', markeredgewidth=1),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='gray',
               markersize=10, label='Groin', markeredgecolor='black', markeredgewidth=1),
        Line2D([0], [0], marker='^', color='w', markerfacecolor='gray',
               markersize=10, label='Axilla', markeredgecolor='black', markeredgewidth=1)
    ]

    # Week border legend
    week_legend = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=10, label='Week 1', markeredgecolor='black', markeredgewidth=0.5),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markersize=10, label='Week 3', markeredgecolor='black', markeredgewidth=2.0)
    ]

    # Add two legends
    legend1 = ax.legend(handles=body_site_legend, title='Body Site (shape)',
                       loc='upper left', fontsize=10, title_fontsize=11, framealpha=0.9)
    ax.add_artist(legend1)

    legend2 = ax.legend(handles=week_legend, title='Timepoint (border)',
                       loc='lower left', fontsize=10, title_fontsize=11, framealpha=0.9)

    ax.set_xlabel('t-SNE Component 1', fontsize=14, fontweight='bold')
    ax.set_ylabel('t-SNE Component 2', fontsize=14, fontweight='bold')
    ax.set_title(f'{organism.replace(".", " ")}{title_suffix}\n' +
                f'Shape=Body Site, Color=Abundance, Border=Week',
                fontsize=15, fontweight='bold')
    ax.grid(alpha=0.3)

    plt.tight_layout()
    return fig

# ============================================================================
# OPTION 2: Faceted by body site, color by abundance, marker by week
# ============================================================================

def create_multilayer_plot_option2(plot_data, data, organism, location=None):
    """
    Create faceted plot:
    - Panels: Body site (3 columns)
    - Color intensity: Organism abundance
    - Marker type: Week (circle vs star)
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Filter by location if specified
    if location:
        subset = plot_data[plot_data['Location'] == location].copy()
        title_suffix = f" - {location}"
    else:
        subset = plot_data.copy()
        title_suffix = ""

    # Get abundance data
    subset['abundance'] = data.loc[subset.index, organism]

    # Define markers for weeks
    week_markers = {
        'Week 1': 'o',   # circle
        'Week 3': '*'     # star
    }

    # Colormap
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('white_red', ['#ffffff', '#d73027'])

    # Get global min/max for consistent coloring
    vmin, vmax = subset['abundance'].min(), subset['abundance'].max()

    # Plot each body site in its own panel
    body_sites = ['Stool', 'Groin', 'Axilla']
    for idx, body_site in enumerate(body_sites):
        ax = axes[idx]
        site_data = subset[subset['SampleType'] == body_site]

        # Plot each week group
        for week_group in ['Week 1', 'Week 3']:
            week_data = site_data[site_data['WeekGroup'] == week_group]

            if len(week_data) > 0:
                scatter = ax.scatter(
                    week_data['tSNE1'],
                    week_data['tSNE2'],
                    c=week_data['abundance'],
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    marker=week_markers[week_group],
                    s=150 if week_group == 'Week 3' else 80,
                    alpha=0.8,
                    edgecolors='black',
                    linewidths=0.5,
                    label=week_group
                )

        ax.set_xlabel('t-SNE 1', fontsize=12, fontweight='bold')
        if idx == 0:
            ax.set_ylabel('t-SNE 2', fontsize=12, fontweight='bold')
        ax.set_title(f'{body_site}\n(n={len(site_data)})', fontsize=13, fontweight='bold')
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9, loc='best')

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, pad=0.02, aspect=30)
    cbar.set_label(f'{organism.replace(".", " ")} (CLR)', fontsize=11, fontweight='bold')

    fig.suptitle(f'{organism.replace(".", " ")}{title_suffix}',
                fontsize=16, fontweight='bold', y=1.02)

    plt.tight_layout()
    return fig

# ============================================================================
# OPTION 3: Faceted by location (columns) and body site (rows)
# ============================================================================

def create_multilayer_plot_option3(plot_data, data, organism):
    """
    Create faceted plot:
    - Columns: Location (Cincinnati, Hangzhou)
    - Rows: Body site (Axilla, Groin, Stool)
    - Color: Organism-specific color (distinct color per organism)
    - Marker type: Week (circle vs star)
    """
    fig, axes = plt.subplots(3, 2, figsize=(12, 16))

    # Get abundance data
    subset = plot_data.copy()
    subset['abundance'] = data.loc[subset.index, organism]

    # Define organism colors (using a colorblind-friendly palette)
    organism_colors = {
        "Staphylococcus.aureus": '#e41a1c',      # red
        "Klebsiella.pneumoniae": '#377eb8',      # blue
        "Klebsiella.oxytoca": '#4daf4a',         # green
        "Enterococcus.faecium": '#6a0dad',       # dark purple
        "Enterococcus.faecalis": '#ff7f00',      # orange
        "Serratia.marcescens": '#11a579',        # teal
        "Escherichia.coli": '#a65628',           # brown
        "Streptococcus.pyogenes": '#f781bf',     # pink
        "Candida.parapsilosis": '#00ced1'        # dark cyan
    }

    # Get the color for this organism
    organism_color = organism_colors.get(organism, '#808080')  # default to gray if not found

    # Create colormap from white to the organism's color
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('white_to_org', ['#ffffff', organism_color])

    # Get global min/max for consistent coloring
    vmin, vmax = subset['abundance'].min(), subset['abundance'].max()

    # Use power normalization to emphasize high-abundance samples
    # gamma < 1 makes colors darker for high values (emphasizes top end)
    # gamma = 0.6 provides moderate emphasis on high-abundance samples
    # This makes lower ~70% of samples lighter and top ~30% progressively darker
    from matplotlib.colors import PowerNorm
    norm = PowerNorm(gamma=0.6, vmin=vmin, vmax=vmax)

    # Define markers for weeks
    week_markers = {
        'Week 1': 'o',   # circle
        'Week 3': '*'     # star
    }

    # Plot each combination
    body_sites = ['Axilla', 'Groin', 'Stool']
    locations = ['Cincinnati', 'Hangzhou']

    for row_idx, body_site in enumerate(body_sites):
        for col_idx, location in enumerate(locations):
            ax = axes[row_idx, col_idx]
            site_loc_data = subset[(subset['SampleType'] == body_site) &
                                   (subset['Location'] == location)]

            # Plot each week group
            for week_group in ['Week 1', 'Week 3']:
                week_data = site_loc_data[site_loc_data['WeekGroup'] == week_group]

                if len(week_data) > 0:
                    scatter = ax.scatter(
                        week_data['tSNE1'],
                        week_data['tSNE2'],
                        c=week_data['abundance'],
                        cmap=cmap,
                        norm=norm,
                        marker=week_markers[week_group],
                        s=150 if week_group == 'Week 3' else 80,
                        alpha=0.8,
                        edgecolors='none',
                        label=week_group
                    )

            # Set labels
            if row_idx == 2:  # Bottom row
                ax.set_xlabel('t-SNE 1', fontsize=11, fontweight='bold')
            if col_idx == 0:  # Left column
                ax.set_ylabel('t-SNE 2', fontsize=11, fontweight='bold')

            # Title showing location on top row, body site on left
            if row_idx == 0:
                ax.set_title(f'{location}\n{body_site} (n={len(site_loc_data)})',
                           fontsize=12, fontweight='bold')
            else:
                ax.set_title(f'{body_site} (n={len(site_loc_data)})',
                           fontsize=12, fontweight='bold')

            ax.grid(alpha=0.3)

            # Only add legend to top-right panel
            if row_idx == 0 and col_idx == 1:
                ax.legend(fontsize=9, loc='upper left', title='Week', framealpha=0.9)

    # Add colorbar with more padding to avoid overlap
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Adjust subplots to make room for colorbar on the right
    plt.tight_layout(rect=[0, 0, 0.92, 0.98])

    # Add colorbar to the right with more padding
    cbar = fig.colorbar(sm, ax=axes, pad=0.03, aspect=40, fraction=0.046)
    cbar.set_label(f'{organism.replace(".", " ")} (CLR)', fontsize=12, fontweight='bold')

    fig.suptitle(f'{organism.replace(".", " ")}',
                fontsize=16, fontweight='bold', y=0.995)

    return fig

# ============================================================================
# OPTION 4: Multi-organism colored view with location × body site layout
# ============================================================================

def create_multilayer_plot_option4(plot_data, data, key_organisms):
    """
    Create faceted plot showing all organisms:
    - Columns: Location (Cincinnati, Hangzhou)
    - Rows: Body site (Axilla, Groin, Stool)
    - Color: Different organism (distinct colors for each organism)
    - Marker type: Week (circle vs star)
    """
    fig, axes = plt.subplots(3, 2, figsize=(14, 18))

    # Define organism colors (using a colorblind-friendly palette)
    organism_colors = {
        "Staphylococcus.aureus": '#e41a1c',      # red
        "Klebsiella.pneumoniae": '#377eb8',      # blue
        "Klebsiella.oxytoca": '#4daf4a',         # green
        "Enterococcus.faecium": '#6a0dad',       # dark purple
        "Enterococcus.faecalis": '#ff7f00',      # orange
        "Serratia.marcescens": '#11a579',        # teal
        "Escherichia.coli": '#a65628',           # brown
        "Streptococcus.pyogenes": '#f781bf',     # pink
        "Candida.parapsilosis": '#00ced1'        # dark cyan
    }

    # Define markers for weeks
    week_markers = {
        'Week 1': 'o',   # circle
        'Week 3': '*'     # star
    }

    # Plot each combination
    body_sites = ['Axilla', 'Groin', 'Stool']
    locations = ['Cincinnati', 'Hangzhou']

    for row_idx, body_site in enumerate(body_sites):
        for col_idx, location in enumerate(locations):
            ax = axes[row_idx, col_idx]
            subset = plot_data.copy()
            site_loc_data = subset[(subset['SampleType'] == body_site) &
                                   (subset['Location'] == location)]

            # For each organism, find samples where it's abundant
            # We'll color by the dominant organism at each sample
            site_loc_data['dominant_organism'] = 'Other'
            site_loc_data['max_abundance'] = -np.inf

            # Find which organism is most abundant at each sample
            for organism in key_organisms:
                if organism in data.columns:
                    abundances = data.loc[site_loc_data.index, organism]
                    mask = abundances > site_loc_data['max_abundance']
                    site_loc_data.loc[mask, 'dominant_organism'] = organism
                    site_loc_data.loc[mask, 'max_abundance'] = abundances[mask]

            # Plot each organism separately
            for organism in key_organisms:
                if organism not in data.columns:
                    continue

                org_data = site_loc_data[site_loc_data['dominant_organism'] == organism]

                # Plot each week group
                for week_group in ['Week 1', 'Week 3']:
                    week_data = org_data[org_data['WeekGroup'] == week_group]

                    if len(week_data) > 0:
                        ax.scatter(
                            week_data['tSNE1'],
                            week_data['tSNE2'],
                            c=organism_colors[organism],
                            marker=week_markers[week_group],
                            s=150 if week_group == 'Week 3' else 80,
                            alpha=0.7,
                            edgecolors='black',
                            linewidths=0.5,
                            label=f'{organism.replace(".", " ")} - {week_group}' if row_idx == 0 and col_idx == 0 else ''
                        )

            # Set labels
            if row_idx == 2:  # Bottom row
                ax.set_xlabel('t-SNE 1', fontsize=11, fontweight='bold')
            if col_idx == 0:  # Left column
                ax.set_ylabel('t-SNE 2', fontsize=11, fontweight='bold')

            # Title showing location on top row, body site on left
            if row_idx == 0:
                ax.set_title(f'{location}\n{body_site} (n={len(site_loc_data)})',
                           fontsize=12, fontweight='bold')
            else:
                ax.set_title(f'{body_site} (n={len(site_loc_data)})',
                           fontsize=12, fontweight='bold')

            ax.grid(alpha=0.3)

    # Create custom legend showing organisms and weeks
    from matplotlib.lines import Line2D

    # Organism legend elements
    org_legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=organism_colors[org],
               markersize=10, label=org.replace(".", " "),
               markeredgecolor='black', markeredgewidth=0.5)
        for org in key_organisms if org in organism_colors
    ]

    # Week legend elements
    week_legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor='gray', markersize=8,
               label='Week 1', markeredgecolor='black', markeredgewidth=0.5),
        Line2D([0], [0], marker='*', color='w',
               markerfacecolor='gray', markersize=12,
               label='Week 3', markeredgecolor='black', markeredgewidth=0.5)
    ]

    # Add legends
    legend1 = fig.legend(handles=org_legend_elements,
                        title='Dominant Organism',
                        loc='center left',
                        bbox_to_anchor=(0.92, 0.65),
                        fontsize=9,
                        title_fontsize=10,
                        framealpha=0.9)

    legend2 = fig.legend(handles=week_legend_elements,
                        title='Timepoint',
                        loc='center left',
                        bbox_to_anchor=(0.92, 0.3),
                        fontsize=9,
                        title_fontsize=10,
                        framealpha=0.9)

    fig.suptitle('Key Organisms Distribution Across Locations and Body Sites',
                fontsize=16, fontweight='bold', y=0.995)

    plt.tight_layout(rect=[0, 0, 0.90, 0.98])
    return fig

# ============================================================================
# Generate plots for all organisms
# ============================================================================

print("\n=== Generating multi-layer visualizations ===\n")

# Create output directory
import os
os.makedirs('revision_figures', exist_ok=True)

# Generate Option 1 plots (combined plot with shapes and borders)
print("Option 1: Shape + Color + Border...")
with PdfPages('revision_figures/tSNE_multilayer_option1_combined.pdf') as pdf:
    for organism in key_organisms:
        if organism in data.columns:
            print(f"  Processing {organism}...")

            # Combined (both locations)
            fig = create_multilayer_plot_option1(plot_data, data, organism)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

# Generate Option 1 plots separated by location
print("\nOption 1: Separated by location...")
with PdfPages('revision_figures/tSNE_multilayer_option1_by_location.pdf') as pdf:
    for organism in key_organisms:
        if organism in data.columns:
            print(f"  Processing {organism}...")

            # Cincinnati
            fig = create_multilayer_plot_option1(plot_data, data, organism, location='Cincinnati')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

            # Hangzhou
            fig = create_multilayer_plot_option1(plot_data, data, organism, location='Hangzhou')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

# Generate Option 2 plots (faceted by body site)
print("\nOption 2: Faceted by body site...")
with PdfPages('revision_figures/tSNE_multilayer_option2_faceted.pdf') as pdf:
    for organism in key_organisms:
        if organism in data.columns:
            print(f"  Processing {organism}...")

            # Combined (both locations)
            fig = create_multilayer_plot_option2(plot_data, data, organism)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

# Generate Option 2 plots separated by location
print("\nOption 2: Faceted by body site, separated by location...")
with PdfPages('revision_figures/tSNE_multilayer_option2_by_location.pdf') as pdf:
    for organism in key_organisms:
        if organism in data.columns:
            print(f"  Processing {organism}...")

            # Cincinnati
            fig = create_multilayer_plot_option2(plot_data, data, organism, location='Cincinnati')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

            # Hangzhou
            fig = create_multilayer_plot_option2(plot_data, data, organism, location='Hangzhou')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

# Generate Option 3 plots (location columns x body site rows)
print("\nOption 3: Location (columns) x Body site (rows)...")
with PdfPages('revision_figures/tSNE_multilayer_option3_location_by_bodysite.pdf') as pdf:
    for organism in key_organisms:
        if organism in data.columns:
            print(f"  Processing {organism}...")
            fig = create_multilayer_plot_option3(plot_data, data, organism)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

# Generate Option 4 plot (all organisms with different colors)
print("\nOption 4: Multi-organism colored view...")
fig = create_multilayer_plot_option4(plot_data, data, key_organisms)
plt.savefig('revision_figures/tSNE_multilayer_option4_all_organisms_colored.pdf', bbox_inches='tight')
plt.close(fig)
print("  ✓ Generated multi-organism plot")

print("\n" + "="*70)
print("MULTI-LAYER VISUALIZATION COMPLETE")
print("="*70)
print("\nGenerated files in revision_figures/:")
print("  - tSNE_multilayer_option1_combined.pdf (shape=site, border=week)")
print("  - tSNE_multilayer_option1_by_location.pdf (separated by location)")
print("  - tSNE_multilayer_option2_faceted.pdf (panels=site, marker=week)")
print("  - tSNE_multilayer_option2_by_location.pdf (faceted + separated)")
print("  - tSNE_multilayer_option3_location_by_bodysite.pdf (cols=location, rows=bodysite)")
print("  - tSNE_multilayer_option4_all_organisms_colored.pdf (NEW: multi-organism colored)")
print("\nVisualization features:")
print("  Option 1: Single plot with shape (body site) + border width (week)")
print("  Option 2: Faceted panels by body site with markers for week")
print("  Option 3: Faceted 3x2 grid (rows=body site, cols=location), color=abundance")
print("  Option 4: Faceted 3x2 grid (rows=body site, cols=location), color=organism")
print("="*70)
