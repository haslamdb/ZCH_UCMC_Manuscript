#!/usr/bin/env python3
"""
BSI-Microbiome Similarity Analysis

Computes Bray-Curtis similarity between BSI (bloodstream infection) profiles
and microbiome centroids at different body sites/weeks for ZCH and UCMC.

Creates improved visualizations:
1. Combined heatmap with both institutions
2. Dumbbell plot for direct ZCH vs UCMC comparison
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.spatial.distance import braycurtis
from scipy.cluster.hierarchy import linkage, dendrogram
import os

# Set up matplotlib for publication-quality PDFs
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


def load_data():
    """Load microbiome and BSI data."""
    # Load microbiome species data
    microbiome_df = pd.read_csv(
        "/home/david/projects/ZCH_UCMC_Manuscript/data/NICUSpeciesReduced.csv",
        index_col=0
    )

    # Load metadata
    metadata_df = pd.read_csv(
        "/home/david/projects/ZCH_UCMC_Manuscript/metadata/AllNICUSampleKey20250206.csv",
        index_col=0
    )

    # Load BSI data
    bsi_df = pd.read_csv(
        "/home/david/projects/ZCH_UCMC_Manuscript/metadata/BSIData20250206.csv"
    )
    bsi_df = bsi_df.set_index('Species')

    return microbiome_df, metadata_df, bsi_df


def compute_centroids(microbiome_df, metadata_df):
    """Compute mean abundance centroids for each body site/week/location combination."""
    # Get common samples
    common_samples = list(set(microbiome_df.index) & set(metadata_df.index))
    microbiome_df = microbiome_df.loc[common_samples]
    metadata_df = metadata_df.loc[common_samples]

    # Normalize to relative abundance
    microbiome_rel = microbiome_df.div(microbiome_df.sum(axis=1), axis=0)

    # Add metadata columns
    microbiome_rel['Location'] = metadata_df.loc[microbiome_rel.index, 'Location']
    microbiome_rel['SampleType'] = metadata_df.loc[microbiome_rel.index, 'SampleType']
    microbiome_rel['SampleCollectionWeek'] = metadata_df.loc[microbiome_rel.index, 'SampleCollectionWeek']

    # Filter to body sites of interest
    microbiome_rel = microbiome_rel[
        microbiome_rel['SampleType'].isin(['Stool', 'Groin', 'Axilla'])
    ]

    # Create grouping variable
    microbiome_rel['Group'] = (
        microbiome_rel['SampleType'].astype(str) + '-' +
        microbiome_rel['SampleCollectionWeek'].astype(str) + '-' +
        microbiome_rel['Location'].astype(str)
    )

    # Compute centroids (mean across samples in each group)
    species_cols = [c for c in microbiome_rel.columns if c not in ['Location', 'SampleType', 'SampleCollectionWeek', 'Group']]
    centroids = microbiome_rel.groupby('Group')[species_cols].mean()

    return centroids


def compute_bray_curtis_similarity(centroids, bsi_df):
    """Compute Bray-Curtis similarity between BSI profiles and microbiome centroids."""
    # Normalize BSI to relative abundance
    bsi_rel = bsi_df.div(bsi_df.sum(axis=0), axis=1)

    # Get common species
    common_species = list(set(centroids.columns) & set(bsi_rel.index))

    # Align data
    centroids_aligned = centroids[common_species]
    bsi_aligned = bsi_rel.loc[common_species]

    # Compute Bray-Curtis similarity (1 - distance) for each combination
    similarity_data = []

    for group in centroids_aligned.index:
        centroid_vec = centroids_aligned.loc[group].values

        for bsi_source in bsi_aligned.columns:
            bsi_vec = bsi_aligned[bsi_source].values

            # Bray-Curtis distance, convert to similarity
            bc_dist = braycurtis(centroid_vec, bsi_vec)
            bc_sim = 1 - bc_dist

            # Parse group info
            parts = group.split('-')
            sample_type = parts[0]
            week = parts[1].replace('Week.', 'Week ')
            location = parts[2]

            # Map BSI source to location
            bsi_location = 'ZCH' if bsi_source == 'ZCH' else 'UCMC'

            similarity_data.append({
                'BodySite': sample_type,
                'Week': week,
                'MicrobiomeLocation': location,
                'BSI_Source': bsi_location,
                'Group': f"{sample_type} {week}",
                'Similarity': bc_sim,
                'Distance': bc_dist
            })

    return pd.DataFrame(similarity_data)


def create_combined_heatmap(similarity_df, output_file="bsi_similarity_heatmap.pdf"):
    """
    Create a single heatmap with BSI source (rows) × body site/week (columns).
    Shows only same-institution comparisons (ZCH BSI vs ZCH microbiome, UCMC BSI vs UCMC microbiome).
    Uses same color scale as full heatmap for consistency.
    """
    # Filter to same-location comparisons (ZCH BSI vs ZCH microbiome, UCMC BSI vs UCMC microbiome)
    same_location = similarity_df[
        ((similarity_df['BSI_Source'] == 'ZCH') & (similarity_df['MicrobiomeLocation'] == 'Hangzhou')) |
        ((similarity_df['BSI_Source'] == 'UCMC') & (similarity_df['MicrobiomeLocation'] == 'Cincinnati'))
    ].copy()

    # Create row labels matching the full heatmap style
    same_location['Comparison'] = same_location['BSI_Source'] + ' BSI vs ' + same_location['BSI_Source'] + ' Microbiome'

    # Pivot for heatmap
    pivot_data = same_location.pivot(index='Comparison', columns='Group', values='Similarity')

    # Reorder columns logically
    col_order = [
        'Axilla Week 1', 'Axilla Week 3',
        'Groin Week 1', 'Groin Week 3',
        'Stool Week 1', 'Stool Week 3'
    ]
    pivot_data = pivot_data[[c for c in col_order if c in pivot_data.columns]]

    # Reorder rows
    row_order = ['ZCH BSI vs ZCH Microbiome', 'UCMC BSI vs UCMC Microbiome']
    pivot_data = pivot_data.reindex([r for r in row_order if r in pivot_data.index])

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 4))

    # Create heatmap with same color scale as full heatmap
    sns.heatmap(
        pivot_data,
        annot=True,
        fmt='.3f',
        cmap='RdYlGn_r',  # Same as full heatmap: red = high, green = low
        vmin=0,
        vmax=0.8,  # Same scale as full heatmap
        linewidths=0.5,
        cbar_kws={'label': 'Bray-Curtis Similarity'},
        ax=ax
    )

    ax.set_xlabel('Microbiome Body Site / Week', fontsize=12)
    ax.set_ylabel('Comparison', fontsize=12)
    ax.set_title('BSI vs Microbiome Composition Similarity\n(Same-institution comparisons only)', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved heatmap to: {output_file}")
    return pivot_data


def create_full_heatmap(similarity_df, output_file="bsi_similarity_full_heatmap.pdf"):
    """
    Create a heatmap showing ALL combinations:
    - ZCH BSI vs ZCH microbiome
    - ZCH BSI vs UCMC microbiome
    - UCMC BSI vs ZCH microbiome
    - UCMC BSI vs UCMC microbiome
    """
    # Create combined row label
    similarity_df = similarity_df.copy()
    similarity_df['Comparison'] = similarity_df['BSI_Source'] + ' BSI vs ' + similarity_df['MicrobiomeLocation'].replace({
        'Hangzhou': 'ZCH', 'Cincinnati': 'UCMC'
    }) + ' Microbiome'

    # Pivot for heatmap
    pivot_data = similarity_df.pivot(index='Comparison', columns='Group', values='Similarity')

    # Reorder columns logically
    col_order = [
        'Axilla Week 1', 'Axilla Week 3',
        'Groin Week 1', 'Groin Week 3',
        'Stool Week 1', 'Stool Week 3'
    ]
    pivot_data = pivot_data[[c for c in col_order if c in pivot_data.columns]]

    # Reorder rows logically
    row_order = [
        'ZCH BSI vs ZCH Microbiome',
        'ZCH BSI vs UCMC Microbiome',
        'UCMC BSI vs ZCH Microbiome',
        'UCMC BSI vs UCMC Microbiome'
    ]
    pivot_data = pivot_data.reindex([r for r in row_order if r in pivot_data.index])

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 6))

    # Create heatmap (reversed colormap: red = high similarity, green = low)
    # vmax=0.8 so the 0.77 Groin Week 1 finding really stands out as dark red
    sns.heatmap(
        pivot_data,
        annot=True,
        fmt='.3f',
        cmap='RdYlGn_r',  # Reversed: red = high, green = low
        vmin=0,
        vmax=0.8,  # Cap at 0.8 so high values pop
        linewidths=0.5,
        cbar_kws={'label': 'Bray-Curtis Similarity'},
        ax=ax
    )

    ax.set_xlabel('Microbiome Body Site / Week', fontsize=12)
    ax.set_ylabel('Comparison', fontsize=12)
    ax.set_title('BSI vs Microbiome Composition Similarity\n(All cross-institutional combinations)', fontsize=14)

    # Add horizontal line to separate same-institution from cross-institution
    ax.axhline(y=2, color='black', linewidth=2)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved full heatmap to: {output_file}")
    return pivot_data


def create_faceted_heatmap(similarity_df, output_file="bsi_similarity_faceted_heatmap.pdf"):
    """
    Create a 2x2 faceted heatmap:
    - Columns: Microbiome source (ZCH, UCMC)
    - Rows: BSI source (ZCH, UCMC)
    Each cell shows body site/week similarity values.
    """
    # Map location names
    similarity_df = similarity_df.copy()
    similarity_df['MicrobiomeSite'] = similarity_df['MicrobiomeLocation'].replace({
        'Hangzhou': 'ZCH', 'Cincinnati': 'UCMC'
    })

    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=True, sharey=True)

    bsi_sources = ['ZCH', 'UCMC']
    microbiome_sources = ['ZCH', 'UCMC']

    col_order = [
        'Axilla Week 1', 'Axilla Week 3',
        'Groin Week 1', 'Groin Week 3',
        'Stool Week 1', 'Stool Week 3'
    ]

    for i, bsi in enumerate(bsi_sources):
        for j, micro in enumerate(microbiome_sources):
            ax = axes[i, j]

            # Filter data
            subset = similarity_df[
                (similarity_df['BSI_Source'] == bsi) &
                (similarity_df['MicrobiomeSite'] == micro)
            ]

            if len(subset) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                continue

            # Pivot
            pivot = subset.pivot_table(index='BSI_Source', columns='Group', values='Similarity', aggfunc='first')
            pivot = pivot[[c for c in col_order if c in pivot.columns]]

            # Plot
            sns.heatmap(
                pivot,
                annot=True,
                fmt='.3f',
                cmap='RdYlGn',
                center=0.5,
                vmin=0,
                vmax=1,
                linewidths=0.5,
                cbar=(j == 1),  # Only show colorbar on right column
                cbar_kws={'label': 'Similarity'} if j == 1 else {},
                ax=ax,
                yticklabels=False
            )

            # Labels
            if i == 0:
                ax.set_title(f'{micro} Microbiome', fontsize=14, fontweight='bold')
            if j == 0:
                ax.set_ylabel(f'{bsi} BSI', fontsize=14, fontweight='bold')
            else:
                ax.set_ylabel('')

            if i == 1:
                ax.set_xlabel('Body Site / Week', fontsize=11)
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            else:
                ax.set_xlabel('')

            # Highlight same-institution comparisons
            if bsi == micro:
                for spine in ax.spines.values():
                    spine.set_edgecolor('gold')
                    spine.set_linewidth(3)

    fig.suptitle('BSI vs Microbiome Similarity: All Combinations\n(Gold border = same institution)',
                 fontsize=16, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved faceted heatmap to: {output_file}")


def create_dumbbell_plot(similarity_df, output_file="bsi_similarity_dumbbell.pdf"):
    """
    Create a dumbbell plot showing ZCH and UCMC side-by-side for each body site/week.
    Makes institutional differences immediately apparent.
    """
    # Filter to same-location comparisons
    same_location = similarity_df[
        ((similarity_df['BSI_Source'] == 'ZCH') & (similarity_df['MicrobiomeLocation'] == 'Hangzhou')) |
        ((similarity_df['BSI_Source'] == 'UCMC') & (similarity_df['MicrobiomeLocation'] == 'Cincinnati'))
    ].copy()

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Get unique groups and sort
    groups = sorted(same_location['Group'].unique())
    y_positions = range(len(groups))

    # Colors
    zch_color = '#E74C3C'  # Red for ZCH
    ucmc_color = '#3498DB'  # Blue for UCMC

    for i, group in enumerate(groups):
        group_data = same_location[same_location['Group'] == group]

        zch_sim = group_data[group_data['BSI_Source'] == 'ZCH']['Similarity'].values
        ucmc_sim = group_data[group_data['BSI_Source'] == 'UCMC']['Similarity'].values

        if len(zch_sim) > 0 and len(ucmc_sim) > 0:
            zch_val = zch_sim[0]
            ucmc_val = ucmc_sim[0]

            # Draw connecting line
            ax.plot([zch_val, ucmc_val], [i, i], color='gray', linewidth=1, zorder=1)

            # Draw points
            ax.scatter(zch_val, i, color=zch_color, s=150, zorder=2, label='ZCH' if i == 0 else '')
            ax.scatter(ucmc_val, i, color=ucmc_color, s=150, zorder=2, label='UCMC' if i == 0 else '')

            # Add value labels
            ax.annotate(f'{zch_val:.3f}', (zch_val, i), textcoords="offset points",
                       xytext=(-5, 10), ha='center', fontsize=9, color=zch_color)
            ax.annotate(f'{ucmc_val:.3f}', (ucmc_val, i), textcoords="offset points",
                       xytext=(-5, -15), ha='center', fontsize=9, color=ucmc_color)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(groups, fontsize=11)
    ax.set_xlabel('Bray-Curtis Similarity to BSI Profile', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_title('BSI-Microbiome Similarity: ZCH vs UCMC Comparison', fontsize=14)
    ax.legend(loc='lower right', fontsize=11)
    ax.grid(axis='x', alpha=0.3)

    # Add vertical line at 0.5 for reference
    ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved dumbbell plot to: {output_file}")


def create_clustered_heatmap(similarity_df, output_file="bsi_similarity_clustered_heatmap.pdf"):
    """
    Create a heatmap with hierarchical clustering to group similar body site/week combinations.
    """
    # Filter to same-location comparisons
    same_location = similarity_df[
        ((similarity_df['BSI_Source'] == 'ZCH') & (similarity_df['MicrobiomeLocation'] == 'Hangzhou')) |
        ((similarity_df['BSI_Source'] == 'UCMC') & (similarity_df['MicrobiomeLocation'] == 'Cincinnati'))
    ].copy()

    # Pivot for heatmap - include location in group name
    same_location['FullGroup'] = same_location['Group'] + ' (' + same_location['BSI_Source'] + ')'
    pivot_data = same_location.pivot(index='FullGroup', columns='BSI_Source', values='Similarity')

    # Create wide format for clustering
    wide_data = same_location.pivot(index='Group', columns='BSI_Source', values='Similarity')

    # Create clustermap
    g = sns.clustermap(
        wide_data,
        annot=True,
        fmt='.3f',
        cmap='RdYlGn',
        center=0.5,
        vmin=0,
        vmax=1,
        linewidths=0.5,
        figsize=(8, 10),
        dendrogram_ratio=0.15,
        cbar_kws={'label': 'Bray-Curtis Similarity'}
    )

    g.ax_heatmap.set_xlabel('BSI Source', fontsize=12)
    g.ax_heatmap.set_ylabel('Body Site / Week', fontsize=12)
    g.fig.suptitle('BSI vs Microbiome Similarity with Clustering', fontsize=14, y=1.02)

    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved clustered heatmap to: {output_file}")


def create_difference_plot(similarity_df, output_file="bsi_similarity_difference.pdf"):
    """
    Create a bar plot showing the difference in similarity between ZCH and UCMC
    for each body site/week. Positive = ZCH more similar, Negative = UCMC more similar.
    """
    # Filter to same-location comparisons
    same_location = similarity_df[
        ((similarity_df['BSI_Source'] == 'ZCH') & (similarity_df['MicrobiomeLocation'] == 'Hangzhou')) |
        ((similarity_df['BSI_Source'] == 'UCMC') & (similarity_df['MicrobiomeLocation'] == 'Cincinnati'))
    ].copy()

    # Compute difference for each group
    diff_data = []
    for group in same_location['Group'].unique():
        group_data = same_location[same_location['Group'] == group]
        zch_sim = group_data[group_data['BSI_Source'] == 'ZCH']['Similarity'].values
        ucmc_sim = group_data[group_data['BSI_Source'] == 'UCMC']['Similarity'].values

        if len(zch_sim) > 0 and len(ucmc_sim) > 0:
            diff_data.append({
                'Group': group,
                'Difference': zch_sim[0] - ucmc_sim[0],
                'ZCH': zch_sim[0],
                'UCMC': ucmc_sim[0]
            })

    diff_df = pd.DataFrame(diff_data)
    diff_df = diff_df.sort_values('Difference')

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = ['#E74C3C' if x > 0 else '#3498DB' for x in diff_df['Difference']]

    bars = ax.barh(diff_df['Group'], diff_df['Difference'], color=colors, edgecolor='darkgray')

    ax.axvline(x=0, color='black', linewidth=1)
    ax.set_xlabel('Similarity Difference (ZCH - UCMC)', fontsize=12)
    ax.set_title('Which Institution\'s Microbiome is More Similar to BSI?\n(Red = ZCH, Blue = UCMC)', fontsize=14)
    ax.grid(axis='x', alpha=0.3)

    # Add value labels
    for bar, val in zip(bars, diff_df['Difference']):
        x_pos = val + 0.01 if val > 0 else val - 0.01
        ha = 'left' if val > 0 else 'right'
        ax.text(x_pos, bar.get_y() + bar.get_height()/2, f'{val:+.3f}',
               ha=ha, va='center', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved difference plot to: {output_file}")

    return diff_df


if __name__ == "__main__":
    # Create output directory
    output_dir = "results/bsi_similarity"
    os.makedirs(output_dir, exist_ok=True)

    print("Loading data...")
    microbiome_df, metadata_df, bsi_df = load_data()
    print(f"  Microbiome: {microbiome_df.shape}")
    print(f"  Metadata: {metadata_df.shape}")
    print(f"  BSI: {bsi_df.shape}")

    print("\nComputing centroids...")
    centroids = compute_centroids(microbiome_df, metadata_df)
    print(f"  Centroids: {centroids.shape}")
    print(f"  Groups: {list(centroids.index)}")

    print("\nComputing Bray-Curtis similarities...")
    similarity_df = compute_bray_curtis_similarity(centroids, bsi_df)
    similarity_df.to_csv(os.path.join(output_dir, "bsi_microbiome_similarity.csv"), index=False)
    print(f"  Saved similarity data to {output_dir}/bsi_microbiome_similarity.csv")

    print("\nCreating visualizations...")

    # 1. Combined heatmap (same-institution only)
    pivot_data = create_combined_heatmap(
        similarity_df,
        output_file=os.path.join(output_dir, "bsi_similarity_heatmap.pdf")
    )

    # 2. Full heatmap with ALL combinations (cross-institutional)
    full_pivot = create_full_heatmap(
        similarity_df,
        output_file=os.path.join(output_dir, "bsi_similarity_full_heatmap.pdf")
    )

    # 3. Faceted 2x2 heatmap
    create_faceted_heatmap(
        similarity_df,
        output_file=os.path.join(output_dir, "bsi_similarity_faceted_heatmap.pdf")
    )

    # 4. Dumbbell plot
    create_dumbbell_plot(
        similarity_df,
        output_file=os.path.join(output_dir, "bsi_similarity_dumbbell.pdf")
    )

    # 5. Clustered heatmap
    create_clustered_heatmap(
        similarity_df,
        output_file=os.path.join(output_dir, "bsi_similarity_clustered_heatmap.pdf")
    )

    # 6. Difference plot
    diff_df = create_difference_plot(
        similarity_df,
        output_file=os.path.join(output_dir, "bsi_similarity_difference.pdf")
    )

    print("\nKey findings:")
    print("\nMost similar microbiome site to BSI:")

    # ZCH
    zch_data = similarity_df[
        (similarity_df['BSI_Source'] == 'ZCH') &
        (similarity_df['MicrobiomeLocation'] == 'Hangzhou')
    ].sort_values('Similarity', ascending=False)
    print(f"\n  ZCH: {zch_data.iloc[0]['Group']} (similarity = {zch_data.iloc[0]['Similarity']:.3f})")

    # UCMC
    ucmc_data = similarity_df[
        (similarity_df['BSI_Source'] == 'UCMC') &
        (similarity_df['MicrobiomeLocation'] == 'Cincinnati')
    ].sort_values('Similarity', ascending=False)
    print(f"  UCMC: {ucmc_data.iloc[0]['Group']} (similarity = {ucmc_data.iloc[0]['Similarity']:.3f})")

    print(f"\nAnalysis complete! Results saved to '{output_dir}'")
