#!/usr/bin/env python3
"""
Generate Publication-Ready Figures

Creates high-quality figures for manuscript:
1. Main Figure 1: Study design + PCA
2. Main Figure 2: Body site AMR burden (boxplots + heatmap)
3. Main Figure 3: Longitudinal trajectories (Week 1→3)
4. Main Figure 4: Gene-level changes (volcano plots)
5. Supplementary: Diversity metrics, clustering, model diagnostics
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import wilcoxon
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Define paths
PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
RESULTS_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results"
QC_DIR = RESULTS_DIR / "qc"
EXPLORATORY_DIR = RESULTS_DIR / "exploratory"
LONG_DIR = RESULTS_DIR / "longitudinal"
PUB_FIGURES_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "figures" / "publication"

# Create directory
PUB_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Set publication-quality plotting style
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts
plt.rcParams['ps.fonttype'] = 42

# Color schemes
LOCATION_COLORS = {'UCMC': '#2E86AB', 'ZCH': '#A23B72'}
SITE_COLORS = {'Axilla': '#F18F01', 'Groin': '#C73E1D', 'Stool': '#6A994E'}
WEEK_COLORS = {'Week.1': '#06A77D', 'Week.3': '#D62828'}

def load_data():
    """Load all data needed for figures"""
    print("Loading data...")

    # Metadata
    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")
    metadata = metadata[
        (metadata['has_amr_data']) &
        (metadata['SampleType'].isin(['Axilla', 'Groin', 'Stool']))
    ].copy()

    # AMR matrix
    amr_matrix = pd.read_csv(EXPLORATORY_DIR / "amr_rpm_matrix.tsv", sep="\t", index_col=0)

    # Diversity
    diversity = pd.read_csv(EXPLORATORY_DIR / "diversity_metrics.tsv", sep="\t")

    # Longitudinal results
    long_results = {}
    for site in ['axilla', 'groin', 'stool']:
        file_path = LONG_DIR / f"gene_longitudinal_{site}.tsv"
        if file_path.exists():
            long_results[site] = pd.read_csv(file_path, sep="\t")

    print(f"Loaded {len(metadata)} samples, {amr_matrix.shape[1]} genes")

    return metadata, amr_matrix, diversity, long_results

def figure1_pca_overview(metadata, amr_matrix):
    """Figure 1: PCA showing body site clustering"""
    print("\nCreating Figure 1: PCA overview...")

    # Run PCA
    scaler = StandardScaler()
    amr_scaled = scaler.fit_transform(amr_matrix)
    pca = PCA()
    pca_coords = pca.fit_transform(amr_scaled)

    # Create DataFrame
    pca_df = pd.DataFrame(pca_coords[:, :2], columns=['PC1', 'PC2'], index=amr_matrix.index)
    pca_df = pca_df.merge(metadata[['sample_name', 'Location', 'SampleType', 'SampleCollectionWeek']],
                          left_index=True, right_on='sample_name', how='left')

    # Create figure
    fig = plt.figure(figsize=(12, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.3)

    # Panel A: PCA by body site
    ax1 = fig.add_subplot(gs[0])
    for site in ['Axilla', 'Groin', 'Stool']:
        site_data = pca_df[pca_df['SampleType'] == site]
        ax1.scatter(site_data['PC1'], site_data['PC2'],
                   c=SITE_COLORS[site], label=site,
                   s=60, alpha=0.7, edgecolors='white', linewidths=0.5)

    ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=11, fontweight='bold')
    ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=11, fontweight='bold')
    ax1.set_title('A. Principal Component Analysis', fontsize=12, fontweight='bold', loc='left')
    ax1.legend(title='Body Site', frameon=True, loc='best')
    ax1.grid(True, alpha=0.3, linewidth=0.5)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Panel B: Variance explained
    ax2 = fig.add_subplot(gs[1])
    var_exp = pca.explained_variance_ratio_[:10] * 100
    ax2.bar(range(1, 11), var_exp, color='#2E86AB', alpha=0.7, edgecolor='black', linewidth=0.8)
    ax2.set_xlabel('Principal Component', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Variance Explained (%)', fontsize=11, fontweight='bold')
    ax2.set_title('B. Scree Plot', fontsize=12, fontweight='bold', loc='left')
    ax2.set_xticks(range(1, 11))
    ax2.grid(True, alpha=0.3, axis='y', linewidth=0.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.savefig(PUB_FIGURES_DIR / "Figure1_PCA.pdf", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {PUB_FIGURES_DIR / 'Figure1_PCA.pdf'}")
    plt.close()

def figure2_bodysite_comparison(metadata):
    """Figure 2: Body site AMR burden comparison"""
    print("\nCreating Figure 2: Body site comparison...")

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))

    # Top row: Total AMR burden by body site
    for i, week in enumerate(['Week.1', 'Week.3']):
        ax = axes[0, i]
        week_data = metadata[metadata['SampleCollectionWeek'] == week]

        data_to_plot = [week_data[week_data['SampleType'] == site]['total_amr_rpm'].values
                       for site in ['Axilla', 'Groin', 'Stool']]

        bp = ax.boxplot(data_to_plot, labels=['Axilla', 'Groin', 'Stool'],
                       patch_artist=True, widths=0.6)

        for patch, site in zip(bp['boxes'], ['Axilla', 'Groin', 'Stool']):
            patch.set_facecolor(SITE_COLORS[site])
            patch.set_alpha(0.7)

        for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(bp[element], color='black', linewidth=1.2)

        ax.set_ylabel('Total AMR RPM', fontsize=11, fontweight='bold')
        title = f'{"A" if i == 0 else "B"}. {week.replace(".", " ")}'
        ax.set_title(title, fontsize=12, fontweight='bold', loc='left')
        ax.grid(True, alpha=0.3, axis='y', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(bottom=0)

    # Third panel: Location comparison
    ax = axes[0, 2]
    data_to_plot = [metadata[metadata['Location'] == loc]['total_amr_rpm'].values
                   for loc in ['UCMC', 'ZCH']]

    bp = ax.boxplot(data_to_plot, labels=['UCMC', 'ZCH'],
                   patch_artist=True, widths=0.5)

    for patch, loc in zip(bp['boxes'], ['UCMC', 'ZCH']):
        patch.set_facecolor(LOCATION_COLORS[loc])
        patch.set_alpha(0.7)

    for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color='black', linewidth=1.2)

    ax.set_ylabel('Total AMR RPM', fontsize=11, fontweight='bold')
    ax.set_title('C. By Location', fontsize=12, fontweight='bold', loc='left')
    ax.grid(True, alpha=0.3, axis='y', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(bottom=0)

    # Bottom row: Gene richness
    for i, week in enumerate(['Week.1', 'Week.3']):
        ax = axes[1, i]
        week_data = metadata[metadata['SampleCollectionWeek'] == week]

        data_to_plot = [week_data[week_data['SampleType'] == site]['unique_amr_genes'].values
                       for site in ['Axilla', 'Groin', 'Stool']]

        bp = ax.boxplot(data_to_plot, labels=['Axilla', 'Groin', 'Stool'],
                       patch_artist=True, widths=0.6)

        for patch, site in zip(bp['boxes'], ['Axilla', 'Groin', 'Stool']):
            patch.set_facecolor(SITE_COLORS[site])
            patch.set_alpha(0.7)

        for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(bp[element], color='black', linewidth=1.2)

        ax.set_ylabel('AMR Gene Richness', fontsize=11, fontweight='bold')
        title = f'{"D" if i == 0 else "E"}. {week.replace(".", " ")}'
        ax.set_title(title, fontsize=12, fontweight='bold', loc='left')
        ax.grid(True, alpha=0.3, axis='y', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(bottom=0)

    # Bottom right: Sample size table
    ax = axes[1, 2]
    ax.axis('off')

    # Create sample size summary
    summary_data = []
    for loc in ['UCMC', 'ZCH']:
        for site in ['Axilla', 'Groin', 'Stool']:
            n = len(metadata[(metadata['Location'] == loc) & (metadata['SampleType'] == site)])
            summary_data.append([loc, site, n])

    summary_df = pd.DataFrame(summary_data, columns=['Location', 'Site', 'n'])
    summary_table = summary_df.pivot(index='Site', columns='Location', values='n')

    table = ax.table(cellText=summary_table.values,
                    rowLabels=summary_table.index,
                    colLabels=summary_table.columns,
                    cellLoc='center',
                    loc='center',
                    bbox=[0.1, 0.3, 0.8, 0.5])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Format table cells
    for key, cell in table.get_celld().items():
        cell.set_edgecolor('black')
        cell.set_linewidth(1)

    ax.set_title('F. Sample Sizes', fontsize=12, fontweight='bold', loc='left')

    plt.tight_layout()
    plt.savefig(PUB_FIGURES_DIR / "Figure2_BodySite_Comparison.pdf", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {PUB_FIGURES_DIR / 'Figure2_BodySite_Comparison.pdf'}")
    plt.close()

def figure3_longitudinal_trajectories(metadata):
    """Figure 3: Longitudinal trajectories Week 1→3"""
    print("\nCreating Figure 3: Longitudinal trajectories...")

    complete = metadata[metadata['subject_complete']].copy()

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, site in enumerate(['Axilla', 'Groin', 'Stool']):
        ax = axes[i]
        site_data = complete[complete['SampleType'] == site]

        subjects = site_data['SubjectID'].unique()

        # Plot trajectories
        for subj in subjects:
            subj_data = site_data[site_data['SubjectID'] == subj].sort_values('SampleCollectionWeek')

            if len(subj_data) == 2:
                weeks = [1 if w == 'Week.1' else 3 for w in subj_data['SampleCollectionWeek']]
                amr_values = subj_data['total_amr_rpm'].values

                color = LOCATION_COLORS[subj_data['Location'].iloc[0]]
                ax.plot(weeks, amr_values, 'o-', color=color,
                       alpha=0.4, linewidth=1.5, markersize=5)

        # Add mean trajectory
        week1_mean = site_data[site_data['SampleCollectionWeek'] == 'Week.1']['total_amr_rpm'].mean()
        week3_mean = site_data[site_data['SampleCollectionWeek'] == 'Week.3']['total_amr_rpm'].mean()
        ax.plot([1, 3], [week1_mean, week3_mean], 'k-', linewidth=3, marker='D',
               markersize=8, label='Mean', zorder=100)

        ax.set_xlabel('Week', fontsize=11, fontweight='bold')
        ax.set_ylabel('Total AMR RPM', fontsize=11, fontweight='bold')

        # Calculate p-value
        week1_vals = site_data[site_data['SampleCollectionWeek'] == 'Week.1'].set_index('SubjectID')['total_amr_rpm']
        week3_vals = site_data[site_data['SampleCollectionWeek'] == 'Week.3'].set_index('SubjectID')['total_amr_rpm']
        common = week1_vals.index.intersection(week3_vals.index)

        if len(common) >= 3:
            stat, pval = wilcoxon(week1_vals.loc[common], week3_vals.loc[common])
            pval_text = f'p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'
        else:
            pval_text = 'p = N/A'

        panel_label = chr(65 + i)  # A, B, C
        ax.set_title(f'{panel_label}. {site} (n={len(subjects)})\n{pval_text}',
                    fontsize=12, fontweight='bold', loc='left')
        ax.set_xticks([1, 3])
        ax.set_xticklabels(['Week 1', 'Week 3'])
        ax.grid(True, alpha=0.3, linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Legend
        if i == 0:
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color=LOCATION_COLORS['UCMC'], linewidth=2, label='UCMC'),
                Line2D([0], [0], color=LOCATION_COLORS['ZCH'], linewidth=2, label='ZCH'),
                Line2D([0], [0], color='black', linewidth=3, marker='D', label='Mean')
            ]
            ax.legend(handles=legend_elements, loc='upper left', frameon=True)

    plt.tight_layout()
    plt.savefig(PUB_FIGURES_DIR / "Figure3_Longitudinal_Trajectories.pdf", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {PUB_FIGURES_DIR / 'Figure3_Longitudinal_Trajectories.pdf'}")
    plt.close()

def figure4_volcano_plots(long_results):
    """Figure 4: Volcano plots of gene-level changes"""
    print("\nCreating Figure 4: Volcano plots...")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (site, results) in enumerate([('axilla', long_results.get('axilla')),
                                          ('groin', long_results.get('groin')),
                                          ('stool', long_results.get('stool'))]):
        ax = axes[i]

        if results is None or len(results) == 0:
            ax.text(0.5, 0.5, f'No data for {site.title()}',
                   ha='center', va='center', transform=ax.transAxes)
            continue

        # Prepare data
        results = results.copy()
        results['-log10_p'] = -np.log10(results['pvalue'] + 1e-300)

        # Color by significance
        colors = []
        for _, row in results.iterrows():
            if row['fdr'] < 0.05 and row['log2_fold_change'] > 0:
                colors.append('#D62828')  # Increased
            elif row['fdr'] < 0.05 and row['log2_fold_change'] < 0:
                colors.append('#2E86AB')  # Decreased
            else:
                colors.append('#CCCCCC')  # Not significant

        # Plot
        ax.scatter(results['log2_fold_change'], results['-log10_p'],
                  c=colors, s=30, alpha=0.6, edgecolors='none')

        # Add significance thresholds
        ax.axhline(y=-np.log10(0.05), color='black', linestyle='--',
                  linewidth=1, alpha=0.5, label='p = 0.05')
        ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8, alpha=0.5)

        # Labels
        ax.set_xlabel('log2 Fold Change (Week 3 / Week 1)', fontsize=11, fontweight='bold')
        ax.set_ylabel('-log10(p-value)', fontsize=11, fontweight='bold')

        n_sig = (results['fdr'] < 0.05).sum()
        n_up = ((results['fdr'] < 0.05) & (results['log2_fold_change'] > 0)).sum()
        n_down = ((results['fdr'] < 0.05) & (results['log2_fold_change'] < 0)).sum()

        panel_label = chr(65 + i)
        ax.set_title(f'{panel_label}. {site.title()}\n{n_sig} significant ({n_up}↑, {n_down}↓)',
                    fontsize=12, fontweight='bold', loc='left')
        ax.grid(True, alpha=0.3, linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Annotate top genes
        if n_sig > 0:
            top_genes = results[results['fdr'] < 0.05].nlargest(3, '-log10_p')
            for _, row in top_genes.iterrows():
                ax.annotate(row['gene'], (row['log2_fold_change'], row['-log10_p']),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8, alpha=0.8)

    plt.tight_layout()
    plt.savefig(PUB_FIGURES_DIR / "Figure4_Volcano_Plots.pdf", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {PUB_FIGURES_DIR / 'Figure4_Volcano_Plots.pdf'}")
    plt.close()

def main():
    print("="*60)
    print("GENERATING PUBLICATION-READY FIGURES")
    print("="*60)

    # Load data
    metadata, amr_matrix, diversity, long_results = load_data()

    # Generate figures
    figure1_pca_overview(metadata, amr_matrix)
    figure2_bodysite_comparison(metadata)
    figure3_longitudinal_trajectories(metadata)
    figure4_volcano_plots(long_results)

    print("\n" + "="*60)
    print("PUBLICATION FIGURES COMPLETE")
    print("="*60)
    print("\nGenerated figures:")
    print(f"  - {PUB_FIGURES_DIR / 'Figure1_PCA.pdf'}")
    print(f"  - {PUB_FIGURES_DIR / 'Figure2_BodySite_Comparison.pdf'}")
    print(f"  - {PUB_FIGURES_DIR / 'Figure3_Longitudinal_Trajectories.pdf'}")
    print(f"  - {PUB_FIGURES_DIR / 'Figure4_Volcano_Plots.pdf'}")

if __name__ == "__main__":
    main()
