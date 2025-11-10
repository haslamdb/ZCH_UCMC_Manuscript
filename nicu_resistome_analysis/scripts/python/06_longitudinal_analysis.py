#!/usr/bin/env python3
"""
Longitudinal Analysis: Week 1 → Week 3 changes

Analyses:
1. Paired comparisons within subjects (complete subjects only)
2. Test for changes in total AMR burden over time
3. Identify genes that increase/decrease from Week 1 to Week 3
4. Stratify by location and body site
5. Create trajectory plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import wilcoxon, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Define paths
PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
RESULTS_DIR = PROJECT_ROOT / "nicu_resistome_analysis" / "results"
EXPLORATORY_DIR = RESULTS_DIR / "exploratory"
LONG_DIR = RESULTS_DIR / "longitudinal"
FIGURES_DIR = PROJECT_ROOT / "nicu_resistome_analysis" / "figures" / "longitudinal"
QC_DIR = RESULTS_DIR / "qc"

# Create directories
LONG_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Set plotting style
sns.set_style("whitegrid")
sns.set_context("talk")

def load_data():
    """Load AMR matrix and metadata"""
    print("Loading data...")

    # Load AMR matrix
    amr_matrix = pd.read_csv(EXPLORATORY_DIR / "amr_rpm_matrix.tsv", sep="\t", index_col=0)

    # Load metadata
    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")

    # Filter to samples with AMR data
    metadata_with_amr = metadata[metadata['sample_name'].isin(amr_matrix.index)].copy()

    print(f"AMR matrix: {amr_matrix.shape}")
    print(f"Metadata: {len(metadata_with_amr)} samples")
    print(f"Complete subjects: {metadata_with_amr[metadata_with_amr['subject_complete']]['SubjectID'].nunique()}")

    return amr_matrix, metadata_with_amr

def test_total_amr_change(metadata):
    """Test for changes in total AMR burden Week 1 → Week 3"""
    print("\n" + "="*60)
    print("TOTAL AMR BURDEN: Week 1 → Week 3")
    print("="*60)

    # Filter to complete subjects only (paired data)
    complete = metadata[metadata['subject_complete']].copy()

    results = []

    for site in ['Axilla', 'Groin', 'Stool']:
        site_data = complete[complete['SampleType'] == site]

        # Get Week 1 and Week 3 data
        week1 = site_data[site_data['SampleCollectionWeek'] == 'Week.1'].set_index('SubjectID')['total_amr_rpm']
        week3 = site_data[site_data['SampleCollectionWeek'] == 'Week.3'].set_index('SubjectID')['total_amr_rpm']

        # Only subjects with both timepoints
        common_subjects = week1.index.intersection(week3.index)
        week1_paired = week1.loc[common_subjects]
        week3_paired = week3.loc[common_subjects]

        if len(common_subjects) < 3:
            print(f"\n{site}: Not enough paired samples (n={len(common_subjects)})")
            continue

        # Wilcoxon signed-rank test (paired)
        stat, pval = wilcoxon(week1_paired, week3_paired)

        # Calculate statistics
        mean_week1 = week1_paired.mean()
        mean_week3 = week3_paired.mean()
        median_week1 = week1_paired.median()
        median_week3 = week3_paired.median()
        mean_change = week3_paired.mean() - week1_paired.mean()
        pct_change = (mean_change / mean_week1) * 100 if mean_week1 > 0 else 0

        print(f"\n{site} (n={len(common_subjects)} subjects):")
        print(f"  Week 1: mean={mean_week1:.1f}, median={median_week1:.1f}")
        print(f"  Week 3: mean={mean_week3:.1f}, median={median_week3:.1f}")
        print(f"  Change: {mean_change:+.1f} ({pct_change:+.1f}%)")
        print(f"  Wilcoxon p = {pval:.4f}")

        results.append({
            'body_site': site,
            'n_subjects': len(common_subjects),
            'week1_mean': mean_week1,
            'week1_median': median_week1,
            'week3_mean': mean_week3,
            'week3_median': median_week3,
            'mean_change': mean_change,
            'percent_change': pct_change,
            'wilcoxon_pvalue': pval
        })

    # Save results
    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        results_df['fdr'] = multipletests(results_df['wilcoxon_pvalue'], method='fdr_bh')[1]
        results_df.to_csv(LONG_DIR / "total_amr_longitudinal.tsv", sep="\t", index=False)
        print(f"\n✓ Saved: {LONG_DIR / 'total_amr_longitudinal.tsv'}")

    return results_df

def test_gene_level_changes(amr_matrix, metadata):
    """Test each gene for Week 1 → Week 3 changes"""
    print("\n" + "="*60)
    print("GENE-LEVEL LONGITUDINAL CHANGES")
    print("="*60)

    # Filter to complete subjects
    complete = metadata[metadata['subject_complete']].copy()

    all_results = {}

    for site in ['Axilla', 'Groin', 'Stool']:
        print(f"\nTesting {site}...")

        site_data = complete[complete['SampleType'] == site]

        # Get Week 1 and Week 3 data with subject mapping
        week1_data = site_data[site_data['SampleCollectionWeek'] == 'Week.1'][['sample_name', 'SubjectID']].drop_duplicates('sample_name')
        week3_data = site_data[site_data['SampleCollectionWeek'] == 'Week.3'][['sample_name', 'SubjectID']].drop_duplicates('sample_name')

        # Filter to samples in AMR matrix
        week1_data = week1_data[week1_data['sample_name'].isin(amr_matrix.index)]
        week3_data = week3_data[week3_data['sample_name'].isin(amr_matrix.index)]

        # Test each gene
        results = []

        for gene in amr_matrix.columns:
            # Get Week 1 values and map to subjects
            week1_gene = amr_matrix.loc[week1_data['sample_name'], gene].copy()
            week1_gene.index = week1_data.set_index('sample_name').loc[week1_gene.index, 'SubjectID']

            # Get Week 3 values and map to subjects
            week3_gene = amr_matrix.loc[week3_data['sample_name'], gene].copy()
            week3_gene.index = week3_data.set_index('sample_name').loc[week3_gene.index, 'SubjectID']

            # Only paired subjects
            common_subjects = week1_gene.index.intersection(week3_gene.index)

            if len(common_subjects) < 3:
                continue

            week1_paired = week1_gene.loc[common_subjects]
            week3_paired = week3_gene.loc[common_subjects]

            # Wilcoxon test
            try:
                stat, pval = wilcoxon(week1_paired, week3_paired)
            except:
                pval = 1.0

            # Statistics
            mean_week1 = week1_paired.mean()
            mean_week3 = week3_paired.mean()
            log2fc = np.log2((mean_week3 + 0.1) / (mean_week1 + 0.1))

            results.append({
                'gene': gene,
                'n_subjects': len(common_subjects),
                'week1_mean': mean_week1,
                'week3_mean': mean_week3,
                'log2_fold_change': log2fc,
                'pvalue': pval
            })

        # Create results dataframe
        results_df = pd.DataFrame(results)

        if len(results_df) > 0:
            # FDR correction
            results_df['fdr'] = multipletests(results_df['pvalue'], method='fdr_bh')[1]

            # Sort by p-value
            results_df = results_df.sort_values('pvalue')

            # Save
            output_file = LONG_DIR / f"gene_longitudinal_{site.lower()}.tsv"
            results_df.to_csv(output_file, sep="\t", index=False)

            # Summary
            n_sig = (results_df['fdr'] < 0.05).sum()
            print(f"  Significant genes (FDR < 0.05): {n_sig}/{len(results_df)}")

            if n_sig > 0:
                print(f"  Top changing genes:")
                top_genes = results_df[results_df['fdr'] < 0.05].head(5)
                for _, row in top_genes.iterrows():
                    direction = "increased" if row['log2_fold_change'] > 0 else "decreased"
                    print(f"    {row['gene']}: {direction} (log2FC={row['log2_fold_change']:.2f}, FDR={row['fdr']:.4f})")

            all_results[site] = results_df

    return all_results

def plot_longitudinal_trajectories(metadata):
    """Plot AMR burden trajectories for complete subjects"""
    print("\nCreating trajectory plots...")

    complete = metadata[metadata['subject_complete']].copy()

    # Figure: Trajectories by body site
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    for i, site in enumerate(['Axilla', 'Groin', 'Stool']):
        ax = axes[i]
        site_data = complete[complete['SampleType'] == site]

        # Get unique subjects
        subjects = site_data['SubjectID'].unique()

        # Plot trajectory for each subject
        for subj in subjects:
            subj_data = site_data[site_data['SubjectID'] == subj].sort_values('SampleCollectionWeek')

            if len(subj_data) == 2:  # Has both timepoints
                weeks = [1 if w == 'Week.1' else 3 for w in subj_data['SampleCollectionWeek']]
                amr_values = subj_data['total_amr_rpm'].values

                # Color by location
                color = '#1f77b4' if subj_data['Location'].iloc[0] == 'UCMC' else '#ff7f0e'
                alpha = 0.3

                ax.plot(weeks, amr_values, 'o-', color=color, alpha=alpha, linewidth=1, markersize=6)

        ax.set_xlabel('Week')
        ax.set_ylabel('Total AMR RPM')
        ax.set_title(f'{site}\n(n={len(subjects)} subjects)')
        ax.set_xticks([1, 3])
        ax.set_xticklabels(['Week 1', 'Week 3'])
        ax.grid(True, alpha=0.3)

        # Legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='#1f77b4', linewidth=2, label='UCMC'),
            Line2D([0], [0], color='#ff7f0e', linewidth=2, label='ZCH')
        ]
        ax.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "longitudinal_trajectories.pdf", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'longitudinal_trajectories.pdf'}")
    plt.close()

def plot_paired_boxplots(metadata):
    """Boxplots showing Week 1 vs Week 3 for complete subjects"""
    print("\nCreating paired boxplots...")

    complete = metadata[metadata['subject_complete']].copy()

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Row 1: Total AMR burden
    for i, site in enumerate(['Axilla', 'Groin', 'Stool']):
        ax = axes[0, i]
        site_data = complete[complete['SampleType'] == site]

        data_to_plot = [
            site_data[site_data['SampleCollectionWeek'] == 'Week.1']['total_amr_rpm'],
            site_data[site_data['SampleCollectionWeek'] == 'Week.3']['total_amr_rpm']
        ]

        bp = ax.boxplot(data_to_plot, labels=['Week 1', 'Week 3'], patch_artist=True)
        bp['boxes'][0].set_facecolor('#2ca02c')
        bp['boxes'][1].set_facecolor('#d62728')

        ax.set_ylabel('Total AMR RPM')
        n_subj = len(site_data[site_data['SampleCollectionWeek'] == 'Week.1'])
        ax.set_title(f'{site}\n(n={n_subj} subjects)')
        ax.grid(True, alpha=0.3, axis='y')

        # Add p-value
        week1 = site_data[site_data['SampleCollectionWeek'] == 'Week.1'].set_index('SubjectID')['total_amr_rpm']
        week3 = site_data[site_data['SampleCollectionWeek'] == 'Week.3'].set_index('SubjectID')['total_amr_rpm']
        common = week1.index.intersection(week3.index)

        if len(common) >= 3:
            stat, pval = wilcoxon(week1.loc[common], week3.loc[common])
            ax.text(0.5, 0.95, f'p = {pval:.4f}',
                    transform=ax.transAxes, ha='center', va='top', fontsize=10)

    # Row 2: Gene richness
    for i, site in enumerate(['Axilla', 'Groin', 'Stool']):
        ax = axes[1, i]
        site_data = complete[complete['SampleType'] == site]

        data_to_plot = [
            site_data[site_data['SampleCollectionWeek'] == 'Week.1']['unique_amr_genes'],
            site_data[site_data['SampleCollectionWeek'] == 'Week.3']['unique_amr_genes']
        ]

        bp = ax.boxplot(data_to_plot, labels=['Week 1', 'Week 3'], patch_artist=True)
        bp['boxes'][0].set_facecolor('#2ca02c')
        bp['boxes'][1].set_facecolor('#d62728')

        ax.set_ylabel('Number of Unique AMR Genes')
        ax.set_title(f'{site}')
        ax.grid(True, alpha=0.3, axis='y')

        # Add p-value
        week1 = site_data[site_data['SampleCollectionWeek'] == 'Week.1'].set_index('SubjectID')['unique_amr_genes']
        week3 = site_data[site_data['SampleCollectionWeek'] == 'Week.3'].set_index('SubjectID')['unique_amr_genes']
        common = week1.index.intersection(week3.index)

        if len(common) >= 3:
            stat, pval = wilcoxon(week1.loc[common], week3.loc[common])
            ax.text(0.5, 0.95, f'p = {pval:.4f}',
                    transform=ax.transAxes, ha='center', va='top', fontsize=10)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "longitudinal_boxplots.pdf", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'longitudinal_boxplots.pdf'}")
    plt.close()

def main():
    print("="*60)
    print("LONGITUDINAL ANALYSIS: Week 1 → Week 3")
    print("="*60)

    # Load data
    amr_matrix, metadata = load_data()

    # Test total AMR burden changes
    total_amr_results = test_total_amr_change(metadata)

    # Test gene-level changes
    gene_results = test_gene_level_changes(amr_matrix, metadata)

    # Plot trajectories
    plot_longitudinal_trajectories(metadata)

    # Plot boxplots
    plot_paired_boxplots(metadata)

    print("\n" + "="*60)
    print("LONGITUDINAL ANALYSIS COMPLETE")
    print("="*60)
    print("\nOutputs:")
    print(f"  - {LONG_DIR / 'total_amr_longitudinal.tsv'}")
    for site in ['Axilla', 'Groin', 'Stool']:
        print(f"  - {LONG_DIR / f'gene_longitudinal_{site.lower()}.tsv'}")
    print("\nFigures:")
    print(f"  - {FIGURES_DIR / 'longitudinal_trajectories.pdf'}")
    print(f"  - {FIGURES_DIR / 'longitudinal_boxplots.pdf'}")

if __name__ == "__main__":
    main()
