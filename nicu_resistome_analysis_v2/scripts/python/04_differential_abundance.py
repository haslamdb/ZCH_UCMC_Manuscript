#!/usr/bin/env python3
"""
Differential Abundance Analysis:
Test UCMC vs ZCH separately for each body site

1. Mann-Whitney U test for each gene
2. Calculate fold change and effect sizes
3. FDR correction for multiple testing
4. Generate volcano plots and heatmaps
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Define paths
PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
RESULTS_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results"
EXPLORATORY_DIR = RESULTS_DIR / "exploratory"
DIFF_DIR = RESULTS_DIR / "differential"
FIGURES_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "figures" / "differential"
QC_DIR = RESULTS_DIR / "qc"

# Create directories
DIFF_DIR.mkdir(parents=True, exist_ok=True)
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

    return amr_matrix, metadata_with_amr

def test_differential_abundance(amr_matrix, metadata, body_site):
    """
    Test UCMC vs ZCH for a specific body site
    """
    print(f"\n{'='*60}")
    print(f"DIFFERENTIAL ABUNDANCE: {body_site}")
    print('='*60)

    # Filter to specific body site
    site_metadata = metadata[metadata['SampleType'] == body_site].copy()
    site_samples = site_metadata['sample_name'].values
    site_amr = amr_matrix.loc[site_samples]

    # Get UCMC and ZCH samples
    ucmc_samples = site_metadata[site_metadata['Location'] == 'UCMC']['sample_name'].values
    zch_samples = site_metadata[site_metadata['Location'] == 'ZCH']['sample_name'].values

    print(f"UCMC samples: {len(ucmc_samples)}")
    print(f"ZCH samples: {len(zch_samples)}")

    # Test each gene
    results = []

    for gene in amr_matrix.columns:
        ucmc_values = site_amr.loc[ucmc_samples, gene].values
        zch_values = site_amr.loc[zch_samples, gene].values

        # Calculate statistics
        ucmc_mean = ucmc_values.mean()
        zch_mean = zch_values.mean()
        ucmc_median = np.median(ucmc_values)
        zch_median = np.median(zch_values)

        # Log2 fold change (add pseudocount to avoid log(0))
        pseudocount = 0.1
        log2fc = np.log2((zch_mean + pseudocount) / (ucmc_mean + pseudocount))

        # Mann-Whitney U test
        if len(ucmc_values) > 0 and len(zch_values) > 0:
            try:
                stat, pval = mannwhitneyu(ucmc_values, zch_values, alternative='two-sided')
            except:
                pval = 1.0
                stat = 0
        else:
            pval = 1.0
            stat = 0

        # Prevalence (% of samples with gene detected)
        ucmc_prevalence = (ucmc_values > 0).sum() / len(ucmc_values) * 100
        zch_prevalence = (zch_values > 0).sum() / len(zch_values) * 100

        results.append({
            'gene': gene,
            'ucmc_mean': ucmc_mean,
            'ucmc_median': ucmc_median,
            'ucmc_prevalence': ucmc_prevalence,
            'zch_mean': zch_mean,
            'zch_median': zch_median,
            'zch_prevalence': zch_prevalence,
            'log2_fold_change': log2fc,
            'pvalue': pval,
            'mann_whitney_stat': stat
        })

    # Create results dataframe
    results_df = pd.DataFrame(results)

    # FDR correction
    results_df['fdr'] = multipletests(results_df['pvalue'], method='fdr_bh')[1]

    # Sort by p-value
    results_df = results_df.sort_values('pvalue')

    # Significance flags
    results_df['significant_p005'] = results_df['pvalue'] < 0.05
    results_df['significant_fdr005'] = results_df['fdr'] < 0.05
    results_df['significant_fdr01'] = results_df['fdr'] < 0.1

    # Save results
    output_file = DIFF_DIR / f"differential_abundance_{body_site.lower()}.tsv"
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"✓ Saved: {output_file}")

    # Summary
    n_sig_p = (results_df['pvalue'] < 0.05).sum()
    n_sig_fdr = (results_df['fdr'] < 0.05).sum()
    n_sig_fdr10 = (results_df['fdr'] < 0.1).sum()

    print(f"\nSignificant genes:")
    print(f"  p < 0.05: {n_sig_p}/{len(results_df)}")
    print(f"  FDR < 0.05: {n_sig_fdr}/{len(results_df)}")
    print(f"  FDR < 0.10: {n_sig_fdr10}/{len(results_df)}")

    if n_sig_fdr > 0:
        print(f"\nTop significant genes (FDR < 0.05):")
        top_genes = results_df[results_df['fdr'] < 0.05].head(10)
        for _, row in top_genes.iterrows():
            direction = "higher in ZCH" if row['log2_fold_change'] > 0 else "higher in UCMC"
            print(f"  {row['gene']}: log2FC = {row['log2_fold_change']:.2f} ({direction}), FDR = {row['fdr']:.4f}")

    return results_df

def plot_volcano(results_df, body_site):
    """Create volcano plot"""
    print(f"\nCreating volcano plot for {body_site}...")

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))

    # -log10(p-value)
    results_df['-log10_pvalue'] = -np.log10(results_df['pvalue'])

    # Color by significance
    colors = []
    for _, row in results_df.iterrows():
        if row['fdr'] < 0.05:
            if row['log2_fold_change'] > 0:
                colors.append('#d62728')  # Red: higher in ZCH
            else:
                colors.append('#1f77b4')  # Blue: higher in UCMC
        else:
            colors.append('#7f7f7f')  # Gray: not significant

    # Scatter plot
    ax.scatter(
        results_df['log2_fold_change'],
        results_df['-log10_pvalue'],
        c=colors,
        s=50,
        alpha=0.6,
        edgecolors='black',
        linewidths=0.5
    )

    # Add threshold lines
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1, alpha=0.5, label='p = 0.05')
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.5)

    # Label top genes
    top_genes = results_df.nsmallest(10, 'fdr')
    for _, row in top_genes.iterrows():
        ax.annotate(
            row['gene'],
            (row['log2_fold_change'], row['-log10_pvalue']),
            fontsize=8,
            alpha=0.7
        )

    ax.set_xlabel('log2(Fold Change) [ZCH/UCMC]')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'Volcano Plot: UCMC vs ZCH ({body_site})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d62728', label='Higher in ZCH (FDR < 0.05)'),
        Patch(facecolor='#1f77b4', label='Higher in UCMC (FDR < 0.05)'),
        Patch(facecolor='#7f7f7f', label='Not significant')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / f"volcano_plot_{body_site.lower()}.pdf", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / f'volcano_plot_{body_site.lower()}.pdf'}")
    plt.close()

def plot_heatmap(amr_matrix, metadata, results_df, body_site, top_n=50):
    """Create heatmap of top differential genes"""
    print(f"\nCreating heatmap for {body_site}...")

    # Get top genes by FDR
    top_genes = results_df.nsmallest(top_n, 'fdr')['gene'].values

    # Filter to body site
    site_metadata = metadata[metadata['SampleType'] == body_site].copy()
    site_samples = site_metadata['sample_name'].values
    site_amr = amr_matrix.loc[site_samples, top_genes]

    # Log-transform
    site_amr_log = np.log10(site_amr + 0.1)

    # Sort samples by location
    site_metadata = site_metadata.set_index('sample_name')
    site_metadata = site_metadata.loc[site_samples]
    site_metadata = site_metadata.sort_values('Location')
    site_amr_log = site_amr_log.loc[site_metadata.index]

    # Plot without col_colors to avoid indexing issues
    g = sns.clustermap(
        site_amr_log.T,  # Transpose: genes as rows, samples as columns
        cmap='viridis',
        row_cluster=True,
        col_cluster=True,  # Allow clustering
        figsize=(16, 12),
        yticklabels=True,
        xticklabels=False,
        cbar_kws={'label': 'log10(RPM + 0.1)'}
    )

    g.fig.suptitle(f'Top {top_n} Differential AMR Genes: UCMC vs ZCH ({body_site})\n(Samples clustered)', y=0.98)

    plt.savefig(FIGURES_DIR / f"heatmap_{body_site.lower()}_top{top_n}.pdf", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / f'heatmap_{body_site.lower()}_top{top_n}.pdf'}")
    plt.close()

def create_summary_table(results_all):
    """Create summary table across all body sites"""
    print("\nCreating summary table...")

    summary = []
    for site, results_df in results_all.items():
        n_genes = len(results_df)
        n_sig_p = (results_df['pvalue'] < 0.05).sum()
        n_sig_fdr05 = (results_df['fdr'] < 0.05).sum()
        n_sig_fdr10 = (results_df['fdr'] < 0.1).sum()

        # Genes higher in ZCH vs UCMC
        zch_higher = results_df[(results_df['fdr'] < 0.05) & (results_df['log2_fold_change'] > 0)]
        ucmc_higher = results_df[(results_df['fdr'] < 0.05) & (results_df['log2_fold_change'] < 0)]

        summary.append({
            'body_site': site,
            'total_genes': n_genes,
            'sig_p005': n_sig_p,
            'sig_fdr005': n_sig_fdr05,
            'sig_fdr010': n_sig_fdr10,
            'higher_in_zch': len(zch_higher),
            'higher_in_ucmc': len(ucmc_higher)
        })

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(DIFF_DIR / "differential_abundance_summary.tsv", sep="\t", index=False)
    print(f"✓ Saved: {DIFF_DIR / 'differential_abundance_summary.tsv'}")

    print("\n" + "="*60)
    print("SUMMARY TABLE")
    print("="*60)
    print(summary_df.to_string(index=False))

    return summary_df

def main():
    print("="*60)
    print("DIFFERENTIAL ABUNDANCE ANALYSIS")
    print("="*60)

    # Load data
    amr_matrix, metadata = load_data()

    # Test each body site
    body_sites = ['Axilla', 'Groin', 'Stool']
    results_all = {}

    for site in body_sites:
        # Differential abundance testing
        results_df = test_differential_abundance(amr_matrix, metadata, site)
        results_all[site] = results_df

        # Volcano plot
        plot_volcano(results_df, site)

        # Heatmap
        plot_heatmap(amr_matrix, metadata, results_df, site, top_n=50)

    # Summary table
    summary_df = create_summary_table(results_all)

    print("\n" + "="*60)
    print("DIFFERENTIAL ABUNDANCE ANALYSIS COMPLETE")
    print("="*60)
    print("\nOutputs:")
    for site in body_sites:
        print(f"  - {DIFF_DIR / f'differential_abundance_{site.lower()}.tsv'}")
    print(f"  - {DIFF_DIR / 'differential_abundance_summary.tsv'}")
    print("\nFigures:")
    for site in body_sites:
        print(f"  - {FIGURES_DIR / f'volcano_plot_{site.lower()}.pdf'}")
        print(f"  - {FIGURES_DIR / f'heatmap_{site.lower()}_top50.pdf'}")

if __name__ == "__main__":
    main()
