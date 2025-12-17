#!/usr/bin/env python3
"""
Visualize ESBL-Antibiotic Associations

Creates publication-ready figures showing:
1. Correlation heatmap (antibiotics × ESBL gene families)
2. Scatter plots for key associations (CTX-M vs cephalosporins)
3. Box plots comparing exposed vs unexposed
4. Forest plot of effect sizes
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Define paths
DIAMOND_DIR = Path("/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis/diamond_traditional")
DATA_DIR = DIAMOND_DIR / "data"
RESULTS_DIR = DIAMOND_DIR / "results" / "esbl_analysis_tables"
FIGURES_DIR = DIAMOND_DIR / "figures" / "esbl_analysis"

def load_data():
    """Load ESBL abundance matrix and antibiotic data"""
    print("Loading data...")

    # Load ESBL abundance matrix
    esbl_matrix = pd.read_csv(RESULTS_DIR / "esbl_abundance_matrix.tsv", sep="\t")
    print(f"  ESBL matrix: {len(esbl_matrix)} samples")

    # Load master metadata with antibiotics
    master_meta = pd.read_csv(DATA_DIR / "master_metadata_fixed.tsv", sep="\t")
    print(f"  Master metadata: {len(master_meta)} samples")

    return esbl_matrix, master_meta

def prepare_antibiotic_data(esbl_matrix, master_meta):
    """Merge and calculate cumulative antibiotic exposure"""
    print("\nPreparing antibiotic exposure data...")

    # Fix ZCH sample names
    master_meta = master_meta.copy()
    zch_mask = master_meta['Location'] == 'ZCH'
    master_meta.loc[zch_mask, 'sample_name'] = 'ZJH_' + master_meta.loc[zch_mask, 'sample_name']

    # Get antibiotic columns
    abx_w1_cols = [c for c in master_meta.columns if c.endswith('_w1') and c != 'None_w1']
    abx_w2_cols = [c for c in master_meta.columns if c.endswith('_w2') and c != 'None_w2']

    # Merge
    keep_cols = ['sample_name', 'Location', 'SampleCollectionWeek', 'SampleType'] + abx_w1_cols + abx_w2_cols
    keep_cols = [c for c in keep_cols if c in master_meta.columns]

    merged = esbl_matrix.merge(
        master_meta[keep_cols],
        on='sample_name',
        how='inner',
        suffixes=('', '_meta')
    )
    print(f"  Merged samples: {len(merged)}")

    # Calculate cumulative exposure for beta-lactams
    beta_lactams = ['Ampicillin', 'Nafcillin', 'Penicillin', 'Meropenem',
                    'Cefotaxime', 'Moxalactam', 'Cefazolin']

    week_col = 'SampleCollectionWeek_meta' if 'SampleCollectionWeek_meta' in merged.columns else 'SampleCollectionWeek'

    for abx in beta_lactams:
        w1_col = f'{abx}_w1'
        w2_col = f'{abx}_w2'

        if w1_col not in merged.columns:
            continue

        merged[f'{abx}_cumulative'] = merged[w1_col].fillna(0)

        if w2_col in merged.columns:
            week3_mask = merged[week_col] == 'Week.3'
            merged.loc[week3_mask, f'{abx}_cumulative'] += merged.loc[week3_mask, w2_col].fillna(0)

    return merged, beta_lactams

def plot_correlation_heatmap(merged, beta_lactams):
    """Create heatmap of antibiotic-ESBL correlations"""
    print("\nCreating correlation heatmap...")

    gene_families = ['CTX-M_total', 'SHV_total', 'TEM_total', 'ESBL_total']
    gene_labels = ['CTX-M', 'SHV', 'TEM', 'Total ESBL']

    # Calculate correlations
    corr_matrix = []
    pval_matrix = []

    available_abx = []
    for abx in beta_lactams:
        if f'{abx}_cumulative' in merged.columns:
            abx_vals = merged[f'{abx}_cumulative']
            if abx_vals.std() > 0 and (abx_vals > 0).sum() >= 5:
                available_abx.append(abx)

    for abx in available_abx:
        row_corr = []
        row_pval = []
        abx_vals = merged[f'{abx}_cumulative']

        for gene in gene_families:
            if gene in merged.columns:
                rho, pval = spearmanr(abx_vals, merged[gene])
                row_corr.append(rho)
                row_pval.append(pval)
            else:
                row_corr.append(np.nan)
                row_pval.append(np.nan)

        corr_matrix.append(row_corr)
        pval_matrix.append(row_pval)

    corr_df = pd.DataFrame(corr_matrix, index=available_abx, columns=gene_labels)
    pval_df = pd.DataFrame(pval_matrix, index=available_abx, columns=gene_labels)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(8, 6))

    # Custom colormap centered at 0
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    sns.heatmap(corr_df, annot=True, fmt='.2f', cmap=cmap, center=0,
                vmin=-0.4, vmax=0.4, ax=ax,
                cbar_kws={'label': 'Spearman ρ'})

    # Add significance markers
    for i, abx in enumerate(available_abx):
        for j, gene in enumerate(gene_labels):
            pval = pval_df.loc[abx, gene]
            if pval < 0.001:
                ax.text(j + 0.5, i + 0.75, '***', ha='center', va='center',
                       fontsize=8, color='black')
            elif pval < 0.01:
                ax.text(j + 0.5, i + 0.75, '**', ha='center', va='center',
                       fontsize=8, color='black')
            elif pval < 0.05:
                ax.text(j + 0.5, i + 0.75, '*', ha='center', va='center',
                       fontsize=8, color='black')

    ax.set_title('Antibiotic Exposure vs ESBL Gene Abundance\n(Spearman correlations)',
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('ESBL Gene Family')
    ax.set_ylabel('Antibiotic (cumulative exposure)')

    # Add legend for significance
    ax.text(1.02, -0.15, '* p<0.05  ** p<0.01  *** p<0.001',
            transform=ax.transAxes, fontsize=8, va='top')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'esbl_antibiotic_correlation_heatmap.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'esbl_antibiotic_correlation_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: esbl_antibiotic_correlation_heatmap.pdf/png")

    return corr_df, pval_df

def plot_ctxm_cephalosporin_scatter(merged):
    """Scatter plots of CTX-M vs cephalosporin exposure"""
    print("\nCreating CTX-M vs cephalosporin scatter plots...")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    cephalosporins = [
        ('Moxalactam', 'Moxalactam\n(3rd gen cephalosporin)'),
        ('Cefotaxime', 'Cefotaxime\n(3rd gen cephalosporin)'),
        ('Penicillin', 'Penicillin')
    ]

    for ax, (abx, label) in zip(axes, cephalosporins):
        abx_col = f'{abx}_cumulative'
        if abx_col not in merged.columns:
            continue

        x = merged[abx_col]
        y = merged['CTX-M_total']

        # Add jitter to x for visualization
        x_jitter = x + np.random.normal(0, 0.1, len(x))

        # Color by location
        colors = merged['Location'].map({'UCMC': '#2ecc71', 'ZCH': '#e74c3c'})

        ax.scatter(x_jitter, y + 0.1, c=colors, alpha=0.5, s=30, edgecolors='none')

        # Calculate correlation
        mask = ~(np.isnan(x) | np.isnan(y))
        rho, pval = spearmanr(x[mask], y[mask])

        ax.set_xlabel(f'{label}\n(days of exposure)', fontsize=10)
        ax.set_ylabel('CTX-M Abundance (RPM)' if ax == axes[0] else '')
        ax.set_yscale('symlog', linthresh=1)

        # Add correlation text
        if pval < 0.001:
            pval_str = 'p < 0.001'
        else:
            pval_str = f'p = {pval:.3f}'
        ax.text(0.95, 0.95, f'ρ = {rho:.2f}\n{pval_str}',
                transform=ax.transAxes, ha='right', va='top',
                fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Add trend line
        if (x > 0).sum() >= 5:
            x_pos = x[x > 0]
            y_pos = y[x > 0]
            z = np.polyfit(x_pos, np.log1p(y_pos), 1)
            x_line = np.linspace(0, x.max(), 100)
            y_line = np.expm1(z[0] * x_line + z[1])
            ax.plot(x_line, y_line, 'k--', alpha=0.5, linewidth=2)

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#2ecc71',
               markersize=10, label='UCMC'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c',
               markersize=10, label='ZCH')
    ]
    axes[2].legend(handles=legend_elements, loc='lower right')

    plt.suptitle('CTX-M Gene Abundance vs Beta-Lactam Exposure',
                 fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'ctxm_vs_cephalosporins_scatter.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'ctxm_vs_cephalosporins_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: ctxm_vs_cephalosporins_scatter.pdf/png")

def plot_exposed_vs_unexposed_boxplots(merged):
    """Box plots comparing ESBL abundance in exposed vs unexposed infants"""
    print("\nCreating exposed vs unexposed box plots...")

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))

    antibiotics = ['Moxalactam', 'Cefotaxime', 'Penicillin',
                   'Ampicillin', 'Meropenem', 'Nafcillin']

    for ax, abx in zip(axes.flat, antibiotics):
        abx_col = f'{abx}_cumulative'
        if abx_col not in merged.columns:
            ax.set_visible(False)
            continue

        # Create exposure groups
        plot_data = merged.copy()
        plot_data['Exposure'] = np.where(plot_data[abx_col] > 0, 'Exposed', 'Unexposed')
        plot_data['CTX-M_plot'] = plot_data['CTX-M_total'].replace(0, 0.01)

        n_exposed = (plot_data['Exposure'] == 'Exposed').sum()
        n_unexposed = (plot_data['Exposure'] == 'Unexposed').sum()

        # Box plot
        sns.boxplot(data=plot_data, x='Exposure', y='CTX-M_plot', ax=ax,
                    order=['Unexposed', 'Exposed'], palette=['#3498db', '#e74c3c'])
        sns.stripplot(data=plot_data, x='Exposure', y='CTX-M_plot', ax=ax,
                      order=['Unexposed', 'Exposed'], color='black', alpha=0.3, size=3)

        ax.set_yscale('log')
        ax.set_xlabel('')
        ax.set_ylabel('CTX-M RPM' if ax in axes[:, 0] else '')
        ax.set_title(f'{abx}\n(n={n_unexposed} vs {n_exposed})', fontsize=10, fontweight='bold')

        # Statistical test
        exposed_vals = plot_data[plot_data['Exposure'] == 'Exposed']['CTX-M_total']
        unexposed_vals = plot_data[plot_data['Exposure'] == 'Unexposed']['CTX-M_total']

        if len(exposed_vals) >= 5 and len(unexposed_vals) >= 5:
            _, pval = mannwhitneyu(exposed_vals, unexposed_vals, alternative='two-sided')
            if pval < 0.001:
                pval_str = 'p < 0.001'
            elif pval < 0.05:
                pval_str = f'p = {pval:.3f}'
            else:
                pval_str = f'p = {pval:.2f}'
            ax.text(0.5, 0.95, pval_str, transform=ax.transAxes, ha='center', va='top', fontsize=9)

    plt.suptitle('CTX-M Abundance: Antibiotic Exposed vs Unexposed',
                 fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'ctxm_exposed_vs_unexposed.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'ctxm_exposed_vs_unexposed.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: ctxm_exposed_vs_unexposed.pdf/png")

def plot_forest_plot(merged, beta_lactams):
    """Forest plot showing effect sizes for all antibiotic-ESBL associations"""
    print("\nCreating forest plot...")

    results = []

    for abx in beta_lactams:
        abx_col = f'{abx}_cumulative'
        if abx_col not in merged.columns:
            continue

        abx_vals = merged[abx_col]
        if abx_vals.std() == 0 or (abx_vals > 0).sum() < 5:
            continue

        for gene, label in [('CTX-M_total', 'CTX-M'), ('ESBL_total', 'Total ESBL')]:
            if gene not in merged.columns:
                continue

            rho, pval = spearmanr(abx_vals, merged[gene])

            # Bootstrap CI for correlation
            n_boot = 1000
            boot_rhos = []
            for _ in range(n_boot):
                idx = np.random.choice(len(merged), len(merged), replace=True)
                boot_rho, _ = spearmanr(abx_vals.iloc[idx], merged[gene].iloc[idx])
                if not np.isnan(boot_rho):
                    boot_rhos.append(boot_rho)

            if len(boot_rhos) > 0:
                ci_low = np.percentile(boot_rhos, 2.5)
                ci_high = np.percentile(boot_rhos, 97.5)
            else:
                ci_low = ci_high = rho

            results.append({
                'Antibiotic': abx,
                'Gene': label,
                'rho': rho,
                'ci_low': ci_low,
                'ci_high': ci_high,
                'pval': pval,
                'significant': pval < 0.05
            })

    results_df = pd.DataFrame(results)

    # Create forest plot - CTX-M only
    ctxm_results = results_df[results_df['Gene'] == 'CTX-M'].copy()
    ctxm_results = ctxm_results.sort_values('rho', ascending=True)

    fig, ax = plt.subplots(figsize=(10, 6))

    y_pos = range(len(ctxm_results))

    # Plot points and error bars
    colors = ['#e74c3c' if sig else '#95a5a6' for sig in ctxm_results['significant']]

    ax.errorbar(ctxm_results['rho'], y_pos,
                xerr=[ctxm_results['rho'] - ctxm_results['ci_low'],
                      ctxm_results['ci_high'] - ctxm_results['rho']],
                fmt='o', capsize=4, capthick=2, markersize=8,
                color='black', ecolor='gray')

    # Color the markers
    for i, (idx, row) in enumerate(ctxm_results.iterrows()):
        color = '#e74c3c' if row['significant'] else '#95a5a6'
        ax.scatter(row['rho'], i, s=100, c=color, zorder=5, edgecolors='black')

    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(ctxm_results['Antibiotic'])
    ax.set_xlabel('Spearman ρ (with 95% CI)', fontsize=11)
    ax.set_ylabel('')
    ax.set_title('Association Between Antibiotic Exposure and CTX-M Abundance\n(Forest Plot)',
                 fontsize=12, fontweight='bold')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e74c3c', edgecolor='black', label='p < 0.05'),
        Patch(facecolor='#95a5a6', edgecolor='black', label='Not significant')
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    ax.set_xlim(-0.3, 0.5)
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'ctxm_antibiotic_forest_plot.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'ctxm_antibiotic_forest_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: ctxm_antibiotic_forest_plot.pdf/png")

    return results_df

def plot_combined_figure(merged):
    """Create a combined publication figure"""
    print("\nCreating combined publication figure...")

    fig = plt.figure(figsize=(14, 10))

    # Layout: 2x2
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # Panel A: CTX-M by location with antibiotic context
    ax1 = fig.add_subplot(gs[0, 0])

    plot_data = merged.copy()
    plot_data['CTX-M_plot'] = plot_data['CTX-M_total'].replace(0, 0.01)

    sns.boxplot(data=plot_data, x='Location', y='CTX-M_plot', ax=ax1,
                order=['UCMC', 'ZCH'], palette=['#2ecc71', '#e74c3c'])
    sns.stripplot(data=plot_data, x='Location', y='CTX-M_plot', ax=ax1,
                  order=['UCMC', 'ZCH'], color='black', alpha=0.3, size=3)

    ax1.set_yscale('log')
    ax1.set_ylabel('CTX-M Abundance (RPM)')
    ax1.set_xlabel('')
    ax1.set_title('A. CTX-M by Location', fontsize=11, fontweight='bold', loc='left')

    # Add sample sizes and p-value
    ucmc_n = (plot_data['Location'] == 'UCMC').sum()
    zch_n = (plot_data['Location'] == 'ZCH').sum()
    ax1.set_xticklabels([f'UCMC\n(n={ucmc_n})', f'ZCH\n(n={zch_n})'])

    ucmc_vals = plot_data[plot_data['Location'] == 'UCMC']['CTX-M_total']
    zch_vals = plot_data[plot_data['Location'] == 'ZCH']['CTX-M_total']
    _, pval = mannwhitneyu(ucmc_vals, zch_vals)
    ax1.text(0.5, 0.95, f'p = {pval:.2e}', transform=ax1.transAxes, ha='center', va='top', fontsize=10)

    # Panel B: Moxalactam scatter
    ax2 = fig.add_subplot(gs[0, 1])

    if 'Moxalactam_cumulative' in merged.columns:
        x = merged['Moxalactam_cumulative']
        y = merged['CTX-M_total']
        colors = merged['Location'].map({'UCMC': '#2ecc71', 'ZCH': '#e74c3c'})

        ax2.scatter(x + np.random.normal(0, 0.1, len(x)), y + 0.1,
                   c=colors, alpha=0.5, s=30, edgecolors='none')

        rho, pval = spearmanr(x, y)
        ax2.set_xlabel('Moxalactam Exposure (days)')
        ax2.set_ylabel('CTX-M Abundance (RPM)')
        ax2.set_yscale('symlog', linthresh=1)
        ax2.set_title('B. CTX-M vs Moxalactam', fontsize=11, fontweight='bold', loc='left')
        ax2.text(0.95, 0.95, f'ρ = {rho:.2f}\np < 0.001',
                transform=ax2.transAxes, ha='right', va='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Panel C: Exposed vs unexposed for Moxalactam
    ax3 = fig.add_subplot(gs[1, 0])

    if 'Moxalactam_cumulative' in merged.columns:
        plot_data['Moxalactam_Exposed'] = np.where(
            plot_data['Moxalactam_cumulative'] > 0, 'Exposed', 'Unexposed'
        )

        sns.boxplot(data=plot_data, x='Moxalactam_Exposed', y='CTX-M_plot', ax=ax3,
                    order=['Unexposed', 'Exposed'], palette=['#3498db', '#e74c3c'])
        sns.stripplot(data=plot_data, x='Moxalactam_Exposed', y='CTX-M_plot', ax=ax3,
                      order=['Unexposed', 'Exposed'], color='black', alpha=0.3, size=3)

        ax3.set_yscale('log')
        ax3.set_xlabel('')
        ax3.set_ylabel('CTX-M Abundance (RPM)')
        ax3.set_title('C. CTX-M: Moxalactam Exposed vs Unexposed', fontsize=11, fontweight='bold', loc='left')

        n_unexp = (plot_data['Moxalactam_Exposed'] == 'Unexposed').sum()
        n_exp = (plot_data['Moxalactam_Exposed'] == 'Exposed').sum()
        ax3.set_xticklabels([f'Unexposed\n(n={n_unexp})', f'Exposed\n(n={n_exp})'])

        exp_vals = plot_data[plot_data['Moxalactam_Exposed'] == 'Exposed']['CTX-M_total']
        unexp_vals = plot_data[plot_data['Moxalactam_Exposed'] == 'Unexposed']['CTX-M_total']
        _, pval = mannwhitneyu(exp_vals, unexp_vals)
        ax3.text(0.5, 0.95, f'p = {pval:.2e}', transform=ax3.transAxes, ha='center', va='top', fontsize=10)

    # Panel D: Body site comparison
    ax4 = fig.add_subplot(gs[1, 1])

    site_col = 'SampleType' if 'SampleType' in plot_data.columns else 'SampleType_meta'
    if site_col in plot_data.columns:
        sns.boxplot(data=plot_data, x=site_col, y='CTX-M_plot', ax=ax4,
                    order=['Axilla', 'Groin', 'Stool'], palette='Set2')
        sns.stripplot(data=plot_data, x=site_col, y='CTX-M_plot', ax=ax4,
                      order=['Axilla', 'Groin', 'Stool'], color='black', alpha=0.3, size=3)

        ax4.set_yscale('log')
        ax4.set_xlabel('')
        ax4.set_ylabel('CTX-M Abundance (RPM)')
        ax4.set_title('D. CTX-M by Body Site', fontsize=11, fontweight='bold', loc='left')

        # Kruskal-Wallis
        from scipy.stats import kruskal
        groups = [plot_data[plot_data[site_col] == site]['CTX-M_total']
                  for site in ['Axilla', 'Groin', 'Stool']]
        _, pval = kruskal(*groups)
        ax4.text(0.5, 0.95, f'Kruskal-Wallis p = {pval:.2e}',
                transform=ax4.transAxes, ha='center', va='top', fontsize=10)

    plt.suptitle('CTX-M ESBL Gene Analysis: Location, Antibiotic Exposure, and Body Site',
                 fontsize=13, fontweight='bold', y=0.98)

    plt.savefig(FIGURES_DIR / 'ctxm_combined_figure.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'ctxm_combined_figure.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: ctxm_combined_figure.pdf/png")

def main():
    print("="*60)
    print("ESBL-ANTIBIOTIC ASSOCIATION VISUALIZATION")
    print("="*60)

    # Load data
    esbl_matrix, master_meta = load_data()

    # Prepare antibiotic data
    merged, beta_lactams = prepare_antibiotic_data(esbl_matrix, master_meta)

    # Generate all plots
    corr_df, pval_df = plot_correlation_heatmap(merged, beta_lactams)
    plot_ctxm_cephalosporin_scatter(merged)
    plot_exposed_vs_unexposed_boxplots(merged)
    forest_results = plot_forest_plot(merged, beta_lactams)
    plot_combined_figure(merged)

    # Save correlation summary
    corr_df.to_csv(RESULTS_DIR / 'antibiotic_esbl_correlations.tsv', sep='\t')
    print(f"\n✓ Saved: {RESULTS_DIR / 'antibiotic_esbl_correlations.tsv'}")

    print("\n" + "="*60)
    print("VISUALIZATION COMPLETE")
    print("="*60)
    print(f"\nFigures saved to: {FIGURES_DIR}")

if __name__ == "__main__":
    main()
