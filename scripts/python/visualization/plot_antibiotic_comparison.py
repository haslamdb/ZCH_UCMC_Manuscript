#!/usr/bin/env python3
"""
Compare antibiotic usage between UCMC and ZCH.
Creates boxplots showing days of each antibiotic per patient.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

def load_and_combine_data():
    """Load both antibiotic files and combine with location indicator."""

    print("Loading antibiotic data...")

    # Load UC data
    uc_df = pd.read_csv('UC_antibiotics_parsed.csv')
    uc_df['Location'] = 'UCMC'

    # Load ZCH data
    zch_df = pd.read_csv('ZCH_antibiotics_parsed.csv')
    zch_df['Location'] = 'ZCH'

    # Get antibiotic columns (exclude record_id, mrn, and Location)
    uc_abx_cols = [col for col in uc_df.columns if any(abx in col for abx in
                   ['Ampicillin', 'Penicillin', 'Nafcillin', 'Cefotaxime', 'Moxalactam',
                    'Ceftazidime', 'Ceftriaxone', 'Cefepime', 'Gentamicin', 'Tobramycin',
                    'Vancomycin', 'Fluconazole', 'Amphotericin', 'Meropenem', 'Azithromycin',
                    'Piperacillin'])]

    # Standardize column names between datasets
    # UC has mrn_v2, ZCH has mrn_v2_zju
    uc_df = uc_df.rename(columns={'mrn_v2': 'MRN'})
    zch_df = zch_df.rename(columns={'mrn_v2_zju': 'MRN'})

    # Select relevant columns for combining
    uc_cols = ['MRN', 'Location'] + uc_abx_cols
    zch_cols = ['MRN', 'Location'] + uc_abx_cols  # Use same column list

    # Combine datasets
    combined_df = pd.concat([uc_df[uc_cols], zch_df[zch_cols]], ignore_index=True)

    print(f"Combined data: {len(uc_df)} UCMC patients + {len(zch_df)} ZCH patients = {len(combined_df)} total")

    return combined_df, uc_abx_cols

def calculate_statistics(data, antibiotic_cols):
    """Calculate statistics comparing UCMC vs ZCH for each antibiotic."""

    stats_results = []

    for col in antibiotic_cols:
        ucmc_data = data[data['Location'] == 'UCMC'][col]
        zch_data = data[data['Location'] == 'ZCH'][col]

        # Remove zeros for statistics (only count patients who received the antibiotic)
        ucmc_nonzero = ucmc_data[ucmc_data > 0]
        zch_nonzero = zch_data[zch_data > 0]

        # Calculate prevalence (% of patients who received this antibiotic)
        ucmc_prevalence = (ucmc_data > 0).sum() / len(ucmc_data) * 100
        zch_prevalence = (zch_data > 0).sum() / len(zch_data) * 100

        # Wilcoxon rank-sum test (Mann-Whitney U) - comparing all values including zeros
        if len(ucmc_data) > 0 and len(zch_data) > 0:
            stat, pval = stats.mannwhitneyu(ucmc_data, zch_data, alternative='two-sided')
        else:
            pval = np.nan

        # Mean days among those who received the antibiotic
        ucmc_mean = ucmc_nonzero.mean() if len(ucmc_nonzero) > 0 else 0
        zch_mean = zch_nonzero.mean() if len(zch_nonzero) > 0 else 0

        stats_results.append({
            'Antibiotic': col,
            'UCMC_N': len(ucmc_nonzero),
            'UCMC_Prevalence_%': ucmc_prevalence,
            'UCMC_Mean_Days': ucmc_mean,
            'ZCH_N': len(zch_nonzero),
            'ZCH_Prevalence_%': zch_prevalence,
            'ZCH_Mean_Days': zch_mean,
            'P_value': pval,
            'Significant': 'Yes' if pval < 0.05 else 'No'
        })

    stats_df = pd.DataFrame(stats_results)

    # Sort by p-value
    stats_df = stats_df.sort_values('P_value')

    return stats_df

def create_boxplot_grid(data, antibiotic_cols, stats_df, output_dir='revision_figures'):
    """Create comprehensive boxplot grid for all antibiotics."""

    Path(output_dir).mkdir(exist_ok=True)

    # Filter to antibiotics that have at least some usage
    antibiotics_with_usage = []
    for col in antibiotic_cols:
        if data[col].sum() > 0:  # At least some patients received this
            antibiotics_with_usage.append(col)

    print(f"\nAntibiotics with usage: {len(antibiotics_with_usage)}/{len(antibiotic_cols)}")

    if len(antibiotics_with_usage) == 0:
        print("No antibiotics with usage found!")
        return

    # Calculate grid dimensions
    n_antibiotics = len(antibiotics_with_usage)
    ncols = 4
    nrows = int(np.ceil(n_antibiotics / ncols))

    # Create figure
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 4*nrows))
    axes = axes.flatten() if n_antibiotics > 1 else [axes]

    # Color palette
    colors = {'UCMC': '#E69F00', 'ZCH': '#56B4E9'}

    for idx, antibiotic in enumerate(antibiotics_with_usage):
        ax = axes[idx]

        # Prepare data for plotting
        plot_data = data[['Location', antibiotic]].copy()
        plot_data = plot_data[plot_data[antibiotic] > 0]  # Only show patients who received it

        if len(plot_data) == 0:
            ax.text(0.5, 0.5, 'No usage', ha='center', va='center')
            ax.set_title(antibiotic.replace('_', ' '))
            ax.axis('off')
            continue

        # Create boxplot
        positions = [0, 1]
        box_data = [plot_data[plot_data['Location'] == 'UCMC'][antibiotic],
                    plot_data[plot_data['Location'] == 'ZCH'][antibiotic]]

        bp = ax.boxplot(box_data, positions=positions, widths=0.6,
                        patch_artist=True, showfliers=False)

        # Color boxes
        for patch, location in zip(bp['boxes'], ['UCMC', 'ZCH']):
            patch.set_facecolor(colors[location])
            patch.set_alpha(0.7)

        # Add individual points with jitter
        for pos, location in enumerate(['UCMC', 'ZCH']):
            location_data = plot_data[plot_data['Location'] == location][antibiotic]
            if len(location_data) > 0:
                # Add jitter to x-axis
                x_jitter = np.random.normal(pos, 0.04, size=len(location_data))
                ax.scatter(x_jitter, location_data, alpha=0.6, s=30,
                          color=colors[location], edgecolors='black', linewidth=0.5)

        # Get p-value for this antibiotic
        pval = stats_df[stats_df['Antibiotic'] == antibiotic]['P_value'].values[0]

        # Format p-value
        if pd.notna(pval):
            if pval < 0.001:
                pval_str = 'p<0.001***'
            elif pval < 0.01:
                pval_str = f'p={pval:.3f}**'
            elif pval < 0.05:
                pval_str = f'p={pval:.3f}*'
            else:
                pval_str = f'p={pval:.3f}'
        else:
            pval_str = 'p=N/A'

        # Add p-value as text
        y_max = plot_data[antibiotic].max()
        ax.text(0.5, y_max * 1.15, pval_str, ha='center', va='bottom',
               fontsize=10, fontweight='bold')

        # Styling
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['UCMC', 'ZCH'], fontsize=11)
        ax.set_ylabel('Days', fontsize=11)
        ax.set_title(antibiotic.replace('_', ' '), fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3)

    # Hide unused subplots
    for idx in range(len(antibiotics_with_usage), len(axes)):
        axes[idx].axis('off')

    plt.tight_layout()

    # Save figure
    output_file = f'{output_dir}/Antibiotic_Comparison_UCMC_vs_ZCH.pdf'
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"\n✓ Saved boxplot figure: {output_file}")

    plt.close()

def create_summary_figure(stats_df, output_dir='revision_figures'):
    """Create summary figure showing significant differences."""

    Path(output_dir).mkdir(exist_ok=True)

    # Filter to antibiotics with p < 0.05
    sig_antibiotics = stats_df[stats_df['P_value'] < 0.05].copy()

    if len(sig_antibiotics) == 0:
        print("\nNo significant differences found between UCMC and ZCH")
        return

    print(f"\n{len(sig_antibiotics)} antibiotics with significant differences (p<0.05)")

    # Create bar plot
    fig, ax = plt.subplots(figsize=(10, max(6, len(sig_antibiotics)*0.4)))

    x = np.arange(len(sig_antibiotics))
    width = 0.35

    ucmc_prev = sig_antibiotics['UCMC_Prevalence_%'].values
    zch_prev = sig_antibiotics['ZCH_Prevalence_%'].values

    ax.barh(x - width/2, ucmc_prev, width, label='UCMC', color='#E69F00', alpha=0.8)
    ax.barh(x + width/2, zch_prev, width, label='ZCH', color='#56B4E9', alpha=0.8)

    ax.set_yticks(x)
    ax.set_yticklabels([abx.replace('_', ' ') for abx in sig_antibiotics['Antibiotic']],
                       fontsize=10)
    ax.set_xlabel('Prevalence (%)', fontsize=12)
    ax.set_title('Antibiotics with Significant Differences in Usage\n(p < 0.05)',
                fontsize=14, fontweight='bold')
    ax.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()

    output_file = f'{output_dir}/Significant_Antibiotic_Differences.pdf'
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"✓ Saved significance summary: {output_file}")

    plt.close()

def main():
    """Main function."""

    print("="*70)
    print("ANTIBIOTIC USAGE COMPARISON: UCMC vs ZCH")
    print("="*70)

    # Load and combine data
    combined_df, antibiotic_cols = load_and_combine_data()

    # Calculate statistics
    print("\nCalculating statistics...")
    stats_df = calculate_statistics(combined_df, antibiotic_cols)

    # Save statistics
    stats_file = 'revision_figures/Antibiotic_Statistics_UCMC_vs_ZCH.csv'
    stats_df.to_csv(stats_file, index=False)
    print(f"✓ Saved statistics: {stats_file}")

    # Print summary
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    print(f"\nAntibiotics with significant differences (p<0.05): {(stats_df['P_value'] < 0.05).sum()}")

    if (stats_df['P_value'] < 0.05).any():
        print("\nSignificant antibiotics:")
        sig_df = stats_df[stats_df['P_value'] < 0.05][['Antibiotic', 'UCMC_Prevalence_%',
                                                         'ZCH_Prevalence_%', 'P_value']]
        print(sig_df.to_string(index=False))

    # Create visualizations
    print("\n" + "="*70)
    print("Creating visualizations...")
    print("="*70)

    create_boxplot_grid(combined_df, antibiotic_cols, stats_df)
    create_summary_figure(stats_df)

    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)
    print("\nGenerated files:")
    print("  - revision_figures/Antibiotic_Comparison_UCMC_vs_ZCH.pdf")
    print("  - revision_figures/Antibiotic_Statistics_UCMC_vs_ZCH.csv")
    print("  - revision_figures/Significant_Antibiotic_Differences.pdf")

if __name__ == '__main__':
    main()
