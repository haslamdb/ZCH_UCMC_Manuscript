#!/usr/bin/env python3
"""
Analyze mecA gene prevalence and abundance differences between ZCH (Hangzhou) and UCMC (Cincinnati)

This script:
1. Loads AMR gene abundance data and sample metadata
2. Filters for mecA gene
3. Compares prevalence (presence/absence) and abundance between sites
4. Performs statistical tests (chi-square for prevalence, Mann-Whitney U for abundance)
5. Controls for potential confounders (sample type, gestational age, antibiotics)
6. Generates visualizations and summary statistics
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


def load_data(project_dir: Path, use_diamond: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load AMR data and metadata.

    Args:
        project_dir: Project directory path
        use_diamond: If True, use DIAMOND traditional results (more sensitive)
                    If False, use the stricter nicu_amr_only.tsv
    """
    if use_diamond:
        # Use DIAMOND traditional results (more comprehensive)
        amr_path = project_dir / "nicu_resistome_analysis/diamond_traditional/data/diamond_amr_combined.tsv"
        print(f"Loading DIAMOND AMR data from {amr_path}...")
        amr_df = pd.read_csv(amr_path, sep='\t')
        # Standardize column names if needed
        if 'rpm' not in amr_df.columns and 'RPM' in amr_df.columns:
            amr_df = amr_df.rename(columns={'RPM': 'rpm'})
        print(f"  Loaded {len(amr_df)} AMR gene observations")

        # Load metadata from DIAMOND directory (has matching sample names)
        meta_path = project_dir / "nicu_resistome_analysis/diamond_traditional/data/master_metadata_fixed.tsv"
        print(f"Loading metadata from {meta_path}...")
        meta_df = pd.read_csv(meta_path, sep='\t')
        print(f"  Loaded {len(meta_df)} samples with metadata")
    else:
        # Use original strict filtering
        amr_path = project_dir / "nicu_resistome_analysis/data/nicu_amr_only.tsv"
        print(f"Loading AMR data from {amr_path}...")
        amr_df = pd.read_csv(amr_path, sep='\t')
        print(f"  Loaded {len(amr_df)} AMR gene observations")

        meta_path = project_dir / "metadata/AllNICUSampleKeyRevised20250320_for_HanzhouCincinnatiSamples.csv"
        print(f"Loading metadata from {meta_path}...")
        meta_df = pd.read_csv(meta_path)
        print(f"  Loaded {len(meta_df)} samples with metadata")

    return amr_df, meta_df


def filter_mecA(amr_df: pd.DataFrame) -> pd.DataFrame:
    """Filter AMR data for mecA gene."""
    mecA_df = amr_df[amr_df['gene_name'].str.lower() == 'meca'].copy()
    print(f"Found {len(mecA_df)} mecA observations across {mecA_df['sample_name'].nunique()} samples")
    return mecA_df


def merge_with_metadata(mecA_df: pd.DataFrame, meta_df: pd.DataFrame) -> pd.DataFrame:
    """Merge mecA data with metadata and create presence/absence for all samples."""
    # Determine which columns are available (structure differs between data sources)
    if 'sample_name' in meta_df.columns:
        # DIAMOND metadata format
        sample_col = 'sample_name'
        patient_col = 'PatientID'
    else:
        sample_col = 'SampleID'
        patient_col = 'PatientID'

    # Get available columns
    desired_cols = [sample_col, patient_col, 'SampleType', 'Location',
                    'GestationalAge', 'GestationCohort', 'SampleCollectionWeek',
                    'PostNatalAntibiotics', 'MaternalAntibiotics']
    available_cols = [c for c in desired_cols if c in meta_df.columns]

    all_samples = meta_df[available_cols].copy()
    if sample_col != 'sample_name':
        all_samples = all_samples.rename(columns={sample_col: 'sample_name'})

    # Standardize location names
    if 'Location' in all_samples.columns:
        all_samples['Location'] = all_samples['Location'].replace({
            'Cincinnati': 'UCMC',
            'Hangzhou': 'ZCH'
        })

    # Aggregate mecA by sample (sum if multiple alleles)
    # Handle both DIAMOND format (has rpm but not tpkm) and other formats
    agg_dict = {'read_count': 'sum'}
    if 'percent_coverage' in mecA_df.columns:
        agg_dict['percent_coverage'] = 'max'
    if 'mean_depth' in mecA_df.columns:
        agg_dict['mean_depth'] = 'max'
    if 'tpkm' in mecA_df.columns:
        agg_dict['tpkm'] = 'sum'
    if 'rpm' in mecA_df.columns:
        agg_dict['rpm'] = 'sum'

    mecA_by_sample = mecA_df.groupby('sample_name').agg(agg_dict).reset_index()

    # Merge - left join to keep all samples
    merged = all_samples.merge(mecA_by_sample, on='sample_name', how='left')

    # Fill NaN with 0 for samples without mecA
    abundance_cols = [c for c in ['read_count', 'percent_coverage', 'mean_depth', 'tpkm', 'rpm']
                      if c in merged.columns]
    merged[abundance_cols] = merged[abundance_cols].fillna(0)

    # Create presence/absence column
    merged['mecA_present'] = (merged['read_count'] > 0).astype(int)

    # Create log-transformed abundance (add pseudocount for zeros)
    if 'tpkm' in merged.columns:
        merged['log_tpkm'] = np.log10(merged['tpkm'] + 1)
    if 'rpm' in merged.columns:
        merged['log_rpm'] = np.log10(merged['rpm'] + 1)

    print(f"Merged dataset: {len(merged)} samples")
    print(f"  mecA present in {merged['mecA_present'].sum()} samples ({100*merged['mecA_present'].mean():.1f}%)")

    return merged


def analyze_prevalence(df: pd.DataFrame) -> dict:
    """Analyze mecA prevalence differences between sites."""
    print("\n" + "="*80)
    print("PREVALENCE ANALYSIS (Presence/Absence)")
    print("="*80)

    results = {}

    # Overall prevalence by location
    prevalence = df.groupby('Location').agg({
        'mecA_present': ['sum', 'count', 'mean']
    })
    prevalence.columns = ['n_positive', 'n_total', 'prevalence']
    prevalence['prevalence_pct'] = prevalence['prevalence'] * 100

    print("\nOverall mecA Prevalence by Location:")
    print(prevalence.to_string())
    results['prevalence_by_location'] = prevalence

    # Chi-square test for independence
    contingency = pd.crosstab(df['Location'], df['mecA_present'])
    chi2, pvalue, dof, expected = stats.chi2_contingency(contingency)

    print(f"\nChi-square test: χ² = {chi2:.4f}, p = {pvalue:.4e}, df = {dof}")
    results['chi2_test'] = {'chi2': chi2, 'pvalue': pvalue, 'dof': dof}

    # Fisher's exact test (more appropriate for small counts)
    if contingency.shape == (2, 2):
        odds_ratio, fisher_p = stats.fisher_exact(contingency)
        print(f"Fisher's exact test: OR = {odds_ratio:.4f}, p = {fisher_p:.4e}")
        results['fisher_test'] = {'odds_ratio': odds_ratio, 'pvalue': fisher_p}

    return results


def analyze_abundance(df: pd.DataFrame) -> dict:
    """Analyze mecA abundance differences between sites (among positive samples)."""
    print("\n" + "="*80)
    print("ABUNDANCE ANALYSIS (Among mecA-positive samples)")
    print("="*80)

    results = {}

    # Filter to mecA-positive samples only
    positive_df = df[df['mecA_present'] == 1].copy()

    if len(positive_df) < 2:
        print("Not enough mecA-positive samples for abundance analysis")
        return results

    # Determine which abundance metrics are available
    available_metrics = [m for m in ['rpm', 'tpkm', 'read_count'] if m in positive_df.columns]

    # Summary statistics by location
    agg_funcs = ['count', 'mean', 'std', 'median']
    agg_dict = {m: agg_funcs for m in available_metrics}

    abundance_stats = positive_df.groupby('Location').agg(agg_dict)

    print("\nAbundance statistics (mecA-positive samples only):")
    print(abundance_stats.to_string())
    results['abundance_stats'] = abundance_stats

    # Mann-Whitney U test for each abundance metric
    for metric in available_metrics:
        zch = positive_df[positive_df['Location'] == 'ZCH'][metric]
        ucmc = positive_df[positive_df['Location'] == 'UCMC'][metric]

        if len(zch) > 0 and len(ucmc) > 0:
            stat, pvalue = stats.mannwhitneyu(zch, ucmc, alternative='two-sided')
            print(f"\nMann-Whitney U test for {metric}:")
            print(f"  ZCH (n={len(zch)}): median = {zch.median():.4f}")
            print(f"  UCMC (n={len(ucmc)}): median = {ucmc.median():.4f}")
            print(f"  U = {stat:.1f}, p = {pvalue:.4e}")
            results[f'mannwhitney_{metric}'] = {'statistic': stat, 'pvalue': pvalue}

    return results


def analyze_by_sample_type(df: pd.DataFrame) -> dict:
    """Stratified analysis by sample type (Stool, Axilla, Groin)."""
    print("\n" + "="*80)
    print("STRATIFIED ANALYSIS BY SAMPLE TYPE")
    print("="*80)

    results = {}

    for sample_type in df['SampleType'].unique():
        if pd.isna(sample_type):
            continue

        subset = df[df['SampleType'] == sample_type]
        print(f"\n--- {sample_type} samples ---")

        prevalence = subset.groupby('Location')['mecA_present'].agg(['sum', 'count', 'mean'])
        prevalence.columns = ['n_positive', 'n_total', 'prevalence']
        print(prevalence.to_string())

        # Chi-square test
        contingency = pd.crosstab(subset['Location'], subset['mecA_present'])
        if contingency.shape == (2, 2) and contingency.values.min() >= 0:
            try:
                chi2, pvalue, _, _ = stats.chi2_contingency(contingency)
                odds_ratio, fisher_p = stats.fisher_exact(contingency)
                print(f"  Chi-square: p = {pvalue:.4e}")
                print(f"  Fisher's exact: OR = {odds_ratio:.4f}, p = {fisher_p:.4e}")
                results[sample_type] = {
                    'chi2_pvalue': pvalue,
                    'fisher_or': odds_ratio,
                    'fisher_pvalue': fisher_p
                }
            except Exception as e:
                print(f"  Could not compute test: {e}")

    return results


def analyze_confounders(df: pd.DataFrame) -> dict:
    """Check for confounding factors between sites."""
    print("\n" + "="*80)
    print("CONFOUNDER ANALYSIS")
    print("="*80)

    results = {}

    # Sample type distribution
    print("\nSample type distribution by location:")
    sample_type_dist = pd.crosstab(df['Location'], df['SampleType'], normalize='index') * 100
    print(sample_type_dist.to_string())

    # Gestational age distribution
    print("\nGestational age cohort distribution by location:")
    gest_dist = pd.crosstab(df['Location'], df['GestationCohort'], normalize='index') * 100
    print(gest_dist.to_string())

    # Antibiotic exposure
    print("\nPostnatal antibiotic exposure by location:")
    abx_dist = pd.crosstab(df['Location'], df['PostNatalAntibiotics'], normalize='index') * 100
    print(abx_dist.to_string())

    # Collection week
    print("\nSample collection week by location:")
    week_dist = pd.crosstab(df['Location'], df['SampleCollectionWeek'], normalize='index') * 100
    print(week_dist.to_string())

    return results


def create_visualizations(df: pd.DataFrame, output_dir: Path):
    """Create publication-quality visualizations."""
    print("\n" + "="*80)
    print("GENERATING VISUALIZATIONS")
    print("="*80)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Set style
    try:
        plt.style.use('seaborn-v0_8-whitegrid')
    except:
        plt.style.use('seaborn-whitegrid')

    # Color mapping for locations
    color_map = {'UCMC': '#2ecc71', 'ZCH': '#e74c3c'}

    # 1. Prevalence bar plot
    fig, ax = plt.subplots(figsize=(8, 6))
    prevalence = df.groupby('Location')['mecA_present'].mean() * 100
    colors = [color_map.get(loc, '#3498db') for loc in prevalence.index]
    bars = prevalence.plot(kind='bar', ax=ax, color=colors, edgecolor='black', linewidth=1.5)
    ax.set_ylabel('mecA Prevalence (%)', fontsize=12)
    ax.set_xlabel('Location', fontsize=12)
    ax.set_title('mecA Gene Prevalence by Hospital Location', fontsize=14, fontweight='bold')

    # Add sample counts
    for i, (loc, prev) in enumerate(prevalence.items()):
        n = df[df['Location'] == loc].shape[0]
        n_pos = df[(df['Location'] == loc) & (df['mecA_present'] == 1)].shape[0]
        ax.text(i, prev + 1, f'n={n_pos}/{n}', ha='center', fontsize=10)

    plt.tight_layout()
    fig.savefig(output_dir / 'mecA_prevalence_by_location.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(output_dir / 'mecA_prevalence_by_location.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved prevalence plot")

    # 2. Prevalence by sample type (grouped bar)
    if 'SampleType' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        prevalence_by_type = df.groupby(['SampleType', 'Location'])['mecA_present'].mean().unstack() * 100
        prevalence_by_type.plot(kind='bar', ax=ax, color=[color_map.get(c, '#3498db') for c in prevalence_by_type.columns],
                                edgecolor='black', linewidth=1)
        ax.set_ylabel('mecA Prevalence (%)', fontsize=12)
        ax.set_xlabel('Sample Type', fontsize=12)
        ax.set_title('mecA Prevalence by Sample Type and Location', fontsize=14, fontweight='bold')
        ax.legend(title='Location')
        plt.xticks(rotation=45)
        plt.tight_layout()
        fig.savefig(output_dir / 'mecA_prevalence_by_sample_type.pdf', dpi=300, bbox_inches='tight')
        fig.savefig(output_dir / 'mecA_prevalence_by_sample_type.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved sample type plot")

    # 3. Abundance distribution (violin plot, mecA-positive only)
    positive_df = df[df['mecA_present'] == 1].copy()
    abundance_col = 'log_rpm' if 'log_rpm' in positive_df.columns else 'log_tpkm'

    if len(positive_df) > 5 and abundance_col in positive_df.columns:
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.violinplot(data=positive_df, x='Location', y=abundance_col,
                       palette=color_map, ax=ax)
        sns.stripplot(data=positive_df, x='Location', y=abundance_col,
                      color='black', alpha=0.5, size=4, ax=ax)
        metric_name = 'RPM' if 'rpm' in abundance_col else 'TPKM'
        ax.set_ylabel(f'log10({metric_name} + 1)', fontsize=12)
        ax.set_xlabel('Location', fontsize=12)
        ax.set_title('mecA Abundance Distribution\n(mecA-positive samples only)',
                     fontsize=14, fontweight='bold')
        plt.tight_layout()
        fig.savefig(output_dir / 'mecA_abundance_distribution.pdf', dpi=300, bbox_inches='tight')
        fig.savefig(output_dir / 'mecA_abundance_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved abundance distribution plot")

    # 4. Prevalence by antibiotic exposure
    if 'PostNatalAntibiotics' in df.columns:
        locations = df['Location'].unique()
        fig, axes = plt.subplots(1, len(locations), figsize=(6*len(locations), 5))
        if len(locations) == 1:
            axes = [axes]

        for i, location in enumerate(locations):
            subset = df[df['Location'] == location]
            prev_by_abx = subset.groupby('PostNatalAntibiotics')['mecA_present'].agg(['mean', 'count'])
            prev_by_abx['mean'] = prev_by_abx['mean'] * 100

            bars = prev_by_abx['mean'].plot(kind='bar', ax=axes[i],
                                             color=color_map.get(location, '#3498db'),
                                             edgecolor='black')
            axes[i].set_title(f'{location}', fontsize=12, fontweight='bold')
            axes[i].set_ylabel('mecA Prevalence (%)')
            axes[i].set_xlabel('Postnatal Antibiotics')
            axes[i].tick_params(axis='x', rotation=45)

            # Add counts
            for j, (idx, row) in enumerate(prev_by_abx.iterrows()):
                axes[i].text(j, row['mean'] + 1, f'n={int(row["count"])}', ha='center', fontsize=8)

        plt.suptitle('mecA Prevalence by Antibiotic Exposure', fontsize=14, fontweight='bold')
        plt.tight_layout()
        fig.savefig(output_dir / 'mecA_prevalence_by_antibiotics.pdf', dpi=300, bbox_inches='tight')
        fig.savefig(output_dir / 'mecA_prevalence_by_antibiotics.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved antibiotic exposure plot")


def save_results(df: pd.DataFrame, results: dict, output_dir: Path):
    """Save analysis results to files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save merged data
    df.to_csv(output_dir / 'mecA_merged_data.csv', index=False)
    print(f"\nSaved merged data to {output_dir / 'mecA_merged_data.csv'}")

    # Save summary statistics
    summary = []

    # Determine abundance column
    abundance_col = 'rpm' if 'rpm' in df.columns else 'tpkm'

    # Overall stats
    for loc in df['Location'].unique():
        subset = df[df['Location'] == loc]
        positive_subset = subset[subset['mecA_present'] == 1]
        summary.append({
            'Location': loc,
            'n_samples': len(subset),
            'n_mecA_positive': subset['mecA_present'].sum(),
            'prevalence_pct': subset['mecA_present'].mean() * 100,
            f'mean_{abundance_col}_positive': positive_subset[abundance_col].mean() if len(positive_subset) > 0 else np.nan,
            f'median_{abundance_col}_positive': positive_subset[abundance_col].median() if len(positive_subset) > 0 else np.nan
        })

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_dir / 'mecA_summary_statistics.csv', index=False)
    print(f"Saved summary statistics to {output_dir / 'mecA_summary_statistics.csv'}")


def print_final_summary(df: pd.DataFrame, results: dict):
    """Print a final summary of key findings."""
    print("\n" + "="*80)
    print("SUMMARY OF KEY FINDINGS")
    print("="*80)

    for loc in df['Location'].unique():
        subset = df[df['Location'] == loc]
        prev = subset['mecA_present'].mean() * 100
        n_pos = subset['mecA_present'].sum()
        n_total = len(subset)
        print(f"\n{loc}:")
        print(f"  Total samples: {n_total}")
        print(f"  mecA positive: {n_pos} ({prev:.1f}%)")

    if 'chi2_test' in results:
        print(f"\nStatistical comparison:")
        print(f"  Chi-square p-value: {results['chi2_test']['pvalue']:.4e}")

    if 'fisher_test' in results:
        print(f"  Fisher's exact OR: {results['fisher_test']['odds_ratio']:.2f}")
        print(f"  Fisher's exact p-value: {results['fisher_test']['pvalue']:.4e}")

    # Interpretation
    if 'fisher_test' in results:
        or_val = results['fisher_test']['odds_ratio']
        p_val = results['fisher_test']['pvalue']
        locations = df['Location'].unique()

        if p_val < 0.05:
            if or_val > 1:
                print(f"\n=> mecA is significantly MORE prevalent in {locations[1]} compared to {locations[0]}")
            else:
                print(f"\n=> mecA is significantly MORE prevalent in {locations[0]} compared to {locations[1]}")
        else:
            print(f"\n=> No significant difference in mecA prevalence between sites")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze mecA prevalence differences between ZCH and UCMC"
    )
    parser.add_argument(
        "--project-dir",
        type=Path,
        default=Path("/home/david/projects/ZCH_UCMC_Manuscript"),
        help="Project directory path"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: project_dir/results/mecA_analysis)"
    )
    args = parser.parse_args()

    # Set output directory
    output_dir = args.output_dir or (args.project_dir / "results/mecA_analysis")

    print("="*80)
    print("mecA Gene Prevalence Analysis: ZCH (Hangzhou) vs UCMC (Cincinnati)")
    print("="*80)

    # Load data
    amr_df, meta_df = load_data(args.project_dir)

    # Filter for mecA
    mecA_df = filter_mecA(amr_df)

    # Merge with metadata
    merged_df = merge_with_metadata(mecA_df, meta_df)

    # Run analyses
    all_results = {}

    # Prevalence analysis
    prevalence_results = analyze_prevalence(merged_df)
    all_results.update(prevalence_results)

    # Abundance analysis
    abundance_results = analyze_abundance(merged_df)
    all_results.update(abundance_results)

    # Stratified by sample type
    stratified_results = analyze_by_sample_type(merged_df)
    all_results['stratified'] = stratified_results

    # Check confounders
    analyze_confounders(merged_df)

    # Create visualizations
    create_visualizations(merged_df, output_dir / "figures")

    # Save results
    save_results(merged_df, all_results, output_dir)

    # Print summary
    print_final_summary(merged_df, all_results)

    print("\n" + "="*80)
    print(f"Analysis complete! Results saved to {output_dir}")
    print("="*80)


if __name__ == "__main__":
    main()
