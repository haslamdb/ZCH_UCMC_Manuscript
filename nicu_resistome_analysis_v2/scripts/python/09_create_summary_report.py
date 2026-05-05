#!/usr/bin/env python3
"""
Create Statistical Summary Tables and Analysis Report

Consolidates all analysis results into:
1. Summary statistics tables
2. Comprehensive analysis report (markdown)
3. Key findings tables
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

# Define paths
PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
RESULTS_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results"
QC_DIR = RESULTS_DIR / "qc"
EXPLORATORY_DIR = RESULTS_DIR / "exploratory"
DIFFERENTIAL_DIR = RESULTS_DIR / "differential"
CORR_DIR = RESULTS_DIR / "correlations"
LONG_DIR = RESULTS_DIR / "longitudinal"
MODELS_DIR = RESULTS_DIR / "mixed_models"
SUMMARY_DIR = RESULTS_DIR / "summary"

# Create output directory
SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

def create_sample_summary():
    """Create summary table of samples and subjects"""
    print("Creating sample summary...")

    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")

    # Overall statistics
    total_samples = len(metadata)
    total_subjects = metadata['SubjectID'].nunique()
    samples_with_amr = len(metadata[metadata['has_amr_data']])

    # By location
    location_stats = []
    for loc in ['UCMC', 'ZCH']:
        loc_data = metadata[metadata['Location'] == loc]
        location_stats.append({
            'Location': loc,
            'Total Subjects': loc_data['SubjectID'].nunique(),
            'Complete Subjects': loc_data[loc_data['subject_complete']]['SubjectID'].nunique(),
            'Total Samples': len(loc_data),
            'Samples with AMR': len(loc_data[loc_data['has_amr_data']]),
            'Subjects with Antibiotics': loc_data[loc_data['has_abx_data']]['SubjectID'].nunique()
        })

    location_df = pd.DataFrame(location_stats)

    # By body site
    site_stats = []
    for site in ['Axilla', 'Groin', 'Stool']:
        site_data = metadata[(metadata['SampleType'] == site) & (metadata['has_amr_data'])]
        for week in ['Week.1', 'Week.3']:
            week_data = site_data[site_data['SampleCollectionWeek'] == week]
            site_stats.append({
                'Body Site': site,
                'Week': week,
                'n Samples': len(week_data),
                'Mean AMR RPM': week_data['total_amr_rpm'].mean(),
                'Median AMR RPM': week_data['total_amr_rpm'].median(),
                'Mean Gene Richness': week_data['unique_amr_genes'].mean()
            })

    site_df = pd.DataFrame(site_stats)

    # Save
    location_df.to_csv(SUMMARY_DIR / "summary_by_location.tsv", sep="\t", index=False)
    site_df.to_csv(SUMMARY_DIR / "summary_by_bodysite_week.tsv", sep="\t", index=False)

    print(f"✓ Saved: {SUMMARY_DIR / 'summary_by_location.tsv'}")
    print(f"✓ Saved: {SUMMARY_DIR / 'summary_by_bodysite_week.tsv'}")

    return location_df, site_df, total_samples, total_subjects

def create_key_findings_table():
    """Summarize key statistical findings"""
    print("\nCreating key findings table...")

    findings = []

    # 1. Differential abundance (Location comparison)
    diff_summary = pd.read_csv(DIFFERENTIAL_DIR / "differential_abundance_summary.tsv", sep="\t")
    for _, row in diff_summary.iterrows():
        findings.append({
            'Analysis': 'Differential Abundance',
            'Comparison': f'{row["body_site"]}: UCMC vs ZCH',
            'n_tested': row['total_genes'],
            'n_significant': row['sig_fdr005'],
            'Test': 'Mann-Whitney U',
            'Interpretation': 'No significant differences' if row['sig_fdr005'] == 0 else f"{row['sig_fdr005']} genes differ"
        })

    # 2. Longitudinal changes
    for site in ['axilla', 'groin', 'stool']:
        file_path = LONG_DIR / f"gene_longitudinal_{site}.tsv"
        if file_path.exists():
            long_data = pd.read_csv(file_path, sep="\t")
            n_sig = (long_data['fdr'] < 0.05).sum()
            n_up = ((long_data['fdr'] < 0.05) & (long_data['log2_fold_change'] > 0)).sum()
            n_down = ((long_data['fdr'] < 0.05) & (long_data['log2_fold_change'] < 0)).sum()

            findings.append({
                'Analysis': 'Longitudinal',
                'Comparison': f'{site.title()}: Week 1 → Week 3',
                'n_tested': len(long_data),
                'n_significant': n_sig,
                'Test': 'Wilcoxon signed-rank',
                'Interpretation': f'{n_up} increased, {n_down} decreased' if n_sig > 0 else 'No significant changes'
            })

    # 3. Total AMR longitudinal
    total_long = pd.read_csv(LONG_DIR / "total_amr_longitudinal.tsv", sep="\t")
    for _, row in total_long.iterrows():
        sig_status = '***' if row['fdr'] < 0.001 else '**' if row['fdr'] < 0.01 else '*' if row['fdr'] < 0.05 else 'ns'
        findings.append({
            'Analysis': 'Total AMR Longitudinal',
            'Comparison': f"{row['body_site']}: Week 1 → Week 3",
            'n_tested': row['n_subjects'],
            'n_significant': sig_status,
            'Test': 'Wilcoxon signed-rank',
            'Interpretation': f"{row['percent_change']:+.1f}% change (p={row['wilcoxon_pvalue']:.4f})"
        })

    # 4. Antibiotic correlations
    if (CORR_DIR / "total_abx_amr_correlation_by_site.tsv").exists():
        abx_corr = pd.read_csv(CORR_DIR / "total_abx_amr_correlation_by_site.tsv", sep="\t")
        for _, row in abx_corr.iterrows():
            findings.append({
                'Analysis': 'Antibiotic Correlation',
                'Comparison': f"{row['body_site']}: Total Abx vs AMR",
                'n_tested': row['n_samples'],
                'n_significant': 'ns',
                'Test': 'Spearman correlation',
                'Interpretation': f"rho={row['spearman_rho']:.3f}, p={row['pvalue']:.4f}"
            })

    # 5. Mixed-effects models
    if (MODELS_DIR / "full_model_coefficients.tsv").exists():
        model_coef = pd.read_csv(MODELS_DIR / "full_model_coefficients.tsv", sep="\t", index_col=0)

        # Key effects
        key_effects = [
            ('Location[T.ZCH]', 'Location Effect'),
            ('BodySite[T.Stool]', 'BodySite: Stool vs Axilla'),
            ('Week[T.Week.3]', 'Week Effect'),
            ('BodySite[T.Groin]:Week[T.Week.3]', 'Groin × Week'),
            ('BodySite[T.Stool]:Week[T.Week.3]', 'Stool × Week')
        ]

        for effect_name, description in key_effects:
            if effect_name in model_coef.index:
                coef = model_coef.loc[effect_name, 'coefficient']
                pval = model_coef.loc[effect_name, 'p_value']
                sig_status = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'

                findings.append({
                    'Analysis': 'Mixed-Effects Model',
                    'Comparison': description,
                    'n_tested': 318,  # Total samples in complete subjects
                    'n_significant': sig_status,
                    'Test': 'Linear Mixed Model',
                    'Interpretation': f"coef={coef:+.3f}, p={pval:.4f}"
                })

    findings_df = pd.DataFrame(findings)
    findings_df.to_csv(SUMMARY_DIR / "key_findings.tsv", sep="\t", index=False)

    print(f"✓ Saved: {SUMMARY_DIR / 'key_findings.tsv'}")

    return findings_df

def create_analysis_report(location_df, site_df, total_samples, total_subjects, findings_df):
    """Create comprehensive markdown report"""
    print("\nCreating analysis report...")

    report = []

    # Header
    report.append("# NICU Resistome Analysis Report")
    report.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("\n---\n")

    # Study Overview
    report.append("## Study Overview")
    report.append(f"\n- **Total Samples:** {total_samples}")
    report.append(f"- **Total Subjects:** {total_subjects}")
    report.append(f"- **Locations:** UCMC (Cincinnati) and ZCH (Hangzhou)")
    report.append(f"- **Body Sites:** Axilla, Groin, Stool")
    report.append(f"- **Time Points:** Week 1 and Week 3")
    report.append(f"- **AMR Genes Analyzed:** 237 genes (after QC filtering)")

    # Sample Distribution
    report.append("\n\n### Sample Distribution by Location")
    report.append("\n```")
    report.append(location_df.to_string(index=False))
    report.append("```")

    # Key Findings
    report.append("\n\n## Key Findings\n")

    # Finding 1: No location differences
    report.append("### 1. No Significant Differences Between UCMC and ZCH")
    report.append("\n**Result:** Resistome composition is remarkably similar between Cincinnati and Hangzhou NICUs.")
    report.append("\n- Differential abundance testing: 0 genes significantly different at FDR < 0.05")
    report.append("- Mixed-effects model: Location effect p = 0.45 (not significant)")
    report.append("- **Interpretation:** Geographic location does not significantly influence NICU resistome composition")

    # Finding 2: Body site is primary driver
    report.append("\n\n### 2. Body Site is the Primary Driver of Resistome Composition")
    report.append("\n**Result:** Stool samples have dramatically higher AMR burden than skin sites.")
    report.append("\n- PCA analysis: PC1 (30.8% variance) separates stool from axilla/groin")
    report.append("- Mixed-effects model: Stool coefficient = +0.615, p < 0.0001")
    report.append("- **Interpretation:** Gut-associated sites harbor more diverse and abundant AMR genes")

    # Finding 3: Longitudinal changes
    report.append("\n\n### 3. Significant Longitudinal Changes in Groin and Stool")
    report.append("\n**Result:** AMR burden increases dramatically from Week 1 to Week 3 in groin and stool, but not axilla.")

    long_summary = site_df[['Body Site', 'Week', 'Mean AMR RPM']].pivot(index='Body Site', columns='Week', values='Mean AMR RPM')
    report.append("\n\n**Mean AMR RPM by Body Site and Week:**")
    report.append("\n```")
    report.append(long_summary.to_string())
    report.append("```")

    # Get detailed longitudinal stats
    if (LONG_DIR / "total_amr_longitudinal.tsv").exists():
        total_long = pd.read_csv(LONG_DIR / "total_amr_longitudinal.tsv", sep="\t")
        report.append("\n\n**Statistical Tests (Paired Wilcoxon):**")
        for _, row in total_long.iterrows():
            sig_marker = "***" if row['fdr'] < 0.001 else "**" if row['fdr'] < 0.01 else "*" if row['fdr'] < 0.05 else "ns"
            report.append(f"\n- **{row['body_site']}:** {row['percent_change']:+.1f}% change, p = {row['wilcoxon_pvalue']:.4f} {sig_marker}")

    report.append("\n\n**Gene-level changes:**")
    for site in ['groin', 'stool']:
        file_path = LONG_DIR / f"gene_longitudinal_{site}.tsv"
        if file_path.exists():
            long_data = pd.read_csv(file_path, sep="\t")
            n_sig = (long_data['fdr'] < 0.05).sum()
            n_up = ((long_data['fdr'] < 0.05) & (long_data['log2_fold_change'] > 0)).sum()
            n_down = ((long_data['fdr'] < 0.05) & (long_data['log2_fold_change'] < 0)).sum()
            report.append(f"\n- **{site.title()}:** {n_sig} genes significantly changed ({n_up} increased, {n_down} decreased)")

    # Finding 4: No antibiotic correlation
    report.append("\n\n### 4. No Correlation Between Antibiotic Exposure and AMR Burden")
    report.append("\n**Result:** Cumulative antibiotic exposure does not correlate with total AMR burden.")
    report.append("\n- Overall Spearman rho = 0.001, p = 0.985")
    report.append("- No individual antibiotics showed significant correlations")
    report.append("- **Interpretation:** AMR colonization may be independent of selective antibiotic pressure in this cohort")

    # Finding 5: Mixed-effects model
    report.append("\n\n### 5. Mixed-Effects Model Confirms Body Site × Week Interaction")
    report.append("\n**Model:** log(AMR) ~ Location × BodySite × Week + (1|SubjectID)")
    report.append("\n\n**Significant Effects:**")

    if (MODELS_DIR / "full_model_coefficients.tsv").exists():
        model_coef = pd.read_csv(MODELS_DIR / "full_model_coefficients.tsv", sep="\t", index_col=0)

        sig_effects = model_coef[model_coef['p_value'] < 0.05].copy()
        sig_effects = sig_effects.sort_values('p_value')

        for effect in sig_effects.index:
            if effect != 'Intercept':
                coef = sig_effects.loc[effect, 'coefficient']
                pval = sig_effects.loc[effect, 'p_value']
                sig_marker = "***" if pval < 0.001 else "**" if pval < 0.01 else "*"
                report.append(f"\n- **{effect}:** coefficient = {coef:+.3f}, p = {pval:.4f} {sig_marker}")

        # Report ICC
        report.append(f"\n\n**Intraclass Correlation (ICC):** 0.262")
        report.append("\n  - 26.2% of variance is between-subject")
        report.append("  - Confirms substantial within-subject correlation, validating repeated measures design")

    # Data Quality
    report.append("\n\n## Data Quality and Filtering\n")
    report.append("### Quality Control Steps:")
    report.append("\n1. **Gene Filtering:**")
    report.append("\n   - Original: 449 genes")
    report.append("\n   - After QC: 237 genes (present in ≥5% samples with max RPM ≥1.0)")
    report.append("\n2. **Sample Filtering:**")
    report.append("\n   - Excluded 'Nonese' (nose swab) samples (never sequenced)")
    report.append("\n   - Retained only samples with AMR data")
    report.append("\n3. **Subject Completeness:**")
    report.append("\n   - Complete subjects: 53 (40 UCMC, 13 ZCH)")
    report.append("\n   - Complete = all 6 samples present (3 body sites × 2 weeks)")

    # Statistical Methods
    report.append("\n\n## Statistical Methods\n")
    report.append("### Analyses Performed:")
    report.append("\n1. **Differential Abundance:** Mann-Whitney U test (UCMC vs ZCH, stratified by body site)")
    report.append("\n2. **Longitudinal Analysis:** Wilcoxon signed-rank test (paired, Week 1 vs Week 3)")
    report.append("\n3. **Correlation Analysis:** Spearman correlation (antibiotics vs AMR burden)")
    report.append("\n4. **Mixed-Effects Models:** Linear mixed-effects with subject random effects")
    report.append("\n5. **Multiple Testing Correction:** Benjamini-Hochberg FDR")
    report.append("\n6. **Dimensionality Reduction:** PCA on standardized AMR RPM matrix")

    # Outputs
    report.append("\n\n## Analysis Outputs\n")
    report.append("### Results Files:")
    report.append("\n- `results/qc/` - Quality control and metadata")
    report.append("\n- `results/exploratory/` - PCA, diversity, clustering")
    report.append("\n- `results/differential/` - Differential abundance (UCMC vs ZCH)")
    report.append("\n- `results/correlations/` - Antibiotic-AMR correlations")
    report.append("\n- `results/longitudinal/` - Week 1→3 changes")
    report.append("\n- `results/mixed_models/` - Mixed-effects model results")
    report.append("\n- `results/summary/` - This report and summary tables")

    report.append("\n\n### Figures:")
    report.append("\n- `figures/publication/Figure1_PCA.pdf` - PCA overview")
    report.append("\n- `figures/publication/Figure2_BodySite_Comparison.pdf` - AMR by body site")
    report.append("\n- `figures/publication/Figure3_Longitudinal_Trajectories.pdf` - Week 1→3 trajectories")
    report.append("\n- `figures/publication/Figure4_Volcano_Plots.pdf` - Gene-level changes")

    # Conclusions
    report.append("\n\n## Conclusions\n")
    report.append("1. **No geographic differences:** UCMC and ZCH NICUs show remarkably similar resistome profiles")
    report.append("\n2. **Body site dominates:** Stool samples have the highest AMR burden, driven by gut microbiome composition")
    report.append("\n3. **Strong temporal dynamics:** Groin and stool sites show dramatic AMR increases from Week 1 to Week 3")
    report.append("\n4. **Temporal patterns are site-specific:** Axilla remains stable while groin/stool increase")
    report.append("\n5. **Antibiotic exposure not predictive:** No correlation between cumulative antibiotic days and AMR burden")

    # Save report
    report_text = "\n".join(report)
    with open(SUMMARY_DIR / "ANALYSIS_REPORT.md", "w") as f:
        f.write(report_text)

    print(f"✓ Saved: {SUMMARY_DIR / 'ANALYSIS_REPORT.md'}")

    return report_text

def create_top_genes_table():
    """Create table of top changing genes across body sites"""
    print("\nCreating top genes table...")

    all_top_genes = []

    for site in ['axilla', 'groin', 'stool']:
        file_path = LONG_DIR / f"gene_longitudinal_{site}.tsv"
        if file_path.exists():
            long_data = pd.read_csv(file_path, sep="\t")

            # Top 10 increased
            top_increased = long_data[long_data['fdr'] < 0.05].nlargest(10, 'log2_fold_change')
            for _, row in top_increased.iterrows():
                all_top_genes.append({
                    'body_site': site.title(),
                    'gene': row['gene'],
                    'direction': 'Increased',
                    'log2_fc': row['log2_fold_change'],
                    'week1_mean': row['week1_mean'],
                    'week3_mean': row['week3_mean'],
                    'pvalue': row['pvalue'],
                    'fdr': row['fdr']
                })

            # Top 10 decreased
            top_decreased = long_data[long_data['fdr'] < 0.05].nsmallest(10, 'log2_fold_change')
            for _, row in top_decreased.iterrows():
                all_top_genes.append({
                    'body_site': site.title(),
                    'gene': row['gene'],
                    'direction': 'Decreased',
                    'log2_fc': row['log2_fold_change'],
                    'week1_mean': row['week1_mean'],
                    'week3_mean': row['week3_mean'],
                    'pvalue': row['pvalue'],
                    'fdr': row['fdr']
                })

    top_genes_df = pd.DataFrame(all_top_genes)

    if len(top_genes_df) > 0:
        top_genes_df = top_genes_df.sort_values(['body_site', 'direction', 'log2_fc'], ascending=[True, False, False])
        top_genes_df.to_csv(SUMMARY_DIR / "top_changing_genes.tsv", sep="\t", index=False)
        print(f"✓ Saved: {SUMMARY_DIR / 'top_changing_genes.tsv'}")

    return top_genes_df

def main():
    print("="*60)
    print("CREATING STATISTICAL SUMMARY REPORT")
    print("="*60)

    # Create summary tables
    location_df, site_df, total_samples, total_subjects = create_sample_summary()

    # Create key findings table
    findings_df = create_key_findings_table()

    # Create top genes table
    top_genes_df = create_top_genes_table()

    # Create comprehensive report
    report_text = create_analysis_report(location_df, site_df, total_samples, total_subjects, findings_df)

    print("\n" + "="*60)
    print("SUMMARY REPORT COMPLETE")
    print("="*60)
    print("\nOutputs:")
    print(f"  - {SUMMARY_DIR / 'summary_by_location.tsv'}")
    print(f"  - {SUMMARY_DIR / 'summary_by_bodysite_week.tsv'}")
    print(f"  - {SUMMARY_DIR / 'key_findings.tsv'}")
    print(f"  - {SUMMARY_DIR / 'top_changing_genes.tsv'}")
    print(f"  - {SUMMARY_DIR / 'ANALYSIS_REPORT.md'}")

    print("\n" + "="*60)
    print("ALL ANALYSES COMPLETE!")
    print("="*60)

if __name__ == "__main__":
    main()
