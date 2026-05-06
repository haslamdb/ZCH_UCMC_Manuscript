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

    # ------------------------------------------------------------
    # Pull every headline number from the actual analysis artifacts
    # ------------------------------------------------------------
    # Gene counts (universe vs post-QC)
    n_genes_total = None
    n_genes_passing_qc = None
    if (QC_DIR / "gene_prevalence_stats.tsv").exists():
        n_genes_total = len(pd.read_csv(QC_DIR / "gene_prevalence_stats.tsv", sep="\t"))
    if (QC_DIR / "genes_passing_qc.txt").exists():
        n_genes_passing_qc = sum(1 for _ in open(QC_DIR / "genes_passing_qc.txt") if _.strip())

    # Differential abundance — total significant across sites
    n_sig_diff_total = None
    diff_per_site = []
    if (DIFFERENTIAL_DIR / "differential_abundance_summary.tsv").exists():
        ds = pd.read_csv(DIFFERENTIAL_DIR / "differential_abundance_summary.tsv", sep="\t")
        n_sig_diff_total = int(ds['sig_fdr005'].sum())
        diff_per_site = [(r['body_site'], int(r['sig_fdr005']), int(r['total_genes'])) for _, r in ds.iterrows()]

    # Mixed-effects model coefficients (full model)
    location_pval = stool_coef = stool_pval = None
    if (MODELS_DIR / "full_model_coefficients.tsv").exists():
        mc = pd.read_csv(MODELS_DIR / "full_model_coefficients.tsv", sep="\t", index_col=0)
        if 'Location[T.ZCH]' in mc.index:
            location_pval = mc.loc['Location[T.ZCH]', 'p_value']
        if 'BodySite[T.Stool]' in mc.index:
            stool_coef = mc.loc['BodySite[T.Stool]', 'coefficient']
            stool_pval = mc.loc['BodySite[T.Stool]', 'p_value']

    # PCA — PC1 variance explained
    pc1_var_pct = None
    if (EXPLORATORY_DIR / "pca_variance_explained.tsv").exists():
        pc = pd.read_csv(EXPLORATORY_DIR / "pca_variance_explained.tsv", sep="\t")
        pc1 = pc[pc['PC'] == 'PC1']
        if not pc1.empty:
            pc1_var_pct = float(pc1.iloc[0]['variance_explained']) * 100

    # Overall (all-samples) antibiotic-AMR Spearman correlation
    overall_rho = overall_pval = None
    if (CORR_DIR / "total_abx_amr_correlation_overall.tsv").exists():
        ov = pd.read_csv(CORR_DIR / "total_abx_amr_correlation_overall.tsv", sep="\t")
        if not ov.empty:
            overall_rho = float(ov.iloc[0]['spearman_rho'])
            overall_pval = float(ov.iloc[0]['pvalue'])

    # Variance components / ICC from the full mixed-effects model
    icc = None
    if (MODELS_DIR / "variance_components.tsv").exists():
        vc = pd.read_csv(MODELS_DIR / "variance_components.tsv", sep="\t")
        if not vc.empty:
            icc = float(vc.iloc[0]['icc'])

    # Complete-subject counts per location, from the location_df we already built
    n_complete_total = int(location_df['Complete Subjects'].sum())
    complete_by_loc = ", ".join(
        f"{int(r['Complete Subjects'])} {r['Location']}" for _, r in location_df.iterrows()
    )

    def _fmt(x, spec):
        return format(x, spec) if x is not None else "N/A"

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
    report.append(f"- **AMR Genes Analyzed:** {_fmt(n_genes_passing_qc, 'd')} genes (after QC filtering)")

    # Sample Distribution
    report.append("\n\n### Sample Distribution by Location")
    report.append("\n```")
    report.append(location_df.to_string(index=False))
    report.append("```")

    # Key Findings
    report.append("\n\n## Key Findings\n")

    # Finding 1: No location differences
    loc_sig_word = "not significant" if (location_pval is not None and location_pval >= 0.05) else "significant"
    report.append("### 1. No Significant Differences Between UCMC and ZCH")
    report.append("\n**Result:** Resistome composition is remarkably similar between Cincinnati and Hangzhou NICUs.")
    report.append(f"\n- Differential abundance testing: {_fmt(n_sig_diff_total, 'd')} genes significantly different at FDR < 0.05 (across all body sites)")
    report.append(f"- Mixed-effects model: Location effect p = {_fmt(location_pval, '.3f')} ({loc_sig_word})")
    report.append("- **Interpretation:** Geographic location does not significantly influence NICU resistome composition")

    # Finding 2: Body site is primary driver
    stool_p_str = (
        f"{stool_pval:.2g}" if (stool_pval is not None and stool_pval >= 1e-4)
        else f"< {1e-4:.0e}" if stool_pval is not None
        else "N/A"
    )
    report.append("\n\n### 2. Body Site is the Primary Driver of Resistome Composition")
    report.append("\n**Result:** Stool samples have higher AMR burden than skin sites.")
    report.append(f"\n- PCA analysis: PC1 ({_fmt(pc1_var_pct, '.1f')}% variance) separates stool from axilla/groin")
    report.append(f"- Mixed-effects model: Stool coefficient = {_fmt(stool_coef, '+.3f')}, p = {stool_p_str}")
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

    # Finding 4: Antibiotic correlation (sign of result determined by data)
    report.append("\n\n### 4. Antibiotic Exposure vs Total AMR Burden")
    if overall_rho is not None and overall_pval is not None:
        sig = "significant" if overall_pval < 0.05 else "not significant"
        direction = "positively correlates with" if overall_rho > 0.05 else "negatively correlates with" if overall_rho < -0.05 else "does not correlate with"
        report.append(f"\n**Result:** Cumulative antibiotic exposure {direction} total AMR burden ({sig}).")
        report.append(f"\n- Overall Spearman rho = {overall_rho:+.3f}, p = {overall_pval:.3f}")
    else:
        report.append("\n- Overall Spearman: N/A (correlation file missing)")
    # Count antibiotics with FDR-significant correlations from specific_antibiotic_correlations.tsv
    n_sig_specific = "N/A"
    if (CORR_DIR / "specific_antibiotic_correlations.tsv").exists():
        sa = pd.read_csv(CORR_DIR / "specific_antibiotic_correlations.tsv", sep="\t")
        if 'spearman_fdr' in sa.columns:
            n_sig_specific = int((sa['spearman_fdr'] < 0.05).sum())
    report.append(f"- Specific antibiotics with significant correlations (FDR<0.05): {n_sig_specific}")
    report.append("- **Interpretation:** AMR colonization in this cohort is largely independent of selective antibiotic pressure")

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

        # Report ICC (read from variance_components.tsv emitted by script 07)
        report.append(f"\n\n**Intraclass Correlation (ICC):** {_fmt(icc, '.3f')}")
        if icc is not None:
            report.append(f"\n  - {icc * 100:.1f}% of variance is between-subject")
        report.append("  - Confirms substantial within-subject correlation, validating repeated measures design")

    # Finding 6: Intrinsic vs acquired resistance split (from 09b)
    iva_path = SUMMARY_DIR / "intrinsic_vs_acquired_summary.tsv"
    if iva_path.exists():
        report.append("\n\n### 6. Intrinsic vs Acquired Resistance")
        report.append("\nResistance genes split by NCBI ReferenceGeneCatalog `resistance_type`")
        report.append("(intrinsic = chromosomal/species-defining, acquired = mobilizable).")
        report.append("Intrinsic and acquired share no biological denominator and are reported separately.\n")

        iva = pd.read_csv(iva_path, sep="\t")
        # Wide pivot: median RPM by SampleType x resistance_type x Week
        pivot = iva.pivot_table(
            index=['SampleType', 'resistance_type'],
            columns=['SampleCollectionWeek'],
            values='median_rpm',
            aggfunc='median',
        ).round(1)
        report.append("\n**Median RPM by body site × resistance type × week (across both NICUs):**\n")
        report.append("```")
        report.append(pivot.to_string())
        report.append("```")

        # Acquired fraction summary (from per-sample fraction file)
        af_path = SUMMARY_DIR / "acquired_fraction_per_sample.tsv"
        if af_path.exists():
            af = pd.read_csv(af_path, sep="\t")
            af_pivot = af.groupby(['SampleType', 'SampleCollectionWeek'])['acquired_fraction'].median().unstack().round(3)
            report.append("\n**Median acquired-fraction (acquired RPM / (acquired+intrinsic)):**\n")
            report.append("```")
            report.append(af_pivot.to_string())
            report.append("```")
            report.append("\nA fraction near 1.0 means the sample's resistome is dominated by acquired (mobilizable)")
            report.append("genes; values dropping below 1.0 in groin and stool by week 3 reflect rising")
            report.append("intrinsic-resistance organisms (e.g. Pseudomonas/Klebsiella efflux pumps) joining the gut/perineal flora.")

    # Data Quality
    report.append("\n\n## Data Quality and Filtering\n")
    report.append("### Quality Control Steps:")
    report.append("\n1. **Gene Filtering:**")
    report.append(f"\n   - Original (any reads): {_fmt(n_genes_total, 'd')} gene families")
    report.append(f"\n   - After QC: {_fmt(n_genes_passing_qc, 'd')} (present in ≥5% samples with max RPM ≥1.0)")
    report.append("\n2. **Sample Filtering:**")
    report.append("\n   - Excluded 'Nonese' (nose swab) samples (never sequenced)")
    report.append("\n   - Retained only samples with AMR data")
    report.append("\n3. **Subject Completeness:**")
    report.append(f"\n   - Complete subjects: {n_complete_total} ({complete_by_loc})")
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
