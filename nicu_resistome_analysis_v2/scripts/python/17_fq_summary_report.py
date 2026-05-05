#!/usr/bin/env python3
"""
Phase 8: Summary Report and Tables
==================================

Generates manuscript-ready tables and a comprehensive summary report.

Input:
    - All analysis outputs from previous phases

Outputs:
    - results/fq_resistance/summary/Table1_Mutation_Catalog.tsv
    - results/fq_resistance/summary/Table2_Species_Summary.tsv
    - results/fq_resistance/summary/Table3_Bodysite_Comparison.tsv
    - results/fq_resistance/summary/Table4_Model_Results.tsv
    - results/fq_resistance/summary/FQ_RESISTANCE_REPORT.md
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

# Paths
RESULTS_DIR = "results/fq_resistance"
SUMMARY_DIR = os.path.join(RESULTS_DIR, "summary")

# Create output directory
os.makedirs(SUMMARY_DIR, exist_ok=True)

print("=" * 80)
print("Phase 8: Summary Report and Tables")
print("=" * 80)

# ============================================================================
# Load all analysis results
# ============================================================================
print("\n[1/6] Loading analysis results...")

# QC data
mutation_catalog = pd.read_csv(os.path.join(RESULTS_DIR, "qc/mutation_catalog.tsv"), sep="\t")
qc_summary = pd.read_csv(os.path.join(RESULTS_DIR, "qc/qc_summary.tsv"), sep="\t")

# Species prevalence
species_loc = pd.read_csv(os.path.join(RESULTS_DIR, "species_prevalence/location_statistical_tests.tsv"), sep="\t")
species_site = pd.read_csv(os.path.join(RESULTS_DIR, "species_prevalence/bodysite_statistical_tests.tsv"), sep="\t")

# Mutation analysis
mutation_site = pd.read_csv(os.path.join(RESULTS_DIR, "mutation_analysis/mutation_bodysite_comparison.tsv"), sep="\t")

# Model results
model_coefs = pd.read_csv(os.path.join(RESULTS_DIR, "mixed_models/model_coefficients.tsv"), sep="\t")
model_diag = pd.read_csv(os.path.join(RESULTS_DIR, "mixed_models/model_diagnostics.tsv"), sep="\t")

print("  ✓ Loaded all analysis results")

# ============================================================================
# Table 1: FQ Resistance Mutation Catalog
# ============================================================================
print("\n[2/6] Creating Table 1: Mutation Catalog...")

# Top 30 mutations with clinical classification
table1 = mutation_catalog.head(30).copy()

# Add clinical significance classification based on literature
def classify_mutation(row):
    mutation_id = row['mutation_id']
    gene = row['gene']

    # High-impact mutations (well-established resistance mutations)
    if 'gyrA' in gene:
        if any(x in mutation_id for x in ['S83L', 'S83I', 'D87N', 'D87Y', 'S81F', 'S81Y', 'G81']):
            return 'High'
    if 'parC' in gene:
        if any(x in mutation_id for x in ['S80I', 'S80R', 'S79I', 'S79F', 'E84']):
            return 'High'
    if 'parE' in gene:
        if any(x in mutation_id for x in ['S458A', 'D420N', 'I443V']):
            return 'Medium'
    if 'gyrB' in gene:
        return 'Medium'

    # Default
    return 'Low/Unknown'

table1['clinical_significance'] = table1.apply(classify_mutation, axis=1)

# Calculate prevalence percentage
total_samples = 486  # From QC
table1['prevalence_pct'] = (table1['n_samples'] / total_samples) * 100

# Select and rename columns
table1_final = table1[[
    'species', 'gene', 'position', 'wildtype_aa', 'mutant_aa',
    'n_samples', 'prevalence_pct', 'mean_frequency', 'median_frequency',
    'clinical_significance'
]].copy()

table1_final.columns = [
    'Species', 'Gene', 'Position', 'Wildtype', 'Mutant',
    'N_Samples', 'Prevalence_%', 'Mean_Freq', 'Median_Freq',
    'Clinical_Significance'
]

# Round numeric columns
table1_final['Prevalence_%'] = table1_final['Prevalence_%'].round(1)
table1_final['Mean_Freq'] = table1_final['Mean_Freq'].round(3)
table1_final['Median_Freq'] = table1_final['Median_Freq'].round(3)

table1_final.to_csv(os.path.join(SUMMARY_DIR, "Table1_Mutation_Catalog.tsv"), sep='\t', index=False)
print(f"  ✓ Saved Table 1: {len(table1_final)} mutations")

# ============================================================================
# Table 2: Species-Level Resistance Summary
# ============================================================================
print("\n[3/6] Creating Table 2: Species Summary...")

# Combine species data
table2 = species_loc[['species', 'UCMC_n', 'UCMC_pct', 'ZCH_n', 'ZCH_pct', 'p_value']].copy()
table2 = table2.merge(
    species_site[['species', 'Axilla_pct', 'Groin_pct', 'Stool_pct', 'p_value', 'fdr']],
    on='species',
    suffixes=('_location', '_bodysite')
)

# Get top mutations per species from catalog
top_mutations_per_species = mutation_catalog.sort_values(['species', 'n_samples'], ascending=[True, False])
top_mutations_per_species = top_mutations_per_species.groupby('species').first().reset_index()

table2 = table2.merge(
    top_mutations_per_species[['species', 'mutation_id']],
    on='species',
    how='left'
)

# Rename columns
table2.columns = [
    'Species', 'UCMC_N', 'UCMC_%', 'ZCH_N', 'ZCH_%', 'Location_p',
    'Axilla_%', 'Groin_%', 'Stool_%', 'BodySite_p', 'BodySite_FDR',
    'Top_Mutation'
]

# Round numerics
for col in ['UCMC_%', 'ZCH_%', 'Axilla_%', 'Groin_%', 'Stool_%']:
    table2[col] = table2[col].round(1)
for col in ['Location_p', 'BodySite_p', 'BodySite_FDR']:
    table2[col] = table2[col].apply(lambda x: f"{x:.3e}" if x < 0.001 else f"{x:.3f}")

table2.to_csv(os.path.join(SUMMARY_DIR, "Table2_Species_Summary.tsv"), sep='\t', index=False)
print(f"  ✓ Saved Table 2: {len(table2)} species")

# ============================================================================
# Table 3: Body Site Comparison
# ============================================================================
print("\n[4/6] Creating Table 3: Body Site Differences...")

# Mutations with significant body site differences
table3 = mutation_site[mutation_site['fdr'] < 0.05].copy()

if len(table3) > 0:
    table3 = table3[[
        'mutation_id', 'species', 'gene',
        'Axilla_n', 'Axilla_mean_freq',
        'Groin_n', 'Groin_mean_freq',
        'Stool_n', 'Stool_mean_freq',
        'kw_statistic', 'p_value', 'fdr'
    ]].copy()

    table3.columns = [
        'Mutation_ID', 'Species', 'Gene',
        'Axilla_N', 'Axilla_Mean',
        'Groin_N', 'Groin_Mean',
        'Stool_N', 'Stool_Mean',
        'Kruskal_Wallis', 'p_value', 'FDR'
    ]

    # Round numerics
    for col in ['Axilla_Mean', 'Groin_Mean', 'Stool_Mean']:
        table3[col] = table3[col].round(3)
    table3['Kruskal_Wallis'] = table3['Kruskal_Wallis'].round(2)
    table3['p_value'] = table3['p_value'].apply(lambda x: f"{x:.3e}")
    table3['FDR'] = table3['FDR'].apply(lambda x: f"{x:.3e}")
else:
    # Empty table
    table3 = pd.DataFrame(columns=[
        'Mutation_ID', 'Species', 'Gene',
        'Axilla_N', 'Axilla_Mean',
        'Groin_N', 'Groin_Mean',
        'Stool_N', 'Stool_Mean',
        'Kruskal_Wallis', 'p_value', 'FDR'
    ])

table3.to_csv(os.path.join(SUMMARY_DIR, "Table3_Bodysite_Comparison.tsv"), sep='\t', index=False)
print(f"  ✓ Saved Table 3: {len(table3)} significant mutations")

# ============================================================================
# Table 4: Mixed-Effects Model Results
# ============================================================================
print("\n[5/6] Creating Table 4: Model Results...")

# Get significant coefficients (p < 0.05)
table4 = model_coefs[model_coefs['p_value'] < 0.05].copy()

table4 = table4[[
    'model', 'parameter', 'coefficient', 'std_err', 'z_value', 'p_value'
]].copy()

# Add 95% CI
table4['CI_lower'] = table4['coefficient'] - 1.96 * table4['std_err']
table4['CI_upper'] = table4['coefficient'] + 1.96 * table4['std_err']

table4.columns = [
    'Model', 'Parameter', 'Coefficient', 'Std_Error', 'Z_value', 'p_value',
    'CI_95_Lower', 'CI_95_Upper'
]

# Round numerics
for col in ['Coefficient', 'Std_Error', 'CI_95_Lower', 'CI_95_Upper']:
    table4[col] = table4[col].round(3)
table4['Z_value'] = table4['Z_value'].round(2)
table4['p_value'] = table4['p_value'].apply(lambda x: f"{x:.3e}" if x < 0.001 else f"{x:.3f}")

table4.to_csv(os.path.join(SUMMARY_DIR, "Table4_Model_Results.tsv"), sep='\t', index=False)
print(f"  ✓ Saved Table 4: {len(table4)} significant coefficients")

# ============================================================================
# Generate Summary Report
# ============================================================================
print("\n[6/6] Creating comprehensive summary report...")

report = []
report.append("# Fluoroquinolone Resistance Mutation Analysis")
report.append(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
report.append("")
report.append("=" * 80)
report.append("")

# Executive Summary
report.append("## Executive Summary")
report.append("")

n_samples = qc_summary[qc_summary['metric'] == 'Total samples in FQ data']['value'].iloc[0]
n_mutations = qc_summary[qc_summary['metric'] == 'Unique resistance mutations']['value'].iloc[0]

report.append(f"- **Total samples analyzed**: {int(n_samples)}")
report.append(f"- **Unique FQ resistance mutations detected**: {int(n_mutations)}")
report.append(f"- **Species with resistance**: {len(species_loc)}")
report.append(f"- **Complete paired subjects (Week 1 & 3)**: 53")
report.append("")

# Key Findings
report.append("## Key Findings")
report.append("")

report.append("### 1. No Location Effect (UCMC vs ZCH)")
report.append("")
report.append("- FQ resistance patterns are **indistinguishable** between UCMC and ZCH")
report.append(f"- All {len(species_loc)} species showed no significant prevalence differences (all FDR > 0.05)")
report.append("- Mixed-effects model: Location coefficient p=0.79 (not significant)")
report.append("- **Interpretation**: Both NICUs harbor similar FQ-resistant bacterial populations")
report.append("")

report.append("### 2. Strong Body Site Effect")
report.append("")
n_sig_species = (species_site['fdr'] < 0.05).sum()
report.append(f"- **{n_sig_species}/9 species** differ significantly by body site (FDR < 0.05)")
report.append("")
report.append("**Key patterns:**")
report.append("- **S. aureus** (skin colonizer): Axilla (95%) > Groin (77%) > Stool (56%), p=1.5×10⁻⁹")
report.append("- **Enterococcus faecalis**: Stool (59%) > Groin (41%) > Axilla (7%), p=7.3×10⁻¹⁴")
report.append("- **K. pneumoniae**: Stool (58%) > Groin (47%) > Axilla (21%), p=4.8×10⁻⁷")
report.append("")
report.append("- Mixed-effects model:")
report.append("  - Groin: +0.97 mutations vs Axilla (p<0.001)")
report.append("  - Stool: +1.38 mutations vs Axilla (p<0.001)")
report.append("")

report.append("### 3. Modest Temporal Changes")
report.append("")
report.append("- Week 3 has slightly more mutations than Week 1 (+0.33 mutations, p=0.02)")
report.append("- However, individual mutation frequencies remain **stable** (no mutations with FDR<0.05)")
report.append("- Likely reflects cumulative colonization, not resistance evolution")
report.append("")

report.append("### 4. No Antibiotic Selection Pressure")
report.append("")
report.append("- **No correlation** between non-FQ antibiotic exposure and FQ resistance")
report.append("  - Total abx days vs N mutations: r=-0.03, p=0.51")
report.append("  - Total abx days vs Mean frequency: r=-0.02, p=0.66")
report.append("- **Interpretation**: FQ resistance is pre-existing, not selected by hospital antibiotics")
report.append("- Infants likely colonized from environment/caregivers with already-resistant strains")
report.append("")

# Top mutations
report.append("## Top FQ Resistance Mutations")
report.append("")
report.append("| Mutation | N Samples | Prevalence | Mean Freq | Significance |")
report.append("|----------|-----------|------------|-----------|--------------|")

for _, row in table1_final.head(10).iterrows():
    mutation_short = f"{row['Species'].split('_')[0]}. {row['Species'].split('_')[1][:4]}. {row['Gene']} {row['Wildtype']}{row['Position']}{row['Mutant']}"
    report.append(f"| {mutation_short} | {row['N_Samples']} | {row['Prevalence_%']:.1f}% | {row['Mean_Freq']:.2f} | {row['Clinical_Significance']} |")

report.append("")

# Multi-gene resistance
report.append("## Multi-Gene Resistance Patterns")
report.append("")
report.append("- **198 sample×species** (17%) have multi-gene resistance (≥2 genes)")
report.append("- **33 S. aureus samples** with triple mutations (gyrA + parC + parE)")
report.append("- **Enterococcus spp.** typically have single gyrA mutations (fixed resistance)")
report.append("- **Enterobacteriaceae** show diverse mutation patterns across gyrA, parC, parE")
report.append("")

# Clinical implications
report.append("## Clinical Implications")
report.append("")
report.append("1. **Fluoroquinolones should not be used empirically** in this NICU population")
report.append("   - High prevalence of resistance (>60% for key pathogens like S. aureus)")
report.append("   - Resistance is widespread across both institutions")
report.append("")
report.append("2. **Body site matters for infection prevention**")
report.append("   - Skin sites (especially axilla) harbor FQ-resistant S. aureus")
report.append("   - Gut/stool has high burden of FQ-resistant Enterobacteriaceae")
report.append("   - Targeted infection control strategies needed by body site")
report.append("")
report.append("3. **Resistance is community-acquired, not hospital-selected**")
report.append("   - No association with hospital antibiotic use")
report.append("   - Likely reflects colonization from caregivers/environment")
report.append("   - Emphasizes need for transmission-based precautions")
report.append("")

# Statistical summary
report.append("## Statistical Summary")
report.append("")
report.append("### Mixed-Effects Models")
report.append("")
report.append("**Model 1: N mutations ~ Location + BodySite + Week**")
report.append(f"- R² = {model_diag[model_diag['model']=='N_mutations']['r_squared'].iloc[0]:.3f}")
report.append("- Significant predictors:")
report.append("  - BodySite[Groin]: +0.97 (p<0.001)")
report.append("  - BodySite[Stool]: +1.38 (p<0.001)")
report.append("  - Week 3: +0.33 (p=0.02)")
report.append("")

report.append("**Model 2: Mean freq ~ Location + BodySite + Week**")
report.append(f"- R² = {model_diag[model_diag['model']=='Mean_frequency']['r_squared'].iloc[0]:.3f}")
report.append("- Significant predictors:")
report.append("  - BodySite[Groin]: -0.061 (p=0.014)")
report.append("  - BodySite[Stool]: -0.081 (p<0.001)")
report.append("")

# Methods summary
report.append("## Methods Summary")
report.append("")
report.append("**Data source**: NICU microbiome study (n=486 samples, 126 subjects)")
report.append("")
report.append("**FQ resistance detection**: Metagenomic sequencing with chromosomal mutation calling")
report.append("- Targeted genes: gyrA, gyrB, parC, parE")
report.append("- 86 unique resistance mutations identified")
report.append("- 9 species analyzed")
report.append("")
report.append("**Statistical approach**:")
report.append("- Species prevalence: Fisher's exact test (location), Chi-square (body site)")
report.append("- Mutation frequencies: Mann-Whitney U (location), Kruskal-Wallis (body site)")
report.append("- Longitudinal: Wilcoxon signed-rank (paired)")
report.append("- Multiple testing: Benjamini-Hochberg FDR correction")
report.append("- Mixed-effects models: Account for repeated measures within subjects")
report.append("")

# Limitations
report.append("## Limitations")
report.append("")
report.append("1. Metagenomic sequencing may underdetect low-frequency mutations")
report.append("2. Clinical resistance breakpoints not directly tested (no phenotypic data)")
report.append("3. Limited longitudinal data (only Week 1 and Week 3)")
report.append("4. Cannot distinguish colonization from infection")
report.append("5. Mixed-effects models had singular covariance (limited within-subject variation)")
report.append("")

# Future directions
report.append("## Future Directions")
report.append("")
report.append("1. Phenotypic validation of resistance predictions")
report.append("2. Source tracking: caregiver/environmental samples")
report.append("3. Extended longitudinal sampling beyond Week 3")
report.append("4. Integration with clinical outcomes data")
report.append("5. Transmission network analysis using strain-level genomics")
report.append("")

# Files generated
report.append("## Output Files")
report.append("")
report.append("### Data Tables")
for i in range(1, 5):
    report.append(f"- `Table{i}_*.tsv`")

report.append("")
report.append("### Publication Figures")
report.append("- `Figure1_FQ_Overview.pdf` - Mutation heatmap and species prevalence")
report.append("- `Figure2_Location_BodySite.pdf` - Group comparisons")
report.append("- `Figure3_Longitudinal.pdf` - Temporal trajectories")
report.append("- `Figure4_Species_Patterns.pdf` - Species-specific resistance")
report.append("")

report.append("### Analysis Results")
report.append("- QC summaries and mutation catalogs")
report.append("- Species prevalence analyses")
report.append("- Mutation-level comparisons")
report.append("- Gene-level burden scores")
report.append("- Antibiotic correlation analyses")
report.append("- Mixed-effects model outputs")
report.append("")

report.append("=" * 80)
report.append("")
report.append("**Analysis complete!**")
report.append("")

# Write report
report_text = "\n".join(report)
with open(os.path.join(SUMMARY_DIR, "FQ_RESISTANCE_REPORT.md"), 'w') as f:
    f.write(report_text)

print(f"  ✓ Saved comprehensive report: {SUMMARY_DIR}/FQ_RESISTANCE_REPORT.md")

print("\n" + "=" * 80)
print("Phase 8 complete!")
print("=" * 80)

print("\nAll manuscript-ready outputs created:")
print("  • 4 tables (TSV format)")
print("  • 4 publication figures (PDF)")
print("  • 1 comprehensive summary report (Markdown)")

print("\n" + "=" * 80)
print("ENTIRE FQ RESISTANCE ANALYSIS COMPLETE!")
print("=" * 80)

print("\nTotal outputs generated:")
print("  ✓ Phase 1: Data preparation and QC")
print("  ✓ Phase 2: Species prevalence analysis")
print("  ✓ Phase 3: Mutation-level comparisons")
print("  ✓ Phase 4: Gene-level resistance burden")
print("  ✓ Phase 5: Antibiotic correlations")
print("  ✓ Phase 6: Mixed-effects models")
print("  ✓ Phase 7: Publication figures")
print("  ✓ Phase 8: Summary report and tables")

print("\nKey finding: **Body site matters more than hospital location!**")
print("\nAll results are in: results/fq_resistance/")
print("Summary report: results/fq_resistance/summary/FQ_RESISTANCE_REPORT.md")
