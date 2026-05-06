#!/usr/bin/env python3
"""
Mixed-Effects Models for AMR Analysis

Statistical framework:
- Fixed effects: Location (UCMC vs ZCH) × BodySite (Axilla/Groin/Stool) × Week (1 vs 3)
- Random effects: Subject (to account for repeated measures)
- Response: Total AMR burden (log-transformed for normality)

Uses linear mixed-effects models (LMM) to properly account for:
1. Within-subject correlation (multiple samples per subject)
2. Main effects and interactions
3. Individual subject variation
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import warnings
warnings.filterwarnings('ignore')

# Define paths
PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
RESULTS_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results"
QC_DIR = RESULTS_DIR / "qc"
MODELS_DIR = RESULTS_DIR / "mixed_models"
FIGURES_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "figures" / "mixed_models"

# Create directories
MODELS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Set plotting style
sns.set_style("whitegrid")
sns.set_context("talk")

def load_data():
    """Load metadata with AMR data"""
    print("Loading data...")

    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")

    # Filter to samples with AMR data and valid body sites
    metadata = metadata[
        (metadata['has_amr_data']) &
        (metadata['SampleType'].isin(['Axilla', 'Groin', 'Stool']))
    ].copy()

    print(f"Total samples: {len(metadata)}")
    print(f"Subjects: {metadata['SubjectID'].nunique()}")
    print(f"Complete subjects: {metadata[metadata['subject_complete']]['SubjectID'].nunique()}")

    return metadata

def prepare_model_data(metadata):
    """Prepare data for mixed-effects modeling"""
    print("\nPreparing model data...")

    # Use the full sample set; mixed-effects models handle unbalanced designs
    # via subject random intercepts. Subsetting to "complete" subjects (all 6
    # samples) costs ~80% of the data without buying inferential cleanliness.
    model_data = metadata.copy()

    # Log-transform total AMR (add pseudocount for zeros)
    model_data['log_total_amr'] = np.log10(model_data['total_amr_rpm'] + 1.0)

    # Create numeric week variable
    model_data['Week_num'] = model_data['SampleCollectionWeek'].map({'Week.1': 1, 'Week.3': 3})

    # Ensure categorical variables are properly encoded
    model_data['Location'] = pd.Categorical(model_data['Location'], categories=['UCMC', 'ZCH'])
    model_data['BodySite'] = pd.Categorical(model_data['SampleType'], categories=['Axilla', 'Groin', 'Stool'])
    model_data['Week'] = pd.Categorical(model_data['SampleCollectionWeek'], categories=['Week.1', 'Week.3'])

    print(f"Model data: {len(model_data)} samples from {model_data['SubjectID'].nunique()} subjects")
    print(f"\nSample distribution:")
    print(model_data.groupby(['Location', 'BodySite', 'Week']).size().unstack(fill_value=0))

    return model_data

def fit_full_model(model_data):
    """Fit full mixed-effects model with 3-way interaction"""
    print("\n" + "="*60)
    print("FULL MIXED-EFFECTS MODEL")
    print("="*60)
    print("\nFormula: log_total_amr ~ Location * BodySite * Week")
    print("  Fixed effects: Location, BodySite, Week, and all interactions")
    print("  Random effects: Subject-specific intercepts (via groups parameter)")

    # Fit model using REML
    # Note: Random effects are specified via groups parameter, not in formula
    formula = "log_total_amr ~ Location * BodySite * Week"
    model = smf.mixedlm(formula, model_data, groups=model_data["SubjectID"])
    result = model.fit(reml=True, method='powell')

    print("\n" + "-"*60)
    print("MODEL SUMMARY")
    print("-"*60)
    print(result.summary())

    # Extract key statistics
    params = result.params
    pvalues = result.pvalues

    # Save full results
    results_table = pd.DataFrame({
        'coefficient': params,
        'std_err': result.bse,
        'z_value': result.tvalues,
        'p_value': pvalues,
        'ci_lower': result.conf_int()[0],
        'ci_upper': result.conf_int()[1]
    })

    results_table.to_csv(MODELS_DIR / "full_model_coefficients.tsv", sep="\t")
    print(f"\n✓ Saved: {MODELS_DIR / 'full_model_coefficients.tsv'}")

    # Interpret key effects
    print("\n" + "-"*60)
    print("KEY FINDINGS")
    print("-"*60)

    # Main effects
    print("\nMain Effects:")
    for effect in ['Location[T.ZCH]', 'BodySite[T.Groin]', 'BodySite[T.Stool]', 'Week[T.Week.3]']:
        if effect in pvalues.index:
            p = pvalues[effect]
            coef = params[effect]
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            print(f"  {effect}: coef={coef:+.3f}, p={p:.4f} {sig}")

    # Two-way interactions
    print("\nTwo-way Interactions:")
    two_way = [col for col in pvalues.index if ':' in col and col.count(':') == 1]
    for effect in two_way:
        p = pvalues[effect]
        coef = params[effect]
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"  {effect}: coef={coef:+.3f}, p={p:.4f} {sig}")

    # Three-way interactions
    print("\nThree-way Interactions:")
    three_way = [col for col in pvalues.index if col.count(':') == 2]
    for effect in three_way:
        p = pvalues[effect]
        coef = params[effect]
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"  {effect}: coef={coef:+.3f}, p={p:.4f} {sig}")

    # Random effects
    subject_var = float(result.cov_re.values[0][0])
    residual_var = float(result.scale)
    icc = subject_var / (subject_var + residual_var)
    print(f"\nRandom Effects:")
    print(f"  Subject variance: {subject_var:.4f}")
    print(f"  Residual variance: {residual_var:.4f}")
    print(f"  ICC (intraclass correlation): {icc:.3f}")

    # Persist for the summary report — otherwise 09_create_summary_report has
    # nowhere to read these from and falls back to hard-coded values.
    pd.DataFrame([{
        'subject_variance': subject_var,
        'residual_variance': residual_var,
        'icc': icc,
    }]).to_csv(MODELS_DIR / "variance_components.tsv", sep="\t", index=False)

    return result, results_table

def fit_reduced_models(model_data):
    """Fit reduced models to test specific hypotheses"""
    print("\n" + "="*60)
    print("REDUCED MODELS FOR HYPOTHESIS TESTING")
    print("="*60)

    models = {}

    # Model 1: Location effect (main effect only)
    print("\n1. Location effect (UCMC vs ZCH):")
    formula1 = "log_total_amr ~ Location"
    model1 = smf.mixedlm(formula1, model_data, groups=model_data["SubjectID"])
    result1 = model1.fit(reml=True, method='powell')
    models['location_only'] = result1

    if 'Location[T.ZCH]' in result1.pvalues.index:
        p = result1.pvalues['Location[T.ZCH]']
        coef = result1.params['Location[T.ZCH]']
        print(f"   Coefficient: {coef:+.3f}, p = {p:.4f}")

    # Model 2: BodySite effect
    print("\n2. Body Site effect:")
    formula2 = "log_total_amr ~ BodySite"
    model2 = smf.mixedlm(formula2, model_data, groups=model_data["SubjectID"])
    result2 = model2.fit(reml=True, method='powell')
    models['bodysite_only'] = result2

    for site in ['BodySite[T.Groin]', 'BodySite[T.Stool]']:
        if site in result2.pvalues.index:
            p = result2.pvalues[site]
            coef = result2.params[site]
            print(f"   {site}: {coef:+.3f}, p = {p:.4f}")

    # Model 3: Week effect (longitudinal change)
    print("\n3. Week effect (Week 1 vs Week 3):")
    formula3 = "log_total_amr ~ Week"
    model3 = smf.mixedlm(formula3, model_data, groups=model_data["SubjectID"])
    result3 = model3.fit(reml=True, method='powell')
    models['week_only'] = result3

    if 'Week[T.Week.3]' in result3.pvalues.index:
        p = result3.pvalues['Week[T.Week.3]']
        coef = result3.params['Week[T.Week.3]']
        print(f"   Coefficient: {coef:+.3f}, p = {p:.4f}")

    # Model 4: BodySite × Week (test if longitudinal changes differ by site)
    print("\n4. BodySite × Week interaction:")
    formula4 = "log_total_amr ~ BodySite * Week"
    model4 = smf.mixedlm(formula4, model_data, groups=model_data["SubjectID"])
    result4 = model4.fit(reml=True, method='powell')
    models['bodysite_week'] = result4

    interactions = [col for col in result4.pvalues.index if ':' in col]
    for effect in interactions:
        p = result4.pvalues[effect]
        coef = result4.params[effect]
        print(f"   {effect}: {coef:+.3f}, p = {p:.4f}")

    # Model 5: Location × Week (test if longitudinal changes differ by location)
    print("\n5. Location × Week interaction:")
    formula5 = "log_total_amr ~ Location * Week"
    model5 = smf.mixedlm(formula5, model_data, groups=model_data["SubjectID"])
    result5 = model5.fit(reml=True, method='powell')
    models['location_week'] = result5

    if 'Location[T.ZCH]:Week[T.Week.3]' in result5.pvalues.index:
        p = result5.pvalues['Location[T.ZCH]:Week[T.Week.3]']
        coef = result5.params['Location[T.ZCH]:Week[T.Week.3]']
        print(f"   Location × Week: {coef:+.3f}, p = {p:.4f}")

    return models

def compare_models(models):
    """Compare model fit using AIC/BIC"""
    print("\n" + "="*60)
    print("MODEL COMPARISON")
    print("="*60)

    comparison = []
    for name, result in models.items():
        comparison.append({
            'model': name,
            'AIC': result.aic,
            'BIC': result.bic,
            'log_likelihood': result.llf,
            'n_params': len(result.params)
        })

    comparison_df = pd.DataFrame(comparison).sort_values('AIC')
    print("\n" + comparison_df.to_string(index=False))

    comparison_df.to_csv(MODELS_DIR / "model_comparison.tsv", sep="\t", index=False)
    print(f"\n✓ Saved: {MODELS_DIR / 'model_comparison.tsv'}")

    return comparison_df

def plot_predicted_values(model_data, full_result):
    """Plot predicted values from full model"""
    print("\nCreating predicted value plots...")

    # Get predicted values
    model_data['predicted'] = full_result.fittedvalues
    model_data['residuals'] = full_result.resid

    # Figure 1: Predicted vs Observed
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Overall fit
    ax = axes[0, 0]
    ax.scatter(model_data['predicted'], model_data['log_total_amr'],
               alpha=0.5, s=50, edgecolors='black', linewidths=0.5)
    ax.plot([model_data['log_total_amr'].min(), model_data['log_total_amr'].max()],
            [model_data['log_total_amr'].min(), model_data['log_total_amr'].max()],
            'r--', linewidth=2, label='Perfect fit')
    ax.set_xlabel('Predicted log10(AMR RPM)')
    ax.set_ylabel('Observed log10(AMR RPM)')
    ax.set_title('Model Fit: Predicted vs Observed')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Residuals vs fitted
    ax = axes[0, 1]
    ax.scatter(model_data['predicted'], model_data['residuals'],
               alpha=0.5, s=50, edgecolors='black', linewidths=0.5)
    ax.axhline(y=0, color='r', linestyle='--', linewidth=2)
    ax.set_xlabel('Predicted log10(AMR RPM)')
    ax.set_ylabel('Residuals')
    ax.set_title('Residual Plot')
    ax.grid(True, alpha=0.3)

    # Q-Q plot for residuals
    ax = axes[1, 0]
    from scipy import stats
    stats.probplot(model_data['residuals'], dist="norm", plot=ax)
    ax.set_title('Normal Q-Q Plot of Residuals')
    ax.grid(True, alpha=0.3)

    # Distribution of residuals
    ax = axes[1, 1]
    ax.hist(model_data['residuals'], bins=30, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Residuals')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of Residuals')
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "model_diagnostics.pdf", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'model_diagnostics.pdf'}")
    plt.close()

    # Figure 2: Predicted values by factors
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # By Location
    ax = axes[0]
    location_means = model_data.groupby(['Location', 'Week'])['predicted'].mean().unstack()
    location_means.plot(kind='bar', ax=ax, color=['#2ca02c', '#d62728'])
    ax.set_ylabel('Predicted log10(AMR RPM)')
    ax.set_title('Predicted AMR by Location and Week')
    ax.set_xlabel('Location')
    ax.legend(title='Week')
    ax.grid(True, alpha=0.3, axis='y')

    # By BodySite
    ax = axes[1]
    site_means = model_data.groupby(['BodySite', 'Week'])['predicted'].mean().unstack()
    site_means.plot(kind='bar', ax=ax, color=['#2ca02c', '#d62728'])
    ax.set_ylabel('Predicted log10(AMR RPM)')
    ax.set_title('Predicted AMR by Body Site and Week')
    ax.set_xlabel('Body Site')
    ax.legend(title='Week')
    ax.grid(True, alpha=0.3, axis='y')

    # By Location and BodySite
    ax = axes[2]
    loc_site_means = model_data.groupby(['Location', 'BodySite', 'Week'])['predicted'].mean().unstack()
    loc_site_means.plot(kind='bar', ax=ax, color=['#2ca02c', '#d62728'])
    ax.set_ylabel('Predicted log10(AMR RPM)')
    ax.set_title('Predicted AMR by Location × Body Site')
    ax.set_xlabel('Location - Body Site')
    ax.legend(title='Week')
    ax.grid(True, alpha=0.3, axis='y')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "predicted_values_by_factors.pdf", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'predicted_values_by_factors.pdf'}")
    plt.close()

def main():
    print("="*60)
    print("MIXED-EFFECTS MODELS FOR AMR ANALYSIS")
    print("="*60)

    # Load data
    metadata = load_data()

    # Prepare model data
    model_data = prepare_model_data(metadata)

    # Fit full model
    full_result, results_table = fit_full_model(model_data)

    # Fit reduced models
    reduced_models = fit_reduced_models(model_data)

    # Add full model to comparison
    reduced_models['full_model'] = full_result

    # Compare models
    comparison = compare_models(reduced_models)

    # Plot diagnostics and predicted values
    plot_predicted_values(model_data, full_result)

    print("\n" + "="*60)
    print("MIXED-EFFECTS MODELING COMPLETE")
    print("="*60)
    print("\nOutputs:")
    print(f"  - {MODELS_DIR / 'full_model_coefficients.tsv'}")
    print(f"  - {MODELS_DIR / 'model_comparison.tsv'}")
    print("\nFigures:")
    print(f"  - {FIGURES_DIR / 'model_diagnostics.pdf'}")
    print(f"  - {FIGURES_DIR / 'predicted_values_by_factors.pdf'}")

if __name__ == "__main__":
    main()
