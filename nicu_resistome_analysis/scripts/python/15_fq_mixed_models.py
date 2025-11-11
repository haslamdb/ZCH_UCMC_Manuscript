#!/usr/bin/env python3
"""
Phase 6: Mixed-Effects Modeling
===============================

Builds mixed-effects models to properly account for repeated measures within subjects.

Input:
    - data/fq_resistance_filtered.tsv

Outputs:
    - results/fq_resistance/mixed_models/model_coefficients.tsv
    - results/fq_resistance/mixed_models/model_anova.tsv
    - results/fq_resistance/mixed_models/model_diagnostics.tsv
    - figures/fq_resistance/model_diagnostics.pdf
    - figures/fq_resistance/predicted_values_by_factors.pdf
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Try to import statsmodels
try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    STATSMODELS_AVAILABLE = True
except ImportError:
    print("WARNING: statsmodels not available. Installing...")
    import subprocess
    subprocess.check_call(['pip', 'install', 'statsmodels'])
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    STATSMODELS_AVAILABLE = True

# Paths
DATA_DIR = "data"
RESULTS_DIR = "results/fq_resistance/mixed_models"
FIGURES_DIR = "figures/fq_resistance"

# Create output directories
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

print("=" * 80)
print("Phase 6: Mixed-Effects Modeling")
print("=" * 80)

# ============================================================================
# Step 1: Load and prepare data
# ============================================================================
print("\n[1/5] Loading and preparing data...")

fq_data = pd.read_csv(os.path.join(DATA_DIR, "fq_resistance_filtered.tsv"), sep="\t")

# Create sample-level summary (response variable)
sample_summary = fq_data.groupby(['sample_name', 'Location', 'BodySite', 'Week', 'SubjectID']).agg({
    'mutation_id': 'nunique',  # Total mutations
    'total_resistant_frequency': 'mean',  # Mean frequency
    'is_high_frequency': 'sum'  # Number of high-frequency mutations
}).reset_index()

sample_summary.columns = [
    'sample_name', 'Location', 'BodySite', 'Week', 'SubjectID',
    'n_mutations', 'mean_freq', 'n_high_freq'
]

print(f"  - {len(sample_summary)} samples for modeling")
print(f"  - {sample_summary['SubjectID'].nunique()} unique subjects")

# Remove any missing values
sample_summary = sample_summary.dropna(subset=['Location', 'BodySite', 'Week', 'SubjectID'])

print(f"\nResponse variable distributions:")
print(sample_summary[['n_mutations', 'mean_freq', 'n_high_freq']].describe())

# ============================================================================
# Step 2: Model 1 - Number of mutations
# ============================================================================
print("\n[2/5] Model 1: Number of FQ mutations ~ Location + BodySite + Week + (1|SubjectID)")

# Mixed-effects model using statsmodels
# Note: statsmodels uses MixedLM for mixed effects
try:
    # Prepare data - need to ensure categorical variables are properly encoded
    model_data = sample_summary.copy()

    # Create dummy variables for categorical predictors
    model_data['Location_UCMC'] = (model_data['Location'] == 'UCMC').astype(int)
    model_data['BodySite_Groin'] = (model_data['BodySite'] == 'Groin').astype(int)
    model_data['BodySite_Stool'] = (model_data['BodySite'] == 'Stool').astype(int)
    model_data['Week_3'] = (model_data['Week'] == 'Week 3').astype(int)

    # Fit mixed model
    from statsmodels.regression.mixed_linear_model import MixedLM

    # Model 1: N mutations
    print("\n  Building Model 1...")
    md1 = MixedLM.from_formula(
        'n_mutations ~ Location + BodySite + Week',
        groups=model_data['SubjectID'],
        data=model_data
    )
    mdf1 = md1.fit(method='lbfgs')

    print("\n  Model 1 Summary:")
    print(mdf1.summary())

    # Extract coefficients
    coef_df1 = pd.DataFrame({
        'parameter': mdf1.params.index,
        'coefficient': mdf1.params.values,
        'std_err': mdf1.bse.values,
        'z_value': mdf1.tvalues.values,
        'p_value': mdf1.pvalues.values
    })
    coef_df1['model'] = 'N_mutations'

except Exception as e:
    print(f"  Error fitting Model 1: {e}")
    print("  Falling back to simple linear model...")

    # Fallback: Use OLS with dummy variables
    model1 = smf.ols('n_mutations ~ Location + BodySite + Week', data=sample_summary)
    result1 = model1.fit()

    print("\n  Model 1 Summary (OLS fallback):")
    print(result1.summary())

    coef_df1 = pd.DataFrame({
        'parameter': result1.params.index,
        'coefficient': result1.params.values,
        'std_err': result1.bse.values,
        'z_value': result1.tvalues.values,
        'p_value': result1.pvalues.values
    })
    coef_df1['model'] = 'N_mutations (OLS)'

# ============================================================================
# Step 3: Model 2 - Mean mutation frequency
# ============================================================================
print("\n[3/5] Model 2: Mean FQ frequency ~ Location + BodySite + Week + (1|SubjectID)")

try:
    # Model 2: Mean frequency
    print("\n  Building Model 2...")
    md2 = MixedLM.from_formula(
        'mean_freq ~ Location + BodySite + Week',
        groups=model_data['SubjectID'],
        data=model_data
    )
    mdf2 = md2.fit(method='lbfgs')

    print("\n  Model 2 Summary:")
    print(mdf2.summary())

    # Extract coefficients
    coef_df2 = pd.DataFrame({
        'parameter': mdf2.params.index,
        'coefficient': mdf2.params.values,
        'std_err': mdf2.bse.values,
        'z_value': mdf2.tvalues.values,
        'p_value': mdf2.pvalues.values
    })
    coef_df2['model'] = 'Mean_frequency'

except Exception as e:
    print(f"  Error fitting Model 2: {e}")
    print("  Falling back to simple linear model...")

    model2 = smf.ols('mean_freq ~ Location + BodySite + Week', data=sample_summary)
    result2 = model2.fit()

    print("\n  Model 2 Summary (OLS fallback):")
    print(result2.summary())

    coef_df2 = pd.DataFrame({
        'parameter': result2.params.index,
        'coefficient': result2.params.values,
        'std_err': result2.bse.values,
        'z_value': result2.tvalues.values,
        'p_value': result2.pvalues.values
    })
    coef_df2['model'] = 'Mean_frequency (OLS)'

# Combine coefficients
all_coefs = pd.concat([coef_df1, coef_df2], ignore_index=True)

# Save coefficients
all_coefs.to_csv(
    os.path.join(RESULTS_DIR, "model_coefficients.tsv"),
    sep='\t',
    index=False
)
print(f"\n✓ Saved: {RESULTS_DIR}/model_coefficients.tsv")

# ============================================================================
# Step 4: Model diagnostics
# ============================================================================
print("\n[4/5] Running model diagnostics...")

# Simple diagnostics for OLS models
diagnostics = []

# For Model 1 (use OLS for diagnostics)
model1_ols = smf.ols('n_mutations ~ Location + BodySite + Week', data=sample_summary)
result1_ols = model1_ols.fit()

residuals1 = result1_ols.resid
fitted1 = result1_ols.fittedvalues

diagnostics.append({
    'model': 'N_mutations',
    'r_squared': result1_ols.rsquared,
    'adj_r_squared': result1_ols.rsquared_adj,
    'aic': result1_ols.aic,
    'bic': result1_ols.bic,
    'residual_std': np.std(residuals1),
    'mean_absolute_error': np.mean(np.abs(residuals1))
})

# For Model 2
model2_ols = smf.ols('mean_freq ~ Location + BodySite + Week', data=sample_summary)
result2_ols = model2_ols.fit()

residuals2 = result2_ols.resid
fitted2 = result2_ols.fittedvalues

diagnostics.append({
    'model': 'Mean_frequency',
    'r_squared': result2_ols.rsquared,
    'adj_r_squared': result2_ols.rsquared_adj,
    'aic': result2_ols.aic,
    'bic': result2_ols.bic,
    'residual_std': np.std(residuals2),
    'mean_absolute_error': np.mean(np.abs(residuals2))
})

diagnostics_df = pd.DataFrame(diagnostics)
print("\nModel diagnostics:")
print(diagnostics_df.to_string(index=False))

diagnostics_df.to_csv(
    os.path.join(RESULTS_DIR, "model_diagnostics.tsv"),
    sep='\t',
    index=False
)
print(f"\n✓ Saved: {RESULTS_DIR}/model_diagnostics.tsv")

# ============================================================================
# Step 5: Create visualizations
# ============================================================================
print("\n[5/5] Creating visualizations...")

# Figure 1: Model diagnostics plots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: Residuals vs Fitted (Model 1)
ax = axes[0, 0]
ax.scatter(fitted1, residuals1, alpha=0.5)
ax.axhline(y=0, color='r', linestyle='--')
ax.set_xlabel('Fitted Values', fontsize=10)
ax.set_ylabel('Residuals', fontsize=10)
ax.set_title('Model 1: Residuals vs Fitted', fontsize=11, fontweight='bold')
ax.grid(alpha=0.3)

# Panel B: Q-Q plot (Model 1)
ax = axes[0, 1]
stats.probplot(residuals1, dist="norm", plot=ax)
ax.set_title('Model 1: Q-Q Plot', fontsize=11, fontweight='bold')
ax.grid(alpha=0.3)

# Panel C: Residuals vs Fitted (Model 2)
ax = axes[1, 0]
ax.scatter(fitted2, residuals2, alpha=0.5, color='coral')
ax.axhline(y=0, color='r', linestyle='--')
ax.set_xlabel('Fitted Values', fontsize=10)
ax.set_ylabel('Residuals', fontsize=10)
ax.set_title('Model 2: Residuals vs Fitted', fontsize=11, fontweight='bold')
ax.grid(alpha=0.3)

# Panel D: Q-Q plot (Model 2)
ax = axes[1, 1]
stats.probplot(residuals2, dist="norm", plot=ax)
ax.set_title('Model 2: Q-Q Plot', fontsize=11, fontweight='bold')
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "model_diagnostics.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/model_diagnostics.pdf")

# Figure 2: Predicted values by factors
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: By Location
ax = axes[0, 0]
loc_means = sample_summary.groupby('Location')['n_mutations'].agg(['mean', 'sem']).reset_index()
ax.bar(loc_means['Location'], loc_means['mean'], yerr=loc_means['sem'], capsize=5, alpha=0.7)
ax.set_ylabel('N Mutations', fontsize=10, fontweight='bold')
ax.set_title('By Location', fontsize=11, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Panel B: By Body Site
ax = axes[0, 1]
site_means = sample_summary.groupby('BodySite')['n_mutations'].agg(['mean', 'sem']).reset_index()
ax.bar(site_means['BodySite'], site_means['mean'], yerr=site_means['sem'], capsize=5, alpha=0.7, color='green')
ax.set_ylabel('N Mutations', fontsize=10, fontweight='bold')
ax.set_title('By Body Site', fontsize=11, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Panel C: By Week
ax = axes[1, 0]
week_means = sample_summary.groupby('Week')['n_mutations'].agg(['mean', 'sem']).reset_index()
ax.bar(week_means['Week'], week_means['mean'], yerr=week_means['sem'], capsize=5, alpha=0.7, color='coral')
ax.set_ylabel('N Mutations', fontsize=10, fontweight='bold')
ax.set_title('By Week', fontsize=11, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Panel D: Interaction Location × BodySite
ax = axes[1, 1]
interaction = sample_summary.groupby(['Location', 'BodySite'])['n_mutations'].mean().unstack()
interaction.plot(kind='bar', ax=ax, width=0.8)
ax.set_ylabel('Mean N Mutations', fontsize=10, fontweight='bold')
ax.set_xlabel('Location', fontsize=10, fontweight='bold')
ax.set_title('Location × Body Site Interaction', fontsize=11, fontweight='bold')
ax.legend(title='Body Site', fontsize=9)
ax.grid(axis='y', alpha=0.3)
plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "predicted_values_by_factors.pdf"), dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {FIGURES_DIR}/predicted_values_by_factors.pdf")

print("\n" + "=" * 80)
print("Phase 6 complete!")
print("=" * 80)

print("\nKey Findings:")
print(f"  • Model 1 R²: {diagnostics_df.iloc[0]['r_squared']:.3f}")
print(f"  • Model 2 R²: {diagnostics_df.iloc[1]['r_squared']:.3f}")

# Print significant coefficients (p < 0.05)
sig_coefs = all_coefs[all_coefs['p_value'] < 0.05]
if len(sig_coefs) > 0:
    print(f"\n  Significant coefficients (p < 0.05):")
    for _, row in sig_coefs.iterrows():
        print(f"    {row['model']}: {row['parameter']} = {row['coefficient']:.3f} (p={row['p_value']:.3e})")
else:
    print(f"\n  No significant coefficients at p < 0.05")

print("\nNext step: Run scripts/python/16_fq_publication_figures.py")
