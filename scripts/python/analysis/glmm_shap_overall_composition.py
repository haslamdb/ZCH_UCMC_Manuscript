# glmm_shap_overall_composition.py
# GLMM + SHAP analysis for overall microbiome composition (PC-based)
# Creates combined figures with GLMM on left and SHAP on right

import sys
import os

# Add utils directory to path for microbiome_transform module
script_dir = os.path.dirname(os.path.abspath(__file__))
utils_dir = os.path.join(script_dir, '..', 'utils')
sys.path.insert(0, utils_dir)

import pandas as pd
import numpy as np
import shap
import microbiome_transform as mt
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
import statsmodels.api as sm
from statsmodels.formula.api import mixedlm
import re

# Set font type to TrueType (Type 42) for editable text in PDFs
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# Set Arial font
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica', 'sans-serif']

# Output directory - use absolute path (updated results with absolute values and filtered antibiotics)
OUTPUT_DIR = os.path.join(os.path.dirname(script_dir), '..', '..', 'results_updated')
OUTPUT_DIR = os.path.abspath(OUTPUT_DIR)
os.makedirs(OUTPUT_DIR, exist_ok=True)


def create_comparison_label(feature_name):
    """Convert encoded feature names to readable comparison format"""
    # Handle different feature patterns - for LMM outputs
    if 'C(' in feature_name and '[T.' in feature_name:
        # Handle categorical variable format from statsmodels: C(variable)[T.value]
        match = re.match(r'C\((\w+)\)\[T\.(.+)\]', feature_name)
        if match:
            category = match.group(1)
            value = match.group(2)
        else:
            return feature_name
    elif '_' in feature_name:
        # Handle sklearn dummy variable format: variable_value
        parts = feature_name.split('_')
        # Special handling for variables like BSI_30D_1
        if len(parts) == 3 and parts[1] == '30D':
            category = f"{parts[0]}_{parts[1]}"
            value = parts[2]
        else:
            category = parts[0]
            value = '_'.join(parts[1:])
    else:
        # Simple features or intercept
        if feature_name == 'Intercept':
            return 'Intercept (Baseline)'
        return feature_name

    # Create simplified labels
    label_map = {
        ('SampleType', 'Groin'): 'Groin Sample',
        ('SampleType', 'Stool'): 'Stool Sample',
        ('SampleType', 'Nonese'): 'Nose Sample',
        ('SampleType', 'Nose'): 'Nose Sample',
        ('Location', 'ZCH'): 'ZCH Location',
        ('Location', 'Hangzhou'): 'ZCH Location',
        ('Location', 'UCMC'): 'UCMC Location',
        ('Location', 'Cincinnati'): 'UCMC Location',
        ('GestationCohort', '28-32 Weeks'): '28-32 Weeks Gestation',
        ('GestationCohort', '33-36 Weeks'): '33-36 Weeks Gestation',
        ('SampleCollectionWeek', 'Week.3'): 'Week 3 Sample',
        ('SampleCollectionWeek', 'Week 3'): 'Week 3 Sample',
        ('MaternalAntibiotics', 'None.Mat.Abx'): 'No Maternal Antibiotics',
        ('MaternalAntibiotics', 'No Mat Abx'): 'No Maternal Antibiotics',
        ('PostNatalAbxCohort', 'No.Infant.Abx'): 'No Infant Antibiotics',
        ('PostNatalAbxCohort', 'No Infant Abx'): 'No Infant Antibiotics',
        ('PostNatalAbxCohort', 'Low.Infant.Abx'): 'Low Infant Antibiotics',
        ('PostNatalAbxCohort', 'Low Infant Abx'): 'Low Infant Antibiotics',
        ('BSI_30D', '1'): 'BSI Positive',
        ('BSI30D', '1'): 'BSI Positive',
        ('NEC_30D', '1'): 'NEC Positive',
        ('NEC30D', '1'): 'NEC Positive',
        ('AnyMilk', '1'): "Mother's Milk",
        ('AnyMilk', 'Mother'): "Mother's Milk",
        ('AnyMilk', 'No.Milk'): "No Milk",
        ('PICC', '1'): 'PICC Present',
        ('PICC', 'PICC_UE'): 'PICC Upper Extremity',
        ('PICC', 'PICC_LE'): 'PICC Lower Extremity',
        ('PICC', 'PICC_Neck'): 'PICC Neck',
        ('PICC', 'axillary'): 'PICC Axillary',
        ('PICC', 'peripheral_UE'): 'PICC Peripheral UE',
        ('UVC', '1'): 'UVC Present',
        ('UVC', 'UVC'): 'UVC Present',
        ('Delivery', 'Vaginal'): 'Vaginal Delivery',
    }

    key = (category, value)
    if key in label_map:
        return label_map[key]

    # Default: return a cleaned version
    return f"{category}: {value}"


def select_best_transformation(df):
    """
    Automatically selects the best microbiome transformation based on data properties.
    """
    zero_fraction = (df == 0).sum().sum() / df.size
    total_reads_variation = df.sum(axis=1).std() / df.sum(axis=1).mean()

    if zero_fraction > 0.30:
        print("High zero fraction detected! Using CLR Transformation.")
        return mt.clr_transform(df)
    elif total_reads_variation > 0.5:
        print("Large sequencing depth variation detected! Applying Rarefaction.")
        return mt.rarefaction(df)
    else:
        print("Data appears normalized, applying Total Sum Scaling (TSS).")
        return mt.tss_transform(df)


def create_combined_figure(glmm_coefs, shap_importance, pc_name, variance_explained, output_path,
                           exclude_shap_features=None, exclude_glmm_features=None):
    """
    Create a combined figure with GLMM coefficients on left and SHAP importance on right.
    Matches the style of bsi_pedictors_viz.py (Klebsiella_oxytoca_comparison.pdf style).

    Parameters:
    - glmm_coefs: DataFrame with GLMM coefficients
    - shap_importance: DataFrame with SHAP importance values
    - pc_name: Name of the principal component (e.g., "PC1")
    - variance_explained: Proportion of variance explained by this PC
    - output_path: Path to save the figure
    - exclude_shap_features: List of feature names to exclude from SHAP plot (optional)
    - exclude_glmm_features: List of feature names to exclude from GLMM plot (optional)
    """

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # --- Left panel: GLMM Coefficients (Absolute Effect Size) ---
    # Filter to significant variables and exclude intercept/group variance
    glmm_plot = glmm_coefs[
        (glmm_coefs['Significant'] == True) &
        (~glmm_coefs['Variable'].str.contains('Intercept|Group Var', case=False, na=False))
    ].copy()

    # Exclude specified GLMM features if provided (supports partial matching)
    if exclude_glmm_features:
        for pattern in exclude_glmm_features:
            glmm_plot = glmm_plot[~glmm_plot['Variable'].str.contains(pattern, case=False, na=False)]

    if len(glmm_plot) == 0:
        # If no significant variables, show top 10 by absolute coefficient
        glmm_plot = glmm_coefs[
            ~glmm_coefs['Variable'].str.contains('Intercept|Group Var', case=False, na=False)
        ].copy()
        # Re-apply exclusion filters
        if exclude_glmm_features:
            for pattern in exclude_glmm_features:
                glmm_plot = glmm_plot[~glmm_plot['Variable'].str.contains(pattern, case=False, na=False)]
        glmm_plot['AbsCoef'] = glmm_plot['Coefficient'].abs()
        glmm_plot = glmm_plot.nlargest(10, 'AbsCoef')

    # Convert to absolute values for ranking without directionality
    glmm_plot['AbsCoefficient'] = glmm_plot['Coefficient'].abs()

    # Sort by absolute coefficient value (ascending so largest at top in horizontal bar)
    glmm_plot = glmm_plot.sort_values('AbsCoefficient', ascending=True)

    # Limit to top 10 for readability (matching bsi_pedictors_viz style)
    if len(glmm_plot) > 10:
        glmm_plot = glmm_plot.tail(10)

    # Create readable labels
    glmm_plot['Label'] = glmm_plot['Variable'].apply(create_comparison_label)

    y_pos = np.arange(len(glmm_plot))

    # Use single color (steelblue) since we're showing absolute effect sizes
    bars = ax1.barh(y_pos, glmm_plot['AbsCoefficient'].values,
                    color='steelblue', alpha=0.7, edgecolor='darkgray', linewidth=1)

    # Add coefficient values and significance stars
    for i, (coef, pval) in enumerate(zip(glmm_plot['AbsCoefficient'].values, glmm_plot['P-value'].values)):
        # Add value labels with larger font
        ax1.text(coef + 0.15, i, f'{coef:.2f}', va='center',
                ha='left', fontsize=12, fontweight='bold')

        # Add significance stars
        if pval < 0.001:
            sig_marker = '***'
        elif pval < 0.01:
            sig_marker = '**'
        elif pval < 0.05:
            sig_marker = '*'
        else:
            sig_marker = ''

        if sig_marker:
            # Position stars further away to avoid overlap with larger labels
            ax1.text(coef + 0.7, i, sig_marker,
                    va='center', ha='left',
                    fontsize=13, fontweight='bold')

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(glmm_plot['Label'].values, fontsize=13)
    ax1.set_xlabel('LMM Absolute Effect Size', fontsize=14)
    ax1.set_title(f'Linear Mixed Model\n{pc_name} ({variance_explained:.1%} variance)',
                  fontsize=16, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_xlim(0, glmm_plot['AbsCoefficient'].max() * 1.15)

    # --- Right panel: SHAP Feature Importance (Mean |SHAP Value| style) ---
    shap_plot = shap_importance.copy()

    # Exclude specified features if provided (supports partial matching)
    if exclude_shap_features:
        for pattern in exclude_shap_features:
            shap_plot = shap_plot[~shap_plot['Feature'].str.contains(pattern, case=False, na=False)]

    # Sort by importance (ascending for horizontal bar chart)
    shap_plot = shap_plot.sort_values('SHAP Importance', ascending=True)

    # Limit to top 10 for readability (matching bsi_pedictors_viz style)
    if len(shap_plot) > 10:
        shap_plot = shap_plot.tail(10)

    y_pos = np.arange(len(shap_plot))

    # Create horizontal bar chart with skyblue color
    ax2.barh(y_pos, shap_plot['SHAP Importance'].values,
             color='skyblue', alpha=0.7, edgecolor='darkblue', linewidth=1)

    # Add value labels
    for i, v in enumerate(shap_plot['SHAP Importance'].values):
        ax2.text(v + 0.002, i, f'{v:.3f}', va='center', fontsize=11)

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(shap_plot['Feature'].values, fontsize=13)
    ax2.set_xlabel('Mean |SHAP Value|', fontsize=14)
    ax2.set_title(f'Random Forest (SHAP Importance)\n{pc_name}',
                  fontsize=16, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_xlim(0, shap_plot['SHAP Importance'].max() * 1.15)

    # Overall title
    fig.suptitle(f'Top Predictors for Overall Microbiome Composition — GLMM and Random Forest',
                 fontsize=18, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved combined figure to {output_path}")


def main():
    print("=" * 60)
    print("GLMM + SHAP Analysis for Overall Microbiome Composition")
    print("=" * 60)

    # Define base directory
    base_dir = os.path.abspath(os.path.join(script_dir, '..', '..', '..'))

    # Load microbiome count data
    print("\nLoading microbiome data...")
    microbiome_path = os.path.join(base_dir, "data", "NICUSpeciesReduced.csv")
    microbiome_df = pd.read_csv(microbiome_path, index_col=0)

    print(f"Loaded microbiome data: {microbiome_df.shape}")

    # Transform microbiome data
    microbiome_transformed = select_best_transformation(microbiome_df)

    # Load clinical metadata
    print("\nLoading metadata...")
    metadata_path = os.path.join(base_dir, "metadata", "AllNICUSampleKey20250206.csv")
    metadata_df = pd.read_csv(metadata_path, index_col=0)

    print(f"Loaded metadata: {metadata_df.shape}")

    # Replace location names for consistency
    metadata_df['Location'] = metadata_df['Location'].replace({'Cincinnati': 'UCMC', 'Hangzhou': 'ZCH'})

    # Define categorical features
    categorical_features = ["SampleType", "Location", "GestationCohort", "SampleCollectionWeek",
                            "MaternalAntibiotics", "PostNatalAbxCohort", "BSI_30D", "NEC_30D",
                            "AnyMilk", "PICC", "UVC", "Delivery"]

    subject_id_col = "Subject"

    # Get common samples
    common_samples = list(set(microbiome_transformed.index) & set(metadata_df.index))
    microbiome_transformed = microbiome_transformed.loc[common_samples]
    metadata_df = metadata_df.loc[common_samples]

    print(f"\nUsing {len(common_samples)} samples with both microbiome and metadata")

    # Perform PCA on microbiome data
    print("\nPerforming PCA on microbiome composition...")
    n_components = 5
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(microbiome_transformed)

    # Create DataFrame with PC scores
    pc_df = pd.DataFrame(
        pca_result,
        index=microbiome_transformed.index,
        columns=[f"PC{i+1}" for i in range(n_components)]
    )

    # Report variance explained
    print("\nVariance explained by each PC:")
    for i, var in enumerate(pca.explained_variance_ratio_):
        print(f"  PC{i+1}: {var:.2%}")

    # Merge PC scores with metadata
    combined_df = pd.merge(pc_df, metadata_df, left_index=True, right_index=True)

    # Store all results
    all_glmm_results = []
    all_shap_results = []

    # Analyze each principal component
    for pc_idx, pc_name in enumerate([f"PC{i+1}" for i in range(n_components)]):
        print(f"\n{'='*60}")
        print(f"Analyzing {pc_name} ({pca.explained_variance_ratio_[pc_idx]:.2%} variance)")
        print('='*60)

        variance_explained = pca.explained_variance_ratio_[pc_idx]

        # --- GLMM Analysis ---
        print(f"\nFitting GLMM for {pc_name}...")

        # Prepare data for GLMM
        df_glmm = combined_df[[pc_name] + [col for col in categorical_features if col in combined_df.columns] + [subject_id_col]].copy()
        df_glmm = df_glmm.dropna(subset=[pc_name, subject_id_col])

        # Encode categorical variables
        formula_parts = []
        for col in categorical_features:
            if col in df_glmm.columns:
                df_glmm[col] = df_glmm[col].fillna('missing')
                df_glmm[col] = df_glmm[col].astype("category")
                formula_parts.append(f"C({col})")

        fixed_effects = " + ".join(formula_parts)
        formula = f"{pc_name} ~ {fixed_effects}"

        try:
            # Fit mixed effects model
            model = mixedlm(formula, df_glmm, groups=df_glmm[subject_id_col])
            result = model.fit(method='bfgs', maxiter=1000)

            # Extract coefficients
            glmm_coefs = pd.DataFrame({
                'Variable': result.params.index,
                'Coefficient': result.params.values,
                'Std_Error': result.bse.values,
                'T_stat': result.tvalues.values,
                'P-value': result.pvalues.values,
                'Significant': result.pvalues < 0.05
            })

            # Add confidence intervals
            conf_int = result.conf_int()
            glmm_coefs['CI_Lower'] = conf_int[0].values
            glmm_coefs['CI_Upper'] = conf_int[1].values

            glmm_coefs['PC'] = pc_name
            all_glmm_results.append(glmm_coefs)

            # Save GLMM coefficients
            glmm_coefs.to_csv(os.path.join(OUTPUT_DIR, f"glmm_coefficients_{pc_name}_overall_composition.csv"), index=False)

            print(f"GLMM fitted successfully. {(glmm_coefs['Significant']).sum()} significant variables.")

        except Exception as e:
            print(f"Error fitting GLMM: {e}")
            glmm_coefs = None

        # --- SHAP Analysis ---
        print(f"\nRunning SHAP analysis for {pc_name}...")

        # Prepare data for Random Forest
        # Use only the same categorical features as GLMM (matching individual organism analysis)
        y = combined_df[pc_name]
        available_features = [col for col in categorical_features if col in metadata_df.columns]
        X = metadata_df.loc[combined_df.index, available_features].copy()

        # Fill missing values (keep as strings for proper dummy variable names)
        for col in available_features:
            X[col] = X[col].fillna('missing')

        # Create dummy variables directly from categorical strings (preserves category names)
        X_encoded = pd.get_dummies(X, drop_first=True)
        y = y.loc[X_encoded.index]

        # Train-test split
        X_train, X_test, y_train, y_test = train_test_split(X_encoded, y, test_size=0.2, random_state=42)

        # Train Random Forest
        rf = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
        rf.fit(X_train, y_train)

        # Compute SHAP values
        explainer = shap.TreeExplainer(rf)
        shap_values = explainer(X_test)

        # Calculate mean absolute SHAP values
        shap_values_mean = np.abs(shap_values.values).mean(axis=0)

        # Create readable feature names
        feature_names_readable = [create_comparison_label(col) for col in X_encoded.columns]

        shap_importance = pd.DataFrame({
            'Feature': feature_names_readable,
            'Feature_Original': X_encoded.columns,
            'SHAP Importance': shap_values_mean
        }).sort_values('SHAP Importance', ascending=False)

        shap_importance['PC'] = pc_name
        all_shap_results.append(shap_importance)

        # Save SHAP importance
        shap_importance.to_csv(os.path.join(OUTPUT_DIR, f"shap_importance_{pc_name}_overall_composition.csv"), index=False)

        print(f"SHAP analysis complete. Top feature: {shap_importance.iloc[0]['Feature']}")

        # --- Create SHAP Summary Plot ---
        try:
            plt.figure(figsize=(12, 8))
            shap.summary_plot(shap_values, X_test, feature_names=feature_names_readable, show=False)
            plt.title(f"SHAP Summary: {pc_name} Overall Composition")
            plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, f"shap_summary_{pc_name}_overall_composition.pdf"), bbox_inches="tight")
            plt.close()
        except Exception as e:
            print(f"Could not create SHAP summary plot: {e}")

        # --- Create Combined Figure ---
        if glmm_coefs is not None:
            # Exclude "Low Infant Antibiotics" - keep only "No Infant Antibiotics"
            exclude_abx_patterns = ['Low.Infant.Abx', 'Low Infant Abx', 'Low Infant Antibiotics']
            create_combined_figure(
                glmm_coefs=glmm_coefs,
                shap_importance=shap_importance,
                pc_name=pc_name,
                variance_explained=variance_explained,
                output_path=os.path.join(OUTPUT_DIR, f"combined_glmm_shap_{pc_name}_overall_composition.pdf"),
                exclude_shap_features=exclude_abx_patterns,
                exclude_glmm_features=exclude_abx_patterns
            )

    # Save combined results
    if all_glmm_results:
        combined_glmm = pd.concat(all_glmm_results, ignore_index=True)
        combined_glmm.to_csv(os.path.join(OUTPUT_DIR, "glmm_all_PCs_overall_composition.csv"), index=False)

    if all_shap_results:
        combined_shap = pd.concat(all_shap_results, ignore_index=True)
        combined_shap.to_csv(os.path.join(OUTPUT_DIR, "shap_all_PCs_overall_composition.csv"), index=False)

    # --- Create Summary Heatmap ---
    print("\nCreating summary heatmap...")
    create_summary_heatmap(all_glmm_results, all_shap_results, pca.explained_variance_ratio_)

    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)
    print(f"\nOutput files saved to: {OUTPUT_DIR}")
    print("\nGenerated files:")
    print("  - glmm_coefficients_PC*_overall_composition.csv")
    print("  - shap_importance_PC*_overall_composition.csv")
    print("  - combined_glmm_shap_PC*_overall_composition.pdf")
    print("  - shap_summary_PC*_overall_composition.pdf")
    print("  - glmm_all_PCs_overall_composition.csv")
    print("  - shap_all_PCs_overall_composition.csv")
    print("  - glmm_shap_summary_heatmap_overall_composition.pdf")


def create_summary_heatmap(glmm_results, shap_results, variance_explained):
    """Create a summary heatmap showing significant associations across all PCs (using absolute effect sizes)"""
    if not glmm_results:
        return

    # Combine all GLMM results
    combined_glmm = pd.concat(glmm_results, ignore_index=True)

    # Exclude "Low Infant Antibiotics" patterns
    exclude_patterns = ['Low.Infant.Abx', 'Low Infant Abx', 'Low Infant Antibiotics']
    for pattern in exclude_patterns:
        combined_glmm = combined_glmm[~combined_glmm['Variable'].str.contains(pattern, case=False, na=False)]

    # Filter to significant variables (excluding intercept and group variance)
    sig_vars = combined_glmm[
        (combined_glmm['Significant'] == True) &
        (~combined_glmm['Variable'].str.contains('Intercept|Group Var', case=False, na=False))
    ]['Variable'].unique()

    if len(sig_vars) == 0:
        print("No significant variables found for heatmap")
        return

    # Create coefficient matrix (using absolute values)
    pcs = combined_glmm['PC'].unique()
    coef_matrix = pd.DataFrame(index=sig_vars, columns=pcs)
    pval_matrix = pd.DataFrame(index=sig_vars, columns=pcs)

    for _, row in combined_glmm.iterrows():
        if row['Variable'] in sig_vars:
            # Use absolute value of coefficient
            coef_matrix.loc[row['Variable'], row['PC']] = abs(row['Coefficient'])
            pval_matrix.loc[row['Variable'], row['PC']] = row['P-value']

    coef_matrix = coef_matrix.astype(float)

    # Create readable labels
    readable_labels = [create_comparison_label(v) for v in coef_matrix.index]

    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(8, len(sig_vars) * 0.4)))

    # Create heatmap with single-color scale (Blues) since we're using absolute values
    sns.heatmap(coef_matrix.values,
                xticklabels=[f"{pc}\n({var:.1%})" for pc, var in zip(pcs, variance_explained)],
                yticklabels=readable_labels,
                cmap='Blues',
                annot=True, fmt='.2f',
                cbar_kws={'label': 'GLMM Absolute Effect Size'},
                ax=ax)

    # Add significance markers
    for i in range(len(sig_vars)):
        for j, pc in enumerate(pcs):
            pval = pval_matrix.iloc[i, j]
            if pd.notna(pval) and pval < 0.05:
                marker = '***' if pval < 0.001 else ('**' if pval < 0.01 else '*')
                ax.text(j + 0.5, i + 0.75, marker, ha='center', va='center',
                       fontsize=8, fontweight='bold', color='white')

    ax.set_title('GLMM Absolute Effect Sizes for Overall Microbiome Composition\n(Significant Variables Across PCs)',
                fontsize=14, fontweight='bold')
    ax.set_xlabel('Principal Components (% Variance Explained)', fontsize=12)
    ax.set_ylabel('Clinical Variables', fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "glmm_shap_summary_heatmap_overall_composition.pdf"),
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, "glmm_shap_summary_heatmap_overall_composition.png"),
                dpi=300, bbox_inches='tight')
    plt.close()

    print("Saved summary heatmap")


if __name__ == "__main__":
    main()
