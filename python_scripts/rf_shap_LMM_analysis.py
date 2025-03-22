# microbiome_shap_analysis.py
# Random Forest + SHAP + Mixed Effects Model with subject-level random effects

import pandas as pd
import numpy as np
import shap
import microbiome_transform as mt  # Import the transformation module
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
import statsmodels.api as sm
from statsmodels.formula.api import mixedlm

def select_best_transformation(df):
    """
    Automatically selects the best microbiome transformation based on data properties.
    
    Parameters:
    - df: Pandas DataFrame (raw microbiome count data)

    Returns:
    - Transformed DataFrame
    """
    zero_fraction = (df == 0).sum().sum() / df.size  # % of zero values
    total_reads_variation = df.sum(axis=1).std() / df.sum(axis=1).mean()  # Variability in sequencing depth

    if zero_fraction > 0.30:  
        print("⚠️ High zero fraction detected! Using CLR Transformation.")
        return mt.clr_transform(df)
    elif total_reads_variation > 0.5:  
        print("⚠️ Large sequencing depth variation detected! Applying Rarefaction.")
        return mt.rarefaction(df)
    else:
        print("✅ Data appears normalized, applying Total Sum Scaling (TSS).")
        return mt.tss_transform(df)

# Load microbiome count data
microbiome_df = pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)

# Auto-select and apply transformation
microbiome_transformed = select_best_transformation(microbiome_df)
microbiome_transformed.to_csv("microbiome_abundance_transformed.csv")

# Load clinical metadata
metadata_df = pd.read_csv("../metadata/AllNICUSampleKey20250206.csv", index_col=0)
subject_id_col = "Subject"  

# Check for missing values in metadata
print(metadata_df.isnull().sum())

categorical_features = ["SampleType", "Location", "GestationCohort", "SampleCollectionWeek", 
                        "MaternalAntibiotics", "PostNatalAbxCohort", "BSI_30D", "NEC_30D", "AnyMilk", "PICC", "UVC"]

# Keep only required columns
metadata_df = metadata_df[categorical_features + [subject_id_col]]

# Merge metadata with transformed microbiome data
data = microbiome_transformed.merge(metadata_df, left_index=True, right_index=True, how="inner")
print(f"Data Shape After Merge: {data.shape}")

# List of microbes to analyze
key_organisms = ["Klebsiella.pneumoniae", "Staphylococcus.aureus", "Escherichia.coli", "Klebsiella.oxytoca",
                 "Staphylococcus.epidermidis", "Streptococcus.pyogenes", "Staphylococcs.capitus", 
                 "Enterococcus.faecium", "Enterococcus.faecalis", "Serratia.marcescens", "Listeria monocytogenes"]

# Store results for all microbes
all_shap_importance = []

for target_microbe in key_organisms:
    print(f"\n🦠 Analyzing: {target_microbe}")

    if target_microbe not in data.columns:
        print(f"⚠️ {target_microbe} not found in data. Skipping.")
        continue

    try:
        y = data[target_microbe]

        # Subset X to matching samples
        X = metadata_df.loc[data.index].copy()

        # Label encode categorical variables
        for col in categorical_features:
            if col in X.columns:  # Check if column exists
                # Handle potential NaN values
                X[col] = X[col].fillna('missing')
                try:
                    X[col] = LabelEncoder().fit_transform(X[col])
                except:
                    print(f"⚠️ Error encoding {col}, dropping column")
                    X = X.drop(columns=[col])

        # Create dummy variables without the Subject column
        X_encoded = pd.get_dummies(X.drop(columns=[subject_id_col]), drop_first=True)

        # Ensure X_encoded is not empty
        if X_encoded.empty:
            print(f"⚠️ No features available for {target_microbe} after preprocessing. Skipping.")
            continue

        # Match y to encoded X
        y = y.loc[X_encoded.index]

        # Check if we have enough data points
        if len(y) < 10:  # Arbitrary threshold
            print(f"⚠️ Not enough data points for {target_microbe}. Skipping.")
            continue

        # Train-test split
        X_train, X_test, y_train, y_test = train_test_split(X_encoded, y, test_size=0.2, random_state=42)

        # Train Random Forest Regressor
        rf = RandomForestRegressor(n_estimators=100, random_state=42)
        rf.fit(X_train, y_train)

        # Compute SHAP values using TreeExplainer (more appropriate for Random Forest)
        explainer = shap.TreeExplainer(rf)
        shap_values = explainer(X_test)

        # SHAP Summary Plot
        try:
            plt.figure(figsize=(12, 8))
            shap.summary_plot(shap_values, X_test, feature_names=X_encoded.columns.tolist(), show=False)
            plt.title(f"SHAP Summary: {target_microbe}")
            plt.tight_layout()
            plt.savefig(f"../results/shap_summary_{target_microbe}.pdf", bbox_inches="tight")
            plt.close()
        except Exception as e:
            print(f"⚠️ Error creating SHAP summary plot: {e}")
            try:
                # Try alternative visualization
                plt.figure(figsize=(12, 8))
                shap.plots.beeswarm(shap_values, max_display=20, show=False)
                plt.title(f"SHAP Summary: {target_microbe}")
                plt.tight_layout()
                plt.savefig(f"../results/shap_summary_{target_microbe}_alternative.pdf", bbox_inches="tight")
                plt.close()
            except Exception as e2:
                print(f"⚠️ Error creating alternative SHAP plot: {e2}")

        # SHAP Feature Importance - New Line+Dot Plot (similar to PyCaret style)
        shap_values_mean = np.abs(shap_values.values).mean(axis=0)
        shap_importance = pd.DataFrame({
            "Feature": X_encoded.columns,
            "SHAP Importance": shap_values_mean
        }).sort_values(by="SHAP Importance", ascending=True)
        
        # Save this microbe's SHAP importance to the combined results
        shap_importance["Microbe"] = target_microbe
        all_shap_importance.append(shap_importance)

        # Create line+dot plot (PyCaret style) for feature importance
        plt.figure(figsize=(12, 8))
        top_n = min(15, len(shap_importance))
        top_features = shap_importance.head(top_n)
        
        # Create horizontal line + dot plot
        plt.hlines(y=range(top_n), 
                 xmin=0, 
                 xmax=top_features["SHAP Importance"].values,
                 color="skyblue", 
                 alpha=0.7, 
                 linewidth=2)
        
        plt.plot(top_features["SHAP Importance"].values, 
                range(top_n), 
                "o", 
                markersize=10, 
                color="blue", 
                alpha=0.8)
        
        # Add feature names
        plt.yticks(range(top_n), top_features["Feature"].values)
        plt.xlabel("SHAP Importance Score")
        plt.title(f"Clinical Variables Associated with {target_microbe} Abundance")
        
        # Add values next to dots
        for i, importance in enumerate(top_features["SHAP Importance"].values):
            plt.text(importance + 0.001, i, f"{importance:.4f}", va='center')
            
        plt.tight_layout()
        plt.savefig(f"../results/shap_feature_importance_{target_microbe}.pdf", bbox_inches="tight")
        plt.close()

    except Exception as e:
        print(f"❌ Error processing {target_microbe}: {str(e)}")
        continue

# Save all feature importance results as a CSV
if all_shap_importance:
    combined_shap_importance = pd.concat(all_shap_importance)
    combined_shap_importance.to_csv("shap_feature_importance_all_microbes.csv", index=False)
else:
    print("⚠️ No SHAP results to save.")

# Mixed Effects Model for Subject-Level Random Effects
mixedlm_results = []

for target_microbe in key_organisms:
    print(f"\n📊 Fitting Mixed Effects Model for: {target_microbe}")

    if target_microbe not in data.columns:
        print(f"⚠️ {target_microbe} not found in data. Skipping.")
        continue

    try:
        # Select only the columns we need
        df = data[[target_microbe] + [col for col in categorical_features if col in data.columns] + [subject_id_col]].copy()
        df = df.dropna(subset=[target_microbe, subject_id_col])  # Ensure target and grouping variable aren't missing
        
        if df.empty:
            print(f"⚠️ No data available for {target_microbe} after dropping NAs. Skipping.")
            continue
            
        df["MicrobeAbundance"] = df[target_microbe]
        
        # Check for extremely low variance in target variable
        if df["MicrobeAbundance"].var() < 1e-6:
            print(f"⚠️ Near-zero variance in target for {target_microbe}. Skipping.")
            continue

        # Encode categorical variables for formula-based modeling
        formula_parts = []
        for col in categorical_features:
            if col in df.columns:
                # Count unique values to avoid categories with too few samples
                value_counts = df[col].value_counts()
                
                # Filter out categories with very few samples (can cause convergence issues)
                rare_categories = value_counts[value_counts < 3].index.tolist()
                if rare_categories:
                    print(f"⚠️ Removing rare categories in {col}: {rare_categories}")
                    df.loc[df[col].isin(rare_categories), col] = "Other"
                
                # Fill NAs and convert to category
                df[col] = df[col].fillna('missing')
                df[col] = df[col].astype("category")
                formula_parts.append(f"C({col})")
        
        if not formula_parts:
            print(f"⚠️ No categorical features available for {target_microbe}. Skipping.")
            continue
            
        # First try regular OLS regression
        # This is more stable and will help identify problematic variables
        fixed_effects = " + ".join(formula_parts)
        ols_formula = f"MicrobeAbundance ~ {fixed_effects}"
        
        print(f"First fitting OLS model to identify potential issues")
        ols_model = sm.formula.ols(ols_formula, data=df).fit()
        
        # Examine OLS results for problematic variables (perfect collinearity, etc.)
        problem_vars = []
        for i, pval in enumerate(ols_model.pvalues):
            if np.isnan(pval):
                problem_vars.append(ols_model.model.exog_names[i])
        
        if problem_vars:
            print(f"⚠️ Found problematic variables in OLS model: {problem_vars}")
            # Remove problematic terms from formula
            for var in problem_vars:
                formula_parts = [part for part in formula_parts if var not in part]
            
            if not formula_parts:
                print(f"⚠️ No valid categorical features left for {target_microbe}. Skipping.")
                continue
                
            fixed_effects = " + ".join(formula_parts)
        
        # Try mixed effects model with simplified formula
        mixed_formula = f"MicrobeAbundance ~ {fixed_effects}"
        print(f"Fitting mixed effects model with formula: {mixed_formula}")

        # Ensure the grouping variable has enough groups and sufficient data
        group_counts = df[subject_id_col].value_counts()
        min_group_size = group_counts.min()
        num_groups = len(group_counts)
        
        print(f"Groups: {num_groups}, Min group size: {min_group_size}")
        
        if num_groups < 5 or min_group_size < 2:
            print(f"⚠️ Not enough subjects or samples per subject for mixed model. Fitting OLS instead.")
            model = sm.formula.ols(mixed_formula, data=df).fit()
            model_type = "OLS"
            result_summary = model.summary()
        else:
            # Try with different optimizers if default fails
            try:
                # First try with default settings but more iterations
                model = mixedlm(mixed_formula, df, groups=df[subject_id_col])
                result = model.fit(maxiter=1000)
                model_type = "MixedLM"
                result_summary = result.summary()
            except:
                try:
                    # Try with L-BFGS-B optimizer 
                    print("⚠️ Default optimizer failed, trying with L-BFGS-B")
                    model = mixedlm(mixed_formula, df, groups=df[subject_id_col])
                    result = model.fit(method="lbfgs", maxiter=2000)
                    model_type = "MixedLM (LBFGS)"
                    result_summary = result.summary()
                except:
                    # Fall back to OLS if mixed model fails
                    print("⚠️ Mixed model failed to converge, falling back to OLS")
                    model = sm.formula.ols(mixed_formula, data=df).fit()
                    model_type = "OLS (fallback)"
                    result_summary = model.summary()
        
        print(result_summary)

        # Extract coefficients manually since read_html might fail
        if hasattr(model, 'params'):
            # Create DataFrame directly from model parameters
            coef_df = pd.DataFrame({
                "Variable": model.params.index.tolist(),
                "Coefficient": model.params.values,
                "P-value": model.pvalues.values if hasattr(model, 'pvalues') else [None] * len(model.params),
                "Model_Type": model_type,
                "Microbe": target_microbe
            })
            mixedlm_results.append(coef_df)
        else:
            print(f"⚠️ Could not extract coefficients: model has no params attribute")

        # Save individual model results
        with open(f"../results/model_results_{target_microbe}.txt", "w") as f:
            f.write(str(result_summary))
            f.write(f"\nModel type: {model_type}")

    except Exception as e:
        print(f"❌ Failed to fit model for {target_microbe}: {str(e)}")
        import traceback
        traceback.print_exc()

# Save combined results if available
if mixedlm_results:
    try:
        combined_results = pd.concat(mixedlm_results, ignore_index=True)
        combined_results.to_csv("../results/model_summary_all_microbes.csv", index=False)
        print(f"✅ Saved results for {len(mixedlm_results)} microbes")
    except Exception as e:
        print(f"❌ Failed to save combined results: {str(e)}")
        # Try to save individual dataframes
        for i, df in enumerate(mixedlm_results):
            df.to_csv(f"../results/model_results_part_{i}.csv", index=False)
else:
    print("⚠️ No mixed model results to save.")