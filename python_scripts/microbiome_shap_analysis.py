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
from sklearn.model_selection import GroupKFold
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
subject_id_col = "Subject"  # Ensure this is the correct column name

# Check for missing values in metadata
print(metadata_df.isnull().sum())

categorical_features = ["SampleType", "Location", "GestationCohort", "SampleCollectionWeek", 
                        "MaternalAntibiotics", "PostNatalAbxCohort", "BSI_30D", "NEC_30D", "AnyMilk", "PICC", "UVC"]

# Keep only required columns
metadata_df = metadata_df[categorical_features + [subject_id_col]]

# Merge metadata with transformed microbiome data
data = microbiome_transformed.merge(metadata_df, left_index=True, right_index=True, how="inner")
print(f"Data Shape After Merge: {data.shape}")

# List of organisms to analyze
key_organisms = ["Klebsiella.pneumoniae", "Staphylococcus.aureus", "Escherichia.coli", "Klebsiella.oxytoca",
                 "Staphylococcus.epidermidis", "Streptococcus.pyogenes"]

# Store results
shap_importance_all_folds = []
shap_summary_all_folds = []

for target_microbe in key_organisms:
    print(f"\n🦠 Analyzing: {target_microbe}")

    if target_microbe not in data.columns:
        print(f"⚠️ {target_microbe} not found in data. Skipping.")
        continue

    y = data[target_microbe]
    X = metadata_df.loc[data.index].copy()
    groups = X[subject_id_col]

    # Encode categorical variables
    for col in categorical_features:
        if col in X.columns:  # Check if column exists
            # Handle potential NaN values
            X[col] = X[col].fillna('missing')
            X[col] = LabelEncoder().fit_transform(X[col])
    
    # Create dummy variables
    X_encoded = pd.get_dummies(X.drop(columns=[subject_id_col], errors='ignore'), drop_first=True)

    # Ensure X_encoded is not empty
    if X_encoded.empty:
        print(f"⚠️ No features available for {target_microbe} after preprocessing. Skipping.")
        continue

    y = y.loc[X_encoded.index]
    groups = groups.loc[X_encoded.index]

    # Check if we have enough data points
    if len(y) < 10:  # Arbitrary threshold
        print(f"⚠️ Not enough data points for {target_microbe}. Skipping.")
        continue

    gkf = GroupKFold(n_splits=5)
    fold_importances = []
    fold_shap_values = []
    
    # Collect RF models for later use
    rf_models = []

    try:
        for fold_idx, (train_idx, test_idx) in enumerate(gkf.split(X_encoded, y, groups)):
            print(f"📂 Fold {fold_idx+1}")

            X_train, X_test = X_encoded.iloc[train_idx], X_encoded.iloc[test_idx]
            y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

            rf = RandomForestRegressor(n_estimators=100, random_state=42)
            rf.fit(X_train, y_train)
            rf_models.append(rf)

            # Use TreeExplainer instead of Explainer for RandomForestRegressor
            explainer = shap.TreeExplainer(rf)
            shap_values = explainer(X_test)

            # Compute SHAP values for this fold
            shap_values_mean = np.abs(shap_values.values).mean(axis=0)
            fold_importances.append(shap_values_mean)
            fold_shap_values.append(shap_values.values)

        # Check if we have any SHAP values to work with
        if not fold_shap_values:
            print(f"⚠️ No SHAP values calculated for {target_microbe}. Skipping.")
            continue

        # Create a representative sample of X for plotting
        X_sample = X_encoded.sample(min(100, len(X_encoded)), random_state=42)

        # Create SHAP Explanation Object with properly sized data
        shap_values_concat = np.concatenate(fold_shap_values, axis=0)
        
        # Create data array of the correct size to match shap_values_concat
        total_samples_needed = len(shap_values_concat)
        
        # Generate a data matrix of the right size by sampling or repeating X_encoded
        if len(X_encoded) >= total_samples_needed:
            data_for_exp = X_encoded.sample(total_samples_needed, random_state=42, replace=False).values
        else:
            # If we need more samples than available, use replacement
            data_for_exp = X_encoded.sample(total_samples_needed, random_state=42, replace=True).values
            
        shap_exp = shap.Explanation(
            values=shap_values_concat,
            base_values=np.zeros(len(shap_values_concat)),  # Default base values
            data=data_for_exp,
            feature_names=X_encoded.columns.tolist()
        )

        # Save SHAP Summary Plot
        try:
            plt.figure(figsize=(12, 8))
            # Use an alternative SHAP visualization method for more stability
            shap.plots.beeswarm(shap_exp, show=False, max_display=20)
            plt.title(f"Aggregated SHAP Summary for {target_microbe}")
            plt.tight_layout()
            plt.savefig(f"shap_summary_{target_microbe}_aggregated.pdf", bbox_inches="tight")
            plt.close()
        except Exception as e:
            print(f"⚠️ Could not create SHAP summary plot: {e}")
            # Alternative: create a simpler plot just with the feature importances
            plt.figure(figsize=(12, 8))
            top_features = X_encoded.columns[np.argsort(-np.mean(fold_importances, axis=0))[:min(20, len(X_encoded.columns))]]
            importances = np.mean(fold_importances, axis=0)[np.argsort(-np.mean(fold_importances, axis=0))[:min(20, len(X_encoded.columns))]]
            plt.barh(top_features, importances)
            plt.title(f"Feature Importance for {target_microbe}")
            plt.tight_layout()
            plt.savefig(f"feature_importance_{target_microbe}_fallback.pdf", bbox_inches="tight")
            plt.close()

        # Aggregate feature importance across folds
        avg_importance = np.mean(fold_importances, axis=0)
        
        # Check for NaN values in importance
        if np.isnan(avg_importance).any():
            print(f"⚠️ NaN values in importance for {target_microbe}. Cleaning up.")
            avg_importance = np.nan_to_num(avg_importance)
        
        shap_importance_df = pd.DataFrame({
            "Feature": X_encoded.columns,
            "SHAP Importance": avg_importance
        }).sort_values(by="SHAP Importance", ascending=False)

        shap_importance_all_folds.append(shap_importance_df.assign(Microbe=target_microbe))

        # Plot the aggregated feature importance
        plt.figure(figsize=(10, 6))
        top_n = min(15, len(shap_importance_df))
        sns.barplot(x="SHAP Importance", y="Feature", data=shap_importance_df.head(top_n), palette="viridis")
        plt.xlabel("SHAP Importance Score")
        plt.ylabel("Features")
        plt.title(f"Overall Feature Importance for {target_microbe} (Aggregated Across Folds)")
        plt.tight_layout()
        plt.savefig(f"shap_feature_importance_{target_microbe}_aggregated.pdf", bbox_inches="tight")
        plt.close()
    
    except Exception as e:
        print(f"❌ Error processing {target_microbe}: {str(e)}")
        continue

# Save all feature importance results as a CSV
if shap_importance_all_folds:
    combined_shap_importance = pd.concat(shap_importance_all_folds)
    combined_shap_importance.to_csv("shap_feature_importance_all_microbes_aggregated.csv", index=False)
else:
    print("⚠️ No SHAP results to save.")

# Mixed Effects Model for Subject-Level Random Effects
mixedlm_results = []

for target_microbe in key_organisms:
    print(f"\n📊 Fitting Mixed Effects Model for: {target_microbe}")

    if target_microbe not in data.columns:
        print(f"⚠️ {target_microbe} not found in data. Skipping.")
        continue

    df = data[[target_microbe] + [col for col in categorical_features if col in data.columns] + [subject_id_col]].copy()
    df = df.dropna(subset=[target_microbe, subject_id_col])  # Ensure target and grouping variable aren't missing
    
    if df.empty:
        print(f"⚠️ No data available for {target_microbe} after dropping NAs. Skipping.")
        continue
        
    df["MicrobeAbundance"] = df[target_microbe]

    # Encode categorical variables for formula-based modeling
    formula_parts = []
    for col in categorical_features:
        if col in df.columns:
            # Fill NAs and convert to category
            df[col] = df[col].fillna('missing')
            df[col] = df[col].astype("category")
            formula_parts.append(f"C({col})")
    
    if not formula_parts:
        print(f"⚠️ No categorical features available for {target_microbe}. Skipping.")
        continue
        
    fixed_effects = " + ".join(formula_parts)
    formula = f"MicrobeAbundance ~ {fixed_effects}"

    try:
        # Ensure the grouping variable has enough groups
        group_counts = df[subject_id_col].value_counts()
        if len(group_counts) < 5:  # Arbitrary threshold
            print(f"⚠️ Not enough subjects ({len(group_counts)}) for mixed model. Fitting without random effects.")
            formula = f"MicrobeAbundance ~ {fixed_effects}"
            model = sm.formula.ols(formula, data=df).fit()
            result_summary = model.summary()
        else:
            model = mixedlm(formula, df, groups=df[subject_id_col])
            result = model.fit()
            result_summary = result.summary()
        
        print(result_summary)

        # Extract coefficient table in a more robust way
        try:
            coef_df = pd.read_html(result_summary.tables[1].as_html(), header=0, index_col=0)[0]
            coef_df = coef_df.reset_index()
            coef_df["Microbe"] = target_microbe
            mixedlm_results.append(coef_df)
        except Exception as e:
            print(f"⚠️ Could not extract coefficient table: {e}")
            # Create a simplified DataFrame from the parameters
            coef_df = pd.DataFrame({
                "Variable": model.params.index,
                "Coef.": model.params.values,
                "Microbe": target_microbe
            })
            mixedlm_results.append(coef_df)

        # Save individual model results
        with open(f"mixedlm_{target_microbe}.txt", "w") as f:
            f.write(str(result_summary))

    except Exception as e:
        print(f"❌ Failed to fit model for {target_microbe}: {str(e)}")

# Save combined results if available
if mixedlm_results:
    try:
        combined_results = pd.concat(mixedlm_results, ignore_index=True)
        combined_results.to_csv("mixedlm_summary_all_microbes.csv", index=False)
    except Exception as e:
        print(f"❌ Failed to save combined results: {str(e)}")
        # Try to save individual dataframes
        for i, df in enumerate(mixedlm_results):
            df.to_csv(f"mixedlm_results_part_{i}.csv", index=False)
else:
    print("⚠️ No mixed model results to save.")