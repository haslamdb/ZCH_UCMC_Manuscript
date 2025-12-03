import pandas as pd
import numpy as np
import shap
import microbiome_transform as mt
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
from scipy import stats
import warnings
from tqdm import tqdm
import os
import re
from statsmodels.genmod.families import Poisson, NegativeBinomial, Binomial
from statsmodels.genmod.families.links import log, logit
import pickle

# Silence common convergence warnings
warnings.filterwarnings("ignore", category=ConvergenceWarning)

def select_best_transformation(df, zero_proportion_threshold=0.25):
    """
    Automatically selects the best microbiome transformation based on data properties.
    For highly zero-inflated data (>25% zeros by default), avoid CLR transformation.
    
    Parameters:
    - df: Pandas DataFrame (raw microbiome count data)
    - zero_proportion_threshold: Threshold to determine when zero-inflation is severe

    Returns:
    - Transformed DataFrame and transformation method used
    """
    zero_fraction = (df == 0).sum().sum() / df.size  # % of zero values
    total_reads_variation = df.sum(axis=1).std() / df.sum(axis=1).mean()  # Sequencing depth variation
    
    print(f"Data diagnostics: {zero_fraction:.2%} zeros, {total_reads_variation:.2f} sequencing depth CV")
    
    # For zero-inflated data, avoid CLR which doesn't handle zeros well without pseudocounts
    if zero_fraction > zero_proportion_threshold:
        print("⚠️ High zero fraction detected! Using simple log transform with pseudocount.")
        return mt.log_transform(df, pseudocount=0.5), "log_transform"
    elif total_reads_variation > 0.5:
        print("⚠️ Large sequencing depth variation detected! Applying TSS normalization.")
        return mt.tss_transform(df), "tss_transform"
    else:
        print("✅ Data has moderate zeros, applying CLR transformation.")
        return mt.clr_transform(df, pseudocount=0.5), "clr_transform"

def fit_zero_inflated_glmm(formula, group_var, df, key_microbe, verbose=True):
    """
    Fit a zero-inflated negative binomial GLMM using a two-part model approach:
    1) Zero model: Logistic mixed model for presence/absence
    2) Count model: Negative binomial mixed model for non-zero counts
    
    Parameters:
    - formula: Formula string specifying fixed effects
    - group_var: Variable name for random effects grouping
    - df: DataFrame containing data
    - key_microbe: Name of the target microbe/response variable
    - verbose: Whether to print detailed output
    
    Returns:
    - Dictionary containing model results
    """
    results = {}
    
    # Create binary presence/absence indicator for zero model
    df['presence'] = (df[key_microbe] > 0).astype(int)
    
    # Create non-zero subset for count model
    nonzero_df = df[df[key_microbe] > 0].copy()
    
    # Check if we have enough non-zero samples for modeling
    if len(nonzero_df) < 10:
        print(f"⚠️ Too few non-zero observations ({len(nonzero_df)}) for count model.")
        return None
        
    # Calculate proportion of zeros for diagnostics
    zero_proportion = 1 - (len(nonzero_df) / len(df))
    print(f"Zero proportion: {zero_proportion:.2%} ({len(df) - len(nonzero_df)}/{len(df)})")
    
    try:
        # 1. Zero model (presence/absence)
        zero_formula = formula.replace(key_microbe, "presence")
        try:
            # Try Bayesian mixed logistic regression first (better for small samples)
            zero_model = BinomialBayesMixedGLM.from_formula(
                zero_formula, 
                groups=df[group_var],
                data=df
            )
            zero_fit = zero_model.fit_vb()
            zero_method = "Bayesian mixed logistic"
        except:
            # Fall back to regular mixed logistic regression
            zero_model = smf.mixedlm(zero_formula, df, groups=df[group_var])
            zero_fit = zero_model.fit()
            zero_method = "Mixed logistic"
        
        # 2. Count model (negative binomial for non-zeros)
        # For count data, we need to handle the distribution differently
        try:
            # Add offset term for normalization if requested
            if 'total_reads' in df.columns:
                offset_term = np.log(nonzero_df['total_reads'])
                count_formula = formula.replace(key_microbe, key_microbe) + " + offset(np.log(total_reads))"
                print("Using sequencing depth offset term in count model")
            else:
                count_formula = formula.replace(key_microbe, key_microbe)
            
            # Try GLM with Negative Binomial family
            count_model = smf.glm(
                count_formula, 
                data=nonzero_df, 
                family=sm.families.NegativeBinomial(link=sm.families.links.log())
            )
            count_fit = count_model.fit(method='newton', maxiter=100)
            count_method = "Negative binomial GLM"
            
            # Check if model converged
            if not count_fit.converged:
                # Try with different optimization method
                count_fit = count_model.fit(method='bfgs', maxiter=200)
                if not count_fit.converged:
                    raise ValueError("Negative binomial model failed to converge")
                
        except Exception as e:
            print(f"Count model error: {str(e)}")
            # Try Poisson as fallback
            try:
                count_model = smf.glm(
                    count_formula if 'count_formula' in locals() else formula.replace(key_microbe, key_microbe), 
                    data=nonzero_df, 
                    family=sm.families.Poisson(link=sm.families.links.log())
                )
                count_fit = count_model.fit()
                count_method = "Poisson GLM"
                
                # Check for overdispersion in Poisson model
                pearson_chi2 = count_fit.pearson_chi2 / count_fit.df_resid
                if pearson_chi2 > 2:
                    print(f"Warning: Poisson model shows overdispersion (Pearson χ²/df = {pearson_chi2:.2f})")
                
            except:
                # Fall back to regular linear model on log-transformed counts
                nonzero_df['log_counts'] = np.log(nonzero_df[key_microbe] + 0.5)  # Add small pseudocount
                count_model = smf.mixedlm(
                    formula.replace(key_microbe, "log_counts"), 
                    nonzero_df,
                    groups=nonzero_df[group_var]
                )
                try:
                    count_fit = count_model.fit()
                    count_method = "Linear mixed model on log-counts"
                except:
                    # Final fallback: simple OLS
                    count_model = smf.ols(
                        formula.replace(key_microbe, "log_counts"), 
                        data=nonzero_df
                    )
                    count_fit = count_model.fit()
                    count_method = "Log-linear OLS"
        
        # Store results
        results = {
            'microbe': key_microbe,
            'zero_model': zero_fit,
            'zero_method': zero_method,
            'count_model': count_fit,
            'count_method': count_method,
            'zero_proportion': 1 - (np.sum(df['presence']) / len(df))
        }
        
        if verbose:
            print(f"\n=== Zero-Inflated GLMM for {key_microbe} ===")
            print(f"Zero proportion: {results['zero_proportion']:.2%}")
            print(f"\n--- Part 1: {zero_method} (presence/absence) ---")
            print(zero_fit.summary())
            print(f"\n--- Part 2: {count_method} (abundance when present) ---")
            print(count_fit.summary())
        
        return results
    
    except Exception as e:
        print(f"❌ Error fitting zero-inflated GLMM for {key_microbe}: {str(e)}")
        return None

def extract_model_results(results_dict):
    """
    Extract key coefficients and statistics from model results.
    
    Parameters:
    - results_dict: Dictionary with model results from fit_zero_inflated_glmm
    
    Returns:
    - DataFrame with combined results
    """
    if results_dict is None:
        return None
    
    # Get coefficients from zero model (presence/absence)
    try:
        if hasattr(results_dict['zero_model'], 'params'):
            zero_params = pd.DataFrame({
                'Variable': results_dict['zero_model'].params.index,
                'Coefficient_Zero': results_dict['zero_model'].params.values,
                'P_value_Zero': results_dict['zero_model'].pvalues.values if hasattr(results_dict['zero_model'], 'pvalues') else np.nan,
                'Model_Type_Zero': results_dict['zero_method'],
                'Microbe': results_dict['microbe']
            })
        else:
            # Handle Bayesian model output differently
            zero_params = pd.DataFrame({
                'Variable': results_dict['zero_model'].mean_params.index,
                'Coefficient_Zero': results_dict['zero_model'].mean_params.values,
                'P_value_Zero': np.nan,  # Bayesian models don't have traditional p-values
                'Model_Type_Zero': results_dict['zero_method'],
                'Microbe': results_dict['microbe']
            })
    except:
        zero_params = pd.DataFrame()
    
    # Get coefficients from count model (abundance when present)
    try:
        count_params = pd.DataFrame({
            'Variable': results_dict['count_model'].params.index,
            'Coefficient_Count': results_dict['count_model'].params.values,
            'P_value_Count': results_dict['count_model'].pvalues.values if hasattr(results_dict['count_model'], 'pvalues') else np.nan,
            'Model_Type_Count': results_dict['count_method'],
            'Microbe': results_dict['microbe']
        })
    except:
        count_params = pd.DataFrame()
    
    # Merge zero and count model results
    if not zero_params.empty and not count_params.empty:
        # Full outer join on Variable and Microbe
        results = pd.merge(
            zero_params, count_params, 
            on=['Variable', 'Microbe'], 
            how='outer'
        )
    elif not zero_params.empty:
        results = zero_params
    elif not count_params.empty:
        results = count_params
    else:
        return None
    
    # Add zero proportion as a diagnostic
    results['Zero_Proportion'] = results_dict['zero_proportion']
    
    return results

def run_zero_inflated_analysis(data, key_organisms, categorical_features, subject_id_col="Subject"):
    """
    Run the complete zero-inflated GLMM analysis pipeline for multiple organisms.
    
    Parameters:
    - data: DataFrame with microbiome data and metadata
    - key_organisms: List of target microbes to analyze
    - categorical_features: List of categorical predictors
    - subject_id_col: Column name for subject ID (random effects grouping)
    
    Returns:
    - Combined results as a DataFrame
    """
    # Check inputs
    if not isinstance(key_organisms, list) or len(key_organisms) == 0:
        raise ValueError("key_organisms must be a non-empty list")
    
    all_results = []
    model_results_list = []
    
    for target_microbe in tqdm(key_organisms, desc="Analyzing microbes"):
        print(f"\n🦠 Analyzing: {target_microbe}")
        
        if target_microbe not in data.columns:
            print(f"⚠️ {target_microbe} not found in data. Skipping.")
            continue
        
        try:
            # Subset data for this microbe
            df = data[[target_microbe] + [col for col in categorical_features if col in data.columns] + [subject_id_col]].copy()
            df = df.dropna(subset=[target_microbe, subject_id_col])  # Ensure target and grouping variable aren't missing
            
            if df.empty:
                print(f"⚠️ No data available for {target_microbe} after dropping NAs. Skipping.")
                continue
            
            # Print zero proportion
            zero_prop = (df[target_microbe] == 0).mean()
            print(f"Zero proportion: {zero_prop:.2%}")
            
            # Skip microbes with too few non-zero values
            if (1 - zero_prop) * len(df) < 10:
                print(f"⚠️ Too few non-zero observations ({int((1-zero_prop)*len(df))}) for {target_microbe}. Skipping.")
                continue
            
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
            
            # Build formula
            fixed_effects = " + ".join(formula_parts)
            formula = f"{target_microbe} ~ {fixed_effects}"
            
            # Apply zero-inflated GLMM
            results = fit_zero_inflated_glmm(formula, subject_id_col, df, target_microbe)
            
            if results is not None:
                model_results_list.append(results)
                
                # Extract and save results
                result_df = extract_model_results(results)
                if result_df is not None:
                    all_results.append(result_df)
                
                # Save individual model results as text
                with open(f"zero_inflated_glmm_{target_microbe}.txt", "w") as f:
                    f.write(f"=== Zero-Inflated GLMM for {target_microbe} ===\n")
                    f.write(f"Zero proportion: {results['zero_proportion']:.2%}\n\n")
                    f.write(f"--- Part 1: {results['zero_method']} (presence/absence) ---\n")
                    f.write(str(results['zero_model'].summary()) + "\n\n")
                    f.write(f"--- Part 2: {results['count_method']} (abundance when present) ---\n")
                    f.write(str(results['count_model'].summary()))
        
        except Exception as e:
            print(f"❌ Failed to analyze {target_microbe}: {str(e)}")
            import traceback
            traceback.print_exc()
    
    # Save model objects for later use
    if model_results_list:
        # Create results directory if it doesn't exist
        os.makedirs("model_results", exist_ok=True)
        
        # Save model objects with pickle
        with open("model_results/zero_inflated_model_objects.pkl", "wb") as f:
            pickle.dump(model_results_list, f)
    
    # Combine all results
    if all_results:
        combined_results = pd.concat(all_results, ignore_index=True)
        combined_results.to_csv("zero_inflated_glmm_all_microbes.csv", index=False)
        print(f"✅ Saved results for {len(model_results_list)} microbes")
        return combined_results
    else:
        print("⚠️ No results to save.")
        return None

# Main analysis code

if __name__ == "__main__":
    print("📊 Starting Zero-Inflated GLMM Analysis Pipeline")
    
    # Load microbiome count data
    microbiome_df = pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)
    
    # Load clinical metadata
    metadata_df = pd.read_csv("../metadata/AllNICUSampleKey20250206.csv", index_col=0)
    subject_id_col = "Subject"  
    
    # Print data diagnostics
    print(f"Microbiome data shape: {microbiome_df.shape}")
    print(f"Metadata shape: {metadata_df.shape}")
    
    # Select categorical features for modeling
    categorical_features = ["SampleType", "Location", "GestationCohort", "SampleCollectionWeek", 
                           "MaternalAntibiotics", "PostNatalAbxCohort", "BSI_30D", "NEC_30D", 
                           "AnyMilk", "PICC", "UVC"]
    
    # Add total read counts for offset terms
    microbiome_df['total_reads'] = microbiome_df.sum(axis=1)
    
    # Auto-select and apply transformation for exploratory analysis
    # For the ZINB models, we'll use raw counts, but transformed data can be useful for visualization
    transformed_df, transform_method = select_best_transformation(microbiome_df)
    print(f"Applied {transform_method} to microbiome data")
    
    # Merge both transformed and raw count data with metadata
    data_transformed = transformed_df.merge(metadata_df, left_index=True, right_index=True, how="inner")
    data = microbiome_df.merge(metadata_df, left_index=True, right_index=True, how="inner")
    
    print(f"Merged data shape: {data.shape}")
    print(f"Proportion of zeros across all taxa: {(data.iloc[:, :microbiome_df.shape[1]] == 0).mean().mean():.2%}")
    
    # Key organisms to analyze
    key_organisms = ["Klebsiella.pneumoniae", "Staphylococcus.aureus", "Escherichia.coli", 
                     "Klebsiella.oxytoca", "Staphylococcus.epidermidis", "Streptococcus.pyogenes"]
    
    # Run zero-inflated GLMM analysis
    results = run_zero_inflated_analysis(data, key_organisms, categorical_features, subject_id_col)
    
    # Additional visualizations
    if results is not None:
        # Create heatmap of coefficients across models
        # First, reshape the data for the heatmap
        coef_data = pd.pivot_table(
            results,
            values=["Coefficient_Zero", "Coefficient_Count"],
            index="Microbe", 
            columns="Variable"
        )
        
        # Plot heatmap for zero model coefficients
        plt.figure(figsize=(12, 8))
        sns.heatmap(coef_data["Coefficient_Zero"].transpose(), 
                   cmap="coolwarm", center=0, annot=True, 
                   fmt=".2f", linewidths=.5)
        plt.title("Zero Model Coefficients (Presence/Absence)")
        plt.tight_layout()
        plt.savefig("zero_model_coefficients_heatmap.pdf")
        plt.close()
        
        # Plot heatmap for count model coefficients
        plt.figure(figsize=(12, 8))
        sns.heatmap(coef_data["Coefficient_Count"].transpose(), 
                   cmap="coolwarm", center=0, annot=True, 
                   fmt=".2f", linewidths=.5)
        plt.title("Count Model Coefficients (Abundance When Present)")
        plt.tight_layout()
        plt.savefig("count_model_coefficients_heatmap.pdf")
        plt.close()
    
    print("✅ Analysis complete!")
