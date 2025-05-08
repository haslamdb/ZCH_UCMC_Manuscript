import pickle
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def find_lmm_results_file():
    """Search for the LMM results file."""
    possible_paths = [
        "results/variance_analysis/lmm_microbiome_community_results.pkl",
        "lmm_microbiome_community_results.pkl"
    ]
    
    # Also search subdirectories
    for root, dirs, files in os.walk('.'):
        for file in files:
            if file == 'lmm_microbiome_community_results.pkl':
                possible_paths.append(os.path.join(root, file))
    
    for path in possible_paths:
        if os.path.exists(path):
            print(f"Found LMM results file at: {path}")
            return path
    
    return None

def print_lmm_variance_explained(file_path):
    """Load and print the variance explained from LMM results."""
    try:
        with open(file_path, 'rb') as f:
            lmm_results = pickle.load(f)
        
        print("\n===== LMM RESULTS: VARIANCE EXPLAINED =====\n")
        
        # Create a summary table
        summary_rows = []
        
        # Process each principal component
        for pc, result in lmm_results.items():
            # Extract R-squared if available
            r_squared = result.get('r_squared')
            if r_squared is not None:
                print(f"{pc} explained variance (R²): {r_squared:.4f} ({r_squared*100:.2f}%)")
                
                # Get significant variables
                coeffs = result['coefficients']
                significant = coeffs[coeffs['Significant']]
                
                print(f"  Number of significant predictors: {len(significant)}")
                
                # Print top predictors
                if len(significant) > 0:
                    top_3 = significant.head(3)
                    print("  Top predictors:")
                    for _, row in top_3.iterrows():
                        var_name = row['Variable']
                        coef = row['Coefficient']
                        p_val = row['P-value']
                        print(f"    - {var_name}: Coef={coef:.4f}, p={p_val:.4f}")
                
                # Add to summary for the heatmap
                summary_rows.append({
                    'PC': pc,
                    'R-squared': r_squared,
                    'Significant Variables': len(significant)
                })
                
                # Add additional rows for significant variables
                for _, row in significant.iterrows():
                    var_name = row['Variable']
                    # Skip intercept
                    if var_name == 'Intercept':
                        continue
                    
                    # Extract variable group name
                    if 'C(' in var_name:
                        var_group = var_name.split(')')[0].replace('C(', '')
                    else:
                        var_group = var_name
                    
                    # Add contribution to that PC
                    contribution = {
                        'Variable': var_group,
                        'PC': pc,
                        'Coefficient': abs(row['Coefficient']),
                        'P-value': row['P-value']
                    }
                    summary_rows.append(contribution)
            else:
                print(f"{pc}: R-squared not available")
        
        # Create a heatmap of variable importance across PCs
        create_variable_importance_heatmap(lmm_results)
        
    except Exception as e:
        print(f"Error processing LMM results file: {str(e)}")

def create_variable_importance_heatmap(lmm_results):
    """Create a heatmap showing variable importance across principal components."""
    try:
        # Extract all unique variable groups
        all_vars = set()
        for pc, result in lmm_results.items():
            coeffs = result['coefficients']
            for _, row in coeffs.iterrows():
                if row['Variable'] in ['Intercept', 'Group Var']:
                    continue
                    
                if 'C(' in row['Variable']:
                    var_group = row['Variable'].split(')')[0].replace('C(', '')
                else:
                    var_group = row['Variable']
                    
                all_vars.add(var_group)
        
        # Create a matrix of coefficients
        heatmap_data = pd.DataFrame(0, 
                                   index=sorted(all_vars),
                                   columns=sorted(lmm_results.keys()))
        
        # Fill in coefficient values (only if significant)
        for pc, result in lmm_results.items():
            coeffs = result['coefficients']
            for _, row in coeffs.iterrows():
                if row['Variable'] in ['Intercept', 'Group Var']:
                    continue
                    
                if 'C(' in row['Variable']:
                    var_group = row['Variable'].split(')')[0].replace('C(', '')
                else:
                    var_group = row['Variable']
                
                # Use the absolute value of the coefficient if significant, 0 otherwise
                if row['Significant']:
                    heatmap_data.loc[var_group, pc] = abs(row['Coefficient'])
        
        # Calculate total importance of each variable
        heatmap_data['Total_Importance'] = heatmap_data.sum(axis=1)
        
        # Sort by total importance
        heatmap_data = heatmap_data.sort_values('Total_Importance', ascending=False)
        
        # Remove Total_Importance column for plotting
        plot_data = heatmap_data.drop(columns=['Total_Importance'])
        
        # Create heatmap
        plt.figure(figsize=(12, max(8, len(all_vars) * 0.4)))
        sns.heatmap(plot_data, 
                   cmap='viridis',
                   annot=True,
                   fmt='.2f',
                   cbar_kws={'label': 'Absolute Coefficient Value (if significant)'})
        
        plt.title('Clinical Factors Influencing Microbiome Principal Components')
        plt.tight_layout()
        
        # Save figure
        output_file = "lmm_variance_explained_heatmap.pdf"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"\nHeatmap visualization saved as {output_file}")
        
        # Print total importance of each variable
        print("\nTotal importance of each variable across all PCs:")
        importance_df = heatmap_data[['Total_Importance']].reset_index()
        importance_df.columns = ['Variable', 'Total_Importance']
        importance_df = importance_df.sort_values('Total_Importance', ascending=False)
        
        for _, row in importance_df.iterrows():
            print(f"  {row['Variable']}: {row['Total_Importance']:.4f}")
        
        # Create bar chart of total variable importance
        plt.figure(figsize=(10, max(6, len(all_vars) * 0.3)))
        sns.barplot(x='Total_Importance', y='Variable', data=importance_df)
        plt.title('Total Importance of Clinical Factors Across All Principal Components')
        plt.tight_layout()
        
        # Save figure
        output_file2 = "lmm_total_variable_importance.pdf"
        plt.savefig(output_file2, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Total importance visualization saved as {output_file2}")
        
    except Exception as e:
        print(f"Error creating heatmap: {str(e)}")

if __name__ == "__main__":
    print("Searching for LMM results file...")
    file_path = find_lmm_results_file()
    
    if file_path:
        print_lmm_variance_explained(file_path)
    else:
        print("\nCould not find LMM results file.")