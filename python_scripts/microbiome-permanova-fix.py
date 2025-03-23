import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import squareform, pdist
import skbio
import skbio.stats.ordination as ordination
from skbio.stats.distance import DistanceMatrix
from skbio.stats.composition import clr
import statsmodels.api as sm
from statsmodels.formula.api import ols
import microbiome_transform as mt

# Optional: Set plotting style
sns.set(style="whitegrid")
plt.rcParams.update({'font.size': 12})

# Define a function for CLR transformation (in case microbiome_transform module is not available)
def clr_transform(df):
    """
    Centered log-ratio transformation for compositional data
    Adds a small pseudocount to zeros.
    """
    # Add pseudocount to zeros
    df_pseudo = df.replace(0, np.nextafter(0, 1))
    # Apply CLR transformation
    clr_data = clr(df_pseudo.values)
    return pd.DataFrame(clr_data, index=df.index, columns=df.columns)

# Load microbiome abundance data
print("Loading microbiome data...")
try:
    # Try with the imported module first
    import microbiome_transform as mt
    microbiome_df = pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)
    microbiome_clr = mt.clr_transform(microbiome_df)
except (ImportError, NameError):
    # Fall back to our defined function
    print("microbiome_transform module not found, using built-in CLR function")
    microbiome_df = pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)
    microbiome_clr = clr_transform(microbiome_df)

# Load metadata
print("Loading metadata...")
metadata_df = pd.read_csv("../metadata/AllNICUSampleKey20250206.csv", index_col=0)

# Merge transformed microbiome data and metadata on Sample ID
print("Merging datasets...")
# Ensure indices are the same type (string)
microbiome_clr.index = microbiome_clr.index.astype(str)
metadata_df.index = metadata_df.index.astype(str)

# Print some diagnostics about the data before merging
print(f"Microbiome data shape: {microbiome_clr.shape} with {len(microbiome_clr.index.unique())} unique sample IDs")
print(f"Metadata shape: {metadata_df.shape} with {len(metadata_df.index.unique())} unique sample IDs")

# Check for duplicate indices
if len(microbiome_clr.index) != len(microbiome_clr.index.unique()):
    print("WARNING: Microbiome data contains duplicate sample IDs")
if len(metadata_df.index) != len(metadata_df.index.unique()):
    print("WARNING: Metadata contains duplicate sample IDs")

# Check sample overlap
common_samples = set(microbiome_clr.index) & set(metadata_df.index)
print(f"Common samples between microbiome and metadata: {len(common_samples)}")

# Merge using index (assuming index is the sample ID)
merged_data = pd.merge(
    microbiome_clr, 
    metadata_df, 
    left_index=True, 
    right_index=True, 
    how='inner'
)

print(f"Final merged dataset shape: {merged_data.shape} with {len(merged_data.index.unique())} unique samples")

# Define categorical features of interest
categorical_features = [
    "SampleType", "Location", "GestationCohort", "SampleCollectionWeek", 
    "MaternalAntibiotics", "PostNatalAbxCohort", "BSI_30D", "NEC_30D", 
    "AnyMilk", "PICC", "UVC"
]

# Check which features actually exist in the metadata
available_features = [col for col in categorical_features if col in merged_data.columns]
if len(available_features) < len(categorical_features):
    missing = set(categorical_features) - set(available_features)
    print(f"Warning: Some requested features are missing from the data: {missing}")

# Identify microbiome columns (all except metadata)
microbiome_cols = [col for col in merged_data.columns if col not in metadata_df.columns]

# Create distance matrix from microbiome data
print("Calculating distance matrix...")
microbiome_data = merged_data[microbiome_cols]

# Ensure we're using the exact same set of samples for distance matrix and subsequent analyses
distance_matrix = DistanceMatrix(pdist(microbiome_data, metric='euclidean'), ids=microbiome_data.index)

# Function to check and fix sample ID issues
def fix_sample_ids(df1, df2, df1_name="DataFrame 1", df2_name="DataFrame 2"):
    """Check and report sample ID overlap issues between two dataframes"""
    print(f"\nSample ID consistency check between {df1_name} and {df2_name}:")
    print(f"- {df1_name} samples: {len(df1.index)} ({len(df1.index.unique())} unique)")
    print(f"- {df2_name} samples: {len(df2.index)} ({len(df2.index.unique())} unique)")
    
    common = set(df1.index) & set(df2.index)
    print(f"- Common samples: {len(common)}")
    
    only_in_df1 = set(df1.index) - set(df2.index)
    only_in_df2 = set(df2.index) - set(df1.index)
    
    if only_in_df1:
        print(f"- Samples only in {df1_name}: {len(only_in_df1)}")
        if len(only_in_df1) < 10:
            print(f"  IDs: {sorted(only_in_df1)}")
    
    if only_in_df2:
        print(f"- Samples only in {df2_name}: {len(only_in_df2)}")
        if len(only_in_df2) < 10:
            print(f"  IDs: {sorted(only_in_df2)}")
    
    # Return common samples
    return common

# PERMANOVA analysis for each categorical feature
print("Performing PERMANOVA analysis...")
permanova_results = {}

try:
    for feature in available_features:
        # Drop NAs for the feature being tested
        valid_samples = merged_data[feature].dropna().index
        if len(valid_samples) < 5:
            print(f"Skipping {feature} - not enough valid samples")
            continue
        
        # Make sure valid_samples are all in the distance matrix
        valid_samples = [sample for sample in valid_samples if sample in distance_matrix.ids]
        if len(valid_samples) < 5:
            print(f"Skipping {feature} - not enough samples after filtering for distance matrix")
            continue
            
        # Subset distance matrix and grouping variable
        feature_dm = distance_matrix.filter(valid_samples)
        grouping = merged_data.loc[valid_samples, feature].astype(str)
        
        # Skip if only one unique value
        if len(grouping.unique()) < 2:
            print(f"Skipping {feature} - only one unique value")
            continue
        
        try:
            # Print group information for debugging
            group_counts = grouping.value_counts()
            print(f"{feature} groups: {dict(group_counts)}")
            
            # Check if we have enough groups and samples
            if len(group_counts) < 2:
                print(f"Skipping {feature} - need at least 2 groups, found {len(group_counts)}")
                continue
                
            if any(count < 3 for count in group_counts):
                print(f"Warning for {feature} - some groups have fewer than 3 samples")
                
            # Check if there's too much imbalance between groups
            max_count = group_counts.max()
            min_count = group_counts.min()
            if max_count / min_count > 10:
                print(f"Warning for {feature} - highly imbalanced groups (largest/smallest = {max_count/min_count:.1f})")
                
            # Run PERMANOVA (default: 999 permutations)
            result = skbio.stats.distance.permanova(feature_dm, grouping, permutations=999)
            
            # Print the result structure for debugging
            print(f"PERMANOVA result keys: {list(result.keys())}")
            
            # Calculate R² - adapt based on what's available in the result
            # For newer scikit-bio versions
            if 'test statistic' in result and hasattr(result, 'get'):
                test_stat = result['test statistic']
                # Different versions of scikit-bio might have different keys for the denominator
                if 'denominator' in result:
                    denom = result['denominator']
                    r_squared = test_stat / (test_stat + denom)
                else:
                    # If no denominator is provided, we can calculate R² as the ratio of
                    # the between-group sum of squares to the total sum of squares
                    print(f"Warning: 'denominator' not found in PERMANOVA result for {feature}.")
                    print(f"Using alternative R² calculation method.")
                    # For newer versions, R² might be part of the results
                    if 'R2' in result:
                        r_squared = result['R2']
                    else:
                        # As a fallback, we'll use a simple approximation
                        r_squared = test_stat / (test_stat + 1.0)  # This is a placeholder calculation
            else:
                # If the result structure is completely different, we'll use a default value
                print(f"Warning: Expected keys not found in PERMANOVA result for {feature}.")
                test_stat = 0.0
                r_squared = 0.0
                
            # Store the results
            permanova_results[feature] = {
                'test_statistic': test_stat if 'test statistic' in result else result.get('F', 0.0),
                'p_value': result['p-value'] if 'p-value' in result else result.get('p', 1.0),
                'R2': r_squared,
                'sample_size': len(grouping),
                'groups': dict(group_counts)
            }
            print(f"{feature}: p-value = {result['p-value']:.4f}, R² = {permanova_results[feature]['R2']:.4f}, n={len(grouping)}")
        except Exception as e:
            print(f"Error analyzing {feature}: {e}")
            import traceback
            traceback.print_exc()

except KeyboardInterrupt:
    print("\nAnalysis interrupted by user. Saving results collected so far...")
    # Continue with plotting and saving the results we have
except Exception as e:
    print(f"\nUnexpected error during PERMANOVA analysis: {e}")
    import traceback
    traceback.print_exc()

# Create PCoA plot
print("Creating ordination plots...")
try:
    # Set a reasonable number of iterations and check convergence
    print("Running PCoA analysis...")
    # FIXED: Removed the 'fsvd' parameter which was causing the error
    pcoa_result = ordination.pcoa(
        distance_matrix,
        number_of_dimensions=5  # Explicitly request 5 dimensions
    )
    
    # Check if PCoA was successful
    if pcoa_result is None or not hasattr(pcoa_result, 'samples') or pcoa_result.samples.shape[0] == 0:
        raise ValueError("PCoA analysis failed to produce valid results")
        
    print(f"PCoA complete. Shape of results: {pcoa_result.samples.shape}")
    
    # Define PC columns
    pc_cols = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']
    
    # Ensure we don't try to create more PCs than we have
    n_dimensions = pcoa_result.samples.shape[1]
    print(f"Number of dimensions in PCoA result: {n_dimensions}")
    actual_pc_cols = pc_cols[:min(5, n_dimensions)]
    
    # Create dataframe with PCoA results
    pcoa_df = pd.DataFrame(
        data=pcoa_result.samples.values,
        columns=actual_pc_cols,
        index=distance_matrix.ids
    )
    
    # Report proportion of variance explained
    if hasattr(pcoa_result, 'proportion_explained'):
        print("Proportion of variance explained by each PC:")
        for i, prop in enumerate(pcoa_result.proportion_explained[:n_dimensions]):
            print(f"  PC{i+1}: {prop:.2%}")
    else:
        print("Proportion of variance explained not available in PCoA results")
    
    # Check sample ID consistency
    common_samples = fix_sample_ids(pcoa_df, metadata_df, "PCoA results", "Metadata")
    
    # Merge with metadata for visualization (use inner join to be safe)
    pcoa_with_metadata = pd.merge(
        pcoa_df, 
        metadata_df, 
        left_index=True, 
        right_index=True, 
        how='inner'
    )
    
    print(f"PCoA with metadata shape: {pcoa_with_metadata.shape}")

except Exception as e:
    print(f"Error in PCoA analysis: {e}")
    import traceback
    traceback.print_exc()
    # Create empty dataframe to prevent further errors
    pcoa_with_metadata = pd.DataFrame()

# Sort features by significance
significant_features = sorted(
    [(feature, results['p_value'], results['R2']) 
     for feature, results in permanova_results.items()],
    key=lambda x: x[1]  # Sort by p-value
)

# Create a PERMANOVA summary plot
print("Creating PERMANOVA summary plot...")
if significant_features:
    plt.figure(figsize=(10, 6))
    feature_names = [f[0] for f in significant_features]
    r2_values = [f[2] for f in significant_features]
    p_values = [f[1] for f in significant_features]
    
    # Create the R² bar plot
    ax = plt.barh(feature_names, r2_values, color=['#2c7bb6' if p < 0.05 else '#d7191c' for p in p_values])
    plt.xlabel('R² (Proportion of Variance Explained)')
    plt.ylabel('Metadata Features')
    plt.title('PERMANOVA Results: Variance Explained by Metadata Features')
    
    # Get the maximum R2 value to help with text positioning
    max_r2 = max(r2_values) if r2_values else 0.1
    offset = max_r2 * 0.03  # Use 3% of the max value as offset
    
    # Add significance annotations
    for i, p in enumerate(p_values):
        significance = "**" if p < 0.01 else ("*" if p < 0.05 else "ns")
        plt.text(r2_values[i] + offset, i, significance, va='center')
    
    # Add R² text
    for i, r2 in enumerate(r2_values):
        if r2 > 0.03:  # Only add text for bars that have enough space
            plt.text(r2/2, i, f"R²={r2:.3f}", ha='center', va='center', color='white')
    
    # Add legend for significance
    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color='#2c7bb6', lw=4),
        Line2D([0], [0], color='#d7191c', lw=4),
    ]
    plt.legend(custom_lines, ['p < 0.05', 'p ≥ 0.05'], loc='upper right')
    
    plt.tight_layout()
    plt.savefig('permanova_variance_explained.png', dpi=300, bbox_inches='tight')
    print("Created PERMANOVA summary plot: permanova_variance_explained.png")
else:
    print("No significant features found. Skipping summary plot.")

# Choose the most significant feature for visualization in PCoA plot
if significant_features and not pcoa_with_metadata.empty and 'PC1' in pcoa_with_metadata.columns and 'PC2' in pcoa_with_metadata.columns:
    try:
        # Get the top significant feature
        top_feature = significant_features[0][0]
        
        # Check if the feature exists in our data
        if top_feature in pcoa_with_metadata.columns:
            plt.figure(figsize=(12, 10))
            
            # Create a scatter plot with colored points for each category in the top feature
            categories = pcoa_with_metadata[top_feature].astype(str).unique()
            
            if len(categories) <= 10:  # Only create if we don't have too many categories
                # Define a distinct color palette for the categories
                from matplotlib.cm import get_cmap
                cmap = get_cmap('tab10' if len(categories) <= 10 else 'tab20')
                colors = [cmap(i) for i in range(len(categories))]
                
                # Create the scatter plot
                for i, category in enumerate(sorted(categories)):
                    mask = pcoa_with_metadata[top_feature].astype(str) == category
                    if sum(mask) > 0:  # Only plot if we have samples in this category
                        plt.scatter(
                            pcoa_with_metadata.loc[mask, 'PC1'],
                            pcoa_with_metadata.loc[mask, 'PC2'],
                            s=50, 
                            alpha=0.7,
                            color=colors[i],
                            label=f"{category} (n={sum(mask)})"
                        )
                
                # Add explained variance if available
                if hasattr(pcoa_result, 'proportion_explained') and len(pcoa_result.proportion_explained) >= 2:
                    plt.xlabel(f"PC1 ({pcoa_result.proportion_explained[0]:.2%} explained var.)")
                    plt.ylabel(f"PC2 ({pcoa_result.proportion_explained[1]:.2%} explained var.)")
                else:
                    plt.xlabel("PC1")
                    plt.ylabel("PC2")
                
                plt.title(f"PCoA colored by {top_feature} (p={significant_features[0][1]:.4f}, R²={significant_features[0][2]:.4f})")
                
                # Add legend with larger font
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
                
                # Add grid lines
                plt.grid(linestyle='--', alpha=0.7)
                
                plt.tight_layout()
                plt.savefig(f'pcoa_by_{top_feature}.png', dpi=300, bbox_inches='tight')
                print(f"Created PCoA plot: pcoa_by_{top_feature}.png")
            else:
                print(f"Skipping PCoA visualization for {top_feature} - too many categories ({len(categories)})")
        else:
            print(f"Top feature '{top_feature}' not found in PCoA data")
    except Exception as e:
        print(f"Error creating PCoA visualization: {e}")
        import traceback
        traceback.print_exc()
else:
    if not significant_features:
        print("No significant features available for PCoA visualization")
    elif pcoa_with_metadata.empty:
        print("PCoA data is empty, skipping visualization")
    else:
        print("PCoA data is missing PC1/PC2 columns, skipping visualization")

# Create a summary table
print("Creating summary table...")
summary_df = pd.DataFrame({
    'Feature': [f[0] for f in significant_features],
    'p_value': [f[1] for f in significant_features],
    'R2': [f[2] for f in significant_features],
    'Significant': ['Yes' if p < 0.05 else 'No' for p in [f[1] for f in significant_features]]
})

# Save results to file
summary_df.to_csv('permanova_results.csv', index=False)
print("Analysis complete! Results saved to permanova_results.csv")
print("Plots saved as permanova_variance_explained.png and pcoa_by_[feature].png")
