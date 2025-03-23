import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import squareform, pdist
import skbio.stats.ordination as ordination
from skbio.stats.distance import DistanceMatrix
from skbio.stats.composition import clr
import statsmodels.api as sm
from statsmodels.formula.api import ols

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
    microbiome_df = pd.read_csv("microbiome_abundance.csv", index_col=0)
    microbiome_clr = mt.clr_transform(microbiome_df)
except (ImportError, NameError):
    # Fall back to our defined function
    print("microbiome_transform module not found, using built-in CLR function")
    microbiome_df = pd.read_csv("microbiome_abundance.csv", index_col=0)
    microbiome_clr = clr_transform(microbiome_df)

# Load metadata
print("Loading metadata...")
metadata_df = pd.read_csv("metadata.csv", index_col=0)

# Merge transformed microbiome data and metadata on Sample ID
print("Merging datasets...")
# Ensure indices are the same type (string)
microbiome_clr.index = microbiome_clr.index.astype(str)
metadata_df.index = metadata_df.index.astype(str)

# Merge using index (assuming index is the sample ID)
merged_data = microbiome_clr.copy()
merged_data = pd.merge(
    microbiome_clr, 
    metadata_df, 
    left_index=True, 
    right_index=True, 
    how='inner'
)

print(f"Final dataset shape: {merged_data.shape}")

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
distance_matrix = DistanceMatrix(pdist(microbiome_data, metric='euclidean'))

# PERMANOVA analysis for each categorical feature
print("Performing PERMANOVA analysis...")
permanova_results = {}

for feature in available_features:
    # Drop NAs for the feature being tested
    valid_samples = merged_data[feature].dropna().index
    if len(valid_samples) < 5:
        print(f"Skipping {feature} - not enough valid samples")
        continue
        
    # Subset distance matrix and grouping variable
    feature_dm = distance_matrix.filter(valid_samples)
    grouping = merged_data.loc[valid_samples, feature].astype(str)
    
    # Skip if only one unique value
    if len(grouping.unique()) < 2:
        print(f"Skipping {feature} - only one unique value")
        continue
    
    try:
        # Run PERMANOVA (default: 999 permutations)
        result = skbio.stats.distance.permanova(feature_dm, grouping)
        permanova_results[feature] = {
            'test_statistic': result['test statistic'],
            'p_value': result['p-value'],
            'R2': result['test statistic'] / (result['test statistic'] + result['denominator'])
        }
        print(f"{feature}: p-value = {result['p-value']:.4f}, R² = {permanova_results[feature]['R2']:.4f}")
    except Exception as e:
        print(f"Error analyzing {feature}: {e}")

# Create PCoA plot
print("Creating ordination plots...")
pcoa_result = ordination.pcoa(distance_matrix)
pcoa_df = pd.DataFrame(
    data=pcoa_result.samples.values,
    columns=['PC1', 'PC2', 'PC3', 'PC4', 'PC5'][:pcoa_result.samples.shape[1]],
    index=distance_matrix.ids
)

# Merge with metadata for visualization
pcoa_with_metadata = pd.merge(
    pcoa_df, 
    metadata_df, 
    left_index=True, 
    right_index=True, 
    how='inner'
)

# Sort features by significance
significant_features = sorted(
    [(feature, results['p_value'], results['R2']) 
     for feature, results in permanova_results.items()],
    key=lambda x: x[1]  # Sort by p-value
)

# Create a PERMANOVA summary plot
print("Creating PERMANOVA summary plot...")
plt.figure(figsize=(10, 6))
feature_names = [f[0] for f in significant_features]
r2_values = [f[2] for f in significant_features]
p_values = [f[1] for f in significant_features]

# Create the R² bar plot
ax = plt.barh(feature_names, r2_values, color=['#2c7bb6' if p < 0.05 else '#d7191c' for p in p_values])
plt.xlabel('R² (Proportion of Variance Explained)')
plt.ylabel('Metadata Features')
plt.title('PERMANOVA Results: Variance Explained by Metadata Features')

# Add significance annotations
for i, p in enumerate(p_values):
    significance = "**" if p < 0.01 else ("*" if p < 0.05 else "ns")
    plt.text(r2_values[i] + 0.005, i, significance, va='center')

# Add R² text
for i, r2 in enumerate(r2_values):
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

# Choose the most significant feature for visualization in PCoA plot
if significant_features:
    top_feature = significant_features[0][0]
    
    plt.figure(figsize=(10, 8))
    
    # Create a scatter plot with colored points for each category in the top feature
    categories = pcoa_with_metadata[top_feature].astype(str).unique()
    
    if len(categories) <= 10:  # Only create if we don't have too many categories
        plt.figure(figsize=(12, 10))
        
        for i, category in enumerate(sorted(categories)):
            mask = pcoa_with_metadata[top_feature].astype(str) == category
            plt.scatter(
                pcoa_with_metadata.loc[mask, 'PC1'],
                pcoa_with_metadata.loc[mask, 'PC2'],
                s=50, 
                alpha=0.7,
                label=f"{category} (n={sum(mask)})"
            )
        
        # Add labels and legend
        plt.xlabel(f"PC1 ({pcoa_result.proportion_explained[0]:.2%} explained var.)")
        plt.ylabel(f"PC2 ({pcoa_result.proportion_explained[1]:.2%} explained var.)")
        plt.title(f"PCoA colored by {top_feature} (p={significant_features[0][1]:.4f}, R²={significant_features[0][2]:.4f})")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f'pcoa_by_{top_feature}.png', dpi=300, bbox_inches='tight')

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
