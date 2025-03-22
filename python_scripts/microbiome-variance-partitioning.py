import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import scipy.stats as stats
from skbio.stats.distance import permanova, DistanceMatrix
from skbio.stats.ordination import pcoa
import microbiome_transform as mt
import os
from tqdm import tqdm
import statsmodels.api as sm
import statsmodels.formula.api as smf
import warnings
import scipy.spatial.distance as ssd
from functools import partial
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

# Activate pandas conversion for R
pandas2ri.activate()

# Import vegan through rpy2
vegan = importr('vegan')
import pickle
import re


# Define wrappers for R functions
def adonis2(formula, data, permutations=999, method="bray", by="terms", parallel=1):
    """
    Wrapper for vegan's adonis2 function
    """
    # Convert formula to R formula
    r_formula = robjects.Formula(formula)
    
    # Run adonis2
    result = vegan.adonis2(r_formula, data=data, permutations=permutations, 
                          method=method, by=by, parallel=parallel)
    
    return result

def rda(formula, data, scale=True):
    """
    Wrapper for vegan's rda function
    """
    # Convert formula to R formula
    r_formula = robjects.Formula(formula)
    
    # Run rda
    result = vegan.rda(r_formula, data=data, scale=scale)
    
    return result

def varpart(y, *args):
    """
    Wrapper for vegan's varpart function
    """
    # Run varpart
    result = vegan.varpart(y, *args)
    
    return result

# Silence common warnings
warnings.filterwarnings("ignore", category=UserWarning)

def select_best_normalization(df, verbose=True):
    """
    Selects the most appropriate normalization method based on data characteristics.
    
    Parameters:
    - df: DataFrame with samples as rows, taxa as columns
    - verbose: Whether to print diagnostic messages
    
    Returns:
    - Normalized DataFrame and the name of the method used
    """
    # Calculate data characteristics
    zero_fraction = (df == 0).sum().sum() / df.size
    total_reads = df.sum(axis=1)
    cv_depth = total_reads.std() / total_reads.mean()
    
    if verbose:
        print(f"Data diagnostics:")
        print(f"  - Dimensions: {df.shape[0]} samples × {df.shape[1]} taxa")
        print(f"  - Zero fraction: {zero_fraction:.2%}")
        print(f"  - Mean sequencing depth: {total_reads.mean():.1f} ± {total_reads.std():.1f}")
        print(f"  - Coefficient of variation in sequencing depth: {cv_depth:.2f}")
    
    # Decision logic for normalization method
    if zero_fraction > 0.80:
        if verbose:
            print("⚠️ Data is extremely sparse (>80% zeros). Using presence/absence transformation.")
        # Convert to presence/absence
        normalized_df = (df > 0).astype(int)
        method = "presence_absence"
    elif zero_fraction > 0.50:
        if verbose:
            print("⚠️ Data is very sparse (>50% zeros). Using Hellinger transformation.")
        # Hellinger transformation (square root of relative abundance)
        normalized_df = np.sqrt(df.div(df.sum(axis=1), axis=0))
        method = "hellinger"
    elif cv_depth > 0.5:
        if verbose:
            print("⚠️ Large variation in sequencing depth. Using rarefaction.")
        # Rarefaction to minimum depth
        try:
            min_depth = int(total_reads.min())
            if min_depth < 1000:
                min_depth = 1000  # Set minimum threshold
            normalized_df = mt.rarefaction(df, min_depth=min_depth)
            method = "rarefaction"
        except Exception as e:
            print(f"Error in rarefaction: {e}. Falling back to TSS normalization.")
            normalized_df = mt.tss_transform(df)
            method = "tss"
    else:
        if verbose:
            print("✅ Using Centered Log Ratio (CLR) transformation.")
        # CLR transformation
        normalized_df = mt.clr_transform(df)
        method = "clr"
    
    return normalized_df, method

def calculate_distance_matrix(normalized_df, method="bray", binary=False):
    """
    Calculate distance matrix using specified method.
    
    Parameters:
    - normalized_df: Normalized DataFrame
    - method: Distance metric ('bray', 'jaccard', 'unifrac', etc.)
    - binary: Whether to use presence/absence data
    
    Returns:
    - DistanceMatrix object
    """
    if method == "bray":
        # Bray-Curtis dissimilarity
        distances = ssd.pdist(normalized_df, metric="braycurtis")
        dist_matrix = ssd.squareform(distances)
    elif method == "jaccard":
        # Jaccard distance
        if binary:
            binary_df = (normalized_df > 0).astype(int)
            distances = ssd.pdist(binary_df, metric="jaccard")
        else:
            distances = ssd.pdist(normalized_df, metric="jaccard")
        dist_matrix = ssd.squareform(distances)
    elif method == "euclidean":
        # Euclidean distance (good for CLR-transformed data)
        distances = ssd.pdist(normalized_df, metric="euclidean")
        dist_matrix = ssd.squareform(distances)
    else:
        raise ValueError(f"Unsupported distance method: {method}")
    
    # Create SkBio DistanceMatrix object
    return DistanceMatrix(dist_matrix, ids=normalized_df.index)

def run_permanova(distance_matrix, metadata, formula, permutations=999):
    """
    Run PERMANOVA to determine factors that explain variation in microbial community.
    
    Parameters:
    - distance_matrix: DistanceMatrix object
    - metadata: DataFrame with sample metadata
    - formula: Formula string (e.g., "~ Location + SampleType")
    - permutations: Number of permutations for significance testing
    
    Returns:
    - PERMANOVA results
    """
    # Make sure metadata index matches distance matrix IDs
    metadata = metadata.loc[distance_matrix.ids]
    
    # Run PERMANOVA using adonis2 from vegan
    result = adonis2(formula=formula, data=metadata, permutations=permutations, 
                     method="bray", by="terms", parallel=1)
    
    return result

def run_db_rda(distance_matrix, metadata, formula, permutations=999):
    """
    Run distance-based Redundancy Analysis (db-RDA).
    
    Parameters:
    - distance_matrix: DistanceMatrix object
    - metadata: DataFrame with sample metadata
    - formula: Formula string (e.g., "~ Location + SampleType")
    - permutations: Number of permutations for significance testing
    
    Returns:
    - db-RDA results
    """
    # Perform Principal Coordinates Analysis (PCoA) on distance matrix
    pcoa_result = pcoa(distance_matrix)
    
    # Extract the PCoA axes that explain significant variation
    # Typically, we'd use Kaiser criterion (eigenvalue > 1) or broken stick model
    # For simplicity, we'll use axes that explain at least 1% of variation
    pcoa_coords = pcoa_result.samples.copy()
    
    # Calculate total explained variance
    total_explained = pcoa_result.proportion_explained.sum()
    
    # Filter to PCoA axes explaining at least 1% of variation
    significant_axes = pcoa_result.proportion_explained[pcoa_result.proportion_explained > 0.01]
    significant_indices = significant_axes.index
    
    print(f"Using {len(significant_axes)} PCoA axes that explain {significant_axes.sum():.1%} of variation")
    
    # Extract significant PCo axes for RDA
    pcoa_df = pd.DataFrame(
        pcoa_coords.loc[:, significant_indices],
        index=distance_matrix.ids
    )
    
    # Make sure metadata index matches PCoA IDs
    metadata = metadata.loc[pcoa_df.index]
    
    # Run RDA
    rda_result = rda(formula=formula, data=metadata, scale=True)
    
    return rda_result

def variance_partitioning(distance_matrix, metadata, variables, permutations=999):
    """
    Perform variance partitioning to determine unique and shared contributions.
    
    Parameters:
    - distance_matrix: DistanceMatrix object
    - metadata: DataFrame with sample metadata
    - variables: List of variables to partition
    - permutations: Number of permutations for significance testing
    
    Returns:
    - Variance partitioning results
    """
    # Perform Principal Coordinates Analysis (PCoA)
    pcoa_result = pcoa(distance_matrix)
    
    # Extract significant PCoA axes
    pcoa_coords = pcoa_result.samples.copy()
    significant_axes = pcoa_result.proportion_explained[pcoa_result.proportion_explained > 0.01]
    significant_indices = significant_axes.index
    
    # Extract significant PCo axes
    pcoa_df = pd.DataFrame(
        pcoa_coords.loc[:, significant_indices],
        index=distance_matrix.ids
    )
    
    # For each variable, create a formula
    variable_matrices = {}
    for var in variables:
        # Make sure metadata index matches PCoA IDs
        var_data = metadata.loc[pcoa_df.index, [var]]
        
        # Create design matrix for this variable
        if pd.api.types.is_categorical_dtype(var_data[var]) or var_data[var].dtype == 'object':
            # For categorical variables, create dummy variables
            var_dummies = pd.get_dummies(var_data[var], prefix=var)
            variable_matrices[var] = var_dummies
        else:
            # For numerical variables, just use as is
            variable_matrices[var] = var_data
    
    # Prepare data for varpart function
    y = pcoa_df  # Response (community composition)
    
    # Run variance partitioning
    if len(variables) == 2:
        result = varpart(y, variable_matrices[variables[0]], variable_matrices[variables[1]])
    elif len(variables) == 3:
        result = varpart(y, variable_matrices[variables[0]], 
                         variable_matrices[variables[1]], 
                         variable_matrices[variables[2]])
    elif len(variables) == 4:
        result = varpart(y, variable_matrices[variables[0]], 
                         variable_matrices[variables[1]], 
                         variable_matrices[variables[2]],
                         variable_matrices[variables[3]])
    else:
        raise ValueError("Variance partitioning supports 2-4 variables")
    
    return result

def plot_permanova_results(permanova_result, output_file=None):
    """
    Create a visualization of PERMANOVA results showing explained variance.
    
    Parameters:
    - permanova_result: PERMANOVA result from adonis2
    - output_file: Path to save the plot
    
    Returns:
    - matplotlib figure
    """
    # Extract R-squared values and p-values
    variables = permanova_result.index[:-1]  # Exclude 'Residual'
    r_squared = permanova_result.loc[variables, 'R2']
    p_values = permanova_result.loc[variables, 'Pr(>F)']
    
    # Sort by R-squared in descending order
    sorted_idx = r_squared.argsort()[::-1]
    variables = variables[sorted_idx]
    r_squared = r_squared[sorted_idx]
    p_values = p_values[sorted_idx]
    
    # Create bar plot
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(variables, r_squared * 100, color='skyblue')
    
    # Add significance markers
    for i, (var, p) in enumerate(zip(variables, p_values)):
        significance = ''
        if p < 0.001:
            significance = '***'
        elif p < 0.01:
            significance = '**'
        elif p < 0.05:
            significance = '*'
        
        if significance:
            ax.text(i, r_squared[i] * 100 + 0.5, significance, 
                   ha='center', va='bottom', fontweight='bold')
    
    # Add residual variance
    residual = permanova_result.loc['Residual', 'R2'] * 100
    ax.axhline(y=100 - residual, color='red', linestyle='--', 
               label=f'Explained variance ({100 - residual:.1f}%)')
    
    # Customize the plot
    ax.set_ylabel('Variance Explained (%)')
    ax.set_xlabel('Variables')
    ax.set_title('PERMANOVA Results: Explained Community Variance by Variable')
    ax.set_ylim(0, max(r_squared * 100) * 1.2)
    
    # Add values on top of bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
               f'{height:.1f}%', ha='center', va='bottom')
    
    # Add p-value legend
    ax.text(0.05, 0.05, '* p < 0.05, ** p < 0.01, *** p < 0.001',
           transform=ax.transAxes, fontsize=10, verticalalignment='bottom')
    
    plt.legend()
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    return fig

def plot_db_rda(rda_result, metadata, output_file=None):
    """
    Create triplot for db-RDA results.
    
    Parameters:
    - rda_result: RDA result
    - metadata: DataFrame with sample metadata
    - output_file: Path to save the plot
    
    Returns:
    - matplotlib figure
    """
    # Extract site scores, species scores, and biplot (constraint) scores
    site_scores = rda_result['sites']
    biplot_scores = rda_result['biplot']
    
    # Create biplot of first two RDA axes
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot site scores
    scatter = ax.scatter(site_scores[:, 0], site_scores[:, 1], c='blue', alpha=0.7)
    
    # Color by a categorical variable if available
    color_by = None
    for col in ['SampleType', 'Location', 'PostNatalAntibiotics', 'GestationCohort']:
        if col in metadata.columns:
            color_by = col
            break
    
    if color_by:
        # Color points by the selected categorical variable
        categories = metadata[color_by].astype('category')
        colors = dict(zip(categories.cat.categories, 
                          plt.cm.tab10(np.linspace(0, 1, len(categories.cat.categories)))))
        
        for category, color in colors.items():
            mask = metadata[color_by] == category
            ax.scatter(site_scores[mask, 0], site_scores[mask, 1], 
                      c=[color], label=category, alpha=0.7)
        
        plt.legend(title=color_by)
    
    # Plot biplot scores (environmental variables)
    biplot_scale = 0.5 * max(abs(site_scores[:, 0].max()), abs(site_scores[:, 1].max()))
    for i in range(biplot_scores.shape[0]):
        ax.arrow(0, 0, biplot_scores[i, 0] * biplot_scale, biplot_scores[i, 1] * biplot_scale,
                color='red', width=0.01, head_width=0.05, length_includes_head=True)
        ax.text(biplot_scores[i, 0] * biplot_scale * 1.1, 
                biplot_scores[i, 1] * biplot_scale * 1.1,
                rda_result['biplot_names'][i], color='darkred', fontweight='bold')
    
    # Add axis labels and title
    prop_explained_1 = rda_result['CCA$eig'][0] / sum(rda_result['tot.chi'])
    prop_explained_2 = rda_result['CCA$eig'][1] / sum(rda_result['tot.chi'])
    
    ax.set_xlabel(f'RDA1 ({prop_explained_1:.1%} explained)')
    ax.set_ylabel(f'RDA2 ({prop_explained_2:.1%} explained)')
    ax.set_title('Distance-based Redundancy Analysis (db-RDA)')
    
    # Add grid lines
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.7)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    return fig

def plot_variance_partitioning(varpart_result, variables, output_file=None):
    """
    Create a Venn diagram visualization of variance partitioning results.
    
    Parameters:
    - varpart_result: Variance partitioning result
    - variables: List of variables used in partitioning
    - output_file: Path to save the plot
    
    Returns:
    - matplotlib figure
    """
    from matplotlib_venn import venn2, venn3
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    if len(variables) == 2:
        # For 2 variables, extract the fractions
        # [a, b, ab]
        fractions = [
            varpart_result['part']['indfract.X1'],  # unique to X1
            varpart_result['part']['indfract.X2'],  # unique to X2
            varpart_result['part']['fract'].iloc[2]  # shared
        ]
        
        # Create Venn diagram
        venn2(subsets=(fractions[0], fractions[1], fractions[2]), 
              set_labels=variables, ax=ax)
        
    elif len(variables) == 3:
        # For 3 variables, extract the fractions
        # [a, b, c, ab, ac, bc, abc]
        fractions = [
            varpart_result['part']['indfract.X1'],  # unique to X1
            varpart_result['part']['indfract.X2'],  # unique to X2
            varpart_result['part']['indfract.X3'],  # unique to X3
            varpart_result['part']['fract'].iloc[3],  # X1 and X2
            varpart_result['part']['fract'].iloc[5],  # X1 and X3
            varpart_result['part']['fract'].iloc[6],  # X2 and X3
            varpart_result['part']['fract'].iloc[7]   # X1, X2, and X3
        ]
        
        # Create Venn diagram
        venn3(subsets=(fractions[0], fractions[1], fractions[3], 
                      fractions[2], fractions[4], fractions[5], fractions[6]),
              set_labels=variables, ax=ax)
    
    else:
        # For 4+ variables, create a bar plot instead
        components = []
        values = []
        
        # Individual fractions
        for i, var in enumerate(variables):
            components.append(f"{var} (unique)")
            values.append(varpart_result['part'][f'indfract.X{i+1}'])
        
        # Add joint fractions for simplicity
        if len(variables) == 4:
            joint_vars = [
                (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
                (0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3),
                (0, 1, 2, 3)
            ]
            
            for indices in joint_vars:
                var_names = [variables[i] for i in indices]
                components.append(" + ".join(var_names))
                # This will need customization depending on the varpart output structure
                # Below is a simplified approach
                values.append(0.01)  # Placeholder
        
        # Create bar plot
        ax.bar(range(len(components)), values)
        ax.set_xticks(range(len(components)))
        ax.set_xticklabels(components, rotation=90)
        ax.set_ylabel('Proportion of Variance Explained')
        ax.set_title('Variance Partitioning Results')
    
    plt.title('Variance Partitioning of Microbiome Composition')
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    return fig

def analyze_microbiome_variance(microbiome_df, metadata_df, 
                               distance_method="bray",
                               output_dir="results",
                               key_variables=None):
    """
    Comprehensive analysis of factors explaining microbiome variance.
    
    Parameters:
    - microbiome_df: DataFrame with microbiome abundance data
    - metadata_df: DataFrame with sample metadata
    - distance_method: Method for calculating distances
    - output_dir: Directory to save results
    - key_variables: List of key variables to include (if None, will detect automatically)
    
    Returns:
    - Dictionary with analysis results
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get common samples between microbiome and metadata
    common_samples = list(set(microbiome_df.index) & set(metadata_df.index))
    microbiome_df = microbiome_df.loc[common_samples]
    metadata_df = metadata_df.loc[common_samples]
    
    print(f"Using {len(common_samples)} samples with both microbiome and metadata")
    
    # Normalize microbiome data
    normalized_df, norm_method = select_best_normalization(microbiome_df)
    print(f"Applied {norm_method} normalization to microbiome data")
    
    # Save normalized data
    normalized_df.to_csv(os.path.join(output_dir, f"microbiome_normalized_{norm_method}.csv"))
    
    # Calculate distance matrix
    binary = norm_method == "presence_absence"
    dist_matrix = calculate_distance_matrix(normalized_df, method=distance_method, binary=binary)
    print(f"Calculated {distance_method} distances between samples")
    
    # Identify potentially important variables if not specified
    if key_variables is None:
        key_variables = []
        
        for col in metadata_df.columns:
            # Skip columns with too many missing values
            if metadata_df[col].isnull().mean() > 0.2:
                continue
                
            # Skip columns with unique values for each sample
            if metadata_df[col].nunique() == len(metadata_df):
                continue
                
            # Skip columns with only one value
            if metadata_df[col].nunique() <= 1:
                continue
                
            # For categorical variables, ensure enough samples per category
            if pd.api.types.is_categorical_dtype(metadata_df[col]) or metadata_df[col].dtype == 'object':
                if metadata_df[col].value_counts().min() < 3:
                    continue
            
            key_variables.append(col)
        
        print(f"Identified {len(key_variables)} potentially important variables: {', '.join(key_variables)}")
    
    # Create formula for PERMANOVA
    formula = " + ".join([f"{var}" for var in key_variables])
    formula = f"~{formula}"
    print(f"Using formula: {formula}")
    
    # Run PERMANOVA
    try:
        permanova_result = run_permanova(dist_matrix, metadata_df, formula)
        print("PERMANOVA Results:")
        print(permanova_result)
        
        # Save PERMANOVA results
        with open(os.path.join(output_dir, "permanova_result.txt"), "w") as f:
            f.write(str(permanova_result))
        
        # Plot PERMANOVA results
        fig = plot_permanova_results(permanova_result, 
                                    output_file=os.path.join(output_dir, "permanova_plot.png"))
    except Exception as e:
        print(f"Error running PERMANOVA: {e}")
        permanova_result = None
    
    # Run db-RDA
    try:
        db_rda_result = run_db_rda(dist_matrix, metadata_df, formula)
        print("db-RDA Results:")
        print(db_rda_result)
        
        # Save db-RDA results
        with open(os.path.join(output_dir, "db_rda_result.txt"), "w") as f:
            f.write(str(db_rda_result))
        
        # Plot db-RDA results
        fig = plot_db_rda(db_rda_result, metadata_df, 
                         output_file=os.path.join(output_dir, "db_rda_plot.png"))
    except Exception as e:
        print(f"Error running db-RDA: {e}")
        db_rda_result = None
    
    # Run variance partitioning for top variables
    if permanova_result is not None:
        # Get top variables by R-squared
        variables = permanova_result.index[:-1]  # Exclude 'Residual'
        r_squared = permanova_result.loc[variables, 'R2']
        p_values = permanova_result.loc[variables, 'Pr(>F)']
        
        # Filter to significant variables with highest R-squared
        sig_vars = [var for var, p in zip(variables, p_values) if p < 0.05]
        
        if len(sig_vars) >= 2:
            # Limit to top 3-4 variables for variance partitioning
            top_vars = sig_vars[:min(4, len(sig_vars))]
            
            try:
                varpart_result = variance_partitioning(dist_matrix, metadata_df, top_vars)
                print("Variance Partitioning Results:")
                print(varpart_result)
                
                # Save variance partitioning results
                with open(os.path.join(output_dir, "varpart_result.txt"), "w") as f:
                    f.write(str(varpart_result))
                
                # Plot variance partitioning
                if len(top_vars) <= 3:  # Venn diagram only works for 2-3 variables
                    fig = plot_variance_partitioning(varpart_result, top_vars,
                                                   output_file=os.path.join(output_dir, "varpart_plot.png"))
            except Exception as e:
                print(f"Error running variance partitioning: {e}")
                varpart_result = None
        else:
            print("Not enough significant variables for variance partitioning")
            varpart_result = None
    else:
        varpart_result = None
    
    # Return results
    results = {
        "normalized_df": normalized_df,
        "normalization_method": norm_method,
        "distance_matrix": dist_matrix,
        "permanova_result": permanova_result,
        "db_rda_result": db_rda_result,
        "varpart_result": varpart_result
    }
    
    # Save results object
    with open(os.path.join(output_dir, "variance_analysis_results.pkl"), "wb") as f:
        pickle.dump(results, f)
    
    return results

if __name__ == "__main__":
    # Create output directory
    os.makedirs("results", exist_ok=True)
    
    print("Loading microbiome data...")
    # Load microbiome count data
    microbiome_df = pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)
    
    print("Loading metadata...")
    # Load clinical metadata
    metadata_df = pd.read_csv("../metadata/AllNICUSampleKey20250206.csv", index_col=0)
    
    # Specify key variables based on domain knowledge
    key_variables = [
        "Location", "SampleType", "SampleCollectionWeek", 
        "GestationCohort", "PostNatalAntibiotics", 
        "MaternalAntibiotics", "AnyMilk"
    ]
    
    # Run comprehensive analysis
    results = analyze_microbiome_variance(
        microbiome_df=microbiome_df,
        metadata_df=metadata_df,
        key_variables=key_variables,
        output_dir="results/variance_analysis"
    )
    
    print("Analysis complete! Results saved to 'results/variance_analysis' directory.")
