
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
from sklearn.manifold import TSNE

# Load microbiome count data
microbiome_df = pd.read_csv("../data/NICUSpeciesReduced.csv", index_col=0)

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
microbiome_clr = mt.clr_transform(microbiome_df)

# Merge microbiome data and metadata on Sample ID
# here we're using clr transformed data
data = microbiome_clr.merge(metadata_df, left_index=True, right_index=True)

print(f"Data Shape After Merge: {data.shape}")

# List of microbes to analyze
key_organisms = ["Klebsiella.pneumoniae", "Staphylococcus.aureus", "Escherichia.coli", "Klebsiella.oxytoca",
                 "Staphylococcus.epidermidis", "Streptococcus.pyogenes", "Staphylococcs.capitus", 
                 "Enterococcus.faecium", "Enterococcus.faecalis", "Serratia.marcescens", "Listeria monocytogenes"]


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import seaborn as sns

# Assuming your data loading code is already done:
# microbiome_df, microbiome_clr, metadata_df, and data are already defined
# key_organisms list is also defined

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable

# Assuming your data loading code is already done:
# microbiome_df, microbiome_clr, metadata_df, and data are already defined
# key_organisms list is also defined

# Define a color for each organism
def get_organism_colors(key_organisms):
    """Generate a distinct color for each organism"""
    # Define a set of distinct colors
    base_colors = [
        '#e41a1c',  # red
        '#377eb8',  # blue
        '#4daf4a',  # green
        '#984ea3',  # purple
        '#ff7f00',  # orange
        '#ffff33',  # yellow
        '#a65628',  # brown
        '#f781bf',  # pink
        '#999999',  # grey
        '#8dd3c7',  # teal
        '#bebada'   # light purple
    ]
    
    # If we have more organisms than colors, we'll cycle through the colors
    organism_colors = {}
    for i, org in enumerate(key_organisms):
        organism_colors[org] = base_colors[i % len(base_colors)]
    
    return organism_colors

# Create a custom colormap for each organism (light to dark gradient)
def create_custom_cmap(base_color):
    """Create a custom colormap that goes from white to the base color"""
    return mcolors.LinearSegmentedColormap.from_list('custom_cmap', ['#ffffff', base_color])

# 1. Generate t-SNE on the CLR-transformed microbiome data
def generate_tsne_plot(microbiome_data, key_organisms, perplexity=30, n_iter=1000):
    # Perform t-SNE
    tsne = TSNE(n_components=2, perplexity=perplexity, n_iter=n_iter, random_state=42)
    tsne_results = tsne.fit_transform(microbiome_data)
    
    # Create a DataFrame with t-SNE results
    tsne_df = pd.DataFrame(data=tsne_results, columns=['t-SNE1', 't-SNE2'], index=microbiome_data.index)
    
    # Get colors for each organism
    organism_colors = get_organism_colors(key_organisms)
    
    # Create plots for each key organism
    for organism in key_organisms:
        if organism in microbiome_data.columns:
            # Create figure
            plt.figure(figsize=(10, 8))
            
            # Get base color for this organism
            base_color = organism_colors[organism]
            
            # Create custom colormap from white to the base color
            custom_cmap = create_custom_cmap(base_color)
            
            # Create a scatter plot with t-SNE coordinates
            scatter = plt.scatter(
                tsne_df['t-SNE1'], 
                tsne_df['t-SNE2'],
                c=microbiome_data[organism],
                cmap=custom_cmap,
                alpha=0.9,
                s=50
            )
            
            # Add colorbar
            cbar = plt.colorbar(scatter)
            cbar.set_label(f'{organism} abundance (CLR)')
            
            # Add labels and title
            plt.xlabel('t-SNE Component 1')
            plt.ylabel('t-SNE Component 2')
            plt.title(f't-SNE plot colored by {organism} abundance')
            
            # Add grid
            plt.grid(alpha=0.3)
            
            # Save as PNG
            plt.tight_layout()
            plt.savefig(f'tsne_{organism.replace(".", "_")}.png', dpi=300)
            
            # Save as PDF
            plt.savefig(f'tsne_{organism.replace(".", "_")}.pdf')
            
            plt.show()
        else:
            print(f"Warning: {organism} not found in the dataset")
    
    # Also create a single plot with multiple subplots for comparison
    create_multi_organism_plot(tsne_df, microbiome_data, key_organisms, organism_colors)
    
    return tsne_df

def create_multi_organism_plot(tsne_df, microbiome_data, key_organisms, organism_colors, max_cols=3):
    # Filter to only include organisms that exist in the dataset
    valid_organisms = [org for org in key_organisms if org in microbiome_data.columns]
    
    if not valid_organisms:
        print("None of the specified key organisms found in the dataset")
        return
    
    # Calculate grid dimensions
    n_organisms = len(valid_organisms)
    n_cols = min(max_cols, n_organisms)
    n_rows = (n_organisms + n_cols - 1) // n_cols
    
    # Create figure and axes
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*5, n_rows*4))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Plot each organism
    for i, organism in enumerate(valid_organisms):
        # Get base color for this organism
        base_color = organism_colors[organism]
        
        # Create custom colormap from white to the base color
        custom_cmap = create_custom_cmap(base_color)
        
        scatter = axes[i].scatter(
            tsne_df['t-SNE1'], 
            tsne_df['t-SNE2'],
            c=microbiome_data[organism],
            cmap=custom_cmap,
            alpha=0.9,
            s=30
        )
        
        # Add colorbar
        cbar = fig.colorbar(scatter, ax=axes[i])
        cbar.set_label(f'Abundance (CLR)')
        
        # Add labels and title
        axes[i].set_xlabel('t-SNE1')
        axes[i].set_ylabel('t-SNE2')
        axes[i].set_title(organism)
        axes[i].grid(alpha=0.3)
    
    # Hide any unused subplots
    for j in range(i+1, len(axes)):
        axes[j].set_visible(False)
    
    plt.tight_layout()
    
    # Save as PNG
    plt.savefig('tsne_all_key_organisms.png', dpi=300)
    
    # Save as PDF
    plt.savefig('tsne_all_key_organisms.pdf')
    
    # Save all individual plots to a single multi-page PDF
    with PdfPages('all_key_organisms_multipage.pdf') as pdf:
        for i, organism in enumerate(valid_organisms):
            # Create a new figure for each organism
            fig_single = plt.figure(figsize=(10, 8))
            
            # Get base color for this organism
            base_color = organism_colors[organism]
            
            # Create custom colormap from white to the base color
            custom_cmap = create_custom_cmap(base_color)
            
            scatter = plt.scatter(
                tsne_df['t-SNE1'], 
                tsne_df['t-SNE2'],
                c=microbiome_data[organism],
                cmap=custom_cmap,
                alpha=0.9,
                s=50
            )
            
            cbar = plt.colorbar(scatter)
            cbar.set_label(f'{organism} abundance (CLR)')
            
            plt.xlabel('t-SNE Component 1')
            plt.ylabel('t-SNE Component 2')
            plt.title(f't-SNE plot colored by {organism} abundance')
            plt.grid(alpha=0.3)
            
            plt.tight_layout()
            pdf.savefig(fig_single)
            plt.close(fig_single)
        
        # Also add the composite figure
        pdf.savefig(fig)
    
    plt.show()

# Function for generating t-SNE plots colored by metadata variables
def tsne_with_metadata(microbiome_data, metadata, subject_id_col, categorical_vars=None, numerical_vars=None, perplexity=30, n_iter=1000):
    """
    Generate t-SNE plot with metadata overlays
    
    Args:
        microbiome_data: CLR-transformed microbiome data
        metadata: Metadata DataFrame
        subject_id_col: Column name for subject ID in metadata
        categorical_vars: List of categorical variables to plot
        numerical_vars: List of numerical variables to plot
        perplexity: t-SNE perplexity parameter
        n_iter: Number of iterations for t-SNE
    """
    # Perform t-SNE
    tsne = TSNE(n_components=2, perplexity=perplexity, n_iter=n_iter, random_state=42)
    tsne_results = tsne.fit_transform(microbiome_data)
    
    # Create a DataFrame with t-SNE results
    tsne_df = pd.DataFrame(data=tsne_results, columns=['t-SNE1', 't-SNE2'], index=microbiome_data.index)
    
    # Merge with metadata
    tsne_with_meta = tsne_df.merge(metadata, left_index=True, right_index=True)
    
    # Create multi-page PDF for all metadata visualizations
    with PdfPages('tsne_metadata_visualizations.pdf') as pdf:
        # If categorical variables are provided, create plots for them
        if categorical_vars:
            for var in categorical_vars:
                if var in tsne_with_meta.columns:
                    # Create figure
                    fig = plt.figure(figsize=(12, 10))
                    
                    # Create a categorical plot
                    sns.scatterplot(
                        data=tsne_with_meta,
                        x='t-SNE1',
                        y='t-SNE2',
                        hue=var,
                        palette='Set1',
                        s=100,
                        alpha=0.8
                    )
                    
                    plt.title(f't-SNE plot colored by {var}')
                    plt.xlabel('t-SNE Component 1')
                    plt.ylabel('t-SNE Component 2')
                    plt.grid(alpha=0.3)
                    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                    
                    plt.tight_layout()
                    
                    # Save as PNG
                    plt.savefig(f'tsne_{var}.png', dpi=300)
                    
                    # Save as individual PDF
                    plt.savefig(f'tsne_{var}.pdf')
                    
                    # Add to multi-page PDF
                    pdf.savefig(fig)
                    
                    plt.show()
                    plt.close(fig)
                else:
                    print(f"Warning: {var} not found in the metadata")
        
        # If numerical variables are provided, create plots for them
        if numerical_vars:
            for var in numerical_vars:
                if var in tsne_with_meta.columns:
                    # Create figure
                    fig = plt.figure(figsize=(10, 8))
                    
                    # Create a numerical color-coded plot
                    scatter = plt.scatter(
                        tsne_with_meta['t-SNE1'],
                        tsne_with_meta['t-SNE2'],
                        c=tsne_with_meta[var],
                        cmap='viridis',
                        alpha=0.8,
                        s=80
                    )
                    
                    # Add colorbar
                    cbar = plt.colorbar(scatter)
                    cbar.set_label(var)
                    
                    plt.title(f't-SNE plot colored by {var}')
                    plt.xlabel('t-SNE Component 1')
                    plt.ylabel('t-SNE Component 2')
                    plt.grid(alpha=0.3)
                    
                    plt.tight_layout()
                    
                    # Save as PNG
                    plt.savefig(f'tsne_{var}.png', dpi=300)
                    
                    # Save as individual PDF
                    plt.savefig(f'tsne_{var}.pdf')
                    
                    # Add to multi-page PDF
                    pdf.savefig(fig)
                    
                    plt.show()
                    plt.close(fig)
                else:
                    print(f"Warning: {var} not found in the metadata")
    
    # Create a combined metadata visualization PDF with all variables
    create_combined_metadata_pdf(tsne_with_meta, categorical_vars, numerical_vars)
    
    return tsne_with_meta

def create_combined_metadata_pdf(tsne_with_meta, categorical_vars=None, numerical_vars=None):
    """
    Create a single PDF with all metadata visualizations arranged in a grid
    """
    # Count total variables to visualize
    total_vars = 0
    if categorical_vars:
        total_vars += sum(1 for var in categorical_vars if var in tsne_with_meta.columns)
    if numerical_vars:
        total_vars += sum(1 for var in numerical_vars if var in tsne_with_meta.columns)
    
    if total_vars == 0:
        print("No valid metadata variables to visualize")
        return
    
    # Calculate grid dimensions
    max_cols = 2
    n_cols = min(max_cols, total_vars)
    n_rows = (total_vars + n_cols - 1) // n_cols
    
    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*6, n_rows*5))
    
    # Handle single subplot case
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Track subplot index
    plot_idx = 0
    
    # Plot categorical variables
    if categorical_vars:
        for var in categorical_vars:
            if var in tsne_with_meta.columns and plot_idx < len(axes):
                # Use Seaborn for categorical plots
                sns.scatterplot(
                    data=tsne_with_meta,
                    x='t-SNE1',
                    y='t-SNE2',
                    hue=var,
                    palette='Set1',
                    s=80,
                    alpha=0.8,
                    ax=axes[plot_idx]
                )
                
                axes[plot_idx].set_title(f'{var}')
                axes[plot_idx].set_xlabel('t-SNE1')
                axes[plot_idx].set_ylabel('t-SNE2')
                axes[plot_idx].grid(alpha=0.3)
                
                # Move legend outside for better visibility
                axes[plot_idx].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                
                plot_idx += 1
    
    # Plot numerical variables
    if numerical_vars:
        for var in numerical_vars:
            if var in tsne_with_meta.columns and plot_idx < len(axes):
                # Create numerical color-coded plot
                scatter = axes[plot_idx].scatter(
                    tsne_with_meta['t-SNE1'],
                    tsne_with_meta['t-SNE2'],
                    c=tsne_with_meta[var],
                    cmap='viridis',
                    alpha=0.8,
                    s=80
                )
                
                # Add colorbar
                cbar = fig.colorbar(scatter, ax=axes[plot_idx])
                cbar.set_label(var)
                
                axes[plot_idx].set_title(f'{var}')
                axes[plot_idx].set_xlabel('t-SNE1')
                axes[plot_idx].set_ylabel('t-SNE2')
                axes[plot_idx].grid(alpha=0.3)
                
                plot_idx += 1
    
    # Hide any unused subplots
    for i in range(plot_idx, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('tsne_all_metadata_variables.pdf')
    plt.savefig('tsne_all_metadata_variables.png', dpi=300)
    plt.show()

# Execute the functions
# First generate t-SNE plot colored by organism abundance
tsne_df = generate_tsne_plot(microbiome_clr, key_organisms)

# Then generate t-SNE plots colored by the requested metadata variables
categorical_vars = ['SampleType', 'Location', 'SampleCollectionWeek']

# Run the metadata visualization
tsne_with_meta = tsne_with_metadata(
    microbiome_clr, 
    metadata_df, 
    subject_id_col, 
    categorical_vars=categorical_vars, 
    numerical_vars=numerical_vars
)

# Optional: create a unified function that runs all visualizations
def run_all_visualizations(microbiome_data, metadata, subject_id_col, key_organisms, 
                          categorical_vars=None, numerical_vars=None, perplexity=30):
    """
    Run all t-SNE visualizations in one function
    """
    # Generate organism abundance plots
    tsne_df = generate_tsne_plot(microbiome_data, key_organisms, perplexity=perplexity)
    
    # Generate metadata visualization plots
    tsne_with_meta = tsne_with_metadata(
        microbiome_data, 
        metadata, 
        subject_id_col, 
        categorical_vars=categorical_vars,
        numerical_vars=numerical_vars,
        perplexity=perplexity
    )
    
    return tsne_df, tsne_with_meta

# Example usage of the unified function
# tsne_df, tsne_with_meta = run_all_visualizations(
#     microbiome_clr,
#     metadata_df,
#     subject_id_col,
#     key_organisms,
#     categorical_vars=['SampleType', 'Location'],
#     numerical_vars=['SampleCollection_week']
# )