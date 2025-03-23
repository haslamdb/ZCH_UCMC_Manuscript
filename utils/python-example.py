#!/usr/bin/env python
# Example Python script using the unified configuration

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.config import load_config

def main():
    # Load the configuration
    config = load_config()
    
    # Access configuration values
    print(f"Loading data from: {config['input_files']['microbiome_counts']}")
    print(f"Saving results to: {config['paths']['results_dir']}")
    print(f"Key organisms to analyze: {', '.join(config['key_organisms'])}")
    
    # Set random seed for reproducibility
    np.random.seed(config['statistics']['random_seed'])
    
    # Load data
    try:
        microbiome_df = pd.read_csv(config['input_files']['microbiome_counts'], index_col=0)
        print(f"Loaded microbiome data with shape: {microbiome_df.shape}")
        
        metadata_df = pd.read_csv(config['input_files']['sample_key'], index_col=0)
        print(f"Loaded metadata with shape: {metadata_df.shape}")
        
        # Show a summary of key organisms
        key_organisms = config['key_organisms']
        organisms_found = [org for org in key_organisms if org in microbiome_df.columns]
        
        print(f"\nFound {len(organisms_found)}/{len(key_organisms)} key organisms in the dataset:")
        for organism in organisms_found:
            non_zero = (microbiome_df[organism] > 0).sum()
            percent = non_zero / len(microbiome_df) * 100
            print(f"  - {organism}: Present in {non_zero} samples ({percent:.1f}%)")
        
        # Create a simple visualization using the config's visualization parameters
        fig_size = config['visualization']['figure_sizes']['default']
        dpi = config['visualization']['dpi']
        color_palette = config['visualization']['color_palettes']['categorical']
        
        # Example: Create a bar plot of organism prevalence
        prevalence = [(microbiome_df[org] > 0).mean() * 100 for org in organisms_found]
        
        plt.figure(figsize=fig_size)
        bars = plt.bar(range(len(organisms_found)), prevalence, color=color_palette[:len(organisms_found)])
        plt.xticks(range(len(organisms_found)), organisms_found, rotation=45, ha='right')
        plt.ylabel('Prevalence (%)')
        plt.title('Prevalence of Key Organisms')
        plt.tight_layout()
        
        # Save the figure to the results directory
        figures_dir = config['paths']['figures_dir']
        os.makedirs(figures_dir, exist_ok=True)
        plt.savefig(os.path.join(figures_dir, 'organism_prevalence.png'), dpi=dpi)
        plt.close()
        
        print(f"\nCreated visualization: {os.path.join(figures_dir, 'organism_prevalence.png')}")
        
    except Exception as e:
        print(f"Error processing data: {str(e)}")

if __name__ == "__main__":
    main()
