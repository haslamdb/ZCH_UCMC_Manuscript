#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pickle
import os
from sklearn.preprocessing import StandardScaler

def load_lmm_results(file_path):
    """Load LMM results from a pickle file."""
    with open(file_path, 'rb') as f:
        results = pickle.load(f)
    return results

def create_interaction_network(lmm_results, output_file="clinical_factor_interactions_network.png", threshold=0.01):
    """
    Create a network diagram showing clinical factor interactions based on LMM results.
    
    Parameters:
    - lmm_results: Dictionary with LMM results for each principal component
    - output_file: Path to save the network plot
    - threshold: Significance threshold for including variables
    """
    # Extract significant variables and their effect sizes
    all_significant_vars = {}
    
    for pc, result in lmm_results.items():
        coeffs = result['coefficients']
        significant = coeffs[coeffs['Significant']]
        
        # Extract unique variable groups (e.g., "Location" from "C(Location)[T.Hangzhou]")
        for _, row in significant.iterrows():
            var_name = row['Variable']
            
            # Skip intercept and group var
            if var_name in ['Intercept', 'Group Var']:
                continue
                
            # Extract variable group using regex
            if 'C(' in var_name:
                var_group = var_name.split(')')[0].replace('C(', '')
            else:
                var_group = var_name
                
            # Add to dictionary with absolute coefficient as weight
            if var_group not in all_significant_vars:
                all_significant_vars[var_group] = {'count': 0, 'pcs': [], 'total_effect': 0}
            
            all_significant_vars[var_group]['count'] += 1
            all_significant_vars[var_group]['pcs'].append(pc)
            all_significant_vars[var_group]['total_effect'] += abs(row['Coefficient'])
    
    # Create graph
    G = nx.Graph()
    
    # Add nodes with attributes
    for var, data in all_significant_vars.items():
        G.add_node(var, 
                  count=data['count'],
                  pcs=','.join(data['pcs']),
                  effect=data['total_effect'],
                  size=100 + (data['count'] * 100))  # Node size based on occurrence
    
    # Add edges between variables that appear in the same PC
    for pc, result in lmm_results.items():
        coeffs = result['coefficients']
        significant = coeffs[coeffs['Significant']]
        
        # Get significant variable groups for this PC
        pc_vars = []
        for _, row in significant.iterrows():
            if row['Variable'] in ['Intercept', 'Group Var']:
                continue
                
            if 'C(' in row['Variable']:
                var_group = row['Variable'].split(')')[0].replace('C(', '')
            else:
                var_group = row['Variable']
                
            pc_vars.append(var_group)
        
        # Connect all pairs of variables that co-occur in this PC
        for i, var1 in enumerate(pc_vars):
            for var2 in pc_vars[i+1:]:
                if G.has_edge(var1, var2):
                    # Increment weight
                    G[var1][var2]['weight'] += 1
                    G[var1][var2]['pcs'] += f",{pc}"
                else:
                    # Create new edge
                    G.add_edge(var1, var2, 
                              weight=1, 
                              pcs=pc)
    
    # Create layout
    pos = nx.spring_layout(G, k=0.5, iterations=200, seed=42)
    
    # Draw network
    plt.figure(figsize=(14, 12))
    
    # Get node sizes
    node_sizes = [G.nodes[n]['size'] for n in G.nodes()]
    
    # Draw nodes with varying size based on occurrence
    nx.draw_networkx_nodes(G, pos, 
                          node_size=node_sizes,
                          node_color='lightblue',
                          edgecolors='darkblue',
                          linewidths=2,
                          alpha=0.8)
    
    # Draw edges with varying thickness based on weight
    edge_weights = [G[u][v]['weight'] * 2 for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, 
                          width=edge_weights,
                          alpha=0.5,
                          edge_color='gray')
    
    # Add node labels
    nx.draw_networkx_labels(G, pos, 
                          font_size=12, 
                          font_weight='bold')
    
    # Add edge labels (shared PCs)
    edge_labels = {(u,v): f"{d['weight']}" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, 
                                edge_labels=edge_labels,
                                font_size=10)
    
    plt.title('Clinical Factor Interaction Network\n(Based on LMM Analysis of Microbiome Principal Components)',
             fontsize=16, pad=20)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    return G

def create_heatmap_from_lmm(lmm_results, output_file="clinical_factor_heatmap.png"):
    """
    Create a heatmap showing the significance of clinical factors across principal components.
    
    Parameters:
    - lmm_results: Dictionary with LMM results for each principal component
    - output_file: Path to save the heatmap plot
    """
    import seaborn as sns
    
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
    
    # Fill in coefficient values
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
    
    # Create heatmap
    plt.figure(figsize=(10, 12))
    sns.heatmap(heatmap_data, 
               cmap='viridis',
               annot=True,
               fmt='.2f',
               cbar_kws={'label': 'Absolute Coefficient Value (if significant)'})
    
    plt.title('Clinical Factors Influencing Microbiome Principal Components')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return heatmap_data

if __name__ == "__main__":
    # Create output directory
    os.makedirs("results/network_plots", exist_ok=True)
    
    # Load LMM results - fallback to original file if needed
    # Try the new file first
    try:
        results_file = "results/variance_analysis/lmm_microbiome_community_results.pkl"
        lmm_results = load_lmm_results(results_file)
        
        # Create network visualization
        network = create_interaction_network(
            lmm_results=lmm_results,
            output_file="results/network_plots/clinical_factor_interactions_network.png",
            threshold=0.01
        )
        
        # Create heatmap visualization
        heatmap_data = create_heatmap_from_lmm(
            lmm_results=lmm_results,
            output_file="results/network_plots/clinical_factor_heatmap.png"
        )
    except Exception as e:
        print(f"Error using variance_analysis_results.pkl: {e}")
        print("The variance_analysis_results.pkl file has a different structure than expected.")
        print("This script requires the lmm_microbiome_community_results.pkl file with specific keys.")
        print("Available keys in variance_analysis_results.pkl: normalized_df, normalization_method, distance_matrix, permanova_result, db_rda_result, varpart_result")
        print("\nPlease use the original lmm_microbiome_community_results.pkl file or modify the script to handle the new structure.")
    
    print("Network visualization completed!")