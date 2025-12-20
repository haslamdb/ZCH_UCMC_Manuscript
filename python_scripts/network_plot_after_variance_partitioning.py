#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import pickle
import os
from sklearn.preprocessing import StandardScaler

# Set up matplotlib for publication-quality PDFs
# Use DejaVu Sans (similar to Arial) - Arial not available on this system
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42  # TrueType fonts (editable in Illustrator)
mpl.rcParams['ps.fonttype'] = 42

def load_lmm_results(file_path):
    """Load LMM results from a pickle file."""
    with open(file_path, 'rb') as f:
        results = pickle.load(f)
    return results

def create_interaction_network(lmm_results, output_file="clinical_factor_interactions_network.pdf", p_threshold=0.05):
    """
    Create a network diagram showing clinical factor interactions based on LMM results.

    Parameters:
    - lmm_results: Dictionary with LMM results for each principal component
    - output_file: Path to save the network plot
    - p_threshold: P-value threshold for including variables (default 0.05)
    """
    # Extract variables and their effect sizes (using max absolute coefficient per factor per PC)
    all_vars = {}

    for pc, result in lmm_results.items():
        coeffs = result['coefficients']
        filtered = coeffs[coeffs['P-value'] < p_threshold]

        # Track max effect size per variable group per PC
        pc_var_effects = {}

        for _, row in filtered.iterrows():
            var_name = row['Variable']

            # Skip intercept and group var
            if var_name in ['Intercept', 'Group Var']:
                continue

            # Extract variable group (e.g., "Location" from "C(Location)[T.Hangzhou]")
            if 'C(' in var_name:
                var_group = var_name.split(')')[0].replace('C(', '')
            else:
                var_group = var_name

            # Keep max absolute coefficient for this variable group in this PC
            abs_coef = abs(row['Coefficient'])
            if var_group not in pc_var_effects or abs_coef > pc_var_effects[var_group]:
                pc_var_effects[var_group] = abs_coef

        # Add to overall tracking
        for var_group, effect in pc_var_effects.items():
            if var_group not in all_vars:
                all_vars[var_group] = {'pcs': [], 'effects': [], 'total_effect': 0}

            all_vars[var_group]['pcs'].append(pc)
            all_vars[var_group]['effects'].append(effect)
            all_vars[var_group]['total_effect'] += effect

    # Create graph
    G = nx.Graph()

    # Add nodes with attributes - size based on total effect size
    max_effect = max(data['total_effect'] for data in all_vars.values()) if all_vars else 1
    for var, data in all_vars.items():
        # Normalize node size: min 600, max 2500
        normalized_size = 600 + (data['total_effect'] / max_effect) * 1900
        G.add_node(var,
                  pcs=','.join(data['pcs']),
                  total_effect=data['total_effect'],
                  size=normalized_size)

    # Add edges weighted by combined effect size
    for pc, result in lmm_results.items():
        coeffs = result['coefficients']
        filtered = coeffs[coeffs['P-value'] < p_threshold]

        # Get max effect per variable group for this PC
        pc_var_effects = {}
        for _, row in filtered.iterrows():
            if row['Variable'] in ['Intercept', 'Group Var']:
                continue

            if 'C(' in row['Variable']:
                var_group = row['Variable'].split(')')[0].replace('C(', '')
            else:
                var_group = row['Variable']

            abs_coef = abs(row['Coefficient'])
            if var_group not in pc_var_effects or abs_coef > pc_var_effects[var_group]:
                pc_var_effects[var_group] = abs_coef

        # Connect pairs with edge weight = geometric mean of their effects
        var_list = list(pc_var_effects.keys())
        for i, var1 in enumerate(var_list):
            for var2 in var_list[i+1:]:
                # Edge weight = geometric mean of the two effects
                edge_weight = np.sqrt(pc_var_effects[var1] * pc_var_effects[var2])

                if G.has_edge(var1, var2):
                    G[var1][var2]['weight'] += edge_weight
                    G[var1][var2]['pcs'] += f",{pc}"
                else:
                    G.add_edge(var1, var2, weight=edge_weight, pcs=pc)

    # Create layout
    pos = nx.spring_layout(G, k=0.8, iterations=200, seed=42)
    
    # Draw network
    plt.figure(figsize=(14, 12))

    # Get node sizes
    node_sizes = [G.nodes[n]['size'] for n in G.nodes()]

    # Draw nodes with varying size based on effect size
    nx.draw_networkx_nodes(G, pos,
                          node_size=node_sizes,
                          node_color='lightblue',
                          edgecolors='darkblue',
                          linewidths=2,
                          alpha=0.8)

    # Normalize edge widths for visualization (min 1.5, max 12)
    if G.edges():
        weights = [G[u][v]['weight'] for u, v in G.edges()]
        max_weight = max(weights) if weights else 1
        min_weight = min(weights) if weights else 0
        weight_range = max_weight - min_weight if max_weight != min_weight else 1
        edge_widths = [1.5 + (G[u][v]['weight'] - min_weight) / weight_range * 10.5 for u, v in G.edges()]
    else:
        edge_widths = []

    # Draw edges with varying thickness based on effect size
    nx.draw_networkx_edges(G, pos,
                          width=edge_widths,
                          alpha=0.5,
                          edge_color='gray')

    # Add node labels
    nx.draw_networkx_labels(G, pos,
                          font_size=11,
                          font_weight='bold')

    # Edge labels removed for cleaner visualization
    # edge_labels = {(u,v): f"{d['weight']:.1f}" for u, v, d in G.edges(data=True)}
    # nx.draw_networkx_edge_labels(G, pos,
    #                             edge_labels=edge_labels,
    #                             font_size=8)

    plt.title('Clinical Factor Interaction Network\n(Edge weights = combined effect sizes from LMM)',
             fontsize=16, pad=20)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    return G

def create_heatmap_from_lmm(lmm_results, output_file="clinical_factor_heatmap.pdf", p_threshold=0.05):
    """
    Create a heatmap showing the significance of clinical factors across principal components.

    Parameters:
    - lmm_results: Dictionary with LMM results for each principal component
    - output_file: Path to save the heatmap plot
    - p_threshold: P-value threshold for including variables (default 0.05)
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
            
            # Use the absolute value of the coefficient if p-value below threshold, 0 otherwise
            if row['P-value'] < p_threshold:
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
        
        # Create network visualization (p_threshold=1.0 includes all factors)
        network = create_interaction_network(
            lmm_results=lmm_results,
            output_file="results/network_plots/clinical_factor_interactions_network.pdf",
            p_threshold=1.0
        )

        # Create heatmap visualization
        heatmap_data = create_heatmap_from_lmm(
            lmm_results=lmm_results,
            output_file="results/network_plots/clinical_factor_heatmap.pdf",
            p_threshold=1.0
        )
    except Exception as e:
        print(f"Error using variance_analysis_results.pkl: {e}")
        print("The variance_analysis_results.pkl file has a different structure than expected.")
        print("This script requires the lmm_microbiome_community_results.pkl file with specific keys.")
        print("Available keys in variance_analysis_results.pkl: normalized_df, normalization_method, distance_matrix, permanova_result, db_rda_result, varpart_result")
        print("\nPlease use the original lmm_microbiome_community_results.pkl file or modify the script to handle the new structure.")
    
    print("Network visualization completed!")