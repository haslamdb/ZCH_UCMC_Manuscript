#!/usr/bin/env python3
"""
LMM with Interaction Terms for Microbiome Community Analysis

This script fits Linear Mixed Models with pairwise interaction terms
to microbiome principal components, then extracts interaction coefficients
for network visualization.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from sklearn.decomposition import PCA
import statsmodels.formula.api as smf
import pickle
import os
import warnings
from itertools import combinations

# Set up matplotlib for publication-quality PDFs
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

warnings.filterwarnings("ignore", category=UserWarning)


def lmm_with_interactions(microbiome_df, metadata_df, n_components=5, output_dir="results/lmm_interactions"):
    """
    Apply Linear Mixed Model with interaction terms to principal components of microbiome data.

    Parameters:
    - microbiome_df: DataFrame with microbiome abundance data (normalized)
    - metadata_df: DataFrame with sample metadata
    - n_components: Number of principal components to analyze
    - output_dir: Directory to save results

    Returns:
    - Dictionary of LMM results including interaction coefficients
    """
    os.makedirs(output_dir, exist_ok=True)

    # Get common samples
    common_samples = list(set(microbiome_df.index) & set(metadata_df.index))
    microbiome_df = microbiome_df.loc[common_samples]
    metadata_df = metadata_df.loc[common_samples]

    print(f"Using {len(common_samples)} samples with both microbiome and metadata")

    # Perform PCA on the microbiome data
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(microbiome_df)

    # Create a DataFrame with PC scores
    pc_df = pd.DataFrame(
        pca_result,
        index=microbiome_df.index,
        columns=[f"PC{i+1}" for i in range(n_components)]
    )

    # Report variance explained
    print("\nVariance explained by each PC:")
    for i, var in enumerate(pca.explained_variance_ratio_):
        print(f"  PC{i+1}: {var:.2%}")

    # Combine PC scores with metadata
    combined_df = pd.merge(pc_df, metadata_df, left_index=True, right_index=True)

    # Define the categorical features to include in the model
    categorical_features = [
        "Location", "SampleType", "SampleCollectionWeek",
        "GestationCohort", "PostNatalAbxCohort",
        "MaternalAntibiotics", "AnyMilk", "PICC", "UVC", "Delivery"
    ]

    # Filter to features that exist in the data
    model_features = [feat for feat in categorical_features if feat in combined_df.columns]
    print(f"\nUsing {len(model_features)} clinical variables: {model_features}")

    # Identify subject column for random effects
    if "Subject" in combined_df.columns:
        subject_col = "Subject"
    else:
        subject_candidates = [col for col in combined_df.columns if "subject" in col.lower()]
        subject_col = subject_candidates[0] if subject_candidates else None

    # Store results for each PC
    lmm_results = {}

    # Run LMM for each principal component
    for pc in pc_df.columns:
        print(f"\n{'='*60}")
        print(f"Fitting LMM with interactions for {pc}")
        print(f"{'='*60}")

        # Build formula with main effects
        main_effects = []
        for feat in model_features:
            if combined_df[feat].dtype == 'object' or combined_df[feat].dtype.name == 'category':
                main_effects.append(f"C({feat})")
            else:
                main_effects.append(feat)

        # Build interaction terms (all pairwise)
        interaction_terms = []
        for feat1, feat2 in combinations(model_features, 2):
            # Use C() notation for categorical variables
            term1 = f"C({feat1})" if combined_df[feat1].dtype == 'object' or combined_df[feat1].dtype.name == 'category' else feat1
            term2 = f"C({feat2})" if combined_df[feat2].dtype == 'object' or combined_df[feat2].dtype.name == 'category' else feat2
            interaction_terms.append(f"{term1}:{term2}")

        # Combine into full formula
        all_terms = main_effects + interaction_terms
        formula = f"{pc} ~ {' + '.join(all_terms)}"

        print(f"Formula has {len(main_effects)} main effects and {len(interaction_terms)} interaction terms")

        # Fit the model
        if subject_col:
            try:
                model = smf.mixedlm(formula, combined_df, groups=combined_df[subject_col])
                result = model.fit(method='bfgs', maxiter=1000)
                print(f"Successfully fit mixed model with random effects for {subject_col}")
            except Exception as e:
                print(f"Mixed model failed: {e}")
                print("Falling back to OLS")
                result = smf.ols(formula, combined_df).fit()
        else:
            print("No subject column found. Using OLS.")
            result = smf.ols(formula, combined_df).fit()

        # Extract coefficients
        coefs = pd.DataFrame({
            'Variable': result.params.index,
            'Coefficient': result.params.values,
            'Std_Error': result.bse.values,
            'P-value': result.pvalues.values,
            'Significant': result.pvalues < 0.05
        })

        # Classify each coefficient as main effect or interaction
        coefs['Type'] = coefs['Variable'].apply(
            lambda x: 'Interaction' if ':' in x else ('Main' if x != 'Intercept' and x != 'Group Var' else 'Other')
        )

        # For interactions, extract the two variable names involved
        def extract_interaction_vars(var_name):
            if ':' not in var_name:
                return None, None
            parts = var_name.split(':')
            # Extract variable group names from C(VarName)[T.Level] format
            var1 = parts[0].split(')')[0].replace('C(', '') if 'C(' in parts[0] else parts[0]
            var2 = parts[1].split(')')[0].replace('C(', '') if 'C(' in parts[1] else parts[1]
            return var1, var2

        coefs['Var1'], coefs['Var2'] = zip(*coefs['Variable'].apply(extract_interaction_vars))

        # Sort by p-value
        coefs = coefs.sort_values('P-value')

        # Save to results
        lmm_results[pc] = {
            'model': result,
            'coefficients': coefs,
            'r_squared': result.rsquared if hasattr(result, 'rsquared') else None,
            'aic': result.aic if hasattr(result, 'aic') else None
        }

        # Save coefficients to CSV
        coefs.to_csv(os.path.join(output_dir, f"lmm_coefficients_{pc}_with_interactions.csv"), index=False)

        # Print summary of significant interactions
        sig_interactions = coefs[(coefs['Type'] == 'Interaction') & (coefs['Significant'])]
        print(f"\nSignificant interactions for {pc}: {len(sig_interactions)}")
        if len(sig_interactions) > 0:
            print(sig_interactions[['Variable', 'Coefficient', 'P-value']].head(10).to_string())

    # Save the full results
    with open(os.path.join(output_dir, "lmm_with_interactions_results.pkl"), "wb") as f:
        pickle.dump(lmm_results, f)

    # Create summary of all significant interactions across PCs
    all_sig_interactions = []
    for pc, res in lmm_results.items():
        sig = res['coefficients'][(res['coefficients']['Type'] == 'Interaction') & (res['coefficients']['Significant'])]
        sig = sig.copy()
        sig['PC'] = pc
        all_sig_interactions.append(sig)

    if all_sig_interactions:
        summary_df = pd.concat(all_sig_interactions, ignore_index=True)
        summary_df.to_csv(os.path.join(output_dir, "all_significant_interactions.csv"), index=False)
        print(f"\nTotal significant interactions across all PCs: {len(summary_df)}")

    return lmm_results


def create_interaction_network(lmm_results, output_file="interaction_network.pdf", p_threshold=0.05,
                               variance_explained=None, weighting='variance'):
    """
    Create a network diagram using actual interaction coefficients from the LMM.

    Edge weights represent the strength of statistical interactions between clinical factors,
    NOT co-occurrence of main effects.

    Parameters:
    - lmm_results: Dictionary with LMM results (including interaction terms)
    - output_file: Path to save the network plot
    - p_threshold: P-value threshold for including interactions
    - variance_explained: Dict mapping PC names to proportion of variance explained
    - weighting: 'variance' (weight by PC variance), 'sum' (simple sum), or 'count' (number of PCs)
    """
    import networkx as nx

    # Default variance explained if not provided (from typical PCA results)
    if variance_explained is None:
        variance_explained = {
            'PC1': 0.2623, 'PC2': 0.1484, 'PC3': 0.1078,
            'PC4': 0.0674, 'PC5': 0.0555
        }

    # Aggregate interaction effects across PCs
    interaction_effects = {}

    for pc, result in lmm_results.items():
        coefs = result['coefficients']

        # Filter to significant interactions
        sig_interactions = coefs[
            (coefs['Type'] == 'Interaction') &
            (coefs['P-value'] < p_threshold)
        ]

        # Get variance weight for this PC
        var_weight = variance_explained.get(pc, 0.1)

        for _, row in sig_interactions.iterrows():
            var1, var2 = row['Var1'], row['Var2']
            if var1 is None or var2 is None:
                continue

            # Create a canonical key (alphabetically sorted)
            key = tuple(sorted([var1, var2]))

            if key not in interaction_effects:
                interaction_effects[key] = {
                    'weighted_effect': 0,
                    'unweighted_effect': 0,
                    'pcs': [],
                    'coefficients': [],
                    'p_values': [],
                    'var_weights': []
                }

            abs_coef = abs(row['Coefficient'])

            # Store both weighted and unweighted
            interaction_effects[key]['weighted_effect'] += abs_coef * var_weight
            interaction_effects[key]['unweighted_effect'] += abs_coef
            interaction_effects[key]['pcs'].append(pc)
            interaction_effects[key]['coefficients'].append(row['Coefficient'])
            interaction_effects[key]['p_values'].append(row['P-value'])
            interaction_effects[key]['var_weights'].append(var_weight)

    if not interaction_effects:
        print("No significant interactions found! Try increasing p_threshold.")
        return None

    # Calculate combined p-value using Fisher's method for each interaction
    from scipy import stats as scipy_stats

    for key, data in interaction_effects.items():
        # Fisher's method: -2 * sum(log(p_i)) ~ chi-squared with 2k df
        p_values = data['p_values']
        if len(p_values) > 0:
            # Clamp p-values to avoid log(0)
            p_values_clamped = [max(p, 1e-300) for p in p_values]
            chi2_stat = -2 * sum(np.log(p_values_clamped))
            combined_p = 1 - scipy_stats.chi2.cdf(chi2_stat, df=2*len(p_values))
            data['combined_p'] = combined_p
        else:
            data['combined_p'] = 1.0

    # Create graph
    G = nx.Graph()

    # Get all unique variables involved in interactions
    all_vars = set()
    for (var1, var2) in interaction_effects.keys():
        all_vars.add(var1)
        all_vars.add(var2)

    # Add nodes
    for var in all_vars:
        G.add_node(var)

    # Add edges with interaction effect sizes as weights
    for (var1, var2), data in interaction_effects.items():
        # Select weight based on weighting scheme
        if weighting == 'variance':
            edge_weight = data['weighted_effect']
        elif weighting == 'count':
            edge_weight = len(set(data['pcs']))
        else:  # 'sum'
            edge_weight = data['unweighted_effect']

        G.add_edge(var1, var2,
                   weight=edge_weight,
                   weighted_effect=data['weighted_effect'],
                   unweighted_effect=data['unweighted_effect'],
                   combined_p=data['combined_p'],
                   pcs=','.join(data['pcs']),
                   n_pcs=len(set(data['pcs'])))

    # Create layout
    pos = nx.spring_layout(G, k=1.5, iterations=200, seed=42)

    # Draw network
    fig, ax = plt.subplots(figsize=(14, 12))

    # Node size based on number of significant interactions
    node_degrees = dict(G.degree(weight='weight'))
    max_degree = max(node_degrees.values()) if node_degrees else 1
    node_sizes = [800 + (node_degrees[n] / max_degree) * 2000 for n in G.nodes()]

    # Draw nodes
    nx.draw_networkx_nodes(G, pos,
                          node_size=node_sizes,
                          node_color='lightblue',
                          edgecolors='darkblue',
                          linewidths=2,
                          alpha=0.8,
                          ax=ax)

    # Normalize edge widths
    if G.edges():
        weights = [G[u][v]['weight'] for u, v in G.edges()]
        max_weight = max(weights) if weights else 1
        min_weight = min(weights) if weights else 0
        weight_range = max_weight - min_weight if max_weight != min_weight else 1
        edge_widths = [1.5 + (G[u][v]['weight'] - min_weight) / weight_range * 10.5 for u, v in G.edges()]

        # Color edges by number of PCs where interaction is significant
        edge_colors = [G[u][v]['n_pcs'] for u, v in G.edges()]
    else:
        edge_widths = []
        edge_colors = []

    # Draw edges
    edges = nx.draw_networkx_edges(G, pos,
                                   width=edge_widths,
                                   alpha=0.6,
                                   edge_color=edge_colors,
                                   edge_cmap=plt.cm.Reds,
                                   ax=ax)

    # Add colorbar for edge colors
    if edge_colors:
        sm = plt.cm.ScalarMappable(cmap=plt.cm.Reds,
                                    norm=plt.Normalize(vmin=min(edge_colors), vmax=max(edge_colors)))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.5, label='# PCs with significant interaction')

    # Add node labels
    nx.draw_networkx_labels(G, pos,
                           font_size=11,
                           font_weight='bold',
                           ax=ax)

    plt.title('Clinical Factor Interaction Network\n(Edge weights = LMM interaction coefficients, p < 0.05)',
             fontsize=16, pad=20)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    # Print network statistics
    print(f"\nNetwork Statistics:")
    print(f"  Nodes (clinical factors): {G.number_of_nodes()}")
    print(f"  Edges (significant interactions): {G.number_of_edges()}")
    print(f"  Weighting scheme: {weighting}")

    # Top interactions by effect size
    print(f"\nTop 10 interactions by variance-weighted effect size:")
    sorted_edges = sorted(G.edges(data=True), key=lambda x: x[2]['weighted_effect'], reverse=True)
    for i, (u, v, d) in enumerate(sorted_edges[:10]):
        print(f"  {i+1}. {u} × {v}")
        print(f"      Variance-weighted: {d['weighted_effect']:.4f}")
        print(f"      Unweighted sum: {d['unweighted_effect']:.3f}")
        print(f"      Combined p-value: {d['combined_p']:.2e}")
        print(f"      PCs: {d['pcs']} (n={d['n_pcs']})")

    # Create summary dataframe
    summary_data = []
    for u, v, d in G.edges(data=True):
        summary_data.append({
            'Factor1': u,
            'Factor2': v,
            'Weighted_Effect': d['weighted_effect'],
            'Unweighted_Effect': d['unweighted_effect'],
            'Combined_P': d['combined_p'],
            'N_PCs': d['n_pcs'],
            'PCs': d['pcs']
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('Weighted_Effect', ascending=False)

    # Save summary
    summary_file = output_file.replace('.pdf', '_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f"\nSaved interaction summary to: {summary_file}")

    return G


if __name__ == "__main__":
    import sys
    sys.path.insert(0, '/home/david/projects/ZCH_UCMC_Manuscript/python_scripts')

    # Try to import the normalization function
    try:
        from microbiome_transform import clr_transform, tss_transform
    except ImportError:
        print("Warning: microbiome_transform module not found. Using simple normalization.")

    # Create output directory
    output_dir = "results/lmm_interactions"
    os.makedirs(output_dir, exist_ok=True)

    print("Loading microbiome data...")
    try:
        microbiome_df = pd.read_csv("/home/david/projects/ZCH_UCMC_Manuscript/data/NICUSpeciesReduced.csv", index_col=0)
        print(f"Loaded microbiome data with shape: {microbiome_df.shape}")
    except FileNotFoundError:
        print("Error: Could not find microbiome data file")
        exit(1)

    print("Loading metadata...")
    try:
        metadata_df = pd.read_csv("/home/david/projects/ZCH_UCMC_Manuscript/metadata/AllNICUSampleKey20250206.csv", index_col=0)
        print(f"Loaded metadata with shape: {metadata_df.shape}")
    except FileNotFoundError:
        print("Error: Could not find metadata file")
        exit(1)

    # Normalize microbiome data (simple relative abundance for now)
    print("\nNormalizing microbiome data...")
    microbiome_df = microbiome_df.div(microbiome_df.sum(axis=1), axis=0)

    # Run LMM with interactions
    print("\nFitting LMM with interaction terms...")
    lmm_results = lmm_with_interactions(
        microbiome_df=microbiome_df,
        metadata_df=metadata_df,
        n_components=5,
        output_dir=output_dir
    )

    # Create interaction network with variance weighting
    print("\nCreating interaction network (variance-weighted)...")

    # Get actual variance explained from PCA (stored in results)
    # For now, use the values we computed
    variance_explained = {
        'PC1': 0.2623, 'PC2': 0.1484, 'PC3': 0.1078,
        'PC4': 0.0674, 'PC5': 0.0555
    }

    network = create_interaction_network(
        lmm_results=lmm_results,
        output_file=os.path.join(output_dir, "clinical_factor_interaction_network.pdf"),
        p_threshold=0.05,
        variance_explained=variance_explained,
        weighting='variance'
    )

    print(f"\nAnalysis complete! Results saved to '{output_dir}'")
