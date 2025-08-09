#!/usr/bin/env python3
"""
Create comparative visualization of Random Forest SHAP importance and 
Linear Mixed Model coefficients for BSI pathogens
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# Set style for publication-quality figures
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def create_pathogen_comparison_plot(pathogen_name, shap_df, lmm_df, top_n=10):
    """
    Create side-by-side comparison of SHAP importance and LMM coefficients
    for a single pathogen
    
    Parameters:
    -----------
    pathogen_name : str
        Name of the pathogen (e.g., "Klebsiella pneumoniae")
    shap_df : pd.DataFrame
        DataFrame with columns ['Feature', 'SHAP_Importance']
    lmm_df : pd.DataFrame
        DataFrame with columns ['Variable', 'Coefficient', 'P_value']
    top_n : int
        Number of top features to display
    """
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Process SHAP data - get top N features
    shap_top = shap_df.nlargest(top_n, 'SHAP_Importance').sort_values('SHAP_Importance')
    
    # Process LMM data - filter significant and get top by absolute coefficient
    lmm_sig = lmm_df[lmm_df['P_value'] < 0.05].copy()
    lmm_sig['abs_coef'] = lmm_sig['Coefficient'].abs()
    lmm_top = lmm_sig.nlargest(top_n, 'abs_coef').sort_values('Coefficient')
    
    # Clean feature names for display
    def clean_feature_name(name):
        """Clean feature names for better display"""
        replacements = {
            'C(': '',
            ')[T.': ':', 
            ']': '',
            '_': ' ',
            'PostNatalAbxCohort': 'Postnatal Abx',
            'SampleCollectionWeek': 'Week',
            'GestationCohort': 'Gestation',
            'MaternalAntibiotics': 'Maternal Abx',
            'No.Infant.Abx': 'No Infant Abx',
            'Low.Infant.Abx': 'Low Infant Abx',
            'None.Mat.Abx': 'No Maternal Abx',
            'PICC UE': 'PICC Upper Extremity',
            'PICC LE': 'PICC Lower Extremity'
        }
        for old, new in replacements.items():
            name = name.replace(old, new)
        return name.strip()
    
    # Plot 1 (LEFT): Linear Mixed Model Coefficients
    y_pos = np.arange(len(lmm_top))
    
    # Color based on positive/negative coefficients
    colors = ['salmon' if x > 0 else 'steelblue' for x in lmm_top['Coefficient'].values]
    
    # Create horizontal bars
    bars = ax1.barh(y_pos, lmm_top['Coefficient'].values, 
                    color=colors, alpha=0.7, edgecolor='darkgray', linewidth=1)
    
    # Add significance stars
    for i, (coef, pval) in enumerate(zip(lmm_top['Coefficient'].values, lmm_top['P_value'].values)):
        # Add value labels
        offset = 0.1 if coef > 0 else -0.1
        ax1.text(coef + offset, i, f'{coef:.2f}', va='center', ha='left' if coef > 0 else 'right', fontsize=9)
        
        # Add significance stars
        if pval < 0.001:
            sig_marker = '***'
        elif pval < 0.01:
            sig_marker = '**'
        elif pval < 0.05:
            sig_marker = '*'
        else:
            sig_marker = ''
        
        if sig_marker:
            ax1.text(coef + (0.3 if coef > 0 else -0.3), i, sig_marker, 
                    va='center', ha='center', fontsize=10, fontweight='bold')
    
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([clean_feature_name(v) for v in lmm_top['Variable']], fontsize=11)
    ax1.set_xlabel('LMM Effect Size (log odds ratio)', fontsize=12)
    ax1.set_title('Linear Mixed Model', fontsize=14, fontweight='bold')
    ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # Add legend for LMM plot
    pos_patch = mpatches.Patch(color='salmon', alpha=0.7, label='Positive effect')
    neg_patch = mpatches.Patch(color='steelblue', alpha=0.7, label='Negative effect')
    ax1.legend(handles=[pos_patch, neg_patch], loc='best', fontsize=10)
    
    # Plot 2 (RIGHT): Random Forest SHAP Importance
    y_pos = np.arange(len(shap_top))
    colors = ['skyblue'] * len(shap_top)
    
    ax2.barh(y_pos, shap_top['SHAP_Importance'].values, 
             color=colors, alpha=0.7, edgecolor='darkblue', linewidth=1)
    
    # Add value labels
    for i, v in enumerate(shap_top['SHAP_Importance'].values):
        ax2.text(v + 0.002, i, f'{v:.3f}', va='center', fontsize=9)
    
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([clean_feature_name(f) for f in shap_top['Feature']], fontsize=11)
    ax2.set_xlabel('Mean |SHAP Value|', fontsize=12)
    ax2.set_title('Random Forest (SHAP Importance)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_xlim(0, shap_top['SHAP_Importance'].max() * 1.15)
    
    # Overall title
    fig.suptitle(f'Top Predictors for BSI Pathogens — RF vs. LMM\n{pathogen_name}', 
                 fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    return fig

def create_combined_comparison(organisms_data):
    """
    Create a combined figure showing multiple organisms
    
    Parameters:
    -----------
    organisms_data : dict
        Dictionary with organism names as keys and tuple of (shap_df, lmm_df) as values
    """
    n_organisms = len(organisms_data)
    fig, axes = plt.subplots(n_organisms, 2, figsize=(16, 6 * n_organisms))
    
    if n_organisms == 1:
        axes = axes.reshape(1, -1)
    
    for idx, (organism, (shap_df, lmm_df)) in enumerate(organisms_data.items()):
        ax1, ax2 = axes[idx]
        
        # Process and plot LMM data (LEFT)
        lmm_sig = lmm_df[lmm_df['P_value'] < 0.05].copy()
        lmm_sig['abs_coef'] = lmm_sig['Coefficient'].abs()
        lmm_top = lmm_sig.nlargest(10, 'abs_coef').sort_values('Coefficient')
        
        y_pos = np.arange(len(lmm_top))
        colors = ['salmon' if x > 0 else 'steelblue' for x in lmm_top['Coefficient'].values]
        
        ax1.barh(y_pos, lmm_top['Coefficient'].values, 
                color=colors, alpha=0.7, edgecolor='darkgray', linewidth=1)
        
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(lmm_top['Variable'], fontsize=10)
        ax1.set_xlabel('LMM Effect Size (log odds ratio)', fontsize=11)
        ax1.set_title('Linear Mixed Model', fontsize=12)
        ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
        ax1.grid(True, alpha=0.3)
        
        # Add organism label
        ax1.text(0.02, 0.98, organism, transform=ax1.transAxes, 
                fontsize=11, fontweight='bold', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Process and plot SHAP data (RIGHT)
        shap_top = shap_df.nlargest(10, 'SHAP_Importance').sort_values('SHAP_Importance')
        
        y_pos = np.arange(len(shap_top))
        ax2.barh(y_pos, shap_top['SHAP_Importance'].values, 
                 color='skyblue', alpha=0.7, edgecolor='darkblue', linewidth=1)
        
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(shap_top['Feature'], fontsize=10)
        ax2.set_xlabel('Mean |SHAP Value|', fontsize=11)
        ax2.set_title(f'Random Forest (SHAP Importance)', fontsize=12)
        ax2.grid(True, alpha=0.3)
    
    fig.suptitle('Top Predictors for BSI Pathogens — RF vs. LMM', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    return fig

# Example usage with sample data
if __name__ == "__main__":
    # Sample data for Klebsiella pneumoniae (you'll replace with your actual data)
    shap_data_kleb = pd.DataFrame({
        'Feature': ['SampleType:Groin', 'SampleType:Stool', 'SampleCollectionWeek:Week3',
                    'PostNatalAbxCohort:No Infant Abx', 'PostNatalAbxCohort:Low Infant Abx',
                    'PICC_UE', 'UVC', 'Location:Hangzhou', 'AnyMilk:No Milk', 
                    'GestationCohort:28-32', 'MaternalAntibiotics:None', 'BSI_30D'],
        'SHAP_Importance': [0.18, 0.16, 0.15, 0.12, 0.10, 0.09, 0.09, 
                           0.05, 0.04, 0.03, 0.03, 0.02]
    })
    
    lmm_data_kleb = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Groin]', 'C(SampleType)[T.Stool]', 
                     'C(SampleCollectionWeek)[T.Week.3]', 'C(PostNatalAbxCohort)[T.No.Infant.Abx]',
                     'C(PostNatalAbxCohort)[T.Low.Infant.Abx]', 'C(PICC)[T.PICC_UE]',
                     'C(UVC)[T.UVC]', 'C(Location)[T.Hangzhou]', 'C(AnyMilk)[T.No.Milk]',
                     'C(Delivery)[T.Vaginal]', 'C(GestationCohort)[T.28-32]', 
                     'C(MaternalAntibiotics)[T.None.Mat.Abx]'],
        'Coefficient': [2.54, 4.16, 3.14, -1.48, -2.41, 1.48, -1.40, 0.89, -1.82, -0.93, -0.98, -0.62],
        'P_value': [0.00, 0.00, 0.00, 0.02, 0.00, 0.03, 0.05, 0.18, 0.09, 0.16, 0.26, 0.29]
    })
    
    # Sample data for Staphylococcus aureus
    shap_data_staph = pd.DataFrame({
        'Feature': ['SampleType:Groin', 'SampleType:Stool', 'SampleCollectionWeek:Week3',
                    'PostNatalAbxCohort:No Infant Abx', 'PICC_UE', 'Location:Hangzhou',
                    'UVC', 'AnyMilk:No Milk', 'Delivery:Vaginal', 'BSI_30D', 
                    'GestationCohort:28-32', 'NEC_30D'],
        'SHAP_Importance': [0.20, 0.15, 0.13, 0.10, 0.09, 0.07, 0.06, 
                           0.06, 0.05, 0.04, 0.03, 0.02]
    })
    
    lmm_data_staph = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Groin]', 'C(SampleType)[T.Stool]',
                     'C(SampleCollectionWeek)[T.Week.3]', 'C(PostNatalAbxCohort)[T.No.Infant.Abx]',
                     'C(UVC)[T.UVC]', 'C(PICC)[T.PICC_UE]', 'C(Location)[T.Hangzhou]',
                     'C(Delivery)[T.Vaginal]', 'C(GestationCohort)[T.28-32]',
                     'C(AnyMilk)[T.No.Milk]', 'C(MaternalAntibiotics)[T.None.Mat.Abx]',
                     'C(NEC_30D)[T.No.NEC]'],
        'Coefficient': [-1.16, -1.49, -0.94, 0.74, 0.83, 0.70, -0.54, 0.50, 0.61, 0.23, 0.32, 0.43],
        'P_value': [0.00, 0.00, 0.00, 0.02, 0.02, 0.02, 0.08, 0.11, 0.13, 0.65, 0.24, 0.47]
    })
    
    # Create single organism plot
    fig1 = create_pathogen_comparison_plot('Klebsiella pneumoniae', shap_data_kleb, lmm_data_kleb)
    plt.savefig('klebsiella_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    fig2 = create_pathogen_comparison_plot('Staphylococcus aureus', shap_data_staph, lmm_data_staph)
    plt.savefig('staphylococcus_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create combined plot for both organisms
    organisms_data = {
        'K. pneumoniae': (shap_data_kleb, lmm_data_kleb),
        'S. aureus': (shap_data_staph, lmm_data_staph)
    }
    
    fig3 = create_combined_comparison(organisms_data)
    plt.savefig('combined_pathogen_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.show()