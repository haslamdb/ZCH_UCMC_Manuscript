#!/usr/bin/env python3
"""
Create comparative visualization of Random Forest SHAP importance and 
Linear Mixed Model coefficients for BSI pathogens
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# Set font type to TrueType (Type 42) for editable text in PDFs
matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType fonts
matplotlib.rcParams['ps.fonttype'] = 42   # TrueType fonts for PostScript too

# Set font to Arial for all text
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']

# Set style for publication-quality figures
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def create_comparison_label(feature_name):
    """Convert encoded feature names to comparison format"""
    # Handle different feature patterns - for LMM outputs
    if 'C(' in feature_name and '[T.' in feature_name:
        # Handle categorical variable format from statsmodels: C(variable)[T.value]
        import re
        match = re.match(r'C\((\w+)\)\[T\.(.+)\]', feature_name)
        if match:
            category = match.group(1)
            value = match.group(2)
        else:
            return feature_name
    elif ':' in feature_name:
        # Handle format like "SampleType:Stool"
        parts = feature_name.split(':')
        category = parts[0]
        value = parts[1] if len(parts) > 1 else ''
    elif '_' in feature_name:
        # Handle sklearn dummy variable format: variable_value
        parts = feature_name.split('_')
        # Special handling for variables like BSI_30D_1
        if len(parts) == 3 and parts[1] == '30D':
            category = f"{parts[0]}_{parts[1]}"
            value = parts[2]
        else:
            category = parts[0]
            value = '_'.join(parts[1:])
    else:
        # Simple features or intercept
        if feature_name == 'Intercept':
            return 'Intercept (Baseline)'
        return feature_name
    
    # Create simplified labels without "vs" comparisons
    if category == 'SampleType':
        if value == 'Groin':
            return 'Groin Sample'
        elif value == 'Stool':
            return 'Stool Sample'
        elif value in ['Nonese', 'Nose']:
            return 'Nose Sample'
        elif value == 'Axilla':
            return 'Axilla Sample'
    elif category == 'Location':
        if value in ['ZCH', 'Hangzhou']:
            return 'ZCH Location'
        elif value in ['Cincinnati', 'UCMC']:
            return 'UCMC Location'
    elif category == 'GestationCohort':
        if value == '28-32 Weeks' or value == '28-32':
            return '28-32 Weeks Gestation'
        elif value == '33-36 Weeks' or value == '33-36':
            return '33-36 Weeks Gestation'
        elif value == '25-28' or value == '23-27 Weeks':
            return '23-27 Weeks Gestation'
    elif category == 'SampleCollectionWeek':
        if value in ['Week.3', 'Week 3', 'Week3']:
            return 'Week 3 Sample'
        elif value in ['Week.1', 'Week 1', 'Week1']:
            return 'Week 1 Sample'
    elif category == 'MaternalAntibiotics':
        if value in ['None.Mat.Abx', 'No Mat Abx', 'None']:
            return 'No Maternal Antibiotics'
        elif value in ['Mat.Abx', 'Mat Abx']:
            return 'Maternal Antibiotics'
    elif category == 'PostNatalAbxCohort':
        if value in ['No.Infant.Abx', 'No Infant Abx']:
            return 'No Infant Antibiotics'
        elif value in ['Low.Infant.Abx', 'Low Infant Abx']:
            return 'Low Infant Antibiotics'
        elif value in ['High.Infant.Abx', 'High Infant Abx']:
            return 'High Infant Antibiotics'
    elif category in ['BSI_30D', 'BSI30D']:
        if value == '1':
            return 'BSI Positive'
        elif value == '0':
            return 'BSI Negative'
    elif category in ['NEC_30D', 'NEC30D']:
        if value == '1':
            return 'NEC Positive'
        elif value == '0':
            return 'NEC Negative'
    elif category == 'AnyMilk':
        if value == '1' or value == 'Mother':
            return 'Mother\'s Milk'
        elif value == 'No Milk' or value == 'No.Milk':
            return 'No Milk'
        elif value == 'Formula':
            return 'Formula'
    elif category == 'PICC':
        if value == '1':
            return 'PICC Present'
        elif value == 'PICC_UE' or value == 'UE':
            return 'PICC Upper Extremity'
        elif value == 'PICC_LE' or value == 'LE':
            return 'PICC Lower Extremity'
        elif value == 'PICC_Neck' or value == 'Neck':
            return 'PICC Neck'
        elif value == 'axillary':
            return 'PICC Axillary'
        elif value == 'peripheral_UE':
            return 'Peripheral Upper Extremity'
    elif category == 'UVC':
        if value == '1' or value == 'UVC':
            return 'UVC Present'
        elif value == '0':
            return 'No UVC'
    elif category == 'Delivery':
        if value == 'Vaginal':
            return 'Vaginal Delivery'
        elif value == 'Cesarean':
            return 'Cesarean Delivery'
    
    # Default: return original if no mapping found
    return feature_name

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
    
    # Use create_comparison_label for better display names
    
    # Plot 1 (LEFT): Linear Mixed Model Coefficients
    y_pos = np.arange(len(lmm_top))
    
    # Color based on positive/negative coefficients
    colors = ['salmon' if x > 0 else 'steelblue' for x in lmm_top['Coefficient'].values]
    
    # Create horizontal bars
    bars = ax1.barh(y_pos, lmm_top['Coefficient'].values, 
                    color=colors, alpha=0.7, edgecolor='darkgray', linewidth=1)
    
    # Add significance stars combined with value labels
    for i, (coef, pval) in enumerate(zip(lmm_top['Coefficient'].values, lmm_top['P_value'].values)):
        # Determine significance markers
        if pval < 0.001:
            sig_marker = '***'
        elif pval < 0.01:
            sig_marker = '**'
        elif pval < 0.05:
            sig_marker = '*'
        else:
            sig_marker = ''
        
        # Combine coefficient value with significance markers
        label_text = f'{coef:.2f}{sig_marker}'
        
        # Position label with consistent offset
        offset = 0.25 if coef > 0 else -0.25
        ax1.text(coef + offset, i, label_text, va='center', 
                ha='left' if coef > 0 else 'right', fontsize=13, fontweight='bold')
    
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([create_comparison_label(v) for v in lmm_top['Variable']], fontsize=13)
    ax1.set_xlabel('LMM Effect Size (log odds ratio)', fontsize=16, fontweight='bold')
    ax1.set_title('Linear Mixed Model', fontsize=16, fontweight='bold')
    ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # Increase x-axis tick label size
    ax1.tick_params(axis='x', labelsize=14)
    
    # Add legend for LMM plot
    pos_patch = mpatches.Patch(color='salmon', alpha=0.7, label='Positive effect')
    neg_patch = mpatches.Patch(color='steelblue', alpha=0.7, label='Negative effect')
    ax1.legend(handles=[pos_patch, neg_patch], loc='best', fontsize=12)
    
    # Plot 2 (RIGHT): Random Forest SHAP Importance
    y_pos = np.arange(len(shap_top))
    
    ax2.barh(y_pos, shap_top['SHAP_Importance'].values, 
             color='skyblue', alpha=0.7, edgecolor='darkblue', linewidth=1)
    
    # Add value labels
    for i, v in enumerate(shap_top['SHAP_Importance'].values):
        ax2.text(v + 0.002, i, f'{v:.3f}', va='center', fontsize=13, fontweight='bold')
    
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([create_comparison_label(f) for f in shap_top['Feature']], fontsize=13)
    ax2.set_xlabel('Mean |SHAP Value|', fontsize=16, fontweight='bold')
    ax2.set_title('Random Forest (SHAP Importance)', fontsize=16, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_xlim(0, shap_top['SHAP_Importance'].max() * 1.15)
    
    # Increase x-axis tick label size
    ax2.tick_params(axis='x', labelsize=14)
    
    # Overall title
    fig.suptitle(f'Top Predictors for BSI Pathogens — LMM and Random Forest\n{pathogen_name}', 
                 fontsize=18, fontweight='bold', y=1.02)
    
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
        
        # Add significance stars combined with value labels for combined plots
        for i, (coef, pval) in enumerate(zip(lmm_top['Coefficient'].values, lmm_top['P_value'].values)):
            # Determine significance markers
            if pval < 0.001:
                sig_marker = '***'
            elif pval < 0.01:
                sig_marker = '**'
            elif pval < 0.05:
                sig_marker = '*'
            else:
                sig_marker = ''
            
            # Combine coefficient value with significance markers
            label_text = f'{coef:.2f}{sig_marker}'
            
            # Position label with consistent offset
            offset = 0.2 if coef > 0 else -0.2
            ax1.text(coef + offset, i, label_text, va='center', 
                    ha='left' if coef > 0 else 'right', fontsize=12, fontweight='bold')
        
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels([create_comparison_label(v) for v in lmm_top['Variable']], fontsize=12)
        ax1.set_xlabel('LMM Effect Size (log odds ratio)', fontsize=16, fontweight='bold')
        ax1.set_title('Linear Mixed Model', fontsize=14)
        ax1.tick_params(axis='x', labelsize=14)  # Increase X-axis tick label size
        ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
        ax1.grid(True, alpha=0.3)
        
        # Add organism label
        ax1.text(0.02, 0.98, organism, transform=ax1.transAxes, 
                fontsize=13, fontweight='bold', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Process and plot SHAP data (RIGHT)
        shap_top = shap_df.nlargest(10, 'SHAP_Importance').sort_values('SHAP_Importance')
        
        y_pos = np.arange(len(shap_top))
        
        ax2.barh(y_pos, shap_top['SHAP_Importance'].values, 
                 color='skyblue', alpha=0.7, edgecolor='darkblue', linewidth=1)
        
        # Add value labels for SHAP importance
        for i, v in enumerate(shap_top['SHAP_Importance'].values):
            ax2.text(v + 0.002, i, f'{v:.3f}', va='center', fontsize=12, fontweight='bold')
        
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels([create_comparison_label(f) for f in shap_top['Feature']], fontsize=12)
        ax2.set_xlabel('Mean |SHAP Value|', fontsize=16, fontweight='bold')
        ax2.set_title(f'Random Forest (SHAP Importance)', fontsize=14)
        ax2.tick_params(axis='x', labelsize=14)  # Increase X-axis tick label size
        ax2.grid(True, alpha=0.3)
    
    fig.suptitle('Top Predictors for BSI Pathogens — LMM and Random Forest', 
                 fontsize=18, fontweight='bold')
    plt.tight_layout()
    return fig

# Example usage with sample data
if __name__ == "__main__":
    # Dictionary to store all organism data
    all_organisms_data = {}
    
    # 1. Staphylococcus aureus (Higher at UCMC: 41.1% vs 10.1%)
    shap_data_staph = pd.DataFrame({
        'Feature': ['SampleType:Groin', 'SampleType:Stool', 'SampleCollectionWeek:Week3',
                    'PostNatalAbxCohort:No Infant Abx', 'PICC_UE', 'Location:ZCH',
                    'UVC', 'AnyMilk:No Milk', 'Delivery:Vaginal', 'BSI_30D', 
                    'GestationCohort:28-32', 'NEC_30D'],
        'SHAP_Importance': [0.20, 0.15, 0.13, 0.10, 0.09, 0.07, 0.06, 
                           0.06, 0.05, 0.04, 0.03, 0.02]
    })
    
    lmm_data_staph = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Groin]', 'C(SampleType)[T.Stool]',
                     'C(SampleCollectionWeek)[T.Week.3]', 'C(PostNatalAbxCohort)[T.No.Infant.Abx]',
                     'C(UVC)[T.UVC]', 'C(PICC)[T.PICC_UE]', 'C(Location)[T.ZCH]',
                     'C(Delivery)[T.Vaginal]', 'C(GestationCohort)[T.28-32]',
                     'C(AnyMilk)[T.No.Milk]', 'C(MaternalAntibiotics)[T.None.Mat.Abx]',
                     'C(NEC_30D)[T.No.NEC]'],
        'Coefficient': [-1.16, -1.49, -0.94, 0.74, 0.83, 0.70, -0.54, 0.50, 0.61, 0.23, 0.32, 0.43],
        'P_value': [0.00, 0.00, 0.00, 0.02, 0.02, 0.02, 0.08, 0.11, 0.13, 0.65, 0.24, 0.47]
    })
    all_organisms_data['Staphylococcus aureus'] = (shap_data_staph, lmm_data_staph)
    
    # 2. Klebsiella pneumoniae (Higher at ZCH: 8.4% vs 25.6%)
    shap_data_kleb = pd.DataFrame({
        'Feature': ['SampleType:Groin', 'SampleType:Stool', 'SampleCollectionWeek:Week3',
                    'PostNatalAbxCohort:No Infant Abx', 'PostNatalAbxCohort:Low Infant Abx',
                    'PICC_UE', 'UVC', 'Location:ZCH', 'AnyMilk:No Milk', 
                    'GestationCohort:28-32', 'MaternalAntibiotics:None', 'BSI_30D'],
        'SHAP_Importance': [0.18, 0.16, 0.15, 0.12, 0.10, 0.09, 0.09, 
                           0.05, 0.04, 0.03, 0.03, 0.02]
    })
    
    lmm_data_kleb = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Groin]', 'C(SampleType)[T.Stool]', 
                     'C(SampleCollectionWeek)[T.Week.3]', 'C(PostNatalAbxCohort)[T.No.Infant.Abx]',
                     'C(PostNatalAbxCohort)[T.Low.Infant.Abx]', 'C(PICC)[T.PICC_UE]',
                     'C(UVC)[T.UVC]', 'C(Location)[T.ZCH]', 'C(AnyMilk)[T.No.Milk]',
                     'C(Delivery)[T.Vaginal]', 'C(GestationCohort)[T.28-32]', 
                     'C(MaternalAntibiotics)[T.None.Mat.Abx]'],
        'Coefficient': [2.54, 4.16, 3.14, -1.48, -2.41, 1.48, -1.40, 0.89, -1.82, -0.93, -0.98, -0.62],
        'P_value': [0.00, 0.00, 0.00, 0.02, 0.00, 0.03, 0.05, 0.18, 0.09, 0.16, 0.26, 0.29]
    })
    all_organisms_data['Klebsiella pneumoniae'] = (shap_data_kleb, lmm_data_kleb)
    
    # 3. Klebsiella oxytoca (Higher at UCMC: 6.5% vs 0.0%)
    # You'll need to replace with actual data from your analysis
    shap_data_k_oxy = pd.DataFrame({
        'Feature': ['SampleType:Stool', 'SampleType:Groin', 'SampleCollectionWeek:Week3',
                    'Location:ZCH', 'PostNatalAbxCohort:No Infant Abx', 'PICC_UE',
                    'UVC', 'AnyMilk:No Milk', 'Delivery:Vaginal', 'BSI_30D'],
        'SHAP_Importance': [0.14, 0.12, 0.11, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03]
    })
    
    lmm_data_k_oxy = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Stool]', 'C(SampleType)[T.Groin]',
                     'C(SampleCollectionWeek)[T.Week.3]', 'C(Location)[T.ZCH]',
                     'C(PostNatalAbxCohort)[T.No.Infant.Abx]', 'C(PICC)[T.PICC_UE]',
                     'C(UVC)[T.UVC]', 'C(Delivery)[T.Vaginal]'],
        'Coefficient': [-1.8, -1.2, -0.9, -2.1, 0.6, 0.5, 0.7, 0.3],
        'P_value': [0.01, 0.02, 0.03, 0.001, 0.04, 0.05, 0.02, 0.10]
    })
    all_organisms_data['Klebsiella oxytoca'] = (shap_data_k_oxy, lmm_data_k_oxy)
    
    # 4. Escherichia coli 
    shap_data_e_coli = pd.DataFrame({
        'Feature': ['SampleType:Stool', 'SampleType:Groin', 'SampleCollectionWeek:Week3',
                    'PostNatalAbxCohort:No Infant Abx', 'Location:ZCH', 'UVC',
                    'GestationCohort:28-32', 'AnyMilk:No Milk', 'PICC_LE', 'Delivery:Vaginal'],
        'SHAP_Importance': [0.17, 0.15, 0.14, 0.11, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04]
    })
    
    lmm_data_e_coli = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Stool]', 'C(SampleType)[T.Groin]',
                     'C(SampleCollectionWeek)[T.Week.3]', 'C(PostNatalAbxCohort)[T.No.Infant.Abx]',
                     'C(UVC)[T.UVC]', 'C(Location)[T.ZCH]', 'C(GestationCohort)[T.28-32]',
                     'C(PICC)[T.PICC_LE]', 'C(Delivery)[T.Vaginal]', 'C(AnyMilk)[T.No.Milk]'],
        'Coefficient': [3.75, 2.07, 2.35, -0.71, -2.13, -0.95, -1.68, -1.32, -0.86, -0.38],
        'P_value': [0.001, 0.001, 0.001, 0.22, 0.001, 0.11, 0.03, 0.05, 0.15, 0.70]
    })
    all_organisms_data['Escherichia coli'] = (shap_data_e_coli, lmm_data_e_coli)
    
    # 5. Serratia marcescens (Higher at UCMC: 6.2% vs 0.8%)
    shap_data_serratia = pd.DataFrame({
        'Feature': ['SampleType:Stool', 'SampleType:Groin', 'SampleCollectionWeek:Week3',
                    'Location:ZCH', 'PostNatalAbxCohort:No Infant Abx', 'PICC_UE',
                    'UVC', 'AnyMilk:No Milk', 'MaternalAntibiotics:None', 'Delivery:Vaginal'],
        'SHAP_Importance': [0.15, 0.13, 0.12, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04]
    })
    
    lmm_data_serratia = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Stool]', 'C(SampleType)[T.Groin]',
                     'C(SampleCollectionWeek)[T.Week.3]', 'C(Location)[T.ZCH]',
                     'C(PostNatalAbxCohort)[T.No.Infant.Abx]', 'C(PICC)[T.PICC_UE]',
                     'C(UVC)[T.UVC]', 'C(AnyMilk)[T.No.Milk]'],
        'Coefficient': [-1.5, -1.1, -0.8, -1.9, 0.5, 0.6, 0.8, -0.3],
        'P_value': [0.01, 0.02, 0.04, 0.002, 0.05, 0.03, 0.01, 0.15]
    })
    all_organisms_data['Serratia marcescens'] = (shap_data_serratia, lmm_data_serratia)
    
    # 6. Enterococcus faecalis (Higher at ZCH: 3.6% vs 13.9%)
    shap_data_e_faecalis = pd.DataFrame({
        'Feature': ['SampleType:Stool', 'SampleType:Groin', 'SampleCollectionWeek:Week3',
                    'Location:ZCH', 'PostNatalAbxCohort:Low Infant Abx', 'PICC_LE',
                    'UVC', 'AnyMilk:No Milk', 'Delivery:Vaginal', 'GestationCohort:28-32'],
        'SHAP_Importance': [0.17, 0.14, 0.13, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05]
    })
    
    lmm_data_e_faecalis = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Stool]', 'C(SampleType)[T.Groin]',
                     'C(SampleCollectionWeek)[T.Week.3]', 'C(Location)[T.ZCH]',
                     'C(PostNatalAbxCohort)[T.Low.Infant.Abx]', 'C(PICC)[T.PICC_LE]',
                     'C(UVC)[T.UVC]', 'C(AnyMilk)[T.No.Milk]'],
        'Coefficient': [2.8, 2.2, 1.9, 1.5, -1.8, 1.2, -0.9, -1.1],
        'P_value': [0.001, 0.003, 0.01, 0.02, 0.01, 0.03, 0.04, 0.05]
    })
    all_organisms_data['Enterococcus faecalis'] = (shap_data_e_faecalis, lmm_data_e_faecalis)
    
    # 7. Enterococcus faecium (Higher at ZCH: 0.0% vs 6.7%)
    shap_data_e_faecium = pd.DataFrame({
        'Feature': ['SampleType:Stool', 'SampleType:Groin', 'SampleCollectionWeek:Week3',
                    'Location:ZCH', 'PostNatalAbxCohort:Low Infant Abx', 'PICC_LE',
                    'UVC', 'AnyMilk:No Milk', 'MaternalAntibiotics:None', 'Delivery:Vaginal'],
        'SHAP_Importance': [0.16, 0.15, 0.14, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06]
    })
    
    lmm_data_e_faecium = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Stool]', 'C(SampleType)[T.Groin]',
                     'C(SampleCollectionWeek)[T.Week.3]', 'C(Location)[T.ZCH]',
                     'C(PostNatalAbxCohort)[T.Low.Infant.Abx]', 'C(PICC)[T.PICC_LE]',
                     'C(UVC)[T.UVC]', 'C(Delivery)[T.Vaginal]'],
        'Coefficient': [3.1, 2.5, 2.0, 1.8, -2.0, 1.4, -1.0, -0.8],
        'P_value': [0.001, 0.002, 0.01, 0.01, 0.008, 0.02, 0.03, 0.06]
    })
    all_organisms_data['Enterococcus faecium'] = (shap_data_e_faecium, lmm_data_e_faecium)
    
    # 8. Streptococcus pyogenes (Higher at ZCH: 0.0% vs 4.2%)
    shap_data_s_pyogenes = pd.DataFrame({
        'Feature': ['SampleType:Axilla', 'SampleType:Groin', 'SampleCollectionWeek:Week1',
                    'Location:ZCH', 'PostNatalAbxCohort:Low Infant Abx', 'PICC_UE',
                    'UVC', 'AnyMilk:No Milk', 'Delivery:Cesarean', 'GestationCohort:33-36'],
        'SHAP_Importance': [0.14, 0.13, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05]
    })
    
    lmm_data_s_pyogenes = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Axilla]', 'C(SampleType)[T.Groin]',
                     'C(SampleCollectionWeek)[T.Week.1]', 'C(Location)[T.ZCH]',
                     'C(PostNatalAbxCohort)[T.Low.Infant.Abx]', 'C(PICC)[T.PICC_UE]',
                     'C(UVC)[T.UVC]', 'C(Delivery)[T.Cesarean]'],
        'Coefficient': [2.3, 1.9, 1.6, 2.1, -1.5, -0.7, -0.9, 0.8],
        'P_value': [0.002, 0.005, 0.01, 0.001, 0.02, 0.05, 0.03, 0.08]
    })
    all_organisms_data['Streptococcus pyogenes'] = (shap_data_s_pyogenes, lmm_data_s_pyogenes)
    
    # 9. Candida parapsilosis (Higher at ZCH: 0.0% vs 4.6%)
    # Note: You mentioned in the manuscript that Candida wasn't the focus, but including for completeness
    shap_data_c_para = pd.DataFrame({
        'Feature': ['SampleType:Stool', 'SampleType:Groin', 'SampleCollectionWeek:Week3',
                    'Location:ZCH', 'PostNatalAbxCohort:High Infant Abx', 'PICC_LE',
                    'UVC', 'AnyMilk:Formula', 'Delivery:Vaginal', 'GestationCohort:25-28'],
        'SHAP_Importance': [0.18, 0.16, 0.14, 0.13, 0.12, 0.10, 0.08, 0.06, 0.04, 0.03]
    })
    
    lmm_data_c_para = pd.DataFrame({
        'Variable': ['C(SampleType)[T.Stool]', 'C(SampleType)[T.Groin]',
                     'C(SampleCollectionWeek)[T.Week.3]', 'C(Location)[T.ZCH]',
                     'C(PostNatalAbxCohort)[T.High.Infant.Abx]', 'C(PICC)[T.PICC_LE]',
                     'C(UVC)[T.UVC]', 'C(AnyMilk)[T.Formula]'],
        'Coefficient': [2.7, 2.1, 1.7, 2.5, 1.9, 1.3, -0.8, 1.0],
        'P_value': [0.001, 0.003, 0.01, 0.001, 0.005, 0.02, 0.05, 0.04]
    })
    all_organisms_data['Candida parapsilosis'] = (shap_data_c_para, lmm_data_c_para)
    
    # Generate individual plots for each organism
    for organism_name, (shap_df, lmm_df) in all_organisms_data.items():
        fig = create_pathogen_comparison_plot(organism_name, shap_df, lmm_df)
        # Clean filename - remove spaces and special characters
        clean_name = organism_name.replace(' ', '_').replace('.', '')
        plt.savefig(f'{clean_name}_comparison.pdf', dpi=300, bbox_inches='tight')
        plt.close()  # Close to save memory
        print(f"Saved plot for {organism_name}")
    
    # Create combined plots for selected organisms (e.g., top gram-positive and gram-negative)
    # Gram-positive selection
    gram_positive_organisms = {
        'S. aureus': all_organisms_data['Staphylococcus aureus'],
        'E. faecalis': all_organisms_data['Enterococcus faecalis'],
        'E. faecium': all_organisms_data['Enterococcus faecium'],
        'S. pyogenes': all_organisms_data['Streptococcus pyogenes']
    }
    
    fig_gp = create_combined_comparison(gram_positive_organisms)
    plt.savefig('gram_positive_pathogens_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Gram-negative selection
    gram_negative_organisms = {
        'K. pneumoniae': all_organisms_data['Klebsiella pneumoniae'],
        'K. oxytoca': all_organisms_data['Klebsiella oxytoca'],
        'E. coli': all_organisms_data['Escherichia coli'],
        'S. marcescens': all_organisms_data['Serratia marcescens']
    }
    
    fig_gn = create_combined_comparison(gram_negative_organisms)
    plt.savefig('gram_negative_pathogens_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Show one example plot
    fig_example = create_pathogen_comparison_plot('Klebsiella pneumoniae', 
                                                  all_organisms_data['Klebsiella pneumoniae'][0],
                                                  all_organisms_data['Klebsiella pneumoniae'][1])
    plt.show()
    
    print("\nAll plots have been generated and saved as PDFs.")