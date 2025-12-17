#!/usr/bin/env python3
"""
Diverging bar chart showing organisms with significantly different 
prevalence between UCMC and ZCH hospitals.

Author: Generated with Claude
Date: 2025-01-17
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# =============================================================================
# CONFIGURATION - Edit these values as needed
# =============================================================================
INPUT_FILE = 'Copy_of_2025_4.csv'  # Path to your CSV file
OUTPUT_DIR = '.'                    # Output directory for figures
OUTPUT_FORMAT = 'pdf'               # 'pdf' or 'png'
DPI = 150                           # DPI for PNG output

# Colors (Tableau palette)
COLOR_ZCH = '#4e79a7'  # Blue for ZCH-elevated organisms
COLOR_UC = '#f28e2b'   # Orange for UC-elevated organisms

# Set PDF/PS font type for Illustrator compatibility (TrueType fonts)
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# =============================================================================
# LOAD AND PROCESS DATA
# =============================================================================
df = pd.read_csv(INPUT_FILE)
df.columns = ['Organism', 'UC', 'ZCH']

# Calculate totals
total_uc = df['UC'].sum()
total_zch = df['ZCH'].sum()

print(f"Total isolates - UC: {total_uc}, ZCH: {total_zch}")

# Filter to significant organisms (marked with *) plus E. coli
df_sig = df[
    df['Organism'].str.contains(r'\*', regex=True) | 
    (df['Organism'] == 'Escherichia coli')
].copy()
df_sig['Organism'] = df_sig['Organism'].str.replace('*', '', regex=False)

print(f"\nSelected organisms: {len(df_sig)}")

# Calculate proportions (% of total isolates at each hospital)
df_sig['UC_pct'] = (df_sig['UC'] / total_uc) * 100
df_sig['ZCH_pct'] = (df_sig['ZCH'] / total_zch) * 100

# Log2 fold change of proportions (with pseudocount for zeros)
pseudo = 0.5
df_sig['log2FC'] = (
    np.log2((df_sig['UC'] + pseudo) / total_uc) - 
    np.log2((df_sig['ZCH'] + pseudo) / total_zch)
)

# Simple difference in percentages
df_sig['pct_diff'] = df_sig['UC_pct'] - df_sig['ZCH_pct']

print("\nData with metrics:")
print(df_sig[['Organism', 'UC', 'ZCH', 'UC_pct', 'ZCH_pct', 'log2FC', 'pct_diff']].to_string())

# =============================================================================
# FIGURE 1: Side-by-side comparison (Log2FC and Percentage Difference)
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Sort ascending so UC-elevated organisms appear at TOP of horizontal bar chart
df_plot_log2 = df_sig.sort_values('log2FC', ascending=True)
df_plot_pct = df_sig.sort_values('pct_diff', ascending=True)

# Plot 1: Log2 Fold Change
ax1 = axes[0]
colors1 = [COLOR_ZCH if x < 0 else COLOR_UC for x in df_plot_log2['log2FC']]
ax1.barh(df_plot_log2['Organism'], df_plot_log2['log2FC'], color=colors1, edgecolor='white', linewidth=0.5)
ax1.axvline(x=0, color='black', linewidth=0.8)
ax1.set_xlabel('Log₂ Fold Change (UC vs ZCH)', fontsize=11)
ax1.set_title('Log₂ Fold Change of Proportions', fontsize=12, fontweight='bold')
ax1.text(ax1.get_xlim()[0] + 0.3, len(df_plot_log2) - 0.3, '← Elevated in ZCH', fontsize=9, color=COLOR_ZCH, fontweight='bold')
ax1.text(ax1.get_xlim()[1] - 0.3, len(df_plot_log2) - 0.3, 'Elevated in UC →', fontsize=9, color=COLOR_UC, fontweight='bold', ha='right')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Plot 2: Percentage Difference
ax2 = axes[1]
colors2 = [COLOR_ZCH if x < 0 else COLOR_UC for x in df_plot_pct['pct_diff']]
ax2.barh(df_plot_pct['Organism'], df_plot_pct['pct_diff'], color=colors2, edgecolor='white', linewidth=0.5)
ax2.axvline(x=0, color='black', linewidth=0.8)
ax2.set_xlabel('Difference in Percentage (UC % − ZCH %)', fontsize=11)
ax2.set_title('Difference in Percentage of Total Isolates', fontsize=12, fontweight='bold')
ax2.text(ax2.get_xlim()[0] + 0.5, len(df_plot_pct) - 0.3, '← Elevated in ZCH', fontsize=9, color=COLOR_ZCH, fontweight='bold')
ax2.text(ax2.get_xlim()[1] - 0.5, len(df_plot_pct) - 0.3, 'Elevated in UC →', fontsize=9, color=COLOR_UC, fontweight='bold', ha='right')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.suptitle('Organisms with Significantly Different Prevalence\nBetween UC and ZCH Hospitals (Infant BSI)', 
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/bsi_comparison_both_metrics.{OUTPUT_FORMAT}', format=OUTPUT_FORMAT,
            dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

# =============================================================================
# FIGURE 2: Single panel - Log2 Fold Change with counts
# =============================================================================
fig, ax = plt.subplots(figsize=(10, 6))

# Sort ascending so UC-elevated organisms appear at TOP
df_plot = df_sig.sort_values('log2FC', ascending=True)

colors = [COLOR_ZCH if x < 0 else COLOR_UC for x in df_plot['log2FC']]
bars = ax.barh(df_plot['Organism'], df_plot['log2FC'], color=colors, edgecolor='white', linewidth=0.5, height=0.7)

ax.axvline(x=0, color='black', linewidth=1)
ax.set_xlabel('Log₂ Fold Change (UC vs ZCH)', fontsize=11)
ax.set_title('Organisms with Significantly Different Prevalence\nBetween UC and ZCH Hospitals (Infant BSI)', 
             fontsize=12, fontweight='bold')

# Add count labels
for bar, (_, row) in zip(bars, df_plot.iterrows()):
    width = bar.get_width()
    label = f"UC:{int(row['UC'])} / ZCH:{int(row['ZCH'])}"
    if width >= 0:
        ax.text(width + 0.1, bar.get_y() + bar.get_height()/2, label, 
                va='center', ha='left', fontsize=8, color='gray')
    else:
        ax.text(width - 0.1, bar.get_y() + bar.get_height()/2, label, 
                va='center', ha='right', fontsize=8, color='gray')

xlim = ax.get_xlim()
ax.text(xlim[0] + 0.3, len(df_plot) - 0.3, '← Elevated in ZCH', fontsize=10, color=COLOR_ZCH, fontweight='bold')
ax.text(xlim[1] - 0.3, len(df_plot) - 0.3, 'Elevated in UC →', fontsize=10, color=COLOR_UC, fontweight='bold', ha='right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/bsi_comparison_log2fc.{OUTPUT_FORMAT}', format=OUTPUT_FORMAT,
            dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

# =============================================================================
# FIGURE 3: Single panel - Percentage Difference with counts
# =============================================================================
fig, ax = plt.subplots(figsize=(10, 6))

# Sort ascending so UC-elevated organisms appear at TOP
df_plot = df_sig.sort_values('pct_diff', ascending=True)

colors = [COLOR_ZCH if x < 0 else COLOR_UC for x in df_plot['pct_diff']]
bars = ax.barh(df_plot['Organism'], df_plot['pct_diff'], color=colors, edgecolor='white', linewidth=0.5, height=0.7)

ax.axvline(x=0, color='black', linewidth=1)
ax.set_xlabel('Difference in Percentage of Total Isolates (UC − ZCH)', fontsize=11)
ax.set_title('Organisms with Significantly Different Prevalence\nBetween UC and ZCH Hospitals (Infant BSI)', 
             fontsize=12, fontweight='bold')

# Add count labels
for bar, (_, row) in zip(bars, df_plot.iterrows()):
    width = bar.get_width()
    label = f"UC:{int(row['UC'])} / ZCH:{int(row['ZCH'])}"
    if width >= 0:
        ax.text(width + 0.3, bar.get_y() + bar.get_height()/2, label, 
                va='center', ha='left', fontsize=8, color='gray')
    else:
        ax.text(width - 0.3, bar.get_y() + bar.get_height()/2, label, 
                va='center', ha='right', fontsize=8, color='gray')

xlim = ax.get_xlim()
ax.text(xlim[0] + 1, len(df_plot) - 0.3, '← Elevated in ZCH', fontsize=10, color=COLOR_ZCH, fontweight='bold')
ax.text(xlim[1] - 1, len(df_plot) - 0.3, 'Elevated in UC →', fontsize=10, color=COLOR_UC, fontweight='bold', ha='right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/bsi_comparison_pct_diff.{OUTPUT_FORMAT}', format=OUTPUT_FORMAT,
            dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

print("\n✓ Figures saved!")
print(f"  - {OUTPUT_DIR}/bsi_comparison_both_metrics.{OUTPUT_FORMAT} (side-by-side comparison)")
print(f"  - {OUTPUT_DIR}/bsi_comparison_log2fc.{OUTPUT_FORMAT} (log2 fold change)")
print(f"  - {OUTPUT_DIR}/bsi_comparison_pct_diff.{OUTPUT_FORMAT} (percentage difference)")
