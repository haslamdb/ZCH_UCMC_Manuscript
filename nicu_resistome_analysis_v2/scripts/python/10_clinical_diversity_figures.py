#!/usr/bin/env python3
"""
Clinical figures: AMR diversity/richness/abundance + volcano panels.

Outputs are designed to read cleanly side-by-side with the clinical-LMM
result tables (`results/differential_*/...`). Specifically:

1. Diversity panel: richness, total RPM, Shannon by Body Site stratified
   by Week and by PostNatalAbxCohort.
2. Volcano panel: gene-family-level LMM coefficient vs -log10(FDR), one
   panel per effect of clinical interest (Stool vs Axilla, Groin vs
   Axilla, Week 3 vs Week 1, High abx vs None, Low abx vs None).

Pulls diversity metrics from `results/exploratory/diversity_metrics.tsv`
(produced by 03) and the gene-family LMM coefficients from
`results/differential_gene_family/gene_family_clinical_effects_long.tsv`
(produced by 04d).
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings('ignore')

PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
V2 = PROJECT_ROOT / "nicu_resistome_analysis_v2"
QC_DIR = V2 / "results" / "qc"
EXPLORATORY_DIR = V2 / "results" / "exploratory"
GF_DIR = V2 / "results" / "differential_gene_family"
FIGURES_DIR = V2 / "figures" / "clinical"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.1)

SITE_ORDER = ['Axilla', 'Groin', 'Stool']
WEEK_ORDER = ['Week.1', 'Week.3']
ABX_ORDER  = ['No.Infant.Abx', 'Low.Infant.Abx', 'High.Infant.Abx']
ABX_LABELS = {'No.Infant.Abx': 'No abx', 'Low.Infant.Abx': 'Low abx', 'High.Infant.Abx': 'High abx'}

SITE_PALETTE = {'Axilla': '#1f77b4', 'Groin': '#ff7f0e', 'Stool': '#2ca02c'}
WEEK_PALETTE = {'Week.1': '#9ecae1', 'Week.3': '#3182bd'}
ABX_PALETTE  = {'No.Infant.Abx': '#bdbdbd', 'Low.Infant.Abx': '#fdae6b', 'High.Infant.Abx': '#e6550d'}


# ---------------------------------------------------------------------------
# Diversity panel
# ---------------------------------------------------------------------------

def load_diversity():
    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")
    div = pd.read_csv(EXPLORATORY_DIR / "diversity_metrics.tsv", sep="\t")
    # Bring in PostNatalAbxCohort (and other covariates) from metadata
    md_cols = ['sample_name', 'PostNatalAbxCohort', 'GestationCohort']
    div = div.merge(metadata[md_cols].drop_duplicates('sample_name'),
                    on='sample_name', how='left')
    # Reduce noise: drop the ~50 samples with no PostNatalAbxCohort if needed
    return div


def diversity_panel(div: pd.DataFrame):
    """Two-row panel: row 1 stratified by Week, row 2 stratified by Abx cohort."""
    metrics = [('richness', 'Gene-family richness'),
               ('total_rpm', 'Total AMR RPM (log10)'),
               ('shannon', 'Shannon diversity')]

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    # Row 1: by site × week
    for j, (col, title) in enumerate(metrics):
        ax = axes[0, j]
        plot_df = div.copy()
        if col == 'total_rpm':
            plot_df['plot_y'] = np.log10(plot_df[col] + 1)
        else:
            plot_df['plot_y'] = plot_df[col]
        sns.boxplot(
            data=plot_df, x='SampleType', y='plot_y',
            hue='SampleCollectionWeek',
            order=SITE_ORDER, hue_order=WEEK_ORDER,
            palette=WEEK_PALETTE, showfliers=False, ax=ax,
        )
        sns.stripplot(
            data=plot_df, x='SampleType', y='plot_y',
            hue='SampleCollectionWeek',
            order=SITE_ORDER, hue_order=WEEK_ORDER,
            dodge=True, size=2.5, alpha=0.4, color='black',
            ax=ax, legend=False,
        )
        ax.set_title(title)
        ax.set_xlabel('')
        ax.set_ylabel('')
        if j == 0:
            ax.legend(title='Week', loc='upper left', frameon=True)
        else:
            ax.legend().remove()

    # Row 2: by site × abx cohort
    for j, (col, title) in enumerate(metrics):
        ax = axes[1, j]
        plot_df = div.dropna(subset=['PostNatalAbxCohort']).copy()
        if col == 'total_rpm':
            plot_df['plot_y'] = np.log10(plot_df[col] + 1)
        else:
            plot_df['plot_y'] = plot_df[col]
        sns.boxplot(
            data=plot_df, x='SampleType', y='plot_y',
            hue='PostNatalAbxCohort',
            order=SITE_ORDER, hue_order=ABX_ORDER,
            palette=ABX_PALETTE, showfliers=False, ax=ax,
        )
        ax.set_title(title)
        ax.set_xlabel('Body site')
        ax.set_ylabel('')
        if j == 0:
            handles, _ = ax.get_legend_handles_labels()
            ax.legend(handles, [ABX_LABELS[a] for a in ABX_ORDER],
                      title='Postnatal infant abx', loc='upper left', frameon=True)
        else:
            ax.legend().remove()

    fig.suptitle('AMR diversity, abundance, and richness by clinical strata',
                 fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    out = FIGURES_DIR / "diversity_by_clinical_strata.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {out}")


def diversity_summary_table(div: pd.DataFrame):
    """Median + IQR for each metric × stratum, written to TSV."""
    rows = []
    for stratifier, order in [('SampleType', SITE_ORDER),
                              ('SampleCollectionWeek', WEEK_ORDER),
                              ('PostNatalAbxCohort', ABX_ORDER)]:
        sub = div.dropna(subset=[stratifier])
        for level in order:
            d = sub[sub[stratifier] == level]
            if d.empty:
                continue
            for metric in ['richness', 'total_rpm', 'shannon']:
                rows.append({
                    'stratifier': stratifier,
                    'level': level,
                    'metric': metric,
                    'n': len(d),
                    'median': float(np.median(d[metric])),
                    'q25': float(np.percentile(d[metric], 25)),
                    'q75': float(np.percentile(d[metric], 75)),
                    'mean': float(d[metric].mean()),
                })
    out_df = pd.DataFrame(rows)
    out_path = V2 / "results" / "clinical" / "diversity_by_strata_summary.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"✓ Saved: {out_path}")


# ---------------------------------------------------------------------------
# Volcano panel from gene-family LMM
# ---------------------------------------------------------------------------

VOLCANO_EFFECTS = [
    ('BodySite[T.Stool]',                      'Stool vs Axilla'),
    ('BodySite[T.Groin]',                      'Groin vs Axilla'),
    ('Week[T.Week.3]',                         'Postnatal Week 3 vs Week 1'),
    ('PostNatalAbxCohort[T.High.Infant.Abx]',  'High infant abx vs None'),
    ('PostNatalAbxCohort[T.Low.Infant.Abx]',   'Low infant abx vs None'),
]


def volcano_panel():
    long = pd.read_csv(GF_DIR / "gene_family_clinical_effects_long.tsv", sep="\t")
    print(f"Loaded gene-family LMM: {long['gene_family'].nunique()} families × "
          f"{long['effect'].nunique()} effects")

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()

    for i, (eff_name, eff_title) in enumerate(VOLCANO_EFFECTS):
        ax = axes[i]
        sub = long[long['effect'] == eff_name].copy()
        if sub.empty:
            ax.set_visible(False)
            continue
        sub['nlp'] = -np.log10(sub['fdr'].clip(lower=1e-300))
        # Color by direction × significance
        colors = []
        for _, r in sub.iterrows():
            if r['fdr'] < 0.05:
                colors.append('#d62728' if r['coef'] > 0 else '#1f77b4')
            else:
                colors.append('#bdbdbd')
        ax.scatter(sub['coef'], sub['nlp'], c=colors, s=22,
                   alpha=0.75, edgecolors='black', linewidths=0.3)

        ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.axvline(0, color='black', linewidth=0.5, alpha=0.5)

        # Annotate top 8 by FDR (split positive/negative for readability)
        sig = sub[sub['fdr'] < 0.05].copy()
        if not sig.empty:
            top_pos = sig[sig['coef'] > 0].nsmallest(5, 'fdr')
            top_neg = sig[sig['coef'] < 0].nsmallest(5, 'fdr')
            for _, r in pd.concat([top_pos, top_neg]).iterrows():
                ax.annotate(r['gene_family'],
                            (r['coef'], r['nlp']),
                            fontsize=7, alpha=0.85,
                            xytext=(2, 2), textcoords='offset points')

        n_up = ((sub['fdr'] < 0.05) & (sub['coef'] > 0)).sum()
        n_dn = ((sub['fdr'] < 0.05) & (sub['coef'] < 0)).sum()
        ax.set_title(f"{eff_title}\n↑{n_up}  ↓{n_dn}  (FDR<0.05)", fontsize=10)
        ax.set_xlabel('LMM coefficient (log10 RPM)')
        ax.set_ylabel('-log10(FDR)')

    # Legend in last (unused) axes slot
    legend_ax = axes[-1]
    legend_ax.axis('off')
    legend_ax.text(0.05, 0.85, 'Legend', fontsize=14, fontweight='bold')
    legend_ax.scatter([0.1], [0.65], color='#d62728', s=80, edgecolor='black', linewidth=0.5)
    legend_ax.text(0.18, 0.65, '↑ in test condition (FDR<0.05)', fontsize=11, va='center')
    legend_ax.scatter([0.1], [0.50], color='#1f77b4', s=80, edgecolor='black', linewidth=0.5)
    legend_ax.text(0.18, 0.50, '↓ in test condition (FDR<0.05)', fontsize=11, va='center')
    legend_ax.scatter([0.1], [0.35], color='#bdbdbd', s=80, edgecolor='black', linewidth=0.5)
    legend_ax.text(0.18, 0.35, 'not significant', fontsize=11, va='center')
    legend_ax.text(0.05, 0.18,
                   "x = LMM coefficient on log10(RPM+1)\n"
                   "Mixed-effects model with subject random intercept;\n"
                   "FDR within each effect (BH).",
                   fontsize=9, style='italic')
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)

    fig.suptitle('Gene-family clinical effects — volcano panel',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    out = FIGURES_DIR / "volcano_clinical_effects.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {out}")


def main():
    print("=" * 60)
    print("CLINICAL DIVERSITY + VOLCANO FIGURES")
    print("=" * 60)
    div = load_diversity()
    diversity_panel(div)
    diversity_summary_table(div)
    volcano_panel()
    print("=" * 60)


if __name__ == "__main__":
    main()
