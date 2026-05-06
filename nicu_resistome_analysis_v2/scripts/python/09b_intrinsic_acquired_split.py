#!/usr/bin/env python3
"""
Intrinsic vs acquired resistance — split headline tables.

`resistance_type` (from the NCBI ReferenceGeneCatalog annotation joined
in by the aggregator) labels each gene_family as `intrinsic` (e.g.,
chromosomal Pseudomonas/Klebsiella efflux pumps), `acquired` (e.g.,
mobilizable bla genes), or NaN (not annotated by NCBI). For clinical
interpretation these populations shouldn't share a denominator —
intrinsic resistance is a feature of the organism, acquired resistance
reflects mobilizable selection.

This script produces a per-sample summary stratified by resistance
type, plus group-level tables by Location × SampleType × Week, and a
small panel figure. Outputs feed into the section that 09's headline
report lacks.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
DATA_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "data"
QC_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results" / "qc"
SUMMARY_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results" / "summary"
FIGURES_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "figures" / "summary"

SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
sns.set_context("talk")

RTYPE_ORDER = ['intrinsic', 'acquired', 'unknown']
RTYPE_PALETTE = {'intrinsic': '#7f7f7f', 'acquired': '#d62728', 'unknown': '#bcbcbc'}


def load_data():
    print("Loading data...")
    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")
    amr = pd.read_csv(DATA_DIR / "nicu_amr_only.tsv", sep="\t",
                      usecols=['sample_name', 'gene_family', 'resistance_type', 'rpm'])
    amr['resistance_type'] = amr['resistance_type'].fillna('unknown')

    metadata = metadata[metadata['has_amr_data']].copy()
    metadata = metadata.drop_duplicates(subset='sample_name', keep='first').copy()
    metadata = metadata[metadata['SampleType'].isin(['Axilla', 'Groin', 'Stool'])].copy()

    print(f"  AMR rows: {len(amr):,}")
    print(f"  Samples (deduped, valid sites): {len(metadata):,}")
    print("  resistance_type distribution (rows):")
    print(amr['resistance_type'].value_counts().to_string())
    return metadata, amr


def per_sample_summary(metadata, amr):
    """One row per sample × resistance_type with rpm sum and gene count."""
    print("\nBuilding per-sample × resistance_type summary...")
    samples = metadata['sample_name'].unique()
    amr = amr[amr['sample_name'].isin(samples)]

    g = amr.groupby(['sample_name', 'resistance_type']).agg(
        rpm_sum=('rpm', 'sum'),
        n_gene_families=('gene_family', 'nunique'),
    ).reset_index()

    # Reindex so every (sample, type) combination exists, with zeros
    full_idx = pd.MultiIndex.from_product([samples, RTYPE_ORDER],
                                          names=['sample_name', 'resistance_type'])
    g = g.set_index(['sample_name', 'resistance_type']).reindex(full_idx, fill_value=0).reset_index()

    g = g.merge(
        metadata[['sample_name', 'SubjectID', 'Location', 'SampleType',
                  'SampleCollectionWeek']],
        on='sample_name', how='left',
    )

    g.to_csv(SUMMARY_DIR / "intrinsic_vs_acquired_per_sample.tsv",
             sep="\t", index=False)
    print(f"✓ Saved: {SUMMARY_DIR/'intrinsic_vs_acquired_per_sample.tsv'}")
    return g


def group_summary(per_sample):
    """Group-level summary: median + IQR rpm by Location × SampleType ×
    Week × resistance_type, plus simple counts."""
    print("\nBuilding group-level summary table...")
    grp_cols = ['Location', 'SampleType', 'SampleCollectionWeek', 'resistance_type']

    summary = per_sample.groupby(grp_cols).agg(
        n_samples=('sample_name', 'nunique'),
        median_rpm=('rpm_sum', 'median'),
        q25_rpm=('rpm_sum', lambda x: np.percentile(x, 25)),
        q75_rpm=('rpm_sum', lambda x: np.percentile(x, 75)),
        mean_rpm=('rpm_sum', 'mean'),
        median_gene_families=('n_gene_families', 'median'),
        prevalence_any_detected=('rpm_sum', lambda x: (x > 0).mean() * 100),
    ).reset_index()

    # Order
    summary['resistance_type'] = pd.Categorical(
        summary['resistance_type'], categories=RTYPE_ORDER, ordered=True)
    summary = summary.sort_values(grp_cols)

    summary.to_csv(SUMMARY_DIR / "intrinsic_vs_acquired_summary.tsv",
                   sep="\t", index=False)
    print(f"✓ Saved: {SUMMARY_DIR/'intrinsic_vs_acquired_summary.tsv'}")

    # Compact print
    pivot = summary.pivot_table(
        index=['SampleType', 'resistance_type'],
        columns=['Location', 'SampleCollectionWeek'],
        values='median_rpm',
    )
    print("\nMedian RPM by Site × Type × Location × Week:")
    print(pivot.round(1).to_string())
    return summary


def fraction_acquired(per_sample):
    """Per-sample acquired/(acquired+intrinsic) ratio — useful clinical metric."""
    print("\nComputing per-sample acquired fraction...")
    p = per_sample.pivot_table(index='sample_name',
                               columns='resistance_type',
                               values='rpm_sum',
                               fill_value=0).reset_index()
    for col in RTYPE_ORDER:
        if col not in p.columns:
            p[col] = 0.0
    denom = p['acquired'] + p['intrinsic']
    p['acquired_fraction'] = np.where(denom > 0, p['acquired'] / denom, np.nan)

    md_keep = per_sample[['sample_name', 'Location', 'SampleType',
                          'SampleCollectionWeek', 'SubjectID']].drop_duplicates()
    p = p.merge(md_keep, on='sample_name', how='left')

    p.to_csv(SUMMARY_DIR / "acquired_fraction_per_sample.tsv",
             sep="\t", index=False)
    print(f"✓ Saved: {SUMMARY_DIR/'acquired_fraction_per_sample.tsv'}")

    print("\nMedian acquired fraction by SampleType × Week:")
    pf = p.groupby(['SampleType', 'SampleCollectionWeek'])['acquired_fraction'].median().unstack()
    print(pf.round(3).to_string())
    return p


def plot_panel(per_sample, frac_df):
    print("\nPlotting...")
    fig, axes = plt.subplots(1, 3, figsize=(22, 6))

    # Panel A: stacked bar of median RPM by SampleType × Week × type
    pivot = (per_sample.groupby(['SampleType', 'SampleCollectionWeek',
                                 'resistance_type'])['rpm_sum'].median()
             .unstack('resistance_type'))
    pivot = pivot.reindex(columns=RTYPE_ORDER, fill_value=0)
    pivot.plot(kind='bar', stacked=True, ax=axes[0],
               color=[RTYPE_PALETTE[t] for t in pivot.columns], edgecolor='black')
    axes[0].set_ylabel('Median RPM (per sample)')
    axes[0].set_title('A. Median resistome RPM\nby resistance type')
    axes[0].set_xlabel('SampleType, Week')
    axes[0].legend(title='Resistance type', loc='upper left')

    # Panel B: boxplot of log10 RPM — intrinsic vs acquired by SampleType
    df = per_sample[per_sample['resistance_type'].isin(['intrinsic', 'acquired'])].copy()
    df['log_rpm'] = np.log10(df['rpm_sum'] + 1)
    sns.boxplot(data=df, x='SampleType', y='log_rpm',
                hue='resistance_type', order=['Axilla', 'Groin', 'Stool'],
                hue_order=['intrinsic', 'acquired'],
                palette={'intrinsic': RTYPE_PALETTE['intrinsic'],
                         'acquired':  RTYPE_PALETTE['acquired']},
                ax=axes[1])
    axes[1].set_ylabel('log10(RPM + 1)')
    axes[1].set_xlabel('Body site')
    axes[1].set_title('B. Intrinsic vs acquired distribution\nby body site')

    # Panel C: acquired fraction by SampleType × Week
    sns.boxplot(data=frac_df, x='SampleType', y='acquired_fraction',
                hue='SampleCollectionWeek',
                order=['Axilla', 'Groin', 'Stool'],
                hue_order=['Week.1', 'Week.3'],
                palette={'Week.1': '#2ca02c', 'Week.3': '#d62728'},
                ax=axes[2])
    axes[2].set_ylabel('Acquired / (acquired + intrinsic)')
    axes[2].set_xlabel('Body site')
    axes[2].set_title('C. Acquired fraction\nby body site × week')
    axes[2].set_ylim(0, 1.05)

    plt.tight_layout()
    out = FIGURES_DIR / "intrinsic_vs_acquired_panel.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {out}")


def main():
    print("=" * 60)
    print("INTRINSIC vs ACQUIRED RESISTANCE SPLIT")
    print("=" * 60)
    metadata, amr = load_data()
    per_sample = per_sample_summary(metadata, amr)
    group_summary(per_sample)
    frac_df = fraction_acquired(per_sample)
    plot_panel(per_sample, frac_df)
    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
