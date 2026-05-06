#!/usr/bin/env python3
"""
ARO-level comparison panel — cross-referenceable to AMRFinderPlus / CARD.

The AMR data carry an `aro_id` (CARD ARO accession, joined in by the
aggregator) for the ~50% of detected gene families that CARD annotates.
This script builds a per-ARO comparison panel that an outside reader
can pull up directly in CARD or AMRFinderPlus output:

  - per-ARO median RPM and prevalence stratified by SampleType,
    SampleCollectionWeek, and PostNatalAbxCohort
  - retain `subclass`, `gene_family`, and `aro_match_method` for context
  - heatmap of top ARO entries by clinical contrast

We use the same metadata clean-up as 04b (drop NaNs in covariates,
restrict to Axilla/Groin/Stool, no Location-vs-Location framing).
"""

import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent))
import _clinical_lmm as lmm  # noqa: E402

warnings.filterwarnings('ignore')

PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
DATA_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "data"
QC_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results" / "qc"
ARO_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results" / "aro_panel"
FIGURES_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "figures" / "aro_panel"

ARO_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
sns.set_context("talk")


def load_data():
    print("Loading data...")
    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")
    cols = ['sample_name', 'aro_id', 'aro_match_method',
            'gene_family', 'subclass', 'resistance_type',
            'megares_class', 'rpm']
    amr = pd.read_csv(DATA_DIR / "nicu_amr_only.tsv", sep="\t", usecols=cols)

    # Restrict to rows that have an ARO annotation (~50%).
    n_before = len(amr)
    amr = amr.dropna(subset=['aro_id']).copy()
    print(f"  Rows with ARO ID: {len(amr):,} / {n_before:,}")

    metadata = metadata[metadata['has_amr_data']].copy()
    metadata = metadata.drop_duplicates(subset='sample_name', keep='first').copy()
    metadata = metadata[metadata['SampleType'].isin(['Axilla', 'Groin', 'Stool'])].copy()
    print(f"  Metadata samples: {len(metadata):,}")
    return metadata, amr


def build_aro_matrix(metadata, amr):
    """Sample × ARO matrix of summed RPM."""
    print("\nBuilding sample × ARO matrix...")
    samples = metadata['sample_name'].unique()
    amr = amr[amr['sample_name'].isin(samples)]

    matrix = amr.pivot_table(
        index='sample_name', columns='aro_id',
        values='rpm', aggfunc='sum', fill_value=0,
    )
    matrix = matrix.reindex(samples, fill_value=0)
    print(f"  Matrix: {matrix.shape}")
    matrix.to_csv(ARO_DIR / "aro_rpm_matrix.tsv", sep="\t")
    return matrix


def aro_lookup(amr):
    """One row per ARO with cross-reference fields."""
    lookup = (amr.groupby('aro_id')
              .agg(gene_family=('gene_family', lambda x: ';'.join(sorted(set(x.dropna())))),
                   subclass=('subclass', lambda x: ';'.join(sorted(set(x.dropna())))),
                   resistance_type=('resistance_type', lambda x: ';'.join(sorted(set(x.dropna())))),
                   megares_class=('megares_class', lambda x: ';'.join(sorted(set(x.dropna())))),
                   aro_match_method=('aro_match_method', lambda x: ';'.join(sorted(set(x.dropna())))),
                   n_rows=('rpm', 'size'))
              .reset_index())
    lookup.to_csv(ARO_DIR / "aro_lookup.tsv", sep="\t", index=False)
    print(f"✓ Saved: {ARO_DIR/'aro_lookup.tsv'} ({len(lookup)} ARO entries)")
    return lookup


def stratified_summary(matrix, metadata):
    """Per-ARO median RPM + prevalence, stratified by clinical contrasts."""
    print("\nBuilding stratified panel...")

    # Wide → long for easy groupby
    long = matrix.reset_index().melt(id_vars='sample_name', var_name='aro_id', value_name='rpm')
    md_keep = metadata[['sample_name', 'SampleType', 'SampleCollectionWeek',
                        'PostNatalAbxCohort']].copy()
    long = long.merge(md_keep, on='sample_name', how='left')

    rows = []
    for keys, sub in long.groupby(['aro_id']):
        for site in ['Axilla', 'Groin', 'Stool']:
            for week in ['Week.1', 'Week.3']:
                m = (sub['SampleType'] == site) & (sub['SampleCollectionWeek'] == week)
                vals = sub.loc[m, 'rpm']
                if len(vals) == 0:
                    continue
                rows.append({
                    'aro_id': keys if isinstance(keys, str) else keys[0],
                    'SampleType': site,
                    'SampleCollectionWeek': week,
                    'n_samples': len(vals),
                    'prevalence_pct': float((vals > 0).mean() * 100),
                    'median_rpm': float(vals.median()),
                    'mean_rpm': float(vals.mean()),
                })

    panel = pd.DataFrame(rows)
    panel.to_csv(ARO_DIR / "aro_panel_by_site_week.tsv", sep="\t", index=False)
    print(f"✓ Saved: {ARO_DIR/'aro_panel_by_site_week.tsv'} ({len(panel)} rows)")

    # Same summary stratified by infant antibiotic cohort
    rows_abx = []
    for keys, sub in long.groupby(['aro_id']):
        for cohort in ['No.Infant.Abx', 'Low.Infant.Abx', 'High.Infant.Abx']:
            m = sub['PostNatalAbxCohort'] == cohort
            vals = sub.loc[m, 'rpm']
            if len(vals) == 0:
                continue
            rows_abx.append({
                'aro_id': keys if isinstance(keys, str) else keys[0],
                'PostNatalAbxCohort': cohort,
                'n_samples': len(vals),
                'prevalence_pct': float((vals > 0).mean() * 100),
                'median_rpm': float(vals.median()),
                'mean_rpm': float(vals.mean()),
            })

    panel_abx = pd.DataFrame(rows_abx)
    panel_abx.to_csv(ARO_DIR / "aro_panel_by_abx_cohort.tsv", sep="\t", index=False)
    print(f"✓ Saved: {ARO_DIR/'aro_panel_by_abx_cohort.tsv'} ({len(panel_abx)} rows)")
    return panel, panel_abx


def heatmap_top_aro(matrix, metadata, lookup, top_n=30):
    """Heatmap of top-N ARO entries by overall median RPM, ordered by site."""
    print(f"\nHeatmap of top {top_n} ARO entries...")

    # Rank ARO by median RPM across all samples
    overall_median = matrix.median(axis=0).sort_values(ascending=False)
    top = overall_median.head(top_n).index.tolist()

    # Compose row labels: aro_id + most-common gene_family
    labels = {}
    lookup_idx = lookup.set_index('aro_id')
    for aro in top:
        family = lookup_idx.loc[aro, 'gene_family'] if aro in lookup_idx.index else ''
        # Just keep first gene family if multi
        family = family.split(';')[0] if family else ''
        labels[aro] = f"{aro} ({family})" if family else aro

    # Median per group: SampleType × Week
    md_sub = metadata.set_index('sample_name')[['SampleType', 'SampleCollectionWeek']]
    md_sub = md_sub[md_sub.index.isin(matrix.index)]
    sub_matrix = matrix.loc[md_sub.index, top]

    grouped = sub_matrix.join(md_sub).groupby(['SampleType', 'SampleCollectionWeek'])[top].median()
    grouped = np.log10(grouped + 0.1)

    fig, ax = plt.subplots(figsize=(max(10, 0.4 * top_n + 4), 6))
    sns.heatmap(grouped.T, cmap='viridis', ax=ax,
                cbar_kws={'label': 'log10(median RPM + 0.1)'},
                yticklabels=[labels[a] for a in top],
                linewidths=0.4, linecolor='white')
    ax.set_xlabel('SampleType, SampleCollectionWeek')
    ax.set_ylabel('ARO')
    ax.set_title(f'Top {top_n} ARO entries — median RPM by body site × week')
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    out = FIGURES_DIR / f"aro_top{top_n}_heatmap.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {out}")


def clinical_lmm_panel(matrix, metadata, lookup):
    """Per-ARO mixed-effects model with the same clinical covariates as
    04b/04d. Lets readers ask 'which ARO entries shift with body site /
    week / abx / etc.', not just 'where do they sit on average'."""
    print("\n" + "=" * 60)
    print("ARO clinical-covariate mixed-effects models")
    print("=" * 60)
    model_md = lmm.prep_model_metadata(metadata)
    print(f"  Model-ready metadata: {len(model_md)} samples, "
          f"{model_md.SubjectID.nunique()} subjects")
    aros = lmm.select_features(matrix)
    print(f"  {len(aros)}/{matrix.shape[1]} ARO entries pass prevalence filter")

    long_df, wide_df = lmm.fit_all(matrix, model_md, aros, feature_label='aro_id')

    # Annotate long_df with gene_family + subclass for readability
    annot = lookup.set_index('aro_id')[['gene_family', 'subclass', 'megares_class']]
    long_df = long_df.merge(annot, left_on='aro_id', right_index=True, how='left')
    wide_df = wide_df.merge(annot, left_on='aro_id', right_index=True, how='left')

    long_df.to_csv(ARO_DIR / "aro_clinical_effects_long.tsv", sep="\t", index=False)
    wide_df.to_csv(ARO_DIR / "aro_clinical_effects_wide.tsv", sep="\t", index=False)
    print(f"✓ Saved: {ARO_DIR/'aro_clinical_effects_long.tsv'}")
    print(f"✓ Saved: {ARO_DIR/'aro_clinical_effects_wide.tsv'}")

    lmm.print_summary(long_df)
    lmm.write_topn(long_df, ARO_DIR / "aro_clinical_effects_topN.tsv",
                   feature_label='aro_id', n=30)
    lmm.heatmap(long_df, FIGURES_DIR / "aro_clinical_effects_heatmap.pdf",
                feature_label='aro_id',
                title='Clinical-covariate effects on ARO-level abundance\n'
                      '(per-ARO linear mixed-effects, SubjectID random intercept; '
                      'top 15 per effect)',
                top_features_per_effect=15)


def main():
    print("=" * 60)
    print("ARO-LEVEL COMPARISON PANEL")
    print("=" * 60)
    metadata, amr = load_data()
    matrix = build_aro_matrix(metadata, amr)
    lookup = aro_lookup(amr)
    panel, panel_abx = stratified_summary(matrix, metadata)
    heatmap_top_aro(matrix, metadata, lookup, top_n=30)

    # Clinical-covariate inferential layer (mirrors 04b/04d)
    clinical_lmm_panel(matrix, metadata, lookup)

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
