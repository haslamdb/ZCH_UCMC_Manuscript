#!/usr/bin/env python3
"""
Differential abundance at the GENE_FAMILY level — clinical covariates.

Companion to 04b (MEGARes class) and 04c (ARO id). Same model, same
covariates, different feature axis: the manuscript-level gene_family
resolution (e.g. all blaCTX-M alleles collapse to a single blaCTX-M row).

Replaces the location-only Mann-Whitney framing in 04 with a
mixed-effects model that captures the contrasts the manuscript actually
cares about: body site, postnatal age, infant antibiotic cohort,
maternal antibiotics, feeding, gestational age.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import _clinical_lmm as lmm  # noqa: E402

PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
DATA_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "data"
QC_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results" / "qc"
EXPLORATORY_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results" / "exploratory"
DIFF_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results" / "differential_gene_family"
FIGURES_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "figures" / "differential_gene_family"

DIFF_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    print("Loading data...")
    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")
    # Use the QC-filtered gene-family RPM matrix that 03 produces — same
    # universe of 250 gene families used by 04 and 06.
    matrix = pd.read_csv(EXPLORATORY_DIR / "amr_rpm_matrix.tsv",
                         sep="\t", index_col=0)
    print(f"  AMR matrix: {matrix.shape}")

    metadata = metadata[metadata['has_amr_data']].copy()
    metadata = metadata.drop_duplicates(subset='sample_name', keep='first').copy()
    metadata = metadata[metadata['SampleType'].isin(['Axilla', 'Groin', 'Stool'])].copy()
    print(f"  Metadata samples (deduped, valid sites): {len(metadata):,}")
    return metadata, matrix


def main():
    print("=" * 60)
    print("GENE-FAMILY DIFFERENTIAL ABUNDANCE — clinical covariates")
    print("=" * 60)
    metadata, matrix = load_data()
    model_md = lmm.prep_model_metadata(metadata)
    print(f"  Model-ready metadata: {len(model_md)} samples, "
          f"{model_md.SubjectID.nunique()} subjects")

    families = lmm.select_features(matrix)
    print(f"  {len(families)}/{matrix.shape[1]} gene families pass prevalence filter")

    long_df, wide_df = lmm.fit_all(matrix, model_md, families,
                                   feature_label='gene_family')

    long_df.to_csv(DIFF_DIR / "gene_family_clinical_effects_long.tsv",
                   sep="\t", index=False)
    wide_df.to_csv(DIFF_DIR / "gene_family_clinical_effects_wide.tsv",
                   sep="\t", index=False)
    print(f"✓ Saved: {DIFF_DIR/'gene_family_clinical_effects_long.tsv'}")
    print(f"✓ Saved: {DIFF_DIR/'gene_family_clinical_effects_wide.tsv'}")

    lmm.print_summary(long_df)
    lmm.write_topn(long_df, DIFF_DIR / "gene_family_clinical_effects_topN.tsv",
                   feature_label='gene_family', n=30)
    # Heatmap shows the union of top-15 per effect — too many genes to plot all
    lmm.heatmap(long_df, FIGURES_DIR / "gene_family_clinical_effects_heatmap.pdf",
                feature_label='gene_family',
                title='Clinical-covariate effects on gene-family abundance\n'
                      '(per-family linear mixed-effects, SubjectID random intercept; '
                      'top 15 per effect)',
                top_features_per_effect=15)

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
