#!/usr/bin/env python3
"""
Differential abundance at the MEGARes-CLASS level — clinical covariates.

For each MEGARes class with sufficient prevalence, fit a linear
mixed-effects model (see _clinical_lmm.py for the formula and effect
definitions). FDR-corrected per covariate across classes.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import _clinical_lmm as lmm  # noqa: E402

PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
DATA_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "data"
QC_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results" / "qc"
DIFF_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results" / "differential_megares"
FIGURES_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "figures" / "differential_megares"

DIFF_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    print("Loading data...")
    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")
    amr = pd.read_csv(DATA_DIR / "nicu_amr_only.tsv", sep="\t",
                      usecols=['sample_name', 'megares_class', 'rpm'])
    n_before = len(amr)
    amr = amr.dropna(subset=['megares_class']).copy()
    print(f"  AMR rows after drop NaN class: {len(amr):,} (dropped {n_before-len(amr):,})")

    metadata = metadata[metadata['has_amr_data']].copy()
    metadata = metadata.drop_duplicates(subset='sample_name', keep='first').copy()
    metadata = metadata[metadata['SampleType'].isin(['Axilla', 'Groin', 'Stool'])].copy()
    print(f"  Metadata samples (deduped, valid sites): {len(metadata):,}")
    return metadata, amr


def build_class_matrix(metadata, amr):
    samples = metadata['sample_name'].unique()
    amr = amr[amr['sample_name'].isin(samples)]
    matrix = amr.pivot_table(
        index='sample_name', columns='megares_class',
        values='rpm', aggfunc='sum', fill_value=0,
    )
    matrix = matrix.reindex(samples, fill_value=0)
    print(f"  sample × MEGARes-class matrix: {matrix.shape}")
    matrix.to_csv(DIFF_DIR / "megares_class_rpm_matrix.tsv", sep="\t")
    return matrix


def main():
    print("=" * 60)
    print("MEGARes-CLASS DIFFERENTIAL ABUNDANCE — clinical covariates")
    print("=" * 60)
    metadata, amr = load_data()
    matrix = build_class_matrix(metadata, amr)
    model_md = lmm.prep_model_metadata(metadata)
    print(f"  Model-ready metadata: {len(model_md)} samples, "
          f"{model_md.SubjectID.nunique()} subjects")

    classes = lmm.select_features(matrix)
    print(f"  {len(classes)}/{matrix.shape[1]} classes pass prevalence filter")

    long_df, wide_df = lmm.fit_all(matrix, model_md, classes,
                                   feature_label='megares_class')

    long_df.to_csv(DIFF_DIR / "megares_class_clinical_effects_long.tsv",
                   sep="\t", index=False)
    wide_df.to_csv(DIFF_DIR / "megares_class_clinical_effects_wide.tsv",
                   sep="\t", index=False)
    print(f"✓ Saved: {DIFF_DIR/'megares_class_clinical_effects_long.tsv'}")
    print(f"✓ Saved: {DIFF_DIR/'megares_class_clinical_effects_wide.tsv'}")

    lmm.print_summary(long_df)
    lmm.write_topn(long_df, DIFF_DIR / "megares_class_clinical_effects_topN.tsv",
                   feature_label='megares_class', n=20)
    lmm.heatmap(long_df, FIGURES_DIR / "megares_class_clinical_effects_heatmap.pdf",
                feature_label='megares_class',
                title='Clinical-covariate effects on MEGARes-class abundance\n'
                      '(per-class linear mixed-effects, SubjectID random intercept)')

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
