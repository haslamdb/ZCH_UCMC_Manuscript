#!/usr/bin/env python3
"""
Clinical-category-level mixed-effects models.

Same model as 04b (MEGARes class), 04d (gene_family), 04c (ARO id) — but
the feature axis is now the clinically-defined screening categories
defined in 11_functional_class_screening.py (ESBL, carbapenemase, mecA,
aminoglycoside, fluoroquinolone, MLS, tetracycline, sulfonamide,
trimethoprim, colistin, VRE).

For each category, the per-sample feature is the SUM of rpm across all
genes that match the category. The model regresses log10(sum_rpm + 1)
on body site, week, infant antibiotic cohort, maternal antibiotics,
feeding, and gestational age, with subject random intercept.

Interpretation question: if none of the clinical covariates explain a
category's abundance, that's evidence that the resistome arises from
environmental acquisition (NICU surfaces, staff, etc.) rather than
clinical variables — the Brooks et al. "NICU is the reservoir" model.
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import _clinical_lmm as lmm  # noqa: E402

# Reuse the category definitions from 11_*
from importlib import import_module
_screen = import_module('11_functional_class_screening')
CATEGORIES = _screen.CATEGORIES

warnings.filterwarnings('ignore')

PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
V2 = PROJECT_ROOT / "nicu_resistome_analysis_v2"
DATA_DIR = V2 / "data"
QC_DIR = V2 / "results" / "qc"
OUT_DIR = V2 / "results" / "category_lmm"
FIGURES_DIR = V2 / "figures" / "category_lmm"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def build_category_matrix(metadata: pd.DataFrame, amr: pd.DataFrame) -> pd.DataFrame:
    """Sample × category RPM matrix. For each sample, each cell is the sum
    of rpm across all rows in the AMR file that match the category's
    classifier function. Includes all detected rows (no RPM threshold) —
    threshold-based presence/absence is what 11_* does; here we want the
    continuous abundance signal so the LMM can detect dose-response."""
    print("Building sample × category RPM matrix...")
    samples = metadata['sample_name'].unique()
    amr = amr[amr['sample_name'].isin(samples)].copy()

    cat_columns = {}
    syms = amr['gene_symbol'].values
    subc = amr['subclass'].values
    rpms = amr['rpm'].values
    sample_arr = amr['sample_name'].values

    for cat_name, fn, _ in CATEGORIES:
        mask = np.array([fn(s, c) for s, c in zip(syms, subc)])
        # group rpm by sample within mask
        if mask.sum() == 0:
            cat_columns[cat_name] = pd.Series(0.0, index=samples)
            print(f"  {cat_name:24s} 0 matching rows")
            continue
        sub = pd.DataFrame({'sample_name': sample_arr[mask], 'rpm': rpms[mask]})
        agg = sub.groupby('sample_name')['rpm'].sum()
        cat_columns[cat_name] = agg.reindex(samples, fill_value=0.0)
        print(f"  {cat_name:24s} {mask.sum():>6d} rows, "
              f"{(agg > 0).sum():>4d} samples >0")

    matrix = pd.DataFrame(cat_columns, index=samples)
    matrix.index.name = 'sample_name'
    matrix.to_csv(OUT_DIR / "category_rpm_matrix.tsv", sep="\t")
    print(f"\n  matrix shape: {matrix.shape}")
    return matrix


def main():
    print("=" * 60)
    print("CLINICAL CATEGORY × CLINICAL-COVARIATE LMM")
    print("=" * 60)

    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")
    metadata = metadata[metadata['has_amr_data']].drop_duplicates('sample_name').copy()
    metadata = metadata[metadata['SampleType'].isin(['Axilla', 'Groin', 'Stool'])].copy()

    amr = pd.read_csv(DATA_DIR / "nicu_amr_only.tsv", sep="\t",
                      usecols=['sample_name', 'gene_symbol', 'subclass', 'rpm'])

    matrix = build_category_matrix(metadata, amr)
    model_md = lmm.prep_model_metadata(metadata)
    print(f"\nModel-ready metadata: {len(model_md)} samples, "
          f"{model_md.SubjectID.nunique()} subjects")

    cat_names = [c for c, _, _ in CATEGORIES]
    # Skip categories with too-low prevalence to fit
    selected = lmm.select_features(matrix[cat_names], prevalence_min=0.02, n_min=10)
    skipped = [c for c in cat_names if c not in selected]
    print(f"\n{len(selected)}/{len(cat_names)} categories pass prevalence filter (>=2%, >=10 nonzero)")
    if skipped:
        print(f"  skipped (too rare): {skipped}")

    long_df, wide_df = lmm.fit_all(matrix, model_md, selected, feature_label='category')
    long_df.to_csv(OUT_DIR / "category_clinical_effects_long.tsv", sep="\t", index=False)
    wide_df.to_csv(OUT_DIR / "category_clinical_effects_wide.tsv", sep="\t", index=False)
    print(f"✓ Saved: {OUT_DIR/'category_clinical_effects_long.tsv'}")
    print(f"✓ Saved: {OUT_DIR/'category_clinical_effects_wide.tsv'}")

    lmm.print_summary(long_df)

    # Per-category effect compactness (helps answer the "any clinical signal?" question)
    print("\nPer-category: number of clinical effects with FDR<0.05")
    by_cat = (long_df[long_df['fdr'] < 0.05].groupby('category').size()
              .reindex(selected, fill_value=0).sort_values(ascending=False))
    for cat, n in by_cat.items():
        print(f"  {cat:24s} {n:>2d} significant clinical effects")

    lmm.write_topn(long_df, OUT_DIR / "category_clinical_effects_topN.tsv",
                   feature_label='category', n=11)
    lmm.heatmap(long_df, FIGURES_DIR / "category_clinical_effects_heatmap.pdf",
                feature_label='category',
                title='Clinical-covariate effects on AMR-category abundance\n'
                      '(per-category linear mixed-effects, SubjectID random intercept)')

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
