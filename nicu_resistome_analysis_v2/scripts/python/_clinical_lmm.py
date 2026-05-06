"""
Shared helpers for clinical-covariate per-feature mixed-effects analysis.

Used by 04b (MEGARes class), 04d (gene_family), 04c (ARO id). All three run
the same model — log10(rpm + 1) regressed on the clinical covariates we
care about, with SubjectID random intercept — and only differ in which
feature axis they pivot on.

The fixed effects are the manuscript's headline contrasts: body site,
postnatal age (week), infant antibiotic cohort, maternal antibiotics,
feeding (donor vs mother's milk), and gestational age. Reference levels
are chosen so positive coefficients have a clean clinical interpretation
("more AMR with X" rather than "less AMR without X").
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')


# (statsmodels coefficient name, human-readable label)
EFFECTS: list[tuple[str, str]] = [
    ('BodySite[T.Groin]',                      'Groin (vs Axilla)'),
    ('BodySite[T.Stool]',                      'Stool (vs Axilla)'),
    ('Week[T.Week.3]',                         'Postnatal Week 3 (vs Week 1)'),
    ('PostNatalAbxCohort[T.Low.Infant.Abx]',   'Low infant abx (vs None)'),
    ('PostNatalAbxCohort[T.High.Infant.Abx]',  'High infant abx (vs None)'),
    ('MaternalAntibiotics[T.Mat.Abx]',         'Maternal abx (vs None)'),
    ('AnyMilk[T.Donor]',                       "Donor milk (vs Mother's)"),
    ('GestationTime',                          'Gestational age (per week)'),
]

FORMULA = ("y ~ BodySite + Week + PostNatalAbxCohort"
           " + MaternalAntibiotics + AnyMilk + GestationTime")


def prep_model_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """Filter and categoricalize covariates. Drops rows with any covariate NaN
    so all features see the same N (avoids spurious comparisons across
    features with different missingness patterns)."""
    md = metadata.copy()
    md = md.dropna(subset=[
        'SampleType', 'SampleCollectionWeek', 'PostNatalAbxCohort',
        'MaternalAntibiotics', 'AnyMilk', 'GestationTime', 'SubjectID',
    ]).copy()
    md['BodySite'] = pd.Categorical(md['SampleType'], categories=['Axilla', 'Groin', 'Stool'])
    md['Week']     = pd.Categorical(md['SampleCollectionWeek'], categories=['Week.1', 'Week.3'])
    md['PostNatalAbxCohort'] = pd.Categorical(
        md['PostNatalAbxCohort'],
        categories=['No.Infant.Abx', 'Low.Infant.Abx', 'High.Infant.Abx'],
    )
    md['MaternalAntibiotics'] = pd.Categorical(
        md['MaternalAntibiotics'],
        categories=['None.Mat.Abx', 'Mat.Abx'],
    )
    md['AnyMilk'] = pd.Categorical(md['AnyMilk'], categories=['Mother', 'Donor'])
    md['GestationTime'] = pd.to_numeric(md['GestationTime'], errors='coerce')
    md = md.dropna(subset=['GestationTime'])
    return md


def select_features(matrix: pd.DataFrame, prevalence_min: float = 0.05,
                    n_min: int = 20) -> list[str]:
    """Keep features detected in >= prevalence_min of samples AND with at
    least n_min nonzero values."""
    prev = (matrix > 0).mean()
    n_nonzero = (matrix > 0).sum()
    keep = (prev >= prevalence_min) & (n_nonzero >= n_min)
    return matrix.columns[keep].tolist()


def fit_one(feature: str, matrix: pd.DataFrame, model_md: pd.DataFrame) -> dict | None:
    """Fit one mixed-effects model for one feature. Returns per-effect
    coef/se/pvalue plus convergence + N. Returns None on hard failure."""
    if feature not in matrix.columns:
        return None
    y = np.log10(matrix[feature].reindex(model_md['sample_name']).values + 1.0)
    df = model_md.copy()
    df['y'] = y
    df = df.dropna(subset=['y'])
    if len(df) < 50:
        return None
    try:
        model = smf.mixedlm(FORMULA, df, groups=df['SubjectID'])
        result = model.fit(reml=True, method='powell', disp=False)
    except Exception as exc:
        return {'_error': str(exc)}
    out = {}
    for eff_name, _ in EFFECTS:
        if eff_name in result.params.index:
            out[eff_name] = {
                'coef':   float(result.params[eff_name]),
                'se':     float(result.bse[eff_name]),
                'pvalue': float(result.pvalues[eff_name]),
            }
        else:
            out[eff_name] = {'coef': np.nan, 'se': np.nan, 'pvalue': np.nan}
    out['_n'] = len(df)
    out['_converged'] = bool(result.converged) if hasattr(result, 'converged') else True
    return out


def fit_all(matrix: pd.DataFrame, model_md: pd.DataFrame,
            features: list[str], feature_label: str = 'feature') -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run fit_one over all features, then FDR-correct per effect across
    features. Returns (long_df, wide_df). long_df has one row per
    (feature, effect); wide_df has one row per feature with effect columns."""
    rows_long = []
    rows_wide = []
    n_skipped = 0
    for feat in features:
        res = fit_one(feat, matrix, model_md)
        if res is None or '_error' in res:
            n_skipped += 1
            continue
        wide = {feature_label: feat, 'n': res['_n'], 'converged': res['_converged']}
        for eff_name, eff_label in EFFECTS:
            stats = res[eff_name]
            wide[f'{eff_name}_coef']   = stats['coef']
            wide[f'{eff_name}_pvalue'] = stats['pvalue']
            rows_long.append({
                feature_label: feat,
                'effect': eff_name,
                'effect_label': eff_label,
                'coef':   stats['coef'],
                'se':     stats['se'],
                'pvalue': stats['pvalue'],
            })
        rows_wide.append(wide)

    long_df = pd.DataFrame(rows_long)
    wide_df = pd.DataFrame(rows_wide)
    print(f"  fit {len(rows_wide)} / {len(features)} features ({n_skipped} skipped)")

    # FDR within each effect across features
    long_df['fdr'] = np.nan
    for eff_name, _ in EFFECTS:
        m = long_df['effect'] == eff_name
        pvals = long_df.loc[m, 'pvalue'].values
        valid = ~np.isnan(pvals)
        if valid.sum() == 0:
            continue
        fdr = np.full_like(pvals, np.nan, dtype=float)
        fdr[valid] = multipletests(pvals[valid], method='fdr_bh')[1]
        long_df.loc[m, 'fdr'] = fdr

    long_df = long_df.sort_values(['effect', 'fdr', 'pvalue'])
    return long_df, wide_df


def print_summary(long_df: pd.DataFrame, fdr_threshold: float = 0.05) -> None:
    print(f"\nSignificant per-effect counts (FDR<{fdr_threshold}):")
    for eff_name, eff_label in EFFECTS:
        sig = long_df[(long_df['effect'] == eff_name) & (long_df['fdr'] < fdr_threshold)]
        n_pos = (sig['coef'] > 0).sum()
        n_neg = (sig['coef'] < 0).sum()
        print(f"  {eff_label:38s}  {len(sig):4d}  (+{n_pos} / -{n_neg})")


def write_topn(long_df: pd.DataFrame, output_path: Path,
               feature_label: str = 'feature', n: int = 20) -> None:
    rows = []
    for eff_name, eff_label in EFFECTS:
        sub = long_df[long_df['effect'] == eff_name].sort_values('fdr').head(n)
        for _, r in sub.iterrows():
            rows.append({
                'effect_label': eff_label,
                feature_label:  r[feature_label],
                'coef':   r['coef'],
                'pvalue': r['pvalue'],
                'fdr':    r['fdr'],
            })
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)
    print(f"✓ Saved: {output_path}")


def heatmap(long_df: pd.DataFrame, output_path: Path,
            feature_label: str = 'feature', title: str = '',
            top_features_per_effect: int = 0) -> None:
    """Heatmap: features (rows) × effects (cols), color = signed
    -log10(fdr). If top_features_per_effect > 0, show the top-N features
    per effect (union); else show all features."""
    if long_df.empty:
        return
    pivot_fdr = long_df.pivot_table(index=feature_label, columns='effect_label',
                                    values='fdr', aggfunc='first')
    pivot_coef = long_df.pivot_table(index=feature_label, columns='effect_label',
                                     values='coef', aggfunc='first')
    sig = -np.log10(pivot_fdr.clip(lower=1e-15))
    sig = sig.where(pivot_coef >= 0, -sig)
    sig = sig.clip(-10, 10)

    if top_features_per_effect > 0:
        keep = set()
        for eff_name, eff_label in EFFECTS:
            if eff_label in pivot_fdr.columns:
                top = pivot_fdr[eff_label].dropna().sort_values().head(top_features_per_effect).index
                keep.update(top)
        sig = sig.loc[sig.index.intersection(keep)]

    order = sig.abs().max(axis=1).sort_values(ascending=False).index
    sig = sig.loc[order]
    col_order = [lbl for _, lbl in EFFECTS if lbl in sig.columns]
    sig = sig[col_order]

    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, max(6, 0.28 * len(sig) + 2)))
    sns.heatmap(sig, cmap='RdBu_r', center=0, ax=ax,
                cbar_kws={'label': 'sign(coef) × -log10(FDR)'},
                linewidths=0.3, linecolor='white')
    ax.set_xlabel('')
    ax.set_ylabel(feature_label)
    if title:
        ax.set_title(title)
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {output_path}")
