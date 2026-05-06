#!/usr/bin/env python3
"""
Clinical-screening view: per-subject × per-body-site presence/absence of
functional AMR categories (ESBL, carbapenemases, mecA, etc.).

The framing is "metagenomics as a screening swab" — for each clinical
category, what fraction of infants would test positive if you sampled
this body site? Differs by body site, and that's the point: stool
swabs miss MRSA, skin swabs miss ESBL, groin may be a sweet spot.

Threshold: a sample is "positive" for a category if at least one gene
matching the category has rpm >= 1 in `data/nicu_amr_only.tsv` (which
already passes HIGH-confidence + >=80% coverage + >=90% identity QC).

A subject × body-site is positive if ANY of that subject's samples at
that body site (across both timepoints) is positive.

Outputs:
- `results/clinical_screening/functional_class_definitions.tsv`
- `results/clinical_screening/functional_class_per_subject_site.tsv`
- `results/clinical_screening/functional_class_prevalence_by_site.tsv`
- Figures: prevalence-by-site bar panel; Venn per category; site-only
  vs combined-site capture comparison.
"""

import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn3, venn3_unweighted

warnings.filterwarnings('ignore')

PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
V2 = PROJECT_ROOT / "nicu_resistome_analysis_v2"
DATA_DIR = V2 / "data"
QC_DIR = V2 / "results" / "qc"
OUT_DIR = V2 / "results" / "clinical_screening"
FIGURES_DIR = V2 / "figures" / "clinical_screening"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.05)

SITE_ORDER = ['Axilla', 'Groin', 'Stool']
SITE_PALETTE = {'Axilla': '#1f77b4', 'Groin': '#ff7f0e', 'Stool': '#2ca02c'}

PRESENCE_RPM = 1.0


# ---------------------------------------------------------------------------
# Functional category definitions
# ---------------------------------------------------------------------------
#
# Each category is a function (gene_symbol, subclass) -> bool. Returns True
# if the row qualifies as that category. Multiple categories may match a
# single row (e.g., a CTX-M is both ESBL and CEPHALOSPORIN); that's fine.
#
# When a regex anchor is at the start (^), it requires gene_symbol to start
# with the pattern, otherwise it allows it anywhere in the symbol. The
# subclass field is used as a fallback for broad categories.

def _matches_pattern(s: str, pattern: str) -> bool:
    if not isinstance(s, str):
        return False
    return bool(re.search(pattern, s, flags=re.IGNORECASE))


def is_esbl(sym: str, subclass: str) -> bool:
    if not isinstance(sym, str):
        return False
    # Any blaCTX-M is treated as ESBL.
    if _matches_pattern(sym, r'^blaCTX-?M'):
        return True
    # blaTEM ESBLs: variants other than the parental TEM-1, TEM-1A/B/C/D, TEM-2.
    if _matches_pattern(sym, r'^blaTEM-?\d'):
        m = re.search(r'^blaTEM-?(\d+)', sym, flags=re.IGNORECASE)
        if m and int(m.group(1)) in (1, 2):
            return False
        return True
    # blaSHV ESBLs: exclude SHV-1 (parental) and SHV-11 (also wild-type-equivalent).
    if _matches_pattern(sym, r'^blaSHV-?\d'):
        m = re.search(r'^blaSHV-?(\d+)', sym, flags=re.IGNORECASE)
        if m and int(m.group(1)) in (1, 11):
            return False
        return True
    # Other ESBL families
    if _matches_pattern(sym, r'^bla(GES|PER|VEB|BES|TLA|SFO|OXY)\b'):
        return True
    return False


def is_carbapenemase(sym: str, subclass: str) -> bool:
    if not isinstance(sym, str):
        return False
    # Class B (metallo-): NDM, VIM, IMP, GIM, SIM, SPM
    if _matches_pattern(sym, r'^bla(NDM|VIM|IMP|GIM|SIM|SPM)-?'):
        return True
    # Class A: KPC, GES (some), SME, IMI, NMC-A
    if _matches_pattern(sym, r'^bla(KPC|SME|IMI|NMC-A)-?'):
        return True
    # OXA-48-like: 48, 162, 181, 204, 232, 244, 370, 436, 438, 484, 505, 519, etc.
    # Cleanest filter: the OXA-48-family numbers cluster; treat any OXA in subclass=CARBAPENEM as carbapenemase
    if _matches_pattern(sym, r'^blaOXA-?(48|162|181|204|232|244|370|436|438|484|505|519)\b'):
        return True
    if isinstance(subclass, str) and 'CARBAPENEM' in subclass.upper() and _matches_pattern(sym, r'^blaOXA-?\d'):
        return True
    return False


def is_methicillin(sym: str, subclass: str) -> bool:
    return _matches_pattern(sym, r'^mec[AC](\b|2|-)')


def is_aminoglycoside(sym: str, subclass: str) -> bool:
    if isinstance(subclass, str):
        sc = subclass.upper()
        if any(k in sc for k in ['AMINOGLYCOSIDE', 'AMIKACIN', 'GENTAMICIN', 'KANAMYCIN',
                                  'STREPTOMYCIN', 'TOBRAMYCIN']):
            return True
    # 16S rRNA methyltransferases (pan-aminoglycoside)
    if _matches_pattern(sym, r'^(rmt[A-G]|armA|npmA)\b'):
        return True
    return False


def is_fluoroquinolone(sym: str, subclass: str) -> bool:
    if isinstance(sym, str):
        # Acquired plasmid-mediated FQ resistance markers
        if _matches_pattern(sym, r"^(qnr[A-Z]?\d*|qep[AB]?|oqx[AB])\b"):
            return True
        if _matches_pattern(sym, r"aac\(?6'?\)?-?Ib-cr"):
            return True
    return False


def is_macrolide(sym: str, subclass: str) -> bool:
    if isinstance(subclass, str):
        sc = subclass.upper()
        if any(k in sc for k in ['CLINDAMYCIN', 'ERYTHROMYCIN', 'AZITHROMYCIN',
                                 'STREPTOGRAMIN B', 'TYLOSIN', 'LINCOMYCIN']):
            return True
    if _matches_pattern(sym, r'^(erm|mef|msr|mph|lnu)[A-Z]?'):
        return True
    return False


def is_tetracycline(sym: str, subclass: str) -> bool:
    if isinstance(subclass, str):
        sc = subclass.upper()
        if 'TETRACYCLINE' in sc or 'TIGECYCLINE' in sc:
            return True
    if _matches_pattern(sym, r'^tet[A-Z]'):
        return True
    return False


def is_sulfonamide(sym: str, subclass: str) -> bool:
    if isinstance(subclass, str) and 'SULFONAMIDE' in subclass.upper():
        return True
    return _matches_pattern(sym, r'^sul\d')


def is_trimethoprim(sym: str, subclass: str) -> bool:
    if isinstance(subclass, str) and 'TRIMETHOPRIM' in subclass.upper():
        return True
    return _matches_pattern(sym, r'^dfr[A-Z]?')


def is_colistin(sym: str, subclass: str) -> bool:
    if _matches_pattern(sym, r'^mcr-?\d'):
        return True
    return False


def is_vre(sym: str, subclass: str) -> bool:
    return _matches_pattern(sym, r'^van[ABCDEFGM](\b|-)')


CATEGORIES: list[tuple[str, callable, str]] = [
    ('ESBL',                is_esbl,           'blaCTX-M, ESBL TEM/SHV variants, GES/PER/VEB families'),
    ('Carbapenemase',       is_carbapenemase,  'KPC/NDM/VIM/IMP/OXA-48-family'),
    ('Methicillin (mecA)',  is_methicillin,    'mecA / mecC (MRSA marker)'),
    ('Aminoglycoside',      is_aminoglycoside, 'subclass aminoglycoside or 16S rmt/armA/npmA'),
    ('Fluoroquinolone',     is_fluoroquinolone,"qnr, qep, oqxAB, aac(6')-Ib-cr"),
    ('Macrolide (MLS)',     is_macrolide,      'subclass erythromycin/clindamycin/azithromycin or erm/mef/msr/mph/lnu'),
    ('Tetracycline',        is_tetracycline,   'subclass tetracycline/tigecycline or tet*'),
    ('Sulfonamide',         is_sulfonamide,    'subclass sulfonamide or sul*'),
    ('Trimethoprim',        is_trimethoprim,   'subclass trimethoprim or dfr*'),
    ('Colistin (mcr)',      is_colistin,       'mcr-* family'),
    ('Vancomycin (VRE)',    is_vre,            'vanA/B/C/D/E/G/M'),
]


def write_definitions():
    rows = [{'category': name, 'definition': desc} for name, _, desc in CATEGORIES]
    pd.DataFrame(rows).to_csv(OUT_DIR / "functional_class_definitions.tsv",
                              sep="\t", index=False)


# ---------------------------------------------------------------------------
# Build per-(subject, body site) presence/absence
# ---------------------------------------------------------------------------

def build_presence(metadata: pd.DataFrame, amr: pd.DataFrame):
    """Returns:
    - sample_pres: one row per sample with boolean columns per category
    - subj_site: one row per (SubjectID, SampleType) — pooled across weeks
    - subj_site_week: one row per (SubjectID, SampleType, SampleCollectionWeek)
    """
    print("\nClassifying rows into clinical categories...")
    # Filter to RPM threshold
    qual = amr[amr['rpm'] >= PRESENCE_RPM].copy()
    print(f"  Rows with RPM >= {PRESENCE_RPM}: {len(qual):,} / {len(amr):,}")

    # Tag each row with categories it matches (for performance, vectorize)
    # Actually we'll loop categories and produce a {sample_name: bool} per category
    sample_md = metadata[['sample_name', 'SubjectID', 'Location', 'SampleType',
                          'SampleCollectionWeek', 'PostNatalAbxCohort']].drop_duplicates('sample_name')

    # Restrict to relevant body sites
    sample_md = sample_md[sample_md['SampleType'].isin(SITE_ORDER)].copy()

    presence_per_sample = {'sample_name': sample_md['sample_name'].values}
    for cat_name, fn, _ in CATEGORIES:
        sym = qual['gene_symbol']
        subc = qual.get('subclass', pd.Series([np.nan]*len(qual)))
        mask = pd.Series(
            [fn(s, c) for s, c in zip(sym, subc)], index=qual.index
        )
        cat_samples = set(qual.loc[mask, 'sample_name'].unique())
        presence_per_sample[cat_name] = sample_md['sample_name'].isin(cat_samples).values
        print(f"  {cat_name:24s}  positive samples: {len(cat_samples):4d} / {len(sample_md):4d}")

    sample_pres = pd.DataFrame(presence_per_sample)
    sample_pres = sample_pres.merge(sample_md, on='sample_name', how='left')

    # Aggregate to (SubjectID, SampleType) — positive at site if ANY sample positive
    cat_cols = [c for c, _, _ in CATEGORIES]
    subj_site = (sample_pres
                 .groupby(['SubjectID', 'SampleType', 'Location'])[cat_cols]
                 .any()
                 .reset_index())
    # Same but keeping Week as a dimension (per-week per-site presence)
    subj_site_week = (sample_pres
                      .groupby(['SubjectID', 'SampleType', 'SampleCollectionWeek',
                                'Location'])[cat_cols]
                      .any()
                      .reset_index())
    print(f"\n(subject, body site) rows: {len(subj_site)}")
    print(f"(subject, body site, week) rows: {len(subj_site_week)}")
    return sample_pres, subj_site, subj_site_week


def prevalence_table(subj_site: pd.DataFrame) -> pd.DataFrame:
    """% subjects positive for each category at each body site."""
    cat_cols = [c for c, _, _ in CATEGORIES]
    rows = []
    for site in SITE_ORDER:
        site_data = subj_site[subj_site['SampleType'] == site]
        n_subj = site_data['SubjectID'].nunique()
        for cat in cat_cols:
            pos_subj = site_data[site_data[cat]]['SubjectID'].nunique()
            rows.append({
                'category': cat,
                'SampleType': site,
                'n_subjects_screened': n_subj,
                'n_subjects_positive': pos_subj,
                'pct_positive': 100.0 * pos_subj / n_subj if n_subj > 0 else 0.0,
            })
    # Also "any site" row
    cat_cols_idx = subj_site.set_index(['SubjectID', 'SampleType'])[cat_cols]
    any_site = cat_cols_idx.groupby('SubjectID').any()
    n_total_subj = any_site.shape[0]
    for cat in cat_cols:
        pos = int(any_site[cat].sum())
        rows.append({
            'category': cat,
            'SampleType': 'Any site',
            'n_subjects_screened': n_total_subj,
            'n_subjects_positive': pos,
            'pct_positive': 100.0 * pos / n_total_subj if n_total_subj > 0 else 0.0,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_prevalence_by_site(prev: pd.DataFrame):
    """Grouped bar: % subjects positive by category × body site (+ Any site)."""
    df = prev.copy()
    cat_order = [c for c, _, _ in CATEGORIES]

    # Order categories by overall prevalence (Any site) descending
    overall = df[df['SampleType'] == 'Any site'].set_index('category')['pct_positive']
    cat_order = overall.sort_values(ascending=False).index.tolist()

    fig, ax = plt.subplots(figsize=(13, 7))
    width = 0.20
    site_seq = SITE_ORDER + ['Any site']
    palette = {**SITE_PALETTE, 'Any site': '#444444'}
    x = np.arange(len(cat_order))
    for i, site in enumerate(site_seq):
        sub = df[df['SampleType'] == site].set_index('category').reindex(cat_order)
        ax.bar(x + (i - 1.5) * width, sub['pct_positive'].values,
               width, label=site, color=palette[site], edgecolor='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(cat_order, rotation=30, ha='right')
    ax.set_ylabel('% of infants positive')
    ax.set_xlabel('Functional resistance category')
    ax.set_title(f'Functional-class screening — % infants positive by body site\n'
                 f'(presence threshold: RPM ≥ {PRESENCE_RPM} in strict-filtered AMR; '
                 f'any sample of subject at that site)')
    ax.legend(title='Site', loc='upper right')
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    out = FIGURES_DIR / "functional_class_prevalence_by_site.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {out}")


def plot_venn_per_category(subj_site: pd.DataFrame):
    """Venn diagram per category: which subjects are positive at each
    body site? Helps see whether stool/skin capture the same subjects."""
    cat_cols = [c for c, _, _ in CATEGORIES]
    n_cats = len(cat_cols)

    fig, axes = plt.subplots(3, 4, figsize=(18, 14))
    axes = axes.flatten()

    for i, cat in enumerate(cat_cols):
        ax = axes[i]
        sets = {}
        for site in SITE_ORDER:
            ss = subj_site[(subj_site['SampleType'] == site) & subj_site[cat]]
            sets[site] = set(ss['SubjectID'].unique())
        total_positive = len(sets['Axilla'] | sets['Groin'] | sets['Stool'])
        if total_positive == 0:
            ax.text(0.5, 0.5, f'{cat}\n(no positives)', ha='center', va='center',
                    transform=ax.transAxes)
            ax.set_axis_off()
            continue
        v = venn3([sets['Axilla'], sets['Groin'], sets['Stool']],
                  set_labels=('Axilla', 'Groin', 'Stool'),
                  set_colors=(SITE_PALETTE['Axilla'], SITE_PALETTE['Groin'], SITE_PALETTE['Stool']),
                  alpha=0.55, ax=ax)
        ax.set_title(f"{cat} (n={total_positive} subjects positive)", fontsize=11)

    # Hide unused axes
    for j in range(n_cats, len(axes)):
        axes[j].set_axis_off()

    fig.suptitle('Subject overlap across body sites — each category',
                 fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    out = FIGURES_DIR / "venn_per_category.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {out}")


def screening_capture_table(subj_site: pd.DataFrame) -> pd.DataFrame:
    """For each category: of all positive-anywhere subjects, what % are
    caught by stool alone, groin alone, axilla alone, and pairwise/all-three
    combinations? This is the 'where to swab' answer."""
    cat_cols = [c for c, _, _ in CATEGORIES]
    rows = []
    for cat in cat_cols:
        sets = {}
        for site in SITE_ORDER:
            ss = subj_site[(subj_site['SampleType'] == site) & subj_site[cat]]
            sets[site] = set(ss['SubjectID'].unique())
        any_pos = sets['Axilla'] | sets['Groin'] | sets['Stool']
        n_any = len(any_pos)
        if n_any == 0:
            continue
        rows.append({
            'category': cat,
            'n_positive_total': n_any,
            'pct_caught_by_axilla_alone': 100.0 * len(sets['Axilla']) / n_any,
            'pct_caught_by_groin_alone':  100.0 * len(sets['Groin'])  / n_any,
            'pct_caught_by_stool_alone':  100.0 * len(sets['Stool'])  / n_any,
            'pct_caught_by_axilla_or_groin': 100.0 * len(sets['Axilla'] | sets['Groin']) / n_any,
            'pct_caught_by_axilla_or_stool': 100.0 * len(sets['Axilla'] | sets['Stool']) / n_any,
            'pct_caught_by_groin_or_stool':  100.0 * len(sets['Groin']  | sets['Stool']) / n_any,
            'pct_caught_by_all_three':       100.0,  # by definition
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "screening_site_capture.tsv", sep="\t", index=False)
    print(f"✓ Saved: {OUT_DIR/'screening_site_capture.tsv'}")
    return out


def plot_capture(capture: pd.DataFrame):
    """Heatmap: category × site combination, color = % positives caught."""
    if capture.empty:
        return
    cat_cols = capture['category'].tolist()
    cols = [
        ('Axilla only',         'pct_caught_by_axilla_alone'),
        ('Groin only',          'pct_caught_by_groin_alone'),
        ('Stool only',          'pct_caught_by_stool_alone'),
        ('Axilla + Groin',      'pct_caught_by_axilla_or_groin'),
        ('Axilla + Stool',      'pct_caught_by_axilla_or_stool'),
        ('Groin + Stool',       'pct_caught_by_groin_or_stool'),
        ('All three sites',     'pct_caught_by_all_three'),
    ]
    mat = capture.set_index('category')[[c for _, c in cols]].copy()
    mat.columns = [n for n, _ in cols]

    fig, ax = plt.subplots(figsize=(11, 0.45 * len(cat_cols) + 2))
    sns.heatmap(mat, cmap='YlGnBu', vmin=0, vmax=100, ax=ax,
                annot=True, fmt='.0f', linewidths=0.4, linecolor='white',
                cbar_kws={'label': '% of positive subjects captured'})
    ax.set_xlabel('Screening strategy')
    ax.set_ylabel('Category')
    ax.set_title("Screening capture — % of positive infants caught by each site combination\n"
                 "(of subjects positive at any site)")
    plt.xticks(rotation=20, ha='right')
    plt.tight_layout()
    out = FIGURES_DIR / "screening_site_capture.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {out}")


def prevalence_by_site_week(subj_site_week: pd.DataFrame) -> pd.DataFrame:
    """% subjects positive for each category at each (body site, week)."""
    cat_cols = [c for c, _, _ in CATEGORIES]
    rows = []
    for site in SITE_ORDER:
        for week in ['Week.1', 'Week.3']:
            sub = subj_site_week[(subj_site_week['SampleType'] == site) &
                                 (subj_site_week['SampleCollectionWeek'] == week)]
            n_subj = sub['SubjectID'].nunique()
            for cat in cat_cols:
                pos = sub[sub[cat]]['SubjectID'].nunique()
                rows.append({
                    'category': cat,
                    'SampleType': site,
                    'SampleCollectionWeek': week,
                    'n_subjects_screened': n_subj,
                    'n_subjects_positive': pos,
                    'pct_positive': 100.0 * pos / n_subj if n_subj > 0 else 0.0,
                })
    return pd.DataFrame(rows)


def prevalence_by_site_location(subj_site: pd.DataFrame) -> pd.DataFrame:
    """% subjects positive for each category at each (body site, location)."""
    cat_cols = [c for c, _, _ in CATEGORIES]
    rows = []
    for site in SITE_ORDER:
        for loc in ['UCMC', 'ZCH']:
            sub = subj_site[(subj_site['SampleType'] == site) &
                            (subj_site['Location'] == loc)]
            n_subj = sub['SubjectID'].nunique()
            for cat in cat_cols:
                pos = sub[sub[cat]]['SubjectID'].nunique()
                rows.append({
                    'category': cat,
                    'SampleType': site,
                    'Location': loc,
                    'n_subjects_screened': n_subj,
                    'n_subjects_positive': pos,
                    'pct_positive': 100.0 * pos / n_subj if n_subj > 0 else 0.0,
                })
    # Also "any site" by location
    cat_cols_idx = subj_site.set_index(['SubjectID', 'Location', 'SampleType'])[cat_cols]
    any_site = cat_cols_idx.groupby(['SubjectID', 'Location']).any().reset_index()
    for loc in ['UCMC', 'ZCH']:
        sub = any_site[any_site['Location'] == loc]
        n = sub['SubjectID'].nunique()
        for cat in cat_cols:
            pos = int(sub[cat].sum())
            rows.append({
                'category': cat,
                'SampleType': 'Any site',
                'Location': loc,
                'n_subjects_screened': n,
                'n_subjects_positive': pos,
                'pct_positive': 100.0 * pos / n if n > 0 else 0.0,
            })
    return pd.DataFrame(rows)


def plot_prevalence_by_week(prev_week: pd.DataFrame):
    """3-panel (one per body site): grouped bars showing Week 1 vs Week 3."""
    cat_cols = [c for c, _, _ in CATEGORIES]
    overall = (prev_week[prev_week['SampleType'] == 'Stool']
               .groupby('category')['pct_positive'].max())
    cat_order = overall.sort_values(ascending=False).index.tolist()

    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    for j, site in enumerate(SITE_ORDER):
        ax = axes[j]
        sub = prev_week[prev_week['SampleType'] == site]
        x = np.arange(len(cat_order))
        width = 0.4
        for i, week in enumerate(['Week.1', 'Week.3']):
            d = sub[sub['SampleCollectionWeek'] == week].set_index('category').reindex(cat_order)
            ax.bar(x + (i - 0.5) * width, d['pct_positive'].values,
                   width, label=week,
                   color=('#9ecae1' if week == 'Week.1' else '#3182bd'),
                   edgecolor='black', linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(cat_order, rotation=30, ha='right')
        ax.set_title(f'{site}', fontsize=12, fontweight='bold')
        if j == 0:
            ax.set_ylabel('% of infants positive')
            ax.legend(title='Week', loc='upper right')
        ax.grid(axis='y', alpha=0.3)
    fig.suptitle('Functional-class prevalence: Week 1 vs Week 3 (per body site)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = FIGURES_DIR / "functional_class_prevalence_by_week.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {out}")


def prevalence_by_site_week_location(subj_site_week: pd.DataFrame) -> pd.DataFrame:
    """% subjects positive at each (body site, week, location) combination —
    the full 3-way breakdown for the 6-panel view."""
    cat_cols = [c for c, _, _ in CATEGORIES]
    rows = []
    for site in SITE_ORDER:
        for week in ['Week.1', 'Week.3']:
            for loc in ['UCMC', 'ZCH']:
                sub = subj_site_week[
                    (subj_site_week['SampleType'] == site) &
                    (subj_site_week['SampleCollectionWeek'] == week) &
                    (subj_site_week['Location'] == loc)
                ]
                n = sub['SubjectID'].nunique()
                for cat in cat_cols:
                    pos = sub[sub[cat]]['SubjectID'].nunique()
                    rows.append({
                        'category': cat,
                        'SampleType': site,
                        'SampleCollectionWeek': week,
                        'Location': loc,
                        'n_subjects_screened': n,
                        'n_subjects_positive': pos,
                        'pct_positive': 100.0 * pos / n if n > 0 else 0.0,
                    })
    return pd.DataFrame(rows)


def plot_prevalence_6panel(prev_3way: pd.DataFrame):
    """6-panel grid: rows = Week 1 / Week 3, cols = Axilla / Groin / Stool;
    within each panel, UCMC vs ZCH bars per category."""
    cat_cols = [c for c, _, _ in CATEGORIES]
    # Order categories by overall (any-site) max prevalence, descending
    overall = prev_3way.groupby('category')['pct_positive'].max()
    cat_order = overall.sort_values(ascending=False).index.tolist()

    fig, axes = plt.subplots(2, 3, figsize=(20, 10), sharey=True)
    weeks = ['Week.1', 'Week.3']
    width = 0.4

    for r, week in enumerate(weeks):
        for c, site in enumerate(SITE_ORDER):
            ax = axes[r, c]
            sub = prev_3way[(prev_3way['SampleCollectionWeek'] == week) &
                            (prev_3way['SampleType'] == site)]
            x = np.arange(len(cat_order))
            for i, loc in enumerate(['UCMC', 'ZCH']):
                d = sub[sub['Location'] == loc].set_index('category').reindex(cat_order)
                ax.bar(x + (i - 0.5) * width, d['pct_positive'].values,
                       width, label=loc,
                       color=('#1f77b4' if loc == 'UCMC' else '#ff7f0e'),
                       edgecolor='black', linewidth=0.5)
            ax.set_xticks(x)
            if r == 1:
                ax.set_xticklabels(cat_order, rotation=30, ha='right')
            else:
                ax.set_xticklabels([])
            n_ucmc = int(sub[sub['Location'] == 'UCMC']['n_subjects_screened'].max() or 0)
            n_zch  = int(sub[sub['Location'] == 'ZCH']['n_subjects_screened'].max() or 0)
            ax.set_title(f'{site} — {week}  (UCMC n={n_ucmc}, ZCH n={n_zch})',
                         fontsize=11, fontweight='bold')
            if c == 0:
                ax.set_ylabel('% of infants positive')
            if r == 0 and c == 0:
                ax.legend(title='Location', loc='upper right')
            ax.grid(axis='y', alpha=0.3)
            ax.set_ylim(0, 100)

    fig.suptitle('Functional-class prevalence — Body site × Week × Location',
                 fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    out = FIGURES_DIR / "functional_class_prevalence_6panel.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {out}")


def plot_prevalence_by_location(prev_loc: pd.DataFrame):
    """4-panel (one per site + 'Any site'): UCMC vs ZCH bars."""
    cat_cols = [c for c, _, _ in CATEGORIES]
    overall = (prev_loc[prev_loc['SampleType'] == 'Any site']
               .groupby('category')['pct_positive'].max())
    cat_order = overall.sort_values(ascending=False).index.tolist()

    site_seq = SITE_ORDER + ['Any site']
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), sharey=True)
    axes = axes.flatten()
    for j, site in enumerate(site_seq):
        ax = axes[j]
        sub = prev_loc[prev_loc['SampleType'] == site]
        x = np.arange(len(cat_order))
        width = 0.4
        for i, loc in enumerate(['UCMC', 'ZCH']):
            d = sub[sub['Location'] == loc].set_index('category').reindex(cat_order)
            ax.bar(x + (i - 0.5) * width, d['pct_positive'].values,
                   width, label=loc,
                   color=('#1f77b4' if loc == 'UCMC' else '#ff7f0e'),
                   edgecolor='black', linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(cat_order, rotation=30, ha='right')
        ax.set_title(f'{site}', fontsize=12, fontweight='bold')
        if j % 2 == 0:
            ax.set_ylabel('% of infants positive')
        if j == 0:
            ax.legend(title='Location', loc='upper right')
        ax.grid(axis='y', alpha=0.3)
    fig.suptitle('Functional-class prevalence by Location (UCMC vs ZCH)',
                 fontsize=13, fontweight='bold', y=1.00)
    plt.tight_layout()
    out = FIGURES_DIR / "functional_class_prevalence_by_location.pdf"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: {out}")


def main():
    print("=" * 60)
    print("FUNCTIONAL CLASS SCREENING")
    print("=" * 60)
    write_definitions()

    metadata = pd.read_csv(QC_DIR / "master_metadata_with_qc.tsv", sep="\t")
    metadata = metadata[metadata['has_amr_data']].drop_duplicates('sample_name').copy()

    cols_needed = ['sample_name', 'gene_symbol', 'gene_family', 'subclass', 'rpm']
    amr = pd.read_csv(DATA_DIR / "nicu_amr_only.tsv", sep="\t", usecols=cols_needed)

    sample_pres, subj_site, subj_site_week = build_presence(metadata, amr)
    sample_pres.to_csv(OUT_DIR / "functional_class_per_sample.tsv", sep="\t", index=False)
    subj_site.to_csv(OUT_DIR / "functional_class_per_subject_site.tsv", sep="\t", index=False)
    subj_site_week.to_csv(OUT_DIR / "functional_class_per_subject_site_week.tsv", sep="\t", index=False)
    print(f"✓ Saved per-sample, per-subject-site, per-subject-site-week tables")

    # Headline (pooled across weeks): % subjects positive by site
    prev = prevalence_table(subj_site)
    prev.to_csv(OUT_DIR / "functional_class_prevalence_by_site.tsv", sep="\t", index=False)

    # Stratifications requested
    prev_week = prevalence_by_site_week(subj_site_week)
    prev_week.to_csv(OUT_DIR / "functional_class_prevalence_by_site_week.tsv", sep="\t", index=False)
    print(f"✓ Saved: {OUT_DIR/'functional_class_prevalence_by_site_week.tsv'}")

    prev_loc = prevalence_by_site_location(subj_site)
    prev_loc.to_csv(OUT_DIR / "functional_class_prevalence_by_site_location.tsv", sep="\t", index=False)
    print(f"✓ Saved: {OUT_DIR/'functional_class_prevalence_by_site_location.tsv'}")

    prev_3way = prevalence_by_site_week_location(subj_site_week)
    prev_3way.to_csv(OUT_DIR / "functional_class_prevalence_by_site_week_location.tsv",
                     sep="\t", index=False)
    print(f"✓ Saved: {OUT_DIR/'functional_class_prevalence_by_site_week_location.tsv'}")

    plot_prevalence_by_site(prev)
    plot_venn_per_category(subj_site)
    capture = screening_capture_table(subj_site)
    plot_capture(capture)
    plot_prevalence_by_week(prev_week)
    plot_prevalence_by_location(prev_loc)
    plot_prevalence_6panel(prev_3way)

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
