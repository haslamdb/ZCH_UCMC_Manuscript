#!/usr/bin/env python3
"""
ESBL Gene Analysis for NICU Resistome Study

Analyzes Extended-Spectrum Beta-Lactamase (ESBL) genes including:
- CTX-M type genes (most clinically relevant ESBLs)
- SHV variants (both ESBL and non-ESBL)
- TEM variants (both ESBL and non-ESBL)
- Carbapenemases (KPC, NDM, VIM, IMP, OXA-48-like) for comparison

Comparisons:
1. UCMC vs ZCH (location)
2. Body site differences (Axilla, Groin, Stool)
3. Association with antibiotic exposure (especially beta-lactams)
4. Temporal changes (Week 1 vs Week 3)

Output: Publication-ready figures and statistical tables
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import mannwhitneyu, kruskal, spearmanr, fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Define paths - update to use the diamond_traditional directory
DIAMOND_DIR = Path("/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis/diamond_traditional")
DATA_DIR = DIAMOND_DIR / "data"
RESULTS_DIR = DIAMOND_DIR / "results" / "esbl_analysis_tables"
FIGURES_DIR = DIAMOND_DIR / "figures" / "esbl_analysis"

# Create directories
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Define ESBL gene patterns
ESBL_PATTERNS = {
    'CTX-M': 'blaCTX-M',  # All CTX-M variants are ESBLs
    'SHV': 'blaSHV',       # Some SHV variants are ESBLs
    'TEM': 'blaTEM',       # Some TEM variants are ESBLs
}

# Carbapenemases for comparison
CARBAPENEMASE_PATTERNS = {
    'KPC': 'blaKPC',
    'NDM': 'blaNDM',
    'VIM': 'blaVIM',
    'IMP': 'blaIMP',
    'OXA-48': 'blaOXA-48',  # OXA-48-like carbapenemases
}

# Known ESBL variants (common clinically relevant ones)
KNOWN_ESBL_TEM = ['blaTEM-3', 'blaTEM-4', 'blaTEM-5', 'blaTEM-6', 'blaTEM-10',
                  'blaTEM-12', 'blaTEM-26', 'blaTEM-52', 'blaTEM-116']
KNOWN_ESBL_SHV = ['blaSHV-2', 'blaSHV-3', 'blaSHV-4', 'blaSHV-5', 'blaSHV-12',
                  'blaSHV-18', 'blaSHV-24', 'blaSHV-28', 'blaSHV-38']

def load_data():
    """Load DIAMOND combined data and metadata"""
    print("Loading data...")

    # Load combined DIAMOND results
    combined_data = pd.read_csv(DATA_DIR / "diamond_amr_combined.tsv", sep="\t")
    print(f"  Combined data: {len(combined_data):,} entries")

    # Load metadata
    metadata = pd.read_csv(DATA_DIR / "diamond_metadata_updated.tsv", sep="\t")
    print(f"  Metadata: {len(metadata)} samples")

    # Load master metadata with antibiotic exposure
    master_metadata = pd.read_csv(DATA_DIR / "master_metadata_fixed.tsv", sep="\t")
    print(f"  Master metadata: {len(master_metadata)} samples")

    return combined_data, metadata, master_metadata

def identify_esbl_genes(combined_data):
    """Identify all ESBL-related genes in the dataset"""
    print("\nIdentifying ESBL and carbapenemase genes...")

    all_genes = combined_data['gene_name'].unique()
    print(f"  Total unique genes in dataset: {len(all_genes)}")

    esbl_genes = {}

    # CTX-M genes (all are ESBLs)
    ctxm_genes = [g for g in all_genes if g.startswith('blaCTX-M')]
    esbl_genes['CTX-M'] = ctxm_genes
    print(f"  CTX-M genes found: {len(ctxm_genes)}")

    # SHV genes
    shv_genes = [g for g in all_genes if g.startswith('blaSHV')]
    esbl_genes['SHV'] = shv_genes
    print(f"  SHV genes found: {len(shv_genes)}")

    # TEM genes
    tem_genes = [g for g in all_genes if g.startswith('blaTEM')]
    esbl_genes['TEM'] = tem_genes
    print(f"  TEM genes found: {len(tem_genes)}")

    # Carbapenemases
    carbapenemase_genes = {}
    for name, pattern in CARBAPENEMASE_PATTERNS.items():
        genes = [g for g in all_genes if g.startswith(pattern)]
        if genes:
            carbapenemase_genes[name] = genes
            print(f"  {name} genes found: {len(genes)}")

    return esbl_genes, carbapenemase_genes

def create_esbl_abundance_matrix(combined_data, esbl_genes, carbapenemase_genes, metadata):
    """Create sample × ESBL gene abundance matrix"""
    print("\nCreating ESBL abundance matrix...")

    # Combine all ESBL and carbapenemase genes
    all_target_genes = []
    for gene_list in esbl_genes.values():
        all_target_genes.extend(gene_list)
    for gene_list in carbapenemase_genes.values():
        all_target_genes.extend(gene_list)

    # Filter combined data to target genes
    esbl_data = combined_data[combined_data['gene_name'].isin(all_target_genes)].copy()
    print(f"  Entries for target genes: {len(esbl_data):,}")

    # Create pivot table (samples × genes with RPM values)
    rpm_matrix = esbl_data.pivot_table(
        index='sample_name',
        columns='gene_name',
        values='rpm',
        aggfunc='sum',
        fill_value=0
    )

    print(f"  ESBL RPM matrix: {rpm_matrix.shape[0]} samples × {rpm_matrix.shape[1]} genes")

    # Merge with metadata
    rpm_matrix = rpm_matrix.reset_index()
    rpm_matrix = rpm_matrix.merge(metadata, on='sample_name', how='inner')

    return rpm_matrix, all_target_genes

def calculate_gene_family_totals(rpm_matrix, esbl_genes, carbapenemase_genes):
    """Calculate total RPM for each gene family"""
    print("\nCalculating gene family totals...")

    # CTX-M total
    ctxm_cols = [c for c in rpm_matrix.columns if c in esbl_genes.get('CTX-M', [])]
    if ctxm_cols:
        rpm_matrix['CTX-M_total'] = rpm_matrix[ctxm_cols].sum(axis=1)
        print(f"  CTX-M: {len(ctxm_cols)} variants")
    else:
        rpm_matrix['CTX-M_total'] = 0

    # SHV total
    shv_cols = [c for c in rpm_matrix.columns if c in esbl_genes.get('SHV', [])]
    if shv_cols:
        rpm_matrix['SHV_total'] = rpm_matrix[shv_cols].sum(axis=1)
        print(f"  SHV: {len(shv_cols)} variants")
    else:
        rpm_matrix['SHV_total'] = 0

    # TEM total
    tem_cols = [c for c in rpm_matrix.columns if c in esbl_genes.get('TEM', [])]
    if tem_cols:
        rpm_matrix['TEM_total'] = rpm_matrix[tem_cols].sum(axis=1)
        print(f"  TEM: {len(tem_cols)} variants")
    else:
        rpm_matrix['TEM_total'] = 0

    # Total ESBLs (CTX-M + SHV + TEM)
    rpm_matrix['ESBL_total'] = (rpm_matrix['CTX-M_total'] +
                                rpm_matrix['SHV_total'] +
                                rpm_matrix['TEM_total'])

    # Carbapenemases
    for name, genes in carbapenemase_genes.items():
        carb_cols = [c for c in rpm_matrix.columns if c in genes]
        if carb_cols:
            rpm_matrix[f'{name}_total'] = rpm_matrix[carb_cols].sum(axis=1)
            print(f"  {name}: {len(carb_cols)} variants")
        else:
            rpm_matrix[f'{name}_total'] = 0

    # Total carbapenemases
    carb_total_cols = [f'{name}_total' for name in carbapenemase_genes.keys()]
    carb_total_cols = [c for c in carb_total_cols if c in rpm_matrix.columns]
    rpm_matrix['Carbapenemase_total'] = rpm_matrix[carb_total_cols].sum(axis=1)

    return rpm_matrix

def compare_locations(rpm_matrix):
    """Compare ESBL abundance between UCMC and ZCH"""
    print("\n" + "="*60)
    print("ESBL ABUNDANCE: UCMC vs ZCH")
    print("="*60)

    results = []
    gene_families = ['CTX-M_total', 'SHV_total', 'TEM_total', 'ESBL_total',
                     'KPC_total', 'NDM_total', 'VIM_total', 'Carbapenemase_total']

    ucmc = rpm_matrix[rpm_matrix['Location'] == 'UCMC']
    zch = rpm_matrix[rpm_matrix['Location'] == 'ZCH']

    print(f"\nSamples: UCMC = {len(ucmc)}, ZCH = {len(zch)}")

    for gene_family in gene_families:
        if gene_family not in rpm_matrix.columns:
            continue

        ucmc_vals = ucmc[gene_family]
        zch_vals = zch[gene_family]

        # Statistics
        stat, pval = mannwhitneyu(ucmc_vals, zch_vals, alternative='two-sided')

        # Effect size (log2 fold change of medians)
        ucmc_median = ucmc_vals.median()
        zch_median = zch_vals.median()

        # Avoid log(0)
        if ucmc_median > 0 and zch_median > 0:
            log2fc = np.log2(ucmc_median / zch_median)
        elif ucmc_median > 0:
            log2fc = np.inf
        elif zch_median > 0:
            log2fc = -np.inf
        else:
            log2fc = 0

        # Prevalence (proportion of samples with gene > 0)
        ucmc_prev = (ucmc_vals > 0).mean() * 100
        zch_prev = (zch_vals > 0).mean() * 100

        results.append({
            'Gene_Family': gene_family.replace('_total', ''),
            'UCMC_median_RPM': ucmc_median,
            'ZCH_median_RPM': zch_median,
            'UCMC_mean_RPM': ucmc_vals.mean(),
            'ZCH_mean_RPM': zch_vals.mean(),
            'UCMC_prevalence_%': ucmc_prev,
            'ZCH_prevalence_%': zch_prev,
            'log2FC_UCMC_vs_ZCH': log2fc,
            'Mann_Whitney_U': stat,
            'p_value': pval
        })

        print(f"\n{gene_family.replace('_total', '')}:")
        print(f"  UCMC: median={ucmc_median:.3f}, mean={ucmc_vals.mean():.3f}, prev={ucmc_prev:.1f}%")
        print(f"  ZCH:  median={zch_median:.3f}, mean={zch_vals.mean():.3f}, prev={zch_prev:.1f}%")
        print(f"  log2FC (UCMC/ZCH): {log2fc:.2f}, p={pval:.2e}")

    results_df = pd.DataFrame(results)

    # FDR correction
    _, pvals_corrected, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
    results_df['p_value_FDR'] = pvals_corrected

    return results_df

def compare_body_sites(rpm_matrix):
    """Compare ESBL abundance across body sites"""
    print("\n" + "="*60)
    print("ESBL ABUNDANCE BY BODY SITE")
    print("="*60)

    results = []
    gene_families = ['CTX-M_total', 'SHV_total', 'TEM_total', 'ESBL_total']

    # Map SampleType to body site names
    site_map = {'Axilla': 'Axilla', 'Groin': 'Groin', 'Stool': 'Stool'}
    rpm_matrix['BodySite'] = rpm_matrix['SampleType'].map(site_map)

    sites = ['Axilla', 'Groin', 'Stool']
    site_data = {site: rpm_matrix[rpm_matrix['BodySite'] == site] for site in sites}

    print(f"\nSamples by site:")
    for site, data in site_data.items():
        print(f"  {site}: n={len(data)}")

    for gene_family in gene_families:
        if gene_family not in rpm_matrix.columns:
            continue

        # Kruskal-Wallis test
        groups = [site_data[site][gene_family] for site in sites]
        h_stat, kw_pval = kruskal(*groups)

        # Pairwise Mann-Whitney U tests
        pairwise_results = {}
        for i, site1 in enumerate(sites):
            for j, site2 in enumerate(sites):
                if i < j:
                    stat, pval = mannwhitneyu(
                        site_data[site1][gene_family],
                        site_data[site2][gene_family],
                        alternative='two-sided'
                    )
                    pairwise_results[f'{site1}_vs_{site2}'] = pval

        result = {
            'Gene_Family': gene_family.replace('_total', ''),
            'Kruskal_Wallis_H': h_stat,
            'Kruskal_Wallis_p': kw_pval
        }

        for site in sites:
            result[f'{site}_median_RPM'] = site_data[site][gene_family].median()
            result[f'{site}_prevalence_%'] = (site_data[site][gene_family] > 0).mean() * 100

        for pair, pval in pairwise_results.items():
            result[f'p_{pair}'] = pval

        results.append(result)

        print(f"\n{gene_family.replace('_total', '')}:")
        print(f"  Kruskal-Wallis H={h_stat:.2f}, p={kw_pval:.2e}")
        for site in sites:
            median = site_data[site][gene_family].median()
            prev = (site_data[site][gene_family] > 0).mean() * 100
            print(f"  {site}: median={median:.3f}, prev={prev:.1f}%")

    return pd.DataFrame(results)

def analyze_antibiotic_associations(rpm_matrix, master_metadata):
    """Correlate ESBL genes with antibiotic exposure"""
    print("\n" + "="*60)
    print("ESBL ASSOCIATION WITH ANTIBIOTIC EXPOSURE")
    print("="*60)

    # Fix ZCH sample names in master metadata
    master_meta = master_metadata.copy()
    zch_mask = master_meta['Location'] == 'ZCH'
    master_meta.loc[zch_mask, 'sample_name'] = 'ZJH_' + master_meta.loc[zch_mask, 'sample_name']

    # Merge with ESBL data
    merged = rpm_matrix.merge(master_meta[['sample_name', 'SampleCollectionWeek'] +
                                          [c for c in master_meta.columns if c.endswith('_w1') or c.endswith('_w2')]],
                               on='sample_name', how='inner', suffixes=('', '_meta'))

    print(f"\nSamples with antibiotic data: {len(merged)}")

    # Calculate cumulative antibiotic exposure
    abx_cols_w1 = [c for c in merged.columns if c.endswith('_w1') and c != 'None_w1']
    abx_cols_w2 = [c for c in merged.columns if c.endswith('_w2') and c != 'None_w2']

    # Key beta-lactam antibiotics
    beta_lactams = ['Ampicillin', 'Nafcillin', 'Penicillin', 'Meropenem',
                    'Cefotaxime', 'Moxalactam', 'Cefazolin']

    results = []
    gene_families = ['CTX-M_total', 'SHV_total', 'TEM_total', 'ESBL_total']

    for abx in beta_lactams:
        w1_col = f'{abx}_w1'
        w2_col = f'{abx}_w2'

        if w1_col not in merged.columns:
            continue

        # Calculate cumulative exposure
        # Use the SampleCollectionWeek column with suffix if it exists, otherwise original
        week_col = 'SampleCollectionWeek_meta' if 'SampleCollectionWeek_meta' in merged.columns else 'SampleCollectionWeek'

        merged[f'{abx}_cumulative'] = merged.apply(
            lambda row: row.get(w1_col, 0) if pd.isna(row.get(w1_col)) == False else 0,
            axis=1
        )

        # Add w2 for Week 3 samples
        if w2_col in merged.columns:
            week3_mask = merged[week_col] == 'Week.3'
            merged.loc[week3_mask, f'{abx}_cumulative'] += merged.loc[week3_mask, w2_col].fillna(0)

        abx_exposure = merged[f'{abx}_cumulative'].fillna(0)
        n_exposed = (abx_exposure > 0).sum()

        if n_exposed < 5:
            continue

        for gene_family in gene_families:
            if gene_family not in merged.columns:
                continue

            gene_vals = merged[gene_family]

            # Spearman correlation
            rho, pval = spearmanr(abx_exposure, gene_vals)

            # Compare exposed vs unexposed
            exposed_vals = gene_vals[abx_exposure > 0]
            unexposed_vals = gene_vals[abx_exposure == 0]

            if len(exposed_vals) >= 5 and len(unexposed_vals) >= 5:
                mw_stat, mw_pval = mannwhitneyu(exposed_vals, unexposed_vals, alternative='two-sided')
            else:
                mw_stat, mw_pval = np.nan, np.nan

            results.append({
                'Antibiotic': abx,
                'Gene_Family': gene_family.replace('_total', ''),
                'n_exposed': n_exposed,
                'n_unexposed': len(unexposed_vals),
                'Spearman_rho': rho,
                'Spearman_p': pval,
                'Exposed_median_RPM': exposed_vals.median() if len(exposed_vals) > 0 else 0,
                'Unexposed_median_RPM': unexposed_vals.median() if len(unexposed_vals) > 0 else 0,
                'Mann_Whitney_p': mw_pval
            })

            if gene_family == 'CTX-M_total' or gene_family == 'ESBL_total':
                print(f"\n{abx} vs {gene_family.replace('_total', '')}:")
                print(f"  n_exposed={n_exposed}, rho={rho:.3f}, p={pval:.3e}")
                print(f"  Exposed median: {exposed_vals.median():.3f}")
                print(f"  Unexposed median: {unexposed_vals.median():.3f}")

    return pd.DataFrame(results)

def compare_weeks(rpm_matrix):
    """Compare ESBL abundance between Week 1 and Week 3"""
    print("\n" + "="*60)
    print("ESBL ABUNDANCE: Week 1 vs Week 3")
    print("="*60)

    results = []
    gene_families = ['CTX-M_total', 'SHV_total', 'TEM_total', 'ESBL_total']

    week1 = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.1']
    week3 = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.3']

    print(f"\nSamples: Week 1 = {len(week1)}, Week 3 = {len(week3)}")

    for gene_family in gene_families:
        if gene_family not in rpm_matrix.columns:
            continue

        w1_vals = week1[gene_family]
        w3_vals = week3[gene_family]

        stat, pval = mannwhitneyu(w1_vals, w3_vals, alternative='two-sided')

        w1_median = w1_vals.median()
        w3_median = w3_vals.median()

        if w1_median > 0 and w3_median > 0:
            log2fc = np.log2(w3_median / w1_median)
        elif w3_median > 0:
            log2fc = np.inf
        elif w1_median > 0:
            log2fc = -np.inf
        else:
            log2fc = 0

        results.append({
            'Gene_Family': gene_family.replace('_total', ''),
            'Week1_median_RPM': w1_median,
            'Week3_median_RPM': w3_median,
            'Week1_prevalence_%': (w1_vals > 0).mean() * 100,
            'Week3_prevalence_%': (w3_vals > 0).mean() * 100,
            'log2FC_W3_vs_W1': log2fc,
            'Mann_Whitney_U': stat,
            'p_value': pval
        })

        print(f"\n{gene_family.replace('_total', '')}:")
        print(f"  Week 1: median={w1_median:.3f}, prev={(w1_vals > 0).mean()*100:.1f}%")
        print(f"  Week 3: median={w3_median:.3f}, prev={(w3_vals > 0).mean()*100:.1f}%")
        print(f"  log2FC (W3/W1): {log2fc:.2f}, p={pval:.2e}")

    return pd.DataFrame(results)

def plot_esbl_comparison(rpm_matrix, output_prefix):
    """Create publication-ready figures for ESBL analysis"""
    print("\n" + "="*60)
    print("GENERATING FIGURES")
    print("="*60)

    # Set style
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")

    # Figure 1: ESBL abundance by location
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    gene_families = ['CTX-M_total', 'SHV_total', 'TEM_total', 'ESBL_total']
    titles = ['CTX-M', 'SHV', 'TEM', 'Total ESBL']

    for ax, gene_family, title in zip(axes.flat, gene_families, titles):
        if gene_family not in rpm_matrix.columns:
            continue

        # Add small value to 0s for log scale visualization
        plot_data = rpm_matrix.copy()
        plot_data[gene_family] = plot_data[gene_family].replace(0, 0.001)

        sns.boxplot(data=plot_data, x='Location', y=gene_family, ax=ax,
                    order=['UCMC', 'ZCH'], palette=['#2ecc71', '#e74c3c'])

        # Add individual points
        sns.stripplot(data=plot_data, x='Location', y=gene_family, ax=ax,
                      order=['UCMC', 'ZCH'], color='black', alpha=0.3, size=3)

        ax.set_ylabel('RPM (Reads Per Million)')
        ax.set_xlabel('')
        ax.set_title(f'{title} Genes', fontsize=12, fontweight='bold')
        ax.set_yscale('log')

        # Add p-value
        ucmc_vals = rpm_matrix[rpm_matrix['Location'] == 'UCMC'][gene_family]
        zch_vals = rpm_matrix[rpm_matrix['Location'] == 'ZCH'][gene_family]
        _, pval = mannwhitneyu(ucmc_vals, zch_vals, alternative='two-sided')

        if pval < 0.001:
            pval_str = f'p < 0.001'
        elif pval < 0.05:
            pval_str = f'p = {pval:.3f}'
        else:
            pval_str = f'p = {pval:.2f}'

        ax.text(0.5, 0.95, pval_str, transform=ax.transAxes,
                ha='center', va='top', fontsize=10)

    plt.suptitle('ESBL Gene Abundance: UCMC vs ZCH', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_prefix / 'esbl_by_location.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'esbl_by_location.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: esbl_by_location.pdf/png")

    # Figure 2: ESBL by body site
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for ax, gene_family, title in zip(axes.flat, gene_families, titles):
        if gene_family not in rpm_matrix.columns:
            continue

        plot_data = rpm_matrix.copy()
        plot_data[gene_family] = plot_data[gene_family].replace(0, 0.001)

        sns.boxplot(data=plot_data, x='SampleType', y=gene_family, ax=ax,
                    order=['Axilla', 'Groin', 'Stool'], palette='Set2')

        sns.stripplot(data=plot_data, x='SampleType', y=gene_family, ax=ax,
                      order=['Axilla', 'Groin', 'Stool'], color='black', alpha=0.3, size=3)

        ax.set_ylabel('RPM (Reads Per Million)')
        ax.set_xlabel('')
        ax.set_title(f'{title} Genes', fontsize=12, fontweight='bold')
        ax.set_yscale('log')

        # Kruskal-Wallis p-value
        groups = [plot_data[plot_data['SampleType'] == site][gene_family] for site in ['Axilla', 'Groin', 'Stool']]
        _, pval = kruskal(*groups)

        if pval < 0.001:
            pval_str = f'Kruskal-Wallis p < 0.001'
        elif pval < 0.05:
            pval_str = f'Kruskal-Wallis p = {pval:.3f}'
        else:
            pval_str = f'Kruskal-Wallis p = {pval:.2f}'

        ax.text(0.5, 0.95, pval_str, transform=ax.transAxes,
                ha='center', va='top', fontsize=10)

    plt.suptitle('ESBL Gene Abundance by Body Site', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_prefix / 'esbl_by_bodysite.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'esbl_by_bodysite.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: esbl_by_bodysite.pdf/png")

    # Figure 3: ESBL by location AND body site (faceted)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for ax, site in zip(axes, ['Axilla', 'Groin', 'Stool']):
        site_data = rpm_matrix[rpm_matrix['SampleType'] == site].copy()
        site_data['ESBL_total'] = site_data['ESBL_total'].replace(0, 0.001)

        sns.boxplot(data=site_data, x='Location', y='ESBL_total', ax=ax,
                    order=['UCMC', 'ZCH'], palette=['#2ecc71', '#e74c3c'])
        sns.stripplot(data=site_data, x='Location', y='ESBL_total', ax=ax,
                      order=['UCMC', 'ZCH'], color='black', alpha=0.3, size=3)

        ax.set_ylabel('ESBL RPM' if ax == axes[0] else '')
        ax.set_xlabel('')
        ax.set_title(f'{site}', fontsize=12, fontweight='bold')
        ax.set_yscale('log')

        # P-value
        ucmc_vals = site_data[site_data['Location'] == 'UCMC']['ESBL_total']
        zch_vals = site_data[site_data['Location'] == 'ZCH']['ESBL_total']

        if len(ucmc_vals) > 0 and len(zch_vals) > 0:
            _, pval = mannwhitneyu(ucmc_vals, zch_vals, alternative='two-sided')
            if pval < 0.001:
                pval_str = 'p < 0.001'
            elif pval < 0.05:
                pval_str = f'p = {pval:.3f}'
            else:
                pval_str = f'p = {pval:.2f}'
            ax.text(0.5, 0.95, pval_str, transform=ax.transAxes, ha='center', va='top', fontsize=10)

    plt.suptitle('Total ESBL Abundance: UCMC vs ZCH by Body Site', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_prefix / 'esbl_location_by_site.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'esbl_location_by_site.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: esbl_location_by_site.pdf/png")

    # Figure 4: Week comparison
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for ax, gene_family, title in zip(axes.flat, gene_families, titles):
        if gene_family not in rpm_matrix.columns:
            continue

        plot_data = rpm_matrix.copy()
        plot_data[gene_family] = plot_data[gene_family].replace(0, 0.001)

        sns.boxplot(data=plot_data, x='SampleCollectionWeek', y=gene_family, ax=ax,
                    order=['Week.1', 'Week.3'], palette=['#3498db', '#9b59b6'])

        sns.stripplot(data=plot_data, x='SampleCollectionWeek', y=gene_family, ax=ax,
                      order=['Week.1', 'Week.3'], color='black', alpha=0.3, size=3)

        ax.set_ylabel('RPM (Reads Per Million)')
        ax.set_xlabel('')
        ax.set_xticklabels(['Week 1', 'Week 3'])
        ax.set_title(f'{title} Genes', fontsize=12, fontweight='bold')
        ax.set_yscale('log')

        # P-value
        w1_vals = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.1'][gene_family]
        w3_vals = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.3'][gene_family]
        _, pval = mannwhitneyu(w1_vals, w3_vals, alternative='two-sided')

        if pval < 0.001:
            pval_str = f'p < 0.001'
        elif pval < 0.05:
            pval_str = f'p = {pval:.3f}'
        else:
            pval_str = f'p = {pval:.2f}'

        ax.text(0.5, 0.95, pval_str, transform=ax.transAxes,
                ha='center', va='top', fontsize=10)

    plt.suptitle('ESBL Gene Abundance: Week 1 vs Week 3', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_prefix / 'esbl_by_week.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'esbl_by_week.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: esbl_by_week.pdf/png")

    # Figure 5: Prevalence heatmap
    fig, ax = plt.subplots(figsize=(10, 6))

    # Calculate prevalence by location and site
    prevalence_data = []
    for location in ['UCMC', 'ZCH']:
        for site in ['Axilla', 'Groin', 'Stool']:
            subset = rpm_matrix[(rpm_matrix['Location'] == location) &
                               (rpm_matrix['SampleType'] == site)]
            if len(subset) > 0:
                for gene_family in ['CTX-M_total', 'SHV_total', 'TEM_total']:
                    prev = (subset[gene_family] > 0).mean() * 100
                    prevalence_data.append({
                        'Location': location,
                        'Site': site,
                        'Gene': gene_family.replace('_total', ''),
                        'Prevalence': prev
                    })

    prev_df = pd.DataFrame(prevalence_data)
    prev_pivot = prev_df.pivot_table(index=['Location', 'Site'], columns='Gene', values='Prevalence')

    sns.heatmap(prev_pivot, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': 'Prevalence (%)'})
    ax.set_title('ESBL Gene Prevalence by Location and Body Site', fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')

    plt.tight_layout()
    plt.savefig(output_prefix / 'esbl_prevalence_heatmap.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'esbl_prevalence_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: esbl_prevalence_heatmap.pdf/png")

def analyze_top_ctxm_variants(combined_data, metadata):
    """Analyze individual CTX-M variants"""
    print("\n" + "="*60)
    print("TOP CTX-M VARIANTS")
    print("="*60)

    # Get CTX-M data
    ctxm_data = combined_data[combined_data['gene_name'].str.startswith('blaCTX-M')].copy()

    # Sum RPM per variant
    variant_totals = ctxm_data.groupby('gene_name')['rpm'].sum().sort_values(ascending=False)

    print("\nTop 20 CTX-M variants by total abundance:")
    for i, (variant, total_rpm) in enumerate(variant_totals.head(20).items(), 1):
        # Prevalence
        n_samples = ctxm_data[ctxm_data['gene_name'] == variant]['sample_name'].nunique()
        print(f"  {i}. {variant}: {total_rpm:.2f} total RPM, {n_samples} samples")

    # Compare top variants by location
    print("\nTop CTX-M variants by location:")

    # Merge with metadata
    ctxm_data = ctxm_data.merge(metadata, on='sample_name', how='left')

    top_variants = variant_totals.head(10).index.tolist()

    results = []
    for variant in top_variants:
        var_data = ctxm_data[ctxm_data['gene_name'] == variant]

        ucmc_rpm = var_data[var_data['Location'] == 'UCMC']['rpm'].sum()
        zch_rpm = var_data[var_data['Location'] == 'ZCH']['rpm'].sum()

        ucmc_samples = var_data[var_data['Location'] == 'UCMC']['sample_name'].nunique()
        zch_samples = var_data[var_data['Location'] == 'ZCH']['sample_name'].nunique()

        results.append({
            'Variant': variant,
            'UCMC_total_RPM': ucmc_rpm,
            'ZCH_total_RPM': zch_rpm,
            'UCMC_samples': ucmc_samples,
            'ZCH_samples': zch_samples
        })

    return pd.DataFrame(results)

def main():
    print("="*60)
    print("ESBL GENE ANALYSIS - NICU RESISTOME")
    print("="*60)
    print("\nAnalyzing Extended-Spectrum Beta-Lactamase genes")
    print("CTX-M, SHV, TEM variants and carbapenemases")

    # Load data
    combined_data, metadata, master_metadata = load_data()

    # Identify ESBL genes
    esbl_genes, carbapenemase_genes = identify_esbl_genes(combined_data)

    # Create abundance matrix
    rpm_matrix, all_target_genes = create_esbl_abundance_matrix(
        combined_data, esbl_genes, carbapenemase_genes, metadata
    )

    # Calculate gene family totals
    rpm_matrix = calculate_gene_family_totals(rpm_matrix, esbl_genes, carbapenemase_genes)

    # Save the ESBL abundance matrix
    rpm_matrix.to_csv(RESULTS_DIR / 'esbl_abundance_matrix.tsv', sep='\t', index=False)
    print(f"\n✓ Saved: {RESULTS_DIR / 'esbl_abundance_matrix.tsv'}")

    # Run comparisons
    location_results = compare_locations(rpm_matrix)
    location_results.to_csv(RESULTS_DIR / 'esbl_ucmc_vs_zch.tsv', sep='\t', index=False)
    print(f"✓ Saved: {RESULTS_DIR / 'esbl_ucmc_vs_zch.tsv'}")

    bodysite_results = compare_body_sites(rpm_matrix)
    bodysite_results.to_csv(RESULTS_DIR / 'esbl_by_bodysite.tsv', sep='\t', index=False)
    print(f"✓ Saved: {RESULTS_DIR / 'esbl_by_bodysite.tsv'}")

    week_results = compare_weeks(rpm_matrix)
    week_results.to_csv(RESULTS_DIR / 'esbl_week1_vs_week3.tsv', sep='\t', index=False)
    print(f"✓ Saved: {RESULTS_DIR / 'esbl_week1_vs_week3.tsv'}")

    abx_results = analyze_antibiotic_associations(rpm_matrix, master_metadata)
    if len(abx_results) > 0:
        abx_results.to_csv(RESULTS_DIR / 'esbl_antibiotic_associations.tsv', sep='\t', index=False)
        print(f"✓ Saved: {RESULTS_DIR / 'esbl_antibiotic_associations.tsv'}")

    # Analyze top CTX-M variants
    ctxm_variants = analyze_top_ctxm_variants(combined_data, metadata)
    ctxm_variants.to_csv(RESULTS_DIR / 'top_ctxm_variants.tsv', sep='\t', index=False)
    print(f"✓ Saved: {RESULTS_DIR / 'top_ctxm_variants.tsv'}")

    # Generate figures
    plot_esbl_comparison(rpm_matrix, FIGURES_DIR)

    # Final summary
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nResults saved to: {RESULTS_DIR}")
    print(f"Figures saved to: {FIGURES_DIR}")

    # Key findings summary
    print("\n" + "="*60)
    print("KEY FINDINGS SUMMARY")
    print("="*60)

    print("\n1. CTX-M genes:")
    ctxm_prev_ucmc = (rpm_matrix[rpm_matrix['Location'] == 'UCMC']['CTX-M_total'] > 0).mean() * 100
    ctxm_prev_zch = (rpm_matrix[rpm_matrix['Location'] == 'ZCH']['CTX-M_total'] > 0).mean() * 100
    print(f"   UCMC prevalence: {ctxm_prev_ucmc:.1f}%")
    print(f"   ZCH prevalence: {ctxm_prev_zch:.1f}%")

    print("\n2. Total ESBL genes:")
    esbl_prev_ucmc = (rpm_matrix[rpm_matrix['Location'] == 'UCMC']['ESBL_total'] > 0).mean() * 100
    esbl_prev_zch = (rpm_matrix[rpm_matrix['Location'] == 'ZCH']['ESBL_total'] > 0).mean() * 100
    print(f"   UCMC prevalence: {esbl_prev_ucmc:.1f}%")
    print(f"   ZCH prevalence: {esbl_prev_zch:.1f}%")

if __name__ == "__main__":
    main()
