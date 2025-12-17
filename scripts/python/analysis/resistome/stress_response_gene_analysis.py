#!/usr/bin/env python3
"""
Stress Response Gene Analysis for NICU Resistome Study

Analyzes bacterial stress response genes including:
- Oxidative stress: sodA/B/C, katA/B/E/G, ahpC/F, oxyR, soxR/S, trxA/B, grxA/B
- Heat shock: dnaK, dnaJ, groEL, groES, clpA/B/P/X, hspR
- Cold shock: cspA/B/C/D
- General stress: rpoS, relA, spoT, uspA-G
- Acid stress: gadA/B/C
- DNA repair: recA/B/C/N
- Envelope/membrane stress: rpoE, pspA

Comparisons:
1. UCMC vs ZCH (location)
2. Week 1 vs Week 3 (temporal)
3. Association with antibiotic exposure (amount and type)
4. Body site differences

Output: Publication-ready figures and statistical tables
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import mannwhitneyu, kruskal, spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Define paths
DIAMOND_DIR = Path("/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis/diamond_traditional")
DATA_DIR = DIAMOND_DIR / "data"
RESULTS_DIR = DIAMOND_DIR / "results" / "stress_response_analysis"
FIGURES_DIR = DIAMOND_DIR / "figures" / "stress_response"

# Create directories
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Define stress response gene categories
STRESS_RESPONSE_GENES = {
    'Oxidative_stress': {
        'description': 'Protection against reactive oxygen species',
        'genes': ['sodA', 'sodB', 'sodC', 'sodC1', 'sodF',  # Superoxide dismutase
                  'katA', 'katB', 'katE', 'katG', 'katP',    # Catalase
                  'ahpC', 'ahpF',                            # Alkyl hydroperoxide reductase
                  'oxyR',                                     # H2O2 sensor/regulator
                  'soxR', 'soxS',                            # Superoxide stress regulators
                  'trxA', 'trxB', 'trxC',                    # Thioredoxin system
                  'tpx', 'bcp',                              # Peroxiredoxins
                  'gorA', 'gorB']                            # Glutathione reductase
    },
    'Heat_shock': {
        'description': 'Protein folding/refolding under heat stress',
        'genes': ['dnaK', 'dnaJ',                           # Hsp70 chaperone system
                  'groEL', 'groES',                          # Hsp60 chaperonin
                  'clpA', 'clpB', 'clpK', 'clpL', 'clpP', 'clpX',  # Clp proteases/chaperones
                  'rpoH']                                    # Heat shock sigma factor
    },
    'Cold_shock': {
        'description': 'Cold shock response proteins',
        'genes': ['cspA', 'cspB', 'cspC', 'cspD']           # Cold shock proteins
    },
    'General_stress': {
        'description': 'General stress response and starvation',
        'genes': ['rpoS',                                    # Stationary phase sigma factor
                  'relA', 'relB', 'relE', 'spoT',           # Stringent response
                  'uspA', 'uspB', 'uspC', 'uspD', 'uspE', 'uspF', 'uspG']  # Universal stress proteins
    },
    'Acid_stress': {
        'description': 'Acid resistance and pH homeostasis',
        'genes': ['gadA', 'gadB', 'gadC']                   # Glutamate decarboxylase system
    },
    'DNA_repair': {
        'description': 'DNA damage repair (SOS response)',
        'genes': ['recA', 'recB', 'recC', 'recN']           # Recombination/repair
    },
    'Envelope_stress': {
        'description': 'Cell envelope stress response',
        'genes': ['rpoE', 'rpoD', 'rpoN',                   # Sigma factors
                  'sigA']                                    # Sigma factor
    },
    'Efflux_stress': {
        'description': 'Mar/Sox/Rob regulon (multi-drug efflux)',
        'genes': ['marA', 'marR', 'rob', 'ramR']            # Transcriptional regulators
    }
}

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

def identify_stress_genes(combined_data):
    """Identify all stress response genes in the dataset"""
    print("\nIdentifying stress response genes...")

    all_genes = combined_data['gene_name'].unique()
    print(f"  Total unique genes in dataset: {len(all_genes)}")

    found_genes = {}
    gene_to_category = {}

    for category, info in STRESS_RESPONSE_GENES.items():
        category_genes = []
        for target_gene in info['genes']:
            # Match exact gene name or gene name with variants
            matches = [g for g in all_genes if g.lower() == target_gene.lower() or
                      g.lower().startswith(target_gene.lower() + '-') or
                      g.lower().startswith(target_gene.lower() + '_')]
            category_genes.extend(matches)
            for match in matches:
                gene_to_category[match] = category

        found_genes[category] = list(set(category_genes))
        if found_genes[category]:
            print(f"  {category}: {len(found_genes[category])} genes")
            if len(found_genes[category]) <= 10:
                print(f"    {found_genes[category]}")

    return found_genes, gene_to_category

def create_stress_abundance_matrix(combined_data, found_genes, metadata):
    """Create sample × stress gene abundance matrix"""
    print("\nCreating stress response gene abundance matrix...")

    # Combine all stress genes
    all_stress_genes = []
    for gene_list in found_genes.values():
        all_stress_genes.extend(gene_list)
    all_stress_genes = list(set(all_stress_genes))

    print(f"  Total stress genes to analyze: {len(all_stress_genes)}")

    # Filter combined data to stress genes
    stress_data = combined_data[combined_data['gene_name'].isin(all_stress_genes)].copy()
    print(f"  Entries for stress genes: {len(stress_data):,}")

    # Create pivot table (samples × genes with RPM values)
    # Sum multiple hits per gene per sample
    rpm_matrix = stress_data.pivot_table(
        index='sample_name',
        columns='gene_name',
        values='rpm',
        aggfunc='sum',
        fill_value=0
    )

    print(f"  Stress gene RPM matrix: {rpm_matrix.shape[0]} samples × {rpm_matrix.shape[1]} genes")

    # Merge with metadata
    rpm_matrix = rpm_matrix.reset_index()
    rpm_matrix = rpm_matrix.merge(metadata, on='sample_name', how='inner')

    return rpm_matrix, all_stress_genes

def calculate_category_totals(rpm_matrix, found_genes):
    """Calculate total RPM for each stress response category"""
    print("\nCalculating category totals...")

    for category, genes in found_genes.items():
        category_cols = [c for c in rpm_matrix.columns if c in genes]
        if category_cols:
            rpm_matrix[f'{category}_total'] = rpm_matrix[category_cols].sum(axis=1)
            print(f"  {category}: {len(category_cols)} genes")
        else:
            rpm_matrix[f'{category}_total'] = 0

    # Total stress response genes
    category_totals = [f'{cat}_total' for cat in found_genes.keys()]
    category_totals = [c for c in category_totals if c in rpm_matrix.columns]
    rpm_matrix['Total_stress_response'] = rpm_matrix[category_totals].sum(axis=1)

    return rpm_matrix

def compare_locations(rpm_matrix, found_genes):
    """Compare stress gene abundance between UCMC and ZCH"""
    print("\n" + "="*60)
    print("STRESS RESPONSE GENES: UCMC vs ZCH")
    print("="*60)

    results = []
    categories = list(found_genes.keys()) + ['Total_stress_response']

    ucmc = rpm_matrix[rpm_matrix['Location'] == 'UCMC']
    zch = rpm_matrix[rpm_matrix['Location'] == 'ZCH']

    print(f"\nSamples: UCMC = {len(ucmc)}, ZCH = {len(zch)}")

    for category in categories:
        col = f'{category}_total' if category != 'Total_stress_response' else category

        if col not in rpm_matrix.columns:
            continue

        ucmc_vals = ucmc[col]
        zch_vals = zch[col]

        # Statistics
        stat, pval = mannwhitneyu(ucmc_vals, zch_vals, alternative='two-sided')

        # Effect size (log2 fold change of medians)
        ucmc_median = ucmc_vals.median()
        zch_median = zch_vals.median()

        if ucmc_median > 0 and zch_median > 0:
            log2fc = np.log2(zch_median / ucmc_median)  # Positive = higher in ZCH
        elif zch_median > 0:
            log2fc = np.inf
        elif ucmc_median > 0:
            log2fc = -np.inf
        else:
            log2fc = 0

        # Prevalence
        ucmc_prev = (ucmc_vals > 0).mean() * 100
        zch_prev = (zch_vals > 0).mean() * 100

        results.append({
            'Category': category.replace('_total', ''),
            'UCMC_median_RPM': ucmc_median,
            'ZCH_median_RPM': zch_median,
            'UCMC_mean_RPM': ucmc_vals.mean(),
            'ZCH_mean_RPM': zch_vals.mean(),
            'UCMC_prevalence_%': ucmc_prev,
            'ZCH_prevalence_%': zch_prev,
            'log2FC_ZCH_vs_UCMC': log2fc,
            'Mann_Whitney_U': stat,
            'p_value': pval
        })

        sig = '*' if pval < 0.05 else ''
        print(f"\n{category.replace('_', ' ')}:{sig}")
        print(f"  UCMC: median={ucmc_median:.3f}, mean={ucmc_vals.mean():.3f}")
        print(f"  ZCH:  median={zch_median:.3f}, mean={zch_vals.mean():.3f}")
        print(f"  log2FC (ZCH/UCMC): {log2fc:.2f}, p={pval:.2e}")

    results_df = pd.DataFrame(results)

    # FDR correction
    _, pvals_corrected, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
    results_df['p_value_FDR'] = pvals_corrected

    return results_df

def compare_weeks(rpm_matrix, found_genes):
    """Compare stress gene abundance between Week 1 and Week 3"""
    print("\n" + "="*60)
    print("STRESS RESPONSE GENES: Week 1 vs Week 3")
    print("="*60)

    results = []
    categories = list(found_genes.keys()) + ['Total_stress_response']

    week1 = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.1']
    week3 = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.3']

    print(f"\nSamples: Week 1 = {len(week1)}, Week 3 = {len(week3)}")

    for category in categories:
        col = f'{category}_total' if category != 'Total_stress_response' else category

        if col not in rpm_matrix.columns:
            continue

        w1_vals = week1[col]
        w3_vals = week3[col]

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
            'Category': category.replace('_total', ''),
            'Week1_median_RPM': w1_median,
            'Week3_median_RPM': w3_median,
            'Week1_prevalence_%': (w1_vals > 0).mean() * 100,
            'Week3_prevalence_%': (w3_vals > 0).mean() * 100,
            'log2FC_W3_vs_W1': log2fc,
            'Mann_Whitney_U': stat,
            'p_value': pval
        })

        sig = '*' if pval < 0.05 else ''
        print(f"\n{category.replace('_', ' ')}:{sig}")
        print(f"  Week 1: median={w1_median:.3f}")
        print(f"  Week 3: median={w3_median:.3f}")
        print(f"  log2FC (W3/W1): {log2fc:.2f}, p={pval:.2e}")

    results_df = pd.DataFrame(results)
    _, pvals_corrected, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
    results_df['p_value_FDR'] = pvals_corrected

    return results_df

def analyze_antibiotic_associations(rpm_matrix, master_metadata, found_genes):
    """Correlate stress genes with antibiotic exposure amount and type"""
    print("\n" + "="*60)
    print("STRESS RESPONSE ASSOCIATION WITH ANTIBIOTIC EXPOSURE")
    print("="*60)

    # Fix ZCH sample names in master metadata
    master_meta = master_metadata.copy()
    zch_mask = master_meta['Location'] == 'ZCH'
    master_meta.loc[zch_mask, 'sample_name'] = 'ZJH_' + master_meta.loc[zch_mask, 'sample_name']

    # Get antibiotic columns
    abx_cols = [c for c in master_meta.columns if c.endswith('_w1') or c.endswith('_w2') or c.endswith('_w3')]
    abx_cols = [c for c in abx_cols if not c.startswith('None')]

    # Merge with stress data
    merged = rpm_matrix.merge(
        master_meta[['sample_name', 'SampleCollectionWeek'] + abx_cols],
        on='sample_name', how='inner', suffixes=('', '_meta')
    )

    print(f"\nSamples with antibiotic data: {len(merged)}")

    # Define antibiotic classes
    antibiotics = {
        'Ampicillin': 'Beta-lactam',
        'Nafcillin': 'Beta-lactam',
        'Penicillin': 'Beta-lactam',
        'Meropenem': 'Carbapenem',
        'Cefotaxime': 'Cephalosporin',
        'Moxalactam': 'Cephalosporin',
        'Ceftazidime': 'Cephalosporin',
        'Gentamicin': 'Aminoglycoside',
        'Tobramycin': 'Aminoglycoside',
        'Vancomycin': 'Glycopeptide',
        'Fluconazole': 'Antifungal',
        'Azithromycin': 'Macrolide',
        'Piperacillin_Tazobactam': 'Beta-lactam/Inhibitor'
    }

    # Calculate cumulative antibiotic exposure for each sample
    week_col = 'SampleCollectionWeek_meta' if 'SampleCollectionWeek_meta' in merged.columns else 'SampleCollectionWeek'

    for abx in antibiotics.keys():
        w1_col = f'{abx}_w1'
        w2_col = f'{abx}_w2'
        w3_col = f'{abx}_w3'

        merged[f'{abx}_cumulative'] = 0

        if w1_col in merged.columns:
            merged[f'{abx}_cumulative'] += merged[w1_col].fillna(0)

        # Add week 2 for Week 3 samples
        if w2_col in merged.columns:
            week3_mask = merged[week_col] == 'Week.3'
            merged.loc[week3_mask, f'{abx}_cumulative'] += merged.loc[week3_mask, w2_col].fillna(0)

        # Add week 3 for Week 3 samples (if collected at Week 3)
        if w3_col in merged.columns:
            week3_mask = merged[week_col] == 'Week.3'
            merged.loc[week3_mask, f'{abx}_cumulative'] += merged.loc[week3_mask, w3_col].fillna(0)

    # Calculate total antibiotic exposure
    cumulative_cols = [f'{abx}_cumulative' for abx in antibiotics.keys() if f'{abx}_cumulative' in merged.columns]
    merged['Total_abx_days'] = merged[cumulative_cols].sum(axis=1)

    results = []
    categories = list(found_genes.keys()) + ['Total_stress_response']

    # Correlate total ABX exposure with each stress category
    print("\n--- Correlation with TOTAL antibiotic exposure ---")
    for category in categories:
        col = f'{category}_total' if category != 'Total_stress_response' else category
        if col not in merged.columns:
            continue

        stress_vals = merged[col]
        abx_vals = merged['Total_abx_days']

        # Only include samples with some stress gene expression
        mask = (stress_vals > 0) | (abx_vals > 0)

        if mask.sum() >= 10:
            rho, pval = spearmanr(abx_vals[mask], stress_vals[mask])

            # Compare exposed vs unexposed
            exposed = stress_vals[abx_vals > 0]
            unexposed = stress_vals[abx_vals == 0]

            if len(exposed) >= 5 and len(unexposed) >= 5:
                mw_stat, mw_pval = mannwhitneyu(exposed, unexposed, alternative='two-sided')
            else:
                mw_stat, mw_pval = np.nan, np.nan

            results.append({
                'Antibiotic': 'Total_exposure',
                'Antibiotic_class': 'All',
                'Stress_category': category,
                'n_exposed': (abx_vals > 0).sum(),
                'n_unexposed': (abx_vals == 0).sum(),
                'Spearman_rho': rho,
                'Spearman_p': pval,
                'Exposed_median_RPM': exposed.median() if len(exposed) > 0 else 0,
                'Unexposed_median_RPM': unexposed.median() if len(unexposed) > 0 else 0,
                'Mann_Whitney_p': mw_pval
            })

            sig = '*' if pval < 0.05 else ''
            print(f"  {category}: rho={rho:.3f}, p={pval:.3e}{sig}")

    # Correlate individual antibiotics
    print("\n--- Correlation with SPECIFIC antibiotics ---")
    for abx, abx_class in antibiotics.items():
        cum_col = f'{abx}_cumulative'
        if cum_col not in merged.columns:
            continue

        abx_exposure = merged[cum_col]
        n_exposed = (abx_exposure > 0).sum()

        if n_exposed < 10:
            continue

        print(f"\n{abx} ({abx_class}), n_exposed={n_exposed}:")

        for category in categories:
            col = f'{category}_total' if category != 'Total_stress_response' else category
            if col not in merged.columns:
                continue

            stress_vals = merged[col]

            # Spearman correlation
            rho, pval = spearmanr(abx_exposure, stress_vals)

            # Compare exposed vs unexposed
            exposed = stress_vals[abx_exposure > 0]
            unexposed = stress_vals[abx_exposure == 0]

            if len(exposed) >= 5 and len(unexposed) >= 5:
                mw_stat, mw_pval = mannwhitneyu(exposed, unexposed, alternative='two-sided')
            else:
                mw_stat, mw_pval = np.nan, np.nan

            results.append({
                'Antibiotic': abx,
                'Antibiotic_class': abx_class,
                'Stress_category': category,
                'n_exposed': n_exposed,
                'n_unexposed': len(unexposed),
                'Spearman_rho': rho,
                'Spearman_p': pval,
                'Exposed_median_RPM': exposed.median() if len(exposed) > 0 else 0,
                'Unexposed_median_RPM': unexposed.median() if len(unexposed) > 0 else 0,
                'Mann_Whitney_p': mw_pval
            })

            if category in ['Oxidative_stress', 'Heat_shock', 'Total_stress_response']:
                sig = '*' if pval < 0.05 else ''
                print(f"    {category}: rho={rho:.3f}, p={pval:.3e}{sig}")

    return pd.DataFrame(results), merged

def analyze_week_antibiotic_interaction(merged, found_genes):
    """Analyze stress gene differences by week stratified by antibiotic exposure"""
    print("\n" + "="*60)
    print("STRESS RESPONSE BY WEEK × ANTIBIOTIC EXPOSURE")
    print("="*60)

    results = []
    categories = list(found_genes.keys()) + ['Total_stress_response']

    week_col = 'SampleCollectionWeek_meta' if 'SampleCollectionWeek_meta' in merged.columns else 'SampleCollectionWeek'

    for category in categories:
        col = f'{category}_total' if category != 'Total_stress_response' else category
        if col not in merged.columns:
            continue

        print(f"\n{category.replace('_', ' ')}:")

        # Week 1 exposed vs unexposed
        w1 = merged[merged[week_col] == 'Week.1']
        w1_exposed = w1[w1['Total_abx_days'] > 0][col]
        w1_unexposed = w1[w1['Total_abx_days'] == 0][col]

        # Week 3 exposed vs unexposed
        w3 = merged[merged[week_col] == 'Week.3']
        w3_exposed = w3[w3['Total_abx_days'] > 0][col]
        w3_unexposed = w3[w3['Total_abx_days'] == 0][col]

        print(f"  Week 1: exposed n={len(w1_exposed)}, unexposed n={len(w1_unexposed)}")
        print(f"  Week 3: exposed n={len(w3_exposed)}, unexposed n={len(w3_unexposed)}")

        # Compare exposed vs unexposed at each week
        if len(w1_exposed) >= 5 and len(w1_unexposed) >= 5:
            _, p_w1 = mannwhitneyu(w1_exposed, w1_unexposed, alternative='two-sided')
            print(f"  Week 1 exposed vs unexposed: p={p_w1:.3e}")
        else:
            p_w1 = np.nan

        if len(w3_exposed) >= 5 and len(w3_unexposed) >= 5:
            _, p_w3 = mannwhitneyu(w3_exposed, w3_unexposed, alternative='two-sided')
            print(f"  Week 3 exposed vs unexposed: p={p_w3:.3e}")
        else:
            p_w3 = np.nan

        # Compare temporal change in exposed samples
        if len(w1_exposed) >= 5 and len(w3_exposed) >= 5:
            _, p_exposed_temporal = mannwhitneyu(w1_exposed, w3_exposed, alternative='two-sided')
            print(f"  Exposed: Week 1 vs Week 3: p={p_exposed_temporal:.3e}")
        else:
            p_exposed_temporal = np.nan

        # Compare temporal change in unexposed samples
        if len(w1_unexposed) >= 5 and len(w3_unexposed) >= 5:
            _, p_unexposed_temporal = mannwhitneyu(w1_unexposed, w3_unexposed, alternative='two-sided')
            print(f"  Unexposed: Week 1 vs Week 3: p={p_unexposed_temporal:.3e}")
        else:
            p_unexposed_temporal = np.nan

        results.append({
            'Category': category,
            'Week1_exposed_n': len(w1_exposed),
            'Week1_unexposed_n': len(w1_unexposed),
            'Week3_exposed_n': len(w3_exposed),
            'Week3_unexposed_n': len(w3_unexposed),
            'Week1_exposed_median': w1_exposed.median() if len(w1_exposed) > 0 else np.nan,
            'Week1_unexposed_median': w1_unexposed.median() if len(w1_unexposed) > 0 else np.nan,
            'Week3_exposed_median': w3_exposed.median() if len(w3_exposed) > 0 else np.nan,
            'Week3_unexposed_median': w3_unexposed.median() if len(w3_unexposed) > 0 else np.nan,
            'p_week1_exp_vs_unexp': p_w1,
            'p_week3_exp_vs_unexp': p_w3,
            'p_exposed_w1_vs_w3': p_exposed_temporal,
            'p_unexposed_w1_vs_w3': p_unexposed_temporal
        })

    return pd.DataFrame(results)

def analyze_individual_genes(rpm_matrix, found_genes, master_metadata):
    """Analyze individual stress response genes"""
    print("\n" + "="*60)
    print("TOP DIFFERENTIALLY ABUNDANT INDIVIDUAL GENES")
    print("="*60)

    results = []

    for category, genes in found_genes.items():
        for gene in genes:
            if gene not in rpm_matrix.columns:
                continue

            # Location comparison
            ucmc_vals = rpm_matrix[rpm_matrix['Location'] == 'UCMC'][gene]
            zch_vals = rpm_matrix[rpm_matrix['Location'] == 'ZCH'][gene]

            if (ucmc_vals > 0).sum() < 5 and (zch_vals > 0).sum() < 5:
                continue

            stat, loc_pval = mannwhitneyu(ucmc_vals, zch_vals, alternative='two-sided')

            ucmc_median = ucmc_vals.median()
            zch_median = zch_vals.median()

            if ucmc_median > 0 and zch_median > 0:
                log2fc_loc = np.log2(zch_median / ucmc_median)
            elif zch_median > 0:
                log2fc_loc = 5  # Cap for display
            elif ucmc_median > 0:
                log2fc_loc = -5
            else:
                log2fc_loc = 0

            # Week comparison
            w1_vals = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.1'][gene]
            w3_vals = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.3'][gene]

            _, week_pval = mannwhitneyu(w1_vals, w3_vals, alternative='two-sided')

            w1_median = w1_vals.median()
            w3_median = w3_vals.median()

            if w1_median > 0 and w3_median > 0:
                log2fc_week = np.log2(w3_median / w1_median)
            elif w3_median > 0:
                log2fc_week = 5
            elif w1_median > 0:
                log2fc_week = -5
            else:
                log2fc_week = 0

            results.append({
                'Gene': gene,
                'Category': category,
                'UCMC_median': ucmc_median,
                'ZCH_median': zch_median,
                'UCMC_prev%': (ucmc_vals > 0).mean() * 100,
                'ZCH_prev%': (zch_vals > 0).mean() * 100,
                'log2FC_ZCH_vs_UCMC': log2fc_loc,
                'p_location': loc_pval,
                'Week1_median': w1_median,
                'Week3_median': w3_median,
                'log2FC_W3_vs_W1': log2fc_week,
                'p_week': week_pval
            })

    results_df = pd.DataFrame(results)

    if len(results_df) > 0:
        # FDR correction
        _, loc_fdr, _, _ = multipletests(results_df['p_location'], method='fdr_bh')
        _, week_fdr, _, _ = multipletests(results_df['p_week'], method='fdr_bh')
        results_df['p_location_FDR'] = loc_fdr
        results_df['p_week_FDR'] = week_fdr

        # Print top differentially abundant genes
        print("\nTop genes by location (ZCH vs UCMC):")
        top_loc = results_df.nsmallest(10, 'p_location')
        for _, row in top_loc.iterrows():
            print(f"  {row['Gene']} ({row['Category']}): log2FC={row['log2FC_ZCH_vs_UCMC']:.2f}, p={row['p_location']:.2e}")

        print("\nTop genes by week (Week 3 vs Week 1):")
        top_week = results_df.nsmallest(10, 'p_week')
        for _, row in top_week.iterrows():
            print(f"  {row['Gene']} ({row['Category']}): log2FC={row['log2FC_W3_vs_W1']:.2f}, p={row['p_week']:.2e}")

    return results_df

def plot_stress_comparison(rpm_matrix, found_genes, output_prefix):
    """Create publication-ready figures for stress gene analysis"""
    print("\n" + "="*60)
    print("GENERATING FIGURES")
    print("="*60)

    plt.style.use('seaborn-v0_8-whitegrid')

    # Figure 1: Stress categories by location
    categories = [cat for cat in found_genes.keys() if f'{cat}_total' in rpm_matrix.columns]
    n_cats = len(categories)

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for i, category in enumerate(categories[:8]):
        ax = axes[i]
        col = f'{category}_total'

        plot_data = rpm_matrix.copy()
        plot_data[col] = plot_data[col].replace(0, 0.001)

        sns.boxplot(data=plot_data, x='Location', y=col, ax=ax,
                    order=['UCMC', 'ZCH'], palette=['#2ecc71', '#e74c3c'])
        sns.stripplot(data=plot_data, x='Location', y=col, ax=ax,
                      order=['UCMC', 'ZCH'], color='black', alpha=0.3, size=3)

        ax.set_ylabel('RPM' if i % 4 == 0 else '')
        ax.set_xlabel('')
        ax.set_title(category.replace('_', ' '), fontsize=10, fontweight='bold')
        ax.set_yscale('log')

        # P-value
        ucmc_vals = rpm_matrix[rpm_matrix['Location'] == 'UCMC'][col]
        zch_vals = rpm_matrix[rpm_matrix['Location'] == 'ZCH'][col]
        _, pval = mannwhitneyu(ucmc_vals, zch_vals, alternative='two-sided')

        if pval < 0.001:
            pval_str = 'p < 0.001'
        elif pval < 0.05:
            pval_str = f'p = {pval:.3f}'
        else:
            pval_str = f'p = {pval:.2f}'

        ax.text(0.5, 0.95, pval_str, transform=ax.transAxes,
                ha='center', va='top', fontsize=9)

    # Hide unused axes
    for i in range(len(categories), 8):
        axes[i].set_visible(False)

    plt.suptitle('Stress Response Gene Categories: UCMC vs ZCH', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_prefix / 'stress_by_location.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'stress_by_location.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: stress_by_location.pdf/png")

    # Figure 2: Stress categories by week
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for i, category in enumerate(categories[:8]):
        ax = axes[i]
        col = f'{category}_total'

        plot_data = rpm_matrix.copy()
        plot_data[col] = plot_data[col].replace(0, 0.001)

        sns.boxplot(data=plot_data, x='SampleCollectionWeek', y=col, ax=ax,
                    order=['Week.1', 'Week.3'], palette=['#3498db', '#9b59b6'])
        sns.stripplot(data=plot_data, x='SampleCollectionWeek', y=col, ax=ax,
                      order=['Week.1', 'Week.3'], color='black', alpha=0.3, size=3)

        ax.set_ylabel('RPM' if i % 4 == 0 else '')
        ax.set_xlabel('')
        ax.set_xticklabels(['Week 1', 'Week 3'])
        ax.set_title(category.replace('_', ' '), fontsize=10, fontweight='bold')
        ax.set_yscale('log')

        # P-value
        w1_vals = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.1'][col]
        w3_vals = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.3'][col]
        _, pval = mannwhitneyu(w1_vals, w3_vals, alternative='two-sided')

        if pval < 0.001:
            pval_str = 'p < 0.001'
        elif pval < 0.05:
            pval_str = f'p = {pval:.3f}'
        else:
            pval_str = f'p = {pval:.2f}'

        ax.text(0.5, 0.95, pval_str, transform=ax.transAxes,
                ha='center', va='top', fontsize=9)

    for i in range(len(categories), 8):
        axes[i].set_visible(False)

    plt.suptitle('Stress Response Gene Categories: Week 1 vs Week 3', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_prefix / 'stress_by_week.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'stress_by_week.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: stress_by_week.pdf/png")

    # Figure 3: Heatmap of stress gene category prevalence
    fig, ax = plt.subplots(figsize=(12, 6))

    prevalence_data = []
    for location in ['UCMC', 'ZCH']:
        for week in ['Week.1', 'Week.3']:
            subset = rpm_matrix[(rpm_matrix['Location'] == location) &
                               (rpm_matrix['SampleCollectionWeek'] == week)]
            if len(subset) > 0:
                for category in categories:
                    col = f'{category}_total'
                    prev = (subset[col] > 0).mean() * 100
                    prevalence_data.append({
                        'Location': location,
                        'Week': week.replace('.', ' '),
                        'Category': category.replace('_', '\n'),
                        'Prevalence': prev
                    })

    prev_df = pd.DataFrame(prevalence_data)
    prev_df['Group'] = prev_df['Location'] + '\n' + prev_df['Week']
    prev_pivot = prev_df.pivot_table(index='Group', columns='Category', values='Prevalence')

    sns.heatmap(prev_pivot, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': 'Prevalence (%)'})
    ax.set_title('Stress Response Gene Prevalence', fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')

    plt.tight_layout()
    plt.savefig(output_prefix / 'stress_prevalence_heatmap.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'stress_prevalence_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: stress_prevalence_heatmap.pdf/png")

    # Figure 4: Total stress response by location and week (faceted)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Total by location
    ax = axes[0]
    plot_data = rpm_matrix.copy()
    plot_data['Total_stress_response'] = plot_data['Total_stress_response'].replace(0, 0.001)

    sns.boxplot(data=plot_data, x='Location', y='Total_stress_response', ax=ax,
                order=['UCMC', 'ZCH'], palette=['#2ecc71', '#e74c3c'])
    sns.stripplot(data=plot_data, x='Location', y='Total_stress_response', ax=ax,
                  order=['UCMC', 'ZCH'], color='black', alpha=0.3, size=3)
    ax.set_ylabel('Total Stress Response RPM')
    ax.set_xlabel('')
    ax.set_title('By Location', fontsize=12, fontweight='bold')
    ax.set_yscale('log')

    # Add p-value
    ucmc_vals = rpm_matrix[rpm_matrix['Location'] == 'UCMC']['Total_stress_response']
    zch_vals = rpm_matrix[rpm_matrix['Location'] == 'ZCH']['Total_stress_response']
    _, pval = mannwhitneyu(ucmc_vals, zch_vals, alternative='two-sided')
    pval_str = f'p = {pval:.3e}' if pval < 0.001 else f'p = {pval:.3f}'
    ax.text(0.5, 0.95, pval_str, transform=ax.transAxes, ha='center', va='top', fontsize=10)

    # Total by week
    ax = axes[1]
    sns.boxplot(data=plot_data, x='SampleCollectionWeek', y='Total_stress_response', ax=ax,
                order=['Week.1', 'Week.3'], palette=['#3498db', '#9b59b6'])
    sns.stripplot(data=plot_data, x='SampleCollectionWeek', y='Total_stress_response', ax=ax,
                  order=['Week.1', 'Week.3'], color='black', alpha=0.3, size=3)
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_xticklabels(['Week 1', 'Week 3'])
    ax.set_title('By Week', fontsize=12, fontweight='bold')
    ax.set_yscale('log')

    # Add p-value
    w1_vals = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.1']['Total_stress_response']
    w3_vals = rpm_matrix[rpm_matrix['SampleCollectionWeek'] == 'Week.3']['Total_stress_response']
    _, pval = mannwhitneyu(w1_vals, w3_vals, alternative='two-sided')
    pval_str = f'p = {pval:.3e}' if pval < 0.001 else f'p = {pval:.3f}'
    ax.text(0.5, 0.95, pval_str, transform=ax.transAxes, ha='center', va='top', fontsize=10)

    plt.suptitle('Total Stress Response Gene Abundance', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_prefix / 'total_stress_summary.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'total_stress_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: total_stress_summary.pdf/png")

def plot_antibiotic_correlations(merged, found_genes, output_prefix):
    """Plot stress response correlation with antibiotic exposure"""
    print("\nGenerating antibiotic correlation figures...")

    categories = list(found_genes.keys()) + ['Total_stress_response']

    # Figure: Scatter plot of total antibiotics vs stress response categories
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for i, category in enumerate(categories[:8]):
        ax = axes[i]
        col = f'{category}_total' if category != 'Total_stress_response' else category

        if col not in merged.columns:
            ax.set_visible(False)
            continue

        x = merged['Total_abx_days']
        y = merged[col]

        # Add jitter for visualization
        ax.scatter(x + np.random.normal(0, 0.1, len(x)),
                  y + 0.001, alpha=0.5, s=20)

        # Correlation
        rho, pval = spearmanr(x, y)

        ax.set_xlabel('Total Antibiotic Days')
        ax.set_ylabel('RPM' if i % 4 == 0 else '')
        ax.set_title(f'{category.replace("_", " ")}\nrho={rho:.2f}, p={pval:.2e}',
                    fontsize=9, fontweight='bold')
        ax.set_yscale('log')

        # Color by significance
        if pval < 0.05:
            ax.set_facecolor('#fff8dc')

    for i in range(len(categories), 8):
        if i < len(axes):
            axes[i].set_visible(False)

    plt.suptitle('Stress Response vs Total Antibiotic Exposure', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_prefix / 'stress_vs_antibiotics_scatter.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'stress_vs_antibiotics_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: stress_vs_antibiotics_scatter.pdf/png")

    # Figure: Exposed vs unexposed comparison
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for i, category in enumerate(categories[:8]):
        ax = axes[i]
        col = f'{category}_total' if category != 'Total_stress_response' else category

        if col not in merged.columns:
            ax.set_visible(False)
            continue

        merged['Abx_exposed'] = (merged['Total_abx_days'] > 0).map({True: 'Exposed', False: 'Unexposed'})

        plot_data = merged.copy()
        plot_data[col] = plot_data[col].replace(0, 0.001)

        sns.boxplot(data=plot_data, x='Abx_exposed', y=col, ax=ax,
                    order=['Unexposed', 'Exposed'], palette=['#95a5a6', '#e74c3c'])

        ax.set_ylabel('RPM' if i % 4 == 0 else '')
        ax.set_xlabel('')
        ax.set_title(category.replace('_', ' '), fontsize=10, fontweight='bold')
        ax.set_yscale('log')

        # P-value
        exposed = merged[merged['Total_abx_days'] > 0][col]
        unexposed = merged[merged['Total_abx_days'] == 0][col]

        if len(exposed) >= 5 and len(unexposed) >= 5:
            _, pval = mannwhitneyu(exposed, unexposed, alternative='two-sided')
            pval_str = f'p = {pval:.2e}' if pval < 0.001 else f'p = {pval:.3f}'
            ax.text(0.5, 0.95, pval_str, transform=ax.transAxes, ha='center', va='top', fontsize=9)

    for i in range(len(categories), 8):
        if i < len(axes):
            axes[i].set_visible(False)

    plt.suptitle('Stress Response: Antibiotic Exposed vs Unexposed', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_prefix / 'stress_exposed_vs_unexposed.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_prefix / 'stress_exposed_vs_unexposed.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: stress_exposed_vs_unexposed.pdf/png")

def main():
    print("="*60)
    print("STRESS RESPONSE GENE ANALYSIS - NICU RESISTOME")
    print("="*60)
    print("\nAnalyzing bacterial stress response genes")
    print("Oxidative, heat shock, cold shock, general stress, acid stress")

    # Load data
    combined_data, metadata, master_metadata = load_data()

    # Identify stress genes
    found_genes, gene_to_category = identify_stress_genes(combined_data)

    # Create abundance matrix
    rpm_matrix, all_stress_genes = create_stress_abundance_matrix(
        combined_data, found_genes, metadata
    )

    # Calculate category totals
    rpm_matrix = calculate_category_totals(rpm_matrix, found_genes)

    # Save the stress gene abundance matrix
    rpm_matrix.to_csv(RESULTS_DIR / 'stress_abundance_matrix.tsv', sep='\t', index=False)
    print(f"\n✓ Saved: {RESULTS_DIR / 'stress_abundance_matrix.tsv'}")

    # Run comparisons
    location_results = compare_locations(rpm_matrix, found_genes)
    location_results.to_csv(RESULTS_DIR / 'stress_ucmc_vs_zch.tsv', sep='\t', index=False)
    print(f"✓ Saved: {RESULTS_DIR / 'stress_ucmc_vs_zch.tsv'}")

    week_results = compare_weeks(rpm_matrix, found_genes)
    week_results.to_csv(RESULTS_DIR / 'stress_week1_vs_week3.tsv', sep='\t', index=False)
    print(f"✓ Saved: {RESULTS_DIR / 'stress_week1_vs_week3.tsv'}")

    # Antibiotic associations
    abx_results, merged = analyze_antibiotic_associations(rpm_matrix, master_metadata, found_genes)
    if len(abx_results) > 0:
        abx_results.to_csv(RESULTS_DIR / 'stress_antibiotic_associations.tsv', sep='\t', index=False)
        print(f"✓ Saved: {RESULTS_DIR / 'stress_antibiotic_associations.tsv'}")

    # Week × antibiotic interaction
    week_abx_results = analyze_week_antibiotic_interaction(merged, found_genes)
    week_abx_results.to_csv(RESULTS_DIR / 'stress_week_antibiotic_interaction.tsv', sep='\t', index=False)
    print(f"✓ Saved: {RESULTS_DIR / 'stress_week_antibiotic_interaction.tsv'}")

    # Individual gene analysis
    gene_results = analyze_individual_genes(rpm_matrix, found_genes, master_metadata)
    if len(gene_results) > 0:
        gene_results.to_csv(RESULTS_DIR / 'stress_individual_genes.tsv', sep='\t', index=False)
        print(f"✓ Saved: {RESULTS_DIR / 'stress_individual_genes.tsv'}")

    # Generate figures
    plot_stress_comparison(rpm_matrix, found_genes, FIGURES_DIR)

    if len(merged) > 0:
        plot_antibiotic_correlations(merged, found_genes, FIGURES_DIR)

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

    # Location differences
    print("\n1. Location differences (ZCH vs UCMC):")
    sig_loc = location_results[location_results['p_value'] < 0.05]
    for _, row in sig_loc.iterrows():
        direction = "higher in ZCH" if row['log2FC_ZCH_vs_UCMC'] > 0 else "higher in UCMC"
        print(f"   {row['Category']}: {direction} (log2FC={row['log2FC_ZCH_vs_UCMC']:.2f}, p={row['p_value']:.2e})")

    # Week differences
    print("\n2. Temporal changes (Week 3 vs Week 1):")
    sig_week = week_results[week_results['p_value'] < 0.05]
    for _, row in sig_week.iterrows():
        direction = "increased" if row['log2FC_W3_vs_W1'] > 0 else "decreased"
        print(f"   {row['Category']}: {direction} (log2FC={row['log2FC_W3_vs_W1']:.2f}, p={row['p_value']:.2e})")

    # Antibiotic associations
    print("\n3. Antibiotic associations (p < 0.05):")
    if len(abx_results) > 0:
        sig_abx = abx_results[(abx_results['Spearman_p'] < 0.05) &
                              (abx_results['Stress_category'].isin(['Oxidative_stress', 'Heat_shock', 'General_stress', 'Total_stress_response']))]
        sig_abx = sig_abx.sort_values('Spearman_p')
        for _, row in sig_abx.head(10).iterrows():
            direction = "positive" if row['Spearman_rho'] > 0 else "negative"
            print(f"   {row['Antibiotic']} vs {row['Stress_category']}: {direction} correlation (rho={row['Spearman_rho']:.2f}, p={row['Spearman_p']:.2e})")

if __name__ == "__main__":
    main()
