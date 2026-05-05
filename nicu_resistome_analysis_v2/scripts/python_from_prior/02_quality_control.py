#!/usr/bin/env python3
"""
Quality control analysis:
1. Identify samples with AMR data
2. Check sequencing depth and coverage
3. Filter low-abundance genes
4. Identify complete subjects for paired analysis
5. Generate QC figures
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Define paths
PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
DATA_DIR = PROJECT_ROOT / "nicu_resistome_analysis" / "data"
RESULTS_DIR = PROJECT_ROOT / "nicu_resistome_analysis" / "results" / "qc"
FIGURES_DIR = PROJECT_ROOT / "nicu_resistome_analysis" / "figures" / "exploratory"

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

def load_data():
    """Load master metadata and AMR data"""
    print("Loading data...")

    # Load master metadata
    metadata = pd.read_csv(RESULTS_DIR / "master_metadata.tsv", sep="\t")

    # Load full AMR data
    amr_data = pd.read_csv(DATA_DIR / "nicu_amr_only.tsv", sep="\t")

    # Load sample summary
    sample_summary = pd.read_csv(DATA_DIR / "sample_summary.tsv", sep="\t")

    print(f"Metadata: {len(metadata)} samples")
    print(f"AMR data: {len(amr_data)} entries from {amr_data['sample_name'].nunique()} samples")
    print(f"Sample summary: {len(sample_summary)} samples")

    return metadata, amr_data, sample_summary

def analyze_sequencing_depth(metadata, sample_summary):
    """Analyze sequencing depth and AMR burden"""
    print("\n" + "="*60)
    print("SEQUENCING DEPTH ANALYSIS")
    print("="*60)

    # Merge metadata with sample summary for samples with AMR data
    data = metadata[metadata['has_amr_data']].copy()

    # Summary statistics
    print(f"\nTotal AMR reads per sample:")
    print(data['total_amr_reads'].describe())

    print(f"\nTotal AMR RPM per sample:")
    print(data['total_amr_rpm'].describe())

    print(f"\nNumber of unique AMR genes per sample:")
    print(data['unique_amr_genes'].describe())

    print(f"\nNumber of drug classes per sample:")
    print(data['num_drug_classes'].describe())

    # By location
    print(f"\n\nBy Location:")
    location_stats = data.groupby('Location').agg({
        'total_amr_reads': ['mean', 'median', 'std'],
        'unique_amr_genes': ['mean', 'median'],
        'num_drug_classes': ['mean', 'median']
    })
    print(location_stats)

    # By body site
    print(f"\n\nBy Body Site:")
    site_stats = data.groupby('SampleType').agg({
        'total_amr_reads': ['mean', 'median', 'std'],
        'unique_amr_genes': ['mean', 'median'],
        'num_drug_classes': ['mean', 'median']
    })
    print(site_stats)

    # By week
    print(f"\n\nBy Week:")
    week_stats = data.groupby('SampleCollectionWeek').agg({
        'total_amr_reads': ['mean', 'median', 'std'],
        'unique_amr_genes': ['mean', 'median'],
        'num_drug_classes': ['mean', 'median']
    })
    print(week_stats)

    # Identify low-depth samples
    low_amr_threshold = data['total_amr_rpm'].quantile(0.1)
    low_samples = data[data['total_amr_rpm'] < low_amr_threshold]
    print(f"\n\n⚠ Low AMR abundance samples (bottom 10%, RPM < {low_amr_threshold:.2f}): {len(low_samples)}")

    return data, low_samples

def identify_complete_subjects(metadata):
    """Identify subjects with complete sampling (all 6 samples with AMR data)"""
    print("\n" + "="*60)
    print("SUBJECT COMPLETENESS ANALYSIS")
    print("="*60)

    # Only consider samples with AMR data and valid body sites
    valid_samples = metadata[
        (metadata['has_amr_data']) &
        (metadata['SampleType'].isin(['Axilla', 'Groin', 'Stool']))
    ].copy()

    # Count samples per subject
    subject_counts = valid_samples.groupby('PatientID').agg({
        'Sample': 'count',
        'SampleType': lambda x: x.nunique(),
        'SampleCollectionWeek': lambda x: x.nunique(),
        'Location': 'first'
    }).rename(columns={
        'Sample': 'n_samples',
        'SampleType': 'n_body_sites',
        'SampleCollectionWeek': 'n_weeks'
    })

    # Complete subjects have 6 samples (3 sites × 2 weeks)
    subject_counts['is_complete'] = (
        (subject_counts['n_samples'] == 6) &
        (subject_counts['n_body_sites'] == 3) &
        (subject_counts['n_weeks'] == 2)
    )

    complete_subjects = subject_counts[subject_counts['is_complete']].index.tolist()

    print(f"Total subjects: {len(subject_counts)}")
    print(f"Complete subjects (6 samples, 3 sites, 2 weeks): {len(complete_subjects)}")
    print(f"  - Cincinnati: {subject_counts[(subject_counts['is_complete']) & (subject_counts['Location'] == 'Cincinnati')].shape[0]}")
    print(f"  - Hangzhou: {subject_counts[(subject_counts['is_complete']) & (subject_counts['Location'] == 'Hangzhou')].shape[0]}")

    # Distribution of sample counts
    print(f"\nDistribution of samples per subject:")
    print(subject_counts['n_samples'].value_counts().sort_index())

    # Save complete subjects list
    complete_subjects_df = pd.DataFrame({
        'PatientID': complete_subjects,
        'Location': subject_counts.loc[complete_subjects, 'Location']
    })
    output_path = RESULTS_DIR / "complete_subjects.tsv"
    complete_subjects_df.to_csv(output_path, sep="\t", index=False)
    print(f"\n✓ Complete subjects list saved to: {output_path}")

    # Save subject completeness summary
    subject_counts.to_csv(RESULTS_DIR / "subject_completeness.tsv", sep="\t")

    # Mark complete subjects in metadata
    metadata['subject_is_complete'] = metadata['PatientID'].isin(complete_subjects)

    return metadata, complete_subjects

def filter_low_abundance_genes(amr_data, min_prevalence=0.05, min_rpm=1.0):
    """Filter genes that are low abundance or rare"""
    print("\n" + "="*60)
    print("GENE FILTERING")
    print("="*60)

    n_samples = amr_data['sample_name'].nunique()

    # Count prevalence per gene
    gene_prevalence = amr_data.groupby('gene_name').agg({
        'sample_name': 'nunique',
        'rpm': ['mean', 'median', 'max']
    })
    gene_prevalence.columns = ['n_samples', 'mean_rpm', 'median_rpm', 'max_rpm']
    gene_prevalence['prevalence'] = gene_prevalence['n_samples'] / n_samples

    # Filter criteria
    genes_to_keep = gene_prevalence[
        (gene_prevalence['prevalence'] >= min_prevalence) &
        (gene_prevalence['max_rpm'] >= min_rpm)
    ].index.tolist()

    print(f"Total genes before filtering: {len(gene_prevalence)}")
    print(f"Genes passing filters (prevalence >= {min_prevalence}, max RPM >= {min_rpm}): {len(genes_to_keep)}")
    print(f"Genes removed: {len(gene_prevalence) - len(genes_to_keep)}")

    # Save gene prevalence stats
    gene_prevalence.to_csv(RESULTS_DIR / "gene_prevalence_stats.tsv", sep="\t")

    # Save filtered gene list
    pd.DataFrame({'gene_name': genes_to_keep}).to_csv(
        RESULTS_DIR / "genes_passing_qc.txt",
        sep="\t",
        index=False,
        header=False
    )

    print(f"✓ Gene prevalence stats saved to: {RESULTS_DIR / 'gene_prevalence_stats.tsv'}")
    print(f"✓ Filtered gene list saved to: {RESULTS_DIR / 'genes_passing_qc.txt'}")

    return genes_to_keep, gene_prevalence

def create_qc_figures(metadata, sample_summary):
    """Create QC visualization figures"""
    print("\n" + "="*60)
    print("CREATING QC FIGURES")
    print("="*60)

    data = metadata[metadata['has_amr_data']].copy()

    # Figure 1: Sequencing depth and AMR burden
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Total AMR reads by location and body site
    sns.boxplot(data=data, x='Location', y='total_amr_rpm', hue='SampleType', ax=axes[0, 0])
    axes[0, 0].set_title('Total AMR RPM by Location and Body Site')
    axes[0, 0].set_ylabel('Total AMR RPM')
    axes[0, 0].tick_params(axis='x', rotation=45)

    # Unique AMR genes
    sns.boxplot(data=data, x='Location', y='unique_amr_genes', hue='SampleType', ax=axes[0, 1])
    axes[0, 1].set_title('Unique AMR Genes by Location and Body Site')
    axes[0, 1].set_ylabel('Number of Unique AMR Genes')
    axes[0, 1].tick_params(axis='x', rotation=45)

    # Drug classes
    sns.boxplot(data=data, x='Location', y='num_drug_classes', hue='SampleType', ax=axes[0, 2])
    axes[0, 2].set_title('Drug Classes by Location and Body Site')
    axes[0, 2].set_ylabel('Number of Drug Classes')
    axes[0, 2].tick_params(axis='x', rotation=45)

    # Longitudinal changes (Week 1 vs Week 3)
    sns.boxplot(data=data, x='SampleCollectionWeek', y='total_amr_rpm', hue='Location', ax=axes[1, 0])
    axes[1, 0].set_title('AMR RPM: Week 1 vs Week 3')
    axes[1, 0].set_ylabel('Total AMR RPM')

    sns.boxplot(data=data, x='SampleCollectionWeek', y='unique_amr_genes', hue='Location', ax=axes[1, 1])
    axes[1, 1].set_title('Unique AMR Genes: Week 1 vs Week 3')
    axes[1, 1].set_ylabel('Number of Unique AMR Genes')

    # High confidence AMR genes
    sns.boxplot(data=data, x='Location', y='high_conf_amr', hue='SampleType', ax=axes[1, 2])
    axes[1, 2].set_title('High Confidence AMR Genes')
    axes[1, 2].set_ylabel('Number of High Confidence AMR Genes')
    axes[1, 2].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    output_path = FIGURES_DIR / "qc_amr_burden_summary.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✓ QC figure saved to: {output_path}")
    plt.close()

    return

def main():
    # Load data
    metadata, amr_data, sample_summary = load_data()

    # Analyze sequencing depth
    data_with_amr, low_samples = analyze_sequencing_depth(metadata, sample_summary)

    # Identify complete subjects
    metadata, complete_subjects = identify_complete_subjects(metadata)

    # Filter low abundance genes
    genes_passing_qc, gene_prevalence = filter_low_abundance_genes(amr_data)

    # Create QC figures
    create_qc_figures(metadata, sample_summary)

    # Save updated metadata with completeness flags
    output_path = RESULTS_DIR / "master_metadata_with_qc.tsv"
    metadata.to_csv(output_path, sep="\t", index=False)
    print(f"\n✓ Updated metadata with QC flags saved to: {output_path}")

    # Summary
    print("\n" + "="*60)
    print("QC SUMMARY")
    print("="*60)
    print(f"Total samples in metadata: {len(metadata)}")
    print(f"Samples with AMR data: {metadata['has_amr_data'].sum()}")
    print(f"Complete subjects: {len(complete_subjects)}")
    print(f"  - Total samples from complete subjects: {metadata[metadata['subject_is_complete']].shape[0]}")
    print(f"Genes passing QC filters: {len(genes_passing_qc)}")
    print(f"\nRecommendation: Use complete subjects for paired/longitudinal analysis")
    print(f"                All samples can be used for cross-sectional comparisons")

if __name__ == "__main__":
    main()
