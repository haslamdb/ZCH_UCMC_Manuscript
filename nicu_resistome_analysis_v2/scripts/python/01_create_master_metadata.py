#!/usr/bin/env python3
"""
Create master metadata file by merging:
1. Sample key with clinical metadata
2. Sample summary with AMR/stress statistics
3. Kraken2 bacterial read counts

Output: master_metadata.tsv
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Define paths
PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
DATA_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "data"
METADATA_DIR = PROJECT_ROOT / "metadata"
RESULTS_DIR = PROJECT_ROOT / "nicu_resistome_analysis_v2" / "results" / "qc"

# Create results directory if it doesn't exist
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def main():
    print("Loading data files...")

    # 1. Load sample key (most recent version)
    sample_key = pd.read_csv(
        METADATA_DIR / "AllNICUSampleKeyRevised20250320_for_HanzhouCincinnatiSamples.csv"
    )
    print(f"Sample key: {sample_key.shape[0]} rows, {sample_key.shape[1]} columns")

    # 2. Load sample summary with AMR/stress statistics
    sample_summary = pd.read_csv(
        DATA_DIR / "sample_summary.tsv",
        sep="\t"
    )
    print(f"Sample summary: {sample_summary.shape[0]} rows, {sample_summary.shape[1]} columns")

    # 3. Load antibiotic exposure data
    antibiotics = pd.read_csv(
        PROJECT_ROOT / "antibiotics_combined.csv"
    )
    print(f"Antibiotic exposure: {antibiotics.shape[0]} subjects, {antibiotics.shape[1]} columns")

    # Subject ID convention (post-canonicalize_zjh_names.py): UCMC uses N01..N68
    # (always two-digit) and ZCH uses N1..N9 + N10..N59 (no leading zero for
    # single-digit, matching the FASTQ filename convention). The (Subject,
    # Location) merge key disambiguates the two.

    # Standardize location names to match antibiotics file and rest of analysis
    sample_key['Location'] = sample_key['Location'].map({'Cincinnati': 'UCMC', 'Hangzhou': 'ZCH'})

    # EXCLUDE nose swabs (labeled as "Nonese") - these were not sequenced
    print(f"Samples before filtering: {len(sample_key)}")
    sample_key = sample_key[sample_key['SampleType'] != 'Nonese']
    print(f"Samples after excluding nose swabs (Nonese): {len(sample_key)}")

    # Extract unique SubjectID from SampleID (which includes ZJH_ prefix for Hangzhou)
    # UCMC samples: N01_1_2 -> N01
    # ZCH samples: ZJH_N01_1_2 -> ZJH_N01
    sample_key['SubjectID'] = sample_key['SampleID'].str.extract(r'^((?:ZJH_)?N\d+)')[0]

    # Merge datasets
    print("\nMerging datasets...")

    # Merge 1: sample_key + sample_summary
    # Match on Sample column (the merged data uses sample names without ZJH_ prefix)
    master = sample_key.merge(
        sample_summary,
        left_on="Sample",
        right_on="sample_name",
        how="left"
    )
    print(f"After merging sample_key + sample_summary: {master.shape[0]} rows")

    # Merge 2: + antibiotic exposure data
    # Match on PatientID AND Location (PatientID is reused across sites!)

    master = master.merge(
        antibiotics,
        left_on=["PatientID", "Location"],
        right_on=["Subject", "Location"],
        how="left",
        suffixes=("", "_abx")
    )
    print(f"After merging + antibiotic exposure: {master.shape[0]} rows")
    print("Note: Using RPM normalization (already calculated in sample_summary)")

    # Quality flags
    master['has_metadata'] = ~master['PatientID'].isna()
    master['has_amr_data'] = ~master['sample_name'].isna()
    master['has_abx_data'] = ~master['Subject'].isna()

    # Calculate total antibiotic exposure days for week 1 and week 3
    # Note: In antibiotics file, _w1 = Week.1 and _w2 = Week.3
    # Get all antibiotic columns (exclude Subject, Location, None, Other)
    abx_w1_cols = [col for col in master.columns if col.endswith('_w1') and col != 'None_w1']
    abx_w3_cols = [col for col in master.columns if col.endswith('_w2') and col != 'None_w2']  # _w2 means Week.3!

    master['total_abx_days_week1'] = master[abx_w1_cols].sum(axis=1)
    master['total_abx_days_week3'] = master[abx_w3_cols].sum(axis=1)
    master['any_abx_week1'] = (master['total_abx_days_week1'] > 0).astype(int)
    master['any_abx_week3'] = (master['total_abx_days_week3'] > 0).astype(int)

    # Assign appropriate week's antibiotic exposure to each sample
    master['total_abx_days'] = master.apply(
        lambda row: row['total_abx_days_week1'] if row['SampleCollectionWeek'] == 'Week.1'
                    else row['total_abx_days_week3'] if row['SampleCollectionWeek'] == 'Week.3'
                    else np.nan,
        axis=1
    )
    master['any_abx'] = (master['total_abx_days'] > 0).astype(int)

    # Calculate completeness by subject (only count valid samples with AMR data)
    print("\nCalculating sample completeness per subject...")
    valid_samples = master[
        (master['has_amr_data']) &
        (master['SampleType'].isin(['Axilla', 'Groin', 'Stool']))
    ]

    completeness = valid_samples.groupby('SubjectID').agg({
        'Sample': 'count',
        'SampleType': 'nunique',
        'SampleCollectionWeek': 'nunique'
    }).rename(columns={'Sample': 'n_valid_samples', 'SampleType': 'n_body_sites', 'SampleCollectionWeek': 'n_weeks'})

    master = master.merge(completeness, on='SubjectID', how='left')

    # Complete = 6 valid samples, 3 body sites, 2 weeks
    master['subject_complete'] = (
        (master['n_valid_samples'] == 6) &
        (master['n_body_sites'] == 3) &
        (master['n_weeks'] == 2)
    )
    master['subject_complete'] = master['subject_complete'].fillna(False)

    # Summary statistics
    print("\n" + "="*60)
    print("MASTER METADATA SUMMARY")
    print("="*60)
    print(f"Total samples: {len(master)}")
    print(f"Unique subjects (total): {master['SubjectID'].nunique()}")
    print(f"  - UCMC: {master[master['Location']=='UCMC']['SubjectID'].nunique()}")
    print(f"  - ZCH: {master[master['Location']=='ZCH']['SubjectID'].nunique()}")
    print(f"\nBy Location:")
    print(master.groupby('Location')['Sample'].count())
    print(f"\nBy Body Site:")
    print(master.groupby('SampleType')['Sample'].count())
    print(f"\nBy Week:")
    print(master.groupby('SampleCollectionWeek')['Sample'].count())
    print(f"\nSamples with AMR data: {master['has_amr_data'].sum()}")
    print(f"Subjects with antibiotic data: {master[master['has_abx_data']]['SubjectID'].nunique()}")
    print(f"  - UCMC: {master[(master['has_abx_data']) & (master['Location']=='UCMC')]['SubjectID'].nunique()}")
    print(f"  - ZCH: {master[(master['has_abx_data']) & (master['Location']=='ZCH')]['SubjectID'].nunique()}")
    print(f"\nAntibiotic exposure (subjects with abx data):")
    abx_subjects = master[master['has_abx_data']].drop_duplicates('SubjectID')
    print(f"  - Any antibiotics at Week 1: {abx_subjects['any_abx_week1'].sum()}/{len(abx_subjects)}")
    print(f"  - Any antibiotics at Week 3: {abx_subjects['any_abx_week3'].sum()}/{len(abx_subjects)}")
    print(f"  - Mean days at Week 1: {abx_subjects['total_abx_days_week1'].mean():.2f}")
    print(f"  - Mean days at Week 3: {abx_subjects['total_abx_days_week3'].mean():.2f}")
    print(f"\nComplete subjects (6 samples): {master[master['subject_complete']]['SubjectID'].nunique()}")
    print(f"  - UCMC: {master[(master['subject_complete']) & (master['Location']=='UCMC')]['SubjectID'].nunique()}")
    print(f"  - ZCH: {master[(master['subject_complete']) & (master['Location']=='ZCH')]['SubjectID'].nunique()}")
    print(f"Incomplete subjects: {master[~master['subject_complete']]['SubjectID'].nunique()}")

    # Save master metadata
    output_path = RESULTS_DIR / "master_metadata.tsv"
    master.to_csv(output_path, sep="\t", index=False)
    print(f"\n✓ Master metadata saved to: {output_path}")

    # Save QC report
    qc_report = master[['SampleID', 'Sample', 'SubjectID', 'PatientID', 'Location', 'SampleType',
                         'SampleCollectionWeek', 'PostNatalAntibiotics',
                         'has_metadata', 'has_amr_data', 'has_abx_data',
                         'n_valid_samples', 'n_body_sites', 'n_weeks', 'subject_complete']]
    qc_path = RESULTS_DIR / "sample_qc_report.tsv"
    qc_report.to_csv(qc_path, sep="\t", index=False)
    print(f"✓ QC report saved to: {qc_path}")

    # Identify problematic samples
    incomplete_subjects = master[~master['subject_complete']]['SubjectID'].unique()
    if len(incomplete_subjects) > 0:
        print(f"\n⚠ WARNING: {len(incomplete_subjects)} subjects have incomplete data:")
        for subj in incomplete_subjects[:10]:  # Show first 10
            n_samples = master[master['SubjectID'] == subj]['Sample'].count()
            print(f"  - {subj}: {n_samples} samples")

    return master

if __name__ == "__main__":
    master = main()
