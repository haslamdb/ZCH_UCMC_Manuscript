#!/usr/bin/env python3
"""
One-shot canonicalization of ZCH sample names across the v2 pipeline.

The historical bug: ZCH single-digit subjects (N1..N9) got zero-padded to
N01..N09 in the sample manifest and propagated everywhere downstream.
The raw FASTQs in /bulkpool/sequence_data/mss_data/ use ZJH_N1_*; the
clinical naming follows that convention.

This script standardizes the convention back to ZJH_N{stripped}_* across
flat data files. (Renaming radar_runs/ dirs is handled separately by
canonicalize_radar_dirs.sh — that's a filesystem op.)

Affected files:
  - data/nicu_amr_only.tsv             (sample_name column)
  - data/nicu_amr_only_permissive.tsv  (sample_name column)
  - data/sample_summary.tsv            (sample_name column)
  - data/sample_summary_permissive.tsv (sample_name column)
  - data/sample_manifest.csv           (sample_id column)
  - <project>/metadata/AllNICUSampleKey...csv  (PatientID, SampleID, Sample for ZCH rows)

The fix is idempotent: it ignores already-canonical names.
"""

import re
import shutil
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
V2 = PROJECT_ROOT / "nicu_resistome_analysis_v2"
DATA = V2 / "data"
METADATA_DIR = PROJECT_ROOT / "metadata"

SAMPLE_PATTERN = re.compile(r'^ZJH_N0(\d)(_|$)')
SUBJECT_PATTERN = re.compile(r'^N0(\d)$')

def canon_sample(s):
    """ZJH_N0X_a_b -> ZJH_NX_a_b; ZJH_N0X -> ZJH_NX; otherwise unchanged."""
    if not isinstance(s, str):
        return s
    return SAMPLE_PATTERN.sub(r'ZJH_N\1\2', s)

def canon_subject(s):
    """N0X -> NX (single-digit only). Only applied to ZCH rows."""
    if not isinstance(s, str):
        return s
    return SUBJECT_PATTERN.sub(r'N\1', s)

def backup_once(path: Path) -> Path:
    """Copy file to <path>.pre_canonicalize if not already done."""
    bak = path.with_suffix(path.suffix + '.pre_canonicalize')
    if not bak.exists():
        shutil.copy2(path, bak)
        print(f"  backup -> {bak.name}")
    else:
        print(f"  backup already present: {bak.name}")
    return bak


def fix_sample_name_column(path: Path, col: str = 'sample_name', sep: str = '\t'):
    print(f"\n{path.relative_to(PROJECT_ROOT)}:")
    if not path.exists():
        print("  (missing — skip)")
        return
    backup_once(path)
    df = pd.read_csv(path, sep=sep, low_memory=False)
    before = df[col].copy()
    df[col] = df[col].apply(canon_sample)
    n_changed = (before != df[col]).sum()
    df.to_csv(path, sep=sep, index=False)
    print(f"  rewrote {col}: {n_changed:,} cells changed (of {len(df):,})")


def fix_sample_key():
    """Fix the source sample key CSV — strip leading zero from ZCH single-digit
    subjects in PatientID, SampleID, and Sample. The file used by 01 is the
    20250320 revised key."""
    path = METADATA_DIR / "AllNICUSampleKeyRevised20250320_for_HanzhouCincinnatiSamples.csv"
    if not path.exists():
        print(f"\n[sample_key] {path} not found — skip")
        return
    print(f"\n{path.relative_to(PROJECT_ROOT)}:")
    backup_once(path)
    sk = pd.read_csv(path)
    is_zch = sk['Location'] == 'Hangzhou'

    n_changes = {}
    # SampleID and Sample contain ZJH_N0X_*
    for col in ['SampleID', 'Sample']:
        before = sk[col].copy()
        sk.loc[is_zch, col] = sk.loc[is_zch, col].apply(canon_sample)
        # Sample column has unprefixed names like N01_1_2 — also strip the
        # leading zero AND prepend ZJH_ so it becomes the canonical form.
        if col == 'Sample':
            zch_unprefixed = is_zch & sk[col].str.match(r'^N0\d_', na=False)
            sk.loc[zch_unprefixed, col] = sk.loc[zch_unprefixed, col].str.replace(
                r'^N0(\d)_', r'ZJH_N\1_', regex=True)
            zch_prefix_missing = is_zch & sk[col].str.match(r'^N\d', na=False) & ~sk[col].str.startswith('ZJH_')
            sk.loc[zch_prefix_missing, col] = 'ZJH_' + sk.loc[zch_prefix_missing, col]
        n_changes[col] = (before != sk[col]).sum()

    # PatientID for ZCH single-digit
    before = sk['PatientID'].copy()
    sk.loc[is_zch, 'PatientID'] = sk.loc[is_zch, 'PatientID'].apply(canon_subject)
    n_changes['PatientID'] = (before != sk['PatientID']).sum()

    sk.to_csv(path, index=False)
    print(f"  changes: {n_changes}")


def main():
    print("=" * 70)
    print("CANONICALIZE ZJH SAMPLE NAMES")
    print("=" * 70)

    fix_sample_name_column(DATA / "nicu_amr_only.tsv", col='sample_name', sep='\t')
    fix_sample_name_column(DATA / "nicu_amr_only_permissive.tsv", col='sample_name', sep='\t')
    fix_sample_name_column(DATA / "sample_summary.tsv", col='sample_name', sep='\t')
    fix_sample_name_column(DATA / "sample_summary_permissive.tsv", col='sample_name', sep='\t')
    fix_sample_name_column(DATA / "sample_manifest.csv", col='sample_id', sep=',')
    fix_sample_key()

    print("\nDone.")


if __name__ == '__main__':
    main()
