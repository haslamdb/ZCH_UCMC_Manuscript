#!/usr/bin/env python3
"""
Rename antibiotic columns in parsed CSV files to use antibiotic names instead of indices.
"""

import pandas as pd

# Antibiotic mapping based on the order provided
antibiotic_names = [
    "Ampicillin",              # 1
    "Penicillin",              # 2
    "Nafcillin",               # 3
    "Cefotaxime",              # 4
    "Moxalactam",              # 5
    "Ceftazidime",             # 6
    "Ceftriaxone",             # 7
    "Cefepime",                # 8
    "Gentamicin",              # 9
    "Tobramycin",              # 10
    "Vancomycin",              # 11
    "Fluconazole",             # 12
    "Amphotericin_B",          # 13
    "Meropenem",               # 14
    "Azithromycin",            # 15
    "Piperacillin_Tazobactam", # 16
    "None",                    # 17
    "Other"                    # 18 (only in week 3)
]

def create_column_mapping(columns, site_suffix):
    """
    Create mapping from old column names to new ones.

    Pattern:
    - antibiotics_duration1_v2_zju -> Ampicillin_w1
    - antibiotics_duration1_v3_zju -> Ampicillin_w3
    - antibiotics_duration1_v2 -> Ampicillin_w1 (UC format)
    """
    mapping = {}

    for col in columns:
        # Skip non-antibiotic columns
        if 'antibiotics_duration' not in col:
            continue

        # Extract the antibiotic index and week
        # Format can be: antibiotics_duration1_v2_zju or antibiotics_duration1_v2
        parts = col.split('_')

        # Find the duration part with the number
        duration_part = None
        week_part = None

        for part in parts:
            if part.startswith('duration'):
                # Extract number from 'duration1', 'duration10', etc.
                num_str = part.replace('duration', '')
                if num_str.isdigit():
                    duration_part = int(num_str)
            elif part.startswith('v'):
                # Extract week from 'v2' or 'v3'
                week_num = part.replace('v', '')
                if week_num.isdigit():
                    week_part = int(week_num)

        # Special handling for columns that end with _zju or just a number
        # antibiotics_duration18_zju (no v2/v3) - this is likely "Other" for week 3
        if duration_part and not week_part:
            # Check if this is the special "Other" column (index 18)
            if duration_part == 18:
                new_name = "Other_w3"
                mapping[col] = new_name
                continue

        if duration_part and week_part:
            # Map index to antibiotic name (1-indexed to 0-indexed)
            if 1 <= duration_part <= len(antibiotic_names):
                antibiotic = antibiotic_names[duration_part - 1]

                # Map v2 -> w1, v3 -> w3
                week = "w1" if week_part == 2 else f"w{week_part - 1}"

                new_name = f"{antibiotic}_{week}"
                mapping[col] = new_name

    return mapping

def rename_columns_in_file(input_file, output_file, site_suffix):
    """Rename columns in a parsed antibiotic file."""

    print(f"\nProcessing {input_file}...")

    # Read the file
    df = pd.read_csv(input_file)

    # Create column mapping
    mapping = create_column_mapping(df.columns, site_suffix)

    print(f"Found {len(mapping)} antibiotic columns to rename")

    # Show a few examples
    print("\nExample mappings:")
    for old, new in list(mapping.items())[:5]:
        print(f"  {old} → {new}")

    # Rename columns
    df = df.rename(columns=mapping)

    # Save
    df.to_csv(output_file, index=False)

    print(f"✓ Saved to {output_file}")

    return df

def main():
    """Main function to rename columns in both files."""

    print("="*70)
    print("RENAMING ANTIBIOTIC COLUMNS")
    print("="*70)

    # Process ZCH file
    zch_df = rename_columns_in_file(
        'ZCH_antibiotics_parsed.csv',
        'ZCH_antibiotics_parsed_renamed.csv',
        'zju'
    )

    # Process UC file
    uc_df = rename_columns_in_file(
        'UC_antibiotics_parsed.csv',
        'UC_antibiotics_parsed_renamed.csv',
        ''
    )

    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)
    print("\nGenerated files:")
    print("  - ZCH_antibiotics_parsed_renamed.csv")
    print("  - UC_antibiotics_parsed_renamed.csv")

    # Show the new column structure
    print("\nNew column names (first 10 antibiotic columns):")
    antibiotic_cols = [col for col in zch_df.columns if any(ab in col for ab in antibiotic_names)][:10]
    for col in antibiotic_cols:
        print(f"  - {col}")

if __name__ == '__main__':
    main()
