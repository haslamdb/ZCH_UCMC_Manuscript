#!/usr/bin/env python3
"""
Parse antibiotic duration data from ZCH and UC CSV files.
Converts cryptic notation into number of days antibiotic was received.
"""

import pandas as pd
import re
from datetime import datetime
from typing import Union
from math import ceil

def parse_date_range(date_str: str) -> int:
    """
    Parse date range in format M.D-M.D and return number of days (inclusive).
    Examples: '6.3-6.8' -> 6 days, '5.28-6.2' -> 6 days (May 28 to June 2)
    """
    # Handle date ranges like "6.3-6.8"
    match = re.match(r'(\d+)\.(\d+)-(\d+)\.(\d+)', date_str)
    if match:
        start_month, start_day, end_month, end_day = map(int, match.groups())

        # Validate dates
        try:
            # Create datetime objects (year doesn't matter for day counting)
            start_date = datetime(2023, start_month, start_day)
            end_date = datetime(2023, end_month, end_day)

            # If end date is before start date, it likely crossed into next year
            if end_date < start_date:
                end_date = datetime(2024, end_month, end_day)

            # Calculate days (inclusive)
            days = (end_date - start_date).days + 1
            return days
        except ValueError as e:
            print(f"Warning: Invalid date in range '{date_str}': {e}")
            # For invalid dates, try to estimate based on day numbers
            # This is a fallback for entries like "2.29-3.4" where Feb 29 doesn't exist
            if start_month == end_month:
                # Same month, just subtract days
                return abs(end_day - start_day) + 1
            else:
                # Different months, rough estimate
                # Count remaining days in start month + days in end month
                days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
                if start_month <= 12 and end_month <= 12:
                    days_left_in_start = days_in_month[start_month - 1] - start_day + 1
                    return days_left_in_start + end_day
                else:
                    print(f"Warning: Invalid month in '{date_str}'")
                    return 0

    return 0

def parse_zch_entry(entry: str) -> int:
    """
    Parse a single ZCH antibiotic entry and return total days (rounded up).

    Patterns:
    - Date ranges: '6.3-6.8' -> count days
    - Single dose: '1 dose (5.27)' -> 1 day
    - With dosing info: '40mg/kg/dose (5.28-5.31)' -> extract date range
    - Every N days: '1dose/3days (7.12-7.21)' or '/3d' -> calculate actual dosing days

    All fractional days are rounded UP to the nearest integer.
    """
    # Handle NaN and empty values
    if pd.isna(entry):
        return 0

    # Convert to string and check if empty
    entry = str(entry).strip()
    if entry == '' or entry == 'nan':
        return 0

    # Pattern 1: Special case like "3.27-3.30=2.5" (explicitly shows days)
    equals_match = re.search(r'=\s*([\d.]+)', entry)
    if equals_match:
        return ceil(float(equals_match.group(1)))

    # Pattern 2: Contains "/3days" or "/3d" pattern with date range
    # Example: "1dose/3days (7.12-7.21)" or "10.24-10.31 (/3days)"
    freq_match = re.search(r'/(\d+)\s*d(?:ays?)?', entry, re.IGNORECASE)
    if freq_match:
        date_match = re.search(r'(\d+\.\d+-\d+\.\d+)', entry)
        if date_match:
            span_days = parse_date_range(date_match.group(1))
            interval = int(freq_match.group(1))
            # Calculate number of actual dosing days
            # If 1 dose every 3 days for 10 day span = ceil(10/3) = 4 doses = 4 days
            actual_doses = (span_days + interval - 1) // interval
            return actual_doses

    # Pattern 3: Date range with possible dosing info like "40mg/kg/dose (5.28-5.31)" or "6.3-6.8"
    # Extract any date range BEFORE checking for dose keywords
    date_matches = re.findall(r'(\d+\.\d+-\d+\.\d+)', entry)
    if date_matches:
        # If multiple date ranges, sum them
        total_days = 0
        for date_range in date_matches:
            total_days += parse_date_range(date_range)
        return total_days

    # Pattern 4: Single dose like "1 dose (5.27)" or "4mg* 1dose" or "4mg 6 dose"
    if 'dose' in entry.lower():
        # Check if it's "X dose" pattern (where X is a number)
        dose_match = re.search(r'(\d+)\s*dose', entry, re.IGNORECASE)
        if dose_match:
            num_doses = int(dose_match.group(1))
            return num_doses

        # Check for just "dose" (without number) - single dose
        if re.search(r'(?<!\d)\s*dose', entry, re.IGNORECASE):
            return 1

    # Pattern 5: Just a single date like "6.18" - treat as 1 day
    single_date_match = re.match(r'^(\d+\.\d+)$', entry)
    if single_date_match:
        return 1

    # Pattern 6: Number at end like "4.7 0.5" - use the trailing number (rounded up)
    number_match = re.search(r'\s+([\d.]+)$', entry)
    if number_match:
        return ceil(float(number_match.group(1)))

    # If we can't parse it, return 0 and we'll review
    print(f"Warning: Could not parse ZCH entry: '{entry}'")
    return 0

def parse_uc_entry(entry: str) -> int:
    """
    Parse a single UC antibiotic entry and return total days (rounded up).

    Patterns:
    - '100mg/kg q12h*2d' -> 2 days
    - '25mg/kg*bid*3d' -> 3 days
    - '25mg/kg q6h*3d+25mg/kg tid*2d' -> 3 + 2 = 5 days
    - Complex with multiple segments separated by '+'

    All fractional days are rounded UP to the nearest integer.
    """
    # Handle NaN and empty values
    if pd.isna(entry):
        return 0

    # Convert to string and check if empty
    entry = str(entry).strip()
    if entry == '' or entry == 'nan':
        return 0

    # Skip text-only entries like "nystatin bid topical"
    if not re.search(r'\*\d+', entry):
        # Check if it might still have dosing info we should count
        # For now, if no *Xd pattern, return 0
        return 0

    total_days = 0.0

    # Split by '+' to handle complex entries like "X*3d+Y*2d"
    segments = entry.split('+')

    for segment in segments:
        segment = segment.strip()

        # Look for *Xd pattern (X days)
        # Can be *2d or *4.5d for partial days
        day_matches = re.findall(r'\*\s*([\d.]+)\s*d(?:ay)?', segment, re.IGNORECASE)

        if day_matches:
            # Sum all day patterns found in this segment (usually just one)
            for day_str in day_matches:
                total_days += float(day_str)

    # Round up to nearest integer
    return ceil(total_days) if total_days > 0 else 0

def parse_antibiotic_file(input_file: str, output_file: str, file_type: str = 'ZCH'):
    """
    Parse antibiotic file and create output with day counts.

    Args:
        input_file: Path to input CSV
        output_file: Path to output CSV
        file_type: 'ZCH' or 'UC' to determine parsing method
    """
    print(f"\n{'='*70}")
    print(f"Parsing {file_type} antibiotic data")
    print(f"{'='*70}")

    # Read the CSV
    df = pd.read_csv(input_file)

    print(f"Loaded {len(df)} records")
    print(f"Columns: {len(df.columns)} antibiotic columns")

    # Create output dataframe with same structure
    output_df = df.copy()

    # Identify antibiotic columns (exclude record_id and mrn columns)
    antibiotic_cols = [col for col in df.columns if 'antibiotics_duration' in col]

    print(f"Found {len(antibiotic_cols)} antibiotic duration columns")

    # Parse each antibiotic column
    parser_func = parse_zch_entry if file_type == 'ZCH' else parse_uc_entry

    for col in antibiotic_cols:
        print(f"\nProcessing {col}...")
        output_df[col] = df[col].apply(parser_func)

        # Show some statistics
        non_zero = output_df[col] > 0
        if non_zero.any():
            print(f"  - {non_zero.sum()} patients received this antibiotic")
            print(f"  - Mean days (when given): {output_df.loc[non_zero, col].mean():.2f}")
            print(f"  - Total days across all patients: {output_df[col].sum():.2f}")

    # Save output
    output_df.to_csv(output_file, index=False)
    print(f"\n✓ Saved parsed data to {output_file}")

    # Summary statistics
    print(f"\n{'='*70}")
    print(f"SUMMARY for {file_type}")
    print(f"{'='*70}")
    total_antibiotic_days = output_df[antibiotic_cols].sum().sum()
    print(f"Total antibiotic-days across all patients: {total_antibiotic_days:.2f}")

    # Per patient statistics
    patient_total_days = output_df[antibiotic_cols].sum(axis=1)
    patients_with_abx = (patient_total_days > 0).sum()
    print(f"Patients who received antibiotics: {patients_with_abx}/{len(df)}")
    if patients_with_abx > 0:
        print(f"Mean antibiotic days per patient (among those who received): {patient_total_days[patient_total_days > 0].mean():.2f}")

    return output_df

def main():
    """Main function to parse both antibiotic files."""

    # Parse ZCH data
    zch_output = parse_antibiotic_file(
        'ZCH_antibiotics.csv',
        'ZCH_antibiotics_parsed.csv',
        'ZCH'
    )

    # Parse UC data
    uc_output = parse_antibiotic_file(
        'UC_antibiotics.csv',
        'UC_antibiotics_parsed.csv',
        'UC'
    )

    print("\n" + "="*70)
    print("PARSING COMPLETE")
    print("="*70)
    print("\nGenerated files:")
    print("  - ZCH_antibiotics_parsed.csv")
    print("  - UC_antibiotics_parsed.csv")
    print("\nPlease review the output files and check for any warnings above.")

if __name__ == '__main__':
    main()
