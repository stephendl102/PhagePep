import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import math
from scipy import stats
#from normalizelibrary import *

# Main PHASTpep file code
# created by Lindsey Brinton at the University of Virginia, 2015
# updated by Lindsey Brinton at the University of Virginia, 2016
# Translated to Python by Stephen Lees, 2022

# ---- Constants ----
EXPECTED_COLS = ['Seq', 'Count', 'Peptide', 'Percentage', 'Error']

def standardize_df(df):
    """Standardize first 5 columns and return a clean copy."""
    if df.shape[1] < 5:
        raise ValueError("Expected at least 5 columns in input files")
    df = df.iloc[:, :5].copy()
    df.columns = EXPECTED_COLS
    return df

def EnrichmentRatio(input_xlsx, output_file):
    """
    Ask user which files are Round 1 and Round 2 (required),
    optionally Round 3,
    compute enrichment ratios,
    and write result to a new Excel file.
    """

    # ---- Display file options ----
    print("Available files:")
    for i, f in enumerate(input_xlsx):
        print(f"{i}: {f}")

    # ---- User selection ----
    r1_idx = int(input("Enter the index of ROUND 1 file: "))
    r2_idx = int(input("Enter the index of ROUND 2 file: "))

    r3_input = input("Enter the index of ROUND 3 file (or press Enter to skip): ")
    r3_idx = int(r3_input) if r3_input.strip() != "" else None

    round1_file = input_xlsx[r1_idx]
    round2_file = input_xlsx[r2_idx]
    round3_file = input_xlsx[r3_idx] if r3_idx is not None else None

    # ---- Read and standardize ----
    df1 = standardize_df(pd.read_excel(round1_file, header=None)).rename(columns={'Count': 'Count_R1'})
    df2 = standardize_df(pd.read_excel(round2_file, header=None)).rename(columns={'Count': 'Count_R2'})
    df3 = None
    if round3_file is not None:
        df3 = standardize_df(pd.read_excel(round3_file, header=None)).rename(columns={'Count': 'Count_R3'})

    # ---- Merge on Peptide ----
    merged = pd.merge(
        df1[['Peptide', 'Count_R1']],
        df2[['Peptide', 'Count_R2']],
        on='Peptide',
        how='outer'
    )

    if df3 is not None:
        merged = pd.merge(
            merged,
            df3[['Peptide', 'Count_R3']],
            on='Peptide',
            how='outer'
        )

    # ---- Ensure numeric counts ----
    for col in merged.columns:
        if col.startswith("Count_"):
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    # ---- Normalize by library size ----
    total_r1 = merged['Count_R1'].sum(skipna=True)
    total_r2 = merged['Count_R2'].sum(skipna=True)

    merged['Norm_R1'] = merged['Count_R1'] / total_r1
    merged['Norm_R2'] = merged['Count_R2'] / total_r2

    if df3 is not None:
        total_r3 = merged['Count_R3'].sum(skipna=True)
        merged['Norm_R3'] = merged['Count_R3'] / total_r3

    # ---- Replace missing values with column minima ----
    norm_cols = ['Norm_R1', 'Norm_R2'] + (['Norm_R3'] if df3 is not None else [])
    for col in norm_cols:
        col_min = merged[col].min(skipna=True)
        merged[col] = merged[col].fillna(col_min)

    # ---- Avoid division by zero ----
    merged['Norm_R1'] = merged['Norm_R1'].replace(0, np.nan)
    merged['Norm_R2'] = merged['Norm_R2'].replace(0, np.nan)
    if df3 is not None:
        merged['Norm_R3'] = merged['Norm_R3'].replace(0, np.nan)

    # ---- Compute enrichment ratios ----
    merged['ER'] = np.log2(merged['Norm_R2'] / merged['Norm_R1'])

    if df3 is not None:
        merged['ER'] = (
            np.log2(merged['Norm_R2'] / merged['Norm_R1']) +
            np.log2(merged['Norm_R3'] / merged['Norm_R2'])
        )

    merged = merged.replace([np.inf, -np.inf], np.nan)

    # ---- Write output ----
    merged.to_excel(output_file[0], index=False)
    print(f"\nRatio file written to: {output_file[0]}")

    return



print(snakemake.input)
print(snakemake.output)
#Diversity(snakemake.output[0], snakemake.input)
EnrichmentRatio(snakemake.input,snakemake.output)
