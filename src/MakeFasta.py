""" function to normalize libraries to translate DNA sequences into amino
# acid sequences of peptides
#
# inputs:
#       mer: number of peptides per sequence
#
#       filenamefastq: string containing the filepath, filename, and
#       extension of the raw data file from deep sequencing
#
#       startflank: sequence prior to the beginning of the peptide sequence
#
#       endflank: sequence following the end of the peptide sequence
#
#       filenameoutput: string containing the filepath, filename, and
#       extension of the excel file to be created
#
#       PhD7: boolean identifier of whether the files correspond to screens
#       completed using NEB's PhD library
#
# output:
#       libraryexcel: table of sequences and frequencies to be exported to excel
#
# Created by Lindsey Brinton at the University of Virginia, 2015
"""

import os
import numpy as np
import collections
import pandas as pd
import re
import math
import time
import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from Bio import motifs
# must pip install openpyxl

def createFasta(PepFile, filenameoutput, max_seqs=500):
    df = pd.read_excel(PepFile, header=0).reset_index(drop=True)

    # Identify string columns
    string_cols = df.select_dtypes(include=["object", "string"]).columns.tolist()

    if not string_cols:
        raise ValueError("No string columns found in input table")

    def is_aa_sequence(seq):
        return bool(re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq))

    def is_dna_sequence(seq):
        return bool(re.fullmatch(r"[ACGT]+", seq))

    # Heuristic: column that looks like amino acid sequences
    aa_col = None
    for col in string_cols:
        sample = df[col].dropna().astype(str).head(10)
        if sample.apply(is_aa_sequence).all() and not sample.apply(is_dna_sequence).all():
            aa_col = col
            break

    # Fallback: just take the first string column
    if aa_col is None:
        aa_col = string_cols[0]

    # Rename detected column to AA
    df = df.rename(columns={aa_col: "AA"})

    print(f"Using column '{aa_col}' as AA")
    print(df.head())

    sequences = df["AA"].astype(str).tolist()

    with open(filenameoutput, "w") as ofile:
        for i, seq in enumerate(sequences[:max_seqs]):
            ofile.write(f">seq_{i}\n{seq}\n")

    ofile.close()

    return


print(snakemake.input)
print(snakemake.output)
for i in range(len(snakemake.input)):
    createFasta(snakemake.input[i],snakemake.output[i])
