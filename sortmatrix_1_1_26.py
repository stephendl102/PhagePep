""" function to combine multiple libraries together such that their
# duplicated sequences are combined into a single row
#
# inputs:
#       matrix: table of all libraries
#
#       libnumber: total number of libraries, not including reference
#       libraries
#
#       posnumber: number of libraries corresponding to positive screens
#
#       matrixnum: number of sequences to be included in final matrix
#
# outputs:
#       librarymatrix: cell array of sorted sequences
#
# Created by Lindsey Brinton at the University of Virginia, 2015
"""
import numpy as np
import pandas as pd
from normalizelibrary import *

CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'Q', #TAG = Q in amber strain
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    'CCN':'P', 'GGN':'G', 'ACN':'T', 'CGN':'R',
    'GCN':'A', 'CTN':'L', 'TCN':'S', 'GTN':'V',
}

def translate(seq):
    return ''.join(
        CODON_TABLE.get(seq[i:i+3], 'X')
        for i in range(0, len(seq), 3)
        if len(seq) % 3 == 0
    )


# --------------------------------------------------
# Main sorting / scoring function
# --------------------------------------------------
def sortmatrix(
    matrix,
    posnumber,  # number of positive libraries (normally 1)
    matrixnum,
    ReadDepths, # array of total number of reads for each screen
    Error=None,
    average_negative=False
):

    """
    Sorts a sequence-count matrix by enrichment score.
    """

    # ---- Basic setup ----
    seqs = matrix['Seq']

    counts = matrix.drop(columns='Seq')

    libnumber = counts.shape[1]
    negnumber = libnumber - posnumber # number of non-targeted cell screens


    print('Libnumber',libnumber)
    print('negative_number: ',negnumber)
    print('sorting the matrix')
    pd.set_option("display.precision", 10)

    # ---- Normalize by read depth ----
    depth_ratio = ReadDepths / ReadDepths[0] # array of depths normalized by
    print('Depth Ratio', depth_ratio)

    counts = counts.div(depth_ratio, axis=1) #normalize all reads by depth

    min_vals = counts.min(skipna=True)
    print('raw mins:', min_vals)
    counts = counts.fillna(min_vals)
    print(Error)
      # ---- Bleedthrough correction ----
    if Error is not None and Error < 0.10:
        pos_col = counts.columns[0]
        other_cols = counts.columns[1:]
        counts[other_cols] = counts[other_cols].sub(Error * counts[pos_col],axis=0)
    else:
        print(f"WARNING: Error rate {Error:.3f} ≥ 10%. Bleedthrough correction skipped.")
    print(counts)

    # ---- Compute column-wise minimums (after bleedthrough) ----
    # ---- Fill NaNs with column minimums ----

    counts = counts.fillna(min_vals)
    print('column-wise minimums (post-bleedthrough):', min_vals)
    counts = counts.clip(lower=min_vals, axis=1)

    # ---- Compute enrichment score ----
    avg_pos = counts.iloc[:, :posnumber].mean(axis=1)

    if negnumber > 0:
        neg_block = counts.iloc[:, posnumber:]
        avg_neg = (
            neg_block.mean(axis=1)
            if average_negative
            else neg_block.max(axis=1)
        )
        score = avg_pos / (avg_neg + avg_pos)
    else:
        score = avg_pos / (1 + avg_pos)

 # ---- Sort and truncate ----
    matrixnum = min(matrixnum, len(score))
    order = np.argsort(-score.values)[:matrixnum]

    df = pd.concat(
        [
            seqs.iloc[order].reset_index(drop=True),
            counts.iloc[order].reset_index(drop=True),
            pd.Series(score.iloc[order].values, name='Score')
        ],
        axis=1
    )

    print('Converting to amino acid sequences')
 # ---- Translate to amino acids ----
    df['Peptide'] = (
        df['Seq']
        .apply(translate)
        .str.replace('_', 'X')
    )

    # Format to work with BiLSTM
    peptide = df.pop('Peptide')
    score = df.pop('Score')
    df.insert(0, 'Peptide', peptide)
    df.drop(columns='Seq', inplace=True)

    df.columns = [
        c if c == 'Peptide'
        else pd.Series(c)
             .str.replace(r'^\d{8}', '', regex=True)
             .str.replace(r'_.*$', '', regex=True)
             .iloc[0]
        for c in df.columns
    ]

    score_cols = [c for c in df.columns if c != 'Peptide']

    if len(score_cols) != 2:
        raise ValueError("Expected exactly two non-Peptide columns")

    c1, c2 = score_cols

    # Assign E-scores
    df[f'{c1} E-score'] = score
    df[f'{c2} E-score'] = 1 - score


    # Duplicate original columns (3 copies each) with numbered suffixes and increment values
    c1_dups = pd.concat([df[c1] + i for i in range(3)], axis=1)
    c1_dups.columns = [f"{c1} {i+1}" for i in range(3)]

    c2_dups = pd.concat([df[c2] + i for i in range(3)], axis=1)
    c2_dups.columns = [f"{c2} {i+1}" for i in range(3)]

    # Rebuild DataFrame in desired order
    df = pd.concat(
        [
            df[['Peptide']],
            c1_dups,
            c2_dups,
            df[[f'{c1} E-score', f'{c2} E-score']]
        ],
        axis=1
    )

    return df
#librarymatrix = sortmatrix(matrix,libnumber,posnumber,matrixnum)
#librarymatrix.to_csv(r'/Users/mbp/Desktop/librarymatrix530.csv')
