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


    # determine average fold change across targeted & non-targeted cell screens
    negnumber = libnumber - posnumber # number of non-targeted cell screens
    print('negative_number: ',negnumber)

    # Column definitions
    PosCol = counts.columns[0]                   # first positive screen
    cols_to_adjust = counts.columns[1:libnumber] # columns subject to bleedthrough correction

    # --------------------------------------------------
    # Minimum read count per column (used as a floor)
    # --------------------------------------------------
    mins = counts.min()
    print('raw mins:', mins)

    # Normalize minima by read depth
    ReadDepthRatio = ReadDepths / ReadDepths[0] #(# of reads normalized by first column reads, more reads = larger fraction)
    mins = mins / ReadDepthRatio                #(divides the minumum by 2 if read depth is twice as large)
    print('normalized mins:', mins)
    print('ReadDepthRatio:', ReadDepthRatio)

    # Replace NaNs with 0 before normalization
    counts = counts.fillna(0)

    # Normalize each column by its read depth
    for i in range(counts.shape[1]):
        counts.iloc[:, i] = counts.iloc[:, i] / ReadDepthRatio[i]

    print('Error Rate: ', Error)

    # --------------------------------------------------
    # Bleedthrough correction (only if error rate < 10%)
    # --------------------------------------------------
    if Error is not None and Error < 0.10:
        counts[cols_to_adjust] = counts[cols_to_adjust].subtract(
            Error * counts[PosCol],
            axis=0
        )
    else:
        print(
            f"WARNING: Error rate {Error:.3f} ≥ 10%. "
            "Bleedthrough correction skipped."
        )

    counts = counts.where(counts >= 1, np.nan)

    # Replace values less than the minimum with NaN
    for i in range(counts.shape[1]):
        counts.iloc[:,i] = counts.iloc[:,i].fillna(mins[i]) #set empty indices to normalized min of each column

    avgpos= counts.iloc[:,0:posnumber].sum(axis=1)/posnumber # average across targeted cell screens
    if average_negative:
        avgneg= counts.iloc[:,(posnumber):libnumber].sum(axis=1)/negnumber
    else:
        avgneg=counts.iloc[:,(posnumber):libnumber]
        avgneg = avgneg.max(axis = 1)
    if negnumber != 0:
        bigDiff=avgpos/avgneg # Difference between target and non-target
    else:
        bigDiff=avgpos
    bigDiff = np.array(bigDiff,dtype=object)
    print('NewMins ',counts.min())

    sortInd = (np.argsort(-bigDiff))
    bigDiffSort=[bigDiff[i] for i in sortInd]

    # create sorted matrix
    if matrixnum>counts.shape[0]:   # cannot put in more sequences than exist
        matrixnum=counts.shape[0]
    #counts = counts/(counts.min().min())

    Nums=counts.iloc[sortInd[0:matrixnum],:]
    Nums.reset_index(inplace = True, drop = True)
    Score= pd.DataFrame(bigDiff[sortInd[0:matrixnum]],columns=['Score'])
    Score.reset_index(inplace = True, drop = True)
    seqs2=seqs[sortInd[0:matrixnum]]
    seqs2.reset_index(inplace = True, drop = True)
    df=pd.concat([seqs2,Nums,Score], axis=1)
    print(df)

    print('Converting to amino acid sequences')
    AAarray = []
    text = ""
    idx = 0
    for idx in range(len(df.iloc[:,0])):
        AA = translate(df.iloc[idx,0])
        AA = AA.replace('_', 'Q')
        AAarray.append(AA)
    df['Peptide'] = AAarray

    return df

#librarymatrix = sortmatrix(matrix,libnumber,posnumber,matrixnum)
#librarymatrix.to_csv(r'/Users/mbp/Desktop/librarymatrix530.csv')
