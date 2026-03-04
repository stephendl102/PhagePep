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

def translate(seq):

    table = {
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        'CCN':'P', 'GGN':'G', 'ACN':'T', 'CGN':'R',
        'GCN':'A', 'CTN':'L', 'TCN':'S', 'GTN':'V',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon in table:
                protein+= table[codon]
            else:
                protein+='X'
    return protein

def sortmatrix(matrix,libnumber,posnumber,matrixnum,average_negative=False):

    print('sorting the matrix')
    pd.set_option("display.precision", 10)
###

###
    # determine average fold change across targeted & non-targeted cell screens
    negnumber = libnumber - posnumber # number of non-targeted cell screens
    Seqs = matrix.iloc[:,0]
    Numsdf = matrix.drop(labels=['Seq'], axis=1)
    #Nums= matrix.iloc[:,1:].values.tolist()
    cols_to_adjust = Numsdf.columns[1:libnumber]  # Exclude the first column (index 0)
    PosCol = Numsdf.columns[0]
    ReadDepthRatio = ReadDepthRatio[]
    print(Numsdf.min())
    print('Numsdfmin^')
    #Numsdf = (Numsdf / Numsdf.min())

    #mins = Numsdf.min()
    Numsdf[cols_to_adjust] = Numsdf[cols_to_adjust].fillna(0)
    #Numsdf[PosCol] = Numsdf[PosCol].fillna(0)
    print(Numsdf)
    print(0.006 * Numsdf[PosCol])
    print(Numsdf[cols_to_adjust])

    Numsdf[cols_to_adjust] = Numsdf.apply(lambda row: row[cols_to_adjust] - 0.006*row[PosCol], axis=1)
    print(Numsdf)
    #print(Numsdf[cols_to_adjust].sub(Numsdf[PosCol]))
    Numsdf[cols_to_adjust] = Numsdf[cols_to_adjust].where(Numsdf[cols_to_adjust] >= 1, 1)


    # Replace values less than the minimum with NaN
    for i in range(Numsdf.shape[1]):
        Numsdf.iloc[:,i] = Numsdf.iloc[:,i].fillna(1) #set empty indices to min of each column

    avgpos= Numsdf.iloc[:,0:posnumber].sum(axis=1)/posnumber # average across targeted cell screens
    if average_negative:
        avgneg= Numsdf.iloc[:,(posnumber):libnumber].sum(axis=1)/negnumber
    else:
        avgneg=Numsdf.iloc[:,(posnumber):libnumber]
        avgneg = avgneg.max(axis = 1)
    if negnumber != 0:
        bigDiff=avgpos/avgneg # Difference between target and non-target
    else:
        bigDiff=avgpos
    bigDiff = np.array(bigDiff,dtype=object)

    sortInd = (np.argsort(-bigDiff))
    bigDiffSort=[bigDiff[i] for i in sortInd]

    # create sorted matrix
    if matrixnum>Numsdf.shape[0]:   # cannot put in more sequences than exist
        matrixnum=Numsdf.shape[0]

    Nums=Numsdf.iloc[sortInd[0:matrixnum],:]
    Nums.reset_index(inplace = True, drop = True)
    Score= pd.DataFrame(bigDiff[sortInd[0:matrixnum]],columns=['Score'])
    Score.reset_index(inplace = True, drop = True)
    Seqs2=Seqs[sortInd[0:matrixnum]]
    Seqs2.reset_index(inplace = True, drop = True)
    df=pd.concat([Seqs2,Nums,Score], axis=1)
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
