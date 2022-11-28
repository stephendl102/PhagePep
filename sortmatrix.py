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

def sortmatrix(matrix,libnumber,posnumber,matrixnum,average_negative=False):

    print('sorting the matrix')

    # determine average fold change across targeted & non-targeted cell screens
    negnumber = libnumber - posnumber # number of non-targeted cell screens
    Seqs = matrix.iloc[:,0]
    Numsdf = matrix.drop(labels=['Peptide'], axis=1)
    #Nums= matrix.iloc[:,1:].values.tolist()
    mins = (Numsdf.min(axis=0))
    for i in range(Numsdf.shape[1]):
        Numsdf.iloc[:,i] = Numsdf.iloc[:,i].fillna(mins[i]) #set empty indices to min of each column

    avgpos= Numsdf.iloc[:,0:posnumber].sum(axis=1)/posnumber # average across targeted cell screens
    if average_negative:
        avgneg= Numsdf.iloc[:,(posnumber):libnumber].sum(axis=1)/negnumber
    else:
        avgneg=Numsdf.iloc[:,(posnumber):libnumber]
        avgneg = avgneg.max(axis = 1)
    bigDiff=avgpos/avgneg # Difference between target and non-target
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
    print(df.shape)
    #df = df.drop(df[(df['libOne'] == df.iloc[-1,1])& (df['libTwo'] == df.iloc[-1,2]) & (df['libThree'] == df.iloc[-1,3])].index)  #&(df['libFour'] == df.iloc[-1,4])&(df['libFive'] == df.iloc[-1,5])&(df['libSix'] == df.iloc[-1,6])#&(df['libSeven'] == df.iloc[-1,7])&(df['libEight'] == df.iloc[-1,8])&(df['libNine'] == df.iloc[-1,9])&(df['libTen'] == df.iloc[-1,10])&(df['libEleven'] == df.iloc[-1,11])&(df['libTwelve'] == df.iloc[-1,12])
    return df

#librarymatrix = sortmatrix(matrix,libnumber,posnumber,matrixnum)
#librarymatrix.to_csv(r'/Users/mbp/Desktop/librarymatrix530.csv')
