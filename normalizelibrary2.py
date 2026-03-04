# function to normalize libraries to read depth and amplification bias
#
# inputs:
#       library: string containing the filepath and filename
#       of the excel file containing the translated library. Library
#       should be translated using the function "translatefastq", as set up
#       in the program "NGSanalyzeMain.m".
#
#       reflibrary: string containing the filepath and filename of the
#       excel file containing the translated reference library. Reference
#       library refers to a sequenced amplification of the naive library
#       that has been diluted 1:10. The aliquot should have the same lot
#       number as the other input library. The reference library should be
#       translated using the function "translatefastq.m", as set up in the
#       program "NGSanalyzeMain.m".
#
#       libraryname: string containing a unique name for this library
#
# output:
#       libraryNtable: table containing sequences and normalized read
#       frequencies for a library
#
# Created by Lindsey Brinton at the University of Virginia, 2015


##changes:divide by reference squared. If missing from reference but in screen, replace with =1

import pandas as pd
import numpy as np
import math
import os
from scipy import stats

def normalizelibrary2(library,reflibrary,libraryname,index, posnum, average_negative=False):
    print('starting library normalization')
    # create table for reference library
    tableref1 = pd.read_excel(reflibrary,header=None)
    tableref1 = tableref1.iloc[:,:-1]
    tableref1.columns = ['Seq', 'RefLibrary','AA']

    # create table for library 1
    table1 = pd.read_excel(library,header=None)
    libDep = table1.iloc[0, 1] / table1.iloc[0, 3]
    table1 = table1.iloc[:,:-1]
    table1.columns = ['Seq','Library','AA']

    # join the reference library and library together
    mergelibrary = tableref1.merge(table1,how='outer', on='Seq', sort=True)
    #print(mergelibrary)

    readlengths = mergelibrary.drop('Seq','columns').sum(skipna=True)    # get readlengths before changes

    RefOverlap = len(mergelibrary.dropna(axis=0))
    #print(mergelibrary.dropna(axis=0))
    print("Number of Sequences Overlapping with Reference: ",RefOverlap)

    Seqs = mergelibrary['Seq'] # array of strings (sequences)
    Nums = mergelibrary.drop(['Seq','AA_x','AA_y'],'columns')       # array of doubles (read frequencies)

    lib=Nums.iloc[:,1]
    reflib=Nums.iloc[:,0]

    greaterthan10ind = []
    for i in range(len(lib)):
        if lib[i]>0:
            greaterthan10ind.append(i)


    mergelibrary.drop(mergelibrary.index[greaterthan10ind],inplace=True)
    RefOverlap = len(mergelibrary.dropna(axis=0))
    print("Number of Sequences Overlapping with Reference: ",RefOverlap)

    lib=lib.iloc[greaterthan10ind]
    reflib=reflib.iloc[greaterthan10ind]
    Seqs=Seqs.iloc[greaterthan10ind]
    # normalize to read length (division)
    print('normalizing to read length')
    #nonzeromode = []
    #for i in range(len(reflib)):
    #    if reflib.iloc[i]>0:
    #        nonzeromode.append(i)
    #mode = stats.mode(reflib.iloc[nonzeromode])[0]
    #print("mode ",mode[0])
    #reflib = reflib.fillna(mode[0])
    reflib = reflib.fillna(1) # make NaN = low frequency
    reflib = np.divide(reflib,readlengths[0])                  # normalize reference library
    #lib = np.divide(lib,readlengths[1])                        # normalize library

    # normalize library to reference library
    print('normalizing to reference library') #squared
    #if index<=posnum:                 #index<=posnum
    #    libraryN = np.divide(lib,reflib)
    #else:
    libraryN = np.divide(lib,(min(reflib)))

    # combine sequences and quantities
    libraryNtable = pd.DataFrame({'Seq':Seqs, libraryname:libraryN})
    libraryNtable.reset_index(inplace = True, drop = True)
    print('library normalization complete')

    return libraryNtable, libDep
