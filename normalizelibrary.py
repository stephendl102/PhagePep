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
from scipy import stats

def normalizelibrary(library,reflibrary,libraryname,average_negative=False):
    print('starting library normalization')
    # create table for reference library
    tableref1 = pd.read_excel(reflibrary,header=None)
    tableref1.columns = ['Peptide', 'RefLibrary']

    # create table for library 1
    table1 = pd.read_excel(library,header=None)
    table1.columns = ['Peptide', 'Library']
    print(table1)

    # join the reference library and library together
    mergelibrary = tableref1.merge(table1,how='outer', on='Peptide', sort=True)
    print(mergelibrary)

    readlengths = mergelibrary.drop('Peptide','columns').sum(skipna=True)    # get readlengths before changes
    print(readlengths)

    RefOverlap = len(mergelibrary.dropna(axis=0))
    print("Number of Sequences Overlapping with Reference: ",RefOverlap)

    Seqs = mergelibrary['Peptide'] # array of strings (sequences)
    Nums = mergelibrary.drop('Peptide','columns')       # array of doubles (read frequencies)

    # only keep sequences with frequencies >10 in library
    #print('deleting sequences with frequencies < 11')
    lib=Nums.iloc[:,1]
    reflib=Nums.iloc[:,0]

    greaterthan10ind = []
    for i in range(len(lib)):
        if lib[i]>10:
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
    lib = np.divide(lib,readlengths[1])                        # normalize library

    # normalize library to reference library
    print('normalizing to reference library')
    libraryN = np.divide(lib,reflib)

    # combine sequences and quantities
    libraryNtable = pd.DataFrame({'Peptide':Seqs, libraryname:libraryN})
    libraryNtable.reset_index(inplace = True, drop = True)
    print('library normalization complete')
    #libraryNtable.to_csv(r'/Users/stephen/Desktop/Desktop/NGS_Data/StrepIOMEGA/CxBxR2/CxBxR2.csv')

    return libraryNtable

#normalizelibrary("/Users/mbp/Desktop/PythonPhastPep/exampleFilesForPart1 copy/exampleAforPHASTpep/exampleAforPHASTpep.xlsx","/Users/mbp/Desktop/PythonPhastPep/exampleFilesForPart2 copy/exampleReferencePHASTpep.xlsx","Hello.xlsx")
