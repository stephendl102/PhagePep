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

def normalizelibrary3(library,reflibrary,libraryname,index, posnum, average_negative=False):
    print('starting library normalization')
    # create table for reference library
    tableref1 = pd.read_excel(reflibrary,header=None)
    tableref1 = tableref1.iloc[:,:-1]
    tableref1.columns = ['Seq', 'RefLibrary','Peptide']

    # create table for library 1
    table1 = pd.read_excel(library,header=None)

    # Extract ErrorRate safely (if present)
    Error = table1.iloc[1, -1] if table1.shape[1] > 1 else None
    print('ErrorRate :', Error)

    # Trim metadata columns if present
    if table1.shape[1] > 2:
        table1 = table1.iloc[:, :-2]

        # ---- Normalize column structure FIRST ----
    if table1.shape[1] == 2:
        str_col = table1.select_dtypes(include=['object', 'string']).columns[0]
        other_col = table1.columns.difference([str_col])[0]

        table1 = pd.DataFrame({
            'Seq': table1[str_col],
            'Library': table1[other_col],
            'Peptide': table1[str_col]
        })

    elif table1.shape[1] == 3:
        table1.columns = ['Seq', 'Library', 'Peptide']
    else:
        raise ValueError("Unexpected number of columns in library table")

    table1 = table1[table1['Library'] != 1] #drops all sequences with only a single read
    print(table1.loc[table1['Seq']=='ADARYKS'])
    print(table1.columns)

        # ---- Compute libDep AFTER structure is known ----
    if 'Library' in table1.columns and table1['Library'].notna().any():
        libDep = table1['Library'].sum()
    else:
        libDep = None
    print('libDep :', libDep)

    # join the reference library and library together
    mergelibrary = tableref1.merge(table1,how='outer', on='Seq', sort=True)

    # Total reads per column (before filtering)
    readlengths = mergelibrary.drop(columns='Seq').select_dtypes(include='number').sum()
    print(readlengths)

    # Count overlapping sequences
    ref_overlap = mergelibrary.dropna().shape[0]
    print("Number of Sequences Overlapping with Reference:", ref_overlap)

    # Rename for clarity
    lib = mergelibrary['Library']
    reflib = mergelibrary['RefLibrary']

    # Keep only sequences present in library
    mask = lib > 0
    filtered = mergelibrary.loc[mask].copy()
    print(filtered)

    # Recompute overlap after filtering
    ref_overlap = filtered.dropna().shape[0]
    print("Number of Sequences Overlapping with Reference:", ref_overlap)

    # Normalize to read length
    print('normalizing to read length')

    libraryNtable = (
        filtered[['Seq', 'Library']]
        .rename(columns={'Library': libraryname})
        .reset_index(drop=True)
    )
    print('library normalization complete')
    print(libraryNtable)

    return libraryNtable, libDep, Error
