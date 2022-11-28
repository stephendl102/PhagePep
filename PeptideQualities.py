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
import math
import time
import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from Bio import motifs
# must pip install openpyxl

def createFasta(PepFile, filenameoutput):
    tableref1 = pd.read_excel(PepFile,header=0)
    print(tableref1)
    list_name = tableref1.iloc[-2000:,1].tolist()
    list_seq = tableref1.iloc[-2000:,0].tolist()
    print(list_seq)


    ofile = open(filenameoutput, "w")

    for i in range(len(list_seq)):
        ofile.write(">" + str(list_name[i]) + "\n" + str(list_seq[i]) + "\n")

    ofile.close()

    return

def AddCharge(PepFile,filenameoutput):
    """Reads Excel"""

    print('File is: ' + str(PepFile))
    #tableref1 = pd.read_excel(PepFile,header=0,index_col=0)
    tableref1 = pd.read_excel(PepFile)
    createFasta(tableref1)
    #tableref1.columns = ['Peptide', 'Frequency']
    RawData = tableref1.iloc[:,0].tolist()
    print(RawData[0:5])

    charge = []
    weight = []
    instability = []
    aromaticity = []
    gravy = []
    for idx, s in enumerate(RawData):
        X = ProteinAnalysis(s)
        charge.append(X.charge_at_pH(7.0))
        weight.append(X.molecular_weight())
        instability.append(X.instability_index())
        aromaticity.append(X.aromaticity())
        X.get_amino_acids_percent()['A']
        gravy.append(X.gravy())
    tableref1['Charge'] = charge
    tableref1['Molecular_Weight'] = weight
    tableref1['Instability'] = instability
    tableref1['Aromaticity'] = aromaticity
    tableref1['gravy'] = gravy
    SeqFreqTable = tableref1

    print('Starting export')
    # export to excel
    iterationsXLS = int(math.ceil(len(SeqFreqTable)/1000000))              # Determine if for loops necessary
    p=1                                                                    # initialize counter
    if iterationsXLS==1:
        print('Exporting to excel: one sheet')
        SeqFreqTable.to_excel(filenameoutput, header=True) # write excel file--> only 1e6 rows each sheet
    else:
        for w in range(iterationsXLS):
            print('Exporting to excel: multiple sheets')
            sheetI = str(p)                                                   # determine sheet to use
            ind2 = int((w+1)*1e6)
            ind1 = int((ind2-1e6)+1)                                                # find indexes within Sequence Array
            ind3 = len(SeqFreqTable)
            with pd.ExcelWriter(filenameoutput) as writer:
                if (ind3-ind1)>(1e6-1):
                    df = SeqFreqTable.iloc[ind1:ind2,:]
                    df2 = df.copy()
                    df2.to_excel(writer, sheet_name="Try1", header=True)
                else:
                    df = SeqFreqTable.iloc[ind1:ind3,:]
                    df3 = df.copy()
                    df3.to_excel(writer, sheet_name="Try2", header=True) # write excel file--> only 1e6 rows each sheet
                p=p+1                                                              # count up
    libraryexcel=SeqFreqTable                                              # output
    print('Charge Information Finished')
    #SeqFreqTable.plot(x='col_name_1', y='col_name_2', style='o')
    return

print(snakemake.input)
for i in range(len(snakemake.input)):
    createFasta(snakemake.input[i],snakemake.output[i])
