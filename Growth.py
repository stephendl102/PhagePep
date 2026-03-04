import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import math
from scipy import stats
#from normalizelibrary import *
from normalizelibrary2 import *
from sortmatrix3_11_24 import *

# Main PHASTpep file code
# created by Lindsey Brinton at the University of Virginia, 2015
# updated by Lindsey Brinton at the University of Virginia, 2016
# Translated to Python by Stephen Lees, 2022

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

def ConcatRounds(Files,index):
        print('starting library normalization')
        # create table for reference library
        ExcelHolder = {}
        for i in range(len(Files)):
            tableref = pd.read_excel(Files[i],header=None)
            tableref.replace(0, np.nan, inplace=True)
            tableref = tableref.dropna()
            tableref.columns = ['Seq', str(os.path.basename(Files[i])),'AA',str(os.path.basename(Files[i]))]
            print(tableref)
            ExcelHolder[i] = tableref

        mergelibrary = ExcelHolder[1].merge(ExcelHolder[0],how='outer', on='Seq', sort=True)
        print("one")
        Nums = mergelibrary.drop(['AA_x','AA_y'],'columns')
        print("two")

        Seqs = Nums['Seq'] # array of strings (sequences)

        Numsfilt = Nums[Nums.iloc[:,1] != 0]
        #Nums = Numsfilt.dropna()
        Nums.reset_index(inplace = True, drop = True)
        print(Nums)
        print(len(Files))

        for i in range(len(Files)-2):
            Nums = Nums.merge(ExcelHolder[i+2],how='outer', on='Seq', sort=True)
            Nums = Nums.drop(['AA'],'columns')
            print(Nums)
        Seqs = Nums['Seq'] # array of strings (sequences)
        """
        AAarray = []
        text = ""
        idx = 0
        for idx in range(len(Nums['Seq'])):
            AA = translate(Nums['Seq'][idx])
            AA = AA.replace('_', 'Q')
            AAarray.append(AA)
        """
        Nums['3over1'] = Nums.iloc[:,-1] / Nums.iloc[:,-3]
        Nums['3over2'] = Nums.iloc[:,-2] / Nums.iloc[:,-6]
        #Nums['AA'] = AAarray
        print('library normalization complete')
        #libraryNtable.to_csv(r'/Users/stephen/Desktop/Desktop/NGS_Data/StrepIOMEGA/CxBxR2/CxBxR2.csv')

        return Nums


def Diversity(RoundNum, outputfile2,libfile):
    print('Round Number is: ' + str(RoundNum))
    print(libfile[0])
    index = 0

    # file 1
    print('reading in set 1')
    print(len(libfile))
    alltogether = ConcatRounds(libfile,index)

    # sort matrix
    print(alltogether)

    # write table
    print('creating matrix file')
    alltogether.to_excel(outputfile2)
    #writetable(sortedmatrix,outputfile2)
    print('matrix file completed')
    #print(toc-tic, 'sec Elapsed')
    return

print(snakemake.input)
print(snakemake.output)
Diversity(snakemake.config["Rounds"], snakemake.output[0], snakemake.input)
