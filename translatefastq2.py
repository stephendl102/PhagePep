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
from difflib import SequenceMatcher
import os
import subprocess
import numpy as np
import collections
import pandas as pd
import math
import time
# must pip install openpyxl

def readFastq(filename):
    """Reads FASTQ file and remove the special characters!"""

    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline() # skip name line
            seq = fh.readline().rstrip() # read base sequence
            fh.readline() # skip placeholder line
            qual = fh.readline().rstrip() #base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def diff(string_table, target_idx): #target_idx is the idx being looped through
    target_freq = string_table.iloc[target_idx,0]

    string_list = string_table.index
    target_string = string_list[target_idx]

    filter_percent= .03

    indexes = []
    Corrected = list(string_list)

    for i, string in enumerate(string_list): #finds all sequences that are similar to current sequence in loop
        if len(string) == len(target_string):
            diff_idx = ([j for j in range(len(string)) if string[j] != target_string[j]])
            if (len(diff_idx) == 1) and string_table.iloc[i,0] <= (filter_percent*target_freq) :
                indexes.append(i)
                Corrected[i]=target_string
            if (len(diff_idx) == 2) and string_table.iloc[i,0] <= (pow(filter_percent, 2)*target_freq) :
                indexes.append(i)
                Corrected[i]=target_string
            if (len(diff_idx) == 3) and string_table.iloc[i,0] <= (pow(filter_percent, 3)*target_freq) :
                indexes.append(i)
                Corrected[i]=target_string
    return indexes, Corrected
    
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

def translatefastq(mer,filenamefastq,startflank,endflank,filenameoutput,PhD7,collapse):
    ip = str(filenamefastq)
    #print(ip)
    #result = subprocess.run(["grep -o 'TCTCAC.*TCGGCCGAA' ", 'wc -l'], stdout=subprocess.PIPE,input=ip)
    #result.stdout.decode('utf-8')
    #print(result)
    print("Importing file " + filenamefastq)
    # import fastq file
    RawData, Qual = readFastq(filenamefastq)

    ReadDepth0 = len(RawData)
    print(ReadDepth0)
    print(RawData[0], Qual[0])
    print(RawData[1], Qual[1])

    print('Input file successfully imported')
    basepairs=3*mer                                               # calculate number of basepairs

    # put nucleotide sequences in array
    print('Isolating peptide sequences')
    First_idx = [s.find(endflank) for idx, s in enumerate(RawData)] # find end, PhD7 = 'GGTGGAGGT'
    Second_idx = [s.rfind(endflank) for idx, s in enumerate(RawData)] # find end, PhD7 = 'GGTGGAGGT'

    print('Eliminating misreads with double endflank')
    # misreads
    for j in range(len(First_idx)):
        if First_idx[j] != Second_idx[j]:                                   # deal with misreads that contain endflank twice
            First_idx[j]=0
        elif First_idx[j]<(3*mer+3):                                         # deal with reads that have endflank at beginning
            First_idx[j]=0
    End = First_idx

    tic = time.time()
    print('Eliminating misreads with no endflank')
    # reads without end flanking region
    RawData_df = pd.DataFrame(RawData)
    QualData_df = pd.DataFrame(Qual)
    indices_df = pd.DataFrame(End)

    endvalues = np.array(indices_df,dtype=object)
    searchval = 0
    ii = np.where(endvalues == 0)[0] # reads where missing or multiple endflank

    RawData_df = RawData_df.drop(ii)
    RawData_df = RawData_df.reset_index(drop=True)
    QualData_df = QualData_df.drop(ii)
    QualData_df = QualData_df.reset_index(drop=True)
    indices_df = indices_df.drop(ii)
    indices_df = indices_df.reset_index(drop=True)

    indicesMat = indices_df.iloc[:,0].tolist()
    RawData = RawData_df.iloc[:,0].tolist()
    QualData = QualData_df.iloc[:,0].tolist()
    toc = time.time()
    print(toc-tic)
    print(len(RawData))

    tic = time.time()
    print('Finding indices with startflank')
    # find indices with startflank
    rep1 = np.ones(len(indicesMat))
    indstart1 = [s-((basepairs+len(startflank))*rep1[idx]) for idx, s in enumerate(indicesMat)]        # start of where startflank should be
    indstart2 = [s-(basepairs+1)*rep1[idx] for idx, s in enumerate(indicesMat)]                           # end of where startflank should be
    startMat = np.chararray([len(indstart1),3])
    startMatCell = []
    for idx, s in enumerate(RawData):
        a = RawData[idx]
        startMatCell.append((a[int(indstart1[idx]):int(indstart2[idx])+1]))

    toc = time.time()
    print(toc-tic)
    print(len(startMatCell))

    print('Isolating indices of correct reads and dropping reads with incorrect length')
    tic = time.time()
    # isolate indices of correct reads
    indicesMat_df = pd.DataFrame(indicesMat)
    RawData_df = pd.DataFrame(RawData)
    QualData_df = pd.DataFrame(QualData)
    startMatCell_df = pd.DataFrame(startMatCell)

    values = np.array(startMatCell,dtype=object)
    searchval = startflank
    ii = np.where(values != searchval)[0] #dropping where the reads at location startflank are not equal to TCT
    RawData_df = RawData_df.drop(ii)
    RawData_df = RawData_df.reset_index(drop=True)
    QualData_df = QualData_df.drop(ii)
    QualData_df = QualData_df.reset_index(drop=True)
    indices_df = indices_df.drop(ii)
    indices_df = indices_df.reset_index(drop=True)
    startMatCell_df = startMatCell_df.drop(ii)
    startMatCell_df = startMatCell_df.reset_index(drop=True)

    indicesMat = indices_df.to_numpy()
    RawData = RawData_df.iloc[:,0].tolist()
    QualData = QualData_df.iloc[:,0].tolist()
    startMatCell = startMatCell_df.to_numpy()
    toc = time.time()
    print(toc-tic)
    print(len(startMatCell))
    tic = time.time()
    print('Making array of indices of sequences')
    # make array of indices of sequences
    rep1 = np.ones(len(indicesMat))
    indSeq1 = [s-(basepairs*rep1[idx]) for idx, s in enumerate(indicesMat)]
    indSeq2 = [s-rep1[idx] for idx, s in enumerate(indicesMat)]

    # pull out sequence
    a = []
    b = []
    for idx, s in enumerate(RawData):
        a.append(list(s))
    for idx, z in enumerate(QualData):
        b.append(list(z))
    RawData = a
    QualData = b
    NukeArray= []
    QualArray= []
    i=0
    for i in range(len(RawData)):                                # Store Peptide Sequences as Charcater List
        NukeArray.append(RawData[i][int(indSeq1[i]):(int(indSeq2[i])+1)])
        QualArray.append(QualData[i][int(indSeq1[i]):(int(indSeq2[i])+1)])
    toc = time.time()
    print(len(NukeArray))
    print(toc-tic)

    print('Collapsing Misreads')
    # determine frequencies
    SeqArray = []
    Qual = []
    Seqtext = ""
    Qualtext = ""
    for idx in range(len(NukeArray)):
        SeqArray.append(Seqtext.join(NukeArray[idx]))
        Qual.append(Qualtext.join(QualArray[idx]))
    tableSeq = collections.Counter(SeqArray)                                # calculate frequencies
    df = pd.DataFrame.from_dict(tableSeq, orient='index')
    print(df)
    DNAFreqTable = df.sort_values(by=0, ascending=False)
    print(DNAFreqTable)
    print('NumOfPep :',len(DNAFreqTable))
    ReadDepth = sum(DNAFreqTable.iloc[:, 0])
    print('ReadDepth ',ReadDepth)

    if collapse == True:
        DropFreq = ReadDepth0/100000
        DNAFreqTable = DNAFreqTable[DNAFreqTable.iloc[:, 0] >= DropFreq] # dropping Frequency <11
        print('Length After Dropping <', DropFreq, ':',len(DNAFreqTable))

        idx = 0
        while idx < .1*len(DNAFreqTable):
            Errors, Correctedidx = diff(DNAFreqTable,idx)
            DNAFreqTable.index = Correctedidx
            print(idx, len(Errors))
            idx += 1

    DNAFreqTable = DNAFreqTable.groupby(DNAFreqTable.index).agg('sum')
    RawData2 = DNAFreqTable.iloc[:,0].tolist()
    NukeArray= []
    i=0
    for i in range(len(RawData2)):                                # Store Peptide Sequences as Charcater List
        NukeArray.append(RawData[i][int(indSeq1[i]):(int(indSeq2[i])+1)])
    toc = time.time()
    print(len(NukeArray))
    print(toc-tic)

    #Make Sure NUKEARRAY IS being updated from DNAFREQTABLE
    print('Len NukArray ',len(NukeArray))
    print('Total Reads ', sum(DNAFreqTable.iloc[:, 0]))
    print('Getting rid of codons not used by PhD7 library')
    tic = time.time()
    # If PhD7 library, get rid of codons not used by PhD7 library
    if PhD7 == 1:
        badRead1= [row[2] for row in NukeArray]
        badRead2= [row[5] for row in NukeArray]
        badRead3=[row[8] for row in NukeArray]
        badRead4=[row[11] for row in NukeArray]
        badRead5=[row[14] for row in NukeArray]
        badRead6=[row[17] for row in NukeArray]
        badRead7=[row[20] for row in NukeArray]
        brRow = []
        for idx in range(len(badRead1)):
            if (badRead1[idx] == 'A' or badRead2[idx]=='A' or badRead2[idx]=='A' or badRead3[idx]=='A'
            or badRead4[idx]=='A' or badRead5[idx]=='A' or badRead6[idx]=='A' or badRead7[idx]=='A'
            or badRead1[idx] == 'C' or badRead2[idx]=='C' or badRead2[idx]=='C' or badRead3[idx]=='C'
            or badRead4[idx]=='C' or badRead5[idx]=='C' or badRead6[idx]=='C' or badRead7[idx]=='C'):
                brRow.append(idx)                                # find indices of instances of bad reads
        for i in sorted(brRow, reverse=True):
            #print(NukeArray[i])
            NukeArray.pop(i)
            #QualArray.pop(i)                                   # delete bad codon reads
    toc = time.time()
    print(toc-tic)
    print('Len NukeArray ',len(NukeArray))
    # convert to amino acids

    print('Converting to amino acid sequences')
    AAarray = []
    text = ""
    idx = 0
    for idx in range(len(DNAFreqTable.index)):
        AA = translate(DNAFreqTable.index[idx])
        AA = AA.replace('_', 'Q')
        AAarray.append(AA)
    DNAFreqTable['AA'] = AAarray
    DNAFreqTable['Normalized_Freq'] = DNAFreqTable[0] / ReadDepth0
    #DNAFreqTable['Normalized_FreqInsert'] = DNAFreqTable[0] / ReadDepth


    print('Determining frequencies')
    print(DNAFreqTable)
    # determine frequencies
    #   tableAA = collections.Counter(AAarray)                                # calculate frequencies
    #df = pd.DataFrame.from_dict(tableAA, orient='index')
    SeqFreqTable = DNAFreqTable.sort_values(by=0, ascending=False)                                   # sort table                                       # Sequences

    # display stats
    print("Number of total valid reads: ", sum(DNAFreqTable.iloc[:, 0]))
    print("Number of unique reads: ", len(SeqFreqTable))

    print('Starting export')
    # export to excel
    iterationsXLS = int(math.ceil(len(SeqFreqTable)/1000000))              # Filter out sequences that make up < .0001% of screen
    p=1                                                                    # initialize counter
    if iterationsXLS==1:
        print('Exporting to excel: one sheet')
        SeqFreqTable.to_excel(filenameoutput, header=False) # write excel file--> only 1e6 rows each sheet
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
                    df2.to_excel(writer, sheet_name="Try1", header=False)
                else:
                    df = SeqFreqTable.iloc[ind1:ind3,:]
                    df3 = df.copy()
                    df3.to_excel(writer, sheet_name="Try2", header=False) # write excel file--> only 1e6 rows each sheet
                p=p+1                                                              # count up
    libraryexcel=SeqFreqTable                                              # output
    print('Part 1 of program finished')
    return

print(snakemake.input)
print(snakemake.output)
for i in range(len(snakemake.input)):
    print(i)
    translatefastq(snakemake.config["mer"],snakemake.input[i],snakemake.config["startflank"],snakemake.config["endflank"],snakemake.output[i],snakemake.config["PHD7"],snakemake.config["collapse"])
