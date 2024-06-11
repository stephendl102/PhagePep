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
import statistics
# must pip install openpyxl

def diff(string_table, target_idx,filter_percent): #target_idx is the idx being looped through
    target_freq = string_table.iloc[target_idx,0]

    string_list = string_table.index
    target_string = string_list[target_idx]
    indexes = []
    Corrected = list(string_list)

    for i, string in enumerate(string_list): #finds all sequences that are similar to current sequence in loop
        if len(string) == len(target_string):
            diff_count = sum([1 for j in range(len(string)) if string[j] != target_string[j]])
            if diff_count == 1 and string_table.iloc[i,0] <= (filter_percent*target_freq) :
                indexes.append(i)
                Corrected[i]=target_string
            if diff_count == 2 and string_table.iloc[i,0] <= (pow(filter_percent, 2)*target_freq) :
                indexes.append(i)
                Corrected[i]=target_string
            if diff_count == 3 and string_table.iloc[i,0] <= (pow(filter_percent, 3)*target_freq) :
                indexes.append(i)
                Corrected[i]=target_string
    return indexes, Corrected

# Function to convert ASCII character to Phred score
def ascii_to_phred(char):
    return ord(char) - 33

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
    ticstart = time.time()
    ip = str(filenamefastq)
    print("Importing file " + filenamefastq)
    basepairs=3*mer                                               # calculate number of basepairs

    tic = time.time()
    pattern = str(startflank) + '.*' + str(endflank)
    print(pattern)
    QualDataCommand = ["grep", "-A", "2", "-o", pattern, filenamefastq]
    awk_command = ["awk", "NR % 4 == 3"]

    grep_process = subprocess.Popen(QualDataCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    awk_process = subprocess.Popen(awk_command, stdin=grep_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    grep_process.stdout.close()
    output, errors = awk_process.communicate()
    Qual = output.decode().splitlines()
    print(len(Qual))
    toc = time.time()
    print('Time Taken: ', toc-tic)
    tic = time.time()

    # Run the command and capture the output
    RawDataCommand = ["grep", "-o", pattern, filenamefastq]
    RawDataResult = subprocess.run(RawDataCommand, capture_output=True, text=True)
    # Check for errors
    if RawDataResult.returncode != 0:
        print("Error:", result.stderr)
    else:
        # Split the output into a list of lines
        RawData = RawDataResult.stdout.splitlines()

    toc = time.time()
    print('Time Taken: ', toc-tic)

    ReadDepth0 = len(RawData)
    print('Total Reads: ', ReadDepth0)
    print(RawData[0], Qual[0])
    print(RawData[1], Qual[1])

    print('Input file successfully imported')

    # put nucleotide sequences in array
    print('Ensuring all reads are :', str(basepairs+len(startflank)+len(endflank)))
    filtered_list1 = [item[len(startflank):-(len(endflank))] for item in RawData if len(item) == (basepairs+len(startflank)+len(endflank))]

    #32 is temporary, it must change based on where startflank occurs
    Qual = [Qual[idx][32+len(startflank):32+(basepairs+len(startflank))] for idx, item in enumerate(RawData) if len(item) >= (basepairs+len(startflank)+len(endflank))]
    RawData = filtered_list1

    ReadDepth0 = len(RawData)
    print('Total Reads: ', ReadDepth0)
    print(RawData[0], Qual[0])
    print(RawData[1], Qual[1])


    # Calculate the average quality score for each read
    average_quality_scores = []
    for quality_scores in Qual:
        phred_scores = [ascii_to_phred(char) for char in quality_scores]
        average_quality = sum(phred_scores) / len(phred_scores)
        average_quality_scores.append(average_quality)

    std_dev = statistics.stdev(average_quality_scores)
    print('std_dev :',std_dev)

    # Calculate the overall average of the average quality scores
    overall_average_quality = (sum(average_quality_scores) / len(average_quality_scores))
    Error_Rate = 1/(10**(overall_average_quality/10))

    adjusted_average_quality = overall_average_quality-2*std_dev
    Adj_Error_Rate = 1/(10**(adjusted_average_quality/10))

    filter_percent= Adj_Error_Rate*basepairs

    print("Average error rate per base call: ", Error_Rate)
    print("5th percentile average error rate per base call: ", Adj_Error_Rate)


    print("Length after filtering wrong length reads:", len(RawData))
    toc = time.time()
    print(toc-tic)

    tic = time.time()
    print('Making Nucleotide and Quality Array')

    NukeArray = [list(item) for item in RawData]
    QualData = [list(item) for item in Qual]

    print('NukeArray0: ', NukeArray[0])
    i=0
    toc = time.time()
    print(len(NukeArray))
    print(toc-tic)

    print('Counting Reads')
    # determine frequencies
    SeqArray = [''.join(seq) for seq in NukeArray]
    Qual = [''.join(qual) for qual in QualData]

    tableSeq = collections.Counter(SeqArray)                                # calculate frequencies
    df = pd.DataFrame.from_dict(tableSeq, orient='index')
    print(df)
    DNAFreqTable = df.sort_values(by=0, ascending=False)
    print(DNAFreqTable)
    print('NumOfPep :',len(DNAFreqTable))
    ReadDepth = sum(DNAFreqTable.iloc[:, 0])
    print('ReadDepth ',ReadDepth)

    if collapse == True:
        print('Collapsing Misreads')
        DropFreq = ReadDepth0/100000
        DNAFreqTable = DNAFreqTable[DNAFreqTable.iloc[:, 0] >= DropFreq] # dropping Frequency <11
        print('Length After Dropping <', DropFreq, ':',len(DNAFreqTable))

        idx = 0  # Initialize the index
        print('Filter Percentage: ', filter_percent)
        while idx < 0.1 * len(DNAFreqTable):
            # Check the condition to end the loop
            if DNAFreqTable.iloc[idx,0] * filter_percent < min(DNAFreqTable.iloc[:,0]):
                print('break')
                break
            Errors, Correctedidx = diff(DNAFreqTable,idx,filter_percent)
            DNAFreqTable.index = Correctedidx
            print(idx, len(Errors))
            idx += 1

    print('Filter Percentage: ', filter_percent)
    DNAFreqTable = DNAFreqTable.groupby(DNAFreqTable.index).agg('sum')
    RawData2 = DNAFreqTable.iloc[:,0].tolist()
    NukeArray = [list(item) for item in DNAFreqTable.index.tolist()]
    toc = time.time()
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
    DNAFreqTable['5th_Percentile_Per_base_error'] = Adj_Error_Rate
    #DNAFreqTable['Normalized_FreqInsert'] = DNAFreqTable[0] / ReadDepth
    print('Total Reads ', sum(DNAFreqTable.iloc[:, 0]))


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
    tocend = time.time()
    print('Total Time :',tocend-ticstart)
    return

print(snakemake.input)
print(snakemake.output)
for i in range(len(snakemake.input)):
    print(i)
    translatefastq(snakemake.config["mer"],snakemake.input[i],snakemake.config["startflank"],snakemake.config["endflank"],snakemake.output[i],snakemake.config["PHD7"],snakemake.config["collapse"])
