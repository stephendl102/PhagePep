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


def read_fasta_file(file_path):
    sequences = {}
    current_sequence = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_sequence = line[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line

    return sequences


def calculate_amino_acid_percent(sequences):
    amino_acid_counts = {}
    total_length = 0

    for sequence in sequences.values():
        total_length += len(sequence)
        for amino_acid in sequence:
            amino_acid_counts[amino_acid] = amino_acid_counts.get(amino_acid, 0) + 1

    amino_acid_percentages = {}
    for amino_acid, count in amino_acid_counts.items():
        percentage = (count / total_length) * 100
        amino_acid_percentages[amino_acid] = percentage

    return amino_acid_percentages


print(snakemake.input)
for i in range(len(snakemake.input)):
    amino_acid_percentages = calculate_amino_acid_percent(read_fasta_file(snakemake.input[i]))
    print(amino_acid_percentages)
    print("Amino Acid Percentages:")
    for amino_acid, percentage in amino_acid_percentages.items():
        print(f"{amino_acid}: {percentage:.2f}%")
