import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import math
from scipy import stats
#from normalizelibrary import *
from normalizelibrary3 import *
from sortmatrix_1_1_26 import *

# Main PHASTpep file code
# created by Lindsey Brinton at the University of Virginia, 2015
# updated by Lindsey Brinton at the University of Virginia, 2016
# Translated to Python by Stephen Lees, 2022

def submit2_Callback(posnum,matnum,outputfile2,libfile,reffile,average_negative=False):
    libnum = len(libfile)
    print('libnum is: ', str(libnum))
    print(libfile[0])
    print([reffile[0]])

        # file 1
    libraryTableHolder = {}
    ReadDepths = []

    print('reading in set 1')
    libraryTableHolder[0], libDep, Error = normalizelibrary3(libfile[0],reffile[0],str(os.path.basename(libfile[0])),0, posnum, average_negative)
    ReadDepths.append(libDep)
    print('ReadDepths',ReadDepths)
    alltogether = libraryTableHolder[0]
    print(alltogether)

    for i in range(len(libfile)-1):
        print('reading in set ', str(i+2))
        libraryTableHolder[i+1], libDep, _ = normalizelibrary3(libfile[i+1],reffile[i+1],str(os.path.basename(libfile[i+1])),i+1, posnum, average_negative)
        ReadDepths.append(libDep)
        print(libraryTableHolder[i+1])
        alltogether = alltogether.merge(libraryTableHolder[i+1],how='outer', on='Seq', sort=True)

    print('ReadDepths ',ReadDepths)
    print(alltogether.loc[alltogether['Seq']=='ADARYKS'])
    # sort matrix
    print(alltogether)
    sortedmatrix=sortmatrix(alltogether,posnum,matnum,ReadDepths,Error,average_negative)
    print(sortedmatrix)

    # write table
    print('creating matrix file')
    sortedmatrix.to_excel(outputfile2)
    #writetable(sortedmatrix,outputfile2)
    print('matrix file completed')
    return

print(snakemake.input)
print(snakemake.output)
submit2_Callback(snakemake.config["posnum"], snakemake.config["matnum"], snakemake.output[0], snakemake.input, (snakemake.config["positive_reffile"]+snakemake.config["negative_reffile"]),snakemake.config["average_negative"])
