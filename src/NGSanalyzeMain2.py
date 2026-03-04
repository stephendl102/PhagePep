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

def submit2_Callback(posnum, matnum, outputfile2,libfile,reffile,average_negative=False):
    libnum = len(libfile)
    print('libnum is: ' + str(libnum))
    print(libfile[0])
    print([reffile[0]])
    index = 0

    # file 1
    print('reading in set 1')
    library1table = normalizelibrary2(libfile[0],reffile[0],str(os.path.basename(libfile[0])),index, posnum, average_negative)
    alltogether = library1table
    index = 1

    if libnum>1:
        index = 2
        # file 2
        print('reading in set 2')
        library2table = normalizelibrary2(libfile[1],reffile[1],str(os.path.basename(libfile[1])),index, posnum, average_negative)

        # combine file 1 and file 2
        print('combining sets 1 & 2')
        A = alltogether.merge(library2table,how='outer', on='Seq', sort=True)
        alltogether=A

    # file 3
    if libnum>2:
        index = 3
        print('adding set 3')
        library3table = normalizelibrary2(libfile[2],reffile[2],str(os.path.basename(libfile[2])),index, posnum, average_negative)
        B = alltogether.merge(library3table,how='outer', on='Seq', sort=True) # add file 3
        alltogether=B

        # file 4
    if libnum>3:
        index = 4
        print('adding set 4')
        library4table = normalizelibrary2(libfile[3],reffile[3],str(os.path.basename(libfile[3])),index, posnum, average_negative)
        C = alltogether.merge(library4table,how='outer', on='Seq', sort=True) # add file 4
        alltogether=C

        # file 5
    if libnum>4:
        index = 5
        print('adding set 5')
        library5table = normalizelibrary2(libfile[4],reffile[4],str(os.path.basename(libfile[4])), index, posnum, average_negative)
        D = alltogether.merge(library5table,how='outer', on='Seq', sort=True) # add file 5
        alltogether=D

    # file 6
    if libnum>5:
        index = 6
        print('adding set 6')
        print(libfile[5],reffile[5])
        library6table = normalizelibrary2(libfile[5],reffile[5],str(os.path.basename(libfile[5])),index, posnum, average_negative)
        E = alltogether.merge(library6table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=E

        # file 7
    if libnum>6:
        index = 7
        print('adding set 7')
        library7table = normalizelibrary2(libfile[6],reffile[6],str(os.path.basename(libfile[6])),index, posnum, average_negative)
        F = alltogether.merge(library7table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=F

    # file 8
    if libnum>7:
        index = 8
        print('adding set 8')
        library8table = normalizelibrary2(libfile[7],reffile[7],str(os.path.basename(libfile[7])),index, posnum, average_negative)
        G = alltogether.merge(library8table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=G

    # file 9
    if libnum>8:
        index = 9
        print('adding set 9')
        library9table = normalizelibrary2(libfile[8],reffile[8],str(os.path.basename(libfile[8])),index,posnum,average_negative)
        H = alltogether.merge(library9table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=H

    # file 10
    if libnum>9:
        index = 10
        print('adding set 10')
        library10table = normalizelibrary2(libfile[9],reffile[9],str(os.path.basename(libfile[9])),index, posnum,average_negative)
        I = alltogether.merge(library10table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=I

        # file 11
    if libnum>10:
        index = 11
        print('adding set 11')
        library11table = normalizelibrary2(libfile[10],reffile[10],str(os.path.basename(libfile[10])),index,posnum,average_negative)
        J = alltogether.merge(library11table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=J

    # file 12
    if libnum>11:
        print('adding set 12')
        library12table = normalizelibrary2(libfile[11],reffile[11],str(os.path.basename(libfile[11])),index,posnum,average_negative)
        K = alltogether.merge(library12table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=K

    # file 13
    if libnum>12:
        print('adding set 13')
        library13table = normalizelibrary2(libfile[12],reffile[12],str(os.path.basename(libfile[12])),index,posnum,average_negative)
        L = alltogether.merge(library13table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=L


    # file 14
    if libnum>13:
        print('adding set 14')
        library14table = normalizelibrary2(libfile[13],reffile[13],str(os.path.basename(libfile[13])),index,posnum,average_negative)
        M = alltogether.merge(library14table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=M

    # file 15
    if libnum>14:
        print('adding set 15')
        library15table = normalizelibrary2(libfile[14],reffile[14],str(os.path.basename(libfile[14])),index,posnum,average_negative)
        N = alltogether.merge(library15table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=N

    # file 16
    if libnum>15:
        print('adding set 16')
        library16table = normalizelibrary2(libfile[15],reffile[15],str(os.path.basename(libfile[15])),index,posnum,average_negative)
        O = alltogether.merge(library16table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=O

    # file 17
    if libnum>16:
        print('adding set 16')
        library16table = normalizelibrary2(libfile[16],reffile[16],str(os.path.basename(libfile[16])),index,posnum,average_negative)
        P = alltogether.merge(library16table,how='outer', on='Seq', sort=True) # add file 6
        alltogether=P

    # file 18
    if libnum>17:
        print('adding set 18')
        libfile18 = getappdata(hMainGui,'libfile18')
        reffile18 = getappdata(hMainGui,'reffile18')
        library18table = normalizelibrary(libfile18,reffile18,'libEighteen')
        [Q,iq,jq]=outerjoin(P,library18table,'MergeKeys',true) # add file 18
        clear('iq','jq','P','library18table')
        alltogether=Q

    # file 19
    if libnum>18:
        print('adding set 19')
        libfile19 = getappdata(hMainGui,'libfile19')
        reffile19 = getappdata(hMainGui,'reffile19')
        library19table = normalizelibrary(libfile19,reffile19,'libNineteen')
        [R,ir,jr]=outerjoin(Q,library19table,'MergeKeys',true) # add file 19
        clear('ir','jr','Q','library19table')
        alltogether=R

    # file 20
    if libnum>19:
        print('adding set 20')
        libfile20 = getappdata(hMainGui,'libfile20')
        reffile20 = getappdata(hMainGui,'reffile20')
        library20table = normalizelibrary(libfile20,reffile20,'libTwenty')
        [S,i_s,js]=outerjoin(R,library20table,'MergeKeys',true) # add file 20
        clear('i_s','js','R','library20table')
        alltogether=S

    # sort matrix
    print(alltogether)
    sortedmatrix=sortmatrix(alltogether, posnum,matnum,average_negative)
    print(sortedmatrix)

    # write table
    print('creating matrix file')
    sortedmatrix.to_excel(outputfile2)
    #writetable(sortedmatrix,outputfile2)
    print('matrix file completed')
    #print(toc-tic, 'sec Elapsed')
    return

print(snakemake.input)
print(snakemake.output)
submit2_Callback(snakemake.config["posnum"], snakemake.config["matnum"], snakemake.output[0], snakemake.input, (snakemake.config["positive_reffile"]+snakemake.config["negative_reffile"]),snakemake.config["average_negative"])
