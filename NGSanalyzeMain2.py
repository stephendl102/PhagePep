import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import math
from scipy import stats
from normalizelibrarySL import *
from sortmatrix import *
#from sortmatrix import *

# Main PHASTpep file code
# created by Lindsey Brinton at the University of Virginia, 2015
# updated by Lindsey Brinton at the University of Virginia, 2016
# Translated to Python by Stephen Lees, 2022

def submit2_Callback(libnum, posnum, matnum, outputfile2,libfile,reffile,average_negative=False):
    print('libnum is: ' + str(libnum))
    print(libfile[0])
    print([reffile[0]])

    # file 1
    print('reading in set 1')
    library1table = normalizelibrary(libfile[0],reffile[0],'libOne',average_negative)
    alltogether = library1table

    if libnum>1:
        # file 2
        print('reading in set 2')
        library2table = normalizelibrary(libfile[1],reffile[1],'libTwo',average_negative)

        # combine file 1 and file 2
        print('combining sets 1 & 2')
        A = alltogether.merge(library2table,how='outer', on='Peptide', sort=True)
        alltogether=A

    # file 3
    if libnum>2:
        print('adding set 3')
        library3table = normalizelibrary(libfile[2],reffile[2],'libThree',average_negative)
        B = alltogether.merge(library3table,how='outer', on='Peptide', sort=True) # add file 3
        alltogether=B

        # file 4
    if libnum>3:
        print('adding set 4')
        library4table = normalizelibrary(libfile[3],reffile[3],'libFour',average_negative)
        C = alltogether.merge(library4table,how='outer', on='Peptide', sort=True) # add file 4
        alltogether=C

        # file 5
    if libnum>4:
        print('adding set 5')
        library5table = normalizelibrary(libfile[4],reffile[4],'libFive',average_negative)
        D = alltogether.merge(library5table,how='outer', on='Peptide', sort=True) # add file 5
        alltogether=D

    # file 6
    if libnum>5:
        print('adding set 6')
        library6table = normalizelibrary(libfile[5],reffile[5],'libSix',average_negative)
        E = alltogether.merge(library6table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=E

        # file 7
    if libnum>6:
        print('adding set 7')
        library7table = normalizelibrary(libfile[6],reffile[6],'libSeven',average_negative)
        F = alltogether.merge(library7table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=F

    # file 8
    if libnum>7:
        print('adding set 8')
        library8table = normalizelibrary(libfile[7],reffile[7],'libEight',average_negative)
        G = alltogether.merge(library8table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=G

    # file 9
    if libnum>8:
        print('adding set 9')
        library9table = normalizelibrary(libfile[8],reffile[8],'libNine',average_negative)
        H = alltogether.merge(library9table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=H

    # file 10
    if libnum>9:
        print('adding set 10')
        library10table = normalizelibrary(libfile[9],reffile[9],'libTen',average_negative)
        I = alltogether.merge(library10table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=I

        # file 11
    if libnum>10:
        print('adding set 11')
        library11table = normalizelibrary(libfile[10],reffile[10],'libEleven',average_negative)
        J = alltogether.merge(library11table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=J

    # file 12
    if libnum>11:
        print('adding set 12')
        library12table = normalizelibrary(libfile[11],reffile[11],'libTwelve',average_negative)
        K = alltogether.merge(library12table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=K

    # file 13
    if libnum>12:
        print('adding set 13')
        library13table = normalizelibrary(libfile[12],reffile[12],'libThirteen',average_negative)
        L = alltogether.merge(library13table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=L


    # file 14
    if libnum>13:
        print('adding set 14')
        library14table = normalizelibrary(libfile[13],reffile[13],'libFourteen',average_negative)
        M = alltogether.merge(library14table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=M

    # file 15
    if libnum>14:
        print('adding set 15')
        library15table = normalizelibrary(libfile[14],reffile[14],'libFifteen',average_negative)
        N = alltogether.merge(library15table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=N

    # file 16
    if libnum>15:
        print('adding set 16')
        library16table = normalizelibrary(libfile[15],reffile[15],'libSixteen',average_negative)
        O = alltogether.merge(library16table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=O

    # file 17
    if libnum>16:
        print('adding set 17')
        libfile17 = getappdata(hMainGui,'libfile17')
        reffile17 = getappdata(hMainGui,'reffile17')
        library17table = normalizelibrary(libfile17,reffile17,'libSeventeen')
        [P,ip,jp]=outerjoin(O,library17table,'MergeKeys',true) # add file 17
        clear('ip','jp','O','library17table')
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
    sortedmatrix=sortmatrix(alltogether,libnum,posnum,matnum,average_negative)
    print(sortedmatrix)

    # write table
    print('creating matrix file')
    sortedmatrix.to_excel(outputfile2)
    #writetable(sortedmatrix,outputfile2)
    print('matrix file completed')
    #print(toc-tic, 'sec Elapsed')
    return

print(snakemake.input)
submit2_Callback(snakemake.config["libnum"], snakemake.config["posnum"], snakemake.config["matnum"], snakemake.output[0], snakemake.input, (snakemake.config["positive_reffile"]+snakemake.config["negative_reffile"]),snakemake.config["average_negative"])
