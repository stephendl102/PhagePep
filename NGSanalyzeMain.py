import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import math
from scipy import stats
from normalizelibrary import *
from sortmatrix import *
#from sortmatrix import *
startflank = 'GGTGGAGGT'
flank = 'TCT'
AA = 7

libnum = 3
posnum = 2
negnum = libnum-posnum
matnum = 100000

libfile = [None] * libnum
libfile[0] = r'/Users/mbp/Documents/MATLAB/Kelly/PHASTpep-master/PHASTpep-master/exampleFilesForPart2/PositiveLibraryExampleAforPHASTpep.xlsx'
libfile[1] = r'/Users/mbp/Documents/MATLAB/Kelly/PHASTpep-master/PHASTpep-master/exampleFilesForPart2/PositiveLibraryExampleBforPHASTpep.xlsx'
libfile[2] = r'/Users/mbp/Documents/MATLAB/Kelly/PHASTpep-master/PHASTpep-master/exampleFilesForPart2/NegativeLibraryExampleAforPHASTpep.xlsx'

reffile = [None] * libnum
reffile[0] = r'/Users/mbp/Documents/MATLAB/Kelly/PHASTpep-master/PHASTpep-master/exampleFilesForPart2/exampleReferencePHASTpep.xlsx'
reffile[1] = r'/Users/mbp/Documents/MATLAB/Kelly/PHASTpep-master/PHASTpep-master/exampleFilesForPart2/exampleReferencePHASTpep.xlsx'
reffile[2] = r'/Users/mbp/Documents/MATLAB/Kelly/PHASTpep-master/PHASTpep-master/exampleFilesForPart2/exampleReferencePHASTpep.xlsx'

outputfile2 = r'/Users/mbp/Documents/MATLAB/Kelly/530GUITEST.xlsx'

# Main PHASTpep file code
# created by Lindsey Brinton at the University of Virginia, 2015
# updated by Lindsey Brinton at the University of Virginia, 2016
# Translated to Python by Stephen Lees, 2022

def submit2_Callback(libnum, posnum, matnum, outputfile2,libfile,reffile):
    print('libnum is: ' + str(libnum))

    # file 1
    print('reading in set 1')
    library1table = normalizelibrary(libfile[0],reffile[0],'libOne')
    alltogether = library1table

    if libnum>1:
        # file 2
        print('reading in set 2')
        library2table = normalizelibrary(libfile[1],reffile[1],'libTwo')

        # combine file 1 and file 2
        print('combining sets 1 & 2')
        A = alltogether.merge(library2table,how='outer', on='Peptide', sort=True)
        alltogether=A

    # file 3
    if libnum>2:
        print('adding set 3')
        library3table = normalizelibrary(libfile[2],reffile[2],'libThree')
        B = alltogether.merge(library3table,how='outer', on='Peptide', sort=True) # add file 3
        alltogether=B

        # file 4
    if libnum>3:
        print('adding set 4')
        library4table = normalizelibrary(libfile[3],reffile[3],'libFour')
        C = alltogether.merge(library4table,how='outer', on='Peptide', sort=True) # add file 4
        alltogether=C

        # file 5
    if libnum>4:
        print('adding set 5')
        library5table = normalizelibrary(libfile[4],reffile[4],'libFive')
        D = alltogether.merge(library5table,how='outer', on='Peptide', sort=True) # add file 5
        alltogether=D

    # file 6
    if libnum>5:
        print('adding set 6')
        library6table = normalizelibrary(libfile[5],reffile[5],'libSix')
        E = alltogether.merge(library6table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=E

        # file 7
    if libnum>6:
        print('adding set 7')
        library7table = normalizelibrary(libfile[6],reffile[6],'libSeven')
        E = alltogether.merge(library7table,how='outer', on='Peptide', sort=True) # add file 6
        alltogether=F

    # file 8
    if libnum>7:
        print('adding set 8')
        libfile8 = getappdata(hMainGui,'libfile8')
        reffile8 = getappdata(hMainGui,'reffile8')
        library8table = normalizelibrary(libfile8,reffile8,'libEight')
        [G,ig,jg]=outerjoin(F,library8table,'MergeKeys',true) # add file 8
        clear('ig','jg','F','library8table')
        alltogether=G

    # file 9
    if libnum>8:
        print('adding set 9')
        libfile9 = getappdata(hMainGui,'libfile9')
        reffile9 = getappdata(hMainGui,'reffile9')
        library9table = normalizelibrary(libfile9,reffile9,'libNine')
        [H,ih,jh]=outerjoin(G,library9table,'MergeKeys',true) # add file 9
        clear('ih','jh','G','library9table')
        alltogether=H

    # file 10
    if libnum>9:
        print('adding set 10')
        libfile10 = getappdata(hMainGui,'libfile10')
        reffile10 = getappdata(hMainGui,'reffile10')
        library10table = normalizelibrary(libfile10,reffile10,'libTen')
        [I,ii,ji]=outerjoin(H,library10table,'MergeKeys',true) # add file 10
        clear('ii','ji','H','library10table')
        alltogether=I

        # file 11
    if libnum>10:
        print('adding set 11')
        libfile11 = getappdata(hMainGui,'libfile11')
        reffile11 = getappdata(hMainGui,'reffile11')
        library11table = normalizelibrary(libfile11,reffile11,'libEleven')
        [J,ij,jj]=outerjoin(I,library11table,'MergeKeys',true) # add file 11
        clear('ij','jj','I','library11table')
        alltogether=J

    # file 12
    if libnum>11:
        print('adding set 12')
        libfile12 = getappdata(hMainGui,'libfile12')
        reffile12 = getappdata(hMainGui,'reffile12')
        library12table = normalizelibrary(libfile12,reffile12,'libTweleve')
        [K,ik,jk]=outerjoin(J,library12table,'MergeKeys',true) # add file 12
        clear('ik','jk','J','library12table')
        alltogether=K

    # file 13
    if libnum>12:
        print('adding set 13')
        libfile13 = getappdata(hMainGui,'libfile13')
        reffile13 = getappdata(hMainGui,'reffile13')
        library13table = normalizelibrary(libfile13,reffile13,'libThirteen')
        [L,il,jl]=outerjoin(K,library13table,'MergeKeys',true) # add file 13
        clear('il','jl','K','library13table')
        alltogether=L

    # file 14
    if libnum>13:
        print('adding set 14')
        libfile14 = getappdata(hMainGui,'libfile14')
        reffile14 = getappdata(hMainGui,'reffile14')
        library14table = normalizelibrary(libfile14,reffile14,'libFourteen')
        [M,im,jm]=outerjoin(L,library14table,'MergeKeys',true) # add file 14
        clear('im','jm','L','library14table')
        alltogether=M

    # file 15
    if libnum>14:
        print('adding set 15')
        libfile15 = getappdata(hMainGui,'libfile15')
        reffile15 = getappdata(hMainGui,'reffile15')
        library15table = normalizelibrary(libfile15,reffile15,'libFifteen')
        [N, i_n, jn] = outerjoin(M,library15table,'MergeKeys',true) # add file 15
        clear('i_n','jn','M','library15table')
        alltogether=N

    # file 16
    if libnum>15:
        print('adding set 16')
        libfile16 = getappdata(hMainGui,'libfile16')
        reffile16 = getappdata(hMainGui,'reffile16')
        library16table = normalizelibrary(libfile16,reffile16,'libSixteen')
        [O,io,jo]=outerjoin(N,library16table,'MergeKeys',true) # add file 16
        clear('io','jo','N','library16table')
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
    sortedmatrix=sortmatrix(alltogether,libnum,posnum,matnum)
    print(sortedmatrix)

    # write table
    print('creating matrix file')
    sortedmatrix.to_excel(outputfile2)
    #writetable(sortedmatrix,outputfile2)
    print('matrix file completed')
    #print(toc-tic, 'sec Elapsed')
    return

submit2_Callback(libnum, posnum, matnum, outputfile2,libfile, reffile)
