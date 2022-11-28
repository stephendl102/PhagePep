""" function to combine multiple libraries together such that their
# duplicated sequences are combined into a single row
#
# inputs:
#       matrix: table of all libraries
#
#       libnumber: total number of libraries, not including reference
#       libraries
#
#       posnumber: number of libraries corresponding to positive screens
#
#       matrixnum: number of sequences to be included in final matrix
#
# outputs:
#       librarymatrix: cell array of sorted sequences
#
# Created by Lindsey Brinton at the University of Virginia, 2015
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import sklearn
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from normalizelibrary import *

filename = "../ReferenceLibraries/PhastpepPaper8_8_22_NaiveCharacters.xlsx"
posnum = 4
libnum = 7

norm = plt.Normalize(1,4)

def visualize(filename,libnum,posnum):
    print('File is: ' + str(filename))
    df = pd.read_excel(filename,index_col=0)
    df = df.drop(df[(df['libOne'] == df.iloc[-1,1]) & (df['libTwo'] == df.iloc[-1,2]) & (df['libThree'] == df.iloc[-1,3])&(df['libFour'] == df.iloc[-1,4])].index) #&(df['libFive'] == df.iloc[-1,5])&(df['libSix'] == df.iloc[-1,6])&(df['libSeven'] == df.iloc[-1,7])&(df['libEight'] == df.iloc[-1,8])&(df['libNine'] == df.iloc[-1,9])&(df['libTen'] == df.iloc[-1,10])&(df['libEleven'] == df.iloc[-1,11])&(df['libTwelve'] == df.iloc[-1,12])
    print(df.loc[0:6])
    df['y'] = df.iloc[:,1:posnum+1].mean(1,numeric_only=True)
    df['x'] = df.iloc[:,posnum+1:libnum+1].mean(1,numeric_only=True)
    global names
    names = df['Peptide'].to_numpy()
    print(names)

    a = df.iloc[:,1:posnum+1]
    if posnum>1:
        PosNorm = a.std(1).div(df.y, axis=0)
    else:
        PosNorm = np.ones(len(df['x']))
        print("Only one positive library, no CV")
    b = df.iloc[:,1:posnum+1:libnum+1]
    NegNorm = b.std(1).div(df.x, axis=0)

    global ax
    global fig
    fig, ax = plt.subplots(1,2)
    ax[0].scatter(df['x'],df['y'])
    ax[0].set_xlabel('Streptavidin-negative Score Average')
    ax[0].set_ylabel('Streptavidin-positive Score Average')
    ax[0].set_title('Screen Normalized Peptide Abundance')

    print(df['Score'])
    global sc
    sc = ax[1].scatter(PosNorm,df['Score'])
    ax[1].set_xlabel('Coefficient of Variation (CV)')
    ax[1].set_ylabel('Score')
    ax[1].set_title('Score vs CV')

    global annot
    annot = ax[1].annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    fig.canvas.mpl_connect("motion_notify_event", hover)
    plt.show()

def update_annot(ind):

    pos = sc.get_offsets()[ind["ind"][0]]
    annot.xy = pos
    text = "{}, {}".format(" ".join(list(map(str,ind["ind"]))),
                           " ".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_alpha(0.4)


def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax[1]:
        cont, ind = sc.contains(event)
        if cont:
            update_annot(ind)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()

def Correlation(filename):
    print('File is: ' + str(filename))
    df = pd.read_excel(filename,index_col=0,header=0)
    print(df)
    y = df.iloc[:,1].values
    y = stats.zscore(y)
    print(np.mean(y))
    print(y)
    X = np.absolute(df.iloc[:,2:].values)
    X = StandardScaler().fit_transform(X)
    print(X)

    pls2 = PLSRegression(n_components=2)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(X)
    PLSR = pls2.fit_transform(X, y)
    Y_pred = pls2.predict(X)
    plt.scatter(Y_pred,y)
    principalDF = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])


    print(pca.explained_variance_ratio_)

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('2 component PCA', fontsize = 20)
    print(y)

    ax.scatter(principalDF.iloc[:,0], principalDF.iloc[:,1], c = y, cmap = "gray", s = 3)
    ax.grid()


    absCharge = np.absolute(df['gravy'].values)
    Frequency = df['Instability'].values
    print(scipy.stats.spearmanr(absCharge, Frequency))
    df.plot.scatter(x='gravy', y='Frequency', c='DarkBlue')
    plt.show()



#visualize(filename,libnum,posnum)
Correlation(filename)






#visualize(snakemake.input,snakemake.config["posnum"])

#librarymatrix = sortmatrix(matrix,libnumber,posnumber,matrixnum)
#librarymatrix.to_csv(r'/Users/mbp/Desktop/librarymatrix530.csv')
