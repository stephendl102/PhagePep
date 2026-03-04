import matplotlib.pyplot as plt
import time
import pandas as pd
#from sortmatrix import *
startflank = 'GGTGGAGGT'
endflank = 'TCT'
mer = 7 # amino acids in peptide
libnum = 2 # total number of libraries not including references
posnum = 1 # total number of positive libraries

# Main PHASTpep file code
# created by Lindsey Brinton at the University of Virginia, 2015
# updated by Lindsey Brinton at the University of Virginia, 2016
# Translated to Python by Stephen Lees, 2022

def NGSanalyzeMain(varargin):
    """
    NGSANALYZEMAIN MATLAB code for NGSanalyzeMain.fig
     NGSANALYZEMAIN, by itself, creates a new NGSANALYZEMAIN or raises the existing
      singleton*.

      H = NGSANALYZEMAIN returns the handle to a new NGSANALYZEMAIN or the handle to
      the existing singleton*.

      NGSANALYZEMAIN('CALLBACK',hObject,eventData,handles,...) calls the local
      function named CALLBACK in NGSANALYZEMAIN.M with the given input arguments.

      NGSANALYZEMAIN('Property','Value',...) creates a new NGSANALYZEMAIN or raises the
      existing singleton*.  Starting from the left, property value pairs are
      applied to the GUI before NGSanalyzeMain_OpeningFcn gets called.  An
      unrecognized property name or invalid value makes property application
      stop.  All inputs are passed to NGSanalyzeMain_OpeningFcn via varargin.

      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
      instance to run (singleton)".

      See also: GUIDE, GUIDATA, GUIHANDLES

      Edit the above text to modify the response to help NGSanalyzeMain

      Last Modified by GUIDE v2.5 06-Aug-2015 18:01:37

      Begin initialization code - DO NOT EDIT
      """
    gui_Singleton = 1
    # 33-47? ***
    varargout = 1

    return varargout

# End initialization code - DO NOT EDIT

# --- Executes just before NGSanalyzeMain is made visible.
def NGSanalyzeMain_OpeningFcn(hObject, eventdata, handles, varargin):
    """
    This function has no output args, see OutputFcn.
    hObject    handle to figure
    eventdata  reserved - to be defined in a future version of MATLAB
    handles    structure with handles and user data (see GUIDATA)
    arargin   command line arguments to NGSanalyzeMain (see VARARGIN)

    Choose default command line output for NGSanalyzeMain
    """
    handles.output = hObject

    # Update handles structure
    guidata(hObject, handles)

    setappdata(0,'hMainGui',gcf); # get current figure

    # UIWAIT makes NGSanalyzeMain wait for user response (see UIRESUME)
    # uiwait(handles.figure1);
    return

# --- Outputs from this function are returned to the command line.
def NGSanalyzeMain_OutputFcn(hObject, eventdata, handles):
    """
    varargout  cell array for returning output args (see VARARGOUT);
    hObject    handle to figure
    eventdata  reserved - to be defined in a future version of MATLAB
    handles    structure with handles and user data (see GUIDATA)

    Get default command line output from handles structure
    """
    varargout[1] = handles.output
    return varargout

# --- Executes on button press in submit1.
def submit1_Callback(hObject, eventdata, handles):
    # access inputs
    hMainGui = getappdata(0,'hMainGui');
    mer = getappdata(hMainGui,'mer');
    phdcheck = getappdata(hMainGui,'phdcheck');
    othercheck = getappdata(hMainGui,'othercheck');
    fileinput1 = getappdata(hMainGui,'fileinput1');
    outputfile1 = getappdata(hMainGui,'outputfile1');

    println(startflank);
    println(endflank);

    # start timer
    tic = time.time()

    #call function to read in .fastq file & export translated matrix to excel
    toexcelfile = translatefastq(mer,fileinput1,startflank,endflank,outputfile1,phdcheck); # toexcelfile is what was exported to excel

    # read timer
    toc = time.time()
    println(toc-tic, 'sec Elapsed')
    # hObject    handle to submit1 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    return

# --- Executes on button press in submit2.
def submit2_Callback(hObject, eventdata, handles):
    tic = time.time()
    # access inputs
    hMainGui = getappdata(0,'hMainGui');
    libnum = getappdata(hMainGui,'libnum');
    println(['libnum is: ',num2str(libnum)]);
    posnum = getappdata(hMainGui,'posnum');
    matnum = getappdata(hMainGui,'matnum');
    outputfile2 = getappdata(hMainGui,'outputfile2');

    # file 1
    println('reading in set 1');
    libfile1 = getappdata(hMainGui,'libfile1');
    reffile1 = getappdata(hMainGui,'reffile1');
    library1table = normalizelibrary(libfile1,reffile1,'libOne');

    # file 2
    println('reading in set 2');
    libfile2 = getappdata(hMainGui,'libfile2');
    reffile2 = getappdata(hMainGui,'reffile2');
    library2table = normalizelibrary(libfile2,reffile2,'libTwo');

    # combine file 1 and file 2
    println('combining sets 1 & 2');
    [A,ia,ja]=outerjoin(library1table,library2table,'MergeKeys',true);
    clear('ia','ja','library1table','library2table');
    alltogether=A;

    # file 3
    if libnum>2:
        println('adding set 3');
        libfile3 = getappdata(hMainGui,'libfile3');
        reffile3 = getappdata(hMainGui,'reffile3');
        library3table = normalizelibrary(libfile3,reffile3,'libThree');
        [B,ib,jb]=outerjoin(A,library3table,'MergeKeys',true); # add file 3
        clear('ib','jb','A','library3table');
        alltogether=B;

        # file 4
    if libnum>3:
        prinln('adding set 4');
        libfile4 = getappdata(hMainGui,'libfile4');
        reffile4 = getappdata(hMainGui,'reffile4');
        library4table = normalizelibrary(libfile4,reffile4,'libFour');
        [C,ic,jc]=outerjoin(B,library4table,'MergeKeys',true); # add file 4
        clear('ic','jc','B','library4table');
        alltogether=C;

        # file 5
    if libnum>4:
        println('adding set 5');
        libfile5 = getappdata(hMainGui,'libfile5');
        reffile5 = getappdata(hMainGui,'reffile5');
        library5table = normalizelibrary(libfile5,reffile5,'libFive');
        [D,id,jd]=outerjoin(C,library5table,'MergeKeys',true); # add file 5
        clear('id','jd','C','library5table');
        alltogether=D;

    # file 6
    if libnum>5:
        println('adding set 6');
        libfile6 = getappdata(hMainGui,'libfile6');
        reffile6 = getappdata(hMainGui,'reffile6');
        library6table = normalizelibrary(libfile6,reffile6,'libSix');
        [E,ie,je]=outerjoin(D,library6table,'MergeKeys',true); # add file 6
        clear('ie','je','D','library6table');
        alltogether=E;

        # file 7
    if libnum>6:
        println('adding set 7');
        libfile7 = getappdata(hMainGui,'libfile7');
        reffile7 = getappdata(hMainGui,'reffile7');
        library7table = normalizelibrary(libfile7,reffile7,'libSeven');
        [F,iff,jf]=outerjoin(E,library7table,'MergeKeys',true); # add file 7
        clear('iff','jf','E','library7table');
        alltogether=F;

    # file 8
    if libnum>7:
        println('adding set 8');
        libfile8 = getappdata(hMainGui,'libfile8');
        reffile8 = getappdata(hMainGui,'reffile8');
        library8table = normalizelibrary(libfile8,reffile8,'libEight');
        [G,ig,jg]=outerjoin(F,library8table,'MergeKeys',true); # add file 8
        clear('ig','jg','F','library8table');
        alltogether=G;

    # file 9
    if libnum>8:
        println('adding set 9');
        libfile9 = getappdata(hMainGui,'libfile9');
        reffile9 = getappdata(hMainGui,'reffile9');
        library9table = normalizelibrary(libfile9,reffile9,'libNine');
        [H,ih,jh]=outerjoin(G,library9table,'MergeKeys',true); # add file 9
        clear('ih','jh','G','library9table');
        alltogether=H;

    # file 10
    if libnum>9:
        println('adding set 10');
        libfile10 = getappdata(hMainGui,'libfile10');
        reffile10 = getappdata(hMainGui,'reffile10');
        library10table = normalizelibrary(libfile10,reffile10,'libTen');
        [I,ii,ji]=outerjoin(H,library10table,'MergeKeys',true); # add file 10
        clear('ii','ji','H','library10table')
        alltogether=I;

        # file 11
    if libnum>10:
        println('adding set 11');
        libfile11 = getappdata(hMainGui,'libfile11');
        reffile11 = getappdata(hMainGui,'reffile11');
        library11table = normalizelibrary(libfile11,reffile11,'libEleven');
        [J,ij,jj]=outerjoin(I,library11table,'MergeKeys',true); # add file 11
        clear('ij','jj','I','library11table')
        alltogether=J;

    # file 12
    if libnum>11:
        println('adding set 12');
        libfile12 = getappdata(hMainGui,'libfile12');
        reffile12 = getappdata(hMainGui,'reffile12');
        library12table = normalizelibrary(libfile12,reffile12,'libTweleve');
        [K,ik,jk]=outerjoin(J,library12table,'MergeKeys',true); # add file 12
        clear('ik','jk','J','library12table')
        alltogether=K;

    # file 13
    if libnum>12:
        println('adding set 13');
        libfile13 = getappdata(hMainGui,'libfile13');
        reffile13 = getappdata(hMainGui,'reffile13');
        library13table = normalizelibrary(libfile13,reffile13,'libThirteen');
        [L,il,jl]=outerjoin(K,library13table,'MergeKeys',true); # add file 13
        clear('il','jl','K','library13table')
        alltogether=L;

    # file 14
    if libnum>13:
        println('adding set 14');
        libfile14 = getappdata(hMainGui,'libfile14');
        reffile14 = getappdata(hMainGui,'reffile14');
        library14table = normalizelibrary(libfile14,reffile14,'libFourteen');
        [M,im,jm]=outerjoin(L,library14table,'MergeKeys',true); # add file 14
        clear('im','jm','L','library14table')
        alltogether=M;

    # file 15
    if libnum>14:
        println('adding set 15');
        libfile15 = getappdata(hMainGui,'libfile15');
        reffile15 = getappdata(hMainGui,'reffile15');
        library15table = normalizelibrary(libfile15,reffile15,'libFifteen');
        [N, i_n, jn] = outerjoin(M,library15table,'MergeKeys',true); # add file 15
        clear('i_n','jn','M','library15table')
        alltogether=N;

    # file 16
    if libnum>15:
        println('adding set 16');
        libfile16 = getappdata(hMainGui,'libfile16');
        reffile16 = getappdata(hMainGui,'reffile16');
        library16table = normalizelibrary(libfile16,reffile16,'libSixteen');
        [O,io,jo]=outerjoin(N,library16table,'MergeKeys',true); # add file 16
        clear('io','jo','N','library16table')
        alltogether=O;

    # file 17
    if libnum>16:
        println('adding set 17');
        libfile17 = getappdata(hMainGui,'libfile17');
        reffile17 = getappdata(hMainGui,'reffile17');
        library17table = normalizelibrary(libfile17,reffile17,'libSeventeen');
        [P,ip,jp]=outerjoin(O,library17table,'MergeKeys',true); # add file 17
        clear('ip','jp','O','library17table')
        alltogether=P;

    # file 18
    if libnum>17:
        println('adding set 18');
        libfile18 = getappdata(hMainGui,'libfile18');
        reffile18 = getappdata(hMainGui,'reffile18');
        library18table = normalizelibrary(libfile18,reffile18,'libEighteen');
        [Q,iq,jq]=outerjoin(P,library18table,'MergeKeys',true); # add file 18
        clear('iq','jq','P','library18table')
        alltogether=Q;

    # file 19
    if libnum>18:
        println('adding set 19');
        libfile19 = getappdata(hMainGui,'libfile19');
        reffile19 = getappdata(hMainGui,'reffile19');
        library19table = normalizelibrary(libfile19,reffile19,'libNineteen');
        [R,ir,jr]=outerjoin(Q,library19table,'MergeKeys',true); # add file 19
        clear('ir','jr','Q','library19table')
        alltogether=R;

    # file 20
    if libnum>19:
        println('adding set 20');
        libfile20 = getappdata(hMainGui,'libfile20');
        reffile20 = getappdata(hMainGui,'reffile20');
        library20table = normalizelibrary(libfile20,reffile20,'libTwenty');
        [S,i_s,js]=outerjoin(R,library20table,'MergeKeys',true); # add file 20
        clear('i_s','js','R','library20table')
        alltogether=S;

    # sort matrix
    sortedmatrix=sortmatrix(alltogether,libnum,posnum,matnum);

    # write table
    println('creating matrix file');
    writetable(sortedmatrix,outputfile2);
    #writetable(sortedmatrix,outputfile2);
    println('matrix file completed');
    toc = time.time()
    println(toc-tic, 'sec Elapsed')

    # hObject    handle to submit2 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    return

def outputfile2_Callback(hObject, eventdata, handles):
    """ hObject    handle to outputfile2 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)

    # Hints: get(hObject,'String') returns contents of outputfile2 as text
    #        str2double(get(hObject,'String')) returns contents of outputfile2 as a double
    """
    hMainGui = getappdata(0,'hMainGui');
    outputfile2 = get(hObject,'String');
    setappdata(hMainGui,'outputfile2',outputfile2);
    return

# --- Executes during object creation, after setting all properties.
def outputfile2_CreateFcn(hObject, eventdata, handles):
    """
    # hObject    handle to outputfile2 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    empty - handles not created until after all CreateFcns called

    # Hint: edit controls usually have a white background on Windows.
    #       See ISPC and COMPUTER.
    """
    if ispc and isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor')):
        set(hObject,'BackgroundColor','white');
    return

# --- Executes on button press in browse.
def browse_Callback(hObject, eventdata, handles):
    # hObject    handle to browse (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    [browsename,browsepath,browseindex] = uigetfile('*.fastq');  # browse to find .fastq file
    inputfile1name=[browsepath,browsename];
    println(['file selected: ', inputfile1name]);
    setappdata(hMainGui,'inputfile1name',inputfile1name);
    if length(inputfile1name)>2:
        println('File selected using browse');
    inputfile1_Callback(hObject, eventdata, handles);
    return

# --- Executes on button press in phd.
def phd_Callback(hObject, eventdata, handles):
    # hObject    handle to phd (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    phdcheck = get(hObject,'Value');
    setappdata(hMainGui,'phdcheck',phdcheck);
    # Hint: get(hObject,'Value') returns toggle state of phd
    return

    # --- Executes on button press in other.
def other_Callback(hObject, eventdata, handles):
    # hObject    handle to other (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    othercheck = get(hObject,'Value');
    setappdata(hMainGui,'othercheck',othercheck);
    # Hint: get(hObject,'Value') returns toggle state of other
    return

def start_Callback(hObject, eventdata, handles):
    # hObject    handle to start (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    startflank = get(hObject,'String');
    setappdata(hMainGui,'startflank',startflank);
    # Hints: get(hObject,'String') returns contents of start as text
    #        str2double(get(hObject,'String')) returns contents of start as a double
    return

    # --- Executes during object creation, after setting all properties.
def start_CreateFcn(hObject, eventdata, handles):
    # hObject    handle to start (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    empty - handles not created until after all CreateFcns called

    # Hint: edit controls usually have a white background on Windows.
    #       See ISPC and COMPUTER.
    if ispc and isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor')):
        set(hObject,'BackgroundColor','white');
    return

def finish_Callback(hObject, eventdata, handles):
    # hObject    handle to finish (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    finishflank = get(hObject,'String');
    setappdata(hMainGui,'finishflank',finishflank);
    # Hints: get(hObject,'String') returns contents of finish as text
    #        str2double(get(hObject,'String')) returns contents of finish as a double
    return

# --- Executes during object creation, after setting all properties.
def finish_CreateFcn(hObject, eventdata, handles):
    # hObject    handle to finish (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    empty - handles not created until after all CreateFcns called

    # Hint: edit controls usually have a white background on Windows.
    #       See ISPC and COMPUTER.
    if ispc and isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor')):
        set(hObject,'BackgroundColor','white');
    return

def mer_Callback(hObject, eventdata, handles):
    # hObject    handle to mer (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    mer = str2double(get(hObject,'String'));
    setappdata(hMainGui,'mer',mer);
    # Hints: get(hObject,'String') returns contents of mer as text
    #        str2double(get(hObject,'String')) returns contents of mer as a double
    return


# --- Executes during object creation, after setting all properties.
def mer_CreateFcn(hObject, eventdata, handles):
    # hObject    handle to mer (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    empty - handles not created until after all CreateFcns called

    # Hint: edit controls usually have a white background on Windows.
    #       See ISPC and COMPUTER.
    if ispc and isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor')):
        set(hObject,'BackgroundColor','white');
    return



def inputfile1_Callback(hObject, eventdata, handles):
    # hObject    handle to inputfile1 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    fileinput1 = get(hObject,'String');
    setappdata(hMainGui,'fileinput1',fileinput1);

    fileinputed1 = getappdata(hMainGui,'inputfile1name'); # replace text after select something in browse

    if length(fileinputed1)>2:
        set(handles.inputfile1,'String',fileinputed1);
        setappdata(hMainGui,'fileinput1',fileinputed1);
        # Hints: get(hObject,'String') returns contents of inputfile1 as text
        #        str2double(get(hObject,'String')) returns contents of inputfile1 as a double
    return


# --- Executes during object creation, after setting all properties.
def inputfile1_CreateFcn(hObject, eventdata, handles):
    # hObject    handle to inputfile1 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    empty - handles not created until after all CreateFcns called

    # Hint: edit controls usually have a white background on Windows.
    #       See ISPC and COMPUTER.
    if ispc and isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor')):
        set(hObject,'BackgroundColor','white');
    return


def lib1_Callback(hObject, eventdata, handles):
    # hObject    handle to lib1 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile1 = get(hObject,'String');
    setappdata(hMainGui,'libfile1',libfile1);

    libfiled1 = getappdata(hMainGui,'libfile1name'); # replace text after select something in browse
    if length(libfiled1)>2:
        set(handles.lib1,'String',libfiled1);
        setappdata(hMainGui,'libfile1',libfiled1);
    return

    # Hints: get(hObject,'String') returns contents of lib1 as text
    #        str2double(get(hObject,'String')) returns contents of lib1 as a double

def lib2_Callback(hObject, eventdata, handles):
    # hObject    handle to lib2 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile2 = get(hObject,'String');
    setappdata(hMainGui,'libfile2',libfile2);

    libfiled2 = getappdata(hMainGui,'libfile2name'); # replace text after select something in browse
    if length(libfiled2)>2:
        set(handles.lib2,'String',libfiled2);
        setappdata(hMainGui,'libfile2',libfiled2);
    # Hints: get(hObject,'String') returns contents of lib2 as text
    #        str2double(get(hObject,'String')) returns contents of lib2 as a double
    return


def lib3_Callback(hObject, eventdata, handles):
    # hObject    handle to lib3 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile3 = get(hObject,'String');
    setappdata(hMainGui,'libfile3',libfile3);

    libfiled3 = getappdata(hMainGui,'libfile3name'); # replace text after select something in browse
    if length(libfiled3)>2:
        set(handles.lib3,'String',libfiled3);
        setappdata(hMainGui,'libfile3',libfiled3);
    # Hints: get(hObject,'String') returns contents of lib3 as text
    #        str2double(get(hObject,'String')) returns contents of lib3 as a double


def lib4_Callback(hObject, eventdata, handles):
    # hObject    handle to lib4 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile4 = get(hObject,'String');
    setappdata(hMainGui,'libfile4',libfile4);

    libfiled4 = getappdata(hMainGui,'libfile4name'); # replace text after select something in browse
    if length(libfiled4)>2:
        set(handles.lib4,'String',libfiled4);
        setappdata(hMainGui,'libfile4',libfiled4);
    # Hints: get(hObject,'String') returns contents of lib4 as text
    #        str2double(get(hObject,'String')) returns contents of lib4 as a double
    return


def lib5_Callback(hObject, eventdata, handles):
    # hObject    handle to lib5 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile5 = get(hObject,'String');
    setappdata(hMainGui,'libfile5',libfile5);

    libfiled5 = getappdata(hMainGui,'libfile5name'); # replace text after select something in browse
    if length(libfiled5)>2:
        set(handles.lib5,'String',libfiled5);
        setappdata(hMainGui,'libfile5',libfiled5);
    # Hints: get(hObject,'String') returns contents of lib5 as text
    #        str2double(get(hObject,'String')) returns contents of lib5 as a double
    return

def lib6_Callback(hObject, eventdata, handles):
    # hObject    handle to lib6 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile6 = get(hObject,'String');
    setappdata(hMainGui,'libfile6',libfile6);

    libfiled6 = getappdata(hMainGui,'libfile6name'); # replace text after select something in browse
    if length(libfiled6)>2:
        set(handles.lib6,'String',libfiled6);
        setappdata(hMainGui,'libfile6',libfiled6);
    # Hints: get(hObject,'String') returns contents of lib6 as text
    #        str2double(get(hObject,'String')) returns contents of lib6 as a double
    return


def lib7_Callback(hObject, eventdata, handles):
    # hObject    handle to lib7 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile7 = get(hObject,'String');
    setappdata(hMainGui,'libfile7',libfile7);

    libfiled7 = getappdata(hMainGui,'libfile7name'); # replace text after select something in browse
    if length(libfiled7)>2:
        set(handles.lib7,'String',libfiled7);
        setappdata(hMainGui,'libfile7',libfiled7);
    # Hints: get(hObject,'String') returns contents of lib7 as text
    #        str2double(get(hObject,'String')) returns contents of lib7 as a double
    return


def lib8_Callback(hObject, eventdata, handles):
    # hObject    handle to lib8 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile8 = get(hObject,'String');
    setappdata(hMainGui,'libfile8',libfile8);

    libfiled8 = getappdata(hMainGui,'libfile8name'); # replace text after select something in browse
    if length(libfiled8)>2:
        set(handles.lib8,'String',libfiled8);
        setappdata(hMainGui,'libfile8',libfiled8);
    # Hints: get(hObject,'String') returns contents of lib8 as text
    #        str2double(get(hObject,'String')) returns contents of lib8 as a double

def lib9_Callback(hObject, eventdata, handles):
    # hObject    handle to lib9 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile9 = get(hObject,'String');
    setappdata(hMainGui,'libfile9',libfile9);

    libfiled9 = getappdata(hMainGui,'libfile9name'); # replace text after select something in browse
    if length(libfiled9)>2:
        set(handles.lib9,'String',libfiled9);
        setappdata(hMainGui,'libfile9',libfiled9);
    # Hints: get(hObject,'String') returns contents of lib9 as text
    #        str2double(get(hObject,'String')) returns contents of lib9 as a double
    return

def lib10_Callback(hObject, eventdata, handles):
    # hObject    handle to lib10 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile10 = get(hObject,'String');
    setappdata(hMainGui,'libfile10',libfile10);

    libfiled10 = getappdata(hMainGui,'libfile10name'); # replace text after select something in browse
    if length(libfiled10)>2:
        set(handles.lib10,'String',libfiled10);
        setappdata(hMainGui,'libfile10',libfiled10);
    # Hints: get(hObject,'String') returns contents of lib10 as text
    #        str2double(get(hObject,'String')) returns contents of lib10 as a double
    return


def lib11_Callback(hObject, eventdata, handles):
    # hObject    handle to lib11 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile11 = get(hObject,'String');
    setappdata(hMainGui,'libfile11',libfile11);

    libfiled11 = getappdata(hMainGui,'libfile11name'); # replace text after select something in browse
    if length(libfiled11)>2:
        set(handles.lib11,'String',libfiled11);
        setappdata(hMainGui,'libfile11',libfiled11);
    # Hints: get(hObject,'String') returns contents of lib11 as text
    #        str2double(get(hObject,'String')) returns contents of lib11 as a double
    return

def lib12_Callback(hObject, eventdata, handles):
    # hObject    handle to lib12 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile12 = get(hObject,'String');
    setappdata(hMainGui,'libfile12',libfile12);

    libfiled12 = getappdata(hMainGui,'libfile12name'); # replace text after select something in browse
    if length(libfiled12)>2:
        set(handles.lib12,'String',libfiled12);
        setappdata(hMainGui,'libfile12',libfiled12);
    # Hints: get(hObject,'String') returns contents of lib12 as text
    #        str2double(get(hObject,'String')) returns contents of lib12 as a double
    return

def lib13_Callback(hObject, eventdata, handles):
    # hObject    handle to lib13 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile13 = get(hObject,'String');
    setappdata(hMainGui,'libfile13',libfile13);

    libfiled13 = getappdata(hMainGui,'libfile13name'); # replace text after select something in browse
    if length(libfiled13)>2:
        set(handles.lib13,'String',libfiled13);
        setappdata(hMainGui,'libfile13',libfiled13);
    # Hints: get(hObject,'String') returns contents of lib13 as text
    #        str2double(get(hObject,'String')) returns contents of lib13 as a double
    return

def lib14_Callback(hObject, eventdata, handles):
    # hObject    handle to lib14 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile14 = get(hObject,'String');
    setappdata(hMainGui,'libfile14',libfile14);

    libfiled14 = getappdata(hMainGui,'libfile14name'); # replace text after select something in browse
    if length(libfiled14)>2:
        set(handles.lib14,'String',libfiled14);
        setappdata(hMainGui,'libfile1',libfiled14);
    # Hints: get(hObject,'String') returns contents of lib14 as text
    #        str2double(get(hObject,'String')) returns contents of lib14 as a double
    return

def lib15_Callback(hObject, eventdata, handles):
    # hObject    handle to lib15 (see GCBO)
    # eventdata  reserved - to be defined in a future version of MATLAB
    # handles    structure with handles and user data (see GUIDATA)
    hMainGui = getappdata(0,'hMainGui');
    libfile15 = get(hObject,'String');
    setappdata(hMainGui,'libfile15',libfile15);

    libfiled15 = getappdata(hMainGui,'libfile15name'); # replace text after select something in browse
    if length(libfiled15)>2:
        set(handles.lib15,'String',libfiled15);
        setappdata(hMainGui,'libfile15',libfiled15);
    # Hints: get(hObject,'String') returns contents of lib15 as text
    #        str2double(get(hObject,'String')) returns contents of lib15 as a double
    return
