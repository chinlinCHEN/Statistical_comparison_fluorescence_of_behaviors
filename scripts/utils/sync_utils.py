import os
import sys
import numpy as np
import math
from scipy import interpolate
import pandas as pd
from multiprocessing import Process
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from skimage import io
import itertools
from itertools import groupby 
import more_itertools as mit
import pickle
import h5py
import utils2p
import utils2p.synchronization

from utils import plot_setting
import utils.math_utils as math_utils




freqTS = 30000.0
samplingFact = 1500
windowOpflowSampling = samplingFact/5
#windowOpflowSampling = 1

#### GAIN AFTER NEW GOOD CALIBRATION 02/16
# gain0X = round(1/0.70,2)
# gain0Y = round(1/0.68,2)
# gain1X = round(1/0.68,2)
# gain1Y = round(1/0.66,2)


#### GAIN AFTER NEW GOOD CALIBRATION 06/06
gain0X = round(1/1.52,2)
gain0Y = round(1/1.48,2)
gain1X = round(1/1.48,2)
gain1Y = round(1/1.44,2)

#opflow conversion factor
rot_to_mm_factor = 2*math.pi*5  #10mm ball
rot_to_deg_factor = 360

## Optic flow bout thrshold

#Smooth window =1 => AP_thrsld=0.12 mm, ML_thrsld=0.12 mm, Yaw_thrsld=5 deg
AP_thrsld=0.12 #mm (based on smooth opflow trace)
ML_thrsld=0.12 #mm
Yaw_thrsld=5 #deg.
window_for_smth=1 #second

duration_for_replace=0.2
duration_for_scoring=0.2


CO2puff_color='k'
GC_color='forestgreen'
AP_color='r'
ML_color='b'
Yaw_color='m'









def convert_behList_to_binarizedList(Beh_list, beh_class):



    print('len(Beh_list)', len(Beh_list))

    new_beh_bin_list=[0]*len(Beh_list)

    beh_idx=np.where(np.asarray(Beh_list)==beh_class)[0]




    for idx in beh_idx:
        new_beh_bin_list[idx]=1

    # print('new_beh_bin_list', new_beh_bin_list)
    # print('len new_beh_bin_list', len(new_beh_bin_list))



    return new_beh_bin_list





def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]




def replace_nan_with_interp(dataset_2d_ori):

    dataset_2d_ori_c=dataset_2d_ori.copy()

    new_dataset_2d=[]

    for roi_i, gc_trace in enumerate(dataset_2d_ori_c):

        nans, x= nan_helper(gc_trace)
        print('x', x)
        gc_trace=np.asarray(gc_trace)
        gc_trace[nans]= np.interp(x(nans), x(~nans), gc_trace[~nans])

        new_dataset_2d.append(gc_trace)


    return new_dataset_2d


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx


def interpBetweenSeconds(diffList,maxFreq,Op0X,Op0Y,Op1X,Op1Y):

    Opflow0X = []
    Opflow0Y = []
    Opflow1X = []
    Opflow1Y = []
    ini = 0
    for i in range(len(diffList)):
        lenListToInterp = diffList[i]
        LinspaceTemp=np.linspace(0,lenListToInterp-1,lenListToInterp)
        LinspaceGoodLength=np.linspace(0,lenListToInterp-1,maxFreq)
        Opflow0X.append(np.interp(LinspaceGoodLength, LinspaceTemp, Op0X[ini:ini+diffList[i]]))
        Opflow0Y.append(np.interp(LinspaceGoodLength, LinspaceTemp, Op0Y[ini:ini+diffList[i]]))
        Opflow1X.append(np.interp(LinspaceGoodLength, LinspaceTemp, Op1X[ini:ini+diffList[i]]))
        Opflow1Y.append(np.interp(LinspaceGoodLength, LinspaceTemp, Op1Y[ini:ini+diffList[i]]))
        ini+=diffList[i]

    OpInter0X = np.concatenate(Opflow0X[:])
    OpInter0Y = np.concatenate(Opflow0Y[:])
    OpInter1X = np.concatenate(Opflow1X[:])
    OpInter1Y = np.concatenate(Opflow1Y[:])

    return OpInter0X, OpInter0Y, OpInter1X, OpInter1Y

def upsampleOpflow(diffList,Op0X,Op0Y,Op1X,Op1Y):
    
    Opflow0X = []
    Opflow0Y = []
    Opflow1X = []
    Opflow1Y = []

    for i in range(len(diffList)):
        Op0XTemp = [Op0X[i]]*diffList[i]
        Op0YTemp = [Op0Y[i]]*diffList[i]
        Op1XTemp = [Op1X[i]]*diffList[i]
        Op1YTemp = [Op1Y[i]]*diffList[i]
        #print ("Op0XTemp",Op0XTemp)
        Opflow0X.extend(Op0XTemp)
        #print ("Opflow0X",Opflow0X)
        Opflow0Y.extend(Op0YTemp)
        Opflow1X.extend(Op1XTemp)
        Opflow1Y.extend(Op1YTemp)


    return Opflow0X, Opflow0Y, Opflow1X, Opflow1Y


def upsampleGC(GCset, frameCntr, numFrame):
    GCset_US=[]
    GCset_UStemp1=[]
    GCset_UStemp2=[]



    print('upsample len frameCntr', len(frameCntr))
    print('upsample len GCset[0]',len(GCset[0]))
    print('upsample len(GCset[0]*numFrame)', len(GCset[0]*numFrame))

    for i in range(0, len(frameCntr)):
        if frameCntr[i]==int(numFrame/3):
            GCStartIdx=i
            break
    for j in range(0, len(frameCntr)): # bottle neck
        if frameCntr[j]==len(GCset[0])*numFrame:
            GCStopIdx=j
            break

    print("GCset[0][219]", GCset[0][219])
    print("GCStartIdx",GCStartIdx)
    print("GCStopIdx",GCStopIdx)
    print("len GC[0]", len(GCset[0]))



    LinspaceTemp=np.linspace(0,len(GCset[0])-1,len(GCset[0]))
    LinspaceGoodLength=np.linspace(0,len(GCset[0])-1,len(frameCntr[GCStartIdx:GCStopIdx]))
    for i in range(0, len(GCset)):
        GCset_UStemp1.append([])
        GCset_UStemp1[i]=np.interp(LinspaceGoodLength, LinspaceTemp, GCset[i])

    lengthInitGC=GCStartIdx
    lengthTailGC=len(frameCntr)-GCStopIdx

    nanInitGC=[]
    for a in range(0,lengthInitGC):
        nanInitGC.append(np.nan)

    nanTailGC=[]
    for b in range(0,lengthTailGC):
        nanTailGC.append(np.nan)    

    print("len nanInitGC",len(nanInitGC))
    print("len nanTailGC",len(nanTailGC))

    for i in range(0, len(GCset)):
        GCset_UStemp2.append([])
        GCset_UStemp2[i]=np.append(nanInitGC,GCset_UStemp1[i])
        GCset_US.append([])
        GCset_US[i]=np.append(GCset_UStemp2[i],nanTailGC)

    print("len frameCntr", len(frameCntr))
    print("len GCset_US", len(GCset_US))
    print("len GCset_US[0]", len(GCset_US[0]))
    print("GCset_US[0][15000]",GCset_US[0][15000])

    return GCset_US


def upsample_risingEdgeCam(risingEdgeCam, Cam):

    print('Upsample camera rising edge data ...')

    print('len Cam', len(Cam))

    zeroDet = np.where(Cam==0)
    arrayZD_risingEdgeCam = zeroDet[0]

    print('arrayZD_risingEdgeCam', arrayZD_risingEdgeCam)
    print('len arrayZD_risingEdgeCam', len(arrayZD_risingEdgeCam))

    risingEdgeCam_Stair=[0]*len(Cam)

    print('len risingEdgeCam_Stair', len(risingEdgeCam_Stair))

    stair_idx_start=0
    for idx, val in enumerate(risingEdgeCam):
        #print('val', val)
        stair_idx_end=val
        #print('stair_idx_start', stair_idx_start, 'stair_idx_end', stair_idx_end)
        #print('risingEdgeCam_Stair[stair_idx_start:stair_idx_end]', risingEdgeCam_Stair[stair_idx_start:stair_idx_end])
        #print('len risingEdgeCam_Stair[stair_idx_start:stair_idx_end]', len(risingEdgeCam_Stair[stair_idx_start:stair_idx_end]))
        risingEdgeCam_Stair[stair_idx_start:stair_idx_end]=[val]*(stair_idx_end-stair_idx_start)
        stair_idx_start=val

    print('stair_idx_end', stair_idx_end)
    print('risingEdgeCam[-1]', risingEdgeCam[-1])
    print('risingEdgeCam_Stair[stair_idx_end-1]', risingEdgeCam_Stair[stair_idx_end-1])
    print('risingEdgeCam_Stair[stair_idx_end]', risingEdgeCam_Stair[stair_idx_end])
    print('risingEdgeCam_Stair[-1]', risingEdgeCam_Stair[-1])
    print('len(risingEdgeCam_Stair)-stair_idx_end', len(risingEdgeCam_Stair)-stair_idx_end)
    print('len risingEdgeCam_Stair[stair_idx_end:]', len(risingEdgeCam_Stair[stair_idx_end:]))
    #print('risingEdgeCam_Stair[stair_idx_end:-1]', risingEdgeCam_Stair[stair_idx_end:-1])

    risingEdgeCam_Stair[stair_idx_end:] = [risingEdgeCam[-1]]*(len(risingEdgeCam_Stair)-stair_idx_end)

    print('risingEdgeCam_Stair[-1]', risingEdgeCam_Stair[-1])
    print('len risingEdgeCam_Stair', len(risingEdgeCam_Stair))

    print('Upsample camera rising edge data is done!')



    return risingEdgeCam_Stair


def upsampleBeh_old(trace_from_7cam, CamTraceThorsync, CamPulsePoints):
    
    print('upsampling behavior data...')
    trace_from_7cam_US=len(CamTraceThorsync)*[np.nan]

    print('CamPulsePoints', CamPulsePoints)
    print('CamTraceThorsync', CamTraceThorsync)
    print('len CamPulsePoints', len(CamPulsePoints))
    print('len CamTraceThorsync', len(CamTraceThorsync))
    print('len trace_from_7cam', len(trace_from_7cam))


    CamPulsePoints=CamPulsePoints[0:len(trace_from_7cam_US)]

    print('CamPulsePoints', CamPulsePoints)
    print('CamPulsePoints[0]', CamPulsePoints[0])
    print('CamPulsePoints[-1]', CamPulsePoints[-1])
    print('len CamPulsePoints', len(CamPulsePoints))




    for i in range(1, len(CamTraceThorsync)):
        if i <= CamPulsePoints[-1]:
            if i in CamPulsePoints:
                # print(i)
                # print(list(CamPulsePoints).index(i))

                ## make sure camera pulses within the amounts of actual photos
                if list(CamPulsePoints).index(i) < len(trace_from_7cam):
                    trace_from_7cam_US[i]=trace_from_7cam[list(CamPulsePoints).index(i)]
                else:
                    #the ending
                    trace_from_7cam_US[i]=trace_from_7cam[-1]

            else:

                trace_from_7cam_US[i]=trace_from_7cam_US[i-1]

        else:
            break

    # print('len trace_from_7cam_US', len(trace_from_7cam_US))
    # print('trace_from_7cam_US[CamPulsePoints[0]-1]', trace_from_7cam_US[CamPulsePoints[0]-1])
    # print('trace_from_7cam_US[CamPulsePoints[0]]', trace_from_7cam_US[CamPulsePoints[0]])
    # print('trace_from_7cam_US[CamPulsePoints[0]]', trace_from_7cam_US[CamPulsePoints[1]])
    # print('trace_from_7cam_US[CamPulsePoints[-1]+1]', trace_from_7cam_US[CamPulsePoints[-1]+1])
    # print('trace_from_7cam_US[CamPulsePoints[-1]]', trace_from_7cam_US[CamPulsePoints[-1]])
    # print('trace_from_7cam_US[CamPulsePoints[-2]]', trace_from_7cam_US[CamPulsePoints[-2]])
    print('upsampling behavior data done!')




    return trace_from_7cam_US



def upsampleBeh(trace_from_7cam, CamTraceThorsync, CamPulsePoints):


    #print('upsampling behavior data...')
    trace_from_7cam_US=np.asarray(len(CamTraceThorsync)*[np.nan])

    #print('len trace_from_7cam', len(trace_from_7cam))
    # print('CamPulsePoints', CamPulsePoints)
    # print('CamTraceThorsync', CamTraceThorsync)
    #print('len CamPulsePoints', len(CamPulsePoints))
    # print('CamPulsePoints[0]', CamPulsePoints[0])
    # print('CamPulsePoints[-1]', CamPulsePoints[-1])
    # print('len CamTraceThorsync', len(CamTraceThorsync))
    # print('len trace_from_7cam', len(trace_from_7cam))


    CamPulsePoints=CamPulsePoints[0:len(trace_from_7cam_US)]
    # print('CamPulsePoints[0]', CamPulsePoints[0])
    # print('CamPulsePoints[-1]', CamPulsePoints[-1])
    # print('CamPulsePoints[-1]-CamPulsePoints[0]', CamPulsePoints[-1]-CamPulsePoints[0])

    # for pulseIdx in CamPulsePoints:
    #     trace_from_7cam_US[pulseIdx]=trace_from_7cam[list(CamPulsePoints).index(pulseIdx)]



    for idx, camPulse_idx in enumerate(CamPulsePoints):
        if idx<len(trace_from_7cam):

            # print(idx)
            # print(trace_from_7cam[idx])


            ## For those have no additional tail in campulse points
            #trace_from_7cam_US[camPulse_idx:CamPulsePoints[idx]]=trace_from_7cam[idx]

            ## For those really have additional tail in campulse points
            trace_from_7cam_US[camPulse_idx:CamPulsePoints[idx+1]]=trace_from_7cam[idx]

        elif idx==len(trace_from_7cam):
            # print(idx)
            # print(trace_from_7cam[idx-1])
            # trace_from_7cam_US[camPulse_idx:CamPulsePoints[idx+1]]=trace_from_7cam[idx-1]
            # #trace_from_7cam_US[camPulse_idx]=trace_from_7cam[idx]

            break

    # print('len trace_from_7cam_US', len(trace_from_7cam_US))
    # print('trace_from_7cam_US[CamPulsePoints[0]-1]', trace_from_7cam_US[CamPulsePoints[0]-1])
    # print('trace_from_7cam_US[CamPulsePoints[0]]', trace_from_7cam_US[CamPulsePoints[0]])
    # print('trace_from_7cam_US[CamPulsePoints[0]]', trace_from_7cam_US[CamPulsePoints[1]])
    # print('trace_from_7cam_US[CamPulsePoints[-1]+1]', trace_from_7cam_US[CamPulsePoints[-1]+1])
    # print('trace_from_7cam_US[CamPulsePoints[-1]]', trace_from_7cam_US[CamPulsePoints[-1]])
    # print('trace_from_7cam_US[CamPulsePoints[-2]]', trace_from_7cam_US[CamPulsePoints[-2]])


    return trace_from_7cam_US





def check_and_determine_numframe(Framecounter, RGBstack):

    print('Framecounter[-1]', Framecounter[-1])
    print('len(RGBstack)', len(RGBstack))

    numframe_temp=(Framecounter[-1])/(len(RGBstack))
    print('numframe_temp', numframe_temp)

    if numframe_temp.is_integer():
        numframe=numframe_temp
        print('numFrame=', numframe)
    else:

        print('framecounters does not match the 2Pimg numbers, it is not dividible')
        sys.exit(0)

    return numframe


def getXAxisAndRatio(timeSec,startVidIdx,stopVidIdx):
    #This function returns the x axis min and max values and the ratio of the number of points in the global time list to the number of points in the behavior time list.

    Ratio_HDseriestoCam=float(len(timeSec))/float(stopVidIdx-startVidIdx)

    return Ratio_HDseriestoCam



def OpenOpflowVector(opflowPath):
    #This function opens the optic flow measurements values initially stored in a txt file.
    opFlowCols = ['sens0X','sens0Y','sens1X','sens1Y','date','time']
    #opFlowAll = pd.read_table(opflowPath, sep='\s+', header=None, names=opFlowCols)
    opFlowAll = pd.read_table(opflowPath, sep=',', header=None, names=opFlowCols)

    op0XSensor = pd.Series(opFlowAll['sens0X'][:]).values
    op0YSensor = pd.Series(opFlowAll['sens0Y'][:]).values
    op1XSensor = pd.Series(opFlowAll['sens1X'][:]).values
    op1YSensor = pd.Series(opFlowAll['sens1Y'][:]).values
    TimeO = pd.Series(opFlowAll['time'][:]).values
    print ("TimeO",TimeO)
    print ("op0XSensor",op0XSensor)
    print ("len op0XSensor", len(op0XSensor))



    return op0XSensor, op0YSensor, op1XSensor, op1YSensor, TimeO

def manageCompleteSensors(opFlowSensorVal0X, opFlowSensorVal0Y, opFlowSensorVal1X, opFlowSensorVal1Y):

    if math.isnan(opFlowSensorVal1Y[-1])==True:
        opFlowSensorVal0X=opFlowSensorVal0X[:-1]
        opFlowSensorVal0Y=opFlowSensorVal0Y[:-1]
        opFlowSensorVal1X=opFlowSensorVal1X[:-1]
        opFlowSensorVal1Y=opFlowSensorVal1Y[:-1]

    return opFlowSensorVal0X, opFlowSensorVal0Y, opFlowSensorVal1X, opFlowSensorVal1Y

def manageTypeOpflow(opFlowSensorVal0X, opFlowSensorVal0Y, opFlowSensorVal1X, opFlowSensorVal1Y):

    if type(opFlowSensorVal0X[-1])!="float64":
        opFlowSensorVal0X =  opFlowSensorVal0X.astype(np.float)
    if type(opFlowSensorVal0Y[-1])!="float64":
        opFlowSensorVal0Y =  opFlowSensorVal0Y.astype(np.float)
    if type(opFlowSensorVal1X[-1])!="float64":
        opFlowSensorVal1X =  opFlowSensorVal1X.astype(np.float)
    if type(opFlowSensorVal1Y[-1])!="float64":
        opFlowSensorVal1Y =  opFlowSensorVal1Y.astype(np.float)

    return opFlowSensorVal0X, opFlowSensorVal0Y, opFlowSensorVal1X, opFlowSensorVal1Y

def findStartAndStopIdx(zero_crossings_Pos_Cam,zero_crossings_Pos_Opflow,startFrameCounterIdx,stopFCTRIdx):

    startIdxTSLTemp = max(zero_crossings_Pos_Cam[0],startFrameCounterIdx,zero_crossings_Pos_Opflow[0])
    if startIdxTSLTemp == zero_crossings_Pos_Cam[0]:
        startIdxTSL = startIdxTSLTemp
    else : 
        val, idx = find_nearest(zero_crossings_Pos_Cam,startIdxTSLTemp)
        if val > startIdxTSLTemp:
            startIdxTSL = val
        else : 
            startIdxTSL = zero_crossings_Pos_Cam[idx+1]

    stopIdxTSLTemp = min(zero_crossings_Pos_Cam[-1],zero_crossings_Pos_Opflow[-1],stopFCTRIdx)
    if stopIdxTSLTemp == zero_crossings_Pos_Cam[-1]:
        stopIdxTSL == stopIdxTSLTemps
    else :
        valStop, idxStop = find_nearest(zero_crossings_Pos_Cam,stopIdxTSLTemp)
        if valStop < stopIdxTSLTemp:
            stopIdxTSL = valStop
        else : 
            stopIdxTSL = zero_crossings_Pos_Cam[idxStop-1]
    print ("zero_crossings_Pos_Cam[-1]",zero_crossings_Pos_Cam[-1])
    print ("zero_crossings_Pos_Opflow[-1]",zero_crossings_Pos_Opflow[-10:-1])
    print ("stopFCTRIdx",stopFCTRIdx)
    print ("stopIdxTSL",stopIdxTSL)

    vidStopTime = (stopIdxTSL-startIdxTSL)/freqTS

    return startIdxTSL, stopIdxTSL, vidStopTime

def risingEdge(value,ArrayToDetect):

    TempArray = ArrayToDetect.astype(int)
    DiffArray = TempArray-value
    PosZeroCrossing = np.where(np.diff(np.sign(DiffArray))>0)[0] 

    return PosZeroCrossing

def startStopIdxFCtr(FCTR):

    FCtrTemp = FCTR.astype(int)
    startFCIdx = np.nonzero(FCtrTemp)[0][0]
    stopFCIdx = np.where(FCtrTemp==FCtrTemp[-1])[0][0]

    return startFCIdx, stopFCIdx

def compareOpflowTSReal(risingEdge, Sensor): 

    if len(Sensor)<len(risingEdge):
        print ("OPFLOW CROPPED")
        print ("len(Sensor)",len(Sensor))
        print ("len(risingEdge)",len(risingEdge))
        risingEdgeCropped = risingEdge[:len(Sensor)]
    else:
        risingEdgeCropped = risingEdge[:]
    
    storeDiff = []
    for i in range(1,len(risingEdgeCropped)):
        storeDiff.append(risingEdgeCropped[i]-risingEdgeCropped[i-1])

    return storeDiff, risingEdgeCropped

def cropOpFlow(start,stop,O0XUS, O0YUS, O1XUS, O1YUS,RSOC):




    if start>RSOC[0]:
        newStartOIdx = start-RSOC[0]
    elif start==RSOC[0]:
        newStartOIdx = 0
    if stop < RSOC[-1]:
        newStopOIdx = RSOC[-1] - stop
    elif stop == RSOC[-1]:
        newStopOIdx = 0

    O0XCUS = O0XUS[newStartOIdx:len(O0XUS)-newStopOIdx]
    O0YCUS = O0YUS[newStartOIdx:len(O0YUS)-newStopOIdx]
    O1XCUS = O1XUS[newStartOIdx:len(O1XUS)-newStopOIdx]
    O1YCUS = O1YUS[newStartOIdx:len(O1YUS)-newStopOIdx]

    return O0XCUS, O0YCUS, O1XCUS, O1YCUS

def findStartStopCamera(startIdx,stopIdx,risingEC):
    # the start and stop always correspond to a start and stop of Cam idx so no need to check if nearest is above or below the value in this case

    if startIdx==risingEC[0]:
        startCam = 0
    else : 
        valStart, idxStart = find_nearest(risingEC,startIdx)
        startCam = idxStart

    if stopIdx == risingEC[-1]:
        stopCam = len(risingEC)
    else : 
        valStop, idxStop = find_nearest(risingEC,stopIdx)
        stopCam = idxStop
    return startCam, stopCam




def CropFrameCtr(trace, startIdxFC, stopIdxFC):

    cropped_trace = trace[startIdxFC:stopIdxFC]

    return cropped_trace


def downSampling(cropTrace,stopTime):


    currentLen = len(cropTrace)
    # print("len O0XCUSTD",len(O0XCUSTD))
    # print("currentLen",currentLen)
    # print("currentLen/(stopTime*samplingFact)",currentLen/(stopTime*samplingFact))
    #lengthAdaptor = int(round(currentLen/(stopTime*samplingFact)))
    # print ("stopTime",stopTime)
    # print ("samplingFact",samplingFact)
    
    lengthAdaptor = currentLen/(stopTime*samplingFact)
    # print('lengthAdaptor', lengthAdaptor)

    if lengthAdaptor-int(lengthAdaptor)==0.99:

        lengthAdaptor = int(currentLen/(stopTime*samplingFact))
    else:
        lengthAdaptor = int(round(currentLen/(stopTime*samplingFact)))


    traceCUSDS=[]
    
    for j in range(currentLen):
        if j%lengthAdaptor==0:
            traceCUSDS.append(cropTrace[j])


    return traceCUSDS




def getSensorChannel(windowLength,sensor,gain):
    # This function performs the smoothing window on opflow and returns each signal times its gain

    window = np.ones(int(windowLength))/float(windowLength)
    smooth_series=np.convolve(sensor,window,'same')
    SensorTG = smooth_series*gain
    if len(SensorTG)<len(sensor):
        nanTail = np.empty(len(sensor)-len(SensorTG))
        nanTail.fill(np.nan)
        sensorComplete = np.append(SensorTG,nanTail)
    else:
        sensorComplete=SensorTG
    return sensorComplete

def GetVelocities(sensor0X,sensor0Y,sensor1X,sensor1Y):
    #This function computes the AP, ML and Yaw rot/s from the sensor measurements.

    velForwF = -((sensor0Y + sensor1Y) * np.cos(np.deg2rad(45)));     #(Y1 + Y2) * cos(45) [invert for clarity]
    velSideF = (sensor0Y - sensor1Y) * np.sin(np.deg2rad(45));        #(Y1 - Y2) * sin(45)
    velTurnF = (sensor0X + sensor1X) / float(2);    #(X1 + X2)/2 

    return velForwF, velSideF, velTurnF  





def normalize_Data(data):
    #print("CO2puff", CO2puff)
    maxdata=max(data)
    mindata=min(data)
    norm_data_list=[]

    for i in range(0, len(data)):
        if maxdata-mindata!=0:
            norm_data=(data[i]-mindata)/(maxdata-mindata)
            norm_data_list.append(norm_data)
        else:
            norm_data_list.append(0)


    return np.asarray(norm_data_list)

def convert_opflow(velForwRot, velSideRot, velTurnRot):

    velForw_mm = velForwRot*rot_to_mm_factor 
    velSide_mm = velSideRot*rot_to_mm_factor
    velForw_deg = velTurnRot*rot_to_deg_factor

    return velForw_mm, velSide_mm, velForw_deg




def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx


def trim_period(trace, timeSec, startTime=0, endTime=150):

    # print('len trace', len(trace))

    crop_startTime, crop_startIdx = find_nearest(timeSec, startTime)


    if startTime == 0:
        crop_startTime, crop_startIdx = find_nearest(timeSec, timeSec[0])
        crop_startIdx=0
    else:
        crop_startTime, crop_startIdx = find_nearest(timeSec, startTime)   

    if endTime == -1:
        crop_stopTime, crop_stopIdx = find_nearest(timeSec, timeSec[-1])
        crop_stopIdx=-1
    else:
        crop_stopTime, crop_stopIdx = find_nearest(timeSec, endTime)

    print('crop_startTime', crop_startTime, 'crop_startIdx', crop_startIdx)
    print('crop_stopTime', crop_stopTime, 'crop_stopIdx', crop_stopIdx)


    trace_crop = trace[crop_startIdx:crop_stopIdx]



    return trace_crop




def Plot_smooth_op_with_class(velForw_mm, velSide_mm, velTurn_deg, timesec, ethoTimeDic, ethoColorDic):

    datapoints_to_time_factor=len(velForw_mm)/timesec[-1]
    print('datapoints_to_time_factor',datapoints_to_time_factor)
    trace_dur = np.linspace(0,timeSec[-1]-timeSec[0],len(velForw_mm[:]))

    # plt.plot(velTurn)
    # plt.pause(2)
    Totalnum_row=3
    rowspan = int(Totalnum_row/Totalnum_row)
    row_count=rowspan

    #smth optic flow
    # velForw_mm=smooth_data(velForw_mm, datapoints_to_time_factor*window_for_smth)
    # velSide_mm=smooth_data(velSide_mm, datapoints_to_time_factor*window_for_smth)
    # velTurn_deg=smooth_data(velTurn_deg, datapoints_to_time_factor*window_for_smth)

    
    fig = plt.figure(facecolor='white', figsize=(Totalnum_row,Totalnum_row), dpi=300)
    fig.subplots_adjust(left=0.2, right = 0.9, wspace = 0.3, hspace = 0.3)

    axAP = plt.subplot(Totalnum_row,1,row_count)
    axAP.plot(trace_dur[1000:], velForw_mm[1000:], color=AP_color,linewidth=1)
    #axvelForw.plot(trace_dur, velForw_smth, color=AP_color,linewidth=1)
    row_count+=rowspan

    axML = plt.subplot(Totalnum_row,1,row_count)
    axML.plot(trace_dur[1000:], velSide_mm[1000:], color=ML_color,linewidth=1)
    #axvelSide.plot(trace_dur, velSide_smth, color=ML_color,linewidth=1)
    row_count+=rowspan

    axYaw = plt.subplot(Totalnum_row,1,row_count)
    axYaw.plot(trace_dur[1000:], velTurn_deg[1000:], color=Yaw_color,linewidth=1)
    #axvelForw.plot(trace_dur, velturn_smth, color=Yaw_color,linewidth=1)
    row_count+=rowspan

    i=0
    for key in ethoTimeDic: 
        if key =='walk_evt' or key=='rest_evt':
            for a in range(0,len(ethoTimeDic[key])):

                axAP.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=ethoColorDic[key], linewidth=0)
                axML.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=ethoColorDic[key], linewidth=0)
                axYaw.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=ethoColorDic[key], linewidth=0)  

            i+=1
            if i >3:
                break

    plt.savefig(outDir + 'beh_clsfy_on_opflow_by_replace.svg', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    plt.clf()
    plt.close()


    return








def reclassify_walk_rest_by_ball_rotation(velForw, velSide, velTurn, OtherBeh_merged_trace, timesec):

    print('reclassify_walk_rest_by_ball_rotation....')
    datapoints_to_time_factor=len(velForw)/timesec[-1]
    print('datapoints_to_time_factor',datapoints_to_time_factor)  

    coarse_beh_series=len(velForw)*[0]




    for i in range(0,len(velForw)):
        # classify RESTING
        if abs(OtherBeh_merged_trace[i])<0.5:
            if abs(velForw[i])<AP_thrsld and abs(velSide[i])<ML_thrsld and abs(velTurn[i])<Yaw_thrsld:
                coarse_beh_series[i]='r'

        # classify WALK
            elif abs(velForw[i])>=AP_thrsld or abs(velSide[i])>=ML_thrsld or abs(velTurn[i])>=Yaw_thrsld:
                coarse_beh_series[i]='w'

        # classify Other
        elif abs(OtherBeh_merged_trace[i])>0.5:
            coarse_beh_series[i]='x'


    norm_rest=len(coarse_beh_series)*[0]
    norm_walk=len(coarse_beh_series)*[0]

    for i, beh in enumerate(coarse_beh_series):
        if beh == 'r':  
            norm_rest[i]=1
        
        elif beh == 'w':            
            norm_walk[i]=1


    #print('norm_rest[:500]', norm_rest[:500])



    return norm_rest, norm_walk





def make_bin_trace_from_behLabels(behLabels_series, target_beh):

    behFullName_list=    ['rest', 'walk', 'eye_groom','antennal_groom','foreleg_groom','hindleg_groom','Abd_groom','Push']
    behLabels_list = ['r','w',  'eg','antg','flg','hlg','abdg','p']

    beh_idx = behFullName_list.index(target_beh)
    target_beh_label=behLabels_list[beh_idx]


    beh_bin_trace=len(behLabels_series)*[0]

    for i in range(0, len(beh_bin_trace)):

        if behLabels_series[i]==target_beh_label:
            beh_bin_trace[i]=1

    return beh_bin_trace




def replacing_short_evt_with_nearbyDominantClass(beh_bin_trace_dic, beh_symbol_series, timeSec, dur_thrsld = 0.5, dur_scoring=0.5, dataFreq=1500):


    beh_symbol_series_rplcd = beh_symbol_series.copy()


    for beh_class, beh_bin_trace in beh_bin_trace_dic.items():

    

        beh_bin_trace_cp=beh_bin_trace.copy()

        idx_beh_evt, _ = Calculate_idx_time_for_bin_beh_trace(beh_bin_trace_cp, timeSec)

        for i, evt_idx in enumerate(idx_beh_evt):

            if len(evt_idx)<dur_thrsld*dataFreq:
                init_search_idx=int(evt_idx[0]-dur_scoring*dataFreq)
                end_search_idx=int(evt_idx[-1]+dur_scoring*dataFreq)
                if init_search_idx<0:
                    init_search_idx=0
                if end_search_idx>len(beh_symbol_series):
                    end_search_idx=-1

                r_score = beh_symbol_series[init_search_idx:end_search_idx].count('r')
                w_score = beh_symbol_series[init_search_idx:end_search_idx].count('w')
                eg_score = beh_symbol_series[init_search_idx:end_search_idx].count('eg')
                antg_score = beh_symbol_series[init_search_idx:end_search_idx].count('antg')
                flg_score = beh_symbol_series[init_search_idx:end_search_idx].count('flg')
                hlg_score = beh_symbol_series[init_search_idx:end_search_idx].count('hlg')
                abdg_score = beh_symbol_series[init_search_idx:end_search_idx].count('abdg')
                p_score = beh_symbol_series[init_search_idx:end_search_idx].count('p')


                behFullName_list=    ['rest','walk',  'eye_groom','antennal_groom','foreleg_groom','hindleg_groom','Abd_groom','Push']
                behLabels_list = ['r','w',  'eg','antg','flg','hlg','abdg','p']
                score_list = [r_score, w_score, eg_score,  antg_score,      flg_score,      hlg_score,      abdg_score, p_score]
                max_score = max(score_list)
                max_index = score_list.index(max_score)

                new_behLabel = behLabels_list[max_index]
                #print('new_behLabels', new_behLabels)
                if new_behLabel == beh_class:
                    continue
                else:
                    beh_symbol_series_rplcd[evt_idx[0]:evt_idx[-1]+1]=len(evt_idx)*[new_behLabel]
                    # beh_symbol_series[evt_idx[0]:evt_idx[-1]+1]=len(evt_idx)*[new_behLabels]


    
    return beh_symbol_series_rplcd



def beh_clsfy_by_opflow(velForw, velSide, velTurn, Groomtrace, timesec):

    print('beh_clsfy_by_opflow....')
    
    datapoints_to_time_factor=len(velForw)/timesec[-1]
    print('datapoints_to_time_factor',datapoints_to_time_factor)

    coarse_beh_series=len(velForw)*[0]

    idx_rest=[]
    rest_series=len(velForw)*[0]
    idx_walk=[]
    idx_groom=[]

    

    for i in range(0,len(velForw)):
        # classify RESTING
        if abs(Groomtrace[i])<0.5:
            if abs(velForw[i])<AP_thrsld and abs(velSide[i])<ML_thrsld and abs(velTurn[i])<Yaw_thrsld and abs(Groomtrace[i])<0.5:
                idx_rest.append(i)
                coarse_beh_series[i]='r'

        # classify WALK
            elif abs(velForw[i])>=AP_thrsld or abs(velSide[i])>=ML_thrsld or abs(velTurn[i])>=Yaw_thrsld:
                idx_walk.append(i)
                coarse_beh_series[i]='w'

        # classify GROOM
        elif abs(Groomtrace[i])>0.5:
            idx_groom.append(i)
            coarse_beh_series[i]='g'

    idx_rest_evt=grouping_consecutivePoints_into_evt(idx_rest)
    idx_walk_evt=grouping_consecutivePoints_into_evt(idx_walk)
    idx_groom_evt=grouping_consecutivePoints_into_evt(idx_groom)

    ##filtering out the short event by replacing the nearby major beh class based on scoring the length of nearby beh
    for i, evt_idx in enumerate(idx_rest_evt):
        if len(evt_idx)<duration_for_replace*datapoints_to_time_factor:
            RestScore = coarse_beh_series[int(evt_idx[0]-duration_for_scoring*datapoints_to_time_factor):int(evt_idx[-1]+duration_for_scoring*datapoints_to_time_factor)].count('r')
            WalkScore = coarse_beh_series[int(evt_idx[0]-duration_for_scoring*datapoints_to_time_factor):int(evt_idx[-1]+duration_for_scoring*datapoints_to_time_factor)].count('w')
            GroomScore = coarse_beh_series[int(evt_idx[0]-duration_for_scoring*datapoints_to_time_factor):int(evt_idx[-1]+duration_for_scoring*datapoints_to_time_factor)].count('g')

            if RestScore==max(RestScore,WalkScore, GroomScore):
                continue
            elif WalkScore==max(RestScore,WalkScore, GroomScore):
                coarse_beh_series[evt_idx[0]:evt_idx[-1]+1]=len(evt_idx)*['w']
            elif GroomScore==max(RestScore,WalkScore, GroomScore):
                coarse_beh_series[evt_idx[0]:evt_idx[-1]+1]=len(evt_idx)*['g']

    for i, evt_idx in enumerate(idx_walk_evt):
        if len(evt_idx)<duration_for_replace*datapoints_to_time_factor:
            RestScore = coarse_beh_series[int(evt_idx[0]-duration_for_scoring*datapoints_to_time_factor):int(evt_idx[-1]+duration_for_scoring*datapoints_to_time_factor)].count('r')
            WalkScore = coarse_beh_series[int(evt_idx[0]-duration_for_scoring*datapoints_to_time_factor):int(evt_idx[-1]+duration_for_scoring*datapoints_to_time_factor)].count('w')
            GroomScore = coarse_beh_series[int(evt_idx[0]-duration_for_scoring*datapoints_to_time_factor):int(evt_idx[-1]+duration_for_scoring*datapoints_to_time_factor)].count('g')

            if RestScore==max(RestScore,WalkScore, GroomScore):
                coarse_beh_series[evt_idx[0]:evt_idx[-1]+1]=len(evt_idx)*['r']              
            elif WalkScore==max(RestScore,WalkScore, GroomScore):
                continue
            elif GroomScore==max(RestScore,WalkScore, GroomScore):
                coarse_beh_series[evt_idx[0]:evt_idx[-1]+1]=len(evt_idx)*['g']

    for i, evt_idx in enumerate(idx_groom_evt):
        if len(evt_idx)<duration_for_replace*datapoints_to_time_factor:
            RestScore = coarse_beh_series[int(evt_idx[0]-duration_for_scoring*datapoints_to_time_factor):int(evt_idx[-1]+duration_for_scoring*datapoints_to_time_factor)].count('r')
            WalkScore = coarse_beh_series[int(evt_idx[0]-duration_for_scoring*datapoints_to_time_factor):int(evt_idx[-1]+duration_for_scoring*datapoints_to_time_factor)].count('w')
            GroomScore = coarse_beh_series[int(evt_idx[0]-duration_for_scoring*datapoints_to_time_factor):int(evt_idx[-1]+duration_for_scoring*datapoints_to_time_factor)].count('g')

            if RestScore==max(RestScore,WalkScore, GroomScore):
                coarse_beh_series[evt_idx[0]:evt_idx[-1]+1]=len(evt_idx)*['r']
            elif WalkScore==max(RestScore,WalkScore, GroomScore):
                coarse_beh_series[evt_idx[0]:evt_idx[-1]+1]=len(evt_idx)*['w']
            elif GroomScore==max(RestScore,WalkScore, GroomScore):
                continue
    
    # Reassign coarse beahvior trace into binary trace of each behavior class
    
    norm_rest=len(coarse_beh_series)*[0]
    norm_walk=len(coarse_beh_series)*[0]
    norm_groom=len(coarse_beh_series)*[0]

    for i, beh in enumerate(coarse_beh_series):
        if beh == 'r':  
            norm_rest[i]=1
        
        elif beh == 'w':            
            norm_walk[i]=1
        
        elif beh == 'g':            
            norm_groom[i]=1


    return coarse_beh_series, norm_rest, norm_walk, norm_groom





def grouping_consecutivePoints_into_evt(idx_beh):

    idx_beh_evt=[]
    for evt in mit.consecutive_groups(idx_beh):
        #type(evt) is map, once it is converted into list, it become an empty array
        #print(list(evt))
        idx_beh_evt.append(list(evt))    

    return idx_beh_evt



def convert_idx_to_timepoint(idx_beh_evt, timesec):

    timesec_beh_evt=[]
    for a in range(0,len(idx_beh_evt)):
        timesec_beh_evt.append([])
        for b in range(0, len(idx_beh_evt[a])):
            timesec_beh_evt[a].append(timesec[idx_beh_evt[a][b]])
    

    return timesec_beh_evt


def Calculate_idx_time_for_bin_beh_trace(beh_bin_trace, timesec):


    idx_beh=[]
    for i, val in enumerate(beh_bin_trace):
        if val>0.5:
            idx_beh.append(i)


    idx_beh_evt=grouping_consecutivePoints_into_evt(idx_beh)
    timesec_beh_evt=convert_idx_to_timepoint(idx_beh_evt, timesec)



    return idx_beh_evt, timesec_beh_evt



def make_new_bin_trace_excluding_otherBeh(main_bin_trace, excld_beh_bin_trace=[], timeSec=[]):

    if len(timeSec)==0:
        timeSec=np.arange(len(main_bin_trace))

    if len(main_bin_trace)!=len(excld_beh_bin_trace):
        exclded_main_idx_evts=main_idx_evts

    else:
        
        idx_main_evt, timesec_main_evt=Calculate_idx_time_for_bin_beh_trace(main_bin_trace, timeSec)
        idx_exludBeh_evt, timesec_exludBeh_evt=Calculate_idx_time_for_bin_beh_trace(excld_beh_bin_trace, timeSec)

        flatten_excldBeh_idx_evts = [item for subl in idx_exludBeh_evt for item in subl]

        exclded_main_idx_evts=[]
        for i, evt_idxs in enumerate(idx_main_evt):
            S1=set(evt_idxs)
            S2=set(flatten_excldBeh_idx_evts)
            S_intsct=S1.intersection(S2)
            if len(S_intsct)==0:
                exclded_main_idx_evts.append(evt_idxs)


        new_bin_trace=np.asarray(list(itertools.repeat(0,len(main_bin_trace))))
        print('new_bin_trace',new_bin_trace)
        print('type new_bin_trace', type(new_bin_trace))

        for i, evt_idxs in enumerate(exclded_main_idx_evts):
            new_bin_trace[np.asarray(evt_idxs)]=1



    return new_bin_trace






def Package_beh_bin_trace_into_idx_time(rest_bin_trace, walk_bin_trace, groom_bin_trace, timesec):

    print('len rest_bin_trace', len(rest_bin_trace))

    idx_rest=[]
    for i, val in enumerate(rest_bin_trace):
        if val>0.5:
            idx_rest.append(i)

    idx_walk=[]
    for i, val in enumerate(walk_bin_trace):
        if val>0.5:
            idx_walk.append(i)

    idx_groom=[]
    for i, val in enumerate(groom_bin_trace):
        if val>0.5:
            idx_groom.append(i)



    ##Package behavior binary trace into event index and time
    idx_rest_evt=grouping_consecutivePoints_into_evt(idx_rest)
    idx_walk_evt=grouping_consecutivePoints_into_evt(idx_walk)
    idx_groom_evt=grouping_consecutivePoints_into_evt(idx_groom)

    timesec_rest_evt=convert_idx_to_timepoint(idx_rest_evt, timesec)
    timesec_walk_evt=convert_idx_to_timepoint(idx_walk_evt, timesec)
    timesec_groom_evt=convert_idx_to_timepoint(idx_groom_evt, timesec)
    

    EthoIdxDic={}
    EthoTimeDic={}
    EthoColorCodeDic={}

    EthoIdxDic.update({'rest_evt':idx_rest_evt})
    EthoTimeDic.update({'rest_evt':timesec_rest_evt})
    EthoColorCodeDic.update({'rest_evt':plot_setting.rest_color})

    EthoIdxDic.update({'walk_evt':idx_walk_evt})
    EthoTimeDic.update({'walk_evt':timesec_walk_evt})
    EthoColorCodeDic.update({'walk_evt':plot_setting.walk_color})

    EthoIdxDic.update({'groom_evt':idx_groom_evt})
    EthoTimeDic.update({'groom_evt':timesec_groom_evt})
    EthoColorCodeDic.update({'groom_evt':plot_setting.groom_color})




    return EthoIdxDic, EthoTimeDic, EthoColorCodeDic


def FindLastCamPhotoIdx(cameraPath):

    listFiles = os.listdir(cameraPath)
    print('len listFiles',  len(listFiles))

    if '.DS_Store' in listFiles:
        print('removing .DS_Store')
        listFiles.remove('.DS_Store')

    if 'Thumbs.db' in listFiles:
        print('removing Thumbs.db')
        listFiles.remove('Thumbs.db')

    if 'capture_metadata.json' in listFiles:
        print('removing .json from img list')
        listFiles.remove('capture_metadata.json')

    # print('listFiles',listFiles)
    # print('len listFiles',  len(listFiles))
    
    #print('len(listFiles)%cameraAmount',len(listFiles)%cameraAmount)

    #camera_2_img_3220.jpg

    camera_num=[0,1,2,3,4,5,6]

    camera_0_imgs=[]
    camera_1_imgs=[]
    camera_2_imgs=[]
    camera_3_imgs=[]
    camera_4_imgs=[]
    camera_5_imgs=[]
    camera_6_imgs=[]
    for f in listFiles:
        #print(f)
        if f.startswith('camera_0'):
            camera_0_imgs.append(f)
        elif f.startswith('camera_1'):
            camera_1_imgs.append(f)
        elif f.startswith('camera_2'):
            camera_2_imgs.append(f)
        elif f.startswith('camera_3'):
            camera_3_imgs.append(f)
        elif f.startswith('camera_4'):
            camera_4_imgs.append(f)
        elif f.startswith('camera_5'):
            camera_5_imgs.append(f)
        elif f.startswith('camera_6'):
            camera_6_imgs.append(f)

    photoCounts_each_cam=[0,0,0,0,0,0,0]
    if len(camera_0_imgs)!=0:
        photoCounts_each_cam[0]=len(camera_0_imgs)
    if len(camera_1_imgs)!=0:
        photoCounts_each_cam[1]=len(camera_1_imgs)
    if len(camera_2_imgs)!=0:
        photoCounts_each_cam[2]=len(camera_2_imgs)
    if len(camera_3_imgs)!=0:
        photoCounts_each_cam[3]=len(camera_3_imgs)
    if len(camera_4_imgs)!=0:
        photoCounts_each_cam[4]=len(camera_4_imgs)
    if len(camera_5_imgs)!=0:
        photoCounts_each_cam[5]=len(camera_5_imgs)
    if len(camera_6_imgs)!=0:
        photoCounts_each_cam[6]=len(camera_6_imgs)


    print('photoCounts_each_cam', photoCounts_each_cam)

    print('np.nonzero(photoCounts_each_cam)[0]', np.nonzero(photoCounts_each_cam)[0])

    active_cam=[]
    for i in np.nonzero(photoCounts_each_cam)[0]:
        active_cam.append(i)
    active_cam_counts=len(active_cam)
    print('active_cam_counts', active_cam_counts)

    beh_photo_sum=sum(photoCounts_each_cam)
    print('beh_photo_sum', beh_photo_sum)

    if beh_photo_sum%active_cam_counts==0:
        vidStopIdx=int(beh_photo_sum/active_cam_counts-1)
    else:
        print('Behavioral cameras have acquisition issue. Amounts of photos are different from at least one of the cameras...')
        sys.exit(0)


    print('vidStopIdx', vidStopIdx)


    return vidStopIdx



def smooth_data(data, windowlen=10):

    if windowlen==0:
        smooth_data=data
    else:
        window = np.ones(int(windowlen))/float(windowlen)
        smooth_data=np.convolve(data,window,'same')

    return smooth_data


# def trim_period(trace, timeSec, startTime=0, endTime=120):

#     print('len trace', len(trace))

#     crop_startTime, crop_startIdx = find_nearest(timeSec, startTime)
#     crop_stopTime, crop_stopIdx = find_nearest(timeSec, endTime)

#     print('crop_startTime', crop_startTime, 'crop_startIdx', crop_startIdx)
#     print('crop_stopTime', crop_stopTime, 'crop_stopIdx', crop_stopIdx)


#     trace_crop = trace[crop_startIdx:crop_stopIdx]



#     return trace_crop



def fallingEdge(ArrayTD):

    zeroDet = np.where(ArrayTD==0)
    arrayZD = zeroDet[0]
    print ("arrayZD",arrayZD)
    checkDiff = np.diff(arrayZD.astype(int))
    print ("checkDiff",checkDiff)
    diff1 = np.where(checkDiff!=1)
    #print ("diff1",diff1)
    idxFall = []
    for i in range(len(diff1)):
        idxFall.append(arrayZD[diff1[i]+1])

    return idxFall




def beh_dir(date, genotype, fly, trial):
    return f"/mnt/data/CLC/{date - 20000000}_{genotype}/Fly{fly}/CO2xzGG/behData_{trial:03d}/"


def sync_dir(date, genotype, fly, trial):
    return f"/mnt/data/CLC/2photonData/CLC/{date}/{genotype}-fly{fly}/{genotype}-fly{fly}-sync-{trial:03d}"


def dir_2p(date, genotype, fly, trial):
    return f"/mnt/data/CLC/2photonData/CLC/{date}/{genotype}-fly{fly}/{genotype}-fly{fly}-{trial:03d}"


def clc_output_dir(date, genotype, fly, trial):
    return f"/mnt/data/CLC/Ascending_Project/{genotype.split('-')[0]}/2P/{date}/{genotype}-fly{fly}/{genotype}-fly{fly}-{trial:03d}/output"


def get_processed_sync_lines_Florian(date, genotype, fly, trial):
    dir_beh = beh_dir(date, genotype, fly, trial)
    dir_sync = sync_dir(date, genotype, fly, trial)
    dir2p = dir_2p(date, genotype, fly, trial)
    h5_path = utils2p.find_sync_file(dir_sync)
    co2_line, cam_line, opt_flow_line, frame_counter, capture_on = utils2p.synchronization.get_lines_from_h5_file(h5_path, ["CO2_Stim", "Basler", "OpFlow", "Frame Counter", "Capture On"])
    try:
        capture_json = utils2p.find_seven_camera_metadata_file(dir_beh)
    except FileNotFoundError:
        capture_json = None
    metadata_2p = utils2p.find_metadata_file(dir2p)
    metadata = utils2p.Metadata(metadata_2p)

    cam_line = utils2p.synchronization.process_cam_line(cam_line, capture_json)

    opt_flow_line = utils2p.synchronization.process_optical_flow_line(opt_flow_line)

    n_flyback_frames = metadata.get_n_flyback_frames()
    n_steps = metadata.get_n_z()
    frame_counter = utils2p.synchronization.process_frame_counter(frame_counter, steps_per_frame=n_flyback_frames + n_steps)

    co2_line = utils2p.synchronization.process_stimulus_line(co2_line)

    mask = np.logical_and(capture_on, frame_counter >= 0)
    mask = np.logical_and(mask, cam_line >= 0)

    co2_line, cam_line, opt_flow_line, frame_counter = utils2p.synchronization.crop_lines(mask, [co2_line, cam_line, opt_flow_line, frame_counter])

    optical_flow_path = utils2p.find_optical_flow_file(dir_beh)
    optical_flow = utils2p.load_optical_flow(optical_flow_path, 0, 0, 0, 0)

    # Ensure all optical flow data was saved and remove frames with missing data
    mask = (opt_flow_line < len(optical_flow["time_stamps"]))
    co2_line, cam_line, opt_flow_line, frame_counter = utils2p.synchronization.crop_lines(mask, [co2_line, cam_line, opt_flow_line, frame_counter])

    return co2_line, cam_line, opt_flow_line, frame_counter


def co2_regressors_Florian(frame_counter, co2_line):
    co2_onset = utils2p.synchronization.reduce_during_2p_frame(frame_counter, co2_line, lambda x: np.max(np.diff(x)))
    co2 = utils2p.synchronization.reduce_during_2p_frame(frame_counter, co2_line, np.mean)
    return co2_onset, co2



def get_PER_regressors_Florian(frame_counter, cam_line, date, genotype, fly, trial):

    pickle_file = os.path.join(clc_output_dir(date, genotype, fly, trial), "PER/camera_6/PER_labels_camera_6.p")
    with open(pickle_file, "rb") as f:
        PER = pickle.load(f)
    event = np.array(PER["evt_bin_trace"])
    length = np.array(PER["med_PER_exten_len"])

    event_regressor = utils2p.synchronization.reduce_during_2p_frame(frame_counter, event[cam_line], np.mean)
    length_regressor = utils2p.synchronization.reduce_during_2p_frame(frame_counter, length[cam_line], np.mean)
    return event_regressor, length_regressor
