import numpy as np
import os
import pickle
from itertools import groupby 
import sys
import h5py
import pandas as pd
import nrrd
from skimage import io
import matplotlib.pyplot as plt
plt.switch_backend('agg')


import utils.sync_utils as sync_utils
import utils.math_utils as math_utils



##main##

NAS_Dir='/mnt/data/'
NAS_AN_Proj_Dir=NAS_Dir+'CLC/Ascending_Project/'
workstation_dir='/mnt/internal_hdd/clc/'


NAS_Dir='/mnt/data/'
local_Dir = '/Users/clc/Documents/EPFL/NeLy/Data/ANproj/'
local_AN_proj_Dir = local_Dir + 'Ascending_Project/'
NAS_AN_Proj_Dir = NAS_Dir+'CLC/Ascending_Project/'

AN_Proj_Dir = NAS_AN_Proj_Dir





def desired_min_max_yaxis(trace, desired_min, desired_max):

    traceMin = round(min(np.asarray(trace).flatten()),1)
    traceMax = round(max(np.asarray(trace).flatten()),1)
    if desired_max != False:
        if abs(traceMax)<desired_max:
            traceMax=desired_max
        else:
            traceMax=np.nanmax(trace)
    else:
        traceMax = np.nanmax(trace)

    if desired_min !=False:
        if abs(traceMin)<abs(desired_min):
            traceMin=desired_min    
        else:
            traceMin = np.nanmin(trace)
    else:
        traceMin = np.nanmin(trace)

    return traceMax, traceMin



def group_expList_per_fly(list_exp):


    list_exp.sort() 
    # print('list_exp', list_exp)
    # print('shape list_exp', np.shape(list_exp))

    exp_list_per_fly = [list(i) for j, i in groupby(list_exp, lambda x: x[0]+x[1]+x[2])] 


    # print('exp_list_per_fly', exp_list_per_fly)
    # print('shape exp_list_per_fly', np.shape(exp_list_per_fly))


    
    return exp_list_per_fly


def import_camPhotos_into_stack(FilePath, cam_num=6):

    print('importing camera photos ... ')

    cam_photo_filename_init='camera_'+str(cam_num)+'_img_'

    cam_photos_stacks=[]
    for cnt in os.listdir(FilePath):
        if cnt.startswith(cam_photo_filename_init):
            # print('cnt', cnt)
            Cam_img = io.imread(FilePath+cnt)
            cam_photos_stacks.append(Cam_img)

    print('shape cam_photos_stacks', np.shape(cam_photos_stacks))


    return cam_photos_stacks




def readH5File(FilePath):

    myFile = h5py.File(FilePath,'r')
    Camera = myFile['DI']['Basler'][:].squeeze()
    OpFlow = myFile['DI']['OpFlow'][:].squeeze()
    if bool(myFile['DI']['CO2_Stim'])== True:
        CO2puff = myFile['DI']['CO2_Stim'][:].squeeze()
    FrameCounter = myFile['CI']['Frame Counter'][:].squeeze()
    
    print('FrameCounter[-1]', FrameCounter[-1])


    return Camera, OpFlow, CO2puff, FrameCounter



def readGCcsv(GC_dir, filename):

    print(GC_dir+filename+'.csv')


    if os.path.exists(GC_dir+filename+'.csv')==True:
        print(">>>Reading GC fiji...")
        print(GC_dir+filename+'.csv')
        file_GC = pd.read_csv(GC_dir+filename+'.csv')


    else:
        print('No GC .csv are found...')
        sys.exit(0)



    GC_set = [list(file_GC['Mean'].values.tolist())]

    baseline_len=int(10*2.43+0.5) #10s, 2.43fps

    dFF_set=math_utils.compute_dFF(GC_set, baseline_len)
    print('shape dFF_set', np.shape(dFF_set))

    dFF_series=dFF_set[0]



    return dFF_series






def readGCfile(GC_dir):

#This function opens fluorescence data.
    
    print(GC_dir+'dFF_dic.p')
    print(os.path.exists(GC_dir+'dFF_dic.p'))
    if os.path.exists(GC_dir+'dFF_dic.p')==True:
        print(">>>Reading GC_axoid...")
        print(GC_dir+'dFF_dic.p')
        file_dRR = pd.read_pickle(GC_dir+'dFF_dic.p')

    else:
        print('No GCamP data are found...')
        sys.exit(0)
    

    dRR_set = list(file_dRR.values())
    print('type dRR_set', type(dRR_set))
    print('shape dRR_set', np.shape(dRR_set))

    print('np.mean dRR_set[0]', np.mean(dRR_set[0]))

    dRR_set_c=dRR_set.copy()
    dRR_set_o=sync_utils.replace_nan_with_interp(dRR_set_c)

    print('type dRR_set', type(dRR_set))
    print('shape dRR_set', np.shape(dRR_set))

    print('np.mean dRR_set[0]', np.mean(dRR_set[0]))


    return dRR_set


def save_GC_dic_plot(path, dFF):

    def save_traces(traces, filename):
        """Pickle and save the traces under filename.pkl."""
        data_dic = {"ROI_" + str(i): traces[i] for i in range(len(traces))}
        with open(os.path.join(path, filename + ".p"), "wb") as f:
            pickle.dump(data_dic, f)


    save_traces(dFF, "dFF_dic")


    def plot_traces(traces, filename, ylabel, ylim=None):
            """Make a plot of the fluorescence traces."""

            traces=np.asarray(traces)
            if traces.size == 0:
                io.imsave(os.path.join(path, filename + ".png"),
                          np.ones((100,400,3), np.uint8) * 255)
                return
            
            color="forestgreen"
            if ylim is not None:
                ymin, ymax = ylim
            else:
                ymin = min(0, np.nanmin(traces[:, 1:]) - 0.05 * np.abs(np.nanmin(traces[:, 1:])))
                ymax = np.nanmax(traces[:, 1:]) * 1.05
            
            fig, axes = plt.subplots(len(traces), 1, 'all', 'all',
                                     figsize=(8, 4 * len(traces)), facecolor='white', dpi=300)
            fig.subplots_adjust(left=0.2, right = 0.9, wspace = 0.3, hspace = 0.3)
            
            # If only one ROI, still make an array of axes
            if not isinstance(axes, np.ndarray):
                axes = np.array([axes])
            for i, ax in enumerate(axes):
                ax.axhline(linestyle='dashed', color='gray', linewidth=0.5)
                ax.plot(traces[i], color, linewidth=1)
                ax.set_xlim(0, 1.05*(len(traces[i]) - 1))
                ax.set_ylabel("ROI#%d\n" % i + ylabel, size=10, color=color)
                ax.set_ylim(ymin, ymax)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                if i < len(traces) - 1:
                    ax.spines['bottom'].set_visible(False)
                    ax.get_xaxis().set_visible(False)
                else:
                    ax.set_xlabel("Frame", size=10)
                    
            fig.savefig(os.path.join(path, filename + ".png"),
                        bbox_inches='tight', facecolor=fig.get_facecolor(),
                        edgecolor='none', transparent=True)
        
    plot_traces(dFF, "ROIs_dFF", ylabel="$\Delta$F/F (%)")


    return



def read_absGC_file(GC_dir):

    #This function opens fluorescence abs data.
    
    print(GC_dir+'GC_abs_dic.p')
    print(os.path.exists(GC_dir+'GC_abs_dic.p'))
    if os.path.exists(GC_dir+'GC_abs_dic.p')==True:
        print(">>>Reading GC_axoid...")
        print(GC_dir+'GC_abs_dic.p')
        file_absF = pd.read_pickle(GC_dir+'GC_abs_dic.p')

    else:
        print('No GCamP data are found...')
        sys.exit(0)
    

    absF_set = list(file_absF.values())
    print('type absF_set', type(absF_set))
    print('shape absF_set', np.shape(absF_set))

    print('np.nanmean absF_set[0]', np.nanmean(absF_set[0]))

    absF_set_c=absF_set.copy()
    new_absF_set=[]
    for i, trace in enumerate(absF_set_c):
        new_trace=np.asarray(trace)*10000
        new_absF_set.append(new_trace)

    #absF_set_o=sync_utils.replace_nan_with_interp(new_absF_set)


    return new_absF_set





def read_nrrd(dir_nrrd, filename):

    #print('reading nrrd file in', dir_nrrd+filename)

    nrrd_data, header = nrrd.read(dir_nrrd+filename)

    #print('shape nrrd_data', np.shape(nrrd_data))


    return nrrd_data



def readBehRFTData(dirBehRFT, filename):
    
    if os.path.exists(dirBehRFT+"/"+filename):
        BehDicData = pickle.load( open( dirBehRFT+"/"+filename, "rb" ) )
        #BehDicData = pickle.load( open( outDirBeh_RFT+"/Manual_BehRFT_DicData.p", "rb" ) )
        Beh_RFT_Standlist = BehDicData['listStand']
        Beh_RFT_Walklist = BehDicData['listWalk']
        Beh_RFT_Groomlist = BehDicData['listGroom']
        Beh_RFT_Otherlist = BehDicData['listOther']
        #file_Beh_RFT = pd.read_pickle(dirBehRFT+"/"+filename)
    else:
        print ("File not found - BehRFT_DicData.p not done yet")
        sys.exit(0) 

    #Beh_RFT_set=list(file_Beh_RFT.values())



    return Beh_RFT_Standlist, Beh_RFT_Walklist, Beh_RFT_Groomlist, Beh_RFT_Otherlist
    #return Beh_RFT_set




def read_PER_data(PER_h5_Dir):

    if os.path.exists(PER_h5_Dir+'PER_labels_camera_6.p')==True:

        PER_dic = pickle.load( open( PER_h5_Dir+'PER_labels_camera_6.p', "rb" ) )

        pbsc0_X=PER_dic['pbsc0_X']
        pbsc0_Y=PER_dic['pbsc0_Y']      

        med_pbsc1_X=PER_dic['med_pbsc1_X']
        med_pbsc1_Y=PER_dic['med_pbsc1_Y']      

        evt_bin_trace=PER_dic['evt_bin_trace']
        
        med_PER_exten_len=PER_dic['med_PER_exten_len']
        norm_baseFold_med_PER_exten_len=PER_dic['norm_baseFold_med_PER_exten_len']

        evt_startIdx_list=PER_dic['evt_startIdx_list']
        evt_endIdx_list=PER_dic['evt_endIdx_list']


    else:
        print('DLC PER.p is not done ...')



    return PER_dic



def read_jointPos3d_Beh(outDir_J_Pos_Beh, filename):

    print('Opening', filename, '...')

    if os.path.exists(outDir_J_Pos_Beh+"/"+filename):
        J_Pos_Beh_DicData = pickle.load( open( outDir_J_Pos_Beh+"/"+filename, "rb" ) )


        # for i, v in J_Pos_Beh_DicData.items():

        #     print(i)
        #     print('type i', type(i))
        #     print('len v', len(v))
        #     print('type v', type(v))

        print('len J_Pos_Beh_DicData', len(J_Pos_Beh_DicData.items()))



    else:
        print ("File not found - "+filename + " not done yet")
        sys.exit(0)


    return J_Pos_Beh_DicData


def open_Beh_Jpos_GC_DicData(pathDic, filename):

    print('Opening',  filename ,'...')

    if os.path.exists(pathDic+'/'+filename):

        Beh_Jpos_GC_DicData = pickle.load(open( pathDic+'/'+filename, "rb" ))


        return Beh_Jpos_GC_DicData

    else:

        print('data not found...')
        sys.exit(0)

        return


def open_GCevt_dic(pathDic):

    print('Opening SyncDic_7CamBeh_GC-RES.p ...')

    if os.path.exists(pathDic+"/SyncDic_7CamBeh_GC-RES.p"):

        Beh_Jpos_GC_DicData = pickle.load(open( pathDic+"/SyncDic_7CamBeh_GC-RES.p", "rb" ))


    return Beh_Jpos_GC_DicData



def read_GCevt_dic(pathDic, filename):

    print('Opening GCevt_dic ...', filename)

    if os.path.exists(pathDic+'/'+filename):

        GCevt_dic = pickle.load(open( pathDic+'/'+filename, "rb" ))


    return GCevt_dic




def open_Beh_GC_DicData(pathDic):

    if os.path.exists(pathDic+"/DicDataForMovie_Beh_GC-RES.p"):
        DicData = pickle.load(open( pathDic+"/DicDataForMovie_Beh_GC-RES.p", "rb" ))
        frameCntr = DicData['frameCntr']
        GCset = DicData['GCset']

        CO2puff = DicData['CO2puff']
        Rest = DicData['rest']
        Walk = DicData['walk']
        Groom = DicData['groom']

        Etho_Idx_Dic = DicData['Etho_Idx_Dic']
        Etho_Timesec_Dic = DicData['Etho_Timesec_Dic']
        Etho_Colorcode_Dic = DicData['Etho_Colorcode_Dic']
        
        timeSec = DicData['timeSec']
        velForw_mm = DicData['velForw']
        velSide_mm = DicData['velSide']
        velTurn_deg = DicData['velTurn']


        startVidIdx = DicData['startVidIdx']
        stopVidIdx = DicData['stopVidIdx']
    else:
        print ("DicDataForMovie_Beh_GC-RES.p not found - Data not analysed yet")
        sys.exit(0)

    return frameCntr, GCset, Rest, Walk, Groom, \
            Etho_Idx_Dic, Etho_Timesec_Dic, Etho_Colorcode_Dic, \
            CO2puff, timeSec, velForw_mm, velSide_mm, velTurn_deg, \



def open_GC_beh_PER_sync_DicData():


    if os.path.exists(pathForDic+'/DicDataForMovie_GC_PER_trim_camera6_'+str(startFrame)+'f-'+str(endFrame)+'f-RES.p'):
        
        DicData = pickle.load(open( pathForDic+'/DicDataForMovie_GC_PER_trim_camera6_'+str(startFrame)+'f-'+str(endFrame)+'f-RES.p', 'rb' ))
        
        frameCntr = DicData['frameCntr']
        GCset = DicData['GCset']

        CO2puff = DicData['CO2puff']
        Rest = DicData['rest']
        Walk = DicData['walk']
        Groom = DicData['groom']

        Etho_Idx_Dic = DicData['Etho_Idx_Dic']
        Etho_Timesec_Dic = DicData['Etho_Timesec_Dic']
        Etho_Colorcode_Dic = DicData['Etho_Colorcode_Dic']
        
        timeSec = DicData['timeSec']
        velForw_mm = DicData['velForw']
        velSide_mm = DicData['velSide']
        velTurn_deg = DicData['velTurn']

        PER_bin_trace = DicData['PER_bin_trace']
        PER_extLen = DicData['PER_extLen']
        PER_norm_baseFold_extenLen = DicData['PER_norm_baseFold_extenLen']

        startVidIdx = DicData['startVidIdx']
        stopVidIdx = DicData['stopVidIdx']


    else:
        print(pathForDic)
        print ("DicDataForMovie_GC_PER_trim_camera6.p not found - Data not analysed yet")
        sys.exit(0)

    return frameCntr, GCset, Rest, Walk, Groom, \
            Etho_Idx_Dic, Etho_Timesec_Dic, Etho_Colorcode_Dic, \
            CO2puff, timeSec, velForw_mm, velSide_mm, velTurn_deg, \
            PER_bin_trace, PER_extLen, PER_norm_baseFold_extenLen, \
            startVidIdx, stopVidIdx



def find_corresponding_evt_from_groupIdxs(idxEthoDic, whichBeh, GCset, baseline=0, fps=0):

    # baseline is for including the previous trace as baseline for the corresponding epoch

    idx_beh_evt = idxEthoDic[whichBeh]
    bsl_len=int(baseline*fps)

    if len(idx_beh_evt)>0:
  
        # print('len GCset', len(GCset))
        # print('shape idx_beh_evt',  np.shape(idx_beh_evt))

        GCsetEvt=[]

        for i in range(0,len(GCset)):    
                GCsetEvt.append([])
                for j in range(0,len(idx_beh_evt)): 
                    if idx_beh_evt[j][0]-bsl_len>0:          
                        GCsetEvt[i].append(GCset[i][idx_beh_evt[j][0]-bsl_len:idx_beh_evt[j][-1]+1])

        # print('shape GCsetEvt', np.shape(GCsetEvt))
        # print('len GCsetEvt[0][0]', len(GCsetEvt[0][0]))


    else:
        GCsetEvt=[]
        for i in range(0,len(GCset)):
            GCsetEvt.append([])
            for j in range(0,len(idx_beh_evt)):           
                GCsetEvt[i].append(np.nan)



    # print('find_corresponding_evt_from_groupIdxs', 'shape GCsetEvt', np.shape(GCsetEvt))


    return GCsetEvt, bsl_len



def convert_3dList_into_3dnumpyArray(list3d_ori):

    list3d=list3d_ori.copy()

    np_array_3d=[]
    for ix, x in enumerate(list3d):
        np_array_3d.append([])
        for iy, y in enumerate(x): 
            np_array_3d[ix].append(np.array(y))

    np_array_3d=np.array([np.array(x) for x in np_array_3d])

    return np_array_3d


def fix_randomColorCode_for_eachObject(GCset):
    
    list_of_color=[]
    for i in range(0, len(GCset)):


        # choose a color and remove it from dict by .popitem()
        color_for_each_ROI=colors.popitem()[0]
        #color_for_each_ROI=random.choice(colors.keys())
        
        
        print('color_for_each_ROI',color_for_each_ROI)
        list_of_color.append(color_for_each_ROI)

        list_of_color=list(colors)


    return list_of_color



# def detect_event_by_diff(data, diff_thrsld, thrsld_facotor, timesec):
#     print("Dectecting events...")
#     print('data duration:', timesec[-1])
    
#     data_samplerate=len(data)/timesec[-1]

#     print('data_samplerate', data_samplerate)

     
#     intvl_dif=int(data_samplerate*sec_diff)
    
#     print("len data", len(data))    
#     ddata_series=[]
    
        
#     ##微分螢光序列
#     #1.微分smooth data
# #    wdwsize = np.int(np.floor(intvl_dif))
# #    smthwindow = np.ones(int(wdwsize))/float(wdwsize)
# #    smooth_data = np.convolve(data,smthwindow,'same')
# #    print("len smooth_data = ", len(smooth_data))
# #    '''
# #    nanseries=[]
# #    for i in range(0,len(smthwindow)):
# #        nanseries.append(np.nan)
# #    
# #    smooth_data=np.append(nanseries,smooth_data[0:len(smooth_data)-len(nanseries)])
# #    print("len smooth_data after= ", len(smooth_data))
# #    '''
    
# #    for i in range(0,len(smooth_data)-int(intvl_dif)):
# #        ddata = smooth_data[i+int(intvl_dif)]-smooth_data[i]
# #        ddata_series.append(ddata) 
    
#     #2.微分 original data    
#     k=0
#     for i in range(0,len(data)-intvl_dif):
#         ddatapoint = data[i+intvl_dif]-data[i]
#         ddata_series.append(ddatapoint) 
#         k+=1
    
#     #補缺的elements in 微分序列
#     for i in range(0,len(data)-len(ddata_series)):
#         ddata_series.append(np.nan)
    

    
    
#     ##finding the max in ddata_series and make idx list of them
#     evt_idx_series=[]
#     srch_window=[]   
#     find_closevalue_window=[]
#     for k in range(int(bsl_dur*data_samplerate), len(ddata_series)-int(evt_dur*data_samplerate)):
        
#         if ddata_series[k]>diff_thrsld:
#             #print("Tagging")
#             srch_window.append(ddata_series[k])
#             #print('srch_window',srch_window)
 

#         else:            
#             if len(srch_window)!= 0:
#                 #print("srch_window=",srch_window)
#                 maxddatainwindow = max(srch_window)
#                 print('maxddatainwindow', maxddatainwindow)
#                 initial_ddata = thrsld_facotor*maxddatainwindow
#                 maxddata_idx = k-(len(srch_window)-srch_window.index(maxddatainwindow))
#                 #idx_series.append(maxddata_idx)  
                               
#                 for i in range(int(maxddata_idx),int(maxddata_idx-bsl_dur*data_samplerate),-1):
#                     find_closevalue_window.append(abs(ddata_series[i]-initial_ddata))
                
#                 #print('find_closevalue_window', find_closevalue_window)
#                 #因為find_closevalue_window搜尋是反序的，所以用find_closevalue_window長度扣掉來順idx：
#                 minddata_idx = len(find_closevalue_window)-find_closevalue_window.index(min(find_closevalue_window))-1
#                 initialddata_idx = k-(len(find_closevalue_window)-minddata_idx)
#                 evt_idx_series.append(initialddata_idx)
                
#                 srch_window = []
#                 find_closevalue_window = []

#     #print('evt_idx_series', evt_idx_series)
    
#     return evt_idx_series, ddata_series, data_samplerate




def trim_row_of_list(list2d, startIdx=0, endIdx=-1, fps=1500, window_s=0):

    trim_2dlist=[]
    for row in list2d:
        #print('len row', len(row))
        if startIdx<len(row)-1:
            trim_row=row[startIdx:endIdx]
            trim_row= sync_utils.smooth_data(row[startIdx:endIdx], windowlen=int(fps*window_s))
            #print('len trim_row', len(trim_row))
            trim_2dlist.append(trim_row)

    return trim_2dlist



def downsampling_trace(trace, downsampling_freq):

    # print('trace', trace)
    # print('len trace', len(trace))
    # print('downsampling_freq', downsampling_freq)

    trace_DS=[]
    for i, v in enumerate(trace):
        if i%int(downsampling_freq)==0:
            trace_DS.append(trace[i])

    # print('trace_DS', trace_DS)

    return trace_DS




def flatten_list(l):


    flatten_list = [item for subl in l for item in subl]

    return flatten_list




