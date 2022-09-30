import cv2
import matplotlib.pyplot as plt
plt.switch_backend('agg')

import matplotlib.colors as mcolors
from matplotlib import cm
import numpy as np
import pandas as pd
import sys
import math
import itertools
import scipy
import scikit_posthocs
import matplotlib as mpl
from itertools import groupby
from scipy.stats import kde
import csv


import utils.general_utils as utils
import utils.math_utils as math_utils
import utils.plot_setting as plot_setting



#print(type(colors))
#print(colors)
CO2puff_color=plot_setting.CO2puff_color
GC_color=plot_setting.GC_color
AP_color=plot_setting.AP_color
ML_color=plot_setting.ML_color
Yaw_color=plot_setting.Yaw_color



## define duration of before-event period and event period
bsl_dur=4
evt_dur=12
##

dGC_thrsld = 47 # event threshold 閾值
dGC_thrsldfactor = 0.3    #dGC 最高值的幾％當作event起始點，即相當於AP kinx

dCO2thrsld=0.3
dCO2thrsldfactor = 1


sec_diff=0.01  #微分所使用的區間長度(秒）

evt_simlrty=0.7 #event similarity, for removing repetitive detected event. Range: 0~1. 1=totally overlapped

##
hst_bin=0.01

#opflow conversion factor
rot_to_mm_factor = 2*math.pi*5  #10mm ball
rot_to_deg_factor = 360

## Optic flow bout thrshold
AP_thrsld=1.2 #mm
ML_thrsld=1.2 #mm
Yaw_thrsld=17 #deg.
bout_len_thrsld=0.75 #s

last_event_idx_for_average=4




def defFlowFrameRange(timeSec, range_s=10):
    #This function returns the flow frame range which is the number of datapoint within 10 seconds of videos.

    xaxisRange=range_s # The length of time duration showing in each frame (s) 
    flowFrameRange = xaxisRange*len(timeSec)/(timeSec[-1]-timeSec[0]) #Convert duration to datapoint amount

    return xaxisRange, flowFrameRange


def find_nearest_nice_integer_of_postive_number(number):

    if np.floor(number/10) > 0:
        nearest_nice_integer=int(np.floor(number/10)*10)
    else:
        nearest_nice_integer=int(np.floor(number))

    return nearest_nice_integer



def find_aligned_pos_of_panelLabel(ref_xy_span=[1,1], target_xy_span=[2,2], ref_position=0.9, direction='vertical'):

    # This function  is for finding the position of text labels that put the text in the same place of panels with different size. 
    # It is aiming for the text alignment across different panels. 
    # ref_panel_span should be smaller than target_panel_span. Better use the smallest panel as the ref_panel_span. 



    if direction == 'vertical':

        ref_h = ref_xy_span[1]
        targ_h = target_xy_span[1]

        if targ_h>ref_h:

            ratio_h=targ_h/ref_h

            target_position=ref_position/ratio_h

        else:
            print('Target panel span should > Reference panel span.')
            sys.exit(0)


    elif direction == 'horizontal':

        ref_w = ref_xy_span[0]
        targ_w = target_xy_span[0]

        if targ_w>ref_w:

            ratio_w=targ_w/ref_w

            target_position=ref_position/ratio_w

        else:
            print('Target panel span should > Reference panel span.')
            sys.exit(0)

    #print('target_position', target_position)

    return target_position


def calc_scale_bar_um_to_px(pixelSize=0.5, scale_bar_length_um=10):

    scale_bar_px=scale_bar_length_um/pixelSize

    # print('scale_bar_length_um', scale_bar_length_um)
    # print('scale_bar_px', scale_bar_px)

    return scale_bar_px




def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n == -1:
        n = cmap.N
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(np.linspace(minval, maxval, n)))
    return new_cmap


    


def Plot_traces(trace_2DList, save_dir, filename, subtitle_list=None, plot_mode='row_by_row', xaxis_series=[]):

    if plot_mode=='row_by_row':
        row_num=len(trace_2DList)

        fig = plt.figure(facecolor='white', figsize=(20,row_num*5), dpi=170)

        for i, trace in enumerate(trace_2DList):
            plt.subplot(row_num,1,i+1)
            plt.title(subtitle_list[i])
            if len(xaxis_series)==0:
                plt.plot(trace) 
            elif len(xaxis_series)!=0:
                plt.plot(xaxis_series, trace) 
              
    elif plot_mode=='overlay':
        row_num=1

        fig = plt.figure(facecolor='white', figsize=(10,5), dpi=170)
        plt.title(subtitle_list)
        plt.subplot(1,1,1)

        for i, trace in enumerate(trace_2DList):
            plt.plot(trace)

    plt.tight_layout()
    plt.savefig(save_dir+filename)
    plt.clf()
    plt.close(fig)


    return








def Plot_Evtavg_overlay(timesecEthoDic, GCsetEvt_ori_beh, baseline=0, epoch_len=-1, y_lim=False, whichBeh='walk', filename='filename', filepath=None):

    print('\nPlotting ', whichBeh, ' GC event trace overlay....', )




    timesec_beh_evt = sorted(timesecEthoDic[whichBeh],key=len)

    dataFreq=(len(timesec_beh_evt[-1])/(timesec_beh_evt[-1][-1]-timesec_beh_evt[-1][0]))
    print('Longest event length in time (s)', timesec_beh_evt[-1][-1]-timesec_beh_evt[-1][0])
    print('dataFreq', dataFreq)

    GCsetEvt_ori_beh[0].sort(key=len)

    x_start = -baseline
    x_end = timesec_beh_evt[-1][-1]-timesec_beh_evt[-1][0]
    x_len = len(GCsetEvt_ori_beh[0][-1])
    #x_len = len(timesec_beh_evt[-1])+int(baseline*dataFreq)

    print('x_len', x_len)
    if epoch_len!=-1:
        epoch_len=int(epoch_len*dataFreq+baseline*dataFreq)




    # print('type timesec_beh_evt', type(timesec_beh_evt))
    # print('shape timesec_beh_evt', np.shape(timesec_beh_evt))
    # print('timesec_beh_evt[-1][0]', timesec_beh_evt[-1][0])
    # print('timesec_beh_evt[-1][-1]', timesec_beh_evt[-1][-1])






    Neuron_num=len(GCsetEvt_ori_beh)

    GCsetEvt=GCsetEvt_ori_beh.copy()

    print('shape GCsetEvt', np.shape(GCsetEvt))
    print('shape GCsetEvt_ori_beh', np.shape(GCsetEvt_ori_beh))
    

    GCsetEvt_for_plot=[]
    GCsetEvt_for_average=[]
    


    if len(timesec_beh_evt)>0:
        trace_dur = np.linspace(x_start, x_end, x_len)
    else:
        trace_dur = [np.nan]

    # print('trace_dur', trace_dur)
    # print('len trace_dur', len(trace_dur))

    row_width=4
    col_width=int(row_width*Neuron_num)-2

    fig = plt.figure(facecolor='white', figsize=(row_width,col_width), dpi=170)
    fig.subplots_adjust(left=0.3, right = 0.9, wspace = 0.3, hspace = 0.1)
    fig.suptitle(filename+'\n'+whichBeh, color='g')


    count_neuron=1





    for i in range(0, len(GCsetEvt_ori_beh)):

      
        
        if len(GCsetEvt[i])==0:
            print('Behavior ',whichBeh,'has 0 event in ROI# ', i ,'skip plot for this behavior')
            continue

        GCsetEvt_for_plot.append([])

        #fig.tight_layout(pad=0.4, w_pad=0.5,h_pad=1.0)
        
        #GCsetEvt_ori_beh[i]=sorted(GCsetEvt_ori_beh[i],key=len)
        #GCsetEvt_for_average[i]=sorted(GCsetEvt_for_average[i],key=len)

        GCsetEvt_ori_beh[i].sort(key=len)




            

        # print('len GCsetEvt_ori_beh[i][0]', len(GCsetEvt_ori_beh[i][0]))
        # print('len GCsetEvt_ori_beh[i][-1]', len(GCsetEvt_ori_beh[i][-1]))
        # print('GCsetEvt_ori_beh[i][0][-1]', GCsetEvt_ori_beh[i][0][-1])
        # print('GCsetEvt_ori_beh[i][-1][-1]', GCsetEvt_ori_beh[i][-1][-1])
        # print('len GCsetEvt_for_average[0][0]', len(GCsetEvt_for_average[0][0]))

        for j in range(0, len(GCsetEvt[i])):

         

            GCsetEvt_for_plot[i].append([])

            if len(GCsetEvt_ori_beh[i][j])!=x_len:
                nan_tail= [np.nan]*int(x_len-len(GCsetEvt_ori_beh[i][j]))
                #print('len(nan_tail)', len(nan_tail))
                #nan_tail= [0]*int(x_len-len(GCsetEvt[i][j]))
                GCevt_nantail=np.append(GCsetEvt[i][j],nan_tail)
                GCsetEvt_for_plot[i][j].extend(GCevt_nantail)
                

            else:
                print('event#',j,'equal to the longest length')
                GCsetEvt_for_plot[i][j].extend(GCsetEvt[i][j])

        
        # print('len GCsetEvt[i][0]', len(GCsetEvt[i][0]))
        # print('len GCsetEvt[i][-1]', len(GCsetEvt[i][-1]))
        print('shape GCsetEvt_for_plot', np.shape(GCsetEvt_for_plot))
        # print('shape GCsetEvt_for_average', np.shape(GCsetEvt_for_average))
        # print('GCsetEvt_for_plot[-1][0][-1]', GCsetEvt_for_plot[-1][0][-1])
        # print('GCsetEvt_for_plot[-1][-1][-1]', GCsetEvt_for_plot[-1][-1][-1])



        # looking for the last index with 4 nonnan values 
        find_early_idx_of_nan=[]
        if len(GCsetEvt_for_plot[i])>last_event_idx_for_average:
            for j in range(0 ,len(GCsetEvt_for_plot[i])):
                for k in range(0, len(GCsetEvt_for_plot[i][j])):
                    
                    if np.isnan(GCsetEvt_for_plot[i][j][k]):
                        find_early_idx_of_nan.append(k) 
                        break  
            print('find_early_idx_of_nan', find_early_idx_of_nan)
            print('len(find_early_idx_of_nan)', len(find_early_idx_of_nan))

            trace_has_no_nan=len(GCsetEvt_for_plot[i])-len(find_early_idx_of_nan)

            if len(find_early_idx_of_nan)==0:
                last_idx_for_average=-1
            elif last_event_idx_for_average>len(find_early_idx_of_nan) and trace_has_no_nan<last_event_idx_for_average:
                last_idx_for_average=1
            elif last_event_idx_for_average>len(find_early_idx_of_nan) and trace_has_no_nan>last_event_idx_for_average:
                last_idx_for_average=-1
            else:
                #last_idx_for_average=sorted(find_early_idx_of_nan)[-last_event_idx_for_average]
                last_idx_for_average=len(GCsetEvt[i][-last_event_idx_for_average])-1

        else:
            last_idx_for_average=-1


        print('last_idx_for_average', last_idx_for_average)


  

        # make GCsetEvt with at least 4 value for average 
        GCsetEvt_for_average.append([])
        for j in range(0 ,len(GCsetEvt_for_plot[i])):
            GCsetEvt_for_average[i].append([])
            for k in range(0, len(GCsetEvt_for_plot[i][j])):

                if k<last_idx_for_average or last_idx_for_average==-1:
                    GCsetEvt_for_average[i][j].append(np.array(GCsetEvt_for_plot[i][j][k]))
                else:
                    GCsetEvt_for_average[i][j].append(np.nan)

     

        GCsetEvt_for_average_cp=GCsetEvt_for_average.copy()

        GCsetEvt_for_average_np=utils.convert_3dList_into_3dnumpyArray(GCsetEvt_for_average_cp)
        


        print('shape GCsetEvt_for_average_np', np.shape(GCsetEvt_for_average_np))
        # print('shape GCsetEvt_for_average_np[0]', np.shape(GCsetEvt_for_average_np[0]))
        # print('shape GCsetEvt_for_average_np[0][0]', np.shape(GCsetEvt_for_average_np[0][0]))
        # print('type GCsetEvt_for_average_np', type(GCsetEvt_for_average_np))
        # print('type GCsetEvt_for_average_np[0]', type(GCsetEvt_for_average_np[0]))
        # print('type GCsetEvt_for_average_np[0][0]', type(GCsetEvt_for_average_np[0][0]))

        


        

        avg_GCtrace = np.nanmean(GCsetEvt_for_average_np[i],axis=0)
        print('len avg_GCtrace', len(avg_GCtrace))
        print('avg_GCtrace', avg_GCtrace)



        axGC= fig.add_subplot(Neuron_num,1,count_neuron)
        for j in range(0, len(GCsetEvt_for_plot[i])):
            #print('len GCsetEvt[i][j]', len(GCsetEvt[i][j]))
            axGC.plot(trace_dur[:epoch_len], GCsetEvt_for_plot[i][j][:epoch_len], color='g', label='ROI#'+str(i), linewidth=0.5, alpha = 0.3)
            
            # print('GCsetEvt[i][j]', GCsetEvt[i][j][-1])
            # print('type GCsetEvt[i][j]', type(GCsetEvt[i][j][-1]))
            # print('shape GCsetEvt[i][j]', np.shape(GCsetEvt[i][j][-1]))

            # #labeling the end of each evt with circle
            # for k in range(0, len(GCsetEvt_for_plot[i][j])):
            #     #print(GCsetEvt[i][j][k])
            #     if np.isnan(GCsetEvt_for_plot[i][j][k])==True:                   
            #         axGC.plot(trace_dur[k-1], GCsetEvt_for_plot[i][j][k-1], marker = 'o', markersize = '6', markeredgecolor='none', color='g', label='ROI#'+str(i), linewidth=0, alpha=0.3)
            #         break
            #     elif k==len(GCsetEvt_for_plot[i][j])-1:
            #         print('nan')
            #         axGC.plot(trace_dur[k], GCsetEvt_for_plot[i][j][k], marker = 'o', markersize = '6',  markeredgecolor='none', color='g', label='ROI#'+str(i), linewidth=0, alpha = 0.3)



        # print('len avg_GCtrace', len(avg_GCtrace))

        # print('shape GCsetEvt', np.shape(GCsetEvt))
        # print('shape GCsetEvt_for_average_np', np.shape(GCsetEvt_for_average_np))

        
        

        axGC.plot(trace_dur[:epoch_len],avg_GCtrace[:epoch_len], color='k', label='ROI#'+str(i), linewidth=1.25, alpha = 1)

        axGC.spines['bottom'].set_visible(False)
        axGC.spines['top'].set_visible(False)
        axGC.spines['right'].set_visible(False)
        axGC.get_xaxis().set_visible(False)
        #vars()[axGC].set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
        axGC.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
        axGC.axvline(0, linestyle='dashed',color='gray',linewidth=0.5)
        #axGC.set_xlabel('time (s)',size=10, color='k')
        axGC.set_ylabel('ROI#'+str(i)+'\n'+r'$\Delta$'+'F/F (%)',size=14,color='k')
        axGC.tick_params(axis='y', colors='k',top='off',labelsize=14)
        #axGC.set_ylim(minGC-0.1*(abs(minGC)),maxGC+0.1*(abs(maxGC)))
        if y_lim!=False:
            axGC.set_ylim(y_lim[0],y_lim[1])


        if i == len(GCsetEvt)-1:
            axGC.spines['bottom'].set_visible(True)
            axGC.get_xaxis().set_visible(True)
            axGC.tick_params(axis='x', colors='k',top='off',labelsize=14)
            axGC.set_xlabel('time (s)',size=14, color='k')

        count_neuron+=1

    
        # leg = plt.legend(loc='best', bbox_to_anchor = (1, 0.5), numpoints=1)
        # plt.savefig(outDirEvents + filename + '_ROI_' + str(i) +'.png', facecolor=fig.get_facecolor(), bbox_extra_artists=(leg,), bbox_inches='tight',edgecolor='none', transparent=True) #bbox_inches='tight', 
    plt.savefig(filepath + filename + '-' + whichBeh + '.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight',      
    plt.savefig(filepath + filename + '-' + whichBeh + '.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 

    plt.clf
    plt.close(fig)


    return





def Plot_Evtavg_overlay_err(timesecEthoDic, GCsetEvt_ori_beh, baseline=0, epoch_len=-1, y_lim=False, whichBeh='walk', filename='filename', filepath=None):

    print('y_lim', y_lim)

    print('\nPlotting ', whichBeh, ' GC event trace overlay with fill_between of the error....', )

    timesec_beh_evt = sorted(timesecEthoDic[whichBeh],key=len)
    least_realNum_amount=4

    dataFreq=(len(timesec_beh_evt[-1])/(timesec_beh_evt[-1][-1]-timesec_beh_evt[-1][0]))
    print('Longest event length in time (s)', timesec_beh_evt[-1][-1]-timesec_beh_evt[-1][0])
    print('dataFreq', dataFreq)

    GCsetEvt_ori_beh[0].sort(key=len)

    x_start = -baseline
    x_end = timesec_beh_evt[-1][-1]-timesec_beh_evt[-1][0]
    x_len = len(GCsetEvt_ori_beh[0][-1])
    #x_len = len(timesec_beh_evt[-1])+int(baseline*dataFreq)

    print('x_len', x_len)
    if epoch_len!=-1:
        epoch_len=int(epoch_len*dataFreq+baseline*dataFreq)




    # print('type timesec_beh_evt', type(timesec_beh_evt))
    # print('shape timesec_beh_evt', np.shape(timesec_beh_evt))
    # print('timesec_beh_evt[-1][0]', timesec_beh_evt[-1][0])
    # print('timesec_beh_evt[-1][-1]', timesec_beh_evt[-1][-1])






    Neuron_num=len(GCsetEvt_ori_beh)

    GCsetEvt=GCsetEvt_ori_beh.copy()
    

    GCsetEvt_for_plot=[]
    GCsetEvt_for_average=[]
    


    if len(timesec_beh_evt)>0:
        trace_dur = np.linspace(x_start, x_end, x_len)
    else:
        trace_dur = [np.nan]

    # print('trace_dur', trace_dur)
    # print('len trace_dur', len(trace_dur))

    row_width=4
    col_width=int(row_width*Neuron_num)-2

    fig = plt.figure(facecolor='white', figsize=(row_width,col_width), dpi=170)
    fig.subplots_adjust(left=0.3, right = 0.9, wspace = 0.3, hspace = 0.1)
    fig.suptitle(filename+'\n'+whichBeh, color='g')


    count_neuron=1
    for i in range(0, len(GCsetEvt_ori_beh)):
        
        if len(GCsetEvt[i])==0:
            print('Behavior ',whichBeh,'has 0 event in ROI# ', i ,'skip plot for this behavior')
            continue

        GCsetEvt_for_plot.append([])

        #fig.tight_layout(pad=0.4, w_pad=0.5,h_pad=1.0)
        
        #GCsetEvt_ori_beh[i]=sorted(GCsetEvt_ori_beh[i],key=len)
        #GCsetEvt_for_average[i]=sorted(GCsetEvt_for_average[i],key=len)

        GCsetEvt_ori_beh[i].sort(key=len)

        # if len(GCsetEvt_ori_beh[i])>least_realNum_amount:
        #     non_nan_avg_len=len(GCsetEvt_ori_beh[i][-least_realNum_amount])
        # else:
        #     non_nan_avg_len=epoch_len
        

        # print('len GCsetEvt_ori_beh[i][0]', len(GCsetEvt_ori_beh[i][0]))
        # print('len GCsetEvt_ori_beh[i][-1]', len(GCsetEvt_ori_beh[i][-1]))
        # print('GCsetEvt_ori_beh[i][0][-1]', GCsetEvt_ori_beh[i][0][-1])
        # print('GCsetEvt_ori_beh[i][-1][-1]', GCsetEvt_ori_beh[i][-1][-1])
        # print('len GCsetEvt_for_average[0][0]', len(GCsetEvt_for_average[0][0]))

        for j in range(0, len(GCsetEvt[i])):

            GCsetEvt_for_plot[i].append([])

            if len(GCsetEvt_ori_beh[i][j])!=x_len:
                nan_tail= [np.nan]*int(x_len-len(GCsetEvt_ori_beh[i][j]))
                #print('len(nan_tail)', len(nan_tail))
                #nan_tail= [0]*int(x_len-len(GCsetEvt[i][j]))
                GCevt_nantail=np.append(GCsetEvt[i][j],nan_tail)
                GCsetEvt_for_plot[i][j].extend(GCevt_nantail)
                

            else:
                print('event#',j,'equal to the longest length')
                GCsetEvt_for_plot[i][j].extend(GCsetEvt[i][j])

        
        # print('len GCsetEvt[i][0]', len(GCsetEvt[i][0]))
        # print('len GCsetEvt[i][-1]', len(GCsetEvt[i][-1]))
        print('shape GCsetEvt_for_plot', np.shape(GCsetEvt_for_plot))
        # print('shape GCsetEvt_for_average', np.shape(GCsetEvt_for_average))
        # print('GCsetEvt_for_plot[-1][0][-1]', GCsetEvt_for_plot[-1][0][-1])
        # print('GCsetEvt_for_plot[-1][-1][-1]', GCsetEvt_for_plot[-1][-1][-1])

        # looking for the last index with 4 nonnan values 
        find_early_idx_of_nan=[]
        if len(GCsetEvt_for_plot[i])>least_realNum_amount:
            for j in range(0 ,len(GCsetEvt_for_plot[i])):
                for k in range(0, len(GCsetEvt_for_plot[i][j])):
                    if np.isnan(GCsetEvt_for_plot[i][j][k]):
                        find_early_idx_of_nan.append(k) 
                        break  

            trace_has_no_nan=len(GCsetEvt_for_plot[i])-len(find_early_idx_of_nan)


            if len(find_early_idx_of_nan)==0:
                last_idx_for_average=-1
            elif least_realNum_amount>len(find_early_idx_of_nan) and trace_has_no_nan<least_realNum_amount:
                last_idx_for_average=1
            elif least_realNum_amount>len(find_early_idx_of_nan) and trace_has_no_nan>least_realNum_amount:
                last_idx_for_average=-1
            else:
                #last_idx_for_average=sorted(find_early_idx_of_nan)[-last_event_idx_for_average]
                last_idx_for_average=len(GCsetEvt[i][-least_realNum_amount])-1

        else:
            last_idx_for_average=-1




            





        # make GCsetEvt with at least 4 value for average 
        GCsetEvt_for_average.append([])
        for j in range(0 ,len(GCsetEvt_for_plot[i])):
            GCsetEvt_for_average[i].append([])
            for k in range(0, len(GCsetEvt_for_plot[i][j])):
                if k<last_idx_for_average:
                    GCsetEvt_for_average[i][j].append(np.array(GCsetEvt_for_plot[i][j][k]))
                else:
                    GCsetEvt_for_average[i][j].append(np.nan)

        GCsetEvt_for_average_cp=GCsetEvt_for_average.copy()

        GCsetEvt_for_average_np=utils.convert_3dList_into_3dnumpyArray(GCsetEvt_for_average_cp)
    
        print('shape GCsetEvt_for_average_np', np.shape(GCsetEvt_for_average_np))
        # print('shape GCsetEvt_for_average_np[0]', np.shape(GCsetEvt_for_average_np[0]))
        # print('shape GCsetEvt_for_average_np[0][0]', np.shape(GCsetEvt_for_average_np[0][0]))
        # print('type GCsetEvt_for_average_np', type(GCsetEvt_for_average_np))
        # print('type GCsetEvt_for_average_np[0]', type(GCsetEvt_for_average_np[0]))
        # print('type GCsetEvt_for_average_np[0][0]', type(GCsetEvt_for_average_np[0][0]))

        
        
        avg_GCtrace, downCI_GCtrace, upCI_GCtrace=math_utils.compute_CI_and_mean_trace(GCsetEvt_for_average_np[i], confidence=0.95, least_realNum_amount=least_realNum_amount)

        evt_id_list=[]
        for evt_id, evt_v in enumerate(GCsetEvt_for_average_np[i]):
            evt_id_list.append('evt_'+str(evt_id))




        # print('len avg_GCtrace', len(avg_GCtrace))

        # print('shape GCsetEvt', np.shape(GCsetEvt))
        # print('shape GCsetEvt_for_average_np', np.shape(GCsetEvt_for_average_np))
        # print('downCI_GCtrace', downCI_GCtrace)
        # print('upCI_GCtrace', upCI_GCtrace)

        axGC= fig.add_subplot(Neuron_num,1,count_neuron)
        

        axGC.plot(trace_dur[:epoch_len],avg_GCtrace[:epoch_len], color='k', label='ROI#'+str(i), linewidth=1.25, alpha = 1)
        axGC.fill_between(trace_dur[:epoch_len], downCI_GCtrace[:epoch_len], upCI_GCtrace[:epoch_len], color='k', linewidth=0, alpha=plot_setting.err_shade_alpha)

        axGC.spines['bottom'].set_visible(False)
        axGC.spines['top'].set_visible(False)
        axGC.spines['right'].set_visible(False)
        axGC.get_xaxis().set_visible(False)
        #vars()[axGC].set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
        axGC.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
        axGC.axvline(0, linestyle='dashed',color='gray',linewidth=0.5)
        #axGC.set_xlabel('time (s)',size=10, color='k')
        axGC.set_ylabel('ROI#'+str(i)+'\n'+r'$\Delta$'+'F/F (%)',size=14,color='k')
        axGC.tick_params(axis='y', colors='k',top='off',labelsize=14)
        #axGC.set_ylim(minGC-0.1*(abs(minGC)),maxGC+0.1*(abs(maxGC)))

        print('max(trace_dur[:epoch_len])', max(trace_dur[:epoch_len]))
        if max(trace_dur[:epoch_len])<1:
            axGC.set_xlim(-1, 1)
        else:
            axGC.set_xlim(-1, None)
        
        # axGC.set_ylim(y_lim[0],y_lim[1])
        if y_lim!=False:
            print()
            print('setting y limit ...', y_lim)
            print()
            axGC.set_ylim(y_lim[0],y_lim[1])



        if i == len(GCsetEvt)-1:
            axGC.spines['bottom'].set_visible(True)
            axGC.get_xaxis().set_visible(True)
            axGC.tick_params(axis='x', colors='k',top='off',labelsize=14)
            axGC.set_xlabel('time (s)',size=14, color='k')

        count_neuron+=1

        with open(filepath+filename+'-'+whichBeh+'-'+'ROI_'+str(i)+'.csv', 'w') as f:
      
            # using csv.writer method from CSV package
            write = csv.writer(f)
              
            write.writerow(['ROI_'+str(i)])
            write.writerow(['Time (s)']+evt_id_list)
            
            write.writerows(zip(trace_dur[:epoch_len], *GCsetEvt_for_average_np[i]))
    
        # leg = plt.legend(loc='best', bbox_to_anchor = (1, 0.5), numpoints=1)
        # plt.savefig(outDirEvents + filename + '_ROI_' + str(i) +'.png', facecolor=fig.get_facecolor(), bbox_extra_artists=(leg,), bbox_inches='tight',edgecolor='none', transparent=True) #bbox_inches='tight', 
    plt.savefig(filepath + filename + '-' + whichBeh + '_CI.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight',      
    plt.savefig(filepath + filename + '-' + whichBeh + '_CI.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 

    plt.clf
    plt.close(fig)








    return



def plot_matrix(x_list, y_list, y_value_matr, second_x_list=False, Gal4_x_list_reformat=False, roi_seperation_marker='-ROI#', figsize=(10,2.5), set_max=False, set_min=False, r2=None, savedir=None, title=None, PlotMethod=None, unit='a.u.', cmap='RdBu', hatch=False, hatchcolor='k'):


    print('Plotting matrix ....')

    # print('shape y_value_matr', np.shape(y_value_matr))

    joint_num=36
    dof_joint_leg=6

    leg_color=['r', 'b', 'g', 'y', 'm','c']
    joint_angle_name=['θ', 'Ψ', 'Φ', 'α', 'β', 'γ']

    data_font_size=plot_setting.data_font_size-10

    x_series = []
    for x in range(0,len(x_list)):
        x_series.append(x)

    print('len x_series', len(x_series))

    print('shape y_value_matr before transpose', np.shape(y_value_matr))

    y_series=[]
    #print('len(np.shape(y_value_matr))', len(np.shape(y_value_matr)))
    if len(np.shape(y_value_matr))>1: 
        for y in range(0,len(y_value_matr[0])): 
            y_series.append(y)
        y_value_matr=np.asarray(y_value_matr).transpose()



    else:
        y_value_matr_new=[]
        print('only 1d, reshaping...')
        for y in range(0,len(y_value_matr)):
            y_value_matr_new.append([y_value_matr[y]])
        y_value_matr=np.asarray(y_value_matr_new).transpose()
        # y_value_matr=y_value_matr_new
        y_series=[0]

    print('y_series', y_series)

    print('shape y_value_matr after transpose', np.shape(y_value_matr))


    
    if PlotMethod=='R^2':
        vmin_mthd=0
        vmax_mthd=1
        index='R^2'
        cbar_tick=[0,1]
        color_nthd='Blues'

        figsize=(10, 2.5) 
        dpi=300
        data_font_size=plot_setting.data_font_size-25


    elif PlotMethod=='Correlation_coefficient':
        vmin_mthd=-1
        vmax_mthd=1
        index='Corr. Coef.'
        cbar_tick=[-1,-0.5,0,0.5,1]
        color_nthd='RdBu'

        figsize=(10, 2.5) 
        dpi=300
        data_font_size=plot_setting.data_font_size-25


    elif PlotMethod=='weight_of_linear_regression':
        vmin_mthd=np.nanmin(np.asarray(y_value_matr))*(-1)
        #vmin_mthd=np.nanmin(np.asarray(y_value_matr))
        vmax_mthd=np.nanmax(np.asarray(y_value_matr))
        index='Weight'
        cbar_tick=None
        color_nthd=cmap

        figsize=(10, 2.5) 
        dpi=300
        data_font_size=plot_setting.data_font_size-25

    elif PlotMethod=='p_value':
        # print('y_value_matr\n', y_value_matr)
        # print('type y_value_matr\n', type(y_value_matr))

        # y_value_matr[np.isnan(y_value_matr)] = 4
        y_value_matr=np.where(y_value_matr>0.05, 4, y_value_matr)
        y_value_matr=np.where(y_value_matr<0.001, 1, y_value_matr)
        y_value_matr=np.where(y_value_matr<0.01, 2, y_value_matr)
        y_value_matr=np.where(y_value_matr<0.05, 3, y_value_matr)
        y_value_matr=4-y_value_matr
        
        # print('y_value_matr\n', y_value_matr)

        mat_p_value=pd.DataFrame(y_value_matr)
        mat_p_value=mat_p_value.set_axis(x_list, axis='columns')
        mat_p_value=mat_p_value.set_axis(y_list, axis='index')

        colormap = mcolors.ListedColormap(['#FFAB91', '#E1AAFA', '#C05EED', '#9D00E5'])
        colorbar_ticklabels = ['NS', 'p < 0.05', 'p < 0.01', 'p < 0.001']

        fig= plt.figure(facecolor='white', figsize=(80, 40), dpi=170)
        ax = plt.subplot(1,1,1)
        # ms = ax.matshow(y_value_matr, cmap = colormap, vmin=y_value_matr.min() - 0.5, vmax=y_value_matr.max() + 0.5, origin = 'lower')

        ms = ax.matshow(y_value_matr, cmap = colormap, vmin=np.nanmin(y_value_matr), vmax=np.nanmax(y_value_matr), origin = 'upper')

        ax.xaxis.tick_bottom() 
        ax.set_xticks(x_series, minor = False)
        ax.set_yticks(y_series, minor = False)
        ax.set_xticklabels([str(xx) for xx in x_list], rotation=90, size=data_font_size)
        ax.set_yticklabels([str(yy) for yy in y_list], size=data_font_size)
        cbar = fig.colorbar(ms,ticks = np.arange(0.001, 1))
        cbar.ax.set_yticklabels(colorbar_ticklabels, size=data_font_size)

        plt.savefig(str(savedir+title+'.png'), facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
        plt.savefig(str(savedir+title+'.pdf'), facecolor=fig.get_facecolor(), edgecolor='none', transparent=True)  
        plt.clf()
        plt.close(fig)     

        return
    
    elif PlotMethod=='dFF_01_w_negative':
        vmin_mthd=np.nanmin(np.asarray(y_value_matr))
        vmax_mthd=np.nanmax(np.asarray(y_value_matr))
        # if abs(vmax_mthd)>=abs(vmin_mthd):
        #     vmin_mthd=-vmax_mthd
        # else:
        #     vmax_mthd=-vmin_mthd

        if vmin_mthd>0:
            vmin_mthd=0

        index=unit
        cbar_tick=None

        minColor = 0.15
        maxColor = 0.5
        color_nthd = truncate_colormap(plt.get_cmap("brg"), minColor, maxColor)
        # color_nthd='Greens'

        ##---------------------------------------------
        # if abs(vmax_mthd)>=abs(vmin_mthd):
        #     vmin_mthd=-vmax_mthd
        # else:
        #     vmax_mthd=-vmin_mthd        
        # color_nthd='bwr'
        ##---------------------------------------------

        # color_nthd='gray'
        # color_nthd='brg'
        # color_nthd='plasma'

        
        dpi=300
        data_font_size=plot_setting.data_font_size-25
        

               

    else:
        if set_min!=False:
            vmin_mthd=set_min
        else:
            vmin_mthd=np.nanmin(np.asarray(y_value_matr))

        if set_max!=False:
            vmax_mthd=set_max
        else:
            vmax_mthd=np.nanmax(np.asarray(y_value_matr))

        if vmin_mthd>0:
            vmin_mthd=0
        index=unit
        cbar_tick=None
        color_nthd=cmap

        
        dpi=300
        data_font_size=plot_setting.data_font_size-25


    frame_dim = 18
    fig = plt.figure(facecolor='white', figsize=figsize, dpi=dpi)
    fig.suptitle(title, fontsize=data_font_size)
    fig.subplots_adjust(wspace = 0.01, hspace=0.01, left=0.25, right = 0.99, bottom = 0.33 , top = 0.86)
    #ax = plt.subplots(1,1)
    ax = plt.subplot2grid((frame_dim, frame_dim),(0, 0),rowspan=frame_dim,colspan=frame_dim)


    # cmap_bas = plt.get_cmap(color_nthd)

    default_cbar=plt.get_cmap(cmap,256)
    newcolors = default_cbar(np.linspace(0, 1, 256))
    white = np.array([256/256, 256/256, 256/256, 1])
    pink = np.array([248/256, 24/256, 148/256, 1])
    newcolors[:3, :] = white
    
    ## Choose default colorbar
    #newcmp = default_cbar
    ## Choose customized colorbar
    from matplotlib.colors import ListedColormap
    newcmp = ListedColormap(newcolors)

    minColor = 0
    maxColor = 1    
    color_nthd_bas = truncate_colormap(newcmp, minColor, maxColor)
    cmap_bas = plt.get_cmap(color_nthd_bas) 

    countROI_each_Gal4 = count_ROIeachGal4_from_roiList(x_list, split_mark=roi_seperation_marker)

    if hatch==False:
        x=np.arange(len(x_series)+1)
        y=np.arange(len(y_series)+1)*(-1)

        current_x=0
        x_Gal4_seperation_tick=[0]
        for i, roi_n in enumerate(countROI_each_Gal4):
            current_x+=(roi_n)
            x_Gal4_seperation_tick.append(x[current_x])

        psm = ax.pcolor(x, y, y_value_matr, cmap=cmap_bas, edgecolor=None, vmin=vmin_mthd, vmax=vmax_mthd)
        ax.set_xticks(x+0.5)
        ax.set_xticks(x_Gal4_seperation_tick, minor=True)
        ax.tick_params(axis='x', which='minor', colors='w')
        ax.xaxis.grid(True,  which='minor', color='#C4C4C4', linewidth=0.5)

        ax.set_yticks(y-0.5)

        ax.set_xlim(min(x), max(x))
        ax.set_ylim(min(y), max(y))
        ax.set_xticklabels(x_list, rotation=90, size=data_font_size)
        ax.set_yticklabels(y_list, size=data_font_size, color='k')  

    else:
        x=np.arange(len(x_series)+1)
        y=np.arange(len(y_series)+1)*(-1)
        psm = ax.pcolor(x, y, y_value_matr, cmap=cmap_bas, hatch=hatch, edgecolor=None, vmin=vmin_mthd, vmax=vmax_mthd)
        ax.set_xticks(x+0.5)
        current_x=0
        x_Gal4_seperation_tick=[0]
        for i, roi_n in enumerate(countROI_each_Gal4):
            current_x+=(roi_n)
            x_Gal4_seperation_tick.append(x[current_x])
        

        ax.set_yticks(y-0.5)

        ax.set_xlim(min(x), max(x))
        ax.set_ylim(min(y), max(y))

        
        ax.set_yticklabels(y_list, size=10, color='k')  


    cbar=fig.colorbar(psm, ax=ax, shrink=0.3, ticks=cbar_tick)
    cbar.ax.tick_params(labelsize=data_font_size)
    #cbar.ax.set_yticklabels(['-1','','-0.5','','0', '','0.5','','1'], size=data_font_size)

    cbar.set_label(index, rotation=270, size=data_font_size)


    if Gal4_x_list_reformat==True:
        neuronid_x_list=[]
        Gal4_x_list=[]
        for i, roi in enumerate(x_list):
            # print('roi.split(roi_seperation_marker)', roi.split(roi_seperation_marker))
            Gal4=roi.split(roi_seperation_marker)[0]
            neuron_id=roi.split(roi_seperation_marker)[1]
            Gal4_x_list.append(Gal4)
            neuronid_x_list.append(neuron_id)
        Gal4_x_list=list(dict.fromkeys(Gal4_x_list))

        ax.set_xticks(x_Gal4_seperation_tick, minor=True)
        ax.set_xticklabels(neuronid_x_list, rotation=0, size=4)


        ## Plot Gal4 axis
        print('len x_Gal4_seperation_tick', len(x_Gal4_seperation_tick))
        ax2 = ax.twiny()
        ax2.xaxis.set_ticks_position("bottom")
        ax2.xaxis.set_label_position("bottom")
        ax2.spines["bottom"].set_position(("axes", -0.3))
        ax2.set_xticks(x_Gal4_seperation_tick, minor=True)
        ax2.tick_params(axis='x', which='major', bottom=False, labelbottom=False)
        ax2.tick_params(axis='x', which='minor', colors='k', direction='in')



        ax3 = ax.twiny()
        ax3.xaxis.set_ticks_position("bottom")
        ax3.xaxis.set_label_position("bottom")
        ax3.spines["bottom"].set_position(("axes", -0.31))
        ax3.get_xaxis().set_visible(True)
        ax3.spines['bottom'].set_visible(False)
        ax3.tick_params(axis='x', which='major', bottom=False, labelbottom=False)
        ax3.tick_params(axis='x', which='minor', colors='k', bottom=False, labelbottom=True)

        x_Gal4_tick=[]
        for i in range(0, len(x_Gal4_seperation_tick)-1):
            Gal4_pos=(x_Gal4_seperation_tick[i]+x_Gal4_seperation_tick[i+1])/2
            print('Gal4_pos', Gal4_pos)
            x_Gal4_tick.append(Gal4_pos)
        x_Gal4_tick.append(x_Gal4_seperation_tick[-1])
        ax3.set_xticks(x_Gal4_tick, minor=True) 
        ax3.set_xticklabels(Gal4_x_list, minor=True, rotation=90, size=6)
        print('len x_Gal4_tick', len(x_Gal4_tick))
        





    else:
        ax.set_xticklabels(x_list, rotation=90, size=4)
        ax.set_yticklabels(y_list, size=6, color='k')  



    if second_x_list!=False:
        ax4=ax.twiny()
        ax4.set_xticks(x+0.5)
        ax4.set_xlim(min(x), max(x))
        ax4.set_xticklabels(second_x_list, rotation=90, size=5)


 

    plt.savefig(str(savedir+title+'.png'), facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.savefig(str(savedir+title+'.pdf'), facecolor=fig.get_facecolor(), edgecolor='none', transparent=True)  
    plt.clf()
    plt.close(fig)


    # print('x_list', x_list)
    print('y_value_matr', y_value_matr)




    return



def count_ROIeachGal4_from_roiList(roiList, split_mark='-'):

    new_roi_list=[]
    for i, roi in enumerate(roiList):
        roi_new=roi.split(split_mark)
        new_roi_list.append(roi_new)

    groupedGal4_roi = [list(v) for i, v in groupby(new_roi_list, lambda x: x[0])]

    countROI_each_Gal4=[]
    for i, subgroup in enumerate(groupedGal4_roi):
        roi_count=len(subgroup)
        countROI_each_Gal4.append(roi_count)


    # print('groupedGal4_roi', groupedGal4_roi)
    # print('countROI_each_Gal4', countROI_each_Gal4)


    return countROI_each_Gal4




def plot_overlay_matrix(x_list, y_list, mat_bas, mat_mask, savedir=None, title=None, colorbar_bas='gist_heat', colorbar_mask='binary', unit='a.u.', hatch='///////', hatchcolor='w'):

    # print('shape mat_bas', np.shape(mat_bas))
    # print('shape mat_mask', np.shape(mat_bas))

    print('Plotting overlay matrix ....')



    data_font_size=plot_setting.data_font_size-25

    x_series = []
    for x in range(0,len(x_list)):
        x_series.append(x)

    print('len x_series', len(x_series))

    y_series=[]
    #print('len(np.shape(mat_bas))', len(np.shape(mat_bas)))
    if len(np.shape(mat_bas))>1: 
        for y in range(0,len(mat_bas[0])): 
            y_series.append(y)
        mat_bas=np.asarray(mat_bas).transpose()


    else:
        y_value_matr_new=[]
        print('only 1d, reshaping...')
        for y in range(0,len(mat_bas)):
            y_series.append(y)
            y_value_matr_new.append([y])
        mat_bas=np.asarray(y_value_matr_new).transpose()


    if len(np.shape(mat_mask))>1: 
        mat_mask=np.asarray(mat_mask).transpose()

    else:
        y_value_matr_new=[]
        print('only 1d, reshaping...')
        for y in range(0,len(mat_mask)):
            y_value_matr_new.append([y])
        mat_mask=np.asarray(y_value_matr_new).transpose()



    countROI_each_Gal4 = count_ROIeachGal4_from_roiList(x_list)


    


    # print('y_series', y_series)

    # print('shape mat_bas after transpose', np.shape(mat_bas))
    # print('shape mat_mask after transpose', np.shape(mat_mask))

    # print('mat_mask after transpose',mat_mask)

    # vmin_mthd_bas=np.nanmin(np.asarray(mat_bas))
    vmin_mthd_bas=np.min(mat_bas[np.nonzero(mat_bas)])
    vmax_mthd_bas=np.nanmax(np.asarray(mat_bas))
    # if abs(vmax_mthd)>=abs(vmin_mthd):
    #     vmin_mthd=-vmax_mthd
    # else:
    #     vmax_mthd=-vmin_mthd


    if vmin_mthd_bas>0:
        vmin_mthd_bas=0


    vmin_mthd_mask=np.nanmin(np.asarray(mat_mask))
    vmax_mthd_mask=np.nanmax(np.asarray(mat_mask))

    # print('vmin_mthd_bas', vmin_mthd_bas)
    # print('vmax_mthd_bas', vmax_mthd_bas)

    # print('vmin_mthd_mask', vmin_mthd_mask)
    # print('vmax_mthd_mask', vmax_mthd_mask)



    index=unit
    cbar_tick=None




    frame_dim = 18
    fig = plt.figure(facecolor='white', figsize=(10, 3.5), dpi=300)
    fig.suptitle(title, fontsize=data_font_size+16)
    fig.subplots_adjust(wspace = 0.01, hspace=0.01, left=0.25, right = 0.99, bottom = 0.2 , top = 0.86)
    ax = plt.subplot2grid((frame_dim, frame_dim),(0, 0),rowspan=frame_dim,colspan=frame_dim)

    # import matplotlib.ticker as ticker
    # from mpl_toolkits.axes_grid.parasite_axes import SubplotHost

    # ax = SubplotHost(fig, 111)
    # fig.add_subplot(ax, rowspan=frame_dim,colspan=frame_dim)

    default_cbar=plt.get_cmap(colorbar_bas,256)
    newcolors = default_cbar(np.linspace(0, 1, 256))
    white = np.array([256/256, 256/256, 256/256, 1])
    pink = np.array([248/256, 24/256, 148/256, 1])
    newcolors[:2, :] = white
    
    # newcmp = default_cbar
    from matplotlib.colors import ListedColormap
    newcmp = ListedColormap(newcolors)

    minColor = 0
    maxColor = 1    
    color_nthd_bas = truncate_colormap(newcmp, minColor, maxColor)
    cmap_bas = plt.get_cmap(newcmp)


    # minColor = 0
    # maxColor = 1
    # color_nthd_bas = truncate_colormap(plt.get_cmap(colorbar_bas), minColor, maxColor)
    # cmap_bas = plt.get_cmap(color_nthd_bas)
    # # cmap_bas = plt.get_cmap(colorbar_bas)

    # cmap_mask = cmap_bas.reversed()
    cmap_mask = plt.get_cmap('binary').reversed()

    # cmap_bas = plt.get_cmap(colorbar_bas)
    # cmap_mask = plt.get_cmap(colorbar_mask)
    # hatch_color=cmap_bas(0.25)

    mpl.rc('hatch', linewidth=0.4)
    mpl.rc('hatch', color='gray')

   

    x=np.arange(len(x_series)+1)
    y=np.arange(len(y_series)+1)*(-1)

    current_x=0
    x_roi_seperation_tick=[0]
    for i, roi_n in enumerate(countROI_each_Gal4):
        current_x+=(roi_n)
        x_roi_seperation_tick.append(x[current_x])

    psm_bas = ax.pcolor(x, y, mat_bas, cmap=cmap_bas, edgecolor=None, vmin=vmin_mthd_bas, vmax=vmax_mthd_bas)
    psm_mask = ax.pcolor(x, y, mat_mask, hatch=hatch, cmap=cmap_mask, edgecolor=None, vmin=vmin_mthd_bas, vmax=vmax_mthd_mask)
    ax.set_xticks(x+0.5)
    ax.set_xticks(x_roi_seperation_tick, minor=True)
    ax.tick_params(axis='x', which='major', colors='k', direction='out', length=2)
    ax.tick_params(axis='x', which='minor', colors='w', direction='out', length=6)
    ax.tick_params(axis='y', which='major', colors='k', direction='out', length=2)
    ax.xaxis.grid(True,  which='minor', color='#C4C4C4', linewidth=0.5, alpha=1)

    ax.set_yticks(y-0.5)

    ax.set_xlim(min(x), max(x))
    ax.set_ylim(min(y), max(y))

    neuronid_x_list=[]
    Gal4_x_list=[]
    for i, roi in enumerate(x_list):
        Gal4=roi.split('#')[0][:-4]
        neuron_id=roi.split('#')[1]
        Gal4_x_list.append(Gal4)
        neuronid_x_list.append(neuron_id)
    ax.set_xticklabels(neuronid_x_list, rotation=0, size=5)
    ax.set_yticklabels(y_list, size=6, color='k')  

    cbar=fig.colorbar(psm_bas, ax=ax, shrink=0.3, ticks=cbar_tick)
    cbar.ax.tick_params(labelsize=data_font_size)
    #cbar.ax.set_yticklabels(['-1','','-0.5','','0', '','0.5','','1'], size=data_font_size)

    cbar.set_label(index, rotation=270, size=data_font_size)


    ax2 = ax.twiny()
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax2.spines["bottom"].set_position(("axes", -0.085))
    ax2.set_xticks(x_roi_seperation_tick, minor=True)
    ax2.tick_params(axis='x', which='major', bottom=False, labelbottom=False)
    ax2.tick_params(axis='x', which='minor', colors='k', direction='in')

    ax3 = ax.twiny()
    ax3.xaxis.set_ticks_position("bottom")
    ax3.xaxis.set_label_position("bottom")
    ax3.spines["bottom"].set_position(("axes", -0.08))
    print('len x_roi_seperation_tick', len(x_roi_seperation_tick))
    x_Gal4_seperation_tick=[]
    for i in range(0, len(x_roi_seperation_tick)-1):
        Gal4_pos=(x_roi_seperation_tick[i]+x_roi_seperation_tick[i+1])/2
        # print('Gal4_pos', Gal4_pos)
        x_Gal4_seperation_tick.append(Gal4_pos)
    x_Gal4_seperation_tick.append(x_roi_seperation_tick[-1])
    print('len x_Gal4_seperation_tick', len(x_Gal4_seperation_tick))
    ax3.get_xaxis().set_visible(True)
    ax3.spines['bottom'].set_visible(False)
    ax3.set_xticks(x_Gal4_seperation_tick, minor=True)
    ax3.tick_params(axis='x', which='major', bottom=False, labelbottom=False)
    ax3.tick_params(axis='x', which='minor', colors='k', bottom=False, labelbottom=True)
    Gal4_x_list=list(dict.fromkeys(Gal4_x_list))
    # print(Gal4_x_list)
    ax3.set_xticklabels(Gal4_x_list, minor=True, rotation=90, size=5)







    



    # offset = 0, -25 # Position of the second axis
    # new_axisline = ax2.get_grid_helper().new_fixed_axis
    # ax2.axis["bottom"] = new_axisline(loc="bottom", axes=ax2, offset=offset)
    # ax2.axis["top"].set_visible(False)
    # ax2.set_xticks(x_Gal4_seperation_tick/max(x_Gal4_seperation_tick))
    # ax2.xaxis.set_major_formatter(ticker.NullFormatter())
    # ax2.xaxis.set_minor_locator(ticker.FixedLocator([0.3, 0.8]))
    # ax2.xaxis.set_minor_formatter(ticker.FixedFormatter(['mammal', 'reptiles']))
 

    plt.savefig(str(savedir+title+'.png'), facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.savefig(str(savedir+title+'.pdf'), facecolor=fig.get_facecolor(), edgecolor='none', transparent=True)  
    plt.clf()
    plt.close(fig)




    return



def plot_overlay_innerv_matrix(x_list, y_list, mat_bas, mat_mask1, mat_mask2, savedir=None, title=None, colorbar_bas='gist_heat', colorbar_mask='binary', unit='a.u.', hatch='///////', hatchcolor='w'):

    # print('shape mat_bas', np.shape(mat_bas))
    # print('shape mat_mask', np.shape(mat_bas))

    print('Plotting overlay matrix ....')






    data_font_size=plot_setting.data_font_size-25

    x_series = []
    for x in range(0,len(x_list)):
        x_series.append(x)

    print('len x_series', len(x_series))

    y_series=[]
    #print('len(np.shape(mat_bas))', len(np.shape(mat_bas)))
    if len(np.shape(mat_bas))>1: 
        for y in range(0,len(mat_bas[0])): 
            y_series.append(y)
        mat_bas=np.asarray(mat_bas).transpose()
    else:
        y_value_matr_new=[]
        print('only 1d, reshaping...')
        for y in range(0,len(mat_bas)):
            y_series.append(y)
            y_value_matr_new.append([y])
        mat_bas=np.asarray(y_value_matr_new).transpose()


    if len(np.shape(mat_mask1))>1: 
        mat_mask1=np.asarray(mat_mask1).transpose()
    else:
        y_value_matr_new=[]
        print('only 1d, reshaping...')
        for y in range(0,len(mat_mask1)):
            y_value_matr_new.append([y])
        mat_mask1=np.asarray(y_value_matr_new).transpose()


    if len(np.shape(mat_mask2))>1: 
        mat_mask2=np.asarray(mat_mask2).transpose()
    else:
        y_value_matr_new=[]
        print('only 1d, reshaping...')
        for y in range(0,len(mat_mask2)):
            y_value_matr_new.append([y])
        mat_mask2=np.asarray(y_value_matr_new).transpose()



    countROI_each_Gal4 = count_ROIeachGal4_from_roiList(x_list)


    


    # print('y_series', y_series)

    # print('shape mat_bas after transpose', np.shape(mat_bas))
    # print('shape mat_mask after transpose', np.shape(mat_mask))

    # print('mat_mask after transpose',mat_mask)

    vmin_mthd_bas=np.nanmin(np.asarray(mat_bas))
    vmax_mthd_bas=np.nanmax(np.asarray(mat_bas))
    # if abs(vmax_mthd)>=abs(vmin_mthd):
    #     vmin_mthd=-vmax_mthd
    # else:
    #     vmax_mthd=-vmin_mthd

    if vmin_mthd_bas>0:
        vmin_mthd_bas=0


    vmin_mthd_mask1=np.nanmin(np.asarray(mat_mask1))
    vmax_mthd_mask1=np.nanmax(np.asarray(mat_mask1))
    if vmin_mthd_mask1==vmax_mthd_mask1:
        vmin_mthd_mask1=vmin_mthd_mask1-vmax_mthd_mask1

    vmin_mthd_mask2=np.nanmin(np.asarray(mat_mask2))
    vmax_mthd_mask2=np.nanmax(np.asarray(mat_mask2))
    if vmin_mthd_mask2==vmax_mthd_mask2:
        vmin_mthd_mask2=vmin_mthd_mask2-vmax_mthd_mask2

    # print('vmin_mthd_bas', vmin_mthd_bas)
    # print('vmax_mthd_bas', vmax_mthd_bas)

    # print('vmin_mthd_mask', vmin_mthd_mask)
    # print('vmax_mthd_mask', vmax_mthd_mask)

    # print('vmin_mthd_mask2', vmin_mthd_mask2)
    # print('vmax_mthd_mask2', vmax_mthd_mask2)



    index=unit
    cbar_tick=None



    # figsize=(10, 2.5)
    figsize=(10, 5)
    frame_dim = 18
    fig = plt.figure(facecolor='white', figsize=(10, 10), dpi=300)
    fig.suptitle(title, fontsize=data_font_size+16)
    fig.subplots_adjust(wspace = 0.01, hspace=0.01, left=0.1, right = 0.99, bottom = 0.5 , top = 0.86)
    #ax = plt.subplots(1,1)
    ax = plt.subplot2grid((frame_dim, frame_dim),(0, 0),rowspan=frame_dim,colspan=frame_dim)


    default_cbar=plt.get_cmap(colorbar_bas,256)
    newcolors = default_cbar(np.linspace(0, 1, 256))
    white = np.array([256/256, 256/256, 256/256, 1])
    pink = np.array([248/256, 24/256, 148/256, 1])
    newcolors[:2, :] = white
    
    # newcmp = default_cbar
    from matplotlib.colors import ListedColormap
    newcmp = ListedColormap(newcolors)

    minColor = 0
    maxColor = 1    
    color_nthd_bas = truncate_colormap(newcmp, minColor, maxColor)
    cmap_bas = plt.get_cmap(newcmp)

    cmap_mask1 = mcolors.ListedColormap(['k', 'g'])
    cmap_mask2 = cmap_bas.reversed()

    # cmap_bas = plt.get_cmap(colorbar_bas)
    # cmap_mask = plt.get_cmap(colorbar_mask)
    # hatch_color=cmap_bas(0.25)



    

    x=np.arange(len(x_series)+1)
    y=np.arange(len(y_series)+1)*(-1)

    current_x=0
    x_roi_seperation_tick=[0] 
    for i, roi_n in enumerate(countROI_each_Gal4):
        current_x+=(roi_n)
        x_roi_seperation_tick.append(x[current_x])  

    psm_bas = ax.pcolor(x, y, mat_bas, cmap=cmap_bas, edgecolor=None, vmin=vmin_mthd_bas, vmax=vmax_mthd_bas)
    mpl.rc('hatch', linewidth=0.5)
    mpl.rc('hatch', color='steelblue')
    psm_mask1 = ax.pcolor(x, y, mat_mask1, hatch='.......', cmap=cmap_mask2, edgecolor=None, vmin=vmin_mthd_mask1, vmax=vmax_mthd_mask1)
    mpl.rc('hatch', color='gray')
    psm_mask2 = ax.pcolor(x, y, mat_mask2, hatch=hatch, cmap=cmap_mask2, edgecolor=None, vmin=vmin_mthd_mask2, vmax=vmax_mthd_mask2)

    ax.set_xticks(x+0.5)
    ax.set_xticks(x_roi_seperation_tick, minor=True)
    ax.tick_params(axis='x', which='major', colors='k', direction='out', length=2)
    ax.tick_params(axis='x', which='minor', colors='w', direction='out', length=6)
    ax.tick_params(axis='y', which='major', colors='k', direction='out', length=2)
    ax.xaxis.grid(True,  which='major', color='k', linewidth=0.2, alpha=0)
    ax.xaxis.grid(True,  which='minor', color='#C4C4C4', linewidth=0.5)

    ax.set_yticks(y-0.5)

    ax.set_xlim(min(x), max(x))
    ax.set_ylim(min(y), max(y))

    neuronid_x_list=[]
    Gal4_x_list=[]
    for i, roi in enumerate(x_list):
        Gal4=roi.split('#')[0][:-4]
        neuron_id=roi.split('#')[1]
        Gal4_x_list.append(Gal4)
        neuronid_x_list.append(neuron_id)

    ax.set_xticklabels(neuronid_x_list, rotation=0, size=6)
    ax.set_yticklabels(y_list, size=6, color='k')  

    cbar=fig.colorbar(psm_bas, ax=ax, shrink=0.3, ticks=cbar_tick)
    cbar.ax.tick_params(labelsize=data_font_size)
    #cbar.ax.set_yticklabels(['-1','','-0.5','','0', '','0.5','','1'], size=data_font_size)

    cbar.set_label(index, rotation=270, size=data_font_size)



    ax2 = ax.twiny()
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax2.spines["bottom"].set_position(("axes", -0.21))
    ax2.set_xticks(x_roi_seperation_tick, minor=True)
    ax2.tick_params(axis='x', which='major', bottom=False, labelbottom=False)
    ax2.tick_params(axis='x', which='minor', colors='k', direction='in')

    ax3 = ax.twiny()
    ax3.xaxis.set_ticks_position("bottom")
    ax3.xaxis.set_label_position("bottom")
    ax3.spines["bottom"].set_position(("axes", -0.20))
    print('len x_roi_seperation_tick', len(x_roi_seperation_tick))
    x_Gal4_seperation_tick=[]
    for i in range(0, len(x_roi_seperation_tick)-1):
        Gal4_pos=(x_roi_seperation_tick[i]+x_roi_seperation_tick[i+1])/2
        # print('Gal4_pos', Gal4_pos)
        x_Gal4_seperation_tick.append(Gal4_pos)
    x_Gal4_seperation_tick.append(x_roi_seperation_tick[-1])
    print('len x_Gal4_seperation_tick', len(x_Gal4_seperation_tick))
    ax3.get_xaxis().set_visible(True)
    ax3.spines['bottom'].set_visible(False)
    ax3.set_xticks(x_Gal4_seperation_tick, minor=True)
    ax3.tick_params(axis='x', which='major', bottom=False, labelbottom=False)
    ax3.tick_params(axis='x', which='minor', colors='k', bottom=False, labelbottom=True)
    Gal4_x_list=list(dict.fromkeys(Gal4_x_list))
    # print(Gal4_x_list)
    ax3.set_xticklabels(Gal4_x_list, minor=True, rotation=90, size=5)
     

    plt.savefig(str(savedir+title+'.png'), facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.savefig(str(savedir+title+'.pdf'), facecolor=fig.get_facecolor(), edgecolor='none', transparent=True)  
    plt.clf()
    plt.close(fig)




    return









def Plot_whole_PER_trace(GCset, CO2puff, \
    PER_bin_trace, PER_extLen_px, PER_norm_baseFold_extenLen,\
    timeSec, velForw, velSide, velTurn, \
    ethoTimeDic, ethoColorDic, filename):



  
    Num_ROI = len(GCset)
    print('len(GCset)', len(GCset))
    
    Totalnum_row = int(Num_ROI+5) # + 3 optic flow traces (AP, ML, Yaw) + 1 CO2 + 1 PER
    rowspan = int(Totalnum_row/Totalnum_row)

    

    print('Totalnum_row', Totalnum_row)
    print('rowspan', rowspan)


    trace_dur = np.linspace(0,timeSec[-1]-timeSec[0],len(CO2puff[:]))
    print('len trace_dur',len(trace_dur))

    fig = plt.figure(facecolor='white', figsize=(Totalnum_row+math.ceil(2*Totalnum_row),Totalnum_row), dpi=90)
    fig.subplots_adjust(left=0.15, right = 0.9, wspace = 0.3, hspace = 0.3)


    row_count=rowspan

    ##plot detected event summary##
    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axCO2 = plt.subplot(Totalnum_row,1,row_count)
    #plt.title(key)    
    axCO2.plot(timeSec, CO2puff[:], color=CO2puff_color,linewidth=1)
    
    axCO2.spines['bottom'].set_visible(False)
    axCO2.spines['top'].set_visible(False)
    axCO2.spines['right'].set_visible(False)
    axCO2.get_xaxis().set_visible(False)
    axCO2.get_yaxis().set_label_coords(-0.1,0.5)
    axCO2.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    #ax1.axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)            
    #axCO2.set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
    axCO2.set_ylim(0,1)
    axCO2.yaxis.set_ticks(np.array([0,1]))
    axCO2.set_ylabel(r'$\rm{CO}_\mathrm{2}$' +' puff',size=10,color=CO2puff_color)
    #plt.savefig(outDirbouts + ROI + '_'+'GCsum'+ '.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    #plt.clf
    row_count+=rowspan

    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axPER_px = plt.subplot(Totalnum_row,1,row_count)
    axPER_px.plot(timeSec, PER_extLen_px[:], color=PER_color,linewidth=1)
    axPER_px.spines['bottom'].set_visible(False)
    axPER_px.spines['top'].set_visible(False)
    axPER_px.spines['right'].set_visible(False)
    axPER_px.get_xaxis().set_visible(False)
    axPER_px.get_yaxis().set_label_coords(-0.1,0.5)
    axPER_px.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    axPER_px.set_ylim(0, max_PER_len)
    #axPER_px.yaxis.set_ticks(np.array([min_PER_len,max_PER_len]))
    axPER_px.set_ylabel('PER length\n(px)',size=10,color=PER_color)


    row_count+=rowspan
    

    ##plot GC evt summary##
    
    for i in range(0, Num_ROI):
        axGC = 'axGC' + str(i)

    #vars()[axGC]

    # print('len trace_dur', len(trace_dur))
    # print('trace_dur[-1]',trace_dur[-1])
    # print('len GCset[i]', len(GCset[i]))
        print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
        vars()[axGC] = plt.subplot(Totalnum_row,1,row_count)
        vars()[axGC].plot(timeSec, GCset[i][:], color=GC_color,linewidth=1)
        
        vars()[axGC].spines['bottom'].set_visible(False)
        vars()[axGC].spines['top'].set_visible(False)
        vars()[axGC].spines['right'].set_visible(False)
        vars()[axGC].get_xaxis().set_visible(False)
        vars()[axGC].get_yaxis().set_label_coords(-0.1,0.5)
        #vars()[axGC].set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
        vars()[axGC].set_ylim(minGC,maxGC)
        vars()[axGC].axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
        #vars()[axGC].axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)
        vars()[axGC].set_ylabel('ROI#'+str(i)+'\n'+r'$\Delta$'+'F/F (%)',size=10,color=GC_color)
        row_count+=rowspan




    ##plot event AP flow summary##
    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axAP = plt.subplot(Totalnum_row,1,row_count)
    axAP.plot(timeSec,velForw[:], color=AP_color,linewidth=1)
    
    axAP.spines['bottom'].set_visible(False)
    axAP.spines['top'].set_visible(False)
    axAP.spines['right'].set_visible(False)

    axAP.get_xaxis().set_visible(False)
    axAP.get_yaxis().set_label_coords(-0.1,0.5)
    axAP.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    #axAP.axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)
    axAP.set_ylabel(r'$\rm{V}_\mathrm{forward}$'+'\n'+r'$\rm{(mm s}^\mathrm{-1}$'+')',size=10,color=AP_color)
    #axAP.set_ylabel(r"$\mathrm{\mathsf{v_{forward}}}$",fontsize=15,color='r',rotation=90,horizontalalignment='left')
    #axAP.set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
    axAP.set_ylim(velForwMin,velForwMax)
    tempListCurAP = [velForwMin,velForwMax]
    newYAxisAP = np.array(tempListCurAP)
    axAP.yaxis.set_ticks(newYAxisAP)
    row_count+=rowspan

    #plt.savefig(outDirbouts + ROI + '_'+'APflowsum'+ '.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    #plt.clf

    ##plot event ML flow summary## 
    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axML = plt.subplot(Totalnum_row,1,row_count)
    axML.plot(timeSec,velSide[:], color=ML_color,linewidth=1)
    
    axML.spines['bottom'].set_visible(False)
    axML.spines['top'].set_visible(False)
    axML.spines['right'].set_visible(False)
    axML.get_xaxis().set_visible(False)
    axML.get_yaxis().set_label_coords(-0.1,0.5)
    axML.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    #axML.axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)
    axML.set_ylabel(r'$\rm{V}_\mathrm{side}$'+'\n'+r'$\rm{(mm s}^\mathrm{-1}$'+')',size=10,color=ML_color)
    #axML.set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
    axML.set_ylim(velSideMin,velSideMax)
    tempListCurML = [velSideMin,velSideMax]
    newYAxisML = np.array(tempListCurML)
    axML.yaxis.set_ticks(newYAxisML)
    row_count+=rowspan

       
    #plt.savefig(outDirbouts + ROI + '_'+'MLflowsum'+ '.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    #plt.clf

    ##plot event Yaw flow summary##
    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axYaw = plt.subplot(Totalnum_row,1,row_count)  
    axYaw.plot(timeSec,velTurn[:], color=Yaw_color,linewidth=1)

    # axYaw.spines['bottom'].set_visible(False)
    axYaw.spines['top'].set_visible(False)
    axYaw.spines['right'].set_visible(False)
    axYaw.get_yaxis().set_label_coords(-0.1,0.5)
    # axYaw.get_xaxis().set_visible(False)
    axYaw.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    #axYaw.axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)
    axYaw.set_xlabel('Time (s)',size=10,color='k')
    axYaw.set_ylabel(r'$\rm{V}_\mathrm{turn}$'+'\n'+r'$\rm{(deg. s}^\mathrm{-1}$'+')',size=10,color=Yaw_color)
    #axYaw.set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
    axYaw.set_ylim(velTurnMin,velTurnMax)
    tempListCurYaw = [velTurnMin,velTurnMax]
    newYAxisYaw = np.array(tempListCurYaw)
    axYaw.yaxis.set_ticks(newYAxisYaw)
    row_count+=rowspan


    
    ## put each ethography in each whole trace

    for key in ethoTimeDic: 
        
        print('plotting ', key, ' wholetrace...')    

        for a in range(0,len(ethoTimeDic[key])):

            axCO2.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)

            for i in range(0, Num_ROI):
                axGC = 'axGC' + str(i)
                vars()[axGC].axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)

            axAP.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)
            axML.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)
            axYaw.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)  
            axPER_px.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)  


    plt.savefig(outDirEvents + 'whole_trace_beh.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    plt.savefig(outDirEvents + 'whole_trace_beh.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    plt.clf()
    plt.close()      

    
    
    return
  



def plot_2D_density_contour(pnts_2D, ref_reg_pnts_2D=False, x_dim=False, y_dim=False, filename='density_contour', filepath=None):


    nbins=100

    x=np.asarray(list(zip(*pnts_2D))[0])
    y=np.asarray(list(zip(*pnts_2D))[1])

    x=np.append(x,[[0, x_dim], [0, x_dim]])
    y=np.append(y,[[0, y_dim], [y_dim, 0]])

    print('x.min(), x.max()', x.min(), x.max())
    print('y.min(), y.max()', y.min(), y.max())
    print('x_dim', x_dim)
    print('y_dim', y_dim)



    data=[x,y]
    k = kde.gaussian_kde(data)
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    # xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))


    fig = plt.figure(facecolor='white', figsize=(8.27, 3.88), dpi=300)
    fig.suptitle(filename, fontsize=8, fontname='Arial', x=0.26)

    default_cbar=plt.get_cmap('OrRd',256)
    newcolors = default_cbar(np.linspace(0, 1, 256))
    white = np.array([256/256, 256/256, 256/256, 1])
    pink = np.array([248/256, 24/256, 148/256, 1])
    newcolors[:26, :] = white
    
    # newcmp = default_cbar
    from matplotlib.colors import ListedColormap
    newcmp = ListedColormap(newcolors)

    minColor = 0
    maxColor = 1    
    color_nthd_bas = truncate_colormap(newcmp, minColor, maxColor)
    cmap_bas = plt.get_cmap(newcmp)


    ax0 = plt.subplot(1,1,1)
    ax0.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=cmap_bas)
    ax0.contour(xi, yi, zi.reshape(xi.shape), 5 , linewidths=0.5)



    ax0.set_xlim(0, x_dim)
    ax0.set_ylim(0, y_dim)
    ax0.invert_yaxis()







    plt.savefig(filepath + filename+'.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    plt.savefig(filepath + filename+'.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    # plt.savefig(filepath + filename+'.svg', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 

    plt.clf()
    plt.close()      
   


    return



def Plot_whole_trace_off_ball(GCset, CO2puff, timeSec, ethoTimeDic, filename='whole_trace_off_ball', filepath=None):

    Num_ROI = len(GCset)
    Totalnum_row = int(Num_ROI+1) #  + 1 CO2 
    rowspan = int(Totalnum_row/Totalnum_row)

    maxGC, minGC=utils.desired_min_max_yaxis(np.asarray(GCset).flatten(), desired_min=False, desired_max=100)

    trace_dur = np.linspace(0,timeSec[-1]-timeSec[0],len(CO2puff[:]))
    print('len trace_dur',len(trace_dur))

    fig = plt.figure(facecolor='white', figsize=(Totalnum_row+math.ceil(2*Totalnum_row),Totalnum_row), dpi=150)
    fig.subplots_adjust(left=0.15, right = 0.9, bottom=0.2, wspace = 0.3, hspace = 0.3)


    row_count=rowspan

    ##plot detected event summary##
    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axCO2 = plt.subplot(Totalnum_row,1,row_count)
    #plt.title(key)    
    axCO2.plot(timeSec, CO2puff[:], color=CO2puff_color,linewidth=1)
    
    axCO2.spines['bottom'].set_visible(False)
    axCO2.spines['top'].set_visible(False)
    axCO2.spines['right'].set_visible(False)
    axCO2.get_xaxis().set_visible(False)
    axCO2.get_yaxis().set_label_coords(-0.1,0.5)
    axCO2.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    #ax1.axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)            
    #axCO2.set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
    axCO2.set_ylim(0,1)
    axCO2.yaxis.set_ticks(np.array([0,1]))
    axCO2.set_ylabel(r'$\rm{CO}_\mathrm{2}$' +' puff',size=10,color=CO2puff_color)
    #plt.savefig(outDirbouts + ROI + '_'+'GCsum'+ '.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    #plt.clf
    row_count+=rowspan

    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    

    ##plot GC evt summary##
    
    for i in range(0, Num_ROI):
        axGC = 'axGC' + str(i)

    #vars()[axGC]

    # print('len trace_dur', len(trace_dur))
    # print('trace_dur[-1]',trace_dur[-1])
    # print('len GCset[i]', len(GCset[i]))
        print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
        vars()[axGC] = plt.subplot(Totalnum_row,1,row_count)
        vars()[axGC].plot(timeSec, GCset[i][:], color=GC_color,linewidth=1)
        
        vars()[axGC].spines['bottom'].set_visible(False)
        vars()[axGC].spines['top'].set_visible(False)
        vars()[axGC].spines['right'].set_visible(False)
        vars()[axGC].get_xaxis().set_visible(False)
        vars()[axGC].get_yaxis().set_label_coords(-0.1,0.5)
        #vars()[axGC].set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
        vars()[axGC].set_ylim(minGC,maxGC)
        vars()[axGC].axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
        #vars()[axGC].axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)
        vars()[axGC].set_ylabel('ROI#'+str(i)+'\n'+r'$\Delta$'+'F/F (%)',size=10,color=GC_color)
        row_count+=rowspan

        if i==Num_ROI-1:
            vars()[axGC].spines['bottom'].set_visible(True)
            vars()[axGC].get_xaxis().set_visible(True)
            vars()[axGC].set_xlabel('Time (s)',size=10,color='k')


    for key in ethoTimeDic: 


        if key=='CO2puff_evt':

            continue

        elif key=='move_evt' or key=='rest_evt':

            print('plotting ', key, ' wholetrace...', plot_setting.EthoColorDic[key])    

            for a in range(0,len(ethoTimeDic[key])):

                for i in range(0, Num_ROI):
                    axGC = 'axGC' + str(i)
                    vars()[axGC].axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)

     

    plt.savefig(filepath + filename+'.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    plt.savefig(filepath + filename+'.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    plt.clf()
    plt.close()      



    return




def Plot_whole_trace(GCset, CO2puff,\
    timeSec, velForw, velSide, velTurn, \
    ethoTimeDic, PER_len=False, filename=None, filepath=None):


    print('len CO2puff', len(CO2puff))
    # print('len PER_len', len(PER_len))
    print('len velSide', len(velSide))
    print('len timeSec', len(timeSec))
    # print('CO2puff[-3000:]', CO2puff[-1000:])
    # print('timeSec[-3000:]', timeSec[-1000:])
    # print('PER_len[-3000:]', PER_len[-1000:])
  
    print('timeSec[0]', timeSec[0])
    print('timeSec[-1]', timeSec[-1])
  
    Num_ROI = len(GCset)
    print('len(GCset)', len(GCset))
    
    Totalnum_row = int(Num_ROI+5) # + 3 optic flow traces (AP, ML, Yaw) + 1 CO2 + 1 PER
    rowspan = int(Totalnum_row/Totalnum_row)

    velForwMax, velForwMin=utils.desired_min_max_yaxis(velForw[1000:], desired_min=-2, desired_max=2)
    velSideMax, velSideMin=utils.desired_min_max_yaxis(velSide[1000:], desired_min=-1, desired_max=1)
    velTurnMax, velTurnMin=utils.desired_min_max_yaxis(velTurn[1000:], desired_min=-20, desired_max=20)
    maxGC, minGC=utils.desired_min_max_yaxis(np.asarray(GCset).flatten(), desired_min=False, desired_max=10)
    

    if PER_len==False:

        # Totalnum_row = int(Num_ROI+4) # + 3 optic flow traces (AP, ML, Yaw) + 1 CO2 + 1 PER
        # rowspan = int(Totalnum_row/Totalnum_row)
        PER_len=[np.nan]*len(GCset[0])
        PER_len_Max=1
        PER_len_Min=0


    else:
        PER_len_Max, PER_len_Min=utils.desired_min_max_yaxis(PER_len, desired_min=-10, desired_max=150)

        print('PER_len_Max', PER_len_Max)
        print('PER_len_Min', PER_len_Min)



    print('velForwMin', velForwMin)
    print('velForwMax', velForwMax)
    print('velSideMin', velSideMin)
    print('velSideMax', velSideMax)
    print('velTurnMin', velTurnMin)
    print('velTurnMax', velTurnMax)
    print('minGC', minGC)

    

    print('Totalnum_row', Totalnum_row)
    print('rowspan', rowspan)


    trace_dur = np.linspace(0,timeSec[-1]-timeSec[0],len(CO2puff[:]))
    print('len trace_dur',len(trace_dur))

    fig = plt.figure(facecolor='white', figsize=(Totalnum_row+math.ceil(2*Totalnum_row),Totalnum_row), dpi=150)
    fig.subplots_adjust(left=0.15, right = 0.9, wspace = 0.3, hspace = 0.3)


    row_count=rowspan

    ##plot detected event summary##
    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axPER = plt.subplot(Totalnum_row,1,row_count)
    #plt.title(key)    
    axPER.plot(timeSec, PER_len, color='k',linewidth=1)
    axPER.plot(timeSec[-1], 0, color='w',linewidth=1, alpha=0) # for includeing NaN tail to be plotted
    
    axPER.spines['bottom'].set_visible(True)
    axPER.spines['top'].set_visible(False)
    axPER.spines['right'].set_visible(False)
    axPER.get_xaxis().set_visible(True)
    axPER.get_yaxis().set_label_coords(-0.1,0.5)
    axPER.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    #ax1.axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)            
    #axCO2.set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
    axPER.set_ylim(PER_len_Min,PER_len_Max)
    axPER.set_ylabel('PER length \n(px)',size=10,color='k')
    #plt.savefig(outDirbouts + ROI + '_'+'GCsum'+ '.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    #plt.clf
    row_count+=rowspan

    ##plot detected event summary##
    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axCO2 = plt.subplot(Totalnum_row,1,row_count)
    #plt.title(key)    
    axCO2.plot(timeSec, CO2puff[:], color=CO2puff_color,linewidth=1)
    
    axCO2.spines['bottom'].set_visible(False)
    axCO2.spines['top'].set_visible(False)
    axCO2.spines['right'].set_visible(False)
    axCO2.get_xaxis().set_visible(False)
    axCO2.get_yaxis().set_label_coords(-0.1,0.5)
    axCO2.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    #ax1.axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)            
    #axCO2.set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
    axCO2.set_ylim(0,1)
    axCO2.yaxis.set_ticks(np.array([0,1]))
    axCO2.set_ylabel(r'$\rm{CO}_\mathrm{2}$' +' puff',size=10,color=CO2puff_color)
    #plt.savefig(outDirbouts + ROI + '_'+'GCsum'+ '.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    #plt.clf
    row_count+=rowspan

    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    

    ##plot GC evt summary##
    
    for i in range(0, Num_ROI):
        axGC = 'axGC' + str(i)

    #vars()[axGC]

    # print('len trace_dur', len(trace_dur))
    # print('trace_dur[-1]',trace_dur[-1])
    # print('len GCset[i]', len(GCset[i]))
        print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
        vars()[axGC] = plt.subplot(Totalnum_row,1,row_count)
        vars()[axGC].plot(timeSec, GCset[i][:], color=GC_color,linewidth=1)
        
        vars()[axGC].spines['bottom'].set_visible(False)
        vars()[axGC].spines['top'].set_visible(False)
        vars()[axGC].spines['right'].set_visible(False)
        vars()[axGC].get_xaxis().set_visible(False)
        vars()[axGC].get_yaxis().set_label_coords(-0.1,0.5)
        #vars()[axGC].set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
        vars()[axGC].set_ylim(minGC,maxGC)
        vars()[axGC].axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
        #vars()[axGC].axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)
        vars()[axGC].set_ylabel('ROI#'+str(i)+'\n'+r'$\Delta$'+'F/F (%)',size=10,color=GC_color)
        row_count+=rowspan




    ##plot event AP flow summary##
    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axAP = plt.subplot(Totalnum_row,1,row_count)
    axAP.plot(timeSec,velForw[:], color=AP_color,linewidth=1)
    
    axAP.spines['bottom'].set_visible(False)
    axAP.spines['top'].set_visible(False)
    axAP.spines['right'].set_visible(False)

    axAP.get_xaxis().set_visible(False)
    axAP.get_yaxis().set_label_coords(-0.1,0.5)
    axAP.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    #axAP.axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)
    axAP.set_ylabel(r'$\rm{V}_\mathrm{forward}$'+'\n'+r'$\rm{(mm s}^\mathrm{-1}$'+')',size=10,color=AP_color)
    #axAP.set_ylabel(r"$\mathrm{\mathsf{v_{forward}}}$",fontsize=15,color='r',rotation=90,horizontalalignment='left')
    #axAP.set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
    axAP.set_ylim(velForwMin,velForwMax)
    tempListCurAP = [velForwMin,velForwMax]
    newYAxisAP = np.array(tempListCurAP)
    axAP.yaxis.set_ticks(newYAxisAP)
    row_count+=rowspan

    #plt.savefig(outDirbouts + ROI + '_'+'APflowsum'+ '.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    #plt.clf

    ##plot event ML flow summary## 
    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axML = plt.subplot(Totalnum_row,1,row_count)
    axML.plot(timeSec,velSide[:], color=ML_color,linewidth=1)
    
    axML.spines['bottom'].set_visible(False)
    axML.spines['top'].set_visible(False)
    axML.spines['right'].set_visible(False)
    axML.get_xaxis().set_visible(False)
    axML.get_yaxis().set_label_coords(-0.1,0.5)
    axML.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    #axML.axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)
    axML.set_ylabel(r'$\rm{V}_\mathrm{side}$'+'\n'+r'$\rm{(mm s}^\mathrm{-1}$'+')',size=10,color=ML_color)
    #axML.set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
    axML.set_ylim(velSideMin,velSideMax)
    tempListCurML = [velSideMin,velSideMax]
    newYAxisML = np.array(tempListCurML)
    axML.yaxis.set_ticks(newYAxisML)
    row_count+=rowspan

       
    #plt.savefig(outDirbouts + ROI + '_'+'MLflowsum'+ '.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    #plt.clf

    ##plot event Yaw flow summary##
    print('Totalnum_row,1,row_count', Totalnum_row,1,row_count)
    axYaw = plt.subplot(Totalnum_row,1,row_count)  
    axYaw.plot(timeSec,velTurn[:], color=Yaw_color,linewidth=1)

    # axYaw.spines['bottom'].set_visible(False)
    axYaw.spines['top'].set_visible(False)
    axYaw.spines['right'].set_visible(False)
    axYaw.get_yaxis().set_label_coords(-0.1,0.5)
    # axYaw.get_xaxis().set_visible(False)
    axYaw.axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
    #axYaw.axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)
    axYaw.set_xlabel('Time (s)',size=10,color='k')
    axYaw.set_ylabel(r'$\rm{V}_\mathrm{turn}$'+'\n'+r'$\rm{(deg. s}^\mathrm{-1}$'+')',size=10,color=Yaw_color)
    #axYaw.set_xlim(timeSec[0],timeSec[-1]+0.05*timeSec[-1])
    axYaw.set_ylim(velTurnMin,velTurnMax)
    tempListCurYaw = [velTurnMin,velTurnMax]
    newYAxisYaw = np.array(tempListCurYaw)
    axYaw.yaxis.set_ticks(newYAxisYaw)
    row_count+=rowspan


    
    ## put each ethography in each whole trace

    for key in ethoTimeDic: 

        

        # if key == 'F_groom_evt' or key=='H_groom_evt' or key=='PER_evt' or key=='Push_evt':
        if key == 'F_groom_evt' or key=='H_groom_evt' or key=='CO2puff_evt' or key=='SixLeg_Move_evt' or key=='puff_evt':

            continue



        elif key=='PER_evt':
            print('plotting ', key, ' wholetrace...', plot_setting.EthoColorDic[key])    

            for a in range(0,len(ethoTimeDic[key])):

                axPER.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)
        
        else:
            print('plotting ', key, ' wholetrace...', plot_setting.EthoColorDic[key])    

            for a in range(0,len(ethoTimeDic[key])):

                #axCO2.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)

                axAP.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)
                axML.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)
                axYaw.axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0) 

                for i in range(0, Num_ROI):
                    axGC = 'axGC' + str(i)
                    vars()[axGC].axvspan(ethoTimeDic[key][a][0], ethoTimeDic[key][a][-1], alpha=0.5, color=plot_setting.EthoColorDic[key], linewidth=0)

     
            


    plt.savefig(filepath + filename+'.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    plt.savefig(filepath + filename+'.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
    plt.clf()
    plt.close()      

    
    
    return
  


def Plot_Beh_evt_sortedbyLength(timesecEthoDic, GCsetEvt, whichBeh, filename):

    print('Plotting individual beh evt...')

    timesec_rest_evt = sorted(timesecEthoDic[whichBeh],key=len)

    print('shape GCsetEvt', np.shape(GCsetEvt))

    Num_ROI = len(GCsetEvt[0])
    
    Totalnum_row = int(Num_ROI) # + 3 optic flow traces (AP, ML, Yaw) + 1 CO2
    rowspan = int(Totalnum_row/Totalnum_row)

    print('Totalnum_row',Totalnum_row)

    row_count=rowspan

    

    trace_dur = np.linspace(0,timesec_rest_evt[-1][-1]-timesec_rest_evt[-1][0],len(timesec_rest_evt[-1]))

    ##plot GC evt summary##
    
    for i in range(0, len(GCsetEvt)):

        fig = plt.figure(facecolor='white', figsize=(int(Totalnum_row),int(Totalnum_row)), dpi=300)
        fig.subplots_adjust(left=0.2, right = 0.9, wspace = 0.1, hspace = 0.1)
        fig.suptitle('ROI#'+str(i), color=ColorCode_eachROI_list[i])
        row_count=rowspan

        GCsetEvt[i]=sorted(GCsetEvt[i],key=len)


        
        for j in range(0, len(GCsetEvt[i])):

            axGC = 'axGC' + str(i) + str(j)
            #print('len(GCsetEvt[i][j]) before', len(GCsetEvt[i][j]))


            if len(GCsetEvt[i][j])!=len(timesec_rest_evt[-1]):
                nan_tail= [np.nan]*int(len(timesec_rest_evt[-1])-len(GCsetEvt[i][j]))

                #print('shape GCsetEvt[i][j]', np.shape(GCsetEvt[i][j]))
                GCsetEvt[i][j]=np.append(GCsetEvt[i][j],nan_tail)

                #print('len(GCsetEvt[i][j]) after', len(GCsetEvt[i][j]))
                #print('len(trace_dur) after', len(trace_dur))

            #if j == 0:
                #vars()[axGC].set_title('ROI#'+str(i))

            #print('shape GCsetEvt[i][j]', np.shape(GCsetEvt[i][j]))
            
            vars()[axGC] = plt.subplot(Totalnum_row,1,row_count)

            vars()[axGC].plot(trace_dur, GCsetEvt[i][j], color=ColorCode_eachROI_list[i],linewidth=1.5)            
            vars()[axGC].spines['bottom'].set_visible(False)
            vars()[axGC].spines['top'].set_visible(False)
            vars()[axGC].spines['right'].set_visible(False)
            vars()[axGC].get_xaxis().set_visible(False)
            vars()[axGC].set_xlim(trace_dur[0],trace_dur[-1])
            vars()[axGC].set_ylim(minGCEvt,maxGCEvt)
            vars()[axGC].tick_params(axis='both', which='major', length=2, labelsize=7)
            #vars()[axGC].axhline(0, linestyle='dashed',color='gray',linewidth=0.5)
            #vars()[axGC].axvline(bsl_dur, linestyle='dashed',color='gray',linewidth=0.5)
            vars()[axGC].set_ylabel('Evt#'+str(j)+'\n'+r'$\Delta$'+'R/R (%)',size=4,color='k')
            vars()[axGC].axvspan(0, timesec_rest_evt[j][-1]-timesec_rest_evt[j][0], alpha=0.5, color='lightgray',edgecolor='none')
            
            row_count+=rowspan

        vars()[axGC].get_xaxis().set_visible(True)
        vars()[axGC].set_xlabel('time (s)',size=7, color='k')
        

        plt.savefig(outDirEvents + filename +'_ROI_' + str(i) + '.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) #bbox_inches='tight', 
        plt.clf
        plt.close(fig)
    

    return



def Plot_GCEvt_avg_err(timeBoundary, \
                     GC_mean_set, GC_down_CI_set, GC_up_CI_set,\
                     rest_mean, rest_down_CI, rest_up_CI,\
                     walk_mean, walk_down_CI, walk_up_CI,\
                     eye_groom_mean, eye_groom_down_CI, eye_groom_up_CI,\
                     antennae_groom_mean, antennae_groom_down_CI, antennae_groom_up_CI,\
                     foreleg_groom_mean, foreleg_groom_down_CI, foreleg_groom_up_CI,\
                     hindleg_groom_mean, hindleg_groom_down_CI, hindleg_groom_up_CI,\
                     Abd_groom_mean, Abd_groom_down_CI, Abd_groom_up_CI,\
                     PER_mean, PER_down_CI, PER_up_CI,\
                     CO2puff_mean, CO2puff_down_CI, CO2puff_up_CI,\
                     velForw_mean, velForw_down_CI, velForw_up_CI,\
                     velSide_mean, velSide_down_CI, velSide_up_CI,\
                     velTurn_mean, velTurn_down_CI, velTurn_up_CI,\
                     target_ROI=0, filename=None, savedir=None):


    print("Plot_event_summary...")

    axvline_color=plot_setting.axvline_color
    axvline_width=plot_setting.axvline_width
    axhline_color=plot_setting.axhline_color
    axhline_width=plot_setting.axhline_width
    data_trace_width=plot_setting.data_trace_width
    data_font_size=plot_setting.data_font_size
    data_ylabel_position=plot_setting.data_ylabel_position

    CO2puff_color=plot_setting.CO2puff_color

    GC_color=plot_setting.GC_color
    GC_accompany_color_list=['black','red','gray','darkorange','deepskyblue','blue','blueviolet','gold','magenta']

    rest_color=plot_setting.rest_color
    walk_color=plot_setting.walk_color

    E_groom_color=plot_setting.E_groom_color
    A_groom_color=plot_setting.A_groom_color
    FL_groom_color=plot_setting.FL_groom_color

    HL_groom_color=plot_setting.HL_groom_color
    Abd_groom_color=plot_setting.Abd_groom_color

    PER_color=plot_setting.PER_color

    AP_color=plot_setting.AP_color
    ML_color=plot_setting.ML_color
    Yaw_color=plot_setting.Yaw_color

    axis_width=2
    tick_length=5
    tick_width=2

    print('shape GC_mean_set', np.shape(GC_mean_set))
    

    
    eventdur = np.linspace(timeBoundary[0],timeBoundary[-1],len(GC_mean_set[0]))
    
    print('eventdur', eventdur)
    print('len eventdur', len(eventdur))

    temp_min_GC_downCI=[]
    for i, trace in enumerate(GC_down_CI_set):
        min_trace=np.nanmin(trace)
        temp_min_GC_downCI.append(min_trace)
    min_GC=min(temp_min_GC_downCI)

    temp_max_GC_upCI=[]
    for i, trace in enumerate(GC_up_CI_set):
        max_trace=np.nanmax(trace)
        temp_max_GC_upCI.append(max_trace)
    max_GC=max(temp_max_GC_upCI)

    GC_down_CI=np.nanmin(list(itertools.chain.from_iterable(GC_down_CI_set)))
    GC_up_CI=np.nanmax(list(itertools.chain.from_iterable(GC_up_CI_set)))



    print('min_GC, max_GC', min_GC, max_GC)

    print('GC_down_CI', GC_down_CI)
    print('GC_up_CI', GC_up_CI)

    if GC_down_CI>0:
        GC_min=0
    else:
        GC_min=GC_down_CI


    
    if min(velForw_down_CI)>0:
        AP_min=0
    else:
        AP_min=min(velForw_down_CI)

 
    if min(velSide_down_CI)>0:
        ML_min=0
    else:
        ML_min=min(velSide_down_CI)
    

    if min(velTurn_down_CI)>0:
        yaw_min=0
    else:
        yaw_min=min(velTurn_down_CI)



    ##plot GC event summary##






    fig = plt.figure(facecolor='white', figsize=(40,30), dpi=120)
    fig.subplots_adjust(wspace = 0.25, hspace=0.6, left=0.1, right = 0.99, bottom = 0.06, top = 0.97)
    
    axGC=plt.subplot(3,3,1)
    for i in range(0,len(GC_mean_set)):
        if i ==target_ROI:
            print('Plotting target ROI')
            line_style='-'
            axGC.plot(eventdur, GC_mean_set[i], line_style, color=GC_color, linewidth=data_trace_width+4, label='ROI_'+str(i))
            axGC.fill_between(eventdur, GC_down_CI_set[i], GC_up_CI_set[i], color=GC_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
        else:
            print('Plotting accompany ROI')
            line_style='--'
            axGC.plot(eventdur, GC_mean_set[i], line_style, color=GC_accompany_color_list[i], linewidth=data_trace_width+4, label='ROI_'+str(i))
            axGC.fill_between(eventdur, GC_down_CI_set[i], GC_up_CI_set[i], color=GC_accompany_color_list[i], linewidth=0, alpha=plot_setting.err_shade_alpha)
        
    if round(min_GC)>0:
        min_GC=-1    
    axGC.set_ylim(GC_min,GC_up_CI)
    axGC.set_ylim(-10,160 )
    axGC.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axGC.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axGC.spines['bottom'].set_visible(False)
    axGC.spines['top'].set_visible(False)
    axGC.spines['right'].set_visible(False)
    axGC.spines['left'].set_linewidth(axis_width)
    axGC.get_xaxis().set_visible(False)
    axGC.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axGC.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axGC.set_ylabel(r'$\Delta$'+'F/F (%)',size=data_font_size, color=GC_color, rotation=90, alpha=1)
    axGC.legend(loc="upper right", prop={'size':data_font_size})


    axWalkRest=plt.subplot(3,3,7)
    axWalkRest.fill_between(eventdur, rest_down_CI, rest_up_CI, color=rest_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axWalkRest.fill_between(eventdur, walk_down_CI, walk_up_CI, color=walk_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axWalkRest.plot(eventdur,rest_mean, color=rest_color,linewidth=data_trace_width+4, label='Walk')
    axWalkRest.plot(eventdur,walk_mean, color=walk_color,linewidth=data_trace_width+4, label='Rest')
    axWalkRest.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axWalkRest.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axWalkRest.spines['bottom'].set_visible(True)
    axWalkRest.spines['top'].set_visible(False)
    axWalkRest.spines['right'].set_visible(False)
    axWalkRest.spines['left'].set_linewidth(axis_width)
    axWalkRest.get_xaxis().set_visible(True)
    axWalkRest.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    # axWalkRest.set_ylim(0,1)
    axWalkRest.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axWalkRest.set_ylabel('Behavior relative freq.',size=data_font_size, color='k', rotation=90, alpha=1)
    axWalkRest.set_xlabel('Time (s)',size=data_font_size, color='k', rotation=0, alpha=1)
    axWalkRest.tick_params(axis='x', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axWalkRest.legend(loc="upper right", prop={'size':data_font_size})

    axEyeAntenFlegGroom=plt.subplot(3,3,2)
    axEyeAntenFlegGroom.fill_between(eventdur, eye_groom_down_CI, eye_groom_up_CI, color=E_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axEyeAntenFlegGroom.fill_between(eventdur, antennae_groom_down_CI, antennae_groom_up_CI, color=A_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axEyeAntenFlegGroom.fill_between(eventdur, foreleg_groom_down_CI, foreleg_groom_up_CI, color=FL_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axEyeAntenFlegGroom.plot(eventdur,eye_groom_mean, color=E_groom_color,linewidth=data_trace_width+4, label='E_groom')
    axEyeAntenFlegGroom.plot(eventdur,antennae_groom_mean, color=A_groom_color,linewidth=data_trace_width+4, label='A_groom')
    axEyeAntenFlegGroom.plot(eventdur,foreleg_groom_mean, color=FL_groom_color,linewidth=data_trace_width+4, label='FL_groom')
    axEyeAntenFlegGroom.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axEyeAntenFlegGroom.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axEyeAntenFlegGroom.spines['bottom'].set_visible(False)
    axEyeAntenFlegGroom.spines['top'].set_visible(False)
    axEyeAntenFlegGroom.spines['right'].set_visible(False)
    axEyeAntenFlegGroom.spines['left'].set_linewidth(axis_width)
    axEyeAntenFlegGroom.get_xaxis().set_visible(False)
    axEyeAntenFlegGroom.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    # axEyeAntenFlegGroom.set_ylim(0,1)
    axEyeAntenFlegGroom.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axEyeAntenFlegGroom.set_ylabel('Behavior relative freq.',size=data_font_size, color='k', rotation=90, alpha=1)
    axEyeAntenFlegGroom.legend(loc="upper right", prop={'size':data_font_size})
 

    axAbdHlegGroom=plt.subplot(3,3,5)
    axAbdHlegGroom.fill_between(eventdur, hindleg_groom_down_CI, hindleg_groom_up_CI, color=HL_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axAbdHlegGroom.fill_between(eventdur, Abd_groom_down_CI, Abd_groom_up_CI, color=Abd_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axAbdHlegGroom.plot(eventdur,hindleg_groom_mean, color=HL_groom_color,linewidth=data_trace_width+4, label='HL_groom')
    axAbdHlegGroom.plot(eventdur,Abd_groom_mean, color=Abd_groom_color,linewidth=data_trace_width+4, label='Abd_groom')
    axAbdHlegGroom.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axAbdHlegGroom.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axAbdHlegGroom.spines['bottom'].set_visible(False)
    axAbdHlegGroom.spines['top'].set_visible(False)
    axAbdHlegGroom.spines['right'].set_visible(False)
    axAbdHlegGroom.spines['left'].set_linewidth(axis_width)
    axAbdHlegGroom.get_xaxis().set_visible(False)
    axAbdHlegGroom.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    # axAbdHlegGroom.set_ylim(0,1)
    axAbdHlegGroom.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axAbdHlegGroom.set_ylabel('Behavior relative freq.',size=data_font_size, color='k', rotation=90, alpha=1)
    axAbdHlegGroom.legend(loc="upper right", prop={'size':data_font_size})


    axCO2=plt.subplot(3,3,4)
    axCO2.fill_between(eventdur, CO2puff_down_CI, CO2puff_up_CI, color=CO2puff_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axCO2.plot(eventdur,CO2puff_mean, color=CO2puff_color,linewidth=data_trace_width+4, label='CO2 puff')
    axCO2.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axCO2.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axCO2.spines['bottom'].set_visible(False)
    axCO2.spines['top'].set_visible(False)
    axCO2.spines['right'].set_visible(False)
    axCO2.spines['left'].set_linewidth(axis_width)
    axCO2.get_xaxis().set_visible(False)
    axCO2.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    # axCO2.set_ylim(0,1)
    axCO2.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axCO2.set_ylabel('Behavior relative freq.',size=data_font_size, color='k', rotation=90, alpha=1)
    axCO2.legend(loc="upper right", prop={'size':data_font_size})


    axPER=plt.subplot(3,3,8)
    axPER.fill_between(eventdur, PER_down_CI, PER_up_CI, color=PER_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axPER.plot(eventdur,PER_mean, color=PER_color,linewidth=data_trace_width+4, label='PER')
    axPER.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axPER.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axPER.spines['bottom'].set_visible(True)
    axPER.spines['top'].set_visible(False)
    axPER.spines['right'].set_visible(False)
    axPER.spines['left'].set_linewidth(axis_width)
    axPER.spines['bottom'].set_linewidth(axis_width)
    axPER.get_xaxis().set_visible(True)
    axPER.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    # axPER.set_ylim(0,1)
    axPER.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axPER.set_ylabel('Behavior relative freq.',size=data_font_size, color='k', rotation=90, alpha=1)
    axPER.set_xlabel('Time (s)',size=data_font_size, color='k', rotation=0, alpha=1)
    axPER.tick_params(axis='x', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axPER.legend(loc="upper right", prop={'size':data_font_size})



    axAP=plt.subplot(3,3,3)
    axAP.fill_between(eventdur, velForw_down_CI, velForw_up_CI, color=AP_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axAP.plot(eventdur,velForw_mean, color=AP_color,linewidth=data_trace_width+4, label='Forward')
    axAP.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axAP.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axAP.spines['bottom'].set_visible(False)
    axAP.spines['top'].set_visible(False)
    axAP.spines['right'].set_visible(False)
    axAP.spines['left'].set_linewidth(axis_width)
    axAP.get_xaxis().set_visible(False)
    axAP.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axAP.set_ylim(AP_min,max(velForw_up_CI))
    axAP.set_ylim(-3,6)
    axAP.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axAP.set_ylabel(r'$\rm{V}_\mathrm{forward}$'+'\n'+r'$\rm{(mm s}^\mathrm{-1}$'+')',size=data_font_size, color='k', rotation=90, alpha=1)
    axAP.legend(loc="upper right", prop={'size':data_font_size})


    axML=plt.subplot(3,3,6)
    axML.fill_between(eventdur, velSide_down_CI, velSide_up_CI, color=ML_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axML.plot(eventdur,velSide_mean, color=ML_color,linewidth=data_trace_width+4, label='Sideway')
    axML.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axML.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axML.spines['bottom'].set_visible(False)
    axML.spines['top'].set_visible(False)
    axML.spines['right'].set_visible(False)
    axML.spines['left'].set_linewidth(axis_width)
    axML.get_xaxis().set_visible(False)
    axML.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axML.set_ylim(ML_min, max(velSide_up_CI))
    axML.set_ylim(-5, 5)
    axML.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axML.set_ylabel(r'$\rm{V}_\mathrm{side}$'+'\n'+r'$\rm{(mm s}^\mathrm{-1}$'+')',size=data_font_size, color='k', rotation=90, alpha=1)
    axML.legend(loc="upper right", prop={'size':data_font_size})

    axYaw=plt.subplot(3,3,9)
    axYaw.fill_between(eventdur, velTurn_down_CI, velTurn_up_CI, color=Yaw_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axYaw.plot(eventdur,velTurn_mean, color=Yaw_color,linewidth=data_trace_width+4, label='Turn')
    axYaw.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axYaw.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axYaw.spines['bottom'].set_visible(True)
    axYaw.spines['top'].set_visible(False)
    axYaw.spines['right'].set_visible(False)
    axYaw.spines['left'].set_linewidth(axis_width)
    axYaw.get_xaxis().set_visible(True)
    axYaw.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axYaw.set_ylim(yaw_min,max(velTurn_up_CI))
    axYaw.set_ylim(-180,180)
    axYaw.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axYaw.set_ylabel(r'$\rm{V}_\mathrm{turn}$'+'\n'+r'$\rm{(deg. s}^\mathrm{-1}$'+')',size=data_font_size, color='k', rotation=90, alpha=1)
    axYaw.set_xlabel('Time (s)',size=data_font_size, color='k', rotation=0, alpha=1)
    axYaw.tick_params(axis='x', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axYaw.legend(loc="upper right", prop={'size':data_font_size})

    #fig.tight_layout()
    print(savedir)
    plt.savefig(savedir+filename+'.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True)
    plt.savefig(savedir+filename+'.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.savefig(savedir+filename+'.svg', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.clf()
    plt.close(fig)

    
    
    return



def Plot_CO2Evt_avg_err(timeBoundary, \
                     GC_mean_set, GC_down_CI_set, GC_up_CI_set,\
                     rest_mean, rest_down_CI, rest_up_CI,\
                     f_walk_mean, f_walk_down_CI, f_walk_up_CI,\
                     b_walk_mean, b_walk_down_CI, b_walk_up_CI,\
                     eye_groom_mean, eye_groom_down_CI, eye_groom_up_CI,\
                     antennae_groom_mean, antennae_groom_down_CI, antennae_groom_up_CI,\
                     foreleg_groom_mean, foreleg_groom_down_CI, foreleg_groom_up_CI,\
                     hindleg_groom_mean, hindleg_groom_down_CI, hindleg_groom_up_CI,\
                     Abd_groom_mean, Abd_groom_down_CI, Abd_groom_up_CI,\
                     PER_mean, PER_down_CI, PER_up_CI,\
                     CO2puff_mean, CO2puff_down_CI, CO2puff_up_CI,\
                     velForw_mean, velForw_down_CI, velForw_up_CI,\
                     velSide_mean, velSide_down_CI, velSide_up_CI,\
                     velTurn_mean, velTurn_down_CI, velTurn_up_CI,\
                     filename, savedir):


    print("Plot_event_summary...")

    axvline_color=plot_setting.axvline_color
    axvline_width=plot_setting.axvline_width
    axhline_color=plot_setting.axhline_color
    axhline_width=plot_setting.axhline_width
    data_trace_width=plot_setting.data_trace_width
    data_font_size=plot_setting.data_font_size
    data_ylabel_position=np.asarray(plot_setting.data_ylabel_position)+[0.05,0]

    CO2puff_color=plot_setting.CO2puff_color

    GC_color=plot_setting.GC_color

    rest_color=plot_setting.rest_color
    walk_color=plot_setting.walk_color

    E_groom_color=plot_setting.E_groom_color
    A_groom_color=plot_setting.A_groom_color
    FL_groom_color=plot_setting.FL_groom_color

    HL_groom_color=plot_setting.HL_groom_color
    Abd_groom_color=plot_setting.Abd_groom_color

    PER_color=plot_setting.PER_color

    AP_color=plot_setting.AP_color
    ML_color=plot_setting.ML_color
    Yaw_color=plot_setting.Yaw_color

    axis_width=2
    tick_length=5
    tick_width=2
    
    print('shape CO2puff_mean', np.shape(CO2puff_mean))
    eventdur = np.linspace(timeBoundary[0],timeBoundary[-1],len(CO2puff_mean))
    
    print('eventdur', eventdur)
    print('len eventdur', len(eventdur))
    temp_min_GC_downCI=[]
    for i, trace in enumerate(GC_down_CI_set):
        min_trace=np.nanmin(trace)
        temp_min_GC_downCI.append(min_trace)
    min_GC=min(temp_min_GC_downCI)

    temp_max_GC_upCI=[]
    for i, trace in enumerate(GC_up_CI_set):
        max_trace=np.nanmax(trace)
        temp_max_GC_upCI.append(max_trace)
    max_GC=max(temp_max_GC_upCI)

    print('min_GC, max_GC', min_GC, max_GC)

    print('CO2puff_mean', CO2puff_mean)


    ##plot GC event summary##
    fig = plt.figure(facecolor='white', figsize=(40,30), dpi=120)
    fig.subplots_adjust(wspace = 0.25, hspace=0.6, left=0.1, right = 0.99, bottom = 0.06, top = 0.97)


    axCO2=plt.subplot(3,3,1)
    axCO2.fill_between(eventdur, CO2puff_down_CI, CO2puff_up_CI, color=CO2puff_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axCO2.plot(eventdur,CO2puff_mean, color=CO2puff_color,linewidth=data_trace_width+4, label='CO2 puff')
    axCO2.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axCO2.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axCO2.spines['bottom'].set_visible(False)
    axCO2.spines['top'].set_visible(False)
    axCO2.spines['right'].set_visible(False)
    axCO2.spines['left'].set_linewidth(axis_width)
    axCO2.get_xaxis().set_visible(False)
    axCO2.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axCO2.set_xlim(timeBoundary[0],timeBoundary[-1])
    #axCO2.set_ylim(0,1)
    axCO2.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axCO2.set_ylabel('CO2 relative freq.',size=data_font_size, color='k', rotation=90, alpha=1)
    axCO2.legend(loc="upper right", prop={'size':data_font_size})


    axGC=plt.subplot(3,3,4)
    for i in range(0,len(GC_mean_set)):
        axGC.fill_between(eventdur, GC_down_CI_set[i], GC_up_CI_set[i], color=GC_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
        axGC.plot(eventdur,GC_mean_set[i], color=GC_color,linewidth=data_trace_width+4)
    if round(min_GC)>0:
        min_GC=-1
    axGC.set_ylim(min_GC,round(max_GC))
    axGC.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axGC.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axGC.spines['bottom'].set_visible(False)
    axGC.spines['top'].set_visible(False)
    axGC.spines['right'].set_visible(False)
    axGC.spines['left'].set_linewidth(axis_width)
    axGC.set_xlim(timeBoundary[0],timeBoundary[-1])
    axGC.get_xaxis().set_visible(False)
    axGC.get_yaxis().set_visible(True)
    axGC.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axGC.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    # tempListCurGC = [round(min_GC),round(max_GC)]
    # newYAxisGC = np.array(tempListCurGC)
    # axGC.yaxis.set_ticks(newYAxisGC)
    axGC.set_ylabel(r'$\Delta$'+'F/F (%)',size=data_font_size, color=GC_color, rotation=90, alpha=1)



    axWalkRest=plt.subplot(3,3,7)
    axWalkRest.fill_between(eventdur, rest_down_CI, rest_up_CI, color=rest_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axWalkRest.fill_between(eventdur, walk_down_CI, walk_up_CI, color=walk_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axWalkRest.plot(eventdur,rest_mean, color=rest_color,linewidth=data_trace_width+4, label='rest')
    axWalkRest.plot(eventdur,walk_mean, color=walk_color,linewidth=data_trace_width+4, label='walk')
    axWalkRest.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axWalkRest.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axWalkRest.spines['bottom'].set_visible(True)
    axWalkRest.spines['top'].set_visible(False)
    axWalkRest.spines['right'].set_visible(False)
    axWalkRest.spines['left'].set_linewidth(axis_width)
    axWalkRest.spines['bottom'].set_linewidth(axis_width)
    axWalkRest.set_xlim(timeBoundary[0],timeBoundary[-1])
    axWalkRest.get_xaxis().set_visible(True)
    axWalkRest.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    #axWalkRest.set_ylim(0,1)
    axWalkRest.tick_params(axis='x', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axWalkRest.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axWalkRest.set_xlabel('Time (s)',size=data_font_size, color='k', rotation=0, alpha=1)
    axWalkRest.set_ylabel('Counts',size=data_font_size, color='k', rotation=90, alpha=1)
    axWalkRest.legend(loc="upper right", prop={'size':data_font_size})

    axEyeAntenFlegGroom=plt.subplot(3,3,2)
    axEyeAntenFlegGroom.fill_between(eventdur, eye_groom_down_CI, eye_groom_up_CI, color=E_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axEyeAntenFlegGroom.fill_between(eventdur, antennae_groom_down_CI, antennae_groom_up_CI, color=A_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axEyeAntenFlegGroom.fill_between(eventdur, foreleg_groom_down_CI, foreleg_groom_up_CI, color=FL_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axEyeAntenFlegGroom.plot(eventdur,eye_groom_mean, color=E_groom_color,linewidth=data_trace_width+4, label='E_groom')
    axEyeAntenFlegGroom.plot(eventdur,antennae_groom_mean, color=A_groom_color,linewidth=data_trace_width+4, label='A_groom')
    axEyeAntenFlegGroom.plot(eventdur,foreleg_groom_mean, color=FL_groom_color,linewidth=data_trace_width+4, label='FL_groom')
    axEyeAntenFlegGroom.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axEyeAntenFlegGroom.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axEyeAntenFlegGroom.spines['bottom'].set_visible(False)
    axEyeAntenFlegGroom.spines['top'].set_visible(False)
    axEyeAntenFlegGroom.spines['right'].set_visible(False)
    axEyeAntenFlegGroom.spines['left'].set_linewidth(axis_width)
    axEyeAntenFlegGroom.get_xaxis().set_visible(False)
    axEyeAntenFlegGroom.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axEyeAntenFlegGroom.set_xlim(timeBoundary[0],timeBoundary[-1])
    #axEyeAntenFlegGroom.set_ylim(0,1)
    axEyeAntenFlegGroom.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axEyeAntenFlegGroom.set_ylabel('Counts',size=data_font_size, color='k', rotation=90, alpha=1)
    axEyeAntenFlegGroom.legend(loc="upper right", prop={'size':data_font_size})
 

    axAbdHlegGroom=plt.subplot(3,3,5)
    axAbdHlegGroom.fill_between(eventdur, hindleg_groom_down_CI, hindleg_groom_up_CI, color=HL_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axAbdHlegGroom.fill_between(eventdur, Abd_groom_down_CI, Abd_groom_up_CI, color=Abd_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axAbdHlegGroom.plot(eventdur,hindleg_groom_mean, color=HL_groom_color,linewidth=data_trace_width+4, label='HL_groom')
    axAbdHlegGroom.plot(eventdur,Abd_groom_mean, color=Abd_groom_color,linewidth=data_trace_width+4, label='Abd_groom')
    axAbdHlegGroom.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axAbdHlegGroom.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axAbdHlegGroom.spines['bottom'].set_visible(False)
    axAbdHlegGroom.spines['top'].set_visible(False)
    axAbdHlegGroom.spines['right'].set_visible(False)
    axAbdHlegGroom.spines['left'].set_linewidth(axis_width)
    axAbdHlegGroom.get_xaxis().set_visible(False)
    axAbdHlegGroom.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axAbdHlegGroom.set_xlim(timeBoundary[0],timeBoundary[-1])
    #axAbdHlegGroom.set_ylim(0,1)
    axAbdHlegGroom.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axAbdHlegGroom.set_ylabel('Counts',size=data_font_size, color='k', rotation=90, alpha=1)
    axAbdHlegGroom.legend(loc="upper right", prop={'size':data_font_size})


    axPER=plt.subplot(3,3,8)
    axPER.fill_between(eventdur, PER_down_CI, PER_up_CI, color=PER_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axPER.plot(eventdur,PER_mean, color=PER_color,linewidth=data_trace_width+4, label='PER')
    axPER.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axPER.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axPER.spines['bottom'].set_visible(True)
    axPER.spines['top'].set_visible(False)
    axPER.spines['right'].set_visible(False)
    axPER.spines['left'].set_linewidth(axis_width)
    axPER.spines['bottom'].set_linewidth(axis_width)
    axPER.get_xaxis().set_visible(True)
    axPER.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axPER.set_xlim(timeBoundary[0],timeBoundary[-1])
    #axPER.set_ylim(0,1)
    axPER.tick_params(axis='x', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axPER.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axPER.set_ylabel('Counts',size=data_font_size, color='k', rotation=90, alpha=1)
    axPER.set_xlabel('Time (s)',size=data_font_size, color='k', rotation=0, alpha=1)
    axPER.legend(loc="upper right", prop={'size':data_font_size})


    axAP=plt.subplot(3,3,3)
    axAP.fill_between(eventdur, velForw_down_CI, velForw_up_CI, color=AP_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axAP.plot(eventdur,velForw_mean, color=AP_color,linewidth=data_trace_width+4)
    axAP.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axAP.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axAP.spines['bottom'].set_visible(False)
    axAP.spines['top'].set_visible(False)
    axAP.spines['right'].set_visible(False)
    axAP.spines['left'].set_linewidth(axis_width)
    axAP.get_xaxis().set_visible(False)
    axAP.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axAP.set_xlim(timeBoundary[0],timeBoundary[-1])
    axAP.set_ylim(round(min(velForw_down_CI)),np.ceil(max(velForw_up_CI)))
    axAP.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axAP.set_ylabel(r'$\rm{V}_\mathrm{forward}$'+'\n'+r'$\rm{(mm s}^\mathrm{-1}$'+')',size=data_font_size, color='k', rotation=90, alpha=1)


    axML=plt.subplot(3,3,6)
    axML.fill_between(eventdur, velSide_down_CI, velSide_up_CI, color=ML_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axML.plot(eventdur,velSide_mean, color=ML_color,linewidth=data_trace_width+4)
    axML.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axML.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axML.spines['bottom'].set_visible(False)
    axML.spines['top'].set_visible(False)
    axML.spines['right'].set_visible(False)
    axML.spines['left'].set_linewidth(axis_width)
    axML.get_xaxis().set_visible(False)
    axML.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axML.set_xlim(timeBoundary[0],timeBoundary[-1])
    axML.set_ylim(round(min(velSide_down_CI)),np.ceil(max(velSide_up_CI)))
    axML.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axML.set_ylabel(r'$\rm{V}_\mathrm{side}$'+'\n'+r'$\rm{(mm s}^\mathrm{-1}$'+')',size=data_font_size, color='k', rotation=90, alpha=1)

    axYaw=plt.subplot(3,3,9)
    axYaw.fill_between(eventdur, velTurn_down_CI, velTurn_up_CI, color=Yaw_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axYaw.plot(eventdur,velTurn_mean, color=Yaw_color,linewidth=data_trace_width+4)
    axYaw.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axYaw.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axYaw.spines['bottom'].set_visible(True)
    axYaw.spines['top'].set_visible(False)
    axYaw.spines['right'].set_visible(False)
    axYaw.spines['left'].set_linewidth(axis_width)
    axYaw.spines['bottom'].set_linewidth(axis_width)
    axYaw.get_xaxis().set_visible(True)
    axYaw.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axYaw.set_xlim(timeBoundary[0],timeBoundary[-1])
    axYaw.set_ylim(round(min(velTurn_down_CI)),np.ceil(max(velTurn_up_CI)))
    axYaw.tick_params(axis='x', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axYaw.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axYaw.set_ylabel(r'$\rm{V}_\mathrm{turn}$'+'\n'+r'$\rm{(deg. s}^\mathrm{-1}$'+')',size=data_font_size, color='k', rotation=90, alpha=1)
    axYaw.set_xlabel('Time (s)',size=data_font_size, color='k', rotation=0, alpha=1)
    #fig.tight_layout()
    
    plt.savefig(savedir+filename+'.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.savefig(savedir+filename+'.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.clf()
    plt.close(fig)

    
    
    return





def Plot_BWevt_avg_err(timeBoundary, \
                     main_mean, main_down_CI, main_up_CI,\
                     vice_mean, vice_down_CI, vice_up_CI,\
                     GC_mean_set, GC_down_CI_set, GC_up_CI_set,\
                     rest_mean, rest_down_CI, rest_up_CI,\
                     f_walk_mean, f_walk_down_CI, f_walk_up_CI,\
                     b_walk_mean, b_walk_down_CI, b_walk_up_CI,\
                     eye_groom_mean, eye_groom_down_CI, eye_groom_up_CI,\
                     antennae_groom_mean, antennae_groom_down_CI, antennae_groom_up_CI,\
                     foreleg_groom_mean, foreleg_groom_down_CI, foreleg_groom_up_CI,\
                     hindleg_groom_mean, hindleg_groom_down_CI, hindleg_groom_up_CI,\
                     Abd_groom_mean, Abd_groom_down_CI, Abd_groom_up_CI,\
                     PER_mean, PER_down_CI, PER_up_CI,\
                     CO2puff_mean, CO2puff_down_CI, CO2puff_up_CI,\
                     BWwoCO2puff_mean, BWwoCO2puff_down_CI, BWwoCO2puff_up_CI,\
                     velForw_mean, velForw_down_CI, velForw_up_CI,\
                     velSide_mean, velSide_down_CI, velSide_up_CI,\
                     velTurn_mean, velTurn_down_CI, velTurn_up_CI,\
                     mainColor='k', mainlabel='main', viceColor='r', vicelabel='vice',\
                     filename=None, savedir=None):


    print("Plot_event_summary...")

    axvline_color=plot_setting.axvline_color
    axvline_width=plot_setting.axvline_width
    axhline_color=plot_setting.axhline_color
    axhline_width=plot_setting.axhline_width
    data_trace_width=plot_setting.data_trace_width
    data_font_size=plot_setting.data_font_size
    data_ylabel_position=np.asarray(plot_setting.data_ylabel_position)+[0.05,0]

    CO2puff_color=plot_setting.CO2puff_color

    GC_color=plot_setting.GC_color

    rest_color=plot_setting.rest_color
    f_walk_color=plot_setting.FW_color
    b_walk_color=plot_setting.BW_color

    E_groom_color=plot_setting.E_groom_color
    A_groom_color=plot_setting.A_groom_color
    FL_groom_color=plot_setting.FL_groom_color

    HL_groom_color=plot_setting.HL_groom_color
    Abd_groom_color=plot_setting.Abd_groom_color

    PER_color=plot_setting.PER_color

    AP_color=plot_setting.AP_color
    ML_color=plot_setting.ML_color
    Yaw_color=plot_setting.Yaw_color

    axis_width=2
    tick_length=5
    tick_width=2
    
    print('shape CO2puff_mean', np.shape(CO2puff_mean))
    eventdur = np.linspace(timeBoundary[0],timeBoundary[-1],len(CO2puff_mean))
    
    print('eventdur', eventdur)
    print('len eventdur', len(eventdur))
    temp_min_GC_downCI=[]
    for i, trace in enumerate(GC_down_CI_set):
        min_trace=np.nanmin(trace)
        temp_min_GC_downCI.append(min_trace)
    min_GC=min(temp_min_GC_downCI)

    temp_max_GC_upCI=[]
    for i, trace in enumerate(GC_up_CI_set):
        max_trace=np.nanmax(trace)
        temp_max_GC_upCI.append(max_trace)
    max_GC=max(temp_max_GC_upCI)

    print('min_GC, max_GC', min_GC, max_GC)

    print('CO2puff_mean', CO2puff_mean)


    ##plot GC event summary##
    fig = plt.figure(facecolor='white', figsize=(40,30), dpi=120)
    fig.subplots_adjust(wspace = 0.25, hspace=0.6, left=0.1, right = 0.99, bottom = 0.06, top = 0.97)


    axMain=plt.subplot(3,3,1)
    axMain.fill_between(eventdur, main_down_CI, main_up_CI, color=mainColor, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axMain.plot(eventdur,main_mean, color=mainColor,linewidth=data_trace_width+4, label=mainlabel)
    axMain.fill_between(eventdur, vice_down_CI, vice_up_CI, color=viceColor, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axMain.plot(eventdur,vice_mean, color=viceColor,linewidth=data_trace_width+4, label=vicelabel)
    axMain.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axMain.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axMain.spines['bottom'].set_visible(False)
    axMain.spines['top'].set_visible(False)
    axMain.spines['right'].set_visible(False)
    axMain.spines['left'].set_linewidth(axis_width)
    axMain.get_xaxis().set_visible(False)
    axMain.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axMain.set_xlim(timeBoundary[0],timeBoundary[-1])
    #axMain.set_ylim(0,1)
    axMain.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axMain.set_ylabel(mainlabel+' relative freq.',size=data_font_size, color='k', rotation=90, alpha=1)
    axMain.legend(loc="upper right", prop={'size':data_font_size})


    axGC=plt.subplot(3,3,4)
    for i in range(0,len(GC_mean_set)):
        if i==0:
            linestyle='-'
        elif i==1:
            linestyle='--'
        elif i==2:
            linestyle='-.'
        else:
            linestyle='---'
        axGC.fill_between(eventdur, GC_down_CI_set[i], GC_up_CI_set[i], color=GC_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
        axGC.plot(eventdur,GC_mean_set[i], color=GC_color,linewidth=data_trace_width+4, linestyle=linestyle)
    if round(min_GC)>0:
        min_GC=-1
    print('min_GC,round(max_GC)', min_GC,round(max_GC))
    axGC.set_ylim(min_GC,round(max_GC))
    axGC.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axGC.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axGC.spines['bottom'].set_visible(False)
    axGC.spines['top'].set_visible(False)
    axGC.spines['right'].set_visible(False)
    axGC.spines['left'].set_linewidth(axis_width)
    axGC.set_xlim(timeBoundary[0],timeBoundary[-1])
    axGC.set_ylim(-5,35)#For SS36112 suppfig
    axGC.get_xaxis().set_visible(False)
    axGC.get_yaxis().set_visible(True)
    axGC.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axGC.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    # tempListCurGC = [round(min_GC),round(max_GC)]
    # newYAxisGC = np.array(tempListCurGC)
    # axGC.yaxis.set_ticks(newYAxisGC)
    axGC.set_ylabel(r'$\Delta$'+'F/F (%)',size=data_font_size, color=GC_color, rotation=90, alpha=1)



    axWalkRest=plt.subplot(3,3,7)
    axWalkRest.fill_between(eventdur, rest_down_CI, rest_up_CI, color=rest_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axWalkRest.fill_between(eventdur, f_walk_down_CI, f_walk_up_CI, color=f_walk_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axWalkRest.fill_between(eventdur, b_walk_down_CI, b_walk_up_CI, color=b_walk_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axWalkRest.plot(eventdur,rest_mean, color=rest_color,linewidth=data_trace_width+4, label='rest')
    axWalkRest.plot(eventdur,f_walk_mean, color=f_walk_color,linewidth=data_trace_width+4, label='f_walk')
    axWalkRest.plot(eventdur,b_walk_mean, color=b_walk_color,linewidth=data_trace_width+4, label='b_walk')
    axWalkRest.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axWalkRest.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axWalkRest.spines['bottom'].set_visible(True)
    axWalkRest.spines['top'].set_visible(False)
    axWalkRest.spines['right'].set_visible(False)
    axWalkRest.spines['left'].set_linewidth(axis_width)
    axWalkRest.spines['bottom'].set_linewidth(axis_width)
    axWalkRest.set_xlim(timeBoundary[0],timeBoundary[-1])
    axWalkRest.get_xaxis().set_visible(True)
    axWalkRest.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    #axWalkRest.set_ylim(0,1)
    axWalkRest.tick_params(axis='x', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axWalkRest.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axWalkRest.set_xlabel('Time (s)',size=data_font_size, color='k', rotation=0, alpha=1)
    axWalkRest.set_ylabel('Counts',size=data_font_size, color='k', rotation=90, alpha=1)
    axWalkRest.legend(loc="upper right", prop={'size':data_font_size})

    axEyeAntenFlegGroom=plt.subplot(3,3,2)
    axEyeAntenFlegGroom.fill_between(eventdur, eye_groom_down_CI, eye_groom_up_CI, color=E_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axEyeAntenFlegGroom.fill_between(eventdur, antennae_groom_down_CI, antennae_groom_up_CI, color=A_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axEyeAntenFlegGroom.fill_between(eventdur, foreleg_groom_down_CI, foreleg_groom_up_CI, color=FL_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axEyeAntenFlegGroom.plot(eventdur,eye_groom_mean, color=E_groom_color,linewidth=data_trace_width+4, label='E_groom')
    axEyeAntenFlegGroom.plot(eventdur,antennae_groom_mean, color=A_groom_color,linewidth=data_trace_width+4, label='A_groom')
    axEyeAntenFlegGroom.plot(eventdur,foreleg_groom_mean, color=FL_groom_color,linewidth=data_trace_width+4, label='FL_groom')
    axEyeAntenFlegGroom.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axEyeAntenFlegGroom.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axEyeAntenFlegGroom.spines['bottom'].set_visible(False)
    axEyeAntenFlegGroom.spines['top'].set_visible(False)
    axEyeAntenFlegGroom.spines['right'].set_visible(False)
    axEyeAntenFlegGroom.spines['left'].set_linewidth(axis_width)
    axEyeAntenFlegGroom.get_xaxis().set_visible(False)
    axEyeAntenFlegGroom.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axEyeAntenFlegGroom.set_xlim(timeBoundary[0],timeBoundary[-1])
    #axEyeAntenFlegGroom.set_ylim(0,1)
    axEyeAntenFlegGroom.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axEyeAntenFlegGroom.set_ylabel('Counts',size=data_font_size, color='k', rotation=90, alpha=1)
    axEyeAntenFlegGroom.legend(loc="upper right", prop={'size':data_font_size})
 

    axAbdHlegGroom=plt.subplot(3,3,5)
    axAbdHlegGroom.fill_between(eventdur, hindleg_groom_down_CI, hindleg_groom_up_CI, color=HL_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axAbdHlegGroom.fill_between(eventdur, Abd_groom_down_CI, Abd_groom_up_CI, color=Abd_groom_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axAbdHlegGroom.plot(eventdur,hindleg_groom_mean, color=HL_groom_color,linewidth=data_trace_width+4, label='HL_groom')
    axAbdHlegGroom.plot(eventdur,Abd_groom_mean, color=Abd_groom_color,linewidth=data_trace_width+4, label='Abd_groom')
    axAbdHlegGroom.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axAbdHlegGroom.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axAbdHlegGroom.spines['bottom'].set_visible(False)
    axAbdHlegGroom.spines['top'].set_visible(False)
    axAbdHlegGroom.spines['right'].set_visible(False)
    axAbdHlegGroom.spines['left'].set_linewidth(axis_width)
    axAbdHlegGroom.get_xaxis().set_visible(False)
    axAbdHlegGroom.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axAbdHlegGroom.set_xlim(timeBoundary[0],timeBoundary[-1])
    #axAbdHlegGroom.set_ylim(0,1)
    axAbdHlegGroom.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axAbdHlegGroom.set_ylabel('Counts',size=data_font_size, color='k', rotation=90, alpha=1)
    axAbdHlegGroom.legend(loc="upper right", prop={'size':data_font_size})


    axPERCO2=plt.subplot(3,3,8)
    axPERCO2.fill_between(eventdur, PER_down_CI, PER_up_CI, color=PER_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axPERCO2.plot(eventdur,PER_mean, color=PER_color,linewidth=data_trace_width+4, label='PER')
    axPERCO2.fill_between(eventdur, CO2puff_down_CI, CO2puff_up_CI, color=CO2puff_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axPERCO2.plot(eventdur,CO2puff_mean, color=CO2puff_color,linewidth=data_trace_width+4, label='CO2 puff')
    axPERCO2.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axPERCO2.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axPERCO2.spines['bottom'].set_visible(True)
    axPERCO2.spines['top'].set_visible(False)
    axPERCO2.spines['right'].set_visible(False)
    axPERCO2.spines['left'].set_linewidth(axis_width)
    axPERCO2.spines['bottom'].set_linewidth(axis_width)
    axPERCO2.get_xaxis().set_visible(True)
    axPERCO2.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axPERCO2.set_xlim(timeBoundary[0],timeBoundary[-1])
    #axPERCO2.set_ylim(0,1)
    axPERCO2.tick_params(axis='x', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axPERCO2.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axPERCO2.set_ylabel('Counts',size=data_font_size, color='k', rotation=90, alpha=1)
    axPERCO2.set_xlabel('Time (s)',size=data_font_size, color='k', rotation=0, alpha=1)
    axPERCO2.legend(loc="upper right", prop={'size':data_font_size})


    axAP=plt.subplot(3,3,3)
    axAP.fill_between(eventdur, velForw_down_CI, velForw_up_CI, color=AP_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axAP.plot(eventdur,velForw_mean, color=AP_color,linewidth=data_trace_width+4)
    axAP.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axAP.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axAP.spines['bottom'].set_visible(False)
    axAP.spines['top'].set_visible(False)
    axAP.spines['right'].set_visible(False)
    axAP.spines['left'].set_linewidth(axis_width)
    axAP.get_xaxis().set_visible(False)
    axAP.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axAP.set_xlim(timeBoundary[0],timeBoundary[-1])
    # axAP.set_ylim(round(min(velForw_down_CI)),np.ceil(max(velForw_up_CI)))
    axAP.set_ylim(-4,4)#For SS36112 suppfig
    axAP.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axAP.set_ylabel(r'$\rm{V}_\mathrm{forward}$'+'\n'+r'$\rm{(mm s}^\mathrm{-1}$'+')',size=data_font_size, color='k', rotation=90, alpha=1)


    axML=plt.subplot(3,3,6)
    axML.fill_between(eventdur, velSide_down_CI, velSide_up_CI, color=ML_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axML.plot(eventdur,velSide_mean, color=ML_color,linewidth=data_trace_width+4)
    axML.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axML.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axML.spines['bottom'].set_visible(False)
    axML.spines['top'].set_visible(False)
    axML.spines['right'].set_visible(False)
    axML.spines['left'].set_linewidth(axis_width)
    axML.get_xaxis().set_visible(False)
    axML.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axML.set_xlim(timeBoundary[0],timeBoundary[-1])
    # axML.set_ylim(round(min(velSide_down_CI)),np.ceil(max(velSide_up_CI)))
    axML.set_ylim(-2,2)#For SS36112 suppfig
    axML.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axML.set_ylabel(r'$\rm{V}_\mathrm{side}$'+'\n'+r'$\rm{(mm s}^\mathrm{-1}$'+')',size=data_font_size, color='k', rotation=90, alpha=1)

    axYaw=plt.subplot(3,3,9)
    axYaw.fill_between(eventdur, velTurn_down_CI, velTurn_up_CI, color=Yaw_color, linewidth=0, alpha=plot_setting.err_shade_alpha)
    axYaw.plot(eventdur,velTurn_mean, color=Yaw_color,linewidth=data_trace_width+4)
    axYaw.axhline(0, linestyle='dashed',color=axhline_color,linewidth=axhline_width)
    axYaw.axvline(0, linestyle='dashed',color=axvline_color,linewidth=axvline_width)
    axYaw.spines['bottom'].set_visible(True)
    axYaw.spines['top'].set_visible(False)
    axYaw.spines['right'].set_visible(False)
    axYaw.spines['left'].set_linewidth(axis_width)
    axYaw.spines['bottom'].set_linewidth(axis_width)
    axYaw.get_xaxis().set_visible(True)
    axYaw.get_yaxis().set_label_coords(data_ylabel_position[0],data_ylabel_position[1])
    axYaw.set_xlim(timeBoundary[0],timeBoundary[-1])
    # axYaw.set_ylim(round(min(velTurn_down_CI)),np.ceil(max(velTurn_up_CI)))
    axYaw.set_ylim(-60,90) #for SS36112 suppfig
    axYaw.tick_params(axis='x', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    axYaw.tick_params(axis='y', colors='k',right='off', labelsize=data_font_size, length=tick_length, width=tick_width) 
    #tempListCurGC = [round(min(GC_down_CI)-0.1*min(GC_down_CI)),round(max(GC_up_CI)+0.1*max(GC_up_CI))]
    axYaw.set_ylabel(r'$\rm{V}_\mathrm{turn}$'+'\n'+r'$\rm{(deg. s}^\mathrm{-1}$'+')',size=data_font_size, color='k', rotation=90, alpha=1)
    axYaw.set_xlabel('Time (s)',size=data_font_size, color='k', rotation=0, alpha=1)
    #fig.tight_layout()
    
    plt.savefig(savedir+filename+'.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.savefig(savedir+filename+'.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.clf()
    plt.close(fig)

    
    
    return



def plot_scatterdots_box_allGroups(data2dlist_allgroups, y_lim, x_labels, p_value_list, savedir, filename):



    plt.boxplot(data2dlist_allgroups, labels=x_labels)
    ngroup = len(data2dlist_allgroups)
    clevels = np.linspace(0., 1., ngroup)

    xs=[]
    for i, data2dlist_eachgroups in enumerate(data2dlist_allgroups):
        xs.append(np.random.normal(i+1, 0.11, len(data2dlist_eachgroups)))

    for x, data2dlist_eachgroups in zip(xs, data2dlist_allgroups):
        print('x', x)
        print('data2dlist_eachgroups', data2dlist_eachgroups)
        plt.scatter(x, data2dlist_eachgroups, color='k', alpha=0.4)

    plt.ylim([y_lim[0], y_lim[-1]])


    
    plt.savefig(savedir+filename, edgecolor='none', transparent=True) 
    plt.clf()
    plt.close()


    return




def plot_group_bar_w_scatterPoints(x_labels, GC_ONBALL_set, GC_OFFBALL_set, y_lim, p_value_list, labels, savedir, filename):



    maxGC=y_lim[1]
    minGC=y_lim[0]
    y_range=maxGC-minGC   


    print('maxGC', maxGC)
    print('minGC', minGC)

    star_lab_list=[]
    for p_val in p_value_list:
        print('p_val', p_val)
        if p_val>0.05:
            star_lab='ns'
        elif p_val<0.05 and p_val >0.01:
            star_lab='*'
        elif p_val<0.01 and p_val >0.001:
            star_lab='**'
        elif p_val<0.001:
            star_lab='***'

        star_lab_list.append(star_lab)


    GCmean_OFFBALL_per_ROI=np.nanmean(GC_OFFBALL_set, axis=1)
    GCmean_ONBALL_per_ROI=np.nanmean(GC_ONBALL_set, axis=1)

    print('GCmean_OFFBALL_per_ROI', GCmean_OFFBALL_per_ROI)
    print('GCmean_ONBALL_per_ROI', GCmean_ONBALL_per_ROI)


    x = np.arange(len(x_labels))+1  # the label locations
    width = 0.35  # the width of the bars


    position1=list(np.arange(1, 3*(len(GC_ONBALL_set)),3))
    position2=list(np.arange(1, 3*(len(GC_OFFBALL_set)),3)+1)
    position_p_value=list(np.arange(1, 3*(len(GC_OFFBALL_set)),3)+0.5)

    print('position1', position1)
    print('position2', position2)
    print('position_p_value', position_p_value)

    print('labels', labels)

    print('len GC_ONBALL_set', len(GC_ONBALL_set))
    print('len GC_OFFBALL_set', len(GC_OFFBALL_set))


    fig, ax = plt.subplots()

    # rects1 = ax.bar(x - width/2, GCmean_ONBALL_per_ROI, width, color=(0,0,0,0), edgecolor='k', label=labels[0])
    bp1=ax.boxplot(GC_ONBALL_set, positions=position1, notch=False, showfliers=True, whis=1.5, widths = 0.8)
    for i in range(len(position1)):
        ax.scatter(position1[i] - width + np.random.random(len(GC_ONBALL_set[i])) * width, GC_ONBALL_set[i], color='k', linewidths=0, alpha=0.3)
        for box in bp1['boxes']:
            # change outline color
            box.set( color='k',  linewidth=1,)
            # # change fill color
            # box.set( facecolor = 'w' )
        ## change color and linewidth of the whiskers
        for whisker in bp1['whiskers']:
            whisker.set(color='k', linewidth=1)

        ## change color and linewidth of the caps
        for cap in bp1['caps']:
            cap.set(color='k', linewidth=0)

        ## change color and linewidth of the medians
        for median in bp1['medians']:
            median.set(color='k', linewidth=1)

        ## change the style of fliers and their fill
        for flier in bp1['fliers']:
            flier.set(marker='o', color='k', alpha=1)

    # rects2 = ax.bar(x + width/2, GCmean_OFFBALL_per_ROI, width, color=(0,0,0,0), edgecolor='r', label=labels[1])
    bp2=ax.boxplot(GC_OFFBALL_set, positions=position2, notch=False, showfliers=True, whis=1.5, widths = 0.8)
    for i in range(len(position2)):
        ax.scatter(position2[i] - width + np.random.random(len(GC_OFFBALL_set[i])) * width, GC_OFFBALL_set[i], color='r', linewidths=0, alpha=0.3)
        for box in bp2['boxes']:
            # change outline color
            box.set( color='r')
            # # change fill color
            # box.set( facecolor = 'w' )
            ## change color and linewidth of the whiskers
        for whisker in bp2['whiskers']:
            whisker.set(color='r', linewidth=1)

        ## change color and linewidth of the caps
        for cap in bp2['caps']:
            cap.set(color='r', linewidth=0)

        ## change color and linewidth of the medians
        for median in bp2['medians']:
            median.set(color='r', linewidth=1)

        ## change the style of fliers and their fill
        for flier in bp2['fliers']:
            flier.set(marker='o', color='r', alpha=1)

    for i, star in enumerate(star_lab_list):
        plt.text(position_p_value[i], maxGC+0.13*y_range, star)
        plt.text(position_p_value[i], maxGC+0.1*y_range, "{:.4f}".format(p_value_list[i]))

    


    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_xlim(0, position2[-1]+1)
    ax.set_ylim(minGC-0.1*y_range, maxGC+0.2*y_range)
    ax.set_ylabel(r'$\Delta$'+'F/F(%)')
    ax.set_title(filename)
    ax.set_xticks(position_p_value)
    ax.set_xticklabels(x_labels)
    ax.legend([bp1["boxes"][0], bp2["boxes"][0]], labels, loc='upper right')


    # ax.bar_label(rects1, padding=3)
    # ax.bar_label(rects2, padding=3)

    fig.tight_layout()

    plt.savefig(savedir+filename+'.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.savefig(savedir+filename+'.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.clf()
    plt.close(fig)


    return




def plot_box_box_w_scatterPoints(labels, GC_OFFBALLmoving_set, GC_ONBALLmoving_set, y_lim, p_value_list, savedir, filename):


    # gapfree_GC = utils.flatten_list(GC_OFFBALLmoving_set)+utils.flatten_list(GC_ONBALLmoving_set)
    # maxGC=np.nanmax(gapfree_GC)
    # minGC=np.nanmin(gapfree_GC)
    # y_range=maxGC-minGC

    maxGC=y_lim[1]
    minGC=y_lim[0]
    y_range=maxGC-minGC   

    print('maxGC', maxGC)
    print('minGC', minGC)

    star_lab_list=[]
    for p_val in p_value_list:
        print('p_val', p_val)
        if p_val>0.05:
            star_lab='ns'
        elif p_val<0.05 and p_val >0.01:
            star_lab='*'
        elif p_val<0.01 and p_val >0.001:
            star_lab='**'
        elif p_val<0.001:
            star_lab='***'

        star_lab_list.append(star_lab)


    GCmean_OFFBALLmoving_per_ROI=np.nanmean(GC_OFFBALLmoving_set, axis=1)
    GCmean_ONBALLmoving_per_ROI=np.nanmean(GC_ONBALLmoving_set, axis=1)

    print('GCmean_OFFBALLmoving_per_ROI', GCmean_OFFBALLmoving_per_ROI)
    print('GCmean_ONBALLmoving_per_ROI', GCmean_ONBALLmoving_per_ROI)


    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()

    rects1 = ax.bar(x - width/2, GCmean_ONBALLmoving_per_ROI, width, color=(0,0,0,0), edgecolor='k', label='On-ball')
    for i in range(len(x)):
        ax.scatter(x[i] - width/2 + np.random.random(len(GC_ONBALLmoving_set[i])) * width/2-width/3, GC_ONBALLmoving_set[i], color='k')

    rects2 = ax.bar(x + width/2, GCmean_OFFBALLmoving_per_ROI, width, color=(0,0,0,0), edgecolor='r', label='Off-ball')
    for i in range(len(x)):
        ax.scatter(x[i] + width/2 + np.random.random(len(GC_OFFBALLmoving_set[i])) * width/2-width/3, GC_OFFBALLmoving_set[i], color='r')

    for i, star in enumerate(star_lab_list):
        plt.text(i, maxGC+0.1*y_range, star)



    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylim(minGC-0.05*y_range, maxGC+0.15*maxGC)
    ax.set_ylabel(r'$\Delta$'+'F/F(%)')
    ax.set_title(filename)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    # ax.bar_label(rects1, padding=3)
    # ax.bar_label(rects2, padding=3)

    fig.tight_layout()

    plt.savefig(savedir+filename+'.png', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.savefig(savedir+filename+'.pdf', facecolor=fig.get_facecolor(), edgecolor='none', transparent=True) 
    plt.clf()
    plt.close(fig)


    return




def sorting_roiID_correspondingMat_based_on_an_order(roi_id_list, corrspd_list, base_of_order_list):



    # if len(roi_id_list)!=len(base_of_order_list):
    #   print('ROI numbers does not match. Reordering ROI IDs cannot be done ...')
    #   sys.exit(0)

    ressembled_name_for_baseOrder_list=[]
    for i, id_name in enumerate(base_of_order_list):
        ressembled_name=id_name.split(' ')[0]+'-ROI#'+id_name.split(' ')[1]
        ressembled_name_for_baseOrder_list.append(ressembled_name)


    sorted_new_roi_id_list=[]
    sorted_new_corrspd_list=[]


    numbers_for_old_order=[]
    for i, id_name in enumerate(roi_id_list):
        if id_name in roi_id_list:
            numbers_for_old_order.append(ressembled_name_for_baseOrder_list.index(id_name))

    # print('numbers_for_old_order', numbers_for_old_order)


    zipped_id_lists = zip(numbers_for_old_order, roi_id_list)
    sorted_zipped_id_lists = sorted(zipped_id_lists)
    sorted_new_roi_id_list = [element for _, element in sorted_zipped_id_lists]
    # print('sorted_new_roi_id_list', sorted_new_roi_id_list)


    zipped_corrspd_lists = zip(numbers_for_old_order, corrspd_list)
    sorted_zipped_corrspd_lists = sorted(zipped_corrspd_lists)
    sorted_new_corrspd_list = [element for _, element in sorted_zipped_corrspd_lists]


    return sorted_new_roi_id_list, sorted_new_corrspd_list


