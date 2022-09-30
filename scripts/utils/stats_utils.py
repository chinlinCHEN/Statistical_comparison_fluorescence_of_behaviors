import sys
import pickle
#import sys
import os
import re
import numpy as np
import pandas as pd
import random
from PIL import Image
Image.Image.tostring = Image.Image.tobytes
import matplotlib
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from skimage import io
from multiprocessing import Pool
import math
import more_itertools as mit
import time
from matplotlib import colors as mcolors
import scipy
import scikit_posthocs
from statsmodels.stats import weightstats as stests


def mask_out_dFF_based_on_p_value(dFF_list_ori, ref_mat):



    dFF_list = [x[:] for x in dFF_list_ori]
    mask_list = [y[:] for y in ref_mat]

    # min_dff=np.nanmin(dFF_list_ori)


    # print('type dFF_list_ori', type(dFF_list_ori))
    # print('type dFF_list_ori[0]', type(dFF_list_ori[0]))
    # print('type dFF_list_ori[0][0]', type(dFF_list_ori[0][0]))

    if len(ref_mat)==0:
        print('base_mat is not given. The refrent matrix for masking is obligatory...')
        sys.exit(0)

    for roi_i, p_values_per_roi in enumerate(ref_mat):
        for beh_i, p_value in enumerate(p_values_per_roi):
            if p_value>0.05 or np.isnan(p_value):
                dFF_list[roi_i][beh_i]=0


    for roi_i, p_values_per_roi in enumerate(ref_mat):
        for beh_i, p_value in enumerate(p_values_per_roi):
            if p_value>0.05 or np.isnan(p_value):
                mask_list[roi_i][beh_i]=1
            else:
                mask_list[roi_i][beh_i]=np.nan


    return dFF_list, mask_list



def mask_out_dFF_based_on_normality(dFF_list_ori, normality_array):


    dFF_list = [x[:] for x in dFF_list_ori]
    mask_list = [y[:] for y in ref_mat]





    if len(ref_mat)==0:
        print('base_mat is not given. The refrent matrix for masking is obligatory...')
        sys.exit(0)

    for roi_i, p_values_per_roi in enumerate(normality_array):
        for beh_i, p_value in enumerate(p_values_per_roi):
            if p_value>0.05 or np.isnan(p_value):
                dFF_list[roi_i][beh_i]=0


    for roi_i, p_values_per_roi in enumerate(ref_mat):
        for beh_i, p_value in enumerate(p_values_per_roi):
            if p_value>0.05 or np.isnan(p_value):
                mask_list[roi_i][beh_i]=1
            else:
                mask_list[roi_i][beh_i]=np.nan



    return


def make_mask_based_on_p_value_and_normality(p_value_mat, normality_arr):

    mask_list = [y[:] for y in p_value_mat]


    if len(p_value_mat)==0:
        print('base_mat is not given. The refrent matrix for masking is obligatory...')
        sys.exit(0)



    # for roi_i, normality_per_roi in enumerate(normality_arr):
    #   if normality_per_roi==True:
    #       for beh_i, mask_val in enumerate(p_value_mat[roi_i]):
    #           mask_list[roi_i][beh_i]=1
    #   else:
    #       for beh_i, mask_val in enumerate(p_value_mat[roi_i]):
    #           mask_list[roi_i][beh_i]=np.nan


    for roi_i, pvalues_per_roi in enumerate(p_value_mat):
        for beh_i, p_value in enumerate(pvalues_per_roi):
            if p_value>0.05 or np.isnan(p_value) or normality_arr[roi_i]==True:
                mask_list[roi_i][beh_i]=1
            else:
                mask_list[roi_i][beh_i]=np.nan




    return mask_list





def prep_eachBehEvts_for_stat(name_list, eachBeh_evts_set, bsl_evts_list, bsl_data_1dlist, min_evt_num=3, bsl_s=1, cutting_head_s=0.7, data_freq=1500):

    print('preparing the data for statistics ...')

    # print('len bsl_data_1dlist', len(bsl_data_1dlist))
    # print('len bsl_evts_list', len(bsl_evts_list))

    evt_count_list=[]
    for i, beh_GCevt in enumerate(eachBeh_evts_set):
        evt_count=len(beh_GCevt)
        evt_count_list.append(evt_count)
    print('evt_count_list', evt_count_list)
    mean_evt_count=int(np.nanmean(evt_count_list))

    # resample_bsl_datapoints=np.random.choice(bsl_data_1dlist, size=len(1*bsl_evts_list))
    resample_bsl_datapoints=np.random.choice(bsl_data_1dlist, size=mean_evt_count) #100, 110, 125, 150, 200 are tried so far
    # print('len resample_bsl_datapoints', len(resample_bsl_datapoints))
    bsl_mean_list=math_utils.compute_mean_with_diffrerent_row_length(bsl_evts_list, samplerate=data_freq, cutting_head_s=cutting_head_s)
    bsl_mean=np.nanmean(bsl_mean_list)
    # print('bsl_mean', bsl_mean)

    all_beh_meanlist_list=[]
    all_beh_mean=[]
    all_beh_evtCount=[]
    non_overlapBsl_perc_list=[]



    


    for i, beh_GCevt in enumerate(eachBeh_evts_set):

        print('-',name_list[i], 'has', len(beh_GCevt), 'events')

        all_beh_evtCount.append(len(beh_GCevt))

        if len(beh_GCevt)>min_evt_num:
            GCevt_fly_wo_bsl = trim_smthRow_of_list(beh_GCevt, startIdx=int(data_freq*(bsl_s)))
            print('len GCevt_fly_wo_bsl', len(GCevt_fly_wo_bsl))
            # print(name_list[i], 'shape GCevt_fly_wo_bsl', np.shape(GCevt_fly_wo_bsl))
            evtmean_list=math_utils.compute_mean_with_diffrerent_row_length(GCevt_fly_wo_bsl, samplerate=data_freq, cutting_head_s=cutting_head_s)
            print('len evtmean_list', len(evtmean_list))
            # print('evtmean_list', evtmean_list)

            # resample_beh_datapoints=np.random.choice(general_utils.flatten_list(GCevt_fly_wo_bsl), size=len(1*beh_GCevt))
            # resample_beh_datapoints=np.random.choice(general_utils.flatten_list(GCevt_fly_wo_bsl), size=100) 
            # print('resample_beh_datapoints', resample_beh_datapoints)

            ## ------------------------------------------------------------------
            ## choose the mean of each epoch or resampled datapoint from every epochs in the same length of epoch numbers
            ## 1.
            beh_data_for_stats=evtmean_list 
            ## 2.
            # beh_data_for_stats=resample_beh_datapoints
            ## ------------------------------------------------------------------

    
            # print(name_list[i], 'shape beh_data_for_stats', np.shape(beh_data_for_stats))

            if len(beh_data_for_stats)>min_evt_num:
                mean_evt=np.nanmean(beh_data_for_stats)-bsl_mean
                # print(name_list[i], 'mean_evt', mean_evt)
                f_GCevt_fly_wo_bsl=general_utils.flatten_list(GCevt_fly_wo_bsl)
                overlap_bsl_perc=overlap_btwn_two_list(f_GCevt_fly_wo_bsl, bsl_data_1dlist)
                print(name_list[i], 'overlap_bsl_perc', overlap_bsl_perc)
                if overlap_bsl_perc>0.9:
                    mean_evt=0
                    beh_data_for_stats=bsl_mean_list            
            else:
                mean_evt=0
                beh_data_for_stats=[np.nan] 
                overlap_bsl_perc=1.1
        else:
            mean_evt=0
            beh_data_for_stats=[np.nan] 
            overlap_bsl_perc=1.1

        all_beh_meanlist_list.append(beh_data_for_stats)
        all_beh_mean.append(mean_evt)
        non_overlapBsl_perc_list.append(1-overlap_bsl_perc)

        # print('all_beh_mean', all_beh_mean)


    return all_beh_meanlist_list, all_beh_mean, non_overlapBsl_perc_list, resample_bsl_datapoints




def KW_w_dataset(all_beh_meanlist_list, beh_name_list):

    ## Krukal-Wallis process
    KW_value, p_value_KW=scipy.stats.kruskal(\
        *all_beh_meanlist_list,
        nan_policy='omit'
        )
    print('p_value_KW\n', p_value_KW)


    if p_value_KW<0.05:

        # posthoc_results_conover=scikit_posthocs.posthoc_conover(np.array(all_beh_meanlist_list), p_adjust='holm')
        posthoc_results_dunn=scikit_posthocs.posthoc_dunn(np.array(all_beh_meanlist_list), p_adjust='holm')
        # posthoc_results_Mann=scikit_posthocs.posthoc_mannwhitney(np.array(all_beh_meanlist_list), p_adjust='holm')

        posthoc_results_kw=posthoc_results_dunn
        posthoc_results_kw=posthoc_results_kw.set_axis(beh_name_list, axis='columns')
        posthoc_results_kw=posthoc_results_kw.set_axis(beh_name_list, axis='index')
        # print('posthoc_results_conover\n', posthoc_results_conover)
        # print('posthoc_results_dunn\n', posthoc_results_dunn)
        # print('posthoc_results_Mann\n', posthoc_results_Mann)
        # print('type posthoc_results_conover', type(posthoc_results))

    else:
        row = [1.0]*len(beh_name_list)
        # print('row\n', row)

        null_mat=[row]*len(beh_name_list)
        # print('null_mat\n', null_mat)

        posthoc_results_kw=pd.DataFrame(null_mat)
        posthoc_results_kw=posthoc_results_kw.set_axis(beh_name_list, axis='columns')
        posthoc_results_kw=posthoc_results_kw.set_axis(beh_name_list, axis='index') 


    return posthoc_results_kw, p_value_KW



def anova_w_dataset(all_beh_meanlist_list, beh_name_list):

    print('ANOVA ...')

    # print('all_beh_meanlist_list', all_beh_meanlist_list)

    new_all_beh_meanlist_list_for_anova=[]
    new_beh_name_list_for_anova=[]
    missing_beh_meanlist_list=[]
    missing_beh_name_list=[]
    for i, meanlist in enumerate(all_beh_meanlist_list):
        print('len(meanlist)', len(meanlist))
        if len(meanlist)>1:
            new_all_beh_meanlist_list_for_anova.append(meanlist)
            new_beh_name_list_for_anova.append(beh_name_list[i])
        else:
            missing_beh_meanlist_list.append(meanlist)
            missing_beh_name_list.append(beh_name_list[i])

            bsl_mean_list=all_beh_meanlist_list[beh_name_list.index('Bsl')]
            new_all_beh_meanlist_list_for_anova.append(bsl_mean_list)
            new_beh_name_list_for_anova.append(beh_name_list[i])






    count_nanlist=sum(map(lambda x : x==[np.nan], missing_beh_meanlist_list))

    print('')

    if count_nanlist==len(all_beh_meanlist_list)-1:

        p_value_ANOVA=1
        print('p_value_ANOVA', p_value_ANOVA)

        row = [1.0]*len(beh_name_list)
        # print('row\n', row)

        null_mat=[row]*len(beh_name_list)
        # print('null_mat\n', null_mat)

        posthoc_results_anova=pd.DataFrame(null_mat)
        posthoc_results_anova=posthoc_results_anova.set_axis(beh_name_list, axis='columns')
        posthoc_results_anova=posthoc_results_anova.set_axis(beh_name_list, axis='index')



    else:

        ANOVA_value, p_value_ANOVA=scipy.stats.f_oneway(\
            *new_all_beh_meanlist_list_for_anova
            )
        print('p_value_ANOVA', p_value_ANOVA)
        #print('new_all_beh_meanlist_list_for_anova', new_all_beh_meanlist_list_for_anova)

        if p_value_ANOVA<0.05:

            print('len new_all_beh_meanlist_list_for_anova', len(new_all_beh_meanlist_list_for_anova))
            print('len new_beh_name_list_for_anova', len(new_beh_name_list_for_anova))
            print('new_beh_name_list_for_anova', new_beh_name_list_for_anova)
 

            posthoc_results_tukey=scikit_posthocs.posthoc_tukey(np.array(new_all_beh_meanlist_list_for_anova))
            # posthoc_results_scheffe=scikit_posthocs.posthoc_scheffe(np.array(new_all_beh_meanlist_list_for_anova))
            # posthoc_results_ttest=scikit_posthocs.posthoc_ttest(np.array(new_all_beh_meanlist_list_for_anova))
            # posthoc_results_tamhane=scikit_posthocs.posthoc_tamhane(np.array(new_all_beh_meanlist_list_for_anova))

            posthoc_results_anova=posthoc_results_tukey
            posthoc_results_anova=posthoc_results_anova.set_axis(new_beh_name_list_for_anova, axis='columns')
            posthoc_results_anova=posthoc_results_anova.set_axis(new_beh_name_list_for_anova, axis='index')

            posthoc_results_anova=add_missing_cols_rows_to_posthocDf(missing_beh_name_list, np.nan, beh_name_list, posthoc_results_anova)

            print('posthoc_results_tukey\n', posthoc_results_tukey)
            # print('posthoc_results_scheffe\n', posthoc_results_scheffe)
            # print('posthoc_results_ttest\n', posthoc_results_ttest)
            # print('posthoc_results_tamhane\n', posthoc_results_tamhane)
            # print('type posthoc_results_tamhane', type(posthoc_results_tamhane))
            print('posthoc_results_anova\n', posthoc_results_anova)

        else:
            row = [1.0]*len(beh_name_list)
            # print('row\n', row)

            null_mat=[row]*len(beh_name_list)
            # print('null_mat\n', null_mat)

            posthoc_results_anova=pd.DataFrame(null_mat)
            posthoc_results_anova=posthoc_results_anova.set_axis(beh_name_list, axis='columns')
            posthoc_results_anova=posthoc_results_anova.set_axis(beh_name_list, axis='index')

    print('posthoc_results_anova\n', posthoc_results_anova)



    return posthoc_results_anova, p_value_ANOVA






