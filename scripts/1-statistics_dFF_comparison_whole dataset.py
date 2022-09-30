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


import utils.general_utils as general_utils
import utils.plot_utils as plot_utils
import utils.plot_setting as plot_setting
import utils.math_utils as math_utils
import utils.sync_utils as sync_utils
import utils.list_inputFiles as list_inputFiles
import utils.stats_utils as stats_utils



experiments=list_inputFiles.all_representative_experiments










def add_missing_cols_rows_to_posthocDf(missing_name_list, assigned_val_for_missing_name, original_name_list, posthoc_results_df):

	print('missing_name_list', missing_name_list)
	# print('original_name_list', original_name_list)
	# print('posthoc_results_df', posthoc_results_df)
	# print('type posthoc_results_df', type(posthoc_results_df))


	val_col_temp = [assigned_val_for_missing_name]*len(original_name_list)
	# print('val_col_temp', val_col_temp)

	# print('posthoc_results_df.columns\n', posthoc_results_df.columns)
	# print('posthoc_results_df.index\n', posthoc_results_df.index)

	new_posthoc_results_df = pd.DataFrame(columns = original_name_list)
	for col in posthoc_results_df.columns:
		# print(col)
		for beh in posthoc_results_df.index:
			# print('beh', beh)
			idx_beh=original_name_list.index(beh)
			# print('idx_beh', idx_beh)

			val_col_temp[idx_beh]=posthoc_results_df[col][beh]

		new_posthoc_results_df[col]=val_col_temp
		val_col_temp = [assigned_val_for_missing_name]*len(original_name_list)

	for col in missing_name_list:
		idx_col=original_name_list.index(col)
		val_col_temp[idx_col]=float(1)
		# print('val_col_temp',val_col_temp)
		new_posthoc_results_df[col]=val_col_temp
		val_col_temp = [assigned_val_for_missing_name]*len(original_name_list)

	new_posthoc_results_df=new_posthoc_results_df.set_axis(original_name_list, axis='index')

	# val_diagonal=posthoc_results_df['Bsl']['Bsl']
	# print('val_diagonal', val_diagonal)
	# for col in new_posthoc_results_df.columns:
	# 	new_posthoc_results_df.loc[col,col]=val_diagonal


	# print('new_posthoc_results_df\n', new_posthoc_results_df)
	# print('type new_posthoc_results_df', type(new_posthoc_results_df))


	return new_posthoc_results_df



def trim_smthRow_of_list(list2d, startIdx=0, endIdx=-1, fps=1500, smth_window_s=0):
    trim_2dlist=[]
    for row in list2d:
        #print('len row', len(row))
        if startIdx<len(row)-1:
            trim_row=row[startIdx:endIdx]
            # trim_row= sync_utils.smooth_data(trim_row, windowlen=int(fps*smth_window_s))
            #print('len trim_row', len(trim_row))
            trim_2dlist.append(trim_row)

    return trim_2dlist


def overlap_btwn_two_list(short_list, long_list):

	common_element_list=np.intersect1d(short_list, long_list)

	# common_element_list=list(set(short_list).intersection(long_list))
	overlap_perc=len(common_element_list)/len(short_list)

	# print('overlap_perc', overlap_perc)

	return overlap_perc






def save_data_for_dFF_matrix(filename='data_for_dFF_matrix_dic.p'):

	data_for_dFF_mat_dic={}

	data_for_dFF_mat_dic.update({'y_list_allBeh':y_list_allBeh})
	data_for_dFF_mat_dic.update({'reordered_ROI_ID_list':reordered_ROI_ID_list})
	data_for_dFF_mat_dic.update({'reordered_mean_dFF_list':reordered_mean_dFF_list})
	data_for_dFF_mat_dic.update({'reordered_mean_dFF_01_list':reordered_mean_dFF_01_list})
	data_for_dFF_mat_dic.update({'reordered_masked_mean_dFF_list':reordered_masked_mean_dFF_list})
	data_for_dFF_mat_dic.update({'reordered_masked_mean_dFF_01_list':reordered_masked_mean_dFF_01_list})
	data_for_dFF_mat_dic.update({'reordered_p_value_list':reordered_p_value_list})
	data_for_dFF_mat_dic.update({'bin_p_value_mask':bin_p_value_mask})
	data_for_dFF_mat_dic.update({'reordered_nonoverlap_list':reordered_nonoverlap_list})
	data_for_dFF_mat_dic.update({'reordered_normality_list':reordered_normality_list})
	data_for_dFF_mat_dic.update({'reordered_std_list':reordered_std_list})


	pickle.dump( data_for_dFF_mat_dic, open( group_dFF_summary_dir+'/'+filename, "wb" ) )



	return



##main##

NAS_Dir=general_utils.NAS_Dir
NAS_AN_Proj_Dir=general_utils.NAS_AN_Proj_Dir
workstation_dir=general_utils.workstation_dir



experiments_group_per_fly=general_utils.group_expList_per_fly(experiments)
print('experiments_group_per_fly', experiments_group_per_fly)





ROI_id_list=[] #[Gal4-neuron (1 x n)]
mean_dFF_list=[]
mean_dFF_01_list=[]
mean_dFF_norm_list=[]
mean_dFF_log_list=[]
p_value_list=[]
nonoverlap_list=[]
bsl_mean_eachFly_list=[]
bsl_mean_01_eachFly_list=[]

normality_list=[]
std_list=[]

count_ROI=0
count_fly=0
count_recrd=0
num_ROI_ls=[]

print('\n#counting ROIs#\n')
for exp_lists_per_fly in experiments_group_per_fly:

	print('exp_lists_per_fly', exp_lists_per_fly)

	count_recrd+=len(exp_lists_per_fly)

	for date, genotype, fly, recrd_num in exp_lists_per_fly:


		Gal4=genotype.split('-')[0]
		fly_beh=fly[0].upper()+fly[1:]

		outDir_AN_recrd=NAS_AN_Proj_Dir+Gal4+'/2P/'+date+'/'+genotype+'-'+fly+'/'+genotype+'-'+fly+'-'+recrd_num+'/output'
		print('outDir_AN_recrd', outDir_AN_recrd)


		# Beh_Jpos_GC_DicData=general_utils.open_Beh_Jpos_GC_DicData(outDir_AN_recrd, 'SyncDic_7CamBeh_GC-RES.p')
		# Beh_Jpos_GC_DicData=general_utils.open_Beh_Jpos_GC_DicData(outDir_AN_recrd, 'SyncDic_7CamBeh20210507_GC-RES.p')
		# Beh_Jpos_GC_DicData=general_utils.open_Beh_Jpos_GC_DicData(outDir_AN_recrd, 'SyncDic_7CamBeh20210507_WalkRest_by_ball_GC-RES.p')
		Beh_Jpos_GC_DicData=general_utils.open_Beh_Jpos_GC_DicData(outDir_AN_recrd, 'SyncDic_7CamBeh_BW_20210619_GC-RES.p')


		GC_set = Beh_Jpos_GC_DicData['GCset']
	
		count_ROI+=len(GC_set)
		
		
		
		num_ROI_ls.append(len(GC_set))	

		for neuron_ID in range(len(GC_set)):
			ROI_id_list.append(Gal4+'-ROI#'+str(neuron_ID))


		break



	count_fly+=1


print('There are ', count_fly, 'flys in total...')
print('There are ', count_ROI, 'ROIs in total...')
print('There are ', count_recrd, 'recordings in total...')

print('ROI_id_list', ROI_id_list)



manual_ROI_order_csv=pd.read_csv(workstation_dir+'from_florian/Ascending_analysis/output/row_order_manual.csv')
manual_ROI_order= manual_ROI_order_csv['x'].tolist()
# print('manual_ROI_order\n', manual_ROI_order)
# print('type manual_ROI_order', type(manual_ROI_order))




for exp_lists_per_fly in experiments_group_per_fly:
	print('Processing per fly for ... ', exp_lists_per_fly )


	Etho_time_fly_Dic={}
	Etho_time_fly_Dic['FW_evt']=[]
	Etho_time_fly_Dic['BW_evt']=[]
	Etho_time_fly_Dic['rest_evt']=[]
	Etho_time_fly_Dic['E_groom_evt']=[]
	Etho_time_fly_Dic['A_groom_evt']=[]
	Etho_time_fly_Dic['FL_groom_evt']=[]
	Etho_time_fly_Dic['HL_groom_evt']=[]
	Etho_time_fly_Dic['Abd_groom_evt']=[]
	Etho_time_fly_Dic['PER_evt']=[]
	Etho_time_fly_Dic['Push_evt']=[]
	Etho_time_fly_Dic['CO2puff_evt']=[]

	Etho_time_fly_Dic['F_groom_evt']=[]
	Etho_time_fly_Dic['H_groom_evt']=[]
	Etho_time_fly_Dic['SixLeg_Move_evt']=[]



	F_Walk_GCevt_fly=[]
	B_Walk_GCevt_fly=[]
	Rest_GCevt_fly=[]
	E_groom_GCevt_fly=[]
	A_groom_GCevt_fly=[]
	FL_groom_GCevt_fly=[]
	HL_groom_GCevt_fly=[]
	Abd_groom_GCevt_fly=[]
	PER_GCevt_fly=[]
	Push_GCevt_fly=[]
	CO2puff_GCevt_fly=[]

	F_groom_GCevt_fly=[]
	H_groom_GCevt_fly=[]
	SixLeg_Move_GCevt_fly=[]



	F_Walk_GCevtNorm01_fly=[]
	B_Walk_GCevtNorm01_fly=[]
	Rest_GCevtNorm01_fly=[]
	E_groom_GCevtNorm01_fly=[]
	A_groom_GCevtNorm01_fly=[]
	FL_groom_GCevtNorm01_fly=[]
	HL_groom_GCevtNorm01_fly=[]
	Abd_groom_GCevtNorm01_fly=[]
	PER_GCevtNorm01_fly=[]
	Push_GCevtNorm01_fly=[]
	CO2puff_GCevtNorm01_fly=[]

	F_groom_GCevtNorm01_fly=[]
	H_groom_GCevtNorm01_fly=[]
	SixLeg_Move_GCevtNorm01_fly=[]


	Walk_APevt_fly=[]
	Rest_APevt_fly=[]
	Groom_APevt_fly=[]


	GC_gapfree_fly=[]
	GC_all_fly_perROI=[]
	GCnorm01_all_fly_perROI=[]

	AP_gapfree_fly=[]

	maxGC_fly=[]
	minGC_fly=[]

	recrdcount=0
	i=0
	#todo: append individual GC event into these array for summarize all data

	### Find the max GC per fly
	GC_fly_preSum=[]
	for date, genotype, fly, recrd_num in exp_lists_per_fly:
		Gal4=genotype.split('-')[0]
		dataDir = NAS_AN_Proj_Dir + Gal4 +'/2P/' + date+'/'+genotype+'-'+fly+'/'+genotype+'-'+fly+'-'+recrd_num + '/'
		outDirGC6_axoid = dataDir+'/output/GC6_auto/final/'
		GC_set = general_utils.readGCfile(outDirGC6_axoid)

		data_rawGC_freq=4

		if len(GC_fly_preSum)<len(GC_set):
			for roi_i, trace in enumerate(GC_set):
				GC_fly_preSum.append([])
		for roi_i, trace in enumerate(GC_set):
			GC_trace=math_utils.smooth_data(trace, windowlen=int(data_rawGC_freq*0.7)) # original = 0.7
			GC_fly_preSum[roi_i].extend(GC_trace)


	maxGC_fly=np.nanmax(GC_fly_preSum, axis=1)
	print('maxGC_fly', maxGC_fly)






	for date, genotype, fly, recrd_num in exp_lists_per_fly:

		Gal4=genotype.split('-')[0]

		flyDir = NAS_AN_Proj_Dir + Gal4 +'/2P/'+ date+'/'+genotype+'-'+fly+'/'
		dataDir = NAS_AN_Proj_Dir + Gal4 +'/2P/' + date+'/'+genotype+'-'+fly+'/'+genotype+'-'+fly+'-'+recrd_num + '/'
		pathForDic = dataDir+'/output/'

		print('dataDir', dataDir)

		outDirEvents = pathForDic + 'Events_df3d_and_ball_beh_BW_20210619/'
		if not os.path.exists(outDirEvents):
			os.makedirs(outDirEvents)
		outDirEventsSumFly = flyDir + 'Summary_Fly/Events_df3d_and_ball_beh_BW_20210619/'
		if not os.path.exists(outDirEventsSumFly):
			os.makedirs(outDirEventsSumFly)




		Beh_Jpos_GC_DicData=general_utils.open_Beh_Jpos_GC_DicData(pathForDic, 'SyncDic_7CamBeh_BW_20210619_GC-RES.p')



		GC_set = Beh_Jpos_GC_DicData['GCset']

		rest = Beh_Jpos_GC_DicData['rest']
		f_walk = Beh_Jpos_GC_DicData['forward_walk']
		b_walk = Beh_Jpos_GC_DicData['backward_walk']
		eye_groom = Beh_Jpos_GC_DicData['eye_groom']
		antennae_groom = Beh_Jpos_GC_DicData['antennae_groom']
		foreleg_groom = Beh_Jpos_GC_DicData['foreleg_groom']
		hindleg_groom = Beh_Jpos_GC_DicData['hindleg_groom']
		Abd_groom = Beh_Jpos_GC_DicData['Abd_groom']
		Push = Beh_Jpos_GC_DicData['Push']	
		PER=Beh_Jpos_GC_DicData['PER']

		CO2puff = Beh_Jpos_GC_DicData['CO2puff']
		PER_exten_len = Beh_Jpos_GC_DicData['PER_exten_len']

		F_groom = Beh_Jpos_GC_DicData['F_groom']
		H_groom = Beh_Jpos_GC_DicData['H_groom']
		print('type H_groom', type(H_groom))

		F_groom = np.asarray(eye_groom)+np.asarray(foreleg_groom)
		SixLeg_move = np.asarray(f_walk)+np.asarray(b_walk)+np.asarray(Push)
		print('type SixLeg_move', type(SixLeg_move))
		


		timeSec = Beh_Jpos_GC_DicData['timeSec']
		velForw_mm = Beh_Jpos_GC_DicData['velForw']
		velSide_mm = Beh_Jpos_GC_DicData['velSide']
		velTurn_deg = Beh_Jpos_GC_DicData['velTurn']



		Etho_Timesec_Dic = Beh_Jpos_GC_DicData['Etho_Timesec_Dic']
		Etho_Idx_Dic = Beh_Jpos_GC_DicData['Etho_Idx_Dic']

		data_freq=len(rest)/timeSec[-1]
		print('data_freq', data_freq)


		filter_kernel_s=0.5
		# filter_kernel_s=0.06
		F_Walk_LP=math_utils.hysteresis_filter(f_walk, n=int(data_freq*filter_kernel_s))*1
		B_Walk_LP=math_utils.hysteresis_filter(b_walk, n=int(data_freq*filter_kernel_s))*1
		Rest_LP=math_utils.hysteresis_filter(rest, n=int(data_freq*filter_kernel_s))*1
		E_Groom_LP=math_utils.hysteresis_filter(eye_groom, n=int(data_freq*filter_kernel_s))*1
		A_Groom_LP=math_utils.hysteresis_filter(antennae_groom, n=int(data_freq*filter_kernel_s))*1
		FL_Groom_LP=math_utils.hysteresis_filter(foreleg_groom, n=int(data_freq*filter_kernel_s))*1
		HL_Groom_LP=math_utils.hysteresis_filter(hindleg_groom, n=int(data_freq*filter_kernel_s))*1
		Abd_Groom_LP=math_utils.hysteresis_filter(Abd_groom, n=int(data_freq*filter_kernel_s))*1
		Push_LP=math_utils.hysteresis_filter(Push, n=int(data_freq*filter_kernel_s))*1

		H_Groom_LP=math_utils.hysteresis_filter(H_groom, n=int(data_freq*filter_kernel_s))*1
		F_Groom_LP=math_utils.hysteresis_filter(F_groom, n=int(data_freq*filter_kernel_s))*1	
		SixLeg_move_LP = math_utils.hysteresis_filter(SixLeg_move, n=int(data_freq*filter_kernel_s))*1	





		idx_f_walk_evt, timesec_f_walk_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(F_Walk_LP, timeSec)
		idx_b_walk_evt, timesec_b_walk_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(B_Walk_LP, timeSec)
		idx_rest_evt, timesec_rest_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(Rest_LP, timeSec)
		idx_E_groom_evt, timesec_E_groom_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(E_Groom_LP, timeSec)
		idx_A_groom_evt, timesec_A_groom_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(A_Groom_LP, timeSec)
		idx_FL_groom_evt, timesec_FL_groom_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(FL_Groom_LP, timeSec)
		idx_HL_groom_evt, timesec_HL_groom_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(HL_Groom_LP, timeSec)
		idx_Abd_groom_evt, timesec_Abd_groom_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(Abd_Groom_LP, timeSec)
		idx_Push_evt, timesec_Push_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(Push_LP, timeSec)
		idx_PER_evt, timesec_PER_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(PER, timeSec)
		idx_CO2puff_evt, timesec_CO2puff_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(CO2puff, timeSec)

		idx_H_groom_evt, timesec_H_groom_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(H_Groom_LP, timeSec)	
		idx_F_groom_evt, timesec_F_groom_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(F_Groom_LP, timeSec)
		idx_SixLeg_move_evt, timesec_SixLeg_move_evt=sync_utils.Calculate_idx_time_for_bin_beh_trace(SixLeg_move_LP, timeSec)



		EthoLP_Idx_Dic={}
		EthoLP_Idx_Dic.update({'rest_evt':idx_rest_evt})
		EthoLP_Idx_Dic.update({'FW_evt':idx_f_walk_evt})
		EthoLP_Idx_Dic.update({'BW_evt':idx_b_walk_evt})
		EthoLP_Idx_Dic.update({'E_groom_evt':idx_E_groom_evt})
		EthoLP_Idx_Dic.update({'A_groom_evt':idx_A_groom_evt})
		EthoLP_Idx_Dic.update({'FL_groom_evt':idx_FL_groom_evt})
		EthoLP_Idx_Dic.update({'HL_groom_evt':idx_HL_groom_evt})
		EthoLP_Idx_Dic.update({'Abd_groom_evt':idx_Abd_groom_evt})
		EthoLP_Idx_Dic.update({'PER_evt':idx_PER_evt})
		EthoLP_Idx_Dic.update({'Push_evt':idx_Push_evt})
		EthoLP_Idx_Dic.update({'CO2puff_evt':idx_CO2puff_evt})

		EthoLP_Idx_Dic.update({'F_groom_evt':idx_F_groom_evt})
		EthoLP_Idx_Dic.update({'H_groom_evt':idx_H_groom_evt})	
		EthoLP_Idx_Dic.update({'SixLeg_Move_evt':idx_SixLeg_move_evt})	


		EthoLP_Timesec_Dic={}
		EthoLP_Timesec_Dic.update({'rest_evt':timesec_rest_evt})
		EthoLP_Timesec_Dic.update({'FW_evt':timesec_f_walk_evt})
		EthoLP_Timesec_Dic.update({'BW_evt':timesec_b_walk_evt})
		EthoLP_Timesec_Dic.update({'E_groom_evt':timesec_E_groom_evt})
		EthoLP_Timesec_Dic.update({'A_groom_evt':timesec_A_groom_evt})
		EthoLP_Timesec_Dic.update({'FL_groom_evt':timesec_FL_groom_evt})
		EthoLP_Timesec_Dic.update({'HL_groom_evt':timesec_HL_groom_evt})
		EthoLP_Timesec_Dic.update({'Abd_groom_evt':timesec_Abd_groom_evt})
		EthoLP_Timesec_Dic.update({'PER_evt':timesec_PER_evt})
		EthoLP_Timesec_Dic.update({'Push_evt':timesec_Push_evt})
		EthoLP_Timesec_Dic.update({'CO2puff_evt':timesec_CO2puff_evt})

		EthoLP_Timesec_Dic.update({'F_groom_evt':timesec_F_groom_evt})
		EthoLP_Timesec_Dic.update({'H_groom_evt':timesec_H_groom_evt})
		EthoLP_Timesec_Dic.update({'SixLeg_Move_evt':timesec_SixLeg_move_evt})		


		Plot_overlay_each_recrd=False
		if Plot_overlay_each_recrd==True:
			plot_utils.Plot_whole_trace(GC_set, CO2puff, PER_exten_len, timeSec, velForw_mm, velSide_mm, velTurn_deg, EthoLP_Timesec_Dic, 'whole_trace_7CamBeh_BW_20210619', filepath=outDirEvents)


		GC_set_norm01=[]

		for ROI_i, GC_trace in enumerate(GC_set):
			GC_trace=math_utils.smooth_data(GC_trace, windowlen=int(data_freq*0.7)) # original = 0.7
			#GC_trace=math_utils.norm_to_max(GC_trace, percentile_th_to_norm=100)
			GC_trace=math_utils.norm_to_val(GC_trace, val=maxGC_fly[ROI_i])
			GC_set_norm01.append(GC_trace)

		# GC_set=GC_set_norm01



		bsl_s=1 #s
		event_dur_s=3 #s

		GC_rest_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'rest_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_f_walk_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'FW_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_b_walk_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'BW_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_E_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'E_groom_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_A_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'A_groom_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_FL_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'FL_groom_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_HL_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'HL_groom_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_Abd_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'Abd_groom_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_Push_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'Push_evt', GC_set, baseline=bsl_s, fps=data_freq)

		GC_PER_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'PER_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_CO2puff_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'CO2puff_evt', GC_set, baseline=bsl_s, fps=data_freq)

		GC_F_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'F_groom_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_H_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'H_groom_evt', GC_set, baseline=bsl_s, fps=data_freq)
		GC_SixLeg_Move_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'SixLeg_Move_evt', GC_set, baseline=bsl_s, fps=data_freq)


		GCnorm01_rest_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'rest_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_f_walk_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'FW_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_b_walk_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'BW_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_E_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'E_groom_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_A_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'A_groom_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_FL_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'FL_groom_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_HL_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'HL_groom_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_Abd_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'Abd_groom_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_Push_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'Push_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)

		GCnorm01_PER_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'PER_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_CO2puff_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'CO2puff_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)

		GCnorm01_F_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'F_groom_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_H_groom_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'H_groom_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)
		GCnorm01_SixLeg_Move_evt_set, _ =general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'SixLeg_Move_evt', GC_set_norm01, baseline=bsl_s, fps=data_freq)

		# AP_walk_evt_set, _ = general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'walk_evt', [velForw_mm], baseline=bsl_s, fps=data_freq)
		# AP_rest_evt_set, _ = general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'rest_evt', [velForw_mm], baseline=bsl_s, fps=data_freq)
		# AP_groom_evt_set, _ = general_utils.find_corresponding_evt_from_groupIdxs(EthoLP_Idx_Dic, 'F_groom_evt', [velForw_mm], baseline=bsl_s, fps=data_freq)


		# print('shape GC_f_walk_evt_set', np.shape(GC_f_walk_evt_set))
		# print('shape GC_rest_evt_set', np.shape(GC_rest_evt_set))
		# print('shape GC_E_groom_evt_set', np.shape(GC_E_groom_evt_set))


		for neuron_idx, GC_per_neruon in enumerate(GC_set):
			if len(F_Walk_GCevt_fly)!=len(GC_set): 
				F_Walk_GCevt_fly.append([])
				F_Walk_GCevtNorm01_fly.append([])

			if len(B_Walk_GCevt_fly)!=len(GC_set): 
				B_Walk_GCevt_fly.append([])
				B_Walk_GCevtNorm01_fly.append([])

			if len(Rest_GCevt_fly)!=len(GC_set):
				Rest_GCevt_fly.append([])
				Rest_GCevtNorm01_fly.append([])

			if len(E_groom_GCevt_fly)!=len(GC_set):
				E_groom_GCevt_fly.append([])
				E_groom_GCevtNorm01_fly.append([])

			if len(A_groom_GCevt_fly)!=len(GC_set): 
				A_groom_GCevt_fly.append([])
				A_groom_GCevtNorm01_fly.append([])

			if len(FL_groom_GCevt_fly)!=len(GC_set):
				FL_groom_GCevt_fly.append([])
				FL_groom_GCevtNorm01_fly.append([])

			if len(HL_groom_GCevt_fly)!=len(GC_set):
				HL_groom_GCevt_fly.append([])
				HL_groom_GCevtNorm01_fly.append([])

			if len(Abd_groom_GCevt_fly)!=len(GC_set):
				Abd_groom_GCevt_fly.append([])
				Abd_groom_GCevtNorm01_fly.append([])

			if len(PER_GCevt_fly)!=len(GC_set):
				PER_GCevt_fly.append([])
				PER_GCevtNorm01_fly.append([])

			if len(Push_GCevt_fly)!=len(GC_set):
				Push_GCevt_fly.append([])
				Push_GCevtNorm01_fly.append([])

			if len(CO2puff_GCevt_fly)!=len(GC_set):
				CO2puff_GCevt_fly.append([])
				CO2puff_GCevtNorm01_fly.append([])

			if len(F_groom_GCevt_fly)!=len(GC_set): 
				F_groom_GCevt_fly.append([])
				F_groom_GCevtNorm01_fly.append([])

			if len(H_groom_GCevt_fly)!=len(GC_set):
				H_groom_GCevt_fly.append([])
				H_groom_GCevtNorm01_fly.append([])

			if len(SixLeg_Move_GCevt_fly)!=len(GC_set):
				SixLeg_Move_GCevt_fly.append([])
				SixLeg_Move_GCevtNorm01_fly.append([])


			if len(GC_all_fly_perROI)!=len(GC_set):
				GC_all_fly_perROI.append([])			
				GCnorm01_all_fly_perROI.append([])	

			for evt_idx, GC_evt in enumerate(GC_f_walk_evt_set[neuron_idx]):
				F_Walk_GCevt_fly[neuron_idx].append(GC_evt)
				F_Walk_GCevtNorm01_fly[neuron_idx].append(GCnorm01_f_walk_evt_set[neuron_idx][evt_idx])            

			for evt_idx, GC_evt in enumerate(GC_b_walk_evt_set[neuron_idx]):
				B_Walk_GCevt_fly[neuron_idx].append(GC_evt)  
				B_Walk_GCevtNorm01_fly[neuron_idx].append(GCnorm01_b_walk_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_rest_evt_set[neuron_idx]):
				Rest_GCevt_fly[neuron_idx].append(GC_evt)
				Rest_GCevtNorm01_fly[neuron_idx].append(GCnorm01_rest_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_E_groom_evt_set[neuron_idx]):
				E_groom_GCevt_fly[neuron_idx].append(GC_evt)
				E_groom_GCevtNorm01_fly[neuron_idx].append(GCnorm01_E_groom_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_A_groom_evt_set[neuron_idx]):
				A_groom_GCevt_fly[neuron_idx].append(GC_evt)            
				A_groom_GCevtNorm01_fly[neuron_idx].append(GCnorm01_A_groom_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_FL_groom_evt_set[neuron_idx]):
				FL_groom_GCevt_fly[neuron_idx].append(GC_evt)
				FL_groom_GCevtNorm01_fly[neuron_idx].append(GCnorm01_FL_groom_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_HL_groom_evt_set[neuron_idx]):
				HL_groom_GCevt_fly[neuron_idx].append(GC_evt)
				HL_groom_GCevtNorm01_fly[neuron_idx].append(GCnorm01_HL_groom_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_Abd_groom_evt_set[neuron_idx]):
				Abd_groom_GCevt_fly[neuron_idx].append(GC_evt)
				Abd_groom_GCevtNorm01_fly[neuron_idx].append(GCnorm01_Abd_groom_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_PER_evt_set[neuron_idx]):
				PER_GCevt_fly[neuron_idx].append(GC_evt)    
				PER_GCevtNorm01_fly[neuron_idx].append(GCnorm01_PER_evt_set[neuron_idx][evt_idx])         

			for evt_idx, GC_evt in enumerate(GC_Push_evt_set[neuron_idx]):
				Push_GCevt_fly[neuron_idx].append(GC_evt)
				Push_GCevtNorm01_fly[neuron_idx].append(GCnorm01_Push_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_CO2puff_evt_set[neuron_idx]):
				CO2puff_GCevt_fly[neuron_idx].append(GC_evt)
				CO2puff_GCevtNorm01_fly[neuron_idx].append(GCnorm01_CO2puff_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_F_groom_evt_set[neuron_idx]):
				F_groom_GCevt_fly[neuron_idx].append(GC_evt)
				F_groom_GCevtNorm01_fly[neuron_idx].append(GCnorm01_F_groom_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_H_groom_evt_set[neuron_idx]):
				H_groom_GCevt_fly[neuron_idx].append(GC_evt)
				H_groom_GCevtNorm01_fly[neuron_idx].append(GCnorm01_H_groom_evt_set[neuron_idx][evt_idx]) 

			for evt_idx, GC_evt in enumerate(GC_SixLeg_Move_evt_set[neuron_idx]):
				SixLeg_Move_GCevt_fly[neuron_idx].append(GC_evt)
				SixLeg_Move_GCevtNorm01_fly[neuron_idx].append(GCnorm01_SixLeg_Move_evt_set[neuron_idx]) 



		for ROI_i, GC_trace in enumerate(GC_set):
			GC_trace=math_utils.smooth_data(GC_trace, windowlen=int(data_freq*0.7))
			GC_all_fly_perROI[ROI_i].append(GC_trace)

			GCnorm01_all_fly_perROI[ROI_i].append(GC_set_norm01[ROI_i])



		Etho_time_fly_Dic['FW_evt'].extend(EthoLP_Timesec_Dic['FW_evt'])
		Etho_time_fly_Dic['BW_evt'].extend(EthoLP_Timesec_Dic['BW_evt'])
		Etho_time_fly_Dic['rest_evt'].extend(EthoLP_Timesec_Dic['rest_evt'])
		Etho_time_fly_Dic['E_groom_evt'].extend(EthoLP_Timesec_Dic['E_groom_evt'])
		Etho_time_fly_Dic['A_groom_evt'].extend(EthoLP_Timesec_Dic['A_groom_evt'])
		Etho_time_fly_Dic['FL_groom_evt'].extend(EthoLP_Timesec_Dic['FL_groom_evt'])
		Etho_time_fly_Dic['HL_groom_evt'].extend(EthoLP_Timesec_Dic['HL_groom_evt'])
		Etho_time_fly_Dic['Abd_groom_evt'].extend(EthoLP_Timesec_Dic['Abd_groom_evt'])
		Etho_time_fly_Dic['PER_evt'].extend(EthoLP_Timesec_Dic['PER_evt'])
		Etho_time_fly_Dic['Push_evt'].extend(EthoLP_Timesec_Dic['Push_evt'])
		Etho_time_fly_Dic['CO2puff_evt'].extend(EthoLP_Timesec_Dic['CO2puff_evt'])

		Etho_time_fly_Dic['F_groom_evt'].extend(EthoLP_Timesec_Dic['F_groom_evt'])
		Etho_time_fly_Dic['H_groom_evt'].extend(EthoLP_Timesec_Dic['H_groom_evt'])
		Etho_time_fly_Dic['SixLeg_Move_evt'].extend(EthoLP_Timesec_Dic['SixLeg_Move_evt'])




		GC_gapfree=[]
		for neuron_idx, GC in enumerate(GC_set):
			GC_gapfree.extend(GC)
			GC_gapfree_fly.extend(GC)

		# maxGC=sorted(GC_gapfree)[int(len(GC_gapfree)*0.99)] 
		# minGC=sorted(GC_gapfree)[int(len(GC_gapfree)*0.01)]

		maxGC=np.nanmax(GC_gapfree)
		minGC=np.nanmin(GC_gapfree)

		GC_lim=[minGC, maxGC]

		Plot_overlay_each_recrd=False
		if Plot_overlay_each_recrd==True:

			## Plot event overlay per behavior class for individual recording
			if len(GC_f_walk_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_f_walk_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='FW_evt', filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_b_walk_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_b_walk_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='BW_evt', filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_rest_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_rest_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='rest_evt',filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_E_groom_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_E_groom_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='E_groom_evt',filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_A_groom_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_A_groom_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='A_groom_evt', filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_FL_groom_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_FL_groom_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='FL_groom_evt',filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_HL_groom_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_HL_groom_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='HL_groom_evt',filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_Abd_groom_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_Abd_groom_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='Abd_groom_evt', filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_PER_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_PER_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='PER_evt',filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_Push_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_Push_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='Push_evt',filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_CO2puff_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_CO2puff_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='CO2puff_evt',filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_F_groom_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_F_groom_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='F_groom_evt',filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_H_groom_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_H_groom_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='H_groom_evt',filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)
			if len(GC_SixLeg_Move_evt_set[0])>0:
				plot_utils.Plot_Evtavg_overlay(EthoLP_Timesec_Dic, GC_SixLeg_Move_evt_set, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='SixLeg_Move_evt',filename=date+'-'+genotype+'-'+fly+'-'+recrd_num, filepath=outDirEvents)


		GC_gapfree=[]

		# AP_gapfree_fly.extend(velForw_mm[10000:])




	maxGC=np.nanmax(GC_gapfree_fly)
	minGC=np.nanmin(GC_gapfree_fly)

	# maxGC=sorted(GC_gapfree_fly)[int(len(GC_gapfree_fly)*0.99)]
	# minGC=sorted(GC_gapfree_fly)[int(len(GC_gapfree_fly)*0.01)]

	GC_lim=[minGC, maxGC]
	GC_lim=[minGC, 100]

	# maxAP=np.nanmax(AP_gapfree_fly)
	# minAP=np.nanmin(AP_gapfree_fly)
	# AP_lim=[minAP, maxAP]

	Plot_overlay_each_fly=True
	if Plot_overlay_each_fly==True:

		## Plot event overlay per behavior class for individual fly
		if len(F_Walk_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, F_Walk_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='FW_evt', filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(B_Walk_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, B_Walk_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='BW_evt', filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(Rest_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, Rest_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='rest_evt',filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(E_groom_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, E_groom_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='E_groom_evt',filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(A_groom_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, A_groom_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='A_groom_evt', filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(FL_groom_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, FL_groom_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='FL_groom_evt',filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(HL_groom_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, HL_groom_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='HL_groom_evt',filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(Abd_groom_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, Abd_groom_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='Abd_groom_evt', filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(PER_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, PER_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='PER_evt',filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(Push_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, Push_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='Push_evt',filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(CO2puff_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, CO2puff_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='CO2puff_evt',filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(F_groom_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, F_groom_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='F_groom_evt',filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(H_groom_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, H_groom_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='H_groom_evt',filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)
		if len(SixLeg_Move_GCevt_fly[0])>0:
			plot_utils.Plot_Evtavg_overlay_err(Etho_time_fly_Dic, SixLeg_Move_GCevt_fly, baseline=bsl_s, epoch_len=event_dur_s, y_lim=GC_lim, whichBeh='SixLeg_Move_evt',filename=date+'-'+genotype+'-'+fly, filepath=outDirEventsSumFly)



 

	print('shape GC_all_fly_perROI', np.shape(GC_all_fly_perROI))
	print('shape GCnorm01_all_fly_perROI', np.shape(GCnorm01_all_fly_perROI))
	print('shape GC_gapfree_fly', np.shape(GC_gapfree_fly))

	
	continue


	##Handling the statistics comapring behavioral dFF vs baseline dFF
	for ROI_i, evts in enumerate(F_Walk_GCevt_fly):

		print('ROI_i', ROI_i)

		

		p_value_fold=[[] for i in range(11)] # 11 types of behavior in total
		# print('p_value_fold', p_value_fold)

		stats_posthoc_fold=6
		for fold in range(0,stats_posthoc_fold):

			beh_name_list=[
			'FW',
			'BW',
			'Push',
			'Rest',
			'EG',
			'AG',
			'FLG',
			'AbdG',
			'HLG',
			'PE',
			'CO2',
			# 'FG',
			# 'HG',
			]
			# print('len beh_name_list', len(beh_name_list))


			Beh_Evts_set_perROI=[
			F_Walk_GCevt_fly[ROI_i],
			B_Walk_GCevt_fly[ROI_i],
			Push_GCevt_fly[ROI_i],
			Rest_GCevt_fly[ROI_i],
			E_groom_GCevt_fly[ROI_i],
			A_groom_GCevt_fly[ROI_i],
			FL_groom_GCevt_fly[ROI_i],
			Abd_groom_GCevt_fly[ROI_i],
			HL_groom_GCevt_fly[ROI_i],
			PER_GCevt_fly[ROI_i],
			CO2puff_GCevt_fly[ROI_i],
			# F_groom_GCevt_fly[ROI_i],
			# H_groom_GCevt_fly[ROI_i],
			]

			# print('shape Beh_Evts_set_perROI', np.shape(Beh_Evts_set_perROI[0][0]))

			Beh_Norm01Evts_set_perROI=[
			F_Walk_GCevtNorm01_fly[ROI_i],
			B_Walk_GCevtNorm01_fly[ROI_i],
			Push_GCevtNorm01_fly[ROI_i],
			Rest_GCevtNorm01_fly[ROI_i],
			E_groom_GCevtNorm01_fly[ROI_i],
			A_groom_GCevtNorm01_fly[ROI_i],
			FL_groom_GCevtNorm01_fly[ROI_i],
			Abd_groom_GCevtNorm01_fly[ROI_i],
			HL_groom_GCevtNorm01_fly[ROI_i],
			PER_GCevtNorm01_fly[ROI_i],
			CO2puff_GCevtNorm01_fly[ROI_i],
			# F_groom_GCevtNorm01_fly[ROI_i],
			# H_groom_GCevtNorm01_fly[ROI_i],
			]

			# print('PER_GCevtNorm01_fly[0]', PER_GCevtNorm01_fly[0])
			# print('shape PER_GCevtNorm01_fly[0]', np.shape(PER_GCevtNorm01_fly[0]))




			# bsl_datapoints, bsl_evts=math_utils.find_baseline_of_trace_by_runningWindow(GC_all_fly_perROI[ROI_i], outDirEventsSumFly, 'thresholding_bsl_runWindow'+'ROI_'+str(ROI_i), samplerate=data_freq)
			_, bsl_evts, _=math_utils.find_baseline_of_trace01(GC_all_fly_perROI[ROI_i], GC_all_fly_perROI[ROI_i], outDirEventsSumFly, 'thresholding_bsl_without01'+'ROI_'+str(ROI_i), samplerate=data_freq)
			bsl_datapoints_01, bsl_evts_01, normal=math_utils.find_baseline_of_trace01(GCnorm01_all_fly_perROI[ROI_i], GC_all_fly_perROI[ROI_i], outDirEventsSumFly, 'thresholding_bsl_within01_'+'ROI_'+str(ROI_i), samplerate=data_freq)

			# bsl_datapoints, bsl_evts=math_utils.find_baseline_of_trace_by_thresholding_histogram(GC_gapfree_fly_perROI[ROI_i], outDirEventsSumFly, 'thresholding_bsl'+'ROI_'+str(ROI_i), samplerate=data_freq, descent=False)
			# bsl_datapoints=general_utils.downsampling_trace(bsl_datapoints, data_freq/4.3)
			# bsl_datapoints = sync_utils.smooth_data(bsl_datapoints, windowlen=int(4.3*0.3))

			## Handling the raw dFF baseline mean
			bsl_evtmean_list=math_utils.compute_mean_with_diffrerent_row_length(bsl_evts, samplerate=data_freq, cutting_head_s=1)
			bsl_mean_ROI=np.nanmean(bsl_evtmean_list)
			bsl_mean_eachFly_list.append([bsl_mean_ROI])


			cutting_head_s=0.7
			min_evt_num=10

			print('\n----',date, genotype, fly, recrd_num, 'ROI#', ROI_i, 'fold=', fold, '----\n')
			all_beh_meanlist_01_list, all_beh_mean_01, all_behdatapoint_nonoverlap_w_Bsl_01, resample_bsl_evtmean_01_list=stats_utils.prep_eachBehEvts_for_stat(\
				beh_name_list, Beh_Norm01Evts_set_perROI, bsl_evts_01, bsl_datapoints_01, min_evt_num=min_evt_num, bsl_s=bsl_s, cutting_head_s=cutting_head_s, data_freq=data_freq)
			
			
			## ---------------------------------------------------------------------------
			# ## Take thresholding baseline and its mean for statistics
			bsl_evtmean_01_list=resample_bsl_evtmean_01_list
			bsl_mean_01=np.nanmean(bsl_evtmean_01_list)
			bsl_mean_01_eachFly_list.append([bsl_mean_01])	
			# print('shape bsl_evtmean_01_list', np.shape(bsl_evtmean_01_list))
			# print('bsl_mean_01', bsl_mean_01)	

			all_beh_meanlist_01_list.append(bsl_evtmean_01_list)

			## ---------------------------------------------------------------------------





			beh_name_list.append('Bsl')
			all_beh_meanlist_list=all_beh_meanlist_01_list

			# print('len beh_name_list', len(beh_name_list))
			# print('beh_name_list', beh_name_list)
			print('len all_beh_meanlist_list', len(all_beh_meanlist_list))
			for i, v in enumerate(all_beh_meanlist_list):
				print('len v', len(v))



			##-------------------------------------------------------
			# posthoc_results_kw, p_value_KW=stats_utils.KW_w_dataset(all_beh_meanlist_list, beh_name_list)
			# p_value_whole=p_value_KW
			# stats_type='Krukal-Wallis test'
			# posthoc_results=posthoc_results_kw

			##-------------------------------------------------------
			posthoc_results_anova, p_value_ANOVA=stats_utils.anova_w_dataset(all_beh_meanlist_list, beh_name_list)
			p_value_whole=p_value_ANOVA		
			stats_type='ANOVA'
			posthoc_results=posthoc_results_anova
			##-------------------------------------------------------


			print('posthoc_results\n', posthoc_results)
			# print('p value against bsl\n', posthoc_results['Bsl'])

			for i, name in enumerate(beh_name_list):
				if name!='Bsl':
					p_value_fold[i].append(posthoc_results['Bsl'][beh_name_list[i]])
			print('fold=',fold,'p-values:\n',posthoc_results['Bsl'])

		print('p_value_fold', p_value_fold)

		
		
		##reassign the p-value against baseline after voting
		for i, name in enumerate(beh_name_list):
			if name!='Bsl':
				votes=p_value_fold[i]
				count_p_lessthan_005_posthoc=sum(map(lambda x : x<0.05, votes))
				print('->', count_p_lessthan_005_posthoc, 'out of', len(votes), '<0.05')


				if count_p_lessthan_005_posthoc<=len(votes)*0.5 and not np.isnan(posthoc_results.loc['Bsl',name]):

					print('-->', name, 'is not significantly drifferent from baseline')

					posthoc_results.loc['Bsl',name]=1
					# print('posthoc_results["Bsl"][name]', posthoc_results['Bsl'][name])

				elif count_p_lessthan_005_posthoc>len(votes)*0.5 and not np.isnan(posthoc_results.loc['Bsl',name]):
					print('-->', name, 'is significantly drifferent from baseline')
					# print('posthoc_results["Bsl"][name]', posthoc_results['Bsl'][name])
					p_value_fold=np.asarray(p_value_fold)
					print('p_value_fold', p_value_fold)
					if np.median(p_value_fold[p_value_fold<0.05])<=0.001:
						assign_p=0.000999
					else:
						assign_p=np.median(p_value_fold[p_value_fold<0.05])
					posthoc_results.loc['Bsl',name]=np.median(assign_p)

				else:
					print('-->', name, 'p-value is', posthoc_results.loc['Bsl',name])

		# print('type posthoc_results\n', type(posthoc_results))
		print('after re-assign based on the vote of p-value < 0.05. posthoc_results\n', posthoc_results)

	



		## re-organizing the posthoc_results by voting. If the vote support the no significance, then put the P-value of the behavior to 1 or 0.5 as a marker. 
		#Don't touch NaN where it is originally there.

		fig = plt.figure(facecolor='white', figsize=(8,5), dpi=170)
		cmap = ['0', '#FFAB91', '#9D00E5', '#C05EED',  '#E1AAFA']
		heatmap_args = {'cmap':cmap, 'linewidths': 0.2, 'linecolor': '1', 'clip_on': False, 'square': True, 'cbar_ax_bbox': [0.80, 0.35, 0.04, 0.3]}
		scikit_posthocs.sign_plot(x=posthoc_results, **heatmap_args)
		if p_value_whole>0.01:
			p_value_str=stats_type+', \np < 0.05'

		elif p_value_whole<0.01 and p_value_whole>0.001:
			p_value_str=stats_type+', \np < 0.01'

		elif p_value_whole<0.001:
			p_value_str=stats_type+', \np < 0.001'
			

		plt.text(-0.2,1.2,p_value_str, color='k', size=10)

		plt.savefig(outDirEventsSumFly + 'posthoc_compar_dFFofBeh_ROI_'+str(ROI_i)+'.png', edgecolor='none', transparent=False) #bbox_inches='tight', 
		plt.clf()
		plt.close() 



		

		all_beh_mean=[
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(F_Walk_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(B_Walk_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(Push_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(Rest_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(E_groom_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(A_groom_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(FL_groom_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(Abd_groom_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(HL_groom_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(PER_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(CO2puff_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		# np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(F_groom_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		# np.nanmean(math_utils.compute_mean_with_diffrerent_row_length(H_groom_GCevt_fly[ROI_i], samplerate=data_freq, cutting_head_s=1))-bsl_mean_ROI,
		]

		dFF_list_ROI=all_beh_mean
		dFF_01_list_ROI=all_beh_mean_01

		# dFF_list_ROI=[mean_F_Walk_evt, mean_B_Walk_evt, mean_Rest_evt, mean_E_groom_evt, mean_A_groom_evt, mean_FL_groom_evt, mean_HL_groom_evt, mean_Abd_groom_evt, mean_Push_evt, mean_PER_evt, mean_CO2puff_evt]
		#dFF_norm_list_ROI=[mean_norm_F_Walk_evt, mean_norm_B_Walk_evt, mean_norm_Rest_evt, mean_norm_E_groom_evt, mean_norm_A_groom_evt, mean_norm_FL_groom_evt, mean_norm_HL_groom_evt, mean_norm_Abd_groom_evt, mean_norm_Push_evt, mean_norm_PER_evt, mean_norm_CO2puff_evt]


		p_value_list_ROI=[]
		for i, name in enumerate(beh_name_list):
			if name!='Bsl':
				p_value_list_ROI.append(posthoc_results.loc['Bsl',name])

		


		print('p_value_list_ROI', p_value_list_ROI)



		# print('dFF_log_list_ROI', dFF_log_list_ROI)
		mean_dFF_list.append(dFF_list_ROI)
		mean_dFF_01_list.append(dFF_01_list_ROI)
		p_value_list.append(p_value_list_ROI)
		nonoverlap_list.append(all_behdatapoint_nonoverlap_w_Bsl_01)
		#mean_dFF_norm_list.append(dFF_norm_list_ROI)

		normality_list.append(normal)
		std_list.append(np.std(general_utils.flatten_list(GC_all_fly_perROI[ROI_i])))
		
		# mean_dFF_log_list.append(dFF_log_list_ROI)

		F_walk_evtmean_list=all_beh_meanlist_01_list[0]
		B_walk_evtmean_list=all_beh_meanlist_01_list[1]
		Push_evtmean_list=all_beh_meanlist_01_list[2]
		Rest_evtmean_list=all_beh_meanlist_01_list[3]
		E_groom_evtmean_list=all_beh_meanlist_01_list[4]
		A_groom_evtmean_list=all_beh_meanlist_01_list[5]
		FL_groom_evtmean_list=all_beh_meanlist_01_list[6]
		Abd_groom_evtmean_list=all_beh_meanlist_01_list[7]
		HL_groom_evtmean_list=all_beh_meanlist_01_list[8]
		PER_evtmean_list=all_beh_meanlist_01_list[9]
		CO2puff_evtmean_list=all_beh_meanlist_01_list[10]


		# hist_bsl, bin_edges_bsl = np.histogram(bsl_evtmean_01_list, density=True)
		
		# fig = plt.figure(facecolor='white', figsize=(20,4*5), dpi=170)

		# plt.subplot(5,1,1)
		# plt.title('FW, BW, Rest')
		# plt.hist(bsl_evtmean_01_list, len(bin_edges_bsl)-1, density=True, color='k', alpha=0.75)
		# if len(F_walk_evtmean_list)>min_evt_num:
		# 	hist_fw, bin_edges_fw = np.histogram(F_walk_evtmean_list, density=True)
		# 	plt.hist(F_walk_evtmean_list, len(bin_edges_fw)-1, density=True, color=plot_setting.FW_color, alpha=0.75)
		# if len(B_walk_evtmean_list)>min_evt_num:
		# 	hist_bw, bin_edges_bw = np.histogram(B_walk_evtmean_list, density=True)
		# 	plt.hist(B_walk_evtmean_list, len(bin_edges_bw)-1, density=True, color=plot_setting.BW_color, alpha=0.5)
		# if len(Rest_evtmean_list)>min_evt_num:
		# 	hist_r, bin_edges_r = np.histogram(Rest_evtmean_list, density=True)
		# 	plt.hist(Rest_evtmean_list, len(bin_edges_r)-1, density=True, color=plot_setting.rest_color, alpha=0.75)

		# plt.subplot(5,1,2)
		# plt.title('EG, AG, FLG')
		# plt.hist(bsl_evtmean_01_list, len(bin_edges_bsl)-1, density=True, color='k', alpha=0.75)
		# if len(E_groom_evtmean_list)>min_evt_num:
		# 	hist_eg, bin_edges_eg = np.histogram(E_groom_evtmean_list, density=True)
		# 	plt.hist(E_groom_evtmean_list, len(bin_edges_eg)-1, density=True, color=plot_setting.E_groom_color, alpha=0.75)
		# if len(A_groom_evtmean_list)>min_evt_num:
		# 	hist_ag, bin_edges_ag = np.histogram(A_groom_evtmean_list, density=True)
		# 	plt.hist(A_groom_evtmean_list, len(bin_edges_ag)-1, density=True, color=plot_setting.A_groom_color, alpha=0.75)
		# if len(FL_groom_evtmean_list)>min_evt_num:
		# 	hist_flg, bin_edges_flg = np.histogram(FL_groom_evtmean_list, density=True)
		# 	plt.hist(FL_groom_evtmean_list, len(bin_edges_flg)-1, density=True, color=plot_setting.FL_groom_color, alpha=0.75)

		# plt.subplot(5,1,3)
		# plt.title('Push, PER, CO2')
		# plt.hist(bsl_evtmean_01_list, len(bin_edges_bsl)-1, density=True, color='k', alpha=0.75)
		# if len(Push_evtmean_list)>min_evt_num:
		# 	hist_push, bin_edges_push = np.histogram(Push_evtmean_list, density=True)
		# 	plt.hist(Push_evtmean_list, len(bin_edges_push)-1, density=True, color=plot_setting.Push_color, alpha=0.75)
		# if len(PER_evtmean_list)>min_evt_num:
		# 	hist_per, bin_edges_per = np.histogram(PER_evtmean_list, density=True)
		# 	plt.hist(PER_evtmean_list, len(bin_edges_per)-1, density=True, color=plot_setting.PER_color, alpha=0.75)
		# if len(CO2puff_evtmean_list)>min_evt_num:
		# 	hist_co2puff, bin_edges_co2puff = np.histogram(CO2puff_evtmean_list, density=True)
		# 	plt.hist(CO2puff_evtmean_list, len(bin_edges_co2puff)-1, density=True, color=plot_setting.CO2puff_color, alpha=0.5)

		# plt.subplot(5,1,4)
		# plt.title('HLG, AbdG,')
		# plt.hist(bsl_evtmean_01_list, len(bin_edges_bsl)-1, density=True, color='k', alpha=0.75)
		# if len(HL_groom_evtmean_list)>min_evt_num:
		# 	hist_hlg, bin_edges_hlg = np.histogram(HL_groom_evtmean_list, density=True)
		# 	plt.hist(HL_groom_evtmean_list, len(bin_edges_hlg)-1, density=True, color=plot_setting.HL_groom_color, alpha=0.75)
		# if len(Abd_groom_evtmean_list)>min_evt_num:
		# 	hist_abdg, bin_edges_abdg = np.histogram(Abd_groom_evtmean_list, density=True)
		# 	plt.hist(Abd_groom_evtmean_list, len(bin_edges_abdg)-1, density=True, color=plot_setting.Abd_groom_color, alpha=0.75)

		# plt.subplot(5,1,5)
		# plt.title('Non-overlabp against baseline points')
		# # print('beh_name_list[:-1]', beh_name_list[:-1])
		# # print('all_behdatapoint_nonoverlap_w_Bsl_01', all_behdatapoint_nonoverlap_w_Bsl_01)
		# plt.bar(beh_name_list[:-1], all_behdatapoint_nonoverlap_w_Bsl_01)



		# plt.savefig(outDirEventsSumFly + 'hist_event_bsl_ROI_'+str(ROI_i)+'.png', edgecolor='none', transparent=False) #bbox_inches='tight', 
		# plt.clf()
		# plt.close()  

		
		dFF_list_ROI=[]
		dFF_norm_list_ROI=[]
		dFF_log_list_ROI=[]
		p_value_list_ROI=[]






reordered_ROI_ID_list, reordered_mean_dFF_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_id_list, mean_dFF_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_mean_dFF_01_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_id_list, mean_dFF_01_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_p_value_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_id_list, p_value_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_nonoverlap_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_id_list, nonoverlap_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_bslmean_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_id_list, bsl_mean_eachFly_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_bslmean_01_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_id_list, bsl_mean_01_eachFly_list, manual_ROI_order)

reordered_ROI_ID_list, reordered_normality_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_id_list, normality_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_std_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_id_list, std_list, manual_ROI_order)


reordered_masked_mean_dFF_list, bin_p_value_mask = stats_utils.mask_out_dFF_based_on_p_value(reordered_mean_dFF_list, ref_mat=reordered_p_value_list)
reordered_masked_mean_dFF_01_list, bin_p_value_mask = stats_utils.mask_out_dFF_based_on_p_value(reordered_mean_dFF_01_list, ref_mat=reordered_p_value_list)




y_list_allBeh=[
'F_Walk', 
'B_Walk', 
'Push', 
'Rest', 
'E_groom', 
'A_groom', 
'FL_rub', 
'Abd_groom', 
'HL_rub', 
'PE', 
'CO2puff'
]


group_dFF_summary_dir= NAS_AN_Proj_Dir+'_Summary/group_dFF_summary_20210619_df3d_and_ball_beh_BW-20foldNormTest-30resampBsl-resampBehdatapoint-6foldmultiCompar-alldata-ANOVA/' 
if not os.path.exists(group_dFF_summary_dir):
	os.makedirs(group_dFF_summary_dir)

save_data_for_dFF_matrix(filename='data_for_dFF_matrix_dic.p')
# save_data_for_dFF_matrix(filename='data_for_dFF_matrix_dic-BehEpochMean.p')


