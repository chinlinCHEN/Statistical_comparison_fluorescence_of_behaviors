import sys
import pickle
import pandas as pd
import numpy as np



import utils.general_utils as general_utils
import utils.plot_utils as plot_utils
import utils.plot_setting as plot_setting
import utils.math_utils as math_utils
import utils.sync_utils as sync_utils
import utils.list_inputFiles as list_inputFiles
import utils.stats_utils as stats_utils






NAS_AN_Proj_Dir=general_utils.NAS_AN_Proj_Dir
NAS_Dir=general_utils.NAS_Dir
workstation_dir=general_utils.workstation_dir

group_dFF_summary_dir= NAS_AN_Proj_Dir+'_Summary/group_dFF_summary_20210619_df3d_and_ball_beh_BW-20foldNormTest-30resampBsl-resampBehdatapoint-6foldmultiCompar-alldata-ANOVA/' 


dFF_for_matPlot_DicData=general_utils.open_Beh_Jpos_GC_DicData(group_dFF_summary_dir, 'data_for_dFF_matrix_dic.p')


print(dFF_for_matPlot_DicData.keys())
print(dFF_for_matPlot_DicData['reordered_mean_dFF_list'])


manual_ROI_order_csv=pd.read_csv(workstation_dir+'from_florian/Ascending_analysis/output/row_order_manual.csv')
manual_ROI_order= manual_ROI_order_csv['x'].tolist()

ROI_ID_list=dFF_for_matPlot_DicData['reordered_ROI_ID_list']
y_list_allBeh=dFF_for_matPlot_DicData['y_list_allBeh']

mean_dFF_list=dFF_for_matPlot_DicData['reordered_mean_dFF_list']
mean_dFF_01_list=dFF_for_matPlot_DicData['reordered_mean_dFF_01_list']
masked_mean_dFF_list=dFF_for_matPlot_DicData['reordered_masked_mean_dFF_list']
masked_mean_dFF_01_list=dFF_for_matPlot_DicData['reordered_masked_mean_dFF_01_list']
p_value_list=dFF_for_matPlot_DicData['reordered_p_value_list']
bin_p_value_mask=dFF_for_matPlot_DicData['bin_p_value_mask']
nonoverlap_list=dFF_for_matPlot_DicData['reordered_nonoverlap_list']
normality_list=dFF_for_matPlot_DicData['reordered_normality_list']


combined_mask=stats_utils.make_mask_based_on_p_value_and_normality(p_value_list, normality_list)

print('shape bin_p_value_mask', np.shape(bin_p_value_mask))
print('shape combined_mask', np.shape(combined_mask))







reordered_ROI_ID_list, reordered_mean_dFF_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_ID_list, mean_dFF_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_mean_dFF_01_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_ID_list, mean_dFF_01_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_masked_mean_dFF_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_ID_list, masked_mean_dFF_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_masked_mean_dFF_01_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_ID_list, masked_mean_dFF_01_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_bin_p_value_mask=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_ID_list, bin_p_value_mask, manual_ROI_order)
reordered_ROI_ID_list, reordered_nonoverlap_list=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_ID_list, nonoverlap_list, manual_ROI_order)
reordered_ROI_ID_list, reordered_combined_mask=plot_utils.sorting_roiID_correspondingMat_based_on_an_order(ROI_ID_list, combined_mask, manual_ROI_order)


# plot_utils.plot_matrix(reordered_ROI_ID_list, y_list_allBeh, reordered_mean_dFF_list, savedir=group_dFF_summary_dir, title='mean_dFF', PlotMethod='other', cmap='Greens')
# plot_utils.plot_matrix(reordered_ROI_ID_list, y_list_allBeh, reordered_mean_dFF_01_list, savedir=group_dFF_summary_dir, title='mean_dFF_01', PlotMethod='dFF_01_w_negative')
# plot_utils.plot_matrix(reordered_ROI_ID_list, y_list_allBeh, reordered_masked_mean_dFF_list, savedir=group_dFF_summary_dir, title='mean_dFF_masked_by_p', PlotMethod='dFF_01_w_negative')
# plot_utils.plot_matrix(reordered_ROI_ID_list, y_list_allBeh, reordered_masked_mean_dFF_01_list, savedir=group_dFF_summary_dir, title='mean_dFF_01_masked_by_p', PlotMethod='dFF_01_w_negative')
# plot_utils.plot_matrix(reordered_ROI_ID_list, y_list_allBeh, reordered_p_value_list, savedir=group_dFF_summary_dir, title='p_value', PlotMethod='p_value', cmap='Greens')
# plot_utils.plot_matrix(reordered_ROI_ID_list, y_list_allBeh, bin_p_value_mask, savedir=group_dFF_summary_dir, title='p_value_bin_mask', PlotMethod='other', cmap='gray', hatch='/', hatchcolor='lightgrey')
# plot_utils.plot_matrix(reordered_ROI_ID_list, y_list_allBeh, reordered_nonoverlap_list, savedir=group_dFF_summary_dir, title='Non-overlap_against_bsl', PlotMethod='other', cmap='Greens')

plot_utils.plot_matrix(reordered_ROI_ID_list, y_list_allBeh, reordered_masked_mean_dFF_01_list, savedir=group_dFF_summary_dir, title='mean_dFF_01_masked_by_p', PlotMethod='dFF_01_w_negative', Gal4_x_list_reformat=True, cmap='BuPu')


plot_utils.plot_overlay_matrix(reordered_ROI_ID_list, y_list_allBeh, reordered_masked_mean_dFF_01_list, reordered_combined_mask, hatch=False, colorbar_bas='BuPu', savedir=group_dFF_summary_dir, title='mean_dFF_01_masked_by_p_w_mask')
