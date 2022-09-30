
# baseColor=[[255,255,200]]
# midColor=[[255,0,0]]
# topcolor=[[0,0,0]]

# baseColor=[[230,255,230]]
# midColor=[[0,255,0]]
# topColor=[[0,50,0]]

# #blue-red color
# baseColor=[[0,0,250]]
# midColor=[[125,0,125]] 
# topColor=[[250,0,0]]

baseColor=[[176,224,230]]
midColor=[[88,112,185]] # average between top_+color and base_color
topColor=[[0,0,139]]

restFrame_color=[[0,153,0]]

tsne_dot_size=10
tsne_circle_size=14
rest_circle_thickness=1
current_circle_thickness=6



timedot_color='none'
timedot_edgecolor='orange' #(255,128,0)
timedot_facecolor='orange'
timedot_edgewidth=0
timedot_size=4

axis_width=2
axis_color='k'
axis_color_blackBack='w'
tick_length=5
tick_width=2

CO2puff_color='k'
GC_color='forestgreen'
# GC_color_no_cur_neuron='darkgreen'
AP_color='r'
ML_color='blue'
Yaw_color='m'



walk_color='skyblue'
rest_color='silver'
groom_color='orchid'

FW_color='blue'
BW_color='springgreen'

E_groom_color='hotpink'
A_groom_color='darkorange'
FL_groom_color='darkviolet'

HL_groom_color='olive'
Abd_groom_color='yellowgreen'

PER_color='gold'

Push_color='skyblue'

F_groom_color='orchid'
H_groom_color='mediumaquamarine'





EthoColorDic={}
EthoColorDic.update({'rest_evt':rest_color})
EthoColorDic.update({'FW_evt':FW_color})
EthoColorDic.update({'BW_evt':BW_color})
EthoColorDic.update({'E_groom_evt':E_groom_color})
EthoColorDic.update({'A_groom_evt':A_groom_color})
EthoColorDic.update({'FL_groom_evt':FL_groom_color})
EthoColorDic.update({'HL_groom_evt':HL_groom_color})
EthoColorDic.update({'Abd_groom_evt':Abd_groom_color})
EthoColorDic.update({'PER_evt':PER_color})
EthoColorDic.update({'Push_evt':Push_color})

EthoColorDic.update({'F_groom_evt':F_groom_color})
EthoColorDic.update({'H_groom_evt':H_groom_color})	
EthoColorDic.update({'groom_evt':groom_color})

EthoColorDic.update({'move_evt':FW_color})
EthoColorDic.update({'walk_evt':FW_color})






data_trace_width=1.2
data_font_size=28
data_ylabel_position_for_ball=[-0.21,0.5]
data_ylabel_position=[-0.20,0.5]
err_shade_alpha=0.3

font_color='k'
font_color_blackBack='w'

axvline_color='gray'
axvline_color_blackBack='gray'
axvline_width=0.5
axhline_color='gray'
axhline_color_blackBack='gray'
axhline_width=0.5


curFrame_color=[[255,69,0]]
restFrame_color = [[0,153,0]]

jointangle_color_each = 'dodgerblue'




Beh_labels=[
'Walk',
'Rest',
'Groom',
'CO2 puff',
]

Ball_rota_labels=[
'velForw',
'velSide',
'velTurn',
]



leg_colors_Victor={
'RH_leg':(102/255,178/255,1),
'RM_leg':(0,0.5,1),
'RF_leg':(0,76/255,153/255),
'LH_leg':(1,102/255,102/255),
'LM_leg':(1,51/255,51/255),
'LF_leg':(204/255,0,0)
}


leg_colors={
'RH_leg':(255/255,204/255,204/255),
'RM_leg':(255/255,128/255,128/255),
'RF_leg':(255/255,0/255,0/255),
'LH_leg':(204/255,204/255,255/255),
'LM_leg':(128/255,128/255,255/255),
'LF_leg':(0/255,0/255,255/255)
}


# leg_colors_test={
# 'RH_leg':(255/255,255/255,255/255),
# 'RM_leg':(255/255,255/255,255/255),
# 'RF_leg':(255/255,255/255,255/255),
# 'LH_leg':(204/255,204/255,255/255),
# 'LM_leg':(128/255,128/255,255/255),
# 'LF_leg':(0/255,0/255,255/255)
# }










