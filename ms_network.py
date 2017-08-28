#THIS IS THE PART OF THE NEURAL NETWORK. IT TAKES THE VALUES TO CALCULATE THE POISSONRATES OF THE INPUT NEURONS AND THEN CALCULATES THE NETWORK OUTPUT 
#the unit in space is 1arcmin!
import pylab as pyl
import numpy as np
import cv2
import sys
import matplotlib.pyplot as plt
import nest
import nest.topology as tp
from microsaccades_functions import *

nest.ResetKernel()   # since we run the script multiple times

#necessary paramater definitions
#frames = 100 #replaced by motion detectors only
#synaptic weights for the 
syn_weight = .1

#----------------------------------------------------------------------------INPUT-RATES-FROM-MS_INPUT

extent = 121.
delay = 1.

sim_title = sys.argv[1]
handle_name = sys.argv[2]
#extent = float(sys.argv[3])
extent_x = float(sys.argv[3])
extent_y = float(sys.argv[4])
delay = float(sys.argv[5])
vel = float(sys.argv[6])
weight = float(sys.argv[6])
#extent_x = extent
#extent_y = extent
#center = extent/2.
center_x = extent_x/2.
center_y = extent_y/2.

#here needs to be a part that transfers potentials into poisson rates
#m_file = open('data/midget_values.data','r+')
m_file = open('data/'+sim_title+'/midget_rates_'+handle_name+'.data','r+')
m_data = np.load(m_file)   
m_file.close()

#p_file = open('data/parasolic_values.data','r+')
p_file = open('data/'+sim_title+'/parasolic_rates_'+handle_name+'.data','r+')
p_data = np.load(p_file)  
p_file.close()

midget_rates=poissonRateMidgets(m_data)
parasolic_rates=poissonRateParasols(p_data)

#print midget_rates

#to check for maximum spike rates in order to adopt conversion of film input
maxs = []
for i in range(len(midget_rates)):
    for j in range(len(midget_rates[0])):
       maxs.append(max(midget_rates[i][j]))
print max(maxs)
mmr = 400.#max(maxs)
        
paxs = []
for i in range(len(parasolic_rates)):
    for j in range(len(parasolic_rates[0])):
       paxs.append(max(parasolic_rates[i][j]))
print max(paxs)
pmr = max(paxs)

        
#for rates 
mrs = []
prs=[]

#print parasolic_rates


#-----------------------------------------------------------------------------------------NETWORK-PART
#---------------------------------------------------------------------------INITIALIZE-POISSON-NEURONS

#storage = open('data/'+sim_title+'/network/'+handle_name+'.data','w+')
#storage.close()

#nest.SetKernelStatus({'resolution': 0.001,
#                      'overwrite_files': True,
#                      'data_path': 'data/'+sim_title+'/network/',
#                      'data_prefix': ''})

#set the initial Poisson rate
rate = 100.

#create all newly needed models
#nest.CopyModel("poisson_generator", "var_poisson_generator",{"rate": rate})
#for the output
nest.CopyModel("spike_detector", "my_spike_detector",{"withgid": True, "withtime": True})
nest.CopyModel('multimeter', 'my_multimeter',{'withtime': True, 'interval': 0.1, 'withgid': True,'record_from': ['Rate']})

#nest.CopyModel("static_synapse", "ex", {"weight" : 0.1})
#nest.CopyModel("static_synapse", "inh", {"weight" : -0.1})

#----------------------------------------------------------------------------------------CREATE-LAYERS

#rows = len(midget_rates)
#cols = len(midget_rates[0])
#print rows, cols
#rows = 40
#cols = 40

#get grid data form previous simulation
gm_file = open('data/'+sim_title+'/m_pos_'+handle_name+'.data','r+')
#gm_file = open('data/test_vid/m_pos_test_vid.data','r+')
gm_data = np.load(gm_file)  
gm_file.close()
gm_data = gm_data.tolist()

gp_file = open('data/'+sim_title+'/p_pos_'+handle_name+'.data','r+')
#gp_file = open('data/test_vid/p_pos_test_vid.data','r+')
gp_data = np.load(gp_file)  
gp_file.close()
gp_data = gp_data.tolist()

gm_pos=[]
#positions for the reichardt detectors, which have to be centered
gm_r_0_pos=[]
gm_r_60_pos=[]
#gm_r_120_pos=[]

#gm_min=(0.,0.)
#gm_max=(0.,0.)
print 'length:'
print len(gm_data),len(gm_data[0])
for i in range(len(gm_data)):
    for j in range(len(gm_data[0])):
        mrs+=[1000.*midget_rates[i][j]]
        gm_pos+=[[0.5*gm_data[i][j][0]+0.25,0.5*gm_data[i][j][1]+0.25]]
        gm_r_0_pos+=[[0.5*gm_data[i][j][0]+0.5,0.5*gm_data[i][j][1]+0.25]]
        gm_r_60_pos+=[[0.5*gm_data[i][j][0]+0.5,0.5*gm_data[i][j][1]+0.25+0.433]]
        #gm_r_120_pos+=[[gm_data[i][j][0],gm_data[i][j][1]]+0.25+0.433]
        '''
        #currently, this is known in the simulations
        if gm_data[i][j][0]>gm_max[0]:
            gm_max=(gm_data[i][j][0],gm_min[1])
        if gm_data[i][j][1]>gm_max[1]:
            gm_max=(gm_max[0],gm_data[i][j][1])
        if gm_data[i][j][0]<gm_min[0]:
            gm_min=(gm_data[i][j][0],gm_min[1])
        if gm_data[i][j][1]<gm_min[1]:
            gm_min=(gm_min[0],gm_data[i][j][1])
            #print gm_data[i][j]
        '''
'''
#print gm_pos
print 'minimum midgets'
print gm_min, gm_max[0]+gm_min[0]+10.
print 'maximum midgets'
print gm_max, gm_max[1]+gm_min[1]+10.
'''       
gp_pos=[]
gp_r_0_pos=[]
gp_r_60_pos=[]
#gm_r_120_pos=[]

for i in range(len(gp_data)):
    for j in range(len(gp_data[0])):
        prs+=[400.*parasolic_rates[i][j]]
        gp_pos+=[[0.5*gp_data[i][j][0]+0.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_0_pos+=[[0.5*gp_data[i][j][0]+2.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_60_pos+=[[0.5*gp_data[i][j][0]+1.25,0.5*gp_data[i][j][1]+0.25+1.732]]
        #gp_r_120_pos+=[[gp_data[i][j][0]-1.25,gp_data[i][j][1]+0.25+1.732]]


#gm_pos=[[np.random.uniform(-0.5,0.5), np.random.uniform(-0.5,0.5)] for j in range(50)]
#print gm_pos
nest.SetDefaults('iaf_psc_alpha',{'I_e' : 374.0})

midgets = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'poisson_generator', 'edge_wrap': True})
parasolic = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'poisson_generator', 'edge_wrap': True})

#120 degree commented for the moment since two build a basis

#these are direction dependent
m_reichardt_0_left = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha'})
m_reichardt_0_right = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
m_reichardt_60_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
m_reichardt_60_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#m_reichardt_120_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#m_reichardt_120_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_0_left = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_0_right = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_60_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_60_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#p_reichardt_120_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#p_reichardt_120_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})

#react to motion in both directions
m_reichardt_0 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_0_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
m_reichardt_60 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_r_60_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#m_reichardt_120 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_0 = tp.CreateLayer({'extent' : [125.,125.], 'center' : [center_x,center_y], 'positions' : gp_r_0_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_60 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_r_60_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#p_reichardt_120 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})

'''
motion_left = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
motion_right = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
motion_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
motion_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
'''
#think about distances etc
#tp.PlotLayer(midgets)
#tp.PlotLayer(parasolic)
#tp.PlotLayer(reichardt_left)
#plt.show()

out_m = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'my_spike_detector'})
out_p = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector'})
out_m_multi = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'my_multimeter'})
out_p_multi = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_multimeter'})

out_m_r_0 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_0_pos, 'elements': 'my_spike_detector'})
out_m_r_0_multi = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_0_pos, 'elements': 'my_multimeter'})
out_m_r_60 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_60_pos, 'elements': 'my_spike_detector'})
out_m_r_60_multi = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_60_pos, 'elements': 'my_multimeter'})

out_p_r_0 = tp.CreateLayer({'extent' : [125.,125.], 'center' : [center_x,center_y], 'positions' : gp_r_0_pos, 'elements': 'my_spike_detector'})
out_p_r_0_multi = tp.CreateLayer({'extent' : [125.,125.], 'center' : [center_x,center_y], 'positions' : gp_r_0_pos, 'elements': 'my_multimeter'})
out_p_r_60 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_r_60_pos, 'elements': 'my_spike_detector'})
out_p_r_60_multi = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_r_60_pos, 'elements': 'my_multimeter'})
#out = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'my_spike_detector'})
#out_multi = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'my_multimeter'})



nest.CopyModel('static_synapse', 'exc', {'weight' : 2.0})
#-----------------------------------------------------------------------------------CREATE-CONNECTIONS
#connections to left/right half of Reichardt detector
m_r_0_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.6,-0.1], 'upper_right' : [0.1,0.1]}}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
m_r_0_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.1,0.6]}}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
m_r_60_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,-0.5], 'upper_right' : [0.5,0.5]}, 'anchor' : [0.25,0.433]}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
m_r_60_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,-0.5], 'upper_right' : [0.5,0.5]}, 'anchor' : [-0.25,-0.433]}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
p_r_0_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.1,-0.1], 'upper_right' : [0.1,0.1]}}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
p_r_0_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.1,2.1]}}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
p_r_60_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.,-2.], 'upper_right' : [2.,2.]}, 'anchor' : [1.,1.732]}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
p_r_60_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.,-2.], 'upper_right' : [2.,2.]}, 'anchor' : [-1.,-1.732]}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}

#connections of left/right half to just motion sensitive Reichardt detector
m_r_0_left_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.3,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
m_r_0_right_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.3,-0.1], 'upper_right' : [0.1,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
m_r_60_up_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.3,-0.5], 'upper_right' : [0.1,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
m_r_60_down_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.3,0.5]}}, 'weights' : weight } #'synapse_model' : 'ex'}
p_r_0_left_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [1.1,20.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
p_r_0_right_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-1.1,-0.1], 'upper_right' : [0.1,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
p_r_60_up_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-1.25,-1.], 'upper_right' : [0.1,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
p_r_60_down_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [1.1,1.]}}, 'weights' : weight } #'synapse_model' : 'ex'}

'''
#left commented for the moment
m_r_hor_motion_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [4.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'ex'}
m_r_hor_motion_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-4.,0.], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'ex'}
m_r_ver_motion_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,-4.], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'ex'}
m_r_ver_motion_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [0.1,4.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'ex'}
'''
#what about delays here?
#par_motion_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [0.1,2.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'inh'}
#par_motion_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,-2.], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'inh'}
#par_motion_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.,0.], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'inh'}
#par_motion_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [2.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'inh'}

out_conndict = {'connection_type' : 'convergent', 'mask' : {'rectangular' : {'lower_left' : [-0.2,-0.2], 'upper_right' : [0.2,0.2]}}}

#connect them

tp.ConnectLayers(midgets,m_reichardt_0_left,m_r_0_left_conndict)
tp.ConnectLayers(midgets,m_reichardt_0_right,m_r_0_right_conndict)
tp.ConnectLayers(midgets,m_reichardt_60_up,m_r_60_up_conndict)
tp.ConnectLayers(midgets,m_reichardt_60_down,m_r_60_down_conndict)

tp.ConnectLayers(parasolic,p_reichardt_0_left,p_r_0_left_conndict)
tp.ConnectLayers(parasolic,p_reichardt_0_right,p_r_0_right_conndict)
tp.ConnectLayers(parasolic,p_reichardt_60_up,p_r_60_up_conndict)
tp.ConnectLayers(parasolic,p_reichardt_60_down,p_r_60_down_conndict)

tp.ConnectLayers(m_reichardt_0_left,m_reichardt_0,m_r_0_left_r_0_conndict)
tp.ConnectLayers(m_reichardt_0_right,m_reichardt_0,m_r_0_right_r_0_conndict)
tp.ConnectLayers(m_reichardt_60_up,m_reichardt_60,m_r_60_up_r_60_conndict)
tp.ConnectLayers(m_reichardt_60_down,m_reichardt_60,m_r_60_down_r_60_conndict)

tp.ConnectLayers(p_reichardt_0_left,p_reichardt_0,p_r_0_left_r_0_conndict)
tp.ConnectLayers(p_reichardt_0_right,p_reichardt_0,p_r_0_right_r_0_conndict)
tp.ConnectLayers(p_reichardt_60_up,p_reichardt_60,p_r_60_up_r_60_conndict)
tp.ConnectLayers(p_reichardt_60_down,p_reichardt_60,p_r_60_down_r_60_conndict)

tp.ConnectLayers(midgets,out_m,out_conndict)
tp.ConnectLayers(parasolic,out_p,out_conndict)
tp.ConnectLayers(midgets,out_m_multi,out_conndict)
tp.ConnectLayers(parasolic,out_p_multi,out_conndict)


tp.ConnectLayers(m_reichardt_0,out_m_r_0,out_conndict)
tp.ConnectLayers(p_reichardt_0,out_p_r_0,out_conndict)
#tp.ConnectLayers(out_m_r_0_multi,m_reichardt_0,out_conndict)
#tp.ConnectLayers(p_reichardt_0,out_p_r_0_multi,out_conndict)

#tp.ConnectLayers(m_reichardt_60,out_m_r_60,out_conndict)
#tp.ConnectLayers(p_reichardt_60,out_p_r_60,out_conndict)
#tp.ConnectLayers(m_reichardt_60,out_m_r_60_multi,out_conndict)
#tp.ConnectLayers(p_reichardt_60,out_p_r_60_multi,out_conndict)

#ctr = tp.FindNearestElement(m_reichardt_0,[6.,12.])
#fig = tp.PlotLayer(midgets,nodesize =80)
#tp.PlotTargets(ctr,midgets,fig=fig,mask=m_r_0_right_conndict['mask'],src_size=250, tgt_color='red',tgt_size=20)

#-------------------------------------------------------------------------------------------SIMULATION
MIDs = nest.GetNodes(midgets)
PIDs = nest.GetNodes(parasolic)

for f in range(2):#00,300):#frames):
    print f
    #reset rates
    for n in range(len(MIDs[0])):
        qr = mrs[n][f]
        #spontaneous firing rate
        #if qr < 5.:
        #    qr = 5.
        nest.SetStatus([MIDs[0][n]], {'rate': qr})
    for n in range(len(PIDs[0])):
        qr = prs[n][f]
        #spontaneous firing rate
        if qr < 5.:
            qr = 5.
        nest.SetStatus([PIDs[0][n]], {'rate': qr})
    #run simulation
    nest.Simulate(1)


#-----------------------------------------------------------------------SAVING-AND-PRINTING-THE-OUTPUT

store_sp = []
s_gids = []

ms_all_file = open('data/mo_det_cal/m_max_spikes_all.txt','a+')

for q in range(120):
    s = tp.FindNearestElement(out_m_r_0,[float(q),12.])
    #mult = tp.FindNearestElement(out_m_r_0_multi,[float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., 101., 100.))
        store_sp +=[h]
        ms_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(sp)+'\n')
        ms_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(h)+'\n')
        ms_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(e)+'\n')
        #plt.plot(t, r, color='b')
        #plt.step(e[:-1], h , color='b', where='post')
        #plt.title('PST histogram and firing rates')
        #plt.ylabel('Spikes per second')

        #plt.subplot(223)
        #plt.hist(np.diff(sp), bins=np.arange(0., 1.005, 0.02),
                    #histtype='step', color='b')
        #plt.title('ISI histogram')
        #plt.show()
    else:
        store_sp +=[0]

ms_all_file.close()
#plt.show()
max_spikes = max(store_sp)
if max_spikes == 0:
    max_spikes = [0]
print max_spikes

k = (delay,vel,max_spikes)
ms_file = open('data/mo_det_cal/m_max_spikes.txt','a+')
ms_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+ str(max_spikes[0])+'\n')
ms_file.close()

ps_all_file = open('data/mo_det_cal/p_max_spikes_all.txt','a+')

store_sp=[]

for q in range(30):
    s = tp.FindNearestElement(out_p_r_0,[4.*float(q),12.])
    #mult = tp.FindNearestElement(out_p_r_0_multi,[4.*float(q),12.])
    s_gids += [s]
    
    dSD = nest.GetStatus(s,keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]
    #pyl.figure(2)
    #pyl.plot(ts, evs, ".")
    #pyl.show()
    
    if evs.any():
        #ev = nest.GetStatus(mult)[0]['events']
        #t = ev['times']
        #r = ev['rate']

        sp = nest.GetStatus(s)[0]['events']['times']
        plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., 101., 100.))
        store_sp +=[h]
        ps_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(sp)+'\n')
        ps_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(h)+'\n')
        ps_all_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+str(e)+'\n')
        #plt.plot(t, r, color='b')
        #plt.step(e[:-1], h , color='b', where='post')
        #plt.title('PST histogram and firing rates')
        #plt.ylabel('Spikes per second')

        #plt.subplot(223)
        #plt.hist(np.diff(sp), bins=np.arange(0., 1.005, 0.02),
                    #histtype='step', color='b')
        #plt.title('ISI histogram')
        #plt.show()
    else:
        store_sp +=[0]

ps_all_file.close()
#plt.show()
max_spikes = max(store_sp)
if max_spikes == 0:
    max_spikes = [0]
print max_spikes

k = (delay,vel,max_spikes)
ps_file = open('data/mo_det_cal/p_max_spikes.txt','a+')
ps_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+ str(max_spikes[0])+'\n')
ps_file.close()