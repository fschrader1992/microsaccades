#THIS IS THE PART OF THE NEURAL NETWORK. IT TAKES THE VALUES TO CALCULATE THE POISSONRATES OF THE INPUT NEURONS AND THEN CALCULATES THE NETWORK OUTPUT 
#the unit in space is 1arcmin!
import pylab as pyl
import numpy as np
import cv2
import sys
import matplotlib.pyplot as plt
import nest
import datetime
import nest.raster_plot
import nest.topology as tp
from microsaccades_functions import *

nest.ResetKernel()   # since we run the script multiple times
#to get different poisson outputs
'''
msd = 1000  #master seed
n_vp = nest.GetKernelStatus('total_num_virtual_procs')
msdrange1 = range(msd, msd+n_vp)
pyrngs = [np.random.RandomState(s) for s in msdrange1]
msdrange2 = range(msd+n_vp+1, msd+1+2*n_vp)
nest.SetKernelStatus({'grng_seed': msd+n_vp,
                      'rng_seeds': msdrange2})
'''
#necessary paramater definitions
#frames = 100 #replaced by motion detectors only

#----------------------------------------------------------------------------INPUT-RATES-FROM-MS_INPUT

extent = 121.
delay = 34. # speed of point in poletti 2010 -> maybe increase a bit for more reacton later on


sim_title = sys.argv[1]
handle_name = sys.argv[2]
#extent = float(sys.argv[3])
extent_x = float(sys.argv[3])
extent_y = float(sys.argv[4])
#delay = float(sys.argv[5])
vel = float(sys.argv[6])
weight = float(sys.argv[7])
#extent_x = extent
#extent_y = extent
#center = extent/2.
center_x = extent_x/2.
center_y = extent_y/2.

print weight,delay
#weight = 0.

#here needs to be a part that transfers potentials into poisson rates
#m_file = open('data/midget_values.data','r+')
m_file = open('data/'+sim_title+'/midget_rates_'+handle_name+'_on.data','r+')
m_data = np.load(m_file)   
m_file.close()

#p_file = open('data/parasolic_values.data','r+')
p_file = open('data/'+sim_title+'/parasolic_rates_'+handle_name+'_on.data','r+')
p_data = np.load(p_file)  
p_file.close()

midget_rates=poissonRateMidgets(m_data)
parasolic_rates=poissonRateParasols(p_data)

#print midget_rates


#to check for maximum spike rates in order to adopt conversion of film input
maxs = []
for i in range(len(midget_rates)):
    for j in range(len(midget_rates[0])):
        maxs.append(max(midget_rates[i][j][200:400]))
mmr = 25000.
#print maxs
print 'midget max: ' + str(max(maxs)) + ' ' + str(mmr)
        
paxs = []
for i in range(len(parasolic_rates)):
    for j in range(len(parasolic_rates[0])):
       paxs.append(max(parasolic_rates[i][j][200:400]))
pmr = 40000.
print 'parsolic max: ' + str(max(paxs)) + ' ' + str(pmr)

ps_file = open('data/mo_det_cal/input_max_spikes.txt','a+')
ps_file.write(str(handle_name)+' m: '+ str(max(maxs))+' p: '+ str(max(paxs))+'\n')
ps_file.close()
        
#for rates 
mrs = []
prs=[]

#print midget_rates
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
#iaf_psc_alpha: ({u'vp': 0, u'synaptic_elements': {}, u'V_reset': -70.0, u'supports_precise_spikes': False, u'V_th': -55.0, u'tau_minus': 20.0, u'tau_m': 10.0, u'I_e': 0.0, u'element_type': <SLILiteral: neuron>, u't_spike': -1.0, u'V_m': -70.0, u'local': True, u'recordables': (<SLILiteral: I_syn_ex>, <SLILiteral: I_syn_in>, <SLILiteral: V_m>, <SLILiteral: weighted_spikes_ex>, <SLILiteral: weighted_spikes_in>), u'E_L': -70.0, u'node_uses_wfr': False, u'tau_syn_ex': 2.0, u'parent': 2913, u'tau_minus_triplet': 110.0, u'thread_local_id': 5582, u'beta_Ca': 0.001, u't_ref': 2.0, u'thread': 0, u'frozen': False, u'archiver_length': 0, u'Ca': 0.0, u'C_m': 250.0, u'global_id': 5583, u'local_id': 2670, u'tau_syn_in': 2.0, u'model': <SLILiteral: iaf_psc_alpha>, u'V_min': -inf, u'tau_Ca': 10000.0},)


nest.CopyModel("iaf_psc_alpha", "iaf_psc_alpha_i",{"I_e" : 340.0})
nest.CopyModel("iaf_psc_alpha", "iaf_psc_alpha_r",{"I_e" : 365.0})
nest.CopyModel("spike_detector", "my_spike_detector",{"withgid": True, "withtime": True})
nest.CopyModel('multimeter', 'my_multimeter',{'withtime': True, 'interval': 0.1, 'withgid': True,'record_from': ['Rate']})

#nest.SetDefaults('iaf_psc_alpha_r',{'I_e' : 374.0})
#nest.SetDefaults('iaf_psc_alpha',{'I_e' : 350.0})
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
        mrs+=[mmr*midget_rates[i][j]]
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
        prs+=[pmr*parasolic_rates[i][j]]
        gp_pos+=[[0.5*gp_data[i][j][0]+0.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_0_pos+=[[0.5*gp_data[i][j][0]+2.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_60_pos+=[[0.5*gp_data[i][j][0]+1.25,0.5*gp_data[i][j][1]+0.25+1.732]]
        #gp_r_120_pos+=[[gp_data[i][j][0]-1.25,gp_data[i][j][1]+0.25+1.732]]


#gm_pos=[[np.random.uniform(-0.5,0.5), np.random.uniform(-0.5,0.5)] for j in range(50)]
#print gm_pos

midgets = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'poisson_generator', 'edge_wrap': True})
parasolic = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'poisson_generator', 'edge_wrap': True})

#120 degree commented for the moment since two build a basis

#these are direction dependent
m_reichardt_0_left = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
m_reichardt_0_right = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
m_reichardt_60_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
m_reichardt_60_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
#m_reichardt_120_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
#m_reichardt_120_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
p_reichardt_0_left = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
p_reichardt_0_right = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
p_reichardt_60_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
p_reichardt_60_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
#p_reichardt_120_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
#p_reichardt_120_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})

#react to motion in both directions
m_reichardt_0 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_0_pos, 'elements': 'iaf_psc_alpha_r', 'edge_wrap': True})
m_reichardt_60 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_r_60_pos, 'elements': 'iaf_psc_alpha_r', 'edge_wrap': True})
#m_reichardt_120 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_0 = tp.CreateLayer({'extent' : [125.,125.], 'center' : [center_x,center_y], 'positions' : gp_r_0_pos, 'elements': 'iaf_psc_alpha_r', 'edge_wrap': True})
p_reichardt_60 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_r_60_pos, 'elements': 'iaf_psc_alpha_r', 'edge_wrap': True})
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

out_m_r_0_left = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'my_spike_detector'})
out_m_r_0_right = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'my_spike_detector'})

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

out_p_r_0_left = tp.CreateLayer({'extent' : [125.,125.], 'center' : [center_x,center_y], 'positions' : gp_r_0_pos, 'elements': 'my_spike_detector'})
out_p_r_0_right = tp.CreateLayer({'extent' : [125.,125.], 'center' : [center_x,center_y], 'positions' : gp_r_0_pos, 'elements': 'my_spike_detector'})


#-----------------------------------------------------------------------------------CREATE-CONNECTIONS
#connections to left/right half of Reichardt detector
m_r_0_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.1,0.6]}}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
m_r_0_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.6,-0.1], 'upper_right' : [0.1,0.1]}}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
m_r_60_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,-0.5], 'upper_right' : [0.5,0.5]}, 'anchor' : [0.25,0.433]}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
m_r_60_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,-0.5], 'upper_right' : [0.5,0.5]}, 'anchor' : [-0.25,-0.433]}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
p_r_0_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.1,-0.1], 'upper_right' : [0.1,0.1]}}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : 4.*delay}}}
p_r_0_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.1,2.1]}}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : 4.*delay}}}
p_r_60_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.,-2.], 'upper_right' : [2.,2.]}, 'anchor' : [1.,1.732]}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : 4.*delay}}}
p_r_60_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.,-2.], 'upper_right' : [2.,2.]}, 'anchor' : [-1.,-1.732]}, 'weights' : weight, 'delays' : {'linear' : {'c' : .1, 'a' : 4.*delay}}}

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

tp.ConnectLayers(midgets,m_reichardt_0_right,m_r_0_right_conndict)
tp.ConnectLayers(midgets,m_reichardt_0_left,m_r_0_left_conndict)
tp.ConnectLayers(midgets,m_reichardt_60_up,m_r_60_up_conndict)
tp.ConnectLayers(midgets,m_reichardt_60_down,m_r_60_down_conndict)

tp.ConnectLayers(parasolic,p_reichardt_0_right,p_r_0_right_conndict)
tp.ConnectLayers(parasolic,p_reichardt_0_left,p_r_0_left_conndict)
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

tp.ConnectLayers(m_reichardt_0_left,out_m_r_0_left,out_conndict)
tp.ConnectLayers(m_reichardt_0_right,out_m_r_0_right,out_conndict)
tp.ConnectLayers(p_reichardt_0_left,out_p_r_0_left,out_conndict)
tp.ConnectLayers(p_reichardt_0_right,out_p_r_0_right,out_conndict)

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


t_start = 200
t_end = 400
for f in range(t_start,t_end):#frames):
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
now = datetime.datetime.now()

store_sp = []
s_gids = []

fig = plt.figure()
ax = plt.subplot(111)

firstq = True
#poisson midget cells
for q in range(1,121):
    s = tp.FindNearestElement(out_m,[float(q/2.),7.5])
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
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        #print q,s,sp
        qa = [q for i in range(h)]
        if firstq == True:
            ax.plot(qa,sp,'r.',label='poisson generators')
            firstq = False
        else:
            ax.plot(qa,sp,'r.',label='_nolegend_')
        
    else:
        store_sp +=[0]

m_poisson_max_spikes = max(store_sp)
if m_poisson_max_spikes == 0:
    m_poisson_max_spikes = [0]
print m_poisson_max_spikes


#pyl.figure()
#V_E = nest.GetStatus(s_gids, 'rate')
#pyl.hist(V_E, bins = 100)
#pyl.show()

#midget motion detectors complete
ms_all_file = open('data/mo_det_cal/m_max_spikes_all.txt','a+')
store_sp = []
s_gids = []

#midget motion detectors left
store_sp = []
s_gids = []

firstq = True
for q in range(1,121):
    s = tp.FindNearestElement(out_m_r_0_right,[float(q/2.),7.5])
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
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        qa = [q for i in range(h)]
        if firstq == True:
            ax.plot(qa,sp,'b+',label='leftward-motion detectors')
            firstq = False
        else:
            ax.plot(qa,sp,'b+',label='_nolegend_')

    else:
        store_sp +=[0]

max_spikes_right = max(store_sp)
if max_spikes_right == 0:
    max_spikes_right = [0]
print max_spikes_right

#midget motion detectors left
store_sp = []
s_gids = []

firstq = True
for q in range(1,121):
    s = tp.FindNearestElement(out_m_r_0_left,[float(q/2.),7.5])
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
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        qa = [q for i in range(h)]
        if firstq == True:
            ax.plot(qa,sp,'kx',label='rightward-motion detectors')
            firstq = False
        else:
            ax.plot(qa,sp,'kx',label='_nolegend_')

    else:
        store_sp +=[0]

max_spikes_left = max(store_sp)
if max_spikes_left == 0:
    max_spikes_left = [0]
print max_spikes_left


#midget motion detectors complete
ms_all_file = open('data/mo_det_cal/m_max_spikes_all.txt','a+')
store_sp = []
s_gids = []

firstq = True
for q in range(1,121):
    s = tp.FindNearestElement(out_m_r_0,[float(q/2.),7.5])
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
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        
        qa = [q for i in range(h)]
        if firstq == True:
            ax.plot(qa,sp,'g*',label='general motion detectors')
            firstq = False
        else:
            ax.plot(qa,sp,'g*',label='_nolegend_')
        
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

#plt.title('midget spikes for motion detectors')
#plt.xlabel('position of neuron in arcmin')
#plt.ylabel('time $t$ in ms')
#plt.legend()
plt.xlim([0,120])
plt.ylim([0,200])

ax.set_title('midget spikes for motion detectors')
ax.set_xlabel('position of neuron in arcmin')
ax.set_ylabel('time $t$ in ms')
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])


# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), 
          fancybox=False, fontsize=8, shadow=False, ncol=2)

plt.savefig('img/mo_det_cal/x_t_midgets_'+handle_name+'.pdf')
out= '/home/schrader/Documents/microsaccades/img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/x_t_midgets_' + str(handle_name) + '_' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '.pdf'
plt.savefig(out)
plt.show()
plt.close()



ms_all_file.close()
#plt.show()
max_spikes = max(store_sp)
if max_spikes == 0:
    max_spikes = [0]
print max_spikes



ms_file = open('data/mo_det_cal/m_max_spikes.txt','a+')
ms_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+ str(max_spikes[0])+' '+ str(m_poisson_max_spikes[0])+' '+ str(max_spikes_left[0])+' '+ str(max_spikes_right[0])+'\n')
ms_file.close()




#parasol poisson cells-------------------------------------------------------------------------------------------------------------

store_sp=[]

fig = plt.figure()
ax = plt.subplot(111)

firstq = True
for q in range(1,31):
    s = tp.FindNearestElement(out_p,[4.*float(q/2.),7.5])
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
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        qa = [4.*q for i in range(h)]
        if firstq == True:
            ax.plot(qa,sp,'r.',label='poisson generators')
            firstq = False
        else:
            ax.plot(qa,sp,'r.',label='_nolegend_')
    else:
        store_sp +=[0]


#plt.show()
p_poisson_max_spikes = max(store_sp)
if p_poisson_max_spikes == 0:
    p_poisson_max_spikes = [0]
print p_poisson_max_spikes

#midget motion detectors right
store_sp=[]

firstq = True
for q in range(1,31):
    s = tp.FindNearestElement(out_p_r_0_right,[4.*float(q/2.),7.5])
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
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        qa = [4*q for i in range(h)]
        if firstq == True:
            ax.plot(qa,sp,'b+',label='leftward-motion detectors')
            firstq = False
        else:
            ax.plot(qa,sp,'b+',label='_nolegend_')
    else:
        store_sp +=[0]


max_spikes_right = max(store_sp)
if max_spikes_right == 0:
    max_spikes_right = [0]
print max_spikes_right

#midget motion detectors left
store_sp=[]

firstq = True
for q in range(1,31):
    s = tp.FindNearestElement(out_p_r_0_left,[4.*float(q/2.),7.5])
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
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        qa = [4*q for i in range(h)]
        if firstq == True:
            ax.plot(qa,sp,'kx',label='rightward-motion detectors')
            firstq = False
        else:
            ax.plot(qa,sp,'kx',label='_nolegend_')
    else:
        store_sp +=[0]


max_spikes_left = max(store_sp)
if max_spikes_left == 0:
    max_spikes_left = [0]
print max_spikes_left


#midget motion detectors complete
ps_all_file = open('data/mo_det_cal/p_max_spikes_all.txt','a+')

store_sp=[]

firstq = True
for q in range(1,31):
    s = tp.FindNearestElement(out_p_r_0,[4.*float(q/2.),7.5])
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
        #plt.subplot(221)
        h, e = np.histogram(sp, bins=np.arange(0., float(t_end-t_start)+1., float(t_end-t_start)))
        store_sp +=[h]
        qa = [4*q for i in range(h)]
        if firstq == True:
            ax.plot(qa,sp,'g*',label='general motion detectors')
            firstq = False
        else:
            ax.plot(qa,sp,'g*',label='_nolegend_')
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

print weight,delay

#plt.title('midget spikes for motion detectors')
#plt.xlabel('position of neuron in arcmin')
#plt.ylabel('time $t$ in ms')
#plt.legend()
plt.xlim([0,120])
plt.ylim([0,200])

ax.set_title('parasolic spikes for motion detectors')
ax.set_xlabel('position of neuron in arcmin')
ax.set_ylabel('time $t$ in ms')
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])


# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), 
          fancybox=False, fontsize=8, shadow=False, ncol=2)

plt.savefig('img/mo_det_cal/x_t_parasols_'+handle_name+'.pdf')
out= '/home/schrader/Documents/microsaccades/img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/x_t_parasols_' + str(handle_name) + '_' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '.pdf'
plt.savefig(out)
plt.show()
plt.close()

k = (delay,vel,max_spikes)
ps_file = open('data/mo_det_cal/p_max_spikes.txt','a+')
ps_file.write(str(weight)+' '+str(delay)+' '+ str(vel)+' '+ str(max_spikes[0])+' '+str(p_poisson_max_spikes[0])+' '+ str(max_spikes_left[0])+' '+ str(max_spikes_right[0])+'\n')
ps_file.close()