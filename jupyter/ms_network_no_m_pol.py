#THIS IS THE PART OF THE NEURAL NETWORK. IT TAKES THE VALUES TO CALCULATE THE MEMBRANE POTENTIALS OF THE INPUT NEURONS AND THEN CALCULATES THE NETWORK OUTPUT 
#the unit in space is 1arcmin!
import pylab as pyl
import numpy as np
import cv2
import sys
import os
import matplotlib.pyplot as plt
import nest
import datetime
import nest.raster_plot
import nest.topology as tp
from microsaccades_functions import *

pgf_with_rc_fonts = {"font.family": "serif","font.serif": [],"font.sans-serif": []}
plt.rcParams.update(pgf_with_rc_fonts)

nest.ResetKernel()   # since we run the script multiple times
#to get different poisson outputs

msd = 1000 #int(np.random.normal(5000,1000,1)[0])  #master seed
nest.SetKernelStatus({'local_num_threads' : 4})
n_vp = nest.GetKernelStatus('total_num_virtual_procs')
msdrange1 = range(msd, msd+n_vp)
pyrngs = [np.random.RandomState(s) for s in msdrange1]
msdrange2 = range(msd+n_vp+1, msd+1+2*n_vp)
nest.SetKernelStatus({'grng_seed': msd+n_vp,
                      'rng_seeds': msdrange2})

I_E=355.

def set_I_e_random(layer):
    r = nest.GetNodes(layer)[0]
    node_info=nest.GetStatus(r)
    localnodes=[(ni['global_id'],ni['vp']) for ni in node_info if ni['local']]
    for gid, vp in localnodes:
        nest.SetStatus([gid], {'I_e' : I_E+pyrngs[vp].uniform(-2.,2.)})
                       
#necessary paramater definitions
#frames = 100 #replaced by motion detectors only

#----------------------------------------------------------------------------INPUT-RATES-FROM-MS_INPUT

extent = 121.
delay = 30. #15. #30. # speed of point in poletti 2010 -> maybe increase a bit for more reacton later on

t_start = 0
t_end = 630 #630 #1000

weight = 40. #39.
weight_std = 0.5 #1.5
#I_E = 410.

mdw=4.0

if len(sys.argv)==11:
    sim_title = sys.argv[1]
    sim_title_2 = sys.argv[2]
    sim_nr = sys.argv[3]
    handle_name = sys.argv[4]
    #extent = float(sys.argv[3])
    extent_x = float(sys.argv[5])
    extent_y = float(sys.argv[6])
    exp_nr = sys.argv[7]
    cond_nr = sys.argv[8]
    sim_length = int(float(sys.argv[9]))
else:
    sim_title = sys.argv[1]
    sim_nr = sys.argv[2]
    handle_name = sys.argv[3]
    #extent = float(sys.argv[3])
    extent_x = float(sys.argv[4])
    extent_y = float(sys.argv[5])
    exp_nr = sys.argv[6]
    cond_nr = sys.argv[7]
    sim_length = int(float(sys.argv[8]))
    sim_title_2 = ''
	
t_end=sim_length


#delay = float(sys.argv[5])
vel = 0. #float(sys.argv[6])
#weight = float(sys.argv[7])
#extent_x = extent
#extent_y = extent
#center = extent/2.
center_x = extent_x/2.
center_y = extent_y/2.

print weight,delay
#weight = 0.

#p_file = open('data/parasolic_values.data','r+')
p_file = open('data/'+str(sim_title)+'/'+str(sim_nr)+'/parasolic_rates_'+str(handle_name)+'_on.data','r+')
p_data = np.load(p_file)  
p_file.close()

parasolic_rates=poissonRateParasols(p_data)

#print midget_rates

#to check for maximum spike rates in order to adopt conversion of film input
#maxs = []
#for i in range(len(midget_rates)):
#    for j in range(len(midget_rates[0])):
#        maxs.append(max(midget_rates[i][j][t_start:t_end]))
mmr = 300000.
#print maxs
#print 'midget max: ' + str(max(maxs)) + ' ' + str(mmr)
        
#paxs = []
#for i in range(len(parasolic_rates)):
#    for j in range(len(parasolic_rates[0])):
#       paxs.append(max(parasolic_rates[i][j][t_start:t_end]))
pmr = 800000. #480000 #480000. #480000.
#print 'parsolic max: ' + str(max(paxs)) + ' ' + str(pmr)

#ps_file = open('data/mo_det_cal/input_max_spikes.txt','a+')
#ps_file.write(str(handle_name)+' m: '+ str(max(maxs))+' p: '+ str(max(paxs))+'\n')
#ps_file.close()
        
#for rates 
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


nest.CopyModel("iaf_psc_alpha", "iaf_psc_alpha_mp") #,{"I_e" : I_E})
nest.CopyModel("iaf_psc_alpha", "iaf_psc_alpha_i",{"I_e" : I_E})
#nest.CopyModel("iaf_psc_alpha", "iaf_psc_alpha_r",{"I_e" : 374.0})
nest.CopyModel("spike_detector", "my_spike_detector",{"withgid": True,
                                                      "withtime": True,
                                                      "to_memory": True,})
#nest.CopyModel('multimeter', 'my_multimeter',{'withtime': True, 'interval': 0.1, 'withgid': True,'record_from': ['Rate']})

#{'distribution' : 'uniform', 'low': weight-1., 'high': weight+1.}
nest.CopyModel('static_synapse','ex',{'weight': 1.0})

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
gp_file = open('data/'+str(sim_title)+str(sim_title_2)+'/'+str(sim_nr)+'/p_pos_'+handle_name+'.data','r+')
#gp_file = open('data/test_vid/p_pos_test_vid.data','r+')
gp_data = np.load(gp_file)  
gp_file.close()
gp_data = gp_data.tolist()


'''
#print gm_pos
print 'minimum midgets'
print gm_min, gm_max[0]+gm_min[0]+10.
print 'maximum midgets'
print gm_max, gm_max[1]+gm_min[1]+10.
'''       
gp_pos=[]
gp_r_0_pos=[]

#mean, stdv, timelength
#p_noise=[[np.random.normal(0,10,t_end) for j in range(len(gp_data[0]))] for i in range(len(gp_data))]
parasolic_rates=pmr*parasolic_rates #+p_noise

for i in range(len(gp_data)):
    for j in range(len(gp_data[0])):
        prs+=[parasolic_rates[i][j]]
        gp_pos+=[[0.5*gp_data[i][j][0]+0.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_0_pos+=[[0.5*gp_data[i][j][0]+2.25,0.5*gp_data[i][j][1]+0.25]]
        #gp_r_60_pos+=[[0.5*gp_data[i][j][0]+1.25,0.5*gp_data[i][j][1]+0.25+1.732]]
        #gp_r_120_pos+=[[gp_data[i][j][0]-1.25,gp_data[i][j][1]+0.25+1.732]]


#gm_pos=[[np.random.uniform(-0.5,0.5), np.random.uniform(-0.5,0.5)] for j in range(50)]
#print gm_pos

#midgets_V_input = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha_mp', 'edge_wrap': True})
#parasolic_V_input = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_mp', 'edge_wrap': True})
parasolic = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_mp', 'edge_wrap': True})

#120 degree commented for the moment since two build a basis

#these are direction dependent
p_reichardt_0_left = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
p_reichardt_0_right = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})

p_reichardt_0_left_vel2 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
p_reichardt_0_right_vel2 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
#p_reichardt_60_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
#p_reichardt_60_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
#p_reichardt_120_up = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})
#p_reichardt_120_down = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_i', 'edge_wrap': True})

#react to motion in both directions
#m_reichardt_0 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_0_pos, 'elements': 'iaf_psc_alpha_r', 'edge_wrap': True})
#m_reichardt_60 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_r_60_pos, 'elements': 'iaf_psc_alpha_r', 'edge_wrap': True})
#m_reichardt_120 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#p_reichardt_0 = tp.CreateLayer({'extent' : [125.,125.], 'center' : [center_x,center_y], 'positions' : gp_r_0_pos, 'elements': 'iaf_psc_alpha_r', 'edge_wrap': True})
#p_reichardt_60 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_r_60_pos, 'elements': 'iaf_psc_alpha_r', 'edge_wrap': True})
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

out_p = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector'})
#out_m_multi = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'my_multimeter'})
#out_p_multi = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_multimeter'})

#out_m_r_0 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_0_pos, 'elements': 'my_spike_detector'})
#out_m_r_0_multi = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_0_pos, 'elements': 'my_multimeter'})
#out_m_r_60 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_60_pos, 'elements': 'my_spike_detector'})
#out_m_r_60_multi = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_r_60_pos, 'elements': 'my_multimeter'})

#out_p_r_0 = tp.CreateLayer({'extent' : [125.,125.], 'center' : [center_x,center_y], 'positions' : gp_r_0_pos, 'elements': 'my_spike_detector'})
#out_p_r_0_multi = tp.CreateLayer({'extent' : [125.,125.], 'center' : [center_x,center_y], 'positions' : gp_r_0_pos, 'elements': 'my_multimeter'})
#out_p_r_60 = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_r_60_pos, 'elements': 'my_spike_detector'})
#out_p_r_60_multi = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_r_60_pos, 'elements': 'my_multimeter'})
#out = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'my_spike_detector'})
#out_multi = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'my_multimeter'})

out_p_r_0_left = tp.CreateLayer({'extent' : [extent_x+4.,extent_y+4.], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector'})
out_p_r_0_right = tp.CreateLayer({'extent' : [extent_x+4.,extent_y+4.], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector'})

out_p_r_0_left_vel2 = tp.CreateLayer({'extent' : [extent_x+4.,extent_y+4.], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector'})
out_p_r_0_right_vel2 = tp.CreateLayer({'extent' : [extent_x+4.,extent_y+4.], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector'})

#out_p_r_60_up = tp.CreateLayer({'extent' : [extent_x+4.,extent_y+4.], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector'})
#out_p_r_60_down = tp.CreateLayer({'extent' : [extent_x+4.,extent_y+4.], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector'})

#out_p_r_120_up = tp.CreateLayer({'extent' : [extent_x+4.,extent_y+4.], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector'})
#out_p_r_120_down = tp.CreateLayer({'extent' : [extent_x+4.,extent_y+4.], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector'})

#randomization---------------------------------------------------------------------------------------------------------------------
 
#set_I_e_random(p_reichardt_0_left)
#set_I_e_random(p_reichardt_0_right)

#---------------------------------------------------------------------------------------------------------------CREATE-CONNECTIONS
#connections to left/right half of Reichardt detector
V_input_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.1,0.1]}}}

#p_r_0_left_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': {'uniform': {'min': weight-weight_std, 'max': weight+weight_std}}, 'mask' : {'rectangular' : {'lower_left' : [-1.1,-0.1], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
#p_r_0_right_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': {'uniform': {'min': weight-weight_std, 'max': weight+weight_std}}, 'mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [1.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
p_r_0_left_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': {'uniform': {'min': weight-weight_std, 'max': weight+weight_std}}, 'mask' : {'rectangular' : {'lower_left' : [-1.1,-0.1], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
p_r_0_right_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': {'uniform': {'min': weight-weight_std, 'max': weight+weight_std}}, 'mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [1.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}

p_r_0_left_vel2_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': {'uniform': {'min': weight-weight_std, 'max': weight+weight_std}}, 'mask' : {'rectangular' : {'lower_left' : [-1.1,-0.1], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : 0.5*delay}}}
p_r_0_right_vel2_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': {'uniform': {'min': weight-weight_std, 'max': weight+weight_std}}, 'mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [1.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : 0.5*delay}}}

#p_r_60_up_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': {'uniform': {'min': weight-weight_std, 'max': weight+weight_std}}, 'mask' : {'rectangular' : {'lower_left' : [-1.,-1.], 'upper_right' : [1.,1.]}, 'anchor' : [0.5,0.866]}, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
#p_r_60_down_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': {'uniform': {'min': weight-weight_std, 'max': weight+weight_std}}, 'mask' : {'rectangular' : {'lower_left' : [-1.,-1.], 'upper_right' : [1.,1.]}, 'anchor' : [-0.5,-0.866]}, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
#p_r_120_up_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': {'uniform': {'min': weight-weight_std, 'max': weight+weight_std}}, 'mask' : {'rectangular' : {'lower_left' : [-1.,-1.], 'upper_right' : [1.,1.]}, 'anchor' : [-0.5,0.866]}, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}
#p_r_120_down_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': {'uniform': {'min': weight-weight_std, 'max': weight+weight_std}}, 'mask' : {'rectangular' : {'lower_left' : [-1.,-1.], 'upper_right' : [1.,1.]}, 'anchor' : [0.5,-0.866]}, 'delays' : {'linear' : {'c' : .1, 'a' : delay}}}

#connections of left/right half to just motion sensitive Reichardt detector
#m_r_0_left_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.3,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
#m_r_0_right_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.3,-0.1], 'upper_right' : [0.1,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
#m_r_60_up_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.3,-0.5], 'upper_right' : [0.1,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
#m_r_60_down_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.3,0.5]}}, 'weights' : weight } #'synapse_model' : 'ex'}
#p_r_0_left_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [1.1,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
#p_r_0_right_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-1.1,-0.1], 'upper_right' : [0.1,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
#p_r_60_up_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-1.25,-1.], 'upper_right' : [0.1,0.1]}}, 'weights' : weight } #'synapse_model' : 'ex'}
#p_r_60_down_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [1.25,1.]}}, 'weights' : weight } #'synapse_model' : 'ex'}

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
#tp.ConnectLayers(midgets_V_input,midgets,V_input_conndict)
#tp.ConnectLayers(parasolic_V_input,parasolic,V_input_conndict)

tp.ConnectLayers(parasolic,p_reichardt_0_right,p_r_0_right_conndict)
tp.ConnectLayers(parasolic,p_reichardt_0_left,p_r_0_left_conndict)

tp.ConnectLayers(parasolic,p_reichardt_0_right_vel2,p_r_0_right_vel2_conndict)
tp.ConnectLayers(parasolic,p_reichardt_0_left_vel2,p_r_0_left_vel2_conndict)
#tp.ConnectLayers(parasolic,p_reichardt_60_up,p_r_60_up_conndict)
#tp.ConnectLayers(parasolic,p_reichardt_60_down,p_r_60_down_conndict)
#tp.ConnectLayers(parasolic,p_reichardt_120_up,p_r_120_up_conndict)
#tp.ConnectLayers(parasolic,p_reichardt_120_down,p_r_120_down_conndict)

#tp.ConnectLayers(m_reichardt_0_left,m_reichardt_0,m_r_0_left_r_0_conndict)
#tp.ConnectLayers(m_reichardt_0_right,m_reichardt_0,m_r_0_right_r_0_conndict)
#tp.ConnectLayers(m_reichardt_60_up,m_reichardt_60,m_r_60_up_r_60_conndict)
#tp.ConnectLayers(m_reichardt_60_down,m_reichardt_60,m_r_60_down_r_60_conndict)

#tp.ConnectLayers(p_reichardt_0_left,p_reichardt_0,p_r_0_left_r_0_conndict)
#tp.ConnectLayers(p_reichardt_0_right,p_reichardt_0,p_r_0_right_r_0_conndict)
#tp.ConnectLayers(p_reichardt_60_up,p_reichardt_60,p_r_60_up_r_60_conndict)
#tp.ConnectLayers(p_reichardt_60_down,p_reichardt_60,p_r_60_down_r_60_conndict)

tp.ConnectLayers(parasolic,out_p,out_conndict)
#tp.ConnectLayers(midgets,out_m_multi,out_conndict)
#tp.ConnectLayers(parasolic,out_p_multi,out_conndict)

tp.ConnectLayers(p_reichardt_0_left,out_p_r_0_left,out_conndict)
tp.ConnectLayers(p_reichardt_0_right,out_p_r_0_right,out_conndict)

tp.ConnectLayers(p_reichardt_0_left_vel2,out_p_r_0_left_vel2,out_conndict)
tp.ConnectLayers(p_reichardt_0_right_vel2,out_p_r_0_right_vel2,out_conndict)

#tp.ConnectLayers(p_reichardt_60_up,out_p_r_60_up,out_conndict)
#tp.ConnectLayers(p_reichardt_60_down,out_p_r_60_down,out_conndict)
#tp.ConnectLayers(p_reichardt_120_up,out_p_r_120_up,out_conndict)
#tp.ConnectLayers(p_reichardt_120_down,out_p_r_120_down,out_conndict)

#tp.ConnectLayers(m_reichardt_0,out_m_r_0,out_conndict)
#tp.ConnectLayers(p_reichardt_0,out_p_r_0,out_conndict)
#tp.ConnectLayers(out_m_r_0_multi,m_reichardt_0,out_conndict)
#tp.ConnectLayers(p_reichardt_0,out_p_r_0_multi,out_conndict)

#tp.ConnectLayers(m_reichardt_60,out_m_r_60,out_conndict)
#tp.ConnectLayers(p_reichardt_60,out_p_r_60,out_conndict)
#tp.ConnectLayers(m_reichardt_60,out_m_r_60_multi,out_conndict)
#tp.ConnectLayers(p_reichardt_60,out_p_r_60_multi,out_conndict)

#ctr = tp.FindNearestElement(m_reichardt_0,[6.,12.])
#fig = tp.PlotLayer(midgets,nodesize =80)
#tp.PlotTargets(ctr,midgets,fig=fig,mask=m_r_0_right_conndict['mask'],src_size=250, tgt_color='red',tgt_size=20)

#fig = tp.PlotLayer(m_reichardt_60_up)
#ctr = tp.FindNearestElement(midgets,[12.,7.5])
#tp.PlotTargets(ctr, m_reichardt_60_up ,fig = fig, mask=m_r_60_down_conndict['mask'], src_size=200, src_color = 'red', tgt_color = 'green' , tgt_size=20 )
#plt.show()




#-------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------SIMULATION

#for updates
PIDs = nest.GetNodes(parasolic)
SIE = 374.7
STDI = 0.45

for f in range(t_start,t_end):#frames):
    print f
    #reset rates
    for n in range(len(PIDs[0])):
        qr = prs[n][f]
        nest.SetStatus([PIDs[0][n]], {'I_e': qr})
    #run simulation
    nest.Simulate(1)



#spontaneous firing rate
'''
if qr < 375.:
	qr = SIE+np.random.normal(0,STDI,1)[0]
if f<5:
	qr = 374.+np.random.normal(0,1.7,1)[0]
'''	
'''
MIDs = nest.GetNodes(midgets)
for n in range(len(MIDs[0])):
    qr = I_E+np.random.normal(0,10,1)[0]
    nest.SetStatus([MIDs[0][n]], {'I_e': qr})
nest.Simulate(200)
'''

#-------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------SAVING-AND-PRINTING-THE-OUTPUT

max_evs_list = []
max_evs_list_200 = []
all_spikes = 0
all_spikes_200 = 0
p_all_spikes = 0
p_all_spikes_200 = 0
p_all_times=[]
p_left_times=[]
p_right_times=[]
p_left_vel2_times=[]
p_right_vel2_times=[]
p_ud_times=[]
p_all_times_200=[]
p_left_times_200=[]
p_right_times_200=[]
p_ud_times_200=[]

def save_spikes(layer_name,layer,mel,asl,tl):
    directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
    sp_file = open(directory+'/'+layer_name+'_'+handle_name+'.data','w+')
    n_evs_l=[0]
    times = []
    for n in range(len(layer[0])):
        n_evs = nest.GetStatus([layer[0][n]],"n_events")[0]
        if n_evs>0:
            #print layer_name+' '+str(layer[0][n])+' '+str(n_evs)+' '+str(nest.GetStatus([layer[0][n]],"events")[0]["times"])
            t_evs = nest.GetStatus([layer[0][n]],"events")[0]["times"]
            sp_file.write(str(layer[0][n])+'\t'+str(nest.GetStatus([layer[0][n]],"events")[0]["times"])+'\n')
            for t in t_evs:
                times+=[t]
            n_evs_l+=[n_evs]
    mel+=[max(n_evs_l)]
    #if layer_name != ('spikes_midgets' or 'spikes_parasols'):
    asl+=sum(n_evs_l) 
    tl=set(times)
    sp_file.close()
    print layer_name
    return 0

def save_spikes_200(layer_name,layer,mel200,asl200,tl200):
    directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
    sp_file = open(directory+'/'+layer_name+'_200_'+handle_name+'.data','w+')
    n_evs_200_l=[0]
    times=[]
    for n in range(len(layer[0])):
        n_evs = nest.GetStatus([layer[0][n]],"n_events")[0]
        if n_evs>0:
            t_evs = nest.GetStatus([layer[0][n]],"events")[0]["times"]
            t_evs_200 = [t_ev for t_ev in t_evs if t_ev>200.]
            if len(t_evs_200)>0:
                #print layer[0][n],t_evs_200
                sp_file.write(str(layer[0][n])+'\t'+str(t_evs_200)+'\n')
                n_evs_200_l+=[len(t_evs_200)]
                for t in t_evs:
                    times+=[t]
    mel200+=[max(n_evs_200_l)]
    #if layer_name != ('spikes_midgets' or 'spikes_parasols'):
    asl200+=sum(n_evs_200_l) 
    tl200=set(times)
    sp_file.close()
    print layer_name + '200'
    
    return 0

GID_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/GID_info.txt','w+')

#parasolic
layerIDs = nest.GetNodes(out_p)
GID_file.write('parasols \t'+str(layerIDs[0][0])+'\t'+str(layerIDs[0][len(layerIDs[0])-1])+'\n')
save_spikes('spikes_parasols',layerIDs,max_evs_list,all_spikes,p_all_times)
save_spikes_200('spikes_parasols',layerIDs,max_evs_list_200,all_spikes_200,p_all_times_200)

#parasolic rightward motion detectors
layerIDs = nest.GetNodes(out_p_r_0_left)
GID_file.write('p_0_left \t'+str(layerIDs[0][0])+'\t'+str(layerIDs[0][len(layerIDs[0])-1])+'\n')
save_spikes('spikes_p_0_left',layerIDs,max_evs_list,p_all_spikes,p_left_times)
save_spikes_200('spikes_p_0_left',layerIDs,max_evs_list_200,p_all_spikes_200,p_left_times_200)
#parasolic leftward motion detectors
layerIDs = nest.GetNodes(out_p_r_0_right)
GID_file.write('p_0_right \t'+str(layerIDs[0][0])+'\t'+str(layerIDs[0][len(layerIDs[0])-1])+'\n')
save_spikes('spikes_p_0_right',layerIDs,max_evs_list,p_all_spikes,p_right_times)
save_spikes_200('spikes_p_0_right',layerIDs,max_evs_list_200,p_all_spikes_200,p_right_times_200)
#parasolic rightward motion detectors vel2
layerIDs = nest.GetNodes(out_p_r_0_left_vel2)
GID_file.write('p_0_left_vel2 \t'+str(layerIDs[0][0])+'\t'+str(layerIDs[0][len(layerIDs[0])-1])+'\n')
save_spikes('spikes_p_0_left_vel2',layerIDs,max_evs_list,p_all_spikes,p_left_times)
save_spikes_200('spikes_p_0_left_vel2',layerIDs,max_evs_list_200,p_all_spikes_200,p_left_times_200)
#parasolic leftward motion detectors vel2
layerIDs = nest.GetNodes(out_p_r_0_right_vel2)
GID_file.write('p_0_right_vel2 \t'+str(layerIDs[0][0])+'\t'+str(layerIDs[0][len(layerIDs[0])-1])+'\n')
save_spikes('spikes_p_0_right_vel2',layerIDs,max_evs_list,p_all_spikes,p_right_times)
save_spikes_200('spikes_p_0_right_vel2',layerIDs,max_evs_list_200,p_all_spikes_200,p_right_times_200)

GID_file.close()
#print max_evs_list
#print max_evs_list_200

#ps_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_max_spikes.txt','a+')
#ps_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+str(max_evs_list[0])+'\t'+str(max_evs_list[1])+'\t'+str(max_evs_list[2])+'\t'+str(max_evs_list[3])+'\t'+str(max_evs_list[4])+'\t'+str(max_evs_list[5])+'\t'+str(max_evs_list[6])+'\t'+str(p_all_spikes)+'\n')
#ps_file.close()

ps_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_spikes_lr.txt','a+')
ps_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+str(len(p_left_times))+'\t'+str(len(p_right_times))+'\t'+str(len(p_left_times)+len(p_right_times))+'\t'+'\t'+str(len(p_left_vel2_times))+'\t'+str(len(p_right_vel2_times))+'\t'+str(len(p_left_vel2_times)+len(p_right_vel2_times))+'\t'+'\n')
ps_file.close()
ps_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_spike_times_lr.txt','a+')
ps_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+str(p_left_times)+'\t'+str(p_right_times)+'\t'+str(p_left_vel2_times)+'\t'+str(p_right_vel2_times)+'\n')
ps_file.close()


#ps_file_200 = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_max_spikes_200.txt','a+')
#ps_file_200.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+str(max_evs_list_200[1])+'\t'+str(max_evs_list_200[8])+'\t'+str(max_evs_list_200[9])+'\t'+str(max_evs_list_200[10])+'\t'+str(max_evs_list_200[11])+'\t'+str(max_evs_list_200[12])+'\t'+str(max_evs_list_200[13])+'\t'+str(p_all_spikes_200)+'\n')
#ps_file_200.close()
'''
ps_200_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_spikes_lr_200.txt','a+')
ps_200_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+str(len(p_left_times))+'\t'+str(len(p_right_times))+'\t'+str(len(p_left_times)+len(p_right_times))+'\t'+'\n')
ps_200_file.close()
ps_200_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_spike_times_lr_200.txt','a+')
ps_200_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+str(p_left_times)+'\t'+str(p_right_times)+'\n')
ps_200_file.close()
'''
print 'finished'

parasols = 0
p_reichardt_0_left = 0
p_reichardt_0_right = 0

out_p = 0
out_p_r_0_left = 0
out_p_r_0_right = 0

sys.exit()