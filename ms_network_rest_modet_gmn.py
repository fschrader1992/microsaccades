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

msd = int(np.random.normal(5000,1000,1)[0])  #master seed
nest.SetKernelStatus({'local_num_threads' : 4})
n_vp = nest.GetKernelStatus('total_num_virtual_procs')
msdrange1 = range(msd, msd+n_vp)
pyrngs = [np.random.RandomState(s) for s in msdrange1]
msdrange2 = range(msd+n_vp+1, msd+1+2*n_vp)
nest.SetKernelStatus({'grng_seed': msd+n_vp,
                      'rng_seeds': msdrange2})
                       
#necessary paramater definitions
#frames = 100 #replaced by motion detectors only

#----------------------------------------------------------------------------INPUT-RATES-FROM-MS_INPUT

extent = 121.
delay = 30. # speed of point in poletti 2010 -> maybe increase a bit for more reacton later on

t_start = 0
t_end = 1000

weight = 40.
weight_std = 1.5
#I_E = 410.

if len(sys.argv)==9:
    sim_title = sys.argv[1]
    sim_title_2 = sys.argv[2]
    sim_nr = sys.argv[3]
    handle_name = sys.argv[4]
    #extent = float(sys.argv[3])
    extent_x = float(sys.argv[5])/2.
    extent_y = float(sys.argv[6])/2.-1.
    exp_nr = sys.argv[7]
    cond_nr = sys.argv[8]
    sar=True
else:
    sim_title = sys.argv[1]
    sim_nr = sys.argv[2]
    handle_name = sys.argv[3]
    #extent = float(sys.argv[3])
    extent_x = float(sys.argv[4])/2.
    extent_y = float(sys.argv[5])/2.-1.
    exp_nr = sys.argv[6]
    cond_nr = sys.argv[7]
    sim_title_2 = ''
    sar=False
    
#delay = float(sys.argv[5])
vel = 0. #float(sys.argv[6])
#weight = float(sys.argv[7])
#extent_x = extent
#extent_y = extent
#center = extent/2.
#extent_x = extent_x/2.
#extent_y = extent_y/2.
center_x = extent_x/2.
center_y = extent_y/2.
#extent_x = extent_x+1.
#extent_y = extent_y+1.

print weight,delay
#weight = 0.

#p_file = open('data/parasolic_values.data','r+')
p_file = open('data/'+str(sim_title)+str(sim_title_2)+'/parasolic_rates_'+str(handle_name)+'_on.data','r+') 
p_data = np.load(p_file)  
p_file.close()

parasolic_rates=poissonRateParasols(p_data)

#print maxs
#print 'midget max: ' + str(max(maxs)) + ' ' + str(mmr)
        
#paxs = []
#for i in range(len(parasolic_rates)):
#    for j in range(len(parasolic_rates[0])):
#       paxs.append(max(parasolic_rates[i][j][t_start:t_end]))
pmr = 480000.
#print 'parsolic max: ' + str(max(paxs)) + ' ' + str(pmr)

#ps_file = open('data/mo_det_cal/input_max_spikes.txt','a+')
#ps_file.write(str(handle_name)+' m: '+ str(max(maxs))+' p: '+ str(max(paxs))+'\n')
#ps_file.close()
        
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

nest.CopyModel("iaf_cond_alpha", "std_iaf_cond_alpha",{'tau_syn_ex': 1.0, 'V_reset': -70.0})
nest.CopyModel("iaf_cond_alpha", "tr8_iaf_cond_alpha",{'tau_syn_ex': 1.0, 'V_reset': -70.0, 't_ref': 8.})
nest.CopyModel("spike_detector", "my_spike_detector",{'withtime': True, 'withgid': True})
nest.CopyModel('static_synapse','ex') #,{'weight': {'distribution' : 'uniform', 'low': .1, 'high': .2}})
nest.CopyModel('static_synapse','inh_400', {'weight': -400.,}) #-400.,}) 
nest.CopyModel('static_synapse','ex_50_delay', {'weight': 50., 'delay': 12.,}) #50
#stimes=[]
#stimes+=[i*20.0+10. for i in range(10)]
nest.CopyModel('spike_generator', 'my_spike_generator')

#nest.SetDefaults('iaf_psc_alpha_r',{'I_e' : 374.0})
#nest.SetDefaults('iaf_psc_alpha',{'I_e' : 350.0})
#nest.CopyModel("static_synapse", "ex", {"weight" : 0.1})
#nest.CopyModel("static_synapse", "inh", {"weight" : -0.1})

#----------------------------------------------------------------------------------------CREATE-LAYERS

gp_file = open('data/'+str(sim_title)+str(sim_title_2)+'/p_pos_'+handle_name+'.data','r+')
#gp_file = open('data/test_vid/p_pos_test_vid.data','r+')
gp_data = np.load(gp_file)  
gp_file.close()
gp_data = gp_data.tolist()
   
gp_pos=[]
gpgm_pos=[]
gp_r_0_pos=[]
gp_r_60_pos=[]
#gm_r_120_pos=[]

for i in range(len(gp_data)):
    for j in range(len(gp_data[0])):
        gp_pos+=[[0.5*gp_data[i][j][0]+0.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_0_pos+=[[0.5*gp_data[i][j][0]+2.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_60_pos+=[[0.5*gp_data[i][j][0]+1.25,0.5*gp_data[i][j][1]+0.25+1.732]]
        if i < len(gp_data)/5 and j < len(gp_data[0])/5:
            gpgm_pos+=[[5*0.5*gp_data[i][j][0]+0.25,5*0.5*gp_data[i][j][1]+0.25]]
        #gp_r_120_pos+=[[gp_data[i][j][0]-1.25,gp_data[i][j][1]+0.25+1.732]]


#gm_pos=[[np.random.uniform(-0.5,0.5), np.random.uniform(-0.5,0.5)] for j in range(50)]
#print gm_pos

#midgets_V_input = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha_mp', 'edge_wrap': True})
#parasolic_V_input = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_mp', 'edge_wrap': True})

#Test, WHETHER THAT WORKS!
OGIDs = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/GID_info.txt','r+')

layers=[] #list of layers
layers2=[]
layers3=[]
gmdet_layers=[]
out_gmdet_layers=[]
gmn_layers=[]
out_gmn_layers=[]
gmc_layers=[] #list of global motion corrected layers
out_layers=[]
layer_names = ['parasols','p_0_left','p_0_right','p_60_up','p_60_down','p_120_up','p_120_down']


for ln in layer_names:
    lr=tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_generator', 'edge_wrap': True})
    lr2=tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_generator', 'edge_wrap': True})
    lr3=tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_generator', 'edge_wrap': True})
    if sar==True:
        l_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(handle_name)+'.data','r+')
    else:
        l_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(handle_name)+'.data','r+')
    LIDs = nest.GetNodes(lr)
    LIDs2 = nest.GetNodes(lr2)
    LIDs3 = nest.GetNodes(lr3)
    deltaID = 0
    #print '---NEW-FILE---'
    OGIDs.seek(0)
    for ogid in OGIDs:
        ogid=ogid.split('\t')
        #print str(ogid[0].strip())
        #print str(ln.strip())
        #print '---'
        if str(ogid[0].strip())==str(ln.strip()):
            deltaID = int(ogid[1])
    lq=0
    for line in l_file:
        nid = int(line.split('[')[0])
        st = line.split('[')[1]
        st = st.split(']')[0]
        st = filter(None, st.split(' '))
        sta=[]
        for s in st: 
            sta+=[float(s.strip())]
        #print ln
        #print np.array(sta)
        #print nid-deltaID
        nest.SetStatus([LIDs[0][nid-deltaID]], {'spike_times': np.array(sta)})
        nest.SetStatus([LIDs2[0][nid-deltaID]], {'spike_times': np.array(sta)})
        nest.SetStatus([LIDs3[0][nid-deltaID]], {'spike_times': np.array(sta)})
        lq+=1        
    l_file.close()
    if ln != 'parasols':
        layers+=[lr]
        layers2+=[lr2]
        layers3+=[lr3]
        gmc_lr=tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'std_iaf_cond_alpha', 'edge_wrap': True})
        gmc_layers+=[gmc_lr]
        olr=tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'my_spike_detector', 'edge_wrap': True})
        out_layers+=[olr]
        gmdet=tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gpgm_pos, 'elements': 'std_iaf_cond_alpha', 'edge_wrap': True})
        gmdet_layers+=[gmdet]
        gmn=tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : [[center_x, center_y]], 'elements': 'std_iaf_cond_alpha', 'edge_wrap': True})
        gmn_layers+=[gmn]
        out_gmn=tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : [[center_x, center_y]], 'elements': 'my_spike_detector', 'edge_wrap': True})
        out_gmn_layers+=[out_gmn]
        out_gmdet=tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gpgm_pos, 'elements': 'my_spike_detector', 'edge_wrap': True})
        out_gmdet_layers+=[out_gmdet]
    else:
        parasolic = lr
        psd = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'tr8_iaf_cond_alpha', 'edge_wrap': True})
    #print lr
            
OGIDs.close()

#---------------------------------------------------------------------------------------------------------------CREATE-CONNECTIONS
#connections to left/right half of Reichardt detector

parasolic_to_psd_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': 40., 'mask' : {'rectangular' : {'lower_left' : [-0.2,-0.2], 'upper_right' : [0.2,0.2]}}}
modet_to_gmdet_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': 10., 'mask' : {'rectangular' : {'lower_left' : [-5.1,-5.1*0.866], 'upper_right' : [4.9,4.9*0.866]}}}
#gmdet_to_gmc_conndict = {'connection_type' : 'divergent', 'synapse_model': 'inh_400', 'mask' : {'rectangular' : {'lower_left' : [-4.9,-4.9*0.866], 'upper_right' : [5.1,5.1*0.866]}}} #{'lower_left' : [-1.2,-2.4*0.866], 'upper_right' : [1.3,2.6*0.866]}}}
modet_to_gmc_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex_50_delay', 'mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.1,0.1]}}}

gmdet_to_gmn_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': 40.}
gmn_to_gmc_conndict = {'connection_type' : 'divergent', 'synapse_model': 'inh_400'} #{'lower_left' : [-1.2,-2.4*0.866], 'upper_right' : [1.3,2.6*0.866]}}}

out_conndict = {'connection_type' : 'convergent', 'mask' : {'rectangular' : {'lower_left' : [-0.2,-0.2], 'upper_right' : [0.2,0.2]}}}

#connect them
#tp.ConnectLayers(midgets_V_input,midgets,V_input_conndict)
#tp.ConnectLayers(parasolic_V_input,parasolic,V_input_conndict)

#tp.ConnectLayers(parasolic,psd,parasolic_to_psd_conndict)
#tp.ConnectLayers(layers[6],gmdet_layers[5],modet_to_gmdet_conndict)
#tp.ConnectLayers(gmdet,out_gmdet,out_conndict)
#tp.ConnectLayers(psd,gmdet,psd_to_gmn_conndict)
#tp.ConnectLayers(gmn,out_gmn,out_conndict)


for i in range(len(layers)):
    tp.ConnectLayers(layers[i],gmc_layers[i],modet_to_gmc_conndict)
    
    tp.ConnectLayers(layers2[i],gmdet_layers[i],modet_to_gmdet_conndict)
    tp.ConnectLayers(gmdet_layers[i],out_gmdet_layers[i],out_conndict)

    tp.ConnectLayers(layers3[i],gmn_layers[i],modet_to_gmn_conndict)
    tp.ConnectLayers(gmn_layers[i],out_gmn_layers[i],out_conndict)
    
    tp.ConnectLayers(gmdet_layers[i],gmc_layers[i],gmdet_to_gmc_conndict)
    tp.ConnectLayers(gmn_layers[i],gmc_layers[i],gmn_to_gmc_conndict)
    
    tp.ConnectLayers(gmc_layers[i],out_layers[i],out_conndict)
   
#ctr = tp.FindNearestElement(parasolic,[120.,120.])
#fig = tp.PlotLayer(parasolic,nodesize=80)
#tp.PlotTargets(ctr,psd,fig=fig,mask=out_conndict['mask'],src_size=250, tgt_color='red',tgt_size=20)

#fig = tp.PlotLayer(psd)
#ctr = tp.FindNearestElement(parasolic,[120.,120.])
#tp.PlotTargets(ctr, psd ,fig = fig, mask=out_conndict['mask'], src_size=200, src_color = 'red', tgt_color = 'green' , tgt_size=20 )
#plt.show()


#-------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------SIMULATION
       
#run simulation
nest.Simulate(t_end-t_start)
  
'''
MIDs = nest.GetNodes(midgets)
for n in range(len(MIDs[0])):
    qr = I_E+np.random.normal(0,10,1)[0]
    nest.SetStatus([MIDs[0][n]], {'I_e': qr})
nest.Simulate(200)
'''

#-------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------SAVING-AND-PRINTING-THE-OUTPUT

def save_spikes(layer_name,layer,asl,tl):
    directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
    sp_file = open(directory+'/'+layer_name+'_'+handle_name+'_modet.data','w+')
    n_evs_l=[0]
    times = []
    for n in range(len(layer[0])):
        n_evs = nest.GetStatus([layer[0][n]],"n_events")[0]
        if n_evs>0:
            #print layer_name+' '+str(layer[0][n])+' '+str(n_evs)+' '+str(nest.GetStatus([layer[0][n]],"events")[0]["times"])
            t_evs = nest.GetStatus([layer[0][n]],"events")[0]["times"]
            sp_file.write(str(layer[0][n])+'\t'+str(t_evs)+'\n')
            for t in t_evs:
                times+=[int(t+0.500)]
            n_evs_l+=[n_evs]
    #if layer_name != ('spikes_midgets' or 'spikes_parasols'):
    tl=set(times)
    asl=len(tl) 
    sp_file.close()
    print tl, asl
    print layer_name
    return asl, tl

def save_spikes_200(layer_name,layer,asl200,tl200):
    directory = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)
    sp_file = open(directory+'/'+layer_name+'_200_'+handle_name+'_corr_modet.data','w+')
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
    tl200=set(times)
    asl200=len(tl200) 
    sp_file.close()
    print layer_name + '200'
    
    return asl200, tl200

GID_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/GID_info_glob_modet_weight.txt','w+')

gmdet_times = []
gmdet_sum_spikes = []
gmn_times = []
gmn_sum_spikes = []
modet_times = []
modet_sum_spikes = []

gmdet_times200 = []
gmn_sum_spikes200 = []
gmn_times200 = []
gmdet_sum_spikes200 = []
modet_times200 = []
modet_sum_spikes200 = []

#parasolic
print '---EVALUATION---'
'''
layerIDs = nest.GetNodes(out_gmdet)
GID_file.write('gmdet_corr \t'+str(layerIDs[0][0])+'\t'+str(layerIDs[0][len(layerIDs[0])-1])+'\n')
gmdet_sum_spikes,gmdet_times = save_spikes('spikes_gmdet_corr',layerIDs,gmdet_sum_spikes,gmdet_times)
gmdet_sum_spikes200,gmdet_times200 = save_spikes_200('spikes_gmdet_corr',layerIDs,gmdet_sum_spikes200,gmdet_times200)

layerIDs = nest.GetNodes(out_gmn)
GID_file.write('gmn_corr \t'+str(layerIDs[0][0])+'\t'+str(layerIDs[0][len(layerIDs[0])-1])+'\n')
gmn_sum_spikes,gmn_times = save_spikes('spikes_gmn_corr',layerIDs,gmn_sum_spikes,gmn_times)
gmn_sum_spikes200,gmn_times200 = save_spikes_200('spikes_gmn_corr',layerIDs,gmn_sum_spikes200,gmn_times200)
'''

#parasolic rightward motion detectors
for i in range(len(layers)):
    modet_times+=[0]
    modet_sum_spikes+=[0]
    modet_times200+=[0]
    modet_sum_spikes200+=[0]
    layerIDs = nest.GetNodes(out_layers[i])
    GID_file.write(str(layer_names[i+1])+'_corr_modet_weight \t'+str(layerIDs[0][0])+'\t'+str(layerIDs[0][len(layerIDs[0])-1])+'\n') #because of parasolic
    modet_sum_spikes[i],modet_times[i] = save_spikes('spikes_'+str(layer_names[i+1])+'_corr_modet_weight',layerIDs,modet_sum_spikes[i],modet_times[i])
    modet_sum_spikes200[i],modet_times200[i] = save_spikes_200('spikes_'+str(layer_names[i+1])+'_corr_modet_weight',layerIDs,modet_sum_spikes200[i],modet_times200[i])

    gmdet_times+=[0]
    gmdet_sum_spikes+=[0]
    gmdet_times200+=[0]
    gmdet_sum_spikes200+=[0]
    layerIDs = nest.GetNodes(out_gmdet_layers[i])
    GID_file.write('gmdet_corr_modet_weight \t'+str(layerIDs[0][0])+'\t'+str(layerIDs[0][len(layerIDs[0])-1])+'\n')
    gmdet_sum_spikes[i],gmdet_times[i] = save_spikes('spikes_gmdet_'+str(layer_names[i+1])+'_corr_modet_weight',layerIDs,gmdet_sum_spikes[i],gmdet_times[i])
    gmdet_sum_spikes200[i],gmdet_times200[i] = save_spikes_200('spikes_gmdet'+str(layer_names[i+1])+'_corr_modet_weight',layerIDs,gmdet_sum_spikes200[i],gmdet_times200[i])

    gmn_times+=[0]
    gmn_sum_spikes+=[0]
    gmn_times200+=[0]
    gmn_sum_spikes200+=[0]
    layerIDs = nest.GetNodes(out_gmn_layers[i])
    GID_file.write('gmn_corr_modet_weight \t'+str(layerIDs[0][0])+'\t'+str(layerIDs[0][len(layerIDs[0])-1])+'\n')
    gmn_sum_spikes[i],gmn_times[i] = save_spikes('spikes_gmn_'+str(layer_names[i+1])+'_corr_modet_weight',layerIDs,gmn_sum_spikes[i],gmn_times[i])
    gmn_sum_spikes200[i],gmn_times200[i] = save_spikes_200('spikes_gmn'+str(layer_names[i+1])+'_corr_modet_weight',layerIDs,gmn_sum_spikes200[i],gmn_times200[i])
        
GID_file.close()

all_spikes = sum(modet_sum_spikes)
all_spikes200 = sum(modet_sum_spikes200)
gmdet_spikes = sum(modet_sum_spikes)
gmdet_spikes200 = sum(modet_sum_spikes200)
gmn_spikes = sum(modet_sum_spikes)
gmn_spikes200 = sum(modet_sum_spikes200)
print 'SUM:'
print all_spikes

ps_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_corr_max_spikes_modet_weight.txt','a+')
ps_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+str(gmn_spikes)+'\t'+str(gmdet_spikes)+'\t'+str(all_spikes)+'\t'+str(modet_sum_spikes[0])+'\t'+str(modet_sum_spikes[1])+'\t'+str(modet_sum_spikes[2])+'\t'+str(modet_sum_spikes[3])+'\t'+str(modet_sum_spikes[4])+'\t'+str(modet_sum_spikes[5])+'\n')
ps_file.close()

ps_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_corr_spike_times_modet_weight.txt','a+')
ps_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\n'+str(modet_times[0])+'\n'+str(modet_times[1])+'\n'+str(modet_times[2])+'\n'+str(modet_times[3])+'\n'+str(modet_times[4])+'\n'+str(modet_times[5])+'\n')
ps_file.close()

ps200_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_corr_max_spikes_modet_weight_200.txt','a+')
ps200_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+str(gmn_spikes200)+'\t'+str(gmdet_spikes200)+'\t'+str(all_spikes200)+'\t'+str(modet_sum_spikes200[0])+'\t'+str(modet_sum_spikes200[1])+'\t'+str(modet_sum_spikes200[2])+'\t'+str(modet_sum_spikes200[3])+'\t'+str(modet_sum_spikes200[4])+'\t'+str(modet_sum_spikes200[5])+'\n')
ps200_file.close()

ps200_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_corr_spike_times_modet_weight_200.txt','a+')
ps200_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\n'+str(modet_times200[0])+'\n'+str(modet_times200[1])+'\n'+str(modet_times200[2])+'\n'+str(modet_times200[3])+'\n'+str(modet_times200[4])+'\n'+str(modet_times200[5])+'\n')
ps200_file.close()

print 'finished'

for lr in layers:
    lr = 0
for lr in out_layers:
    lr = 0
for lr in gmc_layers:
    lr = 0
for lr in gmdet_layers:
    lr = 0
for lr in gmn_layers:
    lr = 0
for lr in out_gmdet_layers:
    lr = 0
for lr in out_gmn_layers:
    lr = 0

parasolic = 0
psd = 0

sys.exit()