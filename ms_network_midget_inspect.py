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

if len(sys.argv)==10:
    sim_title = sys.argv[1]
    sim_title_2 = sys.argv[2]
    sim_nr = sys.argv[3]
    handle_name = sys.argv[4]
    #extent = float(sys.argv[3])
    extent_x = float(sys.argv[5])/2.
    extent_y = float(sys.argv[6])/2.-1.
    exp_nr = sys.argv[7]
    cond_nr = sys.argv[8]
    mdw = float(sys.argv[9])
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
    mdw = float(sys.argv[8])
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
#p_file = open('data/'+str(sim_title)+str(sim_title_2)+'/'+str(sim_nr)+'/parasolic_rates_'+str(handle_name)+'_on.data','r+') 
#p_data = np.load(p_file)  
#p_file.close()

#parasolic_rates=poissonRateParasols(p_data)

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
nest.CopyModel('static_synapse','inh_400_delay', {'weight': -400., 'delay': 15.,}) #-400.,}) 
nest.CopyModel('static_synapse','ex_50', {'weight': 50.,}) #50
#stimes=[]
#stimes+=[i*20.0+10. for i in range(10)]
nest.CopyModel('spike_generator', 'my_spike_generator')

#nest.SetDefaults('iaf_psc_alpha_r',{'I_e' : 374.0})
#nest.SetDefaults('iaf_psc_alpha',{'I_e' : 350.0})
#nest.CopyModel("static_synapse", "ex", {"weight" : 0.1})
#nest.CopyModel("static_synapse", "inh", {"weight" : -0.1})

#----------------------------------------------------------------------------------------CREATE-LAYERS

gm_file = open('data/'+str(sim_title)+str(sim_title_2)+'/'+str(sim_nr)+'/m_pos_'+handle_name+'.data','r+')
#gm_file = open('data/test_vid/p_pos_test_vid.data','r+')
gm_data = np.load(gm_file)  
gm_file.close()
gm_data = gm_data.tolist()
   
gm_pos=[]
gmgm_pos=[]
gm_r_0_pos=[]
gm_r_60_pos=[]
#gm_r_120_pos=[]

'''
for i in range(len(gp_data)):
    for j in range(len(gp_data[0])):
        gp_pos+=[[0.5*gp_data[i][j][0]+0.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_0_pos+=[[0.5*gp_data[i][j][0]+2.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_60_pos+=[[0.5*gp_data[i][j][0]+1.25,0.5*gp_data[i][j][1]+0.25+1.732]]
        if i < len(gp_data)/5 and j < len(gp_data[0])/5:
            gpgm_pos+=[[5*0.5*gp_data[i][j][0]+0.25,5*0.5*gp_data[i][j][1]+0.25]]
        #gp_r_120_pos+=[[gp_data[i][j][0]-1.25,gp_data[i][j][1]+0.25+1.732]]
'''
        
for i in range(len(gm_data)):
    for j in range(len(gm_data[0])):
        gm_pos+=[[0.5*gm_data[i][j][0]+0.25,0.5*gm_data[i][j][1]+0.25]]
        gm_r_0_pos+=[[0.5*gm_data[i][j][0]+0.5,0.5*gm_data[i][j][1]+0.25]]
        gm_r_60_pos+=[[0.5*gm_data[i][j][0]+0.5,0.5*gm_data[i][j][1]+0.25+0.433]]


#gm_pos=[[np.random.uniform(-0.5,0.5), np.random.uniform(-0.5,0.5)] for j in range(50)]
#print gm_pos

#midgets_V_input = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha_mp', 'edge_wrap': True})
#parasolic_V_input = tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha_mp', 'edge_wrap': True})

#Test, WHETHER THAT WORKS!
OGIDs = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/GID_info.txt','r+')
#OGIDs_par = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/GID_info.txt','r+')


layers=[] #list of layers
gmc_layers=[] #list of global motion corrected layers
out_layers=[]
layer_names = ['midgets','m_0_left','m_0_right'] # ,'p_0_left_corr','p_0_right_corr'] #,'p_60_up','p_60_down','p_120_up','p_120_down']

midpos_x=[]
midpos_y=[]
fmdpos_x=[]
fmdpos_y=[]

mid_pos_times_list = []
ml_pos_times_list = []
mr_pos_times_list = []



for ln in layer_names:
    lr=tp.CreateLayer({'extent' : [extent_x,extent_y], 'center' : [center_x,center_y], 'positions' : gm_pos, 'elements': 'my_spike_generator', 'edge_wrap': True})
    #if ln == 'parasols':
    if sar==True:
        l_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(handle_name)+'.data','r+')
    else:
        l_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(handle_name)+'.data','r+')
    '''
    else:
        if sar==True:
            #l_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(handle_name)+'_mdw_'+str(mdw)+'.data','r+')
            l_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(handle_name)+'_modet_gmn_40_30.data','r+')
        else:
            #l_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(handle_name)+'_mdw_'+str(mdw)+'.data','r+')
            l_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(handle_name)+'_modet_gmn_40_30.data','r+')
    '''
    #print ln
    #print '---'
    LIDs = nest.GetNodes(lr)
    deltaID = 0
    #print '---NEW-FILE---'
    '''
    if ln == 'parasols':
        OGIDs_par.seek(0)
        for ogid in OGIDs_par:
            ogid=ogid.split('\t')
            #print str(ogid[0].strip())
            #print str(ln.strip())
            #print '---'
            if str(ogid[0].strip())==str(ln.strip()):
                deltaID = int(ogid[1])
        lq=0
    else:
    '''
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
        if ln == 'm_0_left':
            if len(sta) > 0:
                fmdpos_x += [tp.GetPosition([LIDs[0][nid-deltaID]])[0][1]]
                fmdpos_y += [tp.GetPosition([LIDs[0][nid-deltaID]])[0][0]]
                for s in sta:
                    ml_pos_times_list += [(int(s),tp.GetPosition([LIDs[0][nid-deltaID]])[0][1],tp.GetPosition([LIDs[0][nid-deltaID]])[0][0])]
        if ln == 'm_0_right':
            if len(sta) > 0:
                fmdpos_x += [tp.GetPosition([LIDs[0][nid-deltaID]])[0][1]]
                fmdpos_y += [tp.GetPosition([LIDs[0][nid-deltaID]])[0][0]]
                for s in sta:
                    mr_pos_times_list += [(int(s),tp.GetPosition([LIDs[0][nid-deltaID]])[0][1],tp.GetPosition([LIDs[0][nid-deltaID]])[0][0])]
        if ln == 'midgets':
            if len(sta) > 0:
                midpos_x += [tp.GetPosition([LIDs[0][nid-deltaID]])[0][1]]
                midpos_y += [tp.GetPosition([LIDs[0][nid-deltaID]])[0][0]]
                for s in sta:
                    mid_pos_times_list += [(int(s),tp.GetPosition([LIDs[0][nid-deltaID]])[0][1],tp.GetPosition([LIDs[0][nid-deltaID]])[0][0])]
                
        lq+=1
    l_file.close()
    #print lr
            
OGIDs.close()


ml_times_pos_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/midget_ml_times_pos.txt','w+') #_'+str(mdw)+'
for t in ml_pos_times_list:
    ml_times_pos_file.write(str(t)+'\n')
ml_times_pos_file.close()

mr_times_pos_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/midget_mr_times_pos.txt','w+') #_'+str(mdw)+'
for t in mr_pos_times_list:
    mr_times_pos_file.write(str(t)+'\n')
mr_times_pos_file.close()

mid_times_pos_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/midget_times_pos.txt','w+')
for t in mid_pos_times_list:
    mid_times_pos_file.write(str(t)+'\n')
mid_times_pos_file.close()

#---------------------------------------------------------------------------------------------------------------CREATE-CONNECTIONS
#connections to left/right half of Reichardt detector
'''
parasolic_to_psd_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': 40., 'mask' : {'rectangular' : {'lower_left' : [-0.2,-0.2], 'upper_right' : [0.2,0.2]}}}
psd_to_gmdet_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': 7., 'mask' : {'rectangular' : {'lower_left' : [-5.1,-5.1*0.866], 'upper_right' : [4.9,4.9*0.866]}}}
gmdet_to_gmc_conndict = {'connection_type' : 'divergent', 'synapse_model': 'inh_400_delay', 'mask' : {'rectangular' : {'lower_left' : [-4.9,-4.9*0.866], 'upper_right' : [5.1,5.1*0.866]}}} #{'lower_left' : [-1.2,-2.4*0.866], 'upper_right' : [1.3,2.6*0.866]}}}
modet_to_gmc_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex_50', 'mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.1,0.1]}}}

psd_to_gmn_conndict = {'connection_type' : 'convergent', 'synapse_model': 'ex', 'weights': 7.}
gmn_to_gmc_conndict = {'connection_type' : 'divergent', 'synapse_model': 'inh_400_delay'} #{'lower_left' : [-1.2,-2.4*0.866], 'upper_right' : [1.3,2.6*0.866]}}}

out_conndict = {'connection_type' : 'convergent', 'mask' : {'rectangular' : {'lower_left' : [-0.2,-0.2], 'upper_right' : [0.2,0.2]}}}

#connect them
#tp.ConnectLayers(midgets_V_input,midgets,V_input_conndict)
#tp.ConnectLayers(parasolic_V_input,parasolic,V_input_conndict)

tp.ConnectLayers(parasolic,psd,parasolic_to_psd_conndict)
tp.ConnectLayers(psd,gmdet,psd_to_gmdet_conndict)
tp.ConnectLayers(gmdet,out_gmdet,out_conndict)
tp.ConnectLayers(psd,gmn,psd_to_gmn_conndict)
tp.ConnectLayers(gmn,out_gmn,out_conndict)

for i in range(len(layers)):
    tp.ConnectLayers(layers[i],gmc_layers[i],modet_to_gmc_conndict)
    tp.ConnectLayers(gmdet,gmc_layers[i],gmdet_to_gmc_conndict)
    tp.ConnectLayers(gmn,gmc_layers[i],gmn_to_gmc_conndict)
    tp.ConnectLayers(gmc_layers[i],out_layers[i],out_conndict)
'''

#ctr = tp.FindNearestElement(parasolic,[120.,120.])
#fig = tp.PlotLayer(parasolic,nodesize=80)
#tp.PlotTargets(ctr,psd,fig=fig,mask=out_conndict['mask'],src_size=250, tgt_color='red',tgt_size=20)

#fig = tp.PlotLayer(psd)
#ctr = tp.FindNearestElement(parasolic,[120.,120.])
#tp.PlotTargets(ctr, psd ,fig = fig, mask=out_conndict['mask'], src_size=200, src_color = 'red', tgt_color = 'green' , tgt_size=20 )
#plt.show()


'''
#plt.plot(parpos_y,parpos_x,'b+')
fig = plt.figure(figsize=(4,2))
ax1 = fig.add_subplot(1,1,1, adjustable='box', aspect=1)

ax1.plot(fmdpos_y,fmdpos_x,'g.')
#plt.plot(fmdpos,'rx')
ax1.set_xlabel('position x [arcmin]')
ax1.set_ylabel('position y\n[arcmin]')

plt.savefig('img/murakami/spike_pos_on105_off'+str(cond_nr)+'.pgf')
plt.savefig('img/murakami/spike_pos_on105_off'+str(cond_nr)+'.pdf')
#plt.show()
'''

#-------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------SIMULATION
       
#run simulation
#nest.Simulate(t_end-t_start)
  
'''
MIDs = nest.GetNodes(midgets)
for n in range(len(MIDs[0])):
    qr = I_E+np.random.normal(0,10,1)[0]
    nest.SetStatus([MIDs[0][n]], {'I_e': qr})
nest.Simulate(200)
'''

#ps_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/p_corr_max_spikes.txt','a+')
#ps_file.write(str(handle_name)+'\t'+str(exp_nr)+'\t'+str(cond_nr)+'\t'+str(gmn_sum_spikes)+'\t'+str(gmdet_sum_spikes)+'\t'+str(all_spikes)+'\t'+str(modet_sum_spikes[0])+'\t'+str(modet_sum_spikes[1])+'\t'+str(modet_sum_spikes[2])+'\t'+str(modet_sum_spikes[3])+'\t'+str(modet_sum_spikes[4])+'\t'+str(modet_sum_spikes[5])+'\n')
#ps_file.close()

print 'finished'

for lr in layers:
    lr = 0
midgets = 0

sys.exit()