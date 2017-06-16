#THIS IS THE PART OF THE NEURAL NETWORK. IT TAKES THE VALUES TO CALCULATE THE POISSONRATES OF THE INPUT NEURONS AND THEN CALCULATES THE NETWORK OUTPUT 

import pylab as pyl
import numpy as np
import cv2
import matplotlib.pyplot as plt
import nest
import nest.topology as tp
from microsaccades_functions import *

nest.ResetKernel()   # in case we run the script multiple times from iPython

#necessary paramater definitions
framerate = 30
#synaptic weights for the 
syn_weight = 0.1

#----------------------------------------------------------------------------INPUT-RATES-FROM-MS_INPUT
#here needs to be a part that transfers potentials into poisson rates
#m_file = open('data/midget_values.data','r+')
m_file = open('data/phases/midget_rates_phases_0fr0deg_36px_only_border.data','r+')
m_data = np.load(m_file)   
m_file.close()

#p_file = open('data/parasolic_values.data','r+')
p_file = open('data/phases/parasolic_rates_phases_0fr0deg_36px_only_border.data','r+')
p_data = np.load(p_file)  
p_file.close()

midget_rates=poissonRateMidgets(m_data)
parasolic_rates=poissonRateParasols(p_data)
'''
print midget_rates

#to check for maximum spike rates in order to adopt conversion of film input
maxs = []
for i in range(len(midget_rates)):
    for j in range(len(midget_rates[0])):
       maxs.append(max(midget_rates[i][j]))
print max(maxs)


'''

#-----------------------------------------------------------------------------------------NETWORK-PART
#---------------------------------------------------------------------------INITIALIZE-POISSON-NEURONS

nest.SetKernelStatus({'resolution': 0.01})

#set the initial Poisson rate
rate = 100.

#create all newly needed models
nest.CopyModel("poisson_generator", "var_poisson_generator",{"rate": rate})
#for the output
nest.CopyModel("spike_detector", "my_spike_detector",{"withgid": True, "withtime": True})
#nest.CopyModel('multimeter', 'my_multimeter',{'interval': 0.1, 'withgid': False,'record_from': ['rate']})

#nest.CopyModel("static_synapse", "ex", {"weight" : 0.1})
#nest.CopyModel("static_synapse", "inh", {"weight" : -0.1})

#----------------------------------------------------------------------------------------CREATE-LAYERS

rows = len(midget_rates)
cols = len(midget_rates[0])
print rows, cols
#rows = 40
#cols = 40

#get grid data form previous simulation
gm_file = open('data/test_vid/m_pos_test_vid.data','r+')
gm_data = np.load(gm_file)  
gm_file.close()
gm_data = gm_data.tolist()

gp_file = open('data/test_vid/p_pos_test_vid.data','r+')
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
for i in range(len(gm_data)):
    for j in range(len(gm_data[0])):
        gm_pos+=[[0.5*gm_data[i][j][0]+0.25,0.5*gm_data[i][j][1]]+0.25]
        gm_r_0_pos+=[[0.5*gm_data[i][j][0]+0.5,0.5*gm_data[i][j][1]]+0.25]
        gm_r_60_pos+=[[0.5*gm_data[i][j][0]+0.5,0.5*gm_data[i][j][1]]+0.25+0.433]
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

#print gm_pos
print 'minimum midgets'
print gm_min, gm_max[0]+gm_min[0]+10.
print 'maximum midgets'
print gm_max, gm_max[1]+gm_min[1]+10.
        
gp_pos=[]
gp_r_0_pos=[]
gp_r_60_pos=[]
#gm_r_120_pos=[]

for i in range(len(gp_data)):
    for j in range(len(gp_data[0])):
        gp_pos+=[[0.5*gp_data[i][j][0]+0.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_0_pos+=[[0.5*gp_data[i][j][0]+2.25,0.5*gp_data[i][j][1]+0.25]]
        gp_r_60_pos+=[[0.5*gp_data[i][j][0]+1.25,0.5*gp_data[i][j][1]+0.25+1.732]]
        #gp_r_120_pos+=[[gp_data[i][j][0]-1.25,gp_data[i][j][1]+0.25+1.732]]


#gm_pos=[[np.random.uniform(-0.5,0.5), np.random.uniform(-0.5,0.5)] for j in range(50)]
#print gm_pos
#retink that model
nest.SetDefaults('iaf_psc_alpha',{'I_e' : 374.0})
midgets = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'poisson_generator', 'edge_wrap': True})
parasolic = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'poisson_generator', 'edge_wrap': True})

#120 degree commented for the moment since two build a basis

#these are direction dependent
m_reichardt_0_left = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
m_reichardt_0_right = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
m_reichardt_60_up = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
m_reichardt_60_down = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#m_reichardt_120_up = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#m_reichardt_120_down = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_0_left = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_0_right = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_60_up = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_60_down = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#p_reichardt_120_up = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#p_reichardt_120_down = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})

#react to motion in both directions
m_reichardt_0 = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
m_reichardt_60 = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#m_reichardt_120 = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_0 = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
p_reichardt_60 = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
#p_reichardt_120 = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})

'''
motion_left = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
motion_right = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
motion_up = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
motion_down = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
'''
#think about distances etc
#tp.PlotLayer(midgets)
#tp.PlotLayer(parasolic)
#tp.PlotLayer(reichardt_left)
#plt.show()

out_m = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gm_pos, 'elements': 'my_spike_detector'})
out_p = tp.CreateLayer({'extent' : [121.,121.], 'center' : [60.5,60.5], 'positions' : gp_pos, 'elements': 'my_spike_detector'})
#out = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'my_spike_detector'})
#out_multi = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'my_multimeter'})

#-----------------------------------------------------------------------------------CREATE-CONNECTIONS
#connections to left/right half of Reichardt detector
m_r_0_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.6,-0.1], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}
m_r_0_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.1,0.6]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}
m_r_60_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,-0.5], 'upper_right' : [0.5,0.5]}, 'anchor' : [0.25,0.433]}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}
m_r_60_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,-0.5], 'upper_right' : [0.5,0.5]}, 'anchor' : [-0.25,-0.433]}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}
p_r_0_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.1,-0.1], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}
p_r_0_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.1,-0.1], 'upper_right' : [0.1,2.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}
p_r_60_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.,-2.], 'upper_right' : [2.,2.]}, 'anchor' : [1.,1.732]}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}
p_r_60_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.,-2.], 'upper_right' : [2.,2.]}, 'anchor' : [-1.,-1.732]}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}

#connections of left/right half to just motion sensitive Reichardt detector
m_r_0_left_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,-0.5], 'upper_right' : [0.5,0.5]}} } #'synapse_model' : 'ex'}
m_r_0_right_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,-0.5], 'upper_right' : [0.5,0.5]}} } #'synapse_model' : 'ex'}
m_r_60_up_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,-0.5], 'upper_right' : [0.5,0.5]}} } #'synapse_model' : 'ex'}
m_r_60_down_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,-0.5], 'upper_right' : [0.5,0.5]}} } #'synapse_model' : 'ex'}
p_r_0_left_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-1.25,-2.], 'upper_right' : [1.25,2.]}} } #'synapse_model' : 'ex'}
p_r_0_right_r_0_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-1.25,-2.], 'upper_right' : [1.25,2.]}} } #'synapse_model' : 'ex'}
p_r_60_up_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-1.25,-2.], 'upper_right' : [1.25,2.]}} } #'synapse_model' : 'ex'}
p_r_60_down_r_60_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-1.25,-2.], 'upper_right' : [1.25,2.]}} } #'synapse_model' : 'ex'}

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

out_conndict = {'connection_type' : 'convergent'}

tp.ConnectLayers(midgets,reichardt_left,r_left_conndict)
tp.ConnectLayers(midgets,reichardt_right,r_right_conndict)
tp.ConnectLayers(midgets,reichardt_up,r_up_conndict)
tp.ConnectLayers(midgets,reichardt_down,r_down_conndict)
tp.ConnectLayers(reichardt_left,reichardt_horizontal,r_left_hor_conndict)
tp.ConnectLayers(reichardt_right,reichardt_horizontal,r_right_hor_conndict)
tp.ConnectLayers(reichardt_up,reichardt_vertical,r_up_ver_conndict)
tp.ConnectLayers(reichardt_down,reichardt_vertical,r_down_ver_conndict)
tp.ConnectLayers(reichardt_horizontal,motion_left,r_hor_motion_left_conndict)
tp.ConnectLayers(reichardt_horizontal,motion_right,r_hor_motion_right_conndict)
tp.ConnectLayers(reichardt_vertical,motion_up,r_ver_motion_up_conndict)
tp.ConnectLayers(reichardt_vertical,motion_down,r_ver_motion_down_conndict)
tp.ConnectLayers(parasolic,motion_left,par_motion_left_conndict)
tp.ConnectLayers(parasolic,motion_right,par_motion_right_conndict)
tp.ConnectLayers(parasolic,motion_up,par_motion_up_conndict)
tp.ConnectLayers(parasolic,motion_down,par_motion_down_conndict)
tp.ConnectLayers(midgets,m_out,out_conndict)
tp.ConnectLayers(parasolic,p_out,out_conndict)
#tp.ConnectLayers(motion_left,out,out_conndict)
'''
#-------------------------------------------------------------------------------------------SIMULATION
for f in range(3):#len(midget_rates[0][0])):
    #reset rates
    for row in range(rows):
        for col in range(cols):
            nest.SetStatus(tp.GetElement(midgets,[row,col]), {'rate': midget_rates[row][col][f]})
            nest.SetStatus(tp.GetElement(parasolic,[int(float(row)/4.),int(float(col)/4.)]), {'rate': parasolic_rates[int(float(row)/4.)][int(float(col)/4.)][f]})
    #run simulation
    print f
    nest.Simulate(2)#int(1/framerate))

s = tp.GetElement(out,[0,30])
#mult = tp.GetElement(out_multi,[0,0])

dSD = nest.GetStatus(s,keys="events")[0]
evs = dSD["senders"]
ts = dSD["times"]
pyl.figure(2)
pyl.plot(ts, evs, ".")
pyl.show()



#----------------------------------------------------------------------------------PRINTING-THE-OUTPUT
s = tp.GetElement(out,[0,0])
mult = tp.GetElement(out_multi,[0,0])

dSD = nest.GetStatus(s,keys="events")[0]
evs = dSD["senders"]
ts = dSD["times"]
pyl.figure(2)
pyl.plot(ts, evs, ".")
pyl.show()

ev = nest.GetStatus(mult)[0]['events']
t = ev['times']
r = ev['rate']

sp = nest.GetStatus(s)[0]['events']['times']
plt.subplot(221)
h, e = np.histogram(sp, bins=np.arange(0., 1001., 1.))
plt.plot(t, r, color='b')
plt.step(e[:-1], h , color='b', where='post')
plt.title('PST histogram and firing rates')
plt.ylabel('Spikes per second')

plt.subplot(223)
plt.hist(np.diff(sp), bins=np.arange(0., 1.005, 0.02),
            histtype='step', color='b')
plt.title('ISI histogram')
plt.show()
'''