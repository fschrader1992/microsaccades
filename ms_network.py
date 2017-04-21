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
m_file = open('data/midget_values.data','r+')
m_data = np.load(m_file)   
m_file.close()

p_file = open('data/parasolic_values.data','r+')
p_data = np.load(p_file)  
p_file.close()

midget_rates=poissonRateMidgets(m_data)
parasolic_rates=poissonRateParasols(p_data)

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

#retink that model
nest.SetDefaults('iaf_psc_alpha',{'I_e' : 374.0})
midgets = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2,float(cols)/2], 'elements': 'poisson_generator', 'edge_wrap': True})
parasolic = tp.CreateLayer({'rows': int(float(rows)/4), 'columns': int(float(cols)/4), 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2,float(cols)/2], 'elements': 'poisson_generator', 'edge_wrap': True})

reichardt_left = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2,float(cols)/2], 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
reichardt_right = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2,float(cols)/2], 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
reichardt_up = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2,float(cols)/2], 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
reichardt_down = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2,float(cols)/2], 'elements': 'iaf_psc_alpha', 'edge_wrap': True})

reichardt_horizontal = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2,float(cols)/2 + 0.5], 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
reichardt_vertical = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2 + 0.5,float(cols)/2], 'elements': 'iaf_psc_alpha', 'edge_wrap': True})

motion_left = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2,float(cols)/2 - 1.5], 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
motion_right = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2,float(cols)/2 + 1.5], 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
motion_up = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2 + 1.5,float(cols)/2], 'elements': 'iaf_psc_alpha', 'edge_wrap': True})
motion_down = tp.CreateLayer({'rows': int(float(rows)/4), 'columns': int(float(cols)/4), 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2 - 1.5,float(cols)/2], 'elements': 'iaf_psc_alpha', 'edge_wrap': True})

#think about distances etc
#tp.PlotLayer(midgets)
#tp.PlotLayer(parasolic)
#plt.show()

out = tp.CreateLayer({'rows': int(float(rows)/4), 'columns': int(float(cols)/4), 'extent': [float(cols),float(rows)], 'center' : [float(rows)/2,float(cols)/2], 'elements': 'my_spike_detector'})
#out = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'my_spike_detector'})
#out_multi = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'my_multimeter'})

#-----------------------------------------------------------------------------------CREATE-CONNECTIONS

r_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,-1.], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}
r_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [0.1,1.]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}
r_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [1.,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}
r_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-1.,0.], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}}

r_left_hor_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [0.1,0.5]}}, } #'synapse_model' : 'ex'}
r_right_hor_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,-0.5], 'upper_right' : [0.1,0.1]}}, } #'synapse_model' : 'ex'}
r_up_ver_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-0.5,0.], 'upper_right' : [0.1,0.1]}}, } #'synapse_model' : 'ex'}
r_down_ver_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [0.5,0.1]}}, } #'synapse_model' : 'ex'}

r_hor_motion_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [4.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'ex'}
r_hor_motion_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-4.,0.], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'ex'}
r_ver_motion_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,-4.], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'ex'}
r_ver_motion_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [0.1,4.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'ex'}

#what about delays here?
par_motion_left_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [0.1,2.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'inh'}
par_motion_right_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,-2.], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'inh'}
par_motion_up_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [-2.,0.], 'upper_right' : [0.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'inh'}
par_motion_down_conndict = {'connection_type' : 'convergent','mask' : {'rectangular' : {'lower_left' : [0.,0.], 'upper_right' : [2.1,0.1]}}, 'delays' : {'linear' : {'c' : .1, 'a' : .02}}, } #'synapse_model' : 'inh'}

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
tp.ConnectLayers(midgets,out,out_conndict)
#tp.ConnectLayers(motion_left,out,out_conndict)

#-------------------------------------------------------------------------------------------SIMULATION
for f in range(3):#len(midget_rates[0][0])):
    #reset rates
    for row in range(rows):
        for col in range(cols):
            nest.SetStatus(tp.GetElement(midgets,[row,col]), {'rate': midget_rates[row][col][f]})
            nest.SetStatus(tp.GetElement(parasolic,[int(float(row)/4.),int(float(col)/4.)]), {'rate': parasolic_rates[int(float(row)/4.)][int(float(col)/4.)][f]})
    #run simulation
    print f
    nest.Simulate(10)#int(1/framerate))

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