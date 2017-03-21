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

#----------------------------------------------------------------------------INPUT-RATES-FROM-MS_INPUT
#here needs to be a part that transfers potentials into poisson rates
tf_file = open('data/temp_filter_subtr.data','r+')
input_data = np.load(tf_file)   

midget_rates=poissonRate(input_data)
maxs = []
'''
for i in range(len(midget_rates)):
    for j in range(len(midget_rates[0])):
       maxs.append(max(midget_rates[i][j]))
print max(maxs)
'''
#-----------------------------------------------------------------------------------------NETWORK-PART
#---------------------------------------------------------------------------INITIALIZE-POISSON-NEURONS

nest.SetKernelStatus({'resolution': 0.01})

#set the initial Poisson rate
rate = 1000.


#create all newly needed models
nest.CopyModel("poisson_generator", "var_poisson_generator",{"rate": rate})
#for the output
nest.CopyModel("spike_detector", "my_spike_detector",{"withgid": True, "withtime": True})
#nest.CopyModel('multimeter', 'my_multimeter',{'interval': 0.1, 'withgid': False,'record_from': ['rate']})

#----------------------------------------------------------------------------------------CREATE-LAYERS

rows = len(midget_rates)
cols = len(midget_rates[0])
print rows, cols
#rows = 40
#cols = 40

#retink that model
nest.SetDefaults('iaf_psc_alpha',{'I_e' : 375.0})
midgets = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'poisson_generator', 'edge_wrap': True})
parasolic = tp.CreateLayer({'rows': int(float(rows)/4), 'columns': int(float(cols)/4), 'extent': [float(cols),float(rows)],'elements': 'iaf_psc_alpha', 'edge_wrap': True})

#parasolic = tp.CreateLayer({'rows': int(float(rows)/4), 'columns': int(float(cols)/4), 'extent': [float(cols),float(rows)],'elements': 'iaf_psc_alpha', 'edge_wrap': True})

#think about distances etc
#tp.PlotLayer(midgets)
#tp.PlotLayer(parasolic)
#plt.show()

out = tp.CreateLayer({'rows': int(float(rows)/4), 'columns': int(float(cols)/4), 'extent': [float(cols),float(rows)],'elements': 'my_spike_detector'})
#out = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'my_spike_detector'})
#out_multi = tp.CreateLayer({'rows': rows, 'columns': cols, 'extent': [float(cols),float(rows)],'elements': 'my_multimeter'})

#----------------------------------------------------------------------------------CREATE-CONNECTIONS

mi_par_on_conndict = {'connection_type' : 'convergent','mask' : {'circular' : {'radius' : 2.}}, 'weights' : {'gaussian' : {'p_center' : 1., 'sigma' : 2.}}}
mi_par_off_conndict = {'connection_type' : 'convergent','mask' : {'circular' : {'radius' : 2.}}, 'weights' : {'gaussian' : {'p_center' : -0.5, 'sigma' : 2.}}}
out_conndict = {'connection_type' : 'convergent'}

tp.ConnectLayers(midgets,parasolic,mi_par_on_conndict)
tp.ConnectLayers(midgets,parasolic,mi_par_off_conndict)
#tp.ConnectLayers(midgets,out,out_conndict)
tp.ConnectLayers(parasolic,out,out_conndict)

#simulate and reset
for f in range(len(midget_rates[0][0])):
    #reset rates
    for row in range(rows):
        for col in range(cols):
            nest.SetStatus(tp.GetElement(midgets,[row,col]), {'rate': midget_rates[row][col][f]})
    #run simulation
    print f
    nest.Simulate(4)#int(1/framerate))

s = tp.GetElement(out,[0,0])
#mult = tp.GetElement(out_multi,[0,0])

dSD = nest.GetStatus(s,keys="events")[0]
evs = dSD["senders"]
ts = dSD["times"]
pyl.figure(2)
pyl.plot(ts, evs, ".")
pyl.show()


'''
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