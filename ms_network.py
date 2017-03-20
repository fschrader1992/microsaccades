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
#print data

midget_rates=poissonRate(input_data)

#-----------------------------------------------------------------------------------------NETWORK-PART
#---------------------------------------------------------------------------INITIALIZE-POISSON-NEURONS

nest.SetKernelStatus({'resolution': 0.01})

'''
#midgets = tp.CreateLayer({'rows': len(midget_rates), 'columns': len(midget_rates[0]), 'extent': [float(len(midget_rates[0])),float(len(midget_rates))],'elements': 'poisson_generator', 'edge_wrap': True})

midgets = tp.CreateLayer({'rows': 2, 'columns': 5, 'extent': [5.,2.],'elements': 'poisson_generator', 'edge_wrap': True})

#tp.PlotLayer(midgets)
#plt.show()

#print midgets

#out = tp.CreateLayer({'rows': len(midget_rates), 'columns': len(midget_rates[0]), 'extent': [float(len(midget_rates[0])),float(len(midget_rates))],'elements': 'spike_detector'})
out = tp.CreateLayer({'rows': 2, 'columns': 5, 'extent': [5.,2.],'elements': 'spike_detector'})

conndict = { 'connection_type' : 'convergent'}

tp.ConnectLayers(midgets,out,conndict)

#learning in progress-->
#connect to spike detectors
#spikes = nest.Create("spike_detector", len(midget_rates)*len(midget_rates[0]))

#nest.Connect(midgets,spikes)

m = tp.GetElement(midgets,[1,2])
print m
#simulate and reset
for i in range(10):#len(midget_rates[0][0])/framerate):
    rate = np.random.rand()*100
    tp.SetElement(midgets, params=[{'rate': 100.0}])
    nest.Simulate(25)

#events = nest.GetStatus(out,'n_events')
nest.raster_plot.from_device(out, hist=True)
pylab.show()


g = nest.Create('poisson_generator', n=2,
                params=[{'rate': 100.0},
                        {'rate': 100.0}])
 
m = nest.Create('multimeter', n=2, params={'interval': 0.1, 'withgid': False,
                                           'record_from': ['rate']})
s = nest.Create('spike_detector', n=2, params={'withgid': False})
 
nest.Connect(m, g, 'one_to_one')
nest.Connect(g, s, 'one_to_one')
print(nest.GetStatus(m))
#nest.Simulate(200)
for i in range(5):
    nest.Simulate(50)
    rate = np.random.rand()*100.
    rate2 = np.random.rand()*100.
    nest.SetStatus(g, [{'rate': rate},{'rate': rate2}])

colors = ['b', 'g']
 
for j in range(2):
 
    ev = nest.GetStatus([m[j]])[0]['events']
    t = ev['times']
    r = ev['rate']
 
    sp = nest.GetStatus([s[j]])[0]['events']['times']
    plt.subplot(221)
    h, e = np.histogram(sp, bins=np.arange(0., 51., 1.))
    plt.plot(t, r, color=colors[j])
    plt.step(e[:-1], h * 1000 / 5., color=colors[j], where='post')
    plt.title('PST histogram and firing rates')
    plt.ylabel('Spikes per second')
 
    plt.subplot(223)
    plt.hist(np.diff(sp), bins=np.arange(0., 1.005, 0.02),
             histtype='step', color=colors[j])
    plt.title('ISI histogram')
plt.show()'''