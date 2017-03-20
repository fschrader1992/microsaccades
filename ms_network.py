#THIS IS THE PART OF THE NEURAL NETWORK. IT TAKES THE VALUES TO CALCULATE THE POISSONRATES OF THE INPUT NEURONS AND THEN CALCULATES THE NETWORK OUTPUT 

import pylab as pyl
import numpy as np
import cv2
import matplotlib.pyplot as plt
import nest
import nest.topology as tp
from microsaccades_functions import *

nest.ResetKernel()   # in case we run the script multiple times from iPython

#-----------------------------------------------------------------------------POISSON+SYNAPTIC-WEIGHTS
#here needs to be a part that transfers potentials into poisson rates
tf_file = open('data/temp_filter_subtr.data','r+')
midget_data = np.load(tf_file)   
#print data

input_rates=poissonRate(midget_data)

#-----------------------------------------------------------------------------------------NETWORK-PART
#---------------------------------------------------------------------------INITIALIZE-POISSON-NEURONS

nest.SetKernelStatus({'resolution': 0.01})

g = nest.Create('sinusoidal_poisson_generator', n=2, params=[{'rate': 10000.0, 
                                                              'amplitude': 5000.0,
                                                              'frequency': 10.0, 
                                                              'phase': 0.0},
                                                             {'rate': 0.0, 
                                                              'amplitude': 10000.0,
                                                              'frequency': 5.0, 
                                                              'phase': 90.0}])

m = nest.Create('multimeter', n=2, params={'interval': 0.1, 'withgid': False,
                                           'record_from': ['rate']})
s = nest.Create('spike_detector', n=2, params={'withgid': False})

nest.Connect(m, g, 'one_to_one')
nest.Connect(g, s, 'one_to_one')
print nest.GetStatus(m)
nest.Simulate(200)

colors = ['b', 'g']

for j in range(2):

    ev = nest.GetStatus([m[j]])[0]['events']
    t = ev['times']
    r = ev['rate']

    sp = nest.GetStatus([s[j]])[0]['events']['times']
    plt.subplot(221)
    h, e = np.histogram(sp, bins=np.arange(0., 201., 5.))
    plt.plot(t, r, color=colors[j])
    plt.step(e[:-1], h * 1000 / 5., color=colors[j], where='post')
    plt.title('PST histogram and firing rates')
    plt.ylabel('Spikes per second')

    plt.subplot(223)
    plt.hist(np.diff(sp), bins=np.arange(0., 1.005, 0.02),
             histtype='step', color=colors[j])
    plt.title('ISI histogram')
    
plt.show()

#------------------------
'''midgets = tp.CreateLayer({'rows': 1, 'columns': 5, 'extent': [5.,1.],'elements': 'poisson_generator','edge_wrap': True})

nest.PrintNetwork(depth=2)'''