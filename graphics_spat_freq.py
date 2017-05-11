#!/bin/python
#THIS PART READS THE PIXELS OF A MOVIE AND APPLIES THE TEMPORAL FILTERS TO THE PROPOSED MODEL FOR MICROSACCADES. THE OUTPUT ARE THE POTENTIAL VALUES FOR THE POISSON RATES.
'''
loading the video into an array
-> first get all the pixels in a frame (grayscale)
-> after that calculate the values of each tempral flter at each time and store them in another array which will be the basis for changing poisson rates
'''

import sys
import os
import pylab as pyl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cv2
import datetime
import itertools 
import multiprocessing
from microsaccades_functions import *


now = datetime.datetime.now()
'''
on_tau1 = 5 #.2
on_tau2 = 15 #.1
on_p = .8 #.05
dt=1.
#set values for the spatial filters for the pixel in \mum(?)
alpha = .1
beta = 5 #CSR
sigma = 1 # parasolic/midget cell ratio * sigma (=1)
px_rec_ratio = 3. #pixel to receptor ratio, needed for receptor distance/number
px_dist= 1./px_rec_ratio #set rec distance to one, if changes, also change for parasolic spatial filter
par_m_ratio = 4.

#spatial filter values determined
spat_filter_break_radius = 6 #filter radius in px 

#subsequent values
receptor_dist = 1. #px_rec_ratio*px_dist #for the moment, later maybe hexagonal?

#we actually need to calculate the values of each temporal filter just once
temp_filter_on = tempFilter(40,5.,on_tau1,on_tau2,on_p) 

plt.title('temporal filter with 5ms steps')
plt.xlabel('frame')
plt.ylabel('weight')
plt.plot(temp_filter_on)
plt.savefig('img/temp_filter_5ms.pdf')
plt.show()

'''
#THIS IS THE PART OF THE NEURAL NETWORK. IT TAKES THE VALUES TO CALCULATE THE POISSONRATES OF THE INPUT NEURONS AND THEN CALCULATES THE NETWORK OUTPUT 

import pylab as pyl
import numpy as np
import cv2
import matplotlib.pyplot as plt
from microsaccades_functions import *
from scipy.fftpack import fft


#necessary paramater definitions
framerate = 30
#synaptic weights for the 
syn_weight = 0.1

#----------------------------------------------------------------------------INPUT-RATES-FROM-MS_INPUT
#here needs to be a part that transfers potentials into poisson rates
m1_file = open('data/phases/midget_rates_phases_0fr0deg_24px_only_border.data','r+')
m1_data = np.load(m1_file)   
m1_file.close()

p1_file = open('data/phases/parasolic_rates_phases_0fr0deg_24px_only_border.data','r+')
p1_data = np.load(p1_file)   
p1_file.close()

'''
m2_file = open('data/spat_freq/parasolic_rates_spatfreq_0fr0deg2spat.data','r+')
m2_data = np.load(m2_file)   
m2_file.close()

m3_file = open('data/spat_freq/midget_rates_spatfreq_0fr0deg60spat.data','r+')
m3_data = np.load(m3_file)   
m3_file.close()
'''

out_data = []
for i in range(14,26):
    out_data += [np.sqrt(np.power(m1_data[10,i],2))]
 
d_file = open('data/phase_displacement.data','r+')
disp = np.load(d_file)   
d_file.close()
#plt.plot(np.sqrt(np.power(m1_data[10,i],2)))
#plt.show()

#fig = plt.figure(1)

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')

ax1.set_title('Midget Cells')
ax1.imshow(out_data, interpolation='nearest', cmap='gist_heat', aspect='auto')
ax1.set_ylabel('Relative RF Phase (deg)')

#cax = fig.add_axes([0.,0.,1.,1.])
#cax.patch.set_alpha(0)
#cax.set_frame_on(False)

ax3.plot(disp)
ax3.set_xlabel('Time in ms')
ax3.set_ylabel('Gaze direction (arcmin)')
#for i in range(10,22):
#    plt.plot(m1_data[10,i])
#print sum(m3_data[10,100])
#print sum(m3_data[10,10])
#plt.plot(m2_data[20,20])
#plt.plot(m3_data[20,20])

out_data = []
for i in range(2,10):
    out_data += [np.sqrt(np.power(p1_data[4,i],2))]
 
d_file = open('data/phase_displacement.data','r+')
disp = np.load(d_file)   
d_file.close()
#plt.plot(np.sqrt(np.power(m1_data[10,i],2)))
#plt.show()


#ax3 = fig.add_subplot(222)
ax2.set_title('Parasol Cells')
ax2.imshow(out_data, interpolation='nearest', cmap='gist_heat', aspect='auto')
#cax = fig.add_axes([0.,0.,1.,1.])
#cax.patch.set_alpha(0)
#cax.set_frame_on(False)

#ax4 = fig.add_subplot(224)
ax4.plot(disp)
ax4.set_xlabel('Time in ms')

plt.savefig('img/ana/phases.pdf')
out= 'img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/phases_midget_parasols_' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '.pdf'
plt.savefig(out)
plt.show()

#fourier-analysis

#m3ft=fft(m3_data[5,5])

#plt.plot(m3ft)
#plt.show()
