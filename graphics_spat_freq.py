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
m1_file = open('data/spat_freq/midget_rates_spatfreq_0fr0deg1spat.data','r+')
m1_data = np.load(m1_file)   
m1_file.close()

m2_file = open('data/spat_freq/midget_rates_spatfreq_0fr0deg10spat.data','r+')
m2_data = np.load(m2_file)   
m2_file.close()

m3_file = open('data/spat_freq/midget_rates_spatfreq_0fr0deg60spat.data','r+')
m3_data = np.load(m3_file)   
m3_file.close()

plt.plot(m3_data[10,100,:20])
plt.plot(m3_data[10,10,:20])
print sum(m3_data[10,100,:20])
print sum(m3_data[10,10,:20])
#plt.plot(m2_data[20,20])
#plt.plot(m3_data[20,20])
plt.show()

#fourier-analysis

#m3ft=fft(m3_data[5,5])

#plt.plot(m3ft)
#plt.show()
'''