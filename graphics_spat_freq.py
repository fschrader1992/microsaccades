#!/bin/python
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
import sys
import os
import glob
import pylab as pyl
import numpy as np
import cv2
import matplotlib.pyplot as plt
from microsaccades_functions import *
from scipy.fftpack import fft
import datetime

now = datetime.datetime.now()
'''

#necessary paramater definitions
framerate = 30
#synaptic weights for the 
syn_weight = 0.1

'''

m1_file = open('data/phases/midget_rates_phases_0fr0deg_36px_only_border.data','r+')
m1_data = np.load(m1_file)   
m1_file.close()

p1_file = open('data/phases/parasolic_rates_phases_0fr0deg_36px_only_border.data','r+')
p1_data = np.load(p1_file)   
p1_file.close()

out_data = []
for i in range(2,38):
    out_data += [np.sqrt(np.power(m1_data[10,i],2))]
 
d_file = open('data/phase_displacement.data','r+')
disp = np.load(d_file)   
d_file.close()
#plt.plot(np.sqrt(np.power(m1_data[10,i],2)))
#plt.show()

#fig = plt.figure(1)

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')

ax1.set_title('Midget Cells')
ax1.imshow(m1_data[10,2:38], interpolation='nearest', cmap='gist_heat', aspect='auto')
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
for i in range(0,12):
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

#plt.plot(p1_data[5,5])
#plt.show()


#here for two

m1_file = open('data/phases/midget_rates_phases_0fr0deg_36px_only_border.data','r+')
m1_data = np.load(m1_file)   
m1_file.close()

p1_file = open('data/phases/parasolic_rates_phases_0fr0deg_36px_only_border.data','r+')
p1_data = np.load(p1_file)   
p1_file.close()

p2_file = open('data/phases/parasolic_rates_phases_0fr0deg_36px_only_border_on.data','r+')
p2_data = np.load(p2_file)   
p2_file.close()

p3_file = open('data/phases/parasolic_rates_phases_0fr0deg_36px_only_border_on_off.data','r+')
p3_data = np.load(p3_file)   
p3_file.close()

out_data = []
for i in range(2,38):
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
for i in range(0,12):
    out_data += [np.sqrt(np.power(p1_data[i,4],2))]
 
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

plt.savefig('img/ana/phases_ms2.pdf')
out= 'img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/phases_midget_parasols_ms2_' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '.pdf'
plt.savefig(out)
plt.show()

'''
m1_file = open('data/spat_freq/parasolic_rates_spatfreq_0fr0deg15spat.data','r+')
m1_data = np.load(m1_file)   
m1_file.close()

m2_file = open('data/spat_freq/parasolic_rates_spatfreq_0fr0deg10spat.data','r+')
m2_data = np.load(m2_file)   
m2_file.close()

m3_file = open('data/spat_freq/midget_rates_spatfreq_0fr0deg15spat.data','r+')
m3_data = np.load(m3_file)   
m3_file.close()


x_val=[1,2,5,10,15,20,30,60]

m_first = []
p_first = []
p_second = []
os.chdir("data/spat_freq")
for file in glob.glob("parasolic*spat_ms2.data"):
    p_file = open(file,'r+')
    p_data = np.load(p_file)   
    p_file.close()
    p_ft=fft(p_data[5,5,50:])
    #plt.plot(np.abs(p_ft))
    #plt.show()
    p_ft=np.abs(p_ft)
    p_first+=[p_ft[7].real]
    p_second+=[p_ft[3].real]

for file in glob.glob("midget*spat_ms2.data"):
    m_file = open(file,'r+')
    m_data = np.load(m_file)   
    m_file.close()
    m_ft=fft(m_data[5,5,50:])
    m_ft=fft(np.abs(m_ft))
    m_first+=[m_ft[3].real]
os.chdir('../..')    
#fourier-analysis

plt.subplot(221)

plt.title('Parasol cells')
plt.plot(x_val,p_first)
plt.plot(x_val,p_second)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('spatial frequency (cyc/deg)')
plt.ylabel('Firing rate (spike/s)')

plt.subplot(222)

plt.title('Midget cells')
plt.plot(x_val,m_first)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('spatial frequency (cyc/deg)')
#plt.ylabel('Firing rate (spike/s)')
#m1ft=fft(m3_data[5,5,50:])
#m2ft=fft(m2_data[5,5,50:])
plt.show()
plt.plot(p_ft)
plt.plot(m_ft)

#plt.plot(np.sqrt(np.power(m2_data[5,5],2)))
#plt.plot(m1_data[5,5])
#plt.plot(m3_data[5,5])#(np.sqrt(np.power(m3_data[5,5],2)))
#plt.plot(m1ft)
#plt.plot(m2ft)
plt.show()
'''