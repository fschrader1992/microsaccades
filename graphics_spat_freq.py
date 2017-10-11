#!/bin/python
from __future__ import unicode_literals
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
import datetime

pgf_with_rc_fonts = {"font.family": "serif","font.serif": [],"font.sans-serif": []}
plt.rcParams.update(pgf_with_rc_fonts)

now = datetime.datetime.now()



m1_file = open('data/phases/phases_0fr0deg_32px_only_border/midget_rates_phases_0fr0deg_32px_only_border_on.data','r+')
m1_data = np.load(m1_file)   
m1_file.close()

p1_file = open('data/phases/phases_0fr0deg_32px_only_border/parasolic_rates_phases_0fr0deg_32px_only_border_on.data','r+')
p1_data = np.load(p1_file)   
p1_file.close()
print(p1_data)

#out_data = []
#for i in range(2,38):
#    out_data += [np.sqrt(np.power(m1_data[10,i],2))]
 
d_file = open('data/phase_displacement.data','r+')
disp = np.load(d_file)   
d_file.close()
#plt.plot(np.sqrt(np.power(m1_data[10,i],2)))
#plt.show()


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', gridspec_kw = {'height_ratios':[4, 1]}, figsize=(6,3.5))

ax1.set_title('Midget Cells')
ax1.imshow(m1_data[17,32:65], interpolation='nearest', cmap='gist_heat', aspect='auto', extent=[0,1000,0,360])
#ax1.set_aspect(2)
ax1.set_ylabel('Relative RF Phase (deg)')
ax1.set_xlim(0,1000)
ax1.set_yticks(np.arange(0,361,60))

#cax = fig.add_axes([0.,0.,1.,1.])
#cax.patch.set_alpha(0)
#cax.set_frame_on(False)

ax3.plot(disp)
ax3.set_xlabel('Time in ms')
ax3.set_ylabel('Gaze direction\n (arcmin)')
ax3.set_yticks(np.arange(-15, 16, 15))
#for i in range(10,22):
#    plt.plot(m1_data[10,i])
#print sum(m3_data[10,100])
#print sum(m3_data[10,10])
#plt.plot(m2_data[20,20])
#plt.plot(m3_data[20,20])

out_data = []
for i in range(8,17):
    out_data += [np.sqrt(np.power(p1_data[4,i],2))]
    #out_data += [np.sqrt(np.power(p1_data[5,i],2))]
 
d_file = open('data/phase_displacement.data','r+')
disp = np.load(d_file)   
d_file.close()
#plt.plot(np.sqrt(np.power(m1_data[10,i],2)))
#plt.show()

p_data = []
for i in range(8,17):
    #p_data+=[p1_data[4,i]]
    p_data+=[p1_data[5,i]]
#ax3 = fig.add_subplot(222)
ax2.set_title('Parasol Cells')
ax2.set_xlim(0,1000)
ax2.imshow(p_data, interpolation='nearest', cmap='gist_heat', aspect='auto', extent=[0,1000,0,360])
ax2.set_yticks(np.arange(0,361,60))
#cax = fig.add_axes([0.,0.,1.,1.])
#cax.patch.set_alpha(0)
#cax.set_frame_on(False)

#ax4 = fig.add_subplot(224)
ax4.plot(disp)
ax4.set_xlabel('Time in ms')
ax4.set_yticks(np.arange(-15, 16, 15))

plt.savefig('img/ana/phase.pgf')
plt.savefig('img/ana/phase.pdf')
out= 'img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/phase_' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '.pdf'
plt.savefig(out)
plt.show()

#plt.plot(p1_data[5,5])
#plt.show()





#End of PHASES
'''
#here for two

m1_file = open('data/phases/midget_rates_phases_0fr0deg_36px_only_border_on.data','r+')
m1_data = np.load(m1_file)   
m1_file.close()

p1_file = open('data/phases/parasolic_rates_phases_0fr0deg_36px_only_border.data','r+')
p1_data = np.load(p1_file)   
p1_file.close()

m2_file = open('data/phases/parasolic_rates_phases_0fr0deg_36px_only_border_on.data','r+')
m2_data = np.load(m2_file)   
m2_file.close()

p2_file = open('data/phases/parasolic_rates_phases_0fr0deg_36px_only_border_on_off.data','r+')
p2_data = np.load(p2_file)   
p2_file.close()

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
#plt.show()

'''
#FREQUENCY-ANALYSIS------------------------------------------------------------------------------------------------------------
'''
#x_val=[1,2,5,10,15,20,30,60]
x_val=[1]


m_first = []
m_second = []
p_first = []
p_second = []
os.chdir("data/spat_freq")
for x in x_val:
    for file in glob.glob("para*deg"+str(x)+"spat_on.data"):
        p_file = open(file,'r+')
        p_data = np.load(p_file)   
        p_file.close()
        #print(file)
        plt.plot(p_data[16,16,50:],label='para'+str(x))
        print('p_length: '+str(len(p_data)))
        print('p_height: '+str(len(p_data[1])))
        #print('data' + str(np.max(p_data[15,15,210:260])))
        #p_first+=[np.max(p_data[20,20,210:260])]
        #p_second+=[np.min(p_data[20,20,110:160])]

    for file in glob.glob("midget*deg"+str(x)+"spat.data"):
        m_file = open(file,'r+')
        m_data = np.load(m_file)   
        m_file.close()
        #plt.plot(m_data[20,102,:],label='midg'+str(x))
        print('m_length: '+str(len(m_data)))
        print('m_height: '+str(len(m_data[1])))
        #m_first+=[np.max(m_data[15,15,210:260])]
        #m_second+=[np.min(m_data[15,15,110:160])]
os.chdir('../..')  

plt.ylim(bottom=0)
plt.legend()
plt.show()

for f in range(400):

    fig = plt.figure(1)

    ax = fig.add_subplot(221)
    ax.set_title('midget output')
    plt.imshow(m_data[:,:,f], aspect='auto', interpolation='nearest')
    ax.set_aspect('equal')
    plt.axis('off')

    cax = fig.add_axes([0.,0.,1.,1.])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    #plt.colorbar(orientation='vertical')

    ax = fig.add_subplot(223)
    ax.set_title('parasolic output')
    plt.imshow(p_data[:,:,f], aspect='auto', interpolation='nearest')
    ax.set_aspect('equal')
    plt.axis('off')

    cax = fig.add_axes([0.,0.,1.,.5])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.savefig('video/animation/spat_freq/1/an'+str(f).zfill(3)+'.png')
    plt.close()
#plt.show()
'''
#FREQUENCY-ANALYSIS------------------------------------------------------------------------------------------------------------
'''
#x_val=[1,2,5,10,15,20,30,60]
x_val=[10]


m_first = []
m_second = []
p_first = []
p_second = []
os.chdir("data/mo_det_cal")
for x in x_val:
    for file in glob.glob("para*cal_"+str(x)+"fr.data"):
        p_file = open(file,'r+')
        p_data = np.load(p_file)   
        p_file.close()
        #print(file)
        #plt.plot(p_data[3,15,50:],label='para'+str(x))
        print('p_length: '+str(len(p_data)))
        print('p_height: '+str(len(p_data[1])))
        #print('data' + str(np.max(p_data[15,15,210:260])))
        #p_first+=[np.max(p_data[20,20,210:260])]
        #p_second+=[np.min(p_data[20,20,110:160])]

    for file in glob.glob("midget*cal_"+str(x)+"fr.data"):
        m_file = open(file,'r+')
        m_data = np.load(m_file)   
        m_file.close()
        plt.plot(m_data[10,70,:],label='midg'+str(x))
        print('m_length: '+str(len(m_data)))
        print('m_height: '+str(len(m_data[1])))
        #m_first+=[np.max(m_data[15,15,210:260])]
        #m_second+=[np.min(m_data[15,15,110:160])]
os.chdir('../..')  

#plt.ylim(bottom=0)
plt.legend()
plt.show()

for f in range(300):

    fig = plt.figure(1)

    ax = fig.add_subplot(221)
    ax.set_title('midget output')
    plt.imshow(m_data[:,:,f], aspect='auto', interpolation='nearest')
    ax.set_aspect('equal')
    plt.axis('off')

    cax = fig.add_axes([0.,0.,1.,1.])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    #plt.colorbar(orientation='vertical')

    ax = fig.add_subplot(223)
    ax.set_title('parasolic output')
    plt.imshow(p_data[:,:,f], aspect='auto', interpolation='nearest')
    ax.set_aspect('equal')
    plt.axis('off')

    cax = fig.add_axes([0.,0.,1.,.5])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.savefig('video/animation/spat_freq/mo_det_cal_10/an'+str(f).zfill(3)+'.png')
    plt.close()
'''
'''
f, (ax1, ax3) = plt.subplots(1, 2, sharex='col')

ax1.set_title('Midget Cells')
ax1.plot(x_val,m_first,'rx')
ax1.plot(x_val,m_second,'b+')
ax1.set_ylabel('Response')

#cax = fig.add_axes([0.,0.,1.,1.])
#cax.patch.set_alpha(0)
#cax.set_frame_on(False)

ax3.set_title('Parasolic Cells')
ax3.plot(x_val,p_first,'rx')
ax3.plot(x_val,p_second,'b+')
ax3.set_xlabel('spatial frequency (cyc/deg)')

plt.savefig('img/ana/spatfreq.pdf')
out= 'img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/spatfreq_midget_parasols_' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '.pdf'
plt.savefig(out)

plt.show()
'''
'''
m1_file = open('data/spat_freq/midget_rates_spatfreq_0fr0deg10spat.data','r+')
m1_data = np.load(m1_file)   
m1_file.close()

m2_file = open('data/spat_freq/parasolic_rates_spatfreq_0fr0deg10spat_on_off.data','r+')
m2_data = np.load(m2_file)   
m2_file.close()

m4_file = open('data/spat_freq/parasolic_rates_spatfreq_0fr0deg10spat_on.data','r+')
m4_data = np.load(m4_file)   
m4_file.close()

m3_file = open('data/spat_freq/midget_rates_spatfreq_0fr0deg15spat.data','r+')
m3_data = np.load(m3_file)   
m3_file.close()
'''


'''
x_val=[1,2,10,20,30,60]
m_x=[30,30,30,33,30,31]
p_x=[7,7,7,8,7,8]

m_first = []
p_first = []
p_second = []
m_second = []
os.chdir("data/spat_freq")

i=0
#for file in glob.glob("midget*spat_ms2.data"):
for x in x_val:
    m_file = open('spatfreq_0fr0deg'+str(x)+'spat/midget_rates_'+str(x)+'spat.data','r+')
    m_data = np.load(m_file)   
    m_file.close()
    m_ft=np.fft.fft(m_data[5,m_x[i],200:400])
    #plt.plot(m_data[5,m_x[i],200:400])
    #plt.show()
    m_ft=np.abs(m_ft)
    m_first+=[m_ft[1].real]
    m_second+=[m_ft[2].real]
    i+=1

f_max=max(m_first)
print m_first
for m in range(len(m_first)):
    m_first[m]=m_first[m]/f_max*60.
i=0
for x in x_val:
    p_file = open('spatfreq_0fr0deg'+str(x)+'spat/parasolic_rates_'+str(x)+'spat_on.data','r+')
    p_data = np.load(p_file)   
    p_file.close()
    if i==2 or i==4 or i==6:
        p_ft=np.fft.fft(p_data[5,p_x[i],200:400])
    else:
        p_ft=np.fft.fft(p_data[6,p_x[i],200:400])
    #plt.plot(np.abs(p_ft))
    #plt.show()
    p_ft=np.abs(p_ft)
    p_first+=[p_ft[1].real]
    p_second+=[p_ft[2].real]
    i+=1
 
f_max=max(p_first)
for p in range(len(p_first)):
    p_first[p]=p_first[p]/f_max*60.
for p in range(len(p_second)):
    p_second[p]=p_second[p]/f_max*60. 
 
os.chdir('../..')
    
#plot output

f, (ax1, ax2) = plt.subplots(2, 1, sharex='col', figsize=(3,5))

ax1.set_title('Midget Cells')
ax1.plot(x_val,m_first,'k--',label='First Harmonic')
ax1.plot(x_val,m_first,'kx')
#ax1.plot(x_val,m_second)
ax1.set_ylim(1,100)
ax1.set_xscale('log')
ax1.set_yscale('log')
#ax1.set_xlabel('spatial frequency (cyc/deg)')

ax2.set_title('Parasol Cells')
ax2.plot(x_val,p_first,'k--',label='First Harmonic')
ax2.plot(x_val,p_first,'kx')
ax2.plot(x_val,p_second,'--',color='darkgray',label='Second Harmonic')
ax2.plot(x_val,p_second,'x',color='darkgray')
ax2.set_ylim(1,100)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('spatial frequency (cyc/deg)')
ax2.legend(loc='upper center', bbox_to_anchor=(0.5,1.6), fancybox=False, fontsize=8)

f.text(0.0,0.5, 'Firing Rate (Spikes/s)', va='center', rotation='vertical')

plt.tight_layout()

plt.savefig('img/ana/spat_freq.pdf')
plt.savefig('img/ana/spat_freq.pgf')
out= 'img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/spat_freq_' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '.pdf'
plt.savefig(out)

plt.show()
plt.close()
'''