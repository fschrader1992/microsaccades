# -*- coding: utf-8 -*-
import fileinput
import os.path
from scipy.optimize import curve_fit

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# -*- coding: utf-8 -*-
import fileinput
import os.path
import sys
import colorsys
import datetime
from scipy.optimize import curve_fit

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pgf_with_rc_fonts = {"font.family": "serif","font.serif": [],"font.sans-serif": []}
plt.rcParams.update(pgf_with_rc_fonts)

now = datetime.datetime.now()

def lin_func(x, a, b):
    return a *x + b
def e_func(x, a, b, c):
    return a * np.exp(-b * x) + c
def log_func(x, a, b, c):
    return a * np.log(b * x) + c

#for the I_E curve
data = np.loadtxt('data/mo_det_cal/average_spiking_I_e_no_zeros.txt')
#delay
I_e = data[:,0]
STD_I_e = 1000./data[:,1]
IST = 1000./data[:,2]
#IST_err = data[:,3]

fig = plt.figure(1)	
fig.set_size_inches(5,3)
ax = fig.add_subplot(111)
ax.set_xlabel(r'Constant External Input Current $I_e$ [pA]', fontsize=8)
ax.set_ylabel(r'Spiking Frequency $f_S$ [Hz]', fontsize=8)
#ax.plot(I_e,IST,'k--',label='Spiking Frequency') 
ax.plot(I_e,STD_I_e,'kx') 
ax.tick_params(labelsize=8)

#popt, pcov = curve_fit(lin_func, I_e, IST) #, bounds=(0, [3., 2.])
#ax.plot(I_e, lin_func(I_e, *popt), 'k-', label='fit')
plt.tight_layout()

plt.savefig('img/mo_det_cal/I_e_response_curve_small.pdf')
plt.savefig('img/mo_det_cal/I_e_response_curve_small.pgf')
out= 'img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '_I_e_response_curve.pdf'
plt.savefig(out)
plt.show()
plt.close()
#plt.show()

'''
data = np.loadtxt('data/mo_det_cal/m_max_spikes.txt')
#delay
weights = data[:,0]
ws = set(weights)

for w in ws:
    print w
    data = np.loadtxt('data/mo_det_cal/m_max_spikes.txt')
    delay = data[:,1]
    ds = set(delay)
    #set of maximums
    som = []
    dels = []
    
    for s in ds:
        #print 'w '+str(w)
        #print 's '+str(s)
        vel=[]
        sp=[]
        sp_ps=[]
        sp_r=[]
        sp_l=[]
        for l in data:
            if l[0]==w and l[1]==s:
                vel+=[l[2]]
                sp+=[l[3]]
                sp_ps+=[l[5]]
                sp_r+=[l[5]]
                sp_l+=[l[6]]
            #velocity
            #vel = data[weis==w,2]
            #spike number
            #sp = data[weis==w,3]
        #print 'v: '+str(vel)
        #vel=vel[delay==s]
        #sp=sp[delay==s]
        mf= max(sp)
        v=vel[sp.index(mf)]
        #print v
        dels+=[s]
        som+=[v]
        
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.set_title('midgets, delay '+str(s)+', weight '+str(w))
        ax.set_xlabel('velocity in arcmin/s')
        ax.set_ylabel('maximum firingrate in spikes/s')
        ax.plot(vel,sp,'rx')
        ax.plot(vel,sp_r,'bx')
        ax.plot(vel,sp_l,'gx')
        #print vel
        #print sp
        plt.savefig('img/mo_det_cal/midget_resp_d' + str(int(s)) + '_w'+str(int(w))+'.pdf')
        #plt.show()
        plt.close()
        
    #parasols 

    data = np.loadtxt('data/mo_det_cal/p_max_spikes.txt')
    #delay
    delay = data[:,1]
    ds = set(delay)
    #set of maximums
    som_p = []
    dels_p = []
    
    for s in ds:
        #print 'w '+str(w)
        #print 's '+str(s)
        vel=[]
        sp=[]
        sp_ps=[]
        sp_r=[]
        sp_l=[]
        for l in data:
            if l[0]==w and l[1]==s:
                vel+=[l[2]]
                sp+=[l[3]]
                sp_ps+=[l[5]]
                sp_r+=[l[5]]
                sp_l+=[l[6]]
            #velocity
            #vel = data[weis==w,2]
            #spike number
            #sp = data[weis==w,3]
        #print 'v: '+str(vel)
        #vel=vel[delay==s]
        #sp=sp[delay==s]
        mf= max(sp)
        v=vel[sp.index(mf)]
        #print v
        dels_p+=[s]
        som_p+=[v]
'''
        
'''
        #velocity
        vel = data[weis==w,2]
        #spike number
        sp = data[weis==w,3]

        vel=vel[delay == s]
        sp=sp[delay == s]
        mf= max(sp)
        v=vel[sp==mf]
        #print v
        dels_p+=[s]
        som_p+=[v[0]]
'''
'''
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.set_title('parasols, delay '+str(s)+', weight '+str(w))
        ax.set_xlabel('velocity in arcmin/s')
        ax.set_ylabel('maximum firingrate in spikes/s')
        ax.plot(vel,sp,'rx')
        ax.plot(vel,sp_r,'bx')
        ax.plot(vel,sp_l,'gx')
        plt.savefig('img/mo_det_cal/parasol_resp_d' + str(int(s)) + '_w'+str(int(w))+'.pdf')
        plt.close()
        
    fig = plt.figure(1)
    #print som
    #print dels	
    fig.suptitle('maximum responses for delays with weight '+str(w))
    ax = fig.add_subplot(121)
    ax.set_title('midgets')
    ax.set_xlabel('delay in ms')
    ax.set_ylabel('maximum response velocity in arcmin/s')
    ax.plot(dels,som,'rx')


    #fig = plt.figure(1)
    #print som
    #print dels	
    ax = fig.add_subplot(122)
    ax.set_title('parasols')
    ax.set_xlabel('delay in ms')
    #ax.set_ylabel('maximum response velocity in arcmin/s')
    ax.plot(dels_p,som_p,'bx')


    plt.savefig('img/mo_det_cal/responses_w'+str(int(w))+'.pdf')
    plt.close()
    #plt.show()
'''