# -*- coding: utf-8 -*-
import fileinput
import os.path
from scipy.optimize import curve_fit

import numpy as np
import matplotlib.pyplot as plt

pgf_with_rc_fonts = {"font.family": "serif","font.serif": [],"font.sans-serif": []}
plt.rcParams.update(pgf_with_rc_fonts)

plt.close('all')

'''
bins = [0. for i in range(200)]
for i in range(50):
    bins[i]=1.-0.02*i
bins[2]=1.
#bins[2]=1.
t=[i for i in range(400)]
'''
'''
def pois_func(f,K1,K2,K3,K4):
    return K1/((K2*f+K3)*(K2*f+K3))+K4 #-0.0000019*(f-7.4119)*(f-7.4119)*(f-7.4119)+0.08911 #-K6*f*f #((1.+2*np.pi*f*0.2)*(1.+2*np.pi*f*0.016))+K2

def f_fit(f):
    return 7.53287/((-0.84356*f+1.53989)*(-0.84356*f+1.53989))+0.0261-0.0000019*(f-7.4119)*(f-7.4119)*(f-7.4119)+0.08911 #-K6*f*f #((1.+2*np.pi*f*0.2)*(1.+2*np.pi*f*0.016))+K2
   
def e_func(x, a, b, c):
    return a * np.exp(-b * x) + c
   
def delog(y):
    return 3.33*np.power(10,(y-10)/20)

   
#x = np.array([5.,12.,20.,40.,75.])   
#y = np.array([20.,0.,-10.,-20.,-30.])
x = np.array([5.,12.,20.,40.,50.,65.,80.,100.,110.]) #,125.,140.,150.])   
y = np.array([delog(20.),delog(0.),delog(-10.),delog(-20.),delog(-23.),delog(-22.),delog(-21.),delog(-25.),delog(-28)]) #,delog(-29.),delog(-30.),delog(-30.)])       

print y
#y=delog(y)
#print y
fig = plt.figure(1)   
fig.set_size_inches(6,4)
ax = fig.add_subplot(111)
ax.set_xlabel('Power Spectrum in dB')
ax.set_ylabel('frequency $f$ in Hz')
ax.plot(x,y,'kx',label='data')

xr=np.arange(0,150,0.2)
#z = np.polyfit(x, y, 2)
p = np.poly1d(np.polynomial.polynomial.polyfit(x, y, 6))
#print p.coefficients
#ax.plot(xr,f_fit(xr),'b',label='fit')
ax.plot(xr,p(xr),'b',label='fit')
'''
'''
freqs = []
for k in range(125):
    freqs += [p(k)]
for k in range(125,250):
    freqs += [0.] #-29.-0.01*k]
fd=[]
for k in range(249):
    fd+=[freqs[k+1]-freqs[k]]
#fd.shuffle()
z = np.fft.ifft(freqs)
'''
'''
#popt, pcov = curve_fit(pois_func, x, y) #, bounds=(-1000000.,[1000000., 1000., 1000., 1000., -1.]))
#ax.plot(xr, pois_func(xr, *popt), 'k-', label='fit')
#_ = ax.plot(x, y, '.', xr, p(xr), '-')
#_ = ax.plot(xr, z, 'k-', label='fit')
ax.set_ylim(0,1)
plt.show()
'''
'''
t = np.arange(1000)
n = np.zeros((1000,), dtype=complex)
for f in range(40):
    n[0:400]=K1/((1.+2*np.pi*f*0.2)*(1.+2*np.pi*f*0.016))
#n[40:60] = np.exp(1j*np.random.uniform(0, 2*np.pi, (20,)))
s = np.fft.ifft(n)
plt.plot(t, s.real, 'b-', t, s.imag, 'r--')

#plt.plot(t,np.fft.ifft(bins).real)
plt.show()
'''

#this returns the result of a normal temporal filter (double exponential)
def tempFilterFct(t,tau1,tau2,p):
    return t*t*t/(tau1*tau1*tau1*tau1)*np.exp(-t/tau1)-p*t*t*t/(tau2*tau2*tau2*tau2)*np.exp(-t/tau2)
   
def spatialFilter(x,x0,sigma,alpha,beta):
    return np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))-alpha*np.exp(-(x-x0)*(x-x0)/(2*beta*beta*sigma*sigma))
def spatialFilterCenter(x,x0,sigma):
    return np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))
def spatialFilterSurround(x,x0,sigma,alpha,beta):
    return alpha*np.exp(-(x-x0)*(x-x0)/(2*beta*beta*sigma*sigma))

fig, ax1 = plt.subplots(1, figsize=(3,3))

xr=np.arange(0,160,0.2)
z1=np.zeros(750)
xr2=np.arange(-10,10,0.05)
z2=np.zeros(400)

on_tau1 = 5. #.2
on_tau2 = 15. #.1
on_p = .8 #.05

#set values for the spatial filters for the pixel in \mum(?)
alpha = .1
beta = 3. #CSR
sigma = 1. # parasolic/midget cell ratio * sigma (=1)
px_midget_ratio = 1. #pixel to receptor ratio, needed for receptor distance/number
px_dist= 0.5 #in arcmin, from paper
par_m_ratio = 4. #2arcmin center field sigma from paper

#spatial filter values determined
spat_filter_break_radius = 6 #filter radius in px

#subsequent values
midget_dist = px_midget_ratio#1, rest gone, since unit is 1px #*px_dist #0.5

#we actually need to calculate the values of each temporal filter just once

ax1.plot(xr,tempFilterFct(xr,on_tau1,on_tau2,on_p),'k')
#ax1.plot(xr,z1,'k--')
ax1.set_xlim(0,160)
ax1.set_ylabel('filter amplitude [arbitrary units]', fontsize=8)
ax1.set_xlabel('time [ms]', fontsize=8)
ax1.tick_params(labelsize=8)

plt.savefig('img/temporal_filter_example.pdf')
plt.savefig('img/temporal_filter_example.pgf')
plt.show()
plt.close()

fig2, ax2 = plt.subplots(1, figsize=(3,3))

ax2.plot(xr2,spatialFilterCenter(xr2,0,sigma),c='gray',label='Center')
ax2.plot(xr2,spatialFilterSurround(xr2,0,sigma,alpha,beta),'--',c='gray',label='Surround')
ax2.plot(xr2,spatialFilter(xr2,0,sigma,alpha,beta),'k',label='Filter')
#ax2.plot(xr2,z2,'k--')
ax2.set_ylabel('filter amplitude [arbitrary units]', fontsize=8)
ax2.set_xlabel('spatial position [arcmin]', fontsize=8)
ax2.legend(loc='upper right', fancybox=False, fontsize=8)

ax2.tick_params(labelsize=8)

plt.savefig('img/spatial_filter_example.pdf')
plt.savefig('img/spatial_filter_example.pgf')
plt.show()
