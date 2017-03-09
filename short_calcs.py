import pylab as pl
import numpy as np
import cv2
from itertools import imap
from microsaccades_functions import *

#calculate influence of temporal filter at different times
off_tau1 = .2
off_tau2 = .1
off_p = .5
dt= 0.033
v = tempFilter(200,dt,off_tau1,off_tau2,off_p)
vmax = 0
for i in range(0,200):
    if v[i]*v[i] > vmax*vmax:
        vmax = v[i]
pl.plot(abs(v/vmax))
pl.show()
for i in range(0,200):        
    print i, abs(v[i]/vmax)
