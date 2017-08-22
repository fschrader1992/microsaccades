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
from scipy.optimize import curve_fit

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt('m_max.data')
#delay
delay = data[:,0]
ds = set(delay)
#set of maximums
som = []
dels = []

for s in ds:
	#velocity
	vel = data[:,1]
	#spike number
	sp = data[:,2]

	vel=vel[delay == s]
	sp=sp[delay == s]
	mf= max(sp)
	v=vel[sp==mf]
	print v
	dels+=[s]
	som+=[v[0]]

	'''fig = plt.figure(1)

	ax = fig.add_subplot(221)
	ax.set_title('delay '+str(s))
	ax.set_xlabel('velocity in arcmin/s')
	ax.set_ylabel('maximum firingrate in spikes/s')
	ax.plot(vel,sp, c='b')

	#plt.savefig('.pdf')

	plt.show()
	'''
fig = plt.figure(1)
print som
print dels	
ax = fig.add_subplot(221)
ax.set_title('maximum responses for delays')
ax.set_xlabel('delay')
ax.set_ylabel('maximum response velocity in arcmin/s')
ax.plot(dels,som,c='g')

#plt.savefig('.pdf')

plt.show()