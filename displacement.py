#this creates the displacement on the grating
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from microsaccades_functions import *

#displacement is zero and standard derivation is receptor size (0.5 arcmin)
mu, sigma = 0.,1

dis = np.random.normal(mu, sigma)

x = np.arange(0, 101,1.)
y = np.random.normal(mu, sigma, 101)

tck = interpolate.splrep(x, y, s=0)
xnew = np.arange(0, 100, .1)
ynew = interpolate.splev(xnew, tck, der=0)

ms_disp = -12*np.ones(300)
ms_disp = np.append(ms_disp,30*np.ones(300))
ms_disp = np.append(ms_disp,-18*np.ones(300))
ms_disp = np.append(ms_disp,12*np.ones(100))

ynew = np.add(ynew,ms_disp)

xnewer = np.arange(0, 1000, 1.)
plt.figure()
plt.plot(x, y, 'x', xnewer, ynew)
plt.legend(['Linear', 'Cubic Spline'])
plt.title('Cubic-spline interpolation')
plt.savefig('img/grating_phase/displacements.pdf')
plt.show()

d_data = open('data/phase_displacement.data','w+')
np.save(d_data, ynew)
d_data.close()

'''
x = np.linspace(0, 100, num=101, endpoint=True)
smooth = interp1d(x, dis, kind='cubic')

xnew = np.linspace(0, 100, num=1001, endpoint=True)
plt.plot(x, dis, 'o', xnew, smooth(xnew), '-')
#x=np.arange(0,len(dis),1.)
#res=interp1d(x, dis, kind='cubic')
#print res
#plt.plot(dis)
plt.show()
'''