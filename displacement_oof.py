#this creates the displacement on the grating
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from microsaccades_functions import *

#displacement is zero and standard derivation is receptor size (0.5 arcmin)
mu, sigma = 0.,0.447

disp = np.random.normal(mu,sigma,1000)

d_data = open('data/oof/displacement.data','w+')
np.save(d_data, disp)
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