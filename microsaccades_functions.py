import pylab as pl
import numpy as np

#this returns the result of a normal temporal filter (double exponential)
def tempFilterFct(t,tau1,tau2,p):
    return t*t*t/(tau1*tau1*tau1*tau1)*np.exp(-t/tau1)-p*t*t*t/(tau2*tau2*tau2*tau2)*np.exp(-t/tau2)

def tempFilter(f_max,dt,tau1,tau2,p):
    vals = []
    for f in range(f_max):
        t = f*dt
        vals += [tempFilterFct(t,tau1,tau2,p)]
    return vals

'''def setTempFilterVals(f_max):
    for f in range(f_max):
        #add a new entry to the time list of the pixel i,j
        temp_on = 0
        temp_off = 0
        #assign another counting variable 
        map()
        for k in range(f):
            #later on, if we need to summarize over the already partly summed fields
            #-> pixels4d must be replaced by on/off filed output + be careful about minus!
            temp_on += float(temp_px4d[f-k])*temp_filter_on[k] 
            temp_off += float(temp_px4d[f-k])*temp_filter_off[k]
        temp_filter_vals_on[i][j][f] = temp_on
        temp_filter_vals_off[i][j][f] = temp_off'''