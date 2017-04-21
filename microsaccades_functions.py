import pylab as pyl
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

def spatialFilter(x,x0,sigma,alpha,beta):
    return (np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))-alpha*np.exp(-(x-x0)*(x-x0)/(2*beta*beta*sigma*sigma)))/sigma
#above is Garrett's version
#def spatialFilter(x,x0,sigma,alpha,beta):
#    return np.exp((x-x0)*(x-x0)/(2*sigma*sigma))-alpha*np.exp((x-x0)*(x-x0)/(2*beta*beta*sigma*sigma))

def poissonRateMidgets(pot):
    return abs(pot/518311.*60.) #probably change a bit, but this is result from opposite (maximum) stimulation

def poissonRateParasols(pot):
    return abs(pot/5526366.*60.)
