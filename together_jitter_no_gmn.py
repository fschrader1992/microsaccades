import pylab as pyl
import numpy as np
import cv2
import sys
import os

sim_title = sys.argv[1]
if len(sys.argv)==5:
    sim_title_2 = sys.argv[2]
    sim_nr = sys.argv[3]
    fh = sys.argv[4]
else:
    sim_title_2 = ''
    sim_nr = sys.argv[2]
    fh = sys.argv[3]
    

r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/jitter_p_corr_max_spikes_no_gmn_all_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/murakami/jitter_p_corr_max_spikes_no_gmn_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'network/'+str(sim_nr)+'/p_corr_max_spikes_no_gmn.txt','r+')
for line in f:
    ltw = line
r.write(ltw)
r2.write(ltw)
f.close()
r.close()
r2.close()


r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/jitter_p_corr_max_spikes_no_gmn_200_all_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/murakami/jitter_p_corr_max_spikes_no_gmn_200_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'network/'+str(sim_nr)+'/p_corr_max_spikes_no_gmn_200.txt','r+')
for line in f:
    ltw = line
r.write(ltw)
r2.write(ltw)
f.close()
r.close()
r2.close()