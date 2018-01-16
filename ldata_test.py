import pylab as pyl
import numpy as np
import cv2
import sys
import os

sim_title = sys.argv[1]
if len(sys.argv)==7:
    sim_title_2 = sys.argv[2]
    sim_nr = sys.argv[3]
    fh = sys.argv[4]
    exp = sys.argv[5]
    cond = sys.argv[6]
    sar=True
else:
    sim_nr = sys.argv[2]
    fh = sys.argv[3]
    exp = sys.argv[4]
    cond = sys.argv[5]
    sar=False

layer_names = ['p_0_left','p_0_right','p_60_up','p_60_down','p_120_up','p_120_down']

for ln in layer_names:
    if sar==True:
        l_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(fh)+'_'+str(sim_title_2)+'_'+str(sim_nr)+'.data','r+')
    else:
        l_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(fh)+'.data','r+')
    narr = []
    for line in l_file:
        nid = int(line.split('[')[0])
        arr = line.split('[')[1]
        arr = arr.split(']')[0]
        s = np.array(arr)
        narr+=[(nid,s)]
        
    print narr
    l_file.close()
