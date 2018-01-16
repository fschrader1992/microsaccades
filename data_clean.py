import pylab as pyl
import numpy as np
import cv2
import sys
import os
import datetime

sim_title = sys.argv[1]
if len(sys.argv)==5:
    sim_title_2 = sys.argv[2]
    sim_nr = sys.argv[3]
    fh = sys.argv[4]
    sar=True
else:
    sim_title_2 = ''
    sim_nr = sys.argv[2]
    fh = sys.argv[3]
    sar=False
    
layer_names = ['parasols','p_0_left','p_0_right','p_60_up','p_60_down','p_120_up','p_120_down']


for ln in layer_names:
    if sar==True:
        l_file = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(fh)+'.data'
    else:
        l_file = '/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_'+ln+'_'+str(fh)+'.data'
       
    with open(l_file) as f:
        file_str = f.read()

    # do stuff with file_str
    file_str=file_str.replace(']\n','linebreak')
    file_str=file_str.replace('\n','\t')
    file_str=file_str.replace('linebreak',']\n')

    with open(l_file, "w") as f:
        f.write(file_str)
    f.close()