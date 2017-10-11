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
    
    
directions = ['0_left','0_right','60_up','60_down','120_up','120_down']

sl=[]
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'m_cells_spiking_3_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/poletti2010/m_cells_spiking_3.txt','a+')
cells = 0
for direct in directions:
    if sar==True:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_m_'+direct+'_'+str(fh)+'_'+str(sim_title_2)+'_'+str(sim_nr)+'.data','r+')
    else:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_m_'+direct+'_'+str(fh)+'.data','r+')
    for line in f:
        gid = int(line.split('[')[0])
        if gid%3==0:
            cells+=1
    f.close()
r.write(fh+'\t'+exp+'\t'+cond+'\t'+str(cells)+'\n')
r2.write(fh+'\t'+exp+'\t'+cond+'\t'+str(cells)+'\n')
r.close()
r2.close()

cells = 0
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'p_cells_spiking_3_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/poletti2010/p_cells_spiking_3.txt','a+')
for direct in directions: 
    if sar==True:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_p_'+direct+'_'+str(fh)+'_'+str(sim_title_2)+'_'+str(sim_nr)+'.data','r+')
    else:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_p_'+direct+'_'+str(fh)+'.data','r+')
    for line in f:
        gid = int(line.split('[')[0])
        if gid%3==0:
            cells+=1
    f.close()
r.write(fh+'\t'+exp+'\t'+cond+'\t'+str(cells)+'\n')
r2.write(fh+'\t'+exp+'\t'+cond+'\t'+str(cells)+'\n')
r.close()
r2.close()


cells = 0
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'m_cells_spiking_3_200_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/poletti2010/m_cells_spiking_3_200.txt','a+')
for direct in directions: 
    if sar==True:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_m_'+direct+'_200_'+str(fh)+'_'+str(sim_title_2)+'_'+str(sim_nr)+'.data','r+')
    else:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_m_'+direct+'_200_'+str(fh)+'.data','r+')
    for line in f:
        gid = int(line.split('[')[0])
        if gid%3==0:
            cells+=1
    f.close()
r.write(fh+'\t'+exp+'\t'+cond+'\t'+str(cells)+'\n')
r2.write(fh+'\t'+exp+'\t'+cond+'\t'+str(cells)+'\n')
r.close()
r2.close()

cells = 0
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'p_cells_spiking_3_200_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/poletti2010/p_cells_spiking_3_200.txt','a+')
for direct in directions: 
    if sar==True:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_p_'+direct+'_200_'+str(fh)+'_'+str(sim_title_2)+'_'+str(sim_nr)+'.data','r+')
    else:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_p_'+direct+'_200_'+str(fh)+'.data','r+')
    for line in f:
        gid = int(line.split('[')[0])
        if gid%3==0:
            cells+=1
    f.close()
r.write(fh+'\t'+exp+'\t'+cond+'\t'+str(cells)+'\n')
r2.write(fh+'\t'+exp+'\t'+cond+'\t'+str(cells)+'\n')
r.close()
r2.close()