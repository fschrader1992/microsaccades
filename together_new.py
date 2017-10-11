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
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'m_all_spikes_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/poletti2010/m_all_spikes.txt','a+')
for direct in directions:
    if sar==True:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_m_'+direct+'_'+str(fh)+'_'+str(sim_title_2)+'_'+str(sim_nr)+'.data','r+')
    else:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_m_'+direct+'_'+str(fh)+'.data','r+')
    s = 0
    for line in f:
        line = line.split('[')[1]
        line = line.split(']')[0]
        s += len(line.split(','))
    sl+=[s]
    f.close()
sl+=[sum(sl)]
r.write(fh+'\t'+exp+'\t'+cond+'\t')
r2.write(fh+'\t'+exp+'\t'+cond+'\t')
for ds in sl:
    r.write(str(ds)+'\t')
    r2.write(str(ds)+'\t')
r.write('\n')
r2.write('\n')
r.close()
r2.close()

sl=[]
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'p_all_spikes_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/poletti2010/p_all_spikes.txt','a+')
for direct in directions: 
    if sar==True:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_p_'+direct+'_'+str(fh)+'_'+str(sim_title_2)+'_'+str(sim_nr)+'.data','r+')
    else:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_p_'+direct+'_'+str(fh)+'.data','r+')
    s = 0
    for line in f:
        line = line.split('[')[1]
        line = line.split(']')[0]
        s += len(line.split(','))
    sl+=[s]
    f.close()
sl+=[sum(sl)]
r.write(fh+'\t'+exp+'\t'+cond+'\t')
r2.write(fh+'\t'+exp+'\t'+cond+'\t')
for ds in sl:
    r.write(str(ds)+'\t')
    r2.write(str(ds)+'\t')
r.write('\n')
r2.write('\n')
r.close()
r2.close()



sl=[]
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'m_all_spikes_200_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/poletti2010/m_all_spikes_200.txt','a+')
for direct in directions: 
    if sar==True:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_m_'+direct+'_200_'+str(fh)+'_'+str(sim_title_2)+'_'+str(sim_nr)+'.data','r+')
    else:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_m_'+direct+'_200_'+str(fh)+'.data','r+')
    s = 0
    for line in f:
        line = line.split('[')[1]
        line = line.split(']')[0]
        s += len(line.split(','))
    sl+=[s]
    f.close()
sl+=[sum(sl)]
r.write(fh+'\t'+exp+'\t'+cond+'\t')
r2.write(fh+'\t'+exp+'\t'+cond+'\t')
for ds in sl:
    r.write(str(ds)+'\t')
    r2.write(str(ds)+'\t')
r.write('\n')
r2.write('\n')
r.close()
r2.close()

sl=[]
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'p_all_spikes_200_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/poletti2010/p_all_spikes_200.txt','a+')
for direct in directions: 
    if sar==True:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_p_'+direct+'_200_'+str(fh)+'_'+str(sim_title_2)+'_'+str(sim_nr)+'.data','r+')
    else:
        f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'network/'+str(sim_nr)+'/spikes_p_'+direct+'_200_'+str(fh)+'.data','r+')
    s = 0
    for line in f:
        line = line.split('[')[1]
        line = line.split(']')[0]
        s += len(line.split(','))
    sl+=[s]
    f.close()
sl+=[sum(sl)]
r.write(fh+'\t'+exp+'\t'+cond+'\t')
r2.write(fh+'\t'+exp+'\t'+cond+'\t')
for ds in sl:
    r.write(str(ds)+'\t')
    r2.write(str(ds)+'\t')
r.write('\n')
r2.write('\n')
#print sl
r.close()
r2.close()