import pylab as pyl
import numpy as np
import cv2
import sys
import os

sim_title = sys.argv[1]
sim_title_2 = ''
sim_nr = sys.argv[2]
v = sys.argv[3]
c = sys.argv[4]
sim_nr2 = sys.argv[5]
   
r = open('/home/schrader/Documents/microsaccades/data/mtesting/m_left_spc_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/mtesting/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_m_0_left_'+sim_title+'.data','r+')
s = 0
for line in f:
    s += line.count('.')
    
spc = float(s)/3980./2. #because 200ms simulated
f.close()
r.write(v+'\t'+c+'\t'+sim_nr2+'\t'+str(s)+'\t'+str(spc)+'\n')
f.close()
r.close()

r = open('/home/schrader/Documents/microsaccades/data/mtesting/m_right_spc_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/mtesting/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_m_0_right_'+sim_title+'.data','r+')
s = 0
for line in f:
    s += line.count('.')
    
spc = float(s)/3980./2. #because 200ms simulated
f.close()
r.write(v+'\t'+c+'\t'+sim_nr+'\t'+str(s)+'\t'+str(spc)+'\n')
f.close()
r.close()


r = open('/home/schrader/Documents/microsaccades/data/mtesting/p_left_spc_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/mtesting/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_p_0_left_'+sim_title+'.data','r+')
s = 0
for line in f:
    s += line.count('.')
    
spc = float(s)/240./2. #because 200ms simulated
f.close()
r.write(v+'\t'+c+'\t'+sim_nr+'\t'+str(s)+'\t'+str(spc)+'\n')
f.close()
r.close()

r = open('/home/schrader/Documents/microsaccades/data/mtesting/p_right_spc_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/mtesting/'+str(sim_title)+str(sim_title_2)+'/network/'+str(sim_nr)+'/spikes_p_0_right_'+sim_title+'.data','r+')
s = 0
for line in f:
    s += line.count('.')
    
spc = float(s)/240./2. #because 200ms simulated
f.close()
r.write(v+'\t'+c+'\t'+sim_nr+'\t'+str(s)+'\t'+str(spc)+'\n')
f.close()
r.close()
