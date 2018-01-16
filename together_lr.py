import pylab as pyl
import numpy as np
import cv2
import sys
import os
'''
sim_title = sys.argv[1]
sim_title_2 = ''
sim_nr = sys.argv[2]
v = sys.argv[3]
c = sys.argv[4]
sim_nr2 = sys.argv[5]
'''
if len(sys.argv)==10:
    sim_title = sys.argv[1]
    sim_title_2 = sys.argv[2]
    sim_nr = sys.argv[3]
    handle_name = sys.argv[4]
    #extent = float(sys.argv[3])
    extent_x = float(sys.argv[5])/2.
    extent_y = float(sys.argv[6])/2.-1.
    exp_nr = sys.argv[7]
    cond_nr = sys.argv[8]
    mdw = float(sys.argv[9])
    sar=True
else:
    sim_title = sys.argv[1]
    sim_nr = sys.argv[2]
    handle_name = sys.argv[3]
    #extent = float(sys.argv[3])
    extent_x = float(sys.argv[4])/2.
    extent_y = float(sys.argv[5])/2.-1.
    exp_nr = sys.argv[6]
    cond_nr = sys.argv[7]
    mdw = float(sys.argv[8])
    sim_title_2 = ''
    sar=False

   
r = open('/home/schrader/Documents/microsaccades/data/murakami/midgets_lr_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/spikes_m_0_left_'+handle_name+'.data','r+')
f2 = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/spikes_m_0_right_'+handle_name+'.data','r+')
f3 = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/spikes_midgets_'+handle_name+'.data','r+')
s = 0
for line in f:
    s += line.count('.')
  
s2 = 0
for line in f2:
    s2 += line.count('.')
	
s3 = 0
for line in f3:
    s3 += line.count('.')
    
sall=s+s2
 
#spc = float(s)/3980./2. #because 200ms simulated
f.close()
f2.close()
f3.close()
r.write(sim_nr+'\t'+cond_nr+'\t'+str(s3)+'\t'+str(sall)+'\t'+str(s)+'\t'+str(s2)+'\n')
r.close()
'''
r = open('/home/schrader/Documents/microsaccades/data/mtesting/m_right_spc_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/spikes_m_0_right_'+sim_title+'.data','r+')
s = 0
for line in f:
    s += line.count('.')
    
spc = float(s)/3980./2. #because 200ms simulated
f.close()
r.write(v+'\t'+c+'\t'+sim_nr+'\t'+str(s)+'\t'+str(spc)+'\n')
f.close()
r.close()
'''
'''
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
'''