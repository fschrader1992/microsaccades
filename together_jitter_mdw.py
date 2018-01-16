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
    
if len(sys.argv)==11:
    sim_title = sys.argv[1]
    sim_title_2 = sys.argv[2]
    sim_nr = sys.argv[3]
    handle_name = sys.argv[4]
    #extent = float(sys.argv[3])
    extent_x = float(sys.argv[5])/2.
    extent_y = float(sys.argv[6])/2.-1.
    exp_nr = sys.argv[7]
    cond_nr = sys.argv[8]
    size = float(sys.argv[9])
    mdw = float(sys.argv[10])
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
    size = float(sys.argv[8])
    mdw = float(sys.argv[9])
    sim_title_2 = ''
    sar=False
    
    
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/jitter_p_corr_max_spikes_mdw_4_modet_gmn_40_all_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/jitter/jitter_p_corr_spike_count_mdw_4_modet_gmn_40_all.txt','a+')

par_times_pos_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/par_times_pos.txt','r+')
mr_times_pos_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/mr_times_pos_corr_modet_gmn_40.txt','r+') #_'+str(mdw)+'
ml_times_pos_file = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/network/'+str(sim_nr)+'/ml_times_pos_corr_modet_gmn_40.txt','r+') #_'+str(mdw)+'

par_times = 0.
ml_times = 0.
mr_times = 0.

par_times_pos_file.seek(0)
for line in par_times_pos_file:
    line = line.replace('(','').replace(')','')
    if int(line.split(',')[0]) > 105:
        par_times += 1.
        
ml_times_pos_file.seek(0)
for line in ml_times_pos_file:
    line = line.replace('(','').replace(')','')
    if int(line.split(',')[0]) > 105:
        ml_times += 1.
        
mr_times_pos_file.seek(0)
for line in mr_times_pos_file:
    line = line.replace('(','').replace(')','')
    if int(line.split(',')[0]) > 105:
        mr_times += 1.
 
ltw=str(exp_nr)+'\t'+str(cond_nr)+'\t'+str(size)+'\t'+str(mdw)+'\t'+str(par_times)+'\t'+str(ml_times)+'\t'+str(mr_times)+'\n'
r.write(ltw)
r2.write(ltw)
par_times_pos_file.close()
ml_times_pos_file.close()
mr_times_pos_file.close()
r.close()
r2.close()

'''
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/jitter_p_corr_max_spikes_modet_gmn_40_200_all_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/jitter/jitter_p_corr_max_spikes_modet_gmn_40_200_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'network/'+str(sim_nr)+'/p_corr_max_spikes_modet_gmn_40_200.txt','r+')
for line in f:
    ltw = line
r.write(ltw)
r2.write(ltw)
f.close()
r.close()
r2.close()
'''


'''
r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/jitter_p_corr_max_spikes_no_gmn_all_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/jitter/jitter_p_corr_max_spikes_no_gmn_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'network/'+str(sim_nr)+'/p_corr_max_spikes_no_gmn.txt','r+')
for line in f:
    ltw = line
r.write(ltw)
r2.write(ltw)
f.close()
r.close()
r2.close()


r = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+'/jitter_p_corr_max_spikes_no_gmn_200_all_'+str(fh)+'.txt','w+')
r2 = open('/home/schrader/Documents/microsaccades/data/jitter/jitter_p_corr_max_spikes_no_gmn_200_all.txt','a+')
f = open('/home/schrader/Documents/microsaccades/data/'+str(sim_title)+str(sim_title_2)+'network/'+str(sim_nr)+'/p_corr_max_spikes_no_gmn_200.txt','r+')
for line in f:
    ltw = line
r.write(ltw)
r2.write(ltw)
f.close()
r.close()
r2.close()
'''