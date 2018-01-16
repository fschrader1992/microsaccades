#THIS IS THE PART OF THE NEURAL NETWORK. IT TAKES THE VALUES TO CALCULATE THE MEMBRANE POTENTIALS OF THE INPUT NEURONS AND THEN CALCULATES THE NETWORK OUTPUT 
#the unit in space is 1arcmin!
import pylab as pyl
import numpy as np
import cv2
import sys
import os
import matplotlib.pyplot as plt
#import nest
import datetime
#import nest.raster_plot
#import nest.topology as tp
from microsaccades_functions import *
'''
pgf_with_rc_fonts = {"font.family": "serif","font.serif": [],"font.sans-serif": []}
plt.rcParams.update(pgf_with_rc_fonts)

nest.ResetKernel()   # since we run the script multiple times
#to get different poisson outputs

msd = int(np.random.normal(5000,1000,1)[0])  #master seed
nest.SetKernelStatus({'local_num_threads' : 4})
n_vp = nest.GetKernelStatus('total_num_virtual_procs')
msdrange1 = range(msd, msd+n_vp)
pyrngs = [np.random.RandomState(s) for s in msdrange1]
msdrange2 = range(msd+n_vp+1, msd+1+2*n_vp)
nest.SetKernelStatus({'grng_seed': msd+n_vp,
                      'rng_seeds': msdrange2})

I_E=355.

def set_I_e_random(layer):
    r = nest.GetNodes(layer)[0]
    node_info=nest.GetStatus(r)
    localnodes=[(ni['global_id'],ni['vp']) for ni in node_info if ni['local']]
    for gid, vp in localnodes:
        nest.SetStatus([gid], {'I_e' : I_E+pyrngs[vp].uniform(-2.,2.)})
'''                       
#necessary paramater definitions
#frames = 100 #replaced by motion detectors only

#----------------------------------------------------------------------------INPUT-RATES-FROM-MS_INPUT

extent = 121.
delay = 7. #15. #30. # speed of point in poletti 2010 -> maybe increase a bit for more reacton later on

t_start = 0
t_end = 100 #315 #1000

weight = 40. #40.
weight_std = 1.5
#I_E = 410.

if len(sys.argv)==9:
    sim_title = sys.argv[1]
    sim_title_2 = sys.argv[2]
    sim_nr = sys.argv[3]
    handle_name = sys.argv[4]
    #extent = float(sys.argv[3])
    extent_x = float(sys.argv[5])
    extent_y = float(sys.argv[6])
    exp_nr = sys.argv[7]
    cond_nr = sys.argv[8]
else:
    sim_title = sys.argv[1]
    sim_nr = sys.argv[2]
    handle_name = sys.argv[3]
    #extent = float(sys.argv[3])
    extent_x = float(sys.argv[4])
    extent_y = float(sys.argv[5])
    exp_nr = sys.argv[6]
    cond_nr = sys.argv[7]
    sim_title_2 = ''

#delay = float(sys.argv[5])
vel = 0. #float(sys.argv[6])
#weight = float(sys.argv[7])
#extent_x = extent
#extent_y = extent
#center = extent/2.
center_x = extent_x/2.
center_y = extent_y/2.

print weight,delay
#weight = 0.

#here needs to be a part that transfers potentials into poisson rates
#m_file = open('data/midget_values.data','r+')
#m_file = open('data/'+str(sim_title)+'/'+str(sim_nr)+'/midget_rates_'+str(handle_name)+'_on.data','r+') 
#m_data = np.load(m_file)   
#m_file.close()

#p_file = open('data/parasolic_values.data','r+')

mmr = 300000.
pmr = 480000.

qr=[2,4,8,16,32,64,96]
parasolic_rates=[]
pmins=[]
pmaxs=[]
for q in qr:
    print q
    p_file = open('/home/schrader/Documents/microsaccades/data/jitter/chess_size'+str(q)+'_vel6on105off0/1/parasolic_rates_chess_size'+str(q)+'_vel6on105off0_on.data','r+') 
    p_data =np.load(p_file)  
    p_file.close()
    
    p_noise=[[np.random.normal(0,100,t_end) for j in range(len(p_data[0]))] for i in range(len(p_data))]
    if q==32 or q==64 or q==96:
        p_noise=[[np.random.normal(0,100,200) for j in range(len(p_data[0]))] for i in range(len(p_data))]
    if q==16:
        p_noise=[[np.random.normal(0,100,200) for j in range(len(p_data[0]))] for i in range(len(p_data))]
    pr=pmr*poissonRateParasols(p_data) +p_noise
    parasolic_rates+=[pr]
    
    #get minimum for colors
    pmins+=[pr.min()]
    pmaxs+=[pr.max()]
    
pmin=min(pmins)
pmax=max(pmaxs)

#midget_rates=poissonRateMidgets(m_data)

#print midget_rates

#to check for maximum spike rates in order to adopt conversion of film input
#maxs = []
#for i in range(len(midget_rates)):
#    for j in range(len(midget_rates[0])):
#        maxs.append(max(midget_rates[i][j][t_start:t_end]))
mmr = 300000.
#print maxs
#print 'midget max: ' + str(max(maxs)) + ' ' + str(mmr)
        
#paxs = []
#for i in range(len(parasolic_rates)):
#    for j in range(len(parasolic_rates[0])):
#       paxs.append(max(parasolic_rates[i][j][t_start:t_end]))
pmr = 480000.
#print 'parsolic max: ' + str(max(paxs)) + ' ' + str(pmr)

#ps_file = open('data/mo_det_cal/input_max_spikes.txt','a+')
#ps_file.write(str(handle_name)+' m: '+ str(max(maxs))+' p: '+ str(max(paxs))+'\n')
#ps_file.close()
        
#for rates 
mrs = []
prs=[]

#print midget_rates
#print parasolic_rates


#-----------------------------------------------------------------------------------------NETWORK-PART
#---------------------------------------------------------------------------INITIALIZE-POISSON-NEURONS

#storage = open('data/'+sim_title+'/network/'+handle_name+'.data','w+')
#storage.close()

#nest.SetKernelStatus({'resolution': 0.001,
#                      'overwrite_files': True,
#                      'data_path': 'data/'+sim_title+'/network/',
#                      'data_prefix': ''})

#set the initial Poisson rate
rate = 100.


#----------------------------------------------------------------------------------------CREATE-LAYERS

#rows = len(midget_rates)
#cols = len(midget_rates[0])
#print rows, cols
#rows = 40
#cols = 40

#get grid data form previous simulation
#gp_file = open('data/'+str(sim_title)+str(sim_title_2)+'/'+str(sim_nr)+'/p_pos_'+handle_name+'.data','r+')
#gp_file = open('data/test_vid/p_pos_test_vid.data','r+')
#gp_data = np.load(gp_file)  
#gp_file.close()
#gp_data = gp_data.tolist()


'''
#print gm_pos
print 'minimum midgets'
print gm_min, gm_max[0]+gm_min[0]+10.
print 'maximum midgets'
print gm_max, gm_max[1]+gm_min[1]+10.
'''       
#gp_pos=[]
#gp_r_0_pos=[]
#gp_r_60_pos=[]
#gm_r_120_pos=[]

#mean, stdv, timelength
#p_noise=[[np.random.normal(0,10,t_end) for j in range(len(gp_data[0]))] for i in range(len(gp_data))]
#parasolic_rates=pmr*parasolic_rates #+p_noise


#for f in range(t_start,t_end):#frames):
#print f
#reset rates
#fig = plt.figure()
#fig.set_size_inches(1,1)
#ax = plt.Axes(fig, [0., 0., 1., 1.])
#ax.set_axis_off()
#fig.add_axes(ax)

fig, (axes) = plt.subplots(nrows=3, ncols=2, )

for i in range(6):
    pr= parasolic_rates[i]
    im = axes[int(i/2)][i%2].imshow(pr[:,:,99],cmap=plt.cm.hot,vmin=0.,vmax=pmax,extent=[0,30,0,30])
    axes[int(i/2)][i%2].tick_params(labelsize=8)
    axes[int(i/2)][i%2].set_yticks(np.arange(0,31,10))
    if i%2==1:
        axes[int(i/2)][i%2].set_yticklabels([])
    if i/2<2:
        axes[int(i/2)][i%2].set_xticklabels([])

cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.5)
#plt.tight_layout()
plt.savefig('/home/schrader/Documents/microsaccades/img/jitter/chess_size_noise.pdf',  dpi = 128)
plt.show()
'''
my_image1 = np.linspace(0, 10, 10000).reshape(100,100)
my_image2 = np.sqrt(my_image1.T) + 3
my_image3 = np.sqrt(my_image1.T) + 3
my_image4 = np.sqrt(my_image1.T) + 3

fig, axes = plt.subplots(nrows=2, ncols=2)
im = axes[0][0].imshow(my_image1)
clim=im.properties()['clim']
axes[0][1].imshow(my_image2, clim=clim)
axes[1][0].imshow(my_image2, clim=clim)
axes[1][1].imshow(my_image2, clim=clim)

fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.5)
plt.tight_layout()

plt.show()

#ax.imshow(parasolic_rates[:,:,99], cmap=plt.cm.bwr, interpolation='nearest', animated=True) #gray
#plt.savefig('/home/schrader/Documents/microsaccades/img/animation/jitter/chess_input/'+str(f+1).zfill(3)+'.png',  dpi = 128)
#plt.savefig('/home/schrader/Documents/microsaccades/img/jitter/chess_size'+str(cond_nr)+'.png',  dpi = 128)
plt.close()
''' 