#!/bin/python
#THIS PART READS THE PIXELS OF A MOVIE AND APPLIES THE TEMPORAL FILTERS TO THE PROPOSED MODEL FOR MICROSACCADES. THE OUTPUT ARE THE POTENTIAL VALUES FOR THE POISSON RATES.
'''
loading the video into an array
-> first get all the pixels in a frame (grayscale)
-> after that calculate the values of each tempral flter at each time and store them in another array which will be the basis for changing poisson rates
'''

import sys
import os
import glob
import pylab as pyl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cv2
import datetime
import itertools 
#import multiprocessing
from microsaccades_functions import *

#---------------------------------------------------------------------------------------------------LOAD-FRAMES
now = datetime.datetime.now()
'''
#for parallel computing
try:
    cpus = multiprocessing.cpu_count()
except NotImplementedError:
    cpus = 2   # arbitrary default
	
pool = multiprocessing.Pool(processes=cpus)
'''
sim_title = sys.argv[1]
handle_name = sys.argv[2]


#alternative method for videos
frames = []
frame_number = 400
os.chdir("video/img_input/" + handle_name)
for file in glob.glob("second*.png"):
    #print(file)
    frames+=[file]
frames.sort()
#print frames

f=cv2.imread(frames[0])
height, width = f.shape[:2]
dt = 5.
print height, width

#assign the all time all pixel array/list
pixels4d = [[[] for j in range(width)] for i in range(height)]


for file in frames:

    frame = cv2.imread(file)
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    #store 2D array in pixels4d 
    for i in range(height):
        for j in range(width):
            pixels4d[i][j]+=[float(gray[i,j])]
    
cv2.destroyAllWindows()

os.chdir("../../..")
'''

#load video
cap = cv2.VideoCapture('video/' + str(sim_title) +'/'+ str(handle_name) + '.mp4')
#just to be sure
while not cap.isOpened():
    cap = cv2.VideoCapture('video/' + str(sim_title) +'/'+ str(handle_name) + '.mp4')
    cv2.waitKey(1000)
    print "Wait for the header"
    
#data of video
width = int(cap.get(3))
height = int(cap.get(4))
#dt=1/video_fps
dt = 1. #1/cap.get(5)
frmct = cap.get(7)
print(width, height, dt, frmct)

#assign the all time all pixel array/list
pixels4d = [[[] for j in range(width)] for i in range(height)]


while(cap.isOpened()):
    ret, frame = cap.read()
    #check, whether there's a frame left, if not break loop
    if ret == True:
        frame_number = int(cap.get(1))
        #will probably be abundant later on + what about color vision?
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        #store 2D array in pixels4d 
        for i in range(height):
            for j in range(width):
                pixels4d[i][j]+=[float(gray[i,j])]
        if frame_number == 2:
            break
        cv2.imshow('frame',gray)
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break
    else:
        break
    
cap.release()
cv2.destroyAllWindows()
'''
#ffile =open('pixels4d.data','w+')
#ffile.write(str(pixels4d[9][9]))
#ffile.close

#----------------------------------------------------------------------------------------------------INITIALIZE
#now comes the part of the temporal filters
#later: define list with all different tau1, tau2 and p, for the moment one is enough (for on/off-center) 


#problem here is now, that both create the same output and therefore nothing is produced
#off_tau1 = 5 #.05 --> what is dt in Garrett's code? 0.1ms or 1ms? depending on that, this is the 
#off_tau2 = 15 #.1     right value
#off_p = .8 #.05

on_tau1 = 5. #.2
on_tau2 = 15. #.1
on_p = .8 #.05

#set values for the spatial filters for the pixel in \mum(?)
alpha = .1
beta = 5. #CSR
sigma = 1. # parasolic/midget cell ratio * sigma (=1)
px_midget_ratio = 2. #pixel to receptor ratio, needed for receptor distance/number
px_dist= 0.5 #in arcmin, from paper
par_m_ratio = 4. #2arcmin center field sigma from paper

#spatial filter values determined
spat_filter_break_radius = 6 #filter radius in px 

#subsequent values
midget_dist = px_midget_ratio*px_dist #1.0

#we actually need to calculate the values of each temporal filter just once
temp_filter_on = tempFilter(40,dt,on_tau1,on_tau2,on_p) #mayber 120, compare to paper --> 200 in Garrett's code!      


midget_height = int(np.floor((height-spat_filter_break_radius)/(1.732*px_midget_ratio))) #2*cos(30deg) #2 cells per unit
midget_width = int(np.floor((width-spat_filter_break_radius)/(3.*px_midget_ratio))*4.) #3 is difference between hexagons where pattern repeatesand 4 is the number of cells for that unit

parasol_height = int(np.floor((height-spat_filter_break_radius)/(px_midget_ratio*1.732*1.5))) #2*cos(30deg) #3 cells per unit
parasol_width = int(np.floor((width-spat_filter_break_radius)/(1.5*px_midget_ratio))) #3 is difference between hexagons where pattern repeatesand 4 is the number of cells for that unit

print midget_height,midget_width,parasol_height,parasol_width

midget_pixels4d = np.zeros(shape=(midget_height,midget_width,frame_number))

grid=[[0],[0]]
m_grid=[[0],[0]]
pgrid=[[0],[0]]

midget_grid= [[(0,0) for j in range(midget_width)] for i in range(midget_height)]

#---------------------------------------------------------------------------------------MIDGETS-SPATIAL-FILTERS

def spatFilterPx(mult_list,center,kl,r_break,sigma,alpha,beta):
    add_val = 0
    
    ml_h = len(mult_list)
    ml_w = len(mult_list[0])
    klv_w = 0
    klv_h = 0
    new_kl0,new_kl1 = kl
    
    if kl[0] >= ml_h:
        new_kl0 = kl[0] - ml_h
        klv_h = height/px_midget_ratio
    if kl[0] < 0:
        new_kl0 = kl[0] + ml_h
        klv_h = -height/px_midget_ratio
    if kl[1] >= ml_w:
        new_kl1 = kl[1] - ml_w
        klv_w = width/px_midget_ratio
    if kl[1] < 0:
        new_kl1 = kl[1] + ml_w
        klv_w = -width/px_midget_ratio
        
    kl = (new_kl0,new_kl1)
    grid_shift = (midget_dist*klv_h,midget_dist*klv_w)
    
    mult_val = mult_list[kl[0]][kl[1]]
    kl = float(kl[0])*px_dist,float(kl[1])*px_dist
    dist = np.sqrt(sum(itertools.imap(lambda x,y,z: (x-y-z)*(x-y-z), center, kl, grid_shift)))
    #make it a real circle
    if dist <= r_break:
        add_val = mult_val*spatialFilter(dist,0,sigma,alpha,beta)
    return add_val

def getSpatFilter(ij):
    
    #pos_i: height, pos_j: width
    pos_i = 1.732*ij[0]        #2*cos(30deg)
    pos_j = 3.*np.floor(float(ij[1])/4)
    
    if ij[1]%4 == 0:
        pos_i += 0.866  #cos(30deg)
    elif ij[1]%4 == 1:
        pos_j += .5
    elif ij[1]%4 == 2:
        pos_j += 1.5
    elif ij[1]%4 == 3:
        pos_j += 2.
        pos_i += 0.866  #cos(30deg)
        
    i_low = int(pos_i*px_midget_ratio-spat_filter_break_radius) # -> all j with dist < r
    i_ceil = int(pos_i*px_midget_ratio+spat_filter_break_radius)+1
    j_low = int(pos_j*px_midget_ratio-spat_filter_break_radius) # -> all j with dist < r
    j_ceil = int(pos_j*px_midget_ratio+spat_filter_break_radius)+1
    
    pos_i=midget_dist*pos_i
    pos_j=midget_dist*pos_j
    
    #right, since in opencv (height, width)
    pos=(pos_i,pos_j)
    
    if f==0:
        grid[0] += [ij[0]]
        grid[1] += [ij[1]]
        m_grid[0] += [pos_i]
        m_grid[1] += [pos_j]
        midget_grid[ij[0]][ij[1]] = (pos_i, pos_j) 
    
    #spatial filters can allow negative potentials
    midget_pixels4d[ij[0]][ij[1]][f] = sum(itertools.imap(lambda x: spatFilterPx(pixels3d,pos,x,spat_filter_break_radius,sigma,alpha,beta), itertools.product(range(i_low,i_ceil),range(j_low,j_ceil))))


#apply the spatial filter 
for f in range(frame_number):
    print f
    pixels3d = [[item[f] for item in pxst] for pxst in pixels4d]      
    map(lambda x: getSpatFilter(x), itertools.product(range(midget_height),range(midget_width)))
     
#pyl.figure()
#pyl.subplot(121, aspect='equal')
#pyl.plot(pgrid[0],pgrid[1],'ro')

#plt.plot(m_grid[1],m_grid[0],'bo')
#plt.plot(m_grid[0][:20],m_grid[1][:20],'bo')
#plt.plot(m_grid[0][158:180],m_grid[1][158:180],'bo')
#plt.plot(m_grid[0][316:322],m_grid[1][316:322],'bo')
#plt.show()

#output frame for the filter values at the end
temp_filter_vals  = np.zeros(shape=(midget_height,midget_width,frame_number))
temp_filter_vals_on  = np.zeros(shape=(midget_height,midget_width,frame_number))
temp_filter_vals_off = np.zeros(shape=(midget_height,midget_width,frame_number))

#--------------------------------------------------------------------------------------MIDGETS-TEMPORAL-FILTERS
print '---'
#now apply the temporal filters: output are two lists for on/off fields
for i in range(midget_height):
    print i
    for j in range(midget_width):
        temp_midget_px4d = []
        pop = temp_midget_px4d.pop
        for f in range(frame_number):
            temp_midget_px4d.insert(0, float(midget_pixels4d[i][j][f]))
            if f > 40:
                pop()
            #add a new entry to the time list of the pixel i,j
            trs = sum(itertools.imap(lambda x,y: x*y, temp_midget_px4d, temp_filter_on))
            temp_filter_vals_on[i][j][f] = trs if trs > 0 else 0
            temp_filter_vals[i][j][f] = trs 
            temp_filter_vals_off[i][j][f] = -trs if trs < 0 else 0
            #temp_filter_vals_off[i][j][f] = sum(itertools.imap(lambda x,y: x*y, temp_midget_px4d, temp_filter_off))


#calculate the difference between surround and center fields --> needs to be done? and safe to file
m_output = np.asarray(temp_filter_vals)#np.subtract(temp_filter_vals_on, temp_filter_vals_off), dtype=int)
m_output_on = np.asarray(temp_filter_vals_on)
m_output_off = np.asarray(temp_filter_vals_off)

#--------------------------------------------------------------------------------------PARASOLS-SPATIAL-FILTERS

par_values = [[[0. for f in range(frame_number)] for j in range(parasol_width)] for i in range(parasol_height)] 
par_values_on = [[[0. for f in range(frame_number)] for j in range(parasol_width)] for i in range(parasol_height)] 
par_values_on_off = [[[0. for f in range(frame_number)] for j in range(parasol_width)] for i in range(parasol_height)] 

def spatFilterParasolPx(mult_list,center,kl,r_break,sigma,alpha,beta):
    add_val = 0
    
    ml_h = len(mult_list)
    ml_w = len(mult_list[0])
    klv_w = 0
    klv_h = 0
    new_kl0,new_kl1 = kl
    
    if kl[0] >= ml_h:
        new_kl0 = kl[0] - ml_h
        klv_h = height/px_midget_ratio
    if kl[0] < 0:
        new_kl0 = kl[0] + ml_h
        klv_h = -height/px_midget_ratio
    if kl[1] >= ml_w:
        new_kl1 = kl[1] - ml_w
        klv_w = width/px_midget_ratio
    if kl[1] < 0:
        new_kl1 = kl[1] + ml_w
        klv_w = -width/px_midget_ratio
    kl = (new_kl0,new_kl1)

    grid_shift = (midget_dist*klv_h,midget_dist*klv_w)
    
    #print kl
    mult_val = mult_list[kl[0]][kl[1]]
    kl_val = midget_grid[kl[0]][kl[1]] 
    dist = np.sqrt(sum(itertools.imap(lambda x,y,z: (x-y-z)*(x-y-z), center, kl_val, grid_shift)))
    #make it a real circle
    if dist <= r_break:
        add_val = mult_val*spatialFilter(dist,0,sigma,alpha,beta)
    return add_val


def getSpatFilterParasol(ij):
    
    #ij[0] height, ij[1] width    
    
    
    #pos = [midget_grid[y_disp][x_disp][0], midget_grid[y_disp][x_disp][1] + 1.]
    #print 'ij' + str(ij)
    
    #step size in each direction vice versa compared to midgets
    #x_disp = 3.*ij[1]+1.
    #y_disp = 1.5*1.732*ij[0]#1.5 for average displacement
    
    x_disp = 1.5*ij[1]+1.
    y_disp = 1.5*1.72*ij[0]
    
    if ij[0]%2==0 and ij[1]%2==0: 
        y_disp+=+0.866
    elif ij[0]%2!=0 and ij[1]%2!=0:
        y_disp+=+0.866
    #if ij[0]%4 == 0:
    #pos[0] += 0.866  #cos(30deg)
    #if ij[0]%4 == 1 or ij[0]%4 == 2:
    #    pos[1] += 3.
    #    x_disp += 3.
    
    pos = [y_disp,x_disp]

    x_disp = 2*ij[1]+1
    y_disp = int(ij[0]/2)*3+ij[0]%2 #together +2*ij[0]
    move = par_m_ratio*spat_filter_break_radius/px_midget_ratio
    i_low = int(y_disp-move) # -> all j with dist < r
    i_ceil = int(y_disp+move)+1  
    j_low = int(x_disp-2.*move) # -> all j with dist < r
    j_ceil = int(x_disp+2.*move)+1

    #print j_low,j_ceil
    if f==0:
        pgrid[0] += [pos[0]]
        pgrid[1] += [pos[1]]
    
    pos[0]= midget_dist*pos[0]
    pos[1]= midget_dist*pos[1]
    par_values[ij[0]][ij[1]][f] = sum(itertools.imap(lambda x: spatFilterParasolPx(midgets3d,pos,x,spat_filter_break_radius,par_m_ratio*sigma,alpha,beta), itertools.product(range(i_low,i_ceil),range(j_low,j_ceil))))
    #par_values_on[ij[0]][ij[1]][f] = sum(itertools.imap(lambda x: spatFilterParasolPx(midgets3d_on,pos,x,spat_filter_break_radius,par_m_ratio*sigma,alpha,beta), itertools.product(range(i_low,i_ceil),range(j_low,j_ceil))))
    #par_values_on_off[ij[0]][ij[1]][f] = sum(itertools.imap(lambda x: spatFilterParasolPx(midgets3d_on_off,pos,x,spat_filter_break_radius,par_m_ratio*sigma,alpha,beta), itertools.product(range(i_low,i_ceil),range(j_low,j_ceil))))


print int(midget_width/4.)
#apply the spatial filter 
for f in range(frame_number):
    print f
    midgets3d = [[mi[f] for mi in midgs] for midgs in temp_filter_vals]   
    #midgets3d_on = [[mi[f] for mi in midgs] for midgs in temp_filter_vals_on]
    #midgets3d_on_off = [[mi[f] for mi in midgs] for midgs in temp_filter_vals_off]  
    #for i in range(len(temp_filter_vals)):
    #    for j in range(len(temp_filter_vals[0])):
    #        midgets3d_on_off[i][j] += temp_filter_vals_on[i][j][f]
    map(lambda x: getSpatFilterParasol(x), itertools.product(range(parasol_height), range(parasol_width)))

#plt.plot(pgrid[0][:10],pgrid[1][:10],'go')
#plt.plot(pgrid[0][120:130],pgrid[1][120:130],'go')
#plt.plot(pgrid[1],pgrid[0],'go')
#plt.plot(pgrid[0][55:60],pgrid[1][55:60],'go')
#plt.plot(pgrid[0][105:110],pgrid[1][105:110],'go')
#plt.plot(pgrid[0][75:82],pgrid[1][75:82],'go')
#plt.savefig('img/grid_large.pdf')
#plt.show()

#--------------------------------------------------------------------------------------------------------OUTPUT
p_output = np.asarray(par_values)
p_output_on = np.asarray(par_values_on)
p_output_on_off = np.asarray(par_values_on_off)

ms = len(midget_grid)*len(midget_grid[0])/4.
ps = len(range(int(2*midget_height/3)))*len(range(int(midget_width/4)))
print ms, ps
#pyl.figure()
#pyl.subplot(121, aspect='equal')
#pyl.plot(pgrid[0],pgrid[1],'ro')
#pyl.show()

#plt.plot(pgrid[0],pgrid[1],'ro')
#out= 'img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/grid_mp_' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '.pdf'
#plt.savefig(out)
#plt.savefig('img/grid_mp.pdf')
#plt.show()

#---------------------------------------------------------------------------------------------------SAVE-OUTPUT
m_data = open('data/'+sim_title+'/midget_rates_'+str(handle_name)+'.data','w+')
np.save(m_data, m_output)
m_data.close()

p_data = open('data/'+sim_title+'/parasolic_rates_'+str(handle_name)+'.data','w+')
np.save(p_data, p_output)
p_data.close()

m_data = open('data/'+sim_title+'/midget_rates_'+str(handle_name)+'_on.data','w+')
np.save(m_data, m_output)
m_data.close()

p_data = open('data/'+sim_title+'/parasolic_rates_'+str(handle_name)+'_on.data','w+')
np.save(p_data, p_output)
p_data.close()

m_data = open('data/'+sim_title+'/midget_rates_'+str(handle_name)+'_off.data','w+')
np.save(m_data, m_output)
m_data.close()

p_data = open('data/'+sim_title+'/parasolic_rates_'+str(handle_name)+'_on_off.data','w+')
np.save(p_data, p_output)
p_data.close()

#-----------------------------------------------------------------------------------------------PLOT-SOME-STUFF


fig = plt.figure(1)

ax = fig.add_subplot(221)
ax.set_title('midget output')
plt.imshow(m_output[:,:,250], aspect='auto', interpolation='nearest')
ax.set_aspect('equal')
plt.axis('off')

cax = fig.add_axes([0.,0.,1.,1.])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
#plt.colorbar(orientation='vertical')

ax = fig.add_subplot(223)
ax.set_title('parasolic output')
plt.imshow(p_output[:,:,250], aspect='auto', interpolation='nearest')
ax.set_aspect('equal')
plt.axis('off')

cax = fig.add_axes([0.,0.,1.,.5])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)

#out= 'img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/output_mp_215fr_opposite_1fr0deg_' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '.pdf'
#plt.savefig(out)

out= 'img/video/' + str(sim_title) +'/'+ str(handle_name) + '.pdf'
plt.savefig(out)

#plt.show()
