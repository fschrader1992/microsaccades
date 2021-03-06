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

pgf_with_rc_fonts = {"font.family": "serif","font.serif": [],"font.sans-serif": []}
plt.rcParams.update(pgf_with_rc_fonts)
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
fn = sys.argv[3]


#alternative method for videos
frames = []
frame_number = int(fn)
os.chdir("video/img_input/" + sim_title)
for file in glob.glob("second*.png"):
    #print(file)
    frames+=[file]
frames.sort()
#print frames

f=cv2.imread(frames[0])
height, width = f.shape[:2]
dt = 1.
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
beta = 3. #CSR
sigma = 1. # parasolic/midget cell ratio * sigma (=1)
px_midget_ratio = 1. #pixel to receptor ratio, needed for receptor distance/number
px_dist= 0.5 #in arcmin, from paper
par_m_ratio = 4. #2arcmin center field sigma from paper

#spatial filter values determined
spat_filter_break_radius = 6 #filter radius in px 

#subsequent values
midget_dist = px_midget_ratio#1, rest gone, since unit is 1px #*px_dist #0.5

#we actually need to calculate the values of each temporal filter just once
temp_filter_on = tempFilter(200,dt,on_tau1,on_tau2,on_p) #mayber 120, compare to paper --> 200 in Garrett's code!      


#midget_height = int(np.floor((height-spat_filter_break_radius)/(1.732*px_midget_ratio))) #2*cos(30deg) #2 cells per unit
#midget_width = int(np.floor((width-spat_filter_break_radius)/(3.*px_midget_ratio))*4.) #3 is difference between hexagons where pattern repeatesand 4 is the number of cells for that unit

midget_height = int(height/(px_midget_ratio*0.866))
midget_width = int(width/px_midget_ratio)

#parasol_height = int(np.floor((height-spat_filter_break_radius)/(px_midget_ratio*1.732*1.5))) #2*cos(30deg) #3 cells per unit
#parasol_width = int(np.floor((width-spat_filter_break_radius)/(1.5*px_midget_ratio))) #3 is difference between hexagons where pattern repeatesand 4 is the number of cells for that unit

parasol_height = int(height/(par_m_ratio*0.866))
parasol_width = int(width/par_m_ratio)

print midget_height,midget_width,parasol_height,parasol_width

midget_pixels4d = np.zeros(shape=(midget_height,midget_width,frame_number))
midget_pixels4d = [[[] for j in range(midget_width)] for i in range(midget_height)]

grid=[[0],[0]]
m_grid=[[0],[0]]
pgrid=[[0],[0]]

midget_grid= [[(0,0) for j in range(midget_width)] for i in range(midget_height)]

#---------------------------------------------------------------------------------------MIDGETS-SPATIAL-FILTERS 
print '---'
print 'start of midget part'

def spatFilterPx(center,kl,r_break,sigma,alpha,beta):

    klv_w = 0
    klv_h = 0
    new_kl0,new_kl1 = kl
    kl_pos = float(kl[0]),float(kl[1])
    
    dis_x = (center[0]-kl_pos[0])
    dis_y = (center[1]-kl_pos[1])
    #here, before shift to get value
    dist=np.sqrt(dis_x*dis_x+dis_y*dis_y)

    if dist <= r_break:
        spat_filter_value = spatialFilter(dist,0,sigma,alpha,beta)
    else:
        spat_filter_value = 0.
    
    return spat_filter_value
    
    #if f==1:
    #    print 'dist: '+str(dist) +' dis_x: '+str(dis_x)+' dis_y: '+str(dis_y)+' dis_x_gs: '+str(dis_x_gs)+' dis_y_gs: '+str(dis_y_gs)+ ' kl: ' +str(kl_val)+' center: ' +str(center) + ' gs: ' +str(grid_shift)
        
    #make it a real circle
    '''
    if dist <= r_break:
        #print int(np.rint(1000*dist))
        add_val = mult_val*filter_val_list[int(np.rint(1000*dist))]#*spatialFilter(dist,0,sigma,alpha,beta)
    return add_val
    '''

spat_filt_values = [[[] for j in range(midget_width)] for i in range(midget_height)]
'''
def getSpatFilter(ij):
    
    pos_i = 0.866*ij[0]        #2*cos(30deg)
    pos_j = ij[1]
    
    if ij[0]%2 == 0:
        pos_j += 0.5  #cos(30deg)


    
    i_low = int(pos_i*px_midget_ratio-spat_filter_break_radius) # -> all j with dist < r
    i_ceil = int(pos_i*px_midget_ratio+spat_filter_break_radius)+1
    j_low = int(pos_j*px_midget_ratio-spat_filter_break_radius) # -> all j with dist < r
    j_ceil = int(pos_j*px_midget_ratio+spat_filter_break_radius)+1

    #pos_i=midget_dist*pos_i
    #pos_j=midget_dist*pos_j
    
    #right, since in opencv (height, width)
    pos=(pos_i,pos_j)
    
    if f==0:
        grid[0] += [ij[0]]
        grid[1] += [ij[1]]
        m_grid[0] += [pos_j]
        m_grid[1] += [pos_i]
        midget_grid[ij[0]][ij[1]] = (pos_j, pos_i) 
        
        #create list with numpy arrays here
        #define length + list beforehand
        
        for k in range(i_low,i_ceil):
            for l in range(j_low,j_ceil):
                spat_filt_values[ij[0]][ij[1]] += [spatFilterPx(pos,(k,l),spat_filter_break_radius,sigma,alpha,beta)]
        #print spat_filt_values[ij[0]][ij[1]]

    pxiq = []
    for k in range(i_low,i_ceil):
        for l in range(j_low,j_ceil):
            #print (midget_height+k)%midget_height,(midget_width+l)%midget_width
            pxiq += [pixels3d[(height+k)%height][(width+l)%width]/255.]
   
    
    #exchange the order of that, to make it faster
    #spatial filters can allow negative potentials
    midget_pixels4d[ij[0]][ij[1]][f] = sum(itertools.imap(lambda x,y: x*y, spat_filt_values[ij[0]][ij[1]],pxiq))
'''

def getSpatFilter(ij):
    
    if ij[0]%5==0:
        if ij[1]==0:
            print ij[0]
            
    pos_i = 0.866*ij[0]        #2*cos(30deg)
    pos_j = ij[1]
    
    if ij[0]%2 == 0:
        pos_j += 0.5  #cos(30deg)


    
    i_low = int(pos_i*px_midget_ratio-spat_filter_break_radius) # -> all j with dist < r
    i_ceil = int(pos_i*px_midget_ratio+spat_filter_break_radius)+1
    j_low = int(pos_j*px_midget_ratio-spat_filter_break_radius) # -> all j with dist < r
    j_ceil = int(pos_j*px_midget_ratio+spat_filter_break_radius)+1

    #pos_i=midget_dist*pos_i
    #pos_j=midget_dist*pos_j
    
    #right, since in opencv (height, width)
    pos=(pos_i,pos_j)
    
    #if f==0:
    grid[0] += [ij[0]]
    grid[1] += [ij[1]]
    m_grid[0] += [pos_j]
    m_grid[1] += [pos_i]
    midget_grid[ij[0]][ij[1]] = (pos_j, pos_i) 
    
    #create list with numpy arrays here
    #define length + list beforehand
    
    mp=[0. for g in range(frame_number)]
    q_len = abs((j_ceil-j_low)*(i_ceil-i_low))
    pxiq = []#0 for q in range(q_len)]
    
    for k in range(i_low,i_ceil):
        for l in range(j_low,j_ceil):
            spat_filt_values[ij[0]][ij[1]] += [spatFilterPx(pos,(k,l),spat_filter_break_radius,sigma,alpha,beta)]
            pxiq+=[pixels4d[(height+k)%height][(width+l)%width]]

    #exchange the order of that, to make it faster
    #spatial filters can allow negative potentials
    #print mp
    for q in xrange(q_len):
        w = spat_filt_values[ij[0]][ij[1]][q]
        muli = [o/255.*w for o in pxiq[q]]
        mp = list(itertools.imap(lambda w,z: w+z, muli,mp))
    #print mp
    midget_pixels4d[ij[0]][ij[1]] = mp
    #= sum(itertools.imap(lambda x,y: x*y, spat_filt_values[ij[0]][ij[1]],pxiq))
    

#apply the spatial filter 
#for f in range(frame_number):
#    print f
#    pixels3d = [[item[f] for item in pxst] for pxst in pixels4d]      
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

#plt.plot(midget_pixels4d[20])
#plt.show()
#--------------------------------------------------------------------------------------MIDGETS-TEMPORAL-FILTERS
print '---'
print 'start of temporal part'
#now apply the temporal filters: output are two lists for on/off fields
for i in range(midget_height):
    if i%5==0:
        print i
    for j in range(midget_width):
        temp_midget_px4d = []
        pop = temp_midget_px4d.pop
        for f in range(frame_number):
            temp_midget_px4d.insert(0, float(midget_pixels4d[i][j][f]))
            if f > 200:
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

print '---'
print 'start of parasol part'
parasol_grid= [[(0,0) for j in range(parasol_width)] for i in range(parasol_height)]

par_values = [[[0. for f in range(frame_number)] for j in range(parasol_width)] for i in range(parasol_height)] 
par_values_on = [[[0. for f in range(frame_number)] for j in range(parasol_width)] for i in range(parasol_height)] 
par_values_on_off = [[[0. for f in range(frame_number)] for j in range(parasol_width)] for i in range(parasol_height)] 

p_spat_filt_values = [[[] for j in range(midget_width)] for i in range(midget_height)]

def spatFilterParasol(center,kl,r_break,sigma,alpha,beta):

    klv_w = 0
    klv_h = 0
    new_kl0,new_kl1 = kl
    
    kl_val = midget_grid[(midget_height+kl[0])%midget_height][(midget_width+kl[1])%midget_width]
    n_kl_h=kl_val[0]
    n_kl_w=kl_val[1]
    
    #if that seems strange: midget_grid is transposed in order to show it better in graphics etc.
    if kl[0] >= midget_height:
        n_kl_w=kl_val[1]+height
    if kl[0] < 0:
        n_kl_w=kl_val[1]-height
    if kl[1] >= midget_width:
        n_kl_h=kl_val[0]+width
    if kl[1] < 0:
        n_kl_h=kl_val[0]-width
    
    kl_pos=(n_kl_h,n_kl_w)
    
    dis_x = (center[0]-kl_pos[0])
    dis_y = (center[1]-kl_pos[1])
    #here, before shift to get value
    dist=np.sqrt(dis_x*dis_x+dis_y*dis_y)

    if dist <= r_break:
        spat_filter_value = spatialFilter(dist,0,sigma,alpha,beta)
    else:
        spat_filter_value = 0.
    
    return spat_filter_value

''''
def getSpatFilterParasol(ij):
    
    pos_i = 0.866*par_m_ratio*ij[0]        #2*cos(30deg)
    pos_j = par_m_ratio*ij[1]
    
    if ij[0]%2 == 0:
        pos_j += par_m_ratio*0.5 

    move = (par_m_ratio)*spat_filter_break_radius/px_midget_ratio
    
    i_low = int(par_m_ratio*ij[0]-move/0.866) # -> all j with dist < r
    i_ceil = int(par_m_ratio*ij[0]+move/0.866)+1  
    j_low = int(par_m_ratio*ij[1]-move) # -> all j with dist < r
    j_ceil = int(par_m_ratio*ij[1]+move)+1
    
    pos = (pos_j, pos_i) 
    
    #print j_low,j_ceil
    #if f==0:
    pgrid[0] += [pos_j]
    pgrid[1] += [pos_i]
    parasol_grid[ij[0]][ij[1]] = pos 
        
        
    mp=[0. for g in range(frame_number)]
    q_len = abs((j_ceil-j_low)*(i_ceil-i_low))
    pxiq = []#0 for q in range(q_len)]
    
    for k in range(i_low,i_ceil):
        for l in range(j_low,j_ceil):
            spat_filt_values[ij[0]][ij[1]] += [spatFilterPx(pos,(k,l),spat_filter_break_radius,sigma,alpha,beta)]
            pxiq+=[pixels4d[(height+k)%height][(width+l)%width]]

    print 'creation of spatial filter finished'
    #exchange the order of that, to make it faster
    #spatial filters can allow negative potentials
    #print mp
    for q in xrange(q_len):
        if q % 5==0:
            print q
        w = spat_filt_values[ij[0]][ij[1]][q]
        muli = [o/255.*w for o in pxiq[q]]
        mp = map(lambda w,z: w+z, muli,mp)
    #print mp
    midget_pixels4d[ij[0]][ij[1]] = mp
    
    
        dis = False
        for k in range(i_low,i_ceil):
            for l in range(j_low,j_ceil):
                p_spat_filt_values[ij[0]][ij[1]] += [spatFilterParasol(pos,(k,l),(par_m_ratio)*spat_filter_break_radius,par_m_ratio*sigma,alpha,beta)]
    
    mgiq = []
    for k in range(i_low,i_ceil):
        for l in range(j_low,j_ceil):
            #print (midget_height+k)%midget_height,(midget_width+l)%midget_width
            mgiq += [midgets3d_on[(midget_height+k)%midget_height][(midget_width+l)%midget_width]]
        
    #spatial filters can allow negative potentials
    par_values_on[ij[0]][ij[1]][f] = sum(itertools.imap(lambda x,y: x*y, p_spat_filt_values[ij[0]][ij[1]],mgiq))
'''    
    
def getSpatFilterParasol(ij):
    
    if ij[0]%5==0:
        if ij[1]==0:
            print ij[0]
    
    pos_i = 0.866*par_m_ratio*ij[0]        #2*cos(30deg)
    pos_j = par_m_ratio*ij[1]
    
    if ij[0]%2 == 0:
        pos_j += par_m_ratio*0.5 

    move = (par_m_ratio)*spat_filter_break_radius/px_midget_ratio
    
    i_low = int(par_m_ratio*ij[0]-move/0.866) # -> all j with dist < r
    i_ceil = int(par_m_ratio*ij[0]+move/0.866)+1  
    j_low = int(par_m_ratio*ij[1]-move) # -> all j with dist < r
    j_ceil = int(par_m_ratio*ij[1]+move)+1
    
    pos = (pos_j, pos_i) 
    
    #print j_low,j_ceil
    #if f==0:
    pgrid[0] += [pos_j]
    pgrid[1] += [pos_i]
    parasol_grid[ij[0]][ij[1]] = pos 
        
        
    pp=[0. for g in range(frame_number)]
    q_len = abs((j_ceil-j_low)*(i_ceil-i_low))
    mgiq = []#0 for q in range(q_len)]
    
    #to_graph = [[0. for l in range(j_ceil-j_low)] for k in range(i_ceil-i_low)]
    for k in range(i_low,i_ceil):
        for l in range(j_low,j_ceil):
            p_spat_filt_values[ij[0]][ij[1]] += [spatFilterParasol(pos,(k,l),(par_m_ratio)*spat_filter_break_radius,par_m_ratio*sigma,alpha,beta)]
            mgiq += [temp_filter_vals_on[(midget_height+k)%midget_height][(midget_width+l)%midget_width]]
            #if ij[0]==0:
            #    if ij[1]==16:
            #        to_graph[k-i_low][l-j_low]=spatFilterParasol(pos,(k,l),(par_m_ratio)*spat_filter_break_radius,par_m_ratio*sigma,alpha,beta) 
    
    '''
    #leave this for possible further investigation later on
    if ij[0]==0:
        if ij[1]==16:
            fig = plt.figure(1)

            ax = fig.add_subplot(111)
            ax.set_title('parasolic filter')
            plt.imshow(to_graph, aspect='auto', interpolation='nearest')
            ax.set_aspect('equal')
            plt.axis('off')

            cax = fig.add_axes([0.,0.,1.,1.])
            cax.get_xaxis().set_visible(False)
            cax.get_yaxis().set_visible(False)
            cax.patch.set_alpha(0)
            cax.set_frame_on(False)
            
            plt.show()
            plt.close()
    '''
    #exchange the order of that, to make it faster
    #spatial filters can allow negative potentials
    for q in xrange(q_len):
        w = p_spat_filt_values[ij[0]][ij[1]][q]
        muli = [o*w for o in mgiq[q]]
        pp = list(itertools.imap(lambda w,z: w+z, muli, pp))
    par_values_on[ij[0]][ij[1]] = pp  


print 'midget width: ' + str(int(midget_width))
print 'parasol width: ' + str(int(parasol_width))
#apply the spatial filter 
#for f in range(frame_number):
#print f
#midgets3d = [[mi[f] for mi in midgs] for midgs in temp_filter_vals]   
#midgets3d_on = [[mi[f] for mi in midgs] for midgs in temp_filter_vals_on]
#midgets3d_on_off = [[mi[f] for mi in midgs] for midgs in temp_filter_vals_off]  
#for i in range(len(temp_filter_vals)):
#    for j in range(len(temp_filter_vals[0])):
#        midgets3d_on_off[i][j] += temp_filter_vals_on[i][j][f]
map(lambda x: getSpatFilterParasol(x), itertools.product(range(parasol_height), range(parasol_width)))

print 'parasols done'

m_pos = np.asarray(midget_grid)
m_pos_data = open('/home/schrader/Documents/microsaccades/data/'+sim_title+'/m_pos_'+str(handle_name)+'.data','w+')
np.save(m_pos_data, m_pos)
m_pos_data.close()

p_pos = np.asarray(parasol_grid)
p_pos_data = open('/home/schrader/Documents/microsaccades/data/'+sim_title+'/p_pos_'+str(handle_name)+'.data','w+')
np.save(p_pos_data, p_pos)
p_pos_data.close()
#plt.plot(pgrid[0][:10],pgrid[1][:10],'go')
#plt.plot(pgrid[0][120:130],pgrid[1][120:130],'go')
#plt.plot(pgrid[1],pgrid[0],'go')
#plt.plot(pgrid[0][55:60],pgrid[1][55:60],'go')
#plt.plot(pgrid[0][105:110],pgrid[1][105:110],'go')
#plt.plot(pgrid[0][75:82],pgrid[1][75:82],'go')
#plt.plot(pgrid[1],pgrid[0],'go')
#plt.savefig('img/grid_large_new.pdf')
#plt.show()

#--------------------------------------------------------------------------------------------------------OUTPUT
p_output = np.asarray(par_values)
p_output_on = np.asarray(par_values_on)
#p_output_on_off = np.asarray(par_values_on_off)

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
m_data = open('/home/schrader/Documents/microsaccades/data/'+sim_title+'/midget_rates_'+str(handle_name)+'.data','w+')
np.save(m_data, m_output)
m_data.close()

p_data = open('/home/schrader/Documents/microsaccades/data/'+sim_title+'/parasolic_rates_'+str(handle_name)+'.data','w+')
np.save(p_data, p_output)
p_data.close()

m_data = open('/home/schrader/Documents/microsaccades/data/'+sim_title+'/midget_rates_'+str(handle_name)+'_on.data','w+')
np.save(m_data, m_output_on)
m_data.close()

p_data = open('/home/schrader/Documents/microsaccades/data/'+sim_title+'/parasolic_rates_'+str(handle_name)+'_on.data','w+')
np.save(p_data, p_output_on)
p_data.close()

#not needed for poletti, murakami etc.
#m_data = open('/home/schrader/Documents/microsaccades/data/'+sim_title+'/midget_rates_'+str(handle_name)+'_off.data','w+')
#np.save(m_data, m_output_off)
#m_data.close()

print '--- PROGRAM FINISHED ---'
sys.exit()

#p_data = open('data/'+sim_title+'/parasolic_rates_'+str(handle_name)+'_on_off.data','w+')
#np.save(p_data, p_output_on_off)
#p_data.close()

#-----------------------------------------------------------------------------------------------PLOT-SOME-STUFF

'''
fig = plt.figure(1)

ax = fig.add_subplot(211)
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

ax = fig.add_subplot(212)
ax.set_title('parasolic output')
plt.imshow(p_output_on[:,:,250], aspect='auto', interpolation='nearest')
ax.set_aspect('equal')
plt.axis('off')

cax = fig.add_axes([0.,0.,1.,.5])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)

out= '/home/schrader/Documents/microsaccades/img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '_' + str(handle_name) + '.pdf'
plt.savefig(out)
plt.close()
#plt.show()



#additional for moment, delete later on again
fig = plt.figure(1)

ax = fig.add_subplot(221)
ax.set_title('midget output _on')
plt.imshow(m_output_on[:,:,250], aspect='auto', interpolation='nearest')
ax.set_aspect('equal')
plt.axis('off')

cax = fig.add_axes([0.,0.,1.,1.])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
#plt.colorbar(orientation='vertical')

ax = fig.add_subplot(223)
ax.set_title('parasolic output on')
plt.imshow(p_output_on[:,:,250], aspect='auto', interpolation='nearest')
ax.set_aspect('equal')
plt.axis('off')

cax = fig.add_axes([0.,0.,1.,.5])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)

out= '/home/schrader/Documents/microsaccades/img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '_' + str(handle_name) + '_on.pdf'
plt.savefig(out)

#out= 'img/video/' + str(sim_title) +'/'+ str(handle_name) + '.pdf'
#plt.savefig(out)

#plt.show()
'''