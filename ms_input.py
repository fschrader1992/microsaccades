#THIS PART READS THE PIXELS OF A MOVIE AND APPLIES THE TEMPORAL FILTERS TO THE PROPOSED MODEL FOR MICROSACCADES. THE OUTPUT ARE THE POTENTIAL VALUES FOR THE POISSON RATES.
'''
loading the video into an array
-> first get all the pixels in a frame (grayscale)
-> after that calculate the values of each tempral flter at each time and store them in another array which will be the basis for changing poisson rates
'''

import pylab as pyl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cv2
import datetime
import itertools 
from microsaccades_functions import *

#-------------------------------------------------------------------------------------------LOAD-VIDEO
#load video
cap = cv2.VideoCapture('/video/test1.avi')#opposite_1fr0deg.mp4')
#just to be sure
while not cap.isOpened():
    cap = cv2.VideoCapture('video/test1.avi')#opposite_1fr0deg.mp4')
    cv2.waitKey(1000)
    print "Wait for the header"
    
#data of video
width = int(cap.get(3))
height = int(cap.get(4))
#dt=1/video_fps
dt = 1/cap.get(5)
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

#ffile =open('pixels4d.data','w+')
#ffile.write(str(pixels4d[9][9]))
#ffile.close

#-----------------------------------------------------------------------------------------TIME-FILTERS
#now comes the part of the temporal filters
#later: define list with all different tau1, tau2 and p, for the moment one is enough (for on/off-center) 


#problem here is now, that both create the same output and therefore nothing is produced
off_tau1 = .5 #.05 --> what is dt in Garrett's code? 0.1ms or 1ms? depending on that, this is the 
off_tau2 = 1.5 #.1     right value
off_p = .8 #.05

on_tau1 = .5 #.2
on_tau2 = 1.5 #.1
on_p = .8 #.05

#set values for the spatial filters for the pixel in \mum(?)
alpha = .1
beta = 5 #CSR
sigma = 1 # parasolic/midget cell ratio * sigma (=1)
px_rec_ratio = 3. #pixel to receptor ratio, needed for receptor distance/number
px_dist= 1./px_rec_ratio #set rec distance to one, if changes, also change for parasolic spatial filter
par_m_ratio = 2

#spatial filter values determined
spat_filter_break_radius = 6 #filter radius in px 

#subsequent values
receptor_dist = 1. #px_rec_ratio*px_dist #for the moment, later maybe hexagonal?

#we actually need to calculate the values of each temporal filter just once
temp_filter_on = tempFilter(200,dt,on_tau1,on_tau2,on_p) #mayber 120, compare to paper --> 200 in Garrett's code!      


rec_height = int(np.floor((height-spat_filter_break_radius)/(1.732*px_rec_ratio))) #2*cos(30deg)
rec_width = int(np.floor((width-spat_filter_break_radius)/(3*px_rec_ratio))*4) #3 is difference between hexagons where pattern repeatesand 4 is the number of cells in that distance

print rec_height,rec_width

rec_pixels4d = np.zeros(shape=(rec_height-1,rec_width,frame_number))

grid=[[0],[0]]
pgrid=[[0],[0]]

midget_grid= [[(0,0) for j in range(rec_width)] for i in range(rec_height)]

def spatFilterPx(mult_list,center,kl,r_break,sigma,alpha,beta):
    add_val = 0
    mult_val = mult_list[kl[0]][kl[1]]
    kl = float(kl[0])*px_dist,float(kl[1])*px_dist
    dist = np.sqrt(sum(itertools.imap(lambda x,y: (x-y)*(x-y), center, kl)))
    #make it a real circle
    if dist <= r_break:
        add_val = mult_val*spatialFilter(dist,0,sigma,alpha,beta)
    return add_val

def getSpatFilter(ij):
    
    pos_i = 1.732*(ij[0]+1)        #2*cos(30deg)
    pos_j = np.floor(float(ij[1])/4)*px_rec_ratio + 2
    
    if ij[1]%4 == 0:
        pos_i += 0.866  #cos(30deg)
    elif ij[1]%4 == 1:
        pos_j += .5
    elif ij[1]%4 == 2:
        pos_j += 1.5
    elif ij[1]%4 == 3:
        pos_j += 2
        pos_i += 0.866  #cos(30deg)
        
    i_low = int(pos_i*px_rec_ratio-spat_filter_break_radius+1) # -> all j with dist < r
    i_ceil = i_low + int(2*spat_filter_break_radius)
    j_low = int(pos_j*px_rec_ratio-spat_filter_break_radius+1) # -> all j with dist < r
    j_ceil = j_low + int(2*spat_filter_break_radius)
    
    if f==0:
        grid[0] += [pos_j]
        grid[1] += [pos_i]
        midget_grid[ij[0]][ij[1]] = (pos_i, pos_j)
    
    rec_pixels4d[ij[0]][ij[1]][f] = sum(itertools.imap(lambda x: spatFilterPx(pixels3d,(pos_i,pos_j),x,spat_filter_break_radius,sigma,alpha,beta), itertools.product(range(i_low,i_ceil),range(j_low,j_ceil))))


#apply the spatial filter 
for f in range(frame_number):
    print f
    pixels3d = [[item[f] for item in pxst] for pxst in pixels4d]      
    map(lambda x: getSpatFilter(x), itertools.product(range(rec_height-2),range(rec_width-2)))
     
plt.plot(grid[0],grid[1],'bo')
#plt.show()

#output frame for the filter values at the end
temp_filter_vals_on  = np.zeros(shape=(rec_height,rec_width,frame_number))
#temp_filter_vals_off = np.zeros(shape=(rec_height,rec_width,frame_number))

 
print '---'
#now apply the temporal filters: output are two lists for on/off fields
for i in range(rec_height-2):
    print i
    for j in range(rec_width-2):
        temp_rec_px4d = []
        pop = temp_rec_px4d.pop
        for f in range(frame_number):
            temp_rec_px4d.insert(0, float(rec_pixels4d[i][j][f]))
            if f > 200:
                pop()
            #add a new entry to the time list of the pixel i,j
            temp_filter_vals_on[i][j][f] = rec_pixels4d[i][j][f]#sum(itertools.imap(lambda x,y: x*y, temp_rec_px4d, temp_filter_on))
            #temp_filter_vals_off[i][j][f] = sum(itertools.imap(lambda x,y: x*y, temp_rec_px4d, temp_filter_off))


#calculate the difference between surround and center fields --> needs to be done? and safe to file
m_output = np.asarray(temp_filter_vals_on)#np.subtract(temp_filter_vals_on, temp_filter_vals_off), dtype=int)


#-------------------------------------------------------------------------------------PARASOLIC-OUTPUT

par_values = [[[0. for f in range(frame_number)] for j in range(int(rec_width/4)-8)] for i in range(int(rec_height/2)-4)]
              
#range(int(rec_height/par_m_ratio)-2*par_m_ratio)]

#the spatial filters: use receptor difference + 4*sigma + shift to (.5,.5) in middle of gridcell formed by midgets
par_spat_filter_values = [[spatialFilter(np.sqrt((float(i)-.5)*(float(i)-.5)+(float(j)-.5)*(float(j)-.5)),0,par_m_ratio*sigma,alpha,beta) for j in range(-4,4)] for i in range(-4,4)] #change also, if hexagonal + since receptor distance equals 1

'''
#apply the spatial filter
for f in range(frame_number):
    for i in range(int(rec_height/par_m_ratio)-2*par_m_ratio):
        i_f = par_m_ratio*i+2*par_m_ratio
        for j in range(int(rec_width/par_m_ratio)-2*par_m_ratio):      
            j_f = par_m_ratio*j+2*par_m_ratio
            m_val=0
            for q in range(-4,4): #range(-spat_filter_break_radius,spat_filter_break_radius):
                for p in range(-4,4): #range(-spat_filter_break_radius,spat_filter_break_radius):
                   m_val += temp_filter_vals_on[i_f+q][j_f+p][f]*par_spat_filter_values[q+4][p+4]
            par_values[i][j][f] = m_val
''' 


def spatFilterParasolPx(mult_list,mean,kl,r_break,sigma,alpha,beta):
    add_val = 0
    mult_val = mult_list[kl[0]][kl[1]]
    kl_val = midget_grid[kl[0]][kl[1]]
    dist = np.sqrt(sum(itertools.imap(lambda x,y: (x-y)*(x-y), mean, kl_val)))
    #make it a real circle
    if dist <= r_break:
        add_val = mult_val*spatialFilter(dist,0,sigma,alpha,beta)
    return add_val


def getSpatFilterParasol(ij):
    
    
    pos = [midget_grid[2*ij[0]][4*ij[1]+1][1] + .866, midget_grid[2*ij[0]][4*ij[1]+1][0] + .5]
    #print 'ij' + str(ij)

    i_low = int(2*ij[0]-par_m_ratio*spat_filter_break_radius/px_rec_ratio+1) # -> all j with dist < r
    i_ceil = i_low + int(2*par_m_ratio*spat_filter_break_radius/px_rec_ratio)    
    j_low = int(4*ij[1]+1-par_m_ratio*spat_filter_break_radius/px_rec_ratio) # -> all j with dist < r
    j_ceil = j_low + int(2*par_m_ratio*spat_filter_break_radius/px_rec_ratio)
    
    #print j_low,j_ceil
    if f==0:
        pgrid[0] += [pos[0]]
        pgrid[1] += [pos[1]]
    
    par_values[ij[0]][ij[1]][f] = sum(itertools.imap(lambda x: spatFilterParasolPx(midgets3d,pos,x,spat_filter_break_radius,par_m_ratio*sigma,alpha,beta), itertools.product(range(i_low,i_ceil),range(j_low,j_ceil))))
    
    #par_m_ratio*sigma
 
#for i in itertools.product(range(int(rec_width/4)-8),range(int(rec_height/2)-4)):
#    print i

#apply the spatial filter 
for f in range(frame_number):
    print f
    midgets3d = [[mi[f] for mi in midgs] for midgs in temp_filter_vals_on]      
    map(lambda x: getSpatFilterParasol(x), itertools.product(range(2,int(rec_height/2)-4),range(2,int(rec_width/4)-8)))

p_output = np.asarray(par_values)

plt.plot(pgrid[0],pgrid[1],'ro')
plt.show()

#------------------------------------------------------------------------------------------SAVE-OUTPUT
m_data = open('data/midget_values.data','w+')
np.save(m_data, m_output)
m_data.close()

p_data = open('data/parasolic_values.data','w+')
np.save(p_data, p_output)
p_data.close()

#--------------------------------------------------------------------------------------PLOT-SOME-STUFF

fig = plt.figure(1)

ax = fig.add_subplot(221)
ax.set_title('midget output')
plt.imshow(m_output[:,:,1], interpolation='nearest')
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
plt.imshow(p_output[:,:,1], interpolation='nearest')
ax.set_aspect('equal')
plt.axis('off')

cax = fig.add_axes([0.,0.,1.,1.])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)


now = datetime.datetime.now()

out= 'img/'+ str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/output_mp_test_' + str(now.hour) + '_' + str(now.minute) + '_' + str(now.second) + '.pdf'
plt.savefig(out)

out= 'img/output_mp_test.pdf'
plt.savefig(out)

plt.show()