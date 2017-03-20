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
from itertools import imap
from microsaccades_functions import *

#-------------------------------------------------------------------------------------------LOAD-VIDEO
#load video
cap = cv2.VideoCapture('test1.avi')
#just to be sure
while not cap.isOpened():
    cap = cv2.VideoCapture('test1.avi')
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
                pixels4d[i][j]+=[gray[i,j]]
        if frame_number == 200:
            break
        cv2.imshow('frame',gray)
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break
    else:
        break
    
cap.release()
cv2.destroyAllWindows()

ffile =open('pixels4d.data','w+')
ffile.write(str(pixels4d[9][9]))
ffile.close

#-----------------------------------------------------------------------------------------TIME-FILTERS
#now comes the part of the temporal filters
#later: define list with all different tau1, tau2 and p, for the moment one is enough (for on/off-center) 

off_tau1 = .05
off_tau2 = .1
off_p = .05

on_tau1 = .2
on_tau2 = .1
on_p = .05

#set values for the spatial filters for the pixel in \mum(?)
alpha = .1
beta = 0.5
sigma = .1
px_dist= .1
#subsequent values
diag_px_dist = np.sqrt(2)*px_dist
receptor_dist = 3*px_dist #for the moment, later maybe hexagonal?
max_spat_filter_val = spatialFilter(0,0,sigma,alpha,beta)
next_spat_filter_val = spatialFilter(0,px_dist,sigma,alpha,beta)
diag_spat_filter_val = spatialFilter(0,diag_px_dist,sigma,alpha,beta)

#to get bigger numbers, probably change later on 
#dt = 1000*dt

#we actually need to calculate the values of each temporal filter just once
temp_filter_off = tempFilter(100,dt,off_tau1,off_tau2,off_p)#10 instead of frame_number --> TEST!
temp_filter_on = tempFilter(100,dt,on_tau1,on_tau2,on_p)#mayber 120, compare to paper!      

rec_height = int(height/3)
rec_width = int(width/3)

print rec_height,rec_width

rec_pixels4d = np.zeros(shape=(rec_height,rec_width,frame_number))

for f in range(frame_number):
    for i in range(rec_height):
        i_f = i*3
        for j in range(rec_width):
            j_f = j*3
            rec_pixels4d[i][j][f] = max_spat_filter_val*float(pixels4d[i_f][j_f][f]) + next_spat_filter_val*(float(pixels4d[i_f-1][j_f][f])+float(pixels4d[i_f+1][j_f][f])+float(pixels4d[i_f][j_f-1][f])+float(pixels4d[i_f][j_f+1][f])) + diag_spat_filter_val*(float(pixels4d[i_f-1][j_f-1][f])+float(pixels4d[i_f-1][j_f+1][f])+float(pixels4d[i_f+1][j_f-1][f])+float(pixels4d[i_f+1][j_f+1][f]))

#output frame for the filter values at the end
temp_filter_vals_off = np.zeros(shape=(rec_height,rec_width,frame_number))
temp_filter_vals_on  = np.zeros(shape=(rec_height,rec_width,frame_number))

 
#now apply the temporal filters: output are two lists for on/off fields
for i in range(rec_height):
    #if i%20 == 0:
    print i
    for j in range(rec_width):
        temp_rec_px4d = []
        pop = temp_rec_px4d.pop
        #map(frameCalc,range(0,frame_number))
        for f in range(frame_number):
            temp_rec_px4d.insert(0, float(rec_pixels4d[i][j][f]))
            if f > 100:
                pop()
            #add a new entry to the time list of the pixel i,j
            temp_filter_vals_on[i][j][f] = sum(imap(lambda x,y: x*y, temp_rec_px4d, temp_filter_on))
            temp_filter_vals_off[i][j][f] = sum(imap(lambda x,y: x*y, temp_rec_px4d, temp_filter_off))

#pyl.plot(temp_filter_vals_off[10][10])
#pyl.plot(temp_filter_vals_on[10][10])
#pyl.show()

#calculate the difference between surround and center fields --> needs to be done? and safe to file
temp_filter_subtr = np.asarray(np.subtract(temp_filter_vals_on, temp_filter_vals_off), dtype=int)

#poisson_rates=poissonRate(temp_filter_subtr)

tf_data = open('data/temp_filter_subtr.data','w+')
np.save(tf_data, temp_filter_subtr)
tf_data.close()

fig = plt.figure(figsize=(6, 3.2))

ax = fig.add_subplot(111)
ax.set_title('colorMap')
plt.imshow(temp_filter_subtr[:,:,100])
ax.set_aspect('equal')

cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')
plt.show()