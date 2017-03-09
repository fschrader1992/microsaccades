#THIS PART READS THE PIXELS OF A MOVIE AND APPLIES THE TEMPORAL FILTERS TO THE PROPOSED MODEL FOR MICROSACCADES. THE OUTPUT ARE THE POTENTIAL VALUES FOR THE POISSON RATES.
'''
loading the video into an array
-> first get all the pixels in a frame (grayscale)
-> after that calculate the values of each tempral flter at each time and store them in another array which will be the basis for changing poisson rates
'''

import pylab as pyl
import numpy as np
import matplotlib.pyplot as plt
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
        #will probably be abundant later on
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        #store 2D array in pixels4d 
        for i in range(height):
            for j in range(width):
                pixels4d[i][j]+=[gray[i,j]]
        if frame_number == 100:
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

#to get bigger numbers, probably change later on 
#dt = 1000*dt

#we actually need to calculate the values of each temporal filter just once
temp_filter_off = tempFilter(100,dt,off_tau1,off_tau2,off_p)#100 instead of frame_number --> TEST!
temp_filter_on = tempFilter(100,dt,on_tau1,on_tau2,on_p)#mayber 120, compare to paper!

#output frame for the filter values at the end
#replace through numpy array?
temp_filter_vals_off = np.zeros(shape=(height,width,frame_number))
temp_filter_vals_on  = np.zeros(shape=(height,width,frame_number))


#now apply the temporal filters: output are two lists for on/off fields
for i in range(height):
    #if i%20 == 0:
    print i
    for j in range(width):
        temp_px4d = []
        pop = temp_px4d.pop
        for f in range(frame_number):
            temp_px4d.insert(0, float(pixels4d[i][j][f]))
            if f > 100:
                pop()
            #add a new entry to the time list of the pixel i,j
            temp_filter_vals_on[i][j][f] = sum(imap(lambda x,y: x*y, temp_px4d, temp_filter_on))
            temp_filter_vals_off[i][j][f] = sum(imap(lambda x,y: x*y, temp_px4d, temp_filter_off))

pyl.plot(temp_filter_vals_off[100][100])
pyl.plot(temp_filter_vals_on[100][100])

#pyl.plot(pixels4d[100][100])
pyl.show()

#calculate the difference between surround and center fields --> needs to be done? and safe to file
temp_filter_subtr = np.asarray(np.subtract(temp_filter_vals_on, temp_filter_vals_off), dtype=int)

tf_data = open('data/temp_filter_subtr.data','w+')
np.save(tf_data, temp_filter_subtr)
tf_data.close()

'''
#plot some output
pyl.plot(temp_filter_subtr[100][100])
pyl.plot(pixels4d[100][100])
pyl.show()'''

fig = plt.figure(figsize=(6, 3.2))

ax = fig.add_subplot(111)
ax.set_title('colorMap')
plt.imshow(temp_filter_subtr[:,:,50])
ax.set_aspect('equal')

cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')
plt.show()