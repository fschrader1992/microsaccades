#THIS PART READS THE PIXELS OF A MOVIE AND APPLIES THE TEMPORAL FILTERS TO THE PROPOSED MODEL FOR MICROSACCADES. 
'''
loading the video into an array
-> first get all the pixels in a frame (grayscale)
-> after that calculate the values of each tempral flter at each time and store them in another array which will be the basis for changing poisson rates
'''

import pylab as pl
import numpy as np
import cv2
from itertools import imap
from microsaccades_functions import *

#--------------------------------------------------------------------------------------------LOAD-VIDEO
#load video
cap = cv2.VideoCapture('test1.avi')
#just to be sure
while not cap.isOpened():
    cap = cv2.VideoCapture('test1.avi')
    cv2.waitKey(1000)
    print "Wait for the header"
    
#width and height of the frame
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
        '''for i in range(height):
            for j in range(width):
                pixels4d[i][j]+=[gray[i,j]]'''
        map(pixels4d,gray)
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

#------------------------------------------------------------------------------------------TIME-FILTERS
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
temp_filter_off = tempFilter(100,dt,off_tau1,off_tau2,off_p)#200 instead of frame_number --> TEST!
temp_filter_on = tempFilter(100,dt,on_tau1,on_tau2,on_p)

#pl.plot(temp_filter_off)
#pl.plot(temp_filter_on)
#pl.show()

#output frame for the filter values at the end
temp_filter_vals_off = [[[0.0 for f in range(frame_number)] for j in range(width)] for i in range(height)]
temp_filter_vals_on  = [[[0.0 for f in range(frame_number)] for j in range(width)] for i in range(height)]


#now apply the temporal filters
for i in range(height):
    #if i%20 == 0:
    print i
    for j in range(width):
        #write as function later 
        #setTempFilterVals(frame_number)
        temp_px4d = []
        pop = temp_px4d.pop
        for f in range(frame_number):
            temp_px4d.insert(0, float(pixels4d[i][j][f]))
            if f > 100:
                pop()
            #add a new entry to the time list of the pixel i,j
            temp_filter_vals_on[i][j][f] = sum(imap(lambda x,y: x*y, temp_px4d, temp_filter_on))
            temp_filter_vals_off[i][j][f] = sum(imap(lambda x,y: x*y, temp_px4d, temp_filter_off))
            
#ffile = open('temp_filter_vals_on.data','w+')
#ffile.write(str( temp_filter_vals_on))
#ffile.close

pl.plot(temp_filter_vals_off[100][100])
pl.plot(temp_filter_vals_on[100][100])
#tf = []
'''for f in range(frame_number):
    tf += [temp_filter_vals_on[100][100][f]-temp_filter_vals_off[100][100][f]]
pl.plot(tf)
pl.show()'''
#pl.plot(pixels4d[100][100])
pl.show()