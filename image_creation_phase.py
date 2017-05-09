#!/bin/python
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import cv2
from microsaccades_functions import *

#---------------------------------------------------------------------------------------IMAGE-CREATION

framerate = 3. #get some fancy frequency calculation
stripe_width = 15
gap = 15
image_size = 560 # two times because of rotation
degrees = 0

degrees = sys.argv[1]

file_location = "video/img_input/phases_0fr" + str(degrees) + "deg"

degrees = float(degrees)

film_length = 1000

canvas = np.zeros((image_size, image_size))
current_col = 0
while current_col < image_size:
    if current_col + stripe_width + gap  <= image_size-1:
        canvas[:, current_col:current_col+stripe_width] = 1
        current_col += stripe_width + gap
    elif current_col + stripe_width <= image_size-1:
        canvas[:, current_col:current_col+stripe_width] = 1
        current_col = image_size
    else:
        canvas[:, current_col:] = 1
        current_col = image_size

fig = plt.figure()
fig.set_size_inches(1, 1)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas, cmap='gray')
plt.savefig(file_location + "/first.png",  dpi = 560)

img = cv2.imread(file_location + "/first.png",0)
rows,cols = img.shape

d_file = open('data/phase_displacement.data','r+')
disp = np.load(d_file)   
d_file.close()

for f in range(film_length):
    #------------------------------------------------------------------------------NORMAL-DISPLACEMENT
    fig = plt.figure()
    fig.set_size_inches(2, 2)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(canvas, cmap='gray')
    
    transl = np.float32([[1,0,int(disp[f])],[0,1,0]])
    tlFig = cv2.warpAffine(img,transl,(cols,rows))
    
    '''
    #save to file, comment for additional rotation
    tlFig = tlFig[160:400,160:400]
    fig.set_size_inches(1, 1)
    
    plt.imshow(tlFig,cmap='gray')
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 240)
    plt.close()
    '''
    #-----------------------------------------------------------------------------------------ROTATION
    
    rot = cv2.getRotationMatrix2D((cols/2.,rows/2.),degrees,1)
    rotFig = cv2.warpAffine(tlFig,rot,(cols,rows))
    
    #save to file
    rotFig = rotFig[160:400,160:400]
    fig.set_size_inches(1, 1)
    
    plt.imshow(rotFig,cmap='gray')
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 240)
    plt.close()
    