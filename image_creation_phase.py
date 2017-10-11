#!/bin/python
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import cv2
from microsaccades_functions import *

#---------------------------------------------------------------------------------------IMAGE-CREATION

framerate = 3. #get some fancy frequency calculation
stripe_width = 16
gap = 16
image_height = 30
image_width = 240
degrees = 0

degrees = sys.argv[1]

file_location = "video/img_input/phases_0fr" + str(degrees) + "deg_32px_only_border"

degrees = float(degrees)

film_length = 1000

canvas = np.zeros((image_height, image_width))
current_col = 0
while current_col < image_width:
    #for only border
    if current_col + stripe_width + gap  <= image_width-1:
        canvas[:, current_col] = 1
        current_col += stripe_width + gap
    else:
        canvas[:, current_col:] = 1
        current_col = image_width
    '''
    #for mixxed grating
    if current_col + stripe_width + gap  <= image_size-1:
        canvas[:, current_col:current_col+stripe_width] = 1
        current_col += stripe_width + gap
    elif current_col + stripe_width <= image_size-1:
        canvas[:, current_col:current_col+stripe_width] = 1
        current_col = image_size
    else:
        canvas[:, current_col:] = 1
        current_col = image_size
    '''
fig = plt.figure()
fig.set_size_inches(2.,0.25)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas, cmap='gray')
plt.savefig(file_location + "/first.png",  dpi = 120)

img = cv2.imread(file_location + "/first.png",0)
rows,cols = img.shape

d_file = open('data/phase_displacement.data','r+')
disp = np.load(d_file)   
d_file.close()

for f in range(film_length):
    #------------------------------------------------------------------------------NORMAL-DISPLACEMENT
    fig = plt.figure()
    fig.set_size_inches(2.33,0.25)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(canvas, cmap='gray')
    
    transl = np.float32([[1,0,int(disp[f])],[0,1,0]])
    tlFig = cv2.warpAffine(img,transl,(cols,rows))
    
    
    #save to file, comment for additional rotation
    tlFig = tlFig[0:30,64:184]
    fig.set_size_inches(1,0.25)
    
    plt.imshow(tlFig,cmap='gray')
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 120)
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
    '''