#!/bin/python
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import cv2
from microsaccades_functions import *

#---------------------------------------------------------------------------------------IMAGE-CREATION

#this is to test the responsibility to spatial frequencies

freq = sys.argv[1]
half_wl = 120./(float(freq))

stripe_width = int(half_wl)
gap = int(half_wl)
image_size = 240 # two times because of rotation

degrees = sys.argv[1]

file_location = "video/img_input/spatfreq_0fr0deg" + str(freq) + "spat"

film_length = 300              

even = False
#maybe here a loop to make it last even longer
for f in range(film_length):
    #if f%framerate == 0 and f>0: #not yet
    #    offset += 1
    if even == False:
        offset = int(half_wl)
        even = True
    else:
        offset = 0
        even = False

    canvas = np.zeros((image_size, image_size))
    current_col = int(offset) #there's the problem
    if offset < 0:
        canvas[:, 0:current_col+stripe_width] = 1
        #canvas[:, current_col+stripe_width] = (f%int(framerate))/framerate
        current_col += stripe_width + gap
    while current_col < image_size:
        if current_col + stripe_width + gap  <= image_size-1:
            canvas[:, current_col:current_col+stripe_width] = 1
            #canvas[:, current_col] = 1 - (f%int(framerate))/framerate
            #canvas[:, current_col+stripe_width] = (f%int(framerate))/framerate
            current_col += stripe_width + gap
        elif current_col + stripe_width <= image_size-1:
            canvas[:, current_col:current_col+stripe_width] = 1
            #canvas[:, current_col] = 1 - (f%int(framerate))/framerate
            #canvas[:, current_col+stripe_width] = (f%int(framerate))/framerate
            current_col = image_size
        else:
            #canvas[:, current_col] = 1 - (f%int(framerate))/framerate
            canvas[:, current_col:] = 1
            current_col = image_size

    fig = plt.figure()
    fig.set_size_inches(2,2)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(canvas, cmap='gray')
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 120)
    plt.close()