import numpy as np
import matplotlib.pyplot as plt
import cv2
#from microsaccades_functions import *
import sys
import os

#---------------------------------------------------------------------------------------IMAGE-CREATION

#radius = 4. #dot size of 4arcmin -> 8px
#rect_size = 110. #size of inner radius 
#vel = 0.06 #velocity is 30 arcmin/s -> 60px/1000ms

#direction = float(angle)*np.pi/8. #dot moving direction in arc (get right velocities !)

image_height = 30
image_width = 120

suff = sys.argv[1]
cyc = float(sys.argv[2])
vel = 0.01*float(sys.argv[3])
#film_length = int(sys.argv[3])
film_length = 400

file_location = "/home/schrader/Documents/microsaccades/video/img_input/mtesting/"+str(suff)
#save the different conditions

pos = (image_height/2.+0.5,image_width/2+0.5)
center = image_width/2. #(image_height/2.-0.5,image_width/2.-0.5) 

canvas = np.zeros((image_height, image_width))
current_col = 0
for f in range(film_length):
    center+=vel
    for j in range(image_width):
            #for the center background with sinusoidal grating
            canvas[:,j]=np.sin((cyc/60.*np.pi*float(j)-center))

    
    canvas = cv2.GaussianBlur(canvas,(3,3),0.5)
    fig = plt.figure()
    fig.set_size_inches(1,0.25)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(canvas, cmap='gray')
    #plt.show()

    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 120)
    plt.close()