import numpy as np
import matplotlib.pyplot as plt
import cv2
from microsaccades_functions import *
import sys
import os

#---------------------------------------------------------------------------------------IMAGE-CREATION

image_size = 280
radius = 4. #dot size of 4arcmin -> 8px
rect_size = 110. #size of inner radius 
vel = 0.06 #velocity is 30 arcmin/s -> 60px/1000ms

suff = sys.argv[1]
vel_on = int(sys.argv[2])
angle = sys.argv[3]
rem_on = int(sys.argv[4])
rect = int(sys.argv[5])
cem_on = int(sys.argv[6])

direction = float(angle)*np.pi/8. #dot moving direction in arc (get right velocities !)

#for normal distributed microsaccades
sigma = 0.447

file_location = "/home/schrader/Documents/microsaccades/video/img_input/poletti2010/"+str(suff)
if vel_on==1:
    file_location += "/"+str(angle)+"_8_pi_arc"
#save the different conditions

film_length = 1000 #int(framerate)*(stripe_width+gap)

#get normal 2d distribution
mean = [0, 0]
cov = [[1, 0], [0, 1]]
tl_x, tl_y = sigma*np.random.multivariate_normal(mean, cov, film_length).T


pos = (image_size/2.+0.5,image_size/2+0.5)
center = (image_size/2.-0.5,image_size/2.-0.5) #for rectangle

#loop through the 
for f in range(film_length):
    canvas = np.zeros((image_size, image_size))
    
    #update position
    #for dot moving task-----------------------
    if vel_on==1:
        pos = (pos[0]+vel*np.sin(direction),pos[1]+vel*np.cos(direction))
    #------------------------------------------
	
    #for random eye movements------------------
    if rem_on==1:
        pos = (pos[0]+tl_y[f],pos[1]+tl_x[f])
    #------------------------------------------
	
    #for moving rectangle frame----------------
    if cem_on==1:
        center = (center[0]+tl_y[f],center[1]+tl_x[f])
    #------------------------------------------
    
    for i in range(image_size):
        for j in range(image_size):
            #for the background use sinusoidal grating
            canvas[i,j]=0.5*np.sin((1./120.*float(i)-center[0]))*np.sin((1./120.*float(j)-center[1]))
            
            dist = np.sqrt((float(i)-pos[0])*(float(i)-pos[0])+(float(j)-pos[1])*(float(j)-pos[1]))
            if dist <= 3.5:
                canvas[i,j] = 1
            elif 3.5 < dist and dist  < 4.5:
                canvas[i,j] = (1. -4.5 + dist)*canvas[i,j] + 4.5 - dist
            #rectangle-(exp2)--------------------------
            if rect==1:
                if abs(float(i)-center[0]) > rect_size or abs(float(j)-center[1])> rect_size:
                    canvas[i,j] = 1
            #------------------------------------------
        
    fig = plt.figure()
    fig.set_size_inches(2,2)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(canvas, cmap='gray')
    plt.savefig(file_location + "/first"+str(f+1).zfill(3)+".png",  dpi = 140)
    plt.close()
    
    #--------------------------------------------------------------------------------GAUSSIAN-BLURRING
    img = cv2.imread(file_location + "/first"+str(f+1).zfill(3)+".png",0)
    rows,cols = img.shape
    
    blur = cv2.GaussianBlur(img,(3,3),0.5)
    
    fig = plt.figure()
    fig.set_size_inches(2,2)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(blur, cmap='gray')
    
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 140)
    plt.close()