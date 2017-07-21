import numpy as np
import matplotlib.pyplot as plt
import cv2
from microsaccades_functions import *

#---------------------------------------------------------------------------------------IMAGE-CREATION

image_size = 240
direction = 13.*np.pi/8. #dot moving direction in arc (get right velocities !)
radius = 4. #dot size of 4arcmin -> 8px
rect_size = 110. #size of inner radius 

vel = 0.06 #velocity is 30 arcmin/s -> 60px/1000ms

#for normal distributed microsaccades
sigma = 1

file_location = "video/img_input/poletti2010/exp1/cond3/13_8_pi_arc"
#save the different conditions

film_length = 300 #int(framerate)*(stripe_width+gap)

#get normal 2d distribution
mean = [0, 0]
cov = [[1, 0], [0, 1]]
tl_x, tl_y = sigma*np.random.multivariate_normal(mean, cov, film_length).T


pos = (image_size/2.+0.5,image_size/2+0.5)
center = pos #for rectangle

#loop through the 
for f in range(film_length):
    canvas = np.zeros((image_size, image_size))
    
    #update position
    #for dot moving task-----------------------
    #pos = (pos[0]+vel*np.sin(direction),pos[1]+vel*np.cos(direction))
    #------------------------------------------
	
    #for random eye movements------------------
    #pos = (pos[0]+tl_y[f],pos[1]+tl_x[f])
    #------------------------------------------
	
	#for moving rectangle frame----------------
    #center = (center+tl_y[f],center+tl_x[f])
    #------------------------------------------
    
    for i in range(image_size):
        for j in range(image_size):
            dist = np.sqrt((float(i)-pos[0])*(float(i)-pos[0])+(float(j)-pos[1])*(float(j)-pos[1]))
            if dist <= 3.5:
                canvas[i,j] = 1
            elif 3.5 < dist and dist  < 4.5:
                canvas[i,j] = 4.5 - dist
			#rectangle-(exp2)--------------------------
			#if abs(float(i)-center) > rect_size or abs(float(j)-center)> rect_size:
			#	canvas[i,j] = 1
			#------------------------------------------
        
    fig = plt.figure()
    fig.set_size_inches(2,2)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(canvas, cmap='gray')
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 120)
    plt.close()