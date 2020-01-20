import numpy as np
import matplotlib.pyplot as plt
import cv2
from microsaccades_functions import *
import sys
import os

plt.switch_backend('agg')


#---------------------------------------------------------------------------------------IMAGE-CREATION

image_size = 320 #30 #280
image_height = 30 #280 #120
radius = 4. #dot size of 4arcmin -> 8px
rect_size = 110. #size of inner radius 
vel = 0.06 #velocity is 30 arcmin/s -> 60px/1000ms

suff = sys.argv[1]
vel_on = int(sys.argv[2])
angle = sys.argv[3]
rem_on = int(sys.argv[4])
rect = int(sys.argv[5])
cem_on = int(sys.argv[6])
cem_l_on = int(sys.argv[7])

#for 1D
angle = 0.
direction = float(angle)*np.pi/8. #dot moving direction in arc (get right velocities !) not anymore!

#for normal distributed microsaccades
sigma = 0.447

file_location = "/home/schrader/Documents/microsaccades/video/img_input/poletti2010/"+str(suff)
#if vel_on==1:
#    file_location += "/"+str(angle)+"_8_pi_arc"
#save the different conditions

film_length = 600 #int(framerate)*(stripe_width+gap)

#get normal 2d distribution
#use this to get same statistics as for actual 2D random walk -> just take distance 
mean = [0, 0]
cov = [[1, 0], [0, 1]]
tl_x, tl_y = sigma*np.random.multivariate_normal(mean, cov, film_length).T

rm = tl_x #np.sign(tl_x)*np.sqrt(tl_x*tl_x+tl_y*tl_y)

d_data = open(file_location+'/displacement.data','w+')
np.save(d_data, (tl_x, tl_y))
d_data.close()
pos = (image_height/2.-0.5, image_size/2.-0.5)
center = (image_height/2.-0.5,image_size/2.-0.5) #for rectangle

#loop through the 
for f in range(film_length):
    canvas = np.zeros((image_height, image_size))
    
    #update position
    #for dot moving task-----------------------
    if vel_on==1:
        pos = (pos[0]+vel*np.sin(direction),pos[1]+vel*np.cos(direction))
    #------------------------------------------
	
    #for random eye movements------------------
	#REPLACED FOR 1D!
    if rem_on==1:
        #pos = (pos[0]+tl_y[f],pos[1]+tl_x[f])
        pos = (pos[0],pos[1]+rm[f])
    #------------------------------------------
	
    #for moving rectangle frame----------------
    if cem_l_on==1:
        center = (center[0]+vel*np.sin(direction),center[1]+vel*np.cos(direction))
    #------------------------------------------
    
    #for randomly moving rectangle frame-------
    if cem_on==1:
        #center = (center[0]+tl_y[f],center[1]+tl_x[f])
        center = (center[0],center[1]+rm[f])
    #------------------------------------------
	
    
    for i in range(image_height):
        for j in range(image_size):
            dist = np.sqrt((float(i)-pos[0])*(float(i)-pos[0])+(float(j)-pos[1])*(float(j)-pos[1]))
            if dist <= 3.5:
                canvas[i,j] = 1
            elif 3.5 < dist and dist < 4.5:
                canvas[i,j] = 4.5 - dist
            #rectangle-(exp2)--------------------------
            if rect==1:
                if abs(float(i)-center[0]) > rect_size or abs(float(j)-center[1])> rect_size:
                    canvas[i,j] = 1
            #------------------------------------------
     
    '''
    fig = plt.figure()
    fig.set_size_inches(2,2)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(canvas, cmap='gray')
    plt.savefig(file_location + "/first"+str(f+1).zfill(3)+".png",  dpi = 140)
    plt.close()
    '''
    fig = plt.figure()
    #fig.set_size_inches(0.25,1)
    fig.set_size_inches(2,30./160.)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(canvas, cmap='gray')
    
    plt.savefig(file_location + "/first"+str(f+1).zfill(3)+".png",  dpi = 160) #120)
    plt.close()
    
    img = cv2.imread(file_location + "/first"+str(f+1).zfill(3)+".png",0)
    rows,cols = img.shape
    #--------------------------------------------------------------------------------GAUSSIAN-BLURRING
    img = cv2.imread(file_location + "/first"+str(f+1).zfill(3)+".png",0)
    rows,cols = img.shape
    
    blur = cv2.GaussianBlur(img,(3,3),0.5)
    '''
    fig = plt.figure()
    fig.set_size_inches(2,2)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(blur, cmap='gray')
    '''
    
    fig = plt.figure()
    fig.set_size_inches(2,30./160.)
    #fig.set_size_inches(0.25,1)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(blur, cmap='gray')
    
    #plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 120)
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 160)
    plt.close()