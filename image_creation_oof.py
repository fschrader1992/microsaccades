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
image_width = 300

suff = sys.argv[1]
cyc = float(sys.argv[2])
film_length = int(sys.argv[3])

#for normal distributed microsaccades
mu = 0.
sigma = 0.447

file_location = "/home/schrader/Documents/microsaccades/video/img_input/oof/"+str(suff)
#save the different conditions

#get normal 2d distribution
#mean = [0, 0]
#cov = [[1, 0], [0, 1]]
#tl_x, tl_y = sigma*np.random.multivariate_normal(mean, cov, film_length).T
disp = np.random.normal(mu,sigma,400)
#disp = open('data/oof/displacement.data','w+')
#use displacement here?


pos = (image_height/2.+0.5,image_width/2+0.5)
center = (image_height/2.-0.5,image_width/2.-0.5) 

q = np.random.randint(2, size=(int(image_height/3)+1,int(image_width/3)+1))

canvas = np.zeros((image_height, image_width))
current_col = 0
for i in range(image_height):
	for j in range(image_width):
		#for the center background with sinusoidal grating
		canvas[i,j]=0.5*np.sin((cyc/120.*np.pi*float(i)-center[0]))*np.sin((cyc/120.*np.pi*float(j)-center[1]-60))+0.5
		
		#overwrite the rest with random dots of size 3arcmin = 6px
		#also possible to overwrite this with circle
		if j < 90 or j > 210:
			canvas[i,j]=float(q[int(i/3)][int(j/3)])
			
			#dist = np.sqrt((float(i)-center[0])*(float(i)-center[0])+(float(j)-center[1])*(float(j)-center[1]))
			
			#if dist <= circle_width:
			#	canvas[i,j] = 1
			#elif 3.5 < dist and dist  < 4.5:
			#	canvas[i,j] = (1. -4.5 + dist)*canvas[i,j] + 4.5 - dist
			#rectangle-(exp2)--------------------------
			#if rect==1:
			#	if abs(float(i)-center[0]) > rect_size or abs(float(j)-center[1])> rect_size:
			#		canvas[i,j] = 1
			#------------------------------------------
	

print canvas
	
fig = plt.figure()
fig.set_size_inches(2,0.2)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas, cmap='gray')
#plt.show()
#plt.close()

plt.savefig(file_location + "/first.png",  dpi = 150)

img = cv2.imread(file_location + "/first.png",0)
rows,cols = img.shape

d_file = open('data/phase_displacement.data','r+')
disp = np.load(d_file)   
d_file.close()

for f in range(film_length):
    #------------------------------------------------------------------------------NORMAL-DISPLACEMENT
    fig = plt.figure()
    fig.set_size_inches(2,0.25)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    #ax.imshow(canvas, cmap='gray')
    
    transl = np.float32([[1,0,int(disp[f])],[0,1,0]])
    tlFig = cv2.warpAffine(img,transl,(cols,rows))
    
    
    #save to file, comment for additional rotation
    #fig = plt.figure()
    tlFig = tlFig[:,30:270]
    tlFig = cv2.GaussianBlur(tlFig,(3,3),0.5)
    #fig.set_size_inches(2,0.25)
    
    ax.imshow(tlFig,cmap='gray')
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