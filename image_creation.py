import numpy as np
import matplotlib.pyplot as plt
import cv2
from microsaccades_functions import *

#---------------------------------------------------------------------------------------IMAGE-CREATION

framerate = 33. #get some fancy frequency calculation
stripe_width = 60
gap = 60
offset = -60
image_size = 120
image_height = 30
degrees = 180

#for normal distributed microsaccades
sigma = 1

file_location = "video/img_input/mo_det_cal/mo_det_cal_m33fr"


repeat = int(framerate)*(stripe_width+gap)
offset_start = offset
film_length = 400 #int(framerate)*(stripe_width+gap)

#get normal 2d distribution
mean = [0, 0]
cov = [[1, 0], [0, 1]]
tl_x, tl_y = sigma*np.random.multivariate_normal(mean, cov, film_length).T

#maybe here a loop to make it last even longer
for f in range(film_length):
    if f%framerate == 0 and f>0: #not yet
        offset += 1
        if f%(repeat/1.) == 0:
            offset=offset_start
    
    canvas = np.zeros((image_height, image_size))
    #next line is for no motion
    #canvas[:,int(image_size/2):]=1.
    
    current_col = offset #there's the problem
    if offset < 0:
        canvas[:, 0:current_col+stripe_width] = 1
        canvas[:, current_col+stripe_width] = (f%int(framerate))/framerate
        current_col += stripe_width + gap
    while current_col < image_size:
        if current_col + stripe_width + gap  <= image_size-1:
            canvas[:, current_col+1:current_col+stripe_width] = 1
            canvas[:, current_col] = 1 - (f%int(framerate))/framerate
            canvas[:, current_col+stripe_width] = (f%int(framerate))/framerate
            current_col += stripe_width + gap
        elif current_col + stripe_width <= image_size-1:
            canvas[:, current_col+1:current_col+stripe_width] = 1
            canvas[:, current_col] = 1 - (f%int(framerate))/framerate
            canvas[:, current_col+stripe_width] = (f%int(framerate))/framerate
            current_col = image_size
        else:
            canvas[:, current_col] = 1 - (f%int(framerate))/framerate
            canvas[:, current_col+1:] = 1
            current_col = image_size
    
    fig = plt.figure()
    fig.set_size_inches(1,0.25)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(canvas, cmap='gray')
    
    '''
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 120)
    plt.close()
    '''
    
    #for turn around/dislocation
    plt.savefig(file_location + "/first"+str(f+1).zfill(3)+".png",  dpi = 120)
    plt.close()
    
    img = cv2.imread(file_location + "/first"+str(f+1).zfill(3)+".png",0)
    rows,cols = img.shape
    
    #-----------------------------------------------------------------------------------------ROTATION
    rot = cv2.getRotationMatrix2D((cols/2-1,rows/2),degrees,1)
    rotFig = cv2.warpAffine(img,rot,(cols,rows))
    
    fig = plt.figure()
    #save to file, comment for additional displacement
    #rotFig = rotFig[150:450,150:450]
    fig.set_size_inches(1, 0.25)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    
    ax.imshow(rotFig,cmap='gray')
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 120)
    plt.close()
    
    #------------------------------------------------------------------------------NORMAL-DISPLACEMENT
    
    '''
    transl = np.float32([[1,0,tl_x[f]],[0,1,tl_y[f]]])
    tlFig = cv2.warpAffine(rotFig,transl,(cols,rows))
    
    #save to file
    tlFig = tlFig[150:450,150:450]
    fig.set_size_inches(1, 1)
    
    plt.imshow(tlFig,cmap='gray')
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 300)
    plt.close()
    '''