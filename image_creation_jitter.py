import numpy as np
import matplotlib.pyplot as plt
import math
import cv2
#from microsaccades_functions import *
import sys
import os

def gauss(x,x0,sigma):
    return (np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma)))
            
#---------------------------------------------------------------------------------------IMAGE-CREATION

#radius = 4. #dot size of 4arcmin -> 8px
#rect_size = 110. #size of inner radius
#vel = 0.06 #velocity is 30 arcmin/s -> 60px/1000ms

#direction = float(angle)*np.pi/8. #dot moving direction in arc (get right velocities !)

image_height = 128
image_width = 256
vel = 0.06 #velocity is 30 arcmin/s -> 60px/1000ms

suff = sys.argv[1]
cyc1 = float(sys.argv[2])
cyc2 = float(sys.argv[3])
vel = 0.01*float(sys.argv[4])
film_length = int(sys.argv[5])
rpw = int(sys.argv[6])

#for normal distributed microsaccades
mu = 0.
sigma = 0.447

file_location = "/home/schrader/Documents/microsaccades/video/img_input/jitter/"+str(suff)
#save the different conditions

#get normal 2d distribution
#mean = [0, 0]
#cov = [[1, 0], [0, 1]]
#tl_x, tl_y = sigma*np.random.multivariate_normal(mean, cov, film_length).T
disp = np.random.normal(mu,sigma,500)
#disp = open('data/oof/displacement.data','w+')
#use displacement here?

sigma2 = 16
#rpw = 64 #32
jmin = 120
jmax = 360

pos = (image_height/2.+0.5,image_width/2+0.5)
center = (image_height/2.-0.5,image_width/2.-0.5)


#q = np.random.randint(2, size=(int(image_height/rpw)+1,int(image_width/rpw)+1))
q = np.zeros((int(image_height/rpw)+1,int(image_width/rpw)+1))

for i in range(int(image_height/rpw)+1):
    for j in range(int(image_width/rpw)+1):
        if (i%2==0 and j%2==1) or (j%2==0 and i%2==1):
            q[i][j] = 1.

canvas = np.zeros((image_height, image_width))
canvas2 = np.zeros((image_height, image_width))
current_col = 0

#idea: just in one dimension, on the outer region random square pattern, on the center region randomly dispersed dots, then linear gradient between the two regions

for i in range(image_height):
    for j in range(image_width):
        canvas[i,j]=float(q[int(i/rpw)][int(j/rpw)])
           
        
#print canvas
   
fig = plt.figure()
fig.set_size_inches(2,1)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas, cmap='gray')
plt.savefig(file_location + "/first.png",  dpi = 128)
#plt.show()
plt.close()

fig = plt.figure()
fig.set_size_inches(2,1)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas2, cmap='gray')
plt.savefig(file_location + "/first2.png",  dpi = 128)
#plt.show()
plt.close()


img1 = cv2.imread(file_location + "/first.png",0)
img2 = cv2.imread(file_location + "/first2.png",0)
rows,cols = img1.shape
#d_file = open('data/phase_displacement.data','r+')
#disp = np.load(d_file)  
#d_file.close()

for f in range(film_length):
    if f%(cyc1) >= cyc1-cyc2:
        img=img2
    else:
        img=img1
    #------------------------------------------------------------------------------NORMAL-DISPLACEMENT
    fig = plt.figure()
    fig.set_size_inches(1,1)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    #ax.imshow(canvas, cmap='gray')
    
    #in order to make the tranlation smooth
    img_s=[] #[0 for j in range(len(img[0]))] for i in range(len(img))]
    i=0
    #print img
    img_t=[]
    img_t=img
    img_t = np.insert(img_t, len(img_t[0]), img_t[:,0], axis = 1)
    img_t = np.delete(img_t, 0, 1)
    frac, whole = math.modf(5*vel*float(f))
    img_s = (1.-frac)*img_t + frac*img
    #print '1: '+str(np.amax(img_s))
    #print img_s
    for i in range(len(img_s)):
        for j in range(len(img_s[0])):
            if img_s[i][j] > 255.:
                img_s[i][j] = 255.
            elif img_s[i][j] < 0.:
                img_s[i][j] = 0.
    #print '2: '+str(np.amax(img_s))
    #print img_s
    
    #uncomment for rem
    #transl = np.float32([[1,0,int(disp[f])],[0,1,0]])
    transl = np.float32([[1,0,int(5*vel*float(f))],[0,1,0]])
    tlFig = cv2.warpAffine(img_s,transl,(cols,rows))
   
    #save to file, comment for additional rotation
    #fig = plt.figure()
    tlFig = tlFig[:,128:]
    
    tlFig = cv2.GaussianBlur(tlFig,(3,3),0.5)
    
    #fig.set_size_inches(2,0.25)
    ax.imshow(tlFig,cmap='gray')
    plt.savefig(file_location + "/second"+str(f+1).zfill(3)+".png",  dpi = 128)
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