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

image_height = 120
image_width = 480

suff = sys.argv[1]
cyc1 = float(sys.argv[2])
cyc2 = float(sys.argv[3])
film_length = int(sys.argv[4])

#for normal distributed microsaccades
mu = 0.
sigma = 0.447

file_location = "/home/schrader/Documents/microsaccades/video/img_input/murakami/"+str(suff)
#save the different conditions

#get normal 2d distribution
#mean = [0, 0]
#cov = [[1, 0], [0, 1]]
#tl_x, tl_y = sigma*np.random.multivariate_normal(mean, cov, film_length).T
disp = np.random.normal(mu,sigma,500)
#disp = open('data/oof/displacement.data','w+')
#use displacement here?

sigma2 = 16
rpw = 32
jmin = 120
jmax = 360

pos = (image_height/2.+0.5,image_width/2+0.5)
center = (image_height/2.-0.5,image_width/2.-0.5)

q = np.random.randint(2, size=(int(image_height/rpw)+1,int(image_width/rpw)+1))
qcenter = np.random.rand(8,2)

canvas = np.zeros((image_height, image_width))
canvas2 = np.zeros((image_height, image_width))
current_col = 0

#idea: just in one dimension, on the outer region random square pattern, on the center region randomly dispersed dots, then linear gradient between the two regions

for qc in qcenter:   
    print qc
    for i in range(image_height):
        for j in range(jmin-30,jmax+30):
            dist = np.sqrt((float(i)-image_height*qc[0])*(float(i)-image_height*qc[0])+(float(j)-(jmax-jmin)*qc[1]-jmin)*(float(j)-(jmax-jmin)*qc[1]-jmin))
            val = canvas[i,j]+gauss(dist,0,sigma2)
            if val > 1.:
                val = 1.
            canvas[i,j] = val
            val = canvas2[i,j]+gauss(dist,0,sigma2)
            if val > 1.:
                val = 1.
            canvas2[i,j] = val
            
#here then blurring, not so huge dots (if partners in x/y direction?)
#canvas[i,j]=0.5*np.sin((cyc/120.*np.pi*float(i)-center[0]))*np.sin((cyc/120.*np.pi*float(j)-center[1]-60))+0.5

#canvas[i,j]=float(qcenter[int(i/3)][int(j/3)])
#canvas2[i,j]=float(qcenter[int(i/3)][int(j/3)])

#canvas2 = cv2.GaussianBlur(canvas2,(3,3),0.5)
#overwrite the rest with random dots of size 3arcmin = 6px
#also possible to overwrite this with circle
#90,210 are centers of transition
     
#canvas = cv2.GaussianBlur(canvas,(7,7),20)
for i in range(image_height):
    for j in range(image_width):
        if j <= jmin-30 or j > jmax+30:
            canvas[i,j]=float(q[int(i/rpw)][int(j/rpw)])
            dist = np.sqrt((float(i)-center[0])*(float(i)-center[0])+(float(j)-center[1])*(float(j)-center[1]))
           
        if (j > jmin-30 and j < jmin+30):
            canvas[i,j]=-(jmin-30.-float(j))/60.*canvas[i,j]+(jmin+30.-float(j))/60.*float(q[int(i/rpw)][int(j/rpw)])
        if (j > jmax-30 and j <= jmax+30):
            canvas[i,j]=(jmax+30.-float(j))/60.*canvas[i,j]-(jmax-30.-float(j))/60.*float(q[int(i/rpw)][int(j/rpw)])
            #if dist <= circle_width:
            #    canvas[i,j] = 1
            #elif 3.5 < dist and dist  < 4.5:
            #    canvas[i,j] = (1. -4.5 + dist)*canvas[i,j] + 4.5 - dist
            #rectangle-(exp2)--------------------------
            #if rect==1:
            #    if abs(float(i)-center[0]) > rect_size or abs(float(j)-center[1])> rect_size:
            #        canvas[i,j] = 1
            #------------------------------------------
   

print canvas
   
fig = plt.figure()
fig.set_size_inches(4,1)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas, cmap='gray')
plt.savefig(file_location + "/first.png",  dpi = 120)
#plt.show()
plt.close()

fig = plt.figure()
fig.set_size_inches(4,1)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas2, cmap='gray')
plt.savefig(file_location + "/first2.png",  dpi = 120)
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
    fig.set_size_inches(4,1)
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
    frac, whole = math.modf(disp[f])
    img_s = (1.-frac)*img + frac*img_t
    print '1: '+str(np.amax(img_s))
    #print img_s
    for i in range(len(img_s)):
        for j in range(len(img_s[0])):
            if img_s[i][j] > 255.:
                img_s[i][j] = 255.
            elif img_s[i][j] < 0.:
                img_s[i][j] = 0.
    
    print '2: '+str(np.amax(img_s))
    print img_s
    transl = np.float32([[1,0,int(disp[f])],[0,1,0]])
    tlFig = cv2.warpAffine(img_s,transl,(cols,rows))
   
   
    #save to file, comment for additional rotation
    #fig = plt.figure()
    #tlFig = tlFig[:,30:270]
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