import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cv2
import random
#from microsaccades_functions import *
import sys
import os

#---------------------------------------------------------------------------------------IMAGE-CREATION

#radius = 4. #dot size of 4arcmin -> 8px
#rect_size = 110. #size of inner radius 
#vel = 0.06 #velocity is 30 arcmin/s -> 60px/1000ms

#direction = float(angle)*np.pi/8. #dot moving direction in arc (get right velocities !)
def gauss(x,x0,sigma):
    return (np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma)))

image_height = 600
image_width = 600

suff = 'normal_100_bg'
#suff = sys.argv[1]
#cyc = float(sys.argv[2])
#film_length = int(sys.argv[3])
circle_radius=140

#for normal distributed microsaccades
mu = 0.
sigma = 0.447

#file_location = "/home/schrader/Documents/microsaccades/video/img_input/oof/"+str(suff)
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

sigma2 = 3.
rpw = 3

q = np.random.randint(2, size=(int(image_height/rpw)+1,int(image_width/rpw)+1))
qcenter = []
for k in range(400):
    qcenter += [(np.random.triangular(0, 1, 1, 1)[0],np.random.uniform(0,1,1)[0])] #np.random.rand(100,2)
qall = np.random.rand(2000,2)
#print qcenter

canvas = np.zeros((image_height, image_width))
canvas2 = np.zeros((image_height, image_width))
current_col = 0

'''
#NORMAL ILLUSION            
for i in range(image_height):
    for j in range(image_width):
        dist = np.sqrt((float(i)-center[0])*(float(i)-center[0])+(float(j)-center[1])*(float(j)-center[1]))
        if dist >= circle_radius+30:
            canvas[i,j]=float(q[int(i/rpw)][int(j/rpw)])
            canvas2[i,j]=float(q[int(i/rpw)][int(j/rpw)])
        if dist > circle_radius-30:
            if dist < circle_radius+30:
                canvas[i,j]=(circle_radius-30.-dist)/60.*canvas[i,j]-(circle_radius-30.-dist)/60.*float(q[int(i/rpw)][int(j/rpw)])
                #canvas2[i,j]=(circle_radius-30.-dist)/60.*canvas[i,j]-(circle_radius-30.-dist)/60.*float(q[int(i/rpw)][int(j/rpw)])
                canvas2[i,j]=-(circle_radius-30.-dist)/60.*float(q[int(i/rpw)][int(j/rpw)])

for qc in qcenter:   
    #print qc
    for i in range(int(image_height/2)-circle_radius-15,int(image_height/2)+circle_radius+15):
	for j in range(int(image_width/2)-circle_radius-15,int(image_width/2)+circle_radius+15):
            #print i
            pos=[image_height/2.+circle_radius*qc[0]*np.cos(2.*np.pi*qc[1]),image_width/2.+circle_radius*qc[0]*np.sin(2.*np.pi*qc[1])]
            dist = np.sqrt((pos[0]-float(j))*(pos[0]-float(j))+(pos[1]-float(i))*(pos[1]-float(i)))
            val = canvas[i,j]+gauss(dist,0,sigma2)
            if val > 1.:
                val = 1.
            canvas[i,j] = val
            #val = canvas2[i,j]+gauss(dist,0,sigma2)
            #if val > 1.:
            #    val = 1.
            #canvas2[i,j] = val
''' 

#NORMAL ILLUSION WITH LESS BACKGROUND            
for i in range(image_height):
    for j in range(image_width):
        dist = np.sqrt((float(i)-center[0])*(float(i)-center[0])+(float(j)-center[1])*(float(j)-center[1]))
        if dist >= circle_radius+30:
            canvas[i,j]=(0.4+0.2*float(q[int(i/rpw)][int(j/rpw)]))
            canvas2[i,j]=0.5
            #canvas2[i,j]=(0.1+0.8*float(q[int(i/rpw)][int(j/rpw)]))
        if dist > circle_radius-30:
            if dist < circle_radius+30:
                canvas[i,j]=(circle_radius-30.-dist)/60.*canvas[i,j]-(circle_radius-30.-dist)/60.*(0.4+0.2*float(q[int(i/rpw)][int(j/rpw)]))
                #canvas2[i,j]=(circle_radius-30.-dist)/60.*canvas[i,j]-(circle_radius-30.-dist)/60.*float(q[int(i/rpw)][int(j/rpw)])
                canvas2[i,j]=-(circle_radius-30.-dist)/60.*0.5
               
for qc in qcenter:   
    #print qc
    pos=[image_height/2.+circle_radius*qc[0]*np.cos(2.*np.pi*qc[1]),image_width/2.+circle_radius*qc[0]*np.sin(2.*np.pi*qc[1])]
    dist_c = np.sqrt((pos[0]-center[0])*(pos[0]-center[0])+(pos[1]-center[1])*(pos[1]-center[1]))
    for i in range(int(pos[1])-20,int(pos[1])+20):
	for j in range(int(pos[0])-20,int(pos[0])+20):
            #print i
            dist = np.sqrt((pos[0]-float(j))*(pos[0]-float(j))+(pos[1]-float(i))*(pos[1]-float(i)))
            val = canvas[i,j]+gauss(dist,0,sigma2)
            if val > 1.:
                val = 1.
            canvas[i,j] = val
            val = canvas2[i,j]+gauss(dist,0,sigma2)
            if val > 1.:
                val = 1.
            canvas2[i,j] = val

  
'''
#ILLUSION USING JUST DOTS          
for qc in qall:   
    #print qc
    pos=[image_height*qc[0],image_width*qc[1]]
    dist_c = np.sqrt((pos[0]-center[0])*(pos[0]-center[0])+(pos[1]-center[1])*(pos[1]-center[1]))
    for i in range(int(pos[1])-20,int(pos[1])+20):
	for j in range(int(pos[0])-20,int(pos[0])+20):
            if (i >= 0 and i < image_height) and (j >= 0 and j < image_width):
                #print i
                dist = np.sqrt((pos[0]-float(j))*(pos[0]-float(j))+(pos[1]-float(i))*(pos[1]-float(i)))
                val = canvas[i,j]+gauss(dist,0,sigma2)
                if val > 1.:
                    val = 1.
                canvas[i,j] = val
                val = canvas2[i,j]+gauss(dist,0,sigma2)
                if dist_c <= circle_radius:
                    if val > 1.:
                        val = 1.
                    canvas2[i,j] = val
            
''' 
'''
#ILLUSION WITH LOWER CONTRAST
for i in range(image_height):
    for j in range(image_width):
        dist = np.sqrt((float(i)-center[0])*(float(i)-center[0])+(float(j)-center[1])*(float(j)-center[1]))
        if dist >= circle_radius+30:
            canvas[i,j]=0.8*float(q[int(i/rpw)][int(j/rpw)])
        if dist > circle_radius-30:
            if dist < circle_radius+30:
                canvas[i,j]=(circle_radius-30.-dist)/60.*canvas[i,j]-(circle_radius-30.-dist)/60.*0.8*float(q[int(i/rpw)][int(j/rpw)])

for qc in qcenter:   
    #print qc
    for i in range(int(image_height/2)-circle_radius-15,int(image_height/2)+circle_radius+15):
	for j in range(int(image_width/2)-circle_radius-15,int(image_width/2)+circle_radius+15):
            #print i
            pos=[image_height/2.+circle_radius*qc[0]*np.cos(2.*np.pi*qc[1]),image_width/2.+circle_radius*qc[0]*np.sin(2.*np.pi*qc[1])]
            dist = np.sqrt((pos[0]-float(j))*(pos[0]-float(j))+(pos[1]-float(i))*(pos[1]-float(i)))
            val = canvas[i,j]+gauss(dist,0,sigma2)
            if val > 1.:
                val = 1.
            canvas[i,j] = val

'''
'''
#ILLUSION WITH SAME PATTERN BUT DIFFERENT COLOUR

for i in range(image_height):
    for j in range(image_width):
        dist = np.sqrt((float(i)-center[0])*(float(i)-center[0])+(float(j)-center[1])*(float(j)-center[1]))
        canvas[i,j] = 0.5*float(q[int(i/rpw)][int(j/rpw)])
        canvas2[i,j] = 0.5*float(q[int(i/rpw)][int(j/rpw)])
        #inverse
        #canvas2[i,j]=0.
        if dist >= circle_radius+30:
            canvas2[i,j]=0.
            canvas[i,j]=float(q[int(i/rpw)][int(j/rpw)])
            #inverse
            #canvas2[i,j]=float(q[int(i/rpw)][int(j/rpw)])
        if dist > circle_radius-30:
            if dist < circle_radius+30:
                canvas[i,j]=((-circle_radius+30.+dist)/60.+0.5)*float(q[int(i/rpw)][int(j/rpw)])
                canvas2[i,j]=(circle_radius+30.-dist)/60.*canvas2[i,j]
                #inverse
                #canvas2[i,j]=((-circle_radius+30.+dist)/60.)*float(q[int(i/rpw)][int(j/rpw)])
        #if dist < circle_radius-30:
        #    canvas[i,j]=float(q[int(i/rpw)][int(j/rpw)])
        #    canvas2[i,j]=float(q[int(i/rpw)][int(j/rpw)])
        
        #if dist < circle_radius:
        #    canvas[i,j]=canvas[i,j]*2.
        #    canvas2[i,j]=canvas2[i,j]*2.
        #else:
        #    canvas2[i,j]=0.
'''
'''
#ILLUSION WITH SAME PATTERN BUT SHARP AREAS

for i in range(image_height):
    for j in range(image_width):
        dist = np.sqrt((float(i)-center[0])*(float(i)-center[0])+(float(j)-center[1])*(float(j)-center[1]))
        canvas[i,j] = 0.5*float(q[int(i/rpw)][int(j/rpw)])
        canvas2[i,j] = 0.5*float(q[int(i/rpw)][int(j/rpw)])
        #inverse
        #canvas2[i,j]=0.
        if dist >= circle_radius:
            canvas2[i,j]=0.
            canvas[i,j]=float(q[int(i/rpw)][int(j/rpw)])
            #inverse
            #canvas2[i,j]=float(q[int(i/rpw)][int(j/rpw)])
        #if dist > circle_radius-30:
        #    if dist < circle_radius+30:
        #        canvas[i,j]=((-circle_radius+30.+dist)/60.+0.5)*float(q[int(i/rpw)][int(j/rpw)])
        #        canvas2[i,j]=(circle_radius+30.-dist)/60.*canvas2[i,j]
''' 
'''
#ILLUSION WITH SAME PATTERN AND NO CONTRAST BG

for i in range(image_height):
    for j in range(image_width):
        dist = np.sqrt((float(i)-center[0])*(float(i)-center[0])+(float(j)-center[1])*(float(j)-center[1]))
        canvas[i,j] = 0.3+0.4*float(q[int(i/rpw)][int(j/rpw)])
        canvas2[i,j] = 0.5
        #inverse
        #canvas2[i,j]=0.
        if dist >= circle_radius+30:
            canvas[i,j]=0.5
            #inverse
            #canvas2[i,j]=float(q[int(i/rpw)][int(j/rpw)])
        if dist > circle_radius-30:
            if dist < circle_radius+30:
                canvas[i,j]=((circle_radius+30.-dist)/60.)*(0.3+0.4*float(q[int(i/rpw)][int(j/rpw)]))+0.5*((-circle_radius+30.+dist)/60.)
                #inverse
                #canvas2[i,j]=((-circle_radius+30.+dist)/60.)*float(q[int(i/rpw)][int(j/rpw)])
        #if dist < circle_radius-30:
        #    canvas[i,j]=float(q[int(i/rpw)][int(j/rpw)])
        #    canvas2[i,j]=float(q[int(i/rpw)][int(j/rpw)])
        
        #if dist < circle_radius:
        #    canvas[i,j]=canvas[i,j]*2.
        #    canvas2[i,j]=canvas2[i,j]*2.
        #else:
        #    canvas2[i,j]=0.
'''
'''
#cm = plt.get_cmap('RdBu', 3)
cm = mpl.colors.ListedColormap(['black', 'white', 'chocolate'])

fig = plt.figure()
fig.set_size_inches(3,3)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas, cmap = cm )
plt.savefig("/home/schrader/Documents/microsaccades/img/murakami_illusion/"+str(suff)+".png",  dpi = 200)
plt.close()

fig = plt.figure()
fig.set_size_inches(3,3)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas2, cmap = cm )
plt.savefig("/home/schrader/Documents/microsaccades/img/murakami_illusion/"+str(suff)+"_off.png",  dpi = 200)

'''
#for non-coloured
fig = plt.figure()
fig.set_size_inches(3,3)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas, cmap='gray',norm=mpl.colors.Normalize(vmin=0.,vmax=1.))
plt.savefig("/home/schrader/Documents/microsaccades/img/murakami_illusion/"+str(suff)+".png",  dpi = 200)
plt.close()

fig = plt.figure()
fig.set_size_inches(3,3)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas2, cmap='gray',norm=mpl.colors.Normalize(vmin=0.,vmax=1.))
plt.savefig("/home/schrader/Documents/microsaccades/img/murakami_illusion/"+str(suff)+"_off.png",  dpi = 200)

#plt.show()
#plt.close()

'''
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