#!/bin/python
#THIS PART READS THE PIXELS OF A MOVIE AND APPLIES THE TEMPORAL FILTERS TO THE PROPOSED MODEL FOR MICROSACCADES. THE OUTPUT ARE THE POTENTIAL VALUES FOR THE POISSON RATES.
'''
loading the video into an array
-> first get all the pixels in a frame (grayscale)
-> after that calculate the values of each tempral flter at each time and store them in another array which will be the basis for changing poisson rates
'''

import sys
import os
import glob
import pylab as pyl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cv2
import datetime
import itertools 
#import multiprocessing
from microsaccades_functions import *

#---------------------------------------------------------------------------------------------------LOAD-FRAMES
now = datetime.datetime.now()
'''
#for parallel computing
try:
    cpus = multiprocessing.cpu_count()
except NotImplementedError:
    cpus = 2   # arbitrary default
	
pool = multiprocessing.Pool(processes=cpus)
'''
sim_title = sys.argv[1]
handle_name = sys.argv[2]
fn = sys.argv[3]


#alternative method for videos
frames = []
frame_number = int(fn)
os.chdir("video/img_input/" + sim_title)
for file in glob.glob("second*.png"):
    #print(file)
    frames+=[file]
frames.sort()
#print frames

f=cv2.imread(frames[0])
height, width = f.shape[:2]
dt = 1.
print height, width

#assign the all time all pixel array/list
pixels4d = [[[] for j in range(width)] for i in range(height)]


for file in frames:

    frame = cv2.imread(file)
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    #store 2D array in pixels4d 
    for i in range(height):
        for j in range(width):
            pixels4d[i][j]+=[float(gray[i,j])]
    
cv2.destroyAllWindows()

#os.chdir("../../..")


parasol_grid= [[(0,0) for j in range(3)] for i in range(3)]
p_pos = np.asarray(parasol_grid)
p_pos_data = open('data/phases/p_pos.data','w+')
np.save(p_pos_data, p_pos)
p_pos_data.close()
