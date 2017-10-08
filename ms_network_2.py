#THIS IS THE PART OF THE NEURAL NETWORK. IT TAKES THE VALUES TO CALCULATE THE POISSONRATES OF THE INPUT NEURONS AND THEN CALCULATES THE NETWORK OUTPUT 

import pylab as pyl
import numpy as np
import cv2
#import nest
from microsaccades_functions import *

#-----------------------------------------------------------------------------POISSON+SYNAPTIC-WEIGHTS
#here needs to be a part that transfers potentials into poisson rates
tf_data = open('data/temp_filter_subtr.data','r+')
data = np.load(tf_data)   
print data
print poissonRate(87)

#calculate synaptic weights


#-----------------------------------------------------------------------------------------NETWORK-PART
#---------------------------------------------------------------------------INITIALIZE-POISSON-NEURONS
