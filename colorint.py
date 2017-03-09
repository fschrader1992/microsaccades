import cv2
import numpy as np

out='<name>'.jpg'
img=cv2.imread(out,0)

x1=
x2=
y1=
y2=
redvalue=
greenvalue=
bluevalue=
whitevalue=

#get planes -> check if format is bgr
red=img[:,:,2]
green=img[:,:,1]
blue=img[:,:,0] 

cell=red[y1:y2,x1:x2]
cell=np.sum(abs(cell))
intensity_red=cell/255*redvalue
print('red intensity'+'\t'+str(intensity_red)+'\n')

cell=green[y1:y2,x1:x2]
cell=np.sum(abs(cell))
intensity_green=cell/255*greenvalue
print('green intensity'+'\t'+str(intensity_gree)+'\n')

cell=blue[y1:y2,x1:x2]
cell=np.sum(abs(cell))
intensity_blue=cell/255*bluevalue
print('blue intensity'+'\t'+str(intensity_blue)+'\n')

#fuer weiss auf grayscale schalten, evtl kleineren bildbereich wählen
gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
cell=gray[y1:y2,x1:x2]
cell=cv2.inRange(cell,200,255) #andere intensitäten für weiss setzen
cell=np.sum(abs(cell))
intensity_blue=cell/255*whitevalue
print('white intensity'+'\t'+str(intensity_white)+'\n')

intensity_all=intensity_red+intensity_green+intensity_blue+intensity_white
print('intensity all'+'\t'+str(intensity_all)+'\n')