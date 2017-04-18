import numpy as np
import matplotlib.pyplot as plt
import cv2

#------------------------------------------------------------------------


def rotate_point(center, point, angle):
    """Rotate point around center by angle
    Args:
        center: 
            Center point (tuple, [x,y])
        point: 
            Point to be rotated (tuple, [x,y])
        angle: 
            Angle in radians to rotate by
    Returns:
        New coordinates for the point
    """
    angle = np.radians(angle)
    temp_point = [point[0]-center[0] , point[1]-center[1]]
    temp_point = [ -temp_point[0]*np.cos(angle)+temp_point[1]*np.sin(angle) , temp_point[0]*np.sin(angle)-temp_point[1]*np.cos(angle)]
    temp_point = [temp_point[0]+center[0] , temp_point[1]+center[1]]
    return temp_point

'''
-here should also be a function randomly moving the image and then creating all the frames in a folder
-for moving images in one direction: get input rate and then reset offset + loop
'''


#---------------------------------------------------------------------------------------IMAGE-CREATION
framerate = 3. #get some fancy frequency calculation
film_length = 30

stripe_width = 15
gap = 15
offset = 0
image_size = 292

for f in range(film_length):
    if f%framerate == 0 and f>0: #not yet
        offset += 1

    canvas = np.zeros((image_size, image_size))
    current_col = offset #there's the problem
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
            canvas[:, current_col:] = 1
            current_col = image_size
    '''        
    for i in range(image_size):
        for j in range(image_size):
            canvas2[i][j] = rotate_point(center, canvas[i][j], angle)
    '''
    #canvas2 = np.roll(canvas, 3)

    #plt.imshow(canvas2, cmap='gray')
    fig = plt.figure()
    fig.set_size_inches(1, 1)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(canvas, cmap='gray')
    plt.savefig("video/img_input/opposite2/frame"+str(f+1).zfill(3)+".png",  dpi = 300)
    
    if framerate%10 == 0:
        plt.show()
'''
fig = plt.figure()
fig.set_size_inches(1, 1)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas2, cmap='gray')
plt.savefig("video/img_input/opposite2/test001.png",  dpi = 300)


#rotation
img = cv2.imread('video/img_input/opposite2/test001.png',0)
rows,cols = img.shape

M = cv2.getRotationMatrix2D((cols/2,rows/2),60,1)
dst = cv2.warpAffine(img,M,(cols,rows))

plt.imshow(dst,cmap='gray'),plt.title('Output')
'''

