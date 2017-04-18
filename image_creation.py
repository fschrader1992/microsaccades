import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------

'''
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
    temp_point = ( -temp_point[0]*np.cos(angle)+temp_point[1]*np.sin(angle) , temp_point[0]*np.sin(angle)-temp_point[1]*np.cos(angle))
    temp_point = [temp_point[0]+center[0] , temp_point[1]+center[1]]
    return temp_point
'''
'''
-here should also be a function randomly moving the image and then creating all the frames in a folder
-for moving images in one direction: get input rate and then reset offset + loop



#------------------------------------------------------------------------

stripe_width = 15
gap = 15
offset = 15
image_size = 300

canvas = np.zeros((image_size, image_size))
current_col = 0
while current_col < image_size:
    if current_col + stripe_width + gap <= image_size-1:
        canvas[:, current_col:current_col+stripe_width] = 1
        current_col += stripe_width + gap
    elif current_col + stripe_width <= image_size-1:
        canvas[:, current_col:current_col+stripe_width] = 1
        current_col = image_size
    else:
        canvas[:, current_col:] = 1
        current_col = image_size

#plt.imshow(canvas, cmap='gray')
fig = plt.figure()
fig.set_size_inches(1, 1)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(canvas, cmap='gray')
plt.savefig("test2.png",  dpi = 300)

plt.show()
