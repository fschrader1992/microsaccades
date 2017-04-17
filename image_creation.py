import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
stripe_width = 15
gap = 20
offset = 2
image_size = 250

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

plt.imshow(canvas)

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
        angle = math.radians(angle)
        temp_point = point[0]-center[0] , point[1]-center[1]
        temp_point = ( -temp_point[0]*math.cos(angle)+temp_point[1]*math.sin(angle) , temp_point[0]*math.sin(angle)-temp_point[1]*math.cos(angle))
        temp_point = [temp_point[0]+center[0] , temp_point[1]+center[1]]
        return temp_point