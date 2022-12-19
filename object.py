import numpy as np
from Differentials import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm


class object:

    def __init__(self):
        self.type = None

    def circle(self, radius, center):
        """
            radius: float. circles radius
            center: tuple, containing (x, y)
        """

        self.type = "circle"
        self.radius = radius
        self.center = center


def create_circle(mesh, obj, map):

    center = obj.center 
    radius = obj.radius

    x1 = mesh.xlength[0] #start
    x2 = mesh.xlength[1] #end

    y1 = mesh.ylength[0] #start
    y2 = mesh.ylength[1] #end

    if center[0] < x1 or center[0] > x2:
        ValueError("Circle x position of center should be inside the given domain")
    if center[1] < x1 or center[1] > x2:
        ValueError("Circle y position of center should be inside the given domain")

    c_x1 = center[0] - radius
    c_x2 = center[0] + radius
    c_y1 = center[1] - radius
    c_y2 = center[1] + radius

    if c_x1 < x1 or c_x2 > x2:
        ValueError("Circle should be inside the given domain")
    if c_y1 < y1 or c_y2 > y2:
        ValueError("Circle should be inside the given domain")

    #finding vicinity of the circle:

    below = (mesh.matricies[0] <= c_x1)[0,:]      #A row that includes zeros until c_x1 then ones until end. 
    top = (mesh.matricies[0] <= c_x2)[0,:]        #A row that includes zeros until c_x2 then ones until end. 

    c_x1_index = np.nonzero((below == False)*1)[0][0]                #index of first one apperas in the list
    c_x2_index = np.nonzero((top == False)*1)[0][0]                    #index of first one apperas in the list

    below = (mesh.matricies[1] <= c_y1)[:,1]      #A row that includes zeros until c_x1 then ones until end. 
    top = (mesh.matricies[1] <= c_y2)[:,1]        #A row that includes zeros until c_x2 then ones until end. 

    c_y1_index = np.nonzero((below == False)*1)[0][0]               #index of first one apperas in the list
    c_y2_index = np.nonzero((top == False)*1)[0][0]                  #index of first one apperas in the list

    #finding points lies inside the circle

    circle_matrix = ((mesh.matricies[0][c_y1_index:c_y2_index, c_x1_index:c_x2_index] - center[0])**2 + (mesh.matricies[1][c_y1_index:c_y2_index, c_x1_index:c_x2_index] - center[1])**2 <= radius**2) * -1
    map.area[c_y1_index:c_y2_index, c_x1_index:c_x2_index] += circle_matrix
    # map.area()[c_y1_index:c_y2_index, c_x1_index:c_x2_index]
    
    return map.area