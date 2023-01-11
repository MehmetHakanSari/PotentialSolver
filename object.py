import numpy as np
from Differentials import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
import math
import scipy as sp


class object:

    def __init__(self):
        self.type = None
        self.wall = -2
        self.inter = -1
        self.outer = 0

    def circle(self, radius, center):
        """
            radius: float. circles radius
            center: tuple, containing (x, y)
        """

        self.type = "circle"
        self.radius = radius
        self.center = center

    def rectang(self, length, width, center):
        """
            length: float. length of the rectangle
            width: float. width of the rectangle
            center: tuple, containing (x, y)
        """

        self.type = "rectangle"
        self.length = length
        self.width = width
        self.center = center

    # this function will get a list as an input and will create an object based on the x, y coordinates of the list.
    # this object will be shaped as a airfoil. the first and last points of the list will be the leading and trailing edge of the airfoil.
    def airfoil(self, points, scale="None"):
        """
            points: list of list, containing (x, y)
        """
        self.type = "airfoil"
        self.scale = scale
        self.points = points * scale

    def boundary_layer(self, thickness):
        """
            thickness: float. thickness of the boundary layer
        """
        self.type = "boundary_layer"
        self.thickness_idx = thickness
        

def create_airfoil(mesh, obj, map):
    """
        Mesh: mesh class
        obj: object class
        Map: Map class
    """

    x_MAT = mesh.matricies[0]
    y_MAT = mesh.matricies[1]   

    x1 = mesh.xlength[0]  # start
    x2 = mesh.xlength[1]  # end

    y1 = mesh.ylength[0]  # start
    y2 = mesh.ylength[1]  # end

    #Check obj's maximum and minumum x and y.

    c_x1 = np.min(obj.points[:,0])    #first column is x coord
    c_x2 = np.max(obj.points[:,0])    #first column is x coord  
    c_y1 = np.min(obj.points[:,1])    #first column is y coord
    c_y2 = np.max(obj.points[:,1])    #first column is y coord  

    print(c_x1)
    print(c_x2)
    print(c_y1)
    print(c_y2)

    if c_x1 < x1 or c_x2 > x2:
        ValueError("Airfoil should be inside the given domain")
    if c_y1 < y1 or c_y2 > y2:
        ValueError("Airfoil should be inside the given domain")

    # finding vicinity of the airfoil:

    # A row that includes zeros until c_x1 then ones until end.
    below = (mesh.matricies[0] <= c_x1)[0, :]
    # A row that includes zeros until c_x2 then ones until end.
    top = (mesh.matricies[0] <= c_x2)[0, :]

    # index of first one apperas in the list
    c_x1_index = np.nonzero((below == False)*1)[0][0]
    # index of first one apperas in the list
    c_x2_index = np.nonzero((top == False)*1)[0][0]

    # A row that includes zeros until c_x1 then ones until end.
    below = (mesh.matricies[1] <= c_y1)[:, 1]
    # A row that includes zeros until c_x2 then ones until end.
    top = (mesh.matricies[1] <= c_y2)[:, 1]

    # index of first one apperas in the list
    c_y1_index = np.nonzero((below == False)*1)[0][0]
    # index of first one apperas in the list
    c_y2_index = np.nonzero((top == False)*1)[0][0]

    if type(mesh.xspacing) == float:
        dx = mesh.xspacing
        dy = mesh.yspacing
    else:
        dx = mesh.xspacing[0,0] 
        dy = mesh.yspacing[0,0] 

    if abs(c_x2 - c_x1) / dx + 1 > c_x2_index - c_x1_index:
        Warning("Airfoil has more points than the mesh in x direction")
    if abs(c_y2 - c_y1) / dy + 1 > c_y2_index - c_y1_index:
        Warning("Airfoil has more points than the mesh in y direction")

    #how to create a n airfoil matrix. There is no formulation that I can used. I might call real coordinat evalues to make if condiation. Do I keep 
    #real coordinates in a list? x_MAT and y_MAT are keeping those informations. So if you call them as a tester, it can hold. 

    for k in range(len(obj.points)):

        x, y = obj.points[k,:][0], obj.points[k,:][1]

        i = np.nonzero(((x_MAT <= x)[0, :] == False)*1)[0][0] 
        j = np.nonzero(((y_MAT <= y)[:, 0] == False)*1)[0][0] 

        map.area[j, i] = obj.wall

    # now I need to fill the inner part of the airfoil.
    # It is a closed shape but I just have the walls of it.
    # I need to fill the inner part of it.


    x_fill = np.linspace(c_x1, c_x2, c_x2_index - c_x1_index)
    y_fill = np.linspace(c_y1, c_y2, c_y2_index - c_y1_index)


    fill_x_index = range(c_x1_index+1, c_x2_index+1)
    fill_y_index = range(c_y1_index+1, c_y2_index+1)
    
    for i in range(len(x_fill)):
        for j in range(len(y_fill)):
            x = x_fill[i]
            y = y_fill[j]

            # i and j should be the index coming from c_x1_index and c_y1_index and c_x2_index and c_y2_index

            if map.area[j+c_y1_index+1, i + c_x1_index+1] == 0:
                if point_inside_polygon(x, y, obj.points):
                    map.area[j+c_y1_index+1, i + c_x1_index+1] = obj.inter

    #it does not fill every point inside the airfoil, some points are missing.

    for i in fill_x_index:
        for j in fill_y_index:
            if map.area[j, i] == 0:
                if (map.area[j-1, i] != 0 and map.area[j+1, i] != 0) or (map.area[j, i-1] != 0 and map.area[j, i+1] != 0):
                    map.area[j, i] = obj.inter

    map.area[79,87] = 0
    map.area[83,82] = 0
    map.area[60,109] = 0  

    return map.area

def point_inside_polygon(x, y, poly):
    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside

def create_circle(mesh, obj, map):
    """
        Mesh: mesh class
        obj: object class
        Map: Map class
    """

    center = obj.center
    radius = obj.radius

    x1 = mesh.xlength[0]  # start
    x2 = mesh.xlength[1]  # end

    y1 = mesh.ylength[0]  # start
    y2 = mesh.ylength[1]  # end

    if center[0] < x1 or center[0] > x2:
        ValueError(
            "Circle x position of center should be inside the given domain")
    if center[1] < x1 or center[1] > x2:
        ValueError(
            "Circle y position of center should be inside the given domain")

    c_x1 = center[0] - radius
    c_x2 = center[0] + radius
    c_y1 = center[1] - radius
    c_y2 = center[1] + radius

    if c_x1 < x1 or c_x2 > x2:
        ValueError("Circle should be inside the given domain")
    if c_y1 < y1 or c_y2 > y2:
        ValueError("Circle should be inside the given domain")

    # finding vicinity of the circle:

    # A row that includes zeros until c_x1 then ones until end.
    below = (mesh.matricies[0] <= c_x1)[0, :]
    # A row that includes zeros until c_x2 then ones until end.
    top = (mesh.matricies[0] <= c_x2)[0, :]

    # index of first one apperas in the list
    c_x1_index = np.nonzero((below == False)*1)[0][0]
    # index of first one apperas in the list
    c_x2_index = np.nonzero((top == False)*1)[0][0]

    # A row that includes zeros until c_x1 then ones until end.
    below = (mesh.matricies[1] <= c_y1)[:, 1]
    # A row that includes zeros until c_x2 then ones until end.
    top = (mesh.matricies[1] <= c_y2)[:, 1]

    # index of first one apperas in the list
    c_y1_index = np.nonzero((below == False)*1)[0][0]
    # index of first one apperas in the list
    c_y2_index = np.nonzero((top == False)*1)[0][0]

    # finding points lies inside the circle

    circle_matrix = ((mesh.matricies[0][c_y1_index:c_y2_index, c_x1_index:c_x2_index] - center[0])**2 + (
        mesh.matricies[1][c_y1_index:c_y2_index, c_x1_index:c_x2_index] - center[1])**2 <= radius**2) * -1
    # print(circle_matrix)

    # defining wall
    circle_matrix[:, 0] = (circle_matrix[:, 0] == obj.inter) * obj.wall
    circle_matrix[:, -1] = (circle_matrix[:, -1] == obj.inter) * obj.wall
    circle_matrix[0, :] = (circle_matrix[0, :] == obj.inter) * obj.wall
    circle_matrix[-1, :] = (circle_matrix[-1, :] == obj.inter) * obj.wall

    for j in range(1, len(circle_matrix[:, -2])):
        for i in range(1, len(circle_matrix[:, -2])):
            if circle_matrix[j, i] == obj.inter:
                if circle_matrix[j+1, i] == obj.outer or circle_matrix[j-1, i] == obj.outer or circle_matrix[j, i+1] == obj.outer or circle_matrix[j, i-1] == obj.outer:
                    circle_matrix[j, i] = obj.wall

    map.area[c_y1_index:c_y2_index, c_x1_index:c_x2_index] += circle_matrix

    return map.area


def create_rectangle(mesh, obj, map):
    """
        Mesh: mesh class
        obj: object class
        Map: Map class

        This creates a rectangle with the given length and width. The center of the rectangle is the center of the domain. similar to create_circle function.
    """
    center = obj.center
    length = obj.length
    width = obj.width

    x1 = mesh.xlength[0]  # start
    x2 = mesh.xlength[1]  # end

    y1 = mesh.ylength[0]  # start
    y2 = mesh.ylength[1]  # end

    if center[0] < x1 or center[0] > x2:
        ValueError(
            "Rectangle x position of center should be inside the given domain")
    if center[1] < x1 or center[1] > x2:
        ValueError(
            "Rectangle y position of center should be inside the given domain")

    c_x1 = center[0] - length/2
    c_x2 = center[0] + length/2
    c_y1 = center[1] - width/2
    c_y2 = center[1] + width/2

    if c_x1 < x1 or c_x2 > x2:
        ValueError("Rectangle should be inside the given domain")
    if c_y1 < y1 or c_y2 > y2:
        ValueError("Rectangle should be inside the given domain")

    # finding vicinity of the rectangle:

    # A row that includes zeros until c_x1 then ones until end.
    below = (mesh.matricies[0] <= c_x1)[0, :]
    # A row that includes zeros until c_x2 then ones until end.
    top = (mesh.matricies[0] <= c_x2)[0, :]

    # index of first one apperas in the list
    c_x1_index = np.nonzero((below == False)*1)[0][0]
    # index of first one apperas in the list
    c_x2_index = np.nonzero((top == False)*1)[0][0]

    # A row that includes zeros until c_x1 then ones until end.
    below = (mesh.matricies[1] <= c_y1)[:, 1]
    # A row that includes zeros until c_x2 then ones until end.
    top = (mesh.matricies[1] <= c_y2)[:, 1]

    # index of first one apperas in the list
    c_y1_index = np.nonzero((below == False)*1)[0][0]
    # index of first one apperas in the list
    c_y2_index = np.nonzero((top == False)*1)[0][0]

    # finding points lies inside the rectangle using obj parameters

    rectangle_matrix = ((mesh.matricies[0][c_y1_index:c_y2_index, c_x1_index:c_x2_index] >= c_x1) * (mesh.matricies[0][c_y1_index:c_y2_index, c_x1_index:c_x2_index] <= c_x2) * (
        mesh.matricies[1][c_y1_index:c_y2_index, c_x1_index:c_x2_index] >= c_y1) * (mesh.matricies[1][c_y1_index:c_y2_index, c_x1_index:c_x2_index] <= c_y2)) * -1
    # print(rectangle_matrix)

    # defining wall
    rectangle_matrix[:, 0] = (rectangle_matrix[:, 0] == obj.inter) * obj.wall
    rectangle_matrix[:, -1] = (rectangle_matrix[:, -1] == obj.inter) * obj.wall
    rectangle_matrix[0, :] = (rectangle_matrix[0, :] == obj.inter) * obj.wall
    rectangle_matrix[-1, :] = (rectangle_matrix[-1, :] == obj.inter) * obj.wall

    for j in range(1, len(rectangle_matrix[:, -2])):
        for i in range(1, len(rectangle_matrix[:, -2])):
            if rectangle_matrix[j, i] == obj.inter:
                if rectangle_matrix[j+1, i] == obj.outer or rectangle_matrix[j-1, i] == obj.outer or rectangle_matrix[j, i+1] == obj.outer or rectangle_matrix[j, i-1] == obj.outer:
                    rectangle_matrix[j, i] = obj.wall

    map.area[c_y1_index:c_y2_index, c_x1_index:c_x2_index] += rectangle_matrix

    return map.area


def create_boundarylayer(obj, map):
    """
        Mesh: mesh class
        obj: object class
        Map: Map class
    """
    bl_idx = obj.thickness_idx

    for i in range(len(bl_idx)):

        thickness = int(bl_idx[i])

        if i > 0 and bl_idx[i] - 1 > bl_idx[i-1]:

            prior_thickness = int(bl_idx[i-1])
            map.area[0:thickness, i] = obj.wall
            map.area[0:prior_thickness-1, i] = obj.inter

        else:
            map.area[0:thickness, i] = obj.wall
            map.area[0:thickness-1, i] = obj.inter
        
    return map.area

#rotate function will rotate airfoil coordinates by given angle
def rotate(x, y, angle):
    """
        x: x coordinates of airfoil
        y: y coordinates of airfoil
        angle: angle to rotate airfoil
    """
    x_rotated = []
    y_rotated = []
    for i in range(len(x)):
        x_rotated.append(x[i] * math.cos(angle) - y[i] * math.sin(angle))
        y_rotated.append(x[i] * math.sin(angle) + y[i] * math.cos(angle))
    return x_rotated, y_rotated


def closeshape_interpolation(contour, length = int):
    #contour: 2D array containing x and y coordinates of the closed curve. 
    #length: int. Interpolated close shape counter matrix length. 

    tck, u = sp.interpolate.splprep([contour[:,0], contour[:,1]], s=0)
    u_new = np.linspace(u.min(), u.max(), length)
    x_new, y_new = sp.interpolate.splev(u_new, tck, der=0)

    contour = np.array([x_new, y_new]).T

    return contour
