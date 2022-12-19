import numpy as np
from Differentials import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm


class Map():
    
    def __init__(self, area):
        """
            Area is one matrix that you can show object of different kinds.
        """
        self.area = area

    def __call__(self):
        return self.area

    def create_object(self, physical_domain ,obj):
        """
           self contains the map and the area is updated accordingly. 
           
           physical domain contains the mesh class.
           solution matricies of X and Y can be extracted by using 
           physical_domain.matricies and nodes can be extreacted by physical_domain.nodes.
           
           obj_coordinate: objects coordinates that should be fitted inside the domain. How would object information save into X and Y matricies?
           Map matrix or map object is perfect solution for this. Similar what I have done for fish-fisherman. Map can be translated to to 
           solver class. According to map, it can specify if there is object or not. 

           object coordiante can be list of points. It should enclose itself. Or I can write predefined objects. 
           The think is as the nodes are given already I should automize the object creating. 

           enclosing will take time to implement. 

           obj_type: circle
           obj_coor: radius, center

           obj_type: rectengale
           obj_coor: center, x_length, y_length etc.
        """
        if obj.type == "circle":
            radius = obj.radius
            center = obj.center

        x1 = self.xlength[0] 
        x2 = self.xlength[1]

        y1 = self.ylength[0] 
        y2 = self.ylength[1]

        if center[0] < x1 or center[0] > x2:
            ValueError("Circle x position of center should be inside the given domain")
        if center[1] < x1 or center[1] > x2:
            ValueError("Circle y position of center should be inside the given domain")

        
        """
                c_y2
                ****
               ******
              *********
             ************
       c_x1 ************** c_x2
             ************
              **********
               *******
                *****
                c_y1

                for any x,y that (x - x0)^2 - (y-y0)^2 <= r^2 lay inside the circle. you can start this in the vicinity of the circle.
        """

        c_x1 = center[0] - radius
        c_x2 = center[0] + radius
        c_y1 = center[1] - radius
        c_y2 = center[1] + radius

        if c_x1 < x1 or c_x2 > x2:
            ValueError("Circle should be inside the given domain")
        if c_y1 < y1 or c_y2 > y2:
            ValueError("Circle should be inside the given domain")

        #finding vicinity of the circle:

        below = (self.matricies[0] <= c_x1)[0,:]      #A row that includes zeros until c_x1 then ones until end. 
        top = (self.matricies[0] <= c_x2)[0,:]        #A row that includes zeros until c_x2 then ones until end. 

        c_x1_index = np.nonzero((below == False)*1)[0][0]                #index of first one apperas in the list
        c_x2_index = np.nonzero((top == False)*1)[0][0]                    #index of first one apperas in the list

        below = (self.matricies[1] <= c_y1)[:,1]      #A row that includes zeros until c_x1 then ones until end. 
        top = (self.matricies[1] <= c_y2)[:,1]        #A row that includes zeros until c_x2 then ones until end. 

        c_y1_index = np.nonzero((below == False)*1)[0][0]               #index of first one apperas in the list
        c_y2_index = np.nonzero((top == False)*1)[0][0]                  #index of first one apperas in the list

        #finding points lies inside the circle

        circle_matrix = ((self.matricies[0][c_y1_index:c_y2_index, c_x1_index:c_x2_index] - center[0])**2 + (self.matricies[1][c_y1_index:c_y2_index, c_x1_index:c_x2_index] - center[1])**2 <= radius**2) * -1

        self.map()[c_y1_index:c_y2_index, c_x1_index:c_x2_index] = self.map()[c_y1_index:c_y2_index, c_x1_index:c_x2_index] + circle_matrix

    def show(self):
        
        fig, ax = plt.subplots()

        # z = np.ones((self.nodes[1], self.nodes[0])) * 0.5

        image = ax.pcolormesh(self.area, vmin=-1, vmax=1, edgecolors="black",linewidth=0.1)
        plt.show()
