import numpy as np
from Differentials import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from object import create_circle, create_rectangle, create_airfoil



class Map:
    """
        Creates the map of the given physical domain with its meshed surface. Updates itself by creating object inside the given matrix. 
    """
    
    def __init__(self, physical_domain):
        """
            Area is one matrix that you can show object of different kinds.

            pysical_domain is mesh.
        """
        self.mesh = physical_domain
        self.area = np.zeros((self.mesh.nodes[1], self.mesh.nodes[0]))
        self.interior = -1
        self.wall = -2
        self.outer = 0

    def __call__(self):
        return self.area

    def geometry_guide(self):

        
        """    
                CIRCLE

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
        
        pass

    def create_object(self, obj):
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

        mesh = self.mesh

        if obj.type == "circle":
            self.area = create_circle(mesh, obj ,self)
        if obj.type == "rectangle":
            self.area = create_rectangle(mesh, obj ,self)
        if obj.type == "airfoil":
            self.area = create_airfoil(mesh, obj ,self)
        
                

        # self.area()[c_y1_index:c_y2_index, c_x1_index:c_x2_index] = self.area()[c_y1_index:c_y2_index, c_x1_index:c_x2_index] + circle_matrix

    def show(self):
        
        fig, ax = plt.subplots()

        # z = np.ones((self.nodes[1], self.nodes[0])) * 0.5
        minvalue = np.min(self.area)
        maxvalue = np.max(self.area)
        fig.set_size_inches(15, 15)
        image = ax.pcolormesh(self.area, vmin=minvalue, vmax=maxvalue, edgecolors="none")
        plt.show()
