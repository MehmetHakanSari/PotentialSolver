import numpy as np

class Domain:
    """
        Domain is the phsical domain that problem will be solved. (Formal)
        Unformal:
        In my mind I think to do Block type mesh. However For an airfoil which kind of mesh is used I am unaware of that. 
        C mesh for instance is one way to mesh it right? But could we seperate object boundaries while doing C-mesh. 
        First thing to do is 2D simple 1 block mesh. 

        
        Cartasiean:

        Blocks Physical boundaries:

        For 2D: 

        each element is tuple of 3 elements, specifiying lengths of the block. 

        [(x1, y1, z1),        Left Up Corner
         (x2, y1, z1),        Right Up Corner
         (x2, y2, z1),        Right Down Corner
         (x1, y2, z1)]        Left Down Corner

        Mesh Properties:

        [N_x, N_y, N_z] = Total node numbers in x, y, z direction. 

        
    """

    def __init__(self, pyhsical_domain, nodes):
        #write in PEP8 format.
        #list consist of 4 tuple with 3 elements (2D)
        #list consist of 3 float elements        (2D or 3D)

        self.pyhsical_domain = pyhsical_domain  
        self.nodes = nodes             
        self.matricies = None         

    def uniform_grid_2D(self):
        
        if self.pyhsical_domain[0][0] != self.pyhsical_domain[1][0]:
            #if the given coordinates are in same line for y coordinate pass other coordinate
            x_list = np.linspace(self.pyhsical_domain[0][0], self.pyhsical_domain[1][0], self.nodes[0])
        elif self.pyhsical_domain[1][0] != self.pyhsical_domain[2][0]:
            x_list = np.linspace(self.pyhsical_domain[1][0], self.pyhsical_domain[2][0], self.nodes[0])

        if self.pyhsical_domain[0][1] != self.pyhsical_domain[1][1]:
            #if the given coordinates are in same line for x coordinate pass other coordinate
            y_list = np.linspace(self.pyhsical_domain[0][1], self.pyhsical_domain[1][1], self.nodes[1])
        elif self.pyhsical_domain[1][1] != self.pyhsical_domain[2][1]:
            y_list = np.linspace(self.pyhsical_domain[1][1], self.pyhsical_domain[2][1], self.nodes[1])

        if x_list[-1] < x_list[0]:      #Small value to large value to be consistent for physiscs
            x_list = np.flip(x_list)
        if y_list[-1] < y_list[0]:
            y_list = np.flip(y_list)
            
        x_MAT = np.zeros((self.nodes[1], self.nodes[0]))
        y_MAT = np.zeros((self.nodes[1], self.nodes[0]))

        #for uniform grid one can use np.meshgrid()
        for i in range(len(y_list)):
            x_MAT[i,:] = x_list
        for i in range(len(x_list)):
            y_MAT[:,i] = y_list

        self.matricies = [x_MAT, y_MAT]

    def nonuniform_grid_2D(self, g_x = 1, g_y = 1):
        #g_x is the gradient for delta_x. Float. Default value 1. 
        #g_y is the gradient for delta_Y. Float. Default value 1. 

        #grid size changes for each step. 

        x_list = np.zeros(self.nodes[0])
        y_list = np.zeros(self.nodes[1])

        if self.pyhsical_domain[0][0] != self.pyhsical_domain[1][0]:
            ##calculate first step sizes 
            L_x = abs((self.pyhsical_domain[1][0] - self.pyhsical_domain[0][0])) #length in x direction
            if g_x == 1:
                dx = L_x / (self.nodes[0] - 1)
            else:
                dx = L_x * (1 - g_x) / (1 - g_x**self.nodes[0])  #first step size. The remaining floats should be considered. Rounding error should be expected.
            x_list[0] = self.pyhsical_domain[0][0] 
            for i in range(1,self.nodes[0]):
                x_list[i] = x_list[i-1] + dx * g_x**(i-1)  
            
        elif self.pyhsical_domain[1][0] != self.pyhsical_domain[2][0]:
            L_x = abs((self.pyhsical_domain[2][0] - self.pyhsical_domain[1][0])) #length in x direction
            if g_x == 1:
                dx = L_x / (self.nodes[0] - 1)
            else:
                dx = L_x * (1 - g_x) / (1 - g_x**self.nodes[0])  #first step size. The remaining floats should be considered. Rounding error should be expected.
            x_list[0] = self.pyhsical_domain[1][0] 
            for i in range(1,self.nodes[0]):
                x_list[i] = x_list[i-1] + dx * g_x**(i-1)  
            

        if self.pyhsical_domain[0][1] != self.pyhsical_domain[1][1]:
            L_y = abs((self.pyhsical_domain[1][1] - self.pyhsical_domain[0][1])) #length in x direction
            if g_y == 1:
                dy = L_y / (self.nodes[1] - 1)
            else:
                dy = L_x * (1 - g_x) / (1 - g_y**self.nodes[1])  #first step size. The remaining floats should be considered. Rounding error should be expected.
            y_list[0] = self.pyhsical_domain[0][1] 
            for i in range(1,self.nodes[1]):
                y_list[i] = y_list[i-1] + dy * g_y**(i-1)  
            
        elif self.pyhsical_domain[1][1] != self.pyhsical_domain[2][1]:
            L_y = abs((self.pyhsical_domain[2][1] - self.pyhsical_domain[1][1])) #length in x direction
            if g_y == 1:
                dy = L_y / (self.nodes[1] - 1)
            else:
                dy = L_x * (1 - g_y) / (1 - g_y**self.nodes[1])  #first step size. The remaining floats should be considered. Rounding error should be expected.
            y_list[0] = self.pyhsical_domain[1][1] 
            for i in range(1,self.nodes[1]):
                y_list[i] = y_list[i-1] + dy * g_y**(i-1)  
            
        # x_list = np.array(x_list)
        # y_list = np.array(y_list)

        if x_list[-1] < x_list[0]:      #Small value to large value to be consistent for physiscs
            x_list = np.flip(x_list)
        if y_list[-1] < y_list[0]:
            y_list = np.flip(y_list)
            
        x_MAT = np.zeros((self.nodes[1], self.nodes[0]))
        y_MAT = np.zeros((self.nodes[1], self.nodes[0]))

        #for uniform grid one can use np.meshgrid()
        for i in range(len(y_list)):
            x_MAT[i,:] = x_list
        for i in range(len(x_list)):
            y_MAT[:,i] = y_list

        self.matricies = [x_MAT, y_MAT]
        
        



