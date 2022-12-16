import numpy as np
from Differentials import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from object import object



class Mesh:
    """
        Domain is the phsical domain that problem will be solved. (Formal)
        Unformal:
        In my mind I think to do Block type mesh. However For an airfoil which kind of mesh is used I am unaware of that. 
        C mesh for instance is one way to mesh it right? But could we seperate object boundaries while doing C-mesh. 
        First thing to do is 2D simple 1 block mesh. 

        
        Cartasiean:

        Blocks Physical boundaries:

        For 2D: 

        each element is tuple of 3 elements, specifiying lengths of the block. The order of coordinates should be clockwise. 

        [(x1, y1, z1),        
         (x2, y1, z1),        
         (x2, y2, z1),        
         (x1, y2, z1)]        

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

    def uniform_block_mesh_2D(self):

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

        x_spacing = float(x_list[1] - x_list[0]) 
        y_spacing = float(y_list[1] - y_list[0])   #signs convention might not hold
            
        x_MAT = np.zeros((self.nodes[1], self.nodes[0]))
        y_MAT = np.zeros((self.nodes[1], self.nodes[0]))

        #for uniform grid one can use np.meshgrid()
        for i in range(len(y_list)):
            x_MAT[i,:] = x_list
        for i in range(len(x_list)):
            y_MAT[:,i] = y_list

        self.matricies = [x_MAT, y_MAT]

        return x_spacing, y_spacing

    def nonuniform_block_mesh_2D(self, g_x = 1, g_y = 1):
        #g_x is the gradient for delta_x. Float. Default value 1. 
        #g_y is the gradient for delta_Y. Float. Default value 1. 

        #grid size changes for each step. 

        x_list = np.zeros(self.nodes[0], dtype="float")
        y_list = np.zeros(self.nodes[1], dtype="float")
        x_spacing = np.zeros((1, self.nodes[0] - 1), dtype="float")
        y_spacing = np.zeros((1, self.nodes[1] - 1), dtype="float")

        if self.pyhsical_domain[0][0] != self.pyhsical_domain[1][0]:
            self.xlength = sorted((self.pyhsical_domain[0][0], self.pyhsical_domain[1][0]))
            ##calculate first step sizes 
            L_x = abs((self.pyhsical_domain[1][0] - self.pyhsical_domain[0][0])) #length in x direction
            if g_x == 1:
                dx = L_x / (self.nodes[0] - 1)
            else:
                if g_x > 0:
                    dx = L_x * (1 - g_x) / (1 - g_x**(self.nodes[0] - 1)) 
                    if self.pyhsical_domain[0][0] < self.pyhsical_domain[1][0]:
                        x_list[0] = self.pyhsical_domain[0][0] 
                    elif self.pyhsical_domain[0][0] > self.pyhsical_domain[1][0]:
                        x_list[0] = self.pyhsical_domain[1][0]
                else:
                    g_x = abs(g_x)
                    dx = L_x * (1 - g_x) / (1 - g_x**(self.nodes[0] - 1))
                    dx = -dx
                    if self.pyhsical_domain[0][0] < self.pyhsical_domain[1][0]:
                        x_list[0] = self.pyhsical_domain[0][0] 
                    elif self.pyhsical_domain[0][0] > self.pyhsical_domain[1][0]:
                        x_list[0] = self.pyhsical_domain[1][0]
            for i in range(1,self.nodes[0]):
                x_list[i] = x_list[i-1] + dx * g_x**(i-1)  
                x_spacing[:,i-1] = dx * g_x**(i-1)  
            
        elif self.pyhsical_domain[1][0] != self.pyhsical_domain[2][0]:
            self.xlength = sorted((self.pyhsical_domain[1][0], self.pyhsical_domain[2][0]))
            L_x = abs((self.pyhsical_domain[2][0] - self.pyhsical_domain[1][0])) #length in x direction
            if g_x == 1:
                dx = L_x / (self.nodes[0] - 1)
            else:
                if g_x > 0:
                    dx = L_x * (1 - g_x) / (1 - g_x**(self.nodes[0] - 1)) 
                    if self.pyhsical_domain[1][0] < self.pyhsical_domain[2][0]:
                        x_list[0] = self.pyhsical_domain[1][0] 
                    elif self.pyhsical_domain[1][0] > self.pyhsical_domain[2][0]:
                        x_list[0] = self.pyhsical_domain[2][0]
                else:
                    g_x = abs(g_x)
                    dx = L_x * (1 - g_x) / (1 - g_x**(self.nodes[0] - 1))
                    dx = -dx
                    if self.pyhsical_domain[1][0] < self.pyhsical_domain[2][0]:
                        x_list[0] = self.pyhsical_domain[2][0] 
                    elif self.pyhsical_domain[1][0] > self.pyhsical_domain[2][0]:
                        x_list[0] = self.pyhsical_domain[1][0]
            for i in range(1,self.nodes[0]):
                x_list[i] = x_list[i-1] + dx * g_x**(i-1)
                x_spacing[:,i-1] = dx * g_x**(i-1)    
            

        if self.pyhsical_domain[0][1] != self.pyhsical_domain[1][1]:
            self.ylength = sorted((self.pyhsical_domain[0][1], self.pyhsical_domain[1][1]))
            L_y = abs((self.pyhsical_domain[1][1] - self.pyhsical_domain[0][1])) #length in y direction
            if g_y == 1:
                dy = L_y / (self.nodes[1] - 1)
            else:
                if g_y > 0:
                    dy = L_y * (1 - g_y) / (1 - g_y**(self.nodes[1] - 1)) 
                    if self.pyhsical_domain[0][1] < self.pyhsical_domain[1][1]:
                        y_list[0] = self.pyhsical_domain[0][1] 
                    elif self.pyhsical_domain[0][1] > self.pyhsical_domain[1][1]:
                        y_list[0] = self.pyhsical_domain[1][1]
                else:
                    g_y = abs(g_y)
                    dy = L_y * (1 - g_y) / (1 - g_y**(self.nodes[1] - 1))
                    dy = -dy
                    if self.pyhsical_domain[0][1] < self.pyhsical_domain[1][1]:
                        y_list[0] = self.pyhsical_domain[1][1] 
                    elif self.pyhsical_domain[0][1] > self.pyhsical_domain[1][1]:
                        y_list[0] = self.pyhsical_domain[0][1]
            for i in range(1,self.nodes[1]):
                y_list[i] = y_list[i-1] + dy * g_y**(i-1)
                y_spacing[:,i-1] = dy * g_y**(i-1)   
            
        elif self.pyhsical_domain[1][1] != self.pyhsical_domain[2][1]:
            self.ylength = sorted((self.pyhsical_domain[1][1], self.pyhsical_domain[2][1]))
            L_y = abs((self.pyhsical_domain[1][1] - self.pyhsical_domain[2][1])) #length in y direction
            if g_y == 1:
                dy = L_y / (self.nodes[1] - 1)
            else:
                if g_y > 0:
                    dy = L_y * (1 - g_y) / (1 - g_y**(self.nodes[1] - 1)) 
                    if self.pyhsical_domain[1][1] < self.pyhsical_domain[2][1]:
                        y_list[0] = self.pyhsical_domain[1][1] 
                    elif self.pyhsical_domain[1][1] > self.pyhsical_domain[2][1]:
                        y_list[0] = self.pyhsical_domain[2][1]
                else:
                    g_y = abs(g_y)
                    dy = L_y * (1 - g_y) / (1 - g_y**(self.nodes[1] - 1))
                    dy = -dy
                    if self.pyhsical_domain[1][1] < self.pyhsical_domain[2][1]:
                        y_list[0] = self.pyhsical_domain[2][1] 
                    elif self.pyhsical_domain[1][1] > self.pyhsical_domain[2][1]:
                        y_list[0] = self.pyhsical_domain[1][1]
            for i in range(1,self.nodes[1]):
                y_list[i] = y_list[i-1] + dy * g_y**(i-1) 
                y_spacing[:,i-1] = dy * g_y**(i-1)   

        if x_list[-1] < x_list[0]:      #Small value to large value to be consistent for physiscs
            x_list = np.flip(x_list)
        if y_list[-1] < y_list[0]:
            y_list = np.flip(y_list)
            
        x_MAT = np.zeros((self.nodes[1], self.nodes[0]), dtype="float")
        y_MAT = np.zeros((self.nodes[1], self.nodes[0]), dtype="float")

        #for uniform grid one can use np.meshgrid()
        for i in range(len(y_list)):
            x_MAT[i,:] = x_list
        for i in range(len(x_list)):
            y_MAT[:,i] = y_list

        self.matricies = [x_MAT, y_MAT]

        print(self.xlength)
        print(self.ylength)

        return x_spacing, y_spacing


    def nonuniform_mesh_2D(self, g_x = 1, g_y = 1):

        N_x, N_y = self.nodes[0], self.nodes[1]
        lengths = {}
        for i in range(len(self.pyhsical_domain)):  #or simply write 4. But I want to keep it as parametric if we might change to polygon. 
            lengths["Lx" + str(i)] = self.pyhsical_domain[(i+1) % len(self.pyhsical_domain)][0] - self.pyhsical_domain[i][0]
            lengths["Ly" + str(i)] = self.pyhsical_domain[(i+1) % len(self.pyhsical_domain)][1] - self.pyhsical_domain[i][1]

        x_MAT = np.zeros((N_y, N_x))
        y_MAT = np.zeros((N_y, N_x))
        dx_MAT = np.zeros((N_y, N_x -1))
        dy_MAT = np.zeros((N_y - 1, N_x))

        # Arrange the given coordinates starting from the lowest values of (x , y) to highest in clockwise direction.
        # pysical_domain_matrix = []
        # ordered_physical_corners = []
        # for i in range(len(self.pyhsical_domain)):
        #     pysical_domain_matrix.append(list(self.pyhsical_domain[i]))

        # print(pysical_domain_matrix)

        # pysical_domain_matrix = np.array(pysical_domain_matrix)

        # print(pysical_domain_matrix[:,0])
        # print(pysical_domain_matrix[:,1])

        # temp_list_x = min(pysical_domain_matrix[:,0]) >= pysical_domain_matrix[:,0] #determining west x coordinates
        # x_indixies = []
        # #findone() def
        # for i in range(len(temp_list_x)):
        #     if temp_list_x[i] == 1:
        #         x_indixies.append(i)
        # if sum(temp_list_x) > 1:  #means that there are two same coordinates check their y values to determine north and south
        #     y1 = pysical_domain_matrix[:,1][x_indixies[0]]
        #     y2 = pysical_domain_matrix[:,1][x_indixies[1]]
        #     print(x_indixies)
        #     print(y1)
        #     print(y2)
        #     if y1 < y2:  # if y1 is higher it is in north. check its index. 
        #         ordered_physical_corners.append(pysical_domain_matrix[x_indixies[0],:])
        #         ordered_physical_corners.append(pysical_domain_matrix[x_indixies[1],:])
        #     elif y2 < y1:
        #         ordered_physical_corners.append(pysical_domain_matrix[x_indixies[1],:])
        #         ordered_physical_corners.append(pysical_domain_matrix[x_indixies[0],:])
        #     else:
        #         IndexError("Points lies on same coordinates")
        # else:
        #     ordered_physical_corners.append(pysical_domain_matrix[x_indixies[0],:])
             
        # temp_list_x = min(pysical_domain_matrix[:,0]) <= pysical_domain_matrix[:,0] #determining east x coordinates
        # if sum(temp_list_x) > 1:  #means that there are two same coordinates check their y values to determine north and south
        #     x_indixies = []
        #     #findone() def
        #     for i in range(len(temp_list_x)):
        #         if temp_list_x[i] == 1:
        #             x_indixies.append(i)
        #     y1 = pysical_domain_matrix[:,1][x_indixies[0]]
        #     y2 = pysical_domain_matrix[:,1][x_indixies[1]]

        #     print(x_indixies)
        #     print(y1)
        #     print(y2)

        #     if y1 < y2:  # if y1 is higher it is in north. check its index. 
        #         ordered_physical_corners.append(pysical_domain_matrix[x_indixies[0],:])
        #         ordered_physical_corners.append(pysical_domain_matrix[x_indixies[1],:])
        #     elif y2 < y1:
        #         ordered_physical_corners.append(pysical_domain_matrix[x_indixies[1],:])
        #         ordered_physical_corners.append(pysical_domain_matrix[x_indixies[0],:])
        #     else:
        #         IndexError("Points lies on same coordinates")

        # print(ordered_physical_corners)


          
        for i in range(N_y):                 # i being the indicies in y direction
            dx_MAT[i, :] = (lengths["Lx0"] + (lengths["Lx0"] + lengths["Lx2"]) * i / (N_y - 1)) / (N_x - 1)

        for j in range(N_x):                 # i being the indicies in y direction
            dy_MAT[:, j] = (lengths["Ly0"] + (lengths["Ly3"] + lengths["Ly1"]) * j / (N_x - 1)) / (N_y - 1)            
            
        # print(lengths)
        # print(dx_MAT)
        # print(dy_MAT)

        for i in range(N_y):
            for j in range(N_x):
                x_MAT[i,j] = x_MAT[i,j] + self.pyhsical_domain[0][0] - dx_MAT[i, j -1]
                y_MAT[i,j] = y_MAT[i,j] + self.pyhsical_domain[0][1] - dy_MAT[i - 1, j]
              
            
        # print(x_MAT)
        # print(y_MAT)

        self.matricies = [x_MAT, y_MAT]

    def Jacobi(self, X_spacing, Y_spacing):
        """
            X_spacing: 1D or 2D ndarray. Spacing matrix of the real pysical domain for X matrix 
            Y_spacing: 1D or 2D ndarray. Spacing matrix of the real pysical domain for Y matrix
        """

        if (self.nodes[1]) != X_spacing.shape[0]:                   #element number in y
            IndexError("Matrix sizes at axis 1 do not match for X spacing")
        if (self.nodes[0] - 1) != X_spacing.shape[1]:               #element number in x
            IndexError("Matrix sizes at axis 0 do not match for X spacing")
        if (self.nodes[1] - 1) != Y_spacing.shape[0]:               #element number in y
            IndexError("Matrix sizes at axis 1 do not match for Y spacing")
        if (self.nodes[1]) != Y_spacing.shape[1]:                   #element number in x
            IndexError("Matrix sizes at axis 0 do not match for Y spacing")
            
        #X_spacing nad Y_spacing are 1D vectors for tihs case

        dXdx = OneDcentraldiff(self.matricies[0], X_spacing)
        dXdy = OneDcentraldiff(self.matricies[0], -Y_spacing, axis=1)   #the minus sign put intentionally. 
        dYdx = OneDcentraldiff(self.matricies[1], X_spacing)            #Matrix convention and coordinate convention do not hold   
        dYdy = OneDcentraldiff(self.matricies[1], -Y_spacing, axis=1)   #the minus sign put intentionally.

        #jacobi is 2D at the moment
        J = np.zeros((self.nodes[1], self.nodes[0], 2, 2))

        for j in range(self.nodes[1]):
            for i in range(self.nodes[0]):
                J_page = np.array([[dXdx[j,i], dXdy[j,i] ],[dYdx[j,i], dYdy[j,i] ]])
                J[j,i,:,:] = J_page

        self.Jacobian = J

    def create_object(self, obj):
        """
           self contains physical domain and nodes number for each.
           
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

        

                        

    def plot2D(self):
        """
            plot mesh grid. 
        """

        fig, ax = plt.subplots()

        x_MAT = self.matricies[0]
        y_MAT = self.matricies[1]

        z = np.ones((self.nodes[1], self.nodes[0])) * 0.5

        image = ax.pcolormesh(x_MAT,np.flip(y_MAT), z, vmin=0, vmax=1, edgecolors="black",linewidth=0.1)
        plt.show()
    


def column_TDMA(a_s, a_w, a_n, a_e, phi, y_index, BC_values, i, N_y, N_x, W, E):

    if i == 0:
        W[1:] = -a_s[y_index[0]:y_index[-1], i]                          
        C = 2 * a_s[y_index[0]:, i] + 2 * a_w[y_index[0]:y_index[-1]+1, i - 1]
        E[:-1] = -a_n[y_index[0]:y_index[-1], i]                            
        Q = 2 * a_e[y_index[0]:y_index[-1]+1,i] * phi[y_index[0]:y_index[-1]+1,(i+1)] + 2 * a_e[y_index[0]:y_index[-1]+1,i]**2 * BC_values['W']   #conditions are set for neumann BC.                
        Q[0] += a_n[0,0] * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * a_n[0,i - 1]**2 * BC_values['N']
        Q[-1] += a_s[0,0] * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * a_s[0,i - 1]**2 * BC_values['S']
      
    if i > 0 and i < N_x-1:
        W[1:] = -a_s[y_index[0]:y_index[-1], i]                          
        C = 2 * a_s[y_index[0]:, i] + 2 * a_w[y_index[0]:y_index[-1]+1, i - 1]
        E[:-1] = -a_n[y_index[0]:y_index[-1], i]                            
        Q = a_w[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i-1)] + a_e[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i+1)] #conditions are set for neumann BC.                
        Q[0] += a_n[0,0] * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * a_n[0,i - 1]**2 * BC_values['N']
        Q[-1] += a_s[0,0] * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * a_s[0,i - 1]**2 * BC_values['S']

    if i == N_x - 1:
        W[1:] = -a_s[y_index[0]:y_index[-1], i]                          
        C = 2 * a_s[y_index[0]:, i] + 2 * a_w[y_index[0]:y_index[-1]+1, i - 1]
        E[:-1] = -a_n[y_index[0]:y_index[-1], i]                            
        Q = 2 * a_w[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i-1)] + 2 * a_w[y_index[0]:y_index[-1]+1,i-1]**2 * BC_values['E']             
        Q[0] += a_n[0,0] * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * a_n[0,i - 1]**2 * BC_values['N']
        Q[-1] += a_s[0,0] * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * a_s[0,i - 1]**2 * BC_values['S']
    

        Q = np.flip(Q)                     #The reason of reversing Q is, existing Q is inconsistent with the W and E and C list.

        W[-1] += -a_s[y_index[-1], i] * (y_index[0] == 0)             #Neumann of N-S boundaries. implemented here. 
        E[0] += -a_n[y_index[0], i] * (y_index[-1] == N_y - 1)        #probabaly for different spacing matrixies the east and west should fliped

        if y_index[0] == 0:
            phi[y_index[-1]::-1,i] = TDMA(W,C,E,Q) 
        else:
            phi[y_index[-1]:y_index[0] - 1:-1,i] = TDMA(W,C,E,Q) 


    return phi


class PDE_2D_Solver:
    """
        Solves 2D partial differential equation with given boundary condiations.

        Takes Mesh class.
        BC: dictionary

        Boundary condiations:

        BC = {'W': BC1, 'S': BC2, 'E': BC3, 'N': BC4 }
        Each BC should be string, 'D' or 'N' for Diriclet or Neumann
    """

    def __init__(self, mesh ,BC):
        self.BC = BC
        self.BCvalues = None
        self.solution = None
        self.velocity = None
        self.mesh = mesh

    def solver(self, BC_values, itteration_type  = "column"):
        """
        Construct unknown matrix with its boundary condiations
        
        BC_values: dict 

        BC_values: {'W': BC1, 'S': BC2, 'E': BC3, 'N': BC4}

        BC's are: float or ndarray

        """

        self.BCvalues = BC_values
        mesh = self.mesh

        (N_y, N_x) = np.shape(mesh.matricies[0])
        
        phi = np.zeros((N_y, N_x))       #unknown
        phi_new = np.zeros((N_y, N_x))   #unknown for storing new values
        source = np.zeros((N_y, N_x))

        #Coefficient Matricies 
        #Same if it is uniform-block mesh.
        #Same in row or column if it is non-uniform-block mesh
        #Different for row and coulmn for non-uniform mesh

        #the matrix should include the spacing information. But at the moent I cant configure that if we will apply solver for unstructerde mesh or after Jacobi is defined. Thus I am not 
        #hurrying to create unstrutred mesh and its spacing matrix. 

        a_e = np.zeros((N_y, N_x - 1),dtype="float") + (mesh.matricies[1][1,0] - mesh.matricies[1][0,0])**2
        a_w = np.zeros((N_y, N_x - 1),dtype="float") + (mesh.matricies[1][1,0] - mesh.matricies[1][0,0])**2
        a_n = np.zeros((N_y - 1, N_x),dtype="float") + (mesh.matricies[0][0,1] - mesh.matricies[0][0,0])**2
        a_s = np.zeros((N_y - 1, N_x),dtype="float") + (mesh.matricies[0][0,1] - mesh.matricies[0][0,0])**2

        x_index = list(range(N_x))
        y_index = list(range(N_y))

        #this spacing updates should be more explainable. 
        if self.BC['W'] == "D":
            phi[:,0] = BC_values['W']
            x_index = x_index[1:] 
        if self.BC['S'] == "D":
            phi[-1,:] = BC_values['S']
            y_index = y_index[:-1]
        if self.BC['E'] == "D":
            phi[:,-1] = BC_values['E']
            x_index = x_index[:-1]
        if self.BC['N'] == "D":
            phi[0,:] = BC_values['N']
            y_index = y_index[1:]

        if self.BC['W'] == "N":
            # a_w[:,0] += 1 * a_e[:,-1]
            pass
        if self.BC['S'] == "N":
            a_s =  np.concatenate((a_s, np.array([a_s[-1,:]])), axis=0)
        if self.BC['E'] == "N":
            # a_e[:,-1] += 1 * a_w[:,0]
            pass
        if self.BC['N'] == "N":
            a_n =  np.concatenate((a_n, np.array([a_n[1,:]])), axis=0)
            
        #Column by Column TDMA
        """
            Starting from first column that is unknown. Thus, 
            first we solve phi[:,1]
            than pass to phi[:,2]

            We will solve again this loop. 
        """
        W = np.zeros(len(y_index), dtype="float")
        E = np.zeros(len(y_index), dtype="float")

        for t in range(30):

            # W = np.zeros(len(y_index), dtype="float")
            # E = np.zeros(len(y_index), dtype="float")

            for i in x_index:

                # phi = column_TDMA(a_s, a_w, a_n, a_e, phi, y_index, BC_values, i, N_y, N_x, W, E)

                if i == 0:
                    W[1:] = -a_s[y_index[0]:y_index[-1], i]                          
                    C = 2 * a_s[y_index[0]:, i] + 2 * a_w[y_index[0]:y_index[-1]+1, i - 1]
                    E[:-1] = -a_n[y_index[0]:y_index[-1], i]                            
                    Q = 2 * a_e[y_index[0]:y_index[-1]+1,i] * phi[y_index[0]:y_index[-1]+1,(i+1)] + 2 * a_e[y_index[0]:y_index[-1]+1,i]**2 * BC_values['W']   #conditions are set for neumann BC.                
                    Q[0] += a_n[0,0] * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * a_n[0,i - 1]**2 * BC_values['N']
                    Q[-1] += a_s[0,0] * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * a_s[0,i - 1]**2 * BC_values['S']
                
                if i > 0 and i < N_x-1:
                    W[1:] = -a_s[y_index[0]:y_index[-1], i]                          
                    C = 2 * a_s[y_index[0]:, i] + 2 * a_w[y_index[0]:y_index[-1]+1, i - 1]
                    E[:-1] = -a_n[y_index[0]:y_index[-1], i]                            
                    Q = a_w[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i-1)] + a_e[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i+1)] #conditions are set for neumann BC.                
                    Q[0] += a_n[0,0] * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * a_n[0,i - 1]**2 * BC_values['N']
                    Q[-1] += a_s[0,0] * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * a_s[0,i - 1]**2 * BC_values['S']
                
                if i == N_x - 1:
                    W[1:] = -a_s[y_index[0]:y_index[-1], i]                          
                    C = 2 * a_s[y_index[0]:, i] + 2 * a_w[y_index[0]:y_index[-1]+1, i - 1]
                    E[:-1] = -a_n[y_index[0]:y_index[-1], i]                            
                    Q = 2 * a_w[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i-1)] + 2 * a_w[y_index[0]:y_index[-1]+1,i-1]**2 * BC_values['E']             
                    Q[0] += a_n[0,0] * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * a_n[0,i - 1]**2 * BC_values['N']
                    Q[-1] += a_s[0,0] * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * a_s[0,i - 1]**2 * BC_values['S']
                
                
                Q = np.flip(Q)                     #The reason of reversing Q is, existing Q is inconsistent with the W and E and C list.

                W[-1] += -a_s[y_index[-1], i] * (y_index[0] == 0)             #Neumann of N-S boundaries. implemented here. 
                E[0] += -a_n[y_index[0], i] * (y_index[-1] == N_y - 1)        #probabaly for different spacing matrixies the east and west should fliped

                if y_index[0] == 0:
                    phi[y_index[-1]::-1,i] = TDMA(W,C,E,Q) 
                else:
                    phi[y_index[-1]:y_index[0] - 1:-1,i] = TDMA(W,C,E,Q) 


        self.solution = phi





    def velocityfield(self):
        """
            grad(phi) = d (phi) / dx  + d (phi) / dy =  0
            Continiutiy in fluids. 

            Thus:
            u = d (phi) / dx
            v = d (phi) / dy

        """

        X, Y = self.mesh.matricies[0], self.mesh.matricies[1] 

        (N_y, N_x) = np.shape(X)

        dx = X[0,1] - X[0,0]
        dy = Y[0,0] - Y[1,0]

        (u, v) = TwoDcentraldiff_simple(self.solution, dx, dy)

        W = np.zeros((N_y, N_x, 3))

        W[:,:,0] = u
        W[:,:,1] = v

        self.velocity = W

    def plot2D(self, colormap = "hot"):
        """
            plot Z for domain of the X and Y. 
        """

        fig, ax = plt.subplots()

        twin2 = ax.twinx()
        
        x_MAT = self.mesh.matricies[0]
        y_MAT = self.mesh.matricies[1]

        z_min = self.solution.min()
        z_max = self.solution.max()

        image = ax.pcolormesh(x_MAT,np.flip(y_MAT), self.solution, vmin=z_min, vmax=z_max, edgecolors="black", cmap=colormap, linewidth=0.1)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="6%", pad="2%")

        if self.BC["S"] == "D":
            ax.set_xlabel(('$\phi$ = ' + str(self.BCvalues["S"])) , fontsize=20)
        else:
            ax.set_xlabel(('${\partial \phi}{\partial y}$ = '+ str(self.BCvalues["S"])) , fontsize=20)
        if self.BC["W"] == "D":
            ax.set_ylabel(('$\phi$ = '+ str(self.BCvalues["W"])) , fontsize=20)
        else:
            ax.set_ylabel(('${\partial \phi}{\partial x}$ = ' + str(self.BCvalues["W"])) , fontsize=20)
        if self.BC["N"] == "D":
            ax.set_title(('$\phi$ = ' + str(self.BCvalues["N"])) , fontsize=20)
        else:
            ax.set_title(('$\partial \phi \partial y$ = ' + str(self.BCvalues["N"])) , fontsize=20)
 
# secondary y-axis label
        
        if self.BC["E"] == "D":
            twin2.set_ylabel(('$\phi$ = ' + str(self.BCvalues["E"])) , fontsize=20)
        else:
            twin2.set_ylabel(('$\partial \phi}{\partial y}$ = ' +  str(self.BCvalues["E"])) , fontsize=20)
        

        plt.colorbar(image, cax=cax)
        plt.show()

    def stream(self, streamcolor = "blue"):
        """
            Plots streamplot for velocity of the solution
        """

        x_MAT = self.mesh.matricies[0]
        y_MAT = self.mesh.matricies[1]

        u = self.velocity[:,:,0]
        v = self.velocity[:,:,1]

        plt.streamplot(x_MAT, y_MAT, -u, v, color=streamcolor, linewidth=2, cmap='autumn')
        plt.show()

    def countour(self):
        pass














