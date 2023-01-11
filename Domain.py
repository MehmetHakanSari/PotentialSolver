import numpy as np
from Differentials import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from object import object
from visiual import *
from time import perf_counter



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
            self.xlength = sorted((self.pyhsical_domain[0][0], self.pyhsical_domain[1][0]))
            x_list = np.linspace(self.pyhsical_domain[0][0], self.pyhsical_domain[1][0], self.nodes[0])
        elif self.pyhsical_domain[1][0] != self.pyhsical_domain[2][0]:
            self.xlength = sorted((self.pyhsical_domain[1][0], self.pyhsical_domain[2][0]))
            x_list = np.linspace(self.pyhsical_domain[1][0], self.pyhsical_domain[2][0], self.nodes[0])

        if self.pyhsical_domain[0][1] != self.pyhsical_domain[1][1]:
            #if the given coordinates are in same line for x coordinate pass other coordinate
            self.ylength = sorted((self.pyhsical_domain[0][1], self.pyhsical_domain[1][1]))
            y_list = np.linspace(self.pyhsical_domain[0][1], self.pyhsical_domain[1][1], self.nodes[1])
        elif self.pyhsical_domain[1][1] != self.pyhsical_domain[2][1]:
            self.ylength = sorted((self.pyhsical_domain[1][1], self.pyhsical_domain[2][1]))
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
        self.xspacing = x_spacing
        self.yspacing = y_spacing

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
                x_spacing[:,0:i-1] = dx * g_x**(i-1)  
            
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
        self.xspacing = x_spacing
        self.yspacing = y_spacing

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
    

class ElippticMesh:

    def __init__(self, nodes, GAMA1, GAMA2):
        """
            nodes: 1D ndarray. Number of nodes in zeta and eta direction
            GAMA1: 2D ndarray. GAMA1 matrix contains x and y coordinates of the physical domain curve. Row length equals N_zeta. It is inner boundary. 
            GAMA2: 2D ndarray. GAMA2 matrix contains x and y coordinates of the physical domain curve. Row length equals N_zeta. It is outer boundary.
        """
        self.nodes = nodes
        self.GAMA1 = GAMA1
        self.GAMA2 = GAMA2

    
    def create_elipticmesh(self):
        """
            Create eliptic mesh by gauss-seidel itteration with SOR. Returns self X and Y matrixies inlucding coordinates of the nodal points. 
        """
        #create GAMA3 and GAMA4 boundaries. These boundaries are conntecting inner and outer bondarr togather. 
        #GAMA3 and 4 takes their initial values from the inner and outer boundaries. 

        GAMA1 = self.GAMA1
        GAMA2 = self.GAMA2

        GAMA3 = np.array([[GAMA2[0,0], GAMA2[0,1]],[GAMA1[0,0], GAMA1[0,1]]])
        GAMA4 = np.array([[GAMA2[-1,0], GAMA2[-1,1]],[GAMA1[-1,0], GAMA1[-1,1]]]) # it can be -1 also. Not sure about this.

        GAMA3_x = np.linspace(GAMA3[0,0], GAMA3[1,0], self.nodes[1])
        GAMA3_y = np.linspace(GAMA3[0,1], GAMA3[1,1], self.nodes[1])

        GAMA4_x = np.linspace(GAMA4[0,0], GAMA4[1,0], self.nodes[1])
        GAMA4_y = np.linspace(GAMA4[0,1], GAMA4[1,1], self.nodes[1])

        GAMA3 = np.array([GAMA3_x, GAMA3_y]).T
        GAMA4 = np.array([GAMA4_x, GAMA4_y]).T

        # alpha = np.zeros((self.nodes[1], self.nodes[0]))
        # beta = np.zeros((self.nodes[1], self.nodes[0]))
        # gamma = np.zeros((self.nodes[1], self.nodes[0]))

        #create X and Y matrixes and give boundary conditions from GAMA1, GAMA2, GAMA3 and GAMA4

        X = np.zeros((self.nodes[1], self.nodes[0]))
        Y = np.zeros((self.nodes[1], self.nodes[0]))

        X[0,:] = GAMA1[:,0]
        X[-1,:] = GAMA2[:,0]
        X[:,0] = np.flip(GAMA3[:,0])
        X[:,-1] = np.flip(GAMA4[:,0])
 
        Y[0,:] = GAMA1[:,1]
        Y[-1,:] = GAMA2[:,1]
        Y[:,0] = np.flip(GAMA3[:,1])
        Y[:,-1] = np.flip(GAMA4[:,1])

        error_x, error_y = 1, 1

        maxiteration = 10000
        message = 1000

        N_z = self.nodes[0]
        N_e = self.nodes[1]

        X_new, Y_new = X.T, Y.T

        #interpolate interiour nodes
        for i in range(N_z):
            X[1:-1, i] = np.linspace(X[0, i], X[-1, i], N_e)[1:-1]
            Y[1:-1, i] = np.linspace(Y[0, i], Y[-1, i], N_e)[1:-1]

        # self.X = X
        # self.Y = Y
        # # self.plot_mesh()

        for iteration in range(maxiteration):

            X = X.T
            Y = Y.T

            X_temp = np.append([X[-2, :].copy()], X[0:2, :].copy(), 0) 
            Y_temp = np.append([X[-2, :].copy()], Y[0:2, :].copy(), 0)
            alpha, beta, gamma = Solve_Coeff(X_temp , Y_temp)
            X[0, 1:-1] = SolveEliptic(alpha, beta, gamma, X_temp )
            Y[0, 1:-1] = SolveEliptic(alpha, beta, gamma, Y_temp)

            X[-1, 1:-1] = X[0, 1:-1].copy()
            Y[-1, 1:-1] = Y[0, 1:-1].copy()

            alpha, beta, gamma = Solve_Coeff(X, Y)
            X_new[1:-1, 1:-1] = SolveEliptic(alpha,beta, gamma, X)
            Y_new[1:-1, 1:-1] = SolveEliptic(alpha,beta, gamma, Y)

            error_x = np.max(np.abs(X_new - X))
            error_y = np.max(np.abs(Y_new - Y))
            
            if error_x < 1e-6 and error_y < 1e-6:
                break

            #give message with some interval 
            if iteration % message == 0:
                print("Iteration: ", iteration, "Error is: ", error_x, error_y)

            if iteration == maxiteration - 1:
                print("Max iteration reached. Error is: ", error_x, error_y)

            X = X_new.T.copy()
            Y = Y_new.T.copy()

        self.X = X.T
        self.Y = Y.T
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma


    def plot_mesh(self):
        """
            Plot mesh
        """
        plt.figure(figsize=(15,15))
        
        for i in range(self.nodes[0]):
            plt.plot(self.X[:,i], self.Y[:,i], 'k')
        for i in range(self.nodes[1]):
            plt.plot(self.X[i,:], self.Y[i,:], 'k')
        plt.show()
















