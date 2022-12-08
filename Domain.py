import numpy as np
from Differentials import * 


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
            
        x_MAT = np.zeros((self.nodes[1], self.nodes[0]))
        y_MAT = np.zeros((self.nodes[1], self.nodes[0]))

        #for uniform grid one can use np.meshgrid()
        for i in range(len(y_list)):
            x_MAT[i,:] = x_list
        for i in range(len(x_list)):
            y_MAT[:,i] = y_list

        self.matricies = [x_MAT, y_MAT]

    def nonuniform_block_mesh_2D(self, g_x = 1, g_y = 1):
        #g_x is the gradient for delta_x. Float. Default value 1. 
        #g_y is the gradient for delta_Y. Float. Default value 1. 

        #grid size changes for each step. 

        x_list = np.zeros(self.nodes[0], dtype="float")
        y_list = np.zeros(self.nodes[1], dtype="float")

        if self.pyhsical_domain[0][0] != self.pyhsical_domain[1][0]:
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
            
        elif self.pyhsical_domain[1][0] != self.pyhsical_domain[2][0]:
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
            

        if self.pyhsical_domain[0][1] != self.pyhsical_domain[1][1]:
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
            
        elif self.pyhsical_domain[1][1] != self.pyhsical_domain[2][1]:
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
    

class PDE_2D_Solver:
    """
        Solves 2D partial differential equation with given boundary condiations.

        Takes Mesh class.
        BC: dictionary

        Boundary condiations:

        BC = {'W': BC1, 'S': BC2, 'E': BC3, 'N': BC4 }
        Each BC should be string, 'D' or 'N' for Diriclet or Neumann
    """

    def __init__(self, BC):
        self.BC = BC
        self.solution = None
        self.vector_field = None



    def solver(self, BC_values, mesh):
        """
        Construct unknown matrix with its boundary condiations
        
        BC_values: dict 

        BC_values: {'W': BC1, 'S': BC2, 'E': BC3, 'N': BC4}

        BC's are: float or ndarray

        """

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
            # a_n[0,:] += 1 * a_s[-1,:]
            pass
            
        # print(a_e)
        # print(a_w)
        # print(a_s)
        # print(a_n)
        # print(" ")
        print(phi)
        
        #Column by Column TDMA
        """
            Starting from first column that is unknown. Thus, 
            first we solve phi[:,1]
            than pass to phi[:,2]

            We will solve again this loop. 
        """
        W = np.zeros(len(y_index), dtype="float")
        E = np.zeros(len(y_index), dtype="float")

        
        for t in range(150):

            for i in x_index:
                W[1:] = -a_s[y_index[0]:y_index[-1], i]                          #Neumann BC. implemented here 
                C = 2 * a_s[y_index[0]:, i] + 2 * a_w[y_index[0]:y_index[-1]+1, i - 1]
                E[:-1] = -a_n[y_index[0]:y_index[-1], i]                            #Neumann BC. implemented here
                Q = a_w[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i-1) + 2 * (i-1 == -1)] + a_e[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i+1) - 2 * (i+1 == N_x)] #conditions are set for neumann BC.
           
                Q[0] += a_n[0,0] *  phi[0,i-1] * (y_index[0] == 1) +  (y_index[0] == 0) * 2 * a_w[0,i - 1] * BC_values['N']
                Q[-1] += a_s[0,0] *  phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * a_w[0,i - 1] * BC_values['S']
                Q = np.flip(Q)                     #The reason of reversing Q is, existing Q is inconsistent with the W and E and C list.

                W[-1] +=  -a_s[y_index[-1], i] * (y_index[0] == 0)
                E[0] += -a_n[y_index[0], i] * (y_index[-1] == N_y - 1) 

                if y_index[0] == 0:
                    phi[y_index[-1]::-1,i] = TDMA(W,C,E,Q) 
                else:
                    phi[y_index[-1]:y_index[0] - 1:-1,i] = TDMA(W,C,E,Q) 



        self.solution = phi


    def vector_field(self):
        """
            grad(phi) = d (phi) / dx  + d (phi) / dy =  0
            Continiutiy in fluids. 

            Thus:

            u = d (phi) / dx
            v = d (phi) / dy


        """
        # dx, dy
        # (u, v) = TwoDCentralDiff(self.solution, dx, dy)
















