import numpy as np

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

        x_list = np.zeros(self.nodes[0])
        y_list = np.zeros(self.nodes[1])

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
            
        x_MAT = np.zeros((self.nodes[1], self.nodes[0]))
        y_MAT = np.zeros((self.nodes[1], self.nodes[0]))

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

        for i in range(N_y):                 # i being the indicies in y direction
            dx_MAT[i, :] = (lengths["Lx3"] + (lengths["Lx0"] + lengths["Lx2"]) * i / (N_y - 1)) / (N_x - 1)

        for j in range(N_x):                 # i being the indicies in y direction
            dy_MAT[:, j] = (lengths["Ly0"] + (lengths["Ly3"] + lengths["Ly1"]) * j / (N_x - 1)) / (N_y - 1)            
            
        print(lengths)
        print(dx_MAT)
        print(dy_MAT)

        for i in range(N_y):
            for j in range(N_x):
                x_MAT[i,j] = x_MAT[i,j] + self.pyhsical_domain[0][0] - dx_MAT[i, j -1]
                y_MAT[i,j] = y_MAT[i,j] + self.pyhsical_domain[0][1] - dy_MAT[i - 1, j]
              
            
        print(x_MAT)
        print(y_MAT)

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

    def TDMA(W,C,E,Q):
        """
            W: west
            C: center
            E: east
            Q: source
        """
        X = np.zeros(len(Q))
        C[1:] = C[1:] - E[0:-1] * W[1:] / C[0:-1]
        Q[1:] = Q[1:] - Q[0:-1] * W[1:] / C[0:-1]
        X[-1] = Q[-1] / C[-1]
        for i in range(len(Q),0,-1):
            X[i] = (Q[i] - E[i] * X[i+1]) / C[i]
        
        return X

    def solver(self, BC_values, mesh):
        """
        Construct unknown matrix with its boundary condiations
        
        BC_values: dict 

        BC_values: {'W': BC1, 'S': BC2, 'E': BC3, 'N': BC4}

        BC's are: float or ndarray

        """

        (N_y, N_x) = np.shape(mesh.matricies[0])
        

        phi = np.zeros((N_x, N_y))       #unknown
        phi_new = np.zeros((N_x, N_y))   #unknown for storing new values
        source = np.zeros((N_x, N_y))

        #Coefficient Matricies 
        #Same if it is uniform-block mesh.
        #Same in row or column if it is non-uniform-block mesh
        #Different for row and coulmn for non-uniform mesh
        a_e = np.zeros((N_y, N_x - 1)) + (mesh.matricies[0][0,1] - mesh.matricies[0][0,0])**2
        a_w = np.zeros((N_y, N_x - 1)) + (mesh.matricies[0][0,1] - mesh.matricies[0][0,0])**2
        a_n = np.zeros((N_y - 1, N_x)) + (mesh.matricies[1][1,0] - mesh.matricies[1][0,0])**2
        a_s = np.zeros((N_y - 1, N_x)) + (mesh.matricies[1][1,0] - mesh.matricies[1][0,0])**2


        if self.BC['W'] == "D":
            phi[0,:] = BC_values['W']
        if self.BC['S'] == "D":
            phi[:,-1] = BC_values['S']
        if self.BC['E'] == "D":
            phi[-1,:] = BC_values['E']
        if self.BC['N'] == "D":
            phi[:,0] = BC_values['N']

        if self.BC['W'] == "N":
            a_w[:,0] += 1
        if self.BC['S'] == "N":
            a_s[-1,:] += 1
        if self.BC['E'] == "N":
            a_e[:,-1] += 1
        if self.BC['N'] == "N":
            a_n[0,:] += 1
            
        print(a_e)
        print(a_w)
        print(a_s)
        print(a_n)
        
        #Column by Column TDMA
        """
            Starting from first column that is unknown. Thus, 
            first we solve phi[:,1]
            than pass to phi[:,2]
            than pass to phi[:,3]
            than pass to phi[:,4]

            We will solve again this loop. 
        """
        W = np.zeros(N_y - 2)
        E = np.zeros(N_y - 2)
        for i in range(1, N_x-1):
            W[1:] = a_s[:, i]
            C = 2 * a_s[1:, i] + 2 * a_w[:, i - 1]
            E[:-1] = a_n[:, i]
            Q = a_w[:,i-1] * phi[:,i-1] + a_e[:,i] * phi[:,i+1] 
            Q[0] += a_n[0,0] *  phi[0,i-1]
            Q[-1] += a_s[0,0] *  phi[-1,i-1] 
            phi[1:-1,i] = self.TDMA(W,C,E,Q)

            print(phi)

            if i == 1:
                break















