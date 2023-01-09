import numpy as np
from Differentials import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from object import object
from visiual import *
from time import perf_counter


def nodebynode(x_index, y_index, x_spacing, y_spacing, BCvalues, phi, phi_old, property_map, N_x, N_y, omega,type = "stream"):    
    if type == "stream":
        for i in x_index:

            if i == 0:  #if the west boundary is given neumann

                for j in y_index:
                    
                    if j == 0:   #if the north boundary is given neumann

                        dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2
                        
                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 

                    if j > 0 and j < N_y - 1:
                        
                        dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 

                    if j == N_y - 1:   #if the south boundary is given neumann

                        dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 


            elif i > 0 and i < N_x - 1:  #if the west boundary is given neumann

                for j in y_index:
                    
                    if j == 0:   #if the north boundary is given neumann

                        C_pro = (not(property_map[j,i] == -2)) * 1
                        W_pro = (not(property_map[j,i-1] == -2)) * 1
                        E_pro = (not(property_map[j,i+1] == -2)) * 1
                        S_pro = (not(property_map[j+1,i] == -2)) * 1

                        if C_pro == 0:
                            phi[j,i] = 0
                        else:
                            dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                            dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2  
                            
                            phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] * E_pro + phi[j, i-1] * W_pro) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i] * S_pro) 

                    if j > 0 and j < N_y - 1:

                        C_pro = (not(property_map[j,i] == -2)) * 1
                        W_pro = (not(property_map[j,i-1] == -2)) * 1
                        E_pro = (not(property_map[j,i+1] == -2)) * 1
                        S_pro = (not(property_map[j+1,i] == -2)) * 1
                        N_pro = (not(property_map[j-1,i] == -2)) * 1

                        if C_pro == 0:
                            phi[j,i] = 0
                        else:
                            dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                            dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2  

                            phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] * E_pro + phi[j, i-1] * W_pro) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] * S_pro + phi[j-1, i] * N_pro) 

                    if j == N_y - 1:

                        C_pro = (not(property_map[j,i] == -2)) * 1
                        W_pro = (not(property_map[j,i-1] == -2)) * 1
                        E_pro = (not(property_map[j,i+1] == -2)) * 1
                        N_pro = (not(property_map[j-1,i] == -2)) * 1

                        if C_pro == 0:
                            phi[j,i] = 0
                        else:
                            dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                            dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2

                            phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] * E_pro + phi[j, i-1] * W_pro) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i] * N_pro) 


            elif i == N_x - 1:

                for j in y_index:

                    # property_map(j,i)
                    # property_map(j,i-1)
                    # property_map(j+1,i)
                    # property_map(j-1,i)
                    
                    if j == 0:   #if the north boundary is given neumann
                        
                        dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 

                    if j > 0 and j < N_y - 1:

                        dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 

                    if j == N_y - 1:

                        dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i])
    
    elif type == "potensial":

        for i in x_index:

            if i == 0:  #if the west boundary is given neumann

                for j in y_index:

                    if j == 0:   #if the north boundary is given neumann

                        dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2
                        
                        phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 

                    if j > 0 and j < N_y - 1:
                        
                        dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2  

                        phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 

                    if j == N_y - 1:   #if the south boundary is given neumann

                        dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2

                        phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 


            elif i > 0 and i < N_x - 1:  #interior points   

                for j in y_index:
                    
                    if j == 0:   #if the north boundary is given neumann

                        C_pro = (property_map[j,i] != -2) * 1   #if it is wall check its sides. If the interior is in north or south 
                        #j will be approximated. If interior is west or east i will be approximated

                        if C_pro == 0:
                            phi[j,i] = 0
                        else:
                            dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                            dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2  
                            
                            phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 

                    if j > 0 and j < N_y - 1:

                        dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2  

                        # if (j == 66) and (i == 36):
                            # print(property_map[j,i])
                        if (property_map[j,i] == -2):    #if it is a wall

                            #determine which sides are in the interior
                            
                            east = (property_map[j,i+1] == -1)
                            west = (property_map[j,i-1] == -1)  
                            south = (property_map[j+1,i] == -1) 
                            north = (property_map[j-1,i] == -1) 

                            # if (j == 66) and (i == 36):
                                # print(east, west, south, north)

                            if east:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                                elif north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                            elif west:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                                elif north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                            else:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                                    # if (j == 66) and (i == 36):
                                        # print(phi[j,i])
                                elif north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                        elif (property_map[j,i] == -1):
                            phi[j,i] = 0
                        else:
                            phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                            

                    if j == N_y - 1:

                        C_pro = (not(property_map[j,i] == -2)) * 1

                        if C_pro == 0:
                            phi[j,i] = 0
                        else:
                            dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                            dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2

                            phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 


            elif i == N_x - 1:

                for j in y_index:
                    
                    if j == 0:   #if the north boundary is given neumann
                        
                        dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 

                    if j > 0 and j < N_y - 1:

                        dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 

                    if j == N_y - 1:

                        dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i])

        phi = ((property_map != -1) * 1) * phi 
        phi = (1 - omega) * phi_old + omega * phi
                    
    return phi 
    

def column_TDMA(a_s, a_w, a_n, a_e, phi, y_index, BC_values, x_index, N_y, N_x, W, E):
    for i in x_index:
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

class eliptic_PDE_solver:


    def __init__(self, mesh, BC_values, BC_type = None):

        self.mesh = mesh
        self.BCvalues = BC_values
        self.BC_type = BC_type

    def solver(self):

        N_z = self.mesh.nodes[0]
        N_e = self.mesh.nodes[1]
        alpha = self.mesh.alpha
        beta = self.mesh.beta
        gamma = self.mesh.gamma

        psi = np.zeros((N_e, N_z))  # initial guess for psi
        psi_old = np.zeros((N_e, N_z))  # initial guess for psi

        psi[0,:] = np.ones((1, N_z)) * self.BCvalues['In']
        psi[-1,:] = np.ones((1, N_z)) * self.BCvalues['Out']
        psi[:,0] = np.zeros((N_e)) * self.BCvalues['Cut1']
        psi[:,-1] = np.zeros((N_e)) * self.BCvalues['Cut2']

        maxiteration = 5000
        tolerance = 1e-6
        message = 1000

        for iteration in range(maxiteration):

            psi[1:N_e-1, 1:N_z-1] = ((-0.5) / (alpha[1:N_e-1, 1:N_z-1] + gamma[1:N_e-1, 1:N_z-1] + 1e-9)) * (2 * beta[1:N_e-1, 1:N_z-1] * (psi[2:N_e, 2:N_z] - psi[2:N_e, 0:N_z-2] - psi[0:N_e-2, 2:N_z] + psi[0:N_e-2, 0:N_z-2]) - alpha[1:N_e-1, 1:N_z-1] * (psi[1:N_e-1, 2: N_z] + psi[1:N_e-1 , 0:N_z-2]) - gamma[1:N_e-1, 1:N_z-1] * (psi[2:N_e, 1:N_z-1] + psi[0:N_e-2, 1:N_z-1]))
        
            residual_psi = np.max(np.abs(psi - psi_old))

            if residual_psi < tolerance:
                print("Psi is calculated with residual: ", residual_psi)
                break

            #give message with some interval 
            if iteration % message == 0:
                print("Iteration: ", iteration, "Error is: ", residual_psi)

            if iteration == maxiteration - 1:
                print("Max iteration reached. Error is: ", residual_psi)

            psi_old = psi.copy()

        self.solution = psi

    def velocity_field(self):
        pass


    def plot2D(self):
            
        fig = plt.figure()

        X, Y = self.mesh.X, self.mesh.Y
        #use pcolormesh for better visualization
        plt.pcolormesh(X, Y, self.solution, cmap = 'jet')
        plt.colorbar()
        plt.show()

    def streamplot(self):
        """
            Plots streamplot for velocity of the solution
        """

        x_MAT = self.mesh.X
        y_MAT = self.mesh.Y

        u = self.velocity[:,:,0]
        v = self.velocity[:,:,1]

        speed = np.sqrt(u*u + v*v)
        lw = 1.2*speed/speed.max()
        lw = 1

        plt.streamplot(x_MAT, y_MAT, -u, v, color=streamcolor, density=0.9, linewidth=lw, cmap='winter')
        plt.show()


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

    def solver(self, BC_values, variable, map, omega, toll, itteration_type  = "column"):
        """
        Construct unknown matrix with its boundary condiations
        
        BC_values: dict 

        BC_values: {'W': BC1, 'S': BC2, 'E': BC3, 'N': BC4}

        BC's are: float or ndarray

        variable: Solve according to stream or potensial. Change it after defining wall BC.

        map: map is the properties of the physical domain. It should be called in main script as map.area() 

        """

        self.variable = variable
        self.BCvalues = BC_values
        self.map = map
        mesh = self.mesh
        # property_map = mesh.map()
        property_map = map    #name of the given input can be change like this afterwards 

        # print(property_map)

        (N_y, N_x) = np.shape(mesh.matricies[0])

        if type(mesh.xspacing) == float:
            x_spacing = np.ones((N_y, N_x - 1),dtype="float") * mesh.xspacing
            y_spacing = np.ones((N_y - 1, N_x),dtype="float") * mesh.yspacing
        else:
            x_spacing = np.zeros((N_y, N_x - 1),dtype="float")
            y_spacing = np.zeros((N_y - 1, N_x),dtype="float")

            for i in range(N_y):
                x_spacing[i,:] = mesh.xspacing
            for i in range(N_x):
                y_spacing[:,i] = mesh.yspacing.reshape(N_y-1, )

        
        # phi = np.zeros((N_y, N_x))       #unknown
        phi = np.ones((N_y, N_x)) * np.linspace(BC_values['W'], BC_values['E'], N_x)     #unknown
        phi_old = np.zeros((N_y, N_x))   #unknown for storing new values
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

        # print(y_spacing)
        # print(x_spacing)
        # print(y_index)
        # print(x_index)
        

        message = 100
        mass_old = 1

        start = perf_counter()    
        for t in range(1, 500001):

            if itteration_type == "column":
                phi = column_TDMA(a_s, a_w, a_n, a_e, phi, y_index, BC_values, x_index, N_y, N_x, W, E)
            elif itteration_type == "nodebynode":
                phi = nodebynode(x_index, y_index, x_spacing, y_spacing, self.BCvalues, phi, phi_old, property_map.area, N_x, N_y, omega, type = variable)

            self.solution = phi.copy()

            #calculate the error
            mass = self.mass_conservation()
            residual = np.sum(np.abs(phi - phi_old))
            mass_residaul = abs(mass - mass_old) / mass_old

            if (t) % message == 0:
                # mass = self.mass_conservation()
                print("Residual: ", residual,"Mass Residual:", mass_residaul ,"Mass:", mass, " at ", t, "th iteration", "  Time: ", perf_counter() - start)

            if residual < toll or abs(mass_residaul) < toll or abs(mass) < toll:
                print("Residual: ", residual, "Mass:", mass)
                print("Solution converged at ", t, "th iteration")
                break

            phi_old = phi.copy()
            mass_old = mass
            
 

    def mass_conservation(self):
        """
            This function is for checking the mass conservation. 
            It is not necessary for the solver. 
        """
        self.velocityfield("potensial")
        u = self.velocity[:,:,0]
        v = self.velocity[:,:,1]

        mass = np.sum(u[:,0]) - np.sum(u[:,-1]) + np.sum(u[0,:]) - np.sum(u[-1,:]) + np.sum(v[:,0]) - np.sum(v[:,-1]) + np.sum(v[0,:]) - np.sum(v[-1,:])

        return mass


    def velocityfield(self, type = "potensial"):
        """
            grad(phi) = d (phi) / dx  + d (phi) / dy =  0
            Continiutiy in fluids. 

            Thus for potensial:
            u = d (phi) / dx
            v = d (phi) / dy

            For stream function:
            u = d (psi) / dy
            v = - d (psi) / dx

        """

        X, Y = self.mesh.matricies[0], self.mesh.matricies[1] 

        (N_y, N_x) = np.shape(X)

        dx = X[0,1] - X[0,0]   #normally dx should be taken from the spacing matrix. 
        dy = Y[0,0] - Y[1,0]

        if type == "potensial":
            # (u, v) = TwoDcentraldiff_simple(self.solution, dx, dy)
            u, v  = TwoDcentral_diff_velocity(self)

        elif type == "stream":
            (v, u) = TwoDcentral_diff_velocity(self.solution, dx, dy)


        W = np.zeros((N_y, N_x, 3))

        W[:,:,0] = u
        W[:,:,1] = v

        self.velocity = W

    def plot2D(self, type ="potensial", colormap = "copper"):
        """
            plot Z for domain of the X and Y. 
        """

        fig, ax = plt.subplots()

        twin2 = ax.twinx()
        
        x_MAT = self.mesh.matricies[0]
        y_MAT = self.mesh.matricies[1]

        z_min = self.solution.min()
        z_max = self.solution.max()

        fig.set_size_inches(15, 15)
        image = ax.pcolormesh(x_MAT, np.flip(y_MAT), self.solution, vmin=z_min, vmax=z_max, edgecolors="none", cmap=colormap)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="6%", pad="2%")

        if type == "potensial":
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

        elif type == "stream":
            if self.BC["S"] == "D":
                ax.set_xlabel(('$\psi$ = ' + str(self.BCvalues["S"])) , fontsize=20)
            else:
                ax.set_xlabel(('$\partial \psi \partial y$ = '+ str(self.BCvalues["S"])) , fontsize=20)
            if self.BC["W"] == "D":
                ax.set_ylabel(('$\psi$ = '+ str(self.BCvalues["W"])) , fontsize=20)
            else:
                ax.set_ylabel(('$\partial \psi \partial x}$ = ' + str(self.BCvalues["W"])) , fontsize=20)
            if self.BC["N"] == "D":
                ax.set_title(('$\psi$ = ' + str(self.BCvalues["N"])) , fontsize=20)
            else:
                ax.set_title(('$\partial \psi \partial y$ = ' + str(self.BCvalues["N"])) , fontsize=20)
    
    # secondary y-axis label
            
            if self.BC["E"] == "D":
                twin2.set_ylabel(('$\psi$ = ' + str(self.BCvalues["E"])) , fontsize=20)
            else:
                twin2.set_ylabel(('$\partial \psi \partial y}$ = ' +  str(self.BCvalues["E"])) , fontsize=20)
        

        plt.colorbar(image, cax=cax)
        plt.show()

    def streamplot(self, streamcolor = "blue"):
        """
            Plots streamplot for velocity of the solution
        """

        x_MAT = self.mesh.matricies[0]
        y_MAT = self.mesh.matricies[1]

        u = self.velocity[:,:,0]
        v = self.velocity[:,:,1]

        speed = np.sqrt(u*u + v*v)
        lw = 1.2*speed/speed.max()
        lw = 1

        plt.streamplot(x_MAT, y_MAT, -u, v, color=streamcolor, density=0.9, linewidth=lw, cmap='winter')
        plt.show()

    def countour(self):
        pass

    def quiver(self):
        """
            quiver plot of the solution
        """
        x_MAT = self.mesh.matricies[0]
        y_MAT = self.mesh.matricies[1]

        u = self.velocity[:,:,0]
        v = self.velocity[:,:,1]

        fig, ax = plt.subplots()

        fig.set_size_inches(15, 15)
        plt.quiver(x_MAT, y_MAT, -u, v)
        plt.show()
