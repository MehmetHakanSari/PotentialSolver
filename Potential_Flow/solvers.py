import numpy as np
from Differentials import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from object import object
from visiual import *
from time import perf_counter
from time import sleep
import os 


class eliptic_PDE_solver:


    def __init__(self, mesh, BC_values, BC_type = None):

        self.mesh = mesh
        self.BCvalues = BC_values
        self.BC_type = BC_type

    def solver(self, tolerance, omega):

        N_z = self.mesh.nodes[0]
        N_e = self.mesh.nodes[1]
        alpha = self.mesh.alpha.T[1:-1, 1:-1]
        beta = self.mesh.beta.T[1:-1, 1:-1]
        gamma = self.mesh.gamma.T[1:-1, 1:-1]
        X = self.mesh.X
        Y = self.mesh.Y

        psi = np.zeros((N_e, N_z))  # initial guess for psi
        psi_old = np.zeros((N_e, N_z))  # initial guess for psi

        psi[:,0] = np.ones((N_e)) * self.BCvalues['Cut1']
        psi[:,-1] = np.ones((N_e)) * self.BCvalues['Cut2']
        psi[0,:] = np.ones((1, N_z)) * self.BCvalues['In']
        psi[-1,:] = np.ones((1, N_z)) * self.BCvalues['Out']

        maxiteration = 50000
        message = 1000

        X = X.T
        Y = Y.T

        X_temp = np.append([X[-2, :].copy()], X[0:2, :].copy(), 0) 
        Y_temp = np.append([X[-2, :].copy()], Y[0:2, :].copy(), 0)
        
        alpha_temp, beta_temp, gamma_temp = Solve_Coeff(X_temp , Y_temp)

        psi = psi.T

        for iteration in range(maxiteration):

            psi_old = psi.copy()

            #Periodic BC.
            psi_temp = np.append([psi[-2, :].copy()], psi[0:2, :].copy(), 0)

            psi[0, 1:-1] = SolveEliptic(alpha_temp, beta_temp, gamma_temp, psi_temp)
            psi[-1, :] = psi[0, :].copy()

            psi[:, 0] = psi[0, 1] #Kutta Condiation

            psi[1:-1, 1:-1] = omega * SolveEliptic(alpha, beta, gamma, psi) + (1 - omega) * psi_old
  
            residual_psi = np.max(np.abs(psi - psi_old))

            if residual_psi < tolerance:
                print("Psi is calculated with residual: ", residual_psi, "at itteration: ", iteration)
                break

            #give message with some interval 
            if iteration % message == 0:
                print("Iteration: ", iteration, "Error is: ", residual_psi)

            if iteration == maxiteration - 1:
                print("Max iteration reached. Error is: ", residual_psi)

        self.solution = psi.T


    def velocity_field(self):
        pass
        

    def plot2D(self):
            
        fig = plt.figure(figsize=(10, 8))

        X, Y = self.mesh.X, self.mesh.Y
        #use pcolormesh for better visualization
        plt.pcolormesh(X, Y, self.solution, cmap = 'jet', shading='auto' )
        plt.colorbar()
        plt.show()

    def contour(self):
        plt.figure(figsize=(10, 8), dpi=100)
        cs = plt.contour(self.mesh.X, self.mesh.Y, self.solution, 20, colors='k')
        plt.clabel(cs)
        plt.plot(self.mesh.X[0, :], self.mesh.Y[0, :], 'b-', lw=2) # plot the airfoil
        # plt.xlim((-2., 3.))
        # plt.ylim((-1.5, 1.5))
        plt.show()

    def streamplot(self, streamcolor):
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
        self.continuity = None

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

        (N_y, N_x) = np.shape(mesh.X)

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

       
        x_index = list(range(N_x))
        y_index = list(range(N_y))

        #this spacing updates should be more explainable. 
        if self.BC['W'] == "D":
            phi[:,0] = np.ones((N_y)) * BC_values['W']
            x_index = x_index[1:] 
        if self.BC['S'] == "D":
            phi[-1,:] = np.ones((N_x)) * BC_values['S']
            y_index = y_index[:-1]
        if self.BC['E'] == "D":
            phi[:,-1] = np.ones((N_y)) * BC_values['E']
            x_index = x_index[:-1]
        if self.BC['N'] == "D":
            phi[0,:] = np.ones((N_x)) * BC_values['N']
            y_index = y_index[1:]

        """
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
            """
            
        #Column by Column TDMA
        """
            Starting from first column that is unknown. Thus, 
            first we solve phi[:,1]
            than pass to phi[:,2]

            We will solve again this loop. 
        """
        W = np.zeros(len(y_index), dtype="float")
        E = np.zeros(len(y_index), dtype="float")
        
        message = 1000
        mass_old = 1

        start = perf_counter()    
        for t in range(1, 500001):

            if itteration_type == "column":
                phi = column_TDMA(x_spacing, y_spacing, phi, y_index, x_index, BC_values, property_map.area, N_y, N_x, W, E)
            elif itteration_type == "pointwise":
                phi = pointwise(x_index, y_index, x_spacing, y_spacing, self.BCvalues, phi, phi_old, property_map.area, N_x, N_y, omega, type = variable)

            self.solution = phi.copy()

            #calculate the error
            mass, continuity = self.mass_conservation()
            residual = np.sum(np.abs(phi - phi_old))
            mass_residaul = abs(mass - mass_old) / mass_old

            if (t) % message == 0:
                # mass = self.mass_conservation()
                print("Residual: ", residual,"Mass Residual:", mass_residaul ,"Mass:", mass, " Continuity: ", np.sum(continuity) ," at ", t, "th iteration", "  Time: ", perf_counter() - start)        
                
                #Store data for each message iteration in a seperate directory. First check whether a result directory exists or not.
                #Open a new directory for one solution set. Inside that directory, create a new directory for each message iteration.
            
                if not os.path.exists("results"):
                    os.mkdir("results")
                if not os.path.exists("results/" + self.name):
                    os.mkdir("results/" + self.name)
                if not os.path.exists("results/" + self.name + "/" + str(t)):
                    os.mkdir("results/" + self.name + "/" + str(t))
                
                #store the solution in that directory
                np.save("results/" + self.name + "/" + str(t) + "/solution", phi)
                
                self.velocityfield("potensial")
                
                #Also plot the solution and store it in the same directory
                self.plot2D(phi, variable, "results/" + self.name + "/" + str(t) + "/solution.png")
                self.streamplot("results/" + self.name + "/" + str(t) + "/streamplot.png")
                
                # sleep(1)

            if residual < toll or abs(mass_residaul) < toll:
                print("Residual: ", residual,"Mass Residual:", mass_residaul, "Mass:", mass)
                print("Solution converged at ", t, "th iteration")
                break

            phi_old = phi.copy()
            mass_old = mass.copy()

        plt.pcolormesh(mesh.X, mesh.Y, continuity)
        self.continuity = continuity.copy()
            

    def mass_conservation(self):
        """
            This function is for checking the mass conservation. 
            It is not necessary for the solver. 
        """
        self.velocityfield("potensial")
        u = -self.velocity[:,:,0]
        v = self.velocity[:,:,1]

        X, Y = self.mesh.X, self.mesh.Y 

        dx = float(X[0,1] - X[0,0])   #normally dx should be taken from the spacing matrix. 
        dy = float(Y[0,0] - Y[1,0])

        dudx = OneDcentraldiff(-u, dx, 0)
        dvdy = OneDcentraldiff(v, dy, 1)

        continuity = dudx + dvdy

        mass = np.sum(u[:,0]) - np.sum(u[:,-1]) + np.sum(u[0,:]) - np.sum(u[-1,:]) + np.sum(v[:,0]) - np.sum(v[:,-1]) + np.sum(v[0,:]) - np.sum(v[-1,:])

        return mass, continuity


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

        X, Y = self.mesh.X, self.mesh.Y

        (N_y, N_x) = np.shape(X)

        dx = X[0,1] - X[0,0]   #normally dx should be taken from the spacing matrix. 
        dy = Y[0,0] - Y[1,0]

        if type == "potensial":
            # (u, v) = TwoDcentraldiff_simple(self.solution, dx, dy)
            #u, v  = TwoDcentral_diff_velocity_CD2(self)
            u, v  = TwoDcentral_diff_velocity_CD4(self)

        elif type == "stream":
            (v, u) = TwoDcentral_diff_velocity_CD2(self.solution, dx, dy)


        W = np.zeros((N_y, N_x, 3))

        W[:,:,0] = u
        W[:,:,1] = v

        self.velocity = W

    def plot2D(self, type ="potensial", savepath = "None", colormap = "copper"):
        """
            plot Z for domain of the X and Y. 
        """

        fig, ax = plt.subplots()

        twin2 = ax.twinx()
        
        x_MAT = self.mesh.X
        y_MAT = self.mesh.Y

        z_min = self.solution.min()
        z_max = self.solution.max()

        fig.set_size_inches(15, 15)
        image = ax.pcolormesh(x_MAT, y_MAT, self.solution, vmin=z_min, vmax=z_max, edgecolors="none", cmap=colormap)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="6%", pad="2%")

        if type == "potensial":
            if self.BC["N"] == "D":
                ax.set_xlabel(('$\phi$ = ' + str(self.BCvalues["N"])) , fontsize=20)
            else:
                ax.set_xlabel(('${\partial \phi}{\partial y}$ = '+ str(self.BCvalues["N"])) , fontsize=20)
            if self.BC["W"] == "D":
                ax.set_ylabel(('$\phi$ = '+ str(self.BCvalues["W"])) , fontsize=20)
            else:
                ax.set_ylabel(('${\partial \phi}{\partial x}$ = ' + str(self.BCvalues["W"])) , fontsize=20)
            if self.BC["S"] == "D":
                ax.set_title(('$\phi$ = ' + str(self.BCvalues["S"])) , fontsize=20)
            else:
                ax.set_title(('$\partial \phi \partial y$ = ' + str(self.BCvalues["S"])) , fontsize=20)
    
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
        
        #save the figure for given savepath
        if savepath != None:
            plt.savefig(savepath, dpi=300)
        else:
            plt.show()

    def streamplot(self, savepath = "None", streamcolor = "blue"):
        """
            Plots streamplot for velocity of the solution
        """

        u = self.velocity[:,:,0]
        v = self.velocity[:,:,1]

        speed = np.sqrt(u*u + v*v)
        lw = 10.2*speed/speed.max()
        # lw = 1

        fig, (axs1, axs2) = plt.subplots(1,2)
        fig.set_size_inches(12, 5)
        axs1.streamplot(self.mesh.X, self.mesh.Y, -u, v, color=streamcolor, density=0.9, linewidth=lw, cmap='winter')
        axs2.streamplot(self.mesh.X, self.mesh.Y, -u, v, color=streamcolor, density=0.9, linewidth=1, cmap='winter')
        
        if savepath != None:
            plt.savefig(savepath, dpi=300)
        else:
            plt.show()

    def countour(self):
        pass

    def quiver(self):
        """
            quiver plot of the solution
        """

        u = self.velocity[:,:,0]
        v = self.velocity[:,:,1]

        fig, ax = plt.subplots()

        fig.set_size_inches(15, 15)
        plt.quiver(self.mesh.X, self.mesh.Y, -u, v)
        plt.show()
