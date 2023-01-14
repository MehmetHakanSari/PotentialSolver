
import sympy as sy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
import numpy as np
from Domain import Mesh
from solvers import PDE_2D_Solver
from Differentials import *
from object import *
from visiual import Map


# block_coordinates = [(0,0),(0,5),(5,5),(5,0)]
# block_coordinates = [(0,0),(0,1.4),(1.4,1),(1,0)]
block_coordinates = [(0,0),(0,3),(3,3),(3,0)]
# print(list(block_coordinates[3]))
node_numbers = [100, 100]

space = Mesh(block_coordinates, node_numbers) #it seems like a one block. Build more sopisticated block type structures
# space.nonuniform_block_mesh_2D(-1.2, -1.4)
# space.nonuniform_mesh_2D(-1.2, -1.4)

#|----------------------------------------------
a, b = space.nonuniform_block_mesh_2D(1, 1)

circle = object()
rectangle = object()
bl = object()


circle.circle(0.2, (0.8,0.6))
property_map = Map(space)
# rectangle.rectang(0.9, 0.2, (0.65,0.4))
property_map.create_object(circle)
# property_map.create_object(rectangle)
rectangle.rectang(0.5, 0.3, (2,2.2))
property_map.create_object(rectangle)
# property_map.show()
#|-----------------------------------------------

# property_map.show()
circle.circle(0.5, (1.4,1.4))
property_map.create_object(circle)
# circle.circle(0.4, (2.4,1.4))
# property_map.create_object(circle)

property_map.show()


BCs = {'W': 'D', 'S': 'N', 'E': 'N', 'N': 'D'}
BCs_values = {'W': 1, 'S': 0.6, 'E': 0.1, 'N': 0.2}


solution = PDE_2D_Solver(space,BCs)
solution.solver(BCs_values, "potensial", property_map, 1.4, 1e-7, itteration_type="nodebynode")

# solution.countour()

# compt_dom = Mesh(block_coordinates, node_numbers)
# compt_dom.uniform_block_mesh_2D()

# compt_dom.plot2D()
# space.plot2D()
# compt_dom.Jacobi(a, b)
# print(compt_dom.Jacobian)

# print(OneDcentraldiff(compt_dom.matricies[0], a))
# print(OneDcentraldiff(compt_dom.matricies[0], b, axis=1))
# print(OneDcentraldiff(compt_dom.matricies[1], a))
# print(OneDcentraldiff(compt_dom.matricies[1], -b, axis=1))
# compt_dom.matricies[1] 


#solution.velocityfield("stream")
#solution.plot2D("stream")
#solution.stream()
#solution.quiver()


solution.velocityfield("potensial")
solution.plot2D("potensial")
solution.streamplot()
solution.quiver()


"""
phii = solution.solution
a = solution.velocity[:,:,0]
cont = solution.continuity
Uinf = np.mean(-solution.velocity[:,:,0])
nu = 1.81e-2
x = np.linspace(0.001,1,(node_numbers[0]))
Re = Uinf*x/nu
delta = 5 * x / np.sqrt(Re)
delta_idx = np.round(delta*len(delta))


bl.boundary_layer(delta_idx)
property_map.create_object(bl)
property_map.show()
solution_BL = PDE_2D_Solver(space,BCs)
solution_BL.solver(BCs_values, "potensial", property_map, 1.4, 1e-5, itteration_type="nodebynode")

solution_BL.velocityfield("potensial")
solution_BL.plot2D("potensial")
solution_BL.streamplot()
solution_BL.quiver()
phii = solution.solution
a = solution.velocity[:,:,0]
cont = solution.continuity

areea = property_map.area

#%matplotlib notebook
"""

"""import sympy as sy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
import numpy as np
from Domain import Mesh
from Domain import ElippticMesh
from solvers import PDE_2D_Solver
from Differentials import *
from object import *
from visiual import Map
import scipy as sp
import scipy.interpolate
from importtools import import_geometry
from solvers import eliptic_PDE_solver
from airfoils import *


def circles(x, y, r):
    th = np.linspace(0, 2 * np.pi, 200)
    xunit = r * np.cos(th) + x
    yunit = r * np.sin(th) + y
    circle_list = np.array([xunit, yunit]).T
    return circle_list

def plotMesh(x, y):
    for i in range(len(x[:, 0])):
        plt.plot(x[i, :], y[i, :], 'k-', lw=2)
    for i in range(len(x[0, :])):
        plt.plot(x[:, i], y[:, i], 'k-', lw=2)

##################### Eliptic GRID #####################
naca0012 = import_geometry('naca0012.txt')
naca0012 = closeshape_interpolation(naca0012, 200)
circleee = circles(3, 3, 3)
naca0012[:,:] = naca0012[:,:] * 1.8

naca0012[:,0], naca0012[:,1] = rotate(naca0012[:,0], naca0012[:,1], -np.pi/15)
naca0012[:,0] += 2.4
naca0012[:,1] += 3.0
node_numbers = [200, 30]
space_trial = ElippticMesh(node_numbers, naca0012, circleee) 
space_trial.create_elipticmesh()
space_trial.plot_mesh()
Vinf = 1
AOA = np.pi / 6 * 0
Out_psi = - space_trial.X[-1, :] * Vinf * np.sin(AOA) + space_trial.Y[-1, :] * Vinf * np.cos(AOA)
BCvalues = {'Out': Out_psi, 'In': 0, 'Cut1': 0, 'Cut2': 0}
stream_naca0012 = eliptic_PDE_solver(space_trial, BCvalues)
stream_naca0012.solver()
stream_naca0012.plot2D()
stream_naca0012.contour()
psi = stream_naca0012.solution

##################### Square GRID #####################
# naca0012 = import_geometry('naca0012.txt')
# naca0012 = closeshape_interpolation(naca0012, 500)

# naca0012[-1, :] = naca0012[0, :] 

# naca0012[:,0] += 0.20
# naca0012[:,1] += 0.1

# naca0012[:,0], naca0012[:,1] = rotate(naca0012[:,0], naca0012[:,1], -np.pi/15)

# plt.figure(1)
# plt.plot(naca0012[:,0], naca0012[:,1])



# block_coordinates = [(0,-0.8),(0,0.8),(3.5,0.8),(3.5,-0.8)]
# node_numbers = [150, 180]
# space_trial = Mesh(block_coordinates, node_numbers) 

# a, b = space_trial.uniform_block_mesh_2D()

# airfoil = object()
# airfoil.airfoil(naca0012, scale=2.3)
# airfoil_map = Map(space_trial)
# airfoil_map.create_object(airfoil)
# airfoil_map.show()

# area_matrix = airfoil_map.area
# area_matrix[79,87] = 0
# area_matrix[83,82] = 0
# area_matrix[60,109] = 0  

# plt.figure(2)
# plt.pcolormesh(area_matrix)
# plt.show()

# BCs = {'W': 'D', 'S': 'N', 'E': 'D', 'N': 'N'}
# BCs_values = {'W': 1, 'S': 0, 'E': 0, 'N': 0}

# solution = PDE_2D_Solver(space_trial,BCs)
# solution.solver(BCs_values, "potensial", airfoil_map, 1.5, 1e-6, itteration_type="nodebynode")

# solution.velocityfield("potensial")
# solution.plot2D("potensial")
# solution.streamplot()
# solution.quiver()"""


"""
import numpy as np
import matplotlib.pyplot as plt
from Differentials import *

x = np.linspace(0,1,200)
y = np.linspace(0,1,200)

X, Y = np.meshgrid(x, y)

Z = 2 * X**2 + Y**2 * X

dZdx = 4 * X + Y**2
dZdy = 2 * Y * X

#plt.figure(1)
#plt.streamplot(X, Y, dZdx, dZdy, color=dZdx, cmap="viridis")

dZdX = OneDcentraldiff(Z, 1/200, 0)
dZdY = OneDcentraldiff(Z, 1/200, 1)

dZdX_u, dZdY_u = TwoDcentraldiff_simple(Z, 1/200, 1/200)

dZdX_CD4, dZdY_CD4 = TwoDcentraldiff_CD4(Z, 1/200, 1/200)

plt.figure(2)
plt.pcolormesh(dZdX_u - dZdX)
plt.colorbar()

plt.figure(3)
plt.pcolormesh(dZdY_u - dZdY)
plt.colorbar()

plt.figure(4)
plt.pcolormesh(dZdy + dZdY)
plt.colorbar()

plt.figure(5)
plt.pcolormesh(dZdy - dZdY_CD4)
plt.colorbar()

plt.figure(6)
plt.pcolormesh(dZdx - dZdX_CD4)
plt.colorbar()

plt.figure(7)
plt.pcolormesh(dZdx - dZdX_u)
plt.colorbar()"""



#plt.figure(2)
#plt.streamplot(X, Y, dZdX, -dZdY, color=dZdX, cmap="viridis")












