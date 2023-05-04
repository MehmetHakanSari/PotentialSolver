import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from dataplotter import geometry_plotter

class Airfoil:
    """
    The object, is airfoil class which includes name, geometry, polars, airfoil specific coordiantes. 
    It can call VM-AF to fill up the missing polars
    It can call plotter functions to plot polars and airfoil. 
    It might have self meshing tools and flow graph etc. 
    """
    def __init__(self, name = None, geometry = None, polars = None, data_available = False):
        """
        name: string

        geometry: ndarray (list, tuple) in format [N x 2] or [N x 3]

        polars: dictionary, key is double representing Reynolds Number, value is [N x 7] list. 

        """
        self.name = name
        self.geometry = geometry
        self.polars = polars
        self.data_available = data_available
        self.position_info = {}

    def __call__():
        """
        specifications, important geometrical position vectors, names etc.
        """
        pass

    def plot_geometry(self,  *args):
        """
        Plot geometry of airfoil. 
        """
 
        geometry_plotter(self.geometry, *args)
                

    def plot_polars(self, Re, *args):
        """
        Plot polars for given Re. If Re is not given, plot all polars. 
        """
        if self.polars == None:
            raise FileNotFoundError("polars of airfoil is not aviable!")
        elif Re not in self.polars.keys():
            raise KeyError("Reynolds number is not aviable!")
        else:
            if "lift" in args:
                plt.figure()
                plt.plot(self.polars[Re][:,0], self.polars[Re][:,1])
                plt.xlabel("AoA")
                plt.ylabel("Cl")
                plt.title("Lift Coefficient vs AoA")
                plt.grid()
                plt.show()
            if "drag" in args:
                plt.figure()
                plt.plot(self.polars[Re][:,0], self.polars[Re][:,2])
                plt.xlabel("AoA")
                plt.ylabel("Cd")
                plt.title("Drag Coefficient vs AoA")
                plt.grid()
                plt.show()
            if "lift_drag" in args:
                plt.figure()
                plt.plot(self.polars[Re][:,1], self.polars[Re][:,2])
                plt.xlabel("Cl")
                plt.ylabel("Cd")
                plt.title("Lift vs Drag Coefficient")
                plt.grid()
                plt.show()
            if "moment" in args:
                plt.figure()
                plt.plot(self.polars[Re][:,0], self.polars[Re][:,3])
                plt.xlabel("AoA")
                plt.ylabel("Cm")
                plt.title("Moment Coefficient vs AoA")
                plt.grid()
                plt.show()
            if "liftdrag" in args:
                plt.figure()
                plt.plot(self.polars[Re][:,0], self.polars[Re][:,-1])
                plt.xlabel("AoA")
                plt.ylabel("Cl/Cd")
                plt.title("Lift/Drag Coefficient vs AoA")
                plt.grid()
                plt.show()

    def polar_update(self, Re, polars = None):
        """
        Update the polar for given Re. If Re is not given, update all polars. 
        """
        self.polars[Re] = polars

    def set_position_info(self):
        """
        Set the positon information according to Re = 0 lift and momment coefficient. 
        """

        if self.polars[0] == None:
            raise FileNotFoundError("polars of airfoil is not aviable for Re = 0!")
        
        polar_matrix = self.polars[0]
        dClda = np.gradient(polar_matrix[:,1], polar_matrix[:,0])
        dCm = np.gradient(polar_matrix[:,3], polar_matrix[:,0])

        dCm_mean = np.mean(dCm)
        dClda_mean = np.mean(dClda)

        x_aerodynamic_center = - dCm_mean / dClda_mean + 0.25

        self.position_info["aerodyanmic-center"] = x_aerodynamic_center


    def NACA_geometry_generator(self):
        """
        Generate the NACA geometry from the airfoil name
        """
        if self.name == None:
            raise FileNotFoundError("Airfoil name is not aviable!")
        else:
            #check the name is NACA or not
            if self.name[0:4] == "NACA":
                #check the name is 4 digit or 5 digit
                if len(self.name) == 8:
                    #4 digit
                    self.geometry = four_digit_NACA(int(self.name[5]), int(self.name[6]), int(self.name[7:]))
                elif len(self.name) == 9:
                    #5 digit
                    pass
                else:
                    raise ValueError("Airfoil name is not valid!")
                      

    def call_SM(self, Re):
        pass

    def call_xfoil(self, Re):
        pass

    def calculate_polars(self, Re):
        """
        First check wheter data is exist or not. Then search for it. If it is not exist call VM.
        """
        pass

    def data_search(self):
        """
        If data is exist for given airfoil, do not calculate polars. Search by name!
        It needs a datapath to be defined.
        """
        pass

    def polar_search(self, Re, alpha_range, **kwargs):
        """
        Search for the closest polar for given Re.
        """
        
        if Re not in self.polars:
            return None
        else:
            polar = self.polars[Re]
            #check the AoA range from given kwargs and applied it to polar
            #before check wheter airfoil has datas for given alpha range

            if not(polar[:,0].min() >= alpha_range[0] or polar[:,0].max() <= alpha_range[1]):

                polar = polar[np.logical_and(polar[:,0] >= alpha_range[0], polar[:,0] <= alpha_range[1])]
                lift_statment = True
                drag_statment = True
                moment_statment = True
                pressure_statment = True
                ClCd_statment = True
                dClda_statment = True
                if "Cl" in kwargs:
                    Cl = kwargs["Cl"]
                    lift_statment = polar[:,1].min() >= Cl
                if "Cd" in kwargs:
                    Cd = kwargs["Cd"]
                    drag_statment = polar[:,2].max() <= Cd
                if "Cm_range" in kwargs:
                    Cm_range = kwargs["Cm_range"]
                    moment_statment = np.logical_and(polar[:,3] >= Cm_range[0], polar[:,3] <= Cm_range[1]).any()
                if "Cp_range" in kwargs:
                    Cp_range = kwargs["Cp_range"]
                    pressure_statment = np.logical_and(polar[:,4] >= Cp_range[0], polar[:,4] <= Cp_range[1]).any()
                if "ClCd_range" in kwargs:
                    ClCd_range = kwargs["ClCd_range"]
                    ClCd_statment = np.logical_and(polar[:,-1] >= ClCd_range[0], polar[:,-1] <= ClCd_range[1]).any()
                if "dClda_range" in kwargs:
                    dClda_range = kwargs["dClda_range"]
                    dClda = (polar[0:-1,1] - polar[1:,1]) / (polar[0:-1,0] - polar[1:,0]) 
                    dClda_statment = np.logical_and(dClda >= dClda_range[0], dClda <= dClda_range[1]).any()
                
                if lift_statment and drag_statment and moment_statment and pressure_statment and ClCd_statment and dClda_statment:
                    return True
                
            else:
                return False

        return None


def four_digit_NACA(m, p, t):
    """
    Generate the four digit NACA airfoil geometry
    """
    #m = max camber
    #p = location of max camber
    #t = thickness

    #check the input is valid or not
    if m > 9 or m < 0:
        raise ValueError("m must be between 0 and 9!")
    if p > 9 or p < 0:
        raise ValueError("p must be between 0 and 9!")
    if t > 99 or t < 0:
        raise ValueError("t must be between 0 and 99!")

    #generate the x coordinates
    x = np.linspace(0,1,100)
    #generate the y coordinates
    y = np.zeros(100)

    #calculate the camber line
    for i in range(100):
        if x[i] <= p / 10:
            y[i] = m / (p / 10)**2 * (2 * p / 10 * x[i] - x[i]**2)
        else:
            y[i] = m / (1 - p / 10)**2 * ((1 - 2 * p / 10) + 2 * p / 10 * x[i] - x[i]**2)

    #calculate the thickness
    yt = 5 * t / 100 * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)

    #calculate the upper and lower surface
    xu = x - yt * np.sin(np.arctan(np.gradient(y,x)))
    xl = x + yt * np.sin(np.arctan(np.gradient(y,x)))
    yu = y + yt * np.cos(np.arctan(np.gradient(y,x)))
    yl = y - yt * np.cos(np.arctan(np.gradient(y,x)))

    #generate the geometry
    geometry = np.zeros((200,2))
    geometry[0:100,0] = xu
    geometry[0:100,1] = yu
    geometry[100:,0] = xl[::-1]
    geometry[100:,1] = yl[::-1]

    return geometry

def five_digit_NACA(m, p, t, c):
    """
    Generate the five digit NACA airfoil geometry
    """
    #m = max camber
    #p = location of max camber
    #t = thickness
    #c = chord length

    #check the input is valid or not
    if m > 9 or m < 0:
        raise ValueError("m must be between 0 and 9!")
    if p > 9 or p < 0:
        raise ValueError("p must be between 0 and 9!")
    if t > 99 or t < 0:
        raise ValueError("t must be between 0 and 99!")
    if c < 0:
        raise ValueError("c must be positive!")

    #generate the x coordinates
    x = np.linspace(0,1,100)
    #generate the y coordinates
    y = np.zeros(100)

    #calculate the camber line
    for i in range(100):
        if x[i] <= p / 10:
            y[i] = m / (p / 10)**2 * (2 * p / 10 * x[i] - x[i]**2)
        else:
            y[i] = m / (1 - p / 10)**2 * ((1 - 2 * p / 10) + 2 * p / 10 * x[i] - x[i]**2)

    #calculate the thickness
    yt = 5 * t / 100 * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)

    #calculate the upper and lower surface
    xu = x - yt * np.sin(np.arctan(np.gradient(y,x)))
    xl = x + yt * np.sin(np.arctan(np.gradient(y,x)))
    yu = y + yt * np.cos(np.arctan(np.gradient(y,x)))
    yl = y - yt * np.cos(np.arctan(np.gradient(y,x)))

    #generate the geometry
    geometry = np.zeros((200,2))
    geometry[0:100,0] = xu
    geometry[0:100,1] = yu
    geometry[100:,0] = xl[::-1]
    geometry[100:,1] = yl[::-1]

    #scale the geometry
    geometry[:,0] = geometry[:,0] * c

    return geometry


    


