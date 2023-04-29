import numpy as np
import matplotlib.pyplot as plt

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

    def plot_geometry(self):
        if self.geometry == None:
            raise FileNotFoundError("geometry of airfoil is not aviable!")
        else:
            pass #write geometry plot funcition

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


