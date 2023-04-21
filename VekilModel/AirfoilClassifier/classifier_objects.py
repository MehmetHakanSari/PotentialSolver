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

    def plot_polars(self):
        pass

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
                    ClCd_statment = np.logical_and(polar[:,5] >= ClCd_range[0], polar[:,5] <= ClCd_range[1]).any()
                if "dClda_range" in kwargs:
                    dClda_range = kwargs["dClda_range"]
                    dClda = (polar[0:-1,1] - polar[1:,1]) / (polar[0:-1,0] - polar[1:,0]) 
                    dClda_statment = np.logical_and(dClda >= dClda_range[0], dClda <= dClda_range[1]).any()
                
                if lift_statment and drag_statment and moment_statment and pressure_statment and ClCd_statment and dClda_statment:
                    return True
                
            else:
                return False

        return None


