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


