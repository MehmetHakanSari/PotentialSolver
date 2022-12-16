import numpy as np
from Differentials import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm


class Map():
    
    def __init__(self, area):
        """
            Area is one matrix that you can show object of different kinds.
        """
        self.area = area

    def __call__(self):
        return self.area

    def show(self):
        
        fig, ax = plt.subplots()

        # z = np.ones((self.nodes[1], self.nodes[0])) * 0.5

        image = ax.pcolormesh(self.area, vmin=-1, vmax=1, edgecolors="black",linewidth=0.1)
        plt.show()
