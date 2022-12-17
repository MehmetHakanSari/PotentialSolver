import numpy as np
from Differentials import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm


class object:

    def __init__(self):
        self.type = None

    def circle(self, radius, center):
        """
            radius: float. circles radius
            center: tuple, containing (x, y)
        """

        self.type = "circle"
        self.radius = radius
        self.center = center