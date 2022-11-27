import numpy as np

class Domain:
    """
        Domain is the phsical domain that problem will be solved. (Formal)
        Unformal:
        In my mind I think to do Block type mesh. However For an airfoil which kind of mesh is used I am unaware of that. 
        C mesh for instance is one way to mesh it right? But could we seperate object boundaries while doing C-mesh. 
        First thing to do is 2D simple 1 block mesh. 

        
        Cartasiean:

        Blocks Physical boundaries:

        For 2D: 

        each element is tuple of 3 elements, specifiying lengths of the block. 

        [(x1, y1, z1),        Left Up Corner
         (x2, y1, z1),        Right Up Corner
         (x2, y2, z1),        Right Down Corner
         (x1, y2, z1)]        Left Down Corner

        Mesh Properties:

        [N_x, N_y, N_z] = Total node numbers in x, y, z direction. 

        
    """

    def __init__(self, pyhsical_domain, nodes):
        #write in PEP8 format.
        #list consist of 4 tuple with 3 elements (2D)
        #list consist of 3 float elements        (2D or 3D)

        self.pyhsical_domain = pyhsical_domain  
        self.nodes = nodes                      

    def structured_grid_2D(self):
        x_list = np.linspace(-1, 1, 100)
        y_list = np.linspace(-1, 1, 100)
        


