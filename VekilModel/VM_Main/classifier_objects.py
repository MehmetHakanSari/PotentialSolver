import numpy as np
import matplotlib.pyplot as plt
import os
import glob
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

    def change_name(self, new_name):
        self.name = new_name
        return self


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
                    self.geometry = four_digit_NACA(int(self.name[4]), int(self.name[5]), int(self.name[6:]))
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

    def data_search(self,  datapath = None):
        """
        If data is exist for given airfoil, do not calculate polars. Search by name!
        It needs a datapath to be defined.
        """
        
        airfoil_base_name = self.name[0:4]
        airfoil_directory = os.path.join(os.path.dirname(os.getcwd()), 'data', 'train_data')

        if os.path.exists(airfoil_directory):
            files_and_directories = os.listdir(airfoil_directory)
            Re_folders = [file for file in files_and_directories if file.startswith(airfoil_base_name)]
            for re_folder in Re_folders:
                re_folder_path = os.path.join(airfoil_directory, re_folder)
                Re = int(re_folder_path.split('_')[-1])
                if os.path.isdir(re_folder_path):
                    airfoil_file = os.path.join(re_folder_path, self.name + ".txt")
                    if os.path.isfile(airfoil_file):
                        self.polars[Re] = read_airfoil(airfoil_file, Re)
                                                         
        return self  

        

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


def four_digit_NACA(m, p, t, c=1):
    """
    Generate the four digit NACA airfoil geometry
    """
    #m = max camber
    #p = location of max camber
    #t = thickness
    # print(m, p, t)
    m = m/100
    p = p/10
    t = t/100
    #check the input is valid or not
    if m > 9 or m < 0:
        raise ValueError("m must be between 0 and 9!")
    if p > 9 or p < 0:
        raise ValueError("p must be between 0 and 9!")
    if t > 99 or t < 0:
        raise ValueError("t must be between 0 and 99!")

    #generate the x coordinates
    x = np.linspace(0,c,100)

    Y_c = np.array([m * ((2 * p * x) - np.power(x,2)) / p**2 if x < p else m * ((1 - 2*p) +2*p*x-np.square(x)) / (1 - p)**2 for x in x])
    Y_t = 5*t * (0.2969 * np.sqrt(x) - 0.1260*x - 0.3516 * np.square(x) + 0.2843 * np.power(x,3) - 0.1015 * np.power(x, 4))
    dy_c_dx = np.array([2 * m * (p - x) / p**2 if x < p else 2 * m * (p - x) / (1 - p)**2 for x in x])

    teta = np.arctan(dy_c_dx)
    X_U = np.subtract(x, np.multiply(Y_t,np.sin(teta)))
    Y_U = Y_c + np.multiply(Y_t, np.cos(teta))
    X_L = x + np.multiply(Y_t, np.sin(teta))
    Y_L = Y_c - np.multiply(Y_t, np.cos(teta))
    Y_first = np.flip(Y_U)
    Y = np.round(np.concatenate((Y_first, Y_L)), 6)
    X_first = np.flip(X_U)
    X = np.round(np.concatenate((X_first, X_L)), 6)

    # plt.plot(X_U, Y_U)
    # plt.plot(X_L, Y_L)
    # plt.axis("equal")
    # plt.show()

    X = X.reshape(-1,1)
    Y = Y.reshape(-1,1)

    geometry = np.concatenate((X, Y), axis=1)
    # print(geometry)

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


def readairfoil(foldername, save_path):
    #enter the folder name where the airfoil data is stored
    #folder contains the airfoil data files. This function reads the data form all txt files in the folder
    #returns a dictionary, where the key is the airfoil name and the value is a dictionary of [AOA, Cd, Cl, CDp, Cm].
    #AOD, Cd, Cl, CDp, Cm are dictionary where the key its name and the value is a list of values

    #read  each txt files in the folder, the name of the file is not known in advance
    
    
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    Re = int(foldername.split('_')[-1])
    
    for filename in os.listdir(foldername):
        if filename.endswith(".txt"):
            data_available = False
            #read the file
            f = open(foldername + '/' + filename, 'r')
            #the airfoilname is the file name without the extension
            airfoilname = filename.split('.')[0]
            lines = f.readlines()
        
            #create empty lists to store the AOA, Cd, and Cl, CDp, Cm
            polar_matrix = np.zeros((len(lines[12:]), 8), dtype=float) #xfoil polars starts from line 13

            if len(lines[12:]) > 6:
                data_available = True

            # print(airfoilname)   
            for i, line in enumerate(lines[12:]):
                #split the line into a list
                line = line.split()
                polar_matrix[i,:-1] = np.array(line, dtype=float)
                if float(polar_matrix[i,2] != 0):
                    polar_matrix[i,-1] = float(polar_matrix[i,1]) / float(polar_matrix[i,2]) #ClCd

            polar = {}
            polar[Re] = polar_matrix

            #check if the airfoil exist in savepath directory
            if os.path.exists(save_path + '/' + airfoilname + '.pkl'):
                #if it exists, load the airfoil object
                with open(save_path + '/' + airfoilname + '.pkl', 'rb') as input:
                    airfoil = pickle.load(input)
                
                airfoil.polars.update(polar) 
                airfoil.data_available = data_available

                if type(airfoil.geometry) == None:
                #generate the geometry if it is a NACA airfoil
                    if airfoilname[0:4] == 'NACA':
                        airfoil.NACA_geometry_generator()
                    else:
                        UserWarning("There is no associated name")
                    #update the airfoil object
                    with open(save_path + '/' + airfoilname + '.pkl', 'wb') as output:
                        pickle.dump(airfoil, output, pickle.HIGHEST_PROTOCOL)
                    
                
            else:
                airfoil = Airfoil(name = airfoilname, geometry = None, polars = polar, data_available = data_available)
                if airfoilname[0:4] == 'NACA':
                    airfoil.NACA_geometry_generator()

            #save the airfoil object
            with open(save_path + '/' + airfoilname + '.pkl', 'wb') as output:
                pickle.dump(airfoil, output, pickle.HIGHEST_PROTOCOL)

    return None


def read_airfoil(filename, *args):
    if filename.endswith(".txt"):
        Re = args[0]
        f = open(filename, 'r')
        lines = f.readlines()
    
        #create empty lists to store the AOA, Cd, and Cl, CDp, Cm
        polar_matrix = np.zeros((len(lines[12:]), 8), dtype=float) #xfoil polars starts from line 13
        if len(lines[12:]) > 6:
            data_available = True
        for i, line in enumerate(lines[12:]):
            #split the line into a list
            line = line.split()
            polar_matrix[i,:-1] = np.array(line, dtype=float)
            if float(polar_matrix[i,2] != 0):
                polar_matrix[i,-1] = float(polar_matrix[i,1]) / float(polar_matrix[i,2]) #ClCd

    return polar_matrix



