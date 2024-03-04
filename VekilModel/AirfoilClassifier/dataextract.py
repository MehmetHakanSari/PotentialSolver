import os
import numpy as np
# from classifier_objects import Airfoil
import pickle

#This function reads the AOA, Cd, and Cl from the airfoil data file and returns them as arrays

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
