import os
import numpy as np

#This function reads the AOA, Cd, and Cl from the airfoil data file and returns them as arrays

def readairfoil(foldername):
    #enter the folder name where the airfoil data is stored
    #folder contains the airfoil data files. This function reads the data form all txt files in the folder
    #returns a dictionary, where the key is the airfoil name and the value is a dictionary of [AOA, Cd, Cl, CDp, Cm].
    #AOD, Cd, Cl, CDp, Cm are dictionary where the key its name and the value is a list of values

    #read  each txt files in the folder, the name of the file is not known in advance
    airfoildata = {}
    
    for filename in os.listdir(foldername):
        if filename.endswith(".txt"):
            data_available = False
            #read the file
            f = open(foldername + '/' + filename, 'r')

            #the airfoilname is the file name without the extension
            airfoilname = filename.split('.')[0]

            #read the lines in the file
            lines = f.readlines()
            #the values starts from the 11th line 
            #the first 10 lines are the header

            #create empty lists to store the AOA, Cd, and Cl, CDp, Cm
            aoa = []
            cd = []
            cl = []
            cdp = []
            cm = []
            dcl = []   
             

            for line in lines[12:]:
                #split the line into a list
                line = line.split()
                #append the values to the lists
                aoa.append(float(line[0]))
                cl.append(float(line[1]))
                cd.append(float(line[2]))
                cdp.append(float(line[3]))
                cm.append(float(line[4]))

            if len(aoa) > 12:
                data_available = True
                for i in range(len(aoa)-1):
                    dcl.append((cl[i+1] - cl[i])/(aoa[i+1] - aoa[i]))
                dcl.append(dcl[-1])

            #create a dictionary to store the AOA, Cd, and Cl, CDp, Cm
            airfoil = {'name':(str(airfoilname)),'AOA': np.array(aoa), 'Cd': np.array(cd), 'Cl': np.array(cl), 'CDp': np.array(cdp), 'Cm': np.array(cm), 'dCl': np.array(dcl), 'CmCd': np.array(cl) / np.array(cd), 'data_available': data_available}

            #calculate the rate of change of Cl with respect to AOA
            #create an empty list to store the rate of change of Cl with respect to AOA
  

            #create a dictionary to store the airfoil name and its data
            airfoildata[str(airfoilname)] = airfoil

    return airfoildata