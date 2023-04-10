import os
import numpy as np

#This function reads the AOA, Cd, and Cl from the airfoil data file and returns them as arrays

def readairfoil(foldername):
    #enter the folder name where the airfoil data is stored
    #folder contains the airfoil data files. This function reads the data form all txt files in the folder
    #returns a dictionary, where the key is the airfoil name and the value is a dictionary of [AOA, Cd, Cl, CDp, Cm].
    #AOD, Cd, Cl, CDp, Cm are dictionary where the key its name and the value is a list of values

    #read txt files in the folder

    #create a list from 1 to 9999 but 100 and its  multiples will not exist in the list. Also each number has 4 digits
    #this is to make sure that the airfoil name is 4 digits
    nacanames = [str(i).zfill(4) for i in range(1,10000) if i%100 != 0]
    airfoildata = {}
    for nacaairfoil in nacanames:
        with open(foldername + '/NACA' + str(nacaairfoil) + '.txt') as f:
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

            for i in range(len(aoa)-1):
                dcl.append((cl[i+1] - cl[i])/(aoa[i+1] - aoa[i]))
            dcl.append(dcl[-1])

            #convert the lists to numpy arrays
            aoa = np.array(aoa)
            cd = np.array(cd)
            cl = np.array(cl)
            cdp = np.array(cdp)
            cm = np.array(cm)
            dcl = np.array(dcl)

            #create a dictionary to store the AOA, Cd, and Cl, CDp, Cm
            airfoil = {'name':("NACA" + str(nacaairfoil)),'AOA': aoa, 'Cd': cd, 'Cl': cl, 'CDp': cdp, 'Cm': cm, 'dCl': dcl}

            #calculate the rate of change of Cl with respect to AOA
            #create an empty list to store the rate of change of Cl with respect to AOA
  

            #create a dictionary to store the airfoil name and its data
            airfoildata['NACA' + str(nacaairfoil)] = airfoil

    return airfoildata