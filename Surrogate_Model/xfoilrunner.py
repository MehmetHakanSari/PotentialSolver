"""Runs an XFOIL analysis for a given airfoil and flow conditions"""
import os
import subprocess
import numpy as np

# %% Inputs

def xfoil_analysis(airfoil_name, alpha_i, alpha_f, alpha_step, savepath, Re, n_iter = 100):

    #files will be stored in a seperate folder. First check if the folder exists if not create it.
 
    if not os.path.isdir(savepath):
        os.makedirs(savepath)

    savepath += "/" + str(airfoil_name) + ".txt"
    #file_name = str(airfoil_name) + ".txt"
 
    if os.path.exists(savepath):
        os.remove(savepath)

    input_file = open("input_file.in", 'w')
    input_file.write("PLOP\n")
    input_file.write("G, F\n")
    input_file.write("\n\n")
    input_file.write("LOAD {0}.dat\n".format(airfoil_name))
    input_file.write(airfoil_name + '\n')
    input_file.write("PANE\n")
    input_file.write("OPER\n")
    input_file.write("Visc {0}\n".format(Re))
    input_file.write("PACC\n")
    input_file.write(savepath + "\n\n")
    input_file.write("ITER {0}\n".format(n_iter))
    input_file.write("ASeq {0} {1} {2}\n".format(alpha_i, alpha_f,
                                                alpha_step))
    input_file.write("\n\n")
    input_file.write("quit\n")
    input_file.close()

    subprocess.call("C:/Users/mtaha/OneDrive/Belgeler/GitHub/Panel_Method/XFoil/xfoil.exe < input_file.in", shell=True)

    file_name = np.loadtxt(savepath, skiprows=12)


alpha_i = -5
alpha_f = 18
alpha_step = 1
n_iter = 100

ReS = [1*1e5, 2*1e5,3*1e5, 4*1e5, 5*1e5, 6*1e5, 7*1e5, 8*1e5, 9*1e5, 1e6]
for Re in ReS:
    Re = int(Re)
    
    folder_name = "NACAxxxx" + "_Re_" + str(Re)
    for i in range(6,7000):
        airfoil_name = "NACA" + str(i).zfill(4)
        if 10 < np.mod(i,100) < 40:
            xfoil_analysis(airfoil_name, alpha_i, alpha_f, alpha_step, folder_name, Re, n_iter)
    print("Re = ", Re, " is done")