import numpy as np
import os
import pickle
from classifier_objects import Airfoil

#The function will search for thr airfoil for given specific Cl, Cd, Cm for given Re and AoA

def search_foil(Re, alpha_range, data_path, **kwargs):
    """
        Re: Reynolds number, float
        alpha_range: range of alpha, list
        data_path: path of the data, string
        Cl: minimum lift coefficient, float
        Cd: maximum drag coefficient, float
        Cm: moment coefficient, list
        Cp: pressure coefficient, list
        ClCd: lift to drag ratio, list
    """
    #The user should give the specific Cl, Cd for given range of AoA and Re.
    #The function will search for the airfoil that has the closest Cl, Cd, Cm for given Re and AoA

    condition_satisfied_airfoils = []
    # print(os.listdir(data_path))
    #search for folder names including given Re
    Re_folders = [f for f in os.listdir(data_path) if str(Re) and "pkl" in f]
    #import airfoils using pickle in each Re folder
    for Re_folder in Re_folders:
        for file in os.listdir(data_path + Re_folder):
            if file.endswith(".pkl"):
                with open(data_path + Re_folder + "/" + file, "rb") as f:
                    airfoil = pickle.load(f)
                    if airfoil.data_available:
                        if airfoil.polar_search(Re, alpha_range, **kwargs):
                            condition_satisfied_airfoils.append(airfoil)         
        #search for the airfoil that has the closest Cl, Cd, Cm for given Re and AoA

    return condition_satisfied_airfoils

    


