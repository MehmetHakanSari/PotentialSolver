import matplotlib.pyplot as plt
import numpy as np

#This function plots the airfoil data

def plot_polars(airfoil, Re):
    """
    airfoil: an airfoil object

    Re: Reynolds number, list of Reynolds number, or 'all'
    """

    if Re == 'all':
        Re_list = airfoil.polars.keys()

    if type(Re) == int:
        Re_list = [Re]

    if Re not in airfoil.polars.keys():
        raise KeyError("Reynolds number is not aviable!")

    for Re in Re_list:
        polar = airfoil.polars[Re]
        plt.figure(1, figsize=(6, 10))
        plt.subplot(2, 2, 1)
        #plot AOA vs Cl as a square dots
        plt.plot(polar[:,0], polar[:,1], 's', label = 'Cl-' + airfoil.name + '-' + str(Re))
        #plot AOA vs Cm as a triangle dots
        plt.plot(polar[:,0], polar[:,4], '^', label = 'Cm-' + airfoil.name + '-' + str(Re))
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.xlabel('alpha, degress')
        #plt.title() #change the title to the airfoil name
        plt.legend()
        plt.show()

        plt.subplot(2, 2, 2)
        #plot Cd vs AOA as a circle dots
        plt.plot(polar[:,0], polar[:,2], '-', label = 'Cd-' + airfoil.name + '-' + str(Re))
        #also plot x and y axis
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.xlabel('alpha, degress')
        plt.legend()
        plt.show()

        plt.subplot(2, 2, 3)
        #plot Cd vs Cl as a circle dots
        plt.plot(polar[:,1], polar[:,2], '-', label = 'Cd-' + airfoil.name + '-' + str(Re))
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.xlabel('Cd')
        plt.ylabel('Cl')
        plt.legend()
        plt.show()

        plt.subplot(2, 2, 4)
        plt.plot(polar[:,0], polar[:,-1], '-', label = 'Cl/Cd-' + airfoil.name + '-' + str(Re))
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.xlabel('alpha, degress')
        plt.ylabel('Cl')
        plt.legend()
        plt.show()


def compare_airfoils(airfoil1, airfoil2):
    """
    airfoil1: an airfoil object
    airfoil2: an airfoil object
    """
    #plot values in same graph

    plt.figure(1, figsize=(6, 10))
    #make figure size larger
    #plt.rcParams["figure.figsize"] = (15,15)
    plt.subplot(2, 1, 1)
    #plot AOA vs Cl as a square dots
    plt.plot(airfoil1['AOA'], airfoil1['Cl'], 'bs', label = 'Cl-' + airfoil1['name'])
    plt.plot(airfoil2['AOA'], airfoil2['Cl'], 'rs', fillstyle = 'none', label = 'Cl-' + airfoil2['name'])

    plt.plot(airfoil1['AOA'], airfoil1['Cm'], 'b^', label = 'Cm-' + airfoil1['name'])
    plt.plot(airfoil2['AOA'], airfoil2['Cm'], 'r^', fillstyle = 'none', label = 'Cm-' + airfoil2['name'])

    #also plot x and y axis
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')

    plt.xlabel('alpha, degress')

    plt.legend()
    plt.show()

    plt.subplot(2, 1, 2)
    #plot dCl vs AOA as a circle dots
    plt.plot(airfoil1['AOA'], airfoil1['dCl'], 'o', label = 'dCl-' + airfoil1['name'])
    #plot second airfoil with red color circle dots unfilled
    plt.plot(airfoil2['AOA'], airfoil2['dCl'], 'ro', fillstyle = 'none', label = 'dCl-' + airfoil2['name'])
    #also plot x and y axis
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    plt.xlabel('alpha, degress')
    plt.legend()
    plt.show()

    plt.subplot(3, 1, 3)
    #plot Cd vs Cl as a circle dots
    plt.plot(airfoil1['Cd'], airfoil1['Cl'], '-', label = 'Cd-' + airfoil1['name'])
    plt.plot(airfoil2['Cd'], airfoil2['Cl'], 'r-', label = 'Cd-' + airfoil2['name'])
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    plt.xlabel('Cd')
    plt.ylabel('Cl')
    plt.legend()
    plt.show()
