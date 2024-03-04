import matplotlib.pyplot as plt
import numpy as np
import math
import warnings

#This function plots the airfoil data

def plot_polars(airfoil, Re, **kwargs):
    """
    airfoil: an airfoil object

    Re: Reynolds number, list of Reynolds number, or 'all'
    """

    if Re == 'all':
        Re_list = airfoil.polars.keys()   
    elif Re not in airfoil.polars.keys():
        raise KeyError("Reynolds number is not aviable!")

    if type(Re) == int:
        Re_list = [Re]

    
    LiWi = kwargs.get('LiWi', 1)
    MaSi = kwargs.get('MaSi', 5)

    fig, axs = plt.subplots(2, 2, figsize=(10, 6))
    for Re in Re_list:
        if Re != 0:
            polar = airfoil.polars[Re]
        else:
            continue

        # Plot AOA vs Cl as square dots
        axs[0, 0].plot(polar[:, 0], polar[:, 1], '-s', label= 'Re: ' + str(Re), linewidth=LiWi, markersize=MaSi)
        # Plot AOA vs Cm as triangle dots
        axs[0, 0].plot(polar[:, 0], polar[:, 4], '-^', label= 'Re: ' + str(Re), linewidth=LiWi, markersize=MaSi)
        axs[0, 0].set_title('Cl and Cm')
    
        # Plot Cd vs AOA as circle dots
        axs[0, 1].plot(polar[:, 0], polar[:, 2], '-', label= 'Re: ' + str(Re), linewidth=LiWi, markersize=MaSi)
        axs[0, 1].set_title('Cd')
    
        # Plot Cd vs AOA as circle dots
        axs[1, 0].plot(polar[:,1], polar[:,2], '-', label = 'Re: ' + str(Re),  linewidth=LiWi, markersize=MaSi)
        axs[1, 0].set_title('Cd v Cl')
    
        # Plot Cd vs AOA as circle dots
        axs[1, 1].plot(polar[:,0], polar[:,-1], '-', label = 'Re: ' + str(Re), linewidth=LiWi, markersize=MaSi)
        axs[1, 1].set_title('Cl/Cd')

        
    for i, ax_row in enumerate(axs):
        for j, ax in enumerate(ax_row):
            if i == 0 and j == 0:
                ax.axhline(y=0, color='k')
                ax.axvline(x=0, color='k')
                ax.set_xlabel('alpha, degrees')
                ax.legend()
            if i == 0 and j == 1:
                ax.axhline(y=0, color='k')
                ax.axvline(x=0, color='k')
                ax.set_xlabel('Cd')
                ax.legend()
            if i == 1 and j == 0: 
                ax.axhline(y=0, color='k')
                ax.axvline(x=0, color='k')
                ax.set_xlabel('Cl')
                ax.set_ylabel('Cd')
                ax.legend()
            if i == 1 and j == 1: 
                ax.axhline(y=0, color='k')
                ax.axvline(x=0, color='k')
                ax.set_xlabel('alpha, degrees')
                ax.set_ylabel('Cl/Cd')
                ax.legend()

    # Adjust layout and display the plot
    plt.tight_layout()
    plt.show()

    # for Re in Re_list:
    #     polar = airfoil.polars[Re]
    #     plt.figure(1, figsize=(6, 10))
    #     plt.subplot(2, 2, 1)
    #     #plot AOA vs Cl as a square dots
    #     ax.plot(polar[:,0], polar[:,1], 's', label = 'Cl-' + airfoil.name + '-' + str(Re))
    #     #plot AOA vs Cm as a triangle dots
    #     ax.plot(polar[:,0], polar[:,4], '^', label = 'Cm-' + airfoil.name + '-' + str(Re))
    #     plt.axhline(y=0, color='k')
    #     plt.axvline(x=0, color='k')
    #     plt.xlabel('alpha, degress')
    #     #plt.title() #change the title to the airfoil name
    #     plt.legend()
    #     plt.show()

    #     plt.subplot(2, 2, 2)
    #     #plot Cd vs AOA as a circle dots
    #     plt.plot(polar[:,0], polar[:,2], '-', label = 'Cd-' + airfoil.name + '-' + str(Re))
    #     #also plot x and y axis
    #     plt.axhline(y=0, color='k')
    #     plt.axvline(x=0, color='k')
    #     plt.xlabel('alpha, degress')
    #     plt.legend()
    #     plt.show()

    #     plt.subplot(2, 2, 3)
    #     #plot Cd vs Cl as a circle dots
    #     plt.plot(polar[:,1], polar[:,2], '-', label = 'Cd-' + airfoil.name + '-' + str(Re))
    #     plt.axhline(y=0, color='k')
    #     plt.axvline(x=0, color='k')
    #     plt.xlabel('Cd')
    #     plt.ylabel('Cl')
    #     plt.legend()
    #     plt.show()

    #     plt.subplot(2, 2, 4)
    #     plt.plot(polar[:,0], polar[:,-1], '-', label = 'Cl/Cd-' + airfoil.name + '-' + str(Re))
    #     plt.axhline(y=0, color='k')
    #     plt.axvline(x=0, color='k')
    #     plt.xlabel('alpha, degress')
    #     plt.ylabel('Cl')
    #     plt.legend()
    #     plt.show()


def compare_airfoils(airfoil1, airfoil2, Re, **kwargs):
    """
    airfoil1: an airfoil object
    airfoil2: an airfoil object
    """
    Re1, Re2 = (Re, Re)
    if Re not in airfoil1.polars.keys():
        warnings.warn("Reynolds number is not aviable in airfoil 1. Looking for closest Re")
        Re1_L = airfoil1.polars.keys()
        Re1 = min(Re1_L, key=lambda x: abs(x - Re))
    if Re not in airfoil2.polars.keys():
        warnings.warn("Reynolds number is not aviable in airfoil 2. Looking for closest Re")
        Re2_L = airfoil1.polars.keys()
        Re2 = min(Re2_L, key=lambda x: abs(x - Re))

    if type(Re) != int:
        raise KeyError("Only one Reynolds number is avaibable")
         
    LiWi = kwargs.get('LiWi', 1)
    MaSi = kwargs.get('MaSi', 5)

    fig, axs = plt.subplots(2, 2, figsize=(10, 6))
    if Re != 0:
        polar1 = airfoil1.polars[Re1]
        polar2 = airfoil2.polars[Re2]

    # Plot AOA vs Cl as square dots
    axs[0, 0].plot(polar1[:, 0], polar1[:, 1], '-s', label= 'Re: ' + str(Re1), linewidth=LiWi, markersize=MaSi)
    axs[0, 0].plot(polar2[:, 0], polar2[:, 1], '--s', label= 'Re: ' + str(Re2), linewidth=LiWi, markersize=MaSi)
    # Plot AOA vs Cm as triangle dots
    axs[0, 0].plot(polar1[:, 0], polar1[:, 4], '-^', label= 'Re: ' + str(Re1), linewidth=LiWi, markersize=MaSi)
    axs[0, 0].plot(polar2[:, 0], polar2[:, 4], '--^', label= 'Re: ' + str(Re2), linewidth=LiWi, markersize=MaSi)
    axs[0, 0].set_title('Cl and Cm')

    # Plot Cd vs AOA as circle dots
    axs[0, 1].plot(polar1[:, 0], polar1[:, 2], '-', label= 'Re: ' + str(Re1), linewidth=LiWi, markersize=MaSi)
    axs[0, 1].plot(polar2[:, 0], polar2[:, 2], '--', label= 'Re: ' + str(Re2), linewidth=LiWi, markersize=MaSi)
    axs[0, 1].set_title('Cd')

    # Plot Cd vs AOA as circle dots
    axs[1, 0].plot(polar1[:,1], polar1[:,2], '-', label = 'Re: ' + str(Re1),  linewidth=LiWi, markersize=MaSi)
    axs[1, 0].plot(polar2[:,1], polar2[:,2], '--', label = 'Re: ' + str(Re2),  linewidth=LiWi, markersize=MaSi)
    axs[1, 0].set_title('Cd v Cl')

    # Plot Cd vs AOA as circle dots
    axs[1, 1].plot(polar1[:,0], polar1[:,-1], '-', label = 'Re: ' + str(Re1), linewidth=LiWi, markersize=MaSi)
    axs[1, 1].plot(polar2[:,0], polar2[:,-1], '--', label = 'Re: ' + str(Re2), linewidth=LiWi, markersize=MaSi)
    axs[1, 1].set_title('Cl/Cd')

        
    for i, ax_row in enumerate(axs):
        for j, ax in enumerate(ax_row):
            if i == 0 and j == 0:
                ax.axhline(y=0, color='k')
                ax.axvline(x=0, color='k')
                ax.set_xlabel('alpha, degrees')
                ax.legend()
            if i == 0 and j == 1:
                ax.axhline(y=0, color='k')
                ax.axvline(x=0, color='k')
                ax.set_xlabel('Cd')
                ax.legend()
            if i == 1 and j == 0: 
                ax.axhline(y=0, color='k')
                ax.axvline(x=0, color='k')
                ax.set_xlabel('Cl')
                ax.set_ylabel('Cd')
                ax.legend()
            if i == 1 and j == 1: 
                ax.axhline(y=0, color='k')
                ax.axvline(x=0, color='k')
                ax.set_xlabel('alpha, degrees')
                ax.set_ylabel('Cl/Cd')
                ax.legend()

    # Adjust layout and display the plot
    plt.tight_layout()
    plt.show()


def geometry_plotter(geometry, *args):
    """
    geometry: list, [x, y]
    
    *args:

    alpha_range = [alpha_min, alpha_max]
    """

    # print(*args)
    # print(args)
    # print(geometry[:,0])
    # print(geometry[:,1])

    plt.figure(1, facecolor='#212121', figsize=(15, 5))
    cmap = plt.get_cmap('viridis')
    ax = plt.axes()
    ax.set_facecolor("#424242")
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')

    if len(args) == 0:
        plt.figure(1)
        ax.plot(geometry[:,0], geometry[:,1], color='#00AAAF')
        # print(min(geometry[:,0]) - 0.2 * abs(max(geometry[:,0])), max(geometry[:,0]) + 0.2 * max(geometry[:,0]))
        # print(min(geometry[:,1]) - 0.2 * abs(min(geometry[:,1])), max(geometry[:,1]) + 0.2 * max(geometry[:,1]))
        ax.set_xlim(min(geometry[:,0]) - 0.1 * abs(max(geometry[:,0])), max(geometry[:,0]) + 0.1 * max(geometry[:,0]))
        ax.set_ylim(min(geometry[:,1]) - 0.4 * abs(min(geometry[:,1])), max(geometry[:,1]) + 0.4 * max(geometry[:,1]))
        plt.show()
    else:
        alpha_range = args[0]  
        #plot geometry at different angle of attack
        for alpha in alpha_range:
            #rotate the geometry
            geometry_rotated = rotate(geometry, alpha)
            normalized_angle = (alpha - (-10)) / ((20) - (-10))
            color = cmap(1 - normalized_angle)
            plt.figure(1)
            # plt.plot(geometry_rotated[:,0], geometry_rotated[:,1], '--', color='#00AFFF', linewidth=2)
            plt.plot(geometry_rotated[:,0], geometry_rotated[:,1], '--', color=color, linewidth=2)

        plt.figure(1)
        # plt.axis("equal")
        plt.xlim(-1.1, 0.1)
        plt.ylim(-0.3, 0.3)
        plt.show()

def rotate(geometry, alpha):
    """
        x: x coordinates of airfoil
        y: y coordinates of airfoil
        angle: angle to rotate airfoil
    """
    #The rotation of the airfoil needed to be done rom its trailing edge.
    #Carry translation to the origin

    #The rotation of the airfoil needed to be done rom its trailing edge.
    # 
    #Carry translation to the origin
    geometry[:,0] = geometry[:,0] - geometry[0,0]
    geometry[:,1] = geometry[:,1] - geometry[0,1]
    
    x_rotated = []
    y_rotated = []
    x = geometry[:, 0]
    y = geometry[:, 1]

    angle = -math.radians(alpha)

    for i in range(len(x)):
        x_rotated.append(x[i] * math.cos(angle) - y[i] * math.sin(angle))
        y_rotated.append(x[i] * math.sin(angle) + y[i] * math.cos(angle))

    geometry_rotated = np.array([x_rotated, y_rotated]).T
    return geometry_rotated
    
