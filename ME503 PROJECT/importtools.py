import numpy as np

def import_geometry(file_name):
    """
    The file should contain x and y coordinates in each column. 
    Though for different type of coordinate stroges I can update inside. 
    """
    
    file = []

    with open(file_name) as f:
        for line in f:
            file.append(line.split())

    #convert the data to floats
    for i in range(len(file)):
        for j in range(len(file[i])):
            file[i][j] = float(file[i][j])

    return np.array(file)
