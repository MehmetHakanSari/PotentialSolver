import os
import numpy as np

def prepare_design_matrix(data_path, selection='all', num_samples=None, randomize=False, **kwargs):
    """
    Function to prepare the design matrix from data stored in objects.

    Parameters:
    - data_path (str): Path to the directory containing the data objects.
    - selection (str): Selection method for data. Options: 'all', 'random', 'specified'.
    - num_samples (int): Number of samples to select if selection is 'random'. Ignored if selection is 'all'.
    - randomize (bool): Whether to randomly select samples. Applicable only if selection is 'random'.
    - kwargs (dict): Additional keyword arguments for specifying data selection criteria.

    Returns:
    - design_matrix (numpy.ndarray): Processed design matrix.
    """


    data_objects = []
    for filename in os.listdir(data_path):
        if filename.endswith('.pkl'):
            file_path = os.path.join(data_path, filename)
            data_objects.append(file_path)  
    
    if selection == 'all':
        selected_data = data_objects
    elif selection == 'random':
        if randomize:
            np.random.shuffle(data_objects)  
        selected_data = data_objects
    elif selection == 'specified':
        pass 

    for file_path in selected_data:
        with open(file_path, 'rb') as f:
            airfoil = pickle.load(f)
            m_list = [float(airfoil.name[4])]         #m = max camber
            p_list = [float(airfoil.name[5])]         #p = location of max camber
            t_list = [float(airfoil.name[6:])]        #t = thickness
            keys_list = list(airfoil.polars.keys())
            for key in keys_list:
                if key == 0:
                    pass
                else:
                    data = airfoil.polars[key]
                    AoA = data[:,0] 
                    [float(i) for i in AoA]
                    Cl = data[:,1]
                    Cd = data[:,2]
                    Cm = data[:,3]
                    
                    input_data_matrix = np.array([AoA, [float(key)] * len(AoA), m_list * len(AoA), p_list * len(AoA), t_list * len(AoA)]).T
                    output_data_matrix = np.array([Cl, Cd, Cm]).T
                    num_input_rows, num_input_cols = input_data_matrix.shape
                    num_output_rows, num_output_cols = output_data_matrix.shape
                    if num_input_rows != num_output_rows:
                        UserWarning 

                    with open('input_data.txt', 'a') as f:
                        np.savetxt(f, input_data_matrix, fmt='%f', delimiter='\t')
                    with open('output_data.txt', 'a') as f:
                        np.savetxt(f, output_data_matrix, fmt='%f', delimiter='\t')

    return None


data_path = '../data/airfoil_data/NACAxxxx_pkl/'
prepare_design_matrix(data_path, selection='all', num_samples=None, randomize=False)
