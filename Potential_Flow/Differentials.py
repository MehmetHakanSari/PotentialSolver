import numpy as np

def TwoDcentraldiff_simple(func, dx, dy):
    """
        func: 2D ndarray. Each node has its variables.
        dx  : 1D ndarray. I will update to 2D or not. 
        dy  : 1D ndarray. contains dy for each x  line. 

        If grid changes both in direction of y and x, for dx or dy, I am not sure how to implement those codes. 
    """
    (m, n) = np.shape(func)
    dfuncdx = np.zeros((m, n), dtype="float")
    dfuncdy = np.zeros((m, n), dtype="float")

    if n < 3:
        KeyError("array size is too small")
    if m < 3:
        KeyError("array size is too small")

    dfuncdx[:,0] = (-func[:,2] + 4 * func[:,1] - 3 * func[:,0]) / (2 * dx)
    dfuncdx[:,1:-1] = (func[:,2:] - func[:,0:-2]) / (2 * dx)
    dfuncdx[:,-1] = (func[:,-3] - 4 * func[:, -2] + 3 * func[:,-1]) / (2 * dx)

    dfuncdy[0,:] = (-func[2,:] + 4 * func[1,:] - 3 * func[0,:]) / (2 * dy)
    dfuncdy[1:-1,:] = (func[2:,:] - func[0:-2,:]) / (2 * dy)
    dfuncdy[-1,:] = (func[-3,:] - 4 * func[-2, :] + 3 * func[-1,:]) / (2 * dy)

    return dfuncdx, -dfuncdy

def TwoDcentraldiff_CD4(func, dx, dy):
    """
        func: 2D ndarray. Each node has its variables.
        dx  : 1D ndarray. I will update to 2D or not. 
        dy  : 1D ndarray. contains dy for each x  line. 

        If grid changes both in direction of y and x, for dx or dy, I am not sure how to implement those codes. 
    """
    (m, n) = np.shape(func)
    dfuncdx = np.zeros((m, n), dtype="float")
    dfuncdy = np.zeros((m, n), dtype="float")

    if n < 5:
        raise IndexError("array size is too small")
    if m < 5:
        raise IndexError("array size is too small")
  
    #for j in range(m):
        #for i in range(2,n-2):
           #dfuncdx[j,i] = (-func[j,i+2] + 8 * func[j,i+1] - 8 * func[j,i-1] + func[j,i-2]) / (12 * dx)

    dfuncdx[:,0:2] = (-3 * func[:,4:6] + 16 * func[:,3:5] - 36 * func[:,2:4] + 48 * func[:,1:3] - 25 * func[:,0:2]) / (12 * dx)
    dfuncdx[:,2:-2] = (-func[:,4:] + 8 * func[:,3:-1] - 8 * func[:,1:-3] + func[:,:-4]) / (12 * dx)
    dfuncdx[:,-2:] = (3 * func[:,-6:-4] - 16 * func[:,-5:-3] + 36 * func[:,-4:-2] - 48 * func[:,-3:-1] + 25 * func[:,-2:]) / (12 * dy)

    dfuncdy[0:2,:] = (-3 * func[4:6,:] + 16 * func[3:5,:]  - 36 * func[2:4,:] + 48 * func[1:3,:] - 25 * func[0:2,:]) / (12 * dy)
    dfuncdy[2:-2,:] = (-func[4:,:] + 8 * func[3:-1,:] - 8 * func[1:-3,:] + func[:-4,:]) / (12 * dx)
    dfuncdy[-2:,:] = (3 * func[-6:-4,:] - 16 * func[-5:-3,:]  + 36 * func[-4:-2,:] - 48 * func[-3:-1,:] + 25 * func[-2:,:]) / (12 * dy)

    return dfuncdx, dfuncdy


def OneDcentraldiff(func, dx, axis = 0):
    """
        Function: 1D or 2D array. Functions derivative is taken w.r.t given spacing. be careful that the shape should be (1,6) not (6,). 
        dx      : float or 1D array. Spacing of grid.
        axis    : the direction of the derivative. 0 is row derivative for 2D array. 1 is column derivatives.  

        Default setting is: 

        func: 1 x N
        dx:   1 x (N-1)
    """
    (m, n) = np.shape(func)
    if axis == 0:
        Operation = "row"
    else:
        Operation = "column"
    
    divflag = True
    if type(dx) == float:
        divflag = False
        if m > 1 and n == 1:        #if func is column vector
            Operation = "column"
    else:
        (k,l) = np.shape(dx)
        
        # func: 1D array
        if m > 1 and n == 1:        #if func is column vector
            Operation = "column"
            if l > k:               #but dx is given as row vector
                dx = dx.T
        elif n > 1 and m == 1:  #if func is row vector
            if k > l:           #but dx is given as column vector
                dx = dx.T
    
    if Operation == "row":
        dfuncdx = np.zeros((m,n), dtype="float")

        if n < 3:
            raise IndexError("array size is too small")

        if divflag:
            dfuncdx[:,0] = (-func[:,2] + 4 * func[:,1] - 3 * func[:,0]) / (dx[:,0] + dx[:,1])
            dfuncdx[:,1:-1] = (func[:,2:] - func[:,0:-2]) / (dx[:,0:-1] + dx[:,1:])
            dfuncdx[:,-1] = (func[:,-3] - 4 * func[:, -2] + 3 * func[:,-1]) / (dx[:,-1] + dx[:,-2])
        else:
            dfuncdx[:,0] = (-func[:,2] + 4 * func[:,1] - 3 * func[:,0]) / (2 * dx)
            dfuncdx[:,1:-1] = (func[:,2:] - func[:,:-2]) / (2 * dx)
            dfuncdx[:,-1] = (func[:,-3] - 4 * func[:, -2] + 3 * func[:,-1]) / (2 * dx) 

        return dfuncdx

    elif Operation == "column":
        dfuncdy = np.zeros((m,n), dtype="float")

        if m < 3:
            IndexError("array size is too small")

        if type(dx) != float:
            (k,l) = np.shape(dx)
            if l > k:           #but dx is given as column vector
                dx = dx.T   

        if divflag:
            dy = dx    #last elements should be 0 for y axis.
            dfuncdy[0,:] = (func[2,:] - 4 * func[1,:] + 3 * func[0,:]) / (dy[0,:] + dy[1,:])
            dfuncdy[1:-1,:] = (func[2:,:] - func[:-2,:]) / (dy[0:-1,:] + dy[1:,:])
            dfuncdy[-1,:] = (-func[-3,:] + 4 * func[-2, :] - 3 * func[-1,:]) / (dy[-1,:] + dy[-2,:])
        else:
            dy = dx    #last elements should be 0 for y axis.
            dfuncdy[0,:] = (-func[2,:] + 4 * func[1,:] - 3 * func[0,:]) / (2 * dy)
            dfuncdy[1:-1,:] = (func[2:,:] - func[:-2,:]) / (2 * dy)
            dfuncdy[-1,:] = (func[-3,:] - 4 * func[-2, :] + 3 * func[-1,:]) / (2 * dy)  

        return dfuncdy



def TDMA(W,C,E,Q):
    """
        The inputs should be row array.
        W: west    
        C: center
        E: east
        Q: source
    """
    n = len(Q)
    X = np.zeros(n)
    for i in range(1,n):
        C[i] = C[i] - E[i-1] * W[i] / C[i-1]
        Q[i] = Q[i] - Q[i-1] * W[i] / C[i-1]
    X[-1] = Q[-1] / C[-1]
    for i in range(n-2, -1, -1):
        X[i] = (Q[i] - E[i] * X[i+1]) / C[i]
    
    return X


def TwoDcentral_diff_velocity_CD2(solution):
    
    m, n = np.shape(solution.mesh.X)
    dfuncdx = np.zeros((m,m), dtype="float")
    dfuncdy = np.zeros((m,n), dtype="float")
    func = solution.solution
    
    if type(solution.mesh.xspacing) == float:
        dx = solution.mesh.xspacing
    else:
        dx = solution.mesh.xspacing[0,0]
    if type(solution.mesh.xspacing) == float:
        dy = solution.mesh.yspacing
    else:
        dy = solution.mesh.yspacing[0,0]
        

    dfuncdx[:,0] = (-func[:,2] + 4 * func[:,1] - 3 * func[:,0]) / (2 * dx) #Forward 
    dfuncdx[:,-1] = (func[:,-3] - 4 * func[:,-2] + 3 * func[:, -1]) / (2 * dx) #Backward

    dfuncdy[0,:] = -(-func[2,:] + 4 * func[1,:] - 3 * func[0,:]) / (2 * dy) #Forward 
    dfuncdy[-1,:] = -(func[-3,:] - 4 * func[-2,:] + 3 * func[-1, :]) / (2 * dy) #Backward

    for j in range(m):
         for i in range(1, n-1):
            if solution.map.area[j,i] == -2: #if it is wall
                #check where is the boundary on the east or west side of the wall
                if solution.map.area[j,i+1] == -1:
                    # 3 points backward difference
                    dfuncdx[j,i] = (func[j,i-2] - 4 * func[j,i-1] + 3 * func[j,i]) / (2 * dx) #Backward
                elif solution.map.area[j,i-1] == -1:
                    # 3 points forward difference
                    dfuncdx[j,i] = (-func[j,i+2] + 4 * func[j,i+1] - 3 * func[j,i]) / (2 * dx)
                else:
                    dfuncdx[j,i] = (func[j,i+1] - func[j,i-1]) / (2 * dx)                     #Central
            elif solution.map.area[j,i] == -1:
                dfuncdx[j,i] = 0
            else:
                dfuncdx[j,i] = (func[j,i+1] - func[j,i-1]) / (2 * dx)  #central difference

    for i in range(m):
        for j in range(1, m-1):  #If (y(0)) != (y = 0). the case when matrix index is not coordinate. 
            if solution.map.area[j,i] == -2:
                #check where is the boundary on the north or south side of the wall
                if solution.map.area[j+1,i] == -1: 
                    dfuncdy[j,i] = -(func[j-2,i] - 4 * func[j-1,i] + 3 * func[j,i]) / (2 * dy)  # 3 points backward difference
                elif solution.map.area[j-1,i] == -1:
                    dfuncdy[j,i] = -(-func[j+2,i] + 4 * func[j+1,i] - 3 * func[j,i]) / (2 * dy) # 3 points forward difference
                else:
                    dfuncdy[j,i] = -(func[j+1,i] - func[j-1,i]) / (2 * dy)  #central difference
            elif solution.map.area[j,i] == -1:
                dfuncdy[j,i] = 0
            else:
                dfuncdy[j,i] = -(func[j+1,i] - func[j-1,i]) / (2 * dy) #central difference

    return dfuncdx, dfuncdy

def TwoDcentral_diff_velocity_CD4(solution):
    
    m, n = np.shape(solution.mesh.X)
    dfuncdx = np.zeros((m,m), dtype="float")
    dfuncdy = np.zeros((m,n), dtype="float")
    func = solution.solution
    
    if type(solution.mesh.xspacing) == float:
        dx = solution.mesh.xspacing
    else:
        dx = solution.mesh.xspacing[0,0]
    if type(solution.mesh.xspacing) == float:
        dy = solution.mesh.yspacing
    else:
        dy = solution.mesh.yspacing[0,0]
        

    dfuncdx[:,0:2] = (-3 * func[:,4:6] + 16 * func[:,3:5] - 36 * func[:,2:4] + 48 * func[:,1:3] - 25 * func[:,0:2]) / (12 * dx)
    dfuncdx[:,-2:] = (3 * func[:,-6:-4] - 16 * func[:,-5:-3] + 36 * func[:,-4:-2] - 48 * func[:,-3:-1] + 25 * func[:,-2:]) / (12 * dx)

    dfuncdy[0:2,:] = (-3 * func[4:6,:] + 16 * func[3:5,:]  - 36 * func[2:4,:] + 48 * func[1:3,:] - 25 * func[0:2,:]) / (12 * dy)
    dfuncdy[-2:,:] = (3 * func[-6:-4,:] - 16 * func[-5:-3,:]  + 36 * func[-4:-2,:] - 48 * func[-3:-1,:] + 25 * func[-2:,:]) / (12 * dy)


    for j in range(m):
         for i in range(2, n-2):
            if solution.map.area[j,i] == -2: #if it is wall
                #check where is the boundary on the east or west side of the wall
                if solution.map.area[j,i+1] == -1:
                    # 5 points backward difference
                    dfuncdx[j,i] = (3 * func[j,i-4] - 16 * func[j,i-3] + 36 * func[j,i-2] - 48 * func[j,i-1] + 25 * func[j,i]) / (12 * dx)
                elif solution.map.area[j,i-1] == -1:
                    # 5 points forward difference
                    dfuncdx[j,i] = (-3 * func[j,i+4] + 16 * func[j,i+3] - 36 * func[j,i+2] + 48 * func[j,i+1] - 25 * func[j,i]) / (12 * dx)
                else:
                    dfuncdx[j,i] = (-func[j,i+2] + 8 * func[j,i+1] - 8 * func[j,i-1] + func[j,i-2]) / (12 * dx)
            elif solution.map.area[j,i] == -1:
                dfuncdx[j,i] = 0
            else:
                dfuncdx[j,i] = (-func[j,i+2] + 8 * func[j,i+1] - 8 * func[j,i-1] + func[j,i-2]) / (12 * dx)

    for i in range(m):
        for j in range(2, m-2):  #If (y(0)) != (y = 0). the case when matrix index is not coordinate. 
            if solution.map.area[j,i] == -2:
                #check where is the boundary on the north or south side of the wall
                if solution.map.area[j+1,i] == -1: 
                    dfuncdy[j,i] = (3 * func[j-4,i] - 16 * func[j-3,i]  + 36 * func[j-2,i] - 48 * func[j-1,i] + 25 * func[j,i]) / (12 * dy)  # 5 points backward difference
                elif solution.map.area[j-1,i] == -1:
                    dfuncdy[j,i] = (-3 * func[j+4,i] + 16 * func[j+3,i]  - 36 * func[j+2,i] + 48 * func[j+1,i] - 25 * func[j,i]) / (12 * dy)# 5 points forward difference
                else:
                    dfuncdy[j,i] = (-func[j+2,i] + 8 * func[j+1,i] - 8 * func[j-1,i] + func[j-2,i]) / (12 * dy)
            elif solution.map.area[j,i] == -1:
                dfuncdy[j,i] = 0
            else:
                dfuncdy[j,i] = (-func[j+2,i] + 8 * func[j+1,i] - 8 * func[j-1,i] + func[j-2,i]) / (12 * dy) #central difference

    return dfuncdx, -dfuncdy
             

def Solve_Coeff(x, y):
    '''
    input:
        x: the x coordinate of x_{i=i-1~i+1, j=j-1~j+1}, at least 3x3 array
        y: the y coordinate of y_{i=i-1~i+1, j=j-1~j+1}, at least 3x3 array
    output:
        a, b, c: at least 1x1 float
    '''
    a = 0.25 * (((x[1:-1, 2:] - x[1:-1, :-2])**2) + 
                ((y[1:-1, 2:] - y[1:-1, :-2])**2))
    b = 0.25 * ((x[2:, 1:-1] - x[:-2, 1:-1]) * 
                (x[1:-1, 2:] - x[1:-1, :-2]) + 
                (y[2:, 1:-1] - y[:-2, 1:-1]) * 
                (y[1:-1, 2:] - y[1:-1, :-2]))
    c = 0.25 * (((x[2:, 1:-1] - x[:-2, 1:-1])**2) + 
                ((y[2:, 1:-1] - y[:-2, 1:-1])**2))
    return a, b, c

def SolveEliptic(a, b, c, U):
    '''
    input:
        a, b, c: as described in the content
        U: the result of the last iteration
    output:
        return the result of current iteration
    '''
    return 0.5 * (
                  a * (U[2:, 1:-1] + U[:-2, 1:-1]) + 
                  c * (U[1:-1, 2:] + U[1:-1, :-2]) -
                  b * 0.5 * (U[2:, 2:] - U[2:, :-2] + U[:-2, :-2] - U[:-2, 2:])
                 ) / (a + c)

def pointwise(x_index, y_index, x_spacing, y_spacing, BCvalues, phi, phi_old, property_map, N_x, N_y, omega,type = "stream"):    
    if type == "stream":
        for i in x_index:

            if i == 0:  #if the west boundary is given neumann

                for j in y_index:
                    
                    if j == 0:   #if the north boundary is given neumann

                        dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2
                        
                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 

                    if j > 0 and j < N_y - 1:
                        
                        dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 

                    if j == N_y - 1:   #if the south boundary is given neumann

                        dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 


            elif i > 0 and i < N_x - 1:  #if the west boundary is given neumann

                for j in y_index:
                    
                    if j == 0:   #if the north boundary is given neumann

                        C_pro = (not(property_map[j,i] == -2)) * 1
                        W_pro = (not(property_map[j,i-1] == -2)) * 1
                        E_pro = (not(property_map[j,i+1] == -2)) * 1
                        S_pro = (not(property_map[j+1,i] == -2)) * 1

                        if C_pro == 0:
                            phi[j,i] = 0
                        else:
                            dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                            dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2  
                            
                            phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] * E_pro + phi[j, i-1] * W_pro) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i] * S_pro) 

                    if j > 0 and j < N_y - 1:

                        C_pro = (not(property_map[j,i] == -2)) * 1
                        W_pro = (not(property_map[j,i-1] == -2)) * 1
                        E_pro = (not(property_map[j,i+1] == -2)) * 1
                        S_pro = (not(property_map[j+1,i] == -2)) * 1
                        N_pro = (not(property_map[j-1,i] == -2)) * 1

                        if C_pro == 0:
                            phi[j,i] = 0
                        else:
                            dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                            dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2  

                            phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] * E_pro + phi[j, i-1] * W_pro) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] * S_pro + phi[j-1, i] * N_pro) 

                    if j == N_y - 1:

                        C_pro = (not(property_map[j,i] == -2)) * 1
                        W_pro = (not(property_map[j,i-1] == -2)) * 1
                        E_pro = (not(property_map[j,i+1] == -2)) * 1
                        N_pro = (not(property_map[j-1,i] == -2)) * 1

                        if C_pro == 0:
                            phi[j,i] = 0
                        else:
                            dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                            dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2

                            phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] * E_pro + phi[j, i-1] * W_pro) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i] * N_pro) 


            elif i == N_x - 1:

                for j in y_index:

                    # property_map(j,i)
                    # property_map(j,i-1)
                    # property_map(j+1,i)
                    # property_map(j-1,i)
                    
                    if j == 0:   #if the north boundary is given neumann
                        
                        dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 

                    if j > 0 and j < N_y - 1:

                        dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 

                    if j == N_y - 1:

                        dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i])
    
    elif type == "potensial":

        for i in x_index:

            if i == 0:  #if the west boundary is given neumann

                for j in y_index:

                    if j == 0:   #if the north boundary is given neumann

                        dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2

                        if (property_map[j,i] == -2):    #if it is a wall

                            #determine which sides are in the interior
                            
                            east = (property_map[j,i+1] == -1) 
                            south = (property_map[j+1,i] == -1) 

                            if east:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["N"]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["N"] + 2 * phi[j+1, i]) 
                            else:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["N"]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["N"] + 2 * phi[j+1, i]) 
                        elif (property_map[j,i] == -1):
                            phi[j,i] = 0
                        else:
                            phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["N"] + 2 * phi[j+1, i]) 

                        
                    if j > 0 and j < N_y - 1:
                        
                        dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2 

                        if (property_map[j,i] == -2):    #if it is a wall

                            #determine which sides are in the interior
                            
                            east = (property_map[j,i+1] == -1)
                            south = (property_map[j+1,i] == -1) 
                            north = (property_map[j-1,i] == -1) 

                            if east:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                                elif north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i])  
                            else:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1] ) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                                elif north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1] ) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                        elif (property_map[j,i] == -1):
                            phi[j,i] = 0
                        else:
                            phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 

                        
                    if j == N_y - 1:   #if the south boundary is given neumann

                        dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i]) / 2)**2

                        if (property_map[j,i] == -2):    #if it is a wall

                            #determine which sides are in the interior
                            
                            east = (property_map[j,i+1] == -1)
                            north = (property_map[j-1,i] == -1) 

                            if east:
                                if north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                            else:
                                if north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                        elif (property_map[j,i] == -1):
                            phi[j,i] = 0
                        else:
                            phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["W"] + 2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 

                        


            elif i > 0 and i < N_x - 1:  #interior points   

                for j in y_index:
                    
                    if j == 0:   #if the north boundary is given neumann

                        dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2  

                        if (property_map[j,i] == -2):    #if it is a wall

                            #determine which sides are in the interior     
                            east = (property_map[j,i+1] == -1)
                            west = (property_map[j,i-1] == -1)  
                            south = (property_map[j+1,i] == -1) 
                            
                            if east:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["N"])  
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]+ 2 * np.sqrt(dy2) * BCvalues["N"])    
                            elif west:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["N"])
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i] + 2 * np.sqrt(dy2) * BCvalues["N"])    
                            else:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["N"])
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i] + 2 * np.sqrt(dy2) * BCvalues["N"])    
                        elif (property_map[j,i] == -1):
                            phi[j,i] = 0
                        else:
                            phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i] + 2 * np.sqrt(dy2) * BCvalues["N"])                        
                            
                            
                    if j > 0 and j < N_y - 1:

                        dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2  

                        if (property_map[j,i] == -2):    #if it is a wall

                            #determine which sides are in the interior
                            
                            east = (property_map[j,i+1] == -1)
                            west = (property_map[j,i-1] == -1)  
                            south = (property_map[j+1,i] == -1) 
                            north = (property_map[j-1,i] == -1) 

                            if east:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                                elif north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                            elif west:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                                elif north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                            else:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                                elif north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                        elif (property_map[j,i] == -1):
                            phi[j,i] = 0
                        else:
                            phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                            

                    if j == N_y - 1:  #if the south boundary is given neumann 

                        dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i] + x_spacing[j,i-1]) / 2)**2

                        if (property_map[j,i] == -2):    #if it is a wall

                            #determine which sides are in the interior
                            
                            east = (property_map[j,i+1] == -1)
                            west = (property_map[j,i-1] == -1)  
                            north = (property_map[j-1,i] == -1) 

                            if east:
                                if north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i-1] ) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["S"]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i-1] ) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i] + 2 * np.sqrt(dy2) * BCvalues["S"]) 
                            elif west:
                                if north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i+1]) +  dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["S"]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * phi[j, i+1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i] + 2 * np.sqrt(dy2) * BCvalues["S"]) 
                            else:
                                if north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["S"]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i] + 2 * np.sqrt(dy2) * BCvalues["S"]) 
                        elif (property_map[j,i] == -1):
                            phi[j,i] = 0
                        else:
                            phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (phi[j, i+1] + phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i] + 2 * np.sqrt(dy2) * BCvalues["S"]) 




            elif i == N_x - 1:  #if the east boundary is given neumann

                for j in y_index:
                    
                    if j == 0:   #if the north boundary is given neumann
                        
                        dy2 = ((y_spacing[j,i] + y_spacing[j,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        if (property_map[j,i] == -2):    #if it is a wall
                            #determine which sides are in the interior
                            
                            west = (property_map[j,i-1] == -1)  
                            south = (property_map[j+1,i] == -1) 

                            if west:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["N"]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i] + 2 * np.sqrt(dy2) * BCvalues["N"]) 
                            else:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["N"]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i] + 2 * np.sqrt(dy2) * BCvalues["N"]) 
                        elif (property_map[j,i] == -1):
                            phi[j,i] = 0
                        else:
                            phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i] + 2 * np.sqrt(dy2) * BCvalues["N"]) 
   

                    if j > 0 and j < N_y - 1:

                        dy2 = ((y_spacing[j,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        if (property_map[j,i] == -2):    #if it is a wall

                            #determine which sides are in the interior
                            
                            west = (property_map[j,i-1] == -1)  
                            south = (property_map[j+1,i] == -1) 
                            north = (property_map[j-1,i] == -1) 

                            if west:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] ) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                                elif north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] ) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] ) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                            else:
                                if south:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1] ) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j-1, i]) 
                                elif north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1] ) + dx2 / (2 * dy2 + 2 * dx2) * (2 * phi[j+1, i]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1] ) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 
                        elif (property_map[j,i] == -1):
                            phi[j,i] = 0
                        else:
                            phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (phi[j+1, i] + phi[j-1, i]) 


                    if j == N_y - 1: #if the south boundary is given neumann

                        dy2 = ((y_spacing[j-1,i] + y_spacing[j-1,i]) / 2)**2
                        dx2 = ((x_spacing[j,i-1] + x_spacing[j,i-1]) / 2)**2  

                        if (property_map[j,i] == -2):    #if it is a wall

                            #determine which sides are in the interior
                            
                            west = (property_map[j,i-1] == -1)  
                            north = (property_map[j-1,i] == -1) 

                            if west:
                                if north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["S"]) 
                                else:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["S"] + 2 * phi[j-1, i]) 
                            else:
                                if north:
                                    phi[j,i] =  dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["S"]) 
                                else:
                                    phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["S"] + 2 * phi[j-1, i]) 
                        elif (property_map[j,i] == -1):
                            phi[j,i] = 0
                        else:
                            phi[j,i] = dy2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dx2) * BCvalues["E"] + 2 * phi[j, i-1]) + dx2 / (2 * dy2 + 2 * dx2) * (2 * np.sqrt(dy2) * BCvalues["S"] + 2 * phi[j-1, i])




        phi = ((property_map != -1) * 1) * phi 
        phi = (1 - omega) * phi_old + omega * phi
                    
    return phi 



def column_TDMA(x_spacing, y_spacing, phi, y_index, x_index, BC_values, property_map, N_y, N_x, W, E):
    """
        x & yspacing: dx's and dy's stored in a 2D or 1D ndarrays or a float. 
        phi         : property of the field, 2D ndarray
        y x_index   : indicies for unsolved phi matrix. 
        BC_values   : BC values for neumann Condiation
        property_map: 2D ndarray including objects in the field. 0 for open, -1 for interior, -2 for wall 
        N_y & N_x   : Number of nodes 
        W           : west side of the TDMA
        E           : east side of the TDMA

        Well I will write with constant dx and dy for the time being. If we already use Jacabi Transformation, dx and dy will be constant.
    """

    dx = x_spacing[0,0]
    dy = y_spacing[0,0]


    for i in x_index:

        line = property_map[:, i] 

        #find wall and interior points by checking line. if the element in line -1 it is interior, if it is 2 it is wall:


        inter_indicies = np.nonzero((line == -1))[0]
        wall_indicies = np.nonzero((line == -2))[0] #check wheter nonzero gives a row or column ndarray

        if i == 0:  # if the west wall is neumann
            W[1:] = -dy
                                       
            C = 2 * dy + 2 * dx
            E[:-1] = -dy                            
            Q = 2 * dx * phi[y_index[0]:y_index[-1]+1,(i+1)] + 2 * dx**2 * BC_values['W']   #conditions are set for neumann BC.                
            Q[0] += dy * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * dy**2 * BC_values['N']
            Q[-1] += dy * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * dy**2 * BC_values['S']
        
        if i > 0 and i < N_x-1:
            W[1:] = -dy                          
            C = 2 * dy + 2 * dx
            E[:-1] = -dy                            
            Q = dx * phi[y_index[0]:y_index[-1]+1,(i-1)] + dx * phi[y_index[0]:y_index[-1]+1,(i+1)] #conditions are set for neumann BC.                
            Q[0] += dy * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * dy**2 * BC_values['N']
            Q[-1] += dy * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * dy**2 * BC_values['S']

        if i == N_x - 1: # if the east wall is neumann
            W[1:] = -dy                      
            C = 2 * dy + 2 * dx
            E[:-1] = -dy                           
            Q = 2 * dx * phi[y_index[0]:y_index[-1]+1,(i-1)] + 2 * dx**2 * BC_values['E']             
            Q[0] += dy * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * dy**2 * BC_values['N']
            Q[-1] += dy * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * dy**2 * BC_values['S']

        #we will update east and west by checking lines, if the line is -1 it is interior, if it is -2 it is wall. It the interior is in south, update the east by -dy, if it is in north, update the west by -dy.
        #wall indicies are the indicies of the walls. The index j will be used for updateing E and W.
        for j in wall_indicies:
            if line[j-1] == -1: #north is interior
                W[j] += -dy
            if line[j+1] == -1: #south is interior
                E[j] += -dy

        #the left and right sides will be checked for wall. Update Q accordingly. Use property map to find left and right side of the wall
        for j in wall_indicies:
            if property_map[j, i-1] == -1:
                Q[j] += dx * phi[j,i+1]
            if property_map[j, i+1] == -1:
                Q[j] += dx * phi[j,i-1]

        Q = np.flip(Q)                     #The reason of reversing Q is, existing Q is inconsistent with the W and E and C list.

        W[-1] += -dy * (y_index[0] == 0)             #Neumann of N-S boundaries. implemented here. 
        E[0] += -dy * (y_index[-1] == N_y - 1)        #probabaly for different spacing matrixies the east and west should fliped
        
        if y_index[0] == 0:
            phi[y_index[-1]::-1,i] = TDMA(W,C,E,Q) 
        else:
            phi[y_index[-1]:y_index[0] - 1:-1,i] = TDMA(W,C,E,Q)

        #interior points makes phi zero.
        for j in inter_indicies:
            phi[j,i] = 0 

    return phi



    """
    for i in x_index:
        if i == 0:  # if the west wall is neumann
            W[1:] = -y_spacing[y_index[0]:y_index[-1], i]                          
            C = 2 * y_spacing[y_index[0]:, i] + 2 * y_spacing[y_index[0]:y_index[-1]+1, i - 1]
            E[:-1] = -y_spacing[y_index[0]:y_index[-1], i]                            
            Q = 2 * a_e[y_index[0]:y_index[-1]+1,i] * phi[y_index[0]:y_index[-1]+1,(i+1)] + 2 * a_e[y_index[0]:y_index[-1]+1,i]**2 * BC_values['W']   #conditions are set for neumann BC.                
            Q[0] += a_n[0,0] * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * a_n[0,i - 1]**2 * BC_values['N']
            Q[-1] += a_s[0,0] * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * a_s[0,i - 1]**2 * BC_values['S']
        
        if i > 0 and i < N_x-1:
            W[1:] = -a_s[y_index[0]:y_index[-1], i]                          
            C = 2 * a_s[y_index[0]:, i] + 2 * a_w[y_index[0]:y_index[-1]+1, i - 1]
            E[:-1] = -a_n[y_index[0]:y_index[-1], i]                            
            Q = a_w[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i-1)] + a_e[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i+1)] #conditions are set for neumann BC.                
            Q[0] += a_n[0,0] * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * a_n[0,i - 1]**2 * BC_values['N']
            Q[-1] += a_s[0,0] * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * a_s[0,i - 1]**2 * BC_values['S']

        if i == N_x - 1: # if the east wall is neumann
            W[1:] = -a_s[y_index[0]:y_index[-1], i]                          
            C = 2 * a_s[y_index[0]:, i] + 2 * a_w[y_index[0]:y_index[-1]+1, i - 1]
            E[:-1] = -a_n[y_index[0]:y_index[-1], i]                            
            Q = 2 * a_w[y_index[0]:y_index[-1]+1,i-1] * phi[y_index[0]:y_index[-1]+1,(i-1)] + 2 * a_w[y_index[0]:y_index[-1]+1,i-1]**2 * BC_values['E']             
            Q[0] += a_n[0,0] * phi[0,i-1] * (y_index[0] == 1) + (y_index[0] == 0) * 2 * a_n[0,i - 1]**2 * BC_values['N']
            Q[-1] += a_s[0,0] * phi[-1,i-1] * (y_index[-1] == N_y - 2) + (y_index[-1] == N_y - 1) * 2 * a_s[0,i - 1]**2 * BC_values['S']
        

        Q = np.flip(Q)                     #The reason of reversing Q is, existing Q is inconsistent with the W and E and C list.

        W[-1] += -a_s[y_index[-1], i] * (y_index[0] == 0)             #Neumann of N-S boundaries. implemented here. 
        E[0] += -a_n[y_index[0], i] * (y_index[-1] == N_y - 1)        #probabaly for different spacing matrixies the east and west should fliped
        
        if y_index[0] == 0:
            phi[y_index[-1]::-1,i] = TDMA(W,C,E,Q) 
        else:
            phi[y_index[-1]:y_index[0] - 1:-1,i] = TDMA(W,C,E,Q) 

    return phi
    """




# def Solve_Coeffs(self, X, Y):
#         """
#             input:
#                 x: the x coordinate of x_{i=i-1~i+1, j=j-1~j+1}, at least 3x3 array
#                 y: the y coordinate of y_{i=i-1~i+1, j=j-1~j+1}, at least 3x3 array
#             output:
#                 a, b, c: at least 1x1 float
#         """
#         # X = X.T
#         # Y = Y.T

#         alpha = np.zeros(X.shape)
#         beta = np.zeros(X.shape)
#         gamma = np.zeros(X.shape)

#         # alpha[1:-1,1:-1] = 0.25 * (((X[1:-1, 2:] - X[1:-1, :-2])**2) + 
#         #             ((Y[1:-1, 2:] - Y[1:-1, :-2])**2))
#         # beta[1:-1,1:-1] = 0.25 * ((X[2:, 1:-1] - X[:-2, 1:-1]) * 
#         #             (X[1:-1, 2:] - X[1:-1, :-2]) + 
#         #             (Y[2:, 1:-1] - Y[:-2, 1:-1]) * 
#         #             (Y[1:-1, 2:] - Y[1:-1, :-2]))
#         # gamma[1:-1,1:-1] = 0.25 * (((X[2:, 1:-1] - X[:-2, 1:-1])**2) + 
#         #             ((Y[2:, 1:-1] - Y[:-2, 1:-1])**2))

#         alpha[1:-1,1:-1] = 0.25 * (((X[2:, 1:-1] - X[:-2, 1:-1])**2) + 
#                     ((Y[2:, 1:-1] - Y[:-2, 1:-1])**2))
#         beta[1:-1,1:-1] = 0.25 * ((X[1:-1, 2:] - X[1:-1, :-2]) * 
#                     (X[2:, 1:-1, ] - X[:-2, 1:-1]) + 
#                     (Y[1:-1, 2:] - Y[1:-1, :-2]) * 
#                     (Y[2:, 1:-1] - Y[:-2, 1:-1]))
#         gamma[1:-1,1:-1] = 0.25 * (((X[1:-1, 2:] - X[1:-1, :-2])**2) + 
#                     ((Y[1:-1, 2:] - Y[1:-1, :-2])**2))

#         # alpha = alpha.T
#         # beta = beta.T
#         # gamma = gamma.T

#         return alpha[1:-1,1:-1], beta[1:-1,1:-1], gamma[1:-1,1:-1]

# def Solve_Eliptic(self, alpha, beta, gamma, U):
#     """
#         input:
#             a, b, c: as described in the content
#             U: the result of the last iteration
#         output:
#             return the result of current iteration
#     """
#     # alpha = alpha.T
#     # beta = beta.T
#     # gamma = gamma.T
#     # U = U.T

#     # U[1:-1, 1:-1] = ((-0.5) / (alpha[1:-1, 1:-1] + gamma[1:-1, 1:-1] + 1e-9)) * \
#     # ( 2 * beta[1:-1, 1:-1] * (U[2:, 2:] - U[2:, :-2] - U[:-2, 2:] + U[:-2, :-2]) \
#     # - alpha[1:-1, 1:-1] * (U[1:-1, 2:] + U[1:-1 , :-2])  \
#     # -gamma[1:-1, 1:-1] * (U[2:, 1:-1] + U[:-2, 1:-1]))

#     U = 0.5 * (
#                 alpha * (U[1:-1, 2:] + U[1:-1, :-2]) + 
#                 beta * (U[2:, 1:-1] + U[:-2, 1:-1]) -
#                 gamma * 0.5 * (U[2:, 2:] - U[:-2, 2:] + U[:-2, :-2] - U[2:, :-2])
#                 ) / (alpha + gamma)   #updated to Y,X coordinate

#     return U