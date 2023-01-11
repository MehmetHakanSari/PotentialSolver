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

    dfuncdy[0,:] = (func[2,:] - 4 * func[1,:] + 3 * func[0,:]) / (2 * dy)
    dfuncdy[1:-1,:] = (-func[2:,:] + func[0:-2,:]) / (2 * dy)
    dfuncdy[-1,:] = (-func[-3,:] + 4 * func[-2, :] - 3 * func[-1,:]) / (2 * dy)

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
            dfuncdx[:,1:-1] = (func[:,2:] - func[:,0:-2]) / (2 * dx)
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
            dfuncdy[1:-1,:] = (-func[2:,:] + func[0:-2,:]) / (dy[0:-1,:] + dy[1:,:])
            dfuncdy[-1,:] = (-func[-3,:] + 4 * func[-2, :] - 3 * func[-1,:]) / (dy[-1,:] + dy[-2,:])
        else:
            dy = dx    #last elements should be 0 for y axis.
            dfuncdy[0,:] = (func[2,:] - 4 * func[1,:] + 3 * func[0,:]) / (2 * dy)
            dfuncdy[1:-1,:] = (-func[2:,:] + func[0:-2,:]) / (2 * dy)
            dfuncdy[-1,:] = (-func[-3,:] + 4 * func[-2, :] - 3 * func[-1,:]) / (2 * dy)

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


def TwoDcentral_diff_velocity(solution):
    
    n, m = np.shape(solution.mesh.matricies[0])
    dfuncdx = np.zeros((n,m), dtype="float")
    dfuncdy = np.zeros((n,m), dtype="float")
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

    for j in range(n):
         for i in range(1, m-1):
            if solution.map.area[j,i] == -2:
                #check where is the boundary on the east or west side of the wall
                if solution.map.area[j,i+1] == -1:
                    # 3 points backward difference
                    dfuncdx[j,i] = (func[j,i-2] - 4 * func[j,i-1] + 3 * func[j,i]) / (2 * dx) #Backward
                elif solution.map.area[j,i-1] == -1:
                    # 3 points forward difference
                    dfuncdx[j,i] = (-func[j,i+2] + 4 * func[j,i+1] - 3 * func[j,i]) / (2 * dx)
                else:
                    #central difference
                    dfuncdx[j,i] = (func[j,i+1] - func[j,i-1]) / (2 * dx)                     #Central
            elif solution.map.area[j,i] == -1:
                dfuncdx[j,i] = 0
            else:
                #central difference
                dfuncdx[j,i] = (func[j,i+1] - func[j,i-1]) / (2 * dx)

    for i in range(m):
        for j in range(1, n-1):  #If (y(0)) != (y = 0). the case when matrix index is not coordinate. 
            if solution.map.area[j,i] == -2:
                #check where is the boundary on the north or south side of the wall
                if solution.map.area[j+1,i] == -1: 
                    # 3 points backward difference
                    dfuncdy[j,i] = -(func[j-2,i] - 4 * func[j-1,i] + 3 * func[j,i]) / (2 * dy)
                elif solution.map.area[j-1,i] == -1:
                    # 3 points forward difference
                    dfuncdy[j,i] = -(-func[j+2,i] + 4 * func[j+1,i] - 3 * func[j,i]) / (2 * dy)
                else:
                    #central difference
                    dfuncdy[j,i] = -(func[j+1,i] - func[j-1,i]) / (2 * dy)
            elif solution.map.area[j,i] == -1:
                dfuncdy[j,i] = 0
            else:
                #central difference
                dfuncdy[j,i] = -(func[j+1,i] - func[j-1,i]) / (2 * dy)

    return dfuncdx, dfuncdy
             

def Solve_a_b_c(x, y):
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

def SolveEq(a, b, c, U):
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