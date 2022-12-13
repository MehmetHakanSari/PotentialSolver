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
            IndexError("array size is too small")

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