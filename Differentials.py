import numpy as np

def TwoDcentraldiff(func, dx, dy):
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

def TDMA(W,C,E,Q):
    """
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