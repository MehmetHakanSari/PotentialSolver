import numpy as np


def AFyt(x, t, c):
    '''
    input: 
        x: x coordinate of center line, float
        t: maximum thickness, in fraction of chord length, float
        c: chord lrngth, float
    output:
        half thickness of airfoil at corresponding x coordinate
    '''
    return 5. * t * (0.2969 * ((x/c)**0.5) - 
                     0.126 * (x/c) - 
                     0.3516 * ((x/c)**2) + 
                     0.2843 * ((x/c)**3) - 
                     0.1036 * ((x/c)**4))

def AFyc(x, m, p, c):
    '''
    input:
        x: x coordinate of center line, float
        m: the maximum camber (100 m is the first of the four digits), float
        p: location of maximum camber (10 p is the second digit), float
        c: chord lrngth, float
    output:
        y coordinate of center line at corresponding x coordinate
    '''
    if (x >= 0.0) and (x <= p*c):
        return m * x * (2. * p - (x/c)) / (p**2.)
    elif (x > p*c) and (x <= c):
        return m * (c - x) * (1. + (x/c) - 2. * p) / ((1. - p)**2)
    else:
        raise ValueError


def AFth(x, m, p, c):
    '''
    input:
        x: x coordinate of center line, float
        m: the maximum camber (100 m is the first of the four digits), float
        p: location of maximum camber (10 p is the second digit), float
        c: chord lrngth, float
    output:
        angle between center and horizontal line at corresponding x coordinate
    '''
    if (x >= 0.0) and (x <= p*c):
        return np.arctan(2.0 * m * (p - (x/c)) / (p**2))
    elif (x > p*c) and (x <= c):
        return np.arctan(2.0 * m * (p - (x/c)) / ((1. - p)**2))
    else:
        raise ValueError

def AF(x, t, sign, m, p, c):
    '''
    input:
        x: x coordinate of center line, float
        t: maximum thickness, in fraction of chord length, float
        sign: indicate upper (1) or lower (-1) surface of airfoil
        m: the maximum camber (100 m is the first of the four digits), float
        p: location of maximum camber (10 p is the second digit), float
        c: chord lrngth, float
    output:
        x, y coordinates on airfoil surface at corresponding 
        center line x coordinate
    '''
    if (m == 0.) or (p == 0):
        return x, sign * AFyt(x, t, c)
    else:
        return np.array([x[i] - 
                         sign * AFyt(x[i], t, c) * np.sin(AFth(x[i], m, p, c)) 
                         for i in range(np.size(x))]), \
               np.array([AFyc(x[i], m, p, c) + 
                         sign * AFyt(x[i], t, c) * np.cos(AFth(x[i], m, p, c))
                         for i in range(np.size(x))])


