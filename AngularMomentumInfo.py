'''
Docstring for AngularMomentumInfo

Functions to give the angular momentum matricies and eigen values given a spin
'''


import numpy as np

pi = np.pi

def NumStates(J):
    '''
    Given a value for the angular momentum, returns number of mf values
    '''
    if not(2 * J == int(2 * J)):
        raise ValueError('J must be a half integer')
    return int(2 * J + 1)

def angular_index(J,m):
    '''
    Gives the index of m in the standard basis
    
    :param J: spin number
    :param m: magnetic quantum number
    '''

    index_float = J + m
    index_int = int(index_float)

    assert index_float == index_int

    return index_int

def stateVector(J,m):
    '''
    Returns standard basis state vector given spin and magnetic number
    
    :param J: spin number
    :param m: magnetic number
    '''
    
    N = NumStates(J)
    state = np.zeros(N, dtype=complex)
    state[angular_index(J,m)] = 1
    return state

def LadderOp(J,sign):
    '''
    Given a sign and spin number, returns the matrix representation for the ladder operator associated with it.
    
    :param J: spin number
    :param sign = 1: raising operator
    :param sign = -1: raising operator
    '''
    if sign != 1 and sign != -1:
        raise ValueError('Sign must be +/- 1')
    
    N = NumStates(J)
    matrix = np.zeros((N,N), dtype=complex)
    mjRange = np.arange(-J,J) if sign == 1 else np.arange(-J+1,J+1)

    for m in mjRange:
        matrix[angular_index(J,m+sign),angular_index(J,m)] = np.sqrt(J*(J+1)-m*(m+sign))
    
    return matrix

def idMatrix(J):
    '''
    Identity matrix
    '''
    return np.eye(NumStates(J))

def J_z(J):
    '''
    returns matrix representation for z angular momentum operator
    
    :param J: spin number
    '''
    N = NumStates(J)
    matrix = np.zeros((N,N), dtype=complex)

    for m in np.arange(-J,J+1):
        matrix[angular_index(J,m),angular_index(J,m)] = m

    return matrix

def J_x(J):
    '''
    returns matrix representation for x angular momentum operator
    
    :param J: spin number
    '''
    return 0.5*(LadderOp(J,1)+LadderOp(J,-1))

def J_y(J):
    '''
    returns matrix representation for y angular momentum operator
    
    :param J: spin number
    '''
    return -1j*0.5*(LadderOp(J,1)-LadderOp(J,-1))

def dotProd(I,J):
    return np.kron(J_x(I), J_x(J)) + np.kron(J_y(I), J_y(J)) + np.kron(J_z(I), J_z(J))
