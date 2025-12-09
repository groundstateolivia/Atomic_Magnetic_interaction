# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 16:07:17 2025

@author: Olivia
"""

import numpy as np
# import AngularMomentum as ang

def landeG(J,L1,L2,g1,g2):
    def term(g, a, b):
        return g*(J*(J+1.)+a*(a+1.)-b*(b+1.))/(2.*J*(J+1.))
    return term(g2,L2,L1) + term(g1,L1,L2)


if __name__ == '__main__':
    S = 0.5
    L = 0
    J = 0.5
    I = 1
    F = 1.5
    gs = 2.0023010
    gL = 0.99999587
    gI=-0.0004476540
    gJ = landeG(J,S,L,gs,gL)
    
    print(gJ)
    
    gF = landeG(F,I,J,gI,gJ)
    
    print(gF)