# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 16:46:00 2020

@author: Berni
"""

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj import CEA_Obj

def engine_adiabatic(P_cc, P_e, n_ps, T_cc, F_th, eta_nz, ox, fuel, ofr):
    P, T, rho, kappa, h, u, c, M, r = 1, 2, 3, 4, 5, 6, 7, 8, 9
    cea = CEA_Obj(fuelName=fuel, oxName=ox)
    cea.setupCards(pc_units='bar', Pc=P_cc*10**-5, MR=ofr)
    eng = np.zeros((10, n_ps + 1))
    eng[0,:] = np.arange(0, n_ps + 1)
    eng[P,:] = P_cc - ((P_cc - P_e)/n_ps) * eng[0,:]
    eng[T,0] = cea.get_Temperatures(Pc=P_cc, MR=ofr)[0]
    
    print(eng[T,0])
    
if __name__ == "__main__":
    sample_engine = engine_adiabatic(1000000, 100000, 100, 0,0,0,"N2O","C2H5OH",3)

