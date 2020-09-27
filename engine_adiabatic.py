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
import scipy.constants as const

psi = lambda x : x / const.psi
Pa = lambda x : x * const.psi

def convert_to_CP_string(cea_dict, cutoff=0.01, exclude=["*NO", "*O", "*OH"], pos=1):
    cp_string = ""
    cp_dict = {}
    total = 0
    for key in cea_dict:
        if (cea_dict[key][pos] > cutoff) and not (key in exclude):
            cp_dict[key.strip('*')] = cea_dict[key][pos]
            total += cea_dict[key][pos]
    for key in cp_dict:
        cp_string += "{}[{}]&".format(key, cp_dict[key]/total)
    return cp_string.strip("&")

def engine_adiabatic(P_cc, P_e, n_ps, T_cc, F_th, eta_nz, ox, fuel, ofr):
    step, P, T, rho, kappa, h, u, c, M, r = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    cea = CEA_Obj(fuelName=fuel, oxName=ox)
    eng = np.zeros((10, n_ps + 1))
    eng[step,:] = np.arange(0, n_ps + 1)
    eng[P,:] = P_cc - ((P_cc - P_e)/n_ps) * eng[0,:]
    eng[T,0] = cea.get_Temperatures(Pc=psi(P_cc), MR=ofr)[0]
    mix_cc = convert_to_CP_string(cea.get_SpeciesMoleFractions(Pc=psi(P_cc), MR=ofr)[1])
    
    #print(cea.get_Temperatures(Pc=psi(P_cc), MR=ofr))
    
if __name__ == "__main__":
    sample_engine = engine_adiabatic(1000000, 100000, 100, 0,0,0,"N2O","C2H5OH",3)
