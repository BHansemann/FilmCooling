# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 16:46:00 2020

@author: Bernhard
"""

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj import CEA_Obj
import scipy.constants as const
import thermodynamic_toolbox as thermo

psi = lambda x : x / const.psi
Pa = lambda x : x * const.psi

def prune_cc_products(cea_dict, cutoff=0.01, exclude=["*NO", "*O", "*OH"], pos=1):
    cp_dict = {}
    total = 0
    for key in cea_dict:
        if (cea_dict[key][pos] > cutoff) and not (key in exclude):
            cp_dict[key.strip('*')] = cea_dict[key][pos]
            total += cea_dict[key][pos]
    for key in cp_dict:
        cp_dict[key] = cp_dict[key]/total
    return cp_dict

def engine_adiabatic(P_cc, P_e, n_ps, T_cc, F_th, eta_nz, ox, fuel, ofr, u_cc=1):
    step, P, T, rho, kappa, h, u, c, M, r = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    cea = CEA_Obj(fuelName=fuel, oxName=ox)
    mix_cc = prune_cc_products(cea.get_SpeciesMoleFractions(Pc=psi(P_cc), MR=ofr)[1])
    eng = np.zeros((10, n_ps + 1))
    eng[step,:] = np.arange(0, n_ps + 1)
    eng[P,:] = P_cc - ((P_cc - P_e)/n_ps) * eng[0,:]
    eng[T,0] = cea.get_Temperatures(Pc=psi(P_cc), MR=ofr)[0]
    eng[rho,0] = thermo.get_rho(eng[P,0], eng[T,0], mix_cc)
    eng[kappa,0] = thermo.get_kappa(eng[P,0], eng[T,0], mix_cc)
    eng[h,0] = thermo.mass_mixer(mix_cc, eng[P,0], 'P', eng[T,0], 'T', 'H')
    eng[u,0] = u_cc
    eng[c,0] = thermo.get_speed_of_sound(eng[P,0], eng[T,0], mix_cc)
    eng[M,0] = eng[u,0] / eng[c,0]
    print(eng[:,0])
    for i in eng[step, :]:
        i = int(i)
        if(i != 0):
            eng[T,i] = thermo.isentrop_temp_press(eng[kappa,i-1], eng[T,i-1], eng[P,i-1], eng[P,i])
            eng[rho,i] = thermo.get_rho(eng[P,i], eng[T,i], mix_cc)
            eng[kappa,i] = thermo.get_kappa(eng[P,i], eng[T,i], mix_cc)
            eng[h,i] = thermo.mass_mixer(mix_cc, eng[P,i], 'P', eng[T,i], 'T', 'H')
            eng[u,i] = ((eng[h,i-1] - eng[h,i]) * 2 + eng[u,i-1]**2)**0.5
            eng[c,i] = thermo.get_speed_of_sound(eng[P,i], eng[T,i], mix_cc)
            eng[M,i] = eng[u,i] / eng[c,i]
    m_dot = F_th / eng[u,-1]
    eng[r,:] = np.vectorize(lambda u, rho: (m_dot / (math.pi * u * rho))**0.5)(eng[u,:], eng[rho,:])
    return eng
    
if __name__ == "__main__":
    step, P, T, rho, kappa, h, u, c, M, r = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    sample_engine = engine_adiabatic(1000000, 100000, 100, 0, 500,0,"N2O","C2H5OH",3)
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(212)
    ax1.plot(sample_engine[step,:], sample_engine[M,:])
    ax2.plot(sample_engine[step,:], sample_engine[u,:])
    plt.show()