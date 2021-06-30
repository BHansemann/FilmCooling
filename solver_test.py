# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 16:42:05 2021

@auth"or: Berni
"""

"""import sympy
import math
import numpy as np

#P1, P2, T1, T2, Rs, r1, r2, u1, u2, Qc, hcc, Fr, h1, h2 = sympy.symbols("P1 P2 T1 T2 Rs r1 r2 u1 u2 Qc hcc Fr h1 h2")

#f = (P1 * r1**2 * math.pi * u1) / (Rs * T1) + Qc / hcc - (P2 * r2**2 * math.pi * u2) / (Rs * T2)
#g = (P1 * r1**2 * math.pi * u1 * h1) / (Rs * T1) + Qc * hcc / hcc - Qc - (P2 * r2**2 * math.pi * u2 * h2) / (Rs * T2)
#h = P1/(Rs * T1) * u1**2 * r1**2 * math.pi + P1 * r1**2 * math.pi - Fr - P2/(Rs*T2) * u2**2 * r2**2 * math.pi - P2 * r2**2 * math.pi

def funcmatrix(P_2, T_2, u_2):
    return np.array([P_2 * u_2 * A_2, , ])

def jacobimatrix():
    return np.array([[], [], []], 
                    [[], [], []], 
                    [[], [], []])"""

import thermodynamic_toolbox as thermo

mix = {'CO': 0.22927143797319066,
 'CO2': 0.03355950272079416,
 'H2': 0.21465190951430582,
 'H2O': 0.17864231755549284,
 'N2': 0.3438748322362165}

thermo.get_viscosity_mix(1000000, 2000, mix)