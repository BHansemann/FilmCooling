# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 16:39:05 2021

@author: Bernhard
"""

import nozzle_generator
import numpy as np
import thermodynamic_toolbox as thermo
import math

nozzle = nozzle_generator.Nozzle()
nozzle.sample()

area = lambda r: r**2 * math.pi

step, x, r, P, T, rho, kappa, h, u, c, M, Q, x_cf = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12

P_0 = 1000000
Fuel = "Ethanol"
Oxidizer = "N2O"
OFR = 2.5

Fr=0 #ayayay

#initialize
Rs=1
h_ev=1

data = np.array((12, len(nozzle)))

data[0,:] = nozzle.export_np()[:]

while True:
    for i in range(1, nozzle.get_n_t()):
        s1 = np.array(data[:,i-1])
        s2 = np.array(s1.shape())
        def funcmatrix(P_2, T_2, u_2):
            return np.array([P_2 * u_2 * area(s1[r]) / (T_2 * Rs) - s1[rho] * s1[u] * area(s1[r]) - s1[Q] / h_ev,
                             P_2 * u_2 * area(s1[r]) * s2[h] / (T_2 * Rs) - s1[rho] * s1[u] * area(s1[r]) + s1[h] - mdotf * h_f + s1[Q], #Massenstrom vom film cooling berechnen
                             P_2 * u_2**2 * area(s1[r]) / (T_2 * Rs) + P_2 * area(s1[r]) - s1[rho] * s1[u]**2 * area(s1[r]) - s1[P] * area(s1[r]) + Fr])
        
        def jacobimatrix(P_2, T_2, u_2):
            alpha = u_2 * area(s1[r]) / (T_2 * Rs)
            return np.array([alpha,
                             alpha * -1 * P_2 / T_2,
                             alpha * P_2 / u_2
                             ],
                            
                            [alpha * s2[h],
                             alpha * -1 * P_2 * s2[h] / T_2,
                             alpha * P_2 * s2[h] / u_2
                             ],
                            
                            [alpha * u_2 + area(s2[r]),
                             alpha * -1 * P_2 * u_2,
                             alpha * 2 * P_2
                            ])
        
    if abs(data[M,nozzle.get_n_t()] - 1) < 0.01:
        break
