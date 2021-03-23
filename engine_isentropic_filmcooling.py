# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 16:39:05 2021

@author: Bernhard
"""

import nozzle_generator
import numpy as np
nozzle = nozzle_generator.Nozzle()
nozzle.sample()

step, x, r, P, T, rho, kappa, h, u, c, M, x_cf = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12

P_0 = 1000000
Fuel = "Ethanol"
Oxidizer = "N2O"
OFR = 2.5

data = np.array((12, len(nozzle)))

data[0,:] = nozzle.export_np()[:]

while True:
    for i in range(0, nozzle.get_n_t()):
        pass
    if abs(data[M,nozzle.get_n_t()] - 1) < 0.01:
        break
