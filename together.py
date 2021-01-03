# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 13:16:36 2020

@author: Bernhard
"""

from engine_adiabatic import engine_adiabatic
import contour_generator
import matplotlib.pyplot as plt
import numpy as np
import math
from datetime import datetime

step, P, T, rho, kappa, h, u, c, M, r = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
sample_engine = engine_adiabatic(1000000, 100000, 1000, 0, 500, 0, "N2O", "C2H5OH", 2.5, 0.02)

delta=1000
step_throat=0
for i in sample_engine[step,:]:
    i = int(i)
    #print(abs(1-sample_engine[M,i]))
    if abs(1-sample_engine[M,i]) < delta:
        delta = abs(1-sample_engine[M,i])
        step_throat = i

chamber = 0.02
throat = sample_engine[r,step_throat]
r_exit = sample_engine[r,-1]

nozzle = contour_generator.bell_nozzle(throat, chamber, r_exit, n_steps = 200)
con = math.radians(25)
l_t = (chamber - throat * (1 + 1.5 * (1 - math.cos(con)))) / math.tan(con) + throat * 1.5 * math.sin(con)
i_t = int(200*l_t/nozzle[1,-1])
nozzle = np.append(nozzle, np.zeros((2,201)), axis=0)

for i in nozzle[0,:i_t]:
    i = int(i)
    diff = 1000
    dp = 0
    for j in sample_engine[step,:step_throat]:
        j = int(j)
        if abs(nozzle[2,i]-sample_engine[r,j]) < diff:
            dp=j
            diff=abs(nozzle[2,i]-sample_engine[r,j])
    nozzle[3,i] = sample_engine[P,dp]
    nozzle[4,i] = sample_engine[M,dp]
    
for i in nozzle[0,i_t:]:
    i = int(i)
    diff = 1000
    dp = 0
    for j in sample_engine[step,step_throat+1:]:
        j = int(j)
        if abs(nozzle[2,i]-sample_engine[r,j]) < diff:
            dp=j
            diff=abs(nozzle[2,i]-sample_engine[r,j])
    nozzle[3,i] = sample_engine[P,dp]
    nozzle[4,i] = sample_engine[M,dp]

fig, (ax1, ax2) = plt.subplots(2, sharex=True)
ax1.plot(nozzle[1,:], nozzle[2,:], linewidth=1)
ax2.plot(nozzle[1,:], nozzle[3,:], linewidth=1)
ax1.axis("equal")
plt.title(datetime.now())
#ax1.vlines(l_t, 0, 0.02, linewidth=0.1)
#ax1.hlines(throat, 0, 0.04, linewidth=0.1)
#plt.savefig("graph.svg")
plt.show()