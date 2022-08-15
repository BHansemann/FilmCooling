# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 17:21:32 2022

@author: Bernhard
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter

data_wht = np.array([[-78, -78, -78, -35, -35, -35, 0, 0, 0, 23, 23, 23],
                      [1.5428, 1.7002, 1.8577, 2.2985, 2.0466, 2.4244, 9.6661, 8.4382, 9.6346, 12.4369, 11.3978, 12.2480]])

data_blk = np.array([[-78, -78, -35, -35, 0, 0, 23, 23, 23],
                     [2.2985, 2.3299, 2.6448, 2.4559, 21.9771, 23.2365, 35.8938, 36.1771, 35.7993]])

data_nincl = np.array([[-78, -35, 0],[15.6799, 3.9672, 15.7429]])

avg_wht = np.array([[-78, -35, 0, 23],[1.7002, 2.2565, 9.2463, 12.0276]])
avg_blk = np.array([[-78, -35, 0, 23],[2.3142, 2.5503, 22.6068, 35.9567]])

xvals = np.linspace(-78, 23, 100)

#tck_wht = inter.splrep(avg_wht[0], avg_wht[1], k=3)
#tck_blk = inter.splrep(avg_blk[0], avg_blk[1], k=3)

#spline_wht = inter.splev(xvals, tck_wht)
#spline_blk = inter.splev(xvals, tck_blk)

f_wht = inter.interp1d(avg_wht[0], avg_wht[1], 'slinear')
f_blk = inter.interp1d(avg_blk[0], avg_blk[1], 'slinear')

spline_wht = f_wht(xvals)
spline_blk = f_blk(xvals)

fig1 = plt.figure(figsize=(9,6))
ax1  = fig1.add_subplot(1,1,1)

ax1.plot(data_wht[0], data_wht[1], 'c+', label="Weiße Proben")
ax1.plot(data_blk[0], data_blk[1], 'k+', label="Schwarze Proben")
ax1.plot(data_nincl[0], data_nincl[1], 'rx', label="Ausreißer")
ax1.plot(xvals, spline_wht, 'c-')
ax1.plot(xvals, spline_blk, 'k-')

plt.xticks([-78, -35, 0, 23])
plt.grid(ls='--')
plt.xlabel("Temperatur [°C]")
plt.ylabel("Charpy-Kerbschlagzähigkeit [kJ/m²]")

legend = ax1.legend(loc="upper left")