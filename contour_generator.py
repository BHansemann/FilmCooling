# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 11:45:55 2020

@author: Bernhard
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import bezier
import sympy
from sympy.solvers import solve

def conical_nozzle(r_t, r_ch, r_e, n_steps = 10000, eta = 0, div = 15, con = 25, rf_c = 1, rf_d = 1):
    div = math.radians(div)
    con = math.radians(con)
    if(r_e == 0):
        r_e = eta**0.5 * r_t
    nozzle = np.zeros((3, n_steps + 1))
    nozzle[0,:] = np.arange(0, n_steps + 1)
    nozzle[2,0] = r_ch
    l_1 = (r_ch - r_t * (1 + rf_c * (1 - math.cos(con)))) / math.tan(con)
    l_2 = r_t * rf_c * math.sin(con)
    l_3 = r_t * rf_d * math.sin(div)
    l_4 = (r_e - r_t * (1 + rf_d * (1 - math.cos(div)))) / math.tan(div)
    l_sum = l_1 + l_2 + l_3 + l_4
    print(l_3+l_4)
    print(((r_e-r_t)+rf_d*r_t*((math.cos(math.radians(15))**-1) - 1))/math.tan(math.radians(15)))
    nozzle[1,:] = nozzle[0,:] * (l_sum / n_steps)
    def get_y(x):
        if(x <= l_1):
            return r_ch - x * math.tan(con)
        elif l_1 < x <= (l_1 + l_2):
            return r_t * (1 + rf_c * (1 - math.cos(math.asin((l_2 - x + l_1) / (r_t * rf_c)))))
        elif(l_1 + l_2) < x <= (l_sum - l_4):
            return r_t * (1 + rf_d * (1 - math.cos(math.asin((x - l_1 - l_2) / (r_t * rf_d)))))
        elif(l_sum - l_4) < x <= l_sum:
            return r_t * (1 + rf_d * (1 - math.cos(div))) + (x - (l_1 + l_2 + l_3)) * math.tan(div)
    nozzle[2,:] = np.vectorize(get_y)(nozzle[1,:])
    return nozzle

def bell_nozzle(r_t, r_ch, r_e, n_steps = 10000, l_f = 0.8, eta = 0, div_i = 30, div_f = 10, con = 25, rf_c = 1.5, rf_d = 0.382):
    div_i = math.radians(div_i)
    div_f = math.radians(div_f)
    con = math.radians(con)
    if(eta != 0):
        r_e = eta**0.5 * r_t
    nozzle = np.zeros((3, n_steps + 1))
    nozzle[0,:] = np.arange(0, n_steps + 1)
    nozzle[2,0] = r_ch
    l_1 = (r_ch - r_t * (1 + rf_c * (1 - math.cos(con)))) / math.tan(con)
    l_2 = r_t * rf_c * math.sin(con)
    l_3 = r_t * rf_d * math.sin(div_i)
    l_4 = l_f * (((r_e - r_t) + rf_d * r_t* ((math.cos(math.radians(15)) **-1) - 1)) / math.tan(math.radians(15))) - l_3
    l_sum = l_1 + l_2 + l_3 + l_4
    x_0 = l_sum - l_4
    y_0 = r_t * (1 + rf_d * (1 - math.cos(math.asin((x_0 - l_1 - l_2) / (r_t * rf_d)))))
    x_2 = l_sum
    y_2 = r_e
    x_1 = (y_0 - math.tan(div_i) * x_0 - y_2 + math.tan(div_f) * x_2) / (math.tan(div_f) - math.tan(div_i))
    y_1 = math.tan(div_i) * x_1 + y_0 - math.tan(div_i) * x_0
    bell = bezier.Curve(np.asfortranarray([[x_0, x_1, x_2], [y_0, y_1, y_2]]), degree=2)
    y = sympy.Symbol("y")
    x = sympy.Symbol("x")
    bellexpr = solve(bell.implicitize(), y)
    bellfunc = [sympy.lambdify(x, bellexpr[0], "numpy"), sympy.lambdify(x, bellexpr[1], "numpy")]
    nozzle[1,:] = nozzle[0,:] * (l_sum / n_steps)
    def get_y(x):
        if(x <= l_1):
            return r_ch - x * math.tan(con)
        elif l_1 < x <= (l_1 + l_2):
            return r_t * (1 + rf_c * (1 - math.cos(math.asin((l_2 - x + l_1) / (r_t * rf_c)))))
        elif(l_1 + l_2) < x <= (l_sum - l_4):
            return r_t * (1 + rf_d * (1 - math.cos(math.asin((x - l_1 - l_2) / (r_t * rf_d)))))
        elif(l_sum - l_4) < x <= l_sum:
            if(bellfunc[0](x) > bellfunc[1](x)):
                return bellfunc[0](x)
            else:
                return bellfunc[1](x)
    nozzle[2,:] = np.vectorize(get_y)(nozzle[1,:])
    return nozzle

if __name__ == "__main__":
    #sample_nozzle = conical_nozzle(0.02, 0.03, 0.04)
    sample_nozzle = bell_nozzle(0.02, 0.03, 0.04)
    #print(sample_nozzle)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(sample_nozzle[1,:], sample_nozzle[2,:])
    ax.axis('equal')
    #plt.savefig("graph.svg")
    plt.show()
    
    
