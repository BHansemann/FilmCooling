# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 11:45:55 2020

@author: Berni
"""

import numpy as np
import matplotlib.pyplot as plt
import math

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
    vget_y = np.vectorize(get_y)
    nozzle[2,:] = vget_y(nozzle[1,:])
    
    return nozzle

if __name__ == "__main__":
    sample_nozzle = conical_nozzle(0.02, 0.03, 0.04)
    print(sample_nozzle)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(sample_nozzle[1,:], sample_nozzle[2,:])
    ax.axis('equal')
    plt.show()
