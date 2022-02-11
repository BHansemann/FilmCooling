# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 16:42:05 2021

@author: Berni
"""
import numpy as np
from math import sin
from math import cos

def funcmatrix(vec):
    x, y = tuple(vec)
    return np.array([x**3 + sin(y) + 1, x*y + 1])

def jacobimatrix(vec):
    x, y = tuple(vec)
    return np.array([[3*x**2, cos(y)], [y, x]])

def newton_solver(funcs, jacobi, initguess, cutoff=0.001):
    xn = np.array(initguess)
    i=0
    while(True):
        i += 1
        a = jacobi(xn)
        b = funcs(xn)*-1
        try:
            z = np.linalg.solve(a, b)
        
        except np.linalg.LinAlgError as err:
            #print(err)
            print(jacobi(xn))
            print(funcmatrix(xn))
            print("")
        
        if np.any(abs(z)>cutoff):
            xn += z
            print(z)
        else:
            print(i)
            print(abs(z))
            print(abs(z)<cutoff)
            return xn
        
test = newton_solver(funcmatrix, jacobimatrix, (0.,4.))

co=1e100
for P_test in np.arange(1e6*0.99, 1e6*1.01, 100):
    print(P_test)
    for T_test in np.arange(2300*0.99, 2300*1.01, 100):
        for u_test in np.arange(200*0.99, 200*1.01, 100):
            res = funcmatrix((P_test, T_test, u_test))
            if sum(abs(res)) < co:
                co = sum(abs(res))
                absres = res
                
print(absres)
print(co)
