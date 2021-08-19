# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 16:39:05 2021

@author: Bernhard
"""

import nozzle_generator
import numpy as np
import thermodynamic_toolbox as thermo
import math
from rocketcea.cea_obj import CEA_Obj
import scipy.constants as const
import CoolProp.CoolProp as CP

from engine_adiabatic import engine_adiabatic

area = lambda r: r**2 * math.pi
psi = lambda x : x / const.psi
Pa = lambda x : x * const.psi

columns = 13
step, x, r, P, T, rho, kappa, h, u, c, M, Q, m_dot_c = tuple(range(0,columns))

F_th = 500
P_0 = 1000000
P_e = 100000
fuel = "C2H5OH"
oxidizer = "N2O"
OFR = 2.5

r_cc = 0.02

Fr=0 #ayayay
rec_fac = 1 #recovery Factor

def newton_solver(funcs, jacobi, initguess, cutoff=0.01):
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
            #print(jacobi(xn))
            print(xn)
            print(funcmatrix(xn))
            print("")
        
        if np.any(abs(z)>cutoff):
            xn += z
            #print(z)
        else:
            print(i)
            print(abs(z))
            print(abs(z)<cutoff)
            return xn

#initialize
cea = CEA_Obj(fuelName=fuel, oxName=oxidizer)
mix = thermo.prune_cc_products(cea.get_SpeciesMoleFractions(Pc=psi(P_0), MR=OFR)[1])
Rs = thermo.get_specific_gas_constant(mix)
init_eng = engine_adiabatic(P_0, P_e, 1000, 0, F_th, 0, oxidizer, fuel, OFR, r_cc)

m_dot = F_th/init_eng[6,-1]
m_low = m_dot/2
m_high = m_dot*1.5

r_t = float(init_eng[9].min())
r_e = float(init_eng[9,-1])
print("Chamber Radius: " + str(r_cc))
print("Throat Radius: " + str(r_t))
print("Exit Radius: " + str(r_e))

#alpha_con = float(input("Converging angle (default 35°): ") or "35")
#alpha_divt = float(input("Diverging initial angle (default 21°): ") or "21")
#alpha_dive = float(input("Converging angle (default 11°): ") or "11")
#l_c = float(input("Length of Chamber (default 0.15): ") or "0.15")

alpha_con = 35
alpha_divt = 21
alpha_dive = 11
l_c = 0.15
l_fc = 0.14

nozzle = nozzle_generator.Nozzle()
nozzle.quickset(r_cc, r_t, r_e, alpha_con, alpha_divt, alpha_dive, 1.5, 0.382, l_c, 0, "bell", 0.9, 1000, "m")
nozzle.generate()

data = np.zeros((len(nozzle.export_np()), columns))

data[:,0:3] = nozzle.export_np()[:]

data[0,P] = P_0
data[0,T] = init_eng[2,0]
data[0,rho] = thermo.get_rho(data[0,P], data[0,T], mix)
data[0,kappa] = thermo.get_kappa(data[0,P], data[0,T], mix)
data[0,h] = thermo.mass_mixer(mix, data[0,P], 'P', data[0,T], 'T', 'H')
data[0,c] = thermo.get_speed_of_sound(data[0,P], data[0,T], mix)

while True:
    n=1
    data[0,u] = m_dot / (area(r_cc) * data[0,rho])
    data[0,M] = data[0,u] / data[0,c]
    data[0,Q] = -1 * thermo.get_cone_area(data[0,r], data[0,r], data[1,x]) * thermo.get_heat_transfer_coefficient(data[0,P], data[0,T], data[0,r], data[0,u], mix) * (data[0,T] - thermo.get_adiabatic_wall_temperature(data[0,P], data[0,T], data[0,u], rec_fac, mix))
    data[0,m_dot_c] = data[0,Q] / thermo.get_heat_of_vaporization(data[0,P], "Ethanol")
    try:
        for i in range(1, nozzle.get_n_t()):
            in_cool = data[i][x] < l_fc
            s1 = data[i-1]
            s2 = data[i]
            h_ev = thermo.get_heat_of_vaporization(s1[P], "Ethanol")
            mdotf = (s1[Q] / h_ev) if in_cool else 0
            h_f = CP.PropsSI('H','P',s1[P],'T', s1[T],'Ethanol')
            
            s2[h] = s1[h] #workaround
            
            def funcmatrix(vec):
                P_2, T_2, u_2 = tuple(vec)
                return np.array([P_2 * u_2 * area(s1[r]) / (T_2 * Rs) - s1[rho] * s1[u] * area(s1[r]) - s1[Q] / h_ev,
                                 P_2 * u_2 * area(s1[r]) * s2[h] / (T_2 * Rs) - s1[rho] * s1[u] * area(s1[r]) * s1[h] - mdotf * h_f + s1[Q],
                                 P_2 * u_2**2 * area(s1[r]) / (T_2 * Rs) + P_2 * area(s1[r]) - s1[rho] * s1[u]**2 * area(s1[r]) - s1[P] * area(s1[r]) + Fr])
            
            def jacobimatrix(vec):
                P_2, T_2, u_2 = tuple(vec)
                alpha = u_2 * area(s1[r]) / (T_2 * Rs)
                return np.array([[alpha,
                                 alpha * -1 * P_2 / T_2,
                                 alpha * P_2 / u_2
                                 ],
                                
                                [alpha * s2[h],
                                 alpha * -1 * P_2 * s2[h] / T_2,
                                 alpha * P_2 * s2[h] / u_2
                                 ],
                                
                                [alpha * u_2 + area(s2[r]),
                                 alpha * -1 * P_2 * u_2 / T_2,
                                 alpha * 2 * P_2
                                ]])
            
            P_2, T_2, u_2 = newton_solver(funcmatrix, jacobimatrix, (s1[P]*1.1, s1[T]*1.1, s1[u]*1.1))
            
            print(P_2, T_2, u_2)
            
            s2[P] = P_2
            s2[T] = T_2
            s2[u] = u_2
            
            s2[rho] = thermo.get_rho(s2[P], s2[T], mix)
            s2[kappa] = thermo.get_kappa(s2[P], s2[T], mix)
            s2[h] = thermo.mass_mixer(mix, s2[P], "P", s2[T], "T", "H")
            s2[c] = thermo.get_speed_of_sound(s2[P], s2[T], mix)
            s2[M] = s2[u] / s2[c]
            s2[Q] = -1 * thermo.get_heat_transfer_coefficient(s2[P], s2[T], s2[r], s2[u], mix) if in_cool else 0
            s2[m_dot_c] = s2[Q] / h_ev
            
            data[i] = s2
            
        delta = abs(data[nozzle.get_n_t(),M] - 1)
        print("Loop " + n + ", Delta: " + delta)
        n += 1
            
        if abs(data[nozzle.get_n_t(),M] - 1) < 0.01:
            print("Valid mass flow rate found at " + str(m_dot))
            break
        elif data[nozzle.get_n_t(),M] < 1:
            m_high = m_dot
            m_dot = (m_dot + m_low) / 2
        elif data[nozzle.get_n_t(),M] > 1:
            m_low = m_dot
            m_dot = (m_dot + m_high) / 2
        
    except np.linalg.LinAlgError:
        print("linalg error")
        m_low = m_dot
        m_dot = (m_dot + m_high) / 2
            

print("YAY")
