# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 19:26:05 2022

@author: Bernhard
"""

import thermodynamic_toolbox as thermo
from rocketcea.cea_obj import CEA_Obj
import scipy.constants as const
from scipy.optimize import fsolve
import math
from CoolProp.CoolProp import PropsSI

#Umrechner
area = lambda r: r**2 * math.pi
psi = lambda x : x / const.psi
Pa = lambda x : x * const.psi
K = lambda x : const.convert_temperature(x, "R", "K")

#Eingangsgrößen
m_dot_g = 0.274 #kg/s
P_cc = 10e5 #Pa
r_cc = 0.02 #m
ox = "N2O"
fuel = "C2H5OH"
ofr = 3
u_c = 40 #m/s
l_c = 0.120

cea = CEA_Obj(fuelName=fuel, oxName=ox)
mix_cea = cea.get_SpeciesMoleFractions(Pc=psi(P_cc), MR=ofr)[1] #Gaszusammensetzung in der Brennkammer
mix_cc = thermo.prune_cc_products(mix_cea) #Formatieren und aussortieren
T_cc = K(cea.get_Temperatures(Pc=psi(P_cc), MR=ofr)[0]) #Gastemperatur
print(T_cc)
rho_g = thermo.get_rho(P_cc, T_cc, mix_cc)
u_g = m_dot_g / (area(r_cc) * rho_g)
cp_g = thermo.get_cp(P_cc, T_cc, mix_cc)
mu_g = thermo.get_viscosity_mix(P_cc, T_cc, mix_cc)
print(mu_g)
Pr_g = mu_g * cp_g / thermo.get_thermal_conductivity_mix(P_cc, T_cc, mix_cc)
cp_l = PropsSI("C", "P", 10e5, "Q", 0, "Ethanol")
h_v = thermo.get_heat_of_vaporization(P_cc, "Ethanol")
T_csat = PropsSI('T', 'P', P_cc, 'Q', 0, 'Ethanol')

G_mean = rho_g * u_g * (u_g-u_c) / u_g
Re_g = G_mean * 2 * r_cc / mu_g

def f(x):
    return 1.93 * math.log(Re_g * x**0.5, 10) - 0.537 - x**-0.5
lamb = fsolve(f, 0.01)[0]
f = lamb/4
St_0 = (f/2) / (1.2 + 11.8 * math.sqrt(f/2) * (Pr_g - 1) * Pr_g**(-1/3))
e_t = 0.15
K_t = 1 + 4*e_t


R = (PropsSI("M", thermo.mix_to_CP_string(mix_cc))/PropsSI("M", "Ethanol"))**0.6
def g(x):
    return (math.log(1+R*x[1]/x[0])/(x[1] * R / x[0]) - x[0] / St_0, )#########################