# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 19:26:05 2022

@author: Bernhard Hansemann
"""

import thermodynamic_toolbox as thermo
from rocketcea.cea_obj import CEA_Obj
import scipy.constants as const
from scipy.optimize import fsolve
import math
from CoolProp.CoolProp import PropsSI
from gas_emittance import gas_emittance

#Umrechner
area = lambda r: r**2 * math.pi
psi = lambda x : x / const.psi
Pa = lambda x : x * const.psi
K = lambda x : const.convert_temperature(x, "R", "K")

def filmcooling(m_dot_g=0.274, P_cc = 10e5, r_cc = 0.02, ox = "N2O",
                fuel = "C2H5OH", ofr = 3, u_c = 40, l_c = 0.120, e_t = 0.1):
    
    #Gaszusammensetzung und Temperatur nach NASA CEA
    cea = CEA_Obj(fuelName=fuel, oxName=ox)
    mix_cea = cea.get_SpeciesMoleFractions(Pc=psi(P_cc), MR=ofr)[1] #Gaszusammensetzung in der Brennkammer
    mix_cc = thermo.prune_cc_products(mix_cea) #Formatieren und aussortieren
    
    T_cc = K(cea.get_Temperatures(Pc=psi(P_cc), MR=ofr)[0]) #Gastemperatur
    rho_g = thermo.get_rho(P_cc, T_cc, mix_cc)
    u_g = m_dot_g / (area(r_cc) * rho_g)
    
    #Prandtl Zahl des Gases
    cp_g = thermo.get_cp(P_cc, T_cc, mix_cc)
    mu_g = thermo.get_viscosity_mix(P_cc, T_cc, mix_cc)
    Pr_g = mu_g * cp_g / thermo.get_thermal_conductivity_mix(P_cc, T_cc, mix_cc)
    
    #Korrektur von h*
    cp_l = PropsSI("C", "P", 10e5, "Q", 0, "Ethanol")
    h_v = thermo.get_heat_of_vaporization(P_cc, "Ethanol")
    T_csat = PropsSI('T', 'P', P_cc, 'Q', 0, 'Ethanol')
    h_star = h_v + cp_l * (T_cc - T_csat)
    
    #Reynold Zahl des Gases
    G_mean = rho_g * u_g * (u_g-u_c) / u_g
    Re_g = G_mean * 2 * r_cc / mu_g
    
    #Stanton Zahl nach McKeon et al. und Friend and Metzner
    def f(x):
        return 1.93 * math.log(Re_g * x**0.5, 10) - 0.537 - x**-0.5
    lamb = fsolve(f, 0.01)[0]
    f = lamb/4
    St_0 = (f/2) / (1.2 + 11.8 * math.sqrt(f/2) * (Pr_g - 1) * Pr_g**(-1/3))
    
    #Korrektur fuer RMS-Turbulenz
    K_t = 1 + 4*e_t
    h_0 = G_mean * cp_g * St_0 * K_t
    
    #Waermestrahlung
    eta_g = gas_emittance(r_cc, T_cc, mix_cc["H2O"], mix_cc["CO2"], P_cc)
    q_dot_rad = eta_g * const.Boltzmann * (T_cc**4 - T_csat**4)
    
    #Korrektur von Stanton Zahl
    R = (PropsSI("M", thermo.mix_to_CP_string(mix_cc))/PropsSI("M", "Ethanol"))**0.6
    def g(x):
        (St, F) = x
        eqn_1 = math.log(1+R*F/St)/(F * R / St) - St / St_0
        eqn_2 = ((T_cc - T_csat)+q_dot_rad/h_0) * cp_g/h_star - F/St
        return [eqn_1, eqn_2]
    res = fsolve(g, (0.002, 0.001))
    St = res[0]
    
    
    #Adiabate Wandtemperatu
    r=Pr_g**(1/3)
    T_aw=T_cc*(1+r*(thermo.get_kappa(P_cc, T_cc, mix_cc)-1)/2*(u_g/thermo.get_speed_of_sound(P_cc, T_cc, mix_cc)))
    
    q_dot_conv = St * u_g * rho_g * cp_g * (T_aw - T_csat)* K_t
    
    m_dot_c = (q_dot_conv + q_dot_rad) * 2 * r_cc * math.pi * l_c / h_star
    return m_dot_c

if __name__ == "__main__":
    m_dot_c = filmcooling()
    print(f"Benoetigter Massenstrom [kg/s]: {m_dot_c}")