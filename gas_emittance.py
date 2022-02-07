# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 01:22:47 2021

@author: Simon
"""
#Spaceteam

#-----------------VARAIBALEN
import numpy as np
import scipy.interpolate as interpolate
import scipy.constants as const

D =  0.04 #m
T=2600
N_w = 0.843 #frei gewählt
N_c = 0.157 #frei gewählt
p=10e5
def gas_emittance(r_cc, T, N_w, N_c, p):
    p=p/const.atm #Pascal -> atm
    
    L_eff = 0.95 * r_cc * 2 #m
    
    f_c_H2O = interpolate.interp1d((1000, 2000, 3000), (0.165, 0.9, 2.05), kind="quadratic")
    f_n_H2O = interpolate.interp1d((1000, 2000, 3000), (0.45, 0.65, 0.61), kind="quadratic")
    
    f_c_CO2 = interpolate.interp1d((1000, 1500, 2000), (0.05, 0.075, 0.15), kind="quadratic", fill_value="extrapolate")
    
    c_Co2 = f_c_CO2(T)
    n_Co2 = 0.6 #Konstant ?
    
    c_H2o = f_c_H2O(T)
    n_H2o = f_n_H2O(T)
    
    
    R_H = 461.4 #Wikipedia - spez. Gaskonstante
    R_C = 188.9 # ditto
    
    rho_H = p/(R_H*T) #ideale Gasgleichung - keine Datenblätter zu H20 oder Co2 bei 2300K gefunden
    rho_C = p/(R_C*T)
    
    rho_opt = p * L_eff
    rho_H = rho_opt * N_w
    rho_C = rho_opt * N_c
    
    m = 100 * rho_C
    
    #------------Correction Factor
    
    C1 = 0.26 + 0.74 * np.exp(-2.5 * rho_H)
    C2 = 0.75 + 0.31 * np.exp(-10 * rho_H) #hier auch dichte in atm gemeint ?
    
    Kp_H = 1 + C1 *(1 - np.exp((1-p*(1+N_w))/C2))
    
    Kp_C = 10**(0.036 * rho_C**(-0.483)*(1+(2*np.log(p)**(-m)))**(-1/m))
    
    
    print("Korrekturfaktur Wasser", "  ", Kp_H)
    print("Korrekturfaktor Co2", "  ", Kp_C)
    
    
    n = 5.5 *(1+(1.09*rho_opt)**-3.88)**(-1/3.88)
    Kx = 1 - (abs(((2*N_w)/(N_w + N_c))-1))**n
    
    delta_e = 0.0551 * Kx * (1 - np.exp(-4*rho_opt)) * (1 - np.exp(-12.5 * rho_opt))
    
    #---------------Emittance
    
    
    
    E_f_h = 0.825 #pdf
    E_h = E_f_h * ( 1 + (rho_H/c_H2o)**(-n_H2o))**(-1/n_H2o)
    
    E_f_c = 0.231 #pdf
    E_c = E_f_c*(1+(rho_C/c_Co2)**(-n_Co2))**(-1/n_Co2)
    
    E_H = E_h * Kp_H
    E_C = E_c * Kp_C
    
    E_ges = E_H + E_C + delta_e
    
    print("total gas emmittance", "  ", E_ges, "[Einheit]")
    return E_ges


#-----------------------Heat flux
#sigma = 5.67 * 10**(-8)
#Aw = #Wallabsorptivity - Materialkonstante
#Tv = # saturation temperature of coolant
#Tg = T
#
#Q_rad = sigma * Aw * E_ges*(Tg**(4)-Tv**(4))
#print("Heatflux in ...",Q_rad)