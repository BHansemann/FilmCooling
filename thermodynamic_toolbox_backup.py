# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 12:45:02 2020

@author: Bernhard
"""

import math
import scipy.constants as const
from scipy.optimize import fsolve
import CoolProp.CoolProp as CP
import numpy as np

def isentrop_press_temp(kappa, press, temp1, temp2) -> float:
    return press * (temp2 / temp1)**(kappa / (kappa - 1))

def isentrop_press_vol(kappa, press, vol1, vol2) -> float:
    return press * (vol1 / vol2)**kappa

def isentrop_temp_press(kappa, temp, press1, press2) -> float:
    return temp * (press2 / press1)**((kappa - 1) / kappa)

def isentrop_temp_vol(kappa, temp, vol1, vol2) -> float:
    return temp * (vol2 / vol1)**(1-kappa)

def isentrop_vol_temp(kappa, vol, temp1, temp2) -> float:
    return vol * (temp2 / temp1)**(1 / (1-kappa))

def isentrop_vol_press(kappa, vol, press1, press2) -> float:
    return vol * (press1 / press2)**(1 / kappa)

def mix_to_CP_string(mix: dict) -> str:
    out = ""
    for key in mix:
        out += "{}[{}]&".format(key, mix[key])
    return out.strip("&")

def mass_frac(mix: dict) -> dict:
    massmix = {}
    molarmass_mix = CP.PropsSI('M', mix_to_CP_string(mix))
    for key in mix:
        massmix[key] = mix[key] * CP.PropsSI('M', key) / molarmass_mix
    return massmix

def molar_mixer(mix: dict, input1, type1, input2, type2, output) -> float:
    val = 0
    for key in mix:
        val += mix[key] * CP.PropsSI(output, type1, input1, type2, input2, key)
    return val

def mass_mixer(mix: dict, input1, type1, input2, type2, output) -> float:
    val=0
    massmix = mass_frac(mix)
    for key in massmix:
        val += massmix[key] * CP.PropsSI(output, type1, input1, type2, input2, key)
    return val

def get_kappa(P, T, mix: dict) -> float:
    return molar_mixer(mix, P, 'P', T, 'T', "CPMOLAR") / molar_mixer(mix, P, 'P', T, 'T', "CVMOLAR")

def get_rho(P, T, mix: dict) -> float:
    return P / (T * CP.PropsSI("GAS_CONSTANT", mix_to_CP_string(mix)) / CP.PropsSI("M", mix_to_CP_string(mix)))

def get_specific_gas_constant(mix: dict) -> float:
    return CP.PropsSI("GAS_CONSTANT", mix_to_CP_string(mix)) / CP.PropsSI("M", mix_to_CP_string(mix))

def get_speed_of_sound(P, T, mix: dict) -> float:
    return (get_kappa(P, T, mix) * CP.PropsSI("GAS_CONSTANT", mix_to_CP_string(mix)) * T / CP.PropsSI("M", mix_to_CP_string(mix)))**0.5

def get_cp(P, T, mix: dict):
    return mass_mixer(mix, P, "P", T, "T", "CP0MASS")

def get_friction_factor(Re) -> float: #acc to McKeon et al. [1]
    #if Re < 300E3:
        #pass #error handling
    #if Re >= 18E6:
        #pass #error handling
    def f(x):
        z = 1.93 * math.log(Re * x**0.5, 10) - 0.537 - x**-0.5
        return z
    return fsolve(f, 10**-5)[0]

def get_reynolds_number(P, T, r, u, mix: dict) -> float:
    rho = get_rho(P, T, mix)
    #visc = CP.PropsSI("VISCOSITY", "P", P, "T", T, mix_to_CP_string(mix))
    visc =  get_viscosity_mix(P, T, mix)
    return rho * u * 2*r / visc

def get_viscosity(P, T, fluid):
    if fluid == "CO":
        return 7.3478e-4 * 0.1 #bei 2600K
    else:
        return CP.PropsSI("VISCOSITY", "T", T, "P", P, fluid)#######################################
    
    #Nach VDI, falls CoolProp nicht will
    dipole_moments = {"CO": 0.12,
                      "CO2": 0,
                      "H2O": 1.8,
                      "H2": 0,
                      "N2": 0} #[debye]
    try:
        return CP.PropsSI("VISCOSITY", "T", T, "P", P, fluid)
    except Exception:
        T_c = CP.PropsSI("TCRIT", fluid)
        P_c = CP.PropsSI("PCRIT", fluid)
        T_r = T / T_c
        mu_r = 52.46 * dipole_moments[fluid]**2 * P_c * 1e-5 * T_c**-2
        F=1
        if mu_r <= 0.022:
            F = 1
        elif mu_r <= 0.075:
            F = 1 + 30.55 * (0.292 - CP.PropsSI("Z", fluid))**1.72
        else:
            F = 1 + 30.55 * (0.292 - CP.PropsSI("Z", fluid))**1.72 * abs(0.96 + 0.1 * (T_r-0.7))
        zeta = 0.176 * T_c**(1/6) * CP.PropsSI("M", fluid)**-0.5 * P_c**-(2/3) * 1e-5
        return 1e-7 * F * (0.807 * T_r**0.618 - 0.357 * math.exp(-0.449*T_r) + 0.34 * math.exp(-4.058 * T_r) + 0.018) / zeta
    
    VDI_params = {"CO": (0.01384*10**-5, 0.74306*10**-7, -0.62996*10**-10, 0.03948*10**-12, -0.01032*10**-15),
                  "CO2": (-0.18024*10**-5, 0.65989*10**-7, -0.37108*10**-10, 0.01586*10**-12, -0.00300*10**-15)} #VDI_Waermeatlas S.381, S.450
    if fluid in VDI_params:
        a, b, c, d, e = VDI_params[fluid]
        return a + b*T + c*T**2 + d*T**3 + e*T**4
    else:
        return CP.PropsSI("VISCOSITY", "T", T, "P", P, fluid)

def get_viscosity_mix(P, T, mix: dict):
    #mfmix = get_mass_frac(mix)
    visc = 0
    for i in mix:
        denominator = 0
        for j in mix:
            F = ((1 + (get_viscosity(P, T, i) / get_viscosity(P, T, j))**0.5 * (CP.PropsSI("M", j) / CP.PropsSI("M", i))**0.25)**2
                 / (8 * (1 + CP.PropsSI("M", i) / CP.PropsSI("M", j)))**0.5)
            denominator += F * mix[j]
        visc += mix[i] * get_viscosity(P, T, i) / denominator
    return visc

def get_thermal_conductivity(P, T, fluid):
    if fluid == "CO":
        return 1.8386e-4 * 418.4 #bei 2600K
    else:
        return CP.PropsSI("CONDUCTIVITY", "P", P, "T", T, fluid)

def get_thermal_conductivity_mix(P, T, mix):
    mfmix = get_mass_frac(mix)
    cond = 0
    for i in mix:
        denominator = 0
        for j in mix:
            F = ((1 + (get_thermal_conductivity(P, T, i) / get_thermal_conductivity(P, T, j))**0.5 * (CP.PropsSI("M", j) / CP.PropsSI("M", i))**0.25)**2
                 / (8 * (1 + CP.PropsSI("M", i) / CP.PropsSI("M", j)))**0.5)
            denominator += F * mfmix[j]
        cond += mfmix[i] * get_thermal_conductivity(P, T, i) / denominator
    return cond

def get_prandtl_number(P, T, mix: dict) -> float:
    cp = get_cp(P, T, mix)
    visc = get_viscosity_mix(P, T, mix)
    k = get_thermal_conductivity_mix(P, T, mix)
    return cp * visc / k

def get_stanton_number(P, T, r, u, mix: dict): #acc to Friend and Metzner [2]
    fff = get_friction_factor(get_reynolds_number(P, T, r, u, mix))
    Pr = get_prandtl_number(P, T, mix)
    return (fff/8)/(1.2 + 11.8 * (fff/8)**0.5 * (Pr - 1) * (Pr)**(-1/3))

def get_heat_transfer_coefficient(P, T, r, u, mix: dict):
    St = get_stanton_number(P, T, r, u, mix)
    rho = get_rho(P, T, mix)
    cp = get_cp(P, T, mix)
    return St * u * rho * cp

def get_adiabatic_wall_temperature(P, T, u, recovery_factor, mix: dict):
    M = u / get_speed_of_sound(P, T, mix)
    return T * (1 + recovery_factor * (get_kappa(P, T, mix) - 1) / 2 * M**2)

def get_heat_of_vaporization(P, substance):
    return CP.PropsSI("H", "P", P, "Q", 1, substance) - CP.PropsSI("H", "P", P, "Q", 0, substance)
    
def prune_cc_products(cea_dict, cutoff=0.01, exclude=["*NO", "*O", "*OH", "*H"], pos=1):
    cp_dict = {}
    total = 0
    for key in cea_dict:
        if (cea_dict[key][pos] > cutoff) and not (key in exclude):
            cp_dict[key.strip('*')] = cea_dict[key][pos]
            total += cea_dict[key][pos]
    for key in cp_dict:
        cp_dict[key] = cp_dict[key]/total
    return cp_dict

def get_cone_area(r1, r2, h):
    return (r1 + r2) * math.pi * ((r1-r2)**2 + h**2)**0.5

def get_mass_frac(mix: dict):
    M = 0
    for i in mix:
        M += mix[i] * CP.PropsSI("M", i)
    mdict = {}
    for j in mix:
        mdict[j] = mix[j] * CP.PropsSI("M", j) / M
    return mdict

def newton_solver(funcs, jacobi, initguess, cutoff=0.01):
    xn = np.array(initguess)
    i=0
    while(True):
        i += 1
        a = jacobi(xn)
        b = funcs(xn)*-1
        try:
            z = np.linalg.solve(a, b)
        except np.linalg.LinAlgError:
            z = np.linalg.lstsq(a, b, rcond=None)[0]
        
        if np.any(abs(z) > cutoff):
            xn += z
            #print(z)
        else:
            #print(i)
            #print(abs(z))
            #print(abs(z)<cutoff)
            return xn