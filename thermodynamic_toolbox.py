# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 12:45:02 2020

@author: Bernhard Hansemann
"""

import math
from scipy.optimize import fsolve
import CoolProp.CoolProp as CP

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

def get_speed_of_sound(P, T, mix: dict) -> float:
    return (get_kappa(P, T, mix) * CP.PropsSI("GAS_CONSTANT", mix_to_CP_string(mix)) * T / CP.PropsSI("M", mix_to_CP_string(mix)))**0.5

def get_cp(P, T, mix: dict):
    return mass_mixer(mix, P, "P", T, "T", "CP0MASS")

def get_friction_factor(Re) -> float: #acc to McKeon et al. [1]
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
        return CP.PropsSI("VISCOSITY", "T", T, "P", P, fluid)

def get_viscosity_mix(P, T, mix: dict):
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

def get_mass_frac(mix: dict):
    M = 0
    for i in mix:
        M += mix[i] * CP.PropsSI("M", i)
    mdict = {}
    for j in mix:
        mdict[j] = mix[j] * CP.PropsSI("M", j) / M
    return mdict