# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 12:45:02 2020

@author: Bernhard
"""

import math
import scipy.constants as const
from scipy.optimize import fsolve
import CoolProp.CoolProp as CP

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

def get_speed_of_sound(P, T, mix: dict) -> float:
    return (get_kappa(P, T, mix) * CP.PropsSI("GAS_CONSTANT", mix_to_CP_string(mix)) * T / CP.PropsSI("M", mix_to_CP_string(mix)))**0.5

def get_cp(P, T, mix: dict):
    return mass_mixer(mix, P, "P", T, "T", "CP0MASS")

def get_friction_factor(Re) -> float: #acc to McKeon et al. [1]
    if Re < 300E3:
        pass #error handling
    if Re >= 18E6:
        pass #error handling
    def f(x):
        z = 1.93 * math.log(Re * x**0.5, 10) - 0.537 - x**-0.5
        return z
    return fsolve(f, 10**-5)

def get_reynolds_number(P, T, r, u, mix: dict) -> float:
    rho = get_rho(P, T, mix)
    visc = CP.PropsSI("VISCOSITY", "P", P, "T", T, mix_to_CP_string(mix))
    return rho * u * 2*r / visc

def get_prandtl_number(P, T, mix: dict) -> float:
    cp = get_cp(P, T, mix)
    visc = CP.PropsSI("VISCOSITY", "P", P, "T", T, mix_to_CP_string(mix))
    k = CP.PropsSI("CONDUCTIVITY", "P", P, "T", T, mix_to_CP_string(mix))
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
    