# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 12:45:02 2020

@author: Bernhard
"""

import math
import scipy.constants as const
import CoolProp.CoolProp as CP

def isentrop_press_temp(kappa, press, temp1, temp2):
    return press * (temp2 / temp1)**(kappa / (kappa - 1))

def isentrop_press_vol(kappa, press, vol1, vol2):
    return press * (vol1 / vol2)**kappa

def isentrop_temp_press(kappa, temp, press1, press2):
    return temp * (press2 / press1)**((kappa - 1) / kappa)

def isentrop_temp_vol(kappa, temp, vol1, vol2):
    return temp * (vol2 / vol1)**(1-kappa)

def isentrop_vol_temp(kappa, vol, temp1, temp2):
    return vol * (temp2 / temp1)**(1 / (1-kappa))

def isentrop_vol_press(kappa, vol, press1, press2):
    return vol * (press1 / press2)**(1 / kappa)

def mix_to_CP_string(mix: dict):
    out = ""
    for key in mix:
        out += "{}[{}]&".format(key, mix[key])
    return out.strip("&")

def mass_frac(mix: dict):
    massmix = {}
    molarmass_mix = CP.PropsSI('M', mix_to_CP_string(mix))
    for key in mix:
        massmix[key] = mix[key] * CP.PropsSI('M', key) / molarmass_mix
    return massmix

def molar_mixer(mix: dict, input1, type1, input2, type2, output):
    val = 0
    for key in mix:
        val += mix[key] * CP.PropsSI(output, type1, input1, type2, input2)
    return val

def mass_mixer(mix: dict, input1, type1, input2, type2, output):
    val=0
    massmix = mass_frac(mix)
    for key in massmix:
        val += massmix[key] * CP.PropsSI(output, type1, input1, type2, input2)
    return val