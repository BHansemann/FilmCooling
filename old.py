# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:34:31 2020

@author: Bernhard
"""

# adiabatic approximation of thermodynamic properties inside the engine
import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

# inputs
p_ch = 1000000       # chamber pressure in Pa
p_e = 100000         # exit pressure in Pa
n_ps = 1000          # number of pressure steps
T_ch = 2500         # chamber temperature in Kelvin
F_thrust = 500      # Thrust in N
eta_nozzle = 0.983  # Efficiency of Nozzle
x_n2 = 0.34269      # Mole fraction of exhaust products, via NASA CEA
x_co = 0.23055
x_h2 = 0.21104
x_h2o = 0.17950
x_co2 = 0.03140

components_string = "HEOS::N2[{}]&CO[{}]&H2[{}]&H2O[{}]&CO2[{}]".format(x_n2, x_co, x_h2, x_h2o, x_co2)
thermoprops = ['P', 'T', 'rho', 'kappa', 'h', 'u', 'c', 'M', 'r']
M_mix = CP.PropsSI('M', components_string)
R_mix = CP.PropsSI('GAS_CONSTANT', components_string)
mixture = pd.DataFrame(np.zeros((5, 2)), index=['N2', 'CO', 'H2', 'H2O', 'CO2'], columns=['molar', 'mass'])

mixture.loc['N2', 'molar'] = x_n2
mixture.loc['CO', 'molar'] = x_co
mixture.loc['H2', 'molar'] = x_h2
mixture.loc['H2O', 'molar'] = x_h2o
mixture.loc['CO2', 'molar'] = x_co2


def masscalc(row):
    return row['molar'] * CP.PropsSI('M', row.name) / M_mix


mixture['mass'] = mixture.apply(masscalc, axis=1)

print(mixture)


def mixer(prop, ip1, input1, ip2, input2, mixture=mixture, mode='molar'):
    value = 0
    for index, row in mixture.iterrows():
        value += row[mode]*CP.PropsSI(prop, ip1, input1, ip2, input2, row.name)
    return value


p_step = (p_ch-p_e)/(n_ps-1)
data = pd.DataFrame(np.zeros((n_ps, thermoprops.__len__())), columns=thermoprops)

data['T'][0] = T_ch

for index, row in data.iterrows():
    row['P'] = p_ch-row.name*p_step
    if index != 0:
        row['T'] = data.loc[index-1, 'T'] * (row['P'] / data.loc[index-1, 'P'])**((data.loc[index-1, 'kappa'] - 1) / data.loc[index-1, 'kappa'])
    row['kappa'] = mixer('CP0MOLAR', 'T', row['T'], 'P', row['P']) / mixer('CVMOLAR', 'T', row['T'], 'P', row['T'])
    row['h'] = mixer('H', 'P', row['P'], 'T', row['T'], mode='mass')
    if index != 0:
        row['u'] = (2 * (mixer('H', 'P', data.loc[index-1, 'P'], 'T', data.loc[index-1, 'T'], mode='mass') - row['h'] + (data.loc[index-1, 'u']**2)/2))**0.5
    row['rho'] = row['P'] / (R_mix * row['T'])
    row['c'] = (row['kappa'] * R_mix * row['T'])**0.5
    row['M'] = row['u'] / row['c']

m_dot = F_thrust / (data['u'].iloc[-1]*eta_nozzle)

for index, row in data.iterrows():
    if index != 0:
        row['r'] = (m_dot / (math.pi * row['u'] * row['rho']))**0.5

print(data)
data.plot(subplots=True)
plt.show()