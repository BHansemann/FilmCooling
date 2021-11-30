# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 17:25:04 2021

@author: Bernhard
"""

#Imports
import numpy as np
from CoolProp.CoolProp import PropsSI
import math
from rocketcea.cea_obj import CEA_Obj
import scipy.constants as constants

import nozzle_generator
import thermodynamic_toolbox as thermo
from engine_adiabatic import engine_adiabatic

#Mini-Funktionen
area = lambda r: r**2 * math.pi
psi = lambda x : x / constants.psi
Pa = lambda x : x * constants.psi

#Header für Matrizen
columns = 13
step, x, r, P, T, rho, kappa, h, u, c, M, Q, m_dot_c = tuple(range(0,columns))

#Eingangsvariablen
P_c = 1e6 #10 bar
P_e = 101325 #Normaldruck

F_th = 500 #Schub

Fuel = "C2H5OH"
Fuel_name = "Ethanol"
Oxidizer = "N2O"

OFR = 2.5 #Oxidizer-Fuel Ratio

r_cc = 0.02 #Radius der Brennkammer

rec_fac = 1 #recovery Faktor für adiabatic wall temperature


#Initialisierung
cea = CEA_Obj(fuelName=Fuel, oxName=Oxidizer)
mix = thermo.prune_cc_products(cea.get_SpeciesMoleFractions(Pc=psi(P_c), MR=OFR)[1])
Rs = thermo.get_specific_gas_constant(mix)
init_eng = engine_adiabatic(P_c, P_e, 1000, 0, F_th, 0, Oxidizer, Fuel, OFR, r_cc)

r_t = float(init_eng[9].min()) #throat Radius
r_e = float(init_eng[9,-1]) #exit Radius

#Werte aus Tabellen abgelesen, müssen eventuell korrigiert werden.
alpha_con = 35
alpha_divt = 21
alpha_dive = 11

#Frei angenommen
l_c = 0.15 #Länge der Brennkammer
l_fc = 0.14 #Länge des Kühlfilms

nozzle = nozzle_generator.Nozzle()
nozzle.quickset(r_cc, r_t, r_e, alpha_con, alpha_divt, alpha_dive, 1.5, 0.382, l_c, 0, "bell", 0.9, 1000, "m")
nozzle.generate() #Kontur generieren

data = np.zeros((len(nozzle.export_np()), columns)) #Hauptmatrix generieren
data[:,0:3] = nozzle.export_np()[:] #und mit Kontur füllen

m_dot = F_th/init_eng[6,-1] #erste Annahme
F_r = 0 #for now

while True:
    m_low = m_dot/2 #untere Schranke
    m_high = m_dot*1.5 #obere Schranke
    #Mit Eingangsdaten befüllen
    data[0,P] = P_c
    data[0,T] = init_eng[2,0]
    data[0,rho] = thermo.get_rho(data[0,P], data[0,T], mix)
    data[0,kappa] = thermo.get_kappa(data[0,P], data[0,T], mix)
    data[0,h] = thermo.mass_mixer(mix, data[0,P], 'P', data[0,T], 'T', 'H')
    data[0,c] = thermo.get_speed_of_sound(data[0,P], data[0,T], mix)
    data[0,u] = m_dot / (area(r_cc) * data[0,rho])
    data[0,M] = data[0,u] / data[0,c]
    data[0,Q] = -1 * thermo.get_cone_area(data[0,r], data[0,r], data[1,x]) * thermo.get_heat_transfer_coefficient(data[0,P], data[0,T], data[0,r], data[0,u], mix) * (data[0,T] - thermo.get_adiabatic_wall_temperature(data[0,P], data[0,T], data[0,u], rec_fac, mix))
    data[0,m_dot_c] = data[0,Q] / thermo.get_heat_of_vaporization(data[0,P], "Ethanol")
    
    try:
        for i in range(1, nozzle.get_n_t()):
            s1 = data[i-1] #Slice 1
            s2 = data[i] #Slice 2
            
            in_cool = data[i][x] < l_fc #Liquid film vorhanden?
            h_ev = thermo.get_heat_of_vaporization(s1[P], "Ethanol") #Verdampfungsenthalpie von Ethanol
            mdotf = (s1[Q] / h_ev) if in_cool else 0 #verdampfendes Ethanol
            h_f = PropsSI('H','P',s1[P],'T', s1[T],'Ethanol') #spez Enthalpie vom verdampften Ethanol
            s2[h] = s1[h] #workaround
            
            def funcmatrix(vec): #Funktionsvektor (Massenerhaltung, Energieerhaltung, Impulserhaltung)
                P_2, T_2, u_2 = tuple(vec)
                return np.array([P_2 * u_2 * area(s1[r]) / (T_2 * Rs) - s1[rho] * s1[u] * area(s1[r]) - s1[Q] / h_ev,
                                 P_2 * u_2 * area(s1[r]) * s2[h] / (T_2 * Rs) - s1[rho] * s1[u] * area(s1[r]) * s1[h] - mdotf * h_f + s1[Q],
                                 P_2 * u_2**2 * area(s1[r]) / (T_2 * Rs) + P_2 * area(s1[r]) - s1[rho] * s1[u]**2 * area(s1[r]) - s1[P] * area(s1[r]) + F_r])
            
            def jacobimatrix(vec): #Jacobimatrix des Funktionsvektors
                P_2, T_2, u_2 = tuple(vec)
                alpha = u_2 * area(s1[r]) / (T_2 * Rs) #Gemeinsamer Faktor, zur Vereinfachung
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
            
            P_init = s1[P] * 11
            T_init = s1[T] * 11
            u_init = s1[u] * 11
            P_2, T_2, u_2 = thermo.newton_solver(funcmatrix, jacobimatrix, (P_init, T_init, u_init))
            
            s2[P] = P_2
            s2[T] = T_2
            s2[u] = u_2
            
            s2[rho] = thermo.get_rho(s2[P], s2[T], mix)
            s2[kappa] = thermo.get_kappa(s2[P], s2[T], mix)
            s2[h] = thermo.mass_mixer(mix, s2[P], "P", s2[T], "T", "H")
            s2[c] = thermo.get_speed_of_sound(s2[P], s2[T], mix)
            s2[M] = s2[u] / s2[c]
            if s2[M] > 1:
                raise ValueError(f"M > 1 at step = {s2[step]}")
            s2[Q] = -1 * thermo.get_heat_transfer_coefficient(s2[P], s2[T], s2[r], s2[u], mix) if in_cool else 0
            s2[m_dot_c] = s2[Q] / h_ev
            
            data[i] = s2
        
        
    except np.linalg.LinAlgError as ex:
        print(ex)
        m_dot = m_low
        print(m_dot)
    except ValueError:
        pass
