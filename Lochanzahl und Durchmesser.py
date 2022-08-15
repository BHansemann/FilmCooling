# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 18:01:37 2022

@author: Simon Luk
"""
import math

#----------Durchmesser und Lochanzahl



#gegebene Werte fuer Fuel (Ethanol)

d_1 = 0.0011 #m Durchmesser der OEffnungen im Ringsegment fuer das Ethanol
m_1 = 0.069 #kg/s Massenstrom Fuel
rho_1 = 789.3 #kg/m^3 Dichte Fuel

A_1 = ((d_1/2)**2)*math.pi #m^2 Querschnittsflaeche
v_1 = m_1/(rho_1*A_1) #m/s Geschwindigkeit des Fuels

#-----------Berechnung von A3 (Querschnittflaeche der OEffnung fuer cooling)
#Flaechenverhaeltnis

m_3 = 0.0037736364768785288 #kg/s benoetigter Massenstrom fuer filmcooling


d_2 = 0.003 #mm Durchmesser fuer Fuel+Oxidizer Mischung
m_2 = m_1-m_3 #kg/s Massenstrom der Mischung)
A_2 = (d_2/2)**2*math.pi #m^2 
v_2 = m_2/(rho_1*A_2) #m/s


c = m_3/m_2 #Verhaeltnis zwischen austretenden Massenstroeme

A_3 = A_2*c #m^2
v_3 = m_3/(rho_1*A_3)

#-------- Lochanzahl + Durchmesser fuer A3
i = 6 #Lochanzahl
d_3 = math.sqrt(A_3/(math.pi*i))*2 


print("Lochanzahl:", i)
print("Durchmesser der Bohrungen:", d_3*1000,"mm")