# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 18:01:37 2022

@author: Simon
"""
import math

#----------Durchmesser und Lochanzahl



#gegebene Werte für Fuel (Ethanol)

d_1 = 0.0011 #m Durchmesser der Öffnungen im Ringsegment für das Ethanol
m_1 = 0.069 #kg/s Massenstrom Fuel
rho_1 = 789.3 #kg/m^3 Dichte Fuel

A_1 = ((d_1/2)**2)*math.pi #m^2 Querschnittsfläche
v_1 = m_1/(rho_1*A_1) #m/s Geschwindigkeit des Fuels

#-----------Berechnung von A3 (Querschnittfläche der Öffnung für cooling)
#Flächenverhältnis

m_3 = 0.002624 #kg/s benötigter Massenstrom für filmcooling

d_2 = 0.003 #mm Durchmesser für Fuel+Oxidizer Mischung
m_2 = m_1-m_3 #kg/s Massenstrom der Mischung)
A_2 = (d_2/2)**2*math.pi #m^2 
v_2 = m_2/(rho_1*A_2) #m/s


c = m_3/m_2 #Verhältnis zwischen austretenden Massenströme

A_3 = A_2*c #m^2
v_3 = m_3/(rho_1*A_3)

#-------- Lochanzahl + Durchmesser für A3
i = 4 #Lochanzahl
d_3 = math.sqrt(A_3/(math.pi*i))*2 


print("Lochanzahl:", i)
print("Durchmesser der Bohrungen:", d_3*1000,"mm")