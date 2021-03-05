# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 16:42:05 2021

@author: Berni
"""

import sympy
import math

P1, P2, T1, T2, Rs, r1, r2, u1, u2, Qc, hcc, Fr, h1, h2 = sympy.symbols("P1 P2 T1 T2 Rs r1 r2 u1 u2 Qc hcc Fr h1 h2")

f = (P1 * r1**2 * math.pi * u1) / (Rs * T1) + Qc / hcc - (P2 * r2**2 * math.pi * u2) / (Rs * T2)
g = (P1 * r1**2 * math.pi * u1 * h1) / (Rs * T1) + Qc * hcc / hcc - Qc - (P2 * r2**2 * math.pi * u2 * h2) / (Rs * T2)
h = P1/(Rs * T1) * u1**2 * r1**2 * math.pi + P1 * r1**2 * math.pi - Fr - P2/(Rs*T2) * u2**2 * r2**2 * math.pi - P2 * r2**2 * math.pi

