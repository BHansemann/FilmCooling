# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:14:00 2020

@author: Bernhard
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import bezier
import sympy
from sympy.solvers import solve

class Nozzle:
    '''Placeholder'''
    
    __l1, __l2, __l3, __l4, __l5, __l6 = 0 #needed for calculation
    __step, __x, __y = 0, 1, 2 #key for numpy array columns
    
    def __init__(self):
        self.r_c = 0 #chamber radius
        self.r_t = 0 #throat radius
        self.r_e = 0 #nozzle exit radius
        
        self.rf_cf = 0 #radius fraction of chamber fillet (to chamber radius)
        self.rf_con = 0 #radius fraction of converging throat section
        self.rf_div = 0 #radius fraction of diverging throat section
        
        self.l_c = 0 #length of chamber
        self.l_t = 0 #length to throat
        self.l_e = 0 #length to exit
        self.l_n = 0 #nozzle length
        
        self.lf_n = 0 #nozzle length fraction (to conical 30° nozzle)
        
        self.alpha_con = 0 #convergent half-angle
        self.alpha_divt = 0 #divergent half-angle at throat
        self.alpha_dive = 0 #divergent half-angle at exit
        #all in radians
        
        self.eta = 0 #nozzle efficiency
        self.epsilon #expansion ratio
        
        self.nozzle_type = "bell" #nozzle type (bell or conical) compare with str.casefold()
        self.chamber_type = "cylindrical" #chamber shape (cylindrical), in the future also spherical and near-spherical
        self.chamber = True #with or without chamber
        
        self.l = 0 #overall length
        self.l_t = 0 #length to throat
        self.l_e = 0 #length to exit
        self.l_n = 0 #nozzle length
        
        self.steps = 0 #number of steps
        self.resolution = 0 #number of steps per meter
        
        self.data = np.array(np.zeros((self.steps + 1, 3)))
    
    def __repr__(self):
        pass
    
    def __array__(self):
        pass
    
    def generate(self):
        if self.chamber:
            if (self.chamber_type.casefold() == "cylindrical".casefold() or
                self.chamber_type.casefold() == "c".casefold()):
                    __l1 = self.l_c - self.r_c * self.rf_cf * math.tan(self.alpha_con)
                    __l2 = self.r_c * self.rf_cf * math.sin(self.alpha_con)
            
            elif (self.chamber_type.casefold() == "spherical".casefold() or
                  self.chamber_type.casefold() == "s".casefold()):
                    pass
            elif (self.chamber_type.casefold() == "near spherical".casefold() or
                  self.chamber_type.casefold() == "ns".casefold()):
                    pass
            else:
                pass #raise exception
        else: #if no chamber is selected, l1 and l2 are automatically 0
            __l1 = 0
            __l2 = 0
        
        __l3 = (self.r_c * (1 - self.rf_cf * (1 - math.cos(self.alpha_con))) -
                self.r_t * (1 + self.rf_con * (1 - math.cos(self.alpha_con)))) / (math.tan(self.alpha_con))
        __l4 = self.r_t * self.rf_con * math.sin(self.alpha_con)
        __l5 = self.r_t * self.rf_div * math.sin(self.alpha_divt)
        lf = self.lf_n if (self.nozzle_type.casefold() == "bell".casefold() or self.nozzle_type.casefold() == "b".casefold()) else 1
        __l6 = (self.r_t * (self.rf_div * (math.cos(math.radians(15))**-1 - 1) - 1) + self.r_e) * lf / math.tan(math.radians(15)) - __l5
        
        self.l = __l1 + __l2 + __l3 + __l4 + __l5 + __l6
    
    def export_as(self, ex_type: str):
        pass
    
    def export_csv(self):
        pass
    
    def export_np(self):
        pass
    
    def export_pandas(self):
        pass
    
    def export_graph(self):
        pass