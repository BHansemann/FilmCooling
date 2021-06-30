# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 12:51:47 2021

@author: Bernhard
"""

import numpy as np
from CoolProp.CoolProp import PropsSI
from rocketcea.cea_obj import CEA_Obj
import scipy.constants as const
import thermodynamic_toolbox as thermo
from tqdm import tqdm
import math

__author__ = "Bernhard Hansemann"
__date__ = "2021"
__contact__ = "bernhard.hansemann@spaceteam.at"
__version__ = "1.0.0"

psi = lambda x : x / const.psi
Pa = lambda x : x * const.psi

class AdiabaticEngine:
    
    step, P, T, rho, kappa, h, u, c, M, r = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    
    def __init__(self, fuel, oxidizer, P_cc, P_e, F_t, ofr, r_cc, steps, eta_n=1, stepfactor=10):
        self.fuel = fuel
        self.oxidizer = oxidizer
        self.P_cc = P_cc
        self.P_e = P_e
        self.F_t = F_t
        self.ofr = ofr
        self.r_cc = r_cc
        self.steps = steps
        self.eta_n = eta_n
        self.p_steps = stepfactor*steps
        self.steps = steps
        
        self.p_data = np.zeros((10, self.p_steps+1))
    
    def __repr__(self):
        pass
    
    def __str__(self):
        pass
    
    def generate(self, disable_tqdm=False):
        self.cea = CEA_Obj(fuelName=self.fuel, oxName=self.oxidizer)
        self.mix_cc = self.prune_combustion_products(self.cea.get_SpeciesMoleFractions(Pc=psi(self.P_cc), MR=self.ofr)[1])
        self.p_data[self.step,:] = np.arange(0, self.p_steps+1)
        self.p_data[self.P,:] = self.P_cc - ((self.P_cc - self.P_e)/self.p_steps) * self.p_data[self.step,:]
        self.p_data[self.T,0] = self.cea.get_Temperatures(Pc=psi(self.P_cc), MR=self.ofr)[0] * (5/9)
        self.p_data[self.rho,0] = thermo.get_rho(self.p_data[self.P,0], self.p_data[self.T,0], self.mix_cc)
        self.p_data[self.kappa,0] = thermo.get_kappa(self.p_data[self.P,0], self.p_data[self.T,0], self.mix_cc)
        self.p_data[self.h,0] = thermo.mass_mixer(self.mix_cc, self.p_data[self.P,0], 'P', self.p_data[self.T,0], 'T', 'H')
        self.p_data[self.u,0] = 0
        self.p_data[self.c,0] = thermo.get_speed_of_sound(self.p_data[self.P,0], self.p_data[self.T,0], self.mix_cc)
        self.p_data[self.M,0] = self.p_data[self.u,0] / self.p_data[self.c,0]
        self.p_data[self.r,0] = self.r_cc
        
        for i in tqdm(self.p_data[self.step,1:],desc="Generating slices", disable=disable_tqdm):
            i=int(i)
            self.p_data[1:,i] = self.__next_step(self.p_data[:,i-1], self.p_data[self.P,i])
        
        self.m_dot = self.F_t / self.p_data[self.u,-1]
        
        self.p_data[self.u,0] = self.m_dot / (self.p_data[self.rho,0] * self.p_data[self.r,0]**2 * math.pi)
        
        for i in tqdm(self.p_data[self.step,1:],desc="Calculating velocities", disable=disable_tqdm):
            i = int(i)
            self.p_data[self.u,i] = ((self.p_data[self.h,i-1] - self.p_data[self.h,i]) * 2 + self.p_data[self.u,i-1]**2)**0.5
            self.p_data[self.M,i] = self.p_data[self.u,i] / self.p_data[self.c,i]
        self.p_data[self.r,:] = np.vectorize(lambda u, rho: (self.m_dot / (math.pi * u * rho))**0.5)(self.p_data[self.u,:], self.p_data[self.rho,:])
        
        delta=1000
        self.t_step=0
        for i in self.p_data[self.step,:]:
            i = int(i)
            if abs(self.p_data[self.M,i] - 1) < delta:
                delta = abs(self.p_data[self.M,i] - 1)
                self.t_step = i
        print(self.t_step)
        print(delta)
        
        self.r_e = self.p_data[self.r,-1]
        self.r_t = self.p_data[self.r,self.t_step]
    
    def __next_step(self, prev, P_now):
        T_now = thermo.isentrop_temp_press(prev[self.kappa], prev[self.T], prev[self.P], P_now)
        rho_now = thermo.get_rho(P_now, T_now, self.mix_cc)
        kappa_now = thermo.get_kappa(P_now, T_now, self.mix_cc)
        h_now = thermo.mass_mixer(self.mix_cc, P_now, 'P', T_now, 'T', 'H')
        u_now = ((prev[self.h] - h_now) * 2 + prev[self.u]**2)**0.5
        c_now = thermo.get_speed_of_sound(P_now, T_now, self.mix_cc)
        M_now = u_now / c_now
        
        return np.array([P_now, T_now, rho_now, kappa_now, h_now, u_now, c_now, M_now, 0])
    
    def prune_combustion_products(self, cea_dict, cutoff=0.01, exclude=["*NO", "*O", "*OH", "*H"], pos=1):
        cp_dict = {}
        total = 0
        for key in cea_dict:
            if (cea_dict[key][pos] > cutoff) and not (key in exclude):
                cp_dict[key.strip('*')] = cea_dict[key][pos]
                total += cea_dict[key][pos]
        for key in cp_dict:
            cp_dict[key] = cp_dict[key]/total
        return cp_dict
    
    def matcher(self, nozzle):
        pass
    
    def quickset(self):
        pass
    