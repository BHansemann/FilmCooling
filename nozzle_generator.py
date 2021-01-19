# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:14:00 2020

@author: Bernhard
"""

import numpy as np
import math
import bezier
import sympy
from sympy.solvers import solve

__author__ = "Bernhard Hansemann"
__date__ = "Jan 19, 2021"
__contact__ = "bernhard.hansemann@spaceteam.at"
__version__ = "1.0.0"

class Nozzle:
    '''Placeholder'''
    
    __l1, __l2, __l3, __l4, __l5, __l6 = 0 ,0, 0, 0 ,0 ,0 #needed for calculation
    __step, __cx, __cy = 0, 1, 2 #key for numpy array columns
    
    def __init__(self):
        '''
        Constructor. Declaration of Variables.

        Returns
        -------
        None.

        '''
        self.r_c = 0 #chamber radius
        self.r_t = 0 #throat radius
        self.r_e = 0 #nozzle exit radius
        
        self.rf_cf = 0 #radius fraction of chamber fillet (to chamber radius)
        self.rf_con = 0 #radius fraction of converging throat section
        self.rf_div = 0 #radius fraction of diverging throat section
        
        self.l_c = 0 #length of chamber
        self.l_t = 0 #length to throat
        self.l_n = 0 #nozzle length
        self.l = 0 #overall length
        
        self.n_t = 0 #step number of throat
        
        self.lf_n = 0 #nozzle length fraction (to conical 30° nozzle)
        
        self.alpha_con = 0 #convergent half-angle
        self.alpha_divt = 0 #divergent half-angle at throat
        self.alpha_dive = 0 #divergent half-angle at exit
        #all in radians
        
        self.eta = 0 #nozzle efficiency
        self.epsilon = 0 #expansion ratio
        
        self.nozzle_type = "b" #nozzle type (bell or conical) compare with str.casefold()
        self.chamber_type = "cylindrical" #chamber shape (cylindrical), in the future also spherical and near-spherical
        self.chamber = True #with or without chamber
        
        self.steps = 10000 #number of steps
        self.resolution = 0 #number of steps per meter
        
        self.data = np.array(np.zeros((self.steps + 1, 3))) #representation of the data as a numpy array
        self.data[:,self.__step] = np.arange(0, self.steps + 1)
    
    def __repr__(self):
        csvstring = ""
        for i in self.data[:,self.__step]:
            i=int(i)
            csvstring = csvstring + ','.join(['%.5f' % num for num in self.data[i,:]]) + "\n"
        return csvstring
    
    def __str__(self):
        return self.get_nozzle_parameters()
    
    def __array__(self):
        return self.data
    
    def quickset(self, r_c, r_t, r_e, alpha_con, alpha_divt, alpha_dive, rf_con, rf_div, l_c, rf_cf, nozzle_type, lf_n, steps, lunit, chamber=True):
        '''
        Sets all necessary parameters at once.

        Parameters
        ----------
        r_c : Float
            Chamber radius.
        r_t : Float
            Throat radius.
        r_e : Float
            Exit radius.
        alpha_con : Float
            Converging angle [degrees].
        alpha_divt : Float
            Diverging angle at throat [degrees].
        alpha_dive : Float
            Diverging angle at exit [degrees].
        rf_con : Float
            Radius fraction at throat upstream.
        rf_div : Float
            Radius fraction at throat downstream.
        l_c : Float
            Length of chamber.
        rf_cf : Float
            Chamber-nozzle fillet radius relative to the chamber radius.
        nozzle_type : String
            Type of nozzle ("c" or "b").
        lf_n : Float
            Length fraction of the bell nozzle.
        steps : Integer
            Number of steps.
        lunit : String
            Unit length (e.g. "mm").
        chamber : Boolean, optional
            Should a chamber be generated? The default is True.

        Returns
        -------
        None.

        '''
        self.set_r_c(r_c, unit=lunit)
        self.set_r_t(r_t, unit=lunit)
        self.set_r_e(r_e, unit=lunit)
        self.set_alpha_con(alpha_con)
        self.set_alpha_divt(alpha_divt)
        self.set_alpha_dive(alpha_dive)
        self.set_rf_con(rf_con)
        self.set_rf_div(rf_div)
        self.set_l_c(l_c, unit=lunit)
        self.set_rf_cf(rf_cf)
        self.set_chamber(chamber)
        self.set_chamber_type("cylindrical")
        self.set_nozzle_type(nozzle_type)
        self.set_lf_n(lf_n)
        self.set_steps(steps)
    
    def generate(self):
        '''
        Updates the data. Should be called after new parameters have been set.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        
        #check what chamber type is selected, update relevant lengths and set draw function (tbd)
        if self.chamber:
            if (self.chamber_type.casefold() == "cylindrical".casefold() or
                self.chamber_type.casefold() == "c".casefold()):
                    self.__l1 = self.l_c - self.r_c * self.rf_cf * math.tan(self.alpha_con)
                    self.__l2 = self.r_c * self.rf_cf * math.sin(self.alpha_con)
            
            elif (self.chamber_type.casefold() == "spherical".casefold() or
                  self.chamber_type.casefold() == "s".casefold()):
                    pass
            elif (self.chamber_type.casefold() == "near spherical".casefold() or
                  self.chamber_type.casefold() == "ns".casefold()):
                    pass
            else:
                raise ValueError("Chamber type of \"" + self.chamber_type + "\" invalid!")
        else: #if no chamber is selected, l1 and l2 are automatically 0
            self.__l1 = 0
            self.__l2 = 0
        
        #setting other lengths relevant for nozzle generation
        self.__l3 = (self.r_c * (1 - self.rf_cf * (1 - math.cos(self.alpha_con))) -
                self.r_t * (1 + self.rf_con * (1 - math.cos(self.alpha_con)))) / (math.tan(self.alpha_con))
        self.__l4 = self.r_t * self.rf_con * math.sin(self.alpha_con)
        self.__l5 = self.r_t * self.rf_div * math.sin(self.alpha_divt)
        lf = self.lf_n if (self.nozzle_type.casefold() == "bell".casefold() or self.nozzle_type.casefold() == "b".casefold()) else 1
        self.__l6 = (self.r_t * (self.rf_div * (math.cos(math.radians(15))**-1 - 1) - 1) + self.r_e) * lf / math.tan(math.radians(15)) - self.__l5
        
        #calculating lengths
        self.l = self.__l1 + self.__l2 + self.__l3 + self.__l4 + self.__l5 + self.__l6
        self.l_n = self.l - self.l_c
        self.l_t = self.__l1 + self.__l2 + self.__l3 + self.__l4
        self.n_t = int(self.l_t * (self.steps / self.l))
        
        #calculating and setting number of steps/resolution
        if self.steps == 0 and self.resolution == 0:
            raise ValueError("No resolution or number of steps selected.")
        elif self.steps == 0:
            self.steps = int(self.resolution * self.l)
            self.resolution = self.steps / self.l
        elif self.resolution == 0:
            self.resolution = self.steps / self.l
        else:
            self.resolution == self.steps / self.l
        
        #setting up the numpy array with step number and x length
        self.data = np.array(np.zeros((self.steps + 1, 3)))
        self.data[:,self.__step] = np.arange(0, self.steps + 1)
        self.data[:,self.__cx] = self.data[:,self.__step] * (self.l / self.steps)
        
        #drawing bezier curve for bell nozzle geometry
        x_0 = self.l - self.__l6
        y_0 = self.r_t * (1 + self.rf_div * (1 - math.cos(self.alpha_divt)))
        x_2 = self.l
        y_2 = self.r_e
        x_1 = (y_2 - y_0 + math.tan(self.alpha_divt) * x_0 - math.tan(self.alpha_dive) * x_2) / (math.tan(self.alpha_divt) - math.tan(self.alpha_dive))
        y_1 = math.tan(self.alpha_divt) * x_1 + y_0 - math.tan(self.alpha_divt) * x_0
        bell = bezier.Curve(np.asfortranarray([[x_0, x_1, x_2], [y_0, y_1, y_2]]), degree=2)
        y = sympy.Symbol("y")
        x = sympy.Symbol("x")
        bellexpr = solve(bell.implicitize(), y)
        bellfunc = [sympy.lambdify(x, bellexpr[0], "numpy"), sympy.lambdify(x, bellexpr[1], "numpy")]
        
        #defining the complete function the chamber-nozzle contour follows
        def get_y(x):
            if x <= self.__l1:
                x1 = x
                return self.r_c
            elif self.__l1 < x <= self.__l1 + self.__l2:
                x2 = x - self.__l1
                return self.r_c * (1 - self.rf_cf) + ((self.r_c * self.rf_cf)**2 - x2**2)**0.5
            elif self.__l1 + self.__l2 < x <= self.__l1 + self.__l2 + self.__l3:
                x3 = x - self.__l1 - self.__l2
                return self.r_c * (1 - self.rf_cf * (1 - math.cos(self.alpha_con))) - x3 * math.tan(self.alpha_con)
            elif self.__l1 + self.__l2 + self.__l3 < x <= self.__l1 + self.__l2 + self.__l3 + self.__l4:
                x4 = x - self.__l1 - self.__l2 - self.__l3
                return self.r_t + self.r_t * self.rf_con - ((self.r_t * self.rf_con)**2 - (x4 - self.__l4)**2)**0.5
            elif self.__l1 + self.__l2 + self.__l3 + self.__l4 < x <= self.__l1 + self.__l2 + self.__l3 + self.__l4 + self.__l5:
                x5 = x - self.__l1 - self.__l2 - self.__l3 - self.__l4
                return self.r_t + self.r_t * self.rf_div - ((self.r_t * self.rf_div)**2 - x5**2)**0.5
            elif self.__l1 + self.__l2 + self.__l3 + self.__l4 + self.__l5 < x <= self.l:
                x6 = x - self.__l1 - self.__l2 - self.__l3 - self.__l4 - self.__l5
                if (self.nozzle_type.casefold() == "bell".casefold() or
                    self.nozzle_type.casefold() == "b".casefold()):
                    if(bellfunc[0](x) > bellfunc[1](x)):
                        return bellfunc[0](x)
                    else:
                        return bellfunc[1](x)
                elif (self.nozzle_type.casefold() == "conical".casefold() or
                    self.nozzle_type.casefold() == "c".casefold()):
                    return self.r_t * (1 + self.rf_div * (1 - math.cos(self.alpha_divt))) + x6 * math.tan(self.alpha_divt)
                else:
                    raise ValueError("Nozzle type of\"" + self.nozzle_type + "\" invalid!")
        
        #generating the contour
        self.data[:,self.__cy] = np.vectorize(get_y)(self.data[:,self.__cx])
        
    def sample(self):
        '''
        Sets and generates a sample nozzle and chamber.

        Returns
        -------
        numpy.ndarray
            The data array of the sample nozzle.

        '''
        self.quickset(30, 20, 40, 25, 30, 10, 1.5, 0.382, 80, 0.1, "bell", 0.9, 1000, "mm")
        self.generate()
        return self.data
    
    def set_r_c(self, r_c, unit="m"):
        '''
        Sets the radius of the combustion chamber.

        Parameters
        ----------
        r_c : Float
            Chamber radius.
        unit : String, optional
            Unit of length. The default is "m".

        Returns
        -------
        None.

        '''
        #get conversion factor from unit
        factor = 1
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            raise ValueError("Unknown length unit: " + unit)
        if not(type(r_c) == float or type(r_c) == int):
            raise TypeError("Chamber radius must be a Number.")
        if r_c <= 0:
            raise ValueError("Chamber radius must be greater than 0.")
            
        self.r_c = r_c * factor
        
    def set_r_t(self, r_t, unit="m", fixed_exit=True):
        '''
        Sets throat radius of the nozzle.

        Parameters
        ----------
        r_t : Float
            Throat radius.
        unit : String, optional
            Unit of length. The default is "m".
        fixed_exit : Boolean, optional
            Change epsilon or exit radius? The default is True (change epsilon)

        Returns
        -------
        None.

        '''
        #get conversion factor from unit
        factor = 1
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            raise ValueError("Unknown length unit: " + unit)
        if not(type(r_t) == float or type(r_t) == int):
            raise TypeError("Throat radius must be a Number.")
        if r_t <= 0:
            raise ValueError("Throat radius must be greater than 0.")
            
        self.r_t = r_t * factor
        #chane r_e or epsilon
        if fixed_exit:
            self.epsilon = self.r_e**2 / self.r_t**2
        else:
            self.r_e = (self.epsilon * self.r_t**2)**0.5
        
    def set_r_e(self, r_e, unit="m", fixed_throat=True):
        '''
        Sets exit radius of the nozzle.

        Parameters
        ----------
        r_e : Float
            Exit radius.
        unit : String, optional
            Unit of length. The default is "m".
        fixed_throat : Boolean, optional
            Change epsilon or throat Radius? The default is True. (change epsilon)

        Returns
        -------
        None.

        '''
        #get conversion factor from unit
        factor = 1
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            raise ValueError("Unknown length unit: " + unit)
        if not(type(r_e) == float or type(r_e) == int):
            raise TypeError("Exit radius must be a Number.")
        if r_e <= 0:
            raise ValueError("Exit radius must be greater than 0.")
            
        self.r_e = r_e * factor
        #change r_t or epsilon
        if fixed_throat:
            self.epsilon = self.r_e**2 / self.r_t**2
        else:
            self.r_t = (self.r_e**2 / self.epsilon)**0.5
        
    def set_rf_cf(self, rf_cf):
        '''
        Sets the chamber-nozzle fillet radius as a fraction of the chamber radius.

        Parameters
        ----------
        rf_cf : Float
            Fraction of the chamber radius.

        Returns
        -------
        None.
        
        '''
        if not(type(rf_cf) == float or type(rf_cf) == int):
            raise TypeError("Chamber-nozzle fillet ratio must be a Number.")
        if rf_cf <= 0:
            raise ValueError("Chamber-nozzle fillet ratio must be greater than 0.")
            
        self.rf_cf = rf_cf
        
    def set_rf_con(self, rf_con):
        '''
        Sets the radius of the throat upstream as a fraction of the throat radius.
        Recommendation for bell nozzle: 1.5
        Recommendation for conical Nozzle: 0.5 - 1.5 

        Parameters
        ----------
        rf_con : Float
            Fraction of the throat radius.

        Returns
        -------
        None.

        '''
        if not(type(rf_con) == float or type(rf_con) == int):
            raise TypeError("Radius fraction must be a Number.")
        if rf_con <= 0:
            raise ValueError("Radius fraction must be greater than 0.")
            
        self.rf_con = rf_con
        
    def set_rf_div(self, rf_div):
        '''
        Sets the radius of the throat downstream as a fraction of the throat radius.
        Recommendation for bell nozzle: 0.382
        Recommendation for conical Nozzle: 0.5 - 1.5

        Parameters
        ----------
        rf_div : Float
            Fraction of the throat radius.

        Returns
        -------
        None.

        '''
        if not(type(rf_div) == float or type(rf_div) == int):
            raise TypeError("Radius fraction must be a Number.")
        if rf_div <= 0:
            raise ValueError("Radius fraction must be greater than 0.")
            
        self.rf_div = rf_div
        
    def set_l_c(self, l_c, unit="m"):
        '''
        Sets the length of the combustion chamber.

        Parameters
        ----------
        l_c : Float
            Length of combustion chamber.
        unit : String, optional
            Unit of length. The default is "m".

        Returns
        -------
        None.

        '''
        #get conversion factor from unit
        factor = 1
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            raise ValueError("Unknown length unit: " + unit)
        if not(type(l_c) == float or type(l_c) == int):
            raise TypeError("Chamber length must be a Number.")
        if l_c < 0:
            raise ValueError("Chamber length must be positive.")
            
        self.l_c = l_c * factor
        
    def set_lf_n(self, lf_n):
        '''
        Sets the nozzle length as a fraction of a comparable 30° conical nozzle.
        Only applicable to bell nozzles

        Parameters
        ----------
        lf_n : Float
            Fractional length.

        Returns
        -------
        None.

        '''
        if not(type(lf_n) == float or type(lf_n) == int):
            raise TypeError("Length fraction must be a Number.")
        if not(0 <= lf_n <= 1):
            raise ValueError("Length fraction must be between 0 and 1, is " + str(lf_n))
            
        self.lf_n = lf_n
        
    def set_alpha_con(self, alpha_con, unit="degrees"):
        '''
        Sets the half-angle of the converging section.

        Parameters
        ----------
        alpha_con : Float
            Converging half-angle.
        unit : String, optional
            Unit of angle ("degrees" or "radians"). The default is "degrees".

        Returns
        -------
        None.

        '''
        if not(type(alpha_con) == float or type(alpha_con) == int):
            raise TypeError("Converging angle must be a Number.")
        if alpha_con <= 0:
            raise ValueError("Converging angle must be greater than 0.")
            
        if unit.casefold() == "degrees":
            self.alpha_con = math.radians(alpha_con)
        elif unit.casefold() == "radians":
            self.alpha_con = alpha_con
        else:
            raise ValueError("Unknown angle unit: " + unit)
            
    def set_alpha_divt(self, alpha_divt, unit="degrees"):
        '''
        Sets the half-angle of the diverging section at the throat.

        Parameters
        ----------
        alpha_divt : Float
            Diverging half-angle at throat.
        unit : String, optional
            Unit of angle ("degrees" or "radians"). The default is "degrees".

        Returns
        -------
        None.

        '''
        if not(type(alpha_divt) == float or type(alpha_divt) == int):
            raise TypeError("Diverging throat angle must be a Number.")
        if alpha_divt <= 0:
            raise ValueError("Diverging throat angle must be greater than 0.")
            
        if unit.casefold() == "degrees":
            self.alpha_divt = math.radians(alpha_divt)
        elif unit.casefold() == "radians":
            self.alpha_divt = alpha_divt
        else:
            raise ValueError("Unknown angle unit: " + unit)
            
    def set_alpha_dive(self, alpha_dive, unit="degrees"):
        '''
        Sets the half-angle of the diverging section at the nozzle exit.

        Parameters
        ----------
        alpha_dive : Float
            Diverging half-angle at the exit.
        unit : String, optional
            Unit of angle ("degrees" or "radians"). The default is "degrees".

        Returns
        -------
        None.

        '''
        if not(type(alpha_dive) == float or type(alpha_dive) == int):
            raise TypeError("Diverging exit angle must be a Number.")
        if alpha_dive <= 0:
            raise ValueError("Diverging exit angle must be greater than 0.")
            
        if unit.casefold() == "degrees":
            self.alpha_dive = math.radians(alpha_dive)
        elif unit.casefold() == "radians":
            self.alpha_dive = alpha_dive
        else:
            raise ValueError("Unknown angle unit: " + unit)
    
    def set_epsilon(self, epsilon, fixed_throat = True):
        '''
        Sets the expansion area ratio of the nozzle.
        Overwrites either exit radius or throat radius.

        Parameters
        ----------
        epsilon : Float
            Expansion area ratio.
        fixed_throat : Boolean, optional
            Throat radius is kept, exit radius is adjusted. The default is True.

        Returns
        -------
        None.

        '''
        if not(type(epsilon) == float or type(epsilon) == int):
            raise TypeError("Expansion area ratio must be a Number.")
        if epsilon <= 0:
            raise ValueError("Expansion area ratio must be greater than 0.")
            
        self.epsilon = epsilon
        #change r_t or r_e
        if fixed_throat:
            self.r_e = (epsilon * self.r_t **2)**0.5
        else:
            self.r_t = (self.r_e**2 / epsilon)**0.5
        
    def set_steps(self, steps):
        '''
        Sets the number of discrete calculation steps/sections.
        Overwrites resolution

        Parameters
        ----------
        steps : Integer
            Number of steps.

        Returns
        -------
        None.

        '''
        if type(steps) != int:
            raise TypeError("Number of steps must be an Integer.")
        if steps <= 0:
            raise ValueError("Number of steps must be greater than 0.")
            
        self.steps = int(steps)
        self.resolution = 0
    
    def set_resolution(self, resolution, unit="m"):
        '''
        Sets the resolution of the contour generation (steps per unit of length).
        Overwrites number of steps.

        Parameters
        ----------
        resolution : Float
            Steps per [unit].
        unit : String, optional
            Unit of length. The default is "m".

        Returns
        -------
        None.

        '''
        #get conversion factor from unit
        factor = 1
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            raise ValueError("Unknown length unit: " + unit)
        if not(type(resolution) == float or type(resolution) == int):
            raise TypeError("Resolution must be a Number.")
        if resolution <= 0:
            raise ValueError("Resolution must be greater than 0.")
            
        self.resolution = resolution * factor
        self.steps = 0
    
    def set_nozzle_type(self, nozzle_type):
        '''
        Sets the type of nozzle.
        "bell", "b" for parabolic approximation of a bell nozzle.
        "conical", "c" for conical nozzle.

        Parameters
        ----------
        nozzle_type : String
            Nozzle type.

        Returns
        -------
        None.

        '''
        if type(nozzle_type) != str:
            raise TypeError("Nozzle type must be a String.")
        self.nozzle_type = nozzle_type.casefold()
        
    def set_chamber_type(self, chamber_type):
        '''
        Sets the type of chamber.
        "cylindrical", "c" for cylindrical chamber.
        "spherical", "s" for spherical chamber. NOT YET SUPPORTED
        "near spherical", "ns" for near spherical chamber. NOT YET SUPPORTED

        Parameters
        ----------
        chamber_type : String
            Chamber type.

        Returns
        -------
        None.

        '''
        if type(chamber_type) != str:
            raise TypeError("Chamber type must be a String.")
        self.chamber_type = chamber_type.casefold()
        
    def set_chamber(self, chamber):
        '''
        Sets if a chamber should be generated.

        Parameters
        ----------
        chamber : Boolean
            Should a chamber be generated?

        Returns
        -------
        None.

        '''
        if type(chamber) != bool:
            raise TypeError("Chamber must be a Boolean.")
        self.chamber = chamber
    
    def get_nozzle_parameters(self):
        '''
        Returns all set nozzle and chamber parameters in human legible form. 

        Returns
        -------
        str
            Parameters.

        '''
        params = ("====================================\n" + 
                  "{:20}{:6.2f} mm\n".format("Chamber radius", self.r_c * 1000) +
                  "{:20}{:6.2f} mm\n".format("Throat radius", self.r_t * 1000) + 
                  "{:20}{:6.2f} mm\n".format("Exit radius", self.r_e * 1000) + 
                  "\n" +
                  "{:20}{:6.2f}\n".format("Expansion ratio", self.epsilon) +
                  "\n" +
                  "{:20}{:6.2f} °\n".format("Converging angle", math.degrees(self.alpha_con)) +
                  "{:20}{:6.2f} °\n".format("Diverging angle 1", math.degrees(self.alpha_divt)) +
                  "{:20}{:6.2f} °\n".format("Diverging angle 2", math.degrees(self.alpha_dive)) +
                  "\n" +
                  "{:20}{:6.2f} mm\n".format("Chamber length", self.l_c * 1000) +
                  "{:20}{:6.2f} mm\n".format("Overall length", self.l * 1000) +
                  "\n" +
                  "{:20}{:<6s}\n".format("Nozzle type", self.nozzle_type) +
                  "{:20}{:<6s}\n".format("Chamber type", self.chamber_type) +
                  "\n" +
                  "{:20}{:6d}\n".format("Steps", self.steps) +
                  "{:20}{:6.2f} Steps/mm\n".format("Resolution", self.resolution / 1000) +
                  "===================================="
                  )
        return params
    
    def get_l(self, unit="m"):
        '''
        Returns overall length of nozzle + chamber

        Parameters
        ----------
        unit : String, optional
            Unit of length. The default is "m".

        Returns
        -------
        Float
            Length og nozzle + chamber.

        '''
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            raise ValueError("Unknown length unit: " + unit)
        return self.l / factor
    
    def get_l_n(self, unit="m"):
        '''
        Returns length of nozzle section

        Parameters
        ----------
        unit : String, optional
            Unit of length. The default is "m".

        Returns
        -------
        Float
            length of nozzle.

        '''
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            raise ValueError("Unknown length unit: " + unit)
        return self.l_n / factor
    
    def get_l_t(self, unit="m"):
        '''
        Returns the downstream distance to the throat

        Parameters
        ----------
        unit : String, optional
            Unit of length. The default is "m".

        Returns
        -------
        Float
            length to throat.

        '''
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            raise ValueError("Unknown length unit: " + unit)
        return self.l_t / factor
    
    def get_n_t(self):
        '''
        Returns the step number of the throat

        Returns
        -------
        Integer
            Step number of throat.

        '''
        return self.n_t
    
    def save_csv(self, name, separator=";", dec=","):
        '''
        Saves the data to disc as a .csv file.

        Parameters
        ----------
        name : String
            Name or path of target file. End with .csv or .txt.
        separator : String, optional
            Column separator. The default is ";" for Excel.
        dec : String, optional
            Decimal symbol. The default is "," for european standards. 

        Returns
        -------
        None.

        '''
        frame = self.export_pandas()
        frame.to_csv(name, sep=separator, decimal=dec)
    
    def export_np(self):
        '''
        Returns the data as a numpy array.

        Returns
        -------
        numpy.ndarray
            Data as a numpy array.

        '''
        return self.data
    
    def export_pandas(self):
        '''
        Converts and returns the data as a pandas Dataframe.
        Column names:
            Step: calculation step
            x: length coordinate
            y: radius coordinate

        Returns
        -------
        pandas.DataFrame
            Data frame of data.

        '''
        import pandas as pd
        
        step = np.array(np.zeros(self.data[:,self.__step].shape[0]), dtype=int)
        for i in range(0, len(self.data[:,0]) - 1):
            step[i] = int(self.data[i,self.__step])
        return pd.DataFrame({"x": self.data[:,self.__cx], "y": self.data[:,self.__cy]})
    
    def export_graph(self, name, fformat, date=True, title="Nozzle", dpi=None):
        '''
        Exports the nozzle contour as a picture using matplotlib.

        Parameters
        ----------
        name : String
            Name or path of output file.
        fformat : String
            Desired output file type (e.g. "png", "svg", "pdf"...).
        date : Boolean, optional
            Sets if the current date and time should be added to the titel. The default is True.
        title : String, optional
            Title of the plot. The default is "Nozzle".
        dpi : Float, optional
            DPI of the exported image. The default is 80.   
            
        Returns
        -------
        None.

        '''
        import matplotlib.pyplot as plt
        
        now = ""
        if date:
            from datetime import datetime
            now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.data[:,self.__cx], self.data[:, self.__cy])
        ax.set_xlim(0,self.l)
        ax.set_ylim(bottom=0)
        ax.axis("equal")
        ax.set(xlabel="length [m]", ylabel="radius [m]")
        fulltitle = str(title) + "  " + str(now)
        ax.set_title(fulltitle)
        plt.savefig(name, dpi=dpi, format=fformat)
    
    def draw_contour(self):
        '''
        Draws the Nozzle contour as a Matplotlib plot inline of the console.

        Returns
        -------
        None.

        '''
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.data[:,self.__cx], self.data[:, self.__cy])
        #ax.axis("equal")
        #ax.set_xlim(0,self.l)
        #ax.set_ylim(0,(self.r_c * 1.1) if (self.r_c < self.r_e) else (self.r_e * 1.1) )
        ax.set(xlabel="length [m]", ylabel="radius [m]")
        #plt.axis([0,self.l,0,(self.r_c * 1.1) if (self.r_c < self.r_e) else (self.r_e * 1.1)],aspect = 'equal')
        ax.set_aspect("equal", adjustable="datalim", anchor="SW")
        #ax.set_ylim(0,(self.r_c * 1.1) if (self.r_c < self.r_e) else (self.r_e * 1.1))
        #keine ahnung wie man gleiche Skalierung und eigene Limits macht :/
        plt.show()


if __name__ == "__main__":
    print("This file is not meant to be run alone. It was written to be importet as a custom module.")
    print("Here is a demonstration:")
    test = Nozzle()
    test.sample()
    print("I have generated a sample nozzle and chamber with the following parameters:")
    print(test)
    input()
    print("Here is the contour as a Graph!")
    test.draw_contour()
    input()
    print("Let's change the nozzle type to conical and the expansion ratio to 10. We should also change the diverging angle.")
    test.set_nozzle_type("conical")
    test.set_epsilon(10)
    test.set_alpha_divt(30)
    print("Don't forget to generate!")
    input()
    test.generate()
    print("Here is our new nozzle:")
    print(test)
    test.draw_contour()
    print("For more information contact " + __author__)