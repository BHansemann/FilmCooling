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
    
    __l1, __l2, __l3, __l4, __l5, __l6 = 0 ,0, 0, 0 ,0 ,0 #needed for calculation
    __step, __cx, __cy = 0, 1, 2 #key for numpy array columns
    
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
        self.epsilon = 0 #expansion ratio
        
        self.nozzle_type = "b" #nozzle type (bell or conical) compare with str.casefold()
        self.chamber_type = "cylindrical" #chamber shape (cylindrical), in the future also spherical and near-spherical
        self.chamber = True #with or without chamber
        
        self.l = 0 #overall length
        self.l_t = 0 #length to throat
        self.l_e = 0 #length to exit
        self.l_n = 0 #nozzle length
        
        self.steps = 10000 #number of steps
        self.resolution = 0 #number of steps per meter
        
        self.data = np.array(np.zeros((self.steps + 1, 3)))
        self.data[:,self.__step] = np.arange(0, self.steps + 1)
    
    def __repr__(self):
        pass
    
    def __array__(self):
        pass
    
    def set_values(self):
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
        #print(self.l, __l1, __l2, __l3, __l4, __l5, __l6)
        
        x_0 = self.l - __l6
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
        
        def get_y(x):
            if x <= __l1:
                x1 = x
                return self.r_c
            elif __l1 < x <= __l1 + __l2:
                x2 = x - __l1
                return self.r_c * (1 - self.rf_cf) + ((self.r_c * self.rf_cf)**2 - x2**2)**0.5
            elif __l1 + __l2 < x <= __l1 + __l2 + __l3:
                x3 = x - __l1 - __l2
                return self.r_c * (1 - self.rf_cf * (1 - math.cos(self.alpha_con))) - x3 * math.tan(self.alpha_con)
            elif __l1 + __l2 + __l3 < x <= __l1 + __l2 + __l3 + __l4:
                x4 = x - __l1 - __l2 - __l3
                return self.r_t + self.r_t * self.rf_con - ((self.r_t * self.rf_con)**2 - (x4 - __l4)**2)**0.5
            elif __l1 + __l2 + __l3 + __l4 < x <= __l1 + __l2 + __l3 + __l4 + __l5:
                x5 = x - __l1 - __l2 - __l3 - __l4
                return self.r_t + self.r_t * self.rf_div - ((self.r_t * self.rf_div)**2 - x5**2)**0.5
            elif __l1 + __l2 + __l3 + __l4 + __l5 < x <= self.l:
                x6 = x - __l1 - __l2 - __l3 - __l4 - __l5
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
                    pass #error handling
                    
        self.data[:,self.__cx] = self.data[:,self.__step] * (self.l / self.steps)
        self.data[:,self.__cy] = np.vectorize(get_y)(self.data[:,self.__cx])
        
    def test(self):
        self.r_c = 0.03
        self.r_t = 0.02
        self.r_e = 0.04
        
        self.rf_cf = 0.08
        self.rf_con = 1.5
        self.rf_div = 0.382
        
        self.l_c = 0.08
        
        self.lf_n = 0.8
        
        self.alpha_con = math.radians(25)
        self.alpha_divt = math.radians(30)
        self.alpha_dive = math.radians(10)
        
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
        factor = 1
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            pass #error handling
        self.r_c = r_c * factor
        
    def set_r_t(self, r_t, unit="m"):
        '''
        Sets the radius of the throat.

        Parameters
        ----------
        r_t : Float
            Throat radius.
        unit : String, optional
            Unit of Length. The default is "m".

        Returns
        -------
        None.

        '''
        factor = 1
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            pass #error handling
        self.r_t = r_t * factor
        
    def set_r_e(self, r_e, unit="m"):
        '''
        Sets exit radius of the nozzle.

        Parameters
        ----------
        r_e : Float
            Exit Radius.
        unit : String, optional
            Unit of length. The default is "m".

        Returns
        -------
        None.

        '''
        factor = 1
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            pass #error handling
        self.r_e = r_e * factor
        
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
        factor = 1
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            pass #error handling
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
        if unit.casefold() == "degrees":
            self.alpha_con = math.radians(alpha_con)
        elif unit.casefold() == "radians":
            self.alpha_con = alpha_con
            
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
        if unit.casefold() == "degrees":
            self.alpha_divt = math.radians(alpha_divt)
        elif unit.casefold() == "radians":
            self.alpha_divt = alpha_divt
            
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
        if unit.casefold() == "degrees":
            self.alpha_dive = math.radians(alpha_dive)
        elif unit.casefold() == "radians":
            self.alpha_dive = alpha_dive
    
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
        self.epsilon = epsilon
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
        self.steps = int(steps)
        self.resolution = 0
    
    def set_resolution(self, resolution, unit="m"):
        '''
        Sets the resolution of the contour generation.
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
        factor = 1
        if unit == "m":
            factor = 1
        elif unit == "mm":
            factor = 1E-3
        elif unit == "cm":
            factor = 1E-2
        else:
            pass #error handling
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
        if chamber.casefold() in ["n", "nein", "no", "false", "ohne"]:
            chamber = False
        self.chamber = chamber
    
    def get_nozzle_parameters(self):
        '''
        Returns all set nozzle and chamber parameters in human legible form. 

        Returns
        -------
        str
            Parameters.

        '''
        return ""
    
    def export_as(self, ex_type: str):
        pass
    
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
        step = np.array(np.zeros(1, self.data.shape[0]), dtype=int)
        for i in range(0, len(self.data) - 1):
            step[i] = int(self.data[i,self.__step])
        return pd.DataFrame({"Step": step, "x": self.data[:,self.__cx], "y": self.data[:,self.__cy]})
    
    def export_graph(self, name, fformat, date=True, title="Nozzle"):
        if date:
            from datetime import datetime
            now = datetime.now().strftime("%d/%m/%Y %H:&M:&S")
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.data[:,self.__cx], self.data[:, self.__cy])
        ax.set_xlim(0,self.l)
        ax.set_ylim(bottom=0)
        ax.axis("equal")
        ax.set(xlabel="length [m]", ylabel="radius [m]")
        ax.title(title + now)
        plt.savefig(name, format=fformat)
    
    def draw_contour(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.data[:,self.__cx], self.data[:, self.__cy])
        ax.set_xlim(0,self.l)
        ax.set_ylim(bottom=0)
        ax.axis("equal")
        ax.set(xlabel="length [m]", ylabel="radius [m]")
        plt.show()