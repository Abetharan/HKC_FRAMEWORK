import numpy as np
import impact_norms_py3 as ImNorms
import impact_module_py3 as cf
import os
from scipy import constants
from scipy.interpolate import CubicSpline
import fnmatch
import sys
import re
import matplotlib.pyplot as plt
import HydroRemainingFromImpact as remain


kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")

def findLastIndex(path, var):
    LargestIndex = 0
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, "*_" + var + "_*"):
            k = os.path.splitext(file)[0]
            k = k.split("_")
            timeIndex = k[-1]
            if int(timeIndex) > LargestIndex:
                LargestIndex = int(timeIndex)
            
    return(LargestIndex)
def UnitConversion(varName,var,To,conversionFactor = None, density = None):
    """
    Purpose: Convert SI to appropriate Plasma Units i.e. eV, cm^-3 and Wm^-3
    var = Variable to convert, T, ne, rad(radiation) or x
    To = SI to Plasma if set to True else Plasma to SI
    conversionFactor = if to other units add the conversion factor 
    density = defaults to None but required to convert to Wm^-3 for rad

    Returns: Converted unit
    """
    if varName == "T":
        if To:
            if conversionFactor is not None:
                converted_var = var * conversionFactor
            else:
                converted_var = var * (kb/e)
        else:
            if conversionFactor is not None:
                converted_var = var * (1/conversionFactor)
            else:
                converted_var = var * (e/kb)
    elif varName.lower() == "n":
        if To:
            if conversionFactor is not None:
                converted_var = var * conversionFactor
            else:
                converted_var = var * 1e-6
        else:
            if conversionFactor is not None:
                converted_var = var * (1/conversionFactor)
            else:
                converted_var = var * 1e-6
    elif varName == "rad":
        if To:
            converted_var = var * density 
        else:
            converted_var = var * density
    elif varName == "x":
        if To:
            converted_var = var * 1
        else:
            converted_var = var * 1
    return converted_var

    
BASE_DIR_ = "/media/abetharan/DATADRIVE1/Abetharan/"
RUN_NAME_ = "Test_1"
cycleNo = "/cycle_1"
cycleNoNext =  "/cycle_2"
nextFluidInput = BASE_DIR_ + RUN_NAME_ + cycleNoNext + "/fluid_input"
fluidOutPath = BASE_DIR_ + RUN_NAME_ + cycleNo + "/fluid_output"
kineticInput = BASE_DIR_ + RUN_NAME_ + cycleNo + "/kinetic_input"
kineticOutput = BASE_DIR_ + RUN_NAME_ + cycleNo + "/kinetic_input"

_KINETIC_nx = 30
_KINETIC_ny = 1 
_KINETIC_nv = 300
_KINETIC_np = 1
_FLUID_nx = 100
_CYCLES  = 3

Z = 1#60
Ar = 1#157

avgNe = 1e21 #np.average(fluid_ne) * 1e-6
avgTe = 1000 #np.average(fluid_Te) * (kb / e)
normalised_values = ImNorms.impact_inputs(avgNe, avgTe, Z, Ar, Bz = 0)

def CheckInOut(fluidPath, kineticPath, normalised_values, fluidNx, kineticNx, kineticIN=False, FluidIN=False):
    

    if kineticIN:
        ##Checck kinetic input against fluid output 
        ##required paths are Fluid out put in same CYCLE path as Kinetic in 
        dict_to_match = {'ionden':"NumberDensityI", 'eden':"NumberDensityE", 'laserdep':"InverseBrem", 'rad_to_electron':"Brem",
                        'tmat':'TemperatureE','xf':"Coord", 'zstar': Z}
        
        for file_var_path in os.listdir(kineticPath):
            if fnmatch.fnmatch(file_var_path, "*_xf*"):
                var = cf.load_dict(file_var_path)['mat']
            else:
                var = cf.load_dict(file_var_path)['mat'][1:-1]
            file_var_check = os.path.splitext(file_var_path)[0]
            var_check = file_var_check.split("_")[-1]
            if var_check == 'laserdep' or var_check == 'rad_to_electron':
                continue
            fluid_var_name = dict_to_match[str(var_check)]
            largest_index = findLastIndex(fluidPath, fluid_var_name)
            fluid_x  = np.loadtxt(fluidOutPath + "/Coord_" + str(largest_index) + ".txt") / normalised_values['lambda_mfp']
            fluid_var_list = np.loadtxt(fluidPath + "/" + fluid_var_name + "_" + largest_index + ".txt") 
            
            conv_fluid_var = UnitConversion(fluid_var_name[0], fluid_var_list, To = True) / 
                                (1e6 *  normalised_values['ni'] * 1e21) 
            
            kinetic_x = np.linspace(x_norm[0], x_norm[-1], int(kineticNx + 1))
            fluid_centered_x = [(x_norm[i + 1] + x_norm[i])/2 for i in range(fluidNx)]
            kinetic_centered_x = [(kinetic_x[i + 1] + kinetic_x[i])/2 for i in range(kineticNx)]
            
            cs_ne = CubicSpline(fluid_centered_x, conv_fluid_var)
            kinetic_ne = cs_ne(kinetic_centered_x)
            diff = 

            if any(diff != 0):
                print(fluid_var_name)
                exit(1)
    
    
    
    
    
    
    if FluidIN:
        #Checks kinetic output against fluid input
        #reuqire paths from Kinetic ouput agains the previous cycle agsinst current cycle fluid input path
        print("kappa")


        











