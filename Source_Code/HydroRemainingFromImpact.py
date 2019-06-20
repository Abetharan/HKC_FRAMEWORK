import numpy as np
import os
import fnmatch
import math
from scipy.interpolate import CubicSpline
from scipy import interpolate
protonMass = 1.66e-27
kb = 1.38e-23

def CalculateMass(coord, density, nx):
    mass = np.zeros(nx)
    for i in range(nx):

        mass[i] = density[i]  * (coord[i + 1] - coord[i])
    return(mass)

def ConstCheck(density, ne, ni, Te, Ti, Pe, Pi, Ptot, Ar, Z):
    import sys

    for i in range(len(density)):
        if density[i] != ((ne[i]/ Z[i]) * Ar[i] * protonMass):
            print("DENSITY 1 ARE NOT CALCULATED CONSISTENTLY")
            print("VALUES")
            print(density[i])
            print(((ne[i] * Ar[i] * protonMass) / Z[i]))
            sys.exit(0)
        if density[i] != ni[i] * Ar[i] * protonMass:
            print("DENSITY 2 ARE NOT CALCULATED CONSISTENTLY")
            print("VALUES")
            print(density[i])
            print(ni[i] * Ar[i] * protonMass)
            sys.exit(0)

        if Ptot[i] != (Pe[i] + Pi[i]):
            print("PRESSURES ARE NOT CALCULATED CONSISTENTLY")
            print("VALUES")
            print(Ptot[i])
            print(Pe[i] + Pi[i])
            sys.exit(0)


def InvBrem(coord,ne, nc, wavelength, mZ, coulombLog, TemperatureE,LaserPower, mass, nx, start = "right"):
    mNx = nx
    laserPower = LaserPower
    alpha = np.zeros(mNx)
    CoulombLog = np.zeros(nx) + coulombLog
    CoulombLog.fill(11)
    TransmittedLaser = np.zeros(nx)
    PowerAbsorbed = np.zeros(nx)
    InverseBrem = np.zeros(nx)
    reflection_index = -1
    
    if start == "right":
        i = mNx - 1
    else:
        i = 0
    
    while(True):
        if i < 0 or i >mNx - 1:
            break
        print(i)
        if ne[i] >= nc:
            reflection_index = i
            if i == 0 or i == mNx -1 :
                break        
        if reflection_index > 0:
            print("REFLECTED")
            reflection_index = 0
            if start == "right":
                InverseBrem[0:i] = 0
                i+= 1
            else:
                InverseBrem[i:-1] = 0
                i-= 1

        DistanceTravelled = abs(coord[i + 1] - coord[i])
        beta = ne[i] / nc
        alpha[i] = 1.23E-14 * ne[i] * mZ[i] * CoulombLog[i] * pow(TemperatureE[i], -1.5) * pow(beta, 2) * pow((1 - beta),-0.5)
        
        if reflection_index == 0:
            ReflectedTranmittedLaser = laserPower * np.exp(-alpha[i] * DistanceTravelled)
            reflected_power_absorbed = laserPower - ReflectedTranmittedLaser
            InverseBrem[i] += (1/mass[i]) * reflected_power_absorbed
            PowerAbsorbed[i] += reflected_power_absorbed
            laserPower = ReflectedTranmittedLaser
            TransmittedLaser[i] += ReflectedTranmittedLaser
        else:
            TransmittedLaser[i] += laserPower * np.exp(-alpha[i] * DistanceTravelled)
            PowerAbsorbed[i] += laserPower - TransmittedLaser[i]
            laserPower = TransmittedLaser[i]
            InverseBrem[i] += (1/mass[i]) * PowerAbsorbed[i]

        if reflection_index == 0:
            if start == "right":
                i+= 1
            else:
                i-= 1
        if reflection_index != 0:
            if start == "right":
                i-= 1
            else:
                i+= 1
    
    return(InverseBrem, PowerAbsorbed,TransmittedLaser)

def Brem(ne, Ar, Z, Temperature, nx):
    BremPower = -8.5E-14 * ne * (Temperature ** 0.5) * (Z**2 / Ar)

    return(BremPower)

def Heatflow(electron_thermal_flux, mass, nx):
    HeatConductionE = np.zeros(nx)
    for i in range(nx):   
        HeatConductionE[i] = -(electron_thermal_flux[i + 1] - electron_thermal_flux[i]) / mass[i]
    return(HeatConductionE)

def TextDump(path, coord, velocity, density, ne, ni, Te, Ti, Pe, Pi, Ptot, IntEe, IntEi, DpDTe, DpDTi, Cve,Cvi, mass, invBrem, brem, heatflow,Z ,Ar):
    """ 
    Purpose: Dumps variables to text
    Args: All physical quantities and Path to dump
    Returns: Null
    """
    np.savetxt((path+"coord.txt"), coord)
    np.savetxt((path+"velocity.txt"), velocity)
    np.savetxt((path+"density.txt"), density)
    np.savetxt((path+"ne.txt"), ne)
    np.savetxt((path+"ni.txt"), ni)
    np.savetxt((path+"total_pressure.txt"),  Ptot)
    np.savetxt((path+"ion_pressure.txt"), Pi)
    np.savetxt((path+"electron_pressure.txt"), Pe)
    np.savetxt((path+"electron_temperature.txt"), Te)
    np.savetxt((path+"ion_temperature.txt"), Ti)
    np.savetxt((path+"electron_internal_energy.txt"), IntEe)
    np.savetxt((path+"ion_internal_energy.txt"), IntEi)
    np.savetxt((path+"electron_dp_dt.txt"), DpDTe)
    np.savetxt((path+"ion_dp_dt.txt"), DpDTi)
    np.savetxt((path+"electron_specific_heat.txt"), Cve)
    np.savetxt((path+"ion_specific_heat.txt"), Cvi)
    np.savetxt((path+"mass.txt"), mass)
    np.savetxt((path+"inv_brem.txt"), invBrem)
    np.savetxt((path+"brem.txt"), brem)
    np.savetxt((path+"qe.txt"),heatflow )
    np.savetxt((path+ "Z.txt"), Z)
    np.savetxt((path+ "Ar.txt"),Ar)

def CalculateRemain(kinetic_ne, kinetic_Te, kinetic_qe, normalisationValues, gammaFactor, laserWaveLength, laserPower, fluidNx, NextStepFluidInit, PreviousFluidOutPath, PreviousKineticOutPath, ionfluidDensity = None, ionfluidTemperature = None):
    
    if ionfluidDensity is not None:
        ni = ionfluidDensity
        Ti = ionfluidTemperature
    else:
    ##########################################################################################################################
    ##########################################################################################################################
    ##########################################################################################################################
    ##### Note HERE that if ION FLUIDS are included this NEEDS to be changed.     
    
        LargestIndex = 0
        for file in os.listdir(PreviousFluidOutPath):
            if fnmatch.fnmatch(file, "NumberDensityI_*"):
                k = os.path.splitext(file)[0]
                k = k.split("_")
                timeIndex = k[-1]
                if int(timeIndex) > LargestIndex:
                    LargestIndex = int(timeIndex)
            
        ni = np.loadtxt(PreviousFluidOutPath + "NumberDensityI_" + str(LargestIndex) + ".txt")
        Ti = np.loadtxt(PreviousFluidOutPath + "TemperatureI_" + str(LargestIndex) + ".txt")
        coord = np.loadtxt(PreviousFluidOutPath + "Coord_" + str(LargestIndex) + ".txt")
        velocity = np.loadtxt(PreviousFluidOutPath + "Velocity_" + str(LargestIndex) + ".txt")
        Z = np.loadtxt(PreviousFluidOutPath + "Z_" + str(LargestIndex) + ".txt")
        Ar = np.loadtxt(PreviousFluidOutPath + "Ar_"+ str(LargestIndex) + ".txt")
    #material
    # Ar = normalisationValues['Ar']
    # Z = normalisationValues['Z']
    # Handle the Interpolation back. Cubic 

    kinetic_x = np.loadtxt(PreviousKineticOutPath + "ReturnToHydro_xf.xy", delimiter = "\n") * normalisationValues['lambda_mfp'] 
    kinetic_centered_x = [(kinetic_x[i + 1] + kinetic_x[i])/2 for i in range(int(os.environ["NXM"]))]
    fluid_centered_x = [(coord[i + 1] + coord[i])/2 for i in range(fluidNx)]
    

    #### NOTE HERE IT IS ASSUMED IDEAL GAS THIS SHOULD ALSO BE CHANGED #####
    kinetic_pressure = kinetic_ne * kinetic_Te * kb
    cs_pressure = CubicSpline(kinetic_centered_x, kinetic_pressure)
    cs_ne = CubicSpline(kinetic_centered_x, kinetic_ne)#interpolate.interp1d(kinetic_centered_x, kinetic_ne,fill_value="extrapolate")#
    cs_qe = CubicSpline(kinetic_x, kinetic_qe)
    pressureE = cs_pressure(fluid_centered_x)
    qe = cs_qe(coord)
    ne = cs_ne(fluid_centered_x)
    Te = pressureE/(ne*kb)
    ni = ne/Z

    # cs_Te = interpolate.interp1d(kinetic_centered_x, kinetic_Te,fill_value="extrapolate")
    # cs_qe = interpolate.interp1d(kinetic_x, kinetic_qe,fill_value="extrapolate")

    # ne = cs_ne(fluid_centered_x)
    # Te = cs_Te(fluid_centered_x)
    # qe = cs_qe(coord)
    # import matplotlib.pyplot as plt
    # plt.plot(Te)
    # plt.plot(kinetic_Te)
    # plt.show()

    density = ni * Ar  * protonMass
    ##Noting convention that Specific Heat ha
    specificHeatE =  (ne * kb) / (density * (gammaFactor -1))
    DpDTe = ne *kb
    
    specificHeatI =  (ni * kb) / (density * (gammaFactor -1))
    DpDTi = ni *kb

    #Pressure
    pressureI = ni * kb * Ti
    pressureTotal = pressureE + pressureI

    #Internal Energy
    IntEe = pressureE / ((gammaFactor - 1) * density)
    IntEi = pressureI / ((gammaFactor - 1) * density)

    mass = CalculateMass(coord, density, fluidNx)
    ConstCheck(density, ne, ni, Te, Ti, pressureE, pressureI, pressureTotal, Ar, Z)

    nc = 1.1E15 / pow(laserWaveLength, 2)
    InvBrem_, _, _ = InvBrem(coord, ne, nc, laserWaveLength, Z, 11, Te, laserPower, mass, fluidNx, "left")
    
    Brem_ = Brem(ne, Ar, Z, Te, fluidNx)
    electronheatflow= Heatflow(qe, mass, fluidNx)
    TextDump(path = NextStepFluidInit,
            coord= coord,
            velocity = velocity,
            density = density,
            ne = ne,
            ni = ni ,
            Te = Te,
            Ti = Ti,
            Pe = pressureE,
            Pi = pressureI,
            Ptot = pressureTotal,
            IntEe = IntEe,
            IntEi = IntEi,
            DpDTe = DpDTe,
            DpDTi = DpDTi,
            Cve = specificHeatE,
            Cvi = specificHeatI,
            mass = mass,
            invBrem = InvBrem_,
            brem = Brem_,
            heatflow = electronheatflow, 
            Z = Z,
            Ar = Ar,
        )