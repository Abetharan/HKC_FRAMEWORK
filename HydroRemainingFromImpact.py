import numpy as np
import os
import fnmatch
import math
from scipy.interpolate import CubicSpline

protonMass = 1.66e-27
kb = 1.38e-23

def CalculateMass(coord, density, nx):
    mass = np.zeros(nx)
    for i in range(nx):

        mass[i] = density[i]  * (coord[i + 1] - coord[i])
    return(mass)

def ConstCheck(density, ne, ni, Te, Ti, Pe, Pi, Ptot, massNumber, Z):
    import sys

    if any(density != ni * massNumber * protonMass): #or any(density != ((ne * massNumber * protonMass) / Z)):
        print("DENSITY ARE NOT CALCULATED CONSISTENTLY")
        print("VALUES")
        print(density)
        print(ni * massNumber * protonMass)
        print(((ne * massNumber * protonMass) / Z))
        sys.exit(0)

    if any(Ptot != (Pe + Pi)):
        print("PRESSURES ARE NOT CALCULATED CONSISTENTLY")
        print("VALUES")
        print(Ptot)
        print(Pe + Pi)
        sys.exit(0)


def InvBrem(coord,ne, nc, wavelength, mZ, coulombLog, TemperatureE,LaserPower, mass, nx):
    mNx = nx
    laserPower = LaserPower
    alpha = np.zeros(mNx)
    CoulombLog = np.zeros(nx) + coulombLog
    CoulombLog.fill(11)
    TransmittedLaser = np.zeros(nx)
    PowerAbsorbed = np.zeros(nx)
    InverseBrem = np.zeros(nx)
    reflection_index = 0
    for i in range(mNx - 1, -1, -1):    
        print(i)
        DistanceTravelled = coord[i + 1] - coord[i]
        beta = ne[i] / nc
        alpha[i] = 1.23E-14 * ne[i] * mZ * CoulombLog[i] * pow(TemperatureE[i], -1.5) * pow(beta, 2) * pow((1 - beta),-0.5)
    
        TransmittedLaser[i] = laserPower * np.exp(-alpha[i] * DistanceTravelled)
        PowerAbsorbed[i] = laserPower - TransmittedLaser[i]
        if ne[i] >= nc:
        
            PowerAbsorbed[i] = 0
            reflection_index = i
            reflection_power = laserPower
        
        if reflection_index > 0:
            InverseBrem[i] = 0
        else:
            InverseBrem[i] = (1/mass[i]) * PowerAbsorbed[i]
        if math.isnan(InverseBrem[i]):
            print("kappa")
        laserPower = TransmittedLaser[i]

    if reflection_index > 0:
        print("reflecting")
        for i in range(reflection_index+1, mNx, 1):
    
            DistanceTravelled = coord[i + 1] - coord[i]
            beta = ne[i] / nc
            alpha[i] = 1.23E-14 * ne[i] * mZ * CoulombLog[i] * pow(TemperatureE[i], -1.5) * pow(beta, 2) * pow((1 - beta),-0.5)
        
            reflectedtransmittedLaser = reflection_power * np.exp(-alpha[i] * DistanceTravelled)
            reflectedpowerAbsorbed = reflection_power - reflectedtransmittedLaser
            if ne[i] >= nc:
                reflectedpowerAbsorbed = 0
                reflectedtransmittedLaser = reflection_power
            
            if reflectedpowerAbsorbed == 0:
                InverseBrem[i] += 0
            else:
                InverseBrem[i] += (1/mass[i]) * reflectedpowerAbsorbed
            if math.isnan(InverseBrem[i]):
                print("kappa")
            PowerAbsorbed[i] += reflectedpowerAbsorbed
            reflection_power = reflectedtransmittedLaser
            TransmittedLaser[i] += reflection_power

    return(InverseBrem, PowerAbsorbed, TransmittedLaser)

    
def Brem(ne, massNumber, Z, Temperature, nx):
    BremPower = -8.5E-14 * ne * (Temperature ** 0.5) * (Z**2 / massNumber)

    return(BremPower)

def Heatflow(electron_thermal_flux, mass, nx):
    HeatConductionE = np.zeros(nx)
    for i in range(nx):   
        HeatConductionE[i] = (electron_thermal_flux[i + 1] - electron_thermal_flux[i]) / mass[i]
    return(HeatConductionE)

def TextDump(path, coord, velocity, density, ne, ni, Te, Ti, Pe, Pi, Ptot, IntEe, IntEi, DpDTe, DpDTi, Cve,Cvi, mass, invBrem, brem, heatflow):
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

def CalculateRemain(ne, Te, qe, Z, massNumber, gammaFactor, laserWaveLength, laserPower, fluidNx, NextStepFluidInit, PreviousFluidOutPath, PreviousKineticOutPath, ionfluidDensity = None, ionfluidTemperature = None):
    
    if ionfluidDensity is not None:
        ni = ionfluidDensity
        Ti = ionfluidTemperature
    else:
        
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
    
   
    # Handle the Interpolation back. Cubic 

    kinetic_x = np.loadtxt(PreviousKineticOutPath + "ReturnToHydro_xf.xy", delimiter = "\n")
    kinetic_centered_x = [(kinetic_x[i + 1] + kinetic_x[i])/2 for i in range(int(os.environ["NXM"]))]
    fluid_centered_x = [(coord[i + 1] + coord[i])/2 for i in range(fluidNx)]
    cs_ne = CubicSpline(kinetic_centered_x, ne)
    cs_Te = CubicSpline(kinetic_centered_x, Te)
    cs_qe = CubicSpline(kinetic_x, qe)

    ne = cs_ne(fluid_centered_x)
    Te = cs_Te(fluid_centered_x)
    qe = cs_qe(coord)

    density = ni * massNumber  * protonMass
    ##Noting convention that Specific Heat ha
    specificHeatE =  (ne * kb) / (density * (gammaFactor -1))
    DpDTe = ne *kb
    
    specificHeatI =  (ni * kb) / (density * (gammaFactor -1))
    DpDTi = ni *kb

    #Pressure
    pressureI = ni * kb * Ti
    pressureE = ne * kb * Te
    pressureTotal = pressureE + pressureI

    #Internal Energy
    IntEe = pressureE / ((gammaFactor - 1) * density)
    IntEi = pressureI / ((gammaFactor - 1) * density)

    mass = CalculateMass(coord, density, fluidNx)
    ConstCheck(density, ne, ni, Te, Ti, pressureE, pressureI, pressureTotal, massNumber, Z)

    nc = 1.1E15 / pow(laserWaveLength, 2)
    InvBrem_, _, _ = InvBrem(coord, ne, nc, laserWaveLength, Z, 11, Te, laserPower, mass, fluidNx)
    Brem_ = Brem(ne, massNumber, Z, Te, fluidNx)
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
            heatflow = electronheatflow
        )