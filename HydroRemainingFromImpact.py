import numpy as np
import os
import fnmatch

protonMass = 1.66e-27
kb = 1.38e-23

def CalculateMass(coord, density, nx):
    mass = np.zeros(nx)
    for i in range(nx):

        mass[i] = density[i]  * (coord[i + 1] - coord[i])
    return(mass)

def ConstCheck(density, ne, ni, Te, Ti, Pe, Pi, Ptot, massNumber, Z):
    import sys

    if any(density != ni * massNumber * protonMass) or any(density != ((ne * massNumber * protonMass) / Z)):
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

def InvBrem(coord,ne, nc, wavelength, Z, coulombLog, Temperature,LaserPower, mass, nx):
    mNx = nx
    laserPower = LaserPower
    alpha = np.zeros(mNx)
    CoulombLog = np.zeros(nx) + coulombLog
    CoulombLog.fill(11)
    powerAbsorbed = np.zeros(nx)

    for i in range(mNx): 
        DistanceTravelled = coord[i + 1] - coord[i]
        beta = ne[i] / nc
        alpha[i] = 1.23E-14 * ne[i] * Z * CoulombLog[i] * pow(Temperature[i], -1.5) * pow(beta, 2) * pow((1 - beta),-0.5)
    
        transmittedLaser = laserPower * np.exp(-alpha[i] * DistanceTravelled)
        powerAbsorbed[i] = laserPower - transmittedLaser

        invBremPower = 0.5 * (1/mass[i]) * powerAbsorbed
        laserPower = transmittedLaser
    return(invBremPower,powerAbsorbed)
    
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

def CalculateRemain(ne, Te, qe, Z, massNumber, cycle, gammaFactor, laserWaveLength, laserPower, nx, cycleDumpPath, ionfluidDensity = None, ionfluidTemperature = None):
    
    if ionfluidDensity is not None:
        ni = ionfluidDensity
        Ti = ionfluidTemperature
    else:
        for fileconstituent in cycleDumpPath.split():
            
            if fileconstituent == "cycle_" + "[0-255]":
                new_path  = new_path + "cycle_" + str(cycle)
            else:
                new_path = new_path + fileconstituent
        previous_cycle_fluid_out_path  = new_path + "/fluid_out/"

        LargestIndex = 0
        for file in os.listdir(previous_cycle_fluid_out_path):
            if fnmatch.fnmatch(file, "[a-z]"+ "NumberDensityI"):
                k = os.path.splitext(file)[0]
                k = k.split("_")
                timeIndex = k[-1]
                if timeIndex > LargestIndex:
                    LargestIndex = timeIndex
        ni = np.loadtxt(previous_cycle_fluid_out_path + "NumberDensityI_" + str(LargestIndex) + ".txt")
        Ti = np.loadtxt(previous_cycle_fluid_out_path + "TemperatureI_" + str(LargestIndex) + ".txt")
        coord = np.loadtxt(previous_cycle_fluid_out_path + "Coord_" + str(LargestIndex) + ".txt")
        velocity = np.loadtxt(previous_cycle_fluid_out_path + "Velocity_" + str(LargestIndex) + ".txt")
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

    mass = CalculateMass(coord, density, nx)
    ConstCheck(density, ne, ni, Te, Ti, pressureE, pressureI, pressureTotal, massNumber, Z)

    nc = 1.1E15 / pow(laserWaveLength, 2)
    InvBrem_ = InvBrem(coord, ne, nc, laserWaveLength, Z, 11, Te, laserPower, mass, nx)
    Brem_ = Brem(ne, massNumber, Z, Te, nx)
    electronheatflow= Heatflow(qe, mass, nx)

    TextDump(path = cycleDumpPath + "/fluid_init/",
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