""" Only Valid for idle gases 
    By default initial velocity is always zero
"""
import numpy as np
import impact_profiles_py3 as profiles
import impact_norms_py3 as norm

#global constants
kb = 1.38E-23
e = 1.6E-19
protonMass = 1.67E-27

def LoadTestProblem(testName, nx, gammaFactor, x_l = 0, x_u = 1):
    import TestProblems as tp
    
    if testName is None:
        print("NO TEST NAME GIVEN")
    else:
        coord, density, pressure, energy = tp.set_problems(testName, nx, gammaFactor, x_l, x_u)
    return(coord, density, pressure, energy)

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

def HydroCalculateRemain(density, pressure, Z, massNumber, gamma_factor):

    #Initialise cell centered variables
    numberDensityI = density / (massNumber * protonMass)
    numberDensityE = numberDensityI * Z
    
    ##Derived from .. Total Pressure = Electron pressure + ion pressure = nekbTe + nikbTi  
    # Assuming Te = Ti initally  and ne = ni * Z 
    # Total Pressure = (1 + Z)nekbTe
    temperatureE = pressure / ((1 + Z) * numberDensityE * kb)
    temperatureI = temperatureE

    ##Noting convention that Specific Heat ha
    specificHeatE =  (numberDensityE * kb) / (density * (gamma_factor -1))
    DpDTe = numberDensityE *kb
    
    specificHeatI =  (numberDensityI * kb) / (density * (gamma_factor -1))
    DpDTi = numberDensityI *kb

    #Pressure
    pressureI = numberDensityI * kb * temperatureI
    pressureE = numberDensityE * kb * temperatureE
    pressureTotal = pressureE + pressureI
    #Internal Energy

    IntEe = pressureE / ((gamma_factor - 1) * density)
    IntEi = pressureI / ((gamma_factor - 1) * density)


    ConstCheck(density, numberDensityE, numberDensityI, temperatureE, temperatureI, pressureE, pressureI, pressureTotal, massNumber, Z)

    return(numberDensityE, numberDensityI, temperatureE, temperatureI, specificHeatE, DpDTe, specificHeatI, DpDTi, pressureE, pressureI, pressureTotal, IntEe, IntEi)

def TemperatureCalculateRemain(numberDensityE, temperatureE, Z, massNumber, gammaFactor, temperatureI = None):
    
    numberDensityI = numberDensityE / Z
    density = numberDensityI * massNumber * protonMass

    if temperatureI is None:
        temperatureI = temperatureE
    
    pressureTotal = temperatureE * (1+ Z) * numberDensityE * kb
    ##Noting convention that Specific Heat ha
    specificHeatE =  (numberDensityE * kb) / (density * (gammaFactor -1))
    DpDTe = numberDensityE *kb
    
    specificHeatI =  (numberDensityI * kb) / (density * (gammaFactor -1))
    DpDTi = numberDensityI *kb

    pressureE = pressureTotal / 2
    pressureI = pressureE

    IntEe = pressureE / ((gammaFactor - 1) * density)
    IntEi = pressureI / ((gammaFactor - 1) * density)

    ConstCheck(density, numberDensityE, numberDensityI, temperatureE, temperatureI, pressureE, pressureI, pressureTotal, massNumber, Z)
    return(density, numberDensityE, numberDensityI, temperatureE, temperatureI, specificHeatE, DpDTe, specificHeatI, DpDTi, pressureE, pressureI, pressureTotal, IntEe, IntEi)

def TwoTemperatureBath(x_l, x_u, nx):
    
    
    #Two bath problem
    gammaFactor = 5/3
    Tu = 30
    Td = 3
    Te = (Td**3.5 - (np.linspace(x_l, x_u, nx) / x_u) * (Td**3.5 - Tu**3.5))**(1/3.5)
    Te[0] = Td
    Te[-1] = Tu
    Ti = Te
    nu = 10E21
    ne = nu * (Tu/Te)
    ni = ne 
    density = ni * 1.67E-27
    pe = ne * kb * Te
    pi = pe
    pTotal = pe + pi
    ConstCheck(density, ne, ni, Te,Ti, pe, pi, pTotal, 1, 1)
    specificHeatE =  (ne * kb) / (density * (gammaFactor -1))
    DpDTe = ne *kb
    
    specificHeatI =  (ni * kb) / (density * (gammaFactor -1))
    DpDTi = ni *kb

    IntEe = pe / ((gammaFactor - 1) * density)
    IntEi = pi / ((gammaFactor - 1) * density)


    
    return(density, numberDensityE, numberDensityI, temperatureE, temperatureI, specificHeatE, DpDTe, specificHeatI, DpDTi, pressureE, pressureI, pressureTotal, IntEe, IntEi)
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

def InvBremCheck(ne, nx):
    
    
    #Two bath problem
    gammaFactor = 5/3
    Te = np.zeros(nx) + 3
    Ti = Te
    nu = ne
    ne = np.zeros(nx) + nu
    ni = ne 
    density = ni * 1.67E-27
    pe = ne * kb * Te
    pi = pe
    pTotal = pe + pi
    ConstCheck(density, ne, ni, Te,Ti, pe, pi, pTotal, 1, 1)
    specificHeatE =  (ne * kb) / (density * (gammaFactor -1))
    DpDTe = ne *kb
    
    specificHeatI =  (ni * kb) / (density * (gammaFactor -1))
    DpDTi = ni *kb

    IntEe = pe / ((gammaFactor - 1) * density)
    IntEi = pi / ((gammaFactor - 1) * density)


    
    return(density, ne, ni, Te, Ti, specificHeatE, DpDTe, specificHeatI, DpDTi, pe, pi, pTotal, IntEe, IntEi)

def BremCheck(ne, nx):
    
    
    #Two bath problem
    gammaFactor = 5/3
    Te = np.zeros(nx) + 300 * 11600
    Ti = Te
    nu = ne
    ne = np.zeros(nx) + nu
    ni = ne 
    density = ni * 1.67E-27
    pe = ne * kb * Te
    pi = pe
    pTotal = pe + pi
    ConstCheck(density, ne, ni, Te,Ti, pe, pi, pTotal, 1, 1)
    specificHeatE =  (ne * kb) / (density * (gammaFactor -1))
    DpDTe = ne *kb
    
    specificHeatI =  (ni * kb) / (density * (gammaFactor -1))
    DpDTi = ni *kb

    IntEe = pe / ((gammaFactor - 1) * density)
    IntEi = pi / ((gammaFactor - 1) * density)


    
    return(density, ne, ni, Te, Ti, specificHeatE, DpDTe, specificHeatI, DpDTi, pe, pi, pTotal, IntEe, IntEi)

def TextDump(path, coord, velocity, density, ne, ni, Te, Ti, Pe, Pi, Ptot, IntEe, IntEi, DpDTe, DpDTi, Cve,Cvi, mass, invBrem, absorption, brem):
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
    np.savetxt((path+"alpha.txt"), absorption)
    np.savetxt((path+"brem.txt"), brem)


#### Custom routines
def custom_routine(L, nx, ne, temperature, gamma, Z, massNumber):
    x_l = 0
    x_u = L
    initial_coord = initial_coord = np.linspace(x_l, x_u, nx+1)
    ne = ne + np.zeros(nx)
    ni = ne / Z
    density = massNumber * 1.67E-27 * ni
    np_centered_x = np.array([0.5*(initial_coord[i + 1] + initial_coord[i]) for i in range(len(initial_coord) - 1)])
    dg = np.float(0.9 * (x_u - x_l))
    mid_point = np.float(0.5) * (np.float(x_u - x_l))
    #temperatureE = temperature +  np.float(3000) * np.exp(-(((np_centered_x -  mid_point)**2) / (dg**2)), dtype = float)
    Bz = 0
    normalised_values = norm.impact_inputs(np.average(ne) * 1e-6, temperature, Z, massNumber, Bz)
    temperatureE = 0.5 + profiles.load_profile(nx = nx,
                                                        xmin = 0,
                                                        xmax = 600,
                                                        avg = 0.5,
                                                        amp = 0.0005,
                                                        pos = 0.0,
                                                        nwl = 0.50,
                                                        wid = 1000.0,
                                                        func = '+cos'
                                                        )
                                   
    temperatureE = temperatureE * normalised_values['Te']
   # temperatureE = temperature + np.linspace(11600*10, temperature, nx)
    temperatureI = temperatureE
    pressureE = ne * 1.38E-23 * temperatureE
    pressureI = ni * 1.38E-23 * temperatureI
    pressureTotal  = pressureE + pressureI
    ConstCheck(density, ne, ni, temperatureE,temperatureI, pressureE, pressureI, pressureTotal, massNumber, Z)

    specificHeatE =  (ne * kb) / (density * (gammaFactor -1))
    DpDTe = ne *kb
    
    specificHeatI =  (ni * kb) / (density * (gammaFactor -1))
    DpDTi = ni *kb

    IntEe = pressureE / ((gammaFactor - 1) * density)
    IntEi = pressureI / ((gammaFactor - 1) * density)

    return(initial_coord, density, ne, ni, temperatureE, temperatureI, specificHeatE, DpDTe, specificHeatI, DpDTi, pressureE, pressureI, pressureTotal, IntEe, IntEi)

nx = 100
x_l = 0
x_u = 5.44759e-07  * 600
L = x_u - x_l
massNumber = 157
Z = 64
testName = "hydro_energy_diff"
gammaFactor = 1.4
laserWavelength = 1E-9
LaserPower = 0
coulombLog = 11
#Ev
temperature = 1000

ne = 1E27
nc = 1.1E15 / pow(laserWavelength, 2)

velocity = np.zeros(nx + 1) #+ add any function
#path = "/Users/shiki/Documents/Imperial_College_London/Ph.D./Lagrangian_fluid_code/c_plus_plus_lagrangian/init_data/"
path = "/home/abetharan/HeadlessHydra/init_data/"
coord, density, numberDensityE, numberDensityI, temperatureE, temperatureI, specificHeatE, DpDTe, specificHeatI, DpDTi, pressureE, pressureI, pressureTotal, IntEe, IntEi  = custom_routine(L, nx, ne, temperature, gammaFactor, Z, massNumber)

mass = CalculateMass(coord, density, nx)
InvBrem_, absorption  = InvBrem(coord, numberDensityE, nc, laserWavelength, Z, coulombLog, temperatureE, LaserPower, mass, nx)
Brem = Brem(numberDensityE, massNumber, 0, temperatureE, nx)
TextDump(   path = path,
            coord= coord ,
            velocity = velocity,
            density = density,
            ne = numberDensityE,
            ni = numberDensityI ,
            Te = temperatureE,
            Ti = temperatureI,
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
            absorption = absorption,
            brem = Brem
        )

import matplotlib.pyplot as plt

plt.figure(1)
normalised_values = norm.impact_inputs(np.average(ne)*1e-6, temperature, Z,massNumber, 0)
plotTemperature = (temperatureE * (kb/e))/ normalised_values['Te']
plt.plot(plotTemperature)
#lt.plot(temperatureI)
plt.title("temperature")
plt.figure(2)
plt.plot(density)
plt.title("n density")

plt.figure(3)
plt.plot(pressureE)
plt.plot(pressureI)
plt.title("pressure")

plt.show()
