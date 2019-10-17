""" Only Valid for idle gases 
    By default initial velocity is always zero
"""
import numpy as np
import impact_profiles_py3 as prof
import impact_norms_py3 as norm
import matplotlib.pyplot as plt

#global constants
kb = 1.38064852e-23
protonMass = 1.6726219e-27
electronMass = 9.10938356e-31
epsilon_0 = 8.85418782e-12
electronCharge = 1.60217662e-19
hbar = 1.054571800e-34
c = 299792458;

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

def TextDump(path, coord, velocity, density, ne, ni, Te, Ti, mass, Z, Ar):
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
    #np.savetxt((path+"electron_temperature.txt"), Te)
    #np.savetxt((path+"ion_temperature.txt"), Ti)
    np.savetxt((path+"mass.txt"), mass)
    np.savetxt((path+"Z.txt"), Z)    
    np.savetxt((path+"Ar.txt"), Ar)
    

#### Custom routines
def custom_routine(L, nx, nc, temperature, gamma, Z, massNumber):
    x_l = 0
    x_u = L
    initial_coord = initial_coord = np.linspace(x_l, x_u, nx+1)
    index_for_high_z = [i for i in range(len(initial_coord)) if initial_coord[i] > 0.06e-3]
    Z = np.zeros(nx) + 64
    Ar = np.zeros(nx) + 157
    # Z[:index_for_high_z[0]] = 40
    # Z[index_for_high_z[0]:] = 2
    # Ar[:index_for_high_z[0]] = 157
    # Ar[index_for_high_z[0]:] = 4
    
    # density = np.zeros(nx)
    # density[:index_for_high_z[0]] = 19.3
    # density[index_for_high_z[0]:] = 0.3
    
    ne = np.zeros(nx) + 1e26
    # ne[:index_for_high_z[0]] = 10 * nc
    # ne[index_for_high_z[0]:] = 9e25

    # Z[0:int(nx*0.7)] = 2
    # Z[int(nx*0.7):] = 64
    # Ar = np.zeros(nx)
    # Ar[0:int(nx*0.7)] = 4
    # Ar[int(nx*0.7):] = 157.25
    ni = ne / Z
    density = Ar * ni * protonMass
    temperatureE = np.zeros(nx) + 11604 * 5e-3
    # temperatureE = temperature * prof.load_profile(nx = nx,
    #                                                 xmin = 0,
    #                                                 xmax = 1000,
    #                                                 avg = 1.0,
    #                                                 amp = 0.0005,
    #                                                 pos = 0.0,
    #                                                 nwl = 0.5,
    #                                                 wid = 1000.0,
    #                                                 func = '+cos'
    #                                                     )

   # temperatureE = temperature + np.linspace(11600*10, temperature, nx)



    temperatureI = 0

    return(initial_coord, density, ne, ni, temperatureE, temperatureI, Z, Ar)


nx = 120
x_l = 0
x_u = 600 * 5.75919e-07
L = x_u - x_l
Ar = np.zeros(nx) + 157
Z = np.zeros(nx) + 64
testName = "hydro_energy_diff"
gammaFactor = 1.4
laserWavelength = 188E-9
LaserPower = 0
coulombLog = 11
#Ev
temperature = 300 * 11604


ne = 1E26
nc = ne#1.1E15 / pow(laserWavelength, 2)

velocity = np.zeros(nx + 1) #+ add any function
#path = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/init_data/"
path = "/home/abetharan/HeadlessHydra/init_data/"
coord, density, numberDensityE, numberDensityI, temperatureE, temperatureI, Z, Ar  = custom_routine(L, nx, nc, temperature, gammaFactor, Z, Ar)
mass = CalculateMass(coord, density, nx)

plt.figure(1)
plt.plot(temperatureE)
#plt.plot(temperatureI)
plt.title("temperature")
plt.figure(2)
plt.plot(density)
plt.title("n density")

plt.show()

TextDump(   path = path,
            coord= coord ,
            velocity = velocity,
            density = density,
            ne = numberDensityE,
            ni = numberDensityI ,
            Te = temperatureE,
            Ti = temperatureI,
            mass = mass,
            Z = Z,
            Ar = Ar,
)


