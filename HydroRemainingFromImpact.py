import numpy as np
protonMass = 1.66e-27


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

def CalculateRemain(fluid_out_path, ne, Te, qe):
    
