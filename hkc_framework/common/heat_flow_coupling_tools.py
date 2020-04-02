import numpy as np 
from scipy import constants
kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
epsilon_0 = 8.854188E-12    # Vacuum dielectric constant
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")

def lambda_ei(self, n, T, T_norm = 10, n_norm = 0.25E20, Z_norm = 1.0):
    if T * T_norm < 10.00 * Z_norm ** 2:

        result = 23.00 - math.log(math.sqrt(n * n_norm * 1.00E-6) * Z_norm * (T * T_norm) ** (-3.00/2.00))

    else:

        result = 24.00 - math.log(math.sqrt(n * n_norm * 1.00E-6) / (T * T_norm))   

    return result
def spitzerHarmHeatFlow(cell_centre_coord, temperature, zbar, number_density):
    """
    Purpose: Models Spitzer Harm heat flow 
    Args:
        cell_centre_coord = List of cell centres
        temperature = electron temperature.
        zbar = ionization
        number_density = number density of electrons
    Returns:
        SH_Heat_Flow = Spitzer harm heat flow
    """
    coulomb_log = []
    for i in range(len(zbar)):
        coulomb = self.lambda_ei(1.0, 1.0, T_norm = temperature[i] * (kb/e), n_norm = number_density[i], Z_norm = zbar[i])
        coulomb_log.append(coulomb)

    coulomb_log = np.array(coulomb_log)
    kappaE = 1.843076667547614E-10 * pow(temperature, 2.5) * pow(coulomb_log, -1) *  pow(zbar, -1)
    cell_centre_coord = np.array([(cell_centre_coord[i+1] + cell_centre_coord[i]) / 2 for i in range(len(cell_centre_coord) - 1)]) 
    nx = len(temperature)
    HeatFlowE = np.zeros(nx + 1)

    for i in range(1, nx):
        centered_ke = 0.5 * (kappaE[i] + kappaE[i - 1])
        HeatFlowE[i] = centered_ke *((temperature[i] - temperature[i - 1]) / (cell_centre_coord[i] - cell_centre_coord[i - 1]))
        
    HeatFlowE[0] = 0
    HeatFlowE[-1] = 0
    return(-1*HeatFlowE)

def preHeatModel(coord_non_local, local_heat, non_local_heat, tolerance):
    """ 
    Purpose: Model Pre-heat using an exponential fitting parameter, fit parameter
                is spat out to be used by the fluid code.
    Args:
        coord = cell wall grid
        local_heat = Spitzer-Harm heat flow
        non_local_heat = Kinetic code heat flow
    Returns:
        B = Fitting paramters
        preheat_start = start index for preheat
        preheat_end = end index of preheat
    """
    #interp_local_heat = np.interp(coord_non_local, coord_local, local_heat)
    intersect = local_heat / non_local_heat
    q_q_sh = non_local_heat / local_heat
    intersect[np.isnan(intersect)] = 0 
    preheat_start = np.where(intersect[~np.isnan(intersect)] != 0)[0][-1]
    preheat_end = np.where(abs(non_local_heat[preheat_start:]) < tolerance)[0][0] + preheat_start
    #preheat_end_value = 1e-60
    L = coord_non_local[preheat_end] - coord_non_local[preheat_start] 
    #B =  -1/np.log(preheat_end_value / non_local_heat[preheat_start])
    B = []
    for i in range(preheat_start, preheat_end):
        b = -1* (coord_non_local[i] - coord_non_local[preheat_start]) / (L * np.log(abs(non_local_heat[i])/non_local_heat[preheat_start]))
        B.append(b)
        
    B = np.array(B)
    return(preheat_start, preheat_end, B)

def frontHeatModel(coord_non_local, local_heat, non_local_heat, tolerance):
    """
    Purpose: Model heat wave from top of a bath.
    Args:
        coord = cell wall grid
        local_heat = Spitzer-Harm heat flow
        non_local_heat = Kinetic code heat flow
    Returns:
        B = Fitting paramters
        frontheat_start = start of front heat
        frontheat_end = end of front 
        heat
    """

    #interp_local_heat = np.interp(coord_non_local, coord_local, local_heat)
    intersect = local_heat / non_local_heat
    intersect[np.isnan(intersect)] = 0 
    frontheat_start = np.where(intersect != 0)[0][0]  
    frontheat_end = np.where(abs(non_local_heat[:frontheat_start]) < tolerance)[0][0]
    L = abs(coord_non_local[frontheat_end] - coord_non_local[frontheat_start])
    
    B = []

    for i in range(0, frontheat_start+1):
        b = -1* (coord_non_local[i] - coord_non_local[frontheat_start]) / (L * np.log(abs(non_local_heat[i])/non_local_heat[frontheat_start]))
        B.append(b)
    
    B = np.array(B)
    return(frontheat_start, frontheat_end, B)

def divQHeatFlow(electron_thermal_flux, mass):
    """ Purpose: Find Div.Q
        Args:
            electron_thermal_flux = SOL-KiT heat flux in SI
            mass = areal mass.
            
        Returns: Div.Q
        NOTE: ONLY WORKS FOR SAME GRIDS
    """
    nx = len(electron_thermal_flux)
    HeatConductionE = np.zeros(len(mass))
    #Include first and celll wall for no flow set to 0 per no thermal influx BC = 0 
    if self.boundary_condition =='noflow':
        electron_thermal_flux = np.insert(electron_thermal_flux, 0, 0.)
    #Periodic boundary condition removes last cell wall insert back and set to 0 via thermal influc BC
    electron_thermal_flux = np.append(electron_thermal_flux, 0.)
    #Heat flow from SOL-KiT includes celll walla and cell centre.
    #Thus, step in 2s for 1 to 1 grid. Change if higher resolution grid of SOL-KiT is used.
    # |.|.|.| 
    step = 2
    j = 0
    for i in range(0, nx, step):   
        HeatConductionE[j] = -(electron_thermal_flux[i + 2] - electron_thermal_flux[i]) / mass[j]
        j += 1
    return(HeatConductionE)
