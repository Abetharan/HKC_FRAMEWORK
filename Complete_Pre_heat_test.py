import numpy as np
from scipy import constants
import matplotlib.pyplot as plt  
import math
kB = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")
epsilon_0 = epsilon_0 = 8.854188E-12    # Vacuum dielectric constant

def mCalculateKappaE(TemperatureE, CoulombLogEI, Zbar):

    ke = 1.843076667547614E-10 * pow(TemperatureE, 2.5) * pow(CoulombLogEI, -1) *  pow(Zbar, -1)
    return(ke)

def HeatFlowE(CellCenteredCoord, kappaE, TemperatureE, NumberDensityE):
    nx = len(TemperatureE)
    HeatFlowE = np.zeros(nx + 1)
    for i in range(1, nx):
    
        centered_ke = 0.5 * (kappaE[i] + kappaE[i - 1])
        HeatFlowE[i] = centered_ke *((TemperatureE[i] - TemperatureE[i - 1]) / (CellCenteredCoord[i] - CellCenteredCoord[i - 1]))
        
        #e_flux_lim =   0.05 * pow(me, -0.5) * pow(kB, 1.5) * NumberDensityE[i] * pow(TemperatureE[i], 1.5)
        #HeatFlowE[i] = min([e_flux, e_flux_lim])
    HeatFlowE[0] = 0
    HeatFlowE[-1] = 0
    return(HeatFlowE)

def normlisation(norm_Te, norm_ne, norm_Z):
    def lambda_ei(n, T, T_norm = 10, n_norm = 0.25E20, Z_norm = 1.0):
        if T * T_norm < 10.00 * Z_norm ** 2:

            result = 23.00 - math.log(math.sqrt(n * n_norm * 1.00E-6) * Z_norm * (T * T_norm) ** (-3.00/2.00))

        else:

            result = 24.00 - math.log(math.sqrt(n * n_norm * 1.00E-6) / (T * T_norm))

        return result


    #9.10938E-31            # Electron mass in kg
    el_charge = e #1.602189E-19    # Electron charge
    r_b = bohr_radi#5.29E-11              # Bohr radius in metres

    sigma_0 = math.pi * r_b ** 2  # Cross-section normaliation
    gamma_ee_0 = el_charge ** 4 / (4 * math.pi * (me * epsilon_0) ** 2)

    T_J = norm_Te * el_charge

    gamma_ei_0 = norm_Z ** 2 * gamma_ee_0

    v_t = math.sqrt(2.0 * T_J / me)                  # Velocity normalization
    coulomb = lambda_ei(1.0, 1.0, T_norm = norm_Te, n_norm = norm_ne, Z_norm = norm_Z)
    time_norm = v_t ** 3 / (gamma_ei_0 * (norm_ne/norm_Z) * coulomb) # Time normalization
    x_0 = v_t * time_norm               # Space normalization

    e_0 = me * v_t / (el_charge * time_norm) # Electric field normalization
    q_0 = me * norm_ne * (v_t ** 3)
    print(time_norm)
    dict = {}
    print(coulomb)
    dict['lambda_mfp'] = x_0
    dict['vth'] = v_t
    dict['coulomb_log'] = coulomb
    dict['ne'] = norm_ne
    dict['Te'] = norm_Te
    dict['Z'] = norm_Z
    dict['tau_ei'] = time_norm
    # convert ne to 10**21 cm**-3
    dict['qe'] = q_0
    return(dict)

def PreHeatModel(coord_non_local, local_heat, non_local_heat):
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
    
    """
    #interp_local_heat = np.interp(coord_non_local, coord_local, local_heat)
    intersect = local_heat / non_local_heat
    q_q_sh = non_local_heat / local_heat
    intersect[np.isnan(intersect)] = 0 
    preheat_start = np.where(intersect[~np.isnan(intersect)] != 0)[0][-1]
    tolerance = 0.01 * norms['tau_ei']
    preheat_end = np.where(abs(non_local_heat[preheat_start:]) < tolerance)[0][0] + preheat_start
    #preheat_end_value = 1e-60
    L = coord_non_local[preheat_end] - coord_non_local[preheat_start] 
    #B =  -1/np.log(preheat_end_value / non_local_heat[preheat_start])
    B = []
    for i in range(preheat_start, preheat_end):
       b = -1* (coord_non_local[i] - coord_non_local[preheat_start]) / (L * np.log(non_local_heat[i]/non_local_heat[preheat_start]))
       B.append(b)
    
    B = np.array(B)
    
    #HeatFlow = q_q_sh[preheat_start] *interp_local_heat[preheat_start] * np.exp(-1*(coord_non_local[preheat_start:(preheat_end + 1)] - coord_non_local[preheat_start]) / (B*L))
    HeatFlow = q_q_sh[preheat_start] *local_heat[preheat_start] * np.exp(-1*(((coord_non_local[preheat_start:preheat_end] - coord_non_local[preheat_start]) / (B*L))))
    HeatFlow = np.concatenate((HeatFlow,np.zeros(len(coord_non_local[preheat_start:]) - len(HeatFlow))))
    return(preheat_start,-1, coord_non_local[preheat_start:], HeatFlow)

def FrontHeat(coord_non_local, local_heat, non_local_heat):
    """
    Purpose: Model heat wave from top of a bath.
    Args:
        coord = cell wall grid
        local_heat = Spitzer-Harm heat flow
        non_local_heat = Kinetic code heat flow
    Returns:
        B = Fitting paramters
        frontheat_end = start index for preheat    
    """

    #interp_local_heat = np.interp(coord_non_local, coord_local, local_heat)
    intersect = local_heat / non_local_heat
    intersect[np.isnan(intersect)] = 0 
    q_q_sh = non_local_heat / local_heat
    frontheat_start = np.where(intersect != 0)[0][0]  
    frontheat_end = 0
    L = abs(coord_non_local[frontheat_end] - coord_non_local[frontheat_start])
    
    B = []

    for i in range(0, frontheat_start+1):
       b = -1* (coord_non_local[i] - coord_non_local[frontheat_start]) / (L * np.log(non_local_heat[i]/non_local_heat[frontheat_start]))
       B.append(b)
    
    B = np.array(B)
    HeatFlow = (q_q_sh[frontheat_start] * local_heat[frontheat_start]) * np.exp( -1* (coord_non_local[:(frontheat_start + 1)] - coord_non_local[frontheat_start]) / (B * L))
    
    HeatFlow[0] = 0
    HeatFlow[-1] = 0

    return(frontheat_start, frontheat_end, coord_non_local[:(frontheat_start + 1)], HeatFlow)


def ThermalConduc(HeatFlowE, mass):
    nx = len(HeatFlowE) 
    HeatConductionE = np.zeros(nx)
    if HeatFlowE[0] != 0:
        HeatFlowE[0] = 0
    #HeatFlowE = np.append(HeatFlowE, 0)
    #HeatFlowE = np.insert(HeatFlowE, 0,0)
    j = 0
    for i in range(0, nx - 1):
        HeatConductionE[i] = (HeatFlowE[i + 1] - HeatFlowE[i]) / mass[j];
        j+=1
    return(HeatConductionE)

def LocalHeat(path, index, norms):
    # coord = np.loadtxt(path + '/OUTPUT/GRIDS/X_GRID.txt')
    # Te = np.loadtxt(path + '/OUTPUT/TEMPERATURE/TEMPERATURE_' + index + '.txt')
    # ne = np.loadtxt(path + '/OUTPUT/DENSITY/DENSITY_' + index +'.txt')
    # zbar = np.loadtxt(path + '/INPUT/Z_PROFILE_INPUT.txt')
    # kappaE = mCalculateKappaE(Te, norms['coulomb_log'], zbar)
    # local_heat_flow = HeatFlowE(coord, kappaE, Te, ne)
    
    coord = np.loadtxt(path + '/coord.txt')
    Te = np.loadtxt(path + '/electron_temperature.txt')
    density = np.loadtxt(path + '/density.txt')
    Ar = np.loadtxt(path + '/Ar.txt')
    zbar = np.loadtxt(path + '/Z.txt')
    ni = density/(Ar * mp)
    ne = ni * zbar
    kappaE = mCalculateKappaE(Te, norms['coulomb_log'], zbar)
    local_heat_flow = HeatFlowE(coord, kappaE, Te, ne) / norms['qe']
    
    return(coord, -1* local_heat_flow)


norms = normlisation(3500, 1e27, 36.5)
kinetic_heat_profile = np.loadtxt('/media/abetharan/DATADRIVE2/Abetharan/pre_heat_model_test_problem/OUTPUT/HEAT_FLOW_X/HEAT_FLOW_X_03000.txt')
multiplier = np.loadtxt('/media/abetharan/DATADRIVE2/Abetharan/pre_heat_model_test_problem/OUTPUT/SH_q_ratio/SH_q_ratio_03000.txt')
coord, sh_heat = LocalHeat('/home/abetharan/HYDRO_KINETIC_COUPLING/Non_Linear_Ramp_Investigation/Pre_Heat_Ramp/', '03000', norms)
mass = np.loadtxt('/home/abetharan/HYDRO_KINETIC_COUPLING/Non_Linear_Ramp_Investigation/Pre_Heat_Ramp/mass.txt')
#coord_sol = np.loadtxt('/media/abetharan/DATADRIVE2/Abetharan/pre_heat_model_test_problem/OUTPUT/GRIDS/X_GRID.txt')
#non_local_qe = np.loadtxt('/Users/shiki/Documents/Imperial_College_London/Ph.D./HYDRO_IMPACT_COUPLING/HEAT_FLOW_X_03000.txt')

#coord, sh_heat = LocalHeat('/Users/shiki/Documents/Imperial_College_London/Ph.D./HYDRO_IMPACT_COUPLING/Pre_Heat_Ramp/', '03000',norms )
#coord_sol = np.loadtxt('/Users/shiki/Documents/Imperial_College_London/Ph.D./HYDRO_IMPACT_COUPLING/X_GRID.txt')
#mass = np.loadtxt('/Users/shiki/Documents/Imperial_College_London/Ph.D./HYDRO_IMPACT_COUPLING/Pre_Heat_Ramp/mass.txt')


#multiplier = np.loadtxt('/Users/shiki/Documents/Imperial_College_London/Ph.D./HYDRO_IMPACT_COUPLING/SH_q_ratio_03000.txt')
#inetic_heat_profile  = np.loadtxt('/Users/shiki/Documents/Imperial_College_London/Ph.D./HYDRO_IMPACT_COUPLING/HEAT_FLOW_X_03000.txt')

#Include first and celll wall for no flow set to 0 per no thermal influx BC = 0 

kinetic_heat_profile = np.insert(kinetic_heat_profile, 0, 0.)
multiplier = np.insert(kinetic_heat_profile,0 ,0)
#Periodic boundary condition removes last cell wall insert back and set to 0 via thermal influc BC
kinetic_heat_profile = np.append(kinetic_heat_profile, 0.)
multiplier = np.append(multiplier, 0)
cell_wall_list = coord
##kinetic_grid = np.loadtxt("/Users/shiki/Documents/Imperial_College_London/Ph.D./HYDRO_IMPACT_COUPLING/X_GRID.txt")
#kinetic_grid = np.insert(kinetic_grid, 0, 0)
#kinetic_grid = np.append()

multi_list = []
qe_list = []
#cell_wall_list = []
step = 2

#Get rid of cell-centre quantites

for i in range(0, len(multiplier), step):   
    #add multiplier limiter here if needed
    multi_list.append(multiplier[i])
    qe_list.append(kinetic_heat_profile[i])
    #cell_wall_list.append(kinetic_grid[i])

##Test for pre-heat via looking at NaN outputs expected from q/q_sh
#if nansif(any(np.isnan(multi_list))):

#     multi_list = np.array(multi_list)
#     qe_list = np.array(qe_list)
#     NaN_args = np.argwhere(np.isnan(multi_list))

#     pre_heat_start_index = NaN_args[0]
#     #Determine length to tolerance
#     pre_heat_last_index = np.argwhere(qe_list[pre_heat_start_index:] <= 1e-14)[0] #last point of pre heat

# else:
#     pre_heat_start_index = 0
#     pre_heat_last_index = 0 

multi_list = qe_list / sh_heat
qe_list = np.array(qe_list)


preheat_start, preheat_end, coord_pre_heat, preheat = PreHeatModel(cell_wall_list, sh_heat, qe_list, norms)
frontheat_start, frontheat_end, coord_front_heat, frontheat  = FrontHeat(cell_wall_list, sh_heat, qe_list)

heat_flow = multi_list * sh_heat
combined_flow = np.concatenate((frontheat[:-1],heat_flow[~np.isnan(heat_flow)], preheat[1:]))

combined_thermal_conduc = ThermalConduc(combined_flow, mass)
kinetic_thermal_conduc = ThermalConduc(qe_list, mass)

plt.plot(combined_thermal_conduc, 'k-', label = 'reci_front_heat')
plt.plot(kinetic_thermal_conduc, 'r--')
#plt.plot(qe_list, 'k-', label = 'sol_kit')
#plt.plot(combined_flow, 'r-', label = 'reci_front-heat')
#plt.plot(cell_wall_list, sh_heat, 'b--')
#plt.ylim(0,0.0017)
#plt.plot(cell_wall_list, heat_flow, 'k-')
#plt.plot(coord_front_heat, frontheat, 'r-')
#plt.plot(coord_pre_heat, preheat, 'b-')
plt.legend()
plt.show()