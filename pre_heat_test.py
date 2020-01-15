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

def NonLocalHeatFlow(path, norms):
    norm_qe = norms['vth']**3 * norms['ne'] * me
    qe = np.loadtxt(path) * norm_qe
    return(qe)

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
    local_heat_flow = HeatFlowE(coord, kappaE, Te, ne)
    
    return(coord, local_heat_flow)

def LocalHeatSolKit(path, norms):
    coord = np.loadtxt(path + '/coord.txt')
    Te = np.loadtxt(path + '/electron_temperature.txt')
    density = np.loadtxt(path + '/density.txt')
    zbar = np.loadtxt(path + '/Z.txt')
    Ar = np.loadtxt(path + '/Ar.txt')
    ni = density/(Ar * mp)
    ne = ni * zbar
    kappaE = mCalculateKappaE(Te, norms['coulomb_log'], zbar)
    local_heat_flow = HeatFlowE(coord, kappaE, Te, ne)
    return(coord, local_heat_flow)

def Alt2NonLocalPreHeatModel(coord_local, coord_non_local, local_heat, non_local_heat, norms):
    from scipy.interpolate import CubicSpline, interp1d

    interp_local_heat = np.interp(coord_non_local, coord_local, local_heat)
    intersect = interp_local_heat / non_local_heat
    q_q_sh = non_local_heat / interp_local_heat
    preheat_start = np.where(intersect != 0)[0][-1]
    tolerance = 0.01 * norms['tau_ei']
    preheat_end = np.where(abs(non_local_heat[preheat_start:]) < tolerance)[0][-1] + preheat_start
    preheat_end_value = 1e-60
    L = coord_non_local[preheat_end] - coord_non_local[preheat_start] 
    B =  -1/np.log(preheat_end_value / non_local_heat[preheat_start])
    
    #HeatFlow = q_q_sh[preheat_start] *interp_local_heat[preheat_start] * np.exp(-1*(coord_non_local[preheat_start:(preheat_end + 1)] - coord_non_local[preheat_start]) / (B*L))
    HeatFlow = q_q_sh[preheat_start] *interp_local_heat[preheat_start] * np.exp(-1*(((coord_non_local[preheat_start:(preheat_end + 1)] - coord_non_local[preheat_start]) / (B*L))))
    HeatFlow = np.concatenate((HeatFlow,np.zeros(len(coord_non_local[preheat_start:]) - len(HeatFlow))))
    return(preheat_start,-1, coord_non_local[preheat_start:] * norms['lambda_mfp'], HeatFlow * norms['qe'], interp_local_heat * norms['qe'])

def AltNonLocalPreHeatModel(coord_local, coord_non_local, local_heat, non_local_heat):
    from scipy.interpolate import CubicSpline, interp1d

    interp_local_heat = np.interp(coord_non_local, coord_local, local_heat)
    intersect = interp_local_heat / non_local_heat
    q_q_sh = non_local_qe / interp_local_heat
    preheat_start = np.where(intersect != 0)[0][-1]
#    preheat_end = np.where(abs(non_local_heat[preheat_start:]) < 1e2)[0][-1] + preheat_start
    preheat_end = 1e-13
    non_local_heat_last = non_local_heat[-1]
    L = coord_non_local[-1] - coord_non_local[preheat_start] 
    B =  -1 / (np.log((preheat_end - non_local_heat_last )/ non_local_heat[preheat_start]))
    HeatFlow = q_q_sh[preheat_start] *interp_local_heat[preheat_start] * np.exp(-1*(coord_non_local[preheat_start:] - coord_non_local[preheat_start]) / (B * L)) + non_local_heat_last 
    
    return(preheat_start,-1, coord_non_local[preheat_start:], HeatFlow, interp_local_heat)

def NonLocalPreHeatModel(coord_local, coord_non_local, local_heat, non_local_heat):
    from scipy.interpolate import CubicSpline, interp1d

    interp_local_heat = np.interp(coord_non_local, coord_local, local_heat)
    intersect = interp_local_heat / non_local_heat
    q_q_sh = non_local_qe / interp_local_heat
    preheat_start = np.where(intersect != 0)[0][-1]
#    preheat_end = np.where(abs(non_local_heat[preheat_start:]) < 1e2)[0][-1] + preheat_start
    preheat_end = 1e-13
    non_local_heat_last = non_local_heat[-1]
    L = coord_non_local[-1] - coord_non_local[preheat_start] 
    B =  -1 / (np.log((preheat_end)/ non_local_heat[preheat_start]))
    HeatFlow = q_q_sh[preheat_start] *interp_local_heat[preheat_start] * np.exp(-1*(coord_non_local[preheat_start:] - coord_non_local[preheat_start]) / (B * L)) 
    return(preheat_start,-1, coord_non_local[preheat_start:], HeatFlow, interp_local_heat)

def FrontHeat(coord_local, coord_non_local, local_heat, non_local_heat):

    from scipy.interpolate import CubicSpline, interp1d
    
    interp_local_heat = np.interp(coord_non_local, coord_local, local_heat)
    intersect = interp_local_heat / non_local_heat
    q_q_sh = non_local_qe / interp_local_heat
    frontheat_start = np.where(intersect != 0)[0][0]
    frontheat_end = 0
    L = abs(coord_non_local[frontheat_end] - coord_non_local[frontheat_start])
    B = 1 / (np.log((non_local_heat[frontheat_start] - non_local_heat[0]) / 1e-5))
    HeatFlow = q_q_sh[frontheat_start] * interp_local_heat[frontheat_start] * np.exp((coord_non_local[:(frontheat_start + 1)] - coord_non_local[frontheat_start]) / (B * L)) + non_local_heat[0]

    return(frontheat_start, frontheat_end, coord_non_local[:(frontheat_start + 1)], HeatFlow)

norms = normlisation(3500, 1e27, 36.5)
non_local_qe = NonLocalHeatFlow('/media/abetharan/DATADRIVE2/Abetharan/pre_heat_model_test_problem/OUTPUT/HEAT_FLOW_X/HEAT_FLOW_X_03000.txt', norms)
coord, Local_Heat = LocalHeat('/home/abetharan/HYDRO_KINETIC_COUPLING/Non_Linear_Ramp_Investigation/Pre_Heat_Ramp/', '03000', norms)
coord_sol = np.loadtxt('/media/abetharan/DATADRIVE2/Abetharan/pre_heat_model_test_problem/OUTPUT/GRIDS/X_GRID.txt')

unnorm_local_coord = coord/norms['lambda_mfp']
unnorm_local_qe = -Local_Heat/norms['qe']
unnorm_non_local_qe = non_local_qe/norms['qe']

#pre_heat_start, pre_heat_end, coord_pre_heat, preheat, interpolated_local_heat = NonLocalPreHeatModel(coord,coord_sol*norms['lambda_mfp'], -Local_Heat, non_local_qe)
pre_heat_start, pre_heat_end, coord_pre_heat, preheat, interpolated_local_heat = Alt2NonLocalPreHeatModel(unnorm_local_coord, coord_sol, unnorm_local_qe, unnorm_non_local_qe, norms)
#front_heat_start, front_heat_end, coord_front_heat, frontheat = FrontHeat(coord, coord_sol * norms['lambda_mfp'], -Local_Heat, non_local_qe)

#plt.plot(coord_sol * norms['lambda_mfp'], non_local_qe, 'r--', label = 'non-local')

plt.plot(coord_sol[pre_heat_start:]*norms['lambda_mfp'], interpolated_local_heat[pre_heat_start:], 'k', label = 'local')
plt.plot(coord_pre_heat, preheat, label = 'preheat')
plt.plot(coord_sol[(pre_heat_start):] * norms['lambda_mfp'], non_local_qe[pre_heat_start:], label = 'non-local')
#plt.plot(coord_sol[:front_heat_start]*norms['lambda_mfp'], interpolated_local_heat[:front_heat_start], 'k', label = 'local')
#plt.plot(coord_front_heat, frontheat, label = 'preheat')
#plt.plot(coord_sol[:(front_heat_start + 1)] * norms['lambda_mfp'], non_local_qe[:(front_heat_start  + 1)], label = 'non-local')



#plt.plot(np.linspace(0,coord[-1], 603))
#plt.plot(coord_sol * norms['lambda_mfp'])
#plt.ylim(0, 500)
#plt.plot(coord, non_local_qe/Local_Heat)

plt.ylabel('qe[10E14W/m^2]')
plt.xlabel('x[m]')
plt.legend()
plt.show()