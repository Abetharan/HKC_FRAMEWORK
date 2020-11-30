import itertools
import h5py
import os 
import matplotlib.pyplot as plt 
import random as rand 
import math
import pandas as pd
from scipy.special import comb
from scipy import constants
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy import integrate
import numpy as np 
rand.seed(0)
np.random.seed(0)
#global constant
kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
epsilon_0 = 8.854188E-12    # Vacuum dielectric constant
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")
ev_converter = e/kb

def SpitzerHarm(x, Te, ne, Z):
    def lambda_ei(T_norm , n_norm, Z_norm, return_arg = False):
        coulomb_logs = []
        T_norm = T_norm

        for T,n,Z in zip(T_norm, n_norm, Z_norm):
            if T < 10.00 * Z ** 2:
                result = 23.00 - math.log(math.sqrt(n * 1.00E-6) * Z * (T) ** (-3.00/2.00))
            else:
                result = 24.00 - math.log(math.sqrt(n * 1.00E-6) / (T))   

            if return_arg:
                return result
            else:
                coulomb_logs.append(result)
        return coulomb_logs

    coulomb_log = np.array(lambda_ei(Te * (kb/e), ne , Z))
    kappaE = 1.843076667547614E-10 * pow(Te, 2.5) * pow(coulomb_log, -1) *  pow(Z, -1)
    nx = len(Te)
    HeatFlowE = np.zeros(nx + 1, dtype=np.float64)
    gradTe = np.zeros(nx + 1)
    for i in range(1, nx):
        centered_ke = 0.5 * (kappaE[i] + kappaE[i - 1])
        gradTe[i] = ((Te[i] - Te[i - 1]) / (x[i] - x[i - 1]))
        HeatFlowE[i] = centered_ke * gradTe[i] 
    
    HeatFlowE[0] = 0
    HeatFlowE[-1] = 0
    return(-1 * HeatFlowE)


def free_param_calc(norm_Te, norm_ne, norm_Z):
    def lambda_ei(n, T, T_norm = 10, n_norm = 0.25E20, Z_norm = 1.0):
        if T * T_norm < 10.00 * Z_norm ** 2:

            result = 23.00 - math.log(math.sqrt(n * n_norm * 1.00E-6) * Z_norm * (T * T_norm) ** (-3.00/2.00))

        else:

            result = 24.00 - math.log(math.sqrt(n * n_norm * 1.00E-6) / (T * T_norm))   

        return result
    def lambda_ee(n, T, T_norm, n_norm, Z_norm):
        return (23.5 - math.log(pow(n_norm, 0.5) * pow(T_norm, -5/4)) 
                - pow((1E-5 + pow((math.log(T_norm) - 2),2)/16),0.5))
    

    
    #9.10938E-31            # Electron mass in kg
    el_charge = e #1.602189E-19    # Electron charge
    # gamma_ee_0 = el_charge ** 4 / (4 * math.pi * (me * epsilon_0) ** 2)
    gamma_ee_0 = el_charge ** 4 / pow(4 * math.pi * (epsilon_0), 2)

    T_J = norm_Te * el_charge
    gamma_ei_0 = norm_Z ** 2 * gamma_ee_0
    v_t = math.sqrt(2.0 * T_J / me)                  # Velocity normalization
    time_norm = v_t ** 3 / (gamma_ei_0 * (norm_ne/norm_Z) * lambda_ei(1.0, 1.0, T_norm = norm_Te, n_norm = norm_ne, Z_norm = norm_Z)) # Time normalization
    # time_norm_ee = v_t ** 3 / (gamma_ee_0 * norm_ne * lambda_ee(1.0, 1.0, T_norm = norm_Te, n_norm = norm_ne, Z_norm = norm_Z)) # Time normalization
    # x_0 = v_t * time_norm               # Space normalization
    # x_ee = v_t * time_norm_ee
    # lambda_mfp = np.sqrt(x_0 * x_ee)
    phi = ((np.mean(norm_Z) + 4.2) / (np.mean(norm_Z + 0.24)))
    # lambda_mfp = pow(T_J, 2) / (4 * math.pi * gamma_ee_0 * pow(norm_Z, 0.25) * norm_ne * lambda_ei(1, 1, norm_Te, norm_ne, norm_Z))
    # lambda_ed = 2 * pow(2, 0.5) * lambda_mfp
    lambda_mfp = pow(T_J, 2) / (4 * math.pi * gamma_ee_0 * (norm_Z + 1) * norm_ne * lambda_ei(1, 1, norm_Te, norm_ne, norm_Z))
    lambda_ed = pow((norm_Z + 1), 0.5) * lambda_mfp
    return(lambda_ed)

def new_mfp(Te, ne, Z):
    def lambda_ei(n, T, T_norm = 10, n_norm = 0.25E20, Z_norm = 1.0):
        coulomb_logs = []
        for T,n,Z in zip(T_norm, n_norm, Z_norm):
            if T < 10.00 * Z ** 2:
                result = 23.00 - math.log(math.sqrt(n * 1.00E-6) * Z * (T) ** (-3.00/2.00))
            else:
                result = 24.00 - math.log(math.sqrt(n * 1.00E-6) / (T))   
            coulomb_logs.append(result)
        return(np.array(coulomb_logs))

    def lambda_ee(n, T, T_norm, n_norm, Z_norm):
        return (23.5 - np.log(pow(n_norm, 0.5) * pow(T_norm, -5/4)) 
                - pow((1E-5 + pow((np.log(T_norm) - 2),2)/16),0.5))
    
    gamma_ee_0 = e ** 4 / pow(4 * math.pi * (epsilon_0), 2)
    T_J = Te * e
    ln_ei = lambda_ei(1.0, 1.0, Te, ne, Z)
    ln_ee = lambda_ee(1.0, 1.0, Te, ne, Z)
    Z_star = (np.mean(pow(Z, 2)) * ln_ei / 
                (np.mean(Z) * ln_ee ))
    phi = (Z_star + 4.2) / (Z_star + 0.24)
    lambda_mfp= pow(T_J, 2) / (gamma_ee_0 * ne * pow((Z_star * phi * ln_ei * ln_ee), 0.5))
    return lambda_mfp 

def spatialNorms(Te, ne, Z):
    lambda_ed = np.zeros(len(Te))    
    for idx, (i,j,k) in enumerate(zip(Te, ne, Z)):
        lambda_ed[idx] = free_param_calc(i * (kb/e), j ,k)
    return lambda_ed
# def lambda_ed(Te, ne, Z):
def getWeights(x, ne, lambda_ed, a, z):
    """
    Purpose: Create the weights using 
            Luciani 1983 PRL definition 
    Args: 
        x = Grid 
        z = number of mfp to convolve
        ne = number density
        a = tunable paramter
        lambda_ed = Luciani mfp definition
    """
    counter = 0 
    lambda_mfp = a * lambda_ed 
    search_range = np.nanmin(lambda_mfp) * z
    weight = np.zeros((len(x), len(x))) #Kernel for each x position 
    for i, pos in enumerate(x):
        if i == len(x) - 1 or i == 0:
            continue
        upper_x = pos + search_range 
        lower_x = pos - search_range 
        idx = np.where((x <= upper_x) & (x >= lower_x))
        if lower_x < 0:
            print("k")
        if upper_x > x[-1]:
            print("k2")
        for index in idx[0]:
            if index > len(ne):
                continue 
            density = (ne[i - 1] + ne[index - 1]) / 2 #trapzodial integral rule
            exp_value = (abs(pos - x[index]) * density) / (lambda_mfp[index - 1] * ne[index - 1])
            #print(exp_value)
            if abs(exp_value) < 1:
                #print("below 0 for {} and index {}".format(pos, index))
                counter += 1
            weight[i, index] = (1/(2 * lambda_mfp[index - 1])) * np.exp(-1 * exp_value)
    #print(counter)
    return weight

def convolve(q_sh, weights):
    heat_flow = np.zeros(len(q_sh))
    for idx in range(1, len(q_sh)):
        heat_flow[idx] = np.dot(q_sh, weights[idx, :])
    #enforce no source condition        
    heat_flow[0] = 0
    heat_flow[-1] = 0
    return(heat_flow)
def deleteArrElement(arr,index):
    return  np.delete(arr, index)

def getLuciani(x, ne, lambda_ed, a, z, q_sh, lower_bc, upper_bc):
    """
    Purpose: Create the weights using 
            Luciani 1983 PRL definition 
    Args: 
        x = Grid 
        z = number of mfp to convolve
        ne = number density
        a = tunable paramter
        lambda_ed = Luciani mfp definition
    """
    counter = 0 
    lambda_mfp = a * lambda_ed 
    weight = np.zeros((len(x), len(x))) #Kernel for each x position 

    weights = []
    heat_flow = np.zeros(len(q_sh))


    for i, pos in enumerate(x):
        if i == len(x) - 1 or i == 0:
            continue
        search_range = lambda_mfp[i] * z
        lower_pad_index = 0
        upper_pad_index = 0
        bc_ne = np.copy(ne)
        bc_q_sh = np.copy(q_sh)
        bc_lambda_mfp = np.copy(lambda_mfp)
        bc_x = np.copy(x)
        upper_x = pos + search_range 
        lower_x = pos - search_range 
        idx = np.where((x <= upper_x) & (x >= lower_x))
        if lower_x < 0:
            lower_pad_index = np.where((x <= abs(lower_x)))[0][-1] #valid for reflection not wrap
            if lower_pad_index % 2 == 0:
                pass
            else:
                lower_pad_index +=1

            if lower_bc == 'reflect':
                lower_ne = ne[:lower_pad_index:-1] 
                lower_x_coord = x[:lower_pad_index:-1]
                lower_q_sh = q_sh[:lower_pad_index:-1]
                lower_lambda = lambda_mfp[:lower_pad_index:-1]

            elif lower_bc == 'periodic':
                lower_ne = ne[-1*lower_pad_index::] 
                lower_x_coord = x[-1*lower_pad_index::]
                lower_q_sh = q_sh[-1*lower_pad_index::]
                lower_lambda = lambda_mfp[-1*lower_pad_index::]

            if lower_bc == 'reflect' or lower_bc == 'periodic':
                bc_ne = np.delete(bc_ne, 0)
                bc_q_sh = np.delete(bc_q_sh, 0)
                bc_lambda_mfp = np.delete(bc_lambda_mfp, 0)
                # bc_x = np.delete(bc_x, 0)

                bc_ne = np.concatenate((lower_ne, bc_ne)) 
                bc_q_sh = np.concatenate((lower_q_sh, bc_q_sh)) 
                bc_lambda_mfp = np.concatenate((lower_lambda, bc_lambda_mfp)) 
                dx = np.diff(lower_x_coord)
                cumsum = np.cumsum(dx) * -1 
                bc_x = np.concatenate((cumsum[::-1], bc_x)) 

        if upper_x > x[-1]:
            diff = upper_x - np.max(x)
            upper_pad_index = np.where((x >= (np.max(x) - diff)))[0][0]#valid for periodic not wrap
            if upper_pad_index % 2 == 0:
                pass
            else:
                upper_pad_index -= 1

            if upper_bc == 'reflect':
                upper_ne = ne[:upper_pad_index:-1] 
                upper_x = x[:upper_pad_index:-1]
                upper_q_sh = q_sh[:upper_pad_index:-1]
                upper_lambda = lambda_mfp[:upper_pad_index:-1]
            elif upper_bc == 'periodic':
                upper_ne = ne[:upper_pad_index:1] 
                upper_x_coord = x[:upper_pad_index:1]
                upper_q_sh = q_sh[:upper_pad_index:1]
                upper_lambda = lambda_mfp[:upper_pad_index:1]

            if upper_bc == 'reflect' or upper_bc == 'periodic':
                bc_ne = np.delete(bc_ne, -1)
                bc_q_sh = np.delete(bc_q_sh, -1)
                bc_lambda_mfp = np.delete(bc_lambda_mfp, -1)
                bc_x = np.delete(bc_x, -1)

                bc_ne = np.concatenate((bc_ne, upper_ne)) 
                bc_q_sh = np.concatenate((bc_q_sh, upper_q_sh)) 
                bc_lambda_mfp = np.concatenate((bc_lambda_mfp, upper_lambda)) 
                dx = np.diff(upper_x_coord)
                cumsum = np.cumsum(dx) + bc_x[-1]
                bc_x = np.concatenate((bc_x, cumsum)) 
        
        # plt.plot(bc_x, bc_q_sh, 'rx')
        # plt.show()

        idx = np.where((bc_x <= upper_x) & (bc_x >= lower_x))
        weights = np.zeros(len(bc_ne))
        initial_position_index = np.where(bc_x == pos)[0][0]
        for index in idx[0]:
            dx = np.diff(bc_x)
            density = (bc_ne[i] + bc_ne[index]) / 2 #trapzodial integral rule
            exp_value = (abs(pos - bc_x[index]) * density) / (bc_lambda_mfp[index] * bc_ne[index])
            # if index < initial_position_index:
            #     density = integrate.simps(bc_ne[index:initial_position_index], bc_x[index:initial_position_index]) #simpsons integral rule
            # elif index == initial_position_index:
            #     density = 0 #simpsons integral rule
            # else:
            #     density = integrate.simps(bc_ne[i:index], bc_x[i:index]) #simpsons integral rule
            
            # exp_value = (density) / (bc_lambda_mfp[index] * bc_ne[index])
                
            #print(exp_value)
            if abs(exp_value) < 1:
                #print("below 0 for {} and index {}".format(pos, index))
                counter += 1
            # weight[i, index] = (1/(2 * bc_lambda_mfp[index])) * np.exp(-1 * exp_value)
            weights[index] = (((1/(2 * bc_lambda_mfp[index])) * np.exp(-1 * exp_value)))
        
        heat_flow[i] = np.convolve(bc_q_sh, weights, "valid")
    #print(counter)
    return heat_flow

# def fit_func(x, cell_wall_x, cell_centre_x, cell_wall_Te, Te, cell_wall_ne, ne, cell_wall_Z, Z, true_qe):
#     local_heat = SpitzerHarm(cell_centre_x, Te * (e/kb), ne, Z)
    # lambda_ed = spatialNorms(cell_wall_Te * (e/kb), cell_wall_ne, cell_wall_Z)
    # # lambda_ed = new_mfp(cell_wall_Te, cell_wall_ne, cell_wall_Z)
    # weights = getWeights(cell_wall_x, cell_wall_ne, lambda_ed, x[0], x[1])
    # q = convolve(local_heat, weights)
    # return(true_qe[1:-1] - q[1:-1])

def fit_func(x, full_x_grid, interp_Te, interp_ne, interp_Z, true_qe):
    local_heat = SpitzerHarm(full_x_grid[1::2], interp_Te[1::2] * (e/kb), interp_ne[1::2], interp_Z[1::2])
    lambda_ed = spatialNorms(interp_Te * (e/kb), interp_ne, interp_Z)
    interp_local_heat = np.interp(full_x_grid, full_x_grid[::2], local_heat)
    luciani = getLuciani(full_x_grid, interp_ne, lambda_ed, x[0], x[1], interp_local_heat, 'periodic', 'periodic')
    return(true_qe - luciani[::2])

def OptimizedConvolveHeat(x, cell_wall_x, cell_centre_x, cell_wall_Te, Te, cell_wall_ne, ne, cell_wall_Z, Z):
    # local_heat = SpitzerHarm(cell_centre_x, Te * (e/kb), ne, Z)
    # lambda_ed = spatialNorms(cell_wall_Te * (e/kb), cell_wall_ne, cell_wall_Z)
    # # lambda_ed = new_mfp(cell_wall_Te, cell_wall_ne, cell_wall_Z)
    # weights = getWeights(cell_wall_x, cell_wall_ne, lambda_ed, x[0], x[1])
    # q = convolve(local_heat, weights)
    # return(q)
    local_heat = SpitzerHarm(full_x_grid[1::2], interp_Te[1::2] * (e/kb), interp_ne[1::2], interp_Z[1::2])
    lambda_ed = spatialNorms(interp_Te * (e/kb), interp_ne, interp_Z)
    interp_local_heat = np.interp(full_x_grid, full_x_grid[::2], local_heat)
    luciani = getLuciani(full_x_grid, interp_ne, lambda_ed, x[0], x[1], interp_local_heat, 'periodic', 'periodic')
    return(luciani)

def epperlein_short(nx, L, Z_ = 37.25, ne_ = 1e27, Te_ = 100., perturb = 1e-3, sin = False):

    x_l = 0
    x_u = L
    initial_coord = np.linspace(x_l, x_u, nx+1, dtype=np.float64)
    x_centered = np.array([(initial_coord[i] + initial_coord[i+1]) /2 for i in range(len(initial_coord) -1)], dtype=np.float64)
    Z = np.zeros(nx, dtype = np.float64) + Z_#37.25
    
    ne = np.zeros(nx, dtype = np.float64) + ne_# 1e27
    #temperatureE = np.linspace(3, 300, nx)
    
    ##Relevant for periodic systems .. fluid cannot be
    if sin: 
        Te = np.zeros(nx, dtype=np.float64) + Te_  + perturb * Te_  * np.sin((2*np.pi * x_centered) / np.max(x_centered), dtype=np.float64)
    else:
    ##Reflective
        Te = np.zeros(nx, dtype=np.float64) + Te_  + perturb * Te_  * np.cos((np.pi * x_centered) / np.max(x_centered), dtype=np.float64)


    return(initial_coord, x_centered, Te, ne, Z)

base_path = '/Users/shiki/DATA/HKC_RELATED//FLUID_INPUT/'
# Te = np.loadtxt(os.path.join(base_path, "electron_temperature.txt"), dtype = np.float64)
# rho = np.loadtxt(os.path.join(base_path, "density.txt"), dtype = np.float64)
# Z = np.loadtxt(os.path.join(base_path, "Z.txt"), dtype = np.float64)
# Ar = np.loadtxt(os.path.join(base_path, "Ar.txt"), dtype = np.float64)
# coord = np.loadtxt(os.path.join(base_path, "coord.txt"), dtype = np.float64)
# true_qe = -1 * np.loadtxt(os.path.join('/Users/shiki/DATA/HKC_RELATED//FLUID_OUTPUT/', 'ELECTRON_HEAT_FLOW_X/ELECTRON_HEAT_FLOW_X_0.txt'), dtype = np.float64)
# centered_x = np.array([(coord[i] + coord[i+1]) / 2 for i in range(len(coord) - 1)], dtype = np.float64)
# ne = Z * (rho / (Ar * mp))
coord, centered_x, Te, ne, Z = epperlein_short(100, 769.7652593/2, Z_= 1, ne_ = 1e19)
lambda_ed = spatialNorms(Te, ne, Z)
cell_wall_Te = np.zeros(len(coord) - 2, dtype = np.float64)
cell_wall_ne = np.zeros(len(coord) - 2, dtype = np.float64)
cell_wall_Z = np.zeros(len(coord) - 2, dtype = np.float64)
length_of_sol_kit_grid = len(coord) + len(centered_x)
full_x_grid = np.zeros(length_of_sol_kit_grid)
j = 0
k = 0
for i in range(length_of_sol_kit_grid):
    if i % 2 != 0:
        full_x_grid[i] = centered_x[j]
        j +=1
    else:
        full_x_grid[i] = coord[k]
        k +=1 

# for i in range(len(ne) - 1):
#     cell_wall_ne[i] = (ne[i] + ne[i + 1]) /2 
    # cell_wall_Te[i] = (Te[i] + Te[i + 1]) /2 
    # cell_wall_Z[i] = (Z[i] + Z[i + 1]) /2 

interp_ne = np.interp(full_x_grid, centered_x, ne)
interp_Te = np.interp(full_x_grid, centered_x, Te)
interp_Z = np.interp(full_x_grid, centered_x, Z)
# interp_ne = np.interp(full_x_grid, centered_x, ne)
# cell_wall_Te = np.append(cell_wall_Te, [0,0])
# cell_wall_ne = np.append(cell_wall_ne, [0,0])
# cell_wall_Z = np.append(cell_wall_Z, [0,0])
local_heat = SpitzerHarm(centered_x, Te * (e/kb), ne, Z)
interp_local_heat = np.interp(full_x_grid, coord, local_heat)
lambda_ed = spatialNorms(interp_Te * (e/kb), interp_ne, interp_Z)
# lambda_ed = new_mfp(cell_wall_Te, cell_wall_ne, cell_wall_Z)
# weights = getWeights(coord, cell_wall_ne, lambda_ed, 32, 5)
# convolve_heat = convolve(local_heat, weights)
luciani = getLuciani(full_x_grid, interp_ne, lambda_ed, 0.25,  1, interp_local_heat, 'periodic', 'periodic')
x0 = np.array([1, 1])
res_lsq = least_squares(fit_func, x0, args =(full_x_grid, interp_Te,
                                            interp_ne, interp_Z,  local_heat))
print(res_lsq.x)
y = OptimizedConvolveHeat(res_lsq.x, coord, centered_x,cell_wall_Te,Te, cell_wall_ne, ne, cell_wall_Z, Z)
plt.plot(coord, y[::2], label = 'optimized')
# plt.plot(coord - coord[75], weights[10,:])
plt.plot(coord, local_heat, 'r--', label = 'Spitzer-Harm')
plt.plot(full_x_grid[::2], luciani[::2], 'y--', label = 'Luciani 1983')
# plt.plot(coord[1:-1], (convolve_heat/local_heat)[1:-1])
plt.legend()
plt.show()