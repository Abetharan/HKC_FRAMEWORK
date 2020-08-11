import numpy as np 
from scipy import linalg
import matplotlib.pyplot as plt 
from scipy import constants
import math
kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
epsilon_0 = 8.854188E-12    # Vacuum dielectric constant
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")

def _interpolate(x, x_centered, Te, ne, Z, boundary_condition, normalised_values):
    # NOrmalise SI to Impact norms
    sol_kit_x_grid = x# / normalised_values["lambda_mfp"]
    sol_kit_x_centered_grid = x_centered# / normalised_values["lambda_mfp"]
    length_of_sol_kit_grid = len(sol_kit_x_centered_grid) + len(sol_kit_x_grid)
    full_x_grid = np.zeros(length_of_sol_kit_grid)
    j = 0
    k = 0
    for i in range(length_of_sol_kit_grid):
        if i % 2 != 0:
            full_x_grid[i] = sol_kit_x_centered_grid[j]
            j +=1
        else:
            full_x_grid[i] = sol_kit_x_grid[k]
            k +=1 

    
    #SOL-KiT periodic condition == .\.\.\ i.e. start with cell centre and end with cell wall
    #SOL-KiT noflow == .\.\. i.e. start and end with cell centres
    sol_kit_x_grid = None
    # if boundary_condition == 'periodic':
    #     sol_kit_grid = full_x_grid[:-1]
    
    # elif boundary_condition == "noflow":
    #     sol_kit_grid = full_x_grid[1:-1]

    sol_kit_grid = full_x_grid
    #Conver to SOL-KiT units
    sol_kit_ne = ne #/ (normalised_values['ne'])
    sol_kit_te = Te #* (kb/e)) / normalised_values['Te']
    #SOL_KIT_laser = (f_laser * f_density) / power_density
    #SOL_KIT_brem = (f_brem * f_density) / power_density
    sol_kit_z = Z
    # SOL_KIT_HEATING = SOL_KIT_laser + SOL_KIT_brem
    
    #Require interpolation to get the centre quanties in SOL-KiT this is done via linear interpolations 
    #here we use cubic spline to smooth quanties. 
    sol_kit_inter_ne = np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_ne)
    sol_kit_inter_te = np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_te)
    sol_kit_inter_z =  np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_z)
    #SOL_KIT_inter_ne = spliner_ne(SOL_KIT_grid)
    #SOL_KIT_inter_Te = spliner_Te(SOL_KIT_grid)
    #SOL_KIT_inter_Z = spliner_Z(SOL_KIT_grid)
    return(sol_kit_grid, sol_kit_inter_te, sol_kit_inter_ne, sol_kit_inter_z)
#    return(np.array([(x[i] + x[i+1])/ 2 for i in range(len(x) - 1)]))

def lambda_ei(T_norm = 10, n_norm = 0.25E20, Z_norm = 1.0, return_arg = False):
    coulomb_logs = []
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

def free_param_calc(Te, ne, Z):
    #9.10938E-31            # Electron mass in kg
    el_charge = e #1.602189E-19    # Electron charge
    r_b = bohr_radi #5.29E-11              # Bohr radius in metres

    sigma_0 = math.pi * r_b ** 2  # Cross-section normaliation
    gamma_ee_0 = el_charge ** 4 / (4 * math.pi * (me * epsilon_0) ** 2)

    T_J = Te * el_charge

    gamma_ei_0 = Z ** 2 * gamma_ee_0

    v_t = np.sqrt(2.0 * T_J / me)                  # Velocity normalization
    k = lambda_ei(T_norm = Te, n_norm = ne, Z_norm = Z)
    time_norm = pow(v_t, 3) / (gamma_ei_0 * (ne/Z) * lambda_ei(T_norm = Te, n_norm = ne, Z_norm = Z)) # Time normalization
    x_0 = v_t * time_norm               # Space normalization

    e_0 = me * v_t / (el_charge * time_norm) # Electric field normalization
    q_0 = me * ne * pow(v_t, 3)

    param_dict = {}

    param_dict['ne'] = ne
    param_dict['Te'] = Te
    param_dict['Z'] = Z
    ni = ne / Z
    param_dict['ni'] = ni
    param_dict['vte'] = v_t
    param_dict['tau_ei'] = time_norm
    param_dict['nu_ei'] = 1/time_norm
    param_dict['lambda_mfp'] = x_0
    param_dict['qe'] = q_0
    return(param_dict)

def lambda_star(nue_ei, v_grid, Z):
    lambda_star = np.zeros((len(nue_ei), len(v_grid)))
    for i, local_nue in enumerate(nue_ei):
        zeta = (Z[i] + 0.24) / (Z[i] + 4.2)
        for j, v in enumerate(v_grid):
            lambda_star[i, j] = zeta * (v/local_nue)
    return(lambda_star)

# def lambda_star(ne, Te, Z, coulomb_log, v_grid):
#     lambda_star = np.zeros((len(Te), len(v_grid)))
#     gamma_ee = pow(e, 4) / pow((4 * np.pi * epsilon_0), 2)
#     for i, (n, T, ioniza, lamb_ei) in enumerate(zip(ne, Te, Z, coulomb_log)):
#         zeta = (ioniza + 0.24) / (ioniza + 4.2)
#         for j, v in enumerate(v_grid):
#             lambda_star[i, j] = ((2 * pow(2, 0.5)*(0.5 * me * pow(v, 2))) / 
#                                 (4 * np.pi * pow(ioniza, 0.5) * gamma_ee * lamb_ei * pow(zeta, 0.5)))

    # return(lambda_star)

def electric_field(x, Te, ne, Z):
    n = len(Te)
    E = np.zeros(n)
    for i in range(n - 2):
        dx = x[i + 1] - x[i]
        gamma =  1 + 3 * (Z[i] + 0.477) / (2*(Z[i] + 2.15))
        if i == n - 1:
            print("kapp")
        E[i + 1] = -1 * (kb/(pow((4 * np.pi * epsilon_0),0.5) * e)) * Te[i] * ((ne[i + 1] - ne[i]) / (ne[i] * dx) +
                        gamma * ((Te[i + 1] - Te[i]) /(Te[i] * dx)))
    return(E)

def lambda_E(lambda_star, E, v_grid):
    lambda_E = np.zeros(np.shape(lambda_star))
    for i, e_field in enumerate(E):
        for j, v in enumerate(v_grid):
            k = 1/(lambda_star[i, j]) + abs((e * e_field) 
                                        / (0.5 * me * pow(v, 2))) 
            lambda_E[i, j] = 1/k

    return(lambda_E)

def f_maxwellian(Te, ne, v_grid, v_grid_width):
    """
    Purpose: Init maxwellian f_0 given some parameters
    Args:
        Te = Electron Temperature 
        ne = Number Density
        v_grid = velocity grid being used, taken from previous sol-kit cycle
        v_grid_width = width of velocity cells.
    """
    nv = len(v_grid)
    nx = len(Te)
    f0 = np.zeros((nx, nv))
    f0 = np.zeros((nx, nv))
    n_num = np.zeros(nx)
    #Method from f_init.f90:215::227
    for i in range(nx):
        for v in range(nv):
            v_th = np.sqrt((kb * Te[i])/me)
            f0[i,v] = ne[i] * ((np.pi * 2*v_th) ** (- 3.00/2.00)) * np.exp( - (v_grid[v] ** 2)/ (2*v_th**2))
            n_num[i] = n_num[i] + 4.00 * np.pi * v_grid[v] ** 2 * v_grid_width[v] * f0[i,v]
        for v in range(nv):
            f0[i, v] = f0[i,v] * ne[i]/n_num[i]
    return(f0)

def g_1_maxwellian(lambda_star, f_0_maxwellian, x, Te):
    nx, nv = np.shape(lambda_star)
    g_1 = np.zeros((nx ,nv))
    for i in range(nx - 2):
        for v in range(nv):
            # if v == 0:
            #     g_1[i + 1, v] = 0# -1 * lambda_star[i + 1, v] * f_0_maxwellian[i + 1, v] * (gradient/Te[i + 1]) 
                # continue
            dx = x[i + 1] - x[i]
            gradient = (Te[i + 1] - Te[i]/ dx )
            # centered = (Te[i] + Te[i + 1]) / 2
            g_1[i + 1, v] = -1 * lambda_star[i + 1, v] * f_0_maxwellian[i + 1, v] * (gradient/Te[i]) 
    return(g_1)

def get_delta_f(g_1, nue_ei, Z, lambda_E, v_grid, x):
    nx, nv = np.shape(g_1)
    delta_f0 = np.zeros((nx - 1, nv))
    r = 2
    A = np.zeros(nx)
    B = np.zeros(nx)
    C = np.zeros(nx)
    S = np.zeros(nx - 1)
    for j, v in enumerate(v_grid):
        for i in range(0, nx - 1):
            dx_k = x[i+1] - x[i]
            dx_k_1_2 = x[i+3] - x[i + 1]
            A[i] = lambda_E[i, j] / (3 * dx_k * dx_k_1_2)
            B[i] = -1*((lambda_E[i + 1, j] / (3 * dx_k*dx_k_1_2) +
                lambda_E[i, j] / (3 * dx_k*dx_k_1_2) + 
                ((r * nue_ei[i]) / (v * Z[i]))))
            C[i] = lambda_E[i + 1, j] / (3 * dx_k*dx_k_1_2)
            S[i] = (g_1[i, j] - g_1[i - 1, j]) / (3 * dx_k_1_2)
        mat = createMatrix(A, B, C)
        d = solve(mat, S)
        delta_f0[:, j] = d

    return delta_f0

def integrate(delta_f0, lambda_E, v_grid, v_grid_width, x_grid):

    nx, nv = np.shape(delta_f0)
    dv = np.diff(v_grid)
    dx = np.diff(x_grid)
    dq = np.zeros(nx + 1)
    for i in range(1, nx - 1):
        heat_flow = 0
        for v in range(nv):
            heat_flow +=  pow(v_grid[v], 5) * v_grid_width[v] * lambda_E[i, v] * ((delta_f0[i+1, v] - delta_f0[i, v]) / dx[i])
        dq[i] = ((2 * np.pi * me)/3) * heat_flow
    return dq
def snb_heat_flow(x, v_grid, Te, ne , Z, norms = 0):
    """
    Purpose: Calculates SNB-Heat flow
    Args: 
        x = cell-wall coords nx+1
        v_grid = cell-wall velocity grid nv+1
        Te = cell-centered temperature nx 
        ne = cell-centered density nx
        Z = cell-centered ionisation nx 
        Ar = cell-centered atomix mass nx
    Returns: q_snb
    """

    dv = np.diff(v_grid)
    x_centered = np.array([(x[i+1] + x[i])/ 2 for i in range(len(x) - 1)])
    v_grid_centered = np.array([(v_grid[i+1] + v_grid[i])/ 2 for i in range(len(v_grid) - 1)])
    full_x, interp_Te, interp_ne, interp_Z = _interpolate(x, x_centered, Te, ne, Z, 'noflow', norms)

    free_params= free_param_calc(interp_Te * (kb/e), interp_ne ,interp_Z )
    lamb_star = lambda_star(free_params['nu_ei'], v_grid_centered * free_params['vte'][1], interp_Z)
    # coulomb_logs = lambda_ei(interp_Te,interp_ne, interp_Z)
    # lamb_star = lambda_star(interp_ne, interp_Te, interp_Z,
    #                  coulomb_logs, v_grid_centered * free_params['vte'][1])

    #Electric field only needs to be defined at cell- walls 
    #Every second entry in computer indicies 0 2 etc are cell-walls 
    E_field = electric_field(centered_x, interp_Te[1::2], interp_ne[1::2], interp_Z[1::2])
    #Likewise Lambda_E only needs to be defined at cell-walls 
    lamb_E = lambda_E(lamb_star[::2], E_field, v_grid_centered * free_params['vte'][1])
    lamb_E = lamb_star[::2]
    #f0 should be definedo n the entire grid 
    f0_mb = f_maxwellian(interp_Te, interp_ne, v_grid_centered * free_params['vte'][1], dv * free_params['vte'][1])
    #g_1 only on cell-walls 
    g_1_mb = g_1_maxwellian(lamb_star[::2], f0_mb, x, Te)
     #Lambda_E, g_1 needs to be cell-wall and rest cell-centereed
    delta_f0 = get_delta_f(g_1_mb, free_params['nu_ei'][1::2], Z, lamb_E, v_grid_centered * free_params['vte'][1], full_x) 

    #integrate to get corrections
    dq = integrate(delta_f0, lamb_E, v_grid_centered * free_params['vte'][1], dv * free_params['vte'][1], x_centered)
    q_sh = spitzer_harm_heat(x, Te, ne, Z)
    #q_snb = q_sh - dq
    q_snb = q_sh - dq

    return q_sh, q_snb

def spitzer_harm_heat(x, Te, ne, Z):
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

def createMatrix(A, B, C):
    matrix = np.zeros((len(A) - 1, len(A) - 1))
    for i, (a,b,c) in enumerate(zip(A,B,C)):
        if i == 0:
            matrix[i, 0] = b 
            matrix[i, 1] = c
        elif i == len(A) - 2:
            matrix[-1, -2] = a
            matrix[-1, -1] = b
        elif i < len(A) - 2:
            matrix[i, i - 1] = a 
            matrix[i, i] = b 
            matrix[i, i +1] = c

    return matrix
def solve(A, b):
    """
    Purpose : Solves A*x = b 
    Args:
        A = Matrix (nv) 
        b = source term
    Returns: x 
    """ 
    x = linalg.solve(A, b)
    return(x)

def epperlein_short(nx, L, Z_ = 37.25, ne_ = 1e27, Te_ = 100., perturb = 1e-3, sin = False):

    x_l = 0
    x_u = L
    initial_coord = np.linspace(x_l, x_u, nx+1, dtype=np.float64)
    x_centered = np.array([(initial_coord[i] + initial_coord[i+1]) /2 for i in range(len(initial_coord) -1)], dtype=np.float64)
    Z = np.zeros(nx, dtype = np.float64) + Z_#37.25
    Ar = np.zeros(nx, dtype = np.float64) + 2#37.25
    
    ne = np.zeros(nx, dtype = np.float64) + ne_# 1e27
    #temperatureE = np.linspace(3, 300, nx)
    
    ##Relevant for periodic systems .. fluid cannot be
    if sin: 
        Te = np.zeros(nx, dtype=np.float64) + Te_  + perturb * Te_  * np.sin((2*np.pi * x_centered) / np.max(x_centered), dtype=np.float64)
    else:
    ##Reflective
        Te = np.zeros(nx, dtype=np.float64) + Te_  + perturb * Te_  * np.cos((np.pi * x_centered) / np.max(x_centered), dtype=np.float64)


    return(initial_coord, x_centered, Te, ne, Z, Ar)

nx = 100
nv = 100
coord, centered_x, Te, ne, Z, Ar = epperlein_short(nx, 784.48, Z_= 1, ne_ = 1e19, sin =True)
v_grid = np.linspace(0, 30, nv + 1)

q_sh, q_snb = snb_heat_flow(coord, v_grid, Te *(e/kb), ne, Z, Ar)
plt.plot(coord, q_sh, 'k--')
plt.plot(coord, q_snb)
plt.show()