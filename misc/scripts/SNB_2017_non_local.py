import numpy as np 
from scipy import linalg
import matplotlib.pyplot as plt 
from scipy import constants
import math
import sys
from mpl_toolkits.mplot3d import Axes3D
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

    v_t = np.sqrt(T_J / me)                  # Velocity normalization
    k = lambda_ei(T_norm = Te, n_norm = ne, Z_norm = Z) #0.75 * pow(np.pi/2, 0.5) *
    time_norm =   pow(v_t, 3) / (gamma_ei_0 * (ne/Z) * lambda_ei(T_norm = Te, n_norm = ne, Z_norm = Z)) # Time normalization
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
    param_dict['E'] = e_0
    return(param_dict)

def getLambdaStar(Te, ne, Z, v_grid):
    """
    Purpose: Calculates velocity dependent collision frequency
             and mfp
    Args:
        Te = Temperature
        ne = Numbe r density 
        Z = ionization
        v_grid = centered grid
    Returns: nue_ei, lambda_star
    """
    coulomb_log = lambda_ei(T_norm = Te * (kb/e), n_norm = ne, Z_norm = Z)
    gamma_ee_0 = e ** 4 / pow((4 * math.pi * me * epsilon_0), 2)
    nue_ei = np.zeros((len(Te), len(v_grid)))
    lambda_star = np.zeros((len(Te), len(v_grid)))
    for j,v in enumerate(v_grid):
        nue_ei[:, j] = (4 * np.pi * Z *ne*gamma_ee_0 *coulomb_log) / (pow(v,3)) 
        zeta = (Z + 0.24) / (Z + 4.2)
        lambda_star[:, j] = (v / nue_ei[:, j]) * zeta 
    return nue_ei, lambda_star
    
def lambda_star(nue_ei, v_grid, Z):
    lambda_star = np.zeros((len(nue_ei), len(v_grid)))
    for i, local_nue in enumerate(nue_ei):
        zeta = (Z[i] + 0.24) / (Z[i] + 4.2)
        for j, v in enumerate(v_grid):
            lambda_star[i, j] = zeta * (v/local_nue)
    return(lambda_star)

def electric_field(x, Te, ne, Z):
    """
    Purpose: Calculate Spizter-Harm Electric Field
    Args: 
        x = Cell-centered grid, shape (nx)
        Te = Full grid, shape (2 * nx + 1)
        ne = Full grid, shape (2 * nx + 1)
        Z = Cell-wall grid, shape (nx)
    Returns E_field with 0 field at boundaries. 
    Checked: 13/08/2020
    """
    n = len(x)
    E = np.zeros(n + 1)
    #gamma_ee_0 = e  / pow((4 * math.pi* epsilon_0), 0.5)
    for i in range(1, n):
        dx = x[i] - x[i - 1]
        gamma =  1 + ((3 * (Z[i] + 0.477)) / (2*(Z[i] + 2.15)))
        E[i] = (-1 * (kb/e) * Te[i*2] * 
                    (((ne[i*2 + 1] - ne[i*2 - 1]) / (ne[i*2] * dx)) +
                    gamma * ((Te[i*2 + 1] - Te[i*2 - 1]) /(Te[i*2] * dx))))
    return(E)

def lambda_E(lambda_star, E, v_grid):
    """
    Purpose: Phemenlogical corrections to mfp via addition of E-field
    Args:
        Lambda_star = velocity dependent mfp 
        E = E-field
        v_grid = centered v grid 
    Returns:
        lambda_E = corrected mfp.
    """
    lambda_E = np.zeros(np.shape(lambda_star))
    for i, e_field in enumerate(E):
        for j, v in enumerate(v_grid):
            k = 1/(lambda_star[i, j]) + abs((e * e_field) / (0.5 * me * pow(v, 2))) 
            lambda_E[i, j] = 1/k

    return(lambda_E)

def f_maxwellian(Te, ne, v_grid, v_grid_width):
    """
    Purpose: Init maxwellian f_0 given some parameters
    Args:
        Te = Electron Temperature 
        ne = Number Density
        v_grid = velocity grid being used
        v_grid_width = width of velocity cells.
    """
    nv = len(v_grid)
    nx = len(Te)
    f0 = np.zeros((nx, nv))
    n_num = np.zeros(nx)
    #Method from f_init.f90:215::227
    for i in range(nx):
        v_th = np.sqrt((kb * Te[i])/me)
        for v in range(nv):
            f0[i,v] = ne[i] * ((np.pi * 2*v_th) ** (- 3.00/2.00)) * np.exp( - (v_grid[v] ** 2)/ (2*v_th**2))
            n_num[i] = n_num[i] + 4.00 * np.pi * v_grid[v] ** 2 * v_grid_width[v] * f0[i,v]
        
        #Ensure density is consistent via scaling with numerical density.
        for v in range(nv):
            f0[i, v] = f0[i,v] * ne[i]/n_num[i]
    return(f0)

def g_1_maxwellian(lambda_star, f_0_maxwellian, x, Te):
    """
    Purpose: Calculates g_1
    Args: 
        lambda_star = cell-wall mfp defined in energy groups with collision fix, shape = (nx + 1, nv).
        f_0_maxwellian = maxwell-boltzmann distribution defined on entire grid, shape = (2*nx +1, nv)
        x = cell-centered grid, shape = (nx)
        Te = full grid Te, shape = (2*nx + 1)
    Returns:
        g_1
    Checked: 13/08/2020
    """
    nx, nv = np.shape(lambda_star)
    g_1 = np.zeros((nx ,nv))
    for i in range(1, nx - 1):
        dx = x[i] - x[i - 1]
        gradient = (Te[i*2 + 1] - Te[i*2 - 1])/ dx
        temperature_correction = (gradient/Te[i*2])
        for v in range(nv):
            g_1[i, v] = -1 * lambda_star[i, v] * f_0_maxwellian[2*i, v] * temperature_correction
    return(g_1)

def get_delta_f(g_1, nue_ei, Z, lambda_E, v_grid, x):
    """
    Purpose: Calculates delta_f_0 
    Args:
        g_1: modified SNB f_1 maxwellian, shape (nx + 1, nv)
        nue_ei: Collision Frequency, shape nx
        Z: Ionization, shape nx 
        Lambda_E: SNB modified collision mfp shape (nx + 1, nv) 
        v_grid: Velocity grid/groups, shape (nv)
        x: full X grid has length,  2*nx + 1
    Returns:
        delta_f_0
    Checked: 13/08/2020
    """
    nx, nv = np.shape(g_1)
    delta_f0 = np.zeros((nx - 1, nv))
    r = 2
    for j, v in enumerate(v_grid):
        A = np.zeros(nx - 1)
        B = np.zeros(nx - 1)
        C = np.zeros(nx - 1)
        S = np.zeros(nx - 1)
        for i in range(0, nx - 1):
            if i == 0:
                dx_k = 2 * (x[1] - x[0]) #extrapolate grid unformly
            else:
                dx_k = x[2 * i + 1] - x[2 * i - 1]

            if 2 *i + 3 > len(x) - 1:
                dx_k_1 = 2 * (x[-1] - x[-2]) #extrapolate grid unformly beyond grid 
            else:
                dx_k_1 = x[2 * i + 3] - x[2 * i + 1]

            if i == 0:
                x__1_2 = 0 - (x[1] - x[0])
                dx_k_1_2 = x[2 * i + 3] - x__1_2
            if 2 *i + 3 > len(x) - 1:
                x_k_n_1_2 = x[-1] + (x[-1] - x[-2]) 
                dx_k_1_2 = x_k_n_1_2 - x[2 * i - 1] 
            else:
                dx_k_1_2 = x[2 * i + 3] - x[2 * i - 1]

            A[i] = lambda_E[i, j] / (3 * dx_k * dx_k_1_2)

            B[i] = -1*((lambda_E[i + 1, j] / (3 * dx_k_1*dx_k_1_2) +
                lambda_E[i, j] / (3 * dx_k*dx_k_1_2) + 
                ((r * nue_ei[i,j]) / (v * Z[i]))))
                
            C[i] = lambda_E[i + 1, j] / (3 * dx_k_1*dx_k_1_2)
            
            S[i] = (g_1[i + 1, j] - g_1[i, j]) / (3 * (x[2 * i + 2] - x[2 * i]))
        
        mat = createMatrix(A, B, C)
        d = solve(mat, S)
        delta_f0[:, j] = d

    return delta_f0

def integrate(delta_f0, lambda_E, v_grid, v_grid_width, x_grid):
    """
    Purpose: numerically integrate to get non-local correction
    Args:
        delta_f0 : Pertubed f0 f_mb - non_local_f_0, shape (nx, nv)
        lambda_E : phemenlogically corrected mfp, shape (nx + 1, nv) 
        v_grid : cell-centered velcoity grid, shape (nv) 
        v_grid_width: dv , shape (nv  - 1)
        x_grid : cell centered x grid, shape (nx)
    Returns: 
        dq = non-local correction
    Checked: 13/08/2020
    """
    nx, nv = np.shape(delta_f0)
    dx = np.diff(x_grid)
    dq = np.zeros(nx + 1)
    for i in range(1, nx):
        heat_flow = 0
        grad_delta_f0 = ((delta_f0[i, :] - delta_f0[i - 1, :]) / dx[i - 1]) #Spatial derivative
        #Integrate
        for v in range(nv):
            heat_flow +=  pow(v_grid[v], 5) * v_grid_width[v] * lambda_E[i, v] * grad_delta_f0[v]
        print(heat_flow)
        dq[i] = -1*((2 * np.pi * me)/3) * heat_flow
    return dq
def integrate_g1(g_1, v_grid, v_grid_width):
    nx, nv = np.shape(g_1)
    q = np.zeros(nx)
    for i in range(1, nx):
        heat_flow = 0
        for v in range(nv):
            heat_flow +=  pow(v_grid[v], 5) * v_grid_width[v] * g_1[i, v]
        q[i] = ((2 * np.pi * me)/3) * heat_flow
    return q 

def snb_heat_flow(x, v_grid, Te, ne , Z, norms = 0):
    """
    Purpose: Calculates SNB-Heat flow
    Args: 
        x = cell-wall coords nx+1
        v_grid = cell-wall velocity grid defined in terms of v_th,  nv+1
        Te = cell-centered temperature nx 
        ne = cell-centered density nx
        Z = cell-centered ionisation nx 
        Ar = cell-centered atomix mass nx
    Returns: q_snb
    """
    x_centered = np.array([(x[i+1] + x[i])/ 2 for i in range(len(x) - 1)])
    full_x, interp_Te, interp_ne, interp_Z = _interpolate(x, x_centered, Te, ne, Z, 'noflow', norms)

    free_params= free_param_calc(interp_Te * (kb/e), interp_ne ,interp_Z )
    # print("Vth scaling factor is {}".format(free_params['vte'][1]))
    # vth = np.sqrt(100 * e /me)
    # print("Hand Calculated value is {}".format(vth))
    v_grid_centered = np.array([(v_grid[i+1] + v_grid[i])/ 2 for i in range(len(v_grid) - 1)]) * free_params['vte'][1]
    dv = np.diff(v_grid) * free_params['vte'][1] 
    # print("Wall velocity grid is {} and centered \n {} \n with dv of {}".format(v_grid * free_params['vte'][1] , v_grid_centered, dv))
    
    #Get velocity dependent collision frequency and subsequent mfp on entire grid 
    nu_ei, lamb_star = getLambdaStar(interp_Te, interp_ne, interp_Z, v_grid_centered)
    # print("The velocity dependent collision frequency is {} \n with reference backgroudn collision frequency i.e. vth particle {}".format(nu_ei[0],free_params['nu_ei'][0]))
    # print("The velocity dependent mfp is {} \n with reference backgroudn mfp i.e. vth particle {}".format(lamb_star[0],free_params['lambda_mfp'][0]))

    #Electric field only needs to be defined at cell- walls 
    #Every second entry in computer indicies 0 2 etc are cell-walls 
    #Electric field currently inaccurate
    E_field = electric_field(centered_x, interp_Te, interp_ne, interp_Z[::2])
    # plt.plot(x, E_field, "k--")
    free_params = free_param_calc(np.array([100]), np.array([1e19]), np.array([1]))
    true_E = np.loadtxt("/Users/shiki/DATA/E_FIELD_X_01000.txt") * free_params['E']
    # x = np.loadtxt("/Users/shiki/DATA/GRIDS/X_GRID.txt") * free_params['lambda_mfp']
    # plt.plot(x, true_E)
    # plt.show()
    E_field[:] = 0 #true_E[::2]
    #Likewise Lambda_E only needs to be defined at cell-walls 
    lamb_E = lambda_E(lamb_star[::2], E_field, v_grid_centered)
    #f0 should be defined on the entire grid 
    f0_mb = f_maxwellian(interp_Te, interp_ne, v_grid_centered, dv)
    #g_1 only on cell-walls 
    g_1_mb = g_1_maxwellian(lamb_star[::2], f0_mb, x_centered, interp_Te)
     #Lambda_E, g_1 needs to be cell-wall and rest cell-centereed
    delta_f0 = get_delta_f(g_1_mb, nu_ei[1::2, :], Z, lamb_E, v_grid_centered, full_x) 
    #integrate to get corrections
    q = integrate_g1(g_1_mb, v_grid_centered, dv) #Check if g_1 corresponds to spitzer-harm as it should.
    dq = integrate(delta_f0, lamb_E, v_grid_centered, dv, x_centered) # nonlocal deviation
    q_sh = spitzer_harm_heat(x_centered, Te, ne, Z) # local 

    #Fail if g_1 heat flow error larger than ~1.5 error scales with v-grid resolution% 
    #if this passes, correpsonds to f0_mb, g_1, nue_ei, lamb_star is correctly calculated
    # if any((abs(q - q_sh) / q_sh)[1:-1] * 100 > 1.5):
    #     sys.exit(0)

    q_snb = q - dq

    return q, q_sh, q_snb

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
    """
    Purpose: Creates a tridagional Matrix 
    Args:
        A = off diagonal term defined at n - 1
        B = leading diagonal defined at n 
        C = off diagonal term defined at n + 1
    Returns:
        matrix
    Checked: 13/08/2020
    """
    matrix = np.zeros((len(A), len(A)))
    for i, (a,b,c) in enumerate(zip(A,B,C)):
        if i == 0:
            matrix[i, 0] = b 
            matrix[i, 1] = c
        elif i == len(A) - 1:
            matrix[-1, -2] = a
            matrix[-1, -1] = b
        elif i < len(A) - 1:
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
    Z = np.zeros(nx, dtype = np.float64) + Z_
    Ar = np.zeros(nx, dtype = np.float64) + 2
    
    ne = np.zeros(nx, dtype = np.float64) + ne_
    
    ##Relevant for periodic systems .. fluid cannot be
    if sin: 
        Te = np.zeros(nx, dtype=np.float64) + Te_  + perturb * Te_  * np.sin((2*np.pi * x_centered) / np.max(x_centered), dtype=np.float64)
    else:
    ##Reflective
        Te = np.zeros(nx, dtype=np.float64) + Te_  + perturb * Te_  * np.cos((np.pi * x_centered) / np.max(x_centered), dtype=np.float64)


    return(initial_coord, x_centered, Te, ne, Z, Ar)

nx = 100
nv = 50
#klambda 0.0075
coord, centered_x, Te, ne, Z, Ar = epperlein_short(nx, 393.93740248643060431 * 2.81400056, Z_= 1, ne_ = 1e19, sin = False)
v_grid = np.linspace(0, 10, nv + 1)

q, q_sh, q_snb = snb_heat_flow(coord, v_grid, Te *(e/kb), ne, Z, Ar)
plt.plot(coord, q_sh, 'k-', label = "Spitzer")
plt.plot(coord, q, 'rx', label = "Apprixmated Spitzer")
plt.plot(coord, q_snb, label = "SNB")
# plt.plot(q_snb/q_sh)
plt.legend()
plt.show()