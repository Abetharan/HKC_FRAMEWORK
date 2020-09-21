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

def interpolateThermo(x, x_centered, Te, ne, Z, boundary_condition, normalised_values):
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
    if boundary_condition == 'periodic':
        sol_kit_grid = full_x_grid[:-1]
    
    elif boundary_condition == "reflective":
        sol_kit_grid = full_x_grid[1:-1]

    # sol_kit_grid = full_x_grid
    #Conver to SOL-KiT units
    sol_kit_ne = ne #/ (normalised_values['ne'])
    sol_kit_te = Te #* (kb/e)) / normalised_values['Te']
    sol_kit_z = Z
    
    #Require interpolation to get the centre quanties in SOL-KiT this is done via linear interpolations 
    #here we use cubic spline to smooth quanties. 
    sol_kit_inter_ne = np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_ne)
    sol_kit_inter_te = np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_te)
    sol_kit_inter_z =  np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_z)
    #Extend grid to cell-wall via interpolating using periodic bounds i.e. 
    #end cell-centre and first cell centre, distance from cell wall assumed to be the same
    #thus == 0.5
    if boundary_condition == 'periodic':
        sol_kit_inter_te[0] = 0.5 * (sol_kit_inter_te[1] + sol_kit_inter_te[-1])
        sol_kit_inter_ne[0] = 0.5 * (sol_kit_inter_ne[1] + sol_kit_inter_ne[-1])
        sol_kit_inter_z[0] = 0.5 * (sol_kit_inter_z[1] + sol_kit_inter_z[-1])


    return(sol_kit_grid, sol_kit_inter_te, sol_kit_inter_ne, sol_kit_inter_z)

def getDifferenceX(x, boundary_condition):
    """
    Purpose: Gets dx_k and dx_k_1_2
    Args:
        x = Full X grid 
        boundary_condition = boundary condition
    Returns: dx
    Notes:
    
    Periodic Boundary Grid looks like |.|.|.
    Thus, dx[1] represents dx_1 and goes to dx_3_2 afterwards

    Reflective Boundary Grid looks like .|.|.
    Thus, dx[1] represents dx_3_2 and dx_2 following that. 
   
    in both cases dx[0] and dx[-1] represent the spacing between first/last cell-centre and first/last cell-wall doubled to represent,
    whatever extension required
    """    
    nx = len(x)
    dxc = np.zeros(nx)
    for i in range(1, nx - 1):
        dxc[i] = x[i + 1] -  x[i - 1]
    dxc[0] = 2 * (x[1] - x[0])
    dxc[-1] = 2 * (x[-1] - x[-2])
    return dxc

def interpolateDistribution(full_x, dxc, f_0, f_1, E_field, boundary_condition):
    interpolated_f_0 = np.zeros(np.shape(f_0)) + f_0
    interpolated_f_1 = np.zeros(np.shape(f_1)) + f_1
    nx, nv = np.shape(f_0)
    interpolated_E_field = np.zeros(nx) + E_field
    for i in range(1, nx - 1):
        a_minus = (full_x[i + 1] - full_x[i]) / dxc[i]
        a_plus = (full_x[i] - full_x[i - 1]) / dxc[i]
        for v in range(nv):
            if i % 2 == 0:
                interpolated_f_0[i, v] = a_minus * f_0[i - 1, v] + a_plus * f_0[i + 1, v]
            elif i % 2 == 1:
                interpolated_f_1[i, v] = a_minus * f_1[i - 1, v] + a_plus * f_1[i + 1, v]
                interpolated_E_field[i] = a_minus * E_field[i - 1] + a_plus * E_field[i + 1]

    for i in range(nv):
        if boundary_condition == "periodic":
            interpolated_f_0[-1, i] = 0.5 * (f_0[0,i] + f_0[-2, i])
            interpolated_f_1[0, i] = 0.5 * (f_1[1,i] + f_1[-1,i])
            interpolated_E_field[0, i] = 0.5 * (E_field[1, i] + E_field[-1, i])
        else:
            interpolated_f_1[0, i]= 0.5 * f_1[1,i]
            interpolated_f_1[-1, i] = 0.5 * f_1[-2, i]
            interpolated_E_field[0] = 0.5 * (E_field[1])
            interpolated_E_field[-1] = 0.5 * (E_field[-2])

    return(interpolated_f_0, interpolated_f_1)

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
    
# def lambda_star(nue_ei, v_grid, Z):
#     lambda_star = np.zeros((len(nue_ei), len(v_grid)))
    # for i, local_nue in enumerate(nue_ei):
    #     zeta = (Z[i] + 0.24) / (Z[i] + 4.2)
        # for j, v in enumerate(v_grid):
        #     lambda_star[i, j] = zeta * (v/local_nue)
    # return(lambda_star)

def electric_field(x, Te, ne, Z, boundary_condition):
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
    if boundary_condition == "periodic":
        E = np.zeros(n)
    else:
        E = np.zeros(n - 1)
        n -=1
    #gamma_ee_0 = e  / pow((4 * math.pi* epsilon_0), 0.5)
    
    for i in range(0, n):
        if boundary_condition == "periodic" and i == 0:
            dx = 2 * x[i]
            gamma =  1 + ((3 * (Z[i] + 0.477)) / (2*(Z[i] + 2.15)))
            E[i] = (-1 * (kb/e) * Te[0] * 
                        (((ne[1] - ne[-1]) / (n * dx)) +
                        gamma * ((Te[1] - Te[-1]) /(Te[0] * dx))))
        else:
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
            if 1/(1/(lambda_star[i, j])) - lambda_star[i,j] > 1e-9:
                print("Da fuq {}".format(abs((e * e_field) / (0.5 * me * pow(v, 2)))))
            lambda_E[i, j] = 1/k

    return(lambda_E)

def f_maxwellian(Te, ne, v_grid, v_grid_width, boundary_condition):
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
        n = ne[i]
        for v in range(nv):
            f0[i,v] = n * ((np.pi * 2*v_th) ** (- 3.00/2.00)) * np.exp( - (v_grid[v] ** 2)/ (2*v_th**2))
            n_num[i] = n_num[i] + 4.00 * np.pi * v_grid[v] ** 2 * v_grid_width[v] * f0[i,v]
        
        #Ensure density is consistent via scaling with numerical density.
        for v in range(nv):
            f0[i, v] = f0[i,v] * ne[i]/n_num[i]
    return(f0)

def g_1_maxwellian(lambda_star, f_0_maxwellian, x, Te, boundary_condition):
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
    for i in range(0, nx):
        if boundary_condition == "periodic" and i == 0:
            dx = 2 * x[0]
            gradient = (Te[1] - Te[-1])/ dx
            temperature_correction = gradient/Te[0]
        elif boundary_condition == "periodic":
            dx = x[i] - x[i - 1]
            gradient = (Te[i*2 + 1] - Te[i*2 - 1])/ dx
            temperature_correction = (gradient/Te[i*2])
        else:
            dx = x[i + 1] - x[i]
            gradient = (Te[i*2 + 2] - Te[i*2])/ dx
            temperature_correction = (gradient/Te[i*2 + 1])

        # f_0_maxwellian[0] = 0.5 * (f_0_maxwellian[1, :] +  f_0_maxwellian[-1, :])
        for v in range(nv):
            if boundary_condition == "periodic":
                g_1[i, v] = -1 * lambda_star[i, v] * f_0_maxwellian[2*i, v] * temperature_correction
            else:
                g_1[i, v] = -1 * lambda_star[i, v] * f_0_maxwellian[2*i + 1, v] * temperature_correction
    return(g_1)

def forward_difference(f_3_2, f_1_2, dx, A = None):
    if A is None:
        return (f_3_2 - f_1_2)/dx
    else:
        return A * ((f_3_2 - f_1_2) / dx)
def solve_df0(g_1, nue_ei, Z, lambda_E, v_grid, dx, boundary_condition):
    """
    Purpose: Fill Implicit matrix 
    Args:
        g_1: modified SNB f_1 maxwellian, shape (nx(-1), nv)
        nue_ei: Collision Frequency, shape (nx, nv)
        Z: Ionization, shape nx 
        Lambda_E: SNB modified collision mfp shape (nx(-1), nv) 
        v_grid: Velocity grid/groups, shape (nv)
        x: full X grid has length,  2*nx(-1)
    Returns:
        matrix to be solved.
    """ 
    nx, nv = np.shape(nue_ei)
    delta_f0 = np.zeros((nx, nv))
    r = 2 # tunable parameter
    for j, v in enumerate(v_grid):
        A = np.zeros(nx) #off diagonal n -1
        B = np.zeros(nx) # leading diagonal
        C = np.zeros(nx) # off diagonal n + 1
        S = np.zeros(nx) #Source 
        for i in range(nx):
            S_has_been_solved = False
            if boundary_condition == "periodic":
                if i < nx -1:
                    dx_k = dx[2*i]
                    dx_k_1 = dx[2*i + 2]
                    dx_k_1_2 = dx[2*i + 1]
                    dx_k_diff = dx_k + dx_k_1
                if i == nx - 1:
                    dx_k = dx[2*i]
                    dx_k_1 = dx[-1]
                    dx_k_1_2 = dx[2*i + 1]
                    dx_k_diff = dx_k + dx_k_1
                    A[i] = lambda_E[i, j]/(3*dx_k_diff*dx_k)
                    B[i] = -1 * (lambda_E[0, j]/(3*dx_k_diff*dx_k_1) + lambda_E[i, j]/(3*dx_k_diff*dx_k) + (r * nue_ei[i, j])/(v*Z[i]))
                    C[i] = lambda_E[0, j]/(3*dx_k_diff*dx_k)
                    S[i] = (g_1[0, j] - g_1[i, j]) / (3*dx_k_1_2)
                    if S[i] == 0:
                        S_has_been_solved = True

            if boundary_condition == "reflective":
                if i < nx - 1 and i > 0:
                    dx_k = dx[2*i - 1]
                    dx_k_1 = dx[2*i + 1]
                    dx_k_1_2 = dx[2*i]
                    dx_k_diff = dx_k + dx_k_1

                elif i == 0:
                    dx_k_1_2 = dx[0] #2 * (x[1] - x[0])
                    dx_k = dx[0] # Extrapolating x[1] - x[0] ends up == dx_k_1_2
                    dx_k_1 = dx[1] # k = 1, (3/2 - 1/2) x[2]  - x[0]
                    dx_k_diff = dx_k + dx_k_1
                    A[i] = None #0 #lambda_E[i - 1]/(3*dx_k_diff*dx_k)
                    B[i] = -1 * (lambda_E[i, j]/(3*dx_k_diff*dx_k) + (r * nue_ei[i, j])/(v*Z[i]))
                    C[i] = lambda_E[i, j]/(3*dx_k_diff*dx_k)
                    S[i] = g_1[i, j] / (3*dx_k_1_2)
                if i == nx - 1:
                    dx_k = dx[-2]
                    dx_k_1 = dx[-1]
                    dx_k_1_2 = dx[-1]
                    dx_k_diff = dx_k + dx_k_1
                    A[i] = lambda_E[i - 1, j]/(3*dx_k_diff*dx_k)
                    B[i] = -1 * (lambda_E[i - 1, j]/(3*dx_k_diff*dx_k) + (r * nue_ei[i,j])/(v*Z[i]))
                    C[i] = None #lambda_E[i]/(3*dx_k_diff*dx_k)
                    S[i] = -1*g_1[i - 1, j] / (3*dx_k_1_2)
                    if S[i] == 0:
                        S_has_been_solved = True
                    
            if A[i] == 0:
                if boundary_condition == "periodic":
                    A[i] = lambda_E[i, j]/(3*dx_k_diff*dx_k)  
                else:
                    A[i] = lambda_E[i - 1, j]/(3*dx_k_diff*dx_k)
            if B[i] == 0:
                if boundary_condition == "periodic":
                    B[i] = -1 * (lambda_E[i + 1, j]/(3*dx_k_diff*dx_k_1) + lambda_E[i, j]/(3*dx_k_diff*dx_k) + (r * nue_ei[i,j])/(v*Z[i]))
                else:
                    B[i] = -1 * (lambda_E[i, j]/(3*dx_k_diff*dx_k_1) + lambda_E[i - 1, j]/(3*dx_k_diff*dx_k) + (r * nue_ei[i,j])/(v*Z[i]))

            if C[i] == 0:
                if boundary_condition == "periodic":
                    C[i] = lambda_E[i + 1, j]/(3*dx_k_diff*dx_k_1) 
                else:
                    C[i] = lambda_E[i, j]/(3*dx_k_diff*dx_k)
            
            if S[i] == 0 and not S_has_been_solved:
                if boundary_condition == "periodic":
                    S[i] = forward_difference(g_1[i + 1, j] , g_1[i, j], dx_k_1_2, 1/3) 
                else:
                    S[i] = forward_difference(g_1[i, j], g_1[i - 1, j], dx_k_1_2, 1/3) 

        mat = createMatrix(A, B, C)
        if boundary_condition == "periodic":
            mat[0, -1] = A[0]
            mat[-1, 0] = C[-1]

        d = solve(mat, S)
        delta_f0[:, j] = d
    return delta_f0

def integrate(delta_f0, lambda_E, v_grid, v_grid_width, x_grid, boundary_condition):
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
    dq = np.zeros(np.shape(lambda_E)[0])

    if boundary_condition == "periodic":
        end = nx
    else:
        end = nx - 1
    for i in range(0, end):
        heat_flow = 0
        if boundary_condition == "periodic":
            if i == 0:
                grad_delta_f0 = ((delta_f0[i, :] - delta_f0[-1, :]) / dx[0]) #Spatial derivative
            else:
                grad_delta_f0 = ((delta_f0[i, :] - delta_f0[i - 1, :]) / dx[i - 1]) #Spatial derivative
        else:
            grad_delta_f0 = ((delta_f0[i + 1, :] - delta_f0[i, :]) / dx[i]) #Spatial derivative
        #Integrate
        for v in range(nv):
            if boundary_condition == "periodic":
                heat_flow +=  pow(v_grid[v], 5) * v_grid_width[v] * lambda_E[i, v] * grad_delta_f0[v]
            else:
                heat_flow +=  pow(v_grid[v], 5) * v_grid_width[v] * lambda_E[i, v] * grad_delta_f0[v]
        
        dq[i] = -1*((2 * np.pi * me)/3) * heat_flow

    return dq

def integrate_g1(g_1, v_grid, v_grid_width):
    nx, nv = np.shape(g_1)
    q = np.zeros(nx)
    for i in range(0, nx):
        heat_flow = 0
        for v in range(nv):
            heat_flow +=  pow(v_grid[v], 5) * v_grid_width[v] * g_1[i, v]
        q[i] = ((2 * np.pi * me)/3) * heat_flow
    return q 

def snb_heat_flow(x, v_grid, Te, ne , Z, boundary_condition, norms = 0):
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
    full_x, interp_Te, interp_ne, interp_Z = interpolateThermo(x, x_centered, Te, ne, Z, boundary_condition, norms)
    dxc = getDifferenceX(full_x, boundary_condition)
    free_params= free_param_calc(interp_Te * (kb/e), interp_ne ,interp_Z )
    # print("Vth scaling factor is {}".format(free_params['vte'][1]))
    # vth = np.sqrt(100 * e /me)
    # print("Hand Calculated value is {}".format(vth))
    v_grid_centered = np.array([(v_grid[i+1] + v_grid[i])/ 2 for i in range(len(v_grid) - 1)]) 
    dv = np.diff(v_grid) 
    # print("Wall velocity grid is {} and centered \n {} \n with dv of {}".format(v_grid * free_params['vte'][1] , v_grid_centered, dv))
    
    #Get velocity dependent collision frequency and subsequent mfp on entire grid 
    nu_ei, lamb_star = getLambdaStar(interp_Te, interp_ne, interp_Z, v_grid_centered)

    #Electric field only needs to be defined at cell- walls 
    #Every second entry in computer indicies 0 2 etc are cell-walls 
    #Electric field currently inaccurate
   
    if boundary_condition =="periodic":
        E_field = electric_field(centered_x, interp_Te, interp_ne, interp_Z[::2], boundary_condition)
        #Likewise Lambda_E only needs to be defined at cell-walls 
        lamb_E = lambda_E(lamb_star[::2, :], E_field, v_grid_centered)
        #f0 should be defined on the entire grid 
        f0_mb = f_maxwellian(interp_Te, interp_ne, v_grid_centered, dv, boundary_condition)
        #g_1 only on cell-walls 
        g_1_mb = g_1_maxwellian(lamb_star[::2, :], f0_mb, x_centered, interp_Te, boundary_condition)
        delta_f0 = solve_df0(g_1_mb, nu_ei[1::2, :], Z, lamb_E, v_grid_centered, dxc, boundary_condition) 
    else:
        E_field = electric_field(centered_x, interp_Te, interp_ne, interp_Z[1::2], boundary_condition)
        #Likewise Lambda_E only needs to be defined at cell-walls 
        lamb_E = lambda_E(lamb_star[1::2, :], E_field, v_grid_centered)
        #f0 should be defined on the entire grid 
        f0_mb = f_maxwellian(interp_Te, interp_ne, v_grid_centered, dv, boundary_condition)
        #g_1 only on cell-walls 
        g_1_mb = g_1_maxwellian(lamb_star[1::2, :], f0_mb, x_centered, interp_Te, boundary_condition)

        delta_f0 = solve_df0(g_1_mb, nu_ei[::2, :], Z, lamb_E, v_grid_centered, dxc, boundary_condition) 
        # fig = plt.figure()
        # x_mesh, v_mesh = np.meshgrid(x_centered, v_grid_centered/free_params['vte'][0])
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_surface(x_mesh, v_mesh, delta_f0)
        # ax.set_xlabel('x/m')
        # ax.set_ylabel('v/v_th')
        # ax.set_zlabel('g_1')
        # plt.show()
    #integrate to get corrections
    q = integrate_g1(g_1_mb, v_grid_centered, dv) #Check if g_1 corresponds to spitzer-harm as it should.
    # plt.plot(q)
    dq = integrate(delta_f0, lamb_E, v_grid_centered, dv, x_centered, boundary_condition) # nonlocal deviation
    q_sh = spitzer_harm_heat(x_centered, Te, ne, Z) # local 
    # plt.plot(q_sh)
    # plt.show()
    

    #Fail if g_1 heat flow error larger than ~1.5 error scales with v-grid resolution% 
    #if this passes, correpsonds to f0_mb, g_1, nue_ei, lamb_star is correctly calculated
    # if any((abs(q - q_sh) / q_sh)[1:-1] * 100 > 1.5):
    #     sys.exit(0)
    # q_snb = np.insert(q_snb, 0 ,0)
    # q_snb = np.append(q_snb, 0)
    
    
    if boundary_condition =="periodic":
        return full_x[::2], q, q_sh, q + dq
    else:
        return full_x[1::2], q, q_sh, q + dq

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

    # kappaE = 1.843076667547614E-10 * pow(Te, 2.5) * pow(coulomb_log, -1) *  pow(Z, -1)
    kappaE =  13.6*((Z+0.24)/(Z+4.24))*5.759614586E-11* pow(Te, 2.5) * pow(coulomb_log, -1) *  pow(Z, -1)
    nx = len(Te)
    HeatFlowE = np.zeros(nx + 1, dtype=np.float64)
    gradTe = np.zeros(nx + 1)
    for i in range(1, nx):
        centered_ke = 0.5 * (kappaE[i] + kappaE[i - 1])
        gradTe[i] = ((Te[i] - Te[i - 1]) / (x[i] - x[i - 1]))
        HeatFlowE[i] = centered_ke * gradTe[i] 
    
    # HeatFlowE[0] = 0
    # HeatFlowE[-1] = 0
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
    initial_coord = np.zeros(nx + 1)
    v_grid_dv = 0.607
    for i in range(1, nx + 1):
        v_grid_dv *= 1.07
        initial_coord[i] = initial_coord[i - 1] + v_grid_dv
    # initial_coord1 = np.linspace(x_l, 0.8*x_u, math.ceil((nx + 1)*0.6), dtype=np.float64)
    # initial_coord2 = np.linspace(0.8*x_u, x_u, int((nx+1) *0.4), dtype=np.float64)
    # initial_coord = np.concatenate((initial_coord1,initial_coord2))
    initial_coord = np.linspace(x_l, x_u, nx +1)
    x_centered = np.array([(initial_coord[i] + initial_coord[i+1]) /2 for i in range(len(initial_coord) -1)], dtype=np.float64)
    Z = np.zeros(nx, dtype = np.float64) + Z_
    Ar = np.zeros(nx, dtype = np.float64) + 4
    
    ne = np.zeros(nx, dtype = np.float64) + ne_
    
    ##Relevant for periodic systems .. fluid cannot be
    if sin: 
        Te = np.zeros(nx + 1, dtype=np.float64) + Te_  + perturb * Te_  * np.sin((2*np.pi * initial_coord) / np.max(initial_coord), dtype=np.float64)
    else:
    ##Reflective
        Te = np.zeros(nx + 1, dtype=np.float64) + Te_  + perturb * Te_  * np.cos((np.pi * initial_coord) / np.max(initial_coord), dtype=np.float64)
    Te_centered = np.array([(Te[i] + Te[i+1]) /2 for i in range(len(Te) -1)], dtype=np.float64) 

    return(initial_coord, x_centered, Te_centered, ne, Z, Ar)

def test_finite_difference_scheme(x, x_centered, v_grid, Te, ne, Z, boundary_condition, dt = 1e-15):
    C_v = np.zeros(len(Te)) + 100
    Te_intemediate = np.zeros(len(Te))
    Te_new_time = np.zeros(len(Te))
    spitzer_heat = spitzer_harm_heat(x_centered, Te, ne, Z)
    dx = np.diff(x)
    for i in range(len(Te)):    
        thermal_conduction = (spitzer_heat[i + 1] - spitzer_heat[i]) / dx[i]
        Te_intemediate[i] = Te[i] - (dt/C_v[i]) * thermal_conduction
    
    new_spitzer = spitzer_harm_heat(x_centered, Te_intemediate, ne, Z)
    x_nl, q, q_sh, q_snb = snb_heat_flow(x, v_grid, Te_intemediate, ne, Z, boundary_condition)
    
    if boundary_condition == "reflective":
        q_snb = np.insert(q_snb, 0 ,0)
        q_snb = np.append(q_snb, 0)
    else:
        q_snb = np.append(q_snb, 0)

    for i in range(len(Te)):    
        spitzer_thermal_conduction = (new_spitzer[i + 1] - new_spitzer[i]) / dx[i]
        nl_thermal_conduction =(q_snb[i + 1] - q_snb[i]) / dx[i]
        Te_new_time[i] = Te[i] - (dt/C_v[i]) * (-1 * nl_thermal_conduction + spitzer_thermal_conduction)

    plt.figure(2) 
    plt.plot(x, q_snb, label = "snb")
    plt.plot(x, spitzer_heat, label = "spitzer")
    plt.plot(x_nl, q, label = "g_1")
    plt.legend()
    plt.figure(1)
    plt.plot(x_centered, Te_intemediate)
    plt.plot(x_centered, Te_new_time)
    plt.legend()
    plt.show()





nx = 10
nv = 100

# klambda_dist = [393.7402486430605, 147.65259324114768, 73.82629662057384, 39.37402486430605, 14.765259324114767, 7.382629662057384, 3.9374024864306043, 1.4765259324114768]
klambda_dist = [393.7402486430605, 147.65259324114768, 73.82629662057384, 39.37402486430605]
klambda = [0.0075,0.02, 0.04, 0.075]
sh_q = []
sh_analytical = []
for k, dist in enumerate(klambda_dist):
    coord, centered_x, Te, ne, Z, Ar = epperlein_short(nx, dist*2.81 , Z_= 1, ne_ = 1e19, sin = False)
    # coord = np.loadtxt('/Users/shiki/DATA/spitzer_checks/OUTPUT/GRIDS/X_GRID.txt')
    # Z = np.zeros(len(coord)) + 36.5
    # Te = np.loadtxt('/Users/shiki/DATA/spitzer_checks/OUTPUT/TEMPERATURE/TEMPERATURE_00000.txt')
    # ne = np.zeros(len(coord)) + 1e19
    # coordd = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_Z_interp', skiprows=1)[:, 0] 
    # Z = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_Z_interp', skiprows=1)[:, 1] 
    # Te = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_5ps_TekeV_interp', skiprows=1)[:, 1] * 1E3 
    # ne = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_ne1e20cm3_interp', skiprows=1)[:, 1] * (1e20 * 1e6)
    # snb_brodrick = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_5ps_separatedsnbWcm2', skiprows=1)[:, 1]
    # coord = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_5ps_separatedsnbWcm2', skiprows=1)[:, 0]*1e-6
    # plt.plot(coord, snb_brodrick * 1e4, label = 'true')

    centered_x = np.array([(coord[i+1] + coord[i]) /2 for i in range(len(coord) -1)])
    free_params = free_param_calc(np.array([100]), np.array([1e19]) ,np.array([1]))

    v_grid = np.zeros(nv + 1)
    v_grid_dv = 0.0307
    for j in range(1, nv + 1):
        v_grid_dv *= 1.04
        v_grid[j] = v_grid[j - 1] + v_grid_dv
    # v_grid = np.linspace(0, 30, nv + 1) * free_params['vte'][0]
    v_grid *=free_params['vte'][0] 
    boundary_condition = "reflective"
    # boundary_condition = "periodic"
    # test_finite_difference_scheme(coord, centered_x, v_grid, Te * (e/kb), ne, Z, boundary_condition, dt = 1e-15)

    x, q, q_sh, q_snb = snb_heat_flow(coord, v_grid, Te* (e/kb), ne, Z, boundary_condition)
    sh_q.append(np.average(abs(q_snb/q_sh[1:-1])))
    q_analy = q_sh *(1 - 233.26 * 1 * klambda[k]**2)
    sh_analytical.append(np.average(abs(q_analy[1:-1]/q_sh[1:-1])))
    plt.plot(coord[1:-1], q_sh[1:-1], 'k-', label = "Spitzer")
    plt.plot(x, q, 'r-', label = "Apprixmated Spitzer")
    plt.plot(x, q_snb, label = "SNB")
#     # plt.plot(q_snb/q_sh[1:-1])
    plt.legend()
    plt.show()

# plt.plot(klambda, sh_q, 'k-')