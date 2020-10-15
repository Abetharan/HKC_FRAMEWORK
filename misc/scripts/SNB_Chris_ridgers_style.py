#Only for reflective boundary conditions and nothign else
import numpy as np 
from scipy import linalg
import matplotlib.pyplot as plt 
from scipy import constants
import math
import sys
import klambda_finder as kf
from scipy.optimize import least_squares
from scipy.optimize import minimize
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
        # zeta = Z / (2 * pow((Z+1), 0.5)) #Schurtz
        zeta = (Z + 0.24) / (Z + 4.2) # Brodrick
        lambda_star[:, j] = (v / nue_ei[:, j]) * zeta 
    return lambda_star
    
def electric_field(dx_wall, Te, ne, Z):
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
    n = len(Z)
    E = np.zeros(n - 1)
    for i in range(1, n):
        dx = dx_wall[i - 1]
        gamma =  1 + ((3 * (Z[i] + 0.477)) / (2*(Z[i] + 2.15)))
        E[i - 1] = (-1 * (kb/e) * Te[i*2] * 
                    (((ne[i*2] - ne[i*2 - 2]) / (ne[2*i - 1] * dx)) +
                    gamma * ((Te[i*2] - Te[i*2 - 2]) /(Te[2*i - 1] * dx))))
    return(E)


def lambda_E(lambda_star, E, e_grid):
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
    a = 1.0
    for i, e_field in enumerate(E):
        for j, energy in enumerate(e_grid):
            k = 1/(a * lambda_star[i, j]) + abs((e * e_field) / (kb * energy))
            lambda_E[i, j] = 1/k

    return(lambda_E)

def getLambStar(Te, ne, Z, beta):
    nx,ng = np.shape(beta)
    vt_bx = np.sqrt(kb*Te/me)
    log_lam_ei = lambda_ei(T_norm = Te * (kb/e), n_norm = ne, Z_norm = Z)

    # log_lam_ei = np.zeros(np.shape(log_lam_ei)) + 2.1484
    Y = 4*np.pi*(e**2/(4 * np.pi * epsilon_0 * me))**2
    tau_ei_bx = vt_bx**3/(Y* Z * ne * log_lam_ei)

    lambda_ei_bx = vt_bx*tau_ei_bx
    zeta = (Z + 0.24) / (Z + 4.2) # Brodrick
    # lambda_e_bx = np.sqrt(Z)*lambda_ei_bx
    lambda_e_bx = zeta * lambda_ei_bx
    lambda_g_bx = np.zeros((nx,ng))
    lambda_g_c =  np.zeros((nx,ng))

    for i in range(ng):
        lambda_g_bx[:,i] = 2*pow(beta[:,i], 2)*lambda_e_bx #Defined everywhere in space,  energy centered 
    return lambda_g_bx

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
    # coulomb_log[:] = 2.1484
    # kappaE = 1.843076667547614E-10 * pow(Te, 2.5) * pow(coulomb_log, -1) *  pow(Z, -1)
    kappaE =  13.6*((Z+0.24)/(Z+4.2))*5.759614586E-11* pow(Te, 2.5) * pow(coulomb_log, -1) *  pow(Z, -1)
    nx = len(Te)
    HeatFlowE = np.zeros(nx + 1, dtype=np.float64)
    gradTe = np.zeros(nx + 1)
    for i in range(1, nx):
        centered_ke = 0.5 * (kappaE[i] + kappaE[i - 1])
        gradTe[i] = ((Te[i] - Te[i - 1]) / (x[i] - x[i - 1]))
        HeatFlowE[i] = centered_ke * gradTe[i] 
    
    # HeatFlowE[0] = 0
    # HeatFlowE[-1] = 0
    return(-1 * HeatFlowE[1:-1])

def getEnergyGroups(velocity_grid):
    return(0.5*me*pow(velocity_grid,2))

# def getBeta(wall_energy_groups, Te):
#     #beta = E_g/kT_e
#     beta = np.zeros((len(Te),len(wall_energy_groups)))
#     dbeta = np.zeros((len(Te),len(wall_energy_groups) - 1)) 
#     for i in range(len(Te)):
    #     beta[i,:] = wall_energy_groups / (kb * Te[i])
    #     for j in range(len(wall_energy_groups) - 1):
    #         dbeta[i, j] =  (wall_energy_groups[j + 1] -wall_energy_groups[j]) / (kb * Te[i])
    # return(beta, dbeta)

def getBeta(energy_grid, Te):
    beta = np.zeros((len(Te),len(energy_grid)))
    dbeta = np.zeros((len(Te),len(energy_grid) - 1)) 
    for i in range(len(Te)):
        beta[i, :] = energy_grid / Te[i]
        for j in range(len(energy_grid) - 1):
            dbeta[i, j] = (energy_grid[j + 1] - energy_grid[j]) /  Te[i]
    return beta, dbeta
        
def getSource(q_sh, beta, dbeta):
    nx, ne = np.shape(beta)
    U_g = np.zeros((nx, ne))
    for i in range(nx):
        for e in range(ne):
            U_g[i, e] = (q_sh[i] / 24) * dbeta[i, e] * pow(beta[i, e], 4) * np.exp(-1*beta[i, e])
    return U_g
     
def getH(dx_wall, dx_centered, lambda_star, lambda_E, U, r, Z):
    nx, nv = np.shape(lambda_star)
    H = np.zeros((nx, nv))
    r = np.zeros(nx) + r#  tunable parameter
    for j in range(nv):
        A = np.zeros(nx) #off diagonal n -1
        B = np.zeros(nx) # leading diagonal
        C = np.zeros(nx) # off diagonal n + 1
        S = np.zeros(nx) #Source 
        for i in range(nx):
            if i == 0:
                dx_diff = dx_centered[i]#dx_wall[1] + dx_wall[0] #
                A[i] = None #0 #lambda_E[i - 1]/(3*dx_k_diff*dx_k)
                B[i] = 1 * (lambda_E[i, j]/(3*dx_wall[i] * dx_diff) + r[i] / (lambda_star[i, j] * Z[i]))
                C[i] = -1*lambda_E[i, j]/(3*dx_wall[i] * dx_diff)
                S[i] = -1* U[i, j] / (dx_centered[i])
            elif i == nx - 1:           
                dx_diff = dx_centered[i]#dx_wall[-1] + dx_wall[-2] #
                A[i] = -1*lambda_E[i - 1, j]/(3*dx_wall[i - 1] * dx_diff)
                B[i] = 1 * (lambda_E[i - 1, j]/(3*dx_wall[i - 1] * dx_diff) + r[i] / (lambda_star[i, j] * Z[i]))
                C[i] = None #lambda_E[i]/(3*pow(dx_k)
                S[i] = 1*U[i - 1, j] / (dx_centered[-1])
            else:
                dx_diff = dx_centered[i]# dx_wall[i + 1] + dx_wall[i] 
                A[i] = -1*lambda_E[i - 1, j]/(3*dx_wall[i - 1] * dx_diff)
                B[i] = 1 * (lambda_E[i, j]/(3*dx_wall[i] * dx_diff) + lambda_E[i - 1, j]/(3*dx_wall[i - 1] * dx_diff) + r[i] / (lambda_star[i, j] * Z[i]))
                C[i] = -1*lambda_E[i, j]/(3*dx_wall[i] * dx_diff)
                S[i] = -1 * (U[i, j] - U[i - 1, j]) / dx_centered[i]

        mat = createMatrix(A, B, C)
        d = solve(mat, S)
        H[:, j] = d
    return(H)

def getDqnl(H, lambda_E, dx_wall):
    nx, nv = np.shape(H)
    q_nl = np.zeros(nx - 1)
    for i in range(0, nx - 1):
        heat_flow = 0        
        grad_delta_f0 = ((H[i + 1, :] - H[i, :]) / dx_wall[i])
        for v in range(nv):
            heat_flow += lambda_E[i, v] * grad_delta_f0[v] 
        q_nl[i] = (1/3) * heat_flow
    return q_nl 

def snb_heat_flow(x, e_grid_wall, Te, ne , Z, r):
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
    full_x, interp_Te, interp_ne, interp_Z = interpolateThermo(x, x_centered, Te, ne, Z, 'reflective', None) 
    dx_centered = np.diff(full_x[1::2])
    dx_wall = np.diff(full_x[::2])
    dx_centered = np.insert(dx_centered, 0 , 2 * (full_x[1] - full_x[0]))
    dx_centered = np.append(dx_centered,  2 * (full_x[-1] - full_x[-2]))
    e_grid_centre = np.array([(e_grid_wall[i+1] + e_grid_wall[i])/ 2 for i in range(len(e_grid_wall) - 1)])

    beta, _ = getBeta(e_grid_centre, interp_Te)
    _, dbeta = getBeta(e_grid_wall, interp_Te)

    # lambda_star = getLambdaStar(interp_Te, interp_ne, interp_Z, beta)
    E_field = electric_field(dx_wall, interp_Te, interp_ne, Z)
    lambda_star = getLambStar(interp_Te, interp_ne, interp_Z, beta)
    lamb_E = lambda_E(lambda_star[1::2], E_field, e_grid_centre)
    # lamb_E = lambda_star[1::2] 
    
    q_sh = spitzer_harm_heat(x_centered, Te, ne, Z) # local 
    U = getSource(q_sh, beta[1::2, :], dbeta[1::2, :])
    H = getH(dx_wall, dx_centered, lambda_star[::2], lamb_E, U, r, Z)
    q_nl = getDqnl(H, lamb_E, dx_wall)
    
    q_snb = q_sh  - q_nl
    # q_snb = np.append(q_snb, 0)
    # q_snb = np.insert(q_snb, 0 ,0)
    return q_sh, q_snb  

def fitfunc(x0, x, e_grid_wall, Te, ne , Z, true_qe):
    _, q_snb = snb_heat_flow(x, e_grid_wall, Te, ne , Z, x0)
    multi = q_snb / true_qe 
    return abs(1 - multi)

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
    Ar = np.zeros(nx, dtype = np.float64) + 2
    
    ne = np.zeros(nx, dtype = np.float64) + ne_
    
    ##Relevant for periodic systems .. fluid cannot be
    if sin: 
        Te = np.zeros(nx + 1, dtype=np.float64) + Te_  + perturb * Te_  * np.sin((2*np.pi * initial_coord) / np.max(initial_coord), dtype=np.float64)
    else:
    ##Reflective
        Te = np.zeros(nx + 1, dtype=np.float64) + Te_  + perturb * Te_  * np.cos((np.pi * initial_coord) / np.max(initial_coord), dtype=np.float64)
    Te_centered = np.array([(Te[i] + Te[i+1]) /2 for i in range(len(Te) -1)], dtype=np.float64) 

    return(initial_coord, x_centered, Te_centered, ne, Z, Ar)





nx = 64 
nv = 70
#klambda_dist = [393.7402486430605, 147.65259324114768, 73.82629662057384, 39.37402486430605, 14.765259324114767, 7.382629662057384, 3.9374024864306043, 1.4765259324114768, 0.7382629662057384, 0.49217531080382554]
# klambda_dist = [393.7402486430605, 147.65259324114768, 73.82629662057384, 39.37402486430605, 14.765259324114767, 7.382629662057384, 3.9374024864306043, 1.4765259324114768]
# klambda = [0.0075,0.02, 0.04, 0.075, 0.2, 0.4, 0.75, 2]
#klambda = [0.0075,0.02, 0.04, 0.075, 0.2, 0.4, 0.75, 2, 4, 6]
klambdas = np.geomspace(7.699668410015750308e-03, 2.170061402318336332e+04, 100)
klambda_dist = kf.find_domain_length(klambdas, wavenumber = 0.5)
sh_q = []
for i in range(len(klambda_dist)):#2.81 Z = 1, 1.492 Z = 8
    # coord, centered_x, Te, ne, Z, Ar = epperlein_short(nx, klambda_dist[i]*1.492 , Z_= 8, ne_ = 1e19, sin = False)

    # centered_x = np.array([(coord[i+1] + coord[i]) /2 for i in range(len(coord) -1)])

    Z = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_Z_interp', skiprows=1)[:, 1] 
    Te = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_5ps_TekeV_interp', skiprows=1)[:, 1] * 1E3 *(e/kb)
    ne = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_ne1e20cm3_interp', skiprows=1)[:, 1] * (1e20 * 1e6)
    snb_brodrick = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_5ps_separatedsnbWcm2', skiprows=1)[:, 1]
    impact_brodrick = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_5ps_IMPACTWcm2', skiprows=1)[:, 1]
    coord = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/gdhohlraum_xmic_5ps_separatedsnbWcm2', skiprows=1)[:, 0]*1e-6
    Hykict_snb = np.loadtxt('/Users/shiki/DATA/HyKiCT_OUTPUT/ELECTRON_HEAT_FLOW_X/ELECTRON_HEAT_FLOW_X_0.txt') *-1e-4
    # # v_grid = np.zeros(nv + 1)
    #v_grid_dv = 0.0207
    #v_max = np.sqrt(2 * (200 * e * np.max(Te)) / me)
    # dv = v_max / nv
    # for i in range(1, nv + 1):
        # v_grid_dv *= 1.02
        # v_grid[i] = v_grid[i - 1] + v_grid_dv
    E_max = np.max(Te)* 20
    e_grid = np.linspace(0, E_max, nv + 1)
    # e_grid[0] = 0  
    # v_grid *=  np.sqrt((Te[0] * e)/me)
    x0 = [2]
    q_sh, q_snb = snb_heat_flow(coord, e_grid, Te, ne, Z, 2)
    # res_lsq = least_squares(fitfunc, x0, args =(coord, e_grid, Te* (e/kb), ne, Z, impact_brodrick))
    # print(res_lsq.x)
    # print(q_snb)
    # sh_q.append(np.average(q_snb/q_sh))
    # plt.plot(coord[1:-1], q_sh*1e-4, 'k-', label = "Spitzer")
    plt.plot(coord, snb_brodrick, label = 'brodrick')
    plt.plot(coord[1:-1], impact_brodrick, label = 'impact')
    plt.plot(coord[1:-1], q_snb*1e-4, label = "SNB")
    plt.plot(coord, Hykict_snb, label = "hykict SNB")
    # plt.plot(coord[1:-1], (1e-4*q_snb[1:-1]) / impact_brodrick)
    plt.legend()
    plt.show()

# klambda_br =  np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/linearised_z8_snbr2', skiprows=1)[:, 0] 
# q_q_sh_brodrick = np.loadtxt('/Users/shiki/DATA/Brodrick_2017_data/linearised_z8_snbr2', skiprows=1)[:, 1]
# plt.loglog(klambda_br, q_q_sh_brodrick, label = "Brodrick 2017")
# plt.loglog(klambdas, sh_q, '--', label = "Me")
# plt.xlabel('Klambda')
# plt.ylabel('q/q_sh')
# plt.legend()
# # plt.savefig('/Users/shiki/Documents/Documents – Abetharan’s MacBook Pro/Imperial_College_London/Ph.D./HKC_FRAMEWORK/Documentation/brodrick_vs_me_linearised_z1_snbr2.pdf')
# plt.show()