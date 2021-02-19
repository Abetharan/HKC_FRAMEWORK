import math
import numpy as np
import sys
from scipy import linalg
from scipy import constants
kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
epsilon_0 = 8.854188E-12    # Vacuum dielectric constant
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")
class HeatFlowCouplingTools:
    """
    A set of utility tools that is to be used in the calculation of 
    coupling parameters.
    Tools:
        Spitzer-Harm calculation
        Coulomb Log relevant to Spitzer-harm calculation
        Front-heat parameter finder
        Pre-heat parameter finder
        Div.q calculator 
    """

    def __init__(self):
        self.cell_centered_coord =  np.array([])
        self.cell_wall_coord = np.array([])
        self.electron_temperature = np.array([])
        self.electron_number_density = np.array([])
        self.zbar = np.array([])
        self.mass = np.array([])
        self.coulomb_log = None 
        self.spitzer_harm_heat = np.array([])
        self.vfp_heat = np.array([])
        self.q_vfp_q_sh_multipliers = np.array([])
        self.search_tolerance = 1e-9
        self.snb = False
        self.q_snb = None

    def lambda_ei(self, T_norm , n_norm, Z_norm, return_arg = False, return_array = False):
        coulomb_logs = []
        T_norm = T_norm

        for T,n,Z in zip(T_norm, n_norm, Z_norm):
            if T < 10.00 * Z ** 2:
                result = 23.00 - math.log(math.sqrt(n * 1.00E-6) * Z * (T) ** (-3.00/2.00))
            else:
                result = 24.00 - math.log(math.sqrt(n * 1.00E-6) / (T))   
            if result < 2.5:
                result = 2.5
            if result < 0:
                print("NEGATIVE COULOMB LOGS")
                sys.exit(0)


            if return_arg:
                return result
            else:
                coulomb_logs.append(result)
        if return_array:
            return coulomb_logs
        else:
            self.coulomb_log = np.array(coulomb_logs)

    def spitzerHarmHeatFlow(self):
        """
        Purpose: Models Spitzer Harm heat flow 
        """
        kappaE =  13.6*((self.zbar+0.24)/(self.zbar+4.2))*5.759614586E-11* pow(self.electron_temperature, 2.5) * pow(self.coulomb_log, -1) *  pow(self.zbar, -1)
        nx = len(self.electron_temperature)
        HeatFlowE = np.zeros(nx + 1)

        for i in range(1, nx):
            centered_ke = 2*(kappaE[i]*kappaE[i-1])/(kappaE[i]+kappaE[i-1])
            #centered_ke = 0.5 * (kappaE[i] + kappaE[i - 1])
            HeatFlowE[i] = centered_ke *((self.electron_temperature[i] - self.electron_temperature[i - 1]) 
                            / (self.cell_centered_coord[i] - self.cell_centered_coord[i - 1]))
            
        HeatFlowE[0] = 0
        HeatFlowE[-1] = 0
        self.spitzer_harm_heat = -1 * HeatFlowE

    def divQHeatFlow(self, step = 2,  laser_dir = None):
        """ Purpose: Find Div.Q
            Args:
                electron_thermal_flux = SOL-KiT heat flux in SI
                mass = areal mass.
                step = A singel step represent index to next cell-wall. 
                in the case of SOL-KiT where quantites are defiend in 
                cell-centres and walls, require a step of two. In
                HyKiCT quantites are only defined at one point, thus,
                a step of 1 is sufficient. 
                
            Returns: Div.Q
            NOTE: ONLY WORKS FOR SAME GRIDS
        """
        ##Logic relevant only for Limit Density methods 
        heat_flow = self.vfp_heat

        #Start limit-density
        limit_index = len(heat_flow)
        if laser_dir is not None:
            if laser_dir == "right":
                len_sh = len(self.spitzer_harm_heat)
                append_index = len_sh - limit_index + 1
                heat_flow = heat_flow[1:]
                heat_flow =  np.append(self.spitzer_harm_heat[:append_index], heat_flow)
            else:
                heat_flow = heat_flow[:-1]
                heat_flow = np.append(heat_flow, self.spitzer_harm_heat[limit_index - 1:])
        #End

        nx = len(heat_flow) 
        HeatConductionE = np.zeros(nx - 1)
        for i, m in enumerate(self.mass):
            HeatConductionE[i] = (-(heat_flow[i + 1] - heat_flow[i])) / m #- sign is there because of convention used in HyKiCT 
        return(HeatConductionE)
    
    def getSubtractDq(self, laser_dir = None):
        
        ##Logic relevant only for Limit Density methods 
        heat_flow = self.vfp_heat

        #Start limit-density
        limit_index = len(heat_flow)
        if laser_dir is not None:
            len_sh = len(self.spitzer_harm_heat)
            if laser_dir == "right":
                #Remove first term as its 0 
                heat_flow = heat_flow[1:]
                append_index = len_sh - limit_index + 1
                #To retain spitzer-heat in fluid code remove the subtraction factor. 
                self.spitzer_harm_heat[:append_index] = 0 
                #Pad the VFP heat-flow as it wont be same shape as spitzer-harm
                padding = len_sh - limit_index + 1
                heat_flow = np.pad(heat_flow, (padding,0), 'constant', constant_values = (0,0))
            else:
                #Remove last term as its 0 
                heat_flow = heat_flow[:-1]
                #To retain spitzer-heat in fluid code remove the subtraction factor. 
                self.spitzer_harm_heat[limit_index - 1:] = 0 
                #Pad the VFP heat-flow as it wont be same shape as spitzer-harm
                padding = len_sh - limit_index + 1
                heat_flow = np.pad(heat_flow, (0,padding), 'constant', constant_values = (0,0))
        #End

        if self.snb:
            return heat_flow - self.q_snb
        else:
            return heat_flow - self.spitzer_harm_heat

    def preHeatModel(self, rough = False):
        """ 
        Purpose: Model Pre-heat using an exponential fitting parameter, fit parameter
                    is spat out to be used by the fluid code.
        Args:
        Returns:
            B = Fitting paramters
            preheat_start = start index for preheat
            preheat_end = end index of preheat
        """
        if rough:
            all_index = np.where(np.isinf(self.q_vfp_q_sh_multipliers) == True)[0]
            valid_index = [i for i in all_index if i > int(len(self.vfp_heat)/2)]
            if len(valid_index) == 0:
                return self.pre_heat_start_index ,self.pre_heat_last_index , self.pre_heat_fit_params 
            preheat_start = valid_index[0] - 1
            preheat_end = np.where(abs(self.vfp_heat[preheat_start:]) < self.search_tolerance)
        else:
            preheat_start = np.where(self.q_vfp_q_sh_multipliers[~np.isinf(self.q_vfp_q_sh_multipliers)] != 0)[0][-1]
            preheat_end = np.where(abs(self.vfp_heat[preheat_start:]) < self.search_tolerance)
        #if search fails i.e. valid heat flow in all domain
        if(len(preheat_end[0]) == 0):
            preheat_end = len(self.q_vfp_q_sh_multipliers)
        else:
            preheat_end = preheat_end[0][0] + preheat_start

        L = self.cell_wall_coord[preheat_end] -self.cell_wall_coord[preheat_start] 
        B = []

        for i in range(preheat_start, preheat_end):
            b = (-1* (self.cell_wall_coord[i] - self.cell_wall_coord[preheat_start]) 
                / (L * np.log(abs(self.vfp_heat[i]/self.vfp_heat[preheat_start]))))
            B.append(b)
            
        B = np.array(B)
        return(preheat_start, preheat_end, B)

    def frontHeatModel(self, rough=False):
        """
        Purpose: Model heat wave from top of a bath.
        Args:
        Returns:
            B = Fitting paramters
            frontheat_start = start of front heat
            frontheat_end = end of front 
            heat
        """
        if rough:
            all_index = np.where(np.isinf(self.q_vfp_q_sh_multipliers) == True)[0]
            valid_index = [i for i in all_index if i <= int(len(self.vfp_heat)/2)]
            if len(valid_index) == 0:
                return self.front_heat_start_index ,self.front_heat_last_index , self.front_heat_fit_params
            frontheat_start = valid_index[-1] + 1
            frontheat_end = np.where(abs(self.vfp_heat[:frontheat_start]) < self.search_tolerance)
        else:
            frontheat_start = np.where(self.q_vfp_q_sh_multipliers[~np.isinf(self.q_vfp_q_sh_multipliers)] != 0)[0][0]
            frontheat_end = np.where(abs(self.vfp_heat[:frontheat_start]) < self.search_tolerance)
        #if search fails i.e. valid heat flow in all domain
        if(len(frontheat_end[0]) == 0):
            frontheat_end = np.int64(0)
        else:
            frontheat_end = frontheat_end[0][0]

        L = abs(self.cell_wall_coord[frontheat_end] - self.cell_wall_coord[frontheat_start])
        B = []

        for i in range(0, frontheat_start+1):
            b = (-1* (self.cell_wall_coord[i] - self.cell_wall_coord[frontheat_start]) 
                / (L * np.log(abs(self.vfp_heat[i]/self.vfp_heat[frontheat_start]))))
            B.append(b)

        B = np.array(B)
        return(frontheat_start, frontheat_end, B)

    def _detectAnamalousHeat(self, heat_flow):
        if (all(self.spitzer_harm_heat == 0)):
            return None, None

        start_of_spitzer_harm_heat_flow_index = np.where(abs(heat_flow) > 0)[0][0]
        last_of_spitzer_harm_heat_flow_index = np.where(abs(heat_flow) > 0)[0][-1]
        front = None 
        pre = None
        inf_mulitplier_index = np.where(np.isinf(self.q_vfp_q_sh_multipliers) == True)[0]
        diff_inf_index = np.diff(inf_mulitplier_index)
        if any(diff_inf_index > 1):
            max_diff = np.argmax(diff_inf_index)
            if((inf_mulitplier_index[max_diff] < int(len(self.q_vfp_q_sh_multipliers)/2)) and 
                (inf_mulitplier_index[max_diff + 1] > int(len(self.q_vfp_q_sh_multipliers)/2))):
                pass
            else:
                start_of_spitzer_harm_heat_flow_index = int(len(heat_flow)/2)
                last_of_spitzer_harm_heat_flow_index = int(len(heat_flow)/2)
        if(any(np.isnan(self.q_vfp_q_sh_multipliers[1:-2])) #boundaries are always nan as we have 0 inflow conditions. 
            or any(np.isinf(self.q_vfp_q_sh_multipliers))):
        
            if any(abs(self.vfp_heat[:start_of_spitzer_harm_heat_flow_index]) > self.search_tolerance):
                if(start_of_spitzer_harm_heat_flow_index == int(len(heat_flow)/2)):
                    front = self.frontHeatModel(rough = True)
                else:
                    front = self.frontHeatModel()
            if(any(abs(self.vfp_heat[last_of_spitzer_harm_heat_flow_index:]) > self.search_tolerance)):
                if(start_of_spitzer_harm_heat_flow_index == int(len(heat_flow)/2)):
                    pre = self.preHeatModel(rough = True)
                else:
                    pre = self.preHeatModel()

            self.q_vfp_q_sh_multipliers[np.isnan(self.q_vfp_q_sh_multipliers)] = 0
            self.q_vfp_q_sh_multipliers[np.isinf(self.q_vfp_q_sh_multipliers)] = 0
        return front, pre

    def multiplier(self, laser_dir = None):
        """ Purpose: Find multipliers and exponentially extrapolate for pre-heat
            Returns: Multipliers
        """
        ##Limit density method 
        ##Multipliers -> 1 for region that has been cut off 
        ## If there is heat-flow there going to be using spitzer-harm regardless 
        ## Models currently will work as intended. 
        ##Maybe rework or work around anamalous heat detection function. 
        self.front_heat_start_index = np.int64(0) 
        self.front_heat_last_index = np.int64(0)
        self.pre_heat_start_index = np.int64(0)
        self.pre_heat_last_index = np.int64(0)
        self.front_heat_fit_params = None
        self.pre_heat_fit_params = None

        if self.snb:
            heat_flow = None #self.q_snb 
        else:
            heat_flow = self.spitzer_harm_heat

        limit_index = len(heat_flow)
        if laser_dir is not None:
            len_vfp = len(self.vfp_heat)
            if laser_dir == "right":
                append_index = limit_index - len_vfp# + 1
                heat_flow = heat_flow[append_index:]
            else:
                heat_flow = heat_flow[:len_vfp] 
        #End
        ##Test for pre-heat via looking at NaN outputs expected from q/q_sh
        #if nans
        if self.snb:
            self.q_vfp_q_sh_multipliers = np.array(self.vfp_heat/heat_flow)
        else:
            self.q_vfp_q_sh_multipliers = np.array(self.vfp_heat/heat_flow)


        if laser_dir is not None:
            if laser_dir == "right":
                self.q_vfp_q_sh_multipliers[-1] = 0
            else:
                self.q_vfp_q_sh_multipliers[0] = 0
            
        #Detect if there is Front heat
        #Detect if there is Pre-Heat
        #Modify as bounds will always be nan. 
        if not self.snb:
            front, pre = self._detectAnamalousHeat(heat_flow)
            if front is not None:
                (self.front_heat_start_index, 
                self.front_heat_last_index, 
                self.front_heat_fit_params) = front
            if pre is not None:
                (self.pre_heat_start_index,
                self.pre_heat_last_index,
                self.pre_heat_fit_params) = pre

        if laser_dir is not None:
            padding = limit_index - len_vfp
            if laser_dir == "right":
                self.q_vfp_q_sh_multipliers[0] = 1
                self.q_vfp_q_sh_multipliers = np.pad(self.q_vfp_q_sh_multipliers, (padding, 0), 'constant', constant_values = (1,0))
            else:
                self.q_vfp_q_sh_multipliers[-1] = 1
                self.q_vfp_q_sh_multipliers = np.pad(self.q_vfp_q_sh_multipliers, (0, padding), 'constant', constant_values = (0,1))
            
            if self.front_heat_start_index > 0:
                self.front_heat_start_index += padding
                self.front_heat_last_index += padding
            if self.pre_heat_start_index > 0:
                self.pre_heat_last_index += padding
                self.pre_heat_last_index += padding

        self.q_vfp_q_sh_multipliers[0] = 0
        self.q_vfp_q_sh_multipliers[-1] = 0
        return(self.q_vfp_q_sh_multipliers, self.pre_heat_start_index, 
                self.pre_heat_last_index, self.pre_heat_fit_params, 
                self.front_heat_start_index, self.front_heat_last_index, self.front_heat_fit_params)

    def interpolateThermo(self, boundary_condition, normalised_values):
        # NOrmalise SI to Impact norms
        sol_kit_x_grid = self.cell_wall_coord# / normalised_values["lambda_mfp"]
        sol_kit_x_centered_grid = self.cell_centered_coord# / normalised_values["lambda_mfp"]
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
        sol_kit_ne = self.electron_number_density #/ (normalised_values['ne'])
        sol_kit_te = self.electron_temperature #* (kb/e)) / normalised_values['Te']
        sol_kit_z = self.zbar
        
        #Require interpolation to get the centre quanties in SOL-KiT this is done via linear interpolations 
        #here we use cubic spline to smooth quanties. 
        self.inter_ne = np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_ne)
        self.inter_te = np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_te)
        self.inter_z =  np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_z)
        #Extend grid to cell-wall via interpolating using periodic bounds i.e. 
        #end cell-centre and first cell centre, distance from cell wall assumed to be the same
        #thus == 0.5
        if boundary_condition == 'periodic':
            self.inter_te[0] = 0.5 * (self.inter_te[1] + self.inter_te[-1])
            self.inter_ne[0] = 0.5 * (self.inter_ne[1] + self.inter_ne[-1])
            self.inter_z[0] = 0.5 * (self.inter_z[1] + self.inter_z[-1])
        return(sol_kit_grid)

    def electric_field(self, dx_wall):
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
        n = len(self.zbar)
        self.E = np.zeros(n - 1)
        for i in range(1, n):
            dx = dx_wall[i - 1]
            gamma =  1 + ((3 * (self.zbar[i] + 0.477)) / (2*(self.zbar[i] + 2.15)))
            self.E[i - 1] = (-1 * (kb/e) * self.inter_te[i*2] * 
                    (((self.inter_ne[i*2] - self.inter_ne[i*2 - 2]) / (self.inter_ne[2*i - 1] * dx)) +
                    gamma * ((self.inter_te[i*2] - self.inter_te[i*2 - 2]) /(self.inter_te[2*i - 1] * dx))))
        return(self.E)


    def lambda_E(self, lambda_star, E, e_grid):
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

    def getLambStar(self, beta):
        nx,ng = np.shape(beta)
        vt_bx = np.sqrt(kb*self.inter_te/me)
        log_lam_ei = self.lambda_ei(T_norm = self.inter_te * (kb/e), n_norm = self.inter_ne, Z_norm = self.inter_z, return_array = True)
        Y = 4*np.pi*(e**2/(4 * np.pi * epsilon_0 * me))**2
        tau_ei_bx = vt_bx**3/(Y* self.inter_z * self.inter_ne * log_lam_ei)

        lambda_ei_bx = vt_bx*tau_ei_bx
        zeta = (self.inter_z + 0.24) / (self.inter_z + 4.2) # Brodrick
        # lambda_e_bx = np.sqrt(Z)*lambda_ei_bx
        lambda_e_bx = zeta * lambda_ei_bx
        lambda_g_bx = np.zeros((nx,ng))
        lambda_g_c =  np.zeros((nx,ng))

        for i in range(ng):
            lambda_g_bx[:,i] = 2*pow(beta[:,i], 2)*lambda_e_bx #Defined everywhere in space,  energy centered 
        return lambda_g_bx

    def createMatrix(self, A, B, C):
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
    def solve(self, A, b):
        """
        Purpose : Solves A*x = b 
        Args:
            A = Matrix (nv) 
            b = source term
        Returns: x 
        """ 
        x = linalg.solve(A, b)
        return(x)

    def getBeta(self, energy_grid, Te):
        beta = np.zeros((len(Te),len(energy_grid)))
        dbeta = np.zeros((len(Te),len(energy_grid) - 1)) 
        for i in range(len(Te)):
            beta[i, :] = energy_grid / Te[i]
            for j in range(len(energy_grid) - 1):
                dbeta[i, j] = (energy_grid[j + 1] - energy_grid[j]) /  Te[i]
        return beta, dbeta
            
    def getSource(self, q_sh, beta, dbeta):
        nx, ne = np.shape(beta)
        U_g = np.zeros((nx, ne))
        for i in range(nx):
            for e in range(ne):
                U_g[i, e] = (q_sh[i] / 24) * dbeta[i, e] * pow(beta[i, e], 4) * np.exp(-1*beta[i, e])
        return U_g
        
    def getH(self, dx_wall, dx_centered, lambda_star, lambda_E, U, r, Z):
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

            mat = self.createMatrix(A, B, C)
            d = self.solve(mat, S)
            H[:, j] = d
        return(H)

    def getDqnl(self, H, lambda_E, dx_wall):
        nx, nv = np.shape(H)
        q_nl = np.zeros(nx - 1)
        for i in range(0, nx - 1):
            heat_flow = 0        
            grad_delta_f0 = ((H[i + 1, :] - H[i, :]) / dx_wall[i])
            for v in range(nv):
                heat_flow += lambda_E[i, v] * grad_delta_f0[v] 
            q_nl[i] = (1/3) * heat_flow
        return q_nl 

    def snb_heat_flow(self, ng, max_E, r):
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
        # self.x_centered = np.array([(x[i+1] + x[i])/ 2 for i in range(len(x) - 1)])
        full_x = self.interpolateThermo('reflective', None) 
        dx_centered = np.diff(full_x[1::2])
        dx_wall = np.diff(full_x[::2])
        dx_centered = np.insert(dx_centered, 0 , 2 * (full_x[1] - full_x[0]))
        dx_centered = np.append(dx_centered,  2 * (full_x[-1] - full_x[-2]))

        E_max = np.max(self.electron_temperature) * max_E
        e_grid_wall = np.linspace(0, E_max, ng + 1)
        e_grid_centre = np.array([(e_grid_wall[i+1] + e_grid_wall[i])/ 2 for i in range(len(e_grid_wall) - 1)])

        beta, _ = self.getBeta(e_grid_centre, self.inter_te)
        _, dbeta = self.getBeta(e_grid_wall, self.inter_te)

        E_field = self.electric_field(dx_wall)
        lambda_star = self.getLambStar(beta)
        lamb_E = self.lambda_E(lambda_star[1::2], E_field, e_grid_centre)
        
        # q_sh = spitzer_harm_heat(x_centered, Te, ne, Z) # local 
        U = self.getSource(self.spitzer_harm_heat[1:-1], beta[1::2, :], dbeta[1::2, :])
        H = self.getH(dx_wall, dx_centered, lambda_star[::2], lamb_E, U, r, self.zbar)
        q_nl = self.getDqnl(H, lamb_E, dx_wall)
        
        self.q_snb = self.spitzer_harm_heat[1:-1]  - q_nl
        self.q_snb = np.append(self.q_snb, 0)
        self.q_snb = np.insert(self.q_snb, 0 ,0)
