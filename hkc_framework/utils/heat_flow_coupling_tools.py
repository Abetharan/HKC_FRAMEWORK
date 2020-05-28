import math
import numpy as np 
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

    def lambda_ei(self, T_norm , n_norm, Z_norm, return_arg = False):
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

        self.coulomb_log = np.array(coulomb_logs)

    def spitzerHarmHeatFlow(self):
        """
        Purpose: Models Spitzer Harm heat flow 
        """
        # self.coulomb_log = []
        # for i in range(len(zbar)):
        #     coulomb = self.lambda_ei(T_norm = self.electron_temperature[i] * (kb/e),
        #                              n_norm = self.electron_number_density[i], Z_norm = self.zbar[i])
            # self.coulomb_log.append(coulomb)
# 
        # self.coulomb_log = np.array(self.coulomb_log)
        kappaE = 1.843076667547614E-10 * pow(self.electron_temperature, 2.5) * pow(self.coulomb_log, -1) *  pow(self.zbar, -1)
        nx = len(self.electron_temperature)
        HeatFlowE = np.zeros(nx + 1)

        for i in range(1, nx):
            centered_ke = 0.5 * (kappaE[i] + kappaE[i - 1])
            HeatFlowE[i] = centered_ke *((self.electron_temperature[i] - self.electron_temperature[i - 1]) 
                            / (self.cell_centered_coord[i] - self.cell_centered_coord[i - 1]))
            
        HeatFlowE[0] = 0
        HeatFlowE[-1] = 0
        self.spitzer_harm_heat = -1 * HeatFlowE

    def preHeatModel(self):
        """ 
        Purpose: Model Pre-heat using an exponential fitting parameter, fit parameter
                    is spat out to be used by the fluid code.
        Args:
        Returns:
            B = Fitting paramters
            preheat_start = start index for preheat
            preheat_end = end index of preheat
        """
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
                / (L * np.log(abs(self.vfp_heat[i])/self.vfp_heat[preheat_start])))
            B.append(b)
            
        B = np.array(B)
        return(preheat_start, preheat_end, B)

    def frontHeatModel(self):
        """
        Purpose: Model heat wave from top of a bath.
        Args:
        Returns:
            B = Fitting paramters
            frontheat_start = start of front heat
            frontheat_end = end of front 
            heat
        """

        frontheat_start = np.where(self.q_vfp_q_sh_multipliers[~np.isinf(self.q_vfp_q_sh_multipliers)] != 0)[0][0]
        frontheat_end = np.where(abs(self.vfp_heat[:frontheat_start]) < self.search_tolerance)
        #if search fails i.e. valid heat flow in all domain
        if(len(frontheat_end[0]) == 0):
            frontheat_end = 0
        else:
            frontheat_end = frontheat_end[0][0]

        L = abs(self.cell_wall_coord[frontheat_end] - self.cell_wall_coord[frontheat_start])
        B = []

        for i in range(0, frontheat_start+1):
            b = (-1* (self.cell_wall_coord[i] - self.cell_wall_coord[frontheat_start]) 
                / (L * np.log(abs(self.vfp_heat[i])/self.vfp_heat[frontheat_start])))
            B.append(b)

        B = np.array(B)
        return(frontheat_start, frontheat_end, B)

    def divQHeatFlow(self, step = 2):
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
        nx = len(self.vfp_heat) 
        HeatConductionE = np.zeros(nx - 1)
        for i, m in enumerate(self.mass):
            HeatConductionE[i] = (-(self.vfp_heat[i + 1] - self.vfp_heat[i])) / m
        return(HeatConductionE)
    
    def _detectAnamalousHeat(self):
        if all(self.spitzer_harm_heat == 0):
            return None, None

        start_of_spitzer_harm_heat_flow_index = np.where(self.spitzer_harm_heat > 0)[0][0]
        last_of_spitzer_harm_heat_flow_index = np.where(self.spitzer_harm_heat > 0)[0][-1]
        front = None 
        pre = None
        if(any(np.isnan(self.q_vfp_q_sh_multipliers[1:-2])) #boundaries are always nan as we have 0 inflow conditions. 
            or any(np.isinf(self.q_vfp_q_sh_multipliers))):
            if any(abs(self.vfp_heat[:start_of_spitzer_harm_heat_flow_index]) > self.search_tolerance):
                front = self.frontHeatModel
            if(any(abs(self.vfp_heat[last_of_spitzer_harm_heat_flow_index:]) > self.search_tolerance)):
                pre = self.preHeatModel

            self.q_vfp_q_sh_multipliers[np.isnan(self.q_vfp_q_sh_multipliers)] = 0
            self.q_vfp_q_sh_multipliers[np.isinf(self.q_vfp_q_sh_multipliers)] = 0
        return front, pre

    def multiplier(self):
        """ Purpose: Find multipliers and exponentially extrapolate for pre-heat
            Args: 
                vfp_qe : Kinetic heat flow to be transfered to fluid code.
            Returns: Multipliers
        """
                
        front_heat_start_index = np.int64(0) 
        front_heat_last_index = np.int64(0)
        pre_heat_start_index = np.int64(0)
        pre_heat_last_index = np.int64(0)
        front_heat_fit_params = None
        pre_heat_fit_params = None
        ##Test for pre-heat via looking at NaN outputs expected from q/q_sh
        #if nans
        self.q_vfp_q_sh_multipliers = np.array(self.vfp_heat/self.spitzer_harm_heat)
        self.q_vfp_q_sh_multipliers[0] = 0
        self.q_vfp_q_sh_multipliers[-1] = 0
        #Detect if there is Front heat
        #Detect if there is Pre-Heat
        #Modify as bounds will always be nan. 
        front, pre = self._detectAnamalousHeat()
        if front is not None:
            (front_heat_start_index, 
            front_heat_last_index, 
            front_heat_fit_params) = front()
        if pre is not None:
            (pre_heat_start_index,
            pre_heat_last_index,
            pre_heat_fit_params) = pre()

        return(self.q_vfp_q_sh_multipliers, pre_heat_start_index, 
                pre_heat_last_index, pre_heat_fit_params, 
                front_heat_start_index, front_heat_last_index, front_heat_fit_params)
