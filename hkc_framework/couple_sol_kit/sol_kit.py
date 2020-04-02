"""
Self-consistent class for SOL-KiT
This object will init Sol-KiT run and take and produce the necessary files to run from and to Fluid.  
@author = Abetharan Antony
Date = 25/11/19
"""
import os
import sys
import math
import pprint
import string
import shutil 
import fnmatch
import subprocess
import numpy as np 
from scipy import constants
from Templating import Templating
from distutils.dir_util import copy_tree
import hkc_framework.common.tmpfilecreator as tfc
from hkc_framework.common.kinetic import Kinetic
import hkc_framework.couple_sol_kit.set_sol_kit_param as setparams 

kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
epsilon_0 = 8.854188E-12    # Vacuum dielectric constant
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")

class SOL_KIT(Kinetic):
    
    def __init__(self,IO,CoupleDivQ_, CoupleMulti_, MaintainF0_):

        self._Kinetic = Kinetic()
        self.templater = Templating() 
        self._kinetic_io_obj = IO
    
        self.copy_sol_kit = get_sol_kit_
        self._run_name = IO._RUN_NAME
        self._run_path = IO._RUN_PATH
        self._SOL_KIT_INPUT_PATH = os.path.join(self._run_path, 'INPUT')
        self._SOL_KIT_OUTPUT_PATH = os.path.join(self._run_path, 'OUTPUT')
        self._base_dir = IO._BASE_DIR
        self._src_dir = IO._K_SRC_DIR
        self._fluid_output_path = IO.fluid_output_path
        self._cycle_path = IO.cycle_dump_path
        
        self.boundary_condition = None
        self.maintain_f0 = MaintainF0_
        self.normalised_values = None
        self.CX1 = CX1_ 

        self.CoupleDivQ = CoupleDivQ_
        self.CoupleMulti = CoupleMulti_
        self.sh_heat_flow = 0
        self.mnultipliers = 0
        self.pre_heat_start_index = 0 
        self.pre_heat_last_index = 0 
        self.front_heat_start_index = 0 
        self.front_heat_last_index = 0 

        self.makeTmpFiles()
        self.copySOL_KiT()
        self.normalisation()


    def SOL_KITRun(self):
        """  
        Purpose: Launch the SOL-KiT
        Funcs : Contains the convergance test function which calculates q/q_sh and its results are passed
                into memory of daemon. 
        Args:
        """
        
        def convergance_test(path):
            
            """  
            Purpose: Calculate q/q_sh given a path
            Args: 
                path = path to latest output heat flow
            returns:
                multipliers = q/q_sh for latest heat flow
            """
            
            #Calculate Spitzer-Harm heat flow first time around
            if self.sh_heat_flow == 0:
                largest_fluid_index = self._findLargestIndex(os.path.join(self._fluid_output_path, "ELECTRON_TEMPERATURE"))
                temperature = np.loadtxt(os.path.join(self._kinetic_io_obj.fluid_output_path + "ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_" + str(largest_fluid_index) +".txt"))
                number_density = np.loadtxt(os.path.join(self._kinetic_io_obj.fluid_output_path + "ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_" + str(largest_fluid_index) +".txt"))
                zbar = np.loadtxt(os.path.join(self._kinetic_io_obj.fluid_output_path + "ZBAR/ZBAR_" + str(largest_fluid_index) + ".txt"))
                cell_wall_x = np.loadtxt(os.path.join(self._kinetic_io_obj.fluid_output_path + "CELL_WALL_X/CELL_WALL_X_" + str(largest_fluid_index) + ".txt")) 
                self.sh_heat_flow = self.SpitzerHarmHeatFlow(np.array(cell_wall_x), temperature,zbar,number_density)
            
            kinetic_heat_profile = np.loadtxt(path)
            #Include first and celll wall for no flow set to 0 per no thermal influx BC = 0 
            if self.boundary_condition =='noflow':
                kinetic_heat_profile = np.insert(kinetic_heat_profile, 0, 0.)
            #Periodic boundary condition removes last cell wall insert back and set to 0 via thermal influc BC
            kinetic_heat_profile = np.append(kinetic_heat_profile, 0.) * self.normalised_values['qe']
            self.multipliers = kinetic_heat_profile[::2] / self.sh_heat_flow
    
            return self.multipliers

        os.chdir(self._run_path)
        
        if self.CX1:
            cmd = ["mpiexec", "./SOL-KiT"]    
        else:
            cmd = ["mpirun", "-np", str(self._np), "./SOL-KiT"]

        heat_flow_path = os.path.join(self._run_path, 'OUTPUT/HEAT_FLOW_X/')
        super().Execute(cmd, self._cycle_path, self.maintain_f0, heat_flow_path, convergance_test, self._nx)
    

    def normalisation(self):
        """ Purpose: Calculates SOL-KiT Normalisation
            
            Params: self.norm_{} {ne, Te, Z}
            
            OUT: Saves norm dict self.normlisation_valeues contains (ne, Te,Z, vte, qe, lambda_mfp, tau_ei)
            
            TAKEN FROM: SOL_KiT_tools.py
        """
        #9.10938E-31            # Electron mass in kg
        el_charge = e #1.602189E-19    # Electron charge
        r_b = bohr_radi #5.29E-11              # Bohr radius in metres

        sigma_0 = math.pi * r_b ** 2  # Cross-section normaliation
        gamma_ee_0 = el_charge ** 4 / (4 * math.pi * (me * epsilon_0) ** 2)

        T_J = self._norm_Te * el_charge

        gamma_ei_0 = self._norm_Z ** 2 * gamma_ee_0

        v_t = math.sqrt(2.0 * T_J / me)                  # Velocity normalization
        coulomb_log = self.lambda_ei(1.0, 1.0, T_norm = self._norm_Te, n_norm = self._norm_ne, Z_norm = self._norm_Z)
        time_norm = v_t ** 3 / (gamma_ei_0 * (self._norm_ne/self._norm_Z) * coulomb_log)  # Time normalization
        x_0 = v_t * time_norm               # Space normalization
        e_0 = me * v_t / (el_charge * time_norm) # Electric field normalization
        q_0 = me * self._norm_ne * (v_t ** 3)
    
        dict = {}

        dict['ne'] = self._norm_ne
        dict['Te'] = self._norm_Te
        dict['Z'] = self._norm_Z

        # convert ne to 10**21 cm**-3
        ni = self._norm_ne / self._norm_Z
        dict['ni'] = ni
        dict['vth'] = v_t
        dict['tau_ei'] = time_norm
        dict['nu_ei'] = 1/time_norm
        dict['lambda_mfp'] = x_0
        dict['coulomb_log'] = coulomb_log
        dict['qe'] = q_0
        self.normalised_values = dict

    def makeTmpFiles(self):
        """ Purpose: Create all temporary files defined in tmpFileCreator.py
        """
        tfc.SOL_KIT_GRID_INPUT(self._run_path)
        tfc.SOL_KIT_NORMS(self._run_path)
        tfc.SOL_KIT_SOLVER_PARAMS(self._run_path)
        tfc.SOL_KIT_SWITCHES(self._run_path)

    def copySOL_KiT(self):
        """ Purpose: Copy SOL-KiT and creat necessary file structure for SOL-KiT to run i.e. INPUT/OUTPUT 
            Args: self.run_path 
            LOGIC: TO copy sol-kit or not.
        """
        if self.copy_sol_kit:
            shutil.copy(self._src_dir + '/SOL-KiT', self._run_path)
        
        OUTPUT_path = os.path.join(self._run_path, "OUTPUT")
        if not os.path.exists(OUTPUT_path):
            os.mkdir(OUTPUT_path)
            os.mkdir(os.path.join(OUTPUT_path, 'ATOMIC_EN_TEST'))
            os.mkdir(os.path.join(OUTPUT_path, 'CURRENT_TEST'))
            os.mkdir(os.path.join(OUTPUT_path, 'DEEX_E_RATE'))
            os.mkdir(os.path.join(OUTPUT_path, 'DENSITY'))
            os.mkdir(os.path.join(OUTPUT_path, 'DIST_F'))
            os.mkdir(os.path.join(OUTPUT_path, 'E_FIELD_X'))
            os.mkdir(os.path.join(OUTPUT_path, 'EX_E_RATE'))
            os.mkdir(os.path.join(OUTPUT_path, 'FLOW_VEL_X'))
            os.mkdir(os.path.join(OUTPUT_path, 'GRIDS'))
            os.mkdir(os.path.join(OUTPUT_path, 'HEAT_FLOW_X'))
            os.mkdir(os.path.join(OUTPUT_path, 'HEAT_FLOW_X_TOT'))
            os.mkdir(os.path.join(OUTPUT_path, 'ION_DENS'))
            os.mkdir(os.path.join(OUTPUT_path, 'ION_VEL'))
            os.mkdir(os.path.join(OUTPUT_path, 'NEUTRAL_DENS'))
            os.mkdir(os.path.join(OUTPUT_path, 'NUM_DV_HEATING'))
            os.mkdir(os.path.join(OUTPUT_path, 'QN_TEST'))
            os.mkdir(os.path.join(OUTPUT_path, 'RAD_DEEX_E_RATE'))
            os.mkdir(os.path.join(OUTPUT_path, 'RAD_REC_E_RATE'))
            os.mkdir(os.path.join(OUTPUT_path, 'REC_3B_E_RATE'))
            os.mkdir(os.path.join(OUTPUT_path, 'S_ION_M'))
            os.mkdir(os.path.join(OUTPUT_path, 'SH_q_ratio'))
            os.mkdir(os.path.join(OUTPUT_path, 'SHEAT_DATA'))
            os.mkdir(os.path.join(OUTPUT_path, 'TEMPERATURE'))
            os.mkdir(os.path.join(OUTPUT_path, 'TIMESTEP_DATA'))
            os.mkdir(os.path.join(OUTPUT_path, 'TOT_DENS_DATA'))
            #copy_tree(self._src_dir + '/OUTPUT/', OUTPUT_path)


    def setSOL_KITFiles(self, switch_boundary_condition_, switch_EE_L_ = True, switch_DIAG_EE_ = False, switch_IMPACT_MODE_ = False,
                         switch_COLD_IONS_ = False):
        """ Purpose: Set all neccessary SOL-KiT Files: Switches, norms, tolerances, grid_input
            Args: 
                switch_boundary_condition = string containing periodic, fixed or no flow
                switch_EE_L = arbitary L e-e collision operator
                switch_DIAG_EE = tridiagonal form of ee ... faster
                switch_IMPACT_MODE = start with local f_1
                switch_COLD_IONS = include cold ion motion
            Returns : Nothing ... files are set in INPUT fodler path.            
        """
        INPUT_PATH = os.path.join(self._run_path, "INPUT")
        if not os.path.exists(INPUT_PATH):
            os.makedirs(INPUT_PATH)

        #shutil.move(self._run_path + '/tmpSOL_KIT_SOLVER_PARAMS.txt', INPUT_PATH + '/SOLVER_PARAMS_INPUT.txt')
        shutil.copy(self._src_dir + "/INPUT/SOLVER_PARAMS_INPUT.txt", INPUT_PATH)
        shutil.copy(self._src_dir + "/INPUT/F_INIT_INPUT.txt", INPUT_PATH)
        grid = setParams.GRID(nt = self._nt, prestepnt = self._pre_step_nt, dt = self._dt,
                            predt = self._pre_step_dt, save_freq = self._save_freq, lmax = self._lmax,
                            nv = self._nv, dv = self._dv, v_multi = self._v_multi, nx = self._nx,
                            dx = self._dx, smalldx=0.0) 
        
        norms = setParams.NORMS(Z = self._norm_Z, Ar = self._norm_Ar, Density = np.float(self._norm_ne), 
                                Temp = self._norm_Te)

        self.boundary_condition = switch_boundary_condition_
        switch_fixeddown = False
        switch_fixedup = False
        switch_noflowdown = False
        switch_noflowup = False
        switch_periodic = False
        
        boundaries = ['periodic', 'fixed', 'noflow']
        if 'periodic' == self.boundary_condition:
            switch_periodic = True
        if 'fixed' == self.boundary_condition:
            switch_fixeddown = True
            switch_fixedup = True
        if 'noflow' == self.boundary_condition:
            switch_noflowdown = True
            switch_noflowup = True

        switches = setParams.Switches(EE_L = switch_EE_L_, DIAG_EE = switch_DIAG_EE_, 
                                    PERIODIC = switch_periodic, FIXEDUP = switch_fixedup, 
                                    FIXEDDOWN = switch_fixeddown, NOFLOWUP = switch_noflowup, 
                                    NOFLOWDOWN = switch_noflowdown, IMPACTMODE = switch_IMPACT_MODE_, 
                                    COLDIONS = switch_COLD_IONS_, MAINTAIN = self.maintain_f0)
                                
        self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_GRID.txt',
            writePath=INPUT_PATH, fileName="GRID_INPUT.txt", parameters=grid)
        self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_NORMS.txt',
            writePath=INPUT_PATH, fileName="NORMALIZATION_INPUT.txt", parameters=norms)
        #self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_SOLVER_PARAMS.txt',
        #writePath=INPUT_PATH, fileName="", parameters=)
        self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_SWITCHES.txt',
            writePath=INPUT_PATH, fileName="SWITCHES_INPUT.txt", parameters=switches)


    def _SOL_KITNormsToSI(self, var, var_name):
        """Purpose: Multiply SOL-KiT output with normlise to give back SI unit qunatity
            Args: var = quantity
                  var_name = variable to multiply with 
            Returns : converted values
        """
        return(var * self.normalised_values[var_name])

    def moveSOL_KITFiles(self):
        """Purpose: Move INPUT/OUTPUT into respective KINETIC_INPUT and KINETIC_OUTPUT folders"""
        import shutil
        #move INPUT files
        shutil.move(self._SOL_KIT_INPUT_PATH, self._kinetic_io_obj.kinetic_input_path)
        ###HEATING INPUT FOR FUTURE USE
        # shutil.move(os.path.join(self._SOL_KIT_INPUT_PATH, "NEUT_AND_HEAT_INPUT.txt"))    
        
        #move OUTPUT files
        shutil.move(self._SOL_KIT_OUTPUT_PATH, self._kinetic_io_obj.kinetic_output_path)
        
        
    def InitSOL_KITFromHydro(self):
        """ Purpose: Initilise all SOL-KiT load in files from HYDRO
        """
        last_index = self._findLargestIndex(os.path.join(self._fluid_output_path, "ELECTRON_TEMPERATURE"))
        f_x_grid = np.loadtxt(self._fluid_output_path + "CELL_WALL_X/CELL_WALL_X_" + str(last_index) + ".txt")
        f_x_centered_grid = np.loadtxt(self._fluid_output_path + "CELL_CENTRE_X/CELL_CENTRE_X_" + str(last_index) + ".txt")
        f_v = np.loadtxt(self._fluid_output_path + "VELOCITY/VELOCITY_" + str(last_index) + ".txt")
        f_ne = np.loadtxt(
                self._fluid_output_path + "ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_" + str(last_index) + ".txt")
        f_Te = np.loadtxt(
             self._fluid_output_path + "ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_" + str(last_index) + ".txt")
        f_laser = np.loadtxt(
            self._fluid_output_path + "INVERSE_BREM/INVERSE_BREM_" + str(last_index) + ".txt")
        f_brem = np.loadtxt(self._fluid_output_path + "BREM/BREM_" + str(last_index) + ".txt")
        f_Z = np.loadtxt(self._fluid_output_path + "ZBAR/ZBAR_" + str(last_index) + ".txt")

        # NOrmalise SI to Impact norms
        SOL_KIT_x_grid = f_x_grid / self.normalised_values["lambda_mfp"]
        SOL_KIT_x_centered_grid = f_x_centered_grid / self.normalised_values["lambda_mfp"]
        length_of_sol_kit_grid = len(SOL_KIT_x_centered_grid) + len(SOL_KIT_x_grid)
        full_x_grid = np.zeros(length_of_sol_kit_grid)
        j = 0
        k = 0
        for i in range(length_of_sol_kit_grid):
            if i % 2 != 0:
                full_x_grid[i] = SOL_KIT_x_centered_grid[j]
                j +=1
            else:
                full_x_grid[i] = SOL_KIT_x_grid[k]
                k +=1 
        
        #SOL-KiT periodic condition == .\.\.\ i.e. start with cell centre and end with cell wall
        #SOL-KiT noflow == .\.\. i.e. start and end with cell centres
        SOL_KIT_x_grid = None
        if self.boundary_condition == 'periodic':
            SOL_KIT_grid = full_x_grid[:-1]
        
        elif self.boundary_condition == "noflow":
            SOL_KIT_grid = full_x_grid[1:-1]

        #Conver to SOL-KiT units
        SOL_KIT_ne = f_ne / (self.normalised_values['ne'])
        SOL_KIT_Te = (f_Te * (kb/e)) / self.normalised_values['Te']
        #SOL_KIT_laser = (f_laser * f_density) / power_density
        #SOL_KIT_brem = (f_brem * f_density) / power_density
        SOL_KIT_Z = f_Z
        # SOL_KIT_HEATING = SOL_KIT_laser + SOL_KIT_brem
        
        #Require interpolation to get the centre quanties in SOL-KiT this is done via linear interpolations 
        #here we use cubic spline to smooth quanties. 
        from scipy.interpolate import CubicSpline, interp1d
        SOL_KIT_inter_ne = np.interp(SOL_KIT_grid, SOL_KIT_x_centered_grid, SOL_KIT_ne)
        SOL_KIT_inter_Te = np.interp(SOL_KIT_grid, SOL_KIT_x_centered_grid, SOL_KIT_Te)
        SOL_KIT_inter_Z =  np.interp(SOL_KIT_grid, SOL_KIT_x_centered_grid, SOL_KIT_Z)
        #SOL_KIT_inter_ne = spliner_ne(SOL_KIT_grid)
        #SOL_KIT_inter_Te = spliner_Te(SOL_KIT_grid)
        #SOL_KIT_inter_Z = spliner_Z(SOL_KIT_grid)

        np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "DENS_INPUT.txt"), SOL_KIT_inter_ne, fmt = '%.19f')    
        np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "TEMPERATURE_INPUT.txt"), SOL_KIT_inter_Te, fmt = '%.19f')    
        np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "Z_PROFILE_INPUT.txt"), SOL_KIT_inter_Z, fmt = '%.19f')    
        np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "X_GRID_INPUT.txt"), SOL_KIT_grid, fmt = '%.18f')
        # np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "ION_VEL_INPUT.txt"), SOL_KIT_grid)
        # np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "NEUT_HEAT_INPUT.txt"), SOL_KIT_heating)

    def _Multiplier(self, max_index):
        """ Purpose: Find multipliers and exponentially extrapolate for pre-heat
            Args: 
                max_index : Last index for output files.
            Returns: Multipliers
        """
        
        kinetic_heat_profile  = np.loadtxt(os.path.join(self._run_path, "OUTPUT/HEAT_FLOW_X/HEAT_FLOW_X_" + str(max_index).zfill(5) + '.txt'))
        cell_wall_x = np.loadtxt(os.path.join(self._kinetic_io_obj.fluid_output_path + "CELL_WALL_X/CELL_WALL_X_" + str(largest_fluid_index) + ".txt")) 

        #Include first and celll wall for no flow set to 0 per no thermal influx BC = 0 
        if self.boundary_condition =='noflow':
            kinetic_heat_profile = np.insert(kinetic_heat_profile, 0, 0.)

        #Periodic boundary condition removes last cell wall insert back and set to 0 via thermal influc BC
        kinetic_heat_profile = np.append(kinetic_heat_profile, 0.) * self.normalised_values['qe']
        
        ##Calculate spitzer-harm if not calculated by a backgroudn process.
        if self.sh_heat_flow == 0:
            largest_fluid_index = self._findLargestIndex(os.path.join(self._fluid_output_path, "ELECTRON_TEMPERATURE"))
            temperature = np.loadtxt(os.path.join(self._kinetic_io_obj.fluid_output_path + "ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_" + str(largest_fluid_index) +".txt"))
            number_density = np.loadtxt(os.path.join(self._kinetic_io_obj.fluid_output_path + "ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_" + str(largest_fluid_index) +".txt"))
            zbar = np.loadtxt(os.path.join(self._kinetic_io_obj.fluid_output_path + "ZBAR/ZBAR_" + str(largest_fluid_index) + ".txt"))
            self.sh_heat_flow = self.SpitzerHarmHeatFlow(np.array(cell_wall_x), temperature,
                                             zbar,number_density)

        #Get rid of cell-centre quantites
        qe_list = kinetic_heat_profile[0::2]
        ##Test for pre-heat via looking at NaN outputs expected from q/q_sh
        #if nans
        multi_list = qe_list/self.sh_heat_flow
        #Detect if there is Front heat
        #Detect if there is Pre-Heat
        #Modify as bounds will always be nan. 
        if(any(np.isnan(multi_list)) or any(np.isinf(multi_list))):

            multi_list = np.array(multi_list)
            qe_list = np.array(qe_list)
            NaN_args = np.argwhere(np.isnan(multi_list))
            pre_heat_start_index, pre_heat_last_index, pre_heat_fit_params = self.PreHeatModel(cell_wall_x, self.sh_heat_flow, qe_list)
            front_heat_start_index, front_heat_last_index, front_heat_fit_params = self.FrontHeatModel(cell_wall_x, self.sh_heat_flow, qe_list)
        
        else:
            pre_heat_start_index = 0
            pre_heat_last_index = 0
            front_heat_start_index = 0
            front_heat_last_index = 0
            pre_heat_fit_params = 0 
            front_heat_fit_params = 0 

        return(multi_list, pre_heat_start_index, pre_heat_last_index, pre_heat_fit_params, front_heat_start_index, front_heat_last_index, front_heat_fit_params)



    def SOL_KITInitNextHydroFiles(self):
        """ Purpose: Prepare HYDRO init files.. calculate heat flow qunatities either div.q or multiplier and copy
            last hydro step files.
        """
        qe_path =  os.path.join(self._run_path, "OUTPUT/HEAT_FLOW_X")
        max_index = self._findLargestIndex(qe_path)
        mass = np.loadtxt(self._kinetic_io_obj.fluid_input_path + "/mass.txt")        

        if self.CoupleDivQ:
            max_index_reformated = '{:>05d}'.format(max_index)
            norm_qe = self.normalised_values['vth']**3 * self.normalised_values['ne'] * me
            qe = np.loadtxt(os.path.join(qe_path, "HEAT_FLOW_X_" + max_index_reformated + '.txt')) * norm_qe
            electronheatflow= self._DivQHeatFlow(qe, mass)
        
            np.savetxt(os.path.join(self._kinetic_io_obj.next_fluid_input_path, "qe.txt"), electronheatflow)    
        
        elif self.CoupleMulti:
            multi_list, self.pre_heat_start_index, self.pre_heat_last_index, pre_heat_fit_params, self.front_heat_start_index, self.front_heat_last_index, front_heat_fit_params = self._Multiplier(max_index)
            np.savetxt(os.path.join(self._kinetic_io_obj.next_fluid_input_path, "qe.txt"), multi_list)
            np.savetxt(os.path.join(self._kinetic_io_obj.next_fluid_input_path, "pre_heat_fit_param.txt"), pre_heat_fit_params)
            np.savetxt(os.path.join(self._kinetic_io_obj.next_fluid_input_path, "front_heat_fit_param.txt"), front_heat_fit_params)


        #Coyp last step fluid quants
        largest_fluid_index = self._findLargestIndex(os.path.join(self._fluid_output_path, "ELECTRON_TEMPERATURE"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "CELL_WALL_X/CELL_WALL_X_" + str(largest_fluid_index) +".txt")    
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"coord.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "VELOCITY/VELOCITY_" + str(largest_fluid_index) +".txt") 
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"velocity.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "DENSITY/DENSITY_" + str(largest_fluid_index) +".txt")  
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"density.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_" + str(largest_fluid_index) +".txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"electron_temperature.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "ION_TEMPERATURE/ION_TEMPERATURE_" + str(largest_fluid_index) +".txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"ion_temperature.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "MASS/MASS_" + str(largest_fluid_index) +  ".txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"mass.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_input_path + "Z.txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path, "Z.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_input_path + "Ar.txt")  
                                    ,os.path.join(self._kinetic_io_obj.next_fluid_input_path, "Ar.txt"))





