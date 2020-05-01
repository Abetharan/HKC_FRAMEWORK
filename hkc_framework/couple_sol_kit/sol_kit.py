"""
Self-consistent class for SOL-KiT
This object will init Sol-KiT run and take and produce the necessary files to run from and to Fluid.  
@author = Abetharan Antony
"""
# from __future__ import absolute_import

from distutils.dir_util import copy_tree
import fnmatch
import math
import numpy as np 
import os
import pathlib
import pprint
from scipy import constants
from scipy.interpolate import interp1d
import shutil 
import string
import subprocess
import sys
from utils import templating
from utils import Input
from utils import Kinetic
from utils import findLargestIndex
from utils import HeatFlowCouplingTools
import couple_sol_kit.set_sol_kit_param as setparams
# from . import set_sol_kit_param as setparams
#import set_sol_kit_param as setparams

BOLTZMANN_CONSTANT = constants.value("Boltzmann constant")
ELECTRON_MASS = constants.value("electron mass")
PROTON_MASS = constants.value("proton mass")
ELEMENTARY_CHARGE = constants.value("elementary charge")
VACUUM_PERMITTIVITY = 8.854188E-12    # Vacuum dielectric constant
PLANCK_CONSTANT = constants.value("Planck constant")
BOHR_RADIUS = constants.value("Bohr radius")

class SOL_KIT(Kinetic):
    
    def __init__(self,run_path, k_src_dir, kinetic_input_path, kinetic_output_path,
                k_config_yml_file_path, convergence_monitoring = False, cx1 = False):
        
        # config_yml_file_path = os.path.join(
        #                         pathlib.Path(__file__).parent.absolute(),
                                # 'config.yml')
        self.tmp_input_path = os.path.join(
                                pathlib.Path(__file__).parent.absolute(),
                                'tmp_INPUT')
        self.tmp_output_path = os.path.join(
                                pathlib.Path(__file__).parent.absolute(),
                                'tmp_OUTPUT')

        #objects
        self.heat_flow_tools = HeatFlowCouplingTools() 
        self.init = Input(k_config_yml_file_path)
        
        #paths
        self._run_path = run_path
        self._kinetic_output_path = kinetic_output_path
        self._kinetic_input_path = kinetic_input_path
        self._kinetic_src_dir = k_src_dir
        self.cycle_dump_path = ""
        self.previous_cycle_output_path = ""
        self._sol_kit_input_path = os.path.join(self._run_path, 'INPUT')
        self._sol_kit_output_path = os.path.join(self._run_path, 'OUTPUT')
        
        self.maintain_f0 = self.init.yaml_file['Switches']['f_0_maintain']
        self.load_f1 = False
        self.convergence_monitoring = convergence_monitoring
        self.normalised_values = None
        self.cx1 = cx1
        self._np = self.init.yaml_file['Params']['Np'] 

        self._norm_Te = float(self.init.yaml_file['Norms']['Te'])
        self._norm_Z = float(self.init.yaml_file['Norms']['Z'])
        self._norm_Ar = float(self.init.yaml_file['Norms']['Ar'])
        self._norm_ne = float(self.init.yaml_file['Norms']['Ne'])

        self.copyAndCreateSOL_KiT()
        self.setFiles()
        self.normalisation()

        self.sh_heat_flow = 0

        if(self.init.yaml_file['Boundary']['Periodic']):
            self.boundary_condition = "periodic"
        elif self.init.yaml_file['Boundary']['Fixed_up']:
            self.boundary_condition = "fixed"
        else:
            self.boundary_condition = "noflow"

        if self.cx1:
            cmd = ["mpiexec", "./SOL-KiT"]    
        else:
            cmd = ["mpirun", "-np", str(self._np), "./SOL-KiT"]
        if(convergence_monitoring):
            Kinetic.__init__(self, cmd, convergence_monitoring, self.convergance_test,
                                 self._run_path)
        else:
            Kinetic.__init__(self, cmd)

    def convergance_test(self, path):
        
        """  
        Purpose: Calculate q/q_sh given a path
        Args: 
            path = path to latest output heat flow
        returns:
            multipliers = q/q_sh for latest heat flow
        """
        
        #Calculate Spitzer-Harm heat flow first time around
        kinetic_heat_profile = np.loadtxt(path)
        #Include first and celll wall for no flow set to 0 per no thermal influx BC = 0 
        if self.boundary_condition =='noflow':
            kinetic_heat_profile = np.insert(kinetic_heat_profile, 0, 0.)
        #Periodic boundary condition removes last cell wall insert back and set to 0 via thermal influc BC
        kinetic_heat_profile = np.append(kinetic_heat_profile, 0.) * self.normalised_values['qe']
        multipliers = kinetic_heat_profile[::2] / self.sh_heat_flow

        return multipliers

    def Run(self):
        """  
        Purpose: Launch the SOL-KiT
        Funcs : Contains the convergance test function which calculates q/q_sh and its results are passed
                into memory of daemon. 
        Args:
        """
        os.chdir(self._run_path)
        heat_flow_path = os.path.join(self._run_path, 'OUTPUT/HEAT_FLOW_X/')
        super().Execute(heat_flow_path)
    
    def getPhysicalRunTime(self):
        step = self.init.yaml_file['Params']['Nt'] 
        dt = self.init.yaml_file['Params']['Dt']
        return(dt * step * self.normalised_values['tau_ei']) 

    def normalisation(self):
        """ Purpose: Calculates SOL-KiT Normalisation
            
            Params: self.norm_{} {ne, Te, Z}
            
            OUT: Saves norm dict self.normlisation_valeues contains (ne, Te,Z, vte, qe, lambda_mfp, tau_ei)
            
            TAKEN FROM: SOL_KiT_tools.py
        """
        sigma_0 = math.pi * BOHR_RADIUS ** 2  # Cross-section normaliation
        gamma_ee_0 = ELEMENTARY_CHARGE ** 4 / (4 * math.pi * (ELECTRON_MASS * VACUUM_PERMITTIVITY) ** 2)

        T_J = self.init.yaml_file['Norms']['Te'] * ELEMENTARY_CHARGE

        gamma_ei_0 = self.init.yaml_file['Norms']['Z'] ** 2 * gamma_ee_0

        v_t = math.sqrt(2.0 * T_J / ELECTRON_MASS)                  # Velocity normalization
        coulomb_log = self.heat_flow_tools.lambda_ei(T_norm = [self._norm_Te],
                                                    n_norm = [self._norm_ne],
                                                    Z_norm = [self._norm_Z],
                                                    return_arg = True)

        time_norm = v_t ** 3 / (gamma_ei_0 * (self._norm_ne/self._norm_Z) * coulomb_log)  # Time normalization
        x_0 = v_t * time_norm               # Space normalization
        e_0 = ELECTRON_MASS * v_t / (ELEMENTARY_CHARGE * time_norm) # Electric field normalization
        q_0 = ELECTRON_MASS * self._norm_ne * (v_t ** 3)
    
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

    def copyAndCreateSOL_KiT(self):
        """ Purpose: Copy SOL-KiT and creat necessary file structure for SOL-KiT to run i.e. INPUT/OUTPUT 
            Args: self.run_path 
            LOGIC: TO copy sol-kit or not.
        """
        shutil.copy(self._kinetic_src_dir + '/SOL-KiT', self._run_path)

    def setFiles(self):
        """ Purpose: Set all neccessary SOL-KiT Files: Switches, norms, tolerances, grid_input
            Args: 
                switch_boundary_condition = string containing periodic, fixed or no flow
                switch_EE_L = arbitary L e-e collision operator
                switch_DIAG_EE = tridiagonal form of ee ... faster
                switch_IMPACT_MODE = start with local f_1
                switch_COLD_IONS = include cold ion motion
            Returns : Nothing ... files are set in INPUT fodler path.            
        """
        copy_tree(self.tmp_input_path, self._sol_kit_input_path)
        copy_tree(self.tmp_output_path, self._sol_kit_output_path)

        grid = setparams.GRID(
                            nt = self.init.yaml_file['Params']['Nt'],
                            prestepnt = self.init.yaml_file['Params']['Pre_step_nt'],
                            dt = self.init.yaml_file['Params']['Dt'],
                            predt = self.init.yaml_file['Params']['Pre_step_dt'],
                            save_freq = self.init.yaml_file['Params']['Output_freq'],
                            lmax = self.init.yaml_file['Params']['L_max'],
                            nv = self.init.yaml_file['Params']['Nv'],
                            dv = self.init.yaml_file['Params']['Dv'],
                            v_multi = self.init.yaml_file['Params']['Dv_multi'],
                            nx = self.init.yaml_file['Params']['Nx'],
                            dx = self.init.yaml_file['Params']['Dx'],
                            smalldx=0.0) 
        
        norms = setparams.NORMS(Z = self._norm_Z, Ar = self._norm_Ar, Density = np.float(self._norm_ne), 
                                Temp = self._norm_Te)

        self.switches = setparams.Switches(
                            EE_0 = self.init.yaml_file['Switches']['Coll_ee_0'] ,
                            EE_L = self.init.yaml_file['Switches']['Coll_ee_l'],
                            EI_L = self.init.yaml_file['Switches']['Coll_ei_l'],
                            DIAG_EE = self.init.yaml_file['Switches']['Diag_ee_l'],
                            PERIODIC = self.init.yaml_file['Boundary']['Periodic'],
                            FIXEDUP = self.init.yaml_file['Boundary']['Fixed_up'], 
                            FIXEDDOWN = self.init.yaml_file['Boundary']['Fixed_down'],
                            NOFLOWUP = self.init.yaml_file['Boundary']['No_flow_up'], 
                            NOFLOWDOWN = self.init.yaml_file['Boundary']['No_flow_down'],
                            IMPACTMODE = self.init.yaml_file['Switches']['Local_init'], 
                            COLDIONS = self.init.yaml_file['Switches']['Cold_ion_fluid'],
                            MAINTAIN = self.init.yaml_file['Switches']['f_0_maintain'],
                            OUTPUT_TEMPERATURE = self.init.yaml_file['Output']['Temperature'],
                            OUTPUT_DENSITY= self.init.yaml_file['Output']['Density'], 
                            OUTPUT_VELOCITY= self.init.yaml_file['Output']['Velocity'],
                            OUTPUT_E_FIELD= self.init.yaml_file['Output']['E_field'],
                            OUTPUT_SH_TEST= self.init.yaml_file['Output']['Sh_test'])
        
        templating(tmpfilePath= os.path.join(self._run_path, 'INPUT/tmpSOL_KIT_GRID.txt'),
            writePath=self._sol_kit_input_path, fileName="GRID_INPUT.txt", parameters=grid)
        templating(tmpfilePath= os.path.join(self._run_path, 'INPUT/tmpSOL_KIT_NORMS.txt'),
            writePath=self._sol_kit_input_path, fileName="NORMALIZATION_INPUT.txt", parameters=norms)
        templating(tmpfilePath= os.path.join(self._run_path, 'INPUT/tmpSOL_KIT_SWITCHES.txt'),
            writePath=self._sol_kit_input_path, fileName="SWITCHES_INPUT.txt", parameters=self.switches)
        #remove tmp files
        os.rename(os.path.join(self._sol_kit_input_path, 'tmpSOL_KIT_SOLVER_PARAMS.txt'), os.path.join(self._sol_kit_input_path, 'SOLVER_PARAMS_INPUT.txt'))
        os.remove(os.path.join(self._sol_kit_input_path, 'tmpSOL_KIT_GRID.txt'))
        os.remove(os.path.join(self._sol_kit_input_path, 'tmpSOL_KIT_NORMS.txt'))
        # os.remove(os.path.join(self._sol_kit_input_path, 'tmpSOL_KIT_SWITCHES.txt'))


    def _SOL_KITNormsToSI(self, var, var_name):
        """Purpose: Multiply SOL-KiT output with normlise to give back SI unit qunatity
            Args: var = quantity
                  var_name = variable to multiply with 
            Returns : converted values
        """
        return(var * self.normalised_values[var_name])

    def moveFiles(self):
        """Purpose: Move INPUT/OUTPUT into respective KINETIC_INPUT and KINETIC_OUTPUT folders"""
        #Copy INPUT files
        copy_tree(self._sol_kit_input_path, self._kinetic_input_path)
        # #move OUTPUT files
        # shutil.move(self._sol_kit_output_path, self._kinetic_output_path)

        #create new output folder
        #quicker to move and create new than copy and clearing 
        copy_tree(self._sol_kit_output_path, self._kinetic_output_path)
        with os.scandir(self._sol_kit_output_path) as output_folder:
            for var_folder in output_folder:
                if not var_folder.name.startswith('.') and var_folder.is_dir():
                    for var_txt in os.scandir(var_folder):
                        if not var_txt.name.startswith('.') and var_txt.is_file():
                            os.unlink(var_txt.path)
        
    def initFromHydro(self, f_x_grid, f_x_centered_grid, f_te, f_ne, f_z, f_laser = 0, f_rad = 0):
        """ Purpose: Initilise all SOL-KiT load in files from HYDRO
        Here the SOL-KiT specifics are done, SOL-KiT requires quantities to be defined at all 
        points on the grid. Whereas, the fluid quantities passed are not. 

        NOTE this is where specifics come into place and may need to be modified to
        the specific needs of the code being coupled. Here I assume commmon 
        Lagrangian Practice i.e. quantities are not defined at all points. 
        """

        # for line in proc.stdout:
        # NOrmalise SI to Impact norms
        sol_kit_x_grid = f_x_grid / self.normalised_values["lambda_mfp"]
        sol_kit_x_centered_grid = f_x_centered_grid / self.normalised_values["lambda_mfp"]
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
        if self.boundary_condition == 'periodic':
            sol_kit_grid = full_x_grid[:-1]
        
        elif self.boundary_condition == "noflow":
            sol_kit_grid = full_x_grid[1:-1]

        #Conver to SOL-KiT units
        sol_kit_ne = f_ne / (self.normalised_values['ne'])
        sol_kit_te = (f_te * (BOLTZMANN_CONSTANT/ELEMENTARY_CHARGE)) / self.normalised_values['Te']
        #SOL_KIT_laser = (f_laser * f_density) / power_density
        #SOL_KIT_brem = (f_brem * f_density) / power_density
        sol_kit_z = f_z
        # SOL_KIT_HEATING = SOL_KIT_laser + SOL_KIT_brem
        
        #Require interpolation to get the centre quanties in SOL-KiT this is done via linear interpolations 
        #here we use cubic spline to smooth quanties. 
        sol_kit_inter_ne = np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_ne)
        sol_kit_inter_te = np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_te)
        sol_kit_inter_z =  np.interp(sol_kit_grid, sol_kit_x_centered_grid, sol_kit_z)
        #SOL_KIT_inter_ne = spliner_ne(SOL_KIT_grid)
        #SOL_KIT_inter_Te = spliner_Te(SOL_KIT_grid)
        #SOL_KIT_inter_Z = spliner_Z(SOL_KIT_grid)

        if self.load_f1:
            self.switches['RESTART'] = "T"
            templating(tmpfilePath= os.path.join(self._run_path, 'INPUT/tmpSOL_KIT_SWITCHES.txt'),
            writePath=self._sol_kit_input_path, fileName="SWITCHES_INPUT.txt", parameters=self.switches)
            grid_path = os.path.join(self.previous_cycle_output_path, "".join(["GRIDS"]))
            v_grid = np.loadtxt(os.path.join(grid_path, "V_GRID.txt"))
            v_grid_width = np.loadtxt(os.path.join(grid_path, "V_GRID_WIDTH.txt"))
            self.interpolated_nx = len(sol_kit_inter_te)
            self.__initF_0(sol_kit_inter_te, sol_kit_inter_ne, v_grid, v_grid_width)
            self.__initRestartVector()

        else:
            np.savetxt(os.path.join(self._sol_kit_input_path, "DENS_INPUT.txt"), sol_kit_inter_ne)    
            np.savetxt(os.path.join(self._sol_kit_input_path, "TEMPERATURE_INPUT.txt"), sol_kit_inter_te)    
            np.savetxt(os.path.join(self._sol_kit_input_path, "Z_PROFILE_INPUT.txt"), sol_kit_inter_z)    
            np.savetxt(os.path.join(self._sol_kit_input_path, "X_GRID_INPUT.txt"), sol_kit_grid)
            # np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "ION_VEL_INPUT.txt"), SOL_KIT_grid)
            # np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "NEUT_HEAT_INPUT.txt"), SOL_KIT_heating)
    
    def getLastHeatFlow(self):
        """ Purpose: Gets the last heat flow from SOL-KiT after run finishes.
            Returns: Last heat flow with assumed format of cell-wall only definition of heat flow
        """

        qe_path =  os.path.join(self._run_path, "OUTPUT/HEAT_FLOW_X")
        max_index = findLargestIndex(qe_path)
        kinetic_heat_profile  = np.loadtxt(os.path.join(qe_path, "HEAT_FLOW_X_" + str(max_index).zfill(5) + '.txt'))

        #Include first and celll wall for no flow set to 0 per no thermal influx BC = 0 
        if self.boundary_condition =='noflow':
            kinetic_heat_profile = np.insert(kinetic_heat_profile, 0, 0.)

        #Periodic boundary condition removes last cell wall insert back and set to 0 via thermal influc BC
        self._kinetic_heat_profile = np.append(kinetic_heat_profile, 0.) * self.normalised_values['qe']
        
        return self._kinetic_heat_profile[::2] #here we assume heat flow is only defined on cell-walls 

    def __initRestartVector(self):
        """
        Purpose: To speed up div.q coupling init the next cycle with previous f_1
                done by starting from restart instead of loading in parameters. 
        NOTE:
        if ion things are included this has to be modified            
        """
        num_fields = 1
        num_h = self.init.yaml_file['Params']['L_max'] + 1 #References grid.f90:404
        num_0d = num_h*self.nv  + num_fields #References grid.f90:309
        restart_path =  os.path.join(self._run_path, "".join(["INPUT","/","RESTART","/", "VAR_VEC_INPUT.txt"]))
        restart_vector = np.loadtxt(restart_path)
        for i in range(self.interpolated_nx):
            F_0_POS = i*num_0d + self.nv * (num_h - 1) #References f_init.f90:219
            for v in range(self.nv):
                restart_vector[F_0_POS + v] = self.next_f0[i, v]
        
        np.savetxt(os.path.join(self._run_path, "".join(["INPUT","/","RESTART","/", "VAR_VEC_INPUT.txt"])),
                     restart_vector,
                     fmt ="%22.14e3") #Specific format expected by SOL-KiT
    def __initF_0(self, Te, ne, v_grid, v_grid_width ):
        self.nv = self.init.yaml_file['Params']['Nv'] 
        self.next_f0 = np.zeros((self.interpolated_nx, self.nv))
        f0 = np.zeros((self.interpolated_nx, self.nv))
        n_num = np.zeros(self.interpolated_nx)
        #Method from f_init.f90:215::227
        for i in range(self.interpolated_nx):
            for v in range(self.nv):
                f0[i,v] = ne[i] * ((np.pi * Te[i]) ** (- 3.00/2.00)) * np.exp( - (v_grid[v] ** 2)/ Te[i])
                n_num[i] = n_num[i] + 4.00 * np.pi * v_grid[v] ** 2 * v_grid_width[v] * f0[i,v]
            for v in range(self.nv):
                self.next_f0[i, v] = f0[i,v] * ne[i]/n_num[i]

