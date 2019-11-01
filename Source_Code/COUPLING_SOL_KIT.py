import math
import sys
import pprint
import numpy as np 
import os
import fnmatch
import string
import subprocess
import shutil 
from scipy import constants
import TmpFileCreator as tfc
import SetSOL_KITParams as setParams
from Templating import Templating
from Kinetic import Kinetic
from IO import IO
from distutils.dir_util import copy_tree
kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
epsilon_0 = 8.854188E-12    # Vacuum dielectric constant
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")
class SOL_KIT(Kinetic):
    
    def __init__(self,IO, np_, nv_, nx_, nt_, pre_step_nt_, dx_, dv_, v_multi_, dt_, pre_step_dt_,
                    lmax_, save_freq_, norm_Z_, norm_Ar_, norm_ne_, norm_Te_, get_sol_kit_):

        self._Kinetic = Kinetic()
        self._templater = Templating() 
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
        
        self._np = np_
        self._nv = nv_
        self._nx = nx_
        self._nt = nt_
        self._pre_step_nt = pre_step_nt_
        self._dx = dx_
        self._dv = dv_ 
        self._v_multi = v_multi_
        self._dt = dt_
        self._pre_step_dt = pre_step_dt_
        self._lmax = lmax_
        self._save_freq = save_freq_
        self._norm_Z = norm_Z_
        self._norm_Ar = norm_Ar_
        self._norm_ne = norm_ne_
        self._norm_Te = norm_Te_        
        self.boundary_condition = None
        self.normalised_values = None

        self.makeTmpFiles()
        self.copySOL_KiT()
        self.normalisation()
    def SOL_KITRun(self):
        os.chdir(self._run_path)
        cmd = ["mpirun", "-np", str(self._np), "./SOL-KiT"]    
        super().Execute(cmd, self._cycle_path)
    
    def normalisation(self):

        def lambda_ei(n, T, T_norm = 10, n_norm = 0.25E20, Z_norm = 1.0):
            if T * T_norm < 10.00 * Z_norm ** 2:

                result = 23.00 - math.log(math.sqrt(n * n_norm * 1.00E-6) * Z_norm * (T * T_norm) ** (3.00/2.00))

            else:

                result = 24.00 - math.log(math.sqrt(n * n_norm * 1.00E-6) / (T * T_norm))   

            return result


        #9.10938E-31            # Electron mass in kg
        el_charge = e #1.602189E-19    # Electron charge
        r_b = bohr_radi #5.29E-11              # Bohr radius in metres

        sigma_0 = math.pi * r_b ** 2  # Cross-section normaliation
        gamma_ee_0 = el_charge ** 4 / (4 * math.pi * (me * epsilon_0) ** 2)

        T_J = self._norm_Te * el_charge

        gamma_ei_0 = self._norm_Z ** 2 * gamma_ee_0

        v_t = math.sqrt(2.0 * T_J / me)                  # Velocity normalization
        time_norm = v_t ** 3 / (gamma_ei_0 * self._norm_ne * lambda_ei(1.0, 1.0, T_norm = self._norm_Te, n_norm = self._norm_ne, Z_norm = self._norm_Z)) # Time normalization
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
        dict['vte'] = v_t
        dict['tau_ei'] = time_norm
        dict['nu_ei'] = 1/time_norm
        dict['lambda_mfp'] = x_0
        dict['qe'] = q_0
        self.normalised_values = dict
    def makeTmpFiles(self):
        tfc.SOL_KIT_GRID_INPUT(self._run_path)
        tfc.SOL_KIT_NORMS(self._run_path)
        tfc.SOL_KIT_SOLVER_PARAMS(self._run_path)
        tfc.SOL_KIT_SWITCHES(self._run_path)

    def copySOL_KiT(self):
        if self.copy_sol_kit:
            shutil.copy(self._src_dir + '/SOL-KiT', self._run_path)
        
        OUTPUT_path = os.path.join(self._run_path, "OUTPUT")
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


    def setSOL_KITFiles(self, switch_boundary_condition_, switch_EE_L_ = True, switch_DIAG_EE_ = False, switch_IMPACT_MODE_ = False, switch_COLD_IONS_ = False):
        INPUT_PATH = os.path.join(self._run_path, "INPUT")
        os.makedirs(INPUT_PATH)

        #shutil.move(self._run_path + '/tmpSOL_KIT_SOLVER_PARAMS.txt', INPUT_PATH + '/SOLVER_PARAMS_INPUT.txt')
        shutil.copy(self._src_dir + "/INPUT/SOLVER_PARAMS_INPUT.txt", INPUT_PATH)
        shutil.copy(self._src_dir + "/INPUT/F_INIT_INPUT.txt", INPUT_PATH)
        grid = setParams.GRID(nt = self._nt, prestepnt = self._pre_step_nt, dt = self._dt,
                            predt = self._pre_step_dt, save_freq = self._save_freq, lmax = self._lmax,
                            nv = self._nv, dv = self._dv, v_multi = self._v_multi, nx = self._nx,
                            dx = self._dx, smalldx=0.0) 
        
        norms = setParams.NORMS(Z = self._norm_Z, Ar = self._norm_Ar, Density = self._norm_ne, 
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
                                    COLDIONS = switch_COLD_IONS_)
                                
        self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_GRID.txt',
            writePath=INPUT_PATH, fileName="GRID_INPUT.txt", parameters=grid)
        self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_NORMS.txt',
            writePath=INPUT_PATH, fileName="NORMALIZATION_INPUT.txt", parameters=norms)
        #self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_SOLVER_PARAMS.txt',
        #writePath=INPUT_PATH, fileName="", parameters=)
        self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_SWITCHES.txt',
            writePath=INPUT_PATH, fileName="SWITCHES_INPUT.txt", parameters=switches)


    def _findLargestIndex(self, path):
        files = os.listdir(path)
        from string import ascii_letters
        extracted_indicies = []
        for x in files:
            if x == ".keep":
                continue
            extracted_indicies.append(int(x.strip(ascii_letters + '_.')))    
        #files.rstrip(ascii_letters)
        #files.rstrip('_.')
        #files = int(files)
        return(max(extracted_indicies))

    def _SOL_KITNormsToSI(self, var, var_name):
        return(var * self.normalised_values[var_name])

    def moveSOL_KITFiles(self):
        import shutil
        #move INPUT files
        shutil.move(self._SOL_KIT_INPUT_PATH, self._kinetic_io_obj.kinetic_input_path)
        ###HEATING INPUT FOR FUTURE USE
        # shutil.move(os.path.join(self._SOL_KIT_INPUT_PATH, "NEUT_AND_HEAT_INPUT.txt"))    
        
        #move OUTPUT files
        shutil.move(self._SOL_KIT_OUTPUT_PATH, self._kinetic_io_obj.kinetic_output_path)
        
        
        
    def InitSOL_KITFromHydro(self):
        
        last_index = self._findLargestIndex(os.path.join(self._fluid_output_path, "ELECTRON_TEMPERATURE"))
        f_x_grid = np.loadtxt(self._fluid_output_path + "CELL_WALL_X/Cell_Wall_X_" + str(last_index) + ".txt")
        f_x_centered_grid = np.loadtxt(self._fluid_output_path + "CELL_CENTRE_X/Cell_Centre_X_" + str(last_index) + ".txt")
        f_v = np.loadtxt(self._fluid_output_path + "VELOCITY/Velocity_" + str(last_index) + ".txt")
        f_ne = np.loadtxt(
                self._fluid_output_path + "ELECTRON_NUMBER_DENSITY/Ne_" + str(last_index) + ".txt")
        f_Te = np.loadtxt(
             self._fluid_output_path + "ELECTRON_TEMPERATURE/Te_" + str(last_index) + ".txt")
        f_laser = np.loadtxt(
            self._fluid_output_path + "INVERSE_BREM/Inverse_Brem_" + str(last_index) + ".txt")
        f_brem = np.loadtxt(self._fluid_output_path + "BREM/Brem_" + str(last_index) + ".txt")
        f_Z = np.loadtxt(self._fluid_output_path + "ZBAR/Zbar_" + str(last_index) + ".txt")

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

        if self.boundary_condition == 'periodic':
            SOL_KIT_grid = full_x_grid[:-1]
        
        elif self.boundary_condition == "noflow":
            SOL_KIT_grid = full_x_grid[1:-1]

        SOL_KIT_ne = f_ne / (self.normalised_values['ne'])
        SOL_KIT_Te = (f_Te * (kb/e)) / self.normalised_values['Te']
        #SOL_KIT_laser = (f_laser * f_density) / power_density
        #SOL_KIT_brem = (f_brem * f_density) / power_density
        SOL_KIT_Z = f_Z
        # SOL_KIT_HEATING = SOL_KIT_laser + SOL_KIT_brem
        from scipy.interpolate import CubicSpline, interp1d
        spliner_ne = CubicSpline(SOL_KIT_x_centered_grid, SOL_KIT_ne)
        spliner_Te = CubicSpline(SOL_KIT_x_centered_grid, SOL_KIT_Te)
        spliner_Z = CubicSpline(SOL_KIT_x_centered_grid, SOL_KIT_Z)
        SOL_KIT_inter_ne = spliner_ne(SOL_KIT_grid, extrapolate= True)
        SOL_KIT_inter_Te = spliner_Te(SOL_KIT_grid, extrapolate= True)
        SOL_KIT_inter_Z = spliner_Z(SOL_KIT_grid, extrapolate= True)
        # import matplotlib.pyplot as plt

        np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "DENS_INPUT.txt"), SOL_KIT_inter_ne, fmt = '%.19f')    
        np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "TEMPERATURE_INPUT.txt"), SOL_KIT_inter_Te, fmt = '%.19f')    
        np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "Z_PROFILE_INPUT.txt"), SOL_KIT_inter_Z, fmt = '%.19f')    
        np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "X_GRID_INPUT.txt"), SOL_KIT_grid, fmt = '%.18f')
        # np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "ION_VEL_INPUT.txt"), SOL_KIT_grid)
        # np.savetxt(os.path.join(self._SOL_KIT_INPUT_PATH, "NEUT_HEAT_INPUT.txt"), SOL_KIT_heating)

    def SOL_KITInitNextHydroFiles(self):
        qe_path =  os.path.join(self._run_path, "OUTPUT/HEAT_FLOW_X")
        max_index = self._findLargestIndex(qe_path)
        max_index_reformated = '{:>05d}'.format(max_index)
        qe = np.loadtxt(os.path.join(qe_path, "HEAT_FLOW_X_" + max_index_reformated + '.txt'))
        
        ###
        def Heatflow(electron_thermal_flux, mass):
            nx = len(mass)
            HeatConductionE = np.zeros(nx)
            electron_thermal_flux[0] = 0.
            electron_thermal_flux[-1] = 0.
            for i in range(0, nx, 2):   
                HeatConductionE[i] = -(electron_thermal_flux[i + 1] - electron_thermal_flux[i]) / mass[i]
            return(HeatConductionE)

        mass = np.loadtxt(self._kinetic_io_obj.fluid_input_path + "/mass.txt")        
        electronheatflow= Heatflow(qe, mass)
        
        np.savetxt(os.path.join(self._kinetic_io_obj.next_fluid_input_path, "qe.txt"), electronheatflow)    

        largest_fluid_index = self._findLargestIndex(os.path.join(self._fluid_output_path, "ELECTRON_TEMPERATURE"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "CELL_WALL_X/Cell_Wall_X_" + str(largest_fluid_index) +".txt")    
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"coord.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "VELOCITY/Velocity_" + str(largest_fluid_index) +".txt") 
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"velocity.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "DENSITY/Density_" + str(largest_fluid_index) +".txt")  
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"density.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "ELECTRON_TEMPERATURE/Te_" + str(largest_fluid_index) +".txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"electron_temperature.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "ION_TEMPERATURE/Ti_" + str(largest_fluid_index) +".txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"ion_temperature.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "MASS/Mass_.txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"mass.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "ZBAR/Zbar_" + str(largest_fluid_index) +".txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path, "Z.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "MASS_NUMBER/Ar_" + str(largest_fluid_index) +".txt")  
                                    ,os.path.join(self._kinetic_io_obj.next_fluid_input_path, "Ar.txt"))