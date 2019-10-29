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

kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
epsilon_0 = 8.854188E-12    # Vacuum dielectric constant
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")
class SOL_KIT(Kinetic):
    
    def __init__(self,IO, np_, nv_, nx_, nt_, pre_step_nt_, dx_, dv_, v_multi_, dt_, pre_step_dt_, lmax_, save_freq_, norm_Z_, norm_Ar_, norm_ne_, norm_Te_):

        self._Kinetic = Kinetic()
        self._templater = Templating() 
        self._kinetic_io_obj = IO
    
        self._run_name = IO._RUN_NAME
        self._run_path = IO._RUN_PATH
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

        self.normalised_values = None

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


        me = me #9.10938E-31            # Electron mass in kg
        el_charge = e #1.602189E-19    # Electron charge
        r_b = bohr_radi #5.29E-11              # Bohr radius in metres
        planck_h = planck_h #6.62607004E-34   # Planck constant

        sigma_0 = math.pi * r_b ** 2  # Cross-section normaliation
        gamma_ee_0 = el_charge ** 4 / (4 * math.pi * (me * epsilon_0) ** 2)

        T_J = self._norm_Te * el_charge

        gamma_ei_0 = self._norm_Z ** 2 * gamma_ee_0

        v_t = math.sqrt(2.0 * T_J / me)                  # Velocity normalization
        time_norm = v_t ** 3 / (gamma_ei_0 * self._norm_ne * lambda_ei(1.0, 1.0, T_norm = self._norm_Te, n_norm = self._norm_ne, Z_norm = self._norm_Z)) # Time normalization
        x_0 = v_t * time_norm               # Space normalization

        e_0 = me * v_t / (el_charge * time_norm) # Electric field normalization
        q_0 = me * np * (v_t ** 3)

        dict = {}

        dict['ne'] = self._norm_ne
        dict['te'] = self._norm_Te
        dict['z'] = self._norm_Z

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
        tfc.SOL_KIT_GRID_INPUT(self._run_name)
        tfc.SOL_KIT_NORMS(self._run_name)
        tfc.SOL_KIT_SOLVER_PARAMS(self._run_name)
        tfc.SOL_KIT_SWITCHES(self._run_name)

    def copySOL_KiT(self):

        shutil.copy(self._base_dir + '/SOL-KiT', self._run_path)
        shutil.copytree(self._base_dir + '/OUTPUT', self._run_path)


    def setSOL_KITFiles(self):
        INPUT_PATH = os.path.join(self._run_path, "INPUT")
        os.makedirs(INPUT_PATH)

        shutil.move(self._run_path + '/tmpSOL_KIT_NORMS.txt', INPUT_PATH + '/SOLVER_PARAMS_INPUT.txt')
        
        grid = setParams.GRID(nt = self._nt, prestepnt = self._pre_step_nt, dt = self._dt,
                            predt = self._pre_step_dt, save_freq = self._save_freq, lmax = self._lmax,
                            nv = self._nv, dv = self._dv, v_multi = self._v_multi, nx = self._nx,
                            dx = self._dx, smalldx=0) 
        
        norms = setParams.NORMS(Z = self._norm_Z, Ar = self._norm_Ar, Density = self._norm_ne, 
                                Temp = self._norm_Te)
        
        switches = setParams.Switches(EE_L = True, DIAG_EE = False, 
                                    PERIODIC = False, FIXEDUP = False, 
                                    FIXEDDOWN = False, NOFLOWUP = False, 
                                    NOFLOWDOWN = False, IMPACTMODE = False, 
                                    COLDIONS = False)
                                    
        
        self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_GRID.txt',
            writePath=INPUT_PATH, fileName="GRID_INPUT.txt", parameters=grid)
        self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_NORMS.txt',
            writePath=INPUT_PATH, fileName="NORMALIZATION_INPUT.txt", parameters=norms)
        #self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_SOLVER_PARAMS.txt',
        #writePath=INPUT_PATH, fileName="", parameters=)
        self._templater.templating(tmpfilePath= self._run_path + '/tmpSOL_KIT_SWITCHES.txt',
            writePath=INPUT_PATH, fileName="SWITCHES_INPUT.txt", parameters=switches)

    def fluidToSOL_KIT(self):

        
    def SOL_KITNormsToSI(self):



    def moveSOL_KiTFiles(self):