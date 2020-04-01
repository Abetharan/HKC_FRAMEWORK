""" 
Self-consistently runs the Rad-Hydro code HyKiCT
@author = Abetharan Antony
Last Update Date = 01/04/2020
"""
import numpy as np 
import os
import TmpFileCreator as tfc 
import subprocess
import Templating as temple 
import SetHydroInit as setinit
import yaml
from Fluid import Fluid

class HyKiCT(Fluid):

    def __init__(self, IO, nx_, laser_wave_length_, laser_power_, 
                dur_of_laser_, steps_, fluid_t_max_, initial_dt_,
                dt_global_max_, dt_global_min_, percentage_output_freq_,
                bounday_condition_, CoupleDivQ_, CoupleMulti_, PreHeatStartIndex_, PreHeatLastIndex_, FrontHeatStartIndex_, FrontHeatLastIndex_, initialise_start_files_run_ = True):
        self._templater = temple.Templating()
        self._nx = nx_
        self._f_io_obj = IO
        self._laser_wavelength = laser_wave_length_
        self._laser_power = laser_power_
        self._laser_loc = "left"
        self._dur_of_laser = dur_of_laser_
        self._nt = steps_
        self._fluid_time_max = fluid_t_max_
        self._initial_dt = initial_dt_
        self._dt_global_max = dt_global_max_
        self._dt_global_min = dt_global_min_
        self._output_freq = 1
        self._fluid_src_dir = IO._F_SRC_DIR
        self._cycle_path = IO.cycle_dump_path
        self._init_file_path = IO.fluid_input_path       
        self._base_dir = IO._BASE_DIR
        self._run_path = IO._RUN_PATH
        self._out_file_path = IO.fluid_output_path
        self._feos_material_1 = IO._F_FEOS_1_PATH 
        self._feos_material_2 = IO._F_FEOS_2_PATH
        self._cycle_dump_path = IO.cycle_dump_path
        self._cq = 2
        self._cfl = 0.85
        self._gamma = 5/3
        self._boundary_condition = bounday_condition_
        self._pre_heat_start_index = PreHeatStartIndex_
        self._pre_heat_last_index = PreHeatLastIndex_
        self._front_heat_start_index = FrontHeatStartIndex_
        self._front_heat_last_index = FrontHeatLastIndex_
        if self._fluid_time_max == 0:
            self._output_freq = 100
        else:
            self._output_freq = int(percentage_output_freq_ * 
                        (int(self._fluid_time_max)/int(self._dt_global_min)))
        

        self.makeTmpFiles()
        self.setSwitches(coupleDivQ_ = CoupleDivQ_, coupleMulti_ = CoupleMulti_)    
        if initialise_start_files_run_:
            self.setSwitches(coupleDivQ_= False)
            self._nt = 0
            self._laser_power = 0            
            self.copyELH1()            
        hydroparam = setinit.set_hydro_init(self._nx, self._cq, self._gamma, self._cfl, self._laser_wavelength,  
                                            self._laser_power, self._dur_of_laser, self._laser_loc, self._nt, 
                                            self._fluid_time_max, self._initial_dt, self._dt_global_max, self._dt_global_min, 
                                            self._output_freq, self._boundary_condition, self._init_file_path, 
                                            self._out_file_path, self._switch_path, self._feos_material_1, 
                                            self._feos_material_2, self._pre_heat_start_index, self._pre_heat_last_index,
                                            self._front_heat_start_index, self._front_heat_last_index)

        # Handling templating to create the init file for fluid code
        self._templater.templating(tmpfilePath=self._base_dir + '/tmpHydroParameterInit.txt',
                writePath=self._cycle_dump_path, fileName="HydroParameterInit.txt", parameters=hydroparam)

    def setSwitches(self, viscosity_on_ = False, velocity_on_ = False, 
                    heat_conduction_on_ = True, exchange_on_ = False,
                    bremsstrahlung_on_ = False, inv_brem_on_ = False, 
                    single_temperature_on_ = True, multi_material_ = False,
                    ideal_gas_ = True, fully_ionized_ = True, coupleDivQ_ =True, coupleMulti_ = False):
        """
        Purpose : Set all neccessary switches for ELH-1 for coupling couple{} shoud be set to True.
        Args: 
            viscosity_on_ = include artificial viscosity
            velocity_on_ = calculate velocity
            heat_conduction_on = Calculate heat conduction
            exchange_on = calculate Exchange term between electorn-ion
            bremsstrahlung_on = find bremsstrahlung 
            inv_brem_on = Laser abosprtion via inv_brem
            single_temperature_on = Te==Ti
            multi_material = 2 materials with different Ar
            ideal_gas = EoS is ideal gas
            fully_ionized = run with Zbar=Z
            coupleDivQ = div.q coupoling 
            coupleMulti_ = multipler coupling.
        """              

        self.viscosity = viscosity_on_
        self.velocity = velocity_on_
        self.heat_conduction = heat_conduction_on_
        self.exchange = exchange_on_
        self.bremsstrahlung = bremsstrahlung_on_
        self.inv_brem = inv_brem_on_
        self.single_temp_mode = single_temperature_on_
        self.couple_div_q = coupleDivQ_
        self.couple_multi = coupleMulti_

        self.multi_material = multi_material_
        self.ideal_gas_mode = ideal_gas_
        self.fully_ionized_mode = fully_ionized_
        self.isothermal_mode = False
        self.adibatic_mode = False
        self.p_dv_work_off = False
        self.AdapativeTimeStep = False
        switches = {
            'Viscosity' : str(self.viscosity).lower(),
            'Velocity' : str(self.velocity).lower(),
            'HeatConduction' : str(self.heat_conduction).lower(),
            'Exchange' : str(self.exchange).lower(),
            'Bremsstrahlung' : str(self.bremsstrahlung).lower(),
            'InvBremsstrahlung' : str(self.inv_brem).lower(),
            'IsothermalMode' : str(self.isothermal_mode).lower(),
            'AdiabaticMode' : str(self.adibatic_mode).lower(),
            'pDvWorkOff' : str(self.p_dv_work_off).lower(),
            'CoupleDivQ' : str(self.couple_div_q).lower(),
            'CoupleMulti' : str(self.couple_multi).lower(),
            'SingleTemperature' : str(self.single_temp_mode).lower(),
            'MultiMaterial' : str(self.multi_material).lower(),
            'IdealGas' : str(self.ideal_gas_mode).lower(),
            'FullyIonized':str(self.fully_ionized_mode).lower(),
            'AdapativeTimeStep':str(self.AdapativeTimeStep).lower()
        }
        self._templater.templating(tmpfilePath= self._base_dir + '/tmpFluidSwitch.txt',
                writePath=self._cycle_dump_path, fileName="HydroSwitches.txt", parameters=switches)
        self._switch_path = os.path.join(self._cycle_dump_path,"HydroSwitches.txt")
    
    def copyELH1(self):
        """ Purpose: Copy EHL1 exe. """
        import shutil
        shutil.copy(self._fluid_src_dir + '/HyKiCT', self._run_path)
    
    def makeTmpFiles(self):
        """ Purpose: Create tmpfiles to be templated 
        """
        tfc.ELH1Parameter(self._base_dir)
        tfc.ELH1Switches(self._base_dir)
    
    def ELH1Run(self):
        """ Purpose: Run ELH1 with parameters set previously
        """
        os.chdir(self._run_path)
        cmd = ['./ELH1','-p',
                self._cycle_dump_path+'/HydroParameterInit.txt']
        super().Execute(cmd, self._cycle_path)
    
