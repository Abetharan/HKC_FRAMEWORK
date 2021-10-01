""" 
Hydro-Kinetic Coupling.
Heat-Flow coupling between Hydro and Kinetic codes. 
@author = Abetharan Antony
"""
import argparse
import atexit
import copy
from distutils.dir_util import copy_tree
import logging
import math
import numpy as np
from PIL import Image, ImageFont, ImageDraw
import os
import pathlib
import signal
import shutil
import sys
import time
import utils as util
from couple_sol_kit.sol_kit import SOL_KIT
from couple_hykict.hykict import HyKiCT
from couple_methods import Multiplier, DivQ, Subtract
from scipy import constants
# import time_stepper as time_step 
BOLTZMANN_CONSTANT = constants.value("Boltzmann constant")
ELECTRON_MASS = constants.value("electron mass")
PROTON_MASS = constants.value("proton mass")
ELEMENTARY_CHARGE = constants.value("elementary charge")
VACUUM_PERMITTIVITY = 8.854188E-12    # Vacuum dielectric constant
PLANCK_CONSTANT = constants.value("Planck constant")
BOHR_RADIUS = constants.value("Bohr radius")
class Coupler:
    """
    main class i.e. where all objects are passed into and the coupling is run.
    """
    def __init__(self, init_file_path):
        self.yml_init_file_path = init_file_path
        self.pre_heat_present = False
        self.kinetic_time_taken = []
        self.fluid_time_taken = []
        self.cycle_time_taken = []
        self.init = util.Input(self.yml_init_file_path)
        self.coupling_message = "No Method Chosen"
        self.critical_density = None
        self.laser_dir = None
    def startprint(self, fluid_code, kinetic_code):

        ShowText = ('COUPLING ' + fluid_code + '-' 
                    + kinetic_code)

        font = ImageFont.load_default().font #load the font
        #font = ImageFont.truetype("/home/abetharan/Downloads/arial.ttf",15)
        size = font.getsize(ShowText)  #calc the size of text in pixels
        image = Image.new('1', size, 1)  #create a b/w image
        draw = ImageDraw.Draw(image)
        draw.text((0, 0), ShowText, font=font) #render the text to the bitmap
        for rownum in range(size[1]): 
            line = []
            for colnum in range(size[0]):
                if image.getpixel((colnum, rownum)): line.append(' '),
                else: line.append('#'),
            print(''.join(line))

    def prettyPrint(self, text, color = None):
        length_to_print = 100 + len(text)
        if color is None:
            print("#"*length_to_print)
            print('\033[1m' + '#'*50 + text +'#'*50 + '\033[0m')
        
        else:
            print('\033[31m' + "#"*length_to_print + '\033[0m')
            print('\033[31m' + '#'*50 + text +'#'*50 + '\033[0m')
        
        print('\n')   
    
    def cleanUpHandlers(self):
        for handler in self.logger.handlers:
            handler.close()
            self.logger.removeHandler(handler)  

    def storeToHdf5(self, cycle_no):
        self.io_obj._createHDF5()
        self.logger.info("Store to HDF5")
        self.fluid_obj.storeToHdf5(self.io_obj.hdf5_file, cycle_no)
        self.kin_obj.storeToHdf5(self.io_obj.hdf5_file, cycle_no)
        self.io_obj.hdf5_file.close()
    
    
    def kineticStep(self, cycle_no, conv_heat_flow,
                        fluid_x_grid, fluid_x_centered_grid, 
                                fluid_Te, fluid_ne, fluid_Z):
        """
        Purpose: Init Kinetic, Runs Kinetic.
        Args: 
            cycle_no = Current Cycle Number
            conv_heat_flow = Convergence checking heat-flow typically Spitzer/SNB
            fluid_x_grid = Cell-Wall centered Grid from fluid code 
            fluid_x_centered_grid = Cell centered Grid from fluid code 
            fluid_Te = Temperature from fluid code 
            fluid_ne = Electron Density from fluid code 
            fluid_Z = Zbar from fluid code 
        Returns : Kinetic Heat-Flow 
        """

        tkin_start = time.time()
        self.prettyPrint(' RUNNING ' + self.init.yaml_file['Codes']['Kinetic_code'], color = True)
        self.logger.info("Setting Kinetic Parameters")
        if cycle_no == 0 or self.first_pass:
            #input and output unchanged 
            self.logger.info("Set files... These are UNCHANGED EXCEPT IF LOAD_f1")
            self.kin_obj.setFiles()

        self.kin_obj._kinetic_input_path = self.io_obj.kinetic_input_path
        self.kin_obj._kinetic_output_path = self.io_obj.kinetic_output_path
        self.kin_obj._cycle_dump_path = self.io_obj.cycle_dump_path
        self.logger.info("Kinetic Paths")
        self.logger.info("Kinetic input path {}".format(self.io_obj.kinetic_input_path))
        self.logger.info("Kinetic output path {}".format(self.io_obj.kinetic_output_path))

        if self.init.yaml_file['Mode']['Rescale_init_f0']:
            if cycle_no > 0:
                self.logger.info("Engaging MODE: Rescale f0 and load in previous f_l")
                #all_kinetic_output_paths contain all paths for kinetic output, require last cycle for load F1, thus cycle_no - 1
                self.kin_obj.previous_kinetic_output_path = self.io_obj.all_kinetic_output_path[cycle_no - 1] 
                self.kin_obj.rescale_f0_init = True                                                              
        elif self.init.yaml_file['Mode']['Load_f1']:
            if cycle_no > 0:
                self.logger.info("Engaging MODE: Load F1")
                #all_kinetic_output_paths contain all paths for kinetic output, require last cycle for load F1, thus cycle_no - 1
                self.kin_obj.previous_kinetic_output_path = self.io_obj.all_kinetic_output_path[cycle_no - 1] 
                self.kin_obj.load_f1 = True                                                              
            
        self.logger.info("Initialize Kinetic")
        self.kin_obj.sh_heat_flow = conv_heat_flow

        skipped = self.kin_obj.initFromHydro(fluid_x_grid, fluid_x_centered_grid, 
                            fluid_Te, fluid_ne, fluid_Z, critical_density = self.critical_density, laser_dir = self.laser_dir)
        if skipped and cycle_no > 0:
            self.logger.info("FAILED TO Engage MODE: Load F1 DUE TO CHANGE IN GRID SIZE")


        self.logger.info("Start Kinetic")
        self.kin_obj.Run()
        self.logger.info("End Kinetic")
        #Set hfct to contain the vfp heat flow to do the neccessary coupling calcs.
        self.logger.info("Fetch Last Heat Flow")
        last_heat_flow = self.kin_obj.getLastHeatFlow()
        self.logger.info("Move files Kinetic")
        self.kin_obj.moveFiles()
        tkin_end = time.time()
        self.kinetic_time_taken.append(tkin_end - tkin_start)
        return last_heat_flow

    def fluidStep(self, outpath, next_input_path, no_copy = False, update_time = True):
        """
        Purpose: Sets Fluid files and RUns Fluid code
        Args:
            outpath : Output path to save fluid output. This only 
                        changes when advanced modes such as leap-frog
                        is engaged.
            next_input_path : Path to init next Hydro Cycle
        Returns:
            Last step thermodynamic quantities, ionisation, spatial grid and
            Total physical simulation time. 
        """

        self.prettyPrint(' RUNNING ' + self.init.yaml_file['Codes']['Fluid_code'], color = True)
        self.logger.info("Setting Fluid Parameters")
        tfluid_start = time.time()
        #Set Paths 
        self.fluid_obj._cycle_dump_path = self.io_obj.cycle_dump_path
        self.fluid_obj._fluid_input_path = self.io_obj.fluid_input_path
        self.fluid_obj._fluid_output_path = outpath

        self.fluid_obj.init.yaml_file['Paths']['Init_Path'] = self.io_obj.fluid_input_path
        self.fluid_obj.init.yaml_file['Paths']['Laser_Profile_Path'] = os.path.join(self.io_obj.fluid_input_path, "laser_profile.txt")
        self.fluid_obj.init.yaml_file['Paths']['Out_Path'] = outpath

        self.logger.info("Fluid input path {}".format(self.fluid_obj.init.yaml_file['Paths']['Init_Path']))
        self.logger.info("Fluid output path {}".format(self.fluid_obj.init.yaml_file['Paths']['Out_Path']))

        ################
        #Any other parameter updates that needs to be set can be added here.
        #################
        self.fluid_obj.setFiles()
        self.logger.info("Start Fluid")
        self.fluid_obj.Run()
        self.logger.info("End Fluid")
        self.logger.info("Get Last Fluid Quants")
        (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
            fluid_Z, _, fluid_mass, sim_time, f_specific_heat, f_Ar) = self.fluid_obj.getLastStepQuants(update_time) 
        if no_copy:
            return ((fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
            fluid_Z, _, fluid_mass, f_specific_heat, f_Ar, sim_time))
        self.fluid_obj.setTimes()
        self.logger.info("Copying files to next step Init: {}".format(next_input_path))
        self.fluid_obj.initHydro(next_input_path)
        tfluid_end = time.time()
        self.fluid_time_taken.append(tfluid_end - tfluid_start)

        return (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
            fluid_Z, _, fluid_mass,f_specific_heat,f_Ar, sim_time) 
    
    def returnConvHeatFlow(self,fluid_x_grid,fluid_x_centered_grid, fluid_Te,fluid_ne , fluid_Z,fluid_mass ):
        """
        Purpose: Returns the convergence/coupling heat flows  
        Args:
            Grid + Thermodynamic quants
        Notes:
        """
        self.hfct_obj.electron_temperature = fluid_Te
        self.hfct_obj.electron_number_density = fluid_ne
        self.hfct_obj.zbar = fluid_Z
        self.hfct_obj.cell_wall_coord = fluid_x_grid
        self.hfct_obj.cell_centered_coord = fluid_x_centered_grid
        self.hfct_obj.mass = fluid_mass
        self.hfct_obj.lambda_ei(self.hfct_obj.electron_temperature * (BOLTZMANN_CONSTANT/ELEMENTARY_CHARGE), 
                            self.hfct_obj.electron_number_density,
                            self.hfct_obj.zbar)
        self.logger.info("HFCT Spitzer Calculation")
        self.hfct_obj.spitzerHarmHeatFlow()
        conv_heat_flow = self.hfct_obj.spitzer_harm_heat
        if self.fluid_obj.init.yaml_file['Switches']['SNBHeatFlow']:
            self.hfct_obj.snb = True
            self.hfct_obj.snb_heat_flow(self.fluid_obj.init.yaml_file['FixedParameters']['ng'],self.fluid_obj.init.yaml_file['FixedParameters']['MaxE'], 2)
            conv_heat_flow = self.hfct_obj.q_snb 
        return conv_heat_flow

    def generalCoupling(self, cycle_no, cycles,
                         run_only_fluid = False):
        """
        Purpose: Implements the General Coupling Algorithm
        Args:
            cycle_no = current cycle number
            cycles = Max Cycles
        Notes:
            Algorithm works as follows
            -> Run Fluid coupled and
                -> If Start Kin will be run with 0 steps 
            -> Send last step fluid output 
                to Kinetic and next fluid init. 
            -> Run Kinetic
            -> Init coupling params to next fluid step
            -> Finish 
        """
        no_copy = False
        if cycle_no == cycles - 1: 
            no_copy = True

        (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
        fluid_Z, _, fluid_mass, fluid_specific_heat,fluid_Ar, sim_time) = self.fluidStep(self.io_obj.fluid_output_path,
                                                            self.io_obj.next_fluid_input_path,  
                                                            no_copy = no_copy)
        if self.init.yaml_file['Mode']['Limit_density']:
            self.critical_density = self.init.yaml_file['Coupling_params']['nc_multiplier']* (1114326918632954.5 / pow(self.fluid_obj.init.yaml_file['LaserParams']['Wavelength'], 2)) #Hard-coded limit to 10*nc
            if any(fluid_ne >= self.critical_density):
                self.laser_dir = self.fluid_obj.laser_direction
                self.couple_obj.limit_density = True
            else:
                self.critical_density = None
                self.laser_dir = None

        if not run_only_fluid:
            conv_heat_flow = self.returnConvHeatFlow(fluid_x_grid,fluid_x_centered_grid, 
                                    fluid_Te,fluid_ne , fluid_Z,fluid_mass)

            vfp_heat = self.kineticStep(cycle_no, conv_heat_flow,
                        fluid_x_grid, fluid_x_centered_grid, 
                                fluid_Te, fluid_ne, fluid_Z)

            self.couple_obj.method(conv_heat_flow, vfp_heat, 
                            laser_dir = self.laser_dir, mass = fluid_mass, cell_wall_coord = fluid_x_grid,
                            q_snb = self.hfct_obj.q_snb)
            
            self.logger.info(self.coupling_message)
            self.couple_obj.setCoupleParams(self.io_obj.next_fluid_input_path, fluid_yaml = self.fluid_obj.init.yaml_file, Te = fluid_Te,
                                        no_negative = self.init.yaml_file['Mode']["no_negative_multiplier"], nn_via_material = self.init.yaml_file['Mode']["nn_via_material"],
                                        nn_via_interpolation = self.init.yaml_file['Mode']['nn_via_interpolation'],x = fluid_x_grid, Ar = fluid_Ar, 
                                        material = self.init.yaml_file['Coupling_params']["material_Ar_pacify"])

        if self.init.yaml_file['Mode']['AdaptiveTimeStep']:
            if self.init.yaml_file['Mode']['Couple_divq']:
                dtmax = self.time_stepper.divqtimestepper(fluid_Te, fluid_specific_heat,self.couple_obj.HeatConductionE)
            elif self.init.yaml_file['Mode']['Couple_multi']:
                dtmax = self.time_stepper.multitimestepper(fluid_x_grid,fluid_x_centered_grid, fluid_Te, fluid_ne, fluid_Z, fluid_specific_heat, fluid_mass, self.couple_obj.q_vfp_q_sh_multipliers)
                self.logger.info("t_cycle is {}".format(dtmax))
            self.fluid_obj.updateTime(new_cycle_step = dtmax)
            self.fluid_obj.setTimes()

            
    def leapFrog(self, cycle_no, cycles):
        """
        Purpose: Implements the Leap-Frog Coupling Algorithm
        Args:
            cycle_no = current cycle number
            cycles = Max Cycles
        Notes:
            Algorithm works as follows
            -> Execute General coupling step 
                with HyKiCT run as a normal
                fluid code 
            -> Send coupling params to 
                original F_init
            -> Run HyKiCT in Coupled mode 
            -> Finish 
        """
        self.fluid_obj.setNoCoupleSwitches()
        no_copy = False
        if cycle_no == cycles - 1: 
            no_copy = True

        (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
        fluid_Z, _, fluid_mass, fluid_specific_heat, sim_time) = self.fluidStep(self.io_obj.intermediate_fluid_outpath,
                                                            self.io_obj.next_fluid_input_path,  
                                                            no_copy = no_copy, update_time = False)
        conv_heat_flow = self.returnConvHeatFlow(fluid_x_grid,fluid_x_centered_grid, 
                                fluid_Te,fluid_ne , fluid_Z,fluid_mass)
        if self.init.yaml_file['Mode']['Limit_density']:
            self.critical_density = 10 * (1114326918632954.5 / pow(self.fluid_obj.init.yaml_file['LaserParams']['Wavelength'], 2)) #Hard-coded limit to 10*nc
            if any(fluid_ne >= self.critical_density):
                self.laser_dir = self.fluid_obj.laser_direction
                self.couple_obj.limit_density = True
            else:
                self.critical_density = None
                self.laser_dir = None
        vfp_heat = self.kineticStep(cycle_no, conv_heat_flow,
                    fluid_x_grid, fluid_x_centered_grid, 
                            fluid_Te, fluid_ne, fluid_Z)

        (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
    fluid_Z, _, fluid_mass, sim_time)  = self.fluid_obj.returnInitValues()

        conv_heat_flow = self.returnConvHeatFlow(fluid_x_grid,fluid_x_centered_grid, 
                                fluid_Te,fluid_ne , fluid_Z,fluid_mass)
        
        self.couple_obj.method(conv_heat_flow, vfp_heat, 
                        laser_dir = self.laser_dir, mass = fluid_mass, cell_wall_coord = fluid_x_grid,
                        q_snb = self.hfct_obj.q_snb)
        
        self.logger.info(self.coupling_message)
        self.couple_obj.setCoupleParams(self.io_obj.fluid_input_path, fluid_yaml = self.fluid_obj.init.yaml_file,Te = fluid_Te, no_negative = self.init.yaml_file['Mode']["no_negative_multiplier"])
        
        #Run Fluid in coupled config 
        self.fluid_obj.revertStartKinSwitches()
        _ = self.fluidStep(self.io_obj.fluid_output_path, self.io_obj.next_fluid_input_path, no_copy = no_copy)


    def operatorSplit(self, cycle_no, cycles):
        """
        Purpose: Implements the Operator-Split Coupling Algorithm
        Args:
            cycle_no = current cycle number
            cycles = Max Cycles
        Notes:
            Algorithm works as follows
            -> Execute General coupling step 
                with HyKiCT run as a normal
                fluid code. 
            -> Send end state to Orig F_init 
            -> Send coupling params to 
                original F_init
            -> Run HyKiCT in Coupled mode 
            -> Finish 
        """
        self.fluid_obj.setNoCoupleSwitches()
        no_copy = False
        if cycle_no == cycles - 1: 
            no_copy = True

        (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
        fluid_Z, _, fluid_mass, fluid_specific_heat, sim_time) = self.fluidStep(self.io_obj.intermediate_fluid_outpath,
                                                            self.io_obj.fluid_input_path, False
                                                            )
        
        #Operator-split needs VFP-heat flow from original profile thus initialise kinetic 
        #with initial values
        (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
    fluid_Z, _, fluid_mass, sim_time)  = self.fluid_obj.returnInitValues()
        conv_heat_flow = self.returnConvHeatFlow(fluid_x_grid,fluid_x_centered_grid, 
                                fluid_Te,fluid_ne , fluid_Z,fluid_mass)

        if self.init.yaml_file['Mode']['Limit_density']:
            self.critical_density = 10 * (1114326918632954.5 / pow(self.fluid_obj.init.yaml_file['LaserParams']['Wavelength'], 2)) #Hard-coded limit to 10*nc
            if any(fluid_ne >= self.critical_density):
                self.laser_dir = self.fluid_obj.laser_direction
                self.couple_obj.limit_density = True
            else:
                self.critical_density = None
                self.laser_dir = None

        vfp_heat = self.kineticStep(cycle_no, conv_heat_flow,
                    fluid_x_grid, fluid_x_centered_grid, 
                            fluid_Te, fluid_ne, fluid_Z)

        (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
    fluid_Z, _, fluid_mass, sim_time, _)  = self.fluid_obj.getLastStepQuants()
        self.logger.info("Operator-Split Set HFCT tools")
        conv_heat_flow = self.returnConvHeatFlow(fluid_x_grid,fluid_x_centered_grid, 
                                fluid_Te,fluid_ne , fluid_Z,fluid_mass)
                                
        self.couple_obj.method(conv_heat_flow, vfp_heat, 
                        laser_dir = self.laser_dir, mass = fluid_mass, cell_wall_coord = fluid_x_grid,
                        q_snb = self.hfct_obj.q_snb)
            
        self.logger.info(self.coupling_message)
        self.couple_obj.setCoupleParams(self.io_obj.fluid_input_path, fluid_yaml = self.fluid_obj.init.yaml_file)
        self.fluid_obj.setOperatorSplitSwitch()
        _ = self.fluidStep(self.io_obj.fluid_output_path, self.io_obj.next_fluid_input_path, no_copy = no_copy)

    def main(self):
        RUN_PATH = os.path.join(self.init.yaml_file['Paths']['Base_dir'],
                     self.init.yaml_file['Paths']['Run_name'])        
        # ch = logging.StreamHandler()
        # ch.setLevel(logging.DEBUG)
        if not self.init.yaml_file['Misc']['HPC']:
            self.startprint(self.init.yaml_file["Codes"]["Fluid_code"],
                             self.init.yaml_file["Codes"]["Kinetic_code"])

        #Max Cycles
        cycles = self.init.yaml_file['Coupling_params']['Cycles']
        
        #Admin
        self.continue_step_path = os.path.join(RUN_PATH, 'CONTINUE_STEP.txt')
        overwrite = self.init.yaml_file['Misc']['Overwrite']
        start_cycle = 0
        self.kinetic_time_path = os.path.join(RUN_PATH, 'Kinetic_CPU_time.txt')
        self.fluid_time_path = os.path.join(RUN_PATH, 'Fluid_CPU_time.txt')
        self.cycle_time_path = os.path.join(RUN_PATH, 'Cycle_CPU_time.txt')


        if self.init.yaml_file['Misc']['Continue']:
            if os.path.exists(self.continue_step_path):
                start_cycle = np.loadtxt(self.continue_step_path, dtype=np.int) + 1
                overwrite = False
                self.init.yaml_file['Mode']['Start_from_kinetic'] = False

        self.hfct_obj = util.HeatFlowCouplingTools()
        self.time_stepper = util.TimeStepper(self.init.yaml_file['Coupling_params']['t_cycle'])
        self.io_obj = util.IO(
                        self.init.yaml_file['Paths']['Run_name'],
                        self.init.yaml_file['Paths']['Base_dir'], 
                        self.init.yaml_file['Paths']['K_src_dir'],
                        self.init.yaml_file['Paths']['F_src_dir'], 
                        self.init.yaml_file['Paths']['Init_path'],
                        start_cycle, cycles, overwrite)
        if self.init.yaml_file['Mode']['Couple_leap_frog'] or self.init.yaml_file['Mode']['Couple_operator_split']  :
            self.io_obj.intermediate_folder = True
        #Create Folders
        #Has implicit checks on whether to create folder again or not
        self.io_obj.createDirectoryOfOperation()

        fh = logging.FileHandler(os.path.join(RUN_PATH, 'Coupler.log'))
        fh.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        fh.setFormatter(formatter)
        root = logging.getLogger()
        root.addHandler(fh)
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        self.logger.addHandler(ch)
        atexit.register(self.cleanUpHandlers)

        #Note sizes are not allowed to diverge atm. 
        self.fluid_obj = HyKiCT(
                            self.io_obj._run_path, 
                            self.io_obj._f_src_dir, 
                            self.io_obj._f_init_path,
                            self.init.yaml_file['Paths']['F_config_path'],
                            nx = self.init.yaml_file['Coupling_params']['Nx'],
                            tmax = self.init.yaml_file['Coupling_params']['t_cycle'],
                            start_from_kin= self.init.yaml_file['Mode']['Start_from_kinetic'],
                            CoupleDivQ = self.init.yaml_file['Mode']['Couple_divq'],
                            CoupleMulti = self.init.yaml_file['Mode']['Couple_multi'],
                            CoupleSubtract = self.init.yaml_file['Mode']['Couple_subtract'],
                            )

        if self.fluid_obj.curr_time > 0 and not self.fluid_obj.init.yaml_file['Switches']['LoadInRadEnergy']:
           self.fluid_obj.init.yaml_file['Switches']['LoadInRadEnergy'] = True  
        #Adds the tmax to the init time. 
        if self.init.yaml_file['Misc']['Continue'] and start_cycle > 1:
            last_time = self.fluid_obj.getPhysicalRunTime(self.io_obj.all_fluid_output_path[start_cycle - 1])
            self.fluid_obj.updateTime(curr_time = last_time)
            self.fluid_obj.setTimes()

        if(self.init.yaml_file['Codes']['Kinetic_code'] == 'sol_kit'):
            self.kin_obj = SOL_KIT(
                            self.io_obj._fast_run_path,
                            self.io_obj._k_src_dir,
                            self.io_obj.kinetic_input_path,
                            self.io_obj.kinetic_output_path,
                            self.init.yaml_file['Paths']['K_config_path'],
                            nx = self.init.yaml_file['Coupling_params']['Nx'],
                            cycle = start_cycle,
                            convergence_monitoring = self.init.yaml_file['Misc']['Convergence_monitoring'],
                            cx1 = self.init.yaml_file['Misc']['HPC']
                            )
        if self.init.yaml_file['Mode']['Couple_divq']:
            self.couple_obj = DivQ()
            self.coupling_message = 'Calculating Div.Q'
        elif self.init.yaml_file['Mode']['Couple_multi']:
            self.couple_obj = Multiplier()
            self.coupling_message = 'Calculating Multiplier'
        elif self.init.yaml_file['Mode']['Couple_subtract']:
            self.couple_obj = Subtract()
            self.coupling_message = 'Calculating Subtract'
        else:
            raise Exception('No Valid Base Mode Chosen')
        if self.init.yaml_file['Mode']['no_negative_multiplier']:
            self.couple_obj.lower_limit = self.init.yaml_file['Coupling_params']['nlower_limit']
            self.couple_obj.upper_limit = self.init.yaml_file['Coupling_params']['nupper_limit']

        #Impact not ready
        # else:
            #self.kin_obj = IMPACT()
        if self.init.yaml_file['Misc']['Continue'] and start_cycle > 1:
            self.logger.info("CONTINUING FROM CYCLE")
            self.logger.info(start_cycle)
            cpu_time_path = os.path.join(RUN_PATH, 'Cycle_CPU_time.txt')
            fluid_time_path  = os.path.join(RUN_PATH, 'Fluid_CPU_time.txt')
            kin_time_path  = os.path.join(RUN_PATH, 'Kinetic_CPU_time.txt')
            if os.path.exists(cpu_time_path):
                self.cycle_time_taken = np.loadtxt(cpu_time_path).tolist()
                if isinstance(self.cycle_time_taken, float):
                    self.cycle_time_taken = [self.cycle_time_taken]
            if os.path.exists(fluid_time_path):
                self.fluid_time_taken = np.loadtxt(fluid_time_path).tolist()
                if isinstance(self.fluid_time_taken, float):
                    self.fluid_time_taken = [self.fluid_time_taken]
            if os.path.exists(cpu_time_path):
                self.kinetic_time_taken = np.loadtxt(kin_time_path).tolist()
                if isinstance(self.kinetic_time_taken, float):
                    self.kinetic_time_taken = [self.kinetic_time_taken]

        self.logger.info("Initial Conditions")
        self.logger.info("Run path {}".format(RUN_PATH))
        self.logger.info("Fast Tmp Run path {}".format(self.io_obj._fast_run_path))
        self.logger.info("Start Cycle {}".format(start_cycle)) 
        self.logger.info("Max Cycle {}".format(cycles)) 
        self.logger.info("Overwrite {}".format(overwrite))
        self.logger.info("Number of Process for KINETIC {}".format(self.kin_obj.init.yaml_file['Params']['Np'])) 
        self.logger.info("Div q couple switch is {}".format(self.fluid_obj.init.yaml_file['Switches']['CoupleDivQ']))
        self.logger.info("Multi couple switch is {}".format(self.fluid_obj.init.yaml_file['Switches']['CoupleMulti']))
        self.logger.info("Subtract couple switch is {}".format(self.fluid_obj.init.yaml_file['Switches']['CoupleSubtract']))
        self.logger.info("Leap-Frog couple switch is {}".format(self.init.yaml_file['Mode']['Couple_leap_frog']))
        self.logger.info("Operator-Split couple switch is {}".format(self.init.yaml_file['Mode']['Couple_operator_split']))
        self.logger.info("Steps is {}".format(self.fluid_obj.init.yaml_file['TimeParameters']['steps']))
        self.logger.info("Tmax is {}".format(self.fluid_obj.init.yaml_file['TimeParameters']['t_max']))


        self.logger.info("Starting couple LOOP")
        self.first_pass = True
        run_only_fluid = False
        np.savetxt(os.path.join(RUN_PATH, "status.txt"), np.array([0], dtype=np.int), fmt = '%1.1i')

        for cycle_no in range(start_cycle, cycles, 1):

            t0 = time.time()
            self.prettyPrint(' Initiating CYCLE ' + str(cycle_no)) 
            self.logger.info("Initiating Cycle {}".format(cycle_no))
            #Update Paths
            self.io_obj.cycle_counter = cycle_no
            self.io_obj.nextCyclePathManager()
            self.kin_obj.cycle_dump_path = self.io_obj.cycle_dump_path
            self.fluid_obj.cycle_dump_path = self.io_obj.cycle_dump_path
            self.logger.info("Cycle dump path is {}".format(self.io_obj.cycle_dump_path))
            if cycle_no == 0:
                self.logger.info("Copying init file content")
                copy_tree(self.init.yaml_file['Paths']['Init_path'], self.io_obj.fluid_input_path)
                np.savetxt(os.path.join(RUN_PATH, 'NO_CYCLES.txt'), np.array([cycles - 1]), fmt = '%i')

            if self.init.yaml_file['Mode']['Couple_leap_frog']:
                self.leapFrog(cycle_no, cycles)
            elif self.init.yaml_file['Mode']['Couple_operator_split']:
                self.operatorSplit(cycle_no, cycles)
            else:
                if (cycle_no == 1 and self.init.yaml_file['Mode']['Start_from_kinetic']):
                    self.fluid_obj.revertStartKinSwitches()
                if (cycle_no == cycles - 1):
                    run_only_fluid = True
                self.generalCoupling(cycle_no, cycles,
                            run_only_fluid = run_only_fluid)

            if (self.init.yaml_file['Misc']['HDF5'] and 
                cycle_no % self.init.yaml_file['Coupling_params']['hdf5_output_freq'] == 0
                or cycle_no == cycles - 1):
                self.storeToHdf5(cycle_no)

            self.logger.info('CYCLE COMPLETED UPDATING CONTINUE STEP')
            np.savetxt(self.continue_step_path, np.array([cycle_no]), fmt = '%i')

            t1 = time.time()
            self.cycle_time_taken.append(t1 - t0)
            self.logger.info('CPU TIME FOR CYCLE {} IS {} '.format(cycle_no, t1-t0))
            self.first_pass = False

            np.savetxt(self.kinetic_time_path, np.array([self.kinetic_time_taken]))
            np.savetxt(self.fluid_time_path, np.array([self.fluid_time_taken]))
            np.savetxt(self. cycle_time_path, np.array([self.cycle_time_taken]))

        if self.init.yaml_file['Misc']['HDF5']:
            self.logger.info("Delete All Folders")
            self.io_obj.deleteAll()    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Coupling of Hydrodynamic and Kinetic code. The two codes are specified in the input file.')
    parser.add_argument('-p', '--path', required=True, help = 'Give path to Input yml file.')
    args = vars(parser.parse_args())
    couple = Coupler(args['path'])
    couple.main()