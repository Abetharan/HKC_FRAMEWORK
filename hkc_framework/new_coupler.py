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

    # def multiplier(self):

    # def divq(self): 

    # def subtract(self): 
    # def leapFrog(self):

    # def operator_split(self):

    # def adapative(self): 



    def main(self):
        self.init = util.Input(self.yml_init_file_path)
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
        continue_step_path = os.path.join(RUN_PATH, 'CONTINUE_STEP.txt')
        overwrite = self.init.yaml_file['Misc']['Overwrite']
        pre_heat_fit_params = None
        front_heat_fit_params = None
        start_cycle = 0
        kinetic_time_path = os.path.join(RUN_PATH, 'Kinetic_CPU_time.txt')
        fluid_time_path = os.path.join(RUN_PATH, 'Fluid_CPU_time.txt')
        cycle_time_path = os.path.join(RUN_PATH, 'Cycle_CPU_time.txt')


        if self.init.yaml_file['Misc']['Continue']:
            #Continue in this framework is designed such that
            #a run does not have to fail (for whatever reason) and then continue
            #it will keep going until its objective has been reached which is
            #all cycles have been run through. 
            #Smarter continuing should be implemented such that errors due to
            #codes crashing should also crash the framework. Wheras,
            #HPC facility failures should just make it continue. 
            if os.path.exists(continue_step_path):
                start_cycle = np.loadtxt(continue_step_path, dtype=np.int) + 1
                overwrite = False

        #Create Objects 
        hfct_obj = util.HeatFlowCouplingTools()
        io_obj = util.IO(
                        self.init.yaml_file['Paths']['Run_name'],
                        self.init.yaml_file['Paths']['Base_dir'], 
                        self.init.yaml_file['Paths']['K_src_dir'],
                        self.init.yaml_file['Paths']['F_src_dir'], 
                        self.init.yaml_file['Paths']['Init_path'],
                        start_cycle, cycles, overwrite)
        if self.init.yaml_file['Mode']['Couple_leap_frog'] or self.init.yaml_file['Mode']['Couple_operator_split']  :
            io_obj.intermediate_folder = True
        #Create Folders
        #Has implicit checks on whether to create folder again or not
        io_obj.createDirectoryOfOperation()

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

        
        fluid_obj = HyKiCT(
                            io_obj._run_path, 
                            io_obj._f_src_dir, 
                            io_obj._f_init_path,
                            self.init.yaml_file['Paths']['F_config_path'])
                         
        if(self.init.yaml_file['Codes']['Kinetic_code'] == 'sol_kit'):
            kin_obj = SOL_KIT(
                            io_obj._fast_run_path,
                            io_obj._k_src_dir,
                            io_obj.kinetic_input_path,
                            io_obj.kinetic_output_path,
                            self.init.yaml_file['Paths']['K_config_path'],
                            self.init.yaml_file['Misc']['Convergence_monitoring'],
                            self.init.yaml_file['Misc']['HPC'])
        #Impact not ready
        # else:
            #kin_obj = IMPACT()
        if self.init.yaml_file['Misc']['Continue'] and start_cycle > 1:
                self.logger.info("CONTINUING FROM CYCLE")
                self.logger.info(start_cycle)
        self.logger.info("Initial Conditions")
        self.logger.info("Run path {}".format(RUN_PATH))
        self.logger.info("Fast Tmp Run path {}".format(io_obj._fast_run_path))
        self.logger.info("Start Cycle {}".format(start_cycle)) 
        self.logger.info("Max Cycle {}".format(cycles)) 
        self.logger.info("Overwrite {}".format(overwrite))
        self.logger.info("Number of Process for KINETIC {}".format(kin_obj.init.yaml_file['Params']['Np'])) 

        #Enforce equal size ... Constrain at the moment
        kin_obj.init.yaml_file['Params']['Nx'] = self.init.yaml_file['Coupling_params']['Nx']
        fluid_obj.init.yaml_file['FixedParameters']['nx'] = self.init.yaml_file['Coupling_params']['Nx']
        kin_obj.nx = kin_obj.init.yaml_file['Params']['Nx'] 

        if self.init.yaml_file['Mode']['Couple_divq']:
            couple_obj = DivQ(fluid_obj.init.yaml_file, self.init.yaml_file['Mode']['Start_from_kinetic'])
        elif self.init.yaml_file['Mode']['Couple_multi']:
            couple_obj = Multiplier(fluid_obj.init.yaml_file, self.init.yaml_file['Mode']['Start_from_kinetic'])
        elif self.init.yaml_file['Mode']['Couple_subtract']:
            couple_obj = Subtract(fluid_obj.init.yaml_file, self.init.yaml_file['Mode']['Start_from_kinetic'])
        else:
            raise Exception('No Valid Mode Chosen')

        ##############
        ###START LOOP
        #############
        self.logger.info("Starting couple LOOP")
        self.first_pass = True
        np.savetxt(os.path.join(RUN_PATH, "status.txt"), np.array([0], dtype=np.int), fmt = '%1.1i')
        for cycle_no in range(start_cycle, cycles, 1):

            t0 = time.time()
            self.prettyPrint(' Initiating CYCLE ' + str(cycle_no)) 
            self.logger.info("Initiating Cycle {}".format(cycle_no))
            #Update Paths
            io_obj.cycle_counter = cycle_no
            io_obj.nextCyclePathManager()
            kin_obj.cycle_dump_path = io_obj.cycle_dump_path
            fluid_obj.cycle_dump_path = io_obj.cycle_dump_path
            self.logger.info("Cycle dump path is {}".format(io_obj.cycle_dump_path))
            #For restarting/continue
            #The copying could in principal be removed. 
            #Output directory just self-consistent this way.
            if cycle_no == 0:
                self.logger.info("Copying init file content")
                copy_tree(self.init.yaml_file['Paths']['Init_path'], io_obj.fluid_input_path)
                np.savetxt(os.path.join(RUN_PATH, 'NO_CYCLES.txt'), np.array([cycles - 1]), fmt = '%i')       


            self.prettyPrint(' RUNNING ' + self.init.yaml_file['Codes']['Fluid_code'], color = True)
            self.logger.info("Setting Fluid Parameters")
            tfluid_start = time.time()
            #Set Paths that change
            fluid_obj.init.yaml_file['Paths']['Init_Path'] = io_obj.fluid_input_path
            fluid_obj.init.yaml_file['Paths']['Laser_Profile_Path'] = os.path.join(io_obj.fluid_input_path, "laser_profile.txt")
            
            fluid_obj._cycle_dump_path = io_obj.cycle_dump_path
            fluid_obj._fluid_input_path = io_obj.fluid_input_path
            
            #First pass in Leap frog ... run in default mode with specified flux limiter and 
            #save to a different folder 
            if self.init.yaml_file['Mode']['Couple_leap_frog'] or self.init.yaml_file['Mode']['Couple_operator_split']:
                self.logger.info("Setting Intermediate Fluid Output Path")
                
                if self.init.yaml_file['Mode']['Couple_operator_split']:
                    fluid_obj.init.yaml_file['TimeParameters']['t_max'] = self.init.yaml_file['Coupling_params']['t_fluid_operator'] # REF THIS IS THE TIME MODIFICATION TO OPERATOR-SPLIT
                fluid_obj.init.yaml_file['Paths']['Out_Path'] = io_obj.intermediate_fluid_outpath 
                fluid_obj._fluid_output_path = io_obj.intermediate_fluid_outpath
            else:
                fluid_obj.init.yaml_file['Paths']['Out_Path'] = io_obj.fluid_output_path
                fluid_obj._fluid_output_path = io_obj.fluid_output_path
            
            # if cycle_no >= 1:

            self.logger.info("Fluid input path {}".format(fluid_obj.init.yaml_file['Paths']['Init_Path']))
            self.logger.info("Fluid output path {}".format(fluid_obj.init.yaml_file['Paths']['Out_Path']))

            ################
            #Any other parameter updates that needs to be set can be added here.
            #################
            fluid_obj.setFiles()
            self.logger.info("Start Fluid")
            fluid_obj.Run()
            self.logger.info("End Fluid")

            self.logger.info("Get Last Fluid Quants")
            (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
            fluid_Z, _, fluid_mass, sim_time) = fluid_obj.getLastStepQuants()
         
            tfluid_end = time.time()
            self.fluid_time_taken.append(tfluid_end - tfluid_start)

            ####################
            #Init coupling tools
            ####################
            #self.init heat_flow_tools here 
            self.logger.info("Set HFCT tools")
            hfct_obj.electron_temperature = fluid_Te
            hfct_obj.electron_number_density = fluid_ne
            hfct_obj.zbar = fluid_Z
            hfct_obj.cell_wall_coord = fluid_x_grid
            hfct_obj.cell_centered_coord = fluid_x_centered_grid
            hfct_obj.mass = fluid_mass
            #Calculate spitzer harm from last step fluid quants
            hfct_obj.lambda_ei(hfct_obj.electron_temperature * (BOLTZMANN_CONSTANT/ELEMENTARY_CHARGE), 
                                hfct_obj.electron_number_density,
                                hfct_obj.zbar)
            self.logger.info("HFCT Spitzer Calculation")
            hfct_obj.spitzerHarmHeatFlow()
            if fluid_obj.init.yaml_file['Switches']['SNBHeatFlow']:
                hfct_obj.snb = True
                hfct_obj.snb_heat_flow(fluid_obj.init.yaml_file['FixedParameters']['ng'],fluid_obj.init.yaml_file['FixedParameters']['MaxE'], 2)
         
            #############
            #Kinetic Step
            #############
            #Set any parameter ifles
            #self.init all load in files form hydro
            #RUN
            #Move files if neccessary

            tkin_start = time.time()
            self.prettyPrint(' RUNNING ' + self.init.yaml_file['Codes']['Kinetic_code'], color = True)
            self.logger.info("Setting Kinetic Parameters")
            if cycle_no == 0 or self.first_pass:
                #input and output unchanged 
                self.logger.info("Set files... These are UNCHANGED EXCEPT IF LOAD_f1")
                kin_obj.setFiles()

            kin_obj._kinetic_input_path = io_obj.kinetic_input_path
            kin_obj._kinetic_output_path = io_obj.kinetic_output_path
            kin_obj._cycle_dump_path = io_obj.cycle_dump_path
            self.logger.info("Kinetic Paths")
            self.logger.info("Kinetic input path {}".format(io_obj.kinetic_input_path))
            self.logger.info("Kinetic output path {}".format(io_obj.kinetic_output_path))

            if self.init.yaml_file['Mode']['Load_f1']:
                if cycle_no > 0:
                    self.logger.info("Engaging MODE: Load F1")
                    #all_kinetic_output_paths contain all paths for kinetic output, require last cycle for load F1, thus cycle_no - 1
                    kin_obj.previous_kinetic_output_path = io_obj.all_kinetic_output_path[cycle_no - 1] 
                    kin_obj.load_f1 = True                                                              
            self.logger.info("Initialize Kinetic")
            if not fluid_obj.init.yaml_file['Switches']['SNBHeatFlow']:
                kin_obj.sh_heat_flow = hfct_obj.spitzer_harm_heat
            else:
                kin_obj.sh_heat_flow = hfct_obj.q_snb

            if self.init.yaml_file['Mode']['Limit_density']:
                critical_density = 10 * (1114326918632954.5 / pow(fluid_obj.init.yaml_file['LaserParams']['Wavelength'], 2)) #Hard-coded limit to 10*nc
                laser_dir = "left"#fluid_obj.laser_direction
                # laser_dir = fluid_obj.laser_direction
            else:
                critical_density = None
                laser_dir = None

            kin_obj.initFromHydro(fluid_x_grid, fluid_x_centered_grid, 
                                fluid_Te, fluid_ne, fluid_Z, critical_density = critical_density, laser_dir = laser_dir)
            
            self.logger.info("Start Kinetic")
            kin_obj.Run()
            self.logger.info("End Kinetic")
            #Set hfct to contain the vfp heat flow to do the neccessary coupling calcs.
            self.logger.info("Fetch Last Heat Flow")
            hfct_obj.vfp_heat= kin_obj.getLastHeatFlow()
            self.logger.info("Move files Kinetic")
            kin_obj.moveFiles()
            tkin_end = time.time()
            self.kinetic_time_taken.append(tkin_end - tkin_start)
















if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Coupling of Hydrodynamic and Kinetic code. The two codes are specified in the input file.')
    parser.add_argument('-p', '--path', required=True, help = 'Give path to Input yml file.')
    args = vars(parser.parse_args())
    couple = Coupler(args['path'])
    couple.main()