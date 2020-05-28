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
import signal
import shutil
import sys
import utils as util
from couple_sol_kit.sol_kit import SOL_KIT
from couple_hykict.hykict import HyKiCT
class Coupler:
    """
    main class i.e. where all objects are passed into and the coupling is run.
    """
    def __init__(self, init_file_path):
        self.yml_init_file_path = init_file_path
        self.pre_heat_present = False
        atexit.register(self.cleanUpHandlers)
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
        # self.logger.addHandler(fh)
        self.logger.addHandler(ch)
        if self.init.yaml_file['Misc']['HDF5']:
            io_obj._createHDF5()
        
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
        self.logger.info("Initial Conditions")
        self.logger.info("Run path {}".format(RUN_PATH))
        self.logger.info("Start Cycle {}".format(start_cycle)) 
        self.logger.info("Max Cycle {}".format(cycles)) 
        self.logger.info("Overwrite {}".format(overwrite)) 
        #Enforce equal size ... Constrain at the moment
        kin_obj.init.yaml_file['Params']['Nx'] = self.init.yaml_file['Coupling_params']['Nx']
        fluid_obj.init.yaml_file['FixedParameters']['nx'] = self.init.yaml_file['Coupling_params']['Nx']
        kin_obj.nx = kin_obj.init.yaml_file['Params']['Nx'] 

        fluid_obj.init.yaml_file['Switches']['CoupleDivQ'] = self.init.yaml_file['Coupling_params']['Couple_divq']
        fluid_obj.init.yaml_file['Switches']['CoupleMulti'] = self.init.yaml_file['Coupling_params']['Couple_multi']
        #Store original fluid yaml
        self.original_f_init = copy.deepcopy(fluid_obj.init.yaml_file)
        #modify original to represent run mode

        #initially fluid should always be run in no couple mode.
        #Has checks if continue mode is engaged
        if(start_cycle == 0 and not self.init.yaml_file['Coupling_params']['Start_coupled']):
            self.logger.info("Starting in MODE: No-Coupling")

            fluid_obj.init.yaml_file['Switches']['CoupleDivQ'] = False
            fluid_obj.init.yaml_file['Switches']['CoupleMulti'] = False

            if self.init.yaml_file['Coupling_params']['Start_from_kinetic']:
                self.logger.info("Starting in MODE: Start From Kinetic")

                fluid_obj.init.yaml_file['TimeParameters']['steps'] = 0
                fluid_obj.init.yaml_file['TimeParameters']['t_max'] = 0

                self.logger.info("Steps is {}".format(fluid_obj.init.yaml_file['TimeParameters']['steps']))
                self.logger.info("Tmax is {}".format(fluid_obj.init.yaml_file['TimeParameters']['t_max']))
        
        if (self.init.yaml_file['Coupling_params']['Start_coupled'] and
             self.init.yaml_file['Coupling_params']['Couple_adaptive']):

            self.logger.info("Starting Coupled with Adaptive Switching")

            fluid_obj.init.yaml_file['Switches']['CoupleDivQ'] = True
            fluid_obj.init.yaml_file['Switches']['CoupleMulti'] = False
            fluid_obj.init.yaml_file['TimeParameters']['steps'] = 0

            self.logger.info("Fluid Object initial conditions in MODE: START COUPLED WITH ADAPTIVE SWITCH")
            self.logger.info("Div q couple switch is {}".format(fluid_obj.init.yaml_file['Switches']['CoupleDivQ']))
            self.logger.info("Multi couple switch is {}".format(fluid_obj.init.yaml_file['Switches']['CoupleMulti']))
            self.logger.info("Steps is {}".format(fluid_obj.init.yaml_file['TimeParameters']['steps']))
            self.logger.info("Tmax is {}".format(fluid_obj.init.yaml_file['TimeParameters']['t_max']))

        #Coupling loop
        self.logger.info("Starting couple LOOP")
        for cycle_no in range(start_cycle, cycles, 1):
            self.prettyPrint(' RUNNING CYCLE ' + str(cycle_no)) 
            self.logger.info("Running Cycle {}".format(cycle_no))
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
                np.savetxt(os.path.join(RUN_PATH, 'NO_CYCLES.txt'), np.array([cycles - 1]), fmt = '%i' )       
        
            if cycle_no >= 1:
                #Engage coupling 
                self.logger.info("Engage Coupling if MODE: NO-COUPLING/ Start from Kinetic")
                if(self.init.yaml_file['Coupling_params']['Couple_adaptive']):
                    if self.pre_heat_present or self.init.yaml_file['Coupling_params']['Start_coupled']:
                        self.init.yaml_file['Coupling_params']['Start_coupled'] = False
                        k_physical_time = kin_obj.getPhysicalRunTime()
                        f_run_time = k_physical_time * self.init.yaml_file['Coupling_params']['Eta']
                        if(f_run_time < fluid_obj.init.yaml_file['TimeParameters']['dt']):
                            f_dt = fluid_obj.init.yaml_file['TimeParameters']['dt'] / 100 #100 here an arbitary choice
                        else:
                            f_dt = fluid_obj.init.yaml_file['TimeParameters']['dt']
                        
                        fluid_obj.init.yaml_file['TimeParameters']['steps'] = 0
                        fluid_obj.init.yaml_file['TimeParameters']['t_max'] = f_run_time
                        fluid_obj.init.yaml_file['TimeParameters']['dt'] = f_dt
                        fluid_obj.init.yaml_file['Switches']['CoupleDivQ'] = True 
                        fluid_obj.init.yaml_file['Switches']['CoupleMulti'] = False 
                    else:
                        fluid_obj.init.yaml_file = self.original_f_init
                        fluid_obj.init.yaml_file['Switches']['CoupleDivQ'] = False 
                        fluid_obj.init.yaml_file['Switches']['CoupleMulti'] = True 

                if self.init.yaml_file['Coupling_params']['Start_from_kinetic']:
                    if(self.init.yaml_file['Coupling_params']['Couple_adaptive']):
                        #Configurations already done above
                        pass
                    else:
                        fluid_obj.init.yaml_file = self.original_f_init
          
            ###########
            #Fluid Step
            ###########
            #In this example running HyKiCT
            self.prettyPrint(' RUNNING ' + self.init.yaml_file['Codes']['Fluid_code'], color = True)
            self.logger.info("Setting Fluid Parameters")
            #Set Paths that change
            fluid_obj.init.yaml_file['Paths']['Init_Path'] = io_obj.fluid_input_path
            fluid_obj.init.yaml_file['Paths']['Out_Path'] = io_obj.fluid_output_path
            fluid_obj._cycle_dump_path = io_obj.cycle_dump_path
            fluid_obj._fluid_output_path = io_obj.fluid_output_path
            fluid_obj._fluid_input_path = io_obj.fluid_input_path
            self.logger.info("Fluid input path {}".format(fluid_obj.init.yaml_file['Paths']['Init_Path']))
            self.logger.info("Fluid output path {}".format(fluid_obj.init.yaml_file['Paths']['Out_Path']))
            ################
            #Any other parameter updates that needs to be set can be added here.
            #################
            fluid_obj.setFiles()
            self.logger.info("Start Fluid")
            fluid_obj.Run()
            self.logger.info("End Fluid")

            #Breaks here ... Last cycle allowed to run hydro step
            if cycle_no == cycles - 1:
                if self.init.yaml_file['Misc']['Zip']:
                    io_obj.zipAndDelete()        
                np.savetxt(continue_step_path, np.array([cycle_no]), fmt = '%i')
                self.logger.info("End Coupling")
                break 

            self.logger.info("Get Last Fluid Quants")
            (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
            fluid_Z, _, fluid_mass) = fluid_obj.getLastStepQuants()
         
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
            hfct_obj.lambda_ei(hfct_obj.electron_temperature, 
                                hfct_obj.electron_temperature,
                                hfct_obj.zbar)
            self.logger.info("HFCT Spitzer Calculation")
            hfct_obj.spitzerHarmHeatFlow()
         
            #############
            #Kinetic Step
            #############
            #Set any parameter ifles
            #self.init all load in files form hydro
            #RUN
            #Move files if neccessary
            self.prettyPrint(' RUNNING ' + self.init.yaml_file['Codes']['Kinetic_code'], color = True)
            self.logger.info("Setting Kinetic Parameters")
            if cycle_no == 0:
                #input and output unchanged 
                self.logger.info("Set files... These are UNCHANGED EXCEPT IF LOAD_f1")
                kin_obj.setFiles()

            kin_obj._kinetic_input_path = io_obj.kinetic_input_path
            kin_obj._kinetic_output_path = io_obj.kinetic_output_path
            kin_obj._cycle_dump_path = io_obj.cycle_dump_path
            self.logger.info("Kinetic Paths")
            self.logger.info("Kinetic input path {}".format(io_obj.kinetic_input_path))
            self.logger.info("Kinetic output path {}".format(io_obj.kinetic_output_path))

            if self.init.yaml_file['Coupling_params']['Load_f1']:
                if cycle_no > 0:
                    self.logger.info("Engaging MODE: Load F1")
                    kin_obj.previous_kinetic_output_path = io_obj.preserved_kinetic_output_path[-2] 
                    kin_obj.load_f1 = True
                if self.init.yaml_file['Coupling_params']['Start_coupled']:
                    self.logger.info("Engaging MODE: Load F1 from Start Coupled")
                    #Hack to continue broken sims, assumes RESTART and GRID exists in init folder
                    if os.path.exists(os.path.join(self.init.yaml_file['Paths']['Init_path'], "kinetic_restart")):
                        kin_obj.load_f1 = True
                        kin_obj.previous_kinetic_output_path = os.path.join(self.init.yaml_file['Paths']['Init_path'], "kinetic_restart")

            self.logger.info("Initialize Kinetic")
            kin_obj.initFromHydro(fluid_x_grid, fluid_x_centered_grid, 
                                fluid_Te, fluid_ne, fluid_Z)
            kin_obj.sh_heat_flow = hfct_obj.spitzer_harm_heat
            self.logger.info("Start Kinetic")
            kin_obj.Run()
            self.logger.info("End Kinetic")
            #Set hfct to contain the vfp heat flow to do the neccessary coupling calcs.
            self.logger.info("Fetch Last Heat Flow")
            hfct_obj.vfp_heat= kin_obj.getLastHeatFlow()
            self.logger.info("Move files Kinetic")
            kin_obj.moveFiles()
          
            ##############
            #Coupling Step
            ##############
            #the _ are unused variables fluid_v and fluid_laser, for future
            #use these may be required and thus, left in. 

            self.logger.info("Coupling Heat-Flow")
            if self.init.yaml_file['Coupling_params']['Couple_divq']:
                #qe here is div.q_vfp 
                self.logger.info("Calculate Div.q")
                qe = hfct_obj.divQHeatFlow()
            
            elif (self.init.yaml_file['Coupling_params']['Couple_multi'] or 
                    self.init.yaml_file['Coupling_params']['Couple_adaptive']):
                #qe here is q_vfp/q_sh
                self.logger.info("Calculate Multipliers")
                (qe, pre_heat_start_index, pre_heat_last_index,
                pre_heat_fit_params, front_heat_start_index, 
                front_heat_last_index, front_heat_fit_params)  = hfct_obj.multiplier()
                #If pre heat or front heat is present, in adapativbe coupling
                #We utilise div.q coupling to get rid of these fronts.
                #Otherwise, apply the exponential models. 
                if(self.init.yaml_file['Coupling_params']['Couple_adaptive']):
                    if(pre_heat_start_index > 0 or front_heat_start_index > 0):
                        self.logger.info("Pre-heat present in MODE: Adapative ignore multipliers")
                        self.pre_heat_present = True
                        qe = hfct_obj.divQHeatFlow()
                        pre_heat_fit_params = None
                        front_heat_fit_params = None
                else:
                    self.logger.info("Pre-heat Present in MODE: Multipliers")
                    #Modify yaml for future write
                    fluid_obj.init.yaml_file['FixedParameters']['Preheat_StartIndex'] = pre_heat_start_index.item()
                    fluid_obj.init.yaml_file['FixedParameters']['Preheat_LastIndex'] = pre_heat_last_index.item()
                    fluid_obj.init.yaml_file['FixedParameters']['Frontheat_StartIndex'] = front_heat_start_index.item()
                    fluid_obj.init.yaml_file['FixedParameters']['Frontheat_LastIndex'] = front_heat_last_index.item()
            
            #Finish by fluid init next set of files 
            self.logger.info("Init Fluid")
            fluid_obj.initHydroFromKinetic(io_obj.next_fluid_input_path, qe,
                                            pre_heat_fit_params, front_heat_fit_params)
            if self.init.yaml_file['Misc']['HDF5']:
                self.logger.info("Store to HDF5")
                fluid_obj.storeToHdf5(io_obj.hdf5_file, cycle_no)
                kin_obj.storeToHdf5(io_obj.hdf5_file, cycle_no)
            if self.init.yaml_file['Misc']['Zip']:
                self.logger.info("Store to zip")
                io_obj.zipAndDelete()        
            #update continue file
            np.savetxt(continue_step_path, np.array([cycle_no]), fmt = '%i')
            
        if self.init.yaml_file['Misc']['HDF5']:
            self.logger.info("Delete All Folders")
            io_obj.deleteAll()    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Coupling of Hydrodynamic and Kinetic code. The two codes are specified in the input file.')
    parser.add_argument('-p', '--path', required=True, help = 'Give path to Input yml file.')
    args = vars(parser.parse_args())
    couple = Coupler(args['path'])
    couple.main()








