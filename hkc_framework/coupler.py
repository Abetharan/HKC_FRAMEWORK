""" 
Hydro-Kinetic Coupling.
Heat-Flow coupling between Hydro and Kinetic codes. 
@author = Abetharan Antony
"""
import argparse
from distutils.dir_util import copy_tree
import math
import numpy as np
from PIL import Image, ImageFont, ImageDraw
import os
import signal
import shutil
import utils as util
from couple_sol_kit.sol_kit import SOL_KIT
from couple_hykict.hykict import HyKiCT


class Coupler:
    """
    main class i.e. where all objects are passed into and the coupling is run.
    """
    def __init__(self, init_file_path):
        self.yml_init_file_path = init_file_path
        
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
    
    def main(self):
        self.init = util.Input(self.yml_init_file_path)
        RUN_PATH = os.path.join(self.init.yaml_file['Paths']['Base_dir'],
                     self.init.yaml_file['Paths']['Run_name'])        

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
        
        fluid_obj = HyKiCT(
                            io_obj._run_path, 
                            io_obj._f_src_dir, 
                            io_obj._f_init_path,
                            self.init.yaml_file['Paths']['F_config_path'])
                         
        if(self.init.yaml_file['Codes']['Kinetic_code'] == 'sol_kit'):
            kin_obj = SOL_KIT(
                            io_obj._run_path,
                            io_obj._k_src_dir,
                            io_obj.kinetic_input_path,
                            io_obj.kinetic_output_path,
                            self.init.yaml_file['Paths']['K_config_path'],
                            self.init.yaml_file['Misc']['Convergence_monitoring'],
                            self.init.yaml_file['Misc']['HPC'])
        #Impact not ready
        # else:
            #kin_obj = IMPACT()
        
        #Enforce equal size ... Constrain at the moment
        kin_obj.init.yaml_file['Params']['Nx'] = self.init.yaml_file['Coupling_params']['Nx']
        fluid_obj.init.yaml_file['FixedParameters']['nx'] = self.init.yaml_file['Coupling_params']['Nx']
        kin_obj.nx = kin_obj.init.yaml_file['Params']['Nx'] 

        #Store original fluid yaml
        self.original_f_init = fluid_obj.init.yaml_file
        #modify original to represent run mode
        self.original_f_init['Switches']['CoupleDivQ'] = self.init.yaml_file['Coupling_params']['Couple_divq']
        self.original_f_init['Switches']['CoupleMulti'] = self.init.yaml_file['Coupling_params']['Couple_multi']

        #initially fluid should always be run in no couple mode.
        #Has checks if continue mode is engaged
        if(start_cycle == 0):
            fluid_obj.init.yaml_file['Switches']['CoupleDivQ'] = False
            fluid_obj.init.yaml_file['Switches']['CoupleMulti'] = False
            if self.init.yaml_file['Coupling_params']['Start_from_kinetic']:
                fluid_obj.init.yaml_file['TimeParameters']['steps'] = 0
                fluid_obj.init.yaml_file['TimeParameters']['t_max'] = 0
        
        #Coupling loop
        for cycle_no in range(start_cycle, cycles, 1):
            self.prettyPrint(' RUNNING CYCLE ' + str(cycle_no)) 
            #Update Paths
            io_obj.cycle_counter = cycle_no
            io_obj.nextCyclePathManager()
            kin_obj.cycle_dump_path = io_obj.cycle_dump_path
            fluid_obj.cycle_dump_path = io_obj.cycle_dump_path

            #For restarting/continue
            #The copying could in principal be removed. 
            #Output directory just self-consistent this way.
            if cycle_no == 0:
                copy_tree(self.init.yaml_file['Paths']['Init_path'], io_obj.fluid_input_path)
                np.savetxt(os.path.join(RUN_PATH, 'NO_CYCLES.txt'), np.array([cycles - 1]), fmt = '%i' )       
        
            if cycle_no >= 1:
                #Update paths
                #Engage coupling 
                if(self.init.yaml_file['Coupling_params']['Couple_adaptive']):
                    if(pre_heat_start_index > 0 or front_heat_start_index > 0):
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

                if self.init.yaml_file['Coupling_params']['Start_from_kinetic']:
                    if(self.init.yaml_file['Coupling_params']['Couple_adaptive']):
                        pass
                    else:
                        fluid_obj.init.yaml_file = self.original_f_init
            ###########
            #Fluid Step
            ###########
            #In this example running HyKiCT
            self.prettyPrint(' RUNNING ' + self.init.yaml_file['Codes']['Fluid_code'], color = True)
            #Set Paths that change
            fluid_obj.init.yaml_file['Paths']['Init_Path'] = io_obj.fluid_input_path
            fluid_obj.init.yaml_file['Paths']['Out_Path'] = io_obj.fluid_output_path
            fluid_obj._cycle_dump_path = io_obj.cycle_dump_path
            fluid_obj._fluid_output_path = io_obj.fluid_output_path
            fluid_obj.init_file_path = io_obj.fluid_input_path
            ################
            #Any other parameter updates that needs to be set can be added here.
            #################
            fluid_obj.setFiles()
            fluid_obj.Run()

            #Breaks here ... Last cycle allowed to run hydro step
            if cycle_no == cycles - 1:
                if self.init.yaml_file['Misc']['Zip']:
                    io_obj.zipAndDelete()        
                break 

            (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
            fluid_Z, _, fluid_mass) = fluid_obj.getLastStepQuants()
            ####################
            #Init coupling tools
            ####################
            #self.init heat_flow_tools here 
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
            hfct_obj.spitzerHarmHeatFlow()
            #############
            #Kinetic Step
            #############
            #Set any parameter ifles
            #self.init all load in files form hydro
            #RUN
            #Move files if neccessary
            self.prettyPrint(' RUNNING ' + self.init.yaml_file['Codes']['Kinetic_code'], color = True)
            if cycle_no == 0:
                #input and output unchanged 
                kin_obj.setFiles()

            kin_obj._kinetic_input_path = io_obj.kinetic_input_path
            kin_obj._kinetic_output_path = io_obj.kinetic_output_path
            kin_obj._cycle_dump_path = io_obj.cycle_dump_path
            if self.init.yaml_file['Coupling_params']['Load_f1']:
                if cycle_no > 0:
                    kin_obj.previous_cycle_output_path = io_obj.preserved_kinetic_output_path[-2] 
                    kin_obj.load_f1 = True
            kin_obj.initFromHydro(fluid_x_grid, fluid_x_centered_grid, 
                                fluid_Te, fluid_ne, fluid_Z)
            kin_obj.sh_heat_flow = hfct_obj.spitzer_harm_heat
            kin_obj.Run()
            
            #Set hfct to contain the vfp heat flow to do the neccessary coupling calcs.
            hfct_obj.vfp_heat= kin_obj.getLastHeatFlow()
            kin_obj.moveFiles()
            ##############
            #Coupling Step
            ##############
            #the _ are unused variables fluid_v and fluid_laser, for future
            #use these may be required and thus, left in. 


            if self.init.yaml_file['Coupling_params']['Couple_divq']:
                #qe here is div.q_vfp 
                qe = hfct_obj.divQHeatFlow()
            
            elif (self.init.yaml_file['Coupling_params']['Couple_multi'] or 
                    self.init.yaml_file['Coupling_params']['Couple_adaptive']):
                #qe here is q_vfp/q_sh
                (qe, pre_heat_start_index, pre_heat_last_index,
                pre_heat_fit_params, front_heat_start_index, 
                front_heat_last_index, front_heat_fit_params)  = hfct_obj.multiplier()
                #If pre heat or front heat is present, in adapativbe coupling
                #We utilise div.q coupling to get rid of these fronts.
                #Otherwise, apply the exponential models. 
                if(self.init.yaml_file['Coupling_params']['Couple_adaptive']):
                    if(pre_heat_start_index > 0 or front_heat_start_index > 0):
                        qe = hfct_obj.divQHeatFlow()
                        pre_heat_fit_params = None
                        front_heat_fit_params = None
                else:
                    #Modify yaml for future write
                    fluid_obj.init.yaml_file['FixedParameters']['Preheat_StartIndex'] = pre_heat_start_index.item()
                    fluid_obj.init.yaml_file['FixedParameters']['Preheat_LastIndex'] = pre_heat_last_index.item()
                    fluid_obj.init.yaml_file['FixedParameters']['Frontheat_StartIndex'] = front_heat_start_index.item()
                    fluid_obj.init.yaml_file['FixedParameters']['Frontheat_LastIndex'] = front_heat_last_index.item()
            
            if self.init.yaml_file['Misc']['Zip']:
                io_obj.zipAndDelete()        
            
            #Finish by fluid init next set of files 
            fluid_obj.initHydroFromKinetic(io_obj.next_fluid_input_path, qe,
                                            pre_heat_fit_params, front_heat_fit_params)
            
            #update continue file
            np.savetxt(continue_step_path, np.array([cycle_no]), fmt = '%i')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Coupling of Hydrodynamic and Kinetic code. The two codes are specified in the input file.')
    parser.add_argument('-p', '--path', required=True, help = 'Give path to Input yml file.')
    args = vars(parser.parse_args())
    couple = Coupler(args['path'])
    couple.main()








