""" 
Hydro-Kinetic Coupling.
Heat-Flow coupling between Hydro and Kibnetic codes. 
@author = Abetharan Antony
"""
import argparse
from colorama import Fore
from colorama import Style
from distutils.dir_util import copy_tree
import math
import numpy as np
from PIL import Image, ImageFont, ImageDraw
import os
import shutil
import utils as util
from couple_sol_kit.sol_kit import SOL_KIT
from couple_hykict.hykict import HyKiCT


class Coupler:
    """
    main class i.e. where all objects are passed into and the coupling is run.
    """
    def __init__(self, init_file_path_):
        self.yml_init_file_path = init_file_path_
        
    def start_print(self, fluid_code, kinetic_code):

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

    def pretty_print(self, text, color = None):
        length_to_print = 100 + len(text)
        if color is None:
            print("#"*length_to_print)
            print('\033[1m' + '#'*50 + text +'#'*50 + '\033[0m')
        
        else:
            print('\033[31m' + "#"*length_to_print + '\033[0m')
            print('\033[31m' + '#'*50 + text +'#'*50 + '\033[0m')
        
        print('\n')   
    

    def main(self):
        init = util.Input(self.yml_init_file_path)
        RUN_PATH = os.path.join(init.yaml_file['Paths']['Base_dir'],
                     init.yaml_file['Paths']['Run_name'])        

        if not init.yaml_file['Misc']['HPC']:
            self.start_print(init.yaml_file["Codes"]["Fluid_code"],
                             init.yaml_file["Codes"]["Kinetic_code"])

        overwrite = init.yaml_file['Misc']['Overwrite']
        cycles = init.yaml_file['Coupling_params']['Cycles']
        continue_step_path = os.path.join(RUN_PATH, 'CONTINUE_STEP.txt') 
        initialise_all_folders = False
    
        if init.yaml_file['Misc']['Continue']:
            #Continue in this framework is designed such that
            #a run does not have to fail (for whatever reason) and then continue
            #it will keep going until its objective has been reached which is
            #all cycles have been run through. 
            #Smarter continuing should be implemented such that errors due to
            #codes crashing should also crash the framework. Wheras,
            #HPC facility failures should just make it continue. 
            if os.path.exists(continue_step_path):
                start_cycle = np.loadtxt(continue_step_path, dtype=np.int) + 1
            else:
                start_cycle = 0
                initialise_all_folders = True
                overwrite = True
        else:
            initialise_all_folders = True
            start_cycle = 0
            overwrite = True

        #Create Objects 
        io_obj = util.IO(
                        init.yaml_file['Paths']['Run_name'],
                        init.yaml_file['Paths']['Base_dir'], 
                        init.yaml_file['Paths']['K_src_dir'],
                        init.yaml_file['Paths']['F_src_dir'], 
                        init.yaml_file['Paths']['Init_path'],
                        start_cycle, overwrite, cycles, 
                        initialise_all_folders)
        hfct_obj = util.HeatFlowCouplingTools()
        fluid_obj = HyKiCT(io_obj, init.yaml_file['Paths']['F_config_path'])
                         
        if(init.yaml_file['Codes']['Kinetic_code'] == 'sol_kit'):
            kin_obj = SOL_KIT(io_obj,
                        init.yaml_file['Paths']['K_config_path'],
                        init.yaml_file['Misc']['HPC'])
        #Impact not ready
        # else:
            #kin_obj = IMPACT()
        
        kin_obj.init.yaml_file['Params']['Nx'] = init.yaml_file['Coupling_params']['Nx']
        fluid_obj.init.yaml_file['FixedParameters']['Nx'] = init.yaml_file['Coupling_params']['Nx']
        #For restarting/continue
        #The copying could in principal be removed. 
        #Output directory just self-consistent this way.
        if start_cycle == 0:
            copy_tree(init.yaml_file['Paths']['Init_path'], io_obj.fluid_input_path)
            np.savetxt(os.path.join(RUN_PATH, 'NO_CYCLES.txt'), np.array([cycles - 1]), fmt = '%i' )       
        
        #Specific to HyKiCT
        f_nt = fluid_obj.init.yaml_file['TimeParameters']['steps'] 
        f_tmax = fluid_obj.init.yaml_file['TimeParameters']['t_max']
        if init.yaml_file['Coupling_params']['Start_from_kinetic']:
            fluid_obj.init.yaml_file['TimeParameters']['steps'] = 0
            fluid_obj.init.yaml_file['TimeParameters']['t_max'] = 0
            fluid_obj.init.yaml_file['Switches']['CoupleDivQ'] = False
            fluid_obj.init.yaml_file['Switches']['CoupleMulti'] = False

        for cycle_no in range(start_cycle, cycles, 1):
            self.pretty_print(' RUNNING CYCLE ' + str(cycle_no)) 

            if cycle_no >= 1:
                io_obj.cycle_counter = cycle_no
                io_obj.nextCyclePathManager()
                kin_obj._cycle_path = io_obj.cycle_dump_path
                fluid_obj._cycle_path = io_obj.cycle_dump_path
                if init.yaml_file['Coupling_params']['Start_from_kinetic']:
                    fluid_obj.init.yaml_file['TimeParameters']['steps'] = f_nt
                    fluid_obj.init.yaml_file['TimeParameters']['t_max'] = f_tmax
                    fluid_obj.init.yaml_file['Switches']['CoupleDivQ'] = init.yaml_file['Coupling_params']['Couple_divq']
                    fluid_obj.init.yaml_file['Switches']['CoupleMulti'] = init.yaml_file['Coupling_params']['Couple_multi']

            #Fluid Code runs here 
            #In this example running HyKiCT
            self.pretty_print(' RUNNING ' + init.yaml_file['Codes']['Fluid_code'], color = True)
            
            #Set Paths 
            fluid_obj.init.yaml_file['Paths']['Init_Path'] = io_obj.fluid_input_path
            fluid_obj.init.yaml_file['Paths']['Out_Path'] = io_obj.fluid_output_path
            fluid_obj._cycle_dump_path = io_obj.cycle_dump_path
            ################
            #Any other parameter updates that needs to be set can be added here.
            #################

            fluid_obj.setFiles()
            fluid_obj.Run()

            #the _ are unused variables fluid_v and fluid_laser, for future
            #use these may be required and thus, left in. 
            (fluid_x_grid, fluid_x_centered_grid, _, fluid_ne, fluid_Te,
            fluid_Z, _, fluid_mass) = fluid_obj.getLastStepQuants()

            #Init heat_flow_tools here 
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

            #Breaks here ... Last cycle allowed to run hydro step
            if cycle_no == cycles - 1:
                if init.yaml_file['Misc']['Zip']:
                    io_obj.zipAndDelete()        
                break 

            #Set any parameter ifles
            #Init all load in files form hydro
            #RUN
            #Move files if neccessary

            self.pretty_print(' RUNNING ' + init.yaml_file['Codes']['Kinetic_code'], color = True)
            if cycle_no == 0:
                #input and output unchanged 
                kin_obj.setFiles()

            kin_obj._kinetic_input_path = io_obj.kinetic_input_path
            kin_obj._kinetic_output_path = io_obj.kinetic_output_path
            
            kin_obj.InitFromHydro(fluid_x_grid, fluid_x_centered_grid, 
                                fluid_Te, fluid_ne, fluid_Z)
            kin_obj.Run()
            #Set hfct to contain the vfp heat flow to do the neccessary coupling calcs.
            hfct_obj.vfp_heat= kin_obj.getLastHeatFlow()
            kin_obj.moveFiles()

            #Do Coupling Parameters here
            if init.yaml_file['Coupling_params']['Couple_divq']:
                #qe here is div.q_vfp 
                qe = hfct_obj.divQHeatFlow()
            
            elif init.yaml_file['Coupling_params']['Couple_multi']:
                #qe here is q_vfp/q_sh
                (qe, pre_heat_start_index, pre_heat_last_index,
                pre_heat_fit_params, front_heat_start_index, 
                front_heat_last_index, front_heat_fit_params)  = hfct_obj.multiplier()

            if init.yaml_file['Misc']['Zip']:
                io_obj.zipAndDelete()        
            #update continue file
            np.savetxt(continue_step_path, np.array([cycle_no]), fmt = '%i')
            
            
            #Finish by init next set of files 
            fluid_obj.initHydroFromKinetic(io_obj.next_fluid_input_path, qe,
                                            pre_heat_fit_params, front_heat_fit_params)
            #Modify yaml for future write
            fluid_obj.init.yaml_file['FixedParameters']['Preheat_StartIndex'] = pre_heat_start_index.item()
            fluid_obj.init.yaml_file['FixedParameters']['Preheat_LastIndex'] = pre_heat_last_index.item()
            fluid_obj.init.yaml_file['FixedParameters']['Frontheat_StartIndex'] = front_heat_start_index.item()
            fluid_obj.init.yaml_file['FixedParameters']['Frontheat_LastIndex'] = front_heat_last_index.item()
           
            #To be DELETED
            #OBSOLETE
            #Kinetic Codes run here
            # if _SWITCH_KINETIC_CODE == "IMPACT":
            #     impact_obj = impact.IMPACT(io_obj, k_np, k_nv, k_nx, k_ny, k_dt, k_t_max, k_x_max, k_v_max, k_bc, 100, _SWITCH_CX1)
            #     norms = impact_obj.normalisation(calc = True, ne = ne, Te = Te, Z = Z, Ar = Ar,  Bz = 0)
            #     if cycle_no < 1:
            #         #Find normalisation and prepare IMAPCT for compiling and compile.
            #         self.pretty_print(' COMPILING' + _SWITCH_KINETIC_CODE)
            #         impact_obj.setEnvVar()
            #         impact_obj.setIMPACTParam()
            #         impact_obj.setCustomFuncs()
            #         impact_obj.IMPACTCompile()            
            #     else:
                #     norms = impact_obj.normalisation(calc = False)

                # self.pretty_print(' CONVERT ' + _SWITCH_HYDRO_CODE + 'TO ' + _SWITCH_KINETIC_CODE + 'COMPATIBLE ')
                # #Convert ELH1 output to IMPACT Fortmat
                # impact_obj.InitIMPACTFromHydro()        

                # #RUN IMPACT
                # self.pretty_print(' RUN ' + _SWITCH_KINETIC_CODE, color =True )
                # impact_obj.IMPACTRun()
                
                # self.pretty_print(' CONVERT ' + _SWITCH_KINETIC_CODE + 'TO ' + _SWITCH_HYDRO_CODE + 'COMPATIBLE ')
                #Convert IMPACT FILES To ELH1 readable files
                # if not last_cycle:
                #     impact_obj.IMPACTInitHydroNextFiles()
                #     impact_obj.moveIMPACTFile()
                # else:
                #     impact_obj.moveIMPACTFile()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Coupling of Hydrodynamic and Kinetic code. The two codes are specified in the input file.')
    parser.add_argument('-p', '--path', required=True, help = 'Give path to Input yml file.')
    args = vars(parser.parse_args())
    couple = Coupler(args['path'])
    couple.main()








