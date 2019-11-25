""" 
Hydro-Kinetic Coupling.
Heat-Flow coupling between Hydro and Kibnetic codes. 
@author = Abetharan Antony
Last Update = 25/11/19
"""
import os
import shutil
import math
import argparse
import numpy as np
from Input import Input
import TmpFileCreator as tfc
import COUPLING_ELH1 as elh1
import COUPLING_IMPACT as impact
import COUPLING_SOL_KIT as solkit
import IO as io
import Templating as temple
import Material as material
class Coupler:
    """
    main class i.e. where all objects are passed into and the coupling is run.
    """
    #OPTIONS ARE SOL-KiT(1D NO B-field but e-e anisotropy)/IMPACT(2D with B no e-e anisotropy)
    #_SWITCH_KINETIC_CODE = "SOL_KIT"
    #_SWITCH_KINETIC_CODE = "IMPACT"
    
    def __init__(self, init_file_path_):
        self.yml_init_file_path = init_file_path_
    
    
    def main(self):
        init = Input(self.yml_init_file_path)
        _SWITCH_HYDRO_CODE = "ELH1"
        _SWITCH_KINETIC_CODE = init.KINETIC_CODE
        _SWITCH_CONTINUE = init.CONTINUE
        _SWITCH_CX1 = init.CX1
        
        def prettyprint(text, color = None):
            length_to_print = 100 + len(text)
            if color is None:
                print("#"*length_to_print)
                print('\033[1m' + '#'*50 + text +'#'*50 + '\033[0m')
            
            else:
                print('\033[31m' + "#"*length_to_print + '\033[0m')
                print('\033[31m' + '#'*50 + text +'#'*50 + '\033[0m')
            
            print('\n')   
        
        if not _SWITCH_CX1:
        
            from colorama import Fore
            from colorama import Style

            from PIL import Image, ImageFont, ImageDraw
            ShowText = 'COUPLING ' + _SWITCH_HYDRO_CODE + '-' + _SWITCH_KINETIC_CODE

            font = ImageFont.load_default().font #load the font
            font = ImageFont.truetype("/home/abetharan/Downloads/arial.ttf",15)
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

        base_dir = init.BASE_DIR 
        run_path = init.RUN_PATH
        impact_src_dir = init.IMPACT_SRC_DIR  #"../IMPACT/src"
        sol_kit_src_dir = init.SOL_KIT_SRC_DIR #"../SOL-KiT/"
        run_name = init.RUN_NAME
        
        if _SWITCH_KINETIC_CODE == 'IMPACT':
            k_src_dir = impact_src_dir
            k_nx = init.K_NX 
            k_ny = init.K_NY 
            k_nv = init.K_NV 
            k_np = init.K_NP 
            k_dt = init.K_DT 
            k_t_max = init.K_T_MAX
            k_x_max = init.K_X_MAX
            k_v_max = init.K_V_MAX
            k_bc = init.K_BC
        
        if _SWITCH_KINETIC_CODE == "SOL_KIT":
            k_src_dir = sol_kit_src_dir
            k_nx = init.K_NX 
            k_ny = init.K_NY 
            k_nv = init.K_NV 
            k_np = init.K_NP 
            k_dt = init.K_DT 
            k_nt = init.K_NT
            k_dx = init.K_DX
            k_dv = init.K_DV
            k_v_multi = init.K_DV_MULTI
            k_pre_step_nt = init.K_PRE_STEP_NT
            k_pre_step_dt = init.K_PRE_STEP_DT
            k_save_freq = math.ceil(k_nt * 0.1)        
            k_l_max = init.K_L_MAX
            k_bc = init.K_BC


        # Material Properties
        Z = init.Z
        Ar = init.AR
        # Normalisation
        Te = init.TE
        ne = init.NE
        Bz = init.BZ

        # Fluid initial parameters
        f_init_path = init.F_INIT_PATH
        f_src_dir = init.F_SRC_DIR
        f_nx =k_nx
        f_cq = init.F_CQ 
        f_gamma = init.F_GAMMA
        f_cfl = init.F_CFL
        f_laser_wavelength = init.F_LASER_WAVELENGTH
        f_laser_power = init.F_LASER_POWER
        f_dur_of_laser = init.F_DUR_OF_LASER
        f_laser_loc = init.F_LASER_POWER
        f_steps = init.F_STEPS
        f_fluid_t_max = init.F_FLUID_T_MAX
        f_initial_dt = init.F_INITIAL_DT
        f_dt_global_max = init.F_DT_GLOBAL_MAX
        f_dt_global_min = init.F_DT_GLOBAL_MIN
        f_feos_path_1 = init.F_FEOS_PATH_1
        f_feos_path_2 = init.F_FEOS_PATH_2
        f_output_freq = init.F_OUTPUT_FREQ
        f_boundary_condition = init.F_BOUNDARY_CONDITION
        f_initialise_start_file_run = True
        
        overwrite = False
        last_cycle = False
        cycles = init.CYCLES
        continue_step_path = os.path.join(run_path, 'CONTINUE_STEP.txt') 
        initialise_all_folders = False
        copy_sol_kit = False
    
        if _SWITCH_CONTINUE:
        
            if os.path.exists(continue_step_path):
                start_cycle = np.loadtxt(continue_step_path, dtype=np.int) + 1
            else:
                start_cycle = 0
                initialise_all_folders = True
                copy_sol_kit = True
                overwrite = True
        else:
            initialise_all_folders = True
            copy_sol_kit = True
            start_cycle = 0
            overwrite = True

        if _SWITCH_KINETIC_CODE == "IMPACT":
            k_nx = int(k_nx / k_np)

        for cycle_no in range(start_cycle, cycles, 1):
            prettyprint(' RUNNING CYCLE ' + str(cycle_no)) 

            if cycle_no >= 1:
                f_initialise_start_file_run = False 
                initialise_all_folders = False
                overwrite = False
                copy_sol_kit = False

            if cycle_no == cycles - 1:
                last_cycle = True
            
            print('making stuff') 
            io_obj = io.IO(base_dir, k_src_dir, run_name, f_src_dir, f_init_path,f_feos_path_1, f_feos_path_2, 
                            cycle_no, overwrite, last_cycle, initialise_all_folders, cycles)

            if cycle_no == 0:
                io_obj.copyFluidInit(f_init_path)
                np.savetxt(os.path.join(run_path, 'NO_CYCLES.txt'), np.array([cycles - 1]), fmt = '%i' )       
            
            prettyprint(' RUNNING ' + _SWITCH_HYDRO_CODE, color = True)
            elh1_obj = elh1.ELH1(io_obj, f_nx, f_laser_wavelength, f_laser_power, f_dur_of_laser, 
                                f_steps, f_fluid_t_max, f_initial_dt, f_dt_global_max, f_dt_global_min,
                                f_output_freq, f_boundary_condition, init.COUPLEDIVQ, init.COUPLEMULTI, f_initialise_start_file_run )
            elh1_obj.ELH1Run()
            #Breaks here ... Last cycle allowed to run hydro step
            #kinetic would be useless as qe generated is not used. 
            if cycle_no == cycles - 1:
                if init.ZIP:
                    io_obj.zipAndDelete()        
                break;

            if _SWITCH_KINETIC_CODE == "IMPACT":
                impact_obj = impact.IMPACT(io_obj, k_np, k_nv, k_nx, k_ny, k_dt, k_t_max, k_x_max, k_v_max, k_bc, 100, _SWITCH_CX1)
                norms = impact_obj.normalisation(calc = True, ne = ne, Te = Te, Z = Z, Ar = Ar,  Bz = 0)
                if cycle_no < 1:
                    #Find normalisation and prepare IMAPCT for compiling and compile.
                    prettyprint(' COMPILING' + _SWITCH_KINETIC_CODE)
                    impact_obj.setEnvVar()
                    impact_obj.setIMPACTParam()
                    impact_obj.setCustomFuncs()
                    impact_obj.IMPACTCompile()            
                else:
                    norms = impact_obj.normalisation(calc = False)

                prettyprint(' CONVERT ' + _SWITCH_HYDRO_CODE + 'TO ' + _SWITCH_KINETIC_CODE + 'COMPATIBLE ')
                #Convert ELH1 output to IMPACT Fortmat
                impact_obj.InitIMPACTFromHydro()        

                #RUN IMPACT
                prettyprint(' RUN ' + _SWITCH_KINETIC_CODE, color =True )
                impact_obj.IMPACTRun()
                
                prettyprint(' CONVERT ' + _SWITCH_KINETIC_CODE + 'TO ' + _SWITCH_HYDRO_CODE + 'COMPATIBLE ')
                #Convert IMPACT FILES To ELH1 readable files
                if not last_cycle:
                    impact_obj.IMPACTInitHydroNextFiles()
                    impact_obj.moveIMPACTFile()
                else:
                    impact_obj.moveIMPACTFile()

            if _SWITCH_KINETIC_CODE == "SOL_KIT":
                sol_kit_obj = solkit.SOL_KIT(io_obj, k_np, k_nv, k_nx, k_nt, k_pre_step_nt, k_dx, 
                                    k_dv, k_v_multi, k_dt, k_pre_step_dt, k_l_max, k_save_freq, Z, Ar, ne, Te,
                                    copy_sol_kit, _SWITCH_CX1, init.COUPLEDIVQ, init.COUPLEMULTI) 
                sol_kit_obj.setSOL_KITFiles(switch_boundary_condition_ = k_bc)
                sol_kit_obj.InitSOL_KITFromHydro()
                sol_kit_obj.SOL_KITRun()
                if not last_cycle:
                    sol_kit_obj.SOL_KITInitNextHydroFiles()
                    sol_kit_obj.moveSOL_KITFiles()
                else:    
                    sol_kit_obj.moveSOL_KITFiles()

            np.savetxt(continue_step_path, np.array([cycle_no]), fmt = '%i')
            
            if init.ZIP:
                io_obj.zipAndDelete()        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Coupling of Hydrodynamic and Kinetic code. The two codes are specified in the input file.')
    parser.add_argument('-p', '--path', required=True, help = 'Give path to Input yml file.')
    args = vars(parser.parse_args())
    couple = Coupler(args['path'])
    couple.main()















