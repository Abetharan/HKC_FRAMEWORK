import os
import shutil
import TmpFileCreator as tfc
import COUPLING_ELH1 as elh1
import COUPLING_IMPACT as impact
import IO as io
import Templating as temple
import Material as material
from colorama import Fore
from colorama import Style
def mapBitToChar(im, col, row):
    if im.getpixel((col, row)): return ' '
    else: return '#'

def prettyprint(text, color = None):
    length_to_print = 100 + len(text)
    if color is None:
        print("#"*length_to_print)
        print('\033[1m' + '#'*50 + text +'#'*50 + '\033[0m')
    
    else:
        print('\033[31m' + "#"*length_to_print + '\033[0m')
        print('\033[31m' + '#'*50 + text +'#'*50 + '\033[0m')
    
    print('\n')   

class main:
    """
    main class i.e. where all objects are passed into and the coupling is run.
    """
    #OPTIONS ARE SOL-KiT(1D NO B-field but e-e anisotropy)/IMPACT(2D with B no e-e anisotropy)
    _SWITCH_KINETIC_CODE = "IMPACT"
    _SWITCH_HYDRO_CODE = "ELH1"

    from PIL import Image, ImageFont, ImageDraw
    ShowText = 'COUPLING ' + _SWITCH_HYDRO_CODE + '-' + _SWITCH_KINETIC_CODE

    font = ImageFont.load_default().font #load the font
    font = ImageFont.truetype("/home/abetharan/Downloads/arial.ttf",15)
    size = font.getsize(ShowText)  #calc the size of text in pixels
    image = Image.new('1', size, 1)  #create a b/w image
    draw = ImageDraw.Draw(image)
    draw.text((0, 0), ShowText, font=font) #render the text to the bitmap
    for rownum in range(size[1]): 
    #scan the bitmap:
    # print ' ' for black pixel and 
    # print '#' for white one
        line = []
        for colnum in range(size[0]):
            if image.getpixel((colnum, rownum)): line.append(' '),
            else: line.append('#'),
        print(''.join(line))
    #for r in range(size[1]):
    #    print(''.join([mapBitToChar(image, c, r) for c in range(size[0])]))

    print('\n')   
    base_dir = "/home/abetharan/test/"
    k_src_dir = "/home/abetharan/IMPACT/src"
    run_name = "kap"
    k_nx = 60  
    k_ny = 1
    k_nv = 120
    k_np = 1
    k_dt = 0.9
    k_t_max = 1
    k_x_max = 1000
    k_v_max = 30
    cycles = 50
    
    
    # Material Properties
    Z = 37
    Ar = 157
    # Normalisation
    Te = 10000
    ne = 1e27
    Bz = 0

    # Fluid initial parameters
    f_init_path ="/media/abetharan/DATADRIVE1/Abetharan/Results/fixed_nx/Ncub60/cycle_0/fluid_input/"     
    f_src_dir = "/home/abetharan/HeadlessHydra/Source_Code/release/ELH1"

    f_nx =60
    f_cq = 2
    f_gamma = 1.4
    f_cfl = 0.85
    f_laser_wavelength = 351e-9  # 200e-9
    f_laser_power = 1e15
    f_dur_of_laser = 1e-10
    f_laser_loc = 'left'
    f_steps = 75
    f_fluid_t_max = 0  # 1e-15
    f_initial_dt = 1e-17
    f_dt_global_max = 1e-13
    f_dt_global_min = 1e-16
    f_feos_path_1 = None
    f_feos_path_2 = None
    f_output_freq = 1
    f_boundary_condition = "rigid"
    f_initialise_start_file_run = True
    initialise_all_folders = True

    for cycle_no in range(cycles):
        
        if cycle_no > 1:
            f_initialise_start_file_run = False
            initialise_all_folders = False

        io_obj = io.IO(base_dir, k_src_dir, run_name, f_src_dir, f_init_path,f_feos_path_1, f_feos_path_2, cycle_no, initialise_all_folders, cycles)
        if cycle_no == 0:
            io_obj.copyFluidInit(f_init_path)

        prettyprint(' RUNNING ' + _SWITCH_HYDRO_CODE, color = True)
        elh1_obj = elh1.ELH1(io_obj, f_nx, f_laser_wavelength, f_laser_power, f_dur_of_laser, 
                             f_steps, f_fluid_t_max, f_initial_dt, f_dt_global_max, f_dt_global_min,
                             f_output_freq, f_boundary_condition, f_initialise_start_file_run )
        elh1_obj.ELH1Run()

        impact_obj = impact.IMPACT(io_obj, k_np, k_nv, k_nx, k_ny, k_dt, k_t_max, k_x_max, k_v_max, restart_freq_ = 100)

        if cycle_no < 1:
            #Find normalisation and prepare IMAPCT for compiling and compile.
            prettyprint(' COMPILING' + _SWITCH_KINETIC_CODE)
            norms = impact_obj.normalisation(calc = True, ne = ne, Te = Te, Z = Z, Ar = Ar,  Bz = 0)
            impact_obj.setEnvVar()
            impact_obj.SetIMPACTParam()
            impact_obj.IMPACTCompile()            
        else:
            norms = impact_obj.normalisation(calc = False)

        prettyprint(' CONVERT ' + _SWITCH_HYDRO_CODE + 'TO ' + _SWITCH_KINETIC_CODE + 'COMPATIBLE ')
        #Convert ELH1 output to IMPACT Fortmat
        impact_obj.TxtToImpact()        

        #RUN IMPACT
        prettyprint(' RUN ' + _SWITCH_KINETIC_CODE, color =True )
        impact_obj.IMPACTRun()
        
        prettyprint(' CONVERT ' + _SWITCH_KINETIC_CODE + 'TO ' + _SWITCH_HYDRO_CODE + 'COMPATIBLE ')
        #Convert IMPACT FILES To ELH1 readable files
        impact_obj.ImpactToTxt()



















