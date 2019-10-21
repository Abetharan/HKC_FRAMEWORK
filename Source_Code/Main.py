import os
import shutil
import Coupling as cpl
import TmpFileCreator as tfc
import COUPLING_ELH1 as elh1
import COUPLING_IMPACT as impact
import IO as io
import Templating as temple
import Material as material
import impact_norms_py3 as norms

class main:
    """
    main class i.e. where all objects are passed into and the coupling is run.
    """
    base_dir = "/media/abetharan/DATADRIVE1/Abetharan/"
    k_src_dir = "/home/abetharan/IMPACT/src"
    run_name = "Ncub18"
    k_nx = 60  
    k_ny = 1
    k_nv = 120
    k_np = 5
    k_dt = 300
    k_t_max = 18
    k_x_max = 1000
    k_v_max = 30
    cycles = 50
    
    
    # Material Properties
    Z = 64
    Ar = 157
    # Normalisation
    Te = 1000
    ne = 1e27
    Bz = 0
    normalisation = norms.impact_inputs(ne, Te, Z, Ar, Bz)

    # Fluid initial parameters
    f_init_path ="/media/abetharan/DATADRIVE1/Abetharan/Results/fixed_nx/Ncub60/cycle_0/fluid_input/"     
    f_src_dir = "/home/abetharan/HeadlessHydra/Source_Code/run"

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

    for cycle_no in range(cycles):

        io_obj = io.IO(base_dir, k_src_dir, run_name, f_src_dir, f_init_path,f_feos_path_1, f_feos_path_2, cycle_no)

        impact_obj = impact.IMPACT(io_obj, normalisation, k_np, k_nv, k_nx, k_ny, k_dt, k_t_max, k_x_max, k_v_max)

        elh1_obj = elh1.ELH1(io_obj, f_nx, f_laser_wavelength, f_laser_power, f_dur_of_laser, 
                             f_steps, f_fluid_t_max, f_initial_dt, f_dt_global_max, f_dt_global_min,
                             f_output_freq, f_boundary_condition)

        



















