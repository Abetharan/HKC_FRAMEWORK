import os
import subprocess

os.environ["RUN"]="kappa"
print(os.environ["RUN"])
## laptop
os.environ["SRCDIR"]="/Users/shiki/Documents/Imperial_College_London/Ph.D./IMPACT/src"
#Work station
# os.environ["SRCDIR"]=os.environ["HOME"] +"/IMPACT/src"
print(os.environ["SRCDIR"])

os.environ["BASEDIR"]= os.getcwd()
print(os.environ["BASEDIR"])
  

#-----  Code generation control  -----
os.environ["FPCODE"]     = "fp2df1_fast"
os.environ["FPPROF"]   = ""
os.environ["RUN_ARGS"] = "-log_summary -ksp_monitor -pc_type asm -ksp_type bcgs -ksp_rtol 1e-20 -draw_pause -1 -optionstable"

#-----  Specify behaviour of this run script -os.environ["IS_RESTART
# os.environ["OVERWRITE_RUN"]
# os.environ["DONT_COMIPLE"]
os.environ["DONT_RUN"] = ""
# os.environ["RUN_IN_QUEUE"]  

#------  os.environ["COMPILE-TIME ARRAY DIMENSIONS   ------
os.environ["SYNCHRO_DIMS"] = ""
os.environ["SMT_PK_SPM_ON "] = "1"
os.environ["SMT_PC_COLS_ON"] = "0"

os.environ["NP"]  =   "1"
os.environ["NVM "]= "10"

os.environ["NXM"] =  "1"
os.environ["NYM"] =   "1"
# ...  Text output behaviour  ...
os.environ["TEXTOUTPUT_TO_STDOUT_ROOT "]    = "1"
path = os.environ["BASEDIR"] + "/" + os.environ["RUN"] + "/"
if not os.path.exists(path):
      os.makedirs(path)

file10 = open(path + "fort.10", "w")
file12 = open(path + "fort.12", "w")
file14 = open(path + "fort.14", "w")

fort10 = """ 
 $user_inp
       
       switch_fo_Cee_on          = .true.
       switch_f1_dt_on           = .true.
       switch_fo_dt_on           = .true.
       switch_Cee0_iter_on       = .FALSE. 
       switch_initial_cond_on    = .false.
       switch_heating_cooling_on = .TRUE.
       switch_ei_coll_fix_on     = .true.
       switch_displ_J_on         = .false.
       switch_Cee0_Krook_on      = .false.
       
       switch_disable_force_runtime_size_on = .false.
       
       Bz_implicitness_param       =  1.0d0
       switch_dvfo_centred_diff_on = .true.  
       switch_dv2f1_centred_diff_on = .true.
       switch_dBdt_centred_diff_on = .true.
       
       switch_Ohmic_Fara_consis_on = .true.
       switch_Ohmic_all_on         = .false.
       switch_hydro_on           = .false.
       switch_hydro_fo_on        = .false.
       switch_hydro_f1_on        = .false.
       
       wpe_over_nuei = 11.14576d0
       c_over_vte    = 30.455466d0
       atomic_Z      = 37.664907d0
       atomic_A      = 70.0d0
       
       p_SG      = 2.0d0
       q_SG      = 2.0d0
       p_SG_init = 2.0d0
       
       nv = 300
       nx = 6
       ny = 1
       
       dt   =     0.5d0
       tmax =   750.0d0
       
       x_bc_type     = 1
       xmin          = 0.0d0
       xmax          = 159424.7662d0
       grid_x_type   = 3
       grid_x_ratio  =  1.0d0
       grid_x_alpha  =  0.0d0
       grid_x_offset =  0.0d0
       
       y_bc_type     = 0
       ymin          = -3.0d3
       ymax          = +3.0d3
       grid_y_type   = 0
       grid_y_ratio  = 1.0d0
       grid_y_alpha  = 0.0d0
       grid_y_offset = 0.0d0
       
       vmax          = 30.0d0
       grid_v_type   =  2
       grid_v_ratio  = 30.0d0
       grid_v_alpha  =  0.0d0
       grid_v_offset =  0.0d0
       
       do_user_prof_sub = .TRUE.
       
       prof_Bz_ave    = 0.0d0
       
       do_user_heating_sub = .TRUE.
       
       switch_packed_sparse_on    = .true.
       switch_precomp_mat_cols_on = .false.
       matrix_solver              = 'PETSC'
       matrix_solver_tol          = 1.0d-15
       nonlin_tol                 = 1.0d-12
       nonlin_itmax               = 25
       CCee0_00_its_delta         = 10
       
       initial_cond_rel_dt = 1.0d-2
       initial_cond_nt     = 2
       
       do_out_data_compress = .false.
       op_time_mon_skip     = 10000
       op_restart_freq      = 100
 $end
"""
fort12 =""" 
    7
   0.1d0
   2.0d0
 150.0d0
 300.0d0
 450.0d0
 600.0d0
 750.0d0
   0.0d0
$end
"""
#------  Create graphics output selection file, "fort.14" (NAMELIST)  -----
fort14 = """ 
 $user_gra_sel

    op_save_on(20,1) = 1
    op_save_on(20,3) = 1
    op_save_on(21,3) = 1
    op_save_on(22,3) = 1
    op_save_on(23,3) = 1
    op_save_on(24,3) = 1

    op_save_on(33,3) = 1

 $end 
"""


file10.write(fort10)
file12.write(fort12)
file14.write(fort14)
os.chdir(os.environ['SRCDIR'])
os.system('./fp2df1_compile_run.sh')
#exit()
from scipy.interpolate import CubicSpline
import numpy as np 

init_file_path = ""
out_file_path = " "
fluid_code_path = " "
no_cycles  = 2

for i in range(0, no_cycles, 1):
      os.chdir(fluid_code_path)
      ##need to give some output location 
      os.system('./run')
      fluid_x  = np.readtxt(out_file_path)
      kinetic_x = np.readtxt(kinetic_x_path)
      fluid_ne = np.readtxt(out_file_path + "/NumberDensityE_" + i + ".txt")
      fluid_Te = np.readtxt(out_file_path + "/TemperaturE_" + i + ".txt")
      cs_ne = CubicSpline(fluid_x, fluid_ne)
      cs_Te = CubicSpline(fluid_x, fluid_Te)

      kinetic_ne = cs_ne(kinetic_x)
      kinetic_Te = cs_Te(kinetic_x)

      ##### Launch Impact
      os.chdir(os.environ['BASEDIR'])
      os.system('mpirun -np $NP ./fp2df1_run')









