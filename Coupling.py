import os
import subprocess
import Fort_generator as fort
os.environ["RUN"]="kappa"
print(os.environ["RUN"])
## laptop
#os.environ["SRCDIR"]="/Users/shiki/Documents/Imperial_College_London/Ph.D./IMPACT/src"
#Work station
os.environ["SRCDIR"]=os.environ["HOME"] +"/IMPACT/src"
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
#os.environ["DONT_RUN"] = ""
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

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

#Generates Fort with placeholders
fort12_times = "2 \n 0.1d0 \n 0.2d0 \n $end "
fort.fort_generator(path,fort12_times)

os.chdir(os.environ['BASEDIR'])
#os.system('./hydra_kappa.sh')
# exit()
from scipy.interpolate import CubicSpline
import numpy as np 

init_file_path = ""
out_file_path = " "
fluid_code_path = " "
no_cycles  = 2

# for i in range(0, no_cycles, 1):
#       os.chdir(fluid_code_path)
#       ##need to give some output location 
#       os.system('./run')
#       fluid_x  = np.readtxt(out_file_path)
#       kinetic_x = np.readtxt(kinetic_x_path)
#       fluid_ne = np.readtxt(out_file_path + "/NumberDensityE_" + i + ".txt")
#       fluid_Te = np.readtxt(out_file_path + "/TemperaturE_" + i + ".txt")
#       cs_ne = CubicSpline(fluid_x, fluid_ne)
#       cs_Te = CubicSpline(fluid_x, fluid_Te)

#       kinetic_ne = cs_ne(kinetic_x)
#       kinetic_Te = cs_Te(kinetic_x)

#       ##### Launch Impact
#       os.chdir(os.environ['BASEDIR'])
#       os.system('mpirun -np $NP ./fp2df1_run')









