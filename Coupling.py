import os
from string import Template
import shutil
import subprocess
import FortGenerator as fort
import SetFort10Param as setFort10
import SetFort12Param as SetFort12
import SetEnvVar as SetEnvVar
import scipy.constants as constants
import impact_module_py3 as cf
import HydroImpactIO as io

kb  = 1.38E-23
e = 1.6E-19
# simulation domain sizes and number of processors to use 
nx = 1
ny = 2 
nv = 3
np = 4
dt = 1 #as a ratio of collisional time i.e. 1 is collision time 
tmax = 20 #Number of collision times 
Atomic_Z = 60
Atomic_Ar = 157
#Set Environement variables for compiling
runName = "kappa"
#rcDir = "/Users/shiki/Documents/Imperial_College_London/Ph.D./IMPACT/src"
srcDir = "/home/abetharan/IMPACT/src"
#Return path is Run directory
path  = SetEnvVar.setEnvVar(nx, ny, nv, np, runName,srcDir)
if not os.path.exists(path):
    os.makedirs(path)

#Data location save
init_file_path = path 
out_file_path = path
fluid_code_path = "Source/run"
no_cycles  = 2

#Generate fort files 
fort12TimeStr = SetFort12.createFort12String(["1.0d0","2.0d0","3.0d0","5.0d0"])
fort.fort_generator(path,fort12TimeStr)


#Copy and Rename custom functions to run directory 
shutil.copyfile(os.environ["BASEDIR"] + "/heating.f", path + "/" + runName + "_heating.f")
shutil.copyfile(os.environ["BASEDIR"] + "/user_custom.f", path + "/" + runName + "user_custom.f")
shutil.copyfile(os.environ["BASEDIR"] + "/prof.f", path + "/" + runName + "_prof.f")

#----------------------------------------------------------------#
#Start Coupling sequence
os.chdir(os.environ['BASEDIR'])

os.system('./hydra_kappa.sh')
#os.system('bash ' + os.environ['SRCDIR'] + "/fp2df1_compile_run.sh")
# exit()



for i in range(1, no_cycles+1, 1):

    os.chdir(fluid_code_path)
    cycle_dump_path = path + "/cycle_" + i
    if os.path.exists(cycle_dump_path):
        os.makedirs(cycle_dump_path)

    ##need to give some output location 
    os.system('./run')

    avgNe, avgTe, normalised_values = io.HydroToImpact(cycle_dump_path,Atomic_Z,Atomic_Ar)
    
    #NOTES ON PARAMETERS
    # IF COULOMB PARAMETER < 3 FOKKER-PLANK NOT IDLE I.E. VFP PROB NOT BEST TO USE. 
    # IF VALUES ARE NEGATIVE MOST LIKELY COULOMB PARAMETER IS NOT BEING CALCULATED CORRECTLY
    # I.E. PLASMA TO COLD .. Make sure skin depth is of order < 1 as well.. 

    ne = avgNe #1.0e+19 #Note that NE is in cm^3 i.e. times M^3 by 10^-6 to convert 
    Te = avgTe / (e/kb) #eV
    Z = Atomic_Z 
    Bz = 0.0 #Tesla
    Ar = Atomic_Ar 
    wpe_over_nuei = normalised_values['wpe_over_nu_ei']
    c_over_vte = normalised_values['c_over_vte']
    #on the fly parameter changes
    fort10Param = setFort10.set_fort_10(wpe_over_nuei = wpe_over_nuei, c_over_vte = c_over_vte, 
                                        atomic_Z = Z, atomic_A = Bz, nv = nv, nx = nx ,ny = ny,dt = dt, tmax = tmax)
    filein = open(path +'/tmpfort.10', "r")
    src = Template(filein.read())
    k = src.substitute(fort10Param)
    Fort10 = open(path + "/fort.10", "w")
    Fort10.write(k)

    ##### Launch Impact
    os.chdir(path)
    os.system('mpirun -np $NP ./fp2df1_run')

    for file in os.listdir(path):
        filename, extension = os.path.splitext(file)
        if extension == ".xy" or ".xyz" or ".xyv":
            shutil.move(file, cycle_dump_path)





