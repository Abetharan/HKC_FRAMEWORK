import os
from string import Template
import shutil
import subprocess
import FortGenerator as fort
import SetFort10Param as setFort10
import SetFort12Param as SetFort12
import SetEnvVar as SetEnvVar
import impact_norms_py3 as ImNorms
import scipy.constants as constants
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
srcDir = "/Users/shiki/Documents/Imperial_College_London/Ph.D./IMPACT/src"
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

#os.system('./hydra_kappa.sh')
# exit()
from scipy.interpolate import CubicSpline
import numpy as np 


for i in range(0, no_cycles, 1):

    os.chdir(fluid_code_path)
    cycle_dump_path = path + "/cycle_" + i
    if os.path.exists(cycle_dump_path):
        os.makedirs(cycle_dump_path)

    ##need to give some output location 
    os.system('./run')
 
    fluid_x  = np.loadtxt(cycle_dump_path + "/Coord_" + i + ".txt")
    fluid_ne = np.loadtxt(cycle_dump_path + "/NumberDensityE_" + i + ".txt")
    fluid_ni = np.loadtxt(cycle_dump_path + "/NumberDensityI_" + i + ".txt")
    fluid_Te = np.loadtxt(cycle_dump_path + "/TemperaturE_" + i + ".txt")
    fluid_las_dep = np.loadtxt(cycle_dump_path + "/InvBrem_" + i + ".txt")
    fluid_brem = np.loadtxt(cycle_dump_path + "/Brem_" + i + ".txt")

   ##NOTES ON PARAMETERS
    # IF COULOMB PARAMETER < 3 FOKKER-PLANK NOT IDLE I.E. VFP PROB NOT BEST TO USE. 
    # IF VALUES ARE NEGATIVE MOST LIKELY COULOMB PARAMETER IS NOT BEING CALCULATED CORRECTLY
    # I.E. PLASMA TO COLD .. Make sure skin depth is of order < 1 as well.. 

    ne = np.average(fluid_ne) #1.0e+19 #Note that NE is in cm^3 i.e. times M^3 by 10^-6 to convert 
    Te = np.average(fluid_Te) / (e/kb)  # 0.01 #eV
    Z = Atomic_Z 
    Bz = 0.0 #3.75 #3.75 #Tesla
    Ar = Atomic_Ar 
    normalised_values = ImNorms.impact_inputs(ne,Te,Z,Bz,Ar)
    wpe_over_nuei = normalised_values['wpe_over_nu_ei']
    c_over_vte = normalised_values['c_over_vte']
    #on the fly parameter changes
    fort10Param = setFort10.set_fort_10(wpe_over_nuei = wpe_over_nuei, c_over_vte = c_over_vte, atomic_Z = Z, atomic_A = Bz, nv = nv, nx = nx ,ny = ny,dt = dt, tmax = tmax)
    filein = open(path +'/tmpfort.10', "r")
    src = Template(filein.read())
    k = src.substitute(fort10Param)
    Fort10 = open(path + "/fort.10", "w")
    Fort10.write(k)

    ## NOrmalise SI to Impact norms 
    x_norm = fluid_x / normalised_values["lambda_mfp"]
    ne_norm = fluid_ne / ImNorms.calc_norms("n", normalised_values, forced_power= 0)
    ni_norm = fluid_ne / ImNorms.calc_norms("n", normalised_values, forced_power= 0)
    Te_norm = fluid_Te/ ImNorms.calc_norms("Te", normalised_values, forced_power= 0)


    kinetic_x = np.linspace(fluid_x[0], fluid_x[-1], nx)
               # np.geomspace(fluid_x[0, fluid_x[-1], nx)
               # np.logspace(fluid_x[0, fluid_x[-1], nx)

    cs_ne = CubicSpline(x_norm, ne_norm)
    cs_ni = CubicSpline(x_norm, ni_norm)
    cs_Te = CubicSpline(x_norm, Te_norm)

    kinetic_ne = cs_ne(kinetic_x)
    kinetic_ni = cs_ni(kinetic_x)
    kinetic_Te = cs_Te(kinetic_x)

    impactNeFile = open(path + "/" + runName + "_eden.xy", "w")
    impactNiFile = open(path + "/" + runName + "_ionden.xy", "w")
    impactTeFile = open(path + "/" + runName + "_tmat.xy", "w")
    impactXFile = open(path + "/" + runName + "_xf.xy", "w")

    impactNeFile.write(kinetic_ne)
    impactNiFile.write(kinetic_ni)
    impactTeFile.write(kinetic_Te)
    impactXFile.write(kinetic_x)

    ##### Launch Impact
    for file in os.listdir(path):
        filename, extension = os.path.splitext(file)
        if extension == ".xy" or ".xyz" or ".xyv":
            shutil.move(file, cycle_dump_path)

    os.chdir(path)
    os.system('mpirun -np $NP ./fp2df1_run')

    ##Convert to SI from Impact Norms






