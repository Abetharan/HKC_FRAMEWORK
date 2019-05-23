import os
from string import Template
import shutil
import subprocess
import FortGenerator as fort
import SetFort10Param as setFort10
import SetFort12Param as SetFort12
import SetHydroInit as SetHydro
import SetEnvVar as SetEnvVar
import scipy.constants as constants
import impact_module_py3 as cf
import HydroImpactIO as io

kb  = 1.38E-23
e = 1.6E-19
# simulation domain sizes and number of processors to use 
KINETIC_nx = 30
KINETIC_ny = 1 
KINETIC_nv = 300
KINETIC_np = 1
dt = 1 #as a ratio of collisional time i.e. 1 is collision time 
tmax = 20 #Number of collision times 
Atomic_Z = 1#60
Atomic_Ar = 1#157

#Fluid initial parameters 
Cq = 2
Gamma = 1.4
CFL = 0.85
LaserWavelength = 200e-9
LaserPower = 1e18
durOfLaser = 1e-10
steps = 100
tmax = 10
initialDt = 1e-19
OutputFrequency = 10
#Set Environement variables for compiling
runName = "Test_1"
#rcDir = "/Users/shiki/Documents/Imperial_College_London/Ph.D./IMPACT/src"
IMPACT_SRC_DIR_ = "/home/abetharan/IMPACT/src"
FLUID_SRC_DIR_ = "/home/abetharan/HeadlessHydra/Source_Code/run"
#IMPACT_SRC_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./IMPACT/src"
#FLUID_SRC_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/Source_Code/run"
#INITIAL_CONDITIONS_FILE_PATH = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/init_data/"
INITIAL_CONDITIONS_FILE_PATH = "/home/abetharan/HeadlessHydra/init_data/"
#Return path is Run directory
RunPath  = SetEnvVar.setEnvVar(KINETIC_nx, KINETIC_ny, KINETIC_nv, KINETIC_np, runName, IMPACT_SRC_DIR_)

##If set true any conflicting files will be created 
_OVERWRITE_OK_ = True

if os.path.exists(RunPath):
    if _OVERWRITE_OK_:
        shutil.rmtree(RunPath)
        os.makedirs(RunPath)   

#Data location save
no_cycles  = 2
Fluid_nx = 100

for i in range(1, no_cycles+1, 1):
    #Cycle dump
    cycle_dump_path = RunPath + "cycle_" + str(i)
    if not os.path.exists(cycle_dump_path):
        os.makedirs(cycle_dump_path)

    ##Set on the fly Fluid parameters
    if i == 1:
        init = INITIAL_CONDITIONS_FILE_PATH
        cycle_init_path =  cycle_dump_path + "/fluid_init_data/"
        if not os.path.isdir(cycle_init_path):
            shutil.copytree(init, cycle_init_path)
    else:
        cycle_init_path = cycle_dump_path + "/fluid_init_data/"
        os.makedirs(cycle_init_path)   
        previous_cycle_dump_path = RunPath + "cycle_" + str(i - 1) 
        io.ImpactToHydro(previous_cycle_dump_path, normalised_values, Z, Ar, i, Gamma, LaserWavelength, LaserPower, Fluid_nx)


    #Create Fluid Dump Folder
    fluid_dump_path = cycle_dump_path + "/fluid_out/"
    if not os.path.isdir(fluid_dump_path):
        os.makedirs(fluid_dump_path)
    
    #Create Kinetic Dump folder
    kinetic_dump_path = cycle_dump_path + "/kinetic_out"
    if not os.path.isdir(kinetic_dump_path):
        os.makedirs(kinetic_dump_path)


    #Set Hydro param
    hydroparam = SetHydro.set_hydro_init(Fluid_nx, Atomic_Ar, Atomic_Z, Cq, Gamma, CFL, LaserWavelength,  LaserPower,
    durOfLaser, steps, tmax, initialDt, OutputFrequency, cycle_init_path, fluid_dump_path) 

    # Handling templating to create the init file for fluid code
    filein = open(os.environ['BASEDIR'] +'/tmpHydroParameterInit.txt', "r")
    src = Template(filein.read())
    fluidTemplate = src.substitute(hydroparam)
    HydroInit = open(cycle_dump_path + "/HydroParameterInit.txt", "w")
    HydroInit.write(fluidTemplate)
    HydroInit.close()
    #Call fluid code
    p = subprocess.run([FLUID_SRC_DIR_, '-p', cycle_dump_path+'/HydroParameterInit.txt'])
    
    #Convert all relevant parameters to impact norm and prepare files for load in.
    nc, norm_Te, normalised_values = io.HydroToImpact(fluid_dump_path, cycle_dump_path, Atomic_Z, Atomic_Ar, LaserWavelength, Fluid_nx)
    
    ##Impac Reference plasma. 
    ne = nc #1.0e+19 #Note that NE is in cm^3 i.e. times M^3 by 10^-6 to convert 
    Te = norm_Te / (e/kb) #eV
    Z = Atomic_Z 
    Bz = 0.0 #Tesla
    Ar = Atomic_Ar 
    wpe_over_nuei = normalised_values['wpe_over_nu_ei']
    c_over_vte = normalised_values['c_over_vte']

    if normalised_values["log_lambda"] < 0:
        print("Normalisation values not within good IMPACT PARAMETERS ... change")
        exit(1)

    #Generate fort files     
    #On the fly fort10 changes as requierd of fort 10 files.
    fort10Param = setFort10.set_fort_10(wpe_over_nuei = wpe_over_nuei, c_over_vte = c_over_vte, 
                                        atomic_Z = Z, atomic_A = Bz, nv = KINETIC_nv, nx = KINETIC_nx, ny = KINETIC_ny, dt = dt, tmax = tmax, do_user_prof_sub = ".false.")


    filein = open(os.environ['BASEDIR'] +'/tmpfort.10', "r")
    src = Template(filein.read())
    fort10template = src.substitute(fort10Param)
    Fort10 = open(RunPath + "/fort.10", "w")
    Fort10.write(fort10template)
    Fort10.close()
    fort12TimeStr = SetFort12.createFort12String(["1.0d0","2.0d0","3.0d0","5.0d0"])
    fort.fort_generator(RunPath,fort12TimeStr)

    if i == 1:

        #Copy and Rename custom functions to run directory 
        shutil.copyfile(os.environ["BASEDIR"] + "/heating.f", RunPath + "/" + runName + "_heating.f")
        shutil.copyfile(os.environ["BASEDIR"] + "/control.dat.dummy", RunPath + "/fp2df1_control.dat.dummy")
        filein = open(os.environ['BASEDIR'] +'/user_custom.f', "r")
        src = Template(filein.read())
        UserTemplate = src.substitute({'PATH':"\'" + RunPath + "\'", 'RUNNAME':"\'" + runName + "\'"})
        userCustom = open(RunPath + "/" + runName + "_user_custom.f", "w")
        userCustom.write(UserTemplate)
        userCustom.close()
       # shutil.copyfile(os.environ["BASEDIR"] + "/user_custom.f", RunPath + "/" + runName + "_user_custom.f")
        filein = open(os.environ['BASEDIR'] +'/prof.f', "r")
        src = Template(filein.read())
        #UserTemplate = src.substitute({'PATH':"\"" + RunPath + "\"", 'RUNNAME':"\"" + runName + "\""})
        userProf = open(RunPath + "/" + runName + "_prof.f", "w")
        userProf.write(UserTemplate)
        userProf.close()
        
       # shutil.copyfile(os.environ["BASEDIR"] + "/prof.f", RunPath + "/" + runName + "_prof.f")

        #----------------------------------------------------------------#
        #Start Coupling sequence
        #os.chdir(os.environ['BASEDIR'])
        os.system('./hydra_kappa.sh')
        #os.system('bash ' + os.environ['SRCDIR'] + "/fp2df1_compile_run.sh")
        #exit()

    

    filenames = ['ionden', 'rad_to_electron', 'xf', 'eden', 'laserdep', 'tmat']
    for name in filenames:
        if os.path.splitext(cycle_dump_path + "/" + runName + "_" + name  + ".xy")[-1] != ".xy":
            continue
        shutil.move(cycle_dump_path + "/" + runName + "_" + name  + ".xy", RunPath + "/" + runName + "_" + name + ".xy" )
    
    ##### Launch Impact
    os.chdir(RunPath)
    p = subprocess.run(["mpirun", "-np" , str(KINETIC_np), "./fp2df1_fast"])

    # #Impact to Hydro i/o    
    for file in os.listdir(RunPath):
        filename, extension = os.path.splitext(file)
        if extension == ".xy" or extension == ".xyz" or extension == ".xyv" or extension == ".xyt" or extension == ".dat" or extension == ".t":
            shutil.move(file, kinetic_dump_path)

    print("kappa")




