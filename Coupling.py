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

def Execute(cmd):
    """ 
    Purpose: Runs a command and handles stdout and cerr

    Args: Exe comman that takes the form listed in the documetnation as a str list

    """
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        print(stdout_line )
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

def templating(tmpfilePath, writePath, fileName, parameters):
    """
    Purpose: Handles the creation and templating of files.
    
    Args:
        tmpfilePath = Path to tmp file that contains the necessary styling for templating.
        writePath = Write path for templated file.
        fileName = name of file and extension.
        parameters = substuting dictionary.
    """
       
    filein = open(tmpfilePath, "r")
    src = Template(filein.read())
    subbedtemplate = src.substitute(parameters)
    writeFile = open(writePath + "/" + fileName, "w")
    writeFile.write(subbedtemplate)
    writeFile.close()

def KineticCompile(runPath):
    """  
    Purpose: Handles movement of neccessary template files for custom operators for IMPACT and compiles the code

    Args:
        runPath = Path where IMPACT looks for reference files 
    """
    runName = os.environ['RUN']
    #Copy and Rename custom functions to run directory 
    shutil.copyfile(os.environ["BASEDIR"] + "/heating.f", runPath + "/" + runName + "_heating.f")
    shutil.copyfile(os.environ["BASEDIR"] + "/control.dat.dummy", runPath + "/fp2df1_control.dat.dummy")
    custom_param = {'PATH':"\'" + runPath + "\'", 'RUNNAME':"\'" + runName + "\'"}
    templating(tmpfilePath = os.environ['BASEDIR'] +'/user_custom.f', writePath = runPath, fileName = runName + "_user_custom.f", parameters = custom_param)   
    templating(tmpfilePath = os.environ['BASEDIR'] +'/prof.f', writePath = runPath, fileName = runName + "_prof.f", parameters = custom_param)        
    #Start Coupling sequence
    os.system('./hydra_kappa.sh')


def Kinetic(runPath):
    ##### Launch Impact
    os.chdir(runPath)
    impact_cmd = ["mpirun", "-np" , str(_KINETIC_np), "./fp2df1_fast"]
    Execute(impact_cmd)

def Fluid(parameterPath):
    #Call fluid code
    headless_cmd = [FLUID_SRC_DIR_, '-p', parameterPath+'/HydroParameterInit.txt']
    Execute(headless_cmd)

def SetFluidParam(_FLUID_nx, atomicAr, atomicZ, cq, gamma, cfl, laserWavelength,  laserPower, durOfLaser, 
            steps, fluidTMax, initialDt, dtGlobalMax, dtGlobalMin, outputFrequency, boundaryCondition,mode, cycle_init_path, fluid_dump_path, cycle_dump_path):
        """ 
        Purpose: Handles hydro init and creation of init params textfile for HeadlessHydra.
        Args:
        _FLUID_nx = number of grid points.
        atomicAr  = Mass Number of material
        atomicZ = Ionisation state(fully ionized gas atm so it will just be its atomic number)
        cq = Von-Neumann Richtmyer zonal spread number
        gamma = Gas constant
        cfl = PLACEHOLDER NOT IN USE ANYMORE ... DEFAULT TO 0.85
        laserWavelength = Wavelength of laser in m
        laserPower = Power of laser in W
        durOfLaser = Duration of laser in s
        steps = number of steps.
        fluidTMax = Time maximum
        initialDt = initial step size in t
        dtGlobalMax = Max step size
        dtGlobalMin = Min step size
        outputFrequency = Frequency of dumping data files
        boundaryCondition = Boundary conditions for velocity
        cycle_init_path = Initialisation files locations
        fluid_dump_path = Fluid Dumping path
        cycle_dump_path = Dump path for coupling cycle
        """    
        
        #Set Hydro param
        hydroparam = SetHydro.set_hydro_init(_FLUID_nx, atomicAr, atomicZ, cq, gamma, cfl, laserWavelength,  laserPower,
                                            durOfLaser, steps, fluidTMax, initialDt,dtGlobalMax, dtGlobalMin, outputFrequency, 
                                            boundaryCondition, mode,cycle_init_path, fluid_dump_path) 

        # Handling templating to create the init file for fluid code
        templating(tmpfilePath = os.environ['BASEDIR'] +'/tmpHydroParameterInit.txt', writePath = cycle_dump_path, fileName = "HydroParameterInit.txt", parameters = hydroparam)

def SetKineticParam(normalised_values, _KINETIC_nv, _KINETIC_nx, _KINETIC_ny, dt, kineticTMax, cycle_dump_path, runPath, fort12Output = ["1.0d0","2.0d0","3.0d0","5.0d0"],):
        """
        Purpose: Handles the IO side for launching IMPACT. Involves creation of fort files and moving of files and renaming. Furthermore, handles moving of files 
        when IMPACT completes.
        Args:
            normalised_values = Dict with all reference plasma values.
            _KINETIC_nv = Number of grid points in velocity space
            _KINETIC_nx = Number of grid points in X space
            _KINETIC_ny = Number of grid points in Y space
            dt = step size
            kineticTMax = max time 
            fort12Output = fort12 output times in a str list defaults to = ["1.0d0","2.0d0","3.0d0","5.0d0"]
            cycle_dump_path = path to dump all files after IMPACT launches
            runPath = Path to where IMPACT looks for all its necessary files
        """     
        #Generate fort files     
        #On the fly fort10 changes as requierd of fort 10 files.
        wpe_over_nu_ei = normalised_values["wpe_over_nu_ei"]
        c_over_vte = normalised_values["c_over_vte"]
        Z = normalised_values["Z"]
        A = normalised_values["Ar"]
        fort10Param = setFort10.set_fort_10(wpe_over_nuei = wpe_over_nu_ei, c_over_vte = c_over_vte, 
                                            atomic_Z = Z, atomic_A = A, nv = _KINETIC_nv, nx = _KINETIC_nx, ny = _KINETIC_ny, dt = dt, tmax = kineticTMax,vmax= 13, do_user_prof_sub = ".true.")

        templating(tmpfilePath = os.environ['BASEDIR'] +'/tmpfort.10', writePath = runPath, fileName = "fort.10", parameters = fort10Param)
        fort12TimeStr = SetFort12.createFort12String(fort12Output)
        fort.fort_generator(runPath, fort12TimeStr)

def moveIMPACTFILE(runPath, previousKineticInputPath, previousKineticOutputPath):
    runName = os.environ['RUN']
    filenames = ['ionden', 'rad_to_electron', 'xf', 'eden', 'laserdep', 'tmat', 'zstar']
    
    for name in filenames:
        if os.path.splitext(cycle_dump_path + "/" + runName + "_" + name  + ".xy")[-1] != ".xy":
            continue
        shutil.move(runPath + "/" + runName + "_" + name + ".xy", previousKineticInputPath + "/" + runName + "_" + name  + ".xy", )

    for file in os.listdir(runPath):
        _, extension = os.path.splitext(file)
        if extension == ".xy" or extension == ".xyz" or extension == ".xyv" or extension == ".xyt" or extension == ".dat" or extension == ".t":
            shutil.move(file, previousKineticOutputPath)

def NextCycleFileManager(runPath, CycleStep):
    cycle_dump_path = runPath + "cycle_" + str(CycleStep)
    fluid_input_path = cycle_dump_path + "/fluid_input/"
    fluid_output_path = cycle_dump_path + "/fluid_output/"
    kinetic_input_path = cycle_dump_path + "/kinetic_input/"
    kinetic_output_path = cycle_dump_path + "/kinetic_output/"
    
    if not os.path.exists(cycle_dump_path):
        os.makedirs(cycle_dump_path)
        os.makedirs(fluid_input_path)
        os.makedirs(fluid_output_path)
        os.makedirs(kinetic_input_path)
        os.makedirs(kinetic_output_path)
    if os.path.exists(cycle_dump_path):
        if not os.path.exists(fluid_input_path):
            os.makedirs(kinetic_output_path)
        if not os.path.exists(fluid_output_path):
            os.makedirs(fluid_output_path)
        if not os.path.exists(kinetic_input_path):
            os.makedirs(kinetic_output_path)
        if not os.path.exists(kinetic_output_path):
            os.makedirs(kinetic_output_path)
    
    if CycleStep > 0:
        previous_cycle_dump_path = runPath + "cycle_" + str(CycleStep - 1) + "/"
        previous_fluid_output_path = previous_cycle_dump_path + "/fluid_output/"
        previous_kinetic_output_path = previous_cycle_dump_path + "/kinetic_output/"
        previous_kinetic_input_path =  previous_cycle_dump_path + "/kinetic_input/"
        moveIMPACTFILE(runPath, previous_kinetic_input_path, previous_kinetic_output_path)
       
        return(cycle_dump_path,fluid_input_path, fluid_output_path, kinetic_input_path, kinetic_output_path,
                previous_cycle_dump_path, previous_fluid_output_path, previous_kinetic_input_path, previous_kinetic_output_path)
    else:
        return(cycle_dump_path,fluid_input_path, fluid_output_path, kinetic_input_path, kinetic_output_path,
                    0,0, 0, 0)

# simulation domain sizes and number of processors to use 
_KINETIC_nx = 30
_KINETIC_ny = 1 
_KINETIC_nv = 300
_KINETIC_np = 1
_FLUID_nx = 100
_CYCLES  = 1

#Material Properties
atomicZ = 64
atomicAr = 157

#Kinetic parameters
kineticDt = 0.8 #as a ratio of collisional time i.e. 1 is collision time 
kineticTMax = 10 #Number of collision times 


#Fluid initial parameters 
cq = 2
gamma = 1.4
cfl = 0.85
laserWavelength = 200e-9
laserPower = 0
durOfLaser = 1e-10
steps = 75
fluidTMax =  0 #1e-15
initialDt = 1e-19
dtGlobalMax =1e-15
dtGlobalMin = 1e-17
if fluidTMax == 0:
    outputFrequency = 1
else:
    outputFrequency = int(0.05 * fluidTMax/dtGlobalMin)

boundaryCondition = "rigid" 
path = "/media/abetharan/DATADRIVE1/Abetharan/Test_1/cycle_3/fluid_output/Temperature.txt"
#Set Environement variables for compiling
RUN_NAME_ = "Test_1"
BASE_DIR_ = "/media/abetharan/DATADRIVE1/Abetharan/"
IMPACT_SRC_DIR_ = "/home/abetharan/IMPACT/src"
FLUID_SRC_DIR_ = "/home/abetharan/HeadlessHydra/Source_Code/run"
INITIAL_CONDITIONS_FILE_PATH_ = "/home/abetharan/HeadlessHydra/init_data/"

#IMPACT_SRC_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./IMPACT/src"
#FLUID_SRC_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/Source_Code/run"
#INITIAL_CONDITIONS_FILE_PATH_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/init_data/"
#Return path is Run directory
runPath  = SetEnvVar.setEnvVar(_KINETIC_nx, _KINETIC_ny, _KINETIC_nv, _KINETIC_np, RUN_NAME_, IMPACT_SRC_DIR_, 
                                BASE_DIR_)

##If set true any conflicting files will be created 
_OVERWRITE_OK_ = True

if os.path.exists(runPath):
    if _OVERWRITE_OK_:
        shutil.rmtree(runPath)
        os.makedirs(runPath)   

for i in range(0, _CYCLES, 1):
    cycle_dump_path, fluid_input_path, fluid_output_path, kinetic_input_path, kinetic_output_path, previous_cycle_dump_path, previous_fluid_output_path, previous_kinetic_input_path, previous_kinetic_output_path = NextCycleFileManager(runPath, i)
    
    if i == 0:        
        os.rmdir(fluid_input_path)
        shutil.copytree(INITIAL_CONDITIONS_FILE_PATH_, fluid_input_path)
        mode = "free"
        steps = 1
    if i > 0:
        io.ImpactToHydro(fluid_input_path, previous_fluid_output_path, previous_kinetic_output_path, 
                            normalised_values, gamma, laserWavelength, laserPower, _FLUID_nx)
        mode = "couple"
        steps = 75
    SetFluidParam(_FLUID_nx, atomicAr, atomicZ, cq, gamma, cfl, laserWavelength,  laserPower, durOfLaser, 
                    steps, fluidTMax, initialDt, dtGlobalMax, dtGlobalMin, outputFrequency, boundaryCondition, mode,
                    fluid_input_path, fluid_output_path, cycle_dump_path)

    Fluid(cycle_dump_path)
    normalised_values = io.HydroToImpact(fluid_output_path, kinetic_output_path, runPath, atomicZ, atomicAr, laserWavelength, _FLUID_nx)

    SetKineticParam(normalised_values, _KINETIC_nv, _KINETIC_nx, _KINETIC_ny, kineticDt, kineticTMax, cycle_dump_path, runPath)
    if i == 0:
        KineticCompile(runPath)

    Kinetic(runPath)

    if i == _CYCLES - 1:
        previous_kinetic_output_path = cycle_dump_path + "/kinetic_output/"
        previous_kinetic_input_path =  cycle_dump_path + "/kinetic_input/"
        moveIMPACTFILE(runPath, previous_kinetic_input_path, previous_kinetic_output_path)



