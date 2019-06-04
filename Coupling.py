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

def KineticCompile(RunPath):
    """  
    Purpose: Handles movement of neccessary template files for custom operators for IMPACT and compiles the code

    Args:
        RunPath = Path where IMPACT looks for reference files 
    """
    runName = os.environ['RUN']
    #Copy and Rename custom functions to run directory 
    shutil.copyfile(os.environ["BASEDIR"] + "/heating.f", RunPath + "/" + runName + "_heating.f")
    shutil.copyfile(os.environ["BASEDIR"] + "/control.dat.dummy", RunPath + "/fp2df1_control.dat.dummy")
    custom_param = {'PATH':"\'" + RunPath + "\'", 'RUNNAME':"\'" + runName + "\'"}
    templating(tmpfilePath = os.environ['BASEDIR'] +'/user_custom.f', writePath = RunPath, fileName = runName + "_user_custom.f", parameters = custom_param)   
    templating(tmpfilePath = os.environ['BASEDIR'] +'/prof.f', writePath = RunPath, fileName = runName + "_prof.f", parameters = custom_param)        
    #Start Coupling sequence
    os.system('./hydra_kappa.sh')


def Kinetic(RunPath):
    ##### Launch Impact
    os.chdir(RunPath)
    impact_cmd = ["mpirun", "-np" , str(KINETIC_np), "./fp2df1_fast"]
    Execute(impact_cmd)

def Fluid(RunPath):
    #Call fluid code
    headless_cmd = [FLUID_SRC_DIR_, '-p', RunPath+'/HydroParameterInit.txt']
    Execute(headless_cmd)

def fluidIO(Fluid_nx, Atomic_Ar, Atomic_Z, Cq, Gamma, CFL, LaserWavelength,  LaserPower, durOfLaser, 
            steps, fluid_tmax, initialDt, dtGlobalMax, dtGlobalMin, OutputFrequency, BoundaryCondition, cycle_init_path, fluid_dump_path, cycle_dump_path):
        """ 
        Purpose: Handles hydro init and creation of init params textfile for HeadlessHydra.
        Args:
        Fluid_nx = number of grid points.
        Atomic_Ar  = Mass Number of material
        Atomic_Z = Ionisation state(fully ionized gas atm so it will just be its atomic number)
        Cq = Von-Neumann Richtmyer zonal spread number
        Gamma = Gas constant
        CFL = PLACEHOLDER NOT IN USE ANYMORE ... DEFAULT TO 0.85
        LaserWavelength = Wavelength of laser in m
        LaserPower = Power of laser in W
        durOfLaser = Duration of laser in s
        steps = number of steps.
        fluid_tmax = Time maximum
        initialDt = initial step size in t
        dtGlobalMax = Max step size
        dtGlobalMin = Min step size
        OutputFrequency = Frequency of dumping data files
        BoundaryCondition = Boundary conditions for velocity
        cycle_init_path = Initialisation files locations
        fluid_dump_path = Fluid Dumping path
        cycle_dump_path = Dump path for coupling cycle
        """    
        
        #Set Hydro param
        hydroparam = SetHydro.set_hydro_init(Fluid_nx, Atomic_Ar, Atomic_Z, Cq, Gamma, CFL, LaserWavelength,  LaserPower,
                                            durOfLaser, steps, fluid_tmax, initialDt,dtGlobalMax, dtGlobalMin, OutputFrequency, 
                                            BoundaryCondition, cycle_init_path, fluid_dump_path) 

        # Handling templating to create the init file for fluid code
        templating(tmpfilePath = os.environ['BASEDIR'] +'/tmpHydroParameterInit.txt', writePath = cycle_dump_path, fileName = "HydroParameterInit.txt", parameters = hydroparam)

def kineticIO(wpe_over_nuei, c_over_vte, Z, A, KINETIC_nv, KINETIC_nx, KINETIC_ny, dt, kinetic_tmax, fort12Output = ["1.0d0","2.0d0","3.0d0","5.0d0"],
                cycle_dump_path = False, RunPath = False, file_move = False,):
        """
        Purpose: Handles the IO side for launching IMPACT. Involves creation of fort files and moving of files and renaming. Furthermore, handles moving of files 
        when IMPACT completes.
        Args:
            wpe_over_nuei = IMPACT launch parameter
            c_over_vte = IMPACT launch parameter
            Z = Atomiz number
            A = Mass Number
            KINETIC_nv = Number of grid points in velocity space
            KINETIC_nx = Number of grid points in X space
            KINETIC_ny = Number of grid points in Y space
            dt = step size
            kinetic_tmax = max time 
            fort12Output = fort12 output times in a str list defaults to = ["1.0d0","2.0d0","3.0d0","5.0d0"],
            cycle_dump_path = path to dump all files after IMPACT launches. Defaults to  False
            RunPath = Path to where IMPACT looks for all its necessary files. Defaults to  False
            file_move = Whether IMPACT output files are moved or not. Defaults to False
        """     
        #Generate fort files     
        #On the fly fort10 changes as requierd of fort 10 files.
        fort10Param = setFort10.set_fort_10(wpe_over_nuei = wpe_over_nuei, c_over_vte = c_over_vte, 
                                            atomic_Z = Z, atomic_A = A, nv = KINETIC_nv, nx = KINETIC_nx, ny = KINETIC_ny, dt = dt, tmax = kinetic_tmax, do_user_prof_sub = ".true.")

        templating(tmpfilePath = os.environ['BASEDIR'] +'/tmpfort.10', writePath = RunPath, fileName = "fort.10", parameters = fort10Param)
        fort12TimeStr = SetFort12.createFort12String(fort12Output)
        fort.fort_generator(RunPath, fort12TimeStr)
    
        runName = os.environ['RUN']
        filenames = ['ionden', 'rad_to_electron', 'xf', 'eden', 'laserdep', 'tmat', 'zstar']
        
        for name in filenames:
            if os.path.splitext(cycle_dump_path + "/" + runName + "_" + name  + ".xy")[-1] != ".xy":
                continue
            shutil.move(cycle_dump_path + "/" + runName + "_" + name  + ".xy", RunPath + "/" + runName + "_" + name + ".xy" )
    
        if file_move:

            # #Impact to Hydro i/o    
            for file in os.listdir(RunPath):
                _, extension = os.path.splitext(file)
                if extension == ".xy" or extension == ".xyz" or extension == ".xyv" or extension == ".xyt" or extension == ".dat" or extension == ".t":
                    shutil.move(file, kinetic_dump_path)

def NextCycleFileManager(RunPath, CycleStep):
    cycle_dump_path = RunPath + "cycle_" + str(CycleStep)
    fluid_input_path = cycle_dump_path + "/fluid_input"
    fluid_output_path = cycle_dump_path + "/fluid_output"
    kinetic_input_path = cycle_dump_path + "/kinetic_input"
    kinetic_output_path = cycle_dump_path + "/kinetic_output"
    
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
    
    previous_cycle_dump_path = RunPath + "cycle_" + str(CycleStep - 1)
    previous_fluid_output_path = previous_cycle_dump_path + "/fluid_output"
    previous_kinetic_output_path = previous_cycle_dump_path + "/kinetic_output"

    return(cycle_dump_path,fluid_input_path, fluid_output_path, kinetic_input_path, kinetic_output_path,
             previous_cycle_dump_path, previous_fluid_output_path, previous_kinetic_output_path)




# simulation domain sizes and number of processors to use 
KINETIC_nx = 30
KINETIC_ny = 1 
KINETIC_nv = 300
KINETIC_np = 1
dt = 1 #as a ratio of collisional time i.e. 1 is collision time 
kinetic_tmax = 20 #Number of collision times 
Atomic_Z = 1#60
Atomic_Ar = 1#157

#Fluid initial parameters 
Cq = 2
Gamma = 1.4
CFL = 0.85
LaserWavelength = 200e-9
LaserPower = 1e18
durOfLaser = 1e-10
steps = 0
fluid_tmax = 1e-15
initialDt = 1e-19
OutputFrequency = 1000
dtGlobalMax =1e-15
dtGlobalMin = 1e-17
OutputFrequency = 0.05 * fluid_tmax/dtGlobalMin
BoundaryCondition = "free" 
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
RunPath  = SetEnvVar.setEnvVar(KINETIC_nx, KINETIC_ny, KINETIC_nv, KINETIC_np, runName, IMPACT_SRC_DIR_, 
                                baseDir ="/media/abetharan/DATADRIVE1/Abetharan/HeadlessHydra/" )

##If set true any conflicting files will be created 
_OVERWRITE_OK_ = True

if os.path.exists(RunPath):
    if _OVERWRITE_OK_:
        shutil.rmtree(RunPath)
        os.makedirs(RunPath)   

#Data location save
no_cycles  = 3
Fluid_nx = 100

##Set on the fly Fluid parameters
shutil.copytree(init, cycle_init_path)

normalised_values = io.HydroToImpact(fluid_dump_path, RunPath, Atomic_Z, Atomic_Ar, LaserWavelength, Fluid_nx)



for i in range(0, no_cycles, 1):
   cycle_dump_path,fluid_input_path, fluid_output_path, kinetic_input_path, kinetic_output_path,previous_cycle_dump_path, previous_fluid_output_path, previous_kinetic_output_path = NextCycleFileManager(RunPath, i)
   
   if i > 0:
       io.ImpactToHydro(fluid_input_path, previous_fluid_output_path, previous_kinetic_output_path, normalised_values, Gamma, LaserWavelength, LaserPower, Fluid_nx)

   
   
   
    ##Impact Reference plasma. 
    ne = nc #1.0e+19 #Note that NE is in cm^3 i.e. times M^3 by 10^-6 to convert 
    Te = norm_Te #eV
    Z = Atomic_Z 
    Bz = 0.0 #Tesla
    Ar = Atomic_Ar 
    wpe_over_nuei = normalised_values['wpe_over_nu_ei']
    c_over_vte = normalised_values['c_over_vte']

    if normalised_values["log_lambda"] < 0:
        print("Normalisation values not within good IMPACT PARAMETERS ... change")
        exit(1)
