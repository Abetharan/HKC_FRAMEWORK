import os 
import shutil
import subprocess
import scipy.constants as constants
from string import Template
import FortGenerator as fort
import SetFort10Param as setFort10
import SetFort12Param as SetFort12
import SetHydroInit as SetHydro
import SetEnvVar as SetEnvVar
import impact_module_py3 as cf
import HydroImpactIO as io
kb  = 1.38E-23
e = 1.6E-19

def Execute(cmd):
    """ 
    Purpose: Runs a command and handles stdout and cerr

    Args: Exe comman that takes the form listed in the documetnation as a str list

    """
    try:
        subprocess.run(cmd, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        import sys
        print(e.output)
        sys.exit(1)

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
    #os.system('./fp2df1_compile_run.sh')
    os.system('./hydra_kappa.sh')


def Kinetic(runPath,_KINETIC_np):
    """  
    Purpose: Launches Impact and Sets the number of cores

    Args:
        runPath = Path where IMPACT looks for reference files
        _KINETIC_np = Number of cores being used 
    """
    
    os.chdir(runPath)
    impact_cmd = ["mpirun", "-np", str(_KINETIC_np), "./fp2df1_fast"]
    Execute(impact_cmd)

def Fluid(parameterPath,fluidSrcDir):
    """  
    Purpose: Launches Impact and Sets the number of cores

    Args:
        parameterPath = path where the fluid parameter file is located
        fluidSrcDir = path to fluid exe
    """

    headless_cmd = [fluidSrcDir,'-Vb','-p', parameterPath+'/HydroParameterInit.txt']
    Execute(headless_cmd)

def SetFluidParam(fluidNx, cq, gamma, cfl, laserWavelength,  laserPower, durOfLaser,laserLoc, 
            steps, fluidTMax, initialDt, dtGlobalMax, dtGlobalMin, outputFrequency, boundaryCondition, 
            fluidInitPath, fluidDumpPath, switchPath, FeosPathMaterial1, FeosPathMaterial2, cycleDumpPath):
        """ 
        Purpose: Handles hydro init and creation of init params textfile for HeadlessHydra.
        Args:
        fluidNx = number of grid points.
        atomicAr  = Mass Number of material
        atomicZ = Ionisation state(fully ionized gas atm so it will just be its atomic number)
        cq = Von-Neumann Richtmyer zonal spread number
        gamma = Gas constant
        cfl = PLACEHOLDER NOT IN USE ANYMORE ... DEFAULT TO 0.85
        laserWavelength = Wavelength of laser in m
        laserPower = Power of laser in W
        durOfLaser = Duration of laser in s
        laserLoc = Direction that laser strikes from if left laser goes left -> right vice versa.
        steps = number of steps.
        fluidTMax = Time maximum
        initialDt = initial step size in t
        dtGlobalMax = Max step size
        dtGlobalMin = Min step size
        outputFrequency = Frequency of dumping data files
        boundaryCondition = Boundary conditions for velocity
        fluidInitPath = Fluid Initialisation files locations
        fluidDumpPath = Fluid Dumping path
        switchPath = path to txt file containin switches.
        cycleDumpPath = Dump path for coupling cycle
        """    
        
        #Set Hydro param
        hydroparam = SetHydro.set_hydro_init(fluidNx, cq, gamma, cfl, laserWavelength,  laserPower,
                                            durOfLaser, laserLoc, steps, fluidTMax, initialDt,dtGlobalMax, dtGlobalMin, outputFrequency, 
                                            boundaryCondition, fluidInitPath, fluidDumpPath, switchPath, FeosPathMaterial1, FeosPathMaterial2) 

        # Handling templating to create the init file for fluid code
        templating(tmpfilePath = os.environ['BASEDIR'] +'/tmpHydroParameterInit.txt', writePath = cycleDumpPath, fileName = "HydroParameterInit.txt", parameters = hydroparam)

def SetFluidSwitches(cycleDumpPath, viscosityOn = "true", velocityOn = "true", heatConductionOn = "false", exchangeOn = "false",
                     bremsstrahlungOn = "false", invBremOn = "false", singleTemperatureOn = "false", mode = 'free', MultiMaterial="false", IdealGas="false"):
    """
        Purpose: Sets the fluid swithces.
        Args: 
        cycleDumpPath = Dump path of switch txt file
        viscosityOn = If viscosity is on or not 
        velocityOn = if velocity is on or not 
        heatConductionOn = if heat conduction is on or not 
        exchangeOn = if exchange is on or not 
        bremsstrahlungOn = if brem is on or not 
        invBremOn = if laser is on or not 
        
        Note: accepts ONLY string values of true/false
    """

    switches = {
                'Viscosity': viscosityOn,
                'Velocity':velocityOn ,
                'HeatConduction':heatConductionOn ,
                'Exchange':exchangeOn ,
                'Bremsstrahlung': bremsstrahlungOn ,
                'InvBremsstrahlung':invBremOn ,
                'IsothermalMode' : "false",
                'AdiabaticMode' : "false",
                'pDvWorkOff' : "true",
                'mode':mode,
                'SingleTemperature':singleTemperatureOn,
                'MultiMaterial':MultiMaterial,
                'IdealGas':IdealGas
                }
    templating(tmpfilePath = os.environ['BASEDIR'] + '/tmpFluidSwitch.txt', writePath = cycleDumpPath, fileName="HydroSwitches.txt", parameters = switches)
    
def SetKineticParam(normalised_values, _KINETIC_nv, _KINETIC_nx, _KINETIC_ny, dt, kineticTMax, cycle_dump_path, runPath, fort12Output = ["1.0d0","5.0d0","10.0d0","20.0d0", "30.0d0", "50.0d0", "100.0d0", "200.0d0"],):
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

def moveIMPACTFILE(runPath, cycleDumpPath, previousKineticInputPath, previousKineticOutputPath):
    """ 
    Purpose: Moves IMPACT files and initial parameters files to correct paths.
    Args:
        runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
        cycleDumpPath = cycle path. Path is BASEDIR/runName/cycleNAME
        previsouKineticInputPath = previous cycle kinetic input folder which is located as followeing runPath/previsousCyclePath/kinetic_input
        previousKineticOutPut =  previous cycle kinetic output folder which is located as followeing runPath/previsousCyclePath/kinetic_output
    """
    
    runName = os.environ['RUN']
    filenames = ['ionden', 'rad_to_electron', 'xf', 'eden', 'laserdep', 'tmat', 'zstar']
    
    for name in filenames:
        if os.path.splitext(cycleDumpPath + "/" + runName + "_" + name  + ".xy")[-1] != ".xy":
            continue
        shutil.move(runPath + "/" + runName + "_" + name + ".xy", previousKineticInputPath + "/" + runName + "_" + name  + ".xy", )

    for file in os.listdir(runPath):
        _, extension = os.path.splitext(file)
        if extension == ".xy" or extension == ".xyz" or extension == ".xyv" or extension == ".xyt" or extension == ".dat" or extension == ".t":
            shutil.move(file, previousKineticOutputPath)

def NextCycleFileManager(runPath, cycleStep):
    """
    Purpose: Creates the folders for the next cycle and calls any file moving required before the start of new time step. NAMELY moving IMPACT files.
    Args:
        runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
        cycleStep = cycle number.
    """

    cycle_dump_path = runPath + "cycle_" + str(cycleStep)
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
    
    if cycleStep > 0:
        previous_cycle_dump_path = runPath + "cycle_" + str(cycleStep - 1) + "/"
        previous_fluid_input_path = previous_cycle_dump_path + "/fluid_input/"
        previous_fluid_output_path = previous_cycle_dump_path + "/fluid_output/"
        previous_kinetic_output_path = previous_cycle_dump_path + "/kinetic_output/"
        previous_kinetic_input_path =  previous_cycle_dump_path + "/kinetic_input/"
        moveIMPACTFILE(runPath,cycle_dump_path, previous_kinetic_input_path, previous_kinetic_output_path)
       
        return(cycle_dump_path,fluid_input_path, fluid_output_path, kinetic_input_path, kinetic_output_path,
                previous_fluid_input_path, previous_fluid_output_path, previous_kinetic_input_path, previous_kinetic_output_path)
    else:
        return(cycle_dump_path,fluid_input_path, fluid_output_path, kinetic_input_path, kinetic_output_path,
                    0,0, 0, 0)

