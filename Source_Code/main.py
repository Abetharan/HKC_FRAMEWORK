import os
import shutil
import Coupling as cpl
import TmpFileCreator as tfc
import IMPACT_COUPLING as impactcouple

class main:
    """
    main class i.e. where all objects are passed into and the coupling is run.
    """
    def __init__(self, kinetic_nx, kinetic_ny, kinetic_nx, kinetic_np, fluid_nx, cycles)
    













# simulation domain sizes and number of processors to use
_KINETIC_nx = 60
_KINETIC_ny = 1
_KINETIC_nv = 300
_KINETIC_np = 1
_FLUID_nx = 60
_CYCLES = 51

# Material Properties
atomicZ = 64
atomicAr = 157

# Kinetic parameters
kineticDt = 0.2  # as a ratio of collisional time i.e. 1 is collision time
kineticTMax = 18  # Number of collision times
kineticXmax = 1000.0
kineticVmax = 30.0
# Fluid initial parameters
cq = 2
gamma = 1.4
cfl = 0.85
laserWavelength = 351e-9  # 200e-9
laserPower = 1e15
durOfLaser = 1e-10
laserLoc = 'left'
steps = 75
fluidTMax = 0  # 1e-15
initialDt = 1e-17
dtGlobalMax = 1e-13
dtGlobalMin = 1e-16
if fluidTMax == 0:
    outputFrequency = 1
else:
    outputFrequency = int(0.05 * fluidTMax/dtGlobalMin)

boundaryCondition = "rigid"
# Set Environement variafbles for compiling
interpolation_method = "cubic"
RUN_NAME_ = "Ncub18"
BASE_DIR_ = "/media/abetharan/DATADRIVE1/Abetharan/"
IMPACT_SRC_DIR_ = "/home/abetharan/IMPACT/src"
FLUID_SRC_DIR_ = "/home/abetharan/HeadlessHydra/Source_Code/run"
# FLUID_SRC_DIR_ = "/home/abetharan/HeadlessHydra/run"
INITIAL_CONDITIONS_FILE_PATH_ = "/media/abetharan/DATADRIVE1/Abetharan/Results/fixed_nx/Ncub60/cycle_0/fluid_input/"
#INITIAL_CONDITIONS_FILE_PATH_ = "/home/abetharan/HeadlessHydra/init_data/"

# BASE_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HYDRO_IMPACT_COUPLING/"
# IMPACT_SRC_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./IMPACT/src"
# FLUID_SRC_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/Source_Code/run"
# INITIAL_CONDITIONS_FILE_PATH_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/init_data/"
# Return path is Run directory
runPath = cpl.SetEnvVar.setEnvVar(_KINETIC_nx, _KINETIC_ny, _KINETIC_nv, _KINETIC_np, RUN_NAME_, IMPACT_SRC_DIR_,
                                  BASE_DIR_)
# If set true any conflicting files will be created
_OVERWRITE_OK_ = False

if os.path.exists(runPath):

    if _OVERWRITE_OK_:
        shutil.rmtree(runPath)
        os.makedirs(runPath)

for i in range(0, _CYCLES, 1):
    # Cretes all the directories and returns the neccessary paths
    cycle_dump_path, fluid_input_path, fluid_output_path, kinetic_input_path, kinetic_output_path, previous_fluid_input_path, previous_fluid_output_path, previous_kinetic_input_path, previous_kinetic_output_path = cpl.NextCycleFileManager(
        runPath, i)

    # Initial fluid run and genereates all tmp files that are required
    if i == 0:
        os.rmdir(fluid_input_path)
        shutil.copytree(INITIAL_CONDITIONS_FILE_PATH_, fluid_input_path)
        mode = "free"
        steps = 0
        fluidTMax = 0
        outputFrequency = round(0.05 * fluidTMax/initialDt)

        # Generates all tmp files to be changes
        tfc.hydroParameter(BASE_DIR_)
        tfc.hydroSwitches(BASE_DIR_)
        tfc.impactControlDummy(BASE_DIR_)
        tfc.impactHeating(BASE_DIR_)
        tfc.impactUserCustom(BASE_DIR_)
        tfc.impactProf(BASE_DIR_)
        tfc.impactFort10(BASE_DIR_)

    # Start coupling mode of hydro code sa well as handling IMPACT->Fluid IO
    if i > 0:
        #     cpl.io.ImpactToHydro(fluid_input_path, previous_fluid_output_path, previous_kinetic_output_path,
        #                         normalised_values, gamma, laserWavelength, laserPower, _FLUID_nx)
        cpl.io.ImpactToHydro1(normalised_values, fluid_input_path,  previous_fluid_input_path,
                              previous_fluid_output_path, previous_kinetic_output_path, interpolator=interpolation_method)
        mode = "couple"
        steps = 0
        fluidTMax = 10e-12

        outputFrequency = round(0.05 * fluidTMax/initialDt)

    # Set Switches for fluid run as well as fluid parameters
    cpl.SetFluidSwitches(cycle_dump_path,
                         heatConductionOn="true",
                         exchangeOn="false",
                         bremsstrahlungOn="false",
                         invBremOn="false",
                         singleTemperatureOn="false",
                         mode=mode,
                         MultiMaterial="false",
                         IdealGas="true",
                         FullyIonized="true"
                         )

    switchPath = cycle_dump_path + "/HydroSwitches.txt"
    cpl.SetFluidParam(_FLUID_nx, cq, gamma, cfl, laserWavelength,  laserPower, durOfLaser, laserLoc,
                      steps, fluidTMax, initialDt, dtGlobalMax, dtGlobalMin, outputFrequency, boundaryCondition,
                      fluid_input_path, fluid_output_path, switchPath, "/home/abetharan/HeadlessHydra/feos_tables/He/He", "/home/abetharan/HeadlessHydra/feos_tables/Gd/Gd", cycle_dump_path)
    # Launch hydro
    cpl.Fluid(cycle_dump_path, FLUID_SRC_DIR_)

    # Handle file transfers for hydro to impact
    normalised_values, _ = cpl.io.HydroToImpact(
        fluid_output_path, kinetic_output_path, runPath, laserWavelength, _FLUID_nx, interpolator=interpolation_method,
        normNe=1e20, normTe=300, normZ=64, normAr=157)
    # Creates fort10 file and sets the values
    cpl.SetKineticParam(normalised_values, _KINETIC_nv, _KINETIC_nx,
                        _KINETIC_ny, kineticDt, kineticTMax, kineticXmax, kineticVmax, cycle_dump_path, runPath)

    # Compile on first step this can be modified if we are doing dynamic gridding for IMPACT
    if i == 0:
        cpl.KineticCompile(runPath)

    # Launch Impact
    cpl.Kinetic(runPath, _KINETIC_np)

    # IO for last step
    if i == _CYCLES - 1:
        previous_kinetic_output_path=cycle_dump_path + "/kinetic_output/"
        previous_kinetic_input_path=cycle_dump_path + "/kinetic_input/"
        cpl.moveIMPACTFILE(runPath, cycle_dump_path,
                           previous_kinetic_input_path, previous_kinetic_output_path)
