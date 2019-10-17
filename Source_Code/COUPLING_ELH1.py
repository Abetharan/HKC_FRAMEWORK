import numpy as np 
import TmpFileCreator as tfc 
import subprocess
import Templating as temple 

class Fluid:

    def __init__(self, basedir, nx, laserwavelength, laserpower, 
                duroflaser, steps, fluidtmax, initialdt,
                dtglobalmax, dtglobalmin, percentageoutputfreq,
                boundaycondition, fluidsrcpath, initpath):

        self._nx = nx
        self._laser_wavelength = laserwavelength
        self._laser_power = laserpower
        self._dur_of_laser = duroflaser
        self._nt = steps,
        self._fluid_time_max = fluidtmax
        self._initial_dt = initialdt
        self._dt_global_max = dtglobalmax
        self._dt_global_min = dtglobalmin
        self._output_freq = 0
        self._fluid_run_path = fluidsrcpath
        self._init_file_path = initpath        
        self._base_dir = basedir
 
        if self._fluid_time_max == 0:
            self._output_freq = 1
        else:
            self._output_freq = int(percentageoutputfreq * 
                        self._fluid_time_max/self._dt_global_min)
 
    def setSwitches(self, viscosityon, velocityon, heatconductionon, 
                    exchangeon, bremsstrahlungon, invbremon, 
                    singletemperatureon, mode, multimaterial,
                    idealgas, fullyionized):
        
        self._viscosity = viscosityon
        self._velocity = velocityon
        self._heat_conduction = heatconductionon
        self._exchange = exchangeon
        self._bremstrahlung = bremsstrahlungon
        self._inv_brem = invbremon
        self._single_temp_mode = singletemperatureon
        self._couple_mode = mode
        self._multi_material = multimaterial
        self._ideal_gas_mode = idealgas
        self._fully_ionized_mode = fullyionized

    def makeTmpFiles(self):
        tfc.hydroParameter(self._base_dir)
        tfc.hydroSwitches(self._base_dir)


    def execute(self,parameterPath, fluidSrcDir):
        """  
        Purpose: Launches Impact and Sets the number of cores

        Args:
            parameterPath = path where the fluid parameter file is located
            fluidSrcDir = path to fluid exe
        """
        cmd = [fluidSrcDir,'-p',
                        parameterPath+'/HydroParameterInit.txt']
        try:
            subprocess.run(cmd, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            import sys
            print(e.output)
            sys.exit(1)



    def SetFluidParam(self, fluidInitPath, fluidDumpPath, switchPath, FeosPathMaterial1, FeosPathMaterial2, cycleDumpPath):
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

        # Set Hydro param
        hydroparam = SetHydro.set_hydro_init(fluidNx, cq, gamma, cfl, laserWavelength,  laserPower,
                                            durOfLaser, laserLoc, steps, fluidTMax, initialDt, dtGlobalMax, dtGlobalMin, outputFrequency,
                                            boundaryCondition, fluidInitPath, fluidDumpPath, switchPath, FeosPathMaterial1, FeosPathMaterial2)

        # Handling templating to create the init file for fluid code
        templating(tmpfilePath=os.environ['BASEDIR'] + '/tmpHydroParameterInit.txt',
                writePath=cycleDumpPath, fileName="HydroParameterInit.txt", parameters=hydroparam)

