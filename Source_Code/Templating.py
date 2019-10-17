from string import Template


class Templating:

    def __init__(self)

    def templating(self, tmpfilePath, writePath, fileName, parameters):
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


    def SetFluidSwitches(self, cycleDumpPath, viscosityOn="true", velocityOn="false", heatConductionOn="false", exchangeOn="false",
                        bremsstrahlungOn="false", invBremOn="false", singleTemperatureOn="false", mode='free', MultiMaterial="false", IdealGas="false", FullyIonized="false"):
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
            'Velocity': velocityOn,
            'HeatConduction': heatConductionOn,
            'Exchange': exchangeOn,
            'Bremsstrahlung': bremsstrahlungOn,
            'InvBremsstrahlung': invBremOn,
            'IsothermalMode': "false",
            'AdiabaticMode': "false",
            'pDvWorkOff': "true",
            'mode': mode,
            'SingleTemperature': singleTemperatureOn,
            'MultiMaterial': MultiMaterial,
            'IdealGas': IdealGas,
            'FullyIonized':FullyIonized
        }
        templating(tmpfilePath=os.environ['BASEDIR'] + '/tmpFluidSwitch.txt',
                writePath=cycleDumpPath, fileName="HydroSwitches.txt", parameters=switches)


    def SetKineticParam(self, normalised_values, _KINETIC_nv, _KINETIC_nx, _KINETIC_ny, dt, kineticTMax, kineticXmax, kineticVmax, cycle_dump_path, runPath, fort12Output=["0.2d0", "0.8d0", "1.0d0", "2.0d0", "3.0d0", "5.0d0", "10.0d0", "15.0d0"],):
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
        # Generate fort files
        # On the fly fort10 changes as requierd of fort 10 files.
        wpe_over_nu_ei = normalised_values["wpe_over_nu_ei"]
        c_over_vte = normalised_values["c_over_vte"]
        Z = normalised_values["Z"]
        A = normalised_values["Ar"]
        fort10Param = setFort10.set_fort_10(wpe_over_nuei=wpe_over_nu_ei, c_over_vte=c_over_vte,
                                            atomic_Z=Z, atomic_A=A, nv=_KINETIC_nv, nx=_KINETIC_nx, ny=_KINETIC_ny, dt=dt, tmax=kineticTMax,
                                            xmax = kineticXmax,vmax=kineticVmax, do_user_prof_sub=".true.")

        templating(tmpfilePath=os.environ['BASEDIR'] + '/tmpfort.10',
                writePath=runPath, fileName="fort.10", parameters=fort10Param)
        fort12TimeStr = SetFort12.createFort12String(fort12Output)
        fort.fort_generator(runPath, fort12TimeStr)


