import numpy as np 
import os
import fnmatch
import string
import TmpFileCreator as tfc
import SetFort10Param
import Templating as temple
import subprocess
from Kinetic import Kinetic
from IO import IO
import shutil

class IMPACT(Kinetic):
        
    def __init__(self, IO, normalised_values_, nv_, nx_, ny_, dt_, t_max_, x_max_, v_max_):

        self._run_path = IO._RUN_PATH
        self._base_dir = IO._BASE_DIR
        self._src_dir = IO._SRC_DIR
        self._nv = nv_
        self._nx = nx_
        self._ny = ny_
        self._dt = dt_
        self._t_max = t_max_
        self._x_max = x_max_
               self._v_max = v_max_




    def IMPACTRun(self):
        os.chdir(runPath)
        cmd = ["mpirun", "-np", str(_KINETIC_np), "./fp2df1_fast"]    
        Kinetic.Execute(cmd)

    def makeTmpFiles(self):
        tfc.impactControlDummy()
        tfc.impactFort10()
        tfc.impactHeating()
        tfc.impactOutputformat()
        tfc.impactProf()
        tfc.impactUserCustom()
    
    def SetIMPACTParam(self, normalised_values, _KINETIC_nv, _KINETIC_nx, _KINETIC_ny, dt, kineticTMax, kineticXmax, kineticVmax, cycle_dump_path, runPath, fort12Output=["0.2d0", "0.8d0", "1.0d0", "2.0d0", "3.0d0", "5.0d0", "10.0d0", "15.0d0"],):
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

    def IMPACTCompile(self, runPath):
        """  
        Purpose: Handles movement of neccessary template files for custom operators for IMPACT and compiles the code

        Args:
            runPath = Path where IMPACT looks for reference files 
        """
        runName = os.environ['RUN']
        # Copy and Rename custom functions to run directory
        shutil.copyfile(os.environ["BASEDIR"] + "/heating.f",
                        runPath + "/" + runName + "_heating.f")
        shutil.copyfile(os.environ["BASEDIR"] + "/control.dat.dummy",
                        runPath + "/fp2df1_control.dat.dummy")
        custom_param = {'PATH': "\'" + runPath +
                        "\'", 'RUNNAME': "\'" + runName + "\'"}
        temple.templating(tmpfilePath=os.environ['BASEDIR'] + '/user_custom.f', writePath=runPath,
                fileName=runName + "_user_custom.f", parameters=custom_param)
        temple.templating(tmpfilePath=os.environ['BASEDIR'] + '/prof.f', writePath=runPath,
                fileName=runName + "_prof.f", parameters=custom_param)
        # Start Coupling sequence
        # os.system('./fp2df1_compile_run.sh')
        if(os.path.basename(os.path.normpath(os.getcwd())) == "Source_Code"):
            os.system('./hydra_kappa.sh')
        else:                
            os.chdir(os.getcwd() + "/Source_Code")
            os.system('./hydra_kappa.sh')


    def setEnvVar(self, nx, ny, nv, np, runName, srcDir, baseDir = None):
        os.environ["RUN"]=runName
        print(os.environ["RUN"])
        ## laptop
        os.environ["SRCDIR"]=srcDir
        #Work station
        #os.environ["SRCDIR"]=os.environ["HOME"] +"/IMPACT/src"
        print(os.environ["SRCDIR"])

        if baseDir != None:
            os.environ["BASEDIR"]= baseDir
        else:    
            os.environ["BASEDIR"]= os.getcwd()
        
        print(os.environ["BASEDIR"])


        #-----  Code generation control  -----
        os.environ["FPCODE"]     = "fp2df1_fast"
        os.environ["FPPROF"]   = ""
        os.environ["RUN_ARGS"] = "-log_summary -ksp_monitor -pc_type asm -ksp_type bcgs -ksp_rtol 1e-20 -draw_pause -1 -optionstable"

        #-----  Specify behaviour of this run script -os.environ["IS_RESTART
        # os.environ["OVERWRITE_RUN"]
        # os.environ["DONT_COMIPLE"]
        os.environ["DONT_RUN"] = ""
        # os.environ["RUN_IN_QUEUE"]  

        #------  os.environ["COMPILE-TIME ARRAY DIMENSIONS   ------
        os.environ["SYNCHRO_DIMS"] = ""
        os.environ["SMT_PK_SPM_ON "] = "1"
        os.environ["SMT_PC_COLS_ON"] = "0"

        os.environ["NP"]  =  str(np)
        os.environ["NVM"]= str(nv)

        os.environ["NXM"] =  str(nx)
        os.environ["NYM"] =  str(ny)
        # ...  Text output behaviour  ...
        os.environ["TEXTOUTPUT_TO_STDOUT_ROOT "]    = "1"


    def createFort12String(times):
        timeArray = [len(times)] + times
        strTimeArray = [str(i) for i in timeArray]
        strTimes = """ {} """.format("\n".join(strTimeArray))
        return(strTimes)
    
    def moveFile(self, runPath, cycleDumpPath, previousKineticInputPath, previousKineticOutputPath):
        """ 
        Purpose: Moves IMPACT files and initial parameters files to correct paths.
        Args:
            runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
            cycleDumpPath = cycle path. Path is BASEDIR/runName/cycleNAME
            previsouKineticInputPath = previous cycle kinetic input folder which is located as followeing runPath/previsousCyclePath/kinetic_input
            previousKineticOutPut =  previous cycle kinetic output folder which is located as followeing runPath/previsousCyclePath/kinetic_output
        """

        runName = os.environ['RUN']
        filenames = ['ionden', 'rad_to_electron',
                    'xf', 'eden', 'laserdep', 'tmat', 'zstar']

        for name in filenames:
            if os.path.splitext(cycleDumpPath + "/" + runName + "_" + name + ".xy")[-1] != ".xy":
                continue
            shutil.move(runPath + "/" + runName + "_" + name + ".xy",
                        previousKineticInputPath + "/" + runName + "_" + name + ".xy", )

        for file in os.listdir(runPath):
            _, extension = os.path.splitext(file)
            if extension == ".xy" or extension == ".xyz" or extension == ".xyv" or extension == ".xyt" or extension == ".dat" or extension == ".t":
                shutil.move(file, previousKineticOutputPath)
