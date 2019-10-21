import numpy as np 
import os
import fnmatch
import string
import TmpFileCreator as tfc
import SetFort10Param as sf10p
import FortGenerator as fg
from Templating import Templating
import subprocess
from Kinetic import Kinetic
from IO import IO
import shutil

class IMPACT(Kinetic):
        
    def __init__(self, IO, normalised_values_, np_, nv_, nx_, ny_, dt_, t_max_, x_max_, v_max_,
                    fort12Output=["0.2d0", "0.8d0", "1.0d0", "2.0d0", "3.0d0", "5.0d0", "10.0d0", "15.0d0"]):
        
        self._Kinetic = Kinetic()
        self._templater = Templating(IO) 
        self._run_name = IO._RUN_NAME
        self._run_path = IO._RUN_PATH
        self._base_dir = IO._BASE_DIR
        self._src_dir = IO._K_SRC_DIR
        self._cycle_path = IO._cycle_path
        self._nv = nv_
        self._nx = nx_
        self._ny = ny_
        self._dt = dt_
        self._t_max = t_max_
        self._x_max = x_max_
        self._v_max = v_max_
        self._np = np_
        self._normalised_values = normalised_values_
        self._fort12Output = fort12Output
        self.makeTmpFiles()

    def IMPACTRun(self):
        os.chdir(self._run_path)
        cmd = ["mpirun", "-np", str(self._np), "./fp2df1_fast"]    
        super().Execute(cmd)
    

    def makeTmpFiles(self):
        tfc.impactControlDummy(self._run_path)
        tfc.impactFort10(self._run_path)
        tfc.impactHeating(self._run_path)
        tfc.impactProf(self._run_path)
        tfc.impactUserCustom(self._run_path)
        tfc.impactOutputformat(self._cycle_path)

    def createFort12String(self):
        timeArray = [len(self._fort12Output)] + self._fort12Output
        strTimeArray = [str(i) for i in timeArray]
        strTimes = """ {} """.format("\n".join(strTimeArray))
        return(strTimes)

    def SetIMPACTParam(self):
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
            self._run_path = Path to where IMPACT looks for all its necessary files
        """
        # Generate fort files
        # On the fly fort10 changes as requierd of fort 10 files.
        wpe_over_nu_ei = self._normalised_values["wpe_over_nu_ei"]
        c_over_vte = self._normalised_values["c_over_vte"]
        Z = self._normalised_values["Z"]
        A = self._normalised_values["Ar"]
        fort10Param = sf10p.set_fort_10(wpe_over_nuei=wpe_over_nu_ei, c_over_vte=c_over_vte,
                                            atomic_Z=Z, atomic_A=A, nv= self._nv,nx =self._nx, ny=self._ny, dt=self._dt, tmax=self._t_max,
                                            xmax = self._x_max,vmax=self._v_max, do_user_prof_sub=".true.")

        self._templater.templating(tmpfilePath= self._base_dir + '/tmpfort.10',
                writePath=self._run_path, fileName="fort.10", parameters=fort10Param)

        fort12TimeStr = self.createFort12String()

        fg.fort_generator(self._run_path, fort12TimeStr)

    
    def IMPACTCompile(self):
        """  
        Purpose: Handles movement of neccessary template files for custom operators for IMPACT and compiles the code

        Args:
            self._run_path = Path where IMPACT looks for reference files 
        """
        self._run_path = os.environ['RUN']
        # Copy and Rename custom functions to run directory
        shutil.copyfile(os.environ["BASEDIR"] + "/heating.f",
                        self._run_path + "/" + self._run_name + "_heating.f")
        shutil.copyfile(os.environ["BASEDIR"] + "/control.dat.dummy",
                        self._run_path + "/fp2df1_control.dat.dummy")
        custom_param = {'PATH': "\'" + self._run_path +
                        "\'", 'RUNNAME': "\'" + self._run_path + "\'"}
        self._templater.templating(tmpfilePath=self._base_dir + '/user_custom.f', writePath=self._run_path,
                fileName=self._run_path + "_user_custom.f", parameters=custom_param)
        self._templater.templating(tmpfilePath=self._base_dir + '/prof.f', writePath=self._run_path,
                fileName=self._run_path + "_prof.f", parameters=custom_param)
        # Start Coupling sequence
        # os.system('./fp2df1_compile_run.sh')
        if(os.path.basename(os.path.normpath(os.getcwd())) == "Source_Code"):
            os.system('./hydra_kappa.sh')
        else:                
            os.chdir(os.getcwd() + "/Source_Code")
            os.system('./hydra_kappa.sh')


    def setEnvVar(self):
        os.environ["RUN"]= self._run_name
        print(os.environ["RUN"])
        ## laptop
        os.environ["SRCDIR"]=self._src_dir
        #Work station
        #os.environ["SRCDIR"]=os.environ["HOME"] +"/IMPACT/src"
        print(os.environ["SRCDIR"])

        if self._base_dir != None:
            os.environ["BASEDIR"]= self._base_dir
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

        os.environ["NP"]  =  str(self._np)
        os.environ["NVM"]= str(self._nv)

        os.environ["NXM"] =  str(self._nx)
        os.environ["NYM"] =  str(self._ny)
        # ...  Text output behaviour  ...
        os.environ["TEXTOUTPUT_TO_STDOUT_ROOT "]    = "1"
   
    
