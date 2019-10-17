import numpy as np 
import os
from scipy.interpolate import CubicSpline
from scipy import interpolate
import numpy as np
import os
import fnmatch
import string
from string import Template
from scipy import constants
import impact_norms_py3 as ImNorms
import impact_module_py3 as cf
import HydroRemainingFromImpact as remain
import pickle
import TmpFileCreator as tfc



class IMPACT:
    
    
    def __init__(self, path):

        
                            

def KineticCompile(runPath):
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
    templating(tmpfilePath=os.environ['BASEDIR'] + '/user_custom.f', writePath=runPath,
               fileName=runName + "_user_custom.f", parameters=custom_param)
    templating(tmpfilePath=os.environ['BASEDIR'] + '/prof.f', writePath=runPath,
               fileName=runName + "_prof.f", parameters=custom_param)
    # Start Coupling sequence
    # os.system('./fp2df1_compile_run.sh')
    if(os.path.basename(os.path.normpath(os.getcwd())) == "Source_Code"):
        os.system('./hydra_kappa.sh')
    else:                
        os.chdir(os.getcwd() + "/Source_Code")
        os.system('./hydra_kappa.sh')


def Kinetic(runPath, _KINETIC_np):
    """  
    Purpose: Launches Impact and Sets the number of cores

    Args:
        runPath = Path where IMPACT looks for reference files
        _KINETIC_np = Number of cores being used 
    """

    os.chdir(runPath)
    impact_cmd = ["mpirun", "-np", str(_KINETIC_np), "./fp2df1_fast"]
    Execute(impact_cmd)

    

def setEnvVar(nx, ny, nv, np, runName, srcDir, baseDir = None):
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

    return(path)
