import os

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
    os.environ["NYM"] =   str(ny)
    # ...  Text output behaviour  ...
    os.environ["TEXTOUTPUT_TO_STDOUT_ROOT "]    = "1"
    path = os.environ["BASEDIR"] + "/" + os.environ["RUN"] + "/"

    return(path)
