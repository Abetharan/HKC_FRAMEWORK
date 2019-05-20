#!/bin/tcsh 
#--------------------------------------------------------------------------------
#  31/07/17   hydra_impact_testA3--redo.sh       fp2df1_v4.11.12   (Robert Kingham)
#  --------   ----------------------------       ===============
#
#--------------------------------------------------------------------------------
#  14/07/17   hydra_impact_test01--TOSS3.sh       fp2df1_v4.11.12   (Robert Kingham)
#  --------   -----------------------------       ===============
#  - Update for new TOSS-3 Linux systems  (rztrona, rztopaz)
#--------------------------------------------------------------------------------
#  10/11/16   hydra_impact_test01.sh       fp2df1_v4.11.9   (Robert Kingham)
#  --------   ----------------------       ==============
#--------------------------------------------------------------------------------

#-----  ESSENTIAL variables  -----
#set RUN="kappa"
## laptop
#set SRCDIR="${HOME}/IMPACT/src"
#set BASEDIR="`pwd`"
## cx1
##set SRCDIR="${HOME}/FP/CODE/IMPACT"
###set BASEDIR="${HOME}/FP/RUNS/IMPACT/HYDRA_NL_VALIDATION/BRODRICK_Gd_hohl--spline_interp--full_domain__rerun8"
##set BASEDIR="`pwd`"
## RZ
##set SRCDIR="${HOME}/IMPACT/src"
##set BASEDIR="/nfs/tmp2/marinak"
#set SRCDIR="${HOME}/CODE/IMPACT/src"
##set BASEDIR="/nfs/tmp2/kingham"
#set BASEDIR="${HOME}/RUNS/HYDRA_NL_VALIDATION/HYDRA--IMPACT_COUPLING"
#module load gcc/4.8-redhat mvapich2/2.2 petsc/3.7.2
#module list

echo "SRCDIR  = ${SRCDIR}"
echo "BASEDIR = ${BASEDIR}"
echo

#set RUN="${RUN}"
  
#-----  Code generation control  -----
#set FPCODE     = "fp2df1_fast"
#set FPPROF     = ""
#set RUN_ARGS   = "-log_summary -ksp_monitor -pc_type asm -ksp_type bcgs -ksp_rtol 1e-20 -draw_pause -1 -optionstable"
#set RUN_ARGS   = "-ksp_monitor -pc_type asm -ksp_type bcgs -ksp_rtol 1e-30 -draw_pause -1"


#-----  Specify behaviour of this run script -set IS_RESTART
##set OVERWRITE_RUN
##set DONT_COMIPLE
##set DONT_RUN
##set RUN_IN_QUEUE
  

#------  Set COMPILE-TIME ARRAY DIMENSIONS   ------
#set SYNCHRO_DIMS
#set SMT_PK_SPM_ON  = 1
#set SMT_PC_COLS_ON = 0
#set NP  =  36
#set NP  =  24
#set NP  =   1
#set NVM = 10
##set NXM =   6 ## Set for 36 Cores
#set NXM =   9
#set NXM =  1
#set NYM =   1
# ...  Text output behaviour  ...
#set TEXTOUTPUT_TO_STDOUT_ROOT     = 1

#-------------------------------------------------------------------------------
#
#  Meaning of array subscripts
#  ===========================
#
#  op_save_on( output_quantity ,  output_mode )   =  {  0 ==> OFF
#                                                       1 ==> ON  }
#
#  output_quantity
#  ---------------
#        1)  Z2niX
#        2)  Z2niY
#        3)  ni
#        4)  Z
#        5)  n    (i.e. ne)
#        6)  U    (i.e. Ue)
#        7)  Bz
#  8) - 11)  ExX, ExY, EyY, EyX
# 12) - 15)  jxX, jxY, jyY, jyX
# 16) - 19)  qxX, qxY, qyY, qyX
#       20)  fo
# 21) - 24)  fxX, fxY, fyX, fyY
#       25)  spmat
#       26)  cvect
#       27)  iterr
#       28)  1dvec
#       29)  Te
#       30)  wt
#       31)  Cx
#       32)  Cy
#       33)  se
#
#
#  output_mode
#  -----------
#   1)  At start
#   2)  Time monitor
#   3)  Time dumps
#   4)  Non-linear iters monitor    (not recommended for big runs)
#   5)  Non-linear iters dump       (not recommended for big runs)
#
#-------------------------------------------------------------------------------



#------------------------------------------------------------------------------+
#  COMPILE & RUN CODE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
source ${SRCDIR}/fp2df1_compile_run.sh
#exit
