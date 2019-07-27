import matplotlib.pyplot
matplotlib.rcParams.update({'font.size': 40})
import matplotlib.pyplot as plt
import numpy as np
import os
import impact_module_py3 as cf
import impact_profiles_py3 as prof
import impact_norms_py3 as IN


#BASE_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/"
BASE_DIR_ = "/media/abetharan/DATADRIVE1/Abetharan/Results/fixed_nx"
RUN_NAME_ = "Ncub60"
RUN_DIR = os.path.join(BASE_DIR_, RUN_NAME_)
RUN_DIR1 = os.path.join(BASE_DIR_, "NcubS")
IMPACT_RUN_DIR = "/media/abetharan/DATADRIVE1/Abetharan/Results/non_local"

_NO_CYCLES = 3
Bz = 0
ne = 1e20
Te = 300
Z = 64
Ar = 157

normal_dict = IN.impact_inputs(ne,Te,Z, Ar,Bz) 
lambda_mfp = normal_dict['lambda_mfp']
lambda_mfp_mu = lambda_mfp*1e6
xstep_factor = lambda_mfp_mu  

cycle_path_0 = os.path.join(RUN_DIR, "cycle_0/kinetic_output/default_fxX_07.xyv")
cycle_path_1 = os.path.join(RUN_DIR, "cycle_1/kinetic_output/default_fxX_07.xyv")
cycle_path_2 = os.path.join(RUN_DIR, "cycle_2/kinetic_output/default_fxX_07.xyv")
cycle_path_3 = os.path.join(RUN_DIR, "cycle_3/kinetic_output/default_fxX_04.xyv")
cycle_path_4 = os.path.join(RUN_DIR, "cycle_4/kinetic_output/default_fxX_04.xyv")

# IMPACT_path_0 = os.path.join(RUN_DIR1, "cycle_0/kinetic_output/default_Te_01.xy")
# IMPACT_path_1 = os.path.join(RUN_DIR1, "cycle_1/kinetic_output/default_Te_01.xy")
# IMPACT_path_2 = os.path.join(RUN_DIR1, "cycle_2/kinetic_output/default_Te_01.xy")
# IMPACT_path_3 = os.path.join(RUN_DIR1, "cycle_3/kinetic_output/default_Te_01.xy")
# IMPACT_path_4 = os.path.join(RUN_DIR1, "cycle_4/kinetic_output/default_Te_01.xy")
IMPACT_path_0 = os.path.join(IMPACT_RUN_DIR, "non_local_fxX_06.xyv") # 10ps
IMPACT_path_1 = os.path.join(IMPACT_RUN_DIR, "non_local_fxX_09.xyv") # 20ps
IMPACT_path_2 = os.path.join(IMPACT_RUN_DIR, "non_local_fxX_11.xyv") # 30ps
IMPACT_path_3 = os.path.join(IMPACT_RUN_DIR, "non_local_fxX_13.xyv") # 40ps
IMPACT_path_4 = os.path.join(IMPACT_RUN_DIR, "non_local_fxX_15.xyv") # 50ps

dictionary_of_info = cf.load_dict(cycle_path_0)
dictionary_of_infoI = cf.load_dict(IMPACT_path_0)

xgrid0 = dictionary_of_info['v_grid'][1:-1]
var_mat_dict0 = dictionary_of_info['mat'][1:-1, 0 , 3]
Ixgrid0 = dictionary_of_infoI['v_grid'][1:-1]
Ivar_mat_dict0 = dictionary_of_infoI['mat'][1:-1, 0 , 3]

dictionary_of_info = cf.load_dict(cycle_path_1)
dictionary_of_infoI = cf.load_dict(IMPACT_path_1)

xgrid1 = dictionary_of_info['v_grid'][1:-1]
var_mat_dict1 = dictionary_of_info['mat'][1:-1, 0 , 3]
Ixgrid1 = dictionary_of_infoI['v_grid'][1:-1]
Ivar_mat_dict1 = dictionary_of_infoI['mat'][1:-1, 0 , 3]

dictionary_of_info = cf.load_dict(cycle_path_2)
dictionary_of_infoI = cf.load_dict(IMPACT_path_2)

xgrid2 = dictionary_of_info['v_grid'][1:-1]
var_mat_dict2 = dictionary_of_info['mat'][1:-1, 0 , 3]
Ixgrid2 = dictionary_of_infoI['v_grid'][1:-1]
Ivar_mat_dict2 = dictionary_of_infoI['mat'][1:-1, 0 , 3]

dictionary_of_info = cf.load_dict(cycle_path_3)
dictionary_of_infoI = cf.load_dict(IMPACT_path_3)

xgrid3 = dictionary_of_info['v_grid'][1:-1]
var_mat_dict3 = dictionary_of_info['mat'][1:-1, 0 , 3]
Ixgrid3 = dictionary_of_infoI['v_grid'][1:-1]
Ivar_mat_dict3 = dictionary_of_infoI['mat'][1:-1, 0 , 3]

dictionary_of_info = cf.load_dict(cycle_path_4)
dictionary_of_infoI = cf.load_dict(IMPACT_path_4)

xgrid4 = dictionary_of_info['v_grid'][1:-1]
var_mat_dict4 = dictionary_of_info['mat'][1:-1, 0 , 3]
Ixgrid4 = dictionary_of_infoI['v_grid'][1:-1]
Ivar_mat_dict4 = dictionary_of_infoI['mat'][1:-1, 0 , 3]

# te0 = np.loadtxt(cycle_path_0)
# te1 = np.loadtxt(cycle_path_1) 
# te5 = np.loadtxt(cycle_path_5)
# te9 = np.loadtxt(cycle_path_9)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

f = plt.figure(figsize = (15, 8))
plt.subplots_adjust(hspace =0.3)
plt.plot(xgrid0, var_mat_dict0,'m--', label = r'\textbf{Coupled 10ps}')
#plt.plot(xgrid1, np.log10(abs(var_mat_dict1)),'r--', label = r'\textbf{Coupled 20ps}')
plt.plot(xgrid2, var_mat_dict2,'g--', label = r'\textbf{Coupled 30ps}')
#plt.plot(xgrid3, var_mat_dict3,'c--', label = r'\textbf{Coupled 40ps}')
plt.plot(xgrid4, var_mat_dict4,'y--', label = r'\textbf{Coupled 50ps}')
plt.plot(xgrid0, Ivar_mat_dict0,'m-', label = r'\textbf{IMPACT 10ps}')
#plt.plot(xgrid1, np.log10(abs(Ivar_mat_dict1)),'r-', label = r'\textbf{IMPACT 20ps}' )
plt.plot(xgrid2, Ivar_mat_dict2,'g-', label = r'\textbf{IMPACT 30ps}' )
#plt.plot(xgrid3, Ivar_mat_dict3,'c-', label = r'\textbf{IMPACT 40ps}' )
plt.plot(xgrid4, Ivar_mat_dict4,'y-', label = r'\textbf{IMPACT 50ps}' )





plt.ylabel(r'\textbf{$f_1\vec x$}')
plt.xlabel(r'\textbf{$v / v_{th}$}') 
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.legend(loc=1, prop={'size': 15})
plt.grid(linestyle = '--', linewidth =0.25)
plt.show()
f.savefig("/home/abetharan/HYDRO_IMPACT_COUPLING_/results/3_mfp_17_um_F1_comparison.pdf", bbox_inches='tight')