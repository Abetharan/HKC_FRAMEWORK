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

cycle_path_0 = os.path.join(RUN_DIR, "cycle_0/fluid_output")
cycle_path_1 = os.path.join(RUN_DIR, "cycle_1/fluid_output")
cycle_path_2 = os.path.join(RUN_DIR, "cycle_2/fluid_output")
cycle_path_3 = os.path.join(RUN_DIR, "cycle_3/fluid_output")
cycle_path_4 = os.path.join(RUN_DIR, "cycle_4/fluid_output")
cycle_path_5 = os.path.join(RUN_DIR, "cycle_5/fluid_output")

IMPACT_path_00 = os.path.join(IMPACT_RUN_DIR, "non_local_Te_00.xy") # 0ps
IMPACT_path_0 = os.path.join(IMPACT_RUN_DIR, "non_local_Te_06.xy") # 10ps
IMPACT_path_1 = os.path.join(IMPACT_RUN_DIR, "non_local_Te_09.xy") # 20ps
IMPACT_path_2 = os.path.join(IMPACT_RUN_DIR, "non_local_Te_11.xy") # 30ps
IMPACT_path_3 = os.path.join(IMPACT_RUN_DIR, "non_local_Te_13.xy") # 40ps
IMPACT_path_4 = os.path.join(IMPACT_RUN_DIR, "non_local_Te_15.xy") # 50ps

dictionary_of_infoI = cf.load_dict(IMPACT_path_00)
Ixgrid00 = dictionary_of_infoI['x_grid'][1:-1]
Ivar_mat_dict00 = dictionary_of_infoI['mat'][1:-1, 1] * normal_dict['Te']

dictionary_of_infoI = cf.load_dict(IMPACT_path_0)
Ixgrid0 = dictionary_of_infoI['x_grid'][1:-1] 
Ivar_mat_dict0 = dictionary_of_infoI['mat'][1:-1, 1] * normal_dict['Te']

dictionary_of_infoI = cf.load_dict(IMPACT_path_1)
Ixgrid1 = dictionary_of_infoI['x_grid'][1:-1] 
Ivar_mat_dict1 = dictionary_of_infoI['mat'][1:-1, 1] * normal_dict['Te']

dictionary_of_infoI = cf.load_dict(IMPACT_path_2)
Ixgrid2 = dictionary_of_infoI['x_grid'][1:-1] 
Ivar_mat_dict2 = dictionary_of_infoI['mat'][1:-1, 1] * normal_dict['Te']

dictionary_of_infoI = cf.load_dict(IMPACT_path_3)
Ixgrid3 = dictionary_of_infoI['x_grid'][1:-1]
Ivar_mat_dict3 = dictionary_of_infoI['mat'][1:-1, 1] * normal_dict['Te']

dictionary_of_infoI = cf.load_dict(IMPACT_path_4)
Ixgrid4 = dictionary_of_infoI['x_grid'][1:-1]
Ivar_mat_dict4 = dictionary_of_infoI['mat'][1:-1, 1] * normal_dict['Te']

te0 = np.loadtxt(cycle_path_0  + "/TemperatureE_0.txt") /11604
te1 = np.loadtxt(cycle_path_1  + "/TemperatureE_73.txt")/11604 
te2 = np.loadtxt(cycle_path_2  + "/TemperatureE_73.txt")/11604
te3 = np.loadtxt(cycle_path_3  + "/TemperatureE_73.txt")/11604
te4 = np.loadtxt(cycle_path_4  + "/TemperatureE_73.txt")/11604
te5 = np.loadtxt(cycle_path_5  + "/TemperatureE_73.txt")/11604
xgrid0 = np.loadtxt(cycle_path_0 + "/Coord_0.txt") * 1e6
xgrid0 = [(xgrid0[i + 1] + xgrid0[i]) / 2 for i in range(len(xgrid0) - 1)]
xgrid1 = np.loadtxt(cycle_path_1 + "/Coord_73.txt") * 1e6 
xgrid1 = [(xgrid1[i + 1] + xgrid1[i]) / 2 for i in range(len(xgrid1) - 1)]
xgrid2 = np.loadtxt(cycle_path_2 + "/Coord_73.txt") * 1e6
xgrid2 = [(xgrid2[i + 1] + xgrid2[i]) / 2 for i in range(len(xgrid2) - 1)]
xgrid3 = np.loadtxt(cycle_path_3 + "/Coord_73.txt") * 1e6
xgrid3 = [(xgrid3[i + 1] + xgrid3[i]) / 2 for i in range(len(xgrid3) - 1)]
xgrid4 = np.loadtxt(cycle_path_4 + "/Coord_73.txt") * 1e6
xgrid4 = [(xgrid4[i + 1] + xgrid4[i]) / 2 for i in range(len(xgrid4) - 1)]
xgrid5 = np.loadtxt(cycle_path_5 + "/Coord_73.txt") * 1e6
xgrid5 = [(xgrid5[i + 1] + xgrid5[i]) / 2 for i in range(len(xgrid5) - 1)]

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

f = plt.figure(figsize = (15, 8))
plt.subplots_adjust(hspace =0.3)
plt.plot(xgrid0, te0/Ivar_mat_dict00,'m-', label = r'\textbf{$T_{Coupled}/T_{IMPACT} \quad 0ps$}')
plt.plot(xgrid1, te1/Ivar_mat_dict0, 'r-', label = r'\textbf{$T_{Coupled}/T_{IMPACT} \quad 10ps$}' )
plt.plot(xgrid2, te2/Ivar_mat_dict1, 'g-', label = r'\textbf{$T_{Coupled}/T_{IMPACT} \quad 20ps$}' )
plt.plot(xgrid3, te3/Ivar_mat_dict2, 'c-', label = r'\textbf{$T_{Coupled}/T_{IMPACT} \quad 30ps$}' )
plt.plot(xgrid4, te4/Ivar_mat_dict3, 'y-', label = r'\textbf{$T_{Coupled}/T_{IMPACT} \quad 40ps$}' )
plt.plot(xgrid5, te5/Ivar_mat_dict4, 'b-', label = r'\textbf{$T_{Coupled}/T_{IMPACT} \quad 50ps$}' )
# plt.plot(Ixgrid00 * xstep_factor,Ivar_mat_dict00 * normal_dict['Te'],'m-', label = r'\textbf{IMPACT 0ps}')
# plt.plot(Ixgrid0 * xstep_factor, Ivar_mat_dict0 * normal_dict['Te'],'m-', label = r'\textbf{IMPACT 10ps}')
# plt.plot(Ixgrid1 * xstep_factor, Ivar_mat_dict1 * normal_dict['Te'],'r-', label = r'\textbf{IMPACT 20ps}' )
# plt.plot(Ixgrid2 * xstep_factor, Ivar_mat_dict2 * normal_dict['Te'],'g-', label = r'\textbf{IMPACT 30ps}' )
# plt.plot(Ixgrid3 * xstep_factor, Ivar_mat_dict3 * normal_dict['Te'],'c-', label = r'\textbf{IMPACT 40ps}' )
# plt.plot(Ixgrid4 * xstep_factor, Ivar_mat_dict4 * normal_dict['Te'],'y-', label = r'\textbf{IMPACT 50ps}' )


plt.ylabel(r'\textbf{$T_{Coupled}/T_{IMPACT}$}')
plt.xlabel(r'\textbf{Position/$\mu m$}') 
plt.legend(loc=4, prop={'size': 15})
plt.grid(linestyle = '--', linewidth =0.25)
plt.show()
#f.savefig("/home/abetharan/HYDRO_IMPACT_COUPLING_/results/VFPCoupledTemperatureComparison.pdf", bbox_inches='tight')