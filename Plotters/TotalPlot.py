import matplotlib.pyplot
matplotlib.rcParams.update({'font.size': 40})
import matplotlib.pyplot as plt
import numpy as np
import os
#BASE_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/"
# BASE_DIR_ = '/media/abetharan/DATADRIVE2/Abetharan/Coupled_data/COUPLING_METHODOLOGY/EPPERLEIN_SHORT/EPPERLEIN_SHORT_COUPLED_02/'
BASE_DIR_ = '/media/abetharan/DATADRIVE2/Abetharan/'
RUN_NAME_ = "1_KNT_3_FNT"
RUN_NAME_1 = "1_KNT_0.25_KFNT"
RUN_DIR = os.path.join(BASE_DIR_, RUN_NAME_)
RUN_DIR1= os.path.join(BASE_DIR_, RUN_NAME_1)
Scycle_path_00 = os.path.join(RUN_DIR1,"CYCLE_0/FLUID_OUTPUT/ELECTRON_TEMPERATURE/Te_0.txt")
Scycle_path_0 = os.path.join(RUN_DIR1, "CYCLE_1/FLUID_OUTPUT/ELECTRON_TEMPERATURE/Te_398.txt")
Scycle_path_1 = os.path.join(RUN_DIR1, "CYCLE_2/FLUID_OUTPUT/ELECTRON_TEMPERATURE/Te_398.txt")
Scycle_path_2 = os.path.join(RUN_DIR1, "CYCLE_3/FLUID_OUTPUT/ELECTRON_TEMPERATURE/Te_398.txt")
Scycle_path_3 = os.path.join(RUN_DIR1, "CYCLE_4/FLUID_OUTPUT/ELECTRON_TEMPERATURE/Te_398.txt")
Scycle_path_4 = os.path.join(RUN_DIR1, "CYCLE_4/FLUID_OUTPUT/ELECTRON_TEMPERATURE/Te_398.txt")

cycle_path_00 = os.path.join(RUN_DIR,"CYCLE_0/FLUID_OUTPUT/ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_0.txt")
cycle_path_0 = os.path.join(RUN_DIR, "CYCLE_1/FLUID_OUTPUT/ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_3179.txt")
cycle_path_1 = os.path.join(RUN_DIR, "CYCLE_2/FLUID_OUTPUT/ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_3179.txt")
cycle_path_2 = os.path.join(RUN_DIR, "CYCLE_3/FLUID_OUTPUT/ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_3179.txt")
cycle_path_3 = os.path.join(RUN_DIR, "CYCLE_4/FLUID_OUTPUT/ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_3179.txt")
# cycle_path_4 = os.path.join(RUN_DIR, "CYCLE_6/FLUID_OUTPUT/ELECTRON_TEMPERATURE/Te_3179.txt")
#cycle_path_4 = os.path.join(RUN_DIR, "CYCLE_6/FLUID_OUTPUT/ELECTRON_TEMPERATURE/Te_795.txt")
grid = os.path.join(RUN_DIR, "CYCLE_0/FLUID_OUTPUT/CELL_CENTRE_X/CELL_CENTRE_X_0.txt")

te0 = np.loadtxt(cycle_path_00) / 11604
te10 = np.loadtxt(cycle_path_0) / 11604
te20 = np.loadtxt(cycle_path_1) / 11604 
te30 = np.loadtxt(cycle_path_2) / 11604
te40 = np.loadtxt(cycle_path_3) / 11604
# te50 = np.loadtxt(cycle_path_4) / 11604

# Ste0 = np.loadtxt(Scycle_path_00) / 11604
# Ste10 = np.loadtxt(Scycle_path_0) / 11604
# Ste20 = np.loadtxt(Scycle_path_1) / 11604 
# Ste30 = np.loadtxt(Scycle_path_2) / 11604
# Ste40 = np.loadtxt(Scycle_path_3) / 11604
# Ste50 = np.loadtxt(Scycle_path_4) / 11604




centered_x = np.loadtxt(grid)
#centered_x = [(x[i+1] + x[i])* 1e6 / 2 for i in range(len(x) - 1)] 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

f = plt.figure(figsize = (15, 8))
plt.subplots_adjust(hspace =0.3)
plt.plot(centered_x, te0, label = r'\textbf{0ps}')
plt.plot(centered_x, te10, label = r'\textbf{te10ps}')
plt.plot(centered_x, te20, label = r'\textbf{te20ps}')
plt.plot(centered_x, te30, label = r'\textbf{t30}')
plt.plot(centered_x, te40, label = r'\textbf{te40}')
# plt.plot(centered_x, te50, label = r'\textbf{te50}')
#plt.title(r'\textbf{Epperlein-Short Test}')
plt.ylabel(r'\textbf{$T_{VFP}/T_{sh}$}')
plt.xlabel(r'\textbf{Position/$\mu m$}') 
#plt.xticks([0,100,200,300])
plt.legend(loc=1, prop={'size': 15})
plt.grid(linestyle = '--', linewidth =0.25)
plt.show()
#f.savefig("/home/abetharan/HYDRO_IMPACT_COUPLING_/results/TemperatureSpitzerVFPcomparison.pdf", bbox_inches='tight')