import matplotlib.pyplot
matplotlib.rcParams.update({'font.size': 40})
import matplotlib.pyplot as plt
import numpy as np
import os
#BASE_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/"
BASE_DIR_ = "/media/abetharan/DATADRIVE1/Abetharan/Results/fixed_nx"
RUN_NAME_ = "Ncub60"
RUN_DIR = os.path.join(BASE_DIR_, RUN_NAME_)
RUN_DIR1= os.path.join(BASE_DIR_, "NcubS")
Scycle_path_00 = os.path.join(RUN_DIR1, "cycle_0/fluid_output/TemperatureE_0.txt")
Scycle_path_0 = os.path.join(RUN_DIR1, "cycle_1/fluid_output/TemperatureE_43.txt")
Scycle_path_1 = os.path.join(RUN_DIR1, "cycle_2/fluid_output/TemperatureE_43.txt")
Scycle_path_2 = os.path.join(RUN_DIR1, "cycle_3/fluid_output/TemperatureE_43.txt")
Scycle_path_3 = os.path.join(RUN_DIR1, "cycle_4/fluid_output/TemperatureE_43.txt")
Scycle_path_4 = os.path.join(RUN_DIR1, "cycle_4/fluid_output/TemperatureE_43.txt")

cycle_path_00 = os.path.join(RUN_DIR, "cycle_0/fluid_output/TemperatureE_0.txt")
cycle_path_0 = os.path.join(RUN_DIR, "cycle_1/fluid_output/TemperatureE_73.txt")
cycle_path_1 = os.path.join(RUN_DIR, "cycle_2/fluid_output/TemperatureE_73.txt")
cycle_path_2 = os.path.join(RUN_DIR, "cycle_3/fluid_output/TemperatureE_73.txt")
cycle_path_3 = os.path.join(RUN_DIR, "cycle_4/fluid_output/TemperatureE_73.txt")
cycle_path_4 = os.path.join(RUN_DIR, "cycle_4/fluid_output/TemperatureE_73.txt")

grid = os.path.join(RUN_DIR, "cycle_1/fluid_input/coord.txt")

te0 = np.loadtxt(cycle_path_00) / 11604
te10 = np.loadtxt(cycle_path_0) / 11604
te20 = np.loadtxt(cycle_path_1) / 11604 
te30 = np.loadtxt(cycle_path_2) / 11604
te40 = np.loadtxt(cycle_path_3) / 11604
te50 = np.loadtxt(cycle_path_4) / 11604

Ste0 = np.loadtxt(Scycle_path_00) / 11604
Ste10 = np.loadtxt(Scycle_path_0) / 11604
Ste20 = np.loadtxt(Scycle_path_1) / 11604 
Ste30 = np.loadtxt(Scycle_path_2) / 11604
Ste40 = np.loadtxt(Scycle_path_3) / 11604
Ste50 = np.loadtxt(Scycle_path_4) / 11604




x = np.loadtxt(grid)
centered_x = [(x[i+1] + x[i])* 1e6 / 2 for i in range(len(x) - 1)] 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

f = plt.figure(figsize = (15, 8))
plt.subplots_adjust(hspace =0.3)
plt.plot(centered_x, te0/Ste0, label = r'\textbf{0ps}')
plt.plot(centered_x, te10/Ste10, label = r'\textbf{10ps}')
plt.plot(centered_x, te20/Ste20, label = r'\textbf{20ps}')
plt.plot(centered_x, te30/Ste30, label = r'\textbf{30ps}')
plt.plot(centered_x, te40/Ste40, label = r'\textbf{40ps}')
plt.plot(centered_x, te50/Ste50, label = r'\textbf{50ps}')
#plt.title(r'\textbf{Epperlein-Short Test}')
plt.ylabel(r'\textbf{$T_{VFP}/T_{sh}$}')
plt.xlabel(r'\textbf{Position/$\mu m$}') 
plt.xticks([0,100,200,300])
plt.legend(loc=1, prop={'size': 15})
plt.grid(linestyle = '--', linewidth =0.25)
plt.show()
f.savefig("/home/abetharan/HYDRO_IMPACT_COUPLING_/results/TemperatureSpitzerVFPcomparison.pdf", bbox_inches='tight')