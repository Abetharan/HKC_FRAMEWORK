import matplotlib.pyplot
matplotlib.rcParams.update({'font.size': 20})
import matplotlib.pyplot as plt
import numpy as np
import os
#BASE_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/"
BASE_DIR_ = "/media/abetharan/DATADRIVE1/Abetharan/"
RUN_NAME_ = "nloh5"
RUN_DIR = os.path.join(BASE_DIR_, RUN_NAME_)
_NO_CYCLES = 3
another_path = "/home/abetharan/HeadlessHydra/data_out/"
f = plt.figure(figsize = (20, 20))
total_time = 0
nx = 100
qcycle_path_0 = os.path.join(RUN_DIR, "cycle_1/fluid_input/qe.txt")
qcycle_path_1 = os.path.join(RUN_DIR, "cycle_2/fluid_input/qe.txt")
qcycle_path_5 = os.path.join(RUN_DIR, "cycle_3/fluid_input/qe.txt")
qcycle_path_9 = os.path.join(RUN_DIR, "cycle_4/fluid_input/qe.txt")
cycle_path_0 = os.path.join(RUN_DIR, "cycle_1/fluid_input/electron_temperature.txt")
cycle_path_1 = os.path.join(RUN_DIR, "cycle_2/fluid_input/electron_temperature.txt")
cycle_path_5 = os.path.join(RUN_DIR, "cycle_3/fluid_input/electron_temperature.txt")
cycle_path_9 = os.path.join(RUN_DIR, "cycle_4/fluid_input/electron_temperature.txt")

#cycle_path_9 = os.path.join(RUN_DIR, "cycle_7/fluid_output/TemperatureE_-1.txt")
#path1 = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/data_out/TemperatureE_0.txt"
#path5 = "/Users/shiki/Documents/I                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           mperial_College_London/Ph.D./HeadlessHydra/data_out/TemperatureE_1000.txt"
#path9 = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/data_out/TemperatureE_1999.txt"
te0 = np.loadtxt(cycle_path_0) / 11604
te1 = np.loadtxt(cycle_path_1)  / 11604
te5 = np.loadtxt(cycle_path_5) / 11604
te9 = np.loadtxt(cycle_path_9) / 11604
qte0 = np.loadtxt(qcycle_path_0)
qte1 = np.loadtxt(qcycle_path_1) 
qte5 = np.loadtxt(qcycle_path_5)
qte9 = np.loadtxt(qcycle_path_9)

# te1 = np.loadtxt(path1) 
# te5 = np.loadtxt(path5)
# te9 = np.loadtxt(path9)

f = plt.figure(figsize = (25, 25))
plt.subplots_adjust(hspace =0.3)
ax = f.add_subplot(211)
ax2 = f.add_subplot(212)
ax.plot(te0, label = "cycle 1")
ax.plot(te1, label = "cycle 20")
ax.plot(te5, label = "cycle 30")
ax.plot(te9, label = "cycle 40")
ax.set_title("Temperature")
ax.set_ylabel("Temperature/K")
ax.set_xlabel("Grid Position") 
ax.legend()

ax2.plot(qte0, label = "cycle 1")
ax2.plot(qte1, label = "cycle 20")
ax2.plot(qte5, label = "cycle 30")
ax2.plot(qte9, label = "cycle 40")
ax2.set_title("Heat Flow")
ax2.set_ylabel("Heat flow / Jkg-1")
ax2.set_xlabel("Grid Position") 
ax2.legend()

plt.savefig("/home/abetharan/HYDRO_IMPACT_COUPLING_/results/30_fluid_90_kinetic_Spitzer_cub.png")
#plt.savefig("/home/abetharan/HYDRO_IMPACT_COUPLING_/results/30_fluid_60_kinetic_spitzer_qe.png")
plt.show()