import matplotlib.pyplot
matplotlib.rcParams.update({'font.size': 20})
import matplotlib.pyplot as plt
import numpy as np
import os
#BASE_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/"
BASE_DIR_ = "/media/abetharan/DATADRIVE1/Abetharan/"
RUN_NAME_ = "couple2"
RUN_DIR = os.path.join(BASE_DIR_, RUN_NAME_)
_NO_CYCLES = 3
f = plt.figure(figsize = (20, 20))
total_time = 0
nx = 100
cycle_path_0 = os.path.join(RUN_DIR, "cycle_0/fluid_output/TemperatureE_0.txt")
cycle_path_1 = os.path.join(RUN_DIR, "cycle_9/fluid_output/TemperatureE_100000.txt")
cycle_path_5 = os.path.join(RUN_DIR, "cycle_10/fluid_output/TemperatureE_100000.txt")
cycle_path_9 = os.path.join(RUN_DIR, "cycle_12/fluid_output/TemperatureE_100000.txt")

#cycle_path_9 = os.path.join(RUN_DIR, "cycle_7/fluid_output/TemperatureE_-1.txt")
#path1 = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/data_out/TemperatureE_0.txt"
#path5 = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/data_out/TemperatureE_1000.txt"
#path9 = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/data_out/TemperatureE_1999.txt"
te0 = np.loadtxt(cycle_path_0)
te1 = np.loadtxt(cycle_path_1) 
te5 = np.loadtxt(cycle_path_5)
te9 = np.loadtxt(cycle_path_9)

# te1 = np.loadtxt(path1) 
# te5 = np.loadtxt(path5)
# te9 = np.loadtxt(path9)

plt.figure(1)
plt.plot(te0/11600, label = "cycle 0")
plt.plot(te1/11600, label = "cycle 9")
plt.plot(te5/11600, label = "cycle 19")
plt.plot(te9/11600, label = "cycle 29")
plt.xlabel('Grid position')
plt.ylabel('Temperature/eV')
plt.legend()
plt.show()