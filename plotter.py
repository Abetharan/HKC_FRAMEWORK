import matplotlib.pyplot as plt
import numpy as np
import os
BASE_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/"
RUN_NAME_ = "couple"
RUN_DIR = os.path.join(BASE_DIR_, RUN_NAME_)
_NO_CYCLES = 3
f = plt.figure(figsize = (20, 20))
total_time = 0
nx = 100
cycle_path_0 = os.path.join(RUN_DIR, "cycle_0/fluid_output/TemperatureE_-1.txt")
cycle_path_1 = os.path.join(RUN_DIR, "cycle_1/fluid_output/TemperatureE_-1.txt")
cycle_path_5 = os.path.join(RUN_DIR, "cycle_0/fluid_output/TemperatureE_-1.txt")
cycle_path_9 = os.path.join(RUN_DIR, "cycle_3/fluid_output/TemperatureE_-1.txt")
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
plt.plot(te0)
plt.plot(te1)
#plt.plot(te5)
#plt.plot(te9)
plt.show()