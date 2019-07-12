import matplotlib.pyplot
matplotlib.rcParams.update({'font.size': 20})
import matplotlib.pyplot as plt
import numpy as np
import os
#BASE_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/"
BASE_DIR_ = "/media/abetharan/DATADRIVE1/Abetharan/"
RUN_NAME_ = "Nlin60"
RUN_DIR = os.path.join(BASE_DIR_, RUN_NAME_)
_NO_CYCLES = 3
another_path = "/home/abetharan/HeadlessHydra/data_out/"
f = plt.figure(figsize = (20, 20))
total_time = 0
nx = 100
cycle_path_0 = os.path.join(RUN_DIR, "cycle_1/fluid_input/qe.txt")
cycle_path_1 = os.path.join(RUN_DIR, "cycle_2/fluid_input/qe.txt")
cycle_path_5 = os.path.join(RUN_DIR, "cycle_3/fluid_input/qe.txt")
cycle_path_9 = os.path.join(RUN_DIR, "cycle_4/fluid_input/qe.txt")

#cycle_path_9 = os.path.join(RUN_DIR, "cycle_7/fluid_output/TemperatureE_-1.txt")
#path1 = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/data_out/TemperatureE_0.txt"
#path5 = "/Users/shiki/Documents/I                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           mperial_College_London/Ph.D./HeadlessHydra/data_out/TemperatureE_1000.txt"
#path9 = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/data_out/TemperatureE_1999.txt"
te0 = np.loadtxt(cycle_path_0)
te1 = np.loadtxt(cycle_path_1) 
te5 = np.loadtxt(cycle_path_5)
te9 = np.loadtxt(cycle_path_9)

# te1 = np.loadtxt(path1) 
# te5 = np.loadtxt(path5)
# te9 = np.loadtxt(path9)

plt.figure(1)
plt.plot()



plt.plot(te0, label = "cycle 10")
plt.plot(te1, label = "cycle 15")
plt.plot(te5, label = "cycle 25")
plt.plot(te9, label = "cycle 30")
plt.xlabel('Grid position')
plt.ylabel('Temperature/eV')
plt.legend()
plt.show()