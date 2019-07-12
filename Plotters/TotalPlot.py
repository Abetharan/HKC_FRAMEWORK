import numpy as np
import matplotlib.pyplot as plt



k = np.loadtxt("/media/abetharan/DATADRIVE1/Abetharan/Cub100/cycle_7/fluid_input/electron_temperature.txt")
#k1 = np.loadtxt("/home/abetharan/IMPACT/RUNS/kappa1/kappa1_qxX_01.xy")
k1 = np.loadtxt("/media/abetharan/DATADRIVE1/Abetharan/lin100/cycle_7/fluid_input/electron_temperature.txt")
#k1=  np.loadtxt("/home/abetharan/HeadlessHydra/data_out/TemperatureE_-1.txt")

plt.plot(k, 'r')
plt.plot(k1, 'k--')
plt.show()