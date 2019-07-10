import numpy as np
import matplotlib.pyplot as plt



k = np.loadtxt("/media/abetharan/DATADRIVE1/Abetharan/couple5/cycle_2/fluid_input/qe.txt")
k = np.loadtxt("/home/abetharan/HeadlessHydra/data_out/TemperatureE_2.txt")
k1=  np.loadtxt("/home/abetharan/HeadlessHydra/data_out/TemperatureE_-1.txt")

plt.plot(k, 'r')
plt.plot(k1, 'k--')
plt.show()