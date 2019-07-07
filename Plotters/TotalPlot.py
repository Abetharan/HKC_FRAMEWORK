import numpy as np
import matplotlib.pyplot as plt



k = np.loadtxt("/media/abetharan/DATADRIVE1/Abetharan/coupleS/cycle_2/fluid_input/ne.txt")
k1=  np.loadtxt("/home/abetharan/HeadlessHydra/data_out/NumberDensityE_2.txt")

plt.plot(k, 'r')
plt.plot(k1, 'k--')
plt.show()