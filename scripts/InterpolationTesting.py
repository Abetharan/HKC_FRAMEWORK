from scipy.interpolate import CubicSpline
import numpy as np 
import matplotlib.pyplot as plt
import os

nx = 10000
fluid_centered_x = np.linspace(0, np.pi, nx)
norm_x = fluid_centered_x / 1.11785e-06
ne = np.zeros(nx)
ne[:7500] = 1e25
ne[7500:] = 100e25

ne_norm = ne/(1e6 * 1e19)
kinetic_centered_x = np.linspace(norm_x[0], norm_x[-1], 300)


cs_ne = CubicSpline(norm_x, ne_norm)
kinetic_ne = cs_ne(kinetic_centered_x)
#plt.plot(ne_norm)
plt.plot(kinetic_ne)
plt.show()
