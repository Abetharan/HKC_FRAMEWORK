import numpy as np


colision_time = 5.506e-14

time_steps = np.linspace(80,2000,80)

si_time = time_steps * colision_time
print(si_time)