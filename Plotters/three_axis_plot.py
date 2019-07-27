import matplotlib.pyplot
matplotlib.rcParams.update({'font.size': 40})
import matplotlib.pyplot as plt
import numpy as np
import os

BASE_DIR_ = "/media/abetharan/DATADRIVE1/Abetharan/Results/fixed_nx/"
RUN_NAME_ = "Ncub60"
RUN_DIR = os.path.join(BASE_DIR_, RUN_NAME_)
_NO_CYCLES = 3
cycle_path_0 = os.path.join(RUN_DIR, "cycle_0/fluid_input/electron_temperature.txt")
cycle_path_1 = os.path.join(RUN_DIR, "cycle_0/fluid_input/ne.txt")
cycle_path_5 = os.path.join(RUN_DIR, "cycle_0/fluid_input/Z.txt")
cycle_path_9 = os.path.join(RUN_DIR, "cycle_0/fluid_input/coord.txt")
te0 = np.loadtxt(cycle_path_0)/11604
te1 = np.loadtxt(cycle_path_1)/1e26 
te5 = np.loadtxt(cycle_path_5)
te9 = np.loadtxt(cycle_path_9)
centered_x = [(te9[i+1] + te9[i])* 1e6 / 2 for i in range(len(te9) - 1)] 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

f = plt.figure(figsize = (15, 8))
plt.subplots_adjust(hspace =0.3)


f, ax = plt.subplots(figsize = (15, 8))
color ='tab:blue'
ax.set_xlabel(r'\textbf{Position/$\mu m$}')
ax.set_ylabel(r'\textbf{Temperature/$eV$}')
markers_on = [3]
lns1 = ax.plot(centered_x, te0, color=color, label = r'\textbf{Temperature/$eV$}')
ax.plot(centered_x, te0, color=color, markevery=markers_on, marker ='o', markerfacecolor ="r", markeredgecolor="r",markersize  = 15)
ax.plot(centered_x, te0, color = color, markevery=[19], marker = 'o', markerfacecolor = "g", markeredgecolor = "g",markersize  = 15)

ax.tick_params(axis='y', labelcolor=color)
ax2 = ax.twinx()
color = "tab:black"
lns2= ax2.plot(centered_x, te1, color='k', label = r'\textbf{ne/$cm^{20}$}')
lns3 = ax2.plot(centered_x, te5, color='r', label = r'\textbf{Z}')
ax2.tick_params(axis='y', labelcolor='k')
lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=7,prop={'size': 25})
ax.grid(linestyle = '--', linewidth =0.25)
plt.show()
f.savefig("/home/abetharan/HYDRO_IMPACT_COUPLING_/results/InitialProfile.pdf", bbox_inches='tight')