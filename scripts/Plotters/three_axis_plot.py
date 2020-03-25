import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import constants
import matplotlib.pyplot as plt  
import math
kB = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")
epsilon_0 = epsilon_0 = 8.854188E-12    # Vacuum dielectric constant

BASE_DIR_ = "/home/abetharan/HYDRO_KINETIC_COUPLING/Non_Linear_Ramp_Investigation/"
RUN_NAME_ = "Pre_Heat_Ramp"
RUN_PATH = os.path.join(BASE_DIR_, RUN_NAME_)

coord = np.loadtxt(RUN_PATH + '/coord.txt')
centered_x = np.array([(coord[i+1] + coord[i]) / 2 for i in range(len(coord) - 1)]) * 1e3
Te = np.loadtxt(RUN_PATH + '/electron_temperature.txt') * (1e-3*kB/e)
density = np.loadtxt(RUN_PATH + '/density.txt')
zbar = np.loadtxt(RUN_PATH + '/Z.txt')
Ar = np.loadtxt(RUN_PATH + '/Ar.txt')
ni = density/(Ar * mp)
ne = ni * zbar * 1e-26 

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 25})
plt.rcParams.update({'xtick.major.size':10})
plt.rcParams.update({'xtick.minor.size':6})
plt.rcParams.update({'ytick.major.size':10})
plt.rcParams.update({'ytick.minor.size':6})
plt.rcParams.update({'lines.markeredgewidth':6})
plt.rcParams.update({'lines.markersize':10})
plt.rcParams.update({'legend.fontsize':20})
plt.rcParams.update({'legend.frameon': False})
plt.rcParams.update({'lines.linewidth' : 3})
fig = plt.figure()
fig.set_size_inches(13.5,5.5,forward=True)
ax = fig.add_subplot(111)
ax.set_autoscale_on(True)
#ax.set_ylim(0.002, 0.006)
ax.minorticks_on()
#ax.axis('equal')
ax.tick_params(which='major', length=10, width=2,color='gray', direction='inout', bottom = True, top= True, left= True, right = True)
ax.tick_params(which='minor', length=5, width=2, color='gray', direction='in', bottom = True, top= True, left= True, right = True)
color ='tab:blue'
ax.set_xlabel(r'\textbf{Position/$mm$}')
ax.set_ylabel(r'\textbf{Temperature/$KeV$}')
ax.set_ylim(0,3)
markers_on = [3]
lns1 = ax.plot(centered_x,Te , color=color, label = r'\textbf{Te/$KeV$}')
#ax.plot(centered_x, Te, color=color, markevery=markers_on, marker ='o', markerfacecolor ="r", markeredgecolor="r",markersize  = 15)
#ax.plot(centered_x, Te, color = color, markevery=[19], marker = 'o', markerfacecolor = "g", markeredgecolor = "g",markersize  = 15)

ax.tick_params(axis='y')
ax2 = ax.twinx()
color = "tab:black"
lns2= ax2.plot(centered_x, ne, color='k', label = r'\textbf{ne/$10^{20}cm^{-3}$}')
lns3 = ax2.plot(centered_x, zbar, color='r', label = r'\textbf{Z}')
ax2.tick_params(axis='y', labelcolor='k')
ax2.set_ylabel(r'\textbf{Z $\vert$ ne$(10^{20} cm^{-3})$}')
ax2.set_ylim(0,85)
lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=6,prop={'size': 25})
#ax.grid(linestyle = '--', linewidth =0.25)
plt.tight_layout()
#plt.show()
plt.savefig('/home/abetharan/Documents/MY_CONFERENCES_POSTER_PRESENTATIONS/HPL_CHRISTMAS_MEETING_2019/Pre_heat_init.png')
#f.savefig("/home/abetharan/HYDRO_IMPACT_COUPLING_/results/InitialProfile.pdf", bbox_inches='tight')