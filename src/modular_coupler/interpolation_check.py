from scipy.interpolate import CubicSpline
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import impact_norms_py3 as IN
import impact_module_py3 as cf
from scipy import constants

kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
ne = 1e20
Te = 300
Z = 64
Ar = 157
Bz = 0
normal_dict = IN.impact_inputs(ne, Te, Z, Bz, Ar)
amp = 0.5
nwl = 0.5
# x = np.linspace(0, 600 * 5.75919e-07, 101)
# x1 = np.linspace(0, 600 * 5.75919e-07, 31 )
# x_centerd = np.array([(x[i+1] + x[i]) / 2 for i in range(len(x) - 1)])
# x1_centerd = np.array([(x1[i+1] + x1[i]) / 2 for i in range(len(x1) - 1)])
# l_grid = np.max(x)-np.min(x)
# xo = 0
# temperature = amp * np.sin(2*np.pi*nwl*(x1_centerd - xo)*(1.0/l_grid))
# interpolator = interpolate.interp1d(x1_centerd, temperature ,fill_value="extrapolate") #CubicSpline(x_centerd, temperature)
# interpolated_qe = interpolator(x_centerd)
# plt.plot(x1_centerd, temperature)
# plt.plot(x_centerd, interpolated_qe)
# plt.show()

path = "/media/abetharan/DATADRIVE1/Abetharan/data_results/fluid_fixed_nx_varied_kinetic_nx/klin90"
#Interpolation grid
fluid_x = np.loadtxt(path + "/cycle_0/fluid_input/coord.txt")
#Original Kinetic Grid
kinetic_x = np.loadtxt(path + "/cycle_0/kinetic_output/ReturnToHydro_xf.xy", delimiter = "\n") * normal_dict['lambda_mfp'] 

#Centering for plotting of Thermal Conduction and NOT div.q
kinetic_x_centered = np.array([(kinetic_x[i+1] + kinetic_x[i]) / 2 for i in range(len(kinetic_x) - 1)])
fluid_x_centered = np.array([(fluid_x[i+1] + fluid_x[i]) / 2 for i in range(len(fluid_x) - 1)])

# #Temperature check
# fluid_temperature = np.loadtxt(path + "/cycle_0/fluid_output/TemperatureE_73.txt")

# # interpolater = interpolate.interp1d(kinetic_x, qe, kind="linear")
# # interpolated_values = interpolater(fluid_x)
# path_kinetic = path + "/cycle_0/kinetic_input/klin60_tmat.xy"
# run_name = "klin60"
# dictionary_of_info = cf.load_dict(path_kinetic)
# Te_norm_const = normal_dict['Te'] * (e/kb)
# kinetic_temperature = dictionary_of_info['mat'][1:-1, 1] * Te_norm_const
# plt.plot(fluid_x_centered, fluid_temperature, 'k-', label = 'not interpolated')
# plt.plot(kinetic_x_centered, kinetic_temperature, 'r--', label = 'interpolated')




##Heat flow out from IMPACT 
path_kinetic = path + "/cycle_0/kinetic_output/"
run_name = "default"
dictionary_of_info = cf.load_dict(path_kinetic, run_name, "qxX", "01", None)
norm_const = 9.11E-31 * normal_dict['vte']**3 * normal_dict['ne'] * 1e6 *1e21
Te_norm_const = normal_dict['Te'] * 2
qe = dictionary_of_info['mat'] * norm_const
#Heat flow supplied to next cycle from IMPACT after interpolation/normlisation etc.
fluid_qe  = np.loadtxt(path + "/cycle_1/fluid_input/qe.txt")

interpolater = interpolate.interp1d(kinetic_x, qe, kind="linear")
interpolated_values = interpolater(fluid_x)

#mass for 30 knx
mass30 = np.loadtxt("/media/abetharan/DATADRIVE1/Abetharan/data_results/fixed_nx/lin30_/cycle_0/fluid_input/mass.txt")

#mass for 60 nx
mass60 = np.loadtxt("/media/abetharan/DATADRIVE1/Abetharan/data_results/fluid_varied_nx_fixed_kinetic_nx/lin60_/cycle_0/fluid_input/mass.txt")
knx = len(kinetic_x) - 1 

kinetic_thermal_conduction = np.zeros(knx)
for i in range(knx):
    print(qe[i])
    kinetic_thermal_conduction[i] = -(qe[i+1] - qe[i]) / mass60[i]

fnx = len(fluid_x) - 1
fluid_thermal_conduction = np.zeros(fnx)
delta = np.zeros(fnx)
for i in range(fnx):
    delta[i] = interpolated_values[i+1] - interpolated_values[i]
    fluid_thermal_conduction[i] = -(delta[i]) / mass30[i]

#plot for div.q
plt.plot(kinetic_x_centered, kinetic_thermal_conduction, "k-", label = "not interpolated")
plt.plot(fluid_x_centered, fluid_thermal_conduction, "r--", label = "interpolated")

##plotting for q
# plt.plot(kinetic_x, qe, 'k-', label = "not interpolated")
# plt.plot(fluid_x,interpolated_values, 'r--', label = "interpolated")


plt.legend()
plt.show()
