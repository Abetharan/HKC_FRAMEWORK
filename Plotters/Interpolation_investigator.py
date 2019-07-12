import impact_module_py3 as cf
import impact_profiles_py3 as prof
import impact_norms_py3 as IN
import matplotlib.pyplot as plt
import numpy as np
import os 
from scipy import constants

kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
def interpolated_plotter(path_fluid, path_fluid_coord, path_kinetic, run_name, variable, time_step, normal_dict, heat = False, non_standard = True):
    """ Purpose:
                Plots any variable at specificed time step with SI units. 
        Args:
            Path = path to data
            run_name = Run name i.e. name in front of variables on files
            variable = variable considered
            time_step = time step
     """

     # get dictionary containing reference IMPACT vars
    if non_standard:
        kappa = 0
    else:
        kappa = None
    dictionary_of_info = cf.load_dict(path_kinetic, run_name, variable, time_step, kappa)
    norm_const = 9.11E-31 * normal_dict['vte']**3 * normal_dict['ne'] * 1e6
    var_mat_dict = dictionary_of_info['mat']
    var_array = var_mat_dict[:] * norm_const #* 2 * normal_dict['Te'] * (e/kb) #* norm_const
    xgrid = dictionary_of_info['x_grid'] * normal_dict['lambda_mfp']
    kinetic_grid = [(xgrid[i + 1] + xgrid[i]) / 2 for i in range(len(xgrid) - 1)]
    if not heat:
        search_value = np.loadtxt(path_fluid)
    fluid_coord = np.loadtxt(path_fluid_coord)
    centered_fluid_coord = [(fluid_coord[i+ 1] + fluid_coord[i]) / 2 for i in range(len(fluid_coord) - 1)]

    if heat:
        fluid_HeatConduction = np.loadtxt(os.path.join(path_fluid, "qe.txt"))
        mass = np.loadtxt(os.path.join(path_fluid, "mass.txt"))
        nx = len(var_array) - 1
        
        kinetic_HeatConduction = np.zeros(nx)
        for i in range(nx):
            kinetic_HeatConduction[i] = - (var_array[i + 1] - var_array[i]) / mass[i]
        

        plt.plot(kinetic_grid, kinetic_HeatConduction, "k-", label = "Not Interpolated")
        plt.plot(centered_fluid_coord, fluid_HeatConduction, "r--", label = "Interpolated")

    
    else:
        plt.plot(centered_fluid_coord, search_value, "k-", label = "Not Interpolated")
        plt.plot(xgrid[1:-1], var_array[1:-1, 1], "r--", label = "Interpolated")


    plt.legend()
    plt.show()

def grid_plotter(path_fluid, path_kinetic, run_name, normal_dict):

    dictionary_of_info = cf.load_dict(path_kinetic, run_name, "xf", "00", "0")
    var_mat_dict = dictionary_of_info['mat']
    var_array = var_mat_dict[:] * normal_dict['lambda_mfp']
    
    fluid_coord = np.loadtxt(os.path.join(path_fluid, "coord.txt"))

    plt.plot(var_array, var_array, "ko", )
    plt.plot(fluid_coord, fluid_coord, "b--")
    plt.show()




ne = 1e20
Te = 300
Z = 64
Ar = 157
Bz = 0
path_fluid = "/media/abetharan/DATADRIVE1/Abetharan/fixed_nx/lin30_/cycle_1/fluid_input/"
path_fluid_coord = "/media/abetharan/DATADRIVE1/Abetharan/fixed_nx/lin30_/cycle_1/fluid_input/coord.txt"
path_kinetic = "/media/abetharan/DATADRIVE1/Abetharan/fixed_nx/lin30_/cycle_0/kinetic_output/"
run_name = "default"
variable = "qxX"
time_step = "07"
normal_dict = IN.impact_inputs(ne, Te, Z, Bz, Ar)

interpolated_plotter(path_fluid, path_fluid_coord, path_kinetic, run_name, variable, time_step, normal_dict, True, False)

path_kinetic = "/media/abetharan/DATADRIVE1/Abetharan/fixed_nx/lin60_/cycle_0/kinetic_input/"
path_fluid = "/media/abetharan/DATADRIVE1/Abetharan/fixed_nx/lin60_/cycle_0/fluid_input"
run_name = "lin60"
#grid_plotter(path_fluid, path_kinetic, run_name, normal_dict)