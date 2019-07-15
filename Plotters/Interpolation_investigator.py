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
def interpolated_plotter(path_fluid, path_fluid_coord, path_kinetic, run_name, variable, time_step, norm_const, lambda_mfp, path, heat = False, non_standard = True):
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
    var_mat_dict = dictionary_of_info['mat']
    var_array = var_mat_dict[:] * norm_const #* 2 * normal_dict['Te'] * (e/kb) #* norm_const
    xgrid = dictionary_of_info['x_grid'] * lambda_mfp
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
            kinetic_HeatConduction[i] = - (var_array[i + 1] - var_array[i]) / mass[0]
        

        plt.plot(kinetic_grid, kinetic_HeatConduction, "k-", label = "Not Interpolated")
        plt.plot(centered_fluid_coord, fluid_HeatConduction, "r--", label = "Interpolated")

    
    else:
        plt.plot(centered_fluid_coord, search_value, "k-", label = "Not Interpolated")
        plt.plot(xgrid[1:-1], var_array[1:-1, 1], "r--", label = "Interpolated")

    plt.legend()
    plt.xlabel('grid')
    if heat:
        plt.ylabel('Heat Flow/JKg-1')
    else:
        plt.ylabel('Temperature/K')
    plt.savefig(path)
    plt.show()

def grid_plotter(path_fluid, path_kinetic, run_name, normal_dict,dump_path):

    dictionary_of_info = cf.load_dict(path_kinetic, run_name, "xf", "00", "0")
    var_mat_dict = dictionary_of_info['mat']
    var_array = var_mat_dict[:] * normal_dict['lambda_mfp']
    
    fluid_coord = np.loadtxt(path_fluid)
    
    centered_kinetic_grid = [((var_array[i + 1] + var_array[i]) / 2) for i in range(len(var_array) - 1)]
    centered_fluid_grid = [((fluid_coord[i + 1] + fluid_coord[i]) / 2) for i in range(len(fluid_coord) - 1)]
    

    
    f = plt.figure(figsize = (25, 25))
    plt.subplots_adjust(hspace =0.3)
    ax = f.add_subplot(121)
    ax2 = f.add_subplot(122)
    ax.plot(var_array, var_array, "ko", label = "Kinetic")
    ax.plot(fluid_coord, fluid_coord, "rx", label = "Fluid")
    ax.set_title("Grid position at cell walls")
    ax.set_ylabel('Grid position/m')
    ax.set_xlabel('Grid position/m')
    ax.legend()
    ax2.plot(centered_kinetic_grid, centered_kinetic_grid, "ko", label = "Kinetic")
    ax2.plot(centered_fluid_grid, centered_fluid_grid, "rx", label = "Fluid")
    ax2.set_title("Grid position at cell centers")
    ax2.set_ylabel('Grid position/m')
    ax2.set_xlabel('Grid position/m')
    ax2.legend()
    plt.savefig(dump_path)
    plt.show()




ne = 1e20
Te = 300
Z = 64
Ar = 157
Bz = 0
normal_dict = IN.impact_inputs(ne, Te, Z, Bz, Ar)

#path_fluid = "/media/abetharan/DATADRIVE1/Abetharan/fixed_nx/lin30_/cycle_1/fluid_input/"
#_base_path = "/media/abetharan/DATADRIVE1/Abetharan/data_results/fluid_varied_nx_fixed_kinetic_nx/"
_base_path = "/media/abetharan/DATADRIVE1/Abetharan/data_results/fluid_fixed_nx_varied_kinetic_nx/"
run_name = "klin60"
cycle_no = 40
cycle_name = "cycle_" + str(cycle_no)
fluid_no = 0
# #Temperature
# path_1 = os.path.join(_base_path, run_name, "cycle_" + str(cycle_no))
# path_fluid = os.path.join(path_1, "fluid_output/TemperatureE_" + str(fluid_no) + ".txt")
# path_fluid_coord = os.path.join(path_1, "fluid_output/Coord_" + str(fluid_no) + ".txt")
# path_kinetic = os.path.join(path_1, "kinetic_input/")
# #run_name = "klin60"
# variable = "tmat"
# time_step = "00"
# norm_const = 2 * normal_dict['Te'] * (e/kb)

#Heat Flow
# path_1 = os.path.join(_base_path, run_name, "cycle_" + str(cycle_no + 1))
# path_2 = os.path.join(_base_path, run_name, "cycle_" + str(cycle_no))
# path_fluid = os.path.join(path_1, "fluid_input/")
# path_fluid_coord = os.path.join(path_1, "fluid_input/coord.txt")
# path_kinetic = os.path.join(path_2, "kinetic_output/")
# run_name = "default"
# variable = "qxX"
# time_step = "03"
# norm_const = 9.11E-31 * normal_dict['vte']**3 * normal_dict['ne'] * 1e6 * 1e21

# path_fluid_coord = "/media/abetharan/DATADRIVE1/Abetharan/fixed_nx/lin30_/cycle_1/fluid_input/coord.txt"
# path_kinetic = "/media/abetharan/DATADRIVE1/Abetharan/fixed_nx/lin30_/cycle_0/kinetic_output/"

# 
# dump_path = "/home/abetharan/HYDRO_IMPACT_COUPLING_/results/Te_" + cycle_name + "_cub_fluid_out_kinetic_in_30_fluid_60_kinetic.png"

# interpolated_plotter(path_fluid, path_fluid_coord, path_kinetic, run_name, variable, time_step, norm_const, 
#                     normal_dict['lambda_mfp'], dump_path, False, True)



#grid
#Heat Flow
path_1 = os.path.join(_base_path, run_name, "cycle_" + str(cycle_no))
path_2 = os.path.join(_base_path, run_name, "cycle_" + str(cycle_no))
path_fluid = os.path.join(path_1, "fluid_output/Coord_73.txt")
path_kinetic = os.path.join(path_2, "kinetic_input/")

dump_path = "/home/abetharan/HYDRO_IMPACT_COUPLING_/results/grid_" + cycle_name + "_lin_60_fluid_30_kinetic.png"
grid_plotter(path_fluid, path_kinetic, run_name, normal_dict, dump_path)


