import impact_module_py3 as cf
import impact_profiles_py3 as prof
import impact_norms_py3 as IN
import SOL_KIT_NORMS as solnorms
import matplotlib.pyplot
matplotlib.rcParams.update({'font.size': 20})
import matplotlib.pyplot as plt
import numpy as np
import random
import os 
from scipy import constants

kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")
epsilon_0 = 8.854188E-12    # Vacuum dielectric constant

def IMPACT_plotter(path, run_name, variable, times, ne, Z, Ar, Te, Bz,iter_number = None ):
    """ Purpose:
                Plots any variable at specificed time step with SI units. 
        Args:
            Path = path to data
            run_name = Run name i.e. name in front of variables on files
            variable = variable considered
            time_step = time step
            ne = ele density in per cc
            Z = ionisatoin
            Ar = atomic number
            T = temperature in eV
            Bz = magnetic fields in tesla 
    """
    Bz = 0
    normal_dict = IN.impact_inputs(ne,Te,Z, Ar,Bz) # get dictionary containing reference IMPACT vars
    lambda_mfp = normal_dict['lambda_mfp']
    lambda_mfp_mu = lambda_mfp
    xstep_factor = lambda_mfp_mu
    tstep_factor = normal_dict['tau_ei']*1e12
    #norm_const, title, c_fmt = IN.calc_norms(variable,normal_dict,sample =0.0,forced_power=[])  
    for time_step in times:
        if variable == "Te" or variable == "tmat":
            norm_constant = normal_dict['Te'] * 2
            ylabel = r'\textbf{Temperature/eV}'
            ylabel = "Temperature/eV"
        elif variable == "n":
            norm_constant = normal_dict['ne'] * 1e6
            ylabel = r'\textbf{ne/m^-3}'
            ylabel = "ne/m^-3"

        elif variable == "Z":
            norm_constant = normal_dict["Z"]
            ylabel = r'\textbf{Z}'
            ylabel = "Z"

        elif variable == "qe":
            norm_constant = 9.11E-31 * normal_dict['vte']**3 * normal_dict['ne'] * 1e6 * 1e-6
            ylabel = r'\textbf{q_x/MWm^-2}'
            ylabel = "q_x/Wm^-2"
        elif variable == "Cx":
            norm_constant = normal_dict['vte'] #np.sqrt((2*1.38E-23 * 2 * 11604*10)/ 9.11E-31)
            ylabel = r'\textbf{Cx/ms^-1}'
            ylabel = "Cx/ms^-1"
        elif variable == "ExX":
            norm_constant = (normal_dict['lambda_mfp'] / (normal_dict['tau_ei']**2)) * (9.11E-31/1.6E-19)
            ylabel= "V/m "
        else:
            ylabel = r'\textbf{f_1}'
            norm_constant = 1
            ylabel = "f_1 x"
    
        dictionary_of_info = cf.load_dict(path,run_name,variable, time_step, iter_number)
        var_mat_dict = dictionary_of_info['mat']
        time = dictionary_of_info['time']
        iy  = int(np.shape(var_mat_dict)[-1]/2)

        if len(np.shape(var_mat_dict)) == 1:
            var_array = var_mat_dict[1:-1]  * norm_constant 
            xgrid = dictionary_of_info['x_grid'][1:-1] * xstep_factor
            xlabel = r'\textbf{x/\mum}'
            xlabel = 'x/meter'

        elif len(np.shape(var_mat_dict)) == 2:
            var_array = var_mat_dict[1:-1,1]  * norm_constant 
            xlabel = r'\textbf{x/\mum}'
            xgrid = dictionary_of_info['x_grid'][1:-1] * xstep_factor
            xlabel = 'x/meter'

        else:
            var_array = var_mat_dict[1:-1, 0, 1]  * norm_constant 
            xlabel = "v/v_th"
            xgrid = dictionary_of_info['v_grid'][1:-1]


        plt.plot(xgrid, var_array[:], label = 'IMPACT ' + str(time * tstep_factor) + '/ ps', linestyle = ':')
        plt.ylabel(ylabel)
        plt.xlabel(xlabel) 

def SOL_KIT_PLOTTER(path, variable, times, dt, Te, ne, Z, Ar, x_arr = [0]):

    normal_dict = solnorms.free_param_calc(Te, ne, Z, Ar)
    for x in x_arr:
        plt.figure(x)
        for index in times:
            lambda_mfp = normal_dict['lambda_mfp']
            tstep_factor = normal_dict['tau_ei']*1e12
            grid_x = np.loadtxt(os.path.join(path, 'GRIDS/X_GRID.txt')) *lambda_mfp
            grid_v = np.loadtxt(os.path.join(path, 'GRIDS/V_grid.txt')) * 2 * normal_dict['vte']
            if variable == "Te":
                norm_constant = normal_dict['Te']
                ylabel = r'\textbf{Temperature/eV}'
                ylabel = "Temperature/eV"
                Var = np.loadtxt(os.path.join(path, 'TEMPERATURE/TEMPERATURE_' + str(index).zfill(5) + '.txt')) * norm_constant 
            elif variable == "n":
                norm_constant = normal_dict['ne']
                ylabel = r'\textbf{ne/m^-3}'
                ylabel = "ne/m^-3"
                Var = np.loadtxt(os.path.join(path, 'DENSITY/DENSITY_' + str(index).zfill(5) + '.txt')) * norm_constant        
            elif variable == "qe":
                norm_constant = me * normal_dict['vte']**3 * normal_dict['ne'] * 1e-6
                ylabel = r'\textbf{q_x/MWm^-2}'
                ylabel = "q_x/Wm^-2"
                Var = np.loadtxt(os.path.join(path, 'HEAT_FLOW_X/HEAT_FLOW_X_' + str(index).zfill(5) + '.txt')) * norm_constant        
            elif variable == "FLOW_VEL":
                norm_constant = normal_dict['vte'] 
                ylabel = r'\textbf{Cx/ms^-1}'
                ylabel = "Cx/ms^-1"
                Var = np.loadtxt(os.path.join(path, 'FLOW_VEL_X/FLOW_VEL_X_' + str(index).zfill(5) + '.txt')) * norm_constant 
            elif variable == "ENERGY":
                ylabel = r'\textbf{E/J}'
                ylabel = "E/J"
                v = np.loadtxt(os.path.join(path, 'FLOW_VEL_X/FLOW_VEL_X_' + str(index).zfill(5) + '.txt')) * normal_dict['vte'] 
                T = np.loadtxt(os.path.join(path, 'TEMPERATURE/TEMPERATURE_' + str(index).zfill(5) + '.txt')) * normal_dict['Te'] * (e/kb)
                n = np.loadtxt(os.path.join(path, 'DENSITY/DENSITY_' + str(index).zfill(5) + '.txt')) * normal_dict['ne']        
                Var = n * kb * T + 0.5 * me * pow(v, 2)
            
            elif variable == "ohmic_heating":
                ylabel = r'\textbf{E/J}'
                ylabel = "E_ohmic/W"
                v = np.loadtxt(os.path.join(path, 'FLOW_VEL_X/FLOW_VEL_X_' + str(index).zfill(5) + '.txt')) * normal_dict['vte'] 
                n = np.loadtxt(os.path.join(path, 'DENSITY/DENSITY_' + str(index).zfill(5) + '.txt')) * normal_dict['ne']        
                J = -e * n * v
                Enorm = (normal_dict['lambda_mfp'] / (normal_dict['tau_ei']**2)) * (me/e)
                E = np.loadtxt(os.path.join(path, 'E_FIELD_X/E_FIELD_X_' + str(index).zfill(5) + '.txt')) * Enorm 
                Var = E * J

            elif variable == "sh_q":
                ylabel = r'\textbf{q_{VFP}/q_{SH}}'
                ylabel = "q_{VFP}/q_{SH}"
                Var = np.loadtxt(os.path.join(path, 'SH_q_ratio/SH_q_ratio_' + str(index).zfill(5) + '.txt'))     
            elif variable == "E":
                norm_constant = (normal_dict['lambda_mfp'] / (normal_dict['tau_ei']**2)) * (me/e)
                ylabel= "V/m "
                Var = np.loadtxt(os.path.join(path, 'E_FIELD_X/E_FIELD_X_' + str(index).zfill(5) + '.txt')) * norm_constant 
            
            elif variable == "f0":
                ylabel = r'\textbf{f_1}'
                norm_constant = 1
                Var = np.loadtxt(os.path.join(path, 'DIST_F/F_L' + str(0).rjust(3) +'_' +str(index).zfill(5) + '.txt'))
                Var = Var[:, x - 1]
                ylabel = "f_0 x"

            elif variable == "f1":
                ylabel = r'\textbf{f_1}'
                norm_constant = 1
                Var = np.loadtxt(os.path.join(path, 'DIST_F/F_L' + str(0).rjust(3) +'_' +str(index).zfill(5) + '.txt'))
                Var = Var[:, x - 1]
                ylabel = "f_1 x"

            if variable != "f0" and variable!= "f1":
                xlabel = 'Grid Position'
                plt.plot(grid_x, Var, label = 'SOL-KiT' + str(index * dt * normal_dict['tau_ei'] * 1e12) + '/ ps', linestyle = '--')
                plt.ylabel(ylabel)
                plt.xlabel(xlabel)
            else:
                xlabel = 'Energy/eV'
                plt.semilogy(grid_v, Var, label = 'SOL-KiT at ' + str(x) + 'Time step ' + str(index))
                plt.ylabel(ylabel)
                plt.xlabel(xlabel)
        plt.legend()
def ELH_1_plotter(path, var, times, dt):

    for index in times:
        time = np.loadtxt(os.path.join(path, 'TIME/TIME_' + str(index) + '.txt')) * 1e12
        
        if var not in ['qe','qi', 'v']:
            grid_x = np.loadtxt(os.path.join(path, 'CELL_CENTRE_X/CELL_CENTRE_X_' + str(index) + '.txt'))
        else:
            grid_x = np.loadtxt(os.path.join(path, 'CELL_WALL_X/CELL_WALL_X_' + str(index) + '.txt'))
        
        if var == "Te":
            ylabel = r'\textbf{Te/eV}'
            ylabel = "Temperature/eV"
            Var = np.loadtxt(os.path.join(path, 'ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_' + str(index) + '.txt')) * (kb/e)
        if var == "TI":
            ylabel = r'\textbf{Ti/eV}'
            ylabel = "Temperature/eV"
            Var = np.loadtxt(os.path.join(path, 'ION_TEMPERATURE/ION_TEMPERATURE_' + str(index) + '.txt')) * (kb/e)
        elif var == "rho":
            ylabel = r'\textbf{\rho/Kgm^-3}'
            ylabel = "rho/Kgm^-3"
            Var = np.loadtxt(os.path.join(path, 'DENSITY/DENSITY_' + str(index) + '.txt'))
        elif var == "ne":
            ylabel = r'\textbf{ne/m^-3}'
            ylabel = "ne/m^-3"
            Var = np.loadtxt(os.path.join(path, 'ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_' + str(index) + '.txt'))
        elif var == "ni":
            ylabel = r'\textbf{ni/m^-3}'
            ylabel = "ni/m^-3"
            Var = np.loadtxt(os.path.join(path, 'ION_NUMBER_DENSITY/ION_NUMBER_DENSITY_' + str(index) + '.txt'))
        elif var == "qe":
            ylabel = r'\textbf{q_e/MWm^-2}'
            ylabel = "q_x/Wm^-2"
            Var = np.loadtxt(os.path.join(path, 'ELECTRON_HEAT_FLOW_X/ELECTRON_HEAT_FLOW_X_' + str(index) + '.txt')) * 1e-6 * -1
        elif var == "qi":
            ylabel = r'\textbf{q_i/MWm^-2}'
            ylabel = "q_i/Wm^-2"
            Var = np.loadtxt(os.path.join(path, 'ION_HEAT_FLOW_X/ION_HEAT_FLOW_X_' + str(index) + '.txt')) * 1e-6 * -1
        elif var == "v":
            ylabel = r'\textbf{Velocity/ms^-1}'
            ylabel = "Velocity/ms^-1"
            Var = np.loadtxt(os.path.join(path, 'VELOCITY/VELOCITY_' + str(index) + '.txt'))
        elif var == "InvBrem":
            ylabel = r'\textbf{InvBrem/MW/kg}'
            ylabel= "Inv Brem MW/Kg"
            Var = np.loadtxt(os.path.join(path, 'INVERSE_BREM/INVERSE_BREM_' + str(index) + '.txt')) * 1e-6
        elif var == "Brem":
            ylabel = r'\textbf{Brem/MW/kg}'
            ylabel= "Brem MW/Kg"
            Var = np.loadtxt(os.path.join(path, 'BREM/BREM_' + str(index) + '.txt')) * 1e-6
 
        plt.plot(grid_x, Var, label = 'ELH-1GTR ' + str(time) + '/ps', linestyle ='-')

imp_path = "/home/abetharan/IMPACT/RUNS/high_z_epperlein_short/"
imp_name = "high_z_epperlein_short"#"BrodrickComparison
sol_path = '/Users/shiki/Documents/Imperial_College_London/Ph.D./CX1_DATA/PRE_HEAT_PROPER_INIT_HIGHER_L_OUTPUT/Run_1/OUTPUT/'
elh_1_path = '/home/abetharan/ELH1/data_out'
var = 'Te'
sol_time_step = [0, 1000]
elh_time_step =  [0, 10, 20]
ne = 1.1e27
Te = 3500
Z = 36.5
Ar = 157.0
# IMPACT_plotter(imp_path, imp_name, variable = var, time_step = '11', ne = ne*1e-6, Z = Z, Ar = Ar, Te = Te, Bz = 0)
SOL_KIT_PLOTTER(sol_path, var, sol_time_step, dt = 0.01, Te = Te, ne = ne, Z = Z, Ar = Ar, x_arr= [100])
# ELH_1_plotter(elh_1_path, var, elh_time_step, dt = 1e-15)
plt.legend()
plt.show() 