import matplotlib.pyplot as plt
import numpy as np
import impact_norms_py3 as norm
import impact_profiles_py3 as prof
import impact_module_py3 as cf
import os
import pickle
import fnmatch
kb = 1.38E-23
e = 1.6E-19
protonMass = 1.67E-27
electronMass = 9.11E-31

def kineticFindLastIndex(path, var):
    largest_index = 0
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, "*_" + var + "_*"):
            k = os.path.splitext(file)[0]
            k = k.split("_")
            time_index = k[-1]
            if str.isalpha(time_index):
                continue
            if int(time_index) > largest_index:
                largest_index = int(time_index)
            
    return(largest_index)
def fluidFindLastIndex(path, var):
    largest_index = 0
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, var + "_*"):
            k = os.path.splitext(file)[0]
            k = k.split("_")
            time_index = k[-1]
            if str.isalpha(time_index):
                continue
            if int(time_index) > largest_index:
                largest_index = int(time_index)
            
    return(largest_index)
def FluidLoad(fluidOutPath, timeStep):
    xgrid = np.loadtxt(os.path.join(fluidOutPath, "Coord_" + str(timeStep) + ".txt"))     #m
    time =  np.loadtxt(os.path.join(fluidOutPath, "Time_" + str(timeStep) + ".txt"))     #s
    ne = np.loadtxt(os.path.join(fluidOutPath, "NumberDensityE_" + str(timeStep) + ".txt")) * 1e-6        #cm^-3  
    ni = np.loadtxt(os.path.join(fluidOutPath, "NumberDensityI_" + str(timeStep) + ".txt")) * 1e-6        #cm^-3  
    Te = np.loadtxt(os.path.join(fluidOutPath, "TemperatureE_" + str(timeStep) + ".txt")) * (kb/e)       #ev
    Ti = np.loadtxt(os.path.join(fluidOutPath, "TemperatureI_" + str(timeStep) + ".txt")) * (kb/e)       #ev
    HeatConductionE = np.loadtxt(os.path.join(fluidOutPath, "HeatConductionE_" + str(timeStep) + ".txt"))  #W/kg          
    #e = np.loadtxt(os.path.join(fluidOutPath, "/NumberDensityE_" + i + ".txt"))            
    #ne = np.loadtxt(os.path.join(fluidOutPath, "/NumberDensityE_" + i + ".txt"))            
    #ne = np.loadtxt(os.path.join(fluidOutPath, "/NumberDensityE_" + i + ".txt"))            
    
    return(xgrid,time, ne, ni, Te, Ti, HeatConductionE)

def KineticLoad(kineticOutPath, timeStep):

    #normalisedValues = np.load(os.path.join(kineticOutPath, "normValues.npy" )).item()
    with open(os.path.join(kineticOutPath, "normValues.pkl" ), 'rb') as f:
        normalisedValues = pickle.load(f)
    runName = "default"


    ne_var_arrays = cf.load_dict(kineticOutPath, runName, "n", str(timeStep), iter_number = None)
    ne_var_arrays_without_x = ne_var_arrays['mat'][1:-1, 1]
    normConst = 1e21 * normalisedValues['ne'] * 1e6 #to m**-3
    ne = ne_var_arrays_without_x * normConst

    Te_var_arrays = cf.load_dict(kineticOutPath, runName, "Te", str(timeStep), iter_number = None)['mat'][1:-1,1]
    normConst = normalisedValues['Te'] #To kelvin
    Te = Te_var_arrays * normConst

    qe_var_arrays = cf.load_dict(kineticOutPath, runName, "qxX", str(timeStep), iter_number = None)['mat']
    normConst = 9.11E-31 * normalisedValues['vte']**3*normalisedValues['ne'] 
    qxX = qe_var_arrays * normConst

    lambda_mfp = normalisedValues['lambda_mfp']
    lambda_mfp_mu = lambda_mfp
    xstep_factor = lambda_mfp_mu
    tstep_factor = normalisedValues['tau_ei']
    time = int(timeStep)*tstep_factor
    xgrid = ne_var_arrays['x_grid'][1:-1] * xstep_factor

    return(xgrid, time, ne, Te, qxX)
   

BASE_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HYDRO_IMPACT_COUPLING"
RUN_NAME_ = "Test_1"
RUN_DIR = os.path.join(BASE_DIR_, RUN_NAME_)
_NO_CYCLES = 3
f = plt.figure(figsize = (20, 20))
total_time = 0
for i in range(_NO_CYCLES):
    cycle_path = os.path.join(RUN_DIR, "cycle_" + str(i))
    fluid_out_path = os.path.join(cycle_path, "fluid_output")
    kinetic_out_path = os.path.join(cycle_path, "kinetic_output")
    for fluidKin in range(2):
        if fluidKin == 0:
            largest_time = fluidFindLastIndex(fluid_out_path, "TemperatureE" )
        else:
            largest_time = kineticFindLastIndex(kinetic_out_path, "Te")
        for time_step in range(largest_time):
            if fluidKin == 0:
                xgrid, time, ne, ni, Te, Ti, HeatConductionE = FluidLoad(fluid_out_path, time_step)
                xgrid = [(xgrid[i+1] + xgrid[i]) / 2 for i in range(len(ne))]
            else:
                if time_step < 10:
                    time_step = "0" + str(time_step)
                xgrid, time, ne, Te, qxX = KineticLoad(kinetic_out_path, time_step)

            total_time += time
            # ax = f.add_subplot(211)
            # ax2 = f.add_subplot(212)
            # ax.plot(xgrid, Te, color = 'b', label = "Electron Temperature/eV")
            # #ax.plot(xgrid, Ti, color = 'r', label = "Ion Temperature/eV")
            # ax.set_title("Temperature")    
            # ax2.plot(xgrid, ne, color = 'b' , label = "ne/cm^-3")
            # ax2.set_title("Number Density")    
            # plt.suptitle("current step is: " + str(time))
            plt.plot(xgrid, Te)
            #plt.legend()
            plt.title("Time " + str(time) + "S")
            plt.pause(0.0000000000000002)
            plt.gcf().clear()
        plt.show()
