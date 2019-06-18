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
   

BASE_DIR_ = "/Users/shiki/Documents/Imperial_College_London/Ph.D./HeadlessHydra/"
RUN_NAME_ = "couple"
RUN_DIR = os.path.join(BASE_DIR_, RUN_NAME_)
_NO_CYCLES = 3
f = plt.figure(figsize = (20, 20))
total_time = 0
nx = 100
for i in range(1, _NO_CYCLES):
    cycle_path = os.path.join(RUN_DIR, "cycle_" + str(i))
    fluid_out_path = os.path.join(cycle_path, "fluid_output")
    for i in range(0,10000,500):
        _pathTe = fluid_out_path + "/TemperatureE_" + str(i) + ".txt"
        _pathTi = fluid_out_path + "/TemperatureI_" + str(i) + ".txt"
        _pathx = fluid_out_path + "/TotalPressure_" + str(i) + ".txt"
        _pathNumberDensityE =fluid_out_path + "/NumberDensityE_" + str(i) + ".txt" 
        _pathNumberDensityI =fluid_out_path + "/NumberDensityI_" + str(i) + ".txt" 

        _pathDensity = fluid_out_path + "/Density_" + str(i) + ".txt"
        _pathinvbrem = fluid_out_path + "/InverseBrem_" + str(i) + ".txt"
        _path1 = fluid_out_path + "/Coord_" + str(i) + ".txt"
        _pathTime = fluid_out_path + "/Time_" + str(i) + ".txt"
        _pathViscosity = fluid_out_path + "/Viscosity_" + str(i) + ".txt"
        _pathvelocity = fluid_out_path + "/Velocity_" + str(i) + ".txt"
        _pathHeatConductionE = fluid_out_path + "/HeatConductionE_" + str(i) + ".txt"
        _pathInternalEnergyE = fluid_out_path + "/InternalEnergyE_" + str(i) + ".txt"
        _pathInternalEnergyI = fluid_out_path + "/InternalEnergyI_" + str(i) + ".txt"
        _pathExchange =fluid_out_path + "/Exchange_" + str(i) + ".txt"
        _pathBrem =fluid_out_path + "/Brem_" + str(i) + ".txt"
        _pathPowerAbs = fluid_out_path + "/PowerAbsorbed_" + str(i) + ".txt"
        _pathTransmit = fluid_out_path + "/TransmittedLaser_" + str(i) + ".txt"
        time = np.loadtxt(_pathTime)
        rho = np.loadtxt(_pathDensity)
        velocity = np.loadtxt(_pathvelocity)
        temperatureE = np.loadtxt(_pathTe) #* (1/11600)
        temperatureI = np.loadtxt(_pathTi) #* (1/11600)
        Internal_energyE = np.loadtxt(_pathInternalEnergyE)
        Internal_energyI = np.loadtxt(_pathInternalEnergyI)
        ne = np.loadtxt(_pathNumberDensityE)
        ni = np.loadtxt(_pathNumberDensityI)
        HeatConductionE = np.loadtxt(_pathHeatConductionE)
        invbrem = np.loadtxt(_pathinvbrem)   
        Exchange = np.loadtxt(_pathExchange)
        Brem = np.loadtxt(_pathBrem)
        x = np.loadtxt(_path1)
        viscosity = np.loadtxt(_pathViscosity)
        TransmittedLaser = np.loadtxt(_pathTransmit)
        PowerAbsorbed = np.loadtxt(_pathPowerAbs)
        centered_x = [(x[i + 1] + x[i]) / 2 for i in range(nx)]

        total_time += time
        ax = f.add_subplot(211)
        ax2 = f.add_subplot(212)
        ax.plot(centered_x, temperatureE, color = 'b', label = "Electron Temperature/eV")
        #ax.plot(xgrid, Ti, color = 'r', label = "Ion Temperature/eV")
        ax.set_title("Temperature")    
        ax2.plot(centered_x, ne, color = 'b' , label = "ne/cm^-3")
        ax2.set_title("Number Density")    
        plt.suptitle("current step is: " + str(time))
        plt.title("Time " + str(time) + "S")
        plt.pause(0.0000000000000002)
        plt.gcf().clear()
    plt.show()
