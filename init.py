""" Only Valid for idle gases 
    By default initial velocity is always zero
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
#from ES_Investigation import Kinetic_param_finder as kpf
#global constants
kb = 1.38064852e-23;
protonMass = 1.6726219e-27;
electronMass = 9.10938356e-31;
epsilon_0 = 8.85418782e-12;
electronCharge = 1.60217662e-19;
hbar = 1.054571800e-34;
c = 299792458;
ev_converter = electronCharge/kb
def CalculateMass(coord, density, nx):
    mass = np.zeros(nx)
    for i in range(nx):

        mass[i] = density[i]  * (coord[i + 1] - coord[i])
    return(mass)

def TextDump(path, coord, velocity, density, Te, Ti, mass, Z, Ar, qe = None):
    """ 
    Purpose: Dumps variables to text
    Args: All physical quantities and Path to dump
    Returns: Null
    """
    np.savetxt((path+"/coord.txt"), coord)
    np.savetxt((path+"/velocity.txt"), velocity)
    np.savetxt((path+"/density.txt"), density)
    np.savetxt((path+"/electron_temperature.txt"), Te)
    np.savetxt((path+"/ion_temperature.txt"), Ti)
    np.savetxt((path+"/mass.txt"), mass)
    np.savetxt((path+"/Z.txt"), Z)    
    np.savetxt((path+"/Ar.txt"), Ar)
    if qe is not None:
        np.savetxt((path+"/qe.txt"), qe)

def epperlein_short(nx, L, Z_ = 37.25, Ar_ = 157.25, ne_ = 1e27, Te_ = 100.):

    x_l = 0
    x_u = L
    initial_coord = np.linspace(x_l, x_u, nx+1, dtype=np.float64)
    centered_coord = np.array([(initial_coord[i] + initial_coord[i+1]) /2 for i in range(len(initial_coord) -1)], dtype=np.float64)
    Z = np.zeros(nx) + Z_#37.25
    Ar = np.zeros(nx) + Ar_# 157.25
    
    ne = np.zeros(nx) + ne_# 1e27
    ni = ne / Z
    density =  Ar * 1.67E-27* ni 
    #temperatureE = np.linspace(3, 300, nx)
    
    ##Relevant for periodic systems .. fluid cannot be 
    #temperatureE = np.zeros(nx) + 50000 * ev_converter + 50 * ev_converter * np.sin((2*np.pi * centered_coord) / np.max(centered_coord))
    ##Reflective
    temperatureE = np.zeros(nx, dtype=np.float64) + Te_ * ev_converter + 1e-3 * Te_ * ev_converter * np.cos((np.pi * centered_coord) / np.max(centered_coord), dtype=np.float64)
    temperatureI = np.zeros(nx) + 0

    return(initial_coord, density, ne, ni, temperatureE, temperatureI, Z, Ar)

def SmoothRamp(nx, L ):
    lower_cutof = 0.38
    upper_cutof = 0.47
    x_up = np.linspace(0, lower_cutof*L, int(nx*0.3))
    step = x_up[1] - x_up[0]
    x_ramp = np.linspace(lower_cutof*L + step/10, upper_cutof*L, int(nx*0.4) + 1)
    step = x_ramp[1] - x_ramp[0]
    x_down = np.linspace(upper_cutof*L+step/10, L, int(0.3*nx))
    x = np.concatenate((x_up, x_ramp, x_down))
    # x = np.linspace(0, L, nx + 1)
    #x = np.linspace(0, L, nx + 1)
    x_centered = [(x[i] + x[i + 1]) /2 for i in range(len(x) -  1)]
    
    def smooth(x, lower_value, upper_value, lower_limit, upper_limit):
        Var = np.zeros(len(x))
        counter = 0
        for i in range(len(x)):
            if x[i] < lower_limit:
                Var[i] = lower_value
            elif x[i] > upper_limit:
                Var[i] = upper_value
            else:
                # Var[i] = (upper_value - lower_value) * np.tanh((x[i] - lower_limit) / (upper_limit - lower_limit)) + lower_value
                # Var[i] = (upper_value - lower_value) * (1 /(1 + np.exp((x[i] - lower_limit) / (upper_limit - lower_limit))))+ lower_value
            #     # Var[i] = (upper_value - lower_value) * (3 * ((x[i] - lower_limit) / (upper_limit - lower_limit))**2 
                # - 2 * ((x[i] - lower_limit) / (upper_limit - lower_limit))**3) + lower_value
                Var[i] = (upper_value - lower_value) * (6 * ((x[i] - lower_limit) / (upper_limit - lower_limit))**5 
                - 15 * ((x[i] - lower_limit) / (upper_limit - lower_limit))**4 + 10 * ((x[i] - lower_limit) / (upper_limit - lower_limit))**3) + lower_value
                # print(Var[i])
        return Var         

    Te_Up = 2.5E3 * ev_converter
    Te_Down = 0.27E3 * ev_converter
    lower_limit = lower_cutof*L#x_centered[int(nx*0.4)]
    upper_limit = upper_cutof*L#x_centered[int(nx*0.7)]
    Te = smooth(x_centered, Te_Up, Te_Down, lower_limit, upper_limit)

    ne_Up = 4.65 * 1e26
    ne_Down = 80 * 1e26
    ne = smooth(x_centered, ne_Up, ne_Down, lower_limit, upper_limit) 

    Z_Up = 2
    Z_Down = 36.64
    Z = smooth(x_centered, Z_Up, Z_Down, lower_limit, upper_limit) 
    
    Ar_Up = 4
    Ar_Down = 157.25
    Ar = smooth(x_centered, Ar_Up, Ar_Down, lower_limit, upper_limit) 

    Ti = Te
    ni = ne / Z
    density =  Ar* protonMass* ni 
    plt.figure(1)
    plt.plot(x_centered, Te / (ev_converter), 'xr-', label = 'Te')
    plt.legend()
    plt.figure(2)
    plt.plot(x_centered, ne,'xr-', label = 'ne')
    plt.legend()
    plt.figure(3)
    plt.plot(x_centered, Z,'xr-', label = 'Z')
    plt.legend()
    plt.figure(4)
    plt.plot(x_centered, Ar,'xr-', label = 'Ar')
    plt.legend()
    plt.figure(5)
    plt.plot(x_centered, ni,'xr-', label = 'ni')
    plt.figure(6)
    plt.plot(x_centered, density,'xr-', label = 'rho')
    plt.legend()
    plt.show()
    return(x, density, Te, Ti, Z, Ar)

def ExponentialRamp(nx, L):
    nx_Up = 0.4*nx 
    nx_Down = 0.5*nx
    x_up = np.linspace(0, 0.32*L, int(nx*0.1))
    step = x_up[1] - x_up[0]
    x_ramp = np.linspace(0.4*L + step, 0.45*L, int(nx*0.8) + 1)
    step = x_ramp[1] - x_ramp[0]
    x_down = np.linspace(0.45*L+step, L, int(0.1*nx))
    x = np.concatenate((x_up, x_ramp, x_down))
    #x = np.linspace(0, L, nx + 1)
    x_centered = [(x[i] + x[i + 1]) /2 for i in range(len(x) -  1)]
    L = abs(x_centered[int(nx*0.1)] - x_centered[int(nx*0.9) +1])#abs(x[int(nx_Up)] - x[int(nx_Down)])
    T_Up = 2.56E3 * ev_converter
    T_Down = 0.27E3 * ev_converter
    Te_fitting_param =  (L / (x_centered[int(nx*0.9)+1] - x_centered[int(nx*0.1)]))* np.log(T_Up / T_Down) 
    Te_ramp = T_Up * np.exp(-1 * Te_fitting_param * ((x_centered[int(nx*0.1):int(nx*0.9)+1] - x_centered[int(nx*0.1)]) / L )) 
    Te_begin = np.zeros(int(nx*0.1)) +  T_Up   
    Te_end = np.zeros(int(nx*0.1)) + Te_ramp[-1]
    Te = np.concatenate((Te_begin, Te_ramp, Te_end))
    # plt.plot(x_centered, Te)
    # plt.show()
    ne_Up = 4.65 * 1e26
    ne_Down = 80 * 1e26
    ne_fitting_param =  -1*(L / (x_centered[int(nx*0.9)+1] - x_centered[int(nx*0.1)]))* np.log(ne_Up / ne_Down) 
    # ne_fitting_param = -1 * np.log(ne_Up / ne_Down)
    ne_ramp = ne_Up * np.exp(ne_fitting_param * ((x_centered[int(nx * 0.1):int(nx * 0.9) + 1] - x_centered[int(nx * 0.1)]) / L )) 
    ne_begin = np.zeros(int(nx * 0.1)) +  ne_Up   
    ne_end = np.zeros(int(nx*0.1)) + ne_ramp[-1]
    ne = np.concatenate((ne_begin, ne_ramp, ne_end))
    
    Z_Up = 2
    Z_Down = 36.64
    Z_fitting_param =  -1*(L / (x_centered[int(nx*0.9)+1] - x_centered[int(nx*0.1)]))* np.log(Z_Up / Z_Down) 
    # Z_fitting_param = -1 * np.log(Z_Up / Z_Down)
    Z_ramp = Z_Up * np.exp(ne_fitting_param * ((x_centered[int(nx * 0.1):int(nx * 0.9) + 1] - x_centered[int(nx * 0.1)]) / L )) 
    Z_begin = np.zeros(int(nx*0.1)) +  Z_Up   
    Z_end = np.zeros(int(nx*0.1)) + Z_ramp[-1]
    Z = np.concatenate((Z_begin, Z_ramp, Z_end))

    Ar_Up = 4
    Ar_Down = 157.25
    Ar_fitting_param =  -1*(L / (x_centered[int(nx*0.9)+1] - x_centered[int(nx*0.1)]))* np.log(Ar_Up / Ar_Down) 
    # Ar_fitting_param = -1 * np.log(Ar_Up / Ar_Down)
    Ar_ramp = Ar_Up * np.exp(Ar_fitting_param * ((x_centered[int(nx * 0.1):int(nx * 0.9) + 1] - x_centered[int(nx * 0.1)]) / L )) 
    Ar_begin = np.zeros(int(nx*0.1)) +  Ar_Up   
    Ar_end = np.zeros(int(nx*0.1)) + Ar_ramp[-1]
    Ar = np.concatenate((Ar_begin, Ar_ramp, Ar_end))
    
    Ti = Te
    ni = ne / Z
    density =  Ar* protonMass* ni 
    plt.figure(1)
    plt.plot(x_centered, Te / (ev_converter), 'xr-', label = 'Te')
    plt.legend()
    plt.figure(2)
    plt.plot(x_centered, ne,'xr-', label = 'ne')
    plt.legend()
    plt.figure(3)
    plt.plot(x_centered, Z,'xr-', label = 'Z')
    plt.legend()
    plt.figure(4)
    plt.plot(x_centered, Ar,'xr-', label = 'Ar')
    plt.legend()
    plt.figure(5)
    plt.plot(x_centered, ni,'xr-', label = 'ni')
    plt.figure(6)
    plt.plot(x_centered, density,'xr-', label = 'rho')
    plt.legend()
    #plt.show()
    return(x, density, Te, Ti, Z, Ar)



nx = 200
x_l = 0
#x_u = 1e-19
x_u =  10E-3
L = x_u - x_l
Ar = np.zeros(nx) + 157
Z = np.zeros(nx) + 37.5
testName = "hydro_energy_diff"
gammaFactor = 1.4
laserWavelength = 10e-9
LaserPower = 0
coulombLog = 11
#Ev
temperature = 1e-9

ne = 1E19
nc = ne#1.1E15 / pow(laserWavelength, 2)

velocity = np.zeros(nx + 1) #+ add any function
path = "/home/abetharan/HYDRO_KINETIC_COUPLING/Non_Linear_Ramp_Investigation/"
coord, density, temperatureE, temperatureI, Z, Ar = SmoothRamp(nx, L)
# coord, density, numberDensityE, numberDensityI, temperatureE, temperatureI, specificHeatE, DpDTe, specificHeatI, DpDTi, pressureE, pressureI, pressureTotal, IntEe, IntEi = TwoTemperatureBath(nx = 100)

# coord, density, _,_ temperatureE, temperatureI, Z, Ar = epperlein_short(L,63, 64,157,1e19,100)
# coord, density, temperatureE, temperatureI, Z, Ar = ExponentialRamp(nx, 0.4*nx, 0.5*nx, 2.3E-3, 3500*ev_converter, 0.25*3500*ev_converter)
# coord, density, temperatureE, temperatureI, Z, Ar = ExponentialRamp(nx, L)
# coord, density,ne,ni, temperatureE, temperatureI, Z, Ar = epperlein_short(nx, L)

coord = np.linspace(0, 5e-3, nx + 1)
ne = np.zeros(nx) + 1e25
Z = np.zeros(nx) + 30
ni = ne/Z
Ar = np.zeros(nx) + 63.546
density = Ar * protonMass * ni
temperatureE  = np.zeros(nx) + 100
temperatureI = temperatureE
mass = CalculateMass(coord, density, nx = nx)
qe = np.random.rand(nx)
#new_path ='/home/abetharan/HYDRO_KINETIC_COUPLING/Non_Linear_Ramp_Investigation/Smoothed_Pre_Heat_Ramp' 
# new_path = '/Users/shiki/Documents/Imperial_College_London/Ph.D./HyKiCT/init_data/'
# TextDump(   path = new_path,
#             coord= coord ,
#             velocity = velocity,
#             density = density,
#             Te = temperatureE,
            # Ti = temperatureI,
#             mass = mass,
#             Z = Z,
#             Ar = Ar,
# )
