import numpy as np
import matplotlib.pyplot as plt
Z = 8
CoulombLog = 11
T0 = 30 * 11600
k = 0.5
K0 =1.83E-10 * pow(T0, 2.5) * pow(CoulombLog, -1) *  pow(Z, -1);
nx = 100
nt = 20
cycles = 28
path = "/media/abetharan/DATADRIVE1/Abetharan/couple3/"
time = np.zeros(nt*cycles + 5)
temperatureE = np.zeros((nt*cycles + 5, nx))
analyticalTe = np.zeros((nt*cycles + 5, nx))
q = 0
cycles = 29
for j in range(cycles):
    if j == 0:
        step = 1
        nt = 5
    else:
        step = 5000
        nt = 100000

    for i in range(0, nt,step):
        if i < nt and j > 0 :
            i+=1
        
        _pathTemp = path + "/cycle_" + str(j) + "/fluid_output/TemperatureE_" + str(i) + ".txt"
        _pathCv = path + "/cycle_" + str(j) + "/fluid_output/SpecificHeatE_" + str(i) + ".txt"
        _pathTime = path + "/cycle_" + str(j) + "/fluid_output/Time_" + str(i) + ".txt"
        _pathx = path + "/cycle_" + str(j) + "/fluid_output/Coord_" + str(i) + ".txt"
        _pathrho = path + "/cycle_" + str(j) + "/fluid_output/Density_" + str(i) + ".txt"
        x = np.loadtxt(_pathx) 
        centered_x = [x[k + 1] - x[k] for k in range(nx)]
        specificHeatE = np.loadtxt(_pathCv)
        
        if q == 0:
            time[q] = np.loadtxt(_pathTime)
        else:
            time[q] = time[q - 1] + np.loadtxt(_pathTime)
        print(time[q])
        temperatureE[q, :] = np.loadtxt(_pathTemp)
        Density = np.loadtxt(_pathrho)
        q+=1

gamma = (K0 * pow(k,2) / specificHeatE * Density)
centered_x = np.array(centered_x)
centered_x = np.linspace(0, 1000, 100)
l_grid = centered_x[-1]- centered_x[0]
for i in range(len(time)):
    analyticalTe[i, :] = 0.5 * T0 + 0.0005 * T0 * np.exp(- gamma * time[i])*  np.cos(2*np.pi*k*centered_x*(1.0/l_grid))


for j in range(0, len(time), 1):
    print(j)
    #plt.plot(temperatureE[j,:] , "k")
    plt.plot(analyticalTe[j,:] - analyticalTe[j - 1, :], "r--")
    plt.pause(0.0000000000000002)
    plt.gcf().clear()
plt.show()




