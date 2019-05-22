from scipy.interpolate import CubicSpline
import numpy as np 
import os
import fnmatch
from string import Template
from scipy import constants
import impact_norms_py3 as ImNorms
import impact_module_py3 as cf
import HydroRemainingFromImpact as remain

kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")

def findLastIndex(path, var):
    LargestIndex = 0
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, "[a-z]"+ var):
            k = os.path.splitext(file)[0]
            k = k.split("_")
            timeIndex = k[-1]
            if timeIndex > LargestIndex:
                LargestIndex = timeIndex
            
    return(LargestIndex)

def fileWriteFormat(tmpFileLocation, OpenFile, coord, varData, var):   
    
    
    if var != "coord":
        nx = len(coord) - 1 #becuse index start from 0
        cellCenteredCoord = [(coord[i + 1] + coord[i]) / 2 for i in range(nx)]
        LastDx = cellCenteredCoord[nx - 1] - cellCenteredCoord[nx - 2]
        LastGhostCellCoord = cellCenteredCoord[nx - 1] + LastDx 
        cellCenteredCoordWithGhost = [-cellCenteredCoord[0]] + cellCenteredCoord + [LastGhostCellCoord]
        leadingDim = 3
        varDataWithGhost = [varData[0]] + varData.tolist() + [varData[-1]]
        maxParam = "-6.000000000000000e+00e3 0.000000000000000e+00 6.000000000000000e+00e3"
        length = len(varDataWithGhost)
        TemplateSwapDict = {'leadingDim':leadingDim, 'maxparam':maxParam, 'Len':length ,'xlist':' '.join(str(coo) for coo in cellCenteredCoordWithGhost),
                            'arraylist':'\n'.join((str(var) + "\t" + str(var) + "\t" + str(var)) for var in varDataWithGhost)}

    else:
        leadingDim = 1
        maxParam = "0.000000000000000e+00"
        length = len(coord)
        TemplateSwapDict = {'leadingDim':leadingDim, 'maxparam':maxParam, 'Len':length, 'xlist':' '.join(str(coo) for coo in coord), 
                            'arraylist':'\n'.join(str(var) for var in varData)}
    
    tmpFileIn = open(tmpFileLocation)
    src = Template(tmpFileIn.read())
    result = src.substitute(TemplateSwapDict)
    OpenFile.write(result)
    OpenFile.close()
        
    return(0)


def HydroToImpact(fluidPath, cyclePath, Z, Ar, laserWaveLength, fluidNx):

        
    lastIndex = findLastIndex(fluidPath, "Coord")
    fluid_x  = np.loadtxt(fluidPath + "/Coord_" + str(lastIndex) + ".txt")
    fluid_ne = np.loadtxt(fluidPath + "/NumberDensityE_" + str(lastIndex) + ".txt")
    fluid_ni = np.loadtxt(fluidPath + "/NumberDensityI_" + str(lastIndex) + ".txt")
    fluid_Te = np.loadtxt(fluidPath + "/TemperatureE_" + str(lastIndex) + ".txt")
    fluid_las_dep = np.loadtxt(fluidPath + "/InverseBrem_" + str(lastIndex) + ".txt")
    fluid_brem = np.loadtxt(fluidPath + "/Brem_" + str(lastIndex) + ".txt")
    fluid_density = np.loadtxt(fluidPath + "/Density_" + str(lastIndex) + ".txt")

    nc = 1.1E9 / pow(laserWaveLength,2) 
    avgTe = 1000
    normalised_values = ImNorms.impact_inputs(nc, avgTe, Z, Ar, Bz = 0)

    ## NOrmalise SI to Impact norms 
    x_norm = fluid_x / normalised_values["lambda_mfp"]
    ne_norm = fluid_ne / (1e6 * nc)  # 1e6 there to convert from m^-3 to cm^-3
    ni_norm = fluid_ni / (1e6 *  (nc / Z)) # 1e6 there to convert from m^-3 to cm^-3
    Te_norm = (fluid_Te / ((e/kb) * 2 * avgTe)) 
    laser_norm = (fluid_las_dep * fluid_density) / (normalised_values["vte"] / (me * normalised_values["tau_ei"]))
    brem_norm = (fluid_brem * fluid_density) / (normalised_values["vte"] / (me * normalised_values["tau_ei"]))

    kinetic_x = np.linspace(x_norm[0], x_norm[-1], int(os.environ["NXM"]) + 1)
               # np.geomspace(fluid_x[0, fluid_x[-1], nx)
               # np.logspace(fluid_x[0, fluid_x[-1], nx)
    fluid_centered_x = [(x_norm[i + 1] + x_norm[i])/2 for i in range(fluidNx)]
    kinetic_centered_x = [(kinetic_x[i + 1] + kinetic_x[i])/2 for i in range(int(os.environ["NXM"]))]

    cs_ne = CubicSpline(fluid_centered_x, ne_norm)
    cs_ni = CubicSpline(fluid_centered_x, ni_norm)
    cs_Te = CubicSpline(fluid_centered_x, Te_norm)
    cs_laser = CubicSpline(fluid_centered_x, laser_norm)
    cs_brem = CubicSpline(fluid_centered_x, brem_norm)

    kinetic_ne = cs_ne(kinetic_centered_x)
    kinetic_ni = cs_ni(kinetic_centered_x)
    kinetic_Te = cs_Te(kinetic_centered_x)
    kinetic_laser = cs_laser(kinetic_centered_x)
    kinetic_brem = cs_brem(kinetic_centered_x)

    
    if not os.path.exists(cyclePath + "/tmpWrite.txt"):

        writeStatement = """Version:2
0
2
$leadingDim
$maxparam
$Len     
$xlist\n
$arraylist
        """
        kappa = open(cyclePath + "/tmpWrite.txt", "w")
        kappa.write(writeStatement)
        kappa.close()

    impactNeFile = open(cyclePath + "/" + os.environ["RUN"] + "_eden.xy", "w")
    impactNiFile = open(cyclePath + "/" + os.environ["RUN"] + "_ionden.xy", "w")
    impactTeFile = open(cyclePath + "/" + os.environ["RUN"] + "_tmat.xy", "w")
    impactXFile = open(cyclePath + "/" + os.environ["RUN"] + "_xf.xy", "w")
    impactLaserFile = open(cyclePath + "/" + os.environ["RUN"] + "_laserdep.xy", "w")
    impactRadFile = open(cyclePath + "/" + os.environ["RUN"] + "_rad_to_electron.xy", "w")
    
    fileWriteFormat(cyclePath + "/tmpWrite.txt", impactNeFile, kinetic_x, kinetic_ne, "ne")
    fileWriteFormat(cyclePath + "/tmpWrite.txt", impactNiFile, kinetic_x, kinetic_ni, "ni")
    fileWriteFormat(cyclePath + "/tmpWrite.txt", impactXFile, kinetic_x, kinetic_x, "coord")
    fileWriteFormat(cyclePath + "/tmpWrite.txt", impactTeFile, kinetic_x, kinetic_Te, "Te")
    fileWriteFormat(cyclePath + "/tmpWrite.txt", impactLaserFile, kinetic_x, kinetic_laser, "laser")
    fileWriteFormat(cyclePath + "/tmpWrite.txt", impactRadFile, kinetic_x, kinetic_brem, "Rad")

    return(nc, avgTe, normalised_values)

def ImpactToHydro(cycleDumpPath, normalisedValues, Z, massNumber, cycle, gammaFactor, laserWaveLength, laserPower, fluidNx):

    ##Convert to SI from Impact Norms
    varList = ["n", "Te", "qxX"]
    kineticDumpPath = cycleDumpPath + "/kinetic_out/"
    for var in varList:            
        timeStep = findLastIndex(kineticDumpPath, var)
        varArrays = cf.load_dict(kineticDumpPath,os.environ['RUN'],var, str(timeStep), iter_number = None)
        varList = varArrays['mat'][:, 1]

        if var == "n":
            #outputVar = "ne"
            normConst = normalisedValues['ne'] * 1e6 #to m**-3
            ne = varList * normConst
        elif var == "Te":
            #outputVar = "electron_temperature"
            normConst = normalisedValues['Te'] * 2 * 11600 #To kelvin
            Te = varList * normConst
        elif var == "qxX":
            #outputVar = "electron_heat_flow"
            normConst = 9.11E-31 * normalisedValues['vte']**3*normalisedValues['ne'] * 1e6 
            qxX = varList * normConst

    remain.CalculateRemain(ne, Te, qxX, Z, massNumber, cycle, gammaFactor, laserWaveLength, laserPower, fluidNx, cycleDumpPath)




















 