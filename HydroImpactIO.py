from scipy.interpolate import CubicSpline
import numpy as np 
import os
import fnmatch
from string import Template
from scipy import constants
import impact_norms_py3 as ImNorms

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

def fileWriteFormat(tmpFileLocation, OpenFile, coord, varData):   
    
    if varData != coord:
        nx = len(coord) - 1 #becuse index start from 0
        cellCenteredCoord = [(coord[i + 1] + coord[i]) / 2 for i in range(nx)]
        LastDx = coord[nx - 1] - coord[nx - 2]
        LastGhostCellCoord = cellCenteredCoord[nx - 2] + 0.5 * LastDx 
        cellCenteredCoordWithGhost = [-cellCenteredCoord[0]] + cellCenteredCoord + [LastGhostCellCoord]
        leadingDim = 3
        varDataWithGhost = [varData[0]] + varData + [varData[nx - 1]]
        TemplateSwapDict = {'leadingDim':leadingDim, 'xlist':''.join(cellCenteredCoordWithGhost),
         'arraylist':'\n'.join(varDataWithGhost), 'arraylist1':'\n'.join(varDataWithGhost), 'arraylist2':'\n'.join(varDataWithGhost)}

    else:
        leadingDim = 1
        TemplateSwapDict = {'leadingDim':leadingDim, 'xlist':''.join(coord), 'arraylist':'\n'.join(varData),
            'arraylist1':'', 'arraylist2':''}
    
    tmpFileIn = open(tmpFileLocation)
    src = Template(tmpFileIn.read())
    result = src.substitute(TemplateSwapDict)
    OpenFile.write(result)
        
    return(0)


def HydroToImpact(path, Z, Ar, fluidNx):

        
    lastIndex = findLastIndex(path, "Coord")
    fluid_x  = np.loadtxt(path + "/Coord_" + str(lastIndex) + ".txt")
    fluid_ne = np.loadtxt(path + "/NumberDensityE_" + str(lastIndex) + ".txt")
    fluid_ni = np.loadtxt(path + "/NumberDensityI_" + str(lastIndex) + ".txt")
    fluid_Te = np.loadtxt(path + "/TemperatureE_" + str(lastIndex) + ".txt")
    fluid_las_dep = np.loadtxt(path + "/InverseBrem_" + str(lastIndex) + ".txt")
    fluid_brem = np.loadtxt(path + "/Brem_" + str(lastIndex) + ".txt")

    avgNe = np.average(fluid_ne)
    avgTe = np.average(fluid_Te)
    normalised_values = ImNorms.impact_inputs(avgNe, avgTe, Z, Ar, Bz = 0)

    ## NOrmalise SI to Impact norms 
    x_norm = fluid_x / normalised_values["lambda_mfp"]
    ne_norm = fluid_ne / (1e6* normalised_values['ne'])  # 1e6 there to convert from m^-3 to cm^-3
    ni_norm = fluid_ni / (1e6 * normalised_values['ni']) # 1e6 there to convert from m^-3 to cm^-3
    Te_norm = fluid_Te / ((e/kb)* normalised_values['Te'])


    kinetic_x = np.linspace(x_norm[0], x_norm[-1], os.environ["NXM"])
               # np.geomspace(fluid_x[0, fluid_x[-1], nx)
               # np.logspace(fluid_x[0, fluid_x[-1], nx)
    centered_x = [(x_norm[i + 1] - x_norm[i])/2 for i in range(fluidNx)]
    cs_ne = CubicSpline(centered_x, ne_norm)
    cs_ni = CubicSpline(centered_x, ni_norm)
    cs_Te = CubicSpline(centered_x, Te_norm)

    kinetic_ne = cs_ne(kinetic_x)
    kinetic_ni = cs_ni(kinetic_x)
    kinetic_Te = cs_Te(kinetic_x)

    
    if not os.path.exists(path + "/tmpWrite.txt"):

        writeStatement = """
        Version:2
        0
        2
        $leadingDim     
        $xlist
        \n
        $arraylist $arraylist1 $arraylist2
        """
        kappa = open(path + "/tmpWrite.txt")
        kappa.write(writeStatement)

    impactNeFile = open(path + "/" + os.environ["RUN"] + "_eden.xy", "w")
    impactNiFile = open(path + "/" + os.environ["RUN"] + "_ionden.xy", "w")
    impactTeFile = open(path + "/" + os.environ["RUN"] + "_tmat.xy", "w")
    impactXFile = open(path + "/" + os.environ["RUN"] + "_xf.xy", "w")
    impactLaserFile = open(open(path + "/" + os.environ["RUN"] + "_laserdep.xy", "w"))
    impactRadFile = open(open(path + "/" + os.environ["RUN"] + "_rad_to_electron.xy", "w"))
    
    fileWriteFormat(path + "/tmpWrite.txt", impactNeFile, kinetic_x, kinetic_ne)
    fileWriteFormat(path + "/tmpWrite.txt", impactNiFile, kinetic_x, kinetic_ni)
    fileWriteFormat(path + "/tmpWrite.txt", impactXFile, kinetic_x, kinetic_x)
    fileWriteFormat(path + "/tmpWrite.txt", impactTeFile, kinetic_x, kinetic_Te)
    fileWriteFormat(path + "/tmpWrite.txt", impactLaserFile, kinetic_x, kinetic_Te)
    fileWriteFormat(path + "/tmpWrite.txt", impactRadFile, kinetic_x, kinetic_Te)

    return(np.mean(fluid_ne), np.mean(fluid_Te))

def ImpactToHydro(path):

        ##Convert to SI from Impact Norms
    varList = ["n", "ni", "Te", "qxX"]
    timeStep = findLastIndex(path, 'qxX')

    for var in varList:
        varArrays = cf.load_dict(path,os.environ['RUN'],var, timeStep, iter_number = None)
        varList = varArrays['mat'][:, 1]

        if var == "n":
            outputVar = "ne"
        elif var == "Te":
            outputVar = "electron_temperature"
        else:
            outputVar =var
        fileToWrite = open(path + "/" +outputVar + ".txt" )
