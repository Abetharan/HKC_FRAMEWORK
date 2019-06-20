
def HydroToImpact(path):
     
    fluid_x  = np.loadtxt(path + "/Coord_" + i + ".txt")
    fluid_ne = np.loadtxt(path + "/NumberDensityE_" + i + ".txt")
    fluid_ni = np.loadtxt(path + "/NumberDensityI_" + i + ".txt")
    fluid_Te = np.loadtxt(path + "/TemperaturE_" + i + ".txt")
    fluid_las_dep = np.loadtxt(path + "/InvBrem_" + i + ".txt")
    fluid_brem = np.loadtxt(path + "/Brem_" + i + ".txt")

    ## NOrmalise SI to Impact norms 
    x_norm = fluid_x / normalised_values["lambda_mfp"]
    ne_norm = fluid_ne / (1e6* normalised_values['ne']) 
    ni_norm = fluid_ne / (1e6 * normalised_values['ni'])
    Te_norm = fluid_Te / ((e/kb)* normalised_values['Te'])


    kinetic_x = np.linspace(x_norm[0], x_norm[-1], nx)
               # np.geomspace(fluid_x[0, fluid_x[-1], nx)
               # np.logspace(fluid_x[0, fluid_x[-1], nx)

    cs_ne = CubicSpline(x_norm, ne_norm)
    cs_ni = CubicSpline(x_norm, ni_norm)
    cs_Te = CubicSpline(x_norm, Te_norm)

    kinetic_ne = cs_ne(kinetic_x)
    kinetic_ni = cs_ni(kinetic_x)
    kinetic_Te = cs_Te(kinetic_x)

    impactNeFile = open(path + "/" + runName + "_eden.xy", "w")
    impactNiFile = open(path + "/" + runName + "_ionden.xy", "w")
    impactTeFile = open(path + "/" + runName + "_tmat.xy", "w")
    impactXFile = open(path + "/" + runName + "_xf.xy", "w")

    impactNeFile.write(kinetic_ne)
    impactNiFile.write(kinetic_ni)
    impactTeFile.write(kinetic_Te)
    impactXFile.write(kinetic_x)

    return(np.mean(fluid_ne), np.mean(fluid_Te))

def ImpactToHydro(path):
