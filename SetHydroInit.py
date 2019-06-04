def set_hydro_init(nx, Ar, Z, Cq, Gamma, CFL, LaserWavelength,
 LaserPower, durOfLaser, steps, tmax,initialDt, dtGlobalMax, dtGlobalMin,
OutputFrequency, BoundaryCondition, InitPath, OutPath):

    kappa = {
        'nx':nx,
        'Ar':Ar,
        'Z':Z,
        'Cq':Cq,
        'Gamma':Gamma,
        'CFL':CFL,
        'LaserWavelength':LaserWavelength,
        'LaserPower':LaserPower,
        'durOfLaser':durOfLaser,
        'steps':steps,
        'tmax':tmax,
        'initialDt':initialDt,
        'dtGlobalMax':dtGlobalMax,
        'dtGlobalMin':dtGlobalMin,
        'OutputFrequency':OutputFrequency,
        'BoundaryCondition':BoundaryCondition,
        'InitPath':InitPath,
        'OutPath':OutPath
    }
    return(kappa)