def set_hydro_init(nx, Ar, Z, Cq, Gamma, CFL, laserWavelength,
 laserPower, durOfLaser, laserLoc, steps, tMax,initialDt, dtGlobalMax, dtGlobalMin,
outputFrequency, boundaryCondition,initPath, outPath, switchPath):

    kappa = {
        'nx':nx,
        'Ar':Ar,
        'Z':Z,
        'Cq':Cq,
        'Gamma':Gamma,
        'CFL':CFL,
        'LaserWavelength':laserWavelength,
        'LaserPower':laserPower,
        'durOfLaser':durOfLaser,
        'LaserLocation':laserLoc,
        'steps':steps,
        'tmax':tMax,
        'initialDt':initialDt,
        'dtGlobalMax':dtGlobalMax,
        'dtGlobalMin':dtGlobalMin,
        'OutputFrequency':outputFrequency,
        'BoundaryCondition':boundaryCondition,
        'InitPath':initPath,
        'OutPath':outPath,
        'SwitchPath': switchPath,
    }
    return(kappa)