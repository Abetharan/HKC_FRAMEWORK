def set_hydro_init(nx, Cq, Gamma, CFL, laserWavelength,
 laserPower, durOfLaser, laserLoc, steps, tMax,initialDt, dtGlobalMax, dtGlobalMin,
outputFrequency, boundaryCondition,initPath, outPath, switchPath,FeosPathMaterial1,FeosPathMaterial2):

    kappa = {
        'nx':nx,
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
        'FEOSPathMaterial1':FeosPathMaterial1,
        'FEOSPathMaterial2':FeosPathMaterial2
        
    }
    return(kappa)