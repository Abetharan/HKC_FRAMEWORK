def set_hydro_init(nx, Ar, Z, Cq, Gamma, CFL, LaserWavelength,
 LaserPower, durOfLaser, steps, tmax, initialDt,
OutputFrequency, InitPath, OutPath):

    kappa = {
        'nx':nx,
        'Ar':Ar,
        'Z':Z,
        'Cq':Cq,
        'Gamma':Gamma,
        'LaserWavelength':LaserWavelength,
        'CFL':CFL,
        'LaserPower':LaserPower,
        'durOfLaser':durOfLaser,
        'steps':steps,
        'tmax':tmax,
        'initialDt':initialDt,
        'OutputFrequency':OutputFrequency,
        'InitPath':InitPath,
        'OutPath':OutPath
    }
    return(kappa)