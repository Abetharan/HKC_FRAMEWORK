import numpy as np 
from . import heat_flow_coupling_tools as hfct 
from scipy import constants
BOLTZMANN_CONSTANT = constants.value("Boltzmann constant")
ELECTRON_MASS = constants.value("electron mass")
PROTON_MASS = constants.value("proton mass")
ELEMENTARY_CHARGE = constants.value("elementary charge")
VACUUM_PERMITTIVITY = 8.854188E-12    # Vacuum dielectric constant
PLANCK_CONSTANT = constants.value("Planck constant")
BOHR_RADIUS = constants.value("Bohr radius")

class TimeStepper():

    def __init__(self, tmax):
        self.guess_time = tmax
        self.hfct_obj = hfct.HeatFlowCouplingTools()
    def divqtimestepper(self, Te, cv, divq):
        tmax = self.guess_time
        n_iter_limit = 100
        i = 0
        while(True): 
            evolved_Te = Te + divq*tmax / cv
            reldiff = abs(evolved_Te - Te) / Te
            if any(reldiff > 0.05) or any(np.isnan(reldiff)):
                tmax *= 0.05 
            else:
                break
            if i > n_iter_limit:
                break;
            i += 1
        return tmax

    def conduct(self, heat_flow, mass): 
        nx = len(heat_flow) 
        HeatConductionE = np.zeros(nx - 1)
        for i, m in enumerate(mass):
            HeatConductionE[i] = (-(heat_flow[i + 1] - heat_flow[i])) / m #- sign is there because of convention used in HyKiCT 
        return HeatConductionE

    def multitimestepper(self, x_wall,x_cent, Te, ne, Z, cv, mass, multi):
        tmax = self.guess_time 
        self.hfct_obj.electron_temperature = Te
        self.hfct_obj.electron_number_density = ne
        self.hfct_obj.zbar = Z
        self.hfct_obj.cell_wall_coord = x_wall
        self.hfct_obj.cell_centered_coord = x_cent
        self.hfct_obj.mass = mass
        n_iter_limit = 100
        i = 0
        drop_dec_val = False 
        drop_val = 0.07
        drop_inc_val = False 
        inc_val = 1.02
        while(True):
            evolved_Te = Te 
            dt = tmax * 0.01
            steps = int(tmax / dt)
            for _ in range(steps):
                self.hfct_obj.electron_temperature = evolved_Te
                self.hfct_obj.lambda_ei(self.hfct_obj.electron_temperature * (BOLTZMANN_CONSTANT/ELEMENTARY_CHARGE), 
                                    self.hfct_obj.electron_number_density,
                                    self.hfct_obj.zbar)
                self.hfct_obj.spitzerHarmHeatFlow()
                heat_flow = self.hfct_obj.spitzer_harm_heat  * multi
                heatconduc = self.conduct(heat_flow, mass)
                evolved_Te = evolved_Te + heatconduc*dt / cv
                if any(np.isnan(evolved_Te)):
                    break

            reldiff = abs(evolved_Te - Te) / Te
            if any(reldiff > 0.1) or any(np.isnan(reldiff)):
                if drop_dec_val: 
                    drop_val *=1.2
                    drop_dec_val = False
                    drop_inc_val = True
                tmax *= drop_val
            elif (np.max(reldiff) < 0.09):
                if drop_inc_val: 
                    inc_val -= inc_val*0.001
                    drop_inc_val = False 
                tmax *=inc_val
                drop_dec_val = True
            else:
                break

            if i > n_iter_limit:
                break
            i += 1
        if tmax > self.guess_time:
            tmax = self.guess_time 
        return tmax