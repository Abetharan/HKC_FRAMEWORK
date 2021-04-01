import pytest
import numpy as np 
import os 
import sys
from scipy import constants
mp = constants.value("proton mass")
BOLTZMANN_CONSTANT = constants.value("Boltzmann constant")
ELECTRON_MASS = constants.value("electron mass")
PROTON_MASS = constants.value("proton mass")
ELEMENTARY_CHARGE = constants.value("elementary charge")
VACUUM_PERMITTIVITY = 8.854188E-12    # Vacuum dielectric constant
PLANCK_CONSTANT = constants.value("Planck constant")
BOHR_RADIUS = constants.value("Bohr radius")
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../../../')
from hkc_framework.utils.time_stepper import TimeStepper

class TestTimeStepper():

    def test_multiplier_stepper(self):
        myPath2 = os.path.join(myPath, "test_timestep_files")
        fluid_Te = np.loadtxt(os.path.join(myPath2, "ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_0.txt"))
        fluid_specific_heat = np.loadtxt(os.path.join(myPath2, "ELECTRON_SPECIFIC_HEAT/ELECTRON_SPECIFIC_HEAT_0.txt"))
        fluid_x_centered_grid = np.loadtxt(os.path.join(myPath2, "CELL_CENTRE_X/CELL_CENTRE_X_0.txt"))
        fluid_x_grid = np.loadtxt(os.path.join(myPath2, "CELL_WALL_X/CELL_WALL_X_0.txt"))
        fluid_ne = np.loadtxt(os.path.join(myPath2, "ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_0.txt")) 
        fluid_Z = np.loadtxt(os.path.join(myPath2, "Zbar/Zbar_0.txt")) 
        fluid_mass = np.loadtxt(os.path.join(myPath2, "mass.txt"))
        multi = np.loadtxt(os.path.join(myPath2, "qe.txt")) 
        guess = 1.0e-11
        time_stepper = TimeStepper(guess)
        # tmax = multitimestepper(x, Te, ne, Z, cv, mass, multi, guess)
        tmax = time_stepper.multitimestepper(fluid_x_grid,fluid_x_centered_grid, fluid_Te, fluid_ne, fluid_Z, fluid_specific_heat, fluid_mass, multi)
        assert tmax - 7.04939587297349e-14 < 1e-14

    def test_divq_stepper(self):
        myPath2 = os.path.join(myPath, "test_timestep_files")
        divq = np.loadtxt(os.path.join(myPath2, "ELECTRON_HEAT_CONDUCTION/ELECTRON_HEAT_CONDUCTION_0.txt"))
        fluid_Te = np.loadtxt(os.path.join(myPath2, "ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_0.txt"))
        fluid_specific_heat = np.loadtxt(os.path.join(myPath2, "ELECTRON_SPECIFIC_HEAT/ELECTRON_SPECIFIC_HEAT_0.txt"))
        guess = 1.0e-11
        time_stepper = TimeStepper(guess)
        tmax = time_stepper.divqtimestepper(fluid_Te, fluid_specific_heat, divq)
        assert tmax - 5.0e-13 < 1e-14
