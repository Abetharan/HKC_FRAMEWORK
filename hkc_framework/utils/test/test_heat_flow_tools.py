import pytest
import numpy as np 
import os 
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/..../')
from utils.heat_flow_coupling_tools import HeatFlowCouplingTools


class TestHFCT():

    def test_spitzer_harm(self):
        hfct_obj = HeatFlowCouplingTools()
        fluid_Te = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_1.txt'))
        fluid_ne = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_1.txt'))
        fluid_Z = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/ZBAR/ZBAR_1.txt'))
        fluid_x_grid = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/CELL_WALL_X/CELL_WALL_X_1.txt'))
        fluid_x_centered_grid = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/CELL_CENTRE_X/CELL_CENTRE_X_1.txt'))
        fluid_mass = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/mass.txt'))
        hfct_obj.electron_temperature = fluid_Te
        hfct_obj.electron_number_density = fluid_ne
        hfct_obj.zbar = fluid_Z
        hfct_obj.cell_wall_coord = fluid_x_grid
        hfct_obj.cell_centered_coord = fluid_x_centered_grid
        hfct_obj.mass = fluid_mass
        hfct_obj.lambda_ei(hfct_obj.electron_temperature, 
                            hfct_obj.electron_temperature,
                            hfct_obj.zbar)
        hfct_obj.spitzerHarmHeatFlow()
        true_spitzer_harm = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/ELECTRON_HEAT_FLOW_X/ELECTRON_HEAT_FLOW_X_1.txt'))

        true_spitzer_harm == pytest.approx(hfct_obj.spitzer_harm_heat)

    def test_div_q(self):
        hfct_obj = HeatFlowCouplingTools()
        fluid_Te = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_1.txt'))
        fluid_ne = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_1.txt'))
        fluid_Z = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/ZBAR/ZBAR_1.txt'))
        fluid_x_grid = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/CELL_WALL_X/CELL_WALL_X_1.txt'))
        fluid_x_centered_grid = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/CELL_CENTRE_X/CELL_CENTRE_X_1.txt'))
        fluid_mass = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/mass.txt'))
        hfct_obj.electron_temperature = fluid_Te
        hfct_obj.electron_number_density = fluid_ne
        hfct_obj.zbar = fluid_Z
        hfct_obj.cell_wall_coord = fluid_x_grid
        hfct_obj.cell_centered_coord = fluid_x_centered_grid
        hfct_obj.mass = fluid_mass
        hfct_obj.lambda_ei(hfct_obj.electron_temperature, 
                            hfct_obj.electron_temperature,
                            hfct_obj.zbar)
        hfct_obj.spitzerHarmHeatFlow()
        hfct_obj.vfp_heat = hfct_obj.spitzer_harm_heat
        calc_div_q_heat = hfct_obj.divQHeatFlow()
        true_div_q = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/ELECTRON_HEAT_CONDUCTION/ELECTRON_HEAT_CONDUCTION_1.txt'))

        true_div_q == pytest.approx(calc_div_q_heat)
    
    def test_multiplier(self):

        hfct_obj = HeatFlowCouplingTools()
        true_spitzer_harm = np.zeros(100)
        heat_flow_val = np.linspace(20,100,20)
        
        hfct_obj.cell_wall_coord = np.linspace(0, 100, 101)
        true_spitzer_harm[40:60] = heat_flow_val
        true_vfp_heat = np.linspace(1, 100, 100)
        hfct_obj.spitzer_harm_heat = true_spitzer_harm
        hfct_obj.vfp_heat = true_vfp_heat
    
        
        expected_multipliers = true_vfp_heat / true_spitzer_harm
        expected_multipliers[np.isnan(expected_multipliers)] = 0
        expected_multipliers[np.isinf(expected_multipliers)] = 0

        (q_vfp_q_sh_multipliers, pre_heat_start_index, 
        pre_heat_last_index, pre_heat_fit_params, 
        front_heat_start_index, front_heat_last_index, front_heat_fit_params) = hfct_obj.multiplier()

        assert expected_multipliers == pytest.approx(q_vfp_q_sh_multipliers)
        assert pre_heat_start_index == pytest.approx(59)
        assert pre_heat_last_index == pytest.approx(100)
        assert front_heat_start_index == pytest.approx(40)
        assert front_heat_last_index == pytest.approx(0)

    def test_detection(self):

        hfct_obj = HeatFlowCouplingTools()
        true_spitzer_harm = np.linspace(1,1, 100)
        true_vfp_heat = np.linspace(1, 1, 100)
        
        hfct_obj.spitzer_harm_heat = true_spitzer_harm
        hfct_obj.vfp_heat = true_vfp_heat
        hfct_obj.cell_wall_coord = np.linspace(0, 100, 101)
        
        expected_multipliers = true_vfp_heat / true_spitzer_harm

        (q_vfp_q_sh_multipliers, pre_heat_start_index, 
        pre_heat_last_index, pre_heat_fit_params, 
        front_heat_start_index, front_heat_last_index, front_heat_fit_params) = hfct_obj.multiplier()

        #test detection 
        print(pre_heat_fit_params)
        print(front_heat_fit_params) 
        assert pre_heat_fit_params is None
        assert front_heat_fit_params is None
