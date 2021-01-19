import pytest
import numpy as np 
import os 
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../../../')
from hkc_framework.couple_methods.divq import DivQ
from hkc_framework.couple_methods.multiplier import Multiplier
from hkc_framework.couple_methods.subtract import Subtract
from hkc_framework.utils.heat_flow_coupling_tools import HeatFlowCouplingTools


class TestCoupling():

    def test_div_q(self):
        div_obj = DivQ()
        fluid_mass = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/mass.txt'))
        sh_heat = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/ELECTRON_HEAT_FLOW_X/ELECTRON_HEAT_FLOW_X_0.txt'))

        div_obj.method(sh_heat, sh_heat, mass = fluid_mass)
        true_div_q = np.loadtxt(os.path.join(myPath, 'test_spitzer_harm/ELECTRON_HEAT_CONDUCTION/ELECTRON_HEAT_CONDUCTION_1.txt'))
        assert fluid_mass == pytest.approx(div_obj.mass)
        assert true_div_q == pytest.approx(div_obj.HeatConductionE)
    
    def test_full_grad_multiplier(self):
        
        multi_obj = Multiplier()
        q_sh = np.linspace(100, 1000, 101)
        q_sh[0] = 0
        q_sh[-1] = 0
        q_vfp = 0.05*np.linspace(100, 1000, 101)
        cell_wall_coord = np.linspace(0, 10, 101)
        heat_flow_val = np.linspace(20,100,20)
        q_snb = None 
        multi_obj.method(q_sh, q_vfp,
                    laser_dir = None, cell_wall_coord = cell_wall_coord,
                    q_snb = q_snb)
        expected_multipliers = np.zeros(len(q_vfp)) + 0.05
        assert q_sh == pytest.approx(multi_obj.spitzer_harm_heat)
        assert q_vfp == pytest.approx(multi_obj.vfp_heat)
        assert cell_wall_coord == pytest.approx(multi_obj.cell_wall_coord)
        assert q_snb == multi_obj.q_snb
        assert expected_multipliers == pytest.approx(multi_obj.q_vfp_q_sh_multipliers)
      
    
    def test_zero_grad_areas_multiplier(self):
        
        multi_obj = Multiplier()
        true_spitzer_harm = np.zeros(100)
        heat_flow_val = np.linspace(20,100,20)
        cell_wall_coord = np.linspace(0, 100, 100)
        true_spitzer_harm[40:60] = heat_flow_val
        true_vfp_heat = np.linspace(1, 100, 100)
        true_vfp_heat[0] = 0 
        true_vfp_heat[-1] = 0 
        multi_obj.method(true_spitzer_harm, true_vfp_heat,
                    laser_dir = None, cell_wall_coord = cell_wall_coord,
                    q_snb = None)
    
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

    def test_subtract(self):

        subtract_obj = Subtract()
        true_spitzer_harm = np.zeros(101)
        heat_flow_val = np.linspace(20,100,20)
        true_spitzer_harm[40:60] = heat_flow_val
        true_vfp_heat = np.linspace(1, 100, 101)
        true_vfp_heat[0] = 0 
        true_vfp_heat[-1] = 0 
        subtract_obj.method(true_spitzer_harm, true_vfp_heat, q_snb = None)
        assert true_vfp_heat == pytest.approx(true_spitzer_harm- subtract_obj.subtract_factor)
        

    @pytest.mark.parametrize("method, laser_dir",[("subtract", "right"),
     ("divq", "right"), ("multi", "right")])
    def test_limit_density(self, method, laser_dir):
        hfct_obj = HeatFlowCouplingTools() 
        if method == "subtract":
            couple_obj = Subtract()
        elif method == "divq":
            couple_obj = DivQ()
        elif method == "multi":
            couple_obj = Multiplier()

        fluid_Te = np.loadtxt(os.path.join(myPath, 'test_limit_density/electron_temperature.txt'))
        fluid_ne = np.loadtxt(os.path.join(myPath, 'test_limit_density//')) 
        fluid_Z = np.loadtxt(os.path.join(myPath, 'test_limit_density/zbar.txt'))
        fluid_x_grid = np.loadtxt(os.path.join(myPath, 'test_limit_density/cell_wall_x.txt'))
        fluid_x_centered_grid = np.loadtxt(os.path.join(myPath, 'test_limit_density/cell_centre_x.txt')) 
        fluid_mass = np.loadtxt(os.path.join(myPath, 'test_limit_density/mass.txt'))
        vfp_heat = np.loadtxt(os.path.join(myPath, 'test_limit_density/vfp_heat_flow.txt'))
        sh_heat = np.loadtxt(os.path.join(myPath, 'test_limit_density/sh_heat_flow.txt'))
        vfp_heat_conduction = np.loadtxt(os.path.join(myPath, 'test_limit_density/thermal_conduction.txt'))

        hfct_obj.electron_temperature = fluid_Te
        hfct_obj.electron_number_density = fluid_ne
        hfct_obj.zbar = fluid_Z
        hfct_obj.cell_wall_coord = fluid_x_grid
        hfct_obj.cell_centered_coord = fluid_x_centered_grid
        hfct_obj.mass = fluid_mass
        hfct_obj.lambda_ei(hfct_obj.electron_temperature * (11594), 
                            hfct_obj.electron_number_density,
                            hfct_obj.zbar)
        hfct_obj.spitzerHarmHeatFlow()

        assert sh_heat == pytest.approx(hfct_obj.spitzer_harm_heat)

        couple_obj.method(hfct_obj.spitzer_harm_heat, vfp_heat, 
                        laser_dir = laser_dir, mass = fluid_mass, cell_wall_coord = fluid_x_grid,
                        q_snb = hfct_obj.q_snb)

        if(method == 'divq'):
            assert vfp_heat_conduction == pytest.approx(couple_obj.HeatConductionE)
        if(method == 'multiplier'):
            assert vfp_heat_conduction == pytest.approx(couple_obj.q_vfp_q_sh_multipliers * hfct_obj.spitzer_harm_heat)
        if(method == 'subtract'):
            assert vfp_heat_conduction == pytest.approx(hfct_obj.spitzer_harm_heat- couple_obj.subtract_factor)