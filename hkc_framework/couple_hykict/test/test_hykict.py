import numpy as np 
import pytest
import os
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/..../')
from couple_hykict.hykict import HyKiCT

class TestHyKiCT():

    def test_run(self, tmpdir):
        src_dir = os.environ["F_SRC_DIR"]
        p = tmpdir.mkdir('cycle')
        f_obj = HyKiCT(f_src_dir = src_dir, run_path=p, f_config_yml_file_path= myPath + "/test_run_dir/config.yml")
        f_obj._cycle_dump_path = p
        with pytest.raises(SystemExit):
            f_obj.Run()

    def test_last_step(self, tmpdir):
        src_dir = os.environ["F_SRC_DIR"]
        p = tmpdir.mkdir('cycle')
        f_obj = HyKiCT(f_src_dir = src_dir, run_path=p, f_config_yml_file_path= myPath + "/test_run_dir/config.yml")
        f_obj._cycle_dump_path = p
        output_path = myPath + "/test_last_step/OUTPUT/"
        input_path = myPath + "/test_last_step/INPUT/"
        f_obj._fluid_output_path = output_path
        f_obj._init_file_path = input_path
         
        (f_x_grid, f_x_centered_grid, f_v, 
        f_ne, f_Te, f_Z, f_laser, mass) = f_obj.getLastStepQuants()

        f_x_grid == pytest.approx(np.loadtxt(os.path.join(output_path, 'CELL_WALL_X/CELL_WALL_X_1.txt')))
        f_x_centered_grid == pytest.approx(np.loadtxt(os.path.join(output_path, 'CELL_CENTRE_X/CELL_CENTRE_X_1.txt')))
        f_v == pytest.approx(np.loadtxt(os.path.join(output_path, 'VELOCITY/VELOCITY_1.txt')))
        f_ne == pytest.approx(np.loadtxt(os.path.join(output_path, 'ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_1.txt')))
        f_Te == pytest.approx(np.loadtxt(os.path.join(output_path, 'ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_1.txt')))
        f_Z == pytest.approx(np.loadtxt(os.path.join(output_path, 'ZBAR/ZBAR_1.txt')))
        mass == pytest.approx(np.loadtxt(os.path.join(input_path, 'mass.txt')))

    
    def test_nextInit(self, tmpdir):
        src_dir = os.environ["F_SRC_DIR"]
        p = tmpdir.mkdir('cycle')
        f_obj = HyKiCT(f_src_dir = src_dir, run_path=p, f_config_yml_file_path= myPath + "/test_run_dir/config.yml")
        f_obj._cycle_dump_path = p

        qe = np.random.rand(40)

        output_path = myPath + "/test_last_step/OUTPUT/"
        input_path = myPath + "/test_last_step/INPUT/"
        f_obj._fluid_output_path = output_path
        f_obj._init_file_path = input_path

        f_obj.initHydroFromKinetic(tmpdir, qe, qe, qe)
        true_x = np.loadtxt(os.path.join(output_path, 'CELL_WALL_X/CELL_WALL_X_1.txt'))
        true_v = np.loadtxt(os.path.join(output_path, 'VELOCITY/VELOCITY_1.txt'))
        true_density = np.loadtxt(os.path.join(output_path, 'DENSITY/DENSITY_1.txt'))
        true_Te = np.loadtxt(os.path.join(output_path, 'ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_1.txt'))
        true_Z = np.loadtxt(os.path.join(output_path, 'ZBAR/ZBAR_1.txt'))
        true_mass = np.loadtxt(os.path.join(input_path, 'mass.txt'))
        true_Ar = np.loadtxt(os.path.join(input_path, 'Ar.txt'))
        
        true_x == pytest.approx(np.loadtxt(os.path.join(tmpdir, 'coord.txt')))
        true_v == pytest.approx(np.loadtxt(os.path.join(tmpdir, 'velocity.txt')))
        true_density == pytest.approx(np.loadtxt(os.path.join(tmpdir, 'density.txt')))
        true_Te == pytest.approx(np.loadtxt(os.path.join(tmpdir, 'electron_temperature.txt')))
        true_Z == pytest.approx(np.loadtxt(os.path.join(tmpdir, 'Z.txt')))
        true_mass == pytest.approx(np.loadtxt(os.path.join(tmpdir, 'mass.txt')))
        true_Ar == pytest.approx(np.loadtxt(os.path.join(tmpdir, 'Ar.txt')))
        
        qe == pytest.approx(np.loadtxt(os.path.join(tmpdir, 'qe.txt')))
        qe == pytest.approx(np.loadtxt(os.path.join(tmpdir, 'pre_heat_fit_param.txt')))
        qe == pytest.approx(np.loadtxt(os.path.join(tmpdir, 'front_heat_fit_param.txt')))
