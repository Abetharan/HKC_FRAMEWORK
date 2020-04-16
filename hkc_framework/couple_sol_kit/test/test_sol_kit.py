import numpy as np 
import pytest
import os
from scipy import constants
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')
from couple_sol_kit.sol_kit import SOL_KIT

ELEMENTARY_CHARGE = constants.value("elementary charge")
BOLTZMANN_CONSTANT = constants.value("Boltzmann constant")
class TestSOLKiT():

    def test_run(self, tmpdir):
        src_dir = os.environ["K_SRC_DIR"]
        p = tmpdir.mkdir('cycle')
        k_obj = SOL_KIT(p, p, p, src_dir, p, k_config_yml_file_path= myPath + '/test_run_dir/config.yml')
        with pytest.raises(SystemExit):
            k_obj.Run()

    def test_set(self, tmpdir): 
        src_dir = os.environ["K_SRC_DIR"]
        p = tmpdir.mkdir('cycle')
        k_obj = SOL_KIT(p, p, p, src_dir, p, k_config_yml_file_path= myPath + '/test_run_dir/config.yml')
        k_obj.setFiles()
        assert os.path.exists(os.path.join(p, 'INPUT/GRID_INPUT.txt'))
        assert os.path.exists(os.path.join(p, 'INPUT/NORMALIZATION_INPUT.txt'))
        assert os.path.exists(os.path.join(p, 'INPUT/SWITCHES_INPUT.txt'))
        assert os.path.exists(os.path.join(p, 'INPUT/SOLVER_PARAMS_INPUT.txt'))
        assert os.path.exists(os.path.join(p, 'OUTPUT/'))
        assert len(os.listdir(os.path.join(p, 'OUTPUT'))) == 28

    def test_init(self, tmpdir):
        src_dir = os.environ["K_SRC_DIR"]
        p = tmpdir.mkdir('cycle')
        k_obj = SOL_KIT(p, p, p, src_dir, p, k_config_yml_file_path= myPath + '/test_run_dir/config.yml')
        k_obj.setFiles()
        f_x_grid = np.linspace(0, 100, 101)
        f_x_centered_grid = np.array([(f_x_grid[i + 1] + f_x_grid[i]) / 2 for i in range(len(f_x_grid) - 1)])
        f_te = np.linspace(0, 100, 100) 
        f_ne = np.linspace(0, 100, 100)
        f_z = np.linspace(0, 100, 100)
        f_laser = 0
        f_rad = 0
        k_obj.InitFromHydro(f_x_grid, f_x_centered_grid, f_te, f_ne, f_z, f_laser, f_rad)
        init_density = np.loadtxt(os.path.join(p, 'INPUT/DENS_INPUT.txt')) * k_obj.normalised_values['ne']
        init_Te = np.loadtxt(os.path.join(p, 'INPUT/TEMPERATURE_INPUT.txt'))* k_obj.normalised_values['Te'] * (ELEMENTARY_CHARGE/BOLTZMANN_CONSTANT)
        init_Z = np.loadtxt(os.path.join(p, 'INPUT/Z_PROFILE_INPUT.txt'))
        init_x_grid = np.loadtxt(os.path.join(p, 'INPUT/X_GRID_INPUT.txt'))* k_obj.normalised_values['lambda_mfp']
        assert f_x_grid[1:-1] == pytest.approx(init_x_grid[1::2])
        assert f_te == pytest.approx(init_Te[::2])
        assert f_ne == pytest.approx(init_density[::2])
        assert f_z == pytest.approx(init_Z[::2])

    def test_move(self, tmpdir):
        src_dir = os.environ["K_SRC_DIR"]
        run_path = tmpdir.mkdir('cycle')
        cycle_dump_path = tmpdir.mkdir('kinetic/')
        kinetic_output = cycle_dump_path.mkdir("OUTPUT")
        kinetic_input = cycle_dump_path.mkdir("INPUT")
        k_obj = SOL_KIT(run_path, str(kinetic_input), str(kinetic_output), src_dir,
                         cycle_dump_path, k_config_yml_file_path= myPath + '/test_run_dir/config.yml')
        k_obj.setFiles()
        k_obj.moveFiles()
        assert os.path.exists(os.path.join(cycle_dump_path, 'INPUT'))
        assert os.path.exists(os.path.join(cycle_dump_path, 'INPUT/NORMALIZATION_INPUT.txt'))
        assert os.path.exists(os.path.join(cycle_dump_path, 'INPUT/SWITCHES_INPUT.txt'))
        assert os.path.exists(os.path.join(cycle_dump_path, 'INPUT/SOLVER_PARAMS_INPUT.txt'))
        assert os.path.exists(os.path.join(cycle_dump_path, 'OUTPUT/'))
        assert len(os.listdir(os.path.join(cycle_dump_path, 'OUTPUT'))) == 28

    def test_getLastHeat(self, tmpdir):
        src_dir = os.environ["K_SRC_DIR"]
        run_path = tmpdir.mkdir('cycle')
        cycle_dump_path = tmpdir.mkdir('kinetic/')
        kinetic_output = cycle_dump_path.mkdir("OUTPUT")
        kinetic_input = cycle_dump_path.mkdir("INPUT")
        k_obj = SOL_KIT(run_path, str(kinetic_input), str(kinetic_output), src_dir,
                         cycle_dump_path, k_config_yml_file_path= myPath + '/test_run_dir/config.yml')
        k_obj.setFiles()
        true_qe = np.random.rand(198)
        np.savetxt(os.path.join(k_obj._sol_kit_output_path, 'HEAT_FLOW_X/HEAT_FLOW_X_01000.txt'), true_qe)
        true_qe = np.insert(true_qe, 0,0)
        true_qe = np.append(true_qe, 0) * k_obj.normalised_values['qe']
        qe = k_obj.getLastHeatFlow()
        assert true_qe[::2] == pytest.approx(qe)