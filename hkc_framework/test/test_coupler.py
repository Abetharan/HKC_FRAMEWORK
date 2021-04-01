import pytest 
import numpy as np 
import pytest
import os
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '../')
from coupler import Coupler

class TestCoupler:

    @pytest.mark.skip(reason="no way of currently testing this")
    def test_step(self, tmpdir):
        src_dir = os.environ["F_SRC_DIR"]
        p = tmpdir.mkdir('cycle')
        f_obj = HyKiCT(run_path=p,f_src_dir = src_dir, f_config_yml_file_path= myPath + "/test_run_dir/config.yml")
        f_obj.cycle_dump_path = p
        with pytest.raises(SystemExit):
            f_obj.Run()
            
    @pytest.mark.skip(reason="no way of currently testing this")
    def start_from_kinetic_test(self, tmpdir):
        ksrc_path = os.environ["K_SRC_DIR"]
        fsrc_path = os.environ["F_SRC_DIR"]
        base_path = tmpdir.mkdir('cycle')
        coupler = Coupler(os.path.join(myPath, 'test/OrdinaryTest/new_config.yml'))
        coupler.init.yaml_file['Paths']['Base_dir'] = base_path
        coupler.init.yaml_file['Paths']['K_src_dir'] = ksrc_path
        coupler.init.yaml_file['Paths']['F_src_dir'] = fsrc_path
        coupler.init.yaml_file['Paths']['Init_Path'] = os.path.join(myPath, 'OrdinaryTest/step_problem_100_nx')
        coupler.init.yaml_file['Mode']['Couple_subtract'] = True
        coupler.main()
        fpath = os.path.join(base_path, 'CYCLE_0/FLUID_OUTPUT/')
        variables = ['ELECTRON_HEAT_FLOW_X, ELECTRON_TEMPERATURE, VELOCITY, CELL_CENTRE_X, CELL_WALL_X, ION_TEMPERATURE']
        lengths = []
        for variable in variables:
            lengths.append(len(os.listdir(os.path.join(fpath, variable))))
            
        assert all([a == 1 for a in lengths])


