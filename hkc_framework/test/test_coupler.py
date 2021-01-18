import pytest 
import numpy as np 
import pytest
import os
import sys
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '../')
from new_coupler import Coupler

class TestCoupler:

    def test_multiplier(self, tmpdir):
        ksrc_path = os.environ["K_SRC_DIR"]
        fsrc_path = os.environ["F_SRC_DIR"]
        base_path = tmpdir.mkdir('cycle')
        coupler = Coupler(os.path.join(myPath, 'OrdinaryTest/new_config.yml'))
        coupler.init.yaml_file['Paths']['Base_dir'] = base_path
        coupler.init.yaml_file['Mode']['Start_from_kinetic'] = True
        coupler.init.yaml_file['Paths']['K_src_dir'] = ksrc_path
        coupler.init.yaml_file['Paths']['F_src_dir'] = fsrc_path
        coupler.init.yaml_file['Paths']['F_config_path'] = os.path.join(myPath, 'OrdinaryTest/hykict_config.yml')
        coupler.init.yaml_file['Paths']['K_config_path'] = os.path.join(myPath, 'OrdinaryTest/sol_config.yml')
        coupler.init.yaml_file['Paths']['Init_Path'] = os.path.join(myPath, 'OrdinaryTest/step_problem_100_nx')
        print(coupler.init.yaml_file)
        coupler.init.yaml_file['Mode']['Couple_multi'] = True
        coupler.main()

        kpath = os.path.join(base_path, 'CYCLE_0/KINETIC_OUTPUT/')
        fpath = os.path.join(base_path, 'CYCLE_1/FLUID_OUTPUT/')
        kheat_flow = np.loadtxt(os.path.join(kpath, 'HEAT_FLOW_X/HEAT_FLOW_X_' + str(50).zfill(5) + '.txt')) * 4.300133202205076e+20
        kheat_flow = np.insert(kheat_flow, 0, 0.)
        kheat_flow = np.append(kheat_flow, 0.)
        kheat_flow = kheat_flow[::2]
        fheat_flow = np.loadtxt(os.path.join(fpath, 'ELECTRON_HEAT_FLOW_X/ELECTRON_HEAT_FLOW_X_' + str(0) + '.txt'))* -1
        
        assert len(fheat_flow) == len(kheat_flow)
        assert all([(abs(a - b)/b) < 0.01 for a, b in zip(fheat_flow, kheat_flow)])


    @pytest.mark.skip(reason="no way of currently testing this")
    def test_divq(self, tmpdir):
        ksrc_path = os.environ["K_SRC_DIR"]
        fsrc_path = os.environ["F_SRC_DIR"]
        base_path = tmpdir.mkdir('cycle')
        coupler = Coupler(os.path.join(myPath, 'test/OrdinaryTest/new_config.yml'))
        coupler.init.yaml_file['Paths']['Base_dir'] = base_path
        coupler.init.yaml_file['Paths']['K_src_dir'] = ksrc_path
        coupler.init.yaml_file['Paths']['F_src_dir'] = fsrc_path
        coupler.init.yaml_file['Paths']['Init_Path'] = os.path.join(myPath, 'OrdinaryTest/step_problem_100_nx')
        coupler.init.yaml_file['Mode']['Couple_divq'] = True
        coupler.main()

        kpath = os.path.join(base_path, 'CYCLE_0/KINETIC_OUTPUT/')
        fpath = os.path.join(base_path, 'CYCLE_1/FLUID_OUTPUT/')
        kheat_flow = np.loadtxt(os.path.join(kpath, 'HEAT_FLOW_X/HEAT_FLOW_X_' + str(50).zfill(5) + '.txt')) * 4.300133202205076e+20
        kheat_flow = np.insert(kheat_flow, 0, 0.)
        kheat_flow = np.append(kheat_flow, 0.)
        kheat_flow = kheat_flow[::2]
        fheat_conduc = np.loadtxt(os.path.join(fpath, 'ELECTRON_HEAT_CONDUCTION/ELECTRON_CONDUCTION' + str(0) + '.txt'))
        mass = np.loadtxt(os.path.join(myPath, 'OrdinaryTest/step_problem_100_nx/mass.txt'))

        nx = len(kheat_flow) 
        KHeatConductionE = np.zeros(nx - 1)
        for i, m in enumerate(mass):
            KHeatConductionE[i] = (-(kheat_flow[i + 1] - kheat_flow[i])) / m #- sign is there because of convention used in HyKiCT 

        assert len(fheat_conduc) == len(KHeatConductionE)
        assert all([(abs(a - b)/b) < 0.01 for a, b in zip(fheat_conduc, KHeatConductionE)])
    
    @pytest.mark.skip(reason="no way of currently testing this")    
    def test_subtract(self, tmpdir):
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

        kpath = os.path.join(base_path, 'CYCLE_0/KINETIC_OUTPUT/')
        fpath = os.path.join(base_path, 'CYCLE_1/FLUID_OUTPUT/')
        kheat_flow = np.loadtxt(os.path.join(kpath, 'HEAT_FLOW_X/HEAT_FLOW_X_' + str(50).zfill(5) + '.txt')) * 4.300133202205076e+20
        kheat_flow = np.insert(kheat_flow, 0, 0.)
        kheat_flow = np.append(kheat_flow, 0.)
        kheat_flow = kheat_flow[::2]
        fheat_flow = np.loadtxt(os.path.join(fpath, 'ELECTRON_HEAT_FLOW_X/ELECTRON_HEAT_FLOW_X_' + str(0) + '.txt'))* -1
        
        assert len(fheat_flow) == len(kheat_flow)
        assert all([(abs(a - b)/b) < 0.01 for a, b in zip(fheat_flow, kheat_flow)])

    # def leap_frog_divq(self):

    # def leap_frog_multi_test(self):
    # def leap_frog_subtract_test(self):

    # def os_div_q_test(self):
    # def os_multi_test(self):
    # def os_subtract_test(self):
    
    # def continue_test(self):
    # def overwrite_test(self):
    # def limit_density_test(self, tmpdir):
    #     ksrc_path = os.environ["K_SRC_DIR"]
    #     fsrc_path = os.environ["F_SRC_DIR"]
    #     base_path = tmpdir.mkdir('cycle')
    #     coupler = Coupler(os.path.join(myPath, 'test/Hohlraum_test/new_config.yml'))
    #     coupler.init.yaml_file['Paths']['Base_dir'] = base_path
    #     coupler.init.yaml_file['Paths']['K_src_dir'] = ksrc_path
    #     coupler.init.yaml_file['Paths']['F_src_dir'] = fsrc_path
    #     coupler.init.yaml_file['Paths']['Init_Path'] = os.path.join(myPath, 'OrdinaryTest/step_problem_100_nx')
    #     coupler.init.yaml_file['Mode']['Couple_subtract'] = True
    #     coupler.main()

    #     kpath = os.path.join(base_path, 'CYCLE_0/KINETIC_OUTPUT/')
    #     fpath = os.path.join(base_path, 'CYCLE_1/FLUID_OUTPUT/')
    #     kheat_flow = np.loadtxt(os.path.join(kpath, 'HEAT_FLOW_X/HEAT_FLOW_X_' + str(50).zfill(5) + '.txt')) * 4.300133202205076e+20
    #     kheat_flow = np.insert(kheat_flow, 0, 0.)
        # kheat_flow = np.append(kheat_flow, 0.)
        # kheat_flow = kheat_flow[::2]
        # fheat_flow = np.loadtxt(os.path.join(fpath, 'ELECTRON_HEAT_FLOW_X/ELECTRON_HEAT_FLOW_X_' + str(0) + '.txt'))* -1
        
        # assert len(fheat_flow) == len(kheat_flow)
        # assert all([(abs(a - b)/b) < 0.01 for a, b in zip(fheat_flow, kheat_flow)])
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


