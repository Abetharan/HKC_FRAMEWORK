""" 
Self-consistently runs the Rad-Hydro code HyKiCT
@author = Abetharan Antony
Last Update Date = 01/04/2020
"""
import os
import yaml
import shutil
import pathlib
import subprocess
import numpy as np 
from hkc_framework.common.fluid import Fluid
from hkc_framework.common.input import Input

class HyKiCT(Fluid):

    def __init__(self, io, couple_divq, couple_multi, start_from_kinetic):
        config_yml_file_path = os.path.join(
                                pathlib.Path(__file__).parent.absolute(),
                                'config.yml')
        self.init = Input(config_yml_file_path)

        self._fluid_src_dir = io._F_SRC_DIR
        self._cycle_path = io.cycle_dump_path
        self._init_file_path = io.fluid_input_path       
        self._base_dir = io._BASE_DIR
        self._run_path = io._RUN_PATH
        self._out_file_path = io.fluid_output_path
        self._cycle_dump_path = io.cycle_dump_path
        self._pre_heat_start_index = 0
        self._pre_heat_last_index = 0
        self._front_heat_start_index = 0
        self._front_heat_last_index = 0 
        self._nt = 0 
        self._tmax = 0  
        if not start_from_kinetic:
            self._nt = self.init.yaml_file['TimeParameters']['steps']
            if self._nt == 0:
                self._tmax = self.init.yaml_file['TimeParameters']['t_max']
        
        self.copyHyKiCT()            
    
    def writeYaml(self):
        """ Purpose: Write out config.yml for each cycle"""

        yaml_dump_path = os.path.join(self._cycle_dump_path, 'config.yml')
        yaml.dump(self.init.yaml_file, yaml_dump_path)
        
    def copyHyKiCT(self):
        """ Purpose: Copy HyKiCT exe. """

        if not os.path.exists(os.path.join(self._run_path, 'HyKiCT')):
            shutil.copy(self._fluid_src_dir+ '/HyKiCT', self._run_path)
    
    def hykictRun(self):
        """ Purpose: Run HyKiCT with parameters set previously"""

        os.chdir(self._run_path)
        cmd = ['./ELH1','-p',
                self._cycle_dump_path+'/config.yml']
        super().Execute(cmd, self._cycle_path)
    
