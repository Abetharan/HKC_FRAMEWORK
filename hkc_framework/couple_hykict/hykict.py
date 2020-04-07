""" 
Self-consistently runs the Rad-Hydro code HyKiCT
@author = Abetharan Antony
"""
import os
import yaml
import shutil
import pathlib
import subprocess
import numpy as np 
import yaml
from utils import findLargestIndex
from utils import Fluid
from utils import Input

class HyKiCT(Fluid):

    def __init__(self, io):
        config_yml_file_path = os.path.join(
                                pathlib.Path(__file__).parent.absolute(),
                                'config.yml')
        self.init = Input(config_yml_file_path)

        self._fluid_src_dir = io._f_src_dir
        self._cycle_path = io.cycle_dump_path
        self._init_file_path = io.fluid_input_path       
        self._base_dir = io._base_dir
        self._run_path = io._run_path
        self._fluid_output_path = io.fluid_output_path
        self._cycle_dump_path = io.cycle_dump_path
        self._copyHyKiCT()            
    
    def setFiles(self):
        """ Purpose: Write out config.yml for each cycle"""

        yaml_dump_path = os.path.join(self._cycle_dump_path, 'config.yml')
        with open(yaml_dump_path, 'w') as outfile: 
            yaml.dump(self.init.yaml_file, outfile)
        
    def _copyHyKiCT(self):
        """ Purpose: Copy HyKiCT exe. """

        if not os.path.exists(os.path.join(self._run_path, 'HyKiCT')):
            shutil.copy(self._fluid_src_dir+ '/HyKiCT', self._run_path)
    
    def Run(self):
        """ Purpose: Run HyKiCT with parameters set previously"""

        os.chdir(self._run_path)
        cmd = ['./HyKiCT','-p',
                self._cycle_dump_path+'/config.yml']
        super().Execute(cmd, self._cycle_path)
    
    def getLastStepQuants(self): 

        last_index = findLargestIndex(os.path.join(self._fluid_output_path, "ELECTRON_TEMPERATURE"))
        f_x_grid = np.loadtxt(os.path.join(self._fluid_output_path, 
                            "".join(["CELL_WALL_X/CELL_WALL_X_", str(last_index), ".txt"])), dtype = np.float64)
        f_x_centered_grid = np.loadtxt(os.path.join(self._fluid_output_path,
                            "".join(["CELL_CENTRE_X/CELL_CENTRE_X_", str(last_index), ".txt"])), dtype = np.float64)
        f_v = np.loadtxt(os.path.join(self._fluid_output_path,
                "".join(["VELOCITY/VELOCITY_", str(last_index), ".txt"])), dtype = np.float64)
        f_ne = np.loadtxt(os.path.join(self._fluid_output_path,
                "".join(["ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_", str(last_index), ".txt"])), dtype = np.float64)
        f_Te = np.loadtxt(os.path.join(self._fluid_output_path,
                "".join(["ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_", str(last_index), ".txt"])),dtype = np.float64)
        f_laser = np.loadtxt(os.path.join(self._fluid_output_path,
                    "".join(["INVERSE_BREM/INVERSE_BREM_", str(last_index), ".txt"])),dtype = np.float64)
#legacy        # f_brem = np.loadtxt(self._fluid_output_path + "BREM/BREM_" + str(last_index) + ".txt")
        f_Z = np.loadtxt(os.path.join(self._fluid_output_path,
               "".join(["ZBAR/ZBAR_", str(last_index), ".txt"])), dtype = np.float64)
        mass = np.loadtxt(self._init_file_path + "/mass.txt")        
        
        return(f_x_grid, f_x_centered_grid, f_v, f_ne, f_Te, f_Z, f_laser, mass)
    
    def initHydroFromKinetic(self, next_fluid_input_path, qe, pre_params = None, front_params = None):
        """
        Purpose: 
            Move fluid files
        Args: 
            next_fluid_input_path = path for next cycle init of HyKiCT.
            qe = q_vfp/q_sh or div.q
            pre_params = Pre-Heat fit parameters, default is None.
            front_params = Front-Heat fit parameters, default is None.
        """ 

        largest_fluid_index = findLargestIndex(os.path.join(self._fluid_output_path, "ELECTRON_TEMPERATURE"))
        #Properitary names for HyKiCT. Looks for these files when run in coupled mode
        np.savetxt(os.path.join(next_fluid_input_path,"qe.txt"), qe)

        if pre_params is not None:
            np.savetxt(os.path.join(next_fluid_input_path,"pre_heat_fit_params.txt"), pre_params)
            np.savetxt(os.path.join(next_fluid_input_path,"front_heat_fit_params.txt"), front_params)

        #Standard init files. 
        #These files should correspond to the last state of the fluid step. 
        #In the case of HyKiCT these files are sufficient to init everything else. 
        #For other codes this might not be the case. 
        #NOTE Maybe require radiation density to be loaded in. 

        shutil.copyfile(os.path.join(self._fluid_output_path + "CELL_WALL_X/CELL_WALL_X_" + str(largest_fluid_index) +".txt")    
                                        ,os.path.join(next_fluid_input_path,"coord.txt"))
        shutil.copyfile(os.path.join(self._fluid_output_path + "VELOCITY/VELOCITY_" + str(largest_fluid_index) +".txt") 
                                        ,os.path.join(next_fluid_input_path,"velocity.txt"))
        shutil.copyfile(os.path.join(self._fluid_output_path + "DENSITY/DENSITY_" + str(largest_fluid_index) +".txt")  
                                        ,os.path.join(next_fluid_input_path,"density.txt"))
        shutil.copyfile(os.path.join(self._fluid_output_path + "ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_" + str(largest_fluid_index) +".txt")
                                        ,os.path.join(next_fluid_input_path,"electron_temperature.txt"))
        shutil.copyfile(os.path.join(self._fluid_output_path + "ION_TEMPERATURE/ION_TEMPERATURE_" + str(largest_fluid_index) +".txt")
                                        ,os.path.join(next_fluid_input_path,"ion_temperature.txt"))
        shutil.copyfile(os.path.join(self._init_file_path + "mass.txt")
                                        ,os.path.join(next_fluid_input_path,"mass.txt"))
        shutil.copyfile(os.path.join(self._init_file_path + "Z.txt")
                                        ,os.path.join(next_fluid_input_path, "Z.txt"))
        shutil.copyfile(os.path.join(self._init_file_path + "Ar.txt")  
                                    ,os.path.join(next_fluid_input_path, "Ar.txt"))
