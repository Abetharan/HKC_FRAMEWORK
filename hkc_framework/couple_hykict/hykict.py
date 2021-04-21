""" 
Self-consistently runs the Rad-Hydro code HyKiCT
@author = Abetharan Antony
"""
import h5py
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

    def __init__(self, run_path = None, f_src_dir = None, f_input_path = None,
                 f_config_yml_file_path = None, nx = 0, tmax = 0, start_from_kin = True, **kwargs):

        self.init = Input(f_config_yml_file_path)
        self._fluid_src_dir = f_src_dir
        self._fluid_input_path = f_input_path       
        self._run_path = run_path
        self._fluid_output_path = "" 
        self.cycle_dump_path = ""
        self._copyHyKiCT()
        self.laser_direction = 'right'
        self.tmax = float(tmax)
        self.cycle_step = float(tmax)
        self.init.yaml_file['TimeParameters']['steps'] = 0
        self.init.yaml_file['FixedParameters']['nx'] = nx 
        self.couple_mode = None
        self.cycle_no = 0
        self.curr_time = float(self.init.yaml_file['TimeParameters']['t_init'])
        for couple_mode, logic in kwargs.items():
            if logic:
                self.couple_mode = couple_mode  
                break

        if start_from_kin:
            self.init.yaml_file['TimeParameters']['t_max'] = 0
            self.init.yaml_file['Switches']['CoupleDivQ'] = False
            self.init.yaml_file['Switches']['CoupleMulti'] = False 
            self.init.yaml_file['Switches']['CoupleSubtract'] = False 
        else:
            self.init.yaml_file['TimeParameters']['t_max'] = self.tmax + self.curr_time
            self.init.yaml_file['Switches'][self.couple_mode] = True

    def revertStartKinSwitches(self):
        """
        Purpose: Sets the expected switches to HyKiCT which were
                changed due to Start Kinetic Mode. 
        """
        self.init.yaml_file['Switches'][self.couple_mode] = True
        self.setTimes()

    def setNoCoupleSwitches(self):
        """
        Purpose: Sets switches to no coupling
        """
        self.init.yaml_file['Switches']['CoupleDivQ'] = False
        self.init.yaml_file['Switches']['CoupleMulti'] = False 
        self.init.yaml_file['Switches']['CoupleSubtract'] = False 
        self.init.yaml_file['TimeParameters']['t_max'] = self.tmax

    def setOperatorSplitSwitch(self):
        """
        Purpose: Sets the expected switches to Operator-
                Split coupling
        """
        self.init.yaml_file['Switches'][self.couple_mode] = True
        self.init.yaml_file['Switches']['CoupleOperatorSplit'] = True
        self.init.yaml_file['TimeParameters']['t_max'] = self.tmax

    def setTimes(self):
        """
        Purpose: Updates times to keep track of what time fluid is running.
        """
        self.init.yaml_file['TimeParameters']['t_init'] = float(self.curr_time)
        self.init.yaml_file['TimeParameters']['t_max'] = float(self.tmax)

    def setFiles(self):
        """ Purpose: Write out config.yml for each cycle"""
        yaml_dump_path = os.path.join(self.cycle_dump_path, 'config.yml')
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
                self.cycle_dump_path+'/config.yml']
        super().Execute(cmd, self.cycle_dump_path)

    def getPhysicalRunTime(self, given_path = None):
        if given_path is None:
            last_index = findLargestIndex(os.path.join(self._fluid_output_path, "TIME"))
            return(np.loadtxt(os.path.join(self._fluid_output_path,
                    "".join(['TIME/TIME', str(last_index), ".txt"]))))
        else:
            last_index = findLargestIndex(os.path.join(given_path, "TIME"))
            return(np.loadtxt(os.path.join(given_path,
                    "".join(['TIME/TIME_', str(last_index), ".txt"]))))
    def updateTime(self, curr_time = None, new_cycle_step = None):
        if curr_time is not None:
            self.curr_time = float(curr_time)
        if new_cycle_step is None:
            self.tmax = self.curr_time + float(self.cycle_step)
        else: 
            self.tmax = self.curr_time + float(new_cycle_step)

    def getLastStepQuants(self, update_time = True): 
        """
        Purpose: Get the last output from hykict ouput
        Return: Return numpy arrays of all required files
        """
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
        f_time = np.loadtxt(os.path.join(self._fluid_output_path,
               "".join(["TIME/TIME_", str(last_index), ".txt"])), dtype = np.float64)
        mass = np.loadtxt(self._fluid_input_path + "/mass.txt")        

        f_specific_heat = np.loadtxt(os.path.join(self._fluid_output_path,
                    "".join(["ELECTRON_SPECIFIC_HEAT/ELECTRON_SPECIFIC_HEAT_", str(last_index), ".txt"])),dtype = np.float64)
        
        if update_time:
            self.updateTime(curr_time = f_time)
        
        return(f_x_grid, f_x_centered_grid, f_v, f_ne, f_Te, f_Z, f_laser, mass, f_time, f_specific_heat)
    
    def initHydro(self, next_fluid_input_path): #Remove qe etc
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

        #Standard init files. 
        #These files should correspond to the last state of the fluid step. 
        #In the case of HyKiCT these files are sufficient to init everything else. 
        #For other codes this might not be the case. 

        shutil.copyfile(os.path.join(self._fluid_output_path, "CELL_WALL_X/CELL_WALL_X_" + str(largest_fluid_index) +".txt")    
                                        ,os.path.join(next_fluid_input_path,"coord.txt"))
        shutil.copyfile(os.path.join(self._fluid_output_path,  "VELOCITY/VELOCITY_" + str(largest_fluid_index) +".txt") 
                                        ,os.path.join(next_fluid_input_path,"velocity.txt"))
        shutil.copyfile(os.path.join(self._fluid_output_path,  "DENSITY/DENSITY_" + str(largest_fluid_index) +".txt")  
                                        ,os.path.join(next_fluid_input_path,"density.txt"))
        shutil.copyfile(os.path.join(self._fluid_output_path, "ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_" + str(largest_fluid_index) +".txt")
                                        ,os.path.join(next_fluid_input_path,"electron_temperature.txt"))
        shutil.copyfile(os.path.join(self._fluid_output_path, "ION_TEMPERATURE/ION_TEMPERATURE_" + str(largest_fluid_index) +".txt")
                                        ,os.path.join(next_fluid_input_path,"ion_temperature.txt"))
        #REF IF MASS EVER CHANGES IN FLUID SIMULATION OR WE NOW USING IONISATION MODELS
        #THIS NEEDS TO BE CHANGED AND GET COPIED FROM OUTPUT
        # AR REMAINS CONSTANT AS LONG AS MATERIALS DONT CHANGE
        mass_copy_path  = os.path.join(self._fluid_input_path, "mass.txt")
        mass_new_path = os.path.join(next_fluid_input_path,"mass.txt")
        if not os.path.exists(mass_new_path):
            shutil.copyfile(mass_copy_path, mass_new_path)

        z_copy_path  = os.path.join(self._fluid_input_path, "Z.txt")
        z_new_path = os.path.join(next_fluid_input_path, "Z.txt")
        if not os.path.exists(z_new_path):
            shutil.copyfile(z_copy_path, z_new_path)

        Ar_copy_path  = os.path.join(self._fluid_input_path, "Ar.txt")
        Ar_new_path = os.path.join(next_fluid_input_path, "Ar.txt")
        if not os.path.exists(Ar_new_path):
            shutil.copyfile(Ar_copy_path, Ar_new_path)

        laser_copy_path  = os.path.join(self._fluid_input_path, "laser_profile.txt")
        laser_new_path = os.path.join(next_fluid_input_path, "laser_profile.txt")
        if os.path.exists(laser_copy_path):
            if not os.path.exists(laser_new_path):
                shutil.copyfile(laser_copy_path, laser_new_path)

        if self.init.yaml_file["Switches"]["RadiationTransport"]:
            shutil.copyfile(os.path.join(self._fluid_output_path, "RAD_ENERGY_DENSITY/RAD_ENERGY_DENSITY_" + str(largest_fluid_index) +".txt")
                                            ,os.path.join(next_fluid_input_path,"rad_energy_density.txt"))
            if self.curr_time > 0:
                self.init.yaml_file["Switches"]["LoadInRadEnergy"] = True

            
    def returnInitValues(self):
        last_index = 0
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
        f_time = np.loadtxt(os.path.join(self._fluid_output_path,
               "".join(["TIME/TIME_", str(last_index), ".txt"])), dtype = np.float64)
        mass = np.loadtxt(self._fluid_input_path + "/mass.txt")        

        return(f_x_grid, f_x_centered_grid, f_v, f_ne, f_Te, f_Z, f_laser, mass, f_time)
    def storeToHdf5(self, hdf5_file, cycle):
        """
        Purpose: Store HyKiCT in/output to centralised hdf5
        Args:
            hdf5_file = h5py file 
            cycle = cycle number
        """
        for folders in (self._fluid_input_path, self._fluid_output_path):
            with os.scandir(folders) as subdir:
                for var_dir in subdir:
                    if var_dir.is_dir():
                        with os.scandir(var_dir.path) as all_files:
                            for var_file in all_files:
                                if not var_file.name.startswith("."):
                                    tmp_read_file = np.loadtxt(var_file.path)
                                    if var_dir.name != "TIME":
                                        hdf5_file.create_dataset("".join(["Cycle_", str(cycle), 
                                                                        "/Fluid_Output/", var_dir.name,"/",
                                                                        os.path.splitext(var_file.name)[0]]),
                                                                data = tmp_read_file, compression="gzip")
                                    else:
                                        hdf5_file.create_dataset("".join(["Cycle_", str(cycle), 
                                                                        "/Fluid_Output/", var_dir.name,"/",
                                                                        os.path.splitext(var_file.name)[0]]),
                                                                data = tmp_read_file)
                    if var_dir.is_file():
                        tmp_read_file = np.loadtxt(var_dir.path)
                        hdf5_file.create_dataset("".join(["Cycle_", str(cycle), "/Fluid_Input/", var_dir.name]),
                                                data = tmp_read_file)
