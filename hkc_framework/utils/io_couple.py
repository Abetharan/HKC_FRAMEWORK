"""
Handles anything to do with input/output and paths.
Creates folders and retains paths required to couple.
@author = Abetharan Antony
Last Update = 25/11/19
"""
import atexit
import h5py
import os
import signal
import shutil
import string
import sys
from distutils.dir_util import copy_tree
class IO:

    def __init__(self, run_name, base_dir, k_src_dir, f_src_dir, f_init_path,
                    cycle_counter, max_cycle, overwrite, use_hdf5 = True):

        self._base_dir = base_dir
        self._k_src_dir = k_src_dir
        self._run_name = run_name
        self._f_src_dir = f_src_dir
        self._f_switch_path = None
        self._f_init_path = f_init_path
        self._fast_tmp_base_dir = os.environ["BASE_PATH"] 
        self.max_cycle = max_cycle        
        self._run_path = os.path.join(self._base_dir, self._run_name)
        self._fast_run_path = self._fast_tmp_base_dir
        self.cycle_counter = cycle_counter
        self.cycle_dump_path = None
        self.fluid_input_path = None
        self.fluid_leap_frog_path = None
        self.fluid_output_path = None
        self.kinetic_input_path = None
        self.kinetic_output_path = None
        self.next_fluid_input_path = None
        self.__use_hdf5 = use_hdf5
        self.__overwrite_ok = overwrite
        self.all_cycle_path = []
        self.all_fluid_input_path = []
        self.all_fluid_output_path = []  
        self.all_kinetic_input_path = []
        self.all_kinetic_output_path = []
        self.leap_frog = False 

        #self.createDirectoryOfOperation(self.max_cycle)
    
        #self.nextCyclePathManager()
    
    def zipAndDelete(self):
        """ Purpose: To zip cycle folders to reduce memory overhead and delete the folder"""
        zip_location = os.path.join(self._run_path, 'CYCLE_' + str(self.cycle_counter))
        shutil.make_archive(zip_location, 'zip', self.cycle_dump_path)
        shutil.rmtree(self.cycle_dump_path)


    def setPaths(self):
        """ Purpose: Create all paths based on which cycle is going to run."""
        print("#"*100)
        print('\033[1m' + '#'*50 + ' SET PATHS ' + '#'*50 + '\033[0m')
        print('\n')
        self.cycle_dump_path = os.path.join(self._run_path,
                                "".join(["CYCLE_", str(self.cycle_counter)]))
        self.fluid_input_path = os.path.join(self.cycle_dump_path, "FLUID_INPUT/")
        self.fluid_output_path = os.path.join(self.cycle_dump_path, "FLUID_OUTPUT/")
        self.kinetic_input_path = os.path.join(self.cycle_dump_path, "KINETIC_INPUT/")
        self.kinetic_output_path = os.path.join(self.cycle_dump_path, "KINETIC_OUTPUT/")
        if self.leap_frog:
            self.fluid_leap_frog_path = os.path.join(self.cycle_dump_path, "LEAP_FROG_OUTPUT")
    def createDirectoryOfOperation(self,):
        """ Purpose: Creates all folders required immediately .. reduces overhead later one 
            Args: no_cycles = no of total cycles.
        """
        print("#"*100)
        print('\033[1m' + '#'*50 + ' CREATING ALL FOLDERS ' + '#'*50 + '\033[0m')
        print('\n')
        path_exists = False
        if os.path.exists(self._run_path):
            path_exists = True
            if self.__overwrite_ok:
                shutil.rmtree(self._run_path)
        else:
            #create base folder
            os.makedirs(self._run_path)
        for i in range(self.max_cycle):
            cycle_path = os.path.join(self._run_path, ("CYCLE_" + str(i)))
        
            fluid_input_path = os.path.join(cycle_path, "FLUID_INPUT/")
            fluid_output_path = os.path.join(cycle_path, "FLUID_OUTPUT/")
            kinetic_input_path = os.path.join(cycle_path, "KINETIC_INPUT/")
            kinetic_output_path = os.path.join(cycle_path, "KINETIC_OUTPUT/")
            fluid_leap_frog_path = os.path.join(cycle_path, "LEAP_FROG_OUTPUT/")
            self.all_cycle_path.append(cycle_path)
            self.all_fluid_input_path.append(fluid_input_path)
            self.all_fluid_output_path.append(fluid_output_path)
            self.all_kinetic_input_path.append(kinetic_input_path)
            self.all_kinetic_output_path.append(kinetic_output_path)
            if path_exists:
                if self.leap_frog:

                    if (os.path.exists(fluid_input_path) and
                        os.path.exists(fluid_output_path) and
                        os.path.exists(kinetic_output_path) and
                        os.path.exists(kinetic_output_path) and 
                        os.path.exists(fluid_leap_frog_path)):
                            continue
                        
                elif(os.path.exists(fluid_input_path) and
                    os.path.exists(fluid_output_path) and
                    os.path.exists(kinetic_output_path) and
                    os.path.exists(kinetic_output_path)):
                    continue
	    
            
            os.makedirs(cycle_path)
            os.makedirs(fluid_output_path)
            os.makedirs(fluid_input_path)
            os.makedirs(kinetic_input_path)
            os.makedirs(kinetic_output_path)
            os.makedirs(fluid_leap_frog_path)
    
    def nextCyclePathManager(self):
        """
        Purpose: Creates the folders for the next cycle and calls any file moving required before the start of new time step. NAMELY moving IMPACT files.
        Args:
            runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
            cycleStep = cycle number.
        """

        self.setPaths()
        new_cycle_path = os.path.join(self._run_path, "".join(["CYCLE_", str(self.cycle_counter + 1)]))
        self.next_fluid_input_path = os.path.join(new_cycle_path, "FLUID_INPUT/")

        if not os.path.exists(self.next_fluid_input_path) and not (self.max_cycle - 1  == self.cycle_counter): 
            print("NEXT FLUID PATH HAS NOT BEEN CREATED")
            sys.exit(0)
        

    def returnCurrentPaths(self):
        return(self.cycle_dump_path, self.fluid_input_path, 
        self.fluid_output_path, self.kinetic_input_path,
         self.kinetic_output_path)
    
    def deleteAll(self):
        """
        Purpose: Delete all directories created IF
                hdf5 file storage is used
        NOTE: Slightly redundant general method should
            be changed
        """
        log_path = os.path.join(self._run_path, "Logs")
        os.makedirs(log_path)
        log_cycle_paths = []
        for i in range(self.max_cycle):
            log_cy_path = os.path.join(log_path, "".join(["CYCLE_",str(i)]))
            log_cycle_paths.append(log_cy_path)            
            os.makedirs(log_cy_path)

        for i, cycle_path in enumerate(self.all_cycle_path):
            fluid_log_path = os.path.join(cycle_path, "fluid.log")
            kinetic_log_path = os.path.join(cycle_path, "kinetic.log")
            if os.path.exists(fluid_log_path):
                shutil.move(fluid_log_path, log_cycle_paths[i])
            if os.path.exists(kinetic_log_path):
                shutil.move(kinetic_log_path, log_cycle_paths[i])

        for cycle_path in self.all_cycle_path:
            shutil.rmtree(cycle_path)
        
    def _createHDF5(self):
        """
        Purpose: Create hdf5 file
        """
        self.hdf5_file = h5py.File(os.path.join(self._run_path,
                                 "".join([self._run_name, "_Data.hdf5"])), "a")
        def _closeHDF5():
            """
            Purpose: close the file upon exiting program.
            """
            self.hdf5_file.close()
        atexit.register(_closeHDF5)
