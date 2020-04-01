"""
Handles anything to do with input/output and paths.
Creates folders and retains paths required to couple.
@author = Abetharan Antony
Last Update = 25/11/19
"""
import os
import shutil
import shutil
from distutils.dir_util import copy_tree
class IO:

    def __init__(self, base_dir_, k_src_dir_, run_name_, f_src_dir_, f_init_path_, f_feos_1_, f_feos_2_, cycle_counter_, overwrite_,last_cycle_, initialise_all_folders_, *arg):

        self._BASE_DIR = base_dir_
        self._K_SRC_DIR = k_src_dir_
        self._RUN_NAME = run_name_
        self._F_SRC_DIR = f_src_dir_
        self._F_SWITCH_PATH = None
        self._F_INIT_PATH = f_init_path_
        self._F_FEOS_1_PATH = f_feos_1_
        self._F_FEOS_2_PATH = f_feos_2_

        self.last_cycle = last_cycle_        
        self._RUN_PATH = os.path.join(self._BASE_DIR, self._RUN_NAME)
        self.cycle_counter = cycle_counter_
        self.cycle_dump_path = None
        self.fluid_input_path = None
        self.fluid_output_path = None
        self.kinetic_input_path = None
        self.kinetic_output_path = None
        self.next_fluid_input_path = None

        self.__OVERWRITE_OK = overwrite_
        self.preserved_cycle_path = []
        self.preserved_fluid_input_path = []
        self.preserved_fluid_output_path = []  
        self.preserved_kinetic_input_path = []
        self.preserved_kinetic_output_path = []
        #self.setPaths()
        if initialise_all_folders_:
            if len(arg) > 1:
                print("Only one entry allowed which is Number of cycles. Try again")
                import sys
                sys.exit(1)

            self.createDirectoryOfOperation(int(arg[0]))
    
        self.nextCyclePathManager()
    
    def zipAndDelete(self):
        """ Purpose: To zip cycle folders to reduce memory overhead and delete the folder"""
        zip_location = os.path.join(self._RUN_PATH, 'CYCLE_' + str(self.cycle_counter))
        shutil.make_archive(zip_location, 'zip', self.cycle_dump_path)
        shutil.rmtree(self.cycle_dump_path)


    def setPaths(self):
        """ Purpose: Create all paths based on which cycle is going to run."""
        print("#"*100)
        print('\033[1m' + '#'*50 + ' SET PATHS ' + '#'*50 + '\033[0m')
        print('\n')
        self.cycle_dump_path = os.path.join(self._RUN_PATH, ("CYCLE_" + str(self.cycle_counter)))
        self.fluid_input_path = os.path.join(self.cycle_dump_path + "/FLUID_INPUT/")
        self.fluid_output_path = os.path.join(self.cycle_dump_path + "/FLUID_OUTPUT/")
        self.kinetic_input_path = os.path.join(self.cycle_dump_path + "/KINETIC_INPUT/")
        self.kinetic_output_path = os.path.join(self.cycle_dump_path + "/KINETIC_OUTPUT/")

    def copyFluidInit(self,f_init_path_):
        """ 
        OBSOLTE --- DELETE
        Purpose: Copies fluid input folder
        """
        print("#"*100)
        print('\033[1m' + '#'*50 + ' COPYING INIT FLUID ' + '#'*50 + '\033[0m')
        print('\n')   
        copy_tree(f_init_path_, self.fluid_input_path)
        
    
    def createDirectoryOfOperation(self, no_cycles_):
        """ Purpose: Creates all folders required immediately .. reduces overhead later one 
            Args: no_cycles = no of total cycles.
        """
        print("#"*100)
        print('\033[1m' + '#'*50 + ' CREATING ALL FOLDERS ' + '#'*50 + '\033[0m')
        print('\n')
        path_exists = False
        if os.path.exists(self._RUN_PATH):
            path_exists = True
            if self.__OVERWRITE_OK:
                shutil.rmtree(self._RUN_PATH)
        else:
            #create base folder
            os.makedirs(self._RUN_PATH)
        for i in range(no_cycles_):
            cycle_path = os.path.join(self._RUN_PATH, ("CYCLE_" + str(i)))
        
            fluid_input_path = os.path.join(cycle_path + "/FLUID_INPUT/")
            fluid_output_path = os.path.join(cycle_path + "/FLUID_OUTPUT/")
            kinetic_input_path = os.path.join(cycle_path + "/KINETIC_INPUT/")
            kinetic_output_path = os.path.join(cycle_path + "/KINETIC_OUTPUT/")
            if path_exists:
                if os.path.exists(fluid_input_path) and os.path.exists(fluid_output_path) and os.path.exists(kinetic_output_path) and os.path.exists(kinetic_output_path):
                    continue
            
            os.makedirs(cycle_path)
            os.makedirs(fluid_output_path)
            os.makedirs(fluid_input_path)
            os.makedirs(kinetic_input_path)
            os.makedirs(kinetic_output_path)
    
    def nextCyclePathManager(self):
        """
        Purpose: Creates the folders for the next cycle and calls any file moving required before the start of new time step. NAMELY moving IMPACT files.
        Args:
            runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
            cycleStep = cycle number.
        """

        self.cycle_dump_path = self._RUN_PATH + "/CYCLE_" + str(self.cycle_counter)
        self._F_SWITCH_PATH = os.path.join(self.cycle_dump_path, "/HydroSwitches.txt")
        self.fluid_input_path = os.path.join(self.cycle_dump_path + "/FLUID_INPUT/")
        self.fluid_output_path = os.path.join(self.cycle_dump_path + "/FLUID_OUTPUT/")
        self.kinetic_input_path = os.path.join(self.cycle_dump_path + "/KINETIC_INPUT/")
        self.kinetic_output_path = os.path.join(self.cycle_dump_path + "/KINETIC_OUTPUT/")
        new_cycle_path = self._RUN_PATH + "/CYCLE_" + str(self.cycle_counter + 1)
        self.next_fluid_input_path = os.path.join(new_cycle_path, "FLUID_INPUT/")
        self._F_OUT_PATH = self.fluid_output_path
        if not os.path.exists(self.next_fluid_input_path) and not self.last_cycle:
            import sys
            print("NEXT FLUID PATH HAS NOT BEEN CREATED")
            print(self.preserved_fluid_input_path)
            sys.exit(1)
        
        self.preservePaths()

    #     if self.cycle_counter > 0:
    #         self.previous_cycle_dump_path = self._RUN_PATH + \
    #         "cycle_" + str(self.cycle_counter - 1) + "/"
    #         self.previous_fluid_input_path = os.path.join(self.previous_cycle_dump_path + "/fluid_input/")
    #         self.previous_fluid_output_path = os.path.join(self.previous_cycle_dump_path + "/fluid_output/")
    #         self.previous_kinetic_output_path = os.path.join(self.previous_cycle_dump_path + "/kinetic_output/")
    #         self.previous_kinetic_input_path = os.path.join(self.previous_cycle_dump_path + "/kinetic_input/")
    #         self.moveIMPACTFile()

    def returnCurrentPaths(self):
        return(self.cycle_dump_path, self.fluid_input_path, 
        self.fluid_output_path, self.kinetic_input_path,
         self.kinetic_output_path)
    
    # def returnPreviousPaths(self):
    #     return(self.previous_cycle_dump_path, self.previous_fluid_input_path, 
    #     self.previous_fluid_output_path, self.previous_kinetic_input_path,
    #     self.previous_kinetic_output_path)

    def preservePaths(self, unit_test = False):

        self.preserved_cycle_path.append(self.cycle_dump_path)
        self.preserved_fluid_input_path.append(self.fluid_input_path)
        self.preserved_fluid_output_path.append(self.fluid_output_path)
        self.preserved_kinetic_input_path.append(self.kinetic_input_path)
        self.preserved_kinetic_output_path.append(self.kinetic_output_path)

        if unit_test:
            length_cycle_p = len(self.preserved_cycle_path)
            length_f_i_p = len(self.preserved_fluid_input_path)
            length_f_o_p = len(self.preserved_fluid_output_path)
            length_k_i_p = len(self.preserved_kinetic_input_path)
            length_k_o_p = len(self.preserved_kinetic_output_path)
            length_arrays = [length_cycle_p, length_f_i_p, length_f_o_p, length_k_i_p, length_k_o_p]
            if any(length_arrays != self.cycle_counter):
                import sys
                print("path does not match length of array")
                sys.exit(1)