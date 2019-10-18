import os
import shutil

class IO:

    def __init__(self, base_dir_, init_path_, src_dir_, cycle_counter_ = 0):
        
        self._BASE_PATH = base_dir_
        self._INIT_PATH = init_path_
        self._SRC_DIR = src_dir_
        self._RUN_PATH = os.path.join(self._BASE_PATH, self._RUN_PATH)
        self.cycle_counter = cycle_counter_
        self.cycle_dump_path = None
        self.fluid_input_path = None
        self.fluid_output_path = None
        self.kinetic_input_path = None
        self.kinetic_output_path = None
        self.previous_cycle_dump_path = None
        self.previous_fluid_input_path = None
        self.previous_fluid_output_path = None
        self.previous_kinetic_output_path = None
        self.previous_kinetic_input_path = None
        #os.environ["BASEDIR"] + os.environ["RUN"] + "/"
    
    @classmethod
    def moveKineticFiles(Kinetic, run_path_, cycle_dump_path_,
                        previous_kinetic_input_path_, previous_kinetic_output_path_):
        
        Kinetic.moveFile(run_path_, cycle_dump_path_,
                        previous_kinetic_input_path_, previous_kinetic_output_path_)

    def NextCycleFileManager(self):
        """
        Purpose: Creates the folders for the next cycle and calls any file moving required before the start of new time step. NAMELY moving IMPACT files.
        Args:
            runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
            cycleStep = cycle number.
        """

        self.cycle_dump_path = self._RUN_PATH + "cycle_" + str(self.cycle_counter)
        self.fluid_input_path = os.path.join(self.cycle_dump_path + "/fluid_input/")
        self.fluid_output_path = os.path.join(self.cycle_dump_path + "/fluid_output/")
        self.kinetic_input_path = os.path.join(self.cycle_dump_path + "/kinetic_input/")
        self.kinetic_output_path = os.path.join(self.cycle_dump_path + "/kinetic_output/")

        if not os.path.exists(self.cycle_dump_path):
            os.makedirs(self.cycle_dump_path)
            os.makedirs(self.fluid_input_path)
            os.makedirs(self.fluid_output_path)
            os.makedirs(self.kinetic_input_path)
            os.makedirs(self.kinetic_output_path)
        if os.path.exists(self.cycle_dump_path):
            if not os.path.exists(self.fluid_input_path):
                os.makedirs(self.kinetic_output_path)
            if not os.path.exists(self.fluid_output_path):
                os.makedirs(self.fluid_output_path)
            if not os.path.exists(self.kinetic_input_path):
                os.makedirs(self.kinetic_output_path)
            if not os.path.exists(self.kinetic_output_path):
                os.makedirs(self.kinetic_output_path)

        if self.cycle_counter > 0:
            self.previous_cycle_dump_path = self._RUN_PATH + \
            "cycle_" + str(self.cycle_counter - 1) + "/"
            self.previous_fluid_input_path = os.path.join(self.previous_cycle_dump_path + "/fluid_input/")
            self.previous_fluid_output_path = os.path.join(self.previous_cycle_dump_path + "/fluid_output/")
            self.previous_kinetic_output_path = os.path.join(self.previous_cycle_dump_path + "/kinetic_output/")
            self.previous_kinetic_input_path = os.path.join(self.previous_cycle_dump_path + "/kinetic_input/")
            self.moveKineticFiles(self._RUN_PATH, self.cycle_dump_path,
                        self.previous_kinetic_input_path, self.previous_kinetic_output_path)

    def returnCurrentPaths(self):
        return(self.cycle_dump_path, self.fluid_input_path, 
        self.fluid_output_path, self.kinetic_input_path,
         self.kinetic_output_path)
    
    def returnPreviousPaths(self):
        return(self.previous_cycle_dump_path, self.previous_fluid_input_path, 
        self.previous_fluid_output_path, self.previous_kinetic_input_path,
        self.previous_kinetic_output_path)