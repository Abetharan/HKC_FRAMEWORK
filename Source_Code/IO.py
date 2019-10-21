import os
import shutil

class IO:

    def __init__(self, base_dir_, k_src_dir_, run_name_, f_src_dir_, f_init_path_, f_feos_1_, f_feos_2_, cycle_counter_ = 0):

        self._BASE_PATH = base_dir_
        self._K_SRC_DIR = k_src_dir_
        self._RUN_NAME = run_name_
        self._F_SRC_DIR = f_src_dir_
        self._F_SWITCH_PATH = None
        self._F_INIT_PATH = f_init_path_
        self._F_FEOS_1_PATH = f_feos_1_
        self._F_FEOS_2_PATH = f_feos_2_
        
        self._RUN_PATH = os.path.join(self._BASE_PATH, self._RUN_NAME)
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
 
    def moveIMPACTFile(self):
        """ 
        Purpose: Moves IMPACT files and initial parameters files to correct paths.
        Args:
            runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
            cycleDumpPath = cycle path. Path is BASEDIR/runName/cycleNAME
            previsouKineticInputPath = previous cycle kinetic input folder which is located as followeing runPath/previsousCyclePath/kinetic_input
            previousKineticOutPut =  previous cycle kinetic output folder which is located as followeing runPath/previsousCyclePath/kinetic_output
        """

        filenames = ['ionden', 'rad_to_electron',
                    'xf', 'eden', 'laserdep', 'tmat', 'zstar']

        for name in filenames:
            if os.path.splitext(self.cycle_dump_path + "/" + self._RUN_NAME + "_" + name + ".xy")[-1] != ".xy":
                continue
            shutil.move(self._RUN_PATH + "/" + self._RUN_PATH + "_" + name + ".xy",
                        self.previous_kinetic_input_path + "/" + self._RUN_PATH + "_" + name + ".xy", )

        for file in os.listdir(self._RUN_NAME):
            _, extension = os.path.splitext(file)
            if extension == ".xy" or extension == ".xyz" or extension == ".xyv" or extension == ".xyt" or extension == ".dat" or extension == ".t":
                shutil.move(file, self.previous_kinetic_output_path)
    
    def NextCycleFileManager(self):
        """
        Purpose: Creates the folders for the next cycle and calls any file moving required before the start of new time step. NAMELY moving IMPACT files.
        Args:
            runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
            cycleStep = cycle number.
        """

        self.cycle_dump_path = self._RUN_PATH + "cycle_" + str(self.cycle_counter)
        self._F_SWITCH_PATH = os.path.join(self.cycle_dump_path, "/HydroSwitches.txt")
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
            self.moveIMPACTFile()

    def returnCurrentPaths(self):
        return(self.cycle_dump_path, self.fluid_input_path, 
        self.fluid_output_path, self.kinetic_input_path,
         self.kinetic_output_path)
    
    def returnPreviousPaths(self):
        return(self.previous_cycle_dump_path, self.previous_fluid_input_path, 
        self.previous_fluid_output_path, self.previous_kinetic_input_path,
        self.previous_kinetic_output_path)