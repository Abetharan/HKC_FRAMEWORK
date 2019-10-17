import os
import shutil
class IO:

    def __init__(self, basedir, initpath):
        
        self._BASE_PATH = basedir
        self._INIT_PATH = initpath
        self._RUN_PATH = os.path.join(self._BASE_PATH, self._RUN_PATH)
        #os.environ["BASEDIR"] + os.environ["RUN"] + "/"
    
    def moveIMPACTFILE(self, runPath, cycleDumpPath, previousKineticInputPath, previousKineticOutputPath):
        """ 
        Purpose: Moves IMPACT files and initial parameters files to correct paths.
        Args:
            runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
            cycleDumpPath = cycle path. Path is BASEDIR/runName/cycleNAME
            previsouKineticInputPath = previous cycle kinetic input folder which is located as followeing runPath/previsousCyclePath/kinetic_input
            previousKineticOutPut =  previous cycle kinetic output folder which is located as followeing runPath/previsousCyclePath/kinetic_output
        """

        runName = os.environ['RUN']
        filenames = ['ionden', 'rad_to_electron',
                    'xf', 'eden', 'laserdep', 'tmat', 'zstar']

        for name in filenames:
            if os.path.splitext(cycleDumpPath + "/" + runName + "_" + name + ".xy")[-1] != ".xy":
                continue
            shutil.move(runPath + "/" + runName + "_" + name + ".xy",
                        previousKineticInputPath + "/" + runName + "_" + name + ".xy", )

        for file in os.listdir(runPath):
            _, extension = os.path.splitext(file)
            if extension == ".xy" or extension == ".xyz" or extension == ".xyv" or extension == ".xyt" or extension == ".dat" or extension == ".t":
                shutil.move(file, previousKineticOutputPath)

    def NextCycleFileManager(self, runPath, cycleStep):
        """
        Purpose: Creates the folders for the next cycle and calls any file moving required before the start of new time step. NAMELY moving IMPACT files.
        Args:
            runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
            cycleStep = cycle number.
        """

        cycle_dump_path = runPath + "cycle_" + str(cycleStep)
        fluid_input_path = os.path.join(cycle_dump_path + "/fluid_input/")
        fluid_output_path = os.path.join(cycle_dump_path + "/fluid_output/")
        kinetic_input_path = os.path.join(cycle_dump_path + "/kinetic_input/")
        kinetic_output_path = os.path.join(cycle_dump_path + "/kinetic_output/")

        if not os.path.exists(cycle_dump_path):
            os.makedirs(cycle_dump_path)
            os.makedirs(fluid_input_path)
            os.makedirs(fluid_output_path)
            os.makedirs(kinetic_input_path)
            os.makedirs(kinetic_output_path)
        if os.path.exists(cycle_dump_path):
            if not os.path.exists(fluid_input_path):
                os.makedirs(kinetic_output_path)
            if not os.path.exists(fluid_output_path):
                os.makedirs(fluid_output_path)
            if not os.path.exists(kinetic_input_path):
                os.makedirs(kinetic_output_path)
            if not os.path.exists(kinetic_output_path):
                os.makedirs(kinetic_output_path)

        if cycleStep > 0:
            previous_cycle_dump_path = runPath + \
                "cycle_" + str(cycleStep - 1) + "/"
            previous_fluid_input_path = os.path.join(previous_cycle_dump_path + "/fluid_input/")
            previous_fluid_output_path = os.path.join(previous_cycle_dump_path + "/fluid_output/")
            previous_kinetic_output_path = os.path.join(previous_cycle_dump_path + "/kinetic_output/")
            previous_kinetic_input_path = os.path.join(previous_cycle_dump_path + "/kinetic_input/")

            moveIMPACTFILE(runPath, cycle_dump_path,
                        previous_kinetic_input_path, previous_kinetic_output_path)

            return(cycle_dump_path, fluid_input_path, fluid_output_path, kinetic_input_path, kinetic_output_path,
                previous_fluid_input_path, previous_fluid_output_path, previous_kinetic_input_path, previous_kinetic_output_path)
        else:
            return(cycle_dump_path, fluid_input_path, fluid_output_path, kinetic_input_path, kinetic_output_path,
                0, 0, 0, 0)
