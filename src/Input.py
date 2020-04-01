import os 
import string
import yaml
class Input:

    def __init__(self, path):
        with open(path, 'r') as stream:
            k = yaml.safe_load(stream)

        self.CX1 = k['CX1']
        self.BASE_DIR = k['BASE_DIR']
        self.IMPACT_SRC_DIR = k['IMPACT_SRC_DIR']
        self.SOL_KIT_SRC_DIR = k['SOL_KIT_SRC_DIR']
        self.F_INIT_PATH = k['F_INIT_PATH']
        self.F_SRC_DIR = k['F_SRC_DIR']
        self.RUN_NAME = k['RUN_NAME']
        self.RUN_PATH = os.path.join(self.BASE_DIR, self.RUN_NAME)
        self.CYCLES = k['CYCLES']
        self.KINETIC_CODE = k['KINETIC_CODE']
        self.K_NX = k['K_NX']
        self.K_NY = k['K_NY']
        self.K_NV = k['K_NV']
        self.K_NP = k['K_NP']
        self.K_DT = k['K_DT']
        self.K_T_MAX = k['K_T_MAX']
        self.K_X_MAX = k['K_X_MAX']
        self.K_V_MAX = k['K_V_MAX']
        self.K_DX = k['K_DX']
        self.K_DV = k['K_DV']
        self.K_DV_MULTI = k['K_DV_MULTI']
        self.K_NT = k['K_NT']
        self.K_PRE_STEP_NT = k['K_PRE_STEP_NT']
        self.K_PRE_STEP_DT = k['K_PRE_STEP_DT']
        self.K_L_MAX = k['K_L_MAX']
        self.K_BC = k['K_BC']
        self.K_MAINTAIN = k['K_MAINTAIN']
        self.K_OUTPUT_FREQ = k['K_OUTPUT_FREQ']
        self.TE = float(k['TE'])
        self.NE = float(k['NE'])
        self.Z = float(k['Z'])
        self.AR = float(k['AR'])
        self.F_NX = k['F_NX']
        self.F_CQ= k['F_CQ']
        self.F_GAMMA = k['F_GAMMA']
        self.F_CFL = k['F_CFL']
        self.F_LASER_WAVELENGTH = k['F_LASER_WAVELENGTH']
        self.F_LASER_POWER = k['F_LASER_POWER']
        self.F_DUR_OF_LASER = k['F_DUR_OF_LASER']
        self.F_LASER_LOC = k['F_LASER_LOC']
        self.F_STEPS = k['F_STEPS']
        self.F_FLUID_T_MAX = k['F_FLUID_T_MAX']
        self.F_INITIAL_DT = k['F_INITIAL_DT']
        self.F_DT_GLOBAL_MAX = k['F_DT_GLOBAL_MAX']
        self.F_DT_GLOBAL_MIN = k['F_DT_GLOBAL_MIN']
        self.F_FEOS_PATH_1 = k['F_FEOS_PATH_1']
        self.F_FEOS_PATH_2 = k['F_FEOS_PATH_2']
        self.F_OUTPUT_FREQ = k['F_OUTPUT_FREQ']
        self.F_BOUNDARY_CONDITION = k['F_BOUNDARY_CONDITION']
        self.F_INITIALISE_START_FILE_RUN = k['F_INITIALISE_START_FILE_RUN']
        self.CONTINUE = k['CONTINUE']
        self.ZIP = k['ZIP']
        self.COUPLEDIVQ = k['COUPLEDIVQ']
        self.COUPLEMULTI = k['COUPLEMULTI']