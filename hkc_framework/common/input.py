import os 
import string
import yaml
class Input:

    def __init__(self, path):
        with open(path, 'r') as stream:
            self.yaml_file = yaml.safe_load(stream)

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
        self.CONTINUE = k['CONTINUE']
        self.ZIP = k['ZIP']
        self.COUPLEDIVQ = k['COUPLEDIVQ']
        self.COUPLEMULTI = k['COUPLEMULTI']