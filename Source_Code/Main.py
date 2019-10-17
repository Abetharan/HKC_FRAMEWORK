import os
import shutil
import Coupling as cpl
import TmpFileCreator as tfc
import COUPLING_ELH1 as elh1
import COUPLING_IMPACT as impact
import IO as io
import Templating as temple
import Material as material


class main:
    """
    main class i.e. where all objects are passed into and the coupling is run.
    """
    def __init__(self, kinetic_nx, kinetic_ny, kinetic_nx, kinetic_np, fluid_nx, cycles)
    
