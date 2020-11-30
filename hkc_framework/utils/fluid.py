import atexit
import io
import logging
import numpy as np 
import os
import pathlib
import signal
import subprocess
import sys
class Fluid():
    def Execute(self, cmd, path):
        """  
        Purpose: Launches Rad-Hydro code and Sets the number of cores

        Args:
            runPath = Path to exe looks for reference files
        """
       
        filename = path + '/fluid.log'
        with io.open(filename, 'wb') as writer:
            self.__process = subprocess.Popen(cmd, stdout=writer, stderr = subprocess.PIPE)
            _,err = self.__process.communicate()

            if err:
                path_obj = pathlib.Path(path)
                root_path = path_obj.parent
                np.savetxt(os.path.join(root_path, "status.txt"), np.array([1], dtype=np.int), fmt = '%1.1i')
                sys.exit(1)
       