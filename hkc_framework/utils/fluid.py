import atexit
import io
import logging
import os
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
                sys.exit(1)
       