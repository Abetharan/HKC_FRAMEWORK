import atexit
import io
import logging
import os
import subprocess
import sys
class Fluid():
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
    
    def clean_up(self):
        try:
            for handler in self.logger.handlers:
                handler.close()
                self.logger.removeHandler(handler)
            # logging.shutdown()
        except OSError:
            pass #ignore the error.  The OSError doesn't seem to be documented(?)
                #as such, it *might* be better to process.poll() and check for 
                #`None` (meaning the process is still running), but that 
                #introduces a race condition.  I'm not sure which is better,
                #hopefully someone that knows more about this than I do can 
                #comment.
        self.__pid = self.__process.pid
    def Execute(self, cmd, path):
        """  
        Purpose: Launches Rad-Hydro code and Sets the number of cores

        Args:
            runPath = Path to exe looks for reference files
        """
       
        filename = path + '/fluid.log'
        with io.open(filename, 'wb') as writer:
            atexit.register(self.clean_up)
            self.__process = subprocess.Popen(cmd, stdout=writer, stderr = subprocess.PIPE)
            _,err = self.__process.communicate()

            if err:
                self.logger.warning("Fluid code failed see log")
                if not os.path.exists("HyKiCT"):
                    self.logger.debug("Fluid NO LONGER EXISTS")
                    sys.exit(1)
       