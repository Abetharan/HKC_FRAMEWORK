import subprocess
import sys
import io
import os
class Fluid():

    def Execute(self, cmd, path):
        """  
        Purpose: Launches Rad-Hydro code and Sets the number of cores

        Args:
            runPath = Path to exe looks for reference files
        """
       
        filename = path + '/f_test.log'
        with io.open(filename, 'wb') as writer, io.open(filename, 'rb', 1) as reader:
            process = subprocess.Popen(cmd, stdout=writer)
            while process.poll() is None:
                sys.stdout.write(reader.read().decode('utf-8'))
            # Read the remaining
            sys.stdout.write(reader.read().decode('utf-8'))
       