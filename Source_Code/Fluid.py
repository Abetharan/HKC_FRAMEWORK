import subprocess
import sys

class Fluid:
    def Execute(self, cmd):
        """  
        Purpose: Launches Impact and Sets the number of cores

        Args:
            parameterPath = path where the fluid parameter file is located
            fluidSrcDir = path to fluid exe
        """
        try:
            subprocess.run(cmd, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            import sys
            print(e.output)
            sys.exit(1)


