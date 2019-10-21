import subprocess
import sys
class Kinetic():

    def Execute(self, cmd):
        """  
        Purpose: Launches Impact and Sets the number of cores

        Args:
            runPath = Path where IMPACT looks for reference files
            _KINETIC_np = Number of cores being used 
        """
        try:
            subprocess.run(cmd, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            import sys
            print(e.output)
            sys.exit(1)
    
