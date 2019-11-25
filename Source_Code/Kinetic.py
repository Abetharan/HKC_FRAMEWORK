import subprocess
import sys
import io
import os
class Kinetic():

    def Execute(self, cmd, path):
        """  
        Purpose: Launches Impact and Sets the number of cores

        Args:
            runPath = Path where IMPACT looks for reference files
            _KINETIC_np = Number of cores being used 
        """
       
        filename = path + '/k_test.log'
        with io.open(filename, 'wb') as writer, io.open(filename, 'rb', 1) as reader:
            process = subprocess.Popen(cmd, stdout=writer)
            while process.poll() is None:
                sys.stdout.write(reader.read().decode('utf-8'))
                #time.sleep(0.5)
            # Read the remaining
            sys.stdout.write(reader.read().decode('utf-8'))
       
       
       
       
       # proc = subprocess.Popen(cmd, stdout = subprocess.PIPE)
        #subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=sys.stderr, universal_newlines=True).communicate()
        #for line in proc.stdout:
        #   (key, _, value) = line.partition("=")
        #   os.environ[key] = value
        #   proc.communicate()

#        pprint.pprint(dict(os.environ))
        # try:
        #     subprocess.run(cmd, stderr=subprocess.PIPE)
        # except subprocess.CalledProcessError as e:
        #     import sys
        #     print(e.output)
        #     sys.exit(1)
    
