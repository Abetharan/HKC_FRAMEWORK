import subprocess
import sys
import io
import os
import threading
import time
import numpy as np 
class Kinetic():
    
        
    def Execute(self, cmd, path, output_folder_path, convergence_func, nx):
        """  
        Purpose: Launches Impact and Sets the number of cores

        Args:
            runPath = Path where IMPACT looks for reference files
            _KINETIC_np = Number of cores being used 
        """
       
        def ConvergenceMonitoring(proc, output_folder, convergence_func, nx):
            #Path to heat flow
            convergance = 10
            file_counter = 0
            time.sleep(10)
            multipliers = np.zeros(nx + 1) 
            while True:

                heat_flows = os.listdir(output_folder)
                new_file_counter = len(heat_flows)
                #Last Heat 
                if new_file_counter == 0:
                    latest_heat_flow_path = ''
                else:
                    latest_heat_flow_path = os.path.join(output_folder, heat_flows[-1])
                
                if os.path.exists(latest_heat_flow_path)  and new_file_counter != file_counter:
                    
                    if os.access(latest_heat_flow_path, os.R_OK):
                        #routine
                        curr_multipliers = convergence_func(latest_heat_flow_path)                    
                    else:
                        time.sleep(.2)
                        curr_multipliers = convergence_func(latest_heat_flow_path)                    
                        #routine
                    
                    multipliers = np.vstack((multipliers, curr_multipliers))

                if len(np.shape(multipliers)) >= 2:
                    if np.shape(multipliers)[0] >= 3:
                        multipliers = np.delete(multipliers, 0, 0)
                    if np.shape(multipliers)[0] >= 2:
                        convergance = abs(np.array(multipliers[0, :]) - np.array(multipliers[1, : ]))

                    convergance[np.isnan(convergance)] = 0
                    convergance[np.isinf(convergance)] = 0

                file_counter = new_file_counter

                if np.nanmax(convergance) < 1e-3 and np.nanmax(convergance) != 0:
                    proc.terminate()
                    break
        
        filename = path + '/k_test.log'
        with io.open(filename, 'wb') as writer, io.open(filename, 'rb', 1) as reader:
            process = subprocess.Popen(cmd, stdout=writer)
            monitor = threading.Thread(target=ConvergenceMonitoring, args = (process, output_folder_path, convergence_func, nx))
            monitor.start()

            while process.poll() is None:
                sys.stdout.write(reader.read().decode('utf-8'))

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
    
