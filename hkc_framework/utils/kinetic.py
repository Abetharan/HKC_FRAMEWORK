import subprocess
import sys
import io
import os
import threading
import time
import numpy as np 
class Kinetic():
    
        
    def Execute(self, cmd, log_path, monitor_on,  kinetic_heat_flow_output_folder_path, convergence_func, nx):
        """  
        Purpose: Launch the command relevant to kinetic code specified, if maintain 
                 f_0 is on also spawns a daemon to check for convergence, only relevant
                 for SOL-KIT.
        Args:
            cmd : list of commands of format ['', '', '' ]
            log_path : path to output log file 
            monitor_on : if maintaining f_0 only available with SOL-KiT a daemon is spawned to
                        monitor if q/q_sh is converged and kill SOL-KiT prior to its natural end time.
            kinetic_heat_flow_output_folder_path_path : Path to where heat flow is output by SOL-KiT
            convergence_func : function that calculates convergance values i.e. q/q_sh
            nx : number of cell-centrs.
        """
       
        def ConvergenceMonitoring(proc, kinetic_heat_flow_output_folder_path, convergence_func, nx):
            """  
            Purpose: Function which the daemon i.e. thread is constantly doing.. checking for conergence
            Args:
                proc : subprocess that needs to be terminated once convergence has been achieved. 
                log_path : path to output log file 
                monitor_on : if maintaining f_0 only available with SOL-KiT a daemon is spawned to
                            monitor if q/q_sh is converged and kill SOL-KiT prior to its natural end time.
                kinetic_heat_flow_output_folder_path_path : Path to where heat flow is output by SOL-KiT
                convergence_func : function that calculates convergance values i.e. q/q_sh
                nx : number of cell-centrs.
            """
            #Path to heat flow
            convergance = 10
            file_counter = 0
            time.sleep(10)
            multipliers = np.zeros(nx + 1) 

            while True:
                heat_flows = os.listdir(kinetic_heat_flow_output_folder_path)
                new_file_counter = len(heat_flows)
                #Last Heat 
                if new_file_counter == 0:
                    latest_heat_flow_path = ''
                else:
                    latest_heat_flow_path = os.path.join(kinetic_heat_flow_output_folder_path, heat_flows[-1])
                
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
            
        ##Relevant stuff to create logs 
        filename = log_path + '/k_test.log'
        with io.open(filename, 'wb') as writer, io.open(filename, 'rb', 1) as reader:
            #run command provided
            process = subprocess.Popen(cmd, stdout=writer, stderr = subprocess.PIPE)
            _, err = process.communicate()
            ##Spawn a thread which runs convergencemonitoring function. 
            if monitor_on:
                monitor = threading.Thread(target=ConvergenceMonitoring, args = (process, kinetic_heat_flow_output_folder_path, convergence_func, nx), daemon = True)
                monitor.start()

            while process.poll() is None:
                sys.stdout.write(reader.read().decode('utf-8'))

            # Read the remaining
            sys.stdout.write(reader.read().decode('utf-8'))

            if err:
                print("Kinetic code failed see log")
                sys.exit(0)