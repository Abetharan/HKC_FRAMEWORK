import atexit
import io
import logging
import os
import signal
import subprocess
import sys
import threading
import time
import numpy as np 

class Kinetic():
    def __init__(self, cmd, convergence_monitoring = False, convergence_func = None,
                    thread_log_path = None):
        self.log_path = thread_log_path
        self.cycle = 0
        self.cycle_dump_path = ""
        self.convergence_func = convergence_func
        self.monitor_convergence = convergence_monitoring
        self.nx = 0
        self.cmd = cmd
        self.converged = False

    def convergenceMonitoring(self, kinetic_heat_flow_output_folder_path,  stop_event):
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
        logging.basicConfig(filename=os.path.join(self.log_path, "Monitoring_debug.log"),
                    filemode='a',level=logging.DEBUG,
                    format='(%(threadName)-10s) %(asctime)s  %(message)s',
                    )
        logging.info('Starting Cycle :')
        logging.info(self.cycle)
        #Path to heat flow
        convergance = 10
        file_counter = 0
        time.sleep(10)
        multipliers = np.zeros(self.nx + 1)
        #stop_event can be set internally once convergenced has been reached 
        # or externally to kill the daemon if the alloted time in the 
        # subprocess has been reached. 
        while not stop_event.is_set():
            all_heat_flower_files = os.scandir(kinetic_heat_flow_output_folder_path)
            heat_flow_only = [f for f in all_heat_flower_files if not f.name.startswith('.')]
            sorted_heat_flows = sorted(heat_flow_only, 
                                key = lambda x: int(os.path.splitext(x.name)[0].split('_')[-1]))

            new_file_counter = len(sorted_heat_flows)
            #Last Heat 
            if(os.path.exists(sorted_heat_flows[-1].path)
                and new_file_counter != file_counter):
                logging.info('Latest path') 
                logging.info(sorted_heat_flows[-1].path)
                if os.access(sorted_heat_flows[-1].path, os.R_OK):
                    curr_multipliers = self.convergence_func(sorted_heat_flows[-1].path)
                    #Possible occurence where files is read safe but nothing has been written to it.
                    if(len(curr_multipliers) <= 0):
                        logging.warning("Read safe.. No Content")
                        continue                   
                    multipliers = np.vstack((multipliers, curr_multipliers))
                else:
                    continue

                if len(np.shape(multipliers)) >= 2:
                    if np.shape(multipliers)[0] >= 2:
                        convergance = abs(np.array(multipliers[0, :]) -
                                        np.array(multipliers[1, : ])) / np.array(multipliers[0,:])
                        convergance[np.isnan(convergance)] = 0
                        convergance[np.isinf(convergance)] = 0
                        multipliers = np.delete(multipliers, 0, 0)

                    logging.info("Convergence: ")
                    logging.info(np.nanmax(convergance))

                if np.nanmax(convergance) < 1e-3 and np.nanmax(convergance) != 0:
                    self.clean_up()
                    self.converged = True
                    time.sleep(0.5)
                    logging.info("Converged ....Exiting")
                    stop_event.set()
                #Update file counter
                file_counter = new_file_counter
                logging.info(file_counter)
        
        #Kill thread on exit
        self.cycle+=1
        logging.warning("Killing thread")

    def clean_up(self):
        try:
            self.__process.terminate()
            logging.shutdown()
        except OSError:
            pass #ignore the error.  The OSError doesn't seem to be documented(?)
                #as such, it *might* be better to process.poll() and check for 
                #`None` (meaning the process is still running), but that 
                #introduces a race condition.  I'm not sure which is better,
                #hopefully someone that knows more about this than I do can 
                #comment.
        self.__pid = self.__process.pid
    
    def Execute(self,  kinetic_heat_flow_output_folder_path):
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

        if self.monitor_convergence:
            stop_event = threading.Event()
            monitor = threading.Thread(name = "Convergence_Monitoring", 
            target=self.convergenceMonitoring, args = 
            (kinetic_heat_flow_output_folder_path, stop_event),
            daemon = True)
            monitor.start()

        filename = self.cycle_dump_path + '/kinetic.log'
        with io.open(filename, 'wb') as writer:
        #run command provided
            atexit.register(self.clean_up)
            self.__process = subprocess.Popen(self.cmd, stdout=writer, stderr = subprocess.PIPE)

        _, err = self.__process.communicate()
        if self.__process.poll() is not None and self.monitor_convergence:
            stop_event.set()

        if err and not self.converged:
            print("Kinetic code failed see log")
            sys.exit(0)