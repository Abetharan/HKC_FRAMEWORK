import sys
import pprint
import numpy as np 
import os
import fnmatch
import string
import subprocess
import shutil 
from scipy import constants
import TmpFileCreator as tfc
import SetFort10Param as sf10p
import FortGenerator as fg
from Templating import Templating
from Kinetic import Kinetic
from IO import IO
import impact_norms_py3 as norms
import impact_module_py3 as cf
kb = constants.value("Boltzmann constant")
me = constants.value("electron mass")
mp = constants.value("proton mass")
e = constants.value("elementary charge")


class IMPACT(Kinetic):
        
    def __init__(self, IO, np_, nv_, nx_, ny_, dt_, t_max_, x_max_, v_max_, boundary_condition_, restart_freq_,CX1_, 
                    fort12Output=["0.2d0", "0.8d0", "1.0d0", "2.0d0", "3.0d0", "5.0d0", "10.0d0", "15.0d0"], unit_check = False):
        
        self._Kinetic = Kinetic()
        self._templater = Templating() 
        self._kinetic_io_obj = IO
        self._run_name = IO._RUN_NAME
        self._run_path = IO._RUN_PATH
        self._base_dir = IO._BASE_DIR
        self._src_dir = IO._K_SRC_DIR
        self._kinetic_output_path = IO.kinetic_output_path
        self._kinetic_input_path = IO.kinetic_input_path
        self._fluid_output_path = IO.fluid_output_path
        self._cycle_path = IO.cycle_dump_path
        self._nv = nv_
        self._nx = nx_
        self._ny = ny_
        self._dt = dt_
        self._t_max = t_max_
        self._x_max = x_max_
        self._v_max = v_max_
        self._np = np_
        self._CX1 = CX1_
        self.boundary_condition = boundary_condition_
        self._restart_freq = restart_freq_
        self.normalised_values = None
        self._fort12Output = fort12Output
        self.makeTmpFiles()
        self._run_unit_test = unit_check
        
    def IMPACTRun(self):
        os.chdir(self._run_path)
        
        if self._CX1:
            cmd = ["mpiexec", "./fp2df1_fast"]    
        else:
            cmd = ["mpirun", "-np", str(self._np), "./fp2df1_fast"]    
        
        super().Execute(cmd, self._cycle_path)

    def test_units(self):
        self._run_unit_test = True
        self.InitIMPACTFromHydro()


    def normalisation(self, calc = False , **kwargs):

        if calc:
            #1e6 dicvision as it is assumed that ne is given in m^-3
            self.normalised_values = norms.impact_inputs(kwargs.get("ne")/1e6, kwargs.get("Te"), 
                            kwargs.get("Z"),kwargs.get("Ar"), kwargs.get("Bz"))  
        else:
            return self.normalised_values


    def makeTmpFiles(self):
        tfc.impactControlDummy(self._run_path)
        tfc.impactFort10(self._run_path)
        tfc.impactHeating(self._run_path)
        tfc.impactProf(self._run_path)
        tfc.impactUserCustom(self._run_path)
        tfc.impactOutputformat(self._cycle_path)

    def createFort14(self):
    
        fort14 = """$user_gra_sel
        op_save_on(20,1) = 1
        op_save_on(20,3) = 1
        op_save_on(21,3) = 1
        op_save_on(22,3) = 1
        op_save_on(23,3) = 1
        op_save_on(24,3) = 1
        op_save_on(33,3) = 1
        $end 
        """
        
        file14 = open(os.path.join(self._run_path, 'fort.14'), 'w')
        file14.write(fort14)

    def createFort12(self):
        time_array = [len(self._fort12Output)] + self._fort12Output
        str_time_array = [str(i) for i in time_array]
        str_times = """ {} """.format("\n".join(str_time_array))
        file12 = open(os.path.join(self._run_path, 'fort.12'), 'w')
        file12.write(str_times)

    def setIMPACTParam(self):
        """
        Purpose: Handles the IO side for launching IMPACT. Involves creation of fort files and moving of files and renaming. Furthermore, handles moving of files 
        when IMPACT completes.
        Args:
            self.normalised_values = Dict with all reference plasma values.
            _KINETIC_nv = Number of grid points in velocity space
            _KINETIC_nx = Number of grid points in X space
            _KINETIC_ny = Number of grid points in Y space
            dt = step size
            kineticTMax = max time 
            fort12Output = fort12 output times in a str list defaults to = ["1.0d0","2.0d0","3.0d0","5.0d0"]
            cycle_dump_path = path to dump all files after IMPACT launches
            self._run_path = Path to where IMPACT looks for all its necessary files
        """
        # Generate fort files
        # On the fly fort10 changes as requierd of fort 10 files.
        wpe_over_nu_ei = self.normalised_values["wpe_over_nu_ei"]
        c_over_vte = self.normalised_values["c_over_vte"]
        Z = self.normalised_values["Z"]
        A = self.normalised_values["Ar"]
        fort10Param = sf10p.set_fort_10(wpe_over_nuei=wpe_over_nu_ei, c_over_vte=c_over_vte,x_bc_type = self.boundary_condition,
                                            atomic_Z=Z, atomic_A=A, nv= self._nv,nx =self._nx, ny=self._ny, dt=self._dt, tmax=self._t_max,
                                            xmax = self._x_max,vmax=self._v_max, do_user_prof_sub=".true.", op_restart_freq = self._restart_freq)

        self._templater.templating(tmpfilePath= self._run_path + '/tmpfort.10',
                writePath=self._run_path, fileName="fort.10", parameters=fort10Param)
        self.createFort12()
        self.createFort14()

    def setCustomFuncs(self):

        shutil.copyfile(os.path.join(self._run_path, "heating.f"),
                        os.path.join(self._run_path, self._run_name + "_heating.f"))
        shutil.copyfile(os.path.join(self._run_path, "control.dat.dummy"),
                        os.path.join(self._run_path, "fp2df1_control.dat.dummy"))
        
        custom_param = {'PATH': "\'" + self._run_path +
                        "\'", 'RUNNAME': "\'/" + self._run_name + "\'"}
        
        self._templater.templating(tmpfilePath= os.path.join(self._run_path, 'user_custom.f'), writePath=self._run_path,
                fileName=self._run_name + "_user_custom.f", parameters=custom_param)
        
        self._templater.templating(tmpfilePath=os.path.join(self._run_path, 'prof.f'), writePath=self._run_path,
                fileName=self._run_name + "_prof.f", parameters=custom_param)
        # Start Coupling sequence
        # os.system('./fp2df1_compile_run.sh')
        os.remove(os.path.join(self._run_path, "heating.f"))
        os.remove(os.path.join(self._run_path, "control.dat.dummy"))
        os.remove(os.path.join(self._run_path, 'user_custom.f'))
        os.remove(os.path.join(self._run_path, 'prof.f'))
        os.remove(os.path.join(self._run_path, "tmpfort.10"))

    def IMPACTCompile(self):
        """  
        Purpose: Handles movement of neccessary template files for custom operators for IMPACT and compiles the code

        Args:
            self._run_path = Path where IMPACT looks for reference files 
        """
        cmd = ["tcsh",'-c',"source "  + self._src_dir + "/"+"fp2df1_compile_run.sh"]    
        #cmd = ["bash",'-c',"source "  + os.getcwd() + "/" + "hydra_kappa.sh"]    
        
        subprocess.Popen(cmd, stdout = sys.stdout, stderr = sys.stderr).communicate()
        # for line in proc.stdout:
        #    (key, _, value) = line.partition("=")
        #    os.environ[key] = value
        #    proc.communicate()

        # pprint.pprint(dict(os.environ))
#        subprocess.run(cmd, stderr=subprocess.PIPE)
        # if(os.path.basename(os.path.normpath(os.getcwd())) == "Source_Code"):
        #     os.system('./hydra_kappa.sh')
        # else:                
        #     os.chdir(os.getcwd() + "/Source_Code")
        #     os.system('./hydra_kappa.sh')


    def setEnvVar(self):

        os.environ["RUN"]= self._run_name
        print(os.environ["RUN"])
        ## laptop
        os.environ["SRCDIR"]=self._src_dir
        #Work station
        #os.environ["SRCDIR"]=os.environ["HOME"] +"/IMPACT/src"
        print(os.environ["SRCDIR"])

        if self._base_dir != None:
            os.environ["BASEDIR"]= self._base_dir
        else:    
            os.environ["BASEDIR"]= os.getcwd()
        
        print(os.environ["BASEDIR"])


        #-----  Code generation control  -----
        os.environ["FPCODE"]     = "fp2df1_fast"
        os.environ["FPPROF"]   = ""
        os.environ["RUN_ARGS"] = "-log_summary -ksp_monitor -pc_type asm -ksp_type bcgs -ksp_rtol 1e-20 -draw_pause -1 -optionstable"

        #-----  Specify behaviour of this run script -os.environ["IS_RESTART
        # os.environ["OVERWRITE_RUN"]
        # os.environ["DONT_COMIPLE"]
        os.environ["DONT_RUN"] = ""
        # os.environ["RUN_IN_QUEUE"]  

        #------  os.environ["COMPILE-TIME ARRAY DIMENSIONS   ------
        os.environ["SYNCHRO_DIMS"] = ""
        os.environ["SMT_PK_SPM_ON "] = "1"
        os.environ["SMT_PC_COLS_ON"] = "0"

        os.environ["NP"]  =  str(self._np)
        os.environ["NVM"]= str(self._nv)

        os.environ["NXM"] =  str(self._nx)
        os.environ["NYM"] =  str(self._ny)
        # ...  Text output behaviour  ...
        os.environ["TEXTOUTPUT_TO_STDOUT_ROOT "]    = "1"
   
    def moveIMPACTFile(self):
        """ 
        Purpose: Moves IMPACT files and initial parameters files to correct paths.
        Args:
            runPath = Run path i.e. where all data is located and where fp2df1 is created. Path is BASEDIR/runName
            cycleDumpPath = cycle path. Path is BASEDIR/runName/cycleNAME
            previsouKineticInputPath = previous cycle kinetic input folder which is located as followeing runPath/previsousCyclePath/kinetic_input
            previousKineticOutPut =  previous cycle kinetic output folder which is located as followeing runPath/previsousCyclePath/kinetic_output
        """

        print("#"*100)
        print('\033[1m' + '#'*50 + ' MOVE IMPACT FILES ' + '#'*50 + '\033[0m')
        filenames = ['ionden', 'rad_to_electron',
                    'xf', 'eden', 'laserdep', 'tmat', 'zstar']

        for name in filenames:
            shutil.move(self._run_path + "/" + self._run_name + "_" + name + ".xy",
            self._kinetic_input_path + "/" + self._run_name + "_" + name + ".xy", )

        for file in os.listdir(self._run_path):
            _, extension = os.path.splitext(file)
            if extension == ".xy" or extension == ".xyz" or extension == ".xyv" or extension == ".xyt" or extension == ".dat" or extension == ".t":
                shutil.move(file, self._kinetic_output_path)
        
        import glob
        for file in glob.glob(str(self._run_path) + '/*.out'):
            print(file)
            shutil.move(file, self._cycle_path)    
   
    def IMPACTInitHydroNextFiles(self):
        
        # Convert to SI from Impact Norms
        var = "qxX"
        timeStep = self.findLastIndex(self._run_path, var, "kinetic")
        runName = "default"
        #runName = os.environ['RUN']
        if timeStep < 10:
            time_step = "0" + str(timeStep)

        varArrays = cf.load_dict(self._run_path,
                                runName, var, str(time_step), iter_number=None)
        varList = varArrays['mat'][:]
        #outputVar = "electron_heat_flow"
        normConst = me * \
            pow(self.normalised_values['vte'], 3) * self.normalised_values['ne'] * 1e21 * 1e6
        qe = varList * normConst
        import glob
        LargestIndex = 0
        
        for file in os.listdir(self._kinetic_io_obj.fluid_output_path):
           if fnmatch.fnmatch(file, "Coord_*"):
               k = os.path.splitext(file)[0]
               k = k.split("_")
               timeIndex = k[-1]
               if int(timeIndex) > LargestIndex:
                   LargestIndex = int(timeIndex)
        
        ###
        def Heatflow(electron_thermal_flux, mass):
            nx = len(mass)
            HeatConductionE = np.zeros(nx)
            for i in range(nx):   
                HeatConductionE[i] = -(electron_thermal_flux[i + 1] - electron_thermal_flux[i]) / mass[i]
            return(HeatConductionE)

        mass = np.loadtxt(self._kinetic_io_obj.fluid_input_path + "/mass.txt")        
        electronheatflow= Heatflow(qe, mass)
        
        np.savetxt(os.path.join(self._kinetic_io_obj.next_fluid_input_path, "qe.txt"), electronheatflow)    

        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "CELL_WALL_X/CELL_WALL_X_" + str(LargestIndex) +".txt")    
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"coord.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "VELOCITY/VELOCITY_" + str(LargestIndex) +".txt") 
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"velocity.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "DENSITY/DENSITY_" + str(LargestIndex) +".txt")  
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"density.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_" + str(LargestIndex) +".txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"electron_temperature.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "ION_TEMPERATURE/ION_TEMPERATURE_" + str(LargestIndex) +".txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"ion_temperature.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "MASS/MASS_" + str(LargestIndex) + ".txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path,"mass.txt"))
        shutil.copyfile(os.path.join(self._kinetic_io_obj.fluid_output_path + "ZBAR/ZBAR_" + str(LargestIndex) +".txt")
                                        ,os.path.join(self._kinetic_io_obj.next_fluid_input_path, "Z.txt"))
        
    def findLastIndex(self, path, var, which):
        largest_index = 0
        for file in os.listdir(path):
            if which == "kinetic":
                match = "*_" + var + "_*"
            else:
                match = var + "_*"

            if fnmatch.fnmatch(file, match):
                k = os.path.splitext(file)[0]
                k = k.split("_")
                time_index = k[-1]
                if str.isalpha(time_index):
                    continue
                if int(time_index) > largest_index:
                    largest_index = int(time_index)

        return(largest_index)


    def ImpactFileWriteFormat(self, tmp_file_location_, open_file_, coord_, var_data_, var_, unit_tester = None):
        
        from string import Template
        if var_ != "coord":
            nx = len(coord_) - 1  # becuse index start from 0
            cell_centered_coord = [(coord_[i + 1] + coord_[i]) / 2 for i in range(nx)]
            
            if unit_tester is not None:
                for i in range(nx):
                    error = (cell_centered_coord[i] - unit_tester[i]) / unit_tester[i]
                    if error < 1e-7:
                        import sys
                        print('cell centered coords are different')
                        sys.exit(1)

            last_dx = cell_centered_coord[nx - 1] - cell_centered_coord[nx - 2]
            last_ghost_cell_coord = cell_centered_coord[nx - 1] + last_dx
            cell_centered_coordwith_ghost = [-cell_centered_coord[0]
                                        ] + cell_centered_coord + [last_ghost_cell_coord]
            leadingdim = 3
            var_data_with_ghost = [var_data_[0]] + var_data_.tolist() + [var_data_[-1]]
            maxparam = "-6.000000000000000e+03 0.000000000000000e+00 6.000000000000000e+03"
            length = len(var_data_with_ghost)
            template_swap_dict = {'leadingdim': leadingdim, 'maxparam': maxparam, 'len': length, 'xlist': ' '.join(str(coo) for coo in cell_centered_coordwith_ghost),
                                'arraylist': '\n'.join((str(var) + "\t" + str(var) + "\t" + str(var)) for var in var_data_with_ghost)}

        else:
            leadingdim = 1
            maxparam = "0.000000000000000e+00"
            length = len(coord_)
            template_swap_dict = {'leadingdim': leadingdim, 'maxparam': maxparam, 'len': length, 'xlist': ' '.join(str(coo) for coo in coord_),
                                'arraylist': '\n'.join(str(var) for var in var_data_)}

        tmp_file_in = open(tmp_file_location_)
        src = string.Template(tmp_file_in.read())
        result = src.substitute(template_swap_dict)
        open_file_.write(result)
        open_file_.close()

        return(0)


    def InitIMPACTFromHydro(self):

        lastIndex = self.findLastIndex(os.path.join(self._fluid_output_path, "CELL_WALL_X"), "CELL_WALL_X", "fluid")
        f_x_grid = np.loadtxt(self._fluid_output_path + "CELL_WALL_X/CELL_WALL_X_" + str(lastIndex) + ".txt")
        f_x_centered_grid = np.loadtxt(self._fluid_output_path + "CELL_CENTRE_X/CELL_CENTRE_X_" + str(lastIndex) + ".txt")
        f_v = np.loadtxt(self._fluid_output_path + "VELOCITY/VELOCITY_" + str(lastIndex) + ".txt")
        f_ne = np.loadtxt(
                self._fluid_output_path + "ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_" + str(lastIndex) + ".txt")
        f_ni = np.loadtxt(
                self._fluid_output_path + "ION_NUMBER_DENSITY/ION_NUMBER_DENSITY_" + str(lastIndex) + ".txt")
        f_Te = np.loadtxt(
             self._fluid_output_path + "ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_" + str(lastIndex) + ".txt")
        f_laser = np.loadtxt(
            self._fluid_output_path + "INVERSE_BREM/INVERSE_BREM_" + str(lastIndex) + ".txt")
        f_brem = np.loadtxt(self._fluid_output_path + "BREM/BREM_" + str(lastIndex) + ".txt")
        f_density = np.loadtxt(
            self._fluid_output_path + "DENSITY/DENSITY_" + str(lastIndex) + ".txt")
        f_Z = np.loadtxt(self._fluid_output_path + "ZBAR/ZBAR_" + str(lastIndex) + ".txt")

        # Energy normlisation
        Un = me * self.normalised_values['vte']**2 * \
            self.normalised_values['ne'] * 1e21 * 1e6
        power_density = (Un / (self.normalised_values['tau_ei']))
        # NOrmalise SI to Impact norms
        IMPACT_x_grid = f_x_grid / self.normalised_values["lambda_mfp"]
        IMPACT_x_centered_grid = f_x_centered_grid / self.normalised_values["lambda_mfp"]
        # 1e21 normalisation factor. 1e6 there to convert from m^-3 to cm^-3
        IMPACT_ne = f_ne / (1e6 * self.normalised_values['ne'] * 1e21)
        # 1e6 there to convert from m^-3 to cm^-3
        IMPACT_ni = f_ni / (1e6 * self.normalised_values['ni'] * 1e21)
        IMPACT_Te = (f_Te * (kb/e)) / self.normalised_values['Te']
        IMPACT_laser = (f_laser * f_density) / power_density
            
        IMPACT_brem = (f_brem * f_density) / power_density
            
        IMPACT_Z = f_Z / self.normalised_values['Z']

        if self._run_unit_test:
            for c in range(len(IMPACT_x_centered_grid)):
                tolerance = 1e-8
                x_grid_diff = abs(IMPACT_x_grid[c] * self.normalised_values["lambda_mfp"] - f_x_grid[c])/f_x_grid[c]
                ne_diff = abs(IMPACT_ne[c] * self.normalised_values['ne'] *1e21 *1e6 - f_ne[c]) / f_ne[c]
                ni_diff = abs(IMPACT_ni[c] * self.normalised_values['ni'] * 1e21 * 1e6- f_ni[c]) / f_ni[c]
                Te_diff = abs(IMPACT_Te[c] * self.normalised_values['Te'] - f_Te[c]) / f_Te[c]
                laser_diff = abs(IMPACT_laser[c] * (power_density/f_density[c]) - f_laser[c]) / f_laser[c]
                brem_diff = abs(IMPACT_brem[c] * (power_density/f_density[c])- f_brem[c]) / f_brem[c]
                Z_diff = abs(IMPACT_Z[c] * self.normalised_values['Z'] - f_Z[c]) / f_Z[c]
                diff_array = [x_grid_diff, ne_diff, ni_diff, Te_diff, laser_diff , brem_diff, Z_diff]
                
                if any(diff_array < tolerance):
                    import sys
                    print('Differences in IMPACT and fluid values after normalising')
                    sys.exit(1)


        # if(np.max(ne_norm) / np.min(ne_norm)) > 30 or (np.max(Te_norm) / np.min(Te_norm)) > 5:
        #     kinetic_x = CoordGenerator(x_norm, fluid_v, int(
        #         os.environ["NXM"]) + 1, numberSmoothingCellsRatio=.5)
        # else:
        #     kinetic_x = np.linspace(
        #         x_norm[0], x_norm[-1], int(os.environ["NXM"]) + 1)

        # fluid_centered_x = [(x_norm[i + 1] + x_norm[i])/2 for i in range(fluidNx)]
        # kinetic_centered_x = [(kinetic_x[i + 1] + kinetic_x[i]) /
        #                     2 for i in range(int(os.environ["NXM"]))]

        #if(interpolator == "cubic"):
            #cs_ne = CubicSpline(fluid_centered_x, ne_norm)
            #cs_ni = CubicSpline(fluid_centered_x, ni_norm)
            #cs_Te = CubicSpline(fluid_centered_x, Te_norm)
            #cs_laser = CubicSpline(fluid_centered_x, laser_norm)
            #cs_brem = CubicSpline(fluid_centered_x, brem_norm)
            #cs_Z = CubicSpline(fluid_centered_x, Z_norm)
        #else:
            #cs_ne = interpolate.interp1d(
                #fluid_centered_x, ne_norm, fill_value="extrapolate")
            #cs_ni = interpolate.interp1d(
                #fluid_centered_x, ni_norm, fill_value="extrapolate")
            #cs_Te = interpolate.interp1d(
                #fluid_centered_x, Te_norm, fill_value="extrapolate")
            #cs_laser = interpolate.interp1d(
                #fluid_centered_x, laser_norm, fill_value="extrapolate")
            #cs_brem = interpolate.interp1d(
                #fluid_centered_x, brem_norm, fill_value="extrapolate")
            #cs_Z = interpolate.interp1d(
                #fluid_centered_x, Z_norm, fill_value="extrapolate")

        #kinetic_ne = cs_ne(kinetic_centered_x)
        #kinetic_ni = cs_ni(kinetic_centered_x)
        #kinetic_Te = cs_Te(kinetic_centered_x)
        #kinetic_laser = cs_laser(kinetic_centered_x)
        #kinetic_brem = cs_brem(kinetic_centered_x)
        #kinetic_Z = cs_Z(kinetic_centered_x)
        # import matplotlib.pyplot as plt
        # plt.plot(kinetic_Te)
        # plt.show()

        if not os.path.exists(self._run_path + "/tmpWrite.txt"):
            tfc.impactOutputformat(self._run_path)

        impactNeFile = open(self._run_path + "/" + os.environ["RUN"] + "_eden.xy", "w")
        impactNiFile = open(
            self._run_path + "/" + os.environ["RUN"] + "_ionden.xy", "w")
        impactTeFile = open(self._run_path + "/" + os.environ["RUN"] + "_tmat.xy", "w")
        impactXFile = open(self._run_path + "/" + os.environ["RUN"] + "_xf.xy", "w")
        impactZfile = open(self._run_path + "/" + os.environ["RUN"] + "_zstar.xy", "w")
        impactLaserFile = open(
            self._run_path + "/" + os.environ["RUN"] + "_laserdep.xy", "w")
        impactRadFile = open(
            self._run_path + "/" + os.environ["RUN"] + "_rad_to_electron.xy", "w")
      

        if self._run_unit_test:
            unit_tester = f_x_centered_grid
        else:
            unit_tester = None
        self.ImpactFileWriteFormat(self._run_path + "/tmpWrite.txt",
                        impactNeFile, IMPACT_x_grid, IMPACT_ne, "ne", unit_tester)
        self.ImpactFileWriteFormat(self._run_path + "/tmpWrite.txt",
                        impactNiFile, IMPACT_x_grid, IMPACT_ni, "ni", unit_tester)
        self.ImpactFileWriteFormat(self._run_path + "/tmpWrite.txt", impactXFile,
                        IMPACT_x_grid, IMPACT_x_grid, "coord", unit_tester)
        self.ImpactFileWriteFormat(self._run_path + "/tmpWrite.txt",
                        impactTeFile, IMPACT_x_grid, IMPACT_Te, "Te", unit_tester)
        self.ImpactFileWriteFormat(self._run_path + "/tmpWrite.txt", impactLaserFile,
                        IMPACT_x_grid, IMPACT_laser, "laser", unit_tester)
        self.ImpactFileWriteFormat(self._run_path + "/tmpWrite.txt", impactRadFile,
                        IMPACT_x_grid, IMPACT_brem, "Rad", unit_tester)
        self.ImpactFileWriteFormat(self._run_path + "/tmpWrite.txt", impactZfile,
                       IMPACT_x_grid, IMPACT_Z, "zstar", unit_tester)
       

