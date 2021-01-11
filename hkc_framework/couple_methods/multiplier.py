import numpy as np 


class Multiplier: 

    def __init__(self):
        self.search_tolerance = 1e-9

    def preHeatModel(self, rough = False):
        """ 
        Purpose: Model Pre-heat using an exponential fitting parameter, fit parameter
                    is spat out to be used by the fluid code.
        Args:
        Returns:
            B = Fitting paramters
            preheat_start = start index for preheat
            preheat_end = end index of preheat
        """
        if rough:
            all_index = np.where(np.isinf(self.q_vfp_q_sh_multipliers) == True)[0]
            valid_index = [i for i in all_index if i > int(len(self.vfp_heat)/2)]
            if len(valid_index) == 0:
                return self.pre_heat_start_index ,self.pre_heat_last_index , self.pre_heat_fit_params 
            preheat_start = valid_index[0] - 1
            preheat_end = np.where(abs(self.vfp_heat[preheat_start:]) < self.search_tolerance)
        else:
            preheat_start = np.where(self.q_vfp_q_sh_multipliers[~np.isinf(self.q_vfp_q_sh_multipliers)] != 0)[0][-1]
            preheat_end = np.where(abs(self.vfp_heat[preheat_start:]) < self.search_tolerance)
        #if search fails i.e. valid heat flow in all domain
        if(len(preheat_end[0]) == 0):
            preheat_end = len(self.q_vfp_q_sh_multipliers)
        else:
            preheat_end = preheat_end[0][0] + preheat_start

        L = self.cell_wall_coord[preheat_end] -self.cell_wall_coord[preheat_start] 
        B = []

        for i in range(preheat_start, preheat_end):
            b = (-1* (self.cell_wall_coord[i] - self.cell_wall_coord[preheat_start]) 
                / (L * np.log(abs(self.vfp_heat[i]/self.vfp_heat[preheat_start]))))
            B.append(b)
            
        B = np.array(B)
        return(preheat_start, preheat_end, B)

    def frontHeatModel(self, rough=False):
        """
        Purpose: Model heat wave from top of a bath.
        Args:
        Returns:
            B = Fitting paramters
            frontheat_start = start of front heat
            frontheat_end = end of front 
            heat
        """
        if rough:
            all_index = np.where(np.isinf(self.q_vfp_q_sh_multipliers) == True)[0]
            valid_index = [i for i in all_index if i <= int(len(self.vfp_heat)/2)]
            if len(valid_index) == 0:
                return self.front_heat_start_index ,self.front_heat_last_index , self.front_heat_fit_params
            frontheat_start = valid_index[-1] + 1
            frontheat_end = np.where(abs(self.vfp_heat[:frontheat_start]) < self.search_tolerance)
        else:
            frontheat_start = np.where(self.q_vfp_q_sh_multipliers[~np.isinf(self.q_vfp_q_sh_multipliers)] != 0)[0][0]
            frontheat_end = np.where(abs(self.vfp_heat[:frontheat_start]) < self.search_tolerance)
        #if search fails i.e. valid heat flow in all domain
        if(len(frontheat_end[0]) == 0):
            frontheat_end = 0
        else:
            frontheat_end = frontheat_end[0][0]

        L = abs(self.cell_wall_coord[frontheat_end] - self.cell_wall_coord[frontheat_start])
        B = []

        for i in range(0, frontheat_start+1):
            b = (-1* (self.cell_wall_coord[i] - self.cell_wall_coord[frontheat_start]) 
                / (L * np.log(abs(self.vfp_heat[i]/self.vfp_heat[frontheat_start]))))
            B.append(b)

        B = np.array(B)
        return(frontheat_start, frontheat_end, B)

    def _detectAnamalousHeat(self, heat_flow):
        if (all(self.spitzer_harm_heat == 0)):
            return None, None

        start_of_spitzer_harm_heat_flow_index = np.where(abs(heat_flow) > 0)[0][0]
        last_of_spitzer_harm_heat_flow_index = np.where(abs(heat_flow) > 0)[0][-1]
        front = None 
        pre = None
        inf_mulitplier_index = np.where(np.isinf(self.q_vfp_q_sh_multipliers) == True)[0]
        diff_inf_index = np.diff(inf_mulitplier_index)
        if any(diff_inf_index > 1):
            max_diff = np.argmax(diff_inf_index)
            if((inf_mulitplier_index[max_diff] < int(len(self.q_vfp_q_sh_multipliers)/2)) and 
                (inf_mulitplier_index[max_diff + 1] > int(len(self.q_vfp_q_sh_multipliers)/2))):
                pass
            else:
                start_of_spitzer_harm_heat_flow_index = int(len(heat_flow)/2)
                last_of_spitzer_harm_heat_flow_index = int(len(heat_flow)/2)
        if(any(np.isnan(self.q_vfp_q_sh_multipliers[1:-2])) #boundaries are always nan as we have 0 inflow conditions. 
            or any(np.isinf(self.q_vfp_q_sh_multipliers))):
        
            if any(abs(self.vfp_heat[:start_of_spitzer_harm_heat_flow_index]) > self.search_tolerance):
                if(start_of_spitzer_harm_heat_flow_index == int(len(heat_flow)/2)):
                    front = self.frontHeatModel(rough = True)
                else:
                    front = self.frontHeatModel()
            if(any(abs(self.vfp_heat[last_of_spitzer_harm_heat_flow_index:]) > self.search_tolerance)):
                if(start_of_spitzer_harm_heat_flow_index == int(len(heat_flow)/2)):
                    pre = self.preHeatModel(rough = True)
                else:
                    pre = self.preHeatModel()

            self.q_vfp_q_sh_multipliers[np.isnan(self.q_vfp_q_sh_multipliers)] = 0
            self.q_vfp_q_sh_multipliers[np.isinf(self.q_vfp_q_sh_multipliers)] = 0
        return front, pre

    def method(self, sh_heat_flow, vfp_heat_flow, laser_dir = None, **kwargs):
        """ Purpose: Find multipliers and exponentially extrapolate for pre-heat
            Returns: Multipliers
        """

        self.spitzer_harm_heat = sh_heat_flow
        self.vfp_heat = vfp_heat_flow 
        self.cell_wall_coord = kwargs['cell_wall_coord']
        self.q_snb = kwargs['q_snb'] 
        ##Limit density method 
        ##Multipliers -> 1 for region that has been cut off 
        ## If there is heat-flow there going to be using spitzer-harm regardless 
        ## Models currently will work as intended. 
        ##Maybe rework or work around anamalous heat detection function. 
        self.front_heat_start_index = np.int64(0) 
        self.front_heat_last_index = np.int64(0)
        self.pre_heat_start_index = np.int64(0)
        self.pre_heat_last_index = np.int64(0)
        self.front_heat_fit_params = None
        self.pre_heat_fit_params = None

        if self.q_snb is not None:
            heat_flow = self.q_snb 
        else:
            heat_flow = self.spitzer_harm_heat

        limit_index = len(heat_flow)
        if laser_dir is not None:
            len_vfp = len(self.vfp_heat)
            if laser_dir == "right":
                append_index = limit_index - len_vfp# + 1
                heat_flow = heat_flow[append_index:]
            else:
                heat_flow = heat_flow[:len_vfp] 
        #End
        ##Test for pre-heat via looking at NaN outputs expected from q/q_sh
        #if nans
        # if self.snb:
        #     self.q_vfp_q_sh_multipliers = np.array(self.vfp_heat/heat_flow)
        # else:
        self.q_vfp_q_sh_multipliers = np.array(self.vfp_heat/heat_flow)


        if laser_dir is not None:
            if laser_dir == "right":
                self.q_vfp_q_sh_multipliers[-1] = 0
            else:
                self.q_vfp_q_sh_multipliers[0] = 0
            
        #Detect if there is Front heat
        #Detect if there is Pre-Heat
        #Modify as bounds will always be nan. 
        if self.q_snb is not None:
            front, pre = self._detectAnamalousHeat(heat_flow)
            if front is not None:
                (self.front_heat_start_index, 
                self.front_heat_last_index, 
                self.front_heat_fit_params) = front
            if pre is not None:
                (self.pre_heat_start_index,
                self.pre_heat_last_index,
                self.pre_heat_fit_params) = pre

        if laser_dir is not None:
            padding = limit_index - len_vfp
            if laser_dir == "right":
                self.q_vfp_q_sh_multipliers[0] = 1
                self.q_vfp_q_sh_multipliers = np.pad(self.q_vfp_q_sh_multipliers, (padding, 0), 'constant', constant_values = (1,0))
            else:
                self.q_vfp_q_sh_multipliers[-1] = 1
                self.q_vfp_q_sh_multipliers = np.pad(self.q_vfp_q_sh_multipliers, (0, padding), 'constant', constant_values = (0,1))
            
            if self.front_heat_start_index > 0:
                self.front_heat_start_index += padding
                self.front_heat_last_index += padding
            if self.pre_heat_start_index > 0:
                self.pre_heat_last_index += padding
                self.pre_heat_last_index += padding

        self.q_vfp_q_sh_multipliers[0] = 0
        self.q_vfp_q_sh_multipliers[-1] = 0

        return(self.q_vfp_q_sh_multipliers, self.pre_heat_fit_params, 
                 self.front_heat_fit_params,
                self.front_heat_start_index, 
                self.front_heat_last_index,
                self.pre_heat_start_index,
                self.pre_heat_last_index)

    def setCoupleParams(self, save_path, **kwargs):

        kwargs['fluid_yaml']['FixedParameters']['Preheat_StartIndex'] = pre_heat_start_index.item()
        kwargs['fluid_yaml']['FixedParameters']['Preheat_LastIndex'] = pre_heat_last_index.item()
        kwargs['fluid_yaml']['FixedParameters']['Frontheat_StartIndex'] = front_heat_start_index.item()
        kwargs['fluid_yaml']['FixedParameters']['Frontheat_LastIndex'] = front_heat_last_index.item()
        if pre_heat_last_index > 
        if pre_heat_fit_params is not None:
            np.savetxt(os.path.join(io_obj.fluid_input_path,"pre_heat_fit_param.txt"), pre_heat_fit_params)
            np.savetxt(os.path.join(io_obj.fluid_input_path,"front_heat_fit_param.txt"), front_heat_fit_params)

        np.savetxt(os.path.join(save_path,"qe.txt"), self.q_vfp_q_sh_multipliers)