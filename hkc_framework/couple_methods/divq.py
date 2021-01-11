import numpy as np 



class DivQ:
    def __init__(self,fluid_yaml, start_kin):
        raise Exception('Not Implemented')

    def divQHeatFlow(self, step = 2,  laser_dir = None):
        """ Purpose: Find Div.Q
            Args:
                electron_thermal_flux = SOL-KiT heat flux in SI
                mass = areal mass.
                step = A singel step represent index to next cell-wall. 
                in the case of SOL-KiT where quantites are defiend in 
                cell-centres and walls, require a step of two. In
                HyKiCT quantites are only defined at one point, thus,
                a step of 1 is sufficient. 
                
            Returns: Div.Q
            NOTE: ONLY WORKS FOR SAME GRIDS
        """
        ##Logic relevant only for Limit Density methods 
        heat_flow = self.vfp_heat

        #Start limit-density
        limit_index = len(heat_flow)
        if laser_dir is not None:
            if laser_dir == "right":
                len_sh = len(self.spitzer_harm_heat)
                append_index = len_sh - limit_index + 1
                heat_flow = heat_flow[1:]
                heat_flow =  np.append(self.spitzer_harm_heat[:append_index], heat_flow)
            else:
                heat_flow = heat_flow[:-1]
                heat_flow = np.append(heat_flow, self.spitzer_harm_heat[limit_index - 1:])
        #End

        nx = len(heat_flow) 
        HeatConductionE = np.zeros(nx - 1)
        for i, m in enumerate(self.mass):
            HeatConductionE[i] = (-(heat_flow[i + 1] - heat_flow[i])) / m #- sign is there because of convention used in HyKiCT 
        return(HeatConductionE)