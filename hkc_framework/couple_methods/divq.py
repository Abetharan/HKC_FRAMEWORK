import numpy as np 
import os


class DivQ:
    def __init__(self):
        self.limit_density = False

    def method(self, sh_heat_flow, vfp_heat_flow, laser_dir = None, **kwargs):
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
        heat_flow = vfp_heat_flow
        mass = kwargs['mass']
        #Start limit-density
        limit_index = len(heat_flow)
        if laser_dir is not None:
            if laser_dir == "right":
                len_sh = len(sh_heat_flow)
                append_index = len_sh - limit_index + 1
                heat_flow = heat_flow[1:]
                heat_flow =  np.append(sh_heat_flow[:append_index], heat_flow)
            else:
                heat_flow = heat_flow[:-1]
                heat_flow = np.append(heat_flow, sh_heat_flow[limit_index - 1:])
        #End

        nx = len(heat_flow) 
        self.HeatConductionE = np.zeros(nx - 1)
        for i, m in enumerate(mass):
            self.HeatConductionE[i] = (-(heat_flow[i + 1] - heat_flow[i])) / m #- sign is there because of convention used in HyKiCT 

    def setCoupleParams(self, save_path, **kwargs):
        np.savetxt(os.path.join(save_path,"qe.txt"), self.HeatConductionE)