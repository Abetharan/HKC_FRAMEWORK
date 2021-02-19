import numpy as np 
import os
from .coupling_method import CouplingMethod
class Subtract(CouplingMethod):
    def __init__(self):
        self.limit_density = False
    def method(self,sh_heat_flow, vfp_heat_flow, laser_dir = None, **kwargs):
        """
        Purpose: Calculate coupling parameter relevant to
                the subtract method. 
        Args:
            sh_heat_flow = Spitzer Harm Heat flow
            vfp_heat_flow = VFP heat-flow
            laser_dir = Laser direction for limit density method. 
                    Laser dir acts as a pseudo switch to engage
                    limit density search. Default None.
            **Kwargs = Named optional. 
            Required Optional = q_SNB. 
        Notes:
            Subtract method is a constant factor subtraction
            (or addition however you look at it) to spitzer-harm/SNB
            to get the VFP-heat flow.
        """
        ##Logic relevant only for Limit Density methods 
        heat_flow = vfp_heat_flow
        q_snb = kwargs['q_snb']
        #Start limit-density
        limit_index = len(heat_flow)
        if laser_dir is not None:
            len_sh = len(sh_heat_flow)
            if laser_dir == "right":
                #Remove first term as its 0 
                heat_flow = heat_flow[1:]
                append_index = len_sh - limit_index + 1
                #To retain spitzer-heat in fluid code remove the subtraction factor. 
                sh_heat_flow[:append_index] = 0 
                #Pad the VFP heat-flow as it wont be same shape as spitzer-harm
                padding = len_sh - limit_index + 1
                heat_flow = np.pad(heat_flow, (padding,0), 'constant', constant_values = (0,0))
            else:
                #Remove last term as its 0 
                heat_flow = heat_flow[:-1]
                #To retain spitzer-heat in fluid code remove the subtraction factor. 
                sh_heat_flow[limit_index - 1:] = 0 
                #Pad the VFP heat-flow as it wont be same shape as spitzer-harm
                padding = len_sh - limit_index + 1
                heat_flow = np.pad(heat_flow, (0,padding), 'constant', constant_values = (0,0))
        #End
        
        if q_snb is not None:
            self.subtract_factor = heat_flow - q_snb
        else:
            self.subtract_factor = heat_flow - sh_heat_flow

    def setCoupleParams(self, save_path, **kwargs):
        np.savetxt(os.path.join(save_path,"qe.txt"), self.subtract_factor)