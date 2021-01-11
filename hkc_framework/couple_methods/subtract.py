import numpy as np 



class Subtract:
    def __init__(self,fluid_yaml, start_kin):
        raise Exception('Not Implemented')

    def getSubtractDq(self, laser_dir = None):
        
        ##Logic relevant only for Limit Density methods 
        heat_flow = self.vfp_heat

        #Start limit-density
        limit_index = len(heat_flow)
        if laser_dir is not None:
            len_sh = len(self.spitzer_harm_heat)
            if laser_dir == "right":
                #Remove first term as its 0 
                heat_flow = heat_flow[1:]
                append_index = len_sh - limit_index + 1
                #To retain spitzer-heat in fluid code remove the subtraction factor. 
                self.spitzer_harm_heat[:append_index] = 0 
                #Pad the VFP heat-flow as it wont be same shape as spitzer-harm
                padding = len_sh - limit_index + 1
                heat_flow = np.pad(heat_flow, (padding,0), 'constant', constant_values = (0,0))
            else:
                #Remove last term as its 0 
                heat_flow = heat_flow[:-1]
                #To retain spitzer-heat in fluid code remove the subtraction factor. 
                self.spitzer_harm_heat[limit_index - 1:] = 0 
                #Pad the VFP heat-flow as it wont be same shape as spitzer-harm
                padding = len_sh - limit_index + 1
                heat_flow = np.pad(heat_flow, (0,padding), 'constant', constant_values = (0,0))
        #End

        if self.snb:
            return heat_flow - self.q_snb
        else:
            return heat_flow - self.spitzer_harm_heat
