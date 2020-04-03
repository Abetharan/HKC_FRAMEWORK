
def Switches(
EE_0 = True,
EE_L = True,
EI_L = True,
DIAG_EE = False,
PERIODIC = False,
FIXEDUP = False,
FIXEDDOWN = False,
NOFLOWUP = False,
NOFLOWDOWN = False,
IMPACTMODE = False,
COLDIONS = False,
MAINTAIN = False,
OUTPUT_DENSITY = True,
OUTPUT_TEMPERATURE = True,
OUTPUT_VELOCITY = False,
OUTPUT_E_FIELD = False,
OUTPUT_SH_TEST = False,
):
    params = {
                'EE0' : str(EE_0)[0],
                'EEL' :  str(EE_L)[0],
                'EIL' : str(EI_L)[0],
                'TRI' :  str(DIAG_EE)[0],
                'PERIOD' :  str(PERIODIC)[0],
                'FIXEDUP' :  str(FIXEDUP)[0],
                'FIXEDDOWN' :  str(FIXEDDOWN)[0],
                'NOFLOWUP' :  str(NOFLOWUP)[0],
                'NOFLOWDOWN' :  str(NOFLOWDOWN)[0],
                'IMPACT_MODE' :  str(IMPACTMODE)[0],
                'COLD_ION' :  str(COLDIONS)[0],
                'MAINTAIN' : str(MAINTAIN)[0],
                'OD' : str(OUTPUT_DENSITY)[0],
                'OT' : str(OUTPUT_TEMPERATURE)[0],
                'OV' : str(OUTPUT_VELOCITY)[0],
                'OE' : str(OUTPUT_E_FIELD)[0],
                'OSH': str(OUTPUT_SH_TEST)[0]
                }
    return(params) 

def GRID(
nt = 100,
prestepnt = 4,
dt = 0.01,
predt = 0.00001,
save_freq = 10,
lmax = 1,
nv = 120,
dv = 0.5,
v_multi = 1,
nx = 84,
dx = 0.1,
smalldx = 0
):
    params = {
                'nt' : str(nt),
                'prestepnt' : str(prestepnt),
                'dt' : str(dt),
                'predt' : str(predt),
                'save_freq' : str(save_freq),
                'lmax' : str(lmax),
                'nv' : str(nv),
                'dv' : str(dv),
                'v_multi' : str(v_multi),
                'nx' : str(nx),
                'dx' : str(dx),
                'smalldx' : str(smalldx)
    }
    return(params)


def NORMS(
Z = 1,
Ar = 1,
Density = 1E20,
Temp = 100
):
    params = {
        'Z' : Z,
        'Ar' : Ar,
        'Density' : Density, 
        'Temp' : Temp 
    }

    return(params)

