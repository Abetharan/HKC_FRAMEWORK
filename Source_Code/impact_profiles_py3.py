'''
    This script generates profiles for IMPACT inputs
     
    
    
    
'''

import numpy as np

#---------------------------------------------------------------
def cos_prof(x,xlen,amp,pos,nwl,wid):
    xo = pos
    l_grid = xlen
    cossy = amp * np.cos(2*np.pi*nwl*(x - xo)*(1.0/l_grid))
    return cossy


#---------------------------------------------------------------

def db_tanh_prof(x,xlen,amp,pos,nwl,wid):
    xo = pos
    #l_grid = float(len(x))
    l_grid = xlen
    dubz = amp * ((1.0 + np.tanh(((x-xo)+wid/2)/(nwl*l_grid)))*(1.0 - np.tanh(((x-xo)-wid/2)/(nwl*l_grid)))/2 - 1  )
    return dubz

#---------------------------------------------------------------

def tanh_prof(x,xlen,amp,pos,nwl,wid):
    xo = pos
    ##l_grid = float(len(x))
    sing = amp *(( np.tanh((x-xo)/(nwl*xlen))))

    return sing
#---------------------------------------------------------------
def sh_tanh_prof(x,xlen,amp,pos,nwl,wid):
    #xlen = x[-1]-x[0]
    profile = amp*0.5*(1.0 + np.tanh((x-pos)/(nwl*float(xlen))))

    #         val=val + prof_amp*0.5d0*
    #     +        (1.0d0 + tanh((xx-prof_pos)/(prof_sl*xlen)))
    return profile
#---------------------------------------------------------------

def sh_db_tanh_prof(x,xlen,amp,pos,nwl,wid):
        #xlen = len(x)
        profile= amp*0.25*((np.tanh(((x-pos)+wid/2.0)/(nwl*xlen)) +1.0)*((-1.0)*np.tanh(((x-pos)-wid/2.0)/(nwl*xlen)) + 1.0) )
        return profile
#---------------------------------------------------------------
def load_profile(nx,xmin,xmax,avg,amp,pos,nwl,wid,func):
    '''
    #    takes cos, tanh, db_tanh, custom profiles
    #    amp = 0.5*(9.9e22)
    #    pos = 2500.0
    #    nwl = 0.14
    #    wid = 10000.0
    
    generates IMPACT profile for 'func' and inputs assuming a UNIFORM grid
    '''
    step = (xmax - xmin)/nx
    x  = np.linspace(xmin,xmax,nx)
    fct_type = func[1:]
    p_m = func[0]
    slen = xmax-xmin
    print(' plus or minus: ', p_m )
    print('fct_type:       ', fct_type)
    
    xo = pos
    if fct_type == 'cos':
        y = cos_prof(x,slen,amp,pos,nwl,wid)
    if fct_type == 'tanh':
        y = tanh_prof(x,slen,amp,pos,nwl,wid)
    if fct_type == 'db_tanh':
        y = db_tanh_prof(x,slen,amp,pos,nwl,wid)
    if fct_type == 'sh_db_tanh':
        y = sh_db_tanh_prof(x,slen,amp,pos,nwl,wid)
    if fct_type == 'sh_tanh':
        y = sh_tanh_prof(x,slen,amp,pos,nwl,wid)
    z = np.ones(nx)*avg
    if p_m == '+':
        z = z + y
    if p_m == '*':
        z = z*y
    return z

def load_profile_custom(x,slen,avg,amp,pos,nwl,wid,func):
    '''
    #    takes cos, tanh, db_tanh, custom profiles
    
    generates IMPACT profile for 'func' and inputs for a grid x with length slen= xmax-xmin
    '''
    #step = (xmax - xmin)/nx
    #x  = np.linspace(xmin,xmax,nx)
    fct_type = func[1:]
    p_m = func[0]
    #slen = xmax-xmin
    print(' plus or minus: ', p_m )
    print('fct_type:       ', fct_type)
    
    xo = pos
    if fct_type == 'cos':
        y = cos_prof(x,slen,amp,pos,nwl,wid)
    if fct_type == 'tanh':
        y = tanh_prof(x,slen,amp,pos,nwl,wid)
    if fct_type == 'db_tanh':
        y = db_tanh_prof(x,slen,amp,pos,nwl,wid)
    if fct_type == 'sh_db_tanh':
        y = sh_db_tanh_prof(x,slen,amp,pos,nwl,wid)
    if fct_type == 'sh_tanh':
        y = sh_tanh_prof(x,slen,amp,pos,nwl,wid)
    z = np.ones(len(x))*avg
    if p_m == '+':
        z = z + y
    if p_m == '*':
        z = z*y
    return z
#--------------------------------------------------------------
#--------------------------------------------------------------


