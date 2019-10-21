'''
    Python loaders for IMPACT
    processes header file and arrays
    
    - load_dict(path,fprefix,var,time) -> will load variable array, grids and timestep
                                of a particular time step and return it as dictionary
'''
import numpy as np, re, os, sys
def get_startline(fname):
    '''  
        get the line where the data starts
    '''
    f = open(fname)
    f_list = f.readlines()
    out_line = 9
    for i in range(len(f_list)):
        if re.match(r'\s*\n',f_list[i]):
            out_line = i
            break
    return out_line

#-----------------------------------------------------------------------    
def list_to_float(input):
    list = input.split()
    arr = np.zeros((len(list)),dtype=float)
    arr[:] = list[:]
    return arr

#-----------------------------------------------------------------------    

def construct_fname(path,fprefix,var,time, iter_number):
    
    if var=='fo' or var=='fxX' or var=='fxY' or var=='fyX' or var=='fyY':
        suffix = '.xyv'
    else:
        suffix = '.xy'
        
    #fprefix = 'thydro_hi'
    #time = '00'
    if iter_number is not None:
        fname = path + '/' + fprefix + '_' + var + '.' + suffix
    else:
        fname = path + '/' + fprefix +'_' + var + '_' + time + suffix
    return fname
#-----------------------------------------------------------------------    
def load_dict(path,fprefix = None,var = None,time = None, iter_number = None):
    '''
        Gets the IMPACT header info
        dict = fpg_get_info(fname)
    '''
    if fprefix is None:
        fname = path
    else:
        fname = construct_fname(path,fprefix,var,time, iter_number)
    
    #dict = MyDict()
    #mat = np.loadtxt(fname,skiprows=out_l)
    
    info = open(fname, 'r')
    data = info.readlines()
    mat_info = data[:10]
    out_l = get_startline(fname)
    mat = np.loadtxt(fname,skiprows=out_l)
    dict = {}

    ##print mat_info
    time = float(mat_info[1])
    
    ndims = int(mat_info[2])
    if ndims ==2:
        ny = int(mat_info[3])
        y_grid = list_to_float(mat_info[4])
        #y_grid = np.array(mat_info[4])
        
        nx = int(mat_info[5])
        x_grid = list_to_float(mat_info[6])
        ##print '========'
        ##print 'nx: ',nx,'\nx_grid: \n',x_grid
        ##print 'np.shape: xgrid: ', np.shape(x_grid)
        ##print '============'
        mat = np.loadtxt(fname,skiprows=out_l)
        ##print np.shape(mat)
        grid = [y_grid,x_grid]
        hmat_final = mat
        v_grid = 0.0
    else:
        
        nv = int(mat_info[3])
        v_grid = list_to_float(mat_info[4])
        
        ny = int(mat_info[5])
        y_grid = list_to_float(mat_info[6])
    
        nx = int(mat_info[7])
        x_grid = list_to_float(mat_info[8])
        
        mat = np.loadtxt(fname,skiprows=out_l)    
        dim_array = [nv,ny,nx]
        dim_reverse = [nx,ny,nv]
        
        mat = mat.reshape(dim_reverse)
        hmat_final = np.transpose(mat)
        grid = [v_grid,y_grid,x_grid] #nv,ny,nx
        #print '############################'
        #print '   nv, ny, nx: ',nv, ny, nx
        #print '############################'
        dict['v_grid'] = v_grid
        dict['nv'] = nv
    dict['time'] = time
    dict['y_grid'] = y_grid
    dict['x_grid'] = x_grid
    dict['nx'] = nx
    dict['ny'] = ny
    dict['mat'] = hmat_final    

    return dict    
