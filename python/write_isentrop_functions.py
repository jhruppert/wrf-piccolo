# Functions for isentropic binning.
# 
# James Ruppert  
# jruppert@ou.edu  
# 1/28/25

from netCDF4 import Dataset
import numpy as np
import xarray as xr
from thermo_functions import *

# Drop edges
def drop_edges(var):
    buffer = 80
    try:
        # 3D variable
        return var[:,buffer:-buffer,buffer:-buffer]
    except:
        # 2D variable
        return var[buffer:-buffer,buffer:-buffer]

#### NetCDF variable metadata

def var_regrid_metadata(nt,nz,nbins):

    nbinsm1 = nbins-1

    var_names = [
        'bins',
        'pres',
        'theta_e_mn',
        'pclass_frequency',
        'frequency',
        'tmpk',
        'theta_v_prm',
        'qv',
        'rho',
        'h_diabatic',
        'lw',
        'lwc',
        'sw',
        'swc',
        'w',
        'tmpk_mean',
        'thv_mean',
        'qv_mean',
        'rho_mean',
        'lw_mean',
        'lwc_mean',
        'sw_mean',
        'swc_mean',
        'w_mean',
    ]
    descriptions = [
        'equivalent potential temperature bins',
        'pressure',
        'mean equivalent potential temperature',
        'pclass frequency',
        'frequency',
        'temperature',
        'virtual potential temperature xy-anomaly',
        'water vapor mixing ratio',
        'density',
        'H_DIABATIC',
        'LW heat tendency',
        'LW clear-sky heat tendency',
        'SW heat tendency',
        'SW clear-sky heat tendency',
        'vertical motion',
        'mean temperature',
        'mean virtual potential temperature',
        'mean water vapor mixing ratio',
        'mean density',
        'mean LW heat tendency',
        'mean LW clear-sky heat tendency',
        'mean SW heat tendency',
        'mean SW clear-sky heat tendency',
        'mean vertical motion',
    ]
    units = [
        'K',
        'hPa',
        'K',
        'n-cells',
        'n-cells',
        'K',
        'K',
        'kg/kg',
        'kg/m^3',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'm/s',
        'K',
        'K',
        'kg/kg',
        'kg/m^3',
        'K/s',
        'K/s',
        'K/s',
        'K/s',
        'm/s',
    ]
    dims_all = (nt,nz,nbinsm1)
    dim_names = ('nt','nz','nbinsm1')
    dims_set = [
        [('nbins',),(nbins,)],
        [('nz',),(nz,)],
        [('nt','nz'),(nt,nz)],
        [('nt',),(nt,)],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [dim_names,dims_all],
        [('nt','nz'),(nt,nz)],
        [('nt','nz'),(nt,nz)],
        [('nt','nz'),(nt,nz)],
        [('nt','nz'),(nt,nz)],
        [('nt','nz'),(nt,nz)],
        [('nt','nz'),(nt,nz)],
        [('nt','nz'),(nt,nz)],
        [('nt','nz'),(nt,nz)],
        [('nt','nz'),(nt,nz)],
    ]

    len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dims_set) #len4=len(dim_names)
    if (len1 != len2) or (len1 != len3) or (len1 != len4):
        raise ValueError("Variable info counts are off")

    return var_names, descriptions, units, dims_set

################################

def read_var(datdir, varname, t0, t1):
    # Get date tag for output files
    def get_datetag(datetime):
        string = np.datetime_as_string(datetime, unit='m').replace("-","").replace(" ","").replace(":","")
        return string
    t0_str = get_datetag(t0)
    t1_str = get_datetag(t1)
    tag_extra = '_'+t0_str+'-'+t1_str
    post_proc_file = datdir+varname+tag_extra+'.nc'
    ds = xr.open_dataset(post_proc_file)
    var = ds[varname].sel(XTIME=slice(t0,t1))
    ds.close()
    # Set large values to NaN
    var = var.where(var < 1e10)
    return drop_edges(var)

##########################################
# Use dictionary passing to get multiple
# uses of the same variable wherever possible
##########################################

def var_readcheck(varname, dict_in):
    if varname in dict_in:
    # Check if exists
        return dict_in[varname], dict_in
    else:
    # If not, read in
        # if varname == 'qvapor':
        #     dict_in[varname] = read_var(dict_in['datdir'], "qvapor", dict_in['t0'], dict_in['t1']) # kg/kg
        # elif varname == 'tmpk':
        #     dict_in[varname] = read_var(dict_in['datdir'], "tmpk", dict_in['t0'], dict_in['t1']) # kg/kg
        # else:
        dict_in[varname] = read_var(dict_in['datdir'], varname, dict_in['t0'], dict_in['t1']) # kg/kg
        return dict_in[varname], dict_in

################################

def get_theta_v(dict_in):
    tmpk, dict_in = var_readcheck('tmpk', dict_in) # K
    qv, dict_in = var_readcheck('qv', dict_in) # kg/kg
    pres = tmpk['interp_level']
    pres = pres.expand_dims(dim={'new_dim0':1, 'new_dim2':1, 'new_dim3':1}, axis=(0,2,3))
    theta_v = theta_virtual(tmpk.values, qv.values, pres.values*1e2)
    theta_v.attrs = tmpk.attrs
    theta_v.attrs['units'] = 'K'
    theta_v.attrs['description'] = 'virtual potential temperature'
    return theta_v, dict_in

################################

def get_theta_v_prm(dict_in):
    theta_v, dict_in = get_theta_v(dict_in)
    theta_v_mn = theta_v.mean(dim=['south_north', 'west_east'], keep_attrs=True)
    theta_v -= theta_v_mn.expand_dims(dim={'new_dim2':1, 'new_dim3':1}, axis=(2,3))
    return theta_v, dict_in

################################

def read_all_vars(comm, datdir, t0, t1, proc_var_list):

    varname = 'theta_e'
    theta_e = read_var(datdir, varname, t0, t1) # K
    nz = theta_e.shape[0]

    varname = 'pclass'
    pclass = read_var(datdir, varname, t0, t1)
    # Expand pclass to 3D using xarray
    pclass_z = pclass.expand_dims({'z': np.arange(nz)})

    # Distribute variable processing onto all ranks.
    # Rank[0] then receives all processed results and does write-out.

    dict_pass = {'datdir': datdir, 't0': t0, 't1': t1}

    if comm.rank == 1:
        # Calculating theta-v is memory intensive since it requires
        # tmpk and qv, so do this first
        # invar = get_theta_v(datdir,t0,t1,pres)
        invar, dict_pass = get_theta_v_prm(dict_pass)
    else:
        invar, dict_pass = var_readcheck(proc_var_list[5 + comm.rank], dict_pass) # kg/m^3

    return theta_e, pclass_z, invar, theta_e['interp_level']
    # Pass back numpy variables
    # return theta_e.values, pclass_z.values, invar.values, theta_e['interp_level'].values

################################

def write_ncfile(file_out, var_list, var_names, descriptions, units, dims_set): #, dim_names

    len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dims_set) #len4=len(dim_names)
    if (len1 != len2) or (len1 != len3) or (len1 != len4):
        raise ValueError("Variable info counts are off")

    # dims_val = var_list[0].shape

    ncfile = Dataset(file_out,mode='w', clobber=True)

    # Add dimensions to file
    # Will iterate over entire variable dimension list but will only attempt to add each once
    dim_added = [] # List to track dimensions that have been added already
    for idimset in range(len(dims_set)):
    #     dim = ncfile.createDimension(dim_names[0][idim], dims_val[idim]) # unlimited axis (can be appended to).
        for idim in range(len(dims_set[idimset][0])):
            if dims_set[idimset][0][idim] in dim_added:
                continue
            dim = ncfile.createDimension(dims_set[idimset][0][idim], dims_set[idimset][1][idim]) # unlimited axis (can be appended to).
            dim_added.append(dims_set[idimset][0][idim])

    for ivar in range(len(var_list)):
        print("  Writing var: ",var_names[ivar])
        writevar = ncfile.createVariable(var_names[ivar], np.single, dims_set[ivar][0]) #dim_names[ivar])
        writevar.units = units[ivar]
        writevar.description = descriptions[ivar]
        writevar[...] = var_list[ivar]

    ncfile.close()

    print("Done writing!")
    
    return None