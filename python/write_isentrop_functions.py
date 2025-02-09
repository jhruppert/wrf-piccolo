# Functions for isentropic binning.
# 
# James Ruppert  
# jruppert@ou.edu  
# 1/28/25

import numpy as np
import xarray as xr
from thermo_functions import *

################################
# Theta-e (equivalent potential temperature) bins

def theta_e_bins():
    fmin=305; fmax=365 # K
    nbins = 70
    bins=np.linspace(fmin,fmax,num=nbins)
    return bins

################################
# Drops the edges when reading in variables

def drop_edges(var):
    buffer = 80
    return var.isel(south_north=slice(buffer, -buffer), west_east=slice(buffer, -buffer))

################################
# Main function for reading in variables

def read_var(datdir, varname, t0=None, t1=None, tag_extra=''):
    post_proc_file = datdir+varname+tag_extra+'.nc'
    ds = xr.open_dataset(post_proc_file)
    if t0 is not None:
        var = ds[varname].sel(Time=slice(t0,t1))
    else:
        var = ds[varname]
    # var = ds[varname].isel(Time=slice(20,22), south_north=slice(50,60), west_east=slice(50,60))
    ds.close()
    # Set large values to NaN
    var = var.where(var < 1e10)
    return drop_edges(var)

##########################################
# Use dictionary passing to get multiple
# uses of the same variable wherever possible
##########################################

# def var_readcheck(varname, dict_in, tag_extra=''):
#     if varname in dict_in:
#     # Check if exists
#         return dict_in[varname], dict_in
#     else:
#     # If not, read in
#         # if varname == 'qvapor':
#         #     dict_in[varname] = read_var(dict_in['datdir'], "qvapor", dict_in['t0'], dict_in['t1']) # kg/kg
#         # elif varname == 'tmpk':
#         #     dict_in[varname] = read_var(dict_in['datdir'], "tmpk", dict_in['t0'], dict_in['t1']) # kg/kg
#         # else:
#         dict_in[varname] = read_var(dict_in['datdir'], varname, tag_extra=tag_extra)#, dict_in['t0'], dict_in['t1']) # kg/kg
#         return dict_in[varname], dict_in

################################
# Read in theta-v

def get_theta_v(datdir, tag_extra=''):
    # tmpk, dict_in = var_readcheck('tmpk', dict_in, tag_extra=tag_extra) # K
    # qv, dict_in = var_readcheck('qvapor', dict_in, tag_extra=tag_extra) # kg/kg
    tmpk = read_var(datdir, 'tmpk', tag_extra=tag_extra)#, t0, t1) # K
    qv = read_var(datdir, 'qvapor', tag_extra=tag_extra)#, t0, t1) # K
    pres = tmpk.interp_level.values
    theta_v = theta_virtual(tmpk.to_masked_array(),
                            qv.to_masked_array(),
                            pres[np.newaxis, :, np.newaxis, np.newaxis]*1e2)
    return theta_v, tmpk

################################
# Get theta-v prime

def get_theta_v_prm(datdir, tag_extra=''):
    theta_v, tmpk = get_theta_v(datdir, tag_extra=tag_extra)
    theta_v_mn = np.ma.mean(theta_v, axis=(2,3))
    theta_v -= theta_v_mn[..., np.newaxis, np.newaxis]
    theta_v = xr.DataArray(theta_v, dims=tmpk.dims, coords=tmpk.coords)
    return theta_v

################################
# Read in primary variables

def read_vars_stage1(datdir, tag_postproc, t0=None, t1=None):
    # Read THETA-E
    varname = 'theta_e'
    theta_e = read_var(datdir, varname, tag_extra=tag_postproc)#, t0, t1) # K
    # Read PCLASS
    varname = 'pclass'
    pclass = read_var(datdir, varname, t0=t0, t1=t1)
    return theta_e, pclass

################################
# Read in secondary variables

def read_vars_stage2(datdir, tag_postproc, var_tag):
    if var_tag == 'theta_v':
        invar = get_theta_v_prm(datdir, tag_extra=tag_postproc)
    else:
        invar = read_var(datdir, var_tag, tag_extra=tag_postproc)#, t0, t1) # K
    return invar

# def read_all_vars(comm, datdir, tag_postproc, proc_var_list):

#     varname = 'theta_e'
#     theta_e = read_var(datdir, varname, tag_extra=tag_postproc)#, t0, t1) # K

#     # Read PCLASS
#     varname = 'pclass'
#     pclass = read_var(datdir, varname)#, t0, t1)

#     # Distribute variable processing onto all ranks.
#     # Rank0 then receives all processed results and does write-out.

#     dict_pass = {'datdir': datdir}#, 't0': t0, 't1': t1}

#     for ivar in range(len(proc_var_list)):
#     # if comm.rank == 1:
#         if ivar == 1:
#             # Calculating theta-v is memory intensive since it requires
#             # tmpk and qv, so do this first
#             # invar = get_theta_v(datdir,t0,t1,pres)
#             print("GETTING VAR: theta-v")
#             invar, dict_pass = get_theta_v_prm(dict_pass, tag_extra=tag_postproc)
#         else:
#             # invar, dict_pass = var_readcheck(proc_var_list[comm.rank], dict_pass, tag_extra=tag_postproc) # kg/m^3
#             invar, dict_pass = var_readcheck(proc_var_list[ivar], dict_pass, tag_extra=tag_postproc) # kg/m^3

#     return theta_e, pclass, invar, theta_e['interp_level']

################################

# def write_ncfile(file_out, var_list, var_names, descriptions, units, dims_set): #, dim_names

#     len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dims_set) #len4=len(dim_names)
#     if (len1 != len2) or (len1 != len3) or (len1 != len4):
#         raise ValueError("Variable info counts are off")

#     # dims_val = var_list[0].shape

#     ncfile = Dataset(file_out,mode='w', clobber=True)

#     # Add dimensions to file
#     # Will iterate over entire variable dimension list but will only attempt to add each once
#     dim_added = [] # List to track dimensions that have been added already
#     for idimset in range(len(dims_set)):
#     #     dim = ncfile.createDimension(dim_names[0][idim], dims_val[idim]) # unlimited axis (can be appended to).
#         for idim in range(len(dims_set[idimset][0])):
#             if dims_set[idimset][0][idim] in dim_added:
#                 continue
#             dim = ncfile.createDimension(dims_set[idimset][0][idim], dims_set[idimset][1][idim]) # unlimited axis (can be appended to).
#             dim_added.append(dims_set[idimset][0][idim])

#     for ivar in range(len(var_list)):
#         print("  Writing var: ",var_names[ivar])
#         writevar = ncfile.createVariable(var_names[ivar], np.single, dims_set[ivar][0]) #dim_names[ivar])
#         writevar.units = units[ivar]
#         writevar.description = descriptions[ivar]
#         writevar[...] = var_list[ivar]

#     ncfile.close()

#     print("Done writing!")
    
#     return None

#### NetCDF variable metadata

def get_variable_list():
    proc_var_list = [
        # 'tmpk',
        ### 'theta_v',
        # 'qvapor',
        'rho',
        # 'condh',
        # 'rthratlw',
        # 'rthratlwc',
        # 'rthratsw',
        # 'rthratswc',
        # 'w',
        ]
    pclass_name = [
        'noncloud',
        'deepc',
        'congest',
        'shallowc',
        'strat',
        'anvil',
        'all']
    return proc_var_list, pclass_name

# def var_regrid_metadata(nt,nz,nbins):

#     nbinsm1 = nbins-1

#     var_names = [
#         'bins',
#         'pres',
#         'theta_e_mn',
#         'pclass_frequency',
#         'frequency',
#         'tmpk',
#         'theta_v_prm',
#         'qvapor',
#         'rho',
#         'condh',
#         'rthratlw',
#         'rthratlwc',
#         'rthratsw',
#         'rthratswc',
#         'w',
#         'tmpk_mean',
#         'thv_mean',
#         'qv_mean',
#         'rho_mean',
#         'lw_mean',
#         'lwc_mean',
#         'sw_mean',
#         'swc_mean',
#         'w_mean',
#     ]
#     descriptions = [
#         'equivalent potential temperature bins',
#         'pressure',
#         'mean equivalent potential temperature',
#         'pclass frequency',
#         'frequency',
#         'temperature',
#         'virtual potential temperature xy-anomaly',
#         'water vapor mixing ratio',
#         'density',
#         'H_DIABATIC',
#         'LW heat tendency',
#         'LW clear-sky heat tendency',
#         'SW heat tendency',
#         'SW clear-sky heat tendency',
#         'vertical motion',
#         'mean temperature',
#         'mean virtual potential temperature',
#         'mean water vapor mixing ratio',
#         'mean density',
#         'mean LW heat tendency',
#         'mean LW clear-sky heat tendency',
#         'mean SW heat tendency',
#         'mean SW clear-sky heat tendency',
#         'mean vertical motion',
#     ]
#     units = [
#         'K',
#         'hPa',
#         'K',
#         'n-cells',
#         'n-cells',
#         'K',
#         'K',
#         'kg/kg',
#         'kg/m^3',
#         'K/s',
#         'K/s',
#         'K/s',
#         'K/s',
#         'K/s',
#         'm/s',
#         'K',
#         'K',
#         'kg/kg',
#         'kg/m^3',
#         'K/s',
#         'K/s',
#         'K/s',
#         'K/s',
#         'm/s',
#     ]
#     dims_all = (nt,nz,nbinsm1)
#     dim_names = ('nt','nz','nbinsm1')
#     dims_set = [
#         [('nbins',),(nbins,)],
#         [('nz',),(nz,)],
#         [('nt','nz'),(nt,nz)],
#         [('nt',),(nt,)],
#         [dim_names,dims_all],
#         [dim_names,dims_all],
#         [dim_names,dims_all],
#         [dim_names,dims_all],
#         [dim_names,dims_all],
#         [dim_names,dims_all],
#         [dim_names,dims_all],
#         [dim_names,dims_all],
#         [dim_names,dims_all],
#         [dim_names,dims_all],
#         [dim_names,dims_all],
#         [('nt','nz'),(nt,nz)],
#         [('nt','nz'),(nt,nz)],
#         [('nt','nz'),(nt,nz)],
#         [('nt','nz'),(nt,nz)],
#         [('nt','nz'),(nt,nz)],
#         [('nt','nz'),(nt,nz)],
#         [('nt','nz'),(nt,nz)],
#         [('nt','nz'),(nt,nz)],
#         [('nt','nz'),(nt,nz)],
#     ]

#     len1=len(var_names); len2=len(descriptions); len3=len(units); len4=len(dims_set) #len4=len(dim_names)
#     if (len1 != len2) or (len1 != len3) or (len1 != len4):
#         raise ValueError("Variable info counts are off")

#     return var_names, descriptions, units, dims_set