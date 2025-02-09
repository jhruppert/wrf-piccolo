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

def binvar_set_bins(binvar_tag):
    if binvar_tag == 'pw':
        fmin=35; fmax=80 # mm
        step=1
        # description='Column water vapor [mm]'
    if binvar_tag == 'sf':
        fmin=.3;fmax=1 # %
        step=0.01
        # description='Column saturation fraction'
    bins=np.arange(fmin,fmax+step,step)
    return bins#, description

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
    if varname == 'LWacre':
        ivarname = 'LWUPB'
    else:
        ivarname = varname
    if t0 is not None:
        try:
            var = ds[ivarname].sel(Time=slice(t0,t1))
        except:
            var = ds[ivarname].sel(XTIME=slice(t0,t1))
    else:
        var = ds[ivarname]
    # var = ds[varname].isel(Time=slice(20,22), south_north=slice(50,60), west_east=slice(50,60))
    ds.close()
    # Set large values to NaN
    var = var.where(var < 1e10)
    return drop_edges(var)

################################
# Read in primary variables

def read_vars_stage1(binvar_tag, datdir, t0=None, t1=None):
    # Read binvar
    if binvar_tag == 'pw':
        varname = 'pw'
        binvar = read_var(datdir, varname, t0=t0, t1=t1)
    elif binvar_tag == 'sf':
        binvar = read_var(datdir, 'pw', t0=t0, t1=t1) / read_var(datdir, 'pw_sat', t0=t0, t1=t1)
    # Read key 2D variables
    varname = 'pclass'
    pclass = read_var(datdir, varname, t0=t0, t1=t1)
    varname = 'rainrate'
    rain = read_var(datdir, varname, t0=t0, t1=t1)
    varname = 'LWacre'
    lwacre = read_var(datdir, varname, t0=t0, t1=t1)
    return binvar, pclass, rain, lwacre

#### NetCDF variable metadata

def get_variable_list():
    proc_var_list = [
        'tmpk',
        ### 'theta_v',
        'qvapor',
        'rho',
        'condh',
        'rthratlw',
        'rthratlwc',
        'rthratsw',
        'rthratswc',
        'w',
        ]
    pclass_name = [
        'noncloud',
        'deepc',
        'congest',
        'shallowc',
        'strat',
        'anvil',
        ]
    return proc_var_list, pclass_name
