# Metadata for processing WRF output: reads raw WRF output files and writes to new files
# with a single file per variable, with all time steps in the same file.
# 
# James Ruppert
# jruppert@ou.edu
# 1 Nov 2024

import xarray as xr
from wrf import getvar, ALL_TIMES
from read_wrf_piccolo import *
import numpy as np
from precip_class import *
from thermo_functions import *

##########################################
# Compute vertical integral
##########################################

# def vert_int(var, rho, dz):
#     return np.sum(var*rho*dz, axis=1)
# def vert_int(var, dp):
#     g = 9.81 # m/s2
#     return np.sum(var*dp, axis=1)/g
# Using pressure, Xarray
def vert_int(var, dp):
    g = 9.81 # m/s2
    # return np.sum(var*dp, axis=1)/g
    var_int = (var*dp/g).sum(dim='bottom_top')
    # Setting keep_attrs=True doesn't work, since taking (var*dp) removes them
    # So, copy the attributes over from original var, then modify
    var_int.attrs = var.attrs
    var_int.attrs['units'] = var.attrs['units']+' * kg*m-2'
    var_int.attrs['description'] = var.attrs['description']+' vertically integrated'
    return var_int

##########################################
# Get saturation mixing ratio
##########################################

def get_rv_sat(ds, pwrf, timeidx):
    # Read in temp [K] and use that + pr to calculate saturation mixing ratio
    tmpk = getvar(ds, 'tk', timeidx=timeidx) # K
    rv_sat = rv_saturation(tmpk.values, pwrf.values) # kg/kg
    return rv_sat

##########################################
# Run shell command
##########################################

def runshell(str_command):
    # process = subprocess.Popen(str_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    process = subprocess.Popen(str_command, shell=True, universal_newlines=True)
    process.wait()
    # lines = process.stdout.readlines()
    # for iline in lines:
    #     print(iline)
    return None

##########################################
# Merge WRF variables into nc file
##########################################

def cdo_merge_wrf_variable(outdir, wrffiles, varname_str):

    # First, delete file if exists
    operation_str = 'rm -rf '+outdir+varname_str+'.nc'
    process = subprocess.Popen(operation_str, shell=True, universal_newlines=True)
    # runshell(operation_str)

    # Create string array of file-specific commands
    operation_str1 = 'cdo mergetime'
    # SELLEVIDX dramatically slows down CDO process, so only include it if necessary
    if varname_str == 'REFL_10CM':
        operation_str2 = ' -selname,'+varname_str+' -sellevidx,1 '
    else:
        operation_str2 = ' -selname,'+varname_str+' '
    out_file = outdir+varname_str+'.nc'
    cdo_line = [operation_str1]
    # for ifile in hffiles:
    # for ifile in wrffiles[:-1]:
    for ifile in wrffiles:
        cdo_line.append(operation_str2+ifile)
    # Then join them into one string
    cdo_line_merged = " ".join(cdo_line)

    # Run CDO command
    operation_str = cdo_line_merged+' '+out_file+' > '+outdir+varname_str+'.log 2>&1'
    process = subprocess.Popen(operation_str, shell=True, universal_newlines=True)
    process.wait()

    # Remove duplicate time steps (if necessary)
    remove_duplicate_times_ncfile(out_file)

    print("Done writing merge-file "+varname_str)
    # comm.barrier()
    return None

########################################################
# Use Xarray to remove repeating time steps
########################################################

# For an NC file
def remove_duplicate_times_ncfile(infile):
    outfile = infile.replace('.nc', '_unique.nc')
    ds = xr.open_dataset(infile)
    # Identify unique time records
    try:
        timedim='XTIME'
        _, index = np.unique(ds[timedim], return_index=True)
        ds_unique = ds.isel(XTIME=index)
    except:
        timedim='Time'
        _, index = np.unique(ds[timedim], return_index=True)
        ds_unique = ds.isel(Time=index)
    # create a new dataset which doesn't contain time record duplicates
    # Check if there are repeated time records
    if len(ds[timedim]) == len(ds_unique[timedim]):
        ds.close()
        # print('No repeated time records found.')
    else:
        ds.close()
        ds_unique.to_netcdf(outfile)
        ds_unique.close()
        process = subprocess.Popen('mv '+outfile+' '+infile, shell=True, universal_newlines=True)
        process.wait()
        # print('Repeated time records found and removed.')
    return None

# From an open Xarray dataset or dataarray
def remove_duplicate_times(invar):
    # Identify unique time records
    _, index = np.unique(invar['Time'], return_index=True)
    # remove duplicates
    return invar.isel(Time=index)

##########################################
# Write out variable to a NetCDF file
##########################################

def write_ncfile(outdir, var, var_name):#, dims_set, pres=None, do_pres=False): #, dim_names
    # var = xr.DataArray(var, coords=dummyvar.coords, dims=dummyvar.dims, attrs=dummyvar.attrs, name=var_name)
    file_out = outdir+var_name+'.nc'
    description, units = get_metadata(var_name)
    var.name = var_name
    var.attrs['description'] = description
    var.attrs['units'] = units
    try:
        var.attrs['projection'] = str(var.attrs['projection'])
    except:
        var.attrs['projection'] = ''
    var.to_netcdf(file_out)
    return None

##########################################
# Calculate cloud classification
##########################################

# def wrf_pclass(infile, rho, dz):
#     # Read in and vertically integrate mixing ratios
#     q_list = ['QCLOUD', 'QRAIN','QICE', 'QSNOW', 'QGRAUP']
#     q_var = []
#     for ivar in range(len(q_list)):
#         ivar = wrf_var_read(infile,q_list[ivar]) # kg/kg
#         q_var.append(ivar)
#     q_var = np.stack(q_var, axis=0)
#     q_int = np.sum(q_var*rho*dz, axis=2)
#     return precip_class(q_int)
# def wrf_pclass(infile, dp):
#     # Read in and vertically integrate mixing ratios
#     q_list = ['QCLOUD', 'QRAIN','QICE', 'QSNOW', 'QGRAUP']
#     q_var = []
#     for ivar in range(len(q_list)):
#         ivar = wrf_var_read(infile,q_list[ivar]) # kg/kg
#         q_var.append(ivar)
#     q_var = np.stack(q_var, axis=0)
#     g = 9.81 # m/s2
#     q_int = np.sum(q_var*dp[np.newaxis,...], axis=2)/g
#     return precip_class(q_int)
# Now using Xarray
def wrf_pclass(ds, dp, timeidx=1):
    # Read in and vertically integrate mixing ratios
    q_list = ['QCLOUD', 'QRAIN','QICE', 'QSNOW', 'QGRAUP']
    q_vars = []
    for icld in q_list:
        # ivar = wrf_var_read(infile,q_list[ivar]) # kg/kg
        ivar = getvar(ds, icld, timeidx=timeidx) # kg/kg
        q_vars.append(ivar)
    # q_var = np.stack(q_var, axis=0)
    q_vars = xr.concat(q_vars, dim='QVAR')
    g = 9.81 # m/s2
    # q_int = np.sum(q_var*dp[np.newaxis,...], axis=2)/g
    q_int = (q_vars*dp.expand_dims(dim={'QVAR':5}, axis=0)).sum(dim='bottom_top')/g
    pclass = precip_class(np.array(q_int))
    return xr.DataArray(pclass, coords=ivar[0].coords, dims=ivar[0].dims, attrs=ivar[0].attrs, name='cloud_class')

##########################################
# Full variable lists
##########################################

# 3D variables
def var_list_3d():
    return [
        # Variables that can be directly interpolated
        'qvapor',
        # 'qrain',
        # 'qcloud',
        # 'qice',
        # 'qsnow',
        # 'qgraupel',
        # 'u',
        # 'v',
        'w',
        'condh',
        'rthratlw',
        'rthratsw',
        'rthratlwc',
        'rthratswc',
        # Variables that require special care
        # 't',
        # 'hght',
        # 'rthratlwcrf',
        # 'rthratswcrf',
        # 'theta_e',
        # 'mse',
        # 'rho',
        ]

# 2D variables
def var_list_2d():
    return [
        # Variables that can be processed via CDO
        'RAINNC',
        'hfx',
        'lh',
        'lwupt',
        'swupt',
        'lwuptc',
        'swuptc',
        'lwupb',
        'swupb',
        'lwupbc',
        'swupbc',
        'lwdnt',
        'swdnt',
        'lwdntc',
        'swdntc',
        'lwdnb',
        'swdnb',
        'lwdnbc',
        'swdnbc',
        # Variables that require special care
        # 'lwacre',
        # 'swacre',
        # 'pclass',
        # 'rain',
        # 'pw',
        # 'pw_sat',
        ]

##########################################
# Get variable metadata
##########################################

def get_metadata(var_name):#, nt, nz, nx1, nx2):

    #######################################################
    # Special 2D variables
    #######################################################
    if var_name == 'pclass':
        description = 'cloud classification (0 = nocloud,1=deepc,2=congest,3=shall,4=strat,5=anvil)'
        units = '-'
        # dims = ('nt','pclass','nx1','nx2')
        # dim_set = [dims, (nt,6,nx1,nx2)]
    elif var_name == 'rainrate':
        description = 'rainrate (centered diff)'
        units = 'mm/day'
        # dims = ('nt','nx1','nx2')
        # dim_set = [dims, (nt,nx1,nx2)]
    elif var_name == 'pw':
        description = 'column integrated water vapor'
        units = 'mm'
        # dims = ('nt','nx1','nx2')
        # dim_set = [dims, (nt,nx1,nx2)]
    elif var_name == 'pw_sat':
        description = 'column integrated saturation water vapor'
        units = 'mm'
        # dims = ('nt','nx1','nx2')
        # dim_set = [dims, (nt,nx1,nx2)]
    elif var_name == 'refl_10cm':
        description = 'radar reflectivity at lowest model level'
        units = 'dBZ'
        # dims = ('nt','nx1','nx2')
        # dim_set = [dims, (nt,nx1,nx2)]
    #######################################################
    # Basic 3D variables
    #######################################################
    # elif var_name == 'qvapor':
    #     description = 'water vapor mixing ratio'
    #     units = 'kg/kg'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'qrain':
    #     description = 'rain water mixing ratio'
    #     units = 'kg/kg'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'qcloud':
    #     description = 'cloud water mixing ratio'
    #     units = 'kg/kg'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'qice':
    #     description = 'ice water mixing ratio'
    #     units = 'kg/kg'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'qsnow':
    #     description = 'snow mixing ratio'
    #     units = 'kg/kg'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'qgraupel':
    #     description = 'graupel mixing ratio'
    #     units = 'kg/kg'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'u':
    #     description = 'zonal wind'
    #     units = 'm/s'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'v':
    #     description = 'meridional wind'
    #     units = 'm/s'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'w':
    #     description = 'vertical velocity'
    #     units = 'm/s'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'condh':
    #     description = 'condensation heating'
    #     units = 'K/s'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'rthratlw':
    #     description = 'longwave heating rate'
    #     units = 'K/s'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'rthratsw':
    #     description = 'shortwave heating rate'
    #     units = 'K/s'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'rthratlwc':
    #     description = 'longwave heating rate, clear sky'
    #     units = 'K/s'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'rthratswc':
    #     description = 'shortwave heating rate, clear sky'
    #     units = 'K/s'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    #######################################################
    # Special 3D variables
    #######################################################
    # elif var_name == 't':
    #     description = 'temperature'
    #     units = 'K'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'hght':
    #     description = 'height'
    #     units = 'm'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'rthratlwcrf':
        description = 'longwave cloud-radiation forcing'
        units = 'K/s'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'rthratswcrf':
        description = 'shortwave cloud-radiation forcing'
        units = 'K/s'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'theta_e':
        description = 'equivalent potential temperature'
        units = 'K'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'mse':
        description = 'moist static energy'
        units = 'J/kg'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    # elif var_name == 'rho':
    #     description = 'density'
    #     units = 'kg/m^3'
    #     dims = ('nt','nz','nx1','nx2')
    #     dim_set = [dims, (nt,nz,nx1,nx2)]

    return description, units#, dim_set

##########################################
# 
##########################################
