# Metadata for processing WRF output: reads raw WRF output files and writes to new files
# with a single file per variable, with all time steps in the same file.
# 
# James Ruppert
# jruppert@ou.edu
# 1 Nov 2024

import xarray as xr
from wrf import getvar, ALL_TIMES, vinterp
from read_wrf_piccolo import *
import numpy as np
from precip_class import *
from thermo_functions import *
import os
import sys

##########################################
# Get ensemble member settings
##########################################

def memb_dir_settings(datdir, case, test_process, wrf_dom, memb_dir):
    wrfdir = datdir+case+'/'+memb_dir+'/'+test_process+"/"+wrf_dom+"/"
    outdir = wrfdir+"post_proc/"
    os.makedirs(outdir, exist_ok=True)
    # Get WRF file list, dimensions
    wrffiles = get_wrf_file_list(wrfdir, "wrfout_d01*")
    lat, lon, nx1, nx2, nz, npd = wrf_dims(wrffiles[0])
    nfiles = len(wrffiles)
    # Get XTIME arrays for each wrffile
    # xtime_all = []
    # for file in wrffiles:
    #     ds = Dataset(file)
    #     xtime_min = ds.variables['XTIME'][:]
    #     split = ds.variables['XTIME'].getncattr('units').split(' ')[2:4]
    #     start_date = np.datetime64('T'.join(split))
    #     xtime = np.array([start_date + np.timedelta64(int(x), 'm') for x in xtime_min])
    #     xtime_all.append(xtime)
    # New vertical dimension for pressure levels
    # dp = 25 # hPa
    # pres = np.arange(1000, 25, -dp)
    # nznew = len(pres)
    return outdir, wrffiles, nfiles, npd

##########################################
# Compute vertical integral
##########################################

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
# Use dictionary passing to get multiple
# uses of the same variable wherever possible
##########################################

def var_readcheck(varname, dict_pass):
    if varname in dict_pass:
    # Check if exists
        return dict_pass[varname], dict_pass
    else:
    # If not, create
        if varname == 'dp':
            if 'pwrf' in dict_pass:
                dict_pass[varname] = dict_pass['pwrf'].differentiate('bottom_top')*-1
            else:
                dict_pass['pwrf'] = getvar(dict_pass['ds'], "p", units='Pa', timeidx=dict_pass['timeidx']) # Pa
                dict_pass[varname] = dict_pass['pwrf'].differentiate('bottom_top')*-1 # Pa
        elif varname == 'qvapor':
            dict_pass[varname] = getvar(dict_pass['ds'], varname.upper(), timeidx=dict_pass['timeidx']) # kg/kg
        elif varname == 'pwrf':
            dict_pass[varname] = getvar(dict_pass['ds'], "p", units='Pa', timeidx=dict_pass['timeidx']) # Pa
        elif varname == 'tmpk':
            dict_pass[varname] = getvar(dict_pass['ds'], "tk", timeidx=dict_pass['timeidx']) # K
        return dict_pass[varname], dict_pass

##########################################
# For special 2D variables
##########################################

# Read and reduce variables for a given time step
def get_2d_special_vars_it(ds, it_file, var_list):
    vars_it = {}
    dict_pass = {'ds': ds, 'timeidx': it_file}
    for ivar_str in var_list:
        dp, dict_pass = var_readcheck('dp', dict_pass)
        if ivar_str == "pclass":
            ivar = wrf_pclass(ds, dp, it_file)
        elif ivar_str == "pw":
            qvapor, dict_pass = var_readcheck('qvapor', dict_pass)
            ivar = vert_int(qvapor, dp)
        elif ivar_str == "pw_sat":
            tmpk, dict_pass = var_readcheck('tmpk', dict_pass)
            pwrf, dict_pass = var_readcheck('pwrf', dict_pass)
            qv_s = rv_saturation(tmpk, pwrf) # kg/kg
            qv_s.attrs = tmpk.attrs
            qv_s.attrs['description'] = 'saturation water vapor mixing ratio'
            qv_s.attrs['units'] = 'kg/kg'
            ivar = vert_int(qv_s, dp)
        elif ivar_str == "vmf":
            wa = getvar(ds, "wa", timeidx=it_file)
            ivar = vert_int(wa, dp)
        vars_it[ivar_str] = ivar
    return vars_it

##########################################
# For 3D variables, including vertical interp
##########################################

# Read and vertically interpolate variables for a given time step
# def get_3d_vars_it(ds, it_file, var_list, new_p_levels):
def get_3d_vars_it(ds, it_file, ivar_str, new_p_levels):
    vars_it = {}
    dict_pass = {'ds': ds, 'timeidx': it_file}
    # for ivar_str in var_list:
    if ivar_str == "qvapor":
        ivar_ml, dict_pass = var_readcheck('qvapor', dict_pass)
    elif ivar_str == "w":
        ivar_ml = getvar(ds, "wa", timeidx=it_file) # m/s
    elif ivar_str == "condh":
        ivar_ml = getvar(ds, "H_DIABATIC", timeidx=it_file) # K/s
    elif ivar_str == "theta_e":
        qvapor, dict_pass = var_readcheck('qvapor', dict_pass) # kg/kg
        tmpk, dict_pass = var_readcheck('tmpk', dict_pass) # K
        pwrf, dict_pass = var_readcheck('pwrf', dict_pass) # Pa
        ivar_ml = theta_equiv(tmpk, qvapor, qvapor, pwrf) # K
        ivar_ml.attrs = tmpk.attrs
        ivar_ml.attrs['description'] = 'equivalent potential temperature'
    elif ivar_str == "rho":
        qvapor, dict_pass = var_readcheck('qvapor', dict_pass) # kg/kg
        tmpk, dict_pass = var_readcheck('tmpk', dict_pass) # K
        pwrf, dict_pass = var_readcheck('pwrf', dict_pass) # Pa
        ivar_ml = density_moist(tmpk, qvapor, pwrf) # kg/m^3
        ivar_ml.attrs = tmpk.attrs
        ivar_ml.attrs['description'] = 'density'
        ivar_ml.attrs['units'] = 'kg/m^3'
    elif ivar_str == "tmpk":
        ivar_ml, dict_pass = var_readcheck(ivar_str, dict_pass) # K
    else:
        ivar_ml = getvar(ds, ivar_str.upper(), timeidx=it_file) # K/s
    # Interpolate to new pressure levels
    # vars_it[ivar_str] = vinterp(ds, ivar_ml, 'p', new_p_levels, field_type='p_hpa', extrapolate=False)
    var_it = vinterp(ds, ivar_ml, 'p', new_p_levels, field_type='p_hpa', extrapolate=False)
    return var_it

##########################################
# Loop over WRF input file time steps to
# read and process variables
##########################################

# Using an efficient method of skipping repeat time steps, so
# the job of remove_duplicates is already done here.

# def get_vars_ifile(file, var_list, xtime_read, t0, t1, tag='2D', new_p_levels=None):
def get_vars_ifile(file, ivar_str, xtime_read, t0, t1, tag='2D', new_p_levels=None):
    ds = Dataset(file)
    nt_file = ds.dimensions['Time'].size
    xtime_min = ds.variables['XTIME'][:]
    split = ds.variables['XTIME'].getncattr('units').split(' ')[2:4]
    start_date = np.datetime64('T'.join(split))
    xtime_file = np.array([start_date + np.timedelta64(int(x), 'm') for x in xtime_min])
    # Loop over dataset time steps
    # vars_ifile = {}
    var_ifile = None
    for it_file in range(nt_file):
        # Check if within target time range
        if xtime_file[it_file] < t0 or xtime_file[it_file] > t1:
            continue
        # Next check if this time has already been added to xtime_read
        if xtime_file[it_file] in xtime_read:
            continue
        print('Processing time ',xtime_file[it_file])
        # Otherwise, add new time to xtime_read and process variables
        xtime_read = np.append(xtime_read, xtime_file[it_file])
        # 2D variables or 3D?
        if tag == '2D':
            print('UPDATE THIS FOR 2D')
            sys.exit()
            vars_it = get_2d_special_vars_it(ds, it_file, var_list)
        elif tag == '3D':
            # vars_it = get_3d_vars_it(ds, it_file, var_list, new_p_levels)
            var_it = get_3d_vars_it(ds, it_file, ivar_str, new_p_levels)
        # for ivar_str in var_list:
            # if ivar_str in vars_ifile:
            #     vars_ifile[ivar_str] = xr.concat((vars_ifile[ivar_str], vars_it[ivar_str]), 'Time')
            # else:
            #     vars_ifile[ivar_str] = vars_it[ivar_str]
        try:
            var_ifile = xr.concat((var_ifile, var_it), 'Time')
        except:
            var_ifile = var_it.copy()
    ds.close()
    return var_ifile, xtime_read

##########################################
# Calculate rain rate as centered difference
##########################################

def calculate_rainrate(rainnc_all, npd):
    rainrate = rainnc_all.copy()
    nt_all = rainrate.shape[0]
    rainrate[0] = 0
    rainrate[nt_all-1] = np.nan
    rainrate[1:-1] = (rainnc_all.values[2:] - rainnc_all.values[:-2])*0.5
    rainrate *= npd # Convert to mm/day
    rainrate /= 24 # Convert to mm/hr
    return rainrate

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

##########################################
# Write out variable to a NetCDF file
##########################################

def write_ncfile(outdir, var, var_name, tag_extra=''):#, dims_set, pres=None, do_pres=False): #, dim_names
    # var = xr.DataArray(var, coords=dummyvar.coords, dims=dummyvar.dims, attrs=dummyvar.attrs, name=var_name)
    file_out = outdir+var_name+tag_extra+'.nc'
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
    pclass = precip_class(q_int.values)
    return xr.DataArray(pclass, coords=ivar[0].coords, dims=ivar[0].dims, attrs=ivar[0].attrs, name='cloud_class')

##########################################
# Full variable lists
##########################################

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
       'tmpk',
        # 'hght',
        # 'rthratlwcrf',
        # 'rthratswcrf',
        'theta_e',
        # 'mse',
       'rho',
        ]
    # return [
    #     # Variables that can be directly interpolated
    #    'qvapor',
    #    'w',
    #     ]

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
        # units = 'mm/day'
        units = 'mm/hr'
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
    elif var_name == 'vmf':
        description = 'vertical mass flux'
        units = 'kg/m^2/s'
        # dims = ('nt','nx1','nx2')
        # dim_set = [dims, (nt,nx1,nx2)]
    #######################################################
    # Basic 3D variables
    #######################################################
    elif var_name == 'qvapor':
        description = 'water vapor mixing ratio'
        units = 'kg/kg'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'qrain':
        description = 'rain water mixing ratio'
        units = 'kg/kg'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'qcloud':
        description = 'cloud water mixing ratio'
        units = 'kg/kg'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'qice':
        description = 'ice water mixing ratio'
        units = 'kg/kg'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'qsnow':
        description = 'snow mixing ratio'
        units = 'kg/kg'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'qgraupel':
        description = 'graupel mixing ratio'
        units = 'kg/kg'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'u':
        description = 'zonal wind'
        units = 'm/s'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'v':
        description = 'meridional wind'
        units = 'm/s'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'w':
        description = 'vertical velocity'
        units = 'm/s'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'condh':
        description = 'condensation heating'
        units = 'K/s'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'rthratlw':
        description = 'longwave heating rate'
        units = 'K/s'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'rthratsw':
        description = 'shortwave heating rate'
        units = 'K/s'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'rthratlwc':
        description = 'longwave heating rate, clear sky'
        units = 'K/s'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    elif var_name == 'rthratswc':
        description = 'shortwave heating rate, clear sky'
        units = 'K/s'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
    #######################################################
    # Special 3D variables
    #######################################################
    elif var_name == 'tmpk':
        description = 'temperature'
        units = 'K'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]
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
    elif var_name == 'rho':
        description = 'density'
        units = 'kg/m^3'
        # dims = ('nt','nz','nx1','nx2')
        # dim_set = [dims, (nt,nz,nx1,nx2)]

    return description, units#, dim_set

##########################################
# 
##########################################
