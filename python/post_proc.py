# Script to post-process high-resolution WRF model output.
# 
# Major tasks include computing the following for selected variables:
#   1. domain-averages to produce time series
#   2. vertical integrals
#   3. pressure-level vertical interpolation
# 
# This script leverages the CDO (https://code.mpimet.mpg.de/projects/cdo/) module via
# subprocess (executes as a terminal command) to generate basic time series, which helps
# for checking for RCE. This should be available either by loading as a module or
# installing into the conda/mamba kernel you're running.
# 
# James Ruppert
# jruppert@ou.edu
# 11/2/2024

from netCDF4 import Dataset
import numpy as np
from wrf import getvar, ALL_TIMES#, vinterp
from read_wrf_ideal import *
from post_proc_functions import *
import os
import sys

# If running do_2d_vars:
    # need to use mpirun with 18 CPUs
# If running do_acre:
    # need to use mpirun with 8 CPUs

# Options

test_process = "ctl"

# Basic 2D variables
do_2d_vars = True
# 2D ACRE variables
do_acre = False
# Special 2D variables
do_2d_special = False
# Basic 3D variables
do_3d_vars = False
# Special 3D variables
do_3d_special = False

########################################################
# Directories and model output specs
########################################################

# test_process = "ctl"
test_process = "large2km"

# wrfdir = "/glade/campaign/univ/uokl0049/ideal/largedom2/"+test_process+"/"
# wrfdir = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/wrf-ideal/largedom2/"+test_process+"/"
wrfdir = "/glade/derecho/scratch/ruppert/wrf-ideal/"+test_process+"/"
outdir = wrfdir+"post_proc/"
os.makedirs(outdir, exist_ok=True)

# # Get WRF file list, dimensions
wrffiles = get_wrf_file_list(wrfdir, "wrfout_d01")
hffiles = get_wrf_file_list(wrfdir, "hfout_d01")
lat, lon, nx1, nx2, nz, npd = wrf_dims(wrffiles[0])
nfiles = len(wrffiles)

# New vertical dimension for pressure levels
dp = 25 # hPa
pres = np.arange(1000, 25, -dp)
nznew = len(pres)

# Get variable lists
vars2d = var_list_2d()
# vars3d = var_list_3d()

# CDO command
def runshell(str_command):
    # process = subprocess.Popen(str_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    process = subprocess.Popen(str_command, shell=True, universal_newlines=True)
    process.wait()
    # lines = process.stdout.readlines()
    # for iline in lines:
    #     print(iline)
    return None

########################################################
# Use CDO to process basic 2D variables
########################################################

if do_2d_vars:

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()

    # for ivar in vars2d:
    ivar = comm.rank

    varname_str = vars2d[ivar].upper()

    # First, delete file if exists
    operation_str = 'rm -rf '+outdir+varname_str+'.nc'
    runshell(operation_str)

    operation_str1 = 'cdo mergetime'
    operation_str2 = ' -selname,'+varname_str+' '
    out_file = ' '+outdir+varname_str+'.nc'
    str_out = ' > '+outdir+varname_str+'.log 2>&1'

    # Create string array of file-specific commands
    cdo_line = [operation_str1]
    # for ifile in wrffiles:
    # for ifile in hffiles:
    for ifile in hffiles[:-1]:
        cdo_line.append(operation_str2+ifile)
    # Then join them into one string
    cdo_line_merged = " ".join(cdo_line)

    # Run CDO command
    runshell(cdo_line_merged+out_file+str_out)
    print("Done processing "+varname_str)
    comm.barrier()

    print("Done writing out 2D basic variables")

########################################################
# Use CDO to generate ACRE
########################################################

if do_acre:

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()

    # Remove first if exists
    operation_str = 'rm -rf '+outdir+'LWacre.nc '+outdir+'SWacre.nc'
    runshell(operation_str)

    if comm.rank == 0:
        operation_str = 'cdo sub '+outdir+'LWUPT.nc '+outdir+'LWDNT.nc '+outdir+'lw_t.nc'
    if comm.rank == 1:
        operation_str = 'cdo sub '+outdir+'LWUPB.nc '+outdir+'LWDNB.nc '+outdir+'lw_b.nc'
    if comm.rank == 2:
        operation_str = 'cdo sub '+outdir+'LWUPTC.nc '+outdir+'LWDNTC.nc '+outdir+'lw_tC.nc'
    if comm.rank == 3:
        operation_str = 'cdo sub '+outdir+'LWUPBC.nc '+outdir+'LWDNBC.nc '+outdir+'lw_bC.nc'
    if comm.rank == 4:
        operation_str = 'cdo sub '+outdir+'SWUPT.nc '+outdir+'SWDNT.nc '+outdir+'sw_t.nc'
    if comm.rank == 5:
        operation_str = 'cdo sub '+outdir+'SWUPB.nc '+outdir+'SWDNB.nc '+outdir+'sw_b.nc'
    if comm.rank == 6:
        operation_str = 'cdo sub '+outdir+'SWUPTC.nc '+outdir+'SWDNTC.nc '+outdir+'sw_tC.nc'
    if comm.rank == 7:
        operation_str = 'cdo sub '+outdir+'SWUPBC.nc '+outdir+'SWDNBC.nc '+outdir+'sw_bC.nc'

    if comm.rank < 8:
        runshell(operation_str)
    comm.barrier()

    if comm.rank == 0:
        operation_str = 'cdo sub '+outdir+'lw_b.nc '+outdir+'lw_t.nc '+outdir+'lw_net.nc'
    if comm.rank == 1:
        operation_str = 'cdo sub '+outdir+'lw_bC.nc '+outdir+'lw_tC.nc '+outdir+'lw_netC.nc'
    if comm.rank == 2:
        operation_str = 'cdo sub '+outdir+'sw_b.nc '+outdir+'sw_t.nc '+outdir+'sw_net.nc'
    if comm.rank == 3:
        operation_str = 'cdo sub '+outdir+'sw_bC.nc '+outdir+'sw_tC.nc '+outdir+'sw_netC.nc'

    if comm.rank < 4:
        runshell(operation_str)
    comm.barrier()

    # Calculate the longwave ACRE
    if comm.rank == 0:
        operation_str = 'cdo sub '+outdir+'lw_net.nc '+outdir+'lw_netC.nc '+outdir+'LWacre.nc'
    # Calculate the shortwave ACRE
    if comm.rank == 1:
        operation_str = 'cdo sub '+outdir+'sw_net.nc '+outdir+'sw_netC.nc '+outdir+'SWacre.nc'

    if comm.rank < 2:
        runshell(operation_str)
    comm.barrier()

    # Delete unneeded files
    if comm.rank == 0:
        operation_str = 'rm -rf '+outdir+'LWUPT.nc '+outdir+'LWDNT.nc '+outdir+'LWUPB.nc '+outdir+'LWDNB.nc'
        runshell(operation_str)
        operation_str = 'rm -rf '+outdir+'LWUPTC.nc '+outdir+'LWDNTC.nc '+outdir+'LWUPBC.nc '+outdir+'LWDNBC.nc'
        runshell(operation_str)
        operation_str = 'rm -rf '+outdir+'SWUPT.nc '+outdir+'SWDNT.nc '+outdir+'SWUPB.nc '+outdir+'SWDNB.nc'
        runshell(operation_str)
        operation_str = 'rm -rf '+outdir+'SWUPTC.nc '+outdir+'SWDNTC.nc '+outdir+'SWUPBC.nc '+outdir+'SWDNBC.nc'
        runshell(operation_str)
        operation_str = 'rm -rf '+outdir+'lw_t.nc '+outdir+'lw_b.nc '+outdir+'lw_tC.nc '+outdir+'lw_bC.nc'
        runshell(operation_str)
        operation_str = 'rm -rf '+outdir+'sw_t.nc '+outdir+'sw_b.nc '+outdir+'sw_tC.nc '+outdir+'sw_bC.nc'
        runshell(operation_str)
        operation_str = 'rm -rf '+outdir+'lw_net.nc '+outdir+'lw_netC.nc '+outdir+'sw_net.nc '+outdir+'sw_netC.nc'
        runshell(operation_str)

    print("Done writing out ACRE variables")

########################################################
# Process special 2D variables
########################################################

if do_2d_special:

    # Read in variable from WRF files
    for ifile in range(nfiles):

        # Open the WRF file
        file = wrffiles[ifile]
        print("Processing "+file)
        ds = Dataset(file)

        qv = getvar(ds, "QVAPOR", timeidx=ALL_TIMES)#, cache=cache)
        pwrf = getvar(ds, "p", units='Pa', timeidx=ALL_TIMES)#, cache=cache)
        # hght = getvar(dset, "zstag", units='m', timeidx=ALL_TIMES)#, cache=cache)
        # tmpk = getvar(dset, "tk", timeidx=ALL_TIMES)#, cache=cache)
        # rho = density_moist(tmpk, qv, pwrf)

        # Get dz
        # dz = np.zeros(qv.shape)
        # for iz in range(nz):
        #     dz[:,iz] = hght[:,iz+1] - hght[:,iz]
        # Get dp
        dp = pwrf.differentiate('bottom_top')*-1

        # Read in variables

        # rainrate
        var = getvar(ds, 'RAINNC', timeidx=ALL_TIMES)
        if ifile == 0:
            rainnc_all = var
        else:
            rainnc_all = xr.concat((rainnc_all, var), 'Time')
        # pclass
        var = wrf_pclass(ds, dp)
        if ifile == 0:
            pclass_all = var
        else:
            pclass_all = xr.concat((pclass_all, var), 'Time')
        # pw
        var = vert_int(qv, dp)
        if ifile == 0:
            pw_all = var
        else:
            pw_all = xr.concat((pw_all, var), 'Time')
        # pw_sat
        qvsat = get_rv_sat(ds, pwrf)
        qvsat = xr.DataArray(qvsat, coords=qv.coords, dims=qv.dims, attrs=qv.attrs)
        var = vert_int(qvsat, dp)
        if ifile == 0:
            pw_sat_all = var
        else:
            pw_sat_all = xr.concat((pw_sat_all, var), 'Time')

        ds.close()

    # Calculate rain rate as centered difference
    nt_all = rainnc_all.shape[0]
    rainrate = rainnc_all.copy()
    rainrate[0] = 0
    rainrate[nt_all-1] = np.nan
    rainrate[1:-1] = (rainnc_all.values[2:] - rainnc_all.values[:-2])*0.5
    rainrate *= npd # mm/time step --> mm/day

    # Write out the variables
    var_name='pclass'
    write_ncfile(outdir, pclass_all, var_name)
    var_name='rainrate'
    write_ncfile(outdir, rainrate, var_name)
    var_name='pw'
    write_ncfile(outdir, pw_all, var_name)
    var_name='pw_sat'
    write_ncfile(outdir, pw_sat_all, var_name)

    print("Done writing out special 2D variables")

########################################################
########################################################
