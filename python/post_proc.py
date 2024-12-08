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
from post_proc_functions import *
import os
import xarray as xr

# Post-processing for each test needs to follow the below sequence:

# Basic 2D variables
do_2d_vars = False # Set select=1:ncpus=19:mpiprocs=19:ompthreads=1
# 2D ACRE variables
do_acre = False # Set select=1:ncpus=8:mpiprocs=8:ompthreads=1
# Rainrate
do_rainrate = False # Set select=1:ncpus=1:mpiprocs=1:ompthreads=1
# Reflectivity (lowest model level)
do_refl = False # Set select=1:ncpus=1:mpiprocs=1:ompthreads=1
# Special 2D variables
do_2d_special = True # Set select=1:ncpus=1:mpiprocs=1:ompthreads=1
# # Basic 3D variables
# do_3d_vars = False
# # Special 3D variables
# do_3d_special = False

########################################################
# Directories and test selection
########################################################

datdir = "/glade/derecho/scratch/ruppert/piccolo/"
# datdir = "/glade/campaign/univ/uokl0053/"
# datdir = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/piccolo/"

case = "sept1-4"
test_process = "ctl"
wrf_dom = "wrf_fine"
nmem = 5 # number of ensemble members

# Ens-member string tags (e.g., memb_01, memb_02, etc.)
memb0=1 # Starting member to read
memb_nums_str=np.arange(memb0,nmem+memb0,1).astype(str)
nustr = np.char.zfill(memb_nums_str, 2)
memb_all=np.char.add('memb_',nustr)

########################################################
# Use CDO to process basic 2D variables
########################################################

if do_2d_vars:

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()

    # Get variable list
    vars2d = var_list_2d()

    for memb_dir in memb_all:
    # for memb_dir in memb_all[0:1]:

        # for ivar in vars2d:
        ivar = comm.rank

        varname_str = vars2d[ivar].upper()

        outdir, wrffiles, nfiles, npd = memb_dir_settings(datdir, case, test_process, wrf_dom, memb_dir)
        cdo_merge_wrf_variable(outdir, wrffiles, varname_str)

        # comm.barrier()

    print("Done writing out 2D basic variables")

########################################################
# Get reflectivity (lowest level)
########################################################

if do_refl:

    # memb_dir = memb_all[comm.rank]
    for memb_dir in memb_all:
    # for memb_dir in memb_all[-1:]:

        print("Processing reflectivity for "+memb_dir)

        varname_str = 'REFL_10CM'
        outdir, wrffiles, nfiles, npd = memb_dir_settings(datdir, case, test_process, wrf_dom, memb_dir)
        cdo_merge_wrf_variable(outdir, wrffiles, varname_str)

        # comm.barrier()

    print("Done writing out 2D basic variables")

########################################################
# Use CDO to generate ACRE
########################################################

if do_acre:

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()

    for memb_dir in memb_all:

        outdir, wrffiles, nfiles, npd = memb_dir_settings(datdir, case, test_process, wrf_dom, memb_dir)

        # Remove first if exists
        operation_str = 'rm -rf '+outdir+'LWacre.nc '+outdir+'SWacre.nc'
        process = subprocess.Popen(operation_str, shell=True, universal_newlines=True)
        # runshell(operation_str)

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

        process = subprocess.Popen(operation_str, shell=True, universal_newlines=True)
        process.wait()
        # runshell(operation_str)

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
            process = subprocess.Popen(operation_str, shell=True, universal_newlines=True)
            process.wait()
            # runshell(operation_str)
        comm.barrier()

        # Calculate the longwave ACRE
        if comm.rank == 0:
            operation_str = 'cdo sub '+outdir+'lw_net.nc '+outdir+'lw_netC.nc '+outdir+'LWacre.nc'
        # Calculate the shortwave ACRE
        if comm.rank == 1:
            operation_str = 'cdo sub '+outdir+'sw_net.nc '+outdir+'sw_netC.nc '+outdir+'SWacre.nc'

        if comm.rank < 2:
            process = subprocess.Popen(operation_str, shell=True, universal_newlines=True)
            process.wait()
            # runshell(operation_str)
        comm.barrier()

        # Delete unneeded files
        if comm.rank == 0:
            for operation_str in [
                # 'rm -rf '+outdir+'LWUPT.nc '+outdir+'LWDNT.nc '+outdir+'LWUPB.nc '+outdir+'LWDNB.nc',
                # KEEP OLR
                'rm -rf '+outdir+'LWDNT.nc '+outdir+'LWUPB.nc '+outdir+'LWDNB.nc',
                'rm -rf '+outdir+'LWUPTC.nc '+outdir+'LWDNTC.nc '+outdir+'LWUPBC.nc '+outdir+'LWDNBC.nc',
                'rm -rf '+outdir+'SWUPT.nc '+outdir+'SWDNT.nc '+outdir+'SWUPB.nc '+outdir+'SWDNB.nc',
                'rm -rf '+outdir+'SWUPTC.nc '+outdir+'SWDNTC.nc '+outdir+'SWUPBC.nc '+outdir+'SWDNBC.nc',
                'rm -rf '+outdir+'lw_t.nc '+outdir+'lw_b.nc '+outdir+'lw_tC.nc '+outdir+'lw_bC.nc',
                'rm -rf '+outdir+'sw_t.nc '+outdir+'sw_b.nc '+outdir+'sw_tC.nc '+outdir+'sw_bC.nc',
                'rm -rf '+outdir+'lw_net.nc '+outdir+'lw_netC.nc '+outdir+'sw_net.nc '+outdir+'sw_netC.nc'
            ]:
                process = subprocess.Popen(operation_str, shell=True, universal_newlines=True)

    print("Done writing out ACRE variables")

########################################################
# Loop over ensemble members to calculate rainrate
########################################################

if do_rainrate:

    for memb_dir in memb_all:

        outdir, wrffiles, nfiles, npd = memb_dir_settings(datdir, case, test_process, wrf_dom, memb_dir)

        ds = xr.open_dataset(outdir+'RAINNC.nc')
        rainnc = ds['RAINNC']
        # Get rainrate
        rainrate = calculate_rainrate(rainnc, npd)
        # rainrate = xr.DataArray(rainrate, coords=rainnc.coords, dims=rainnc.dims, attrs=rainnc.attrs)
        # Write out
        var_name='rainrate'
        write_ncfile(outdir, rainrate, var_name)

        print("Done writing out rainrate")

########################################################
# Loop over ensemble members for special 2D variables
########################################################

if do_2d_special:

    var_list_special = ['pclass', 'pw', 'vmf']

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    if comm.rank == 0: print(nproc)
    print("Rank "+str(comm.rank))

    # if comm.rank < 5:

    # for memb_dir in memb_all:
    memb_dir = memb_all[comm.rank]
    outdir, wrffiles, nfiles, npd = memb_dir_settings(datdir, case, test_process, wrf_dom, memb_dir)

    print("Processing special 2D variables for "+memb_dir)

    # Read in variable from WRF files
    for ifile in range(nfiles):
    # for ifile in range(1):

        # Open the WRF file
        wrffile = wrffiles[ifile]
        print("Processing "+wrffile)

        # Get variables for entire file
        vars_ifile = get_2d_special_vars_iwrf(wrffile, var_list_special)
        # pclass_ifile, pw_ifile, pw_sat_ifile, vmf_ifile = get_2d_special_vars_iwrf(wrffile)

        # Concatenate variables
        if ifile == 0:
            # pclass
            pclass_all = pclass_ifile
            # pw
            pw_all = pw_ifile
            # pw_sat
            # pw_sat_all = pw_sat_ifile
            # vertical mass flux
            vmf_all = vmf_ifile
        else:
            pclass_all = xr.concat((pclass_all, pclass_ifile), 'Time')
            pw_all = xr.concat((pw_all, pw_ifile), 'Time')
            # pw_sat_all = xr.concat((pw_sat_all, pw_sat_ifile), 'Time')
            vmf_all = xr.concat((vmf_all, vmf_ifile), 'Time')

    # Remove duplicate time steps
    pclass_all = pclass_all.drop_duplicates(dim="Time", keep='first')
    pw_all     = pw_all.drop_duplicates(dim="Time", keep='first')
    # pw_sat_all = pw_sat_all.drop_duplicates(dim="Time", keep='first')

    # Write out the variables
    var_name='pclass'
    write_ncfile(outdir, pclass_all, var_name)
    var_name='pw'
    write_ncfile(outdir, pw_all, var_name)
    # var_name='pw_sat'
    # write_ncfile(outdir, pw_sat_all, var_name)

print("Done writing out special 2D variables")

########################################################
########################################################
