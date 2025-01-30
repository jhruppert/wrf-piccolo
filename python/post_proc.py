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

import numpy as np
from post_proc_functions import *
import os
import xarray as xr

# Post-processing for each test needs to follow the below sequence:

# Basic 2D variables
do_2d_vars = False # Set select=1:ncpus=19:mpiprocs=19:ompthreads=1
# 2D ACRE variables
do_acre = False # Set select=1:ncpus=8:mpiprocs=8:ompthreads=1
# Rainrate
do_rainrate = False # Set select=5:ncpus=1:mpiprocs=1:ompthreads=1
# Reflectivity (lowest model level)
do_refl = False # Set select=1:ncpus=1:mpiprocs=1:ompthreads=1
# Special 2D variables
do_2d_special = False # Set select=5:ncpus=1:mpiprocs=1:ompthreads=1
vars_2dspecial = ['pw_sat']#['pclass', 'pw', 'vmf', 'pw_sat']
# Basic 3D variables - includes vertical pressure interpolation
do_3d_vars = True
# Special 3D variables
# do_3d_special = False

########################################################
# Directories and test selection
########################################################

case = "sept1-4"
# test_process = "ctl"
# test_process = "ncrf12h"

wrf_dom = "wrf_fine"
nmem = 5 # number of ensemble members

# Define time period to process 3D variables
t0_3d = np.datetime64('2024-09-02T00:00:00')
t1_3d = np.datetime64('2024-09-03T00:00:00')
# t1_3d = np.datetime64('2024-09-02T00:40:00')

# Scratch
# datdir = "/glade/derecho/scratch/ruppert/piccolo/"
# Campaign storage
datdir = "/glade/campaign/univ/uokl0053/"
# OSCER
# datdir = "/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/piccolo/"

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

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    print("Rank "+str(comm.rank))

    # for memb_dir in memb_all:
    memb_dir = memb_all[comm.rank]

    print("Processing rainrate for "+memb_dir)

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

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    print("Rank "+str(comm.rank))

    # for memb_dir in memb_all:
    memb_dir = memb_all[comm.rank]

    print("Processing special 2D variables for "+memb_dir)

    outdir, wrffiles, nfiles, npd = memb_dir_settings(datdir, case, test_process, wrf_dom, memb_dir)

    # Read in variable from WRF files
    vars_alltime = {}
    for ifile in range(nfiles):
    # for ifile in range(1):

        # Open the WRF file
        wrffile = wrffiles[ifile]
        print("Processing "+wrffile)

        # Get variables for entire file
        vars_ifile = get_vars_ifile(wrffile, vars_2dspecial, tag='2D')

        # Concatenate variables
        for ivar_str in vars_2dspecial:
            if ifile == 0:
                vars_alltime[ivar_str] = vars_ifile[ivar_str]
            else:
                vars_alltime[ivar_str] = xr.concat((vars_alltime[ivar_str], vars_ifile[ivar_str]), 'Time')

    for ivar_str in vars_2dspecial:
        # Remove duplicate time steps
        vars_alltime[ivar_str] = vars_alltime[ivar_str].drop_duplicates(dim="Time", keep='first')
        # Write out the variables
        write_ncfile(outdir, vars_alltime[ivar_str], ivar_str)

    print("Done writing out special 2D variables")

########################################################
# Loop over ensemble members to process 3D variables
########################################################

if do_3d_vars:

    for test_process in ["ctl", "ncrf12h"]:

        # Define new output pressure levels
        dp=25
        new_p_levels=np.arange(1000,0,-dp)
        # new_p_levels=np.array([500,400])

        # Get date tag for output files
        def get_datetag(datetime):
            # string = np.datetime_as_string(datetime, unit='h').replace("-","").replace(" ","").replace(":","")
            string = np.datetime_as_string(datetime, unit='m').replace("-","").replace(" ","").replace(":","")
            return string
        t0_str = get_datetag(t0_3d)
        t1_str = get_datetag(t1_3d)
        tag_extra = '_'+t0_str+'-'+t1_str

        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        nproc = comm.Get_size()

        # Get variable list
        vars3d = var_list_3d()

        # for memb_dir in memb_all:
        memb_dir = memb_all[comm.rank]

        print("Processing special 3D variables for "+memb_dir)

        outdir, wrffiles, nfiles, npd = memb_dir_settings(datdir, case, test_process, wrf_dom, memb_dir)

        # for ivar_str in vars3d:
        for ivar_str in vars3d[3:]:
        # for ivar_str in vars3d[0:1]:

            # Read in variable from WRF files
            # vars_alltime = {}
            xtime_read = np.array([], dtype='datetime64[s]')
            var_alltime = None
            for ifile in range(nfiles):

                # Open the WRF file
                wrffile = wrffiles[ifile]
                print()
                # print("Processing "+wrffile)
                print("Processing "+ivar_str+" for "+wrffile)

                # Get variables for entire file
                var_ifile, xtime_read = get_vars_ifile(wrffile, ivar_str, xtime_read, t0_3d, t1_3d, tag='3D', new_p_levels=new_p_levels)
                # Check if dictionary is empty
                # if not var_ifile:
                if var_ifile is None:
                    continue

                # Concatenate variable
                # Check if dectionary key exists
                # if ivar_str in vars_alltime:
                try:
                    var_alltime = xr.concat((var_alltime, var_ifile), 'Time')
                except:
                    var_alltime = var_ifile.copy()

            # Remove duplicate time steps
            # vars_alltime[ivar_str] = vars_alltime[ivar_str].drop_duplicates(dim="Time", keep='first')
            # Write out the variables
            write_ncfile(outdir, var_alltime, ivar_str, tag_extra=tag_extra)

    print("Done writing out 3D variables")