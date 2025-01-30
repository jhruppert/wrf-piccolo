# Script to write out isentropic-binned data from WRF output.
# 
# Assumes output is in a single netcdf file on pressure levels.
# 
# James Ruppert  
# jruppert@ou.edu  
# 1/28/25

import numpy as np
import xarray as xr
from thermo_functions import *
from post_proc_functions import *
from mpi4py import MPI
from write_isentrop_functions import *

# Parallelization notes:
#   Using mpi4py to distribute the work of reading and processing
#   large-dimensional numpy arrays. The processed results are then
#   passed back to the rank-0 node, which does the netcdf write-out.

# Testing mode:
#   all loops shortened to a single iteration and nt = 3
# testing=True
testing=False

comm = MPI.COMM_WORLD
nproc = comm.Get_size()

########################################################
# Main settings
########################################################

# Time bounds of processed 3D variables
t0_3d = np.datetime64('2024-09-02T00:00:00')
t1_3d = np.datetime64('2024-09-03T00:00:00')

proc_var_list = ['tmpk', 'theta_v', 'qvapor', 'rho', 'condh', 'rthratlw', 'rthratlwc', 'rthratsw', 'rthratswc', 'w']
nvars = len(proc_var_list)

pclass_name = ['noncloud','deepc','congest','shallowc','strat','anvil','all']
npclass = len(pclass_name)
npclass_orig = npclass-2

case = "sept1-4"
test_process = "ctl"
# test_process = "ncrf12h"

########################################################
# Directories and test selection
########################################################

wrf_dom = "wrf_fine"
nmem = 5 # number of ensemble members

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

# Kill all loops after single iteration for testing mode
if testing:
    nvars=1
    nmem = 1
    # npclass = 1

########################################################
# Main functions
########################################################

def run_binning(ipclass, bins, theta_e, invar, pclass_z):

    shape = theta_e.shape
    nt = shape[0]
    nz = shape[1]
    nbins = bins.size

    # Mask out based on precipitation class
    # pclass_name = ['noncloud','deepc','congest','shallowc','strat','anvil','all']
    if ipclass <= 5:
        indices = (pclass_z != ipclass)
        theta_e_masked = np.ma.masked_where(indices, theta_e, copy=True)
        invar_masked = np.ma.masked_where(indices, invar, copy=True)
    elif ipclass == 6:
    #     # MCS = mask out NONCLOUD and SHALLOW
    #     indices = ((pclass_z != 0) & (pclass_z != 3))
    #     theta_e_masked = np.ma.masked_where(indices, theta_e, copy=True)
    #     invar_masked = np.ma.masked_where(indices, invar, copy=True)
    # elif ipclass == 7:
        # Unmasked for "ALL" category
        theta_e_masked = theta_e
        invar_masked = invar

    # Frequency of cloud-type vs. time
    pclass_count = np.ndarray(nt, dtype=np.float64)
    for it in range(nt):
        pclass_count[it] = np.ma.count(theta_e_masked[it,2,:,:])

    theta_e_mean = np.ma.mean(theta_e_masked, axis=(2,3))

    # Bin the variables from (x,y) --> (bin)

    dims = (nt,nz,nbins-1)
    invar_binned = np.full(dims, np.nan)
    freq_binned = np.ndarray(dims, dtype=np.float64)

    nmin = 3 # minimum points to average

    for it in range(nt):
        for iz in range(nz):
            for ibin in range(nbins-1):
                indices = ((theta_e_masked[it,iz,:,:] >= bins[ibin]) & (theta_e_masked[it,iz,:,:] < bins[ibin+1])).nonzero()
                binfreq = indices[0].size
                freq_binned[it,iz,ibin] = np.array(binfreq, dtype=np.float64)
                # Take mean across ID'd cells
                if binfreq > nmin:
                    invar_binned[it,iz,ibin] = np.ma.mean(invar_masked[it,iz,indices[0],indices[1]])

    # return freq_binned, invar_binned, theta_e_mean, pclass_count
    return freq_binned.values, invar_binned.values, theta_e_mean.values, pclass_count.values

################################

def driver_loop_write_ncdf(comm, datdir, bins, t0, t1, proc_var_list):

    nt = dims[0]
    nz = dims[1]

    # Read variables
    theta_e, pclass_z, invar, pres = read_all_vars(comm, datdir,t0,t1,proc_var_list)

    for ipclass in range(npclass):
    # for ipclass in range(0,2):

        if comm.rank == 0:
            print()
            print("Running ipclass: ",pclass_name[ipclass])

        freq_binned, invar_binned, theta_e_mean, pclass_count = run_binning(ipclass,bins,theta_e,invar,pclass_z)

        # Consolidate rebinned data onto Rank0 and write to netCDF file

        if comm.rank > 0:

            comm.Send(np.ascontiguousarray(invar_binned, dtype=np.float64), dest=0, tag=comm.rank)

        else:
        # Rank 0 - all data write-out

            var_list_write=[]
            var_list_write.append(bins)
            var_list_write.append(pres)
            var_list_write.append(theta_e_mean)
            var_list_write.append(pclass_count)
            var_list_write.append(freq_binned)
            var_list_write.append(invar_binned)

            for irank in range(1,nvars):
                dims = (nt,nz,nbins-1)
                invar_binned = np.empty(dims, dtype=np.float64)
                comm.Recv(invar_binned, source=irank, tag=irank)
                # check that the unique arrays are appearing on process 0
                # print()
                # print(invar_binned[1,:,30])
                var_list_write.append(invar_binned)

            # Write out to netCDF file

            pclass_tag = pclass_name[ipclass]
            file_out = datdir+'binned_isentrop_'+pclass_tag+'.nc'
            var_names, descriptions, units, dims_set = var_regrid_metadata(nt,nz,nbins)
            write_ncfile(file_out, var_list_write, var_names, descriptions, units, dims_set)

    return

########################################################
# Top-level loop
########################################################

# #### Index aka Bin variable settings

# Theta-e (equivalent potential temperature)
fmin=305; fmax=365 # K
nbins = 70
bins=np.linspace(fmin,fmax,num=nbins)

# Loop over ensemble members

for memb_dir in memb_all[0:1]:

    if comm.rank == 0:
        print()
        print("Processing "+memb_dir)

    # Get dimensions and files
    outdir, postproc_files, nt, nx, ny = get_postproc_dims(datdir, case, test_process, wrf_dom, memb_dir)

    driver_loop_write_ncdf(comm, outdir, bins, t0_3d, t1_3d, proc_var_list)
