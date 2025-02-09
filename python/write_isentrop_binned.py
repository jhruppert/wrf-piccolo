# Script to write out isentropic-binned data from WRF output.
# 
# Assumes output is in a single netcdf file on pressure levels.
# 
# James Ruppert  
# jruppert@ou.edu  
# 1/28/25

import numpy as np
from thermo_functions import *
from post_proc_functions import *
from mpi4py import MPI
from write_isentrop_functions import *
import pickle
from read_wrf_piccolo import *

comm = MPI.COMM_WORLD
nproc = comm.Get_size()
print("Number of processes: ",nproc)

########################################################
# Main settings
########################################################

# Skips processing of variable if set to True and file is found
dont_overwrite=True
dont_overwrite=False

# Time bounds of processed 3D variables
t0_3d = np.datetime64('2024-09-02T00:00:00')
t1_3d = np.datetime64('2024-09-03T00:00:00')

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

# Get date tag for post_proc output files
def get_datetag(datetime):
    string = np.datetime_as_string(datetime, unit='m').replace("-","").replace(" ","").replace(":","")
    return string
t0_str = get_datetag(t0_3d)
t1_str = get_datetag(t1_3d)
tag_postproc = '_'+t0_str+'-'+t1_str

########################################################
# Main functions
########################################################

def run_binning(ipclass, theta_e, invar, pclass):

    bins = theta_e_bins()

    nt = theta_e.shape[0]
    nz = theta_e.shape[1]
    nbins = bins.size

    # Expand PCLASS to 3D
    pclass_z = pclass.expand_dims(dim={theta_e.dims[1]: theta_e.interp_level}, axis=1)
    pclass_z = pclass_z.values

    # Mask out based on precipitation class
    # pclass_name = ['noncloud','deepc','congest','shallowc','strat','anvil','all']
    if ipclass <= 5:
        indices = (pclass_z != ipclass)
        theta_e_masked = np.ma.masked_where(indices, theta_e.to_masked_array())#, copy=True)
        invar_masked = np.ma.masked_where(indices, invar.to_masked_array())#, copy=True)
        # condition = (pclass_z.values == ipclass)
        # theta_e_masked = theta_e.where(condition)
        # invar_masked = invar.where(condition)
    elif ipclass == 6:
    #     # MCS = mask out NONCLOUD and SHALLOW
    #     indices = ((pclass_z != 0) & (pclass_z != 3))
    #     theta_e_masked = np.ma.masked_where(indices, theta_e, copy=True)
    #     invar_masked = np.ma.masked_where(indices, invar, copy=True)
    # elif ipclass == 7:
        # Unmasked for "ALL" category
        theta_e_masked = theta_e.to_masked_array()
        invar_masked = invar.to_masked_array()

    # Mean profiles
    theta_e_mean = np.ma.mean(theta_e_masked, axis=(2,3))
    invar_mean   = np.ma.mean(invar_masked, axis=(2,3))

    # Frequency of pclass vs. time
    pclass_count = np.ndarray(nt, dtype=np.int64)
    for it in range(nt):
        # pclass_count[it] = np.ma.count(theta_e_masked[it,2,:,:])
        pclass_count[it] = np.count_nonzero((pclass[it,...] == ipclass))

    # Transform the variables from (x,y) --> (bin)

    dims = (nt,nz,nbins-1)
    invar_binned = np.full(dims, np.nan)
    freq_binned = np.ndarray(dims, dtype=np.int64)

    nmin = 3 # minimum points to average

    for it in range(nt):
        for iz in range(nz):
            for ibin in range(nbins-1):
                indices = ((theta_e_masked[it,iz,:,:] >= bins[ibin]) & (theta_e_masked[it,iz,:,:] < bins[ibin+1])).nonzero()
                binfreq = indices[0].size
                freq_binned[it,iz,ibin] = np.array(binfreq, dtype=np.int64)
                # Take mean across ID'd cells
                if binfreq > nmin:
                    invar_binned[it,iz,ibin] = np.ma.mean(invar_masked[it,iz,indices[0],indices[1]])

    return bins, freq_binned, invar_binned, theta_e_mean, invar_mean, pclass_count

################################

def driver_loop_write_ncdf(memb_dir, datdir, tag_postproc, t0, t1):

    # Get variable metadata and PCLASS names
    proc_var_list, pclass_name = get_variable_list()
    npclass = len(pclass_name)

    # Read independent variables
    theta_e, pclass = read_vars_stage1(datdir, tag_postproc, t0=t0, t1=t1)
    pres = theta_e.interp_level.values

    for ivar in range(len(proc_var_list)):

        # Check if variable was already done
        if dont_overwrite:
            skip=True
            for ipclass in range(npclass):
                pickle_file_ivar = datdir+'isentrop/'+proc_var_list[ivar]+'_'+pclass_name[ipclass]+'_'+tag_postproc+'.pkl'
                if not os.path.isfile(pickle_file_ivar):
                    skip=False
                    break
            if skip:
                print("Skipping "+proc_var_list[ivar]+" for "+memb_dir)
                continue

        # Read variable to process
        invar = read_vars_stage2(datdir, tag_postproc, proc_var_list[ivar])

        print()
        print("Running variable: ",proc_var_list[ivar]," for "+memb_dir)

        # Loop over precipitation classes
        for ipclass in range(npclass):

            pickle_file_main = datdir+'isentrop/mainvars_'+pclass_name[ipclass]+'_'+tag_postproc+'.pkl'
            pickle_file_ivar = datdir+'isentrop/'+proc_var_list[ivar]+'_'+pclass_name[ipclass]+'_'+tag_postproc+'.pkl'

            # Check if processing of variable/ipclass is completed
            if dont_overwrite:
                if os.path.isfile(pickle_file_ivar):
                    print("Skipping "+proc_var_list[ivar]+'/'+pclass_name[ipclass]+" for "+memb_dir)
                    continue

            print()
            print("Running ipclass: ",pclass_name[ipclass]+" for "+memb_dir)

            bins, freq_binned, invar_binned, theta_e_mean, invar_mean, pclass_count = run_binning(ipclass, theta_e, invar, pclass)
            # run_binning(ipclass, theta_e, invar, pclass)
            # continue

            # Write out pickle files
            if ivar == 0:
                with open(pickle_file_main, 'wb') as file:
                    pickle.dump({'bins':bins,
                                'pres':pres,
                                'theta_e_mean':theta_e_mean,
                                'pclass_count':pclass_count,
                                'freq_binned':freq_binned,},
                                file)
            with open(pickle_file_ivar, 'wb') as file:
                pickle.dump({'invar_binned':invar_binned,
                            'invar_mean':invar_mean,},
                            file)

    return

########################################################
# Top-level loop
########################################################

# Loop over ensemble members

# for memb_dir in memb_all:

memb_dir = memb_all[comm.rank]

# if comm.rank == 0:

# Get dimensions and files
outdir, postproc_files, nt, nx, ny = get_postproc_dims(datdir, case, test_process, wrf_dom, memb_dir)

driver_loop_write_ncdf(memb_dir, outdir, tag_postproc, t0_3d, t1_3d)

print("Finished processing "+memb_dir)
