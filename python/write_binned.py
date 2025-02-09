# Script to write out binned data based on a 2D index variable from WRF output.
# 
# Assumes output is in a single netcdf file on pressure levels.
# 
# James Ruppert  
# jruppert@ou.edu  
# 2/8/25

import numpy as np
from thermo_functions import *
from post_proc_functions import *
from mpi4py import MPI
from write_binned_functions import *
import pickle
from read_wrf_piccolo import *

comm = MPI.COMM_WORLD
nproc = comm.Get_size()
print("Number of processes: ",nproc)

########################################################
# Main settings
########################################################

# Which variable to use as index for binning
# binvar_tag = 'pw'
binvar_tag = 'sf'

# Skips processing of variable if set to True and file is found
dont_overwrite=True
dont_overwrite=False

# Time bounds of processed 3D variables
t0_3d = np.datetime64('2024-09-02T00:00:00')
t1_3d = np.datetime64('2024-09-03T00:00:00')
# t1_3d = np.datetime64('2024-09-02T02:00:00')

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

def run_binning(bins, pclass_name, binvar, invar, pclass, rain, lwacre):

    nbins = bins.size
    nz = invar.shape[1]
    npclass = len(pclass_name)
    invar_reshaped = np.transpose(invar, (1,0,2,3))

    # Transform the variables from (x,y) --> (bin)

    # 2D variables
    freq_binned = np.ndarray((nbins-1), dtype=np.int64)
    rain_binned = np.full((nbins-1), np.nan)
    lwacre_binned = np.full((nbins-1), np.nan)
    pclass_freq = np.full((nbins-1,npclass), np.nan)

    # 3D variable
    invar_binned = np.full((nz,nbins-1), np.nan)

    nmin = 3 # minimum points to average

    for ibin in range(nbins-1):
        indices = ((binvar >= bins[ibin]) & (binvar < bins[ibin+1])).nonzero()
        binfreq = indices[0].size
        freq_binned[ibin] = np.array(binfreq, dtype=np.int64)
        # Take mean across ID'd cells
        if binfreq > nmin:
            rain_binned[ibin] = np.ma.mean(rain[indices[0],indices[1],indices[2]])
            lwacre_binned[ibin] = np.ma.mean(lwacre[indices[0],indices[1],indices[2]])
            for iz in range(nz):
                invar_binned[iz,ibin] = np.ma.mean(invar_reshaped[iz,indices[0],indices[1],indices[2]])
        # Get p-class frequency
        for ipclass in range(npclass):
            indices_pclass = ((binvar >= bins[ibin]) & (binvar < bins[ibin+1]) & 
                (pclass == ipclass)).nonzero()
            pclass_freq[ibin, ipclass] = indices_pclass[0].shape[0]

    return bins, freq_binned, invar_binned, pclass_freq, rain_binned, lwacre_binned

################################

def driver_loop_write_ncdf(binvar_tag, memb_dir, datdir, tag_postproc, t0, t1):

    # Get variable metadata and PCLASS names
    proc_var_list, pclass_name = get_variable_list()

    bins = binvar_set_bins(binvar_tag)

    # Read independent variables
    binvar, pclass, rain, lwacre = read_vars_stage1(binvar_tag, datdir, t0=t0, t1=t1)

    for ivar in range(len(proc_var_list)):

        # Check if variable was already done
        pickle_file = datdir+'binned_'+binvar_tag+'/'+proc_var_list[ivar]+tag_postproc+'.pkl'
        if dont_overwrite:
            if os.path.isfile(pickle_file):
                print("Skipping "+proc_var_list[ivar]+" for "+memb_dir)
                continue

        # Read variable to process
        invar = read_var(datdir, proc_var_list[ivar], tag_extra=tag_postproc)#, t0, t1)
        pres = invar.interp_level.values

        print()
        print("Running variable: ",proc_var_list[ivar]," for "+memb_dir)

        bins, freq_binned, invar_binned, pclass_freq, rain_binned, lwacre_binned = \
            run_binning(bins, pclass_name, binvar.to_masked_array(), invar.to_masked_array(), pclass.to_masked_array(), rain.to_masked_array(), lwacre.to_masked_array())

        # Write out pickle file
        with open(pickle_file, 'wb') as file:
            pickle.dump({'bins':bins,
                         'pres':pres,
                         'freq_binned':freq_binned,
                         'pclass_freq':pclass_freq,
                         'invar_binned':invar_binned,
                         'rain_binned':rain_binned,
                         'lwacre_binned':lwacre_binned},
                         file)

    return

########################################################
# Top-level loop
########################################################

# Loop over ensemble members

for memb_dir in memb_all[1:]:

# memb_dir = memb_all[comm.rank]

# if comm.rank == 0:

# Get dimensions and files
    outdir, postproc_files, nt, nx, ny = get_postproc_dims(datdir, case, test_process, wrf_dom, memb_dir)

    driver_loop_write_ncdf(binvar_tag, memb_dir, outdir, tag_postproc, t0_3d, t1_3d)

    print("Finished processing "+memb_dir)

print("Done rebinning data!!")
