#!/bin/bash
# 
# Bash script to prepare and execute WRF for ensemble simulations with downscaling
# using ndown.exe for a two-domain setup.
# 
# WPS needs to be run before this step to generate the initial and boundary
# condition "met_em*" files. For that step, use scripts "download_gefs.py" and
# then "batch_wps_gefs.job" located in the "run_scripts" directory.
# 
# Once WPS is completed, the below code assumes an existing directory tree as follows:
#   in your scratch directory <scratch>:
#      <scratch>/case_name/memb_01/met_em/ (for each ensemble member)
#   where "case_name" is the name of the case (e.g., sept1-4) and "met_em/" contains
#   the "met_em*" files corresponding to the ensemble member.
# 
# Stages when using NDOWN:
#   1. Prepare coarse domain (REAL)
#     - alias: "real1"; subdirectory: wrf_coarse
#     - run real.exe for coarse domain
#   2. Run coarse domain (WRF)
#     - alias: "wrf1"; subdirectory: wrf_coarse
#     - run wrf.exe for coarse domain
#   3. Get new ICs in prep for ndown (REAL)
#     - alias: "real2"; subdirectory: wrf_fine
#     - run real.exe for coarse and fine domain
#   4. Ndown to prepare fine domain (REAL)
#     - alias: "real3"; subdirectory: wrf_fine
#     - run ndown.exe
#   5. Run fine domain (WRF)
#     - alias: "wrf2"; subdirectory: wrf_fine
#     - run wrf.exe for fine domain
# 
# James Ruppert
# 18 Nov 2024

###################################################
# Settings
###################################################

# Selection to run the different stages of REAL, NDOWN, WRF
# ...which should go in the order as listed
# stage="real1" # Run real for coarse domain
# stage="wrf1" # Run wrf for coarse domain
# stage="real2" # Get new ICs in prep for ndown
  # -- can run real2 concurrently with wrf1
# stage="real3" # Run ndown to prepare fine domain
# stage="wrf2" # Run wrf for fine domain
stage="wrfrst" # Restart wrf for fine domain

# Select case name
case_name="sept1-4"

# Select test (e.g., "ctl", "ncrf", etc.)
test_name="ctl"

# Number of ensemble members
nens=5

# Restart
restart_t_stamp="2024-09-02_18:00:00"

###################################################
# Supercomputer environment-specific settings

# Current ongoing projects (totals as of as of 11/16/24):
#  - XX UFSU0031: Allison's exploratory allocation (remaining: 465k)
#  - UOKL0053: James's PICCOLO large allocation (20M)
#  - UOKL0049: James's TC-CRF large allocation (22M)
#  - UOKL0056: Frederick's TC small allocation (500k)

  system='derecho'
  if [[ ${system} == 'derecho' ]]; then
    queue="main"
    if [[ $stage == *"real"* ]]; then
      bigN=5
    elif [[ $stage == *"wrf"* ]]; then
      bigN=33
    fi
    project_code="UFSU0031" # Project to charge core hours against
    node_line="select=${bigN}:ncpus=128:mpiprocs=128:ompthreads=1" # Batch script node line
    submit_command="qsub" # Job submission command
    mpi_command="mpiexec" # MPI executable command
    work_dir=${work}/wrf-piccolo # parent directory containing a bunch of stuff
    ensemb_dir=${scratch}/piccolo/${case_name} # where each ensemble simulation is run
    sourc_file=${work_dir}/bashrc_wrf_der # source file for setting environment variables
    wrf_run_dir=$work_dir/tests_compiled/$test_name # directory with all WRF run code for selected test
  elif [[ ${system} == 'oscer' ]]; then
    echo "NEED TO UPDATE THIS"
    exit 0
  fi

###################################################
# Case-specific WRF job settings

  if [[ ${case_name} == 'sept1-4' ]]; then
  # CTL job settings
    if [[ $stage == *"real"* ]]; then
      run_time='02:00' # HH:MM Job run time
    elif [[ $stage == *"wrf"* ]]; then
      run_time='12:00' # HH:MM Job run time
    fi
    # Mechanism-denial tests
    if [[ ${test_name} == 'ncrf' ]]; then
      restart_t_stamp="2024-09-02_00:00:00"
      run_time='05:00' # HH:MM Job run time
      ndays=1
      restart_base='ctl'
    fi
  fi

###################################################
# Set working subdirectory "wrf_dir"

  jobname="${case_name}_${test_name}"
  if [[ $stage == "real1" ]] || [[ $stage == "wrf1" ]]; then
    wrf_dir="wrf_coarse"
  elif [[ $stage == "real2" ]] || [[ $stage == "real3" ]] || [[ $stage == "wrf2" ]]; then
    wrf_dir="wrf_fine"
  fi


###################################################
# Start parent loop for each ensemble member
###################################################

for em in $(seq -w 01 $nens); do # Ensemble member
# for em in $(seq -w 01 04); do # Ensemble member
# for em in 05; do # Ensemble member

  cd $ensemb_dir

  # Create directory tree (e.g., as "memb_02/ctl/wrf") for ensemble member
  memb_dir="$ensemb_dir/memb_${em}"
  mkdir -p $memb_dir # -p ignores if directory already exists
  test_dir=$memb_dir/$test_name
  mkdir -p $test_dir
  mkdir -p $test_dir/$wrf_dir
  cd $test_dir/$wrf_dir

  # Delete rsl-out from previous steps
  /bin/rm -f rsl.*

  echo "Running ${stage} in $test_dir/$wrf_dir"


###################################################
# Running REAL and NDOWN stages

  if [[ $stage == *"real"* ]]; then

    # Link met_em* files to current directory
    ln -sf $memb_dir/met_em/met_em* .
    # Copy WRF run directory contents for selected test to current directory
    /bin/cp -rafL ${wrf_run_dir}/* .
    # List of additional output variables
    /bin/cp $work_dir/namelists/var_extra_output .
    # Source file for environmental modules
    /bin/cp $sourc_file ./bashrc_wrf
    # Remove extraneous namelist if exists and grab the one needed for the test
    /bin/rm -f namelist.input

    # Prepare start data for REAL, NDOWN
    if [[ $stage == "real1" ]]; then
    # Real step for coarse domain
      # Coarse namelist
      namelist_file=${work_dir}/namelists/namelist.input.wrf.${case_name}.${test_name}.ndown_pt1
      /bin/cp $namelist_file ./namelist.input
      exec_name="real.exe"
    elif [[ $stage == "real2" ]]; then
    # Real step: get new ICs in prep for ndown
      # Need 2-domain setup to run real, get new ICs at start-time of d02
      namelist_file=${work_dir}/namelists/namelist.input.wrf.${case_name}.${test_name}.ndown_pt1
      /bin/cp $namelist_file ./namelist.input
      # Modify namelist
      sed -i "s/start_hour.*/start_hour = 12, 12,/" namelist.input
      sed -i "s/max_dom.*/max_dom = 2,/" namelist.input
      exec_name="real.exe"
    elif [[ $stage == "real3" ]]; then
    # Ndown step
      # Need 2-domain setup to run ndown
      namelist_file=${work_dir}/namelists/namelist.input.wrf.${case_name}.${test_name}.ndown_pt1
      /bin/cp $namelist_file ./namelist.input
      # Modify namelist to get new ICs at start-time of d02 for ndown
      sed -i "s/start_hour.*/start_hour = 12, 12,/" namelist.input
      sed -i "s/interval_seconds.*/interval_seconds = 3600/" namelist.input
      sed -i "s/max_dom.*/max_dom = 2,/" namelist.input
      mv wrfinput_d02 wrfndi_d02
      mv wrflowinp_d02 wrflowinp_d01
      ln -sf ../wrf_coarse/wrfout_d* .
      exec_name="ndown.exe"
    fi

    # Create REAL batch script from header base script
    cat $work_dir/run_scripts/header_${system}.txt > batch_real.job

# Append the below lines to batch file [remove indentation]
echo "

source bashrc_wrf

# Run REAL
${mpi_command} ./${exec_name}
" >> batch_real.job

    # Modify batch script by replacing placeholders
    sed -i "s/PROJECT/${project_code}/g" batch_real.job
    sed -i "s/QUEUE/${queue}/g" batch_real.job
    sed -i "s/JOBNAME/${jobname}/g" batch_real.job
    sed -i "s/EMM/${em}/g" batch_real.job
    sed -i "s/NODELINE/${node_line}/g" batch_real.job
    # if [[ ${system} == 'oscer' ]]; then
    #   sed -i "s/NNODES/${smn}/g" batch_real.job
    # fi
    sed -i "s/TIMSTR/${run_time}/g" batch_real.job

    # Submit REAL job
    # if [[ `grep SUCCESS rsl.error.0000 | wc -l` -eq 0 ]]; then
      ${submit_command} batch_real.job > submit_real_out.txt
    # fi

###################################################
# Running WRF

  elif [[ $stage == *"wrf"* ]]; then

  #  JOBID=$(grep Submitted submit_real_out.txt | cut -d' ' -f 4)

    # Remove and replace for specific test
    /bin/rm -f namelist.input

    # Prepare start data for WRF
    if [[ $stage == "wrf1" ]]; then
    # Run wrf for coarse domain
      namelist_file=${work_dir}/namelists/namelist.input.wrf.${case_name}.${test_name}.ndown_pt1
    elif [[ $stage == "wrf2" ]]; then
    # Run wrf for fine domain
      namelist_file=${work_dir}/namelists/namelist.input.wrf.${case_name}.${test_name}.ndown_pt2
      # overwrite d01 with d02 output from ndown
      mv wrfinput_d02 wrfinput_d01
      mv wrfbdy_d02 wrfbdy_d01
      # Delete symbolic links to wrf_coarse output
      /bin/rm -f wrfout_d*
    elif [[ $stage == "wrfrst" ]]; then
    # In case of restart, grab new copy of corresponding namelist
      namelist_file=${work_dir}/namelists/namelist.input.wrf.${case_name}.${test_name}.ndown_pt2
    fi
    /bin/cp $namelist_file ./namelist.input

    # Restarts/mech-denial tests: link to restart files, BCs from restart-base
    if [[ $stage == "wrfrst" ]]; then
      sed -i "s/restart.*/restart = .true.,/" namelist.input
      MM=`echo $restart_t_stamp | cut -d'-' -f 2`
      DD=`echo $restart_t_stamp | cut -d'-' -f 3 | cut -d'_' -f 1`
      HH=`echo $restart_t_stamp | cut -d'_' -f 2 | cut -d':' -f 1`
      sed -i "s/restart.*/restart = .true.,/" namelist.input
      sed -i "s/start_month.*/start_day = ${MM},/" namelist.input
      sed -i "s/start_day.*/start_day = ${DD},/" namelist.input
      sed -i "s/start_hour.*/start_hour = ${HH},/" namelist.input
    fi
    if [[ ${test_name} == *'crf'* ]] || [[ ${test_name} == *'STRAT'* ]]; then
      ln -sf "$memb_dir/${restart_base}/wrfrst_d01_${restart_t_stamp}" .
      ln -sf "$memb_dir/${restart_base}/wrfrst_d02_${restart_t_stamp}" .
      ln -sf "$memb_dir/${restart_base}/wrfbdy_d01" .
      ln -sf "$memb_dir/${restart_base}/wrflowinp_d01" .
      ln -sf "$memb_dir/${restart_base}/wrflowinp_d02" .
    fi

    # Create WRF batch script from header base script
    cat $work_dir/run_scripts/header_${system}.txt > batch_wrf_${test_name}.job

# Append the below lines to batch file [remove indentation]
echo "

source bashrc_wrf

# Run WRF
${mpi_command} ./wrf.exe

mkdir -p rsl_out
mv rsl.* rsl_out/
" >> batch_wrf_${test_name}.job

    # Modify batch script by replacing placeholders
    sed -i "s/PROJECT/${project_code}/g" batch_wrf_${test_name}.job
    sed -i "s/QUEUE/${queue}/g" batch_wrf_${test_name}.job
    sed -i "s/JOBNAME/${jobname}/g" batch_wrf_${test_name}.job
    sed -i "s/EMM/${em}/g" batch_wrf_${test_name}.job
    sed -i "s/NODELINE/${node_line}/g" batch_wrf_${test_name}.job
    # if [[ ${system} == 'oscer' ]]; then
    #   sed -i "s/NNODES/${nnodes}/g" batch_wrf_${test_name}.job
    # fi
    sed -i "s/TIMSTR/${run_time}/g" batch_wrf_${test_name}.job

    # Submit WRF job
    # if [[ `grep SUCCESS rsl.error.0000 | wc -l` -eq 0 ]] then
      ${submit_command} batch_wrf_${test_name}.job > submit_wrf_out.txt
    # fi
    # tail submit_wrf_out.txt

  fi

  done # Ensemble member loop

exit
