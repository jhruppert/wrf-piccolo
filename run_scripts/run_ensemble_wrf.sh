#!/bin/bash
# 
# Bash script to prepare and execute WRF for ensemble simulations.
# 
# WPS needs to be run before this step to generate the initial and boundary
# condition "met_em*" files. For that step, use scripts "download_gefs.py" and
# then "batch_wps_gefs.job" located in the "run_scripts" directory.
# 
# This step expects there to be a directory structure as follows:
#   in your scratch directory <scratch>:
#      <scratch>/case_name/memb_01/met_em/ (etc. for each ensemble member)
#   where "case_name" is the name of the case (e.g., sept1-4) and "met_em/" contains
#   the "met_em*" files corresponding to the ensemble member.
# 
# James Ruppert
# 18 Nov 2024

###################################################
# Settings
###################################################

# Selection to run REAL (initial & boundary conditions) or WRF
run_type="real"
run_type="wrf"

# Select case name
case_name="sept1-4"

# Select test (e.g., "ctl", "ncrf", etc.)
test_name="ctl"

# Number of ensemble members
nens=5

# WRF simulation details
  jobname="${case_name}_${test_name}"
  # Restart
    irestart=0
#    run_time='04:00' # HH:MM

###################################################
# Supercomputer environment-specific settings

# Current ongoing projects (totals as of as of 11/16/24):
#  - UFSU0031: Allison's exploratory allocation (remaining: 465k)
#  - UOKL0053: James's PICCOLO large allocation (20M)
#  - UOKL0049: James's TC-CRF large allocation (22M)
#  - UOKL0056: Frederick's TC small allocation (500k)

  system='derecho'
  if [[ ${system} == 'derecho' ]]; then
    queue="main"
    if [[ $run_type == "real" ]]; then
      bigN=10
    elif [[ $run_type == "wrf" ]]; then
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
    if [[ $run_type == "real" ]]; then
      run_time='02:00' # HH:MM Job run time
    elif [[ $run_type == "wrf" ]]; then
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
# Start parent loop for each ensemble member
###################################################

cd $ensemb_dir

# for em in 0{1..${nens}}; do # Ensemble member
for em in 01; do # Ensemble member

  # Create directory tree (e.g., as "memb_02/ctl/wrf") for ensemble member
  memb_dir="$ensemb_dir/memb_${em}"
  mkdir -p $memb_dir # -p ignores if directory already exists
  test_dir=$memb_dir/$test_name
  mkdir -p $test_dir
  mkdir -p $test_dir/wrf
  cd $test_dir/wrf

  echo "Running: $test_dir"


###################################################
# First step: REAL

  if [[ $run_type == "real" ]]; then

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
    namelist_file=${work_dir}/namelists/namelist.input.wrf.${case_name}.${test_name}
    /bin/cp $namelist_file ./namelist.input

    # Create REAL batch script from header base script
    cat $work_dir/run_scripts/header_${system}.txt > batch_real.job

# Append the below lines to batch file [remove indentation]
echo "

source bashrc_wrf

# Run REAL
${mpi_command} ./real.exe
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
    #   ${submit} batch_real.job > submit_real_out.txt
    # fi

###################################################
# Second step: WRF

  elif [[ $run_type == "wrf" ]]; then

  #  JOBID=$(grep Submitted submit_real_out.txt | cut -d' ' -f 4)

    # Remove and replace for specific test
    /bin/rm -f namelist.input

    # Delete text-out from REAL
    /bin/rm -f rsl.*
    /bin/rm -f namelist.output

    # In case of restart, grab new copy of corresponding namelist
    if [ ${irestart} -eq 1 ]; then
      namelist_file=${work_dir}/namelists/namelist.input.wrf.${case_name}.${test_name}.restart
    else
      namelist_file=${work_dir}/namelists/namelist.input.wrf.${case_name}.${test_name}
    fi
    /bin/cp $namelist_file ./namelist.input

    # Modify NAMELIST for nproc specs
    # sed -i '/nproc_x/c\ nproc_x = 20,' namelist.input
    # sed -i '/nproc_y/c\ nproc_y = 19,' namelist.input

    # Create symbolic links to restart files and BCs from restart-base for mechanism-denial test
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

# mkdir -p ../text_out
# mkdir -p ../post
# mkdir -p ../post/d01
# mkdir -p ../post/d02
# mv wrfout* wrfrst* ../
# mv namelist.out* out_wrf.* ../text_out/
# cp namelist.input ../text_out/
# if [[ ${test_name} == 'ctl' ]]; then
#   mv wrfinput* wrfbdy* wrflow* ../
# fi
" >> batch_wrf_${test_name}.job

    # Modify batch script by replacing placeholders
    sed -i "s/PROJECT/${project_code}/g" batch_real.job
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
      # ${submit} batch_wrf_${test_name}.job > submit_wrf_out.txt
    # fi
    # tail submit_wrf_out.txt

  fi

  cd ..

  done # Ensemble member loop

exit
