#! /bin/bash 

'''
description:
-------------
This shell script is used to:

1) replay raw data file  generate a .root file which contains the relevant
   leaf variables for checking the calibration of each detector
   (./replay_deut_checkCalib.sh)

2) analyze the raw .root file (checkCalib.C) to generate histogram objects
   to check whether the relevant detectors are properly calibrated

The shell script takes as input a list of runs the user would like to check

'''


# USER SET RUN LIST (here are some examples)
#--------------------
#filename='calib_test_runs.txt' 
#filename='h_singles_aug08.txt'
filename='optics_aug08.txt'
#filename='hms_xem.txt'

for run in $(cat $filename) ; do  
    
    #run=16036
    evt=200000
    
    # generic filename to be read by analzye_script
    root_file="ROOTfiles/deut_replay_checkCalib_${run}_${evt}.root"
    
    # analys_script arguments 
    hms_pid="p"
    shms_pid="e"
    
    # define scripts
    replay_script="replay_deut_checkCalib.sh"
    analyze_script="checkCalib.C"
    
    # define commands to run the scripts
    replay_cmd="./${replay_script} ${run} ${evt}"
    analyze_cmd="root -l -q -b \"${analyze_script}(\\\"${root_file}\\\", ${run}, \\\"${hms_pid}\\\", \\\"${shms_pid}\\\")\""  
    
    # ecevute the script
    echo "${replay_cmd}"
    eval ${replay_cmd}
    echo "${analyze_cmd}" 
    eval ${analyze_cmd}

done
