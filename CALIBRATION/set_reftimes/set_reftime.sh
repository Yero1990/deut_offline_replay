#!/bin/bash

#user input
runNum=$1    # run number
evtNum=$2    # event number
replay=$3    # replay data (default no replay)
singles=$4   # is it a singles run?   
#daq_mode=$4  # shms(singles) or coin (default coin) always coin for deut23

if [ -z "$1" ] || [ -z "$2" ]; then 
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=::=:=:"
    echo "" 
    echo "Usage:  ./set_reftime.sh <run_number> <evt_number> <replay(opt)> <singles?(opt)>"  
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=::=:=:"
    
    exit 0
fi

daq_mode="coin"
set_refTimes=1
debug=0
#singles=1

if [ -z "$4" ]; then
    singles=0
    echo ""
    echo "Analyzing coincidence run ${runNum}"
    echo ""
else
    singles=1
    echo ""
    echo "Analyzing singles run ${runNum}"
    echo ""
fi        

# Which analysis file type are we doing? "reftime" or "timewin"?
ana_type="reftime"

# this rootfile name pattern assumes pattern defined in replay_deut.C script (please do NOT modify replay script) 
filename="../../ROOTfiles/${ana_type}/deut_replay_${ana_type}_${runNum}_${evtNum}.root"

replay_script="SCRIPTS/COIN/PRODUCTION/replay_deut.C"
analysis_script="scripts/set_reftimes.C"

runHcana="./hcana -q \"${replay_script}(${runNum}, ${evtNum}, \\\"${ana_type}\\\")\""
runAna="root -l -q -b \"${analysis_script}(\\\"${filename}\\\", ${runNum}, \\\"${daq_mode}\\\", ${set_refTimes}, ${debug},${singles})\""  

# change to top direcotry and run analyzer to produce specified ROOTfile
if [ "${replay}" == "replay" ]; then 
    cd ../../
    eval $runHcana
    # change back to original directory and execute analysis script
    cd CALIBRATION/set_reftimes
    eval $runAna
else
    eval $runAna
fi
