#!/bin/bash

# Author: C. Yero
# Date: Sep 25, 2022

# Brief: shell script for automating hodoscopes calibration.

#user input
runNum=$1    # run number
evtNum=$2    # event number

if [ -z "$1" ] || [ -z "$2" ]; then 
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo "" 
    echo "Usage:  ./do_hodo_calib.sh <run_number> <evt_number>"  
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    
    exit 0
fi

daq_mode="coin"
debug=0
  
# Which analysis file type are we doing? "hod_calib"
ana_type="hodcalib"

# this rootfile name pattern assumes pattern defined in replay_cafe.C script (please do NOT modify replay script) 
filename="ROOTfiles/${ana_type}/cafe_replay_${ana_type}_${runNum}_${evtNum}.root"

# re-define other file names for later use (in re-naming)
filename_uncalib="ROOTfiles/${ana_type}/cafe_replay_hodoUnCalib_${runNum}_${evtNum}.root"  # uncalibrated file (initially replayed)
filename_twcalib="ROOTfiles/${ana_type}/cafe_replay_hodoTWCalib_${runNum}_${evtNum}.root"  # intermediate calibrated (time-walk calibrated)
filename_calib="ROOTfiles/${ana_type}/cafe_replay_hodoCalib_${runNum}_${evtNum}.root"      # fully calibrated (fit hodo matrix)


replay_script="SCRIPTS/COIN/PRODUCTION/replay_cafe.C"

analysis_script_1="timeWalkHistos.C"
analysis_script_2="timeWalkCalib.C"
analysis_script_3="fitHodoCalib.C"

# run hcana command
runHcana="./hcana -q \"${replay_script}(${runNum}, ${evtNum}, \\\"${ana_type}\\\")\""

# change to top direcotry and run analyzer to produce specified ROOTfile
cd ../../
eval $runHcana

# rename file
mv "${filename}" "${filename_uncalib}"

# change back to original directory and execute analysis script
cd CALIBRATION/shms_hodo_calib

# run hodoscope calibration analysis script 1 ( generates a root file with histogram objects for time-walk calibration)
runAna="root -l -q -b \"${analysis_script_1}(\\\"${filename_uncalib}\\\", ${runNum}, \\\"${daq_mode}\\\")\""
eval $runAna

# run hodoscope calibration analysis script 2 ( performs time-walk correction on the histogram object file generated by script 1)
# this code generates a param file:  ./phodo_TWcalib_${runNum}.param
runAna="root -l -q -b \"${analysis_script_2}(${runNum})\""
eval $runAna

# move param file to appropiate location and make symbolic link
mv "phodo_TWcalib_${runNum}.param" "../../PARAM/SHMS/HODO/cafe2022/calib/"
cd ../../PARAM/SHMS/HODO/
ln -sf "cafe2022/calib/phodo_TWcalib_${runNUM}.param phodo_TWcalib.param" 

# replay data a 2nd time (with updated phodo_TWcalib)
cd ../../../
eval $runHcana

# rename file
mv "${filename}" "${filename_twcalib}"

# change back to original directory and execute analysis script

cd CALIBRATION/shms_hodo_calib

# run hodoscope calibration analysis script 3 (Determine the the effective propagation speeed in the paddle, 
#the time difference between the positive and negative PMTs and then the relative time difference of all paddles 
#compared to paddle 7 in plane S1X. 

# this code generates a param file:  ./phodo_Vpcalib_${runNum}.param
runAna="root -l -q -b \"${analysis_script_3}(\\\"${filename_twcalib}\\\", ${runNum})\""
eval $runAna

# move param file to appropiate location and make symbolic link
mv "phodo_Vpcalib_${runNum}.param" "../../PARAM/SHMS/HODO/cafe2022/calib/"
cd ../../PARAM/SHMS/HODO/
ln -sf "cafe2022/calib/phodo_Vpcalib_${runNUM}.param phodo_Vpcalib.param" 

# replay data a 3rd time (with updated phodo_TWcalib and phodo_Vpcalib)
# but this time, reset evtNum to only 100k, since we just want to check calibration was done
evtNum=100000 
cd ../../../
eval $runHcana

# rename file
mv "${filename}" "${filename_calib}"
