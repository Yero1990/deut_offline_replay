#!/bin/bash


#--------------------------------------------
# Author: C. Yero
# Date: Oct 20, 2022
# email: cyero@jlab.org, cyero002@gmail.com
#--------------------------------------------

# Brief: This shell script replays the raw (.dat) files 
#        and outputs the raw (.root) files. The leaf
#        variables to be written are specified under 
#        DEF-files/deut_{replay_type}.def and may be modified.
#        Where {replay_type} can be either production where only
#        necessary leaf variables are written for doing physics analysis,
#        or calibration (i.e., "optics", "hodcalib", "dccalib", "calcalib", "scalers",
#        "reftime", "timewin")

# What type of input is {replay_type} ?
# Answer: {replay_type} is an input based on the suffix of a symbolic link made to this shell script.
# For examlpe, if a symbolic link is made:  ln -sf UTILS_DEUT/offline_scripts/replay_deut.sh replay_deut_suffix.sh
# then the {replay_type} becomes "suffix".  This way, there exists ONLY this shell script, from which different
# shell scripts can be linked to, and based on the "suffix", then a different part of this shell script would be executed

# example 1: create a production replay script (replay raw data file with all leaf variables)
# ln -sf  UTILS_DEUT/offline_scripts/replay_deut.sh replay_deut_prod.sh

# example 2:  create calibration replay script (replay raw data file with only leaf variabled relevant to calibration)
# ln -sf  UTILS_DEUT/offline_scripts/replay_deut.sh replay_deut_hodcalib.sh

# Which replay type are we doing? physics analysis ("prod"), or calibration ("hodcalib", "dccalib", "calcalib", "scalers", "reftime", "timewin")
replay_type=${0##*_}
replay_type=${replay_type%%.sh}     

HCREPLAY="/work/hallc/c-deuteron/$USER/deut_offline_replay"
echo "HCREPLAY=${HCREPLAY}"

# change to top-level directory
echo ""
echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
echo ""
echo "changing to the top-level directory . . ."
echo "cd ${HCREPLAY}"
cd ${HCREPLAY}
echo ""
echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
echo ""

# replay script
replay_script="${HCREPLAY}/SCRIPTS/COIN/PRODUCTION/replay_deut.C" 


# ==========================
# replay production
# ==========================
if [ "${replay_type}" = "prod" ]; then
    
    # Display help output if no argumenr specified
    if [  $# -eq 0 ]; then
	echo "" 
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
	echo ""
	echo "Brief: This shell script replays the raw (.dat) files "
	echo "       and outputs the raw (.root) files. The leaf "
	echo "       variables to be written are specified under "
	echo "       DEF-files/deut_${replay_type}.def and may be modified. "
	echo ""
	echo "-----------------------------------------"
	echo "Usage 1):  ./replay_deut_${replay_type}.sh run evt "
	echo "-----------------------------------------"
	echo ""
	echo "run: run number"
	echo ""
	echo "evt: event number; defaults to -1 (all events)"
	echo "      if no argument is given"
	echo ""
	echo "example 1: ./replay_deut_${replay_type}.sh 3288 100000"
	echo ""
	echo "------------------------------------------------"
	echo "Usage 2):  ./replay_deut_${replay_type}.sh target kin evt "
	echo "------------------------------------------------"
	echo ""
	echo "target:  dummy, h2, d2, c12 "
	echo "target runlist is read from: UTILS_DEUT/runlist/<target>_<kin>.txt"
	echo ""
	echo "kin: singles, coin "
	echo ""
	echo "evt: event number defaults to -1 (all events), "
	echo "unless explicitly specified as 3rd argument"
	echo ""
	echo "example 2: ./replay_deut_${replay_type}.sh d2 deep"
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 
	echo "" 
	exit 0  
    fi


    run=$1
    evt=$2 

    # check if 1st argument is an integer (i.e., run number, else attempt to read from runlist)
    if [ $run -eq $run ]; then
	
	
	# check if event number is specified
	if [ -z $evt ]; then
	    evt=-1
	    echo "No event number spedified, defaulting to evt=${evt} (all events)"
	fi
	
	# hcana command
	run_hcana="./hcana -q \"${replay_script}(${run}, ${evt}, \\\"${replay_type}\\\")\""
	

	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	echo "" 
	date
	echo ""
	echo ""
	echo "Running HCANA Deut Replay on the run ${run}:"
	echo " -> SCRIPT:  ${replay_script}"
	echo " -> RUN:     ${run}"
	echo " -> NEVENTS: ${evt}"
	echo " -> COMMAND: ${run_hcana}"
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	
	sleep 2
	eval ${run_hcana}
        
	
    else
	
	# order of arguments now is as follows: target, kin, evt (run number excluded)
	# read 1st and 2nd arguments
	target=$1
	kin=$2
	
	evt=$3
	# check if event number is specified
	if [ -z $evt ]; then
	    $evt=-1
	    echo "No event number spedified, defaulting to evt=${evt} (all events)"
	fi
	
	# runlist
	filename="UTILS_DEUT/runlist/${target}_${kin}.txt"
	
	for run in $(cat $filename) ; do    
	    
	    # hcana command
	    run_hcana="./hcana -q \"${replay_script}(${run}, ${evt}, \\\"${replay_type}\\\")\""
	    
	    {
		echo ""
		echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
		echo "" 
		date
		echo ""
		echo ""
		echo "Running HCANA Deut Replay on the run ${runNum}:"
		echo " -> RUNLIST: ${filename}"
		echo " -> SCRIPT:  ${replay_script}"
		echo " -> RUN:     ${run}"
		echo " -> NEVENTS: ${evt}"
		echo " -> COMMAND: ${run_hcana}"
		echo ""
		echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
		
		sleep 2
		eval ${run_hcana}
	    }
	    
	    
	done
    fi
    

else
    
    # {replay_type} for calibration may be one of these:  ("hodcalib", "dccalib", "calcalib", "scalers", "reftime", "timewin" or "checkCalib")
    # Display help output if no argumenr specified
    if [  $# -eq 0 ]; then
	echo "" 
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
	echo ""
	echo "Brief: This shell script replays the raw (.dat) files "
	echo "       and outputs the raw (.root) files for calibration"
	echo "       purposes. The leaf variables to be written are "
	echo "       specified under: DEF-files/deut_${replay_type}.def "
	echo "       and may be modified. "
	echo ""
	echo "-----------------------------------------"
	echo "Usage 1):  ./replay_deut_${replay_type}.sh run evt "
	echo "-----------------------------------------"
	echo ""
	echo "run: run number"
	echo ""
	echo "evt: event number defaults to -1 (all events), "
	echo "unless explicitly specified as 3rd argument"
	echo ""
	echo "example 1: ./replay_deut_${replay_type}.sh 3288 100000"
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 
	echo "" 
	exit 0  

	
    else

	# user input (run number)
	run=$1
	# user input (event number)
	evt=$2
	
	# check if event number is specified
	if [ -z $evt ]; then
	    evt=-1
	    echo "No event number spedified, defaulting to evt=${evt} (all events)"
	fi
	
	# hcana command                                                                       
        run_hcana="./hcana -q \"${replay_script}(${run}, ${evt}, \\\"${replay_type}\\\")\""   

	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	echo "" 
	date
	echo ""
	echo ""
	echo "Running HCANA Deut Replay on the run ${run}:"
	echo " -> SCRIPT:  ${replay_script}"
	echo " -> RUN:     ${run}"
	echo " -> NEVENTS: ${evt}"
	echo " -> COMMAND: ${run_hcana}"
	echo ""
	echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	
	sleep 2
	eval ${run_hcana}	   
	
    fi
    
fi
