import os

run_nums = open("deut_offline_replay/run-list.dat","r")
spec = 'shms'       # shms, hms, or coin


for run in run_nums:
    r = int(run)
    tapeDIR = f'/mss/hallc/c-deuteron/raw/{spec}_all_{r}.dat'
    jcacheCMD = f'jcache get {tapeDIR} -e gvill@jlab.org -D 60'

    os.system(jcacheCMD)
run_nums.close()

