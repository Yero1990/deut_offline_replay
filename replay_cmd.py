
replay_DIR = './deut_offline_replay/'
replay_script = 'replay_deut_prod.sh'
nevents = '-1'

with open('runs_to_replay.txt','r') as readf:
   runs_to_replay =[]
   for line in readf:
        runs_to_replay.append(line.strip())
#print (runs_to_replay)

with open('replay_cmd.txt','w') as writef:
    for run in runs_to_replay:
        replay_cmd = replay_DIR + replay_script + ' ' + f'{run}' + ' ' + nevents + '\n'
        #print(replay_cmd)
        writef.write(replay_cmd)