
import numpy as np 
import os


def mutation_counter_launch(logfile,count_dir= './count/', 
                dir_launch= '..',main_dir= './', outlog= 'muted.log'):
    '''
    launch mutation counter.
    - read mut.log to know which have not been yet processed.
    - launch process_chromosomes.py using simulation name. 
    '''
    with open(logfile,'r') as fp:
        lines= fp.readlines()
    
    
    sims= [x.strip() for x in lines]
    sims= [x for x in sims if x]
    sims= [x.split()[0] for x in sims]
    chroms= [x.split('.')[0].split('C')[-1].strip('chr') for x in sims]
    
    job= 'python process_chromosomes.py -c {} -r {} -s {} -v {}_ -q {} -d {}'
    
    sims= [job.format(chroms[x],*[sims[x]]*4,dir_launch) for x in range(len(sims))]
    
    os.chdir(count_dir)
    for sim in sims:
        
        os.system(sim)
    
    os.chdir(main_dir)

    with open(outlog,'a') as fp:
        fp.write('\n' + ''.join(lines))

    open(logfile,'w').close()


#from tools.SLiM_pipe_tools import mutation_counter_launch

mutlog= 'regions.log'

## directories
main_dir= '/home/jgarc235/Sim_compare/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'

print('launch mutation counter.')
mutation_counter_launch(mutlog,count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir)

