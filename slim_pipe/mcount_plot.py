
from tools.mcounter_tools import (
    read_vcf_allel, ind_assignment_scatter_v1, MC_sample_matrix_v1,
    heatmap_v2
)

#from tools.SLiM_pipe_tools import mutation_counter_launch
import re
import pandas as pd
import os
import numpy as np

## directories
main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
muted_dir= main_dir + 'mutation_counter/data/mutation_count/'
sims_dir= main_dir + 'mutation_counter/data/sims_dem/'
diffs= False

mutlog= 'toMut.log'
min_size= 70
sampling= [5,100,10]

data = MC_sample_matrix_v1(min_size= min_size, samp= sampling, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,
                          muted_dir= muted_dir, diffs= diffs,
                       exclude= False)


### new - pair reference sims and subsetted populations.
### extract kmer comparisons (proportions or pvals) using heatmap_v2.
### make function. 

from tools.mcounter_tools import mcounter_deploy

p_value= 1e-5
test_m= 'fisher'
individually= False
exclude= False
frequency_range= [0,1]
extract= 'pval'

pop_asso, count_data= mcounter_deploy(data,p_value= p_value, test_m= test_m, individually= individually,
                                        exclude= exclude, frequency_range= frequency_range, extract= extract,
                                     muted_dir= muted_dir)

##############
############## process grids


sims_dir= main_dir + 'mutation_counter/data/sims/'

available= list(count_data.keys())
subsamp= len(count_data)
avail_sub= np.random.choice(available,subsamp,replace= False)

### 1. extract grids
grids= [count_data[s]['grids'] for s in avail_sub]
grid_shape= grids[0].shape()
grid_total= np.prod(grid_shape)


props= [count_data[s]['prop'] for s in avail_sub]

#grids= list(it.chain(*grids))

## mask infinite values and compute std.
#grids= [np.ma.masked_where(a == np.inf, a) for a in grids]
grid_mean= [np.mean(x) for x in grids] 
grid_std= [np.std(x) for x in grids]
prop_mean= [np.mean(x) for x in props]

### 2. calculate proportions across smulations
pop_proportions= [count_data[s]['sizes'][1] / count_data[s]['sizes'][0] for s in avail_sub]

### 3. batch names
batch_names= [count_data[s]['batch'] for s in avail_sub]
batch_dict= {
    z:[x for x in range(len(avail_sub)) if batch_names[x] == z] for z in list(set(batch_names))
}



###############
############### plots


fig_dir= 'Figures/'


plot_data= {
	'pval': grid_mean,
	'prop': prop_mean
}

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

for strata in plot_data.keys():

	ydata= plot_data[strata]

	xlab= 'relative sampling'
	ylab= 'mean matrix p-val'

	d=0

	colors= ['ro', 'bs', 'g^']

	append_all= []
	for i in batch_dict.keys():
	    
	    y= [ydata[x] for x in batch_dict[i]]
	    x= [pop_proportions[x] for x in batch_dict[i]]
	    
	    plt.figure(figsize=(10, 10))
	    
	    plt.plot(x,y,colors[d])
	    plt.xlabel(xlab + ' {} comparisons'.format(len(batch_dict[i])))
	    plt.ylim(0,1.03)
	    plt.ylabel(ylab)
	    plt.title(i)
	    
	    plt.savefig(fig_dir + '{}_{}.png'.format(i, strata),bbox_inches='tight')
	    
	    append_all.extend([x,y,colors[d]])
	    #plt.show()
	    
	    d += 1


	plt.figure(figsize=(10, 10))

	plt.plot(*append_all)
	plt.xlabel(xlab)
	plt.ylim(0,1.03)
	plt.ylabel(ylab)
	plt.title('combined stats')

	plt.savefig(fig_dir + 'combined_{}.png'.format(strata),bbox_inches='tight')


