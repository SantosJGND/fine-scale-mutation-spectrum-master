
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
bases= 'ATCG'
ksize= 3


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


############## mutation grids
############## mutation grids
from functools import reduce  # forward compatibility for Python 3
import operator

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
)

bases= 'ATCG'
ksize= 3

mutations= get_mutations(bases= bases,ksize= ksize)
kmers, kmer_idx= kmer_comp_index(mutations)

mut_lib= kmer_mut_index(mutations)
labels= [kmer_idx[x][0] for x in sorted(kmer_idx.keys())]
grid_labels= np.array(labels).reshape(24,4)
list_labels= grid_labels.reshape(1,np.prod(grid_labels.shape))[0]
##############
############## process grids

sims_dir= main_dir + 'mutation_counter/data/sims/'

available= list(count_data.keys())
subsamp= len(count_data)
avail_sub= np.random.choice(available,subsamp,replace= False)

### 1. extract grids
grids= [count_data[s]['grids'] for s in avail_sub]
grid_shape= grids[0].shape
grid_total= np.prod(grid_shape)

props= [count_data[s]['prop'] for s in avail_sub]

## extract statistics per mutation.
mut_grid= {}
mut_prop= {}

for row in range(grid_shape[0]):
    for col in range(grid_shape[1]):
        mut= grid_labels[row,col]
        mut_grid[mut]= []
        mut_prop[mut]= []

        for idx in range(len(avail_sub)):
            mut_grid[mut].append(grids[idx][row,col])
            mut_prop[mut].append(props[idx][row,col])

## mask infinite values and compute std.
#grids= [np.ma.masked_where(a == np.inf, a) for a in grids]
#grid_mean= [np.mean(x) for x in grids] 
#grid_std= [np.std(x) for x in grids]
#prop_mean= [np.mean(x) for x in props]

### 2. calculate proportions across smulations
pop_proportions= [count_data[s]['sizes'][1] / count_data[s]['sizes'][0] for s in avail_sub]
pop_proportions= [round(x,3) for x in pop_proportions]
### 3. batch names
batch_names= [count_data[s]['batch'] for s in avail_sub]
batch_dict= {
    z:[x for x in range(len(avail_sub)) if batch_names[x] == z] for z in list(set(batch_names))
}


###############
############### plots

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt


fig_dir= 'Figures/kmers'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'

## dictionary to hold values across mutation contexts. 
compound_kmer= {
        y: {
                g: {
                    z: [] for z in list(set(pop_proportions))
                } for g in batch_dict.keys()
            } for y in ['pval','prop']
        }

## plot first for every mutation context.
for kmer in mut_grid.keys():
    fig_kmer= fig_dir + '/' + kmer
    os.makedirs(fig_kmer, exist_ok=True)
    fig_kmer= fig_kmer + '/'
    
    plot_data= {
        'pval': mut_grid[kmer],
        'prop': mut_prop[kmer]
    }
    
    for strata in plot_data.keys():
        batch_hold= {}
        ydata= plot_data[strata]
        
        xlab= 'relative sampling'
        ylab= 'mean matrix p-val'

        d=0

        colors= ['ro', 'bs', 'g^']

        append_all= []
        for i in batch_dict.keys():

            xprep= [pop_proportions[x] for x in batch_dict[i]]
            #print(xprep)
            xprep= {
                z: [x for x in range(len(xprep)) if xprep[x] == z] for z in list(set(xprep))
            }
            

            yprep= [ydata[x] for x in batch_dict[i]]
            yprep= {
                z: [yprep[x] for x in xprep[z]] for z in xprep.keys()
            }
            for prop_z in yprep.keys():
                compound_kmer[strata][i][prop_z].extend(yprep[prop_z])

            x= sorted(xprep.keys())
            y= [np.mean(yprep[z]) for z in x]
            error= [np.std(yprep[z]) for z in x]

            batch_hold[i]= {
                'x': x,
                'y': y,
                'error': error
            }

            plt.figure(figsize=(10, 10))

            plt.errorbar(x,y,yerr=error)
            plt.xlabel(xlab + ' {} comparisons'.format(len(batch_dict[i])))
            plt.ylim(0,1.03)
            plt.ylabel(ylab)
            plt.title(i)

            plt.savefig(fig_kmer + '{}_{}_{}.png'.format(kmer,i, strata),bbox_inches='tight')
            plt.close()
            
            #append_all.extend([x,y,colors[d]])
            #plt.show()

            d += 1



        plt.figure(figsize=(10, 10))

        for i in batch_hold.keys():
            plt.errorbar(batch_hold[i]['x'],batch_hold[i]['y'],yerr=batch_hold[i]['error'],label= i)

        plt.xlabel(xlab)
        plt.ylim(0,1.03)
        plt.ylabel(ylab)
        plt.title('combined stats')
        
        plt.legend()
        plt.savefig(fig_kmer + 'combined_{}_{}.png'.format(kmer,strata),bbox_inches='tight')
        plt.close()


##########################################
##########################################



xlab= 'relative sampling'
ylab= 'mean matrix p-val'

for i in batch_dict.keys(): 
    plt.figure(figsize=(10, 10))

    for kmer in mut_grid.keys():

        plot_data= {
            #'pval': mut_grid[kmer],
            'prop': mut_prop[kmer]
        }

        for strata in plot_data.keys():
            ydata= plot_data[strata]

            d=0

            colors= ['ro', 'bs', 'g^']

            append_all= []
            xprep= [pop_proportions[x] for x in batch_dict[i]]
            #print(xprep)
            xprep= {
                z: [x for x in range(len(xprep)) if xprep[x] == z] for z in list(set(xprep))
            }


            yprep= [ydata[x] for x in batch_dict[i]]
            yprep= {
                z: [yprep[x] for x in xprep[z]] for z in xprep.keys()
            }
            for prop_z in yprep.keys():
                compound_kmer[strata][i][prop_z].extend(yprep[prop_z])

            x= sorted(xprep.keys())
            y= [np.mean(yprep[z]) for z in x]
            error= [np.std(yprep[z]) for z in x]

            batch_hold[i]= {
                'x': x,
                'y': y,
                'error': error
            }
            
            plt.errorbar(x,y,yerr=error,label= i)

    plt.xlabel(xlab)
    plt.ylim(0,1.03)
    plt.ylabel(ylab)
    plt.title('combined stats')

    #plt.legend()
    plt.savefig(fig_dir + 'combined_{}_{}.png'.format(i,strata),bbox_inches='tight')
    plt.close()



#############################################################
#############################################################


for strata in ['prop','pval']:

    plt.figure(figsize=(10,10))
    
    for batch in batch_dict.keys():
        global_x= sorted(list(compound_kmer[strata][batch].keys()))
        global_y= [np.mean(compound_kmer[strata][batch][x]) for x in global_x]
        global_error= [np.std(compound_kmer[strata][batch][x]) for x in global_x]
        
        plt.errorbar(global_x,global_y,yerr=global_error)
    plt.xlabel(xlab)
    plt.ylim(0,1.03)
    plt.ylabel(ylab)
    plt.title('combined stats')

    plt.legend()
    plt.savefig(fig_dir + 'combined_{}_{}.png'.format('kmers',strata),bbox_inches='tight')
    plt.close()


##############################################################
##############################################################