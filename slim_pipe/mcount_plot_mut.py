
from tools.mcounter_tools import (
    read_vcf_allel, ind_assignment_scatter_v1, MC_sample_matrix_v1,
    heatmap_v2, read_args
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
sample_sim= 0

data, data_freqs = MC_sample_matrix_v1(min_size= min_size, samp= sampling, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,
                          muted_dir= muted_dir, diffs= diffs,sample_sim= sample_sim,
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
                                     muted_dir= muted_dir, data_freqs= data_freqs)


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
grid_diffs= [count_data[s]['diffs'] for s in avail_sub]
## extract statistics per mutation.
mut_grid= {}
mut_prop= {}
mut_diffs= {}

for row in range(grid_shape[0]):
    for col in range(grid_shape[1]):
        mut= grid_labels[row,col]
        mut_grid[mut]= []
        mut_prop[mut]= []
        mut_diffs[mut]= []

        for idx in range(len(avail_sub)):
            mut_grid[mut].append(grids[idx][row,col])
            mut_prop[mut].append(props[idx][row,col])
            mut_diffs[mut].append(grid_diffs[idx][row,col]**2)

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



##################################################### KMER
## dictionary to hold values across mutation contexts. 
compound_kmer= {
        y: {
                g: {
                    z: [] for z in list(set(pop_proportions))
                } for g in batch_dict.keys()
            } for y in ['pval','prop','diffs']
        }


## plot first for every mutation context.
for kmer in mut_grid.keys():
    fig_kmer= fig_dir + '/' + kmer
    os.makedirs(fig_kmer, exist_ok=True)
    fig_kmer= fig_kmer + '/'
    
    plot_data= {
        'pval': mut_grid[kmer],
        'prop': mut_prop[kmer],
        'diffs': mut_diffs[kmer]
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
            #plt.ylim(0,1.5)
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
        #plt.ylim(0,1.5)
        plt.ylabel(ylab)
        plt.title('combined stats')
        plt.legend()

        plt.savefig(fig_kmer + 'combined_{}_{}.png'.format(kmer,strata),bbox_inches='tight')
        plt.close()


##########################################
########################################## GEN

xlab= 'relative sampling'
ylab= 'mean matrix p-val'

for i in batch_dict.keys(): 
    plt.figure(figsize=(20, 10))

    for kmer in mut_grid.keys():

        plot_data= {
            #'pval': mut_grid[kmer],
            #'prop': mut_prop[kmer],
            'diffs': mut_diffs[kmer]
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
            
            plt.errorbar(x,y,yerr=error)

    plt.xlabel(xlab)
    #plt.ylim(0,1.5)
    plt.ylabel(ylab)
    plt.title('combined stats')

    plt.savefig(fig_dir + 'combined_{}_{}.png'.format(i,strata),bbox_inches='tight')
    plt.close()

####################################################
#################################################### grid SSD

xlab= 'relative sampling'
ylab= 'mean matrix p-val'

grid_whole= {}

for i in batch_dict.keys(): 
    plt.figure(figsize=(20, 10))

    xprep= [pop_proportions[x] for x in batch_dict[i]]
    xprep= {
         z: [x for x in range(len(xprep)) if xprep[x] == z] for z in list(set(xprep))
    }


    batch_grids= [grid_diffs[x] for x in batch_dict[i]]
    y_prep= {
        z: [batch_grids[x] for x in xprep[z]] for z in xprep.keys()
    }

    y_prep= {
        z: [np.sum(x**2) for x in y_prep[z]] for z in y_prep.keys()
    }


    surface= sorted(xprep.keys())
    y= [np.mean(y_prep[x]) for x in surface]
    error= [np.std(y_prep[x]) for x in surface]

    grid_whole[i]= [surface,y,error]

    plt.errorbar(surface,y,yerr=error)    

    plt.xlabel(xlab)
    #plt.ylim(0,1.5)
    plt.ylabel(ylab)
    plt.title('grid SSD / sample proportion')

    plt.savefig(fig_dir + 'gridSSD_{}.png'.format(i),bbox_inches='tight')
    plt.close()

plt.figure(figsize=(20, 10))

for i in grid_whole.keys():
    plt.errorbar(grid_whole[i][0],grid_whole[i][1],yerr=grid_whole[i][2],label= i)    

plt.xlabel(xlab)
#plt.ylim(0,1.5)
plt.ylabel(ylab)
plt.title('grid SSD / sample proportion')
plt.legend()
plt.savefig(fig_dir + 'gridSSD_combined_.png',bbox_inches='tight')
plt.close()

############################################################# 
############################################################# STRATA


for strata in ['prop','pval','diffs']:

    plt.figure(figsize=(20,10))
    
    for batch in batch_dict.keys():
        global_x= sorted(list(compound_kmer[strata][batch].keys()))
        global_y= [np.mean(compound_kmer[strata][batch][x]) for x in global_x]
        global_error= [np.std(compound_kmer[strata][batch][x]) for x in global_x]
        
        plt.errorbar(global_x,global_y,yerr=global_error,label= batch)
    plt.xlabel(xlab)
    plt.ylim(0,1.5)
    plt.ylabel(ylab)
    plt.title('combined stats')

    plt.legend()
    plt.savefig(fig_dir + 'combined_{}_{}.png'.format('kmers',strata),bbox_inches='tight')
    plt.close()


    plt.figure(figsize=(20,10))
    
    for batch in batch_dict.keys():
        global_x= sorted(list(compound_kmer[strata][batch].keys()))
        
        global_error= [np.std(compound_kmer[strata][batch][x]) for x in global_x]
        plt.plot(global_x,global_error,label= batch)
    
    plt.xlabel(xlab)
    #plt.ylim(0,1.5)
    plt.ylabel('variance')
    plt.title('combined stats')
    
    plt.legend()
    plt.savefig(fig_dir + 'combined_{}_{}.png'.format('variance',strata),bbox_inches='tight')
    plt.close()



##############################################################
############################################################## SFS


fig_dir= 'Figures/kmers/'
mu= 1e-8
N_inds= 1000
xlab= 'sampling proportion'
ylab= 'sum of sqared differences.'

batch_sfs= {}

for i in batch_dict.keys():
    sims= [avail_sub[x] for x in batch_dict[i]]
    print(len(sims))
    counts= []
    props= [pop_proportions[x] for x in batch_dict[i]] 

    for sim_idx in sims:
        freqs_dict= count_data[sim_idx]['freqs']
        sfs= []
        sizes_sim= count_data[sim_idx]['sizes']
        N_inds= min(sizes_sim)
        for pop in freqs_dict.keys():

            #print(gen_time)
            freqs= freqs_dict[pop]

            sizeN= [x[1] for x in freqs if x[1] > 0]
            freqs= np.repeat([x[0] for x in freqs if x[1] > 0],sizeN) / sizes_sim[pop]
            
            bin_count= np.histogram(freqs,bins= N_inds,range= [0,1])[0]
            bin_count= bin_count / sum(bin_count)
            sfs.append(bin_count)

        dist_vec= sfs[0] - sfs[1] #/ np.sum(indian + x)
        dist_vec= dist_vec**2
        dist_vec= np.sqrt(np.sum(dist_vec)) / N_inds
        
        counts.append(dist_vec)

    props_dict= {
        z: [counts[x] for x in range(len(props)) if props[x] == z] for z in list(set(props))
    }

    props_sorted= sorted(props_dict.keys())
    props_means= [np.mean(props_dict[x]) for x in props_sorted]
    props_std= [np.std(props_dict[x]) for x in props_sorted]

    batch_sfs[i]= [props_sorted,props_means,props_std]

    plt.figure(figsize=(20,10))

    plt.xlabel(xlab)
    #plt.ylim(0,1.03)
    plt.ylabel(ylab)
    plt.title('sfs_differences.')
    
    plt.errorbar(props_sorted,props_means,yerr=props_std)
    
    plt.savefig(fig_dir + 'SFS_{}.png'.format(i),bbox_inches='tight')
    plt.close()


plt.figure(figsize=(15,10))

plt.xlabel(xlab)
#plt.ylim(0,1.03)
plt.ylabel(ylab)
plt.title('sfs_differences')

for i in batch_sfs.keys():
    plt.errorbar(batch_sfs[i][0],batch_sfs[i][1],yerr=batch_sfs[i][2],label= i)

plt.legend()
plt.savefig(fig_dir + 'SFS_combined.png',bbox_inches='tight')
plt.close() 






