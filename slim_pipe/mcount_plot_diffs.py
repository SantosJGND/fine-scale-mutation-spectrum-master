
from tools.mcounter_tools import (
    read_vcf_allel, ind_assignment_scatter_v1, MC_sample_matrix_v1,
    heatmap_v2, read_args
)

#from tools.SLiM_pipe_tools import mutation_counter_launch
import re
import pandas as pd
import os
import numpy as np
import itertools as it
import collections

def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


## directories
main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
muted_dir= main_dir + 'mutation_counter/data/mutation_count/'
sims_dir= main_dir + 'mutation_counter/data/sims/'
diffs= False

mutlog= 'toMut.log'
min_size= 70
sampling= [5,100,10]
bases= 'ATCG'
ksize= 3
sample_sim= 0
stepup= ""

data, data_freqs = MC_sample_matrix_v1(min_size= min_size, samp= sampling, stepup= stepup, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,
                          muted_dir= muted_dir, diffs= diffs,sample_sim= sample_sim,
                       exclude= False)


###

###

def run_stats(ref_sim,ref_pair,data,data_freqs= {}):
    '''
    co-factor function to md counter comparisons, deploy heatmap and calculate kmer proportion differences 
    between pairs of population.
    - ref pair: list of tuples. can't be dictionary because of repeated pops / reference tags. 
    '''
    batch= ref_sim.split('C')[0]
    sizes= [data[x[0]]['sizes'][x[1]] for x in ref_pair]
    #

    chromosomes= [ref_sim.split('.')[0].split('C')[1]]

    pop_counts= {
        g: data[g[0]]['counts'][g[1]] for g in ref_pair
    }

    num_variants= {
        g: data[g[0]]['Nvars'][g[1]] for g in ref_pair
    }

    ratio_grid, sig_cells= heatmap_v2(chromosomes,pop_counts,num_variants,
                                      {},frequency_range, exclude, p_value, muted_dir,tag= '',
                                      test= test_m,output= 'pval')
    
    pop_counts= {
        z: s / np.sum(s) for z,s in pop_counts.items()
    }

    grid_diffs= pop_counts[ref_pair[0]] - pop_counts[ref_pair[1]]

    comb_stats= {
        'grids': ratio_grid,
        'sigs': sig_cells,
        'sizes': sizes,
        'batch': batch,
        'diffs': grid_diffs
    }

    if data_freqs:
        comb_stats['freqs']= {
            ref_pair.index(x): data_freqs[x[0]][x[1]] for x in ref_pair
        }
    
    return comb_stats




def mcounter_deploy_v2(data,p_value= 1e-5, test_m= 'fisher', individually= False,
                            exclude= False, frequency_range= [0,1], data_freqs= {}, extract= 'pval',
                            muted_dir= '', tag_ref= '_ss'):
    '''
    Parse data dictionary.
        data: {sim: {counts:{pop:g}, Nvars:{pop:g}, sizes:{pop:g}}}
    i: use sim and pop IDs to create dictionary connecting original populations to 
    subset populations created using ind_assignment_scatter_v1.
    ii: for each pair of reference/subset populations, launch heatmapv2. return grid pvals or proportions,
    and proportion of mutations in subset population. allows for fisher or chi2 test for pval.
    - v2: compares sub pops to ref full pops other than its own
    '''
    
    avail= list(data.keys())
    ref_idx= [int(tag_ref in avail[x]) for x in range(len(avail) )]
    categ= {
        z: [x for x in range(len(avail)) if ref_idx[x] == z] for z in [0,1]
    }

    pop_asso= {avail[x]:recursively_default_dict() for x in categ[0]}

    for av in categ[1]:
        dat= [x for x in data[avail[av]]['counts'].keys() if tag_ref in x]
        ref_sim= avail[av].split(tag_ref)[0]
        ref_pop= [x.split('.')[0].strip(tag_ref) for x in dat]
        for p in range(len(dat)):
            pop_asso[ref_sim][ref_pop[p]][avail[av]]= dat[p]

    d= 0
    count_data= recursively_default_dict()

    for ref in pop_asso.keys():
        batch= ref.split('C')[0]
        
        for pop in pop_asso[ref].keys():
            for sub in pop_asso[ref][pop].keys():
                
                ref_pair= [(ref, pop),(sub, pop_asso[ref][pop][sub])]
                
                count_data[d]= run_stats(ref,ref_pair,data,data_freqs= data_freqs)
                
                count_data[d]['pop']= pop
                
                count_data[d]['other']= [] 
                
                
                for ref2 in pop_asso.keys():
                    for pop2 in pop_asso[ref2].keys():
                        if [ref,pop] == [ref2,pop2]:
                            continue
                        if ref2.split('C')[0] != batch: 
                            continue
                        ##
                        pop_dict= {
                            ref2: pop2,
                            sub: pop_asso[ref][pop][sub]
                        }
                        ref_pair= [(ref2, pop2),(sub, pop_asso[ref][pop][sub])]
                        
                        pair_stats= run_stats(ref,ref_pair,data,data_freqs= data_freqs)
                                                
                        count_data[d]['other'].append(pair_stats['diffs'])
                
                d += 1
    
    return pop_asso, count_data




###########
##########


from tools.mcounter_tools import mcounter_deploy

p_value= 1e-5
test_m= 'fisher'
individually= False
exclude= False
frequency_range= [0,1]
extract= 'pval'

pop_asso, count_data= mcounter_deploy_v2(data,p_value= p_value, test_m= test_m, individually= individually,
                                        exclude= exclude, frequency_range= frequency_range, extract= extract,
                                     muted_dir= muted_dir, data_freqs= data_freqs)


############## mutation grids
############## mutation grids
from functools import reduce  # forward compatibility for Python 3
import operator

from tools.fasta_utilities import (
    get_mutations, kmer_comp_index, kmer_mut_index
)

bases= 'ACGT'
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


grid_diffs= [count_data[s]['diffs'] for s in avail_sub]
## extract statistics per mutation.
mut_grid= {}
mut_diffs= {}

for row in range(grid_shape[0]):
    for col in range(grid_shape[1]):
        mut= grid_labels[row,col]
        mut_grid[mut]= []
        mut_diffs[mut]= []

        for idx in range(len(avail_sub)):
            mut_grid[mut].append(grids[idx][row,col])
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
            } for y in ['pval','diffs']
        }


## plot first for every mutation context.
for kmer in mut_grid.keys():
    fig_kmer= fig_dir + '/' + kmer
    os.makedirs(fig_kmer, exist_ok=True)
    fig_kmer= fig_kmer + '/'
    
    plot_data= {
        'pval': mut_grid[kmer],
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



####################################################
#################################################### grid SSD II
Nbins= 30
bins= np.linspace(0,1,Nbins)
bins= np.round(bins,4)
bins= [(bins[x-1],bins[x]) for x in range(1,len(bins))]

other_diffs= [count_data[s]['other'] for s in avail_sub]
pop_vector= [count_data[s]['pop'] for s in avail_sub]
pop_set= list(set(pop_vector))

pop_batch_dict= {
    ba: {
        pop: [x for x in batch_dict[ba] if pop_vector[x] == pop] for pop in pop_set
    } for ba in batch_dict.keys()
}

xlab= 'relative sampling'
ylab= 'mean matrix p-val'

view_sets= ['ref','anti']
grid_whole= {
    pop: {    
        view:{} for view in view_sets
    } for pop in pop_set
}

for i in batch_dict.keys():
    for pop in pop_batch_dict[i].keys():
        plt.figure(figsize=(20, 10))

        xprep= [pop_proportions[x] for x in pop_batch_dict[i][pop]]
        
        xprep= {
            sum(bi) / 2: [x for x in range(len(xprep)) if xprep[x] > bi[0] and xprep[x] <= bi[1]] for bi in bins
        }
        #xprep= {
        #     z: [x for x in range(len(xprep)) if xprep[x] == z] for z in list(set(xprep))
        #}

        ### grids
        batch_grids= [grid_diffs[x] for x in pop_batch_dict[i][pop]]
        y_prep= {
            z: [batch_grids[x] for x in xprep[z]] for z in xprep.keys()
        }

        y_prep= {
            z: [np.sqrt(np.sum(x**2)) for x in y_prep[z]] for z in y_prep.keys()
        }

        surface= sorted(xprep.keys())
        y= [np.mean(y_prep[x]) for x in surface]
        error= [np.std(y_prep[x]) for x in surface]

        grid_whole[pop]['ref'][i]= [surface,y,error]

        plt.errorbar(surface,y,yerr=error,label= 'ref')

        ###
        batch_grids= [other_diffs[x] for x in pop_batch_dict[i][pop]]
        xprep= [pop_proportions[x] for x in pop_batch_dict[i][pop]]
        xprep= np.repeat(xprep,[len(x) for x in batch_grids])
        batch_grids= list(it.chain(*batch_grids))

        xprep= {
            sum(bi) / 2: [x for x in range(len(xprep)) if xprep[x] > bi[0] and xprep[x] <= bi[1]] for bi in bins
        }
        
        y_prep= {
            z: [batch_grids[x] for x in xprep[z]] for z in xprep.keys()
        }

        y_prep= {
            z: [np.sqrt(np.sum(x**2)) for x in y_prep[z]] for z in y_prep.keys()
        }

        surface= sorted(xprep.keys())
        y= [np.mean(y_prep[x]) for x in surface]
        error= [np.std(y_prep[x]) for x in surface]

        grid_whole[pop]['anti'][i]= [surface,y,error]


        plt.errorbar(surface,y,yerr=error,label= 'control')    

        plt.xlabel(xlab)
        #plt.ylim(0,1.5)
        plt.ylabel(ylab)
        plt.title('grid SSD / sample proportion - control')

        plt.legend()
        plt.savefig(fig_dir + 'gridSSD_{}_control.png'.format(pop),bbox_inches='tight')
        plt.close()

############################################################# 
############################################################# STRATA


for strata in ['pval','diffs']:

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


