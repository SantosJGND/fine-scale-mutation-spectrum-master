
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
import gzip

def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)


## directories
main_dir= os.getcwd() + '/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'
muted_dir= main_dir + 'mutation_counter/data/mutation_count/'
sims_dir= main_dir + 'data/phase1_100kb/'
indfile= 'integrated_call_samples.20101123.ALL.panel_regions.txt'
diffs= True

mutlog= 'toMut.log'
min_size= 40
sampling= [5,100,5]
bases= 'ATCG'
ksize= 3
sample_sim= 0
freq_extract= True
stepup= ''

data, data_freqs = MC_sample_matrix_v1(min_size= min_size, samp= sampling, stepup= stepup,count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir, indfile= indfile,
                          muted_dir= muted_dir, diffs= diffs,sample_sim= sample_sim, freq_extract= freq_extract,
                       exclude= False)




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
                            muted_dir= '', tag_ref= '_ss',pool = False):
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

                if not pool:
                    continue
                
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


available= list(count_data.keys())
subsamp= len(count_data)
avail_sub= np.random.choice(available,subsamp,replace= False)


## extract statistics per mutation.

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


fig_dir= 'Figures/SFS'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'


##############################################################
############################################################## SFS

pop_vector= [count_data[s]['pop'] for s in avail_sub]
pop_set= list(set(pop_vector))

pop_batch_dict= {
    ba: {
        pop: [x for x in batch_dict[ba] if pop_vector[x] == pop] for pop in pop_set
    } for ba in batch_dict.keys()
}



mu= 1e-8
N_inds= 1000
xlab= 'sampling proportion'
ylab= 'sum of sqared differences.'




batch_sfs= {ba: {} for ba in pop_batch_dict.keys()}

for i in batch_dict.keys():
    for pop_i in pop_batch_dict[i].keys():
    
        sims= [avail_sub[x] for x in pop_batch_dict[i][pop_i]]

        counts= []
        props= [pop_proportions[x] for x in pop_batch_dict[i][pop_i]] 

        for sim_idx in sims:
            freqs_dict= count_data[sim_idx]['freqs']
            sfs= []
            sizes_sim= count_data[sim_idx]['sizes']
            N_inds= min(sizes_sim) #- int(min(sizes_sim) % 2 == 0)
            #N_inds= int(np.mean(sizes_sim))

            for pop in freqs_dict.keys():

                #print(gen_time)
                freqs= freqs_dict[pop]

                sizeN= [x[1] for x in freqs if x[1] > 0]
                freqs= np.repeat([x[0] for x in freqs if x[1] > 0],sizeN) / sizes_sim[pop]

                #freqs= [x for x in freqs if x > 0]
                bin_count= np.histogram(freqs,bins= N_inds,range= [0,1])[0]
                bin_count= bin_count / np.sum(bin_count)
                #bin_count= np.array(bin_count)
                sfs.append(bin_count)

            dist_vec= sfs[0] - sfs[1] 

            dist_vec= dist_vec**2
            #dist_vec= np.mean(dist_vec)
            dist_vec= np.sqrt(np.sum(dist_vec)) / N_inds
            #dist_vec= np.sqrt(sum(dist_vec)) / N_inds
            
            counts.append(dist_vec)

        props_dict= {
            z: [counts[x] for x in range(len(props)) if props[x] == z] for z in list(set(props))
        }

        props_sorted= sorted(props_dict.keys())
        props_means= [np.mean(props_dict[x]) for x in props_sorted]
        props_std= [np.std(props_dict[x]) for x in props_sorted]

        batch_sfs[i][pop_i]= [props_sorted,props_means,props_std]


    plt.figure(figsize=(20,10))

    plt.xlabel(xlab)
    #plt.ylim(0,1.03)
    plt.ylabel(ylab)
    plt.title('sfs_differences.')
    
    for pop in batch_sfs[i].keys():
        props_sorted,props_means,props_std = batch_sfs[i][pop]
        plt.errorbar(props_sorted,props_means,yerr=props_std,label= pop)
    
    plt.savefig(fig_dir + 'SFS_{}.png'.format(i),bbox_inches='tight')
    plt.close()


