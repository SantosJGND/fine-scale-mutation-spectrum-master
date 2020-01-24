
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
sims_dir= main_dir + 'mutation_counter/data/marVar100kb/'
diffs= False


fig_dir= 'Figures/Samp_increment/'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'


mutlog= 'toMut.log'
min_size= 5
sampling= [50,50,5]
stepup= 'increment'
sample_sim= 0

data, data_freqs = MC_sample_matrix_v1(min_size= min_size, samp= sampling, stepup= stepup, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,
                          muted_dir= muted_dir, diffs= diffs,
                       exclude= False,sample_sim= sample_sim)




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
            x: data_freqs[x[0]][x[1]] for x in ref_pair
        }
    
    return comb_stats





def md_increment_compare(data, muted_dir= '', tag_ref= '_ss'):
    '''
    Parse data dictionary.
        data: {sim: {counts:{pop:g}, Nvars:{pop:g}, sizes:{pop:g}}}
    i: use sim and pop IDs to create dictionary connecting original populations to 
    subset populations created using ind_assignment_scatter_v1.
    ii. for each population for each reference, organise first by popuation sizes (not proportions).
    iii. calculate pairwise differences between sets of counts are contiguous sample sizes. 
    '''
    avail= list(data.keys())
    ref_idx= [int(tag_ref in avail[x]) for x in range(len(avail) )]
    categ= {
        z: [x for x in range(len(avail)) if ref_idx[x] == z] for z in [0,1]
    }
    
    #### population size diffs per population per simulation
    pop_asso= {avail[x]:recursively_default_dict() for x in categ[0]}
    
    for av in categ[1]:
        dat= [x for x in data[avail[av]]['counts'].keys() if tag_ref in x]
        dat_size= [data[avail[av]]['sizes'][x] for x in dat]
        
        ref_sim= avail[av].split(tag_ref)[0]
        ref_pop= [x.split('.')[0].strip(tag_ref) for x in dat]
        dat_size= [dat_size[x] for x in range(len(dat))]
        dat_size= [round(x,3) for x in dat_size]
        
        for p in range(len(dat)):
            pop_asso[ref_sim][ref_pop[p]][dat_size[p]][avail[av]]= dat[p]
    
    d= 0
    ### combine simulation combination and population size ranges.
    stats_dict= recursively_default_dict()
    for ref_sim in pop_asso.keys():
        batch= ref_sim.split('C')[0]
        
        for pop in data[ref_sim]['counts'].keys():
            
            available_sizes= sorted(list(pop_asso[ref_sim][pop].keys()))
            available_sizes= [x for x in available_sizes if x < 50]
            size_counts= {}
            
            for si in available_sizes:
                t= [(v,g) for v,g in pop_asso[ref_sim][pop][si].items()]
                t= [data[v[0]]['counts'][v[1]] for v in t]
                t= [x.reshape(1,np.prod(x.shape)) / np.sum(x) for x in t if np.sum(x)]
                
                t= [x[0] for x in t]
                
                t= np.array(t)
                
                size_counts[si]= t
            
            stats_dict[ref_sim][pop]= {
                'means': [],
                'stds': []
            }
            
            for idx in range(1,len(available_sizes)):
                set1= size_counts[available_sizes[idx]]
                set2= size_counts[available_sizes[idx - 1]]
                
                dists= set_SSD(set1,set2)
                
                stats_dict[ref_sim][pop]['means'].append(np.mean(dists))
                stats_dict[ref_sim][pop]['stds'].append(np.std(dists))
            
            stats_dict[ref_sim][pop]['sizes']= available_sizes[1:]
                
    
    return stats_dict



def set_SSD(set1,set2):
    '''
    return sum of squared differences between every pair of vectors across two sets.
    '''
    dists= []
    
    for indian in set1:
        
        dist_vec= [(x - indian) for x in set2] #/ np.sum(indian + x)
        dist_vec= [z**2 for z in dist_vec]
        dist_vec= [np.sum(x) for x in dist_vec]
        dists.extend(dist_vec)
    
    return dists





#####################################
#####################################

stats_dict= md_increment_compare(data)

batch_keys= list(stats_dict.keys())
batch_list= [x.split('C')[0] for x in batch_keys]
batch_dict= {
    z: [batch_keys[x] for x in range(len(batch_keys)) if batch_list[x] == z] for z in list(set(batch_list))
}

## average results by pops within batches:
## assumes all simulations within the same batch have the same populations. 

batch_stats= {}

for ba in batch_dict.keys():
    pops= list(stats_dict[batch_dict[ba][0]].keys())
    
    pop_dict= {
        pop: {
            z: np.array([stats_dict[ref][pop][z] for ref in batch_dict[ba]]) for z in stats_dict[batch_dict[ba][0]][pops[0]].keys()
        } for pop in pops
    }
    
    pop_dict= {
        pop: {
            z: np.mean(g,axis= 0) for z,g in pop_dict[pop].items()
        } for pop in pops
    }
    
    batch_stats[ba]= pop_dict


########
######## ttest between comparisons.

from scipy import stats
thet= recursively_default_dict()

for ba in batch_stats.keys():
    for pop in batch_stats[ba].keys():
        ts= []
        sizes_avail= batch_stats[ba][pop]['sizes']
        
        for idx in range(1,len(sizes_avail)):
            rvs= {
                z: [batch_stats[ba][pop]['means'][z],batch_stats[ba][pop]['stds'][z]] for z in [idx-1,idx]
            } 
            rvs= {
                z: stats.norm.rvs(loc=g[0],scale=g[1],size=1000) for z,g in rvs.items()
            }
            
            tstat= stats.ttest_ind(rvs[idx],rvs[idx-1],equal_var = False)
            
            ts.append(tstat[1])
        
        thet[ba][pop]= {
            'ts': ts,
            'sizes': sizes_avail[1:]
        }



#####################
##################### PLOTS

#### plot convergence in distance

for ep in batch_stats.keys():
    for pop in batch_stats[ep].keys():

        plt.figure(figsize=(15,10))

        plt.xlabel('relative sample size')
        #plt.ylim(0,1.03)
        plt.ylabel('SSD')
        plt.title('kmer spectrum, sample size convergence.')

        
        plt.errorbar(batch_stats[ep][pop]['sizes'],batch_stats[ep][pop]['means'],yerr=batch_stats[ep][pop]['stds'],label= '-'.join([ep,pop]))

        plt.legend()
        plt.savefig(fig_dir + 'kmerSpectrum_convergence_{}.png'.format(pop),bbox_inches='tight')
        plt.close() 


    plt.figure(figsize=(15,10))

    plt.xlabel('relative sample size')
    #plt.ylim(0,1.03)
    plt.ylabel('SSD')
    plt.title('kmer spectrum, sample size convergence.')

    for pop in batch_stats[ep].keys():
        plt.errorbar(batch_stats[ep][pop]['sizes'],batch_stats[ep][pop]['means'],yerr=batch_stats[ep][pop]['stds'],label= '-'.join([ep,pop]))

    plt.legend()
    plt.savefig(fig_dir + 'kmerSpectrum_convergence_total.png',bbox_inches='tight')
    plt.close() 



## II. plot batch stats
#


for ep in thet.keys():
    for pop in thet[ep].keys():

        plt.figure(figsize=(15,10))

        plt.xlabel('relative sample size')
        #plt.ylim(0,1.03)
        plt.ylabel('log pval')
        plt.title('kmer spectrum, sample size convergence, pval.')

        
        plt.scatter(thet[ep][pop]['sizes'],[np.log(x) for x in thet[ep][pop]['ts']],label= '-'.join([ep,pop]))

        plt.legend()
        plt.savefig(fig_dir + 'kmerSpectrum_convergencePval_{}.png'.format(pop),bbox_inches='tight')
        plt.close() 


    plt.figure(figsize=(15,10))

    plt.xlabel('relative sample size')
    #plt.ylim(0,1.03)
    plt.ylabel('log pval')
    plt.title('kmer spectrum, sample size convergencen pval.')

    for pop in batch_stats[ep].keys():
        plt.scatter(thet[ep][pop]['sizes'],[np.log(x) for x in thet[ep][pop]['ts']],label= '-'.join([ep,pop]))

    plt.legend()
    plt.savefig(fig_dir + 'kmerSpectrum_convergencePval_total.png',bbox_inches='tight')
    plt.close() 

