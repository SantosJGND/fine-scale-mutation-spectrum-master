
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



def md_reference_comp(data,p_value= 1e-5, test_m= 'fisher', individually= False, Nbins= 10,
                            exclude= False, frequency_range= [0,1], data_freqs= {}, extract= 'pval',
                            muted_dir= '', tag_ref= '_ss'):
    '''
    Parse data dictionary.
        data: {sim: {counts:{pop:g}, Nvars:{pop:g}, sizes:{pop:g}}}
    i: use sim and pop IDs to create dictionary connecting original populations to 
    subset populations created using ind_assignment_scatter_v1.
    ii: for each pair of reference populations, launch heatmapv2. return grid pvals or proportions,
    and proportion of mutations in subset population. allows for fisher or chi2 test for pval.
    '''
    
    bins= np.linspace(0,1,Nbins)
    bins= np.round(bins,4)
    bins= [(bins[x-1],bins[x]) for x in range(1,len(bins))]
            
    avail= list(data.keys())
    ref_idx= [int(tag_ref in avail[x]) for x in range(len(avail) )]
    categ= {
        z: [x for x in range(len(avail)) if ref_idx[x] == z] for z in [0,1]
    }
    
    print([len(categ[x]) for x in [0,1]])
    
    ### possible combinations per simulation.
    ref_combos= {}
       
    for idx in categ[0]:
        ref= avail[idx]
        ref_combs= list(data[ref]['counts'].keys())
        ref_combs= it.combinations(ref_combs,2)
        ref_combs= list(ref_combs)
        
        comb_dict= {
            x: {} for x in ref_combs
        }
        
        comb_stats= {}
        
        for pair in ref_combs:
            batch= ref.split('C')[0]
            pop1, pop2= pair
            
            sizes= [data[ref]['sizes'][x] for x in pair]
            #
            
            chromosomes= [ref.split('.')[0].split('C')[1]]
            
            pop_counts= {
                x: data[ref]['counts'][x] for x in pair
            }

            num_variants= {
                z: data[ref]['Nvars'][z] for z in pair
            }

            ratio_grid, sig_cells= heatmap_v2(chromosomes,pop_counts,num_variants,
                                              {},frequency_range, exclude, 
                                                p_value, muted_dir,tag= '',test= test_m,output= 'pval')

            pop_counts= {
                z: pop_counts[z] / np.sum(pop_counts[z]) for z in pop_counts.keys()
            }
            
            dist_prop= pop_counts[pop1] / pop_counts[pop2]
            dist_prop= np.nan_to_num(dist_prop)

            grid_diffs= pop_counts[pop1] - pop_counts[pop2]

            comb_stats[pair]= {
                'grids': ratio_grid,
                'sigs': sig_cells,
                'sizes': sizes,
                'batch': batch,
                'diffs': grid_diffs
            }

            if data_freqs[ref]:
                comb_stats[pair]['freqs']= {
                    pop1: data_freqs[ref][pop1],
                    pop1: data_freqs[ref][pop2]
                }
        
        ref_combos[ref]= {
            'combs': comb_dict,
            'sizes': data[ref]['sizes'],
            'stats': comb_stats
        }
    
    #### population size diffs per population per simulation
    pop_asso= {avail[x]:recursively_default_dict() for x in categ[0]}
    
    for av in categ[1]:
        dat= [x for x in data[avail[av]]['counts'].keys() if tag_ref in x]
        dat_size= [data[avail[av]]['sizes'][x] for x in dat]
        ref_sim= avail[av].split(tag_ref)[0]
        ref_pop= [x.split('.')[0].strip(tag_ref) for x in dat]
        dat_size= [dat_size[x] / data[ref_sim]['sizes'][ref_pop[x]] for x in range(len(dat))]
        dat_size= [round(x,3) for x in dat_size]
        for p in range(len(dat)):
            pop_asso[ref_sim][ref_pop[p]][dat_size[p]][avail[av]]= dat[p]
    
    d= 0
    ### combine simulation combination and population size ranges.
    
    for ref_sim in pop_asso.keys():
        print(ref_sim)
        batch= ref.split('C')[0]
        
        for combo in ref_combos[ref_sim]['combs'].keys():
            
            pop1, pop2= combo
            
            available_sizes= {
                z: sorted(list(pop_asso[ref_sim][z].keys())) for z in combo
            }
            #available_sizes= {
            #    z: [round(x / ref_combos[ref_sim]['sizes'][z], 3) for x in available_sizes[z]] for z in combo
            #}
            bins_dict= {
                b: {
                    z: [x for x in available_sizes[z] if x > b[0] and x <= b[1]] for z in combo
                } for b in bins
            }
            
            print([len(bins_dict[b][pop1]) for b in bins])
            print([len(bins_dict[b][pop2]) for b in bins])
            
            bins_combs= {
                b: [(x,y) for x in bins_dict[b][pop1] for y in bins_dict[b][pop2]] for b in bins
            }

            print([len(bins_combs[c]) for c in bins_combs.keys()])
            
            for bend in bins_combs.keys():
                ref_combos[ref_sim]['combs'][combo][bend]= []
                
                for size_combo in bins_combs[bend]:
                    i,j= size_combo
                    
                    for sub1 in pop_asso[ref_sim][pop1][i].keys():
                        for sub2 in pop_asso[ref_sim][pop2][j].keys():
                            
                            ref_pair= {
                                sub1: pop_asso[ref_sim][pop1][i][sub1],
                                sub2: pop_asso[ref_sim][pop2][j][sub2]
                            }
                            
                            sizes= [data[x]['sizes'][g] for x,g in ref_pair.items()]
                            #

                            chromosomes= [ref_sim.split('.')[0].split('C')[1]]

                            pop_counts= {
                                x: data[g]['counts'][x] for g,x in ref_pair.items()
                            }

                            num_variants= {
                                x: data[z]['Nvars'][x] for z,x in ref_pair.items()
                            }

                            ratio_grid, sig_cells= heatmap_v2(chromosomes,pop_counts,num_variants,
                                                              {},frequency_range, exclude, 
                                                                p_value, muted_dir,tag= '',test= test_m,output= 'pval')

                            pop_counts= {
                                z: s / np.sum(s) for z,s in pop_counts.items()
                            }

                            dist_prop= pop_counts[ref_pair[sub1]] / pop_counts[ref_pair[sub2]]
                            dist_prop= np.nan_to_num(dist_prop)

                            grid_diffs= pop_counts[ref_pair[sub1]] - pop_counts[ref_pair[sub2]]

                            comb_stats= {
                                'grids': ratio_grid,
                                'sigs': sig_cells,
                                'sizes': sizes,
                                'batch': batch,
                                'diffs': grid_diffs
                            }

                            if data_freqs:
                                comb_stats['freqs']= {
                                    pop1: data_freqs[sub1][ref_pair[sub1]],
                                    pop2: data_freqs[sub2][ref_pair[sub2]]
                                }

                            ref_combos[ref_sim]['combs'][pair][bend].append(comb_stats)
    
    return pop_asso, ref_combos




p_value= 1e-5
test_m= 'chi2'
individually= False
exclude= False
frequency_range= [0,1]
extract= 'pval'
Nbins= 20

pop_asso, ref_combos= md_reference_comp(data,p_value= p_value, test_m= test_m, individually= individually, Nbins= Nbins,
                                        exclude= exclude, frequency_range= frequency_range, extract= extract,
                                     muted_dir= muted_dir, data_freqs= data_freqs)




###########



bins= np.linspace(0,1,Nbins)
bins= np.round(bins,4)
bins= [(bins[x-1],bins[x]) for x in range(1,len(bins))]

requested= 'diffs'

sim_stats= list(ref_combos.keys())

combo_grids= {
    ref: {comb: {
            z: [] for z in bins
        } for comb in ref_combos[ref]['combs'].keys()
    } for ref in ref_combos.keys()
}

combo_ref= {x: [] for x in ref_combos[sim_stats[0]]['combs'].keys()}

for sim in sim_stats:
    for comb in ref_combos[sim]['combs'].keys():
        combo_ref[comb].append(ref_combos[sim]['stats'][comb][requested])
        
        for bi in ref_combos[sim]['combs'][comb].keys():
            for idx in ref_combos[sim]['combs'][comb][bi]:
                combo_grids[sim][comb][bi].append(idx[requested])



grad_comb= {
    ref: {
        comb: {
            bi: [np.sqrt(np.sum((x - ref_combos[ref]['stats'][comb][requested])**2)) for x in combo_grids[ref][comb][bi]] for bi in combo_grids[ref][comb].keys()
        } for comb in combo_grids[ref].keys()
    } for ref in combo_grids.keys()
}

grad_comb= {
    comb: {
        bi: list(it.chain(*[grad_comb[r][comb][bi] for r in grad_comb.keys()])) for bi in bins
    } for comb in ref_combos[sim_stats[0]]['combs'].keys()
}


#################################################
#################################################

ydict= {
    comb: {
        round(sum(bi)/ 2,3): np.mean(g) for bi,g in grad_comb[comb].items()
    } for comb in grad_comb.keys()
}

yerror= {
    comb: {
        round(sum(bi)/ 2,3): np.std(g) for bi,g in grad_comb[comb].items()
    } for comb in grad_comb.keys()
}


xdict= {
    comb: sorted(ydict[comb].keys()) for comb in ydict.keys()
}

xlab= "relative size"
ylab= "norm Euc distance"

plt.figure(figsize=(15,10))

plt.xlabel(xlab)
#plt.ylim(0,1.03)
plt.ylabel(ylab)
plt.title('Differences to pairwise full count comp')

for i in xdict.keys():
    plt.errorbar(xdict[i],[ydict[i][x] for x in xdict[i]],yerr=[yerror[i][x] for x in xdict[i]],label= i)

plt.legend()
plt.savefig(fig_dir + 'Pair_wise difference.png',bbox_inches='tight')
plt.close() 


fig= [
    go.Scatter(
        x= xdict[i],
        y= [ydict[i][x] for x in xdict[i]],
        error_y= dict(
            array= [yerror[i][x] for x in xdict[i]],
            type= 'data',
            #symmetric= True,
            visible=True
        ),
        name= '-'.join(i)
    ) for i in xdict.keys()
]

layout= go.Layout()

Figure= go.Figure(data=fig, layout= layout)
iplot(Figure)