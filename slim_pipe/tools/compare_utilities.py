import numpy as np
import itertools as it
from itertools import product
import os

from tools.plot_utilities import Population, frequency_breakdown, heatmap


def heatmap_mutation_labels():
    
    comp = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
    }
    ypos, ylabel = [], []

    mut_index = {}
    row, col = 0, 0

    labels= []

    for b2, d in [('A', 'T'), ('A', 'C'), ('A', 'G'),
                  ('C', 'T'), ('C', 'G'), ('C', 'A')]:

        for b1 in 'ACGT':
            row_lab= []
            col = 0
            ypos.append(row+0.5)
            if b1 == 'T' and b2 == 'C' and d == 'A':
                ylabel.append('5\'-'+b1)
            elif b1 == 'C':
                ylabel.append(b2+r'$\to$'+d+r'  '+b1)
            else:
                ylabel.append(b1)
            for b3 in 'ACGT':
                mut_index[(b1+b2+b3, d)] = (row, col)

                mut_index[(comp[b3]+comp[b2]+comp[b1], comp[d])] = (row, col)
                row_lab.append('_'.join([b1+b2+b3, d]))

                col += 1
            labels.append(row_lab)
            row += 1
    
    return labels



def get_available_muts(muted_log):
    ''' read log of mutation counts '''
    
    with open(muted_log,'r') as fp:
        available= fp.readlines()
    
    available= [x.strip() for x in available]
    
    return available


def pops_from_sim(sim,sims_dir= './mutation_counter/data/sims/',pop_set= True):
    '''read sim specific int to pop assignment, return pops.'''
    sim_dir= sims_dir + '{}/'.format(sim)
    ID_file= sim_dir + "ind_assignments.txt"

    pops= []
    with open(ID_file,'r') as sample_id_lines:
        for line in sample_id_lines:
            line= str.encode(line)
            sample_id, population = line.split()[:2]
            pops.append(population.decode())
    
    if pop_set:
        return list(set(pops))
    else:    
        return pops



def count_compare(sim, frequency_range= [0,1], p_value= 1e-5,muted_dir= './mutation_counter/data/mutation_count/',
                  sims_dir= './mutation_counter/data/sims/', exclude= False):
    
    ''' perform pairwise population comparison of mutation counts for particular simulation'''
    pops= pops_from_sim(sim,sims_dir= sims_dir)

    ### change this 
    focus= pops[0]

    ## chromosome 
    chromosomes= [sim.split('.')[0].split('C')[1]]
    chromosome_groups = [chromosomes]

    ## get population pairs:
    population_pairs = [[(i), (i + 1) % len(pops)] for i in range(len(pops))] 
    population_pairs= [[pops[x] for x in y] for y in population_pairs]
    population_pairs= list(it.chain(*population_pairs))
    #print(population_pairs)
    pop_pair_names= zip(population_pairs[::2],
                               population_pairs[1::2])

    pop_pair_names= ['-'.join(list(x)) for x in pop_pair_names]

    population_pairs= zip(population_pairs[::2],
                               population_pairs[1::2])
    
    ### 
    chrom_pop= list(product(chromosome_groups,list(population_pairs)))

    heatmaps = [
        heatmap(
            chromosomes, population_pair, frequency_range, exclude, 
            p_value, sim, muted_dir
        ) for chromosomes, population_pair in chrom_pop
    ]

    ratio_grids, significant_indices = zip(*heatmaps)
    
    return ratio_grids, significant_indices


def deploy_count(available, frequency_range= [0,1], p_value= 1e-5,muted_dir= './mutation_counter/data/mutation_count/',
                  sims_dir= './mutation_counter/data/sims/'):
    
    ''' deploy count_compare() across simulations read from. '''
    data= {}
    
    for sim in available:
        
        ratio_grids, significant_indices= count_compare(sim, frequency_range= frequency_range, p_value= p_value,
                                               muted_dir= muted_dir, sims_dir= sims_dir)
        
        data[sim] ={
            'grids':ratio_grids,
            'sigs': significant_indices
        }
        
    return data
