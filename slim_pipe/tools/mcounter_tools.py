import re
import numpy as np
import pandas as pd 
import gzip
from itertools import product
import os

from tools.compare_utilities import (
    get_available_muts, count_compare, deploy_count, pops_from_sim, check_availability, clean_empty
)


import time

import collections

def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)



###########################################################
########################################################### data management

def process_log(muted, sims_dir= ''):
    '''
    verify directories indicated in log exist in given directory.
    '''
    available= get_available_muts(muted)

    ### cleaning data set 
    ### i.e. accounting for aborted runs.
    available, miss_data= check_availability(available, dir_check=sims_dir)
    available, empty= clean_empty(available,str_format= '',dir_check= sims_dir,requested= ['.vcf.gz'])

    return available



def process_dir(sims_dir= ''):
    '''
    verify directories indicated in log exist in given directory.
    '''
    available= [name for name in os.listdir(sims_dir)]
    
    ### cleaning data set 
    ### i.e. accounting for aborted runs.
    available, miss_data= check_availability(available, dir_check=sims_dir)
    available, empty= clean_empty(available,str_format= '',dir_check= sims_dir,requested= ['.vcf.gz','fa.gz'])
    
    print('missing: {}, no vcf: {}'.format(len(miss_data),len(empty)))
    return available



def dict_write(new_dict,inds,outemp= 'ind_assignments{}.txt',dir_sim= '',tag= ''):
    '''
    cofactor to ind_assignment_scatter()
    '''
    inds= np.array(inds)
    new_array= [[(inds[x,0],v) for x in new_dict[v]] for v in new_dict.keys()]
    new_array= list(it.chain(*new_array))
    new_array= ['\t'.join(x) for x in new_array]
    
    out= dir_sim + outemp.format(tag)
    
    with open(out,'w') as f:
        f.write('\n'.join(new_array))



##### combinations

def ind_assignment_scatter(reference,main_dir= '',pops= 'ind_assignments.txt',
                          min_size= 80, samp= [5,20,10], outemp= 'ind_assignments{}.txt'):
    '''
    read ind assignments for a given window; 
    chose one population;
    subset that pop in some way. 
    '''
    dir_sim= dir_launch + '/data/sims/' + reference + '/'
    ind_assignments= dir_sim + pops
    
    with open(ind_assignments,'r') as f:
        inds= f.readlines()
    
    inds= [x.split() for x in inds]
    pops= np.array(inds)[:,1]
    pop_dict= {
        z: [x for x in range(len(pops)) if pops[x] == z] for z in list(set(pops))
    }
    
    tag_list= []
    
    ## criterium of choice. chose only one pop.
    pop_chose= [x for x in pop_dict.keys() if len(pop_dict[x]) >= min_size]
    if len(pop_chose):
        pop_chose= pop_chose[0]
        N= len(pop_dict[pop_chose])
        pop_list= pop_dict[pop_chose]

        for each in np.linspace(samp[0],N,samp[1]):  
            each= int(each)
            for perm in range(samp[2]):
                tag= '_' + '.'.join([pop_chose,str(each),str(perm)])
                
                smaller= np.random.choice(pop_list,each,replace= False)
                smaller= [int(x in smaller) for x in pop_list]

                new_pop= {
                    tag + '.s' + str(z): [pop_list[x] for x in range(len(smaller)) if smaller[x] == z] for z in list(set(smaller))
                }
                
                new_dict= {v:g for v,g in pop_dict.items() if v != pop_chose}
                new_dict.update(new_pop)

                dict_write(new_dict,inds,outemp= outemp, dir_sim= dir_sim, tag= tag)

                tag_list.append(tag)

    return tag_list




def ind_assignment_scatter_v1(reference,dir_sim= '',indfile= 'ind_assignments.txt',
                          min_size= 80, samp= [5,20,10], stepup= "increment",outemp= 'ind_assignments{}.txt',write_out= False):
    '''
    read ind assignments for a given window; 
    chose one population;
    subset that pop in some way.
    - v1: instead of writting new pop_assignment files, return them. 
    '''
    
    ind_assignments= dir_sim + reference + '/' + indfile
    
    with open(ind_assignments,'r') as f:
        inds= f.readlines()
    
    inds= [x.split() for x in inds]
    pops= np.array(inds)[:,1]
    pop_dict= {
        z: [x for x in range(len(pops)) if pops[x] == z] for z in list(set(pops))
    }
    
    tag_list= []
    tag_dict= {}
    
    ## criterium of choice. chose only one pop.
    pop_avail= [x for x in pop_dict.keys() if len(pop_dict[x]) >= min_size]
    for pop_chose in pop_avail:
        
        N= len(pop_dict[pop_chose])
        pop_list= pop_dict[pop_chose]

        if stepup== 'increment':
            timetable= np.linspace(1,samp[0],samp[1])
        else:
            timetable= np.linspace(samp[0],N,samp[1])

        for each in timetable:  
            each= int(each)
            for perm in range(samp[2]):
                tag= '_ss' + '.'.join([pop_chose,str(each),str(perm)])
                
                smaller= np.random.choice(pop_list,each,replace= False)
                smaller= [int(x in smaller) for x in pop_list]
                
                new_pop= {
                    tag + '.s' + str(z): [pop_list[x] for x in range(len(smaller)) if smaller[x] == z] for z in [1]
                }
                
                #new_dict= {v:g for v,g in pop_dict.items() if v != pop_chose}
                #new_dict.update(new_pop)
                new_dict= new_pop

                if write_out:
                    dict_write(new_dict,inds,outemp= outemp, dir_sim= dir_sim, tag= tag)
                else:
                    tag_dict[tag]= new_dict
                tag_list.append(tag)

    if write_out:
        return tag_list
    else: 
        return tag_list, tag_dict, pop_dict


#######################
#######################  Allele count 

from tools.fasta_utilities import (
	geno_muts_v2, get_mutations, get_by_path, vcf_muts_matrix,
	kmer_comp_index, kmer_mut_index
	)


def MC_sample_matrix(logfile, min_size= 80, samp= [5,20,10], pops= 'ind_assignments.txt', outemp= 'ind_assignments{}.txt',
                    count_dir= './count/', dir_launch= '..',main_dir= './', sim_dir= 'mutation_counter/data/sims/', muted_dir= 'mutation_counter/data/mutation_count/',
                    outlog= 'indy.log', row= 24,col= 4,exclude= False):
    '''
    launch mutation counter pipeline on manipulated population assignments.
    Use matrix multiplication to extract counts. 
    '''
    
    sims= process_dir(sims_dir= main_dir+sim_dir)
    print(len(sims))
    tags= []
    sim_extend= []
    chroms= []
    
    data= {}
    
    for sim in sims:
        
        ## chromosome
        chrom= sim.split('.')[0].split('C')[-1].strip('chr')
        chromosomes= [sim.split('.')[0].split('C')[1]]
        chromosome_groups = [chromosomes]

        if exclude:
            files= read_exclude()
        else:
            files= {}

        ### read vcf

        row_info= 6
        header_info= 9
        phased= False
        vcf_dir= sim_dir + sim + '/'
        vcf_file= vcf_dir + sim + '_' + 'chr' + chrom + '.vcf.gz'

        genotype, summary, Names= read_geno_nanumv3(vcf_file, header_info= header_info,phased= phased)
        
        
        ## read fasta
        fasta_file= vcf_dir + 'chr{}_{}.fa.gz'.format(chrom,sim)

        with gzip.open(fasta_file,'r') as f:
            lines= f.readlines()
            lines= [x.decode() for x in lines]

        refseq= lines[1].strip()

        ###
        positions= [int(x) for x in summary.POS]
        wstart= int(min(positions))
        wend= int(max(positions))
        
        Wlen= wend - wstart
        ksize= 3 # odd.
        bases = 'ACGT'
        collapsed= True
        
        
        genotype_parse= [x for x in range(summary.shape[0]) if int(summary.POS[x])-1 >= wstart and int(summary.POS[x])-1 <= wend]
        Window= genotype[:,genotype_parse]
        subset_summary= summary.loc[genotype_parse,:].reset_index()
        
        ##
        mut_matrix, flag_reverse= vcf_muts_matrix_v1(refseq,subset_summary,start= wstart,end= wend,ksize= ksize,bases=bases, collapse= collapsed)
        if flag_reverse:
            Window[:,flag_reverse]= 2 - Window[:,flag_reverse]
        
        ind_collapsed_mat= geno_muts_v2(np.array(Window), mut_matrix)
        
        tag_list, tag_dict, pop_dict= ind_assignment_scatter_v1(sim,main_dir= main_dir,
                          min_size= min_size, samp= samp, outemp= outemp)
        #print(tag_list)
        
        ## counts for no tag sim:
        pop_counts= {
            z: np.sum(ind_collapsed_mat[pop_dict[z],:],axis= 0) for z in pop_dict.keys()
        }
        
        pop_counts= {
            z:g.reshape(row,col) for z,g in pop_counts.items()
        }
        
        num_variants= {
            z: np.sum(ind_collapsed_mat[pop_dict[z],:]) for z in pop_dict.keys()
        }
        
        data[sim]= {
            'counts': pop_counts,
            'Nvars': num_variants,
            'sizes': {z:len(g) for z,g in pop_dict.items()}
        }
        
        
        if len(tag_list):
            ###
            sim_extend.append(sim)
            tags.append('')
            chroms.append(chrom)
            ###
            
            for idx in range(len(tag_list)):
                
                sim_extend.extend([sim]*len(tag_list))
                tags.extend(tag_list)
                chroms.extend([chrom]*len(tag_list))
                
                ##
                tag= tag_list[idx]
                ind_file= outemp.format(tags[idx])
                new_sim= sim + tag

                pop_dict= tag_dict[tag]

                pop_sizes= {
                    z: len(g) for z,g in pop_dict.items()
                }
                
                pops= list(set(pop_dict.keys()))
                
                ###
                pop_counts= {
                    z: np.sum(ind_collapsed_mat[pop_dict[z],:],axis= 0) for z in pop_dict.keys()
                }
                
                pop_counts= {
                    z:g.reshape(row,col) for z,g in pop_counts.items()
                }
                num_variants= {
                    z: np.sum(ind_collapsed_mat[pop_dict[z],:]) for z in pop_dict.keys()
                }
                
                data[new_sim]= {
                    'counts': pop_counts,
                    'Nvars': num_variants,
                    'sizes': {z:len(g) for z,g in pop_dict.items()}
                }
    
    return data




### New, adapt to using matrices instead of written files. 


def count_popKmers(Window, mut_matrix, pop_dict, single= True, frequency_range= [0,1],row=24,col=4):
    '''
    Extract population mutation counts from _ind x kmer_ mutation matrix. 
    '''
    
    pop_counts= {}
    num_variants= {}
    
    for pop in pop_dict.keys():
        pop_gen= Window[pop_dict[pop],:]
        freqs= np.sum(pop_gen,axis= 0) / pop_gen.shape[0]
        ## discount alleles outside freq range.
        in_out= (freqs < frequency_range[0]) | (freqs > frequency_range[1])
        
        pop_gen[:,in_out]= 0
        
        if single: 
            pop_gen= np.sum(pop_gen,axis= 0) > 0
            pop_gen= np.array(pop_gen,dtype= int).reshape(1,len(pop_gen))
        
        pop_collapsed_mat= geno_muts_v2(pop_gen, mut_matrix)
        pop_summed= np.sum(pop_collapsed_mat,axis= 0)
        
        pop_counts[pop]= pop_summed.reshape(row,col)

        num_variants[pop]= np.sum(pop_collapsed_mat)

    return {
        'counts': pop_counts,
        'Nvars': num_variants,
        'sizes': {z:len(g) for z,g in pop_dict.items()}
    }



def read_diffs(tag,diff_dir= '',start= 0):
    '''
    read file of differences to ancestral sequence.
    '''

    filename= diff_dir + tag + '_diffs.txt.gz'

    with gzip.open(filename,'r') as f:
    	snps = f.readlines()
        
    snps= [x.decode() for x in snps[1:] if len(x)]
    snps= [x.split() for x in snps if 'SNP' in x]

    snps= {
    	str(int(x[0]) - start): [x[2],x[3]] for x in snps
    }

    return snps




def MC_sample_matrix_v1(min_size= 80, samp= [5,20,10], stepup= "increment", diffs= False, frequency_range= [0,1],indfile= 'ind_assignments.txt', outemp= 'ind_assignments{}.txt',
                    count_dir= './count/', dir_launch= '..',main_dir= './', sim_dir= 'mutation_counter/data/sims/', muted_dir= 'mutation_counter/data/mutation_count/',
                    outlog= 'indy.log', row= 24,col= 4, single= True, exclude= False, print_summ= False, sample_sim= 0,collapsed= True,bases= 'ACGT',ksize= 3,ploidy= 2, freq_extract= False):
    '''
    launch mutation counter pipeline on manipulated population assignments.
    Use matrix multiplication to extract counts. 
    - v1 relies on count_popKmers() function to count mutations per pop. allows freq. filter and single mutaiton count.  
    '''
    
    ti= time.time()
    sims= process_dir(sims_dir= sim_dir)
    print('available {}'.format(len(sims)))

    tags= []
    sim_extend= []
    chroms= []
    
    data_kmer= {}
    data_freqs= {}
    #sim_sample= np.random.choice(sims,8,replace= False)
    if sample_sim == 0:
        sample_sim= len(sims)

    print('sample {}'.format(sample_sim))
    sim_sub= np.random.choice(sims,sample_sim,replace= False)
    
    for sim in sim_sub:
        
        ## chromosome
        chrom= sim.split('.')[0].split('C')[-1].strip('chr')
        chromosomes= [sim.split('.')[0].split('C')[1]]
        chromosome_groups = [chromosomes]

        if exclude:
            files= read_exclude()
        else:
            files= {}

        ### read vcf

        vcf_dir= sim_dir + sim + '/'
        vcf_file= vcf_dir + sim + '_' + 'chr' + chrom + '.vcf.gz'
        
        t0= time.time()
        print(sim)

        genotype, summary, Names= read_vcf_allel(vcf_file)
        t1= time.time()
        
        read_time= t1- t0

        if len(genotype) == 0:
            continue
        
        print(genotype.shape)
        
        ## read fasta
        fasta_file= vcf_dir + 'chr{}_{}.fa.gz'.format(chrom,sim)

        with gzip.open(fasta_file,'r') as f:
            lines= f.readlines()
            lines= [x.decode() for x in lines]

        refseq= lines[1].strip()

        ### 
        positions= [int(x) for x in summary.POS]
        wstart= int(min(positions))
        wend= int(max(positions))
        
        Wlen= wend - wstart
        
        genotype_parse= [x for x in range(summary.shape[0]) if int(summary.POS[x])-1 >= wstart and int(summary.POS[x])-1 <= wend]
        Window= genotype[:,genotype_parse]
        subset_summary= summary.loc[genotype_parse,:].reset_index()
        
        ##
        t0= time.time()
        mut_matrix, flag_reverse, flag_remove= vcf_muts_matrix_v1(refseq,subset_summary,start= wstart,end= wend,ksize= ksize,
        													bases=bases, collapse= collapsed)
        
        retain= [x for x in range(Window.shape[1]) if x not in flag_remove]
        Window= Window[:,retain]
        subset_summary= subset_summary.loc[retain,:].reset_index()

        t1= time.time()
        time_mut= t1 - t0

        if diffs:
        	sim_start= sim.split('.')[-1]
        	diff_snps= read_diffs(sim,diff_dir= vcf_dir, start= int(sim_start))

        	summary_diff= [x for x in range(subset_summary.shape[0]) if subset_summary.POS[x] in diff_snps.keys()]

        	flag_reverse.extend(summary_diff)
        	flag_reverse= list(set(flag_reverse))
        
        
        if flag_reverse:
            Window[:,flag_reverse]= ploidy - Window[:,flag_reverse]
        
        ind_collapsed_mat= geno_muts_v2(np.array(Window), mut_matrix)
        
        tag_list, tag_dict, pop_dict= ind_assignment_scatter_v1(sim,dir_sim= sim_dir,
                          min_size= min_size, samp= samp, stepup= stepup, outemp= outemp,indfile= indfile)
        #print(tag_list)
        total_inds= sum([len(x) for x in pop_dict.values()])
        if Window.shape[0] < total_inds:
            continue
        ## counts for no tag sim:
        s0= time.time()
        data_kmer[sim]= count_popKmers(Window, mut_matrix, pop_dict, single= single, 
                                  frequency_range= frequency_range,row=row,col=col)

        if freq_extract:
            pop_freqs= pop_dict_SFS(Window,pop_dict)
            data_freqs[sim]= pop_freqs
        
        t1= time.time()
        count_time= t1- t0
        
        if len(tag_list):
            ###
            sim_extend.append(sim)
            tags.append('')
            chroms.append(chrom)
            ###
            
            for idx in range(len(tag_list)):
                
                sim_extend.extend([sim]*len(tag_list))
                tags.extend(tag_list)
                chroms.extend([chrom]*len(tag_list))
                
                ##
                tag= tag_list[idx]
                ind_file= outemp.format(tags[idx])
                new_sim= sim + tag

                pop_dict= tag_dict[tag]
                
                data_kmer[new_sim]= count_popKmers(Window, mut_matrix, pop_dict, single= single, 
                                  frequency_range= frequency_range,row=row,col=col)

                if freq_extract:
                    pop_freqs= pop_dict_SFS(Window,pop_dict)
                    data_freqs[new_sim]= pop_freqs
                

        if print_summ:
            print('mut_matrix time: {} s'.format(time_mut / 60))
            print('count time: {} s'.format(count_time / 60))
            print('est total count time: {} s'.format(count_time*len(tag_list) / 60))
            print('replicates: {}'.format(len(tag_list)))
            print('read time: {} s'.format(read_time / 60))

    tf= time.time()
    time_elapsed= tf - ti
    
    print('time elapsed: {}s'.format(time_elapsed))
    

    return data_kmer, data_freqs


##################################################################
##### mutation matrix construction - from mutation counter utilities. 

### new, control for switch in alt- ref allele.
### attention: this shouldn"t be necessary at this point, as the data were simulated from
### the fasta file being read. 

def vcf_muts_matrix_v1(refseq,summary,start= 0,end= 0,ksize= 3,bases='ATCG', collapse= True):
    ''' 
    Return matrix of mutation contexts by SNP in genotype array
    Each mutation is mapped to list of possible mutations as a binary vector.
    - v1 determines if alternative allele = reference allele in fasta. 
        if so, allele is switched, position idx is flagged. 
    '''
    
    mutations= get_mutations(bases= bases,ksize= ksize)
    kmers, kmer_idx= kmer_comp_index(mutations)
    
    mut_lib= kmer_mut_index(mutations)
    
    if end == 0:
        end= max(summary.POS)
    
    k5= int(ksize/2)
    k3= ksize - k5
    pos_mut= []
    flag_reverse= []
    flag_remove= []
    
    for x in range(summary.shape[0]):
        pos= int(summary.POS[x]) - 1
        if pos >=  start and pos <= end:
            kmer= refseq[pos-k5: pos + k3]
            if 'N' in kmer:
                flag_remove.append(x)
                continue
            mut= kmer + summary.ALT[x]
            
            if kmer[1] == summary.ALT[x]:
                flag_reverse.append(x)
                mut= kmer+summary.REF[x]
            
            if len(mut) != 4: 
                print(kmer)
                print(summary.REF[x],summary.ALT[x])
                print(x,pos)
                print(len(refseq),summary.shape[0])
                if collapse:
                    mut_array=np.zeros(len(kmer_idx))
                    pos_mut.append(mut_array)
                    continue
                else:
                    mut_array=np.zeros(len(mutations))
                    pos_mut.append(mut_array)
                    continue
            if collapse:
                mut_index= kmers[mut]
                mut_array=np.zeros(len(kmer_idx))
            else:
                mut_index= get_by_path(mut_lib, list(mut))
                mut_array=np.zeros(len(mutations))
            
            mut_array[mut_index]= 1
            pos_mut.append(mut_array)
    
    pos_mut= np.array(pos_mut).T
    
    return pos_mut, flag_reverse, flag_remove


#####################################################################
#####################################################################
### v1 in plot_utilities module.
from scipy.stats import chi2_contingency
from fisher import pvalue

def heatmap_v2(chromosomes,pop_counts, num_variants, population_dict,frequency_range, exclude, 
                p_value, muted_dir,tag= '',output= 'pval',row= 24, col= 4, test= 'fisher'):

    '''
    pairwise comparison of count matrices. Chi2 applied cell-wise. 
    p-value or proportion - output argument. 
    - v2: count matrices are provided in pop_counts dictionary. 
    '''
    if exclude:
        files= read_exclude()
    else:
        files= {}
    
    refpop, pop = list(pop_counts.keys())

    ratio_grid = np.zeros((row, col))
    sig_x, sig_y = [], []
    
    for i in range(row):
        for j in range(col):
            chi_array= np.array([
                    [pop_counts[pop][i][j], num_variants[pop]],
                    [pop_counts[refpop][i][j], num_variants[refpop]]
                ])

            chi_0= np.sum(chi_array,axis= 1)
            chi_1= np.sum(chi_array,axis= 0)
            
            if chi_0[0] == 0 or chi_0[1] == 0:
                ratio_grid[i][j] = np.nan
                sig_x.append(j+0.5)
                sig_y.append(i+0.5)
            
            elif chi_1[0] == 0 or chi_1[1] == 0:
                ratio_grid[i][j] = 1
            
            else:
                ##
                if test == 'chi2':
                    _, this_pval, _, _ = chi2_contingency(
                        chi_array
                    )
                else:
                    p= pvalue(pop_counts[pop][i][j], num_variants[pop],
                        pop_counts[refpop][i][j], num_variants[refpop])
                    this_pval= p.two_tail
                    
                if output == 'pval':
                    ratio_grid[i][j] = this_pval
                else:
                    ratio_grid[i][j] = (pop_counts[pop][i][j] * num_variants[refpop] /
                                        (num_variants[pop] * pop_counts[refpop][i][j]))
                if this_pval < p_value:
                    sig_x.append(j+0.5)
                    sig_y.append(i+0.5)

    return ratio_grid, (sig_x, sig_y)


#############################################################
#############################################################

import allel


def read_vcf_allel(file_vcf):
    '''
    Use scikit allel to read vcf file. Organise variant information into summary pandas df. 
    '''
    
    vcf_ori= allel.read_vcf(file_vcf)
    
    if not vcf_ori:
        print('empty vcf.')
        return {}, {}, {}

    #print(vcf_ori.keys())
    ### get genotype array
    geno= vcf_ori['calldata/GT']

    mult_alt= []
    indel= []
    single= []

    for idx in range(geno.shape[0]):
            ## eliminate +1 segregating mutations.
            if vcf_ori['variants/ALT'][idx][1]:
                gen_t= geno[idx]
                gen_t[gen_t > 1] = 0
                geno[idx]= gen_t


            if len(vcf_ori['variants/REF'][idx]) != 1 or len(vcf_ori['variants/ALT'][idx][0]) != 1:
                indel.append(idx)
            else:
                single.append(idx)

    

    
    geno= allel.GenotypeArray(geno)
    geno= geno.to_n_alt().T
    
    ## setup summary
    column_names= ['CHROM','POS','ID','REF','ALT','QUAL','FILTER']

    alts= [vcf_ori['variants/ALT'][x][0] for x in range(geno.shape[1])]
    PASS= [['.','PASS'][int(vcf_ori['variants/FILTER_PASS'][x])] for x in range(geno.shape[1])]

    summary= [
        vcf_ori['variants/CHROM'],
        vcf_ori['variants/POS'],
        vcf_ori['variants/ID'],
        vcf_ori['variants/REF'],
        alts,
        vcf_ori['variants/QUAL'],
        PASS,

    ]
    
    summary= np.array(summary).T
    
    if len(indel):
        #print('mutliple ref loci: {}'.format(geno.shape[1] - len(indel)))
        geno= geno[:,single]
        summary= summary[single,:]
    
    summary= pd.DataFrame(summary,columns= column_names)
    
    return geno, summary, vcf_ori['samples']


################################################################
################################################################
################################################################


def mcounter_deploy(data,p_value= 1e-5, test_m= 'fisher', individually= False,
                            exclude= False, frequency_range= [0,1], data_freqs= {}, extract= 'pval',
                            muted_dir= '', tag_ref= '_ss'):
    '''
    Parse data dictionary.
        data: {sim: {counts:{pop:g}, Nvars:{pop:g}, sizes:{pop:g}}}
    i: use sim and pop IDs to create dictionary connecting original populations to 
    subset populations created using ind_assignment_scatter_v1.
    ii: for each pair of reference/subset populations, launch heatmapv2. return grid pvals or proportions,
    and proportion of mutations in subset population. allows for fisher or chi2 test for pval.
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

        for pop in pop_asso[ref].keys():
            for sub in pop_asso[ref][pop].keys():

                batch= ref.split('C')[0]

                pop_dict= {
                    ref: pop,
                    sub: pop_asso[ref][pop][sub]
                }

                sizes= [data[ref]['sizes'][pop], data[sub]['sizes'][pop_asso[ref][pop][sub]]]
                #print(sizes)

                chromosomes= [x.split('.')[0].split('C')[1] for x in pop_dict.keys()]

                pop_counts= {
                    x: data[x]['counts'][z] for x,z in pop_dict.items() 
                }

                num_variants= {
                    x: data[x]['Nvars'][z] for x,z in pop_dict.items() 
                }

                ratio_grid, sig_cells= heatmap_v2(chromosomes,pop_counts,num_variants,
                                                  pop_dict,frequency_range, exclude, 
                                                    p_value, muted_dir,tag= '',test= test_m,output= 'pval')

                pop_counts[sub]= pop_counts[sub] / np.sum(pop_counts[sub])
                pop_counts[ref]= pop_counts[ref] / np.sum(pop_counts[ref])

                dist_prop= pop_counts[sub] / pop_counts[ref]
                dist_prop= np.nan_to_num(dist_prop)

                grid_diffs= pop_counts[sub] - pop_counts[ref]

                count_data[d]= {
                    'grids': ratio_grid,
                    'sigs': sig_cells,
                    'sizes': sizes,
                    'batch': batch,
                    'prop': dist_prop,
                    'pop': pop,
                    'diffs': grid_diffs
                }

                if data_freqs:
                    count_data[d]['freqs']= {
                        0: data_freqs[ref][pop],
                        1: data_freqs[sub][pop_asso[ref][pop][sub]]
                    }


                d += 1
    
    return pop_asso, count_data




def read_windows_SFS(diffs= False, frequency_range= [0,1],indfile= 'ind_assignments.txt', outemp= 'ind_assignments{}.txt',
                    sim_dir= 'mutation_counter/data/sims/', muted_dir= 'mutation_counter/data/mutation_count/',
                    outlog= 'indy.log', row= 24,col= 4, single= True, exclude= False,args= True):
    '''
    
    '''
    
    ti= time.time()
    sims= process_dir(sims_dir= sim_dir)
    print(len(sims))
    tags= []
    sim_extend= []
    chroms= []
    
    data_kmer= {}
    data= {}
    for sim in sims:
        
        ## chromosome
        chrom= sim.split('.')[0].split('C')[-1].strip('chr')
        chromosomes= [sim.split('.')[0].split('C')[1]]
        chromosome_groups = [chromosomes]

        if exclude:
            files= read_exclude()
        else:
            files= {}

        ### read vcf

        vcf_dir= sim_dir + sim + '/'
        vcf_file= vcf_dir + sim + '_' + 'chr' + chrom + '.vcf.gz'
        
        t0= time.time()
        genotype, summary, Names= read_vcf_allel(vcf_file)
        t1= time.time()
        
        read_time= t1- t0

        if len(genotype) == 0:
            continue
        
        #print(genotype.shape, sim)
        
        ## read fasta
        fasta_file= vcf_dir + 'chr{}_{}.fa.gz'.format(chrom,sim)

        with gzip.open(fasta_file,'r') as f:
            lines= f.readlines()
            lines= [x.decode() for x in lines]

        refseq= lines[1].strip()

        ### 
        positions= [int(x) for x in summary.POS]
        wstart= int(min(positions))
        wend= int(max(positions))
        
        Wlen= wend - wstart
        ksize= 3 # odd.
        bases = 'ACGT'
        collapsed= True
        ploidy= 2
        
        genotype_parse= [x for x in range(summary.shape[0]) if int(summary.POS[x])-1 >= wstart and int(summary.POS[x])-1 <= wend]
        Window= genotype[:,genotype_parse]
        subset_summary= summary.loc[genotype_parse,:].reset_index()
        
        ##
        t0= time.time()
        mut_matrix, flag_reverse, flag_remove= vcf_muts_matrix_v1(refseq,subset_summary,start= wstart,end= wend,ksize= ksize,
                                                                bases=bases, collapse= collapsed)

        retain= [x for x in range(Window.shape[1]) if x not in flag_remove]
        Window= Window[:,retain]
        subset_summary= subset_summary.loc[retain,:].reset_index()

        t1= time.time()
        time_mut= t1 - t0

        if diffs:
            sim_start= sim.split('.')[-1]
            diff_snps= read_diffs(sim,diff_dir= vcf_dir, start= int(sim_start))

            summary_diff= [x for x in range(subset_summary.shape[0]) if subset_summary.POS[x] in diff_snps.keys()]

            flag_reverse.extend(summary_diff)
            flag_reverse= list(set(flag_reverse))
        
        
        if flag_reverse:
            Window[:,flag_reverse]= ploidy - Window[:,flag_reverse]
        
        
        pop_dict, pop_freqs= ind_assignment_SFS(sim,Window,dir_sim= sim_dir,indfile= indfile)
        #print(tag_list)
        total_inds= sum([len(x) for x in pop_dict.values()])
        if Window.shape[0] < total_inds:
            continue
        ## counts for no tag sim:
        s0= time.time()
        
        data_kmer[sim]= count_popKmers(Window, mut_matrix, pop_dict, single= single, 
                                  frequency_range= frequency_range,row=row,col=col)
        
        data[sim]= {
            'freqs': pop_freqs,
            'inds': pop_dict,
            'geno': Window,
            'summary': subset_summary
        }

        if args:
            data[sim]['args']= read_args(sim,vcf_dir)


    
    tf= time.time()
    time_elapsed= tf - ti
    
    print('time elapsed: {}s'.format(time_elapsed))
    
    return data_kmer, data



def read_args(reference,sim_dir= ''):

    filename= sim_dir + reference + '_args.txt'

    with open(filename,'r') as fp:
        args= fp.readlines()

    args= [x.split() for x in args]
    args= {x[0]:x[1] for x in args}

    return args


def ind_assignment_SFS(reference,Window, dir_sim= '',indfile= 'ind_assignments.txt'):
    '''
    read ind assignments for a given window; 
    chose one population;
    subset that pop in some way.
    - v1: instead of writting new pop_assignment files, return them. 
    '''
    
    ind_assignments= dir_sim + reference + '/' + indfile
    
    with open(ind_assignments,'r') as f:
        inds= f.readlines()
    
    inds= [x.split() for x in inds]
    pops= np.array(inds)[:,1]
    pop_dict= {
        z: [x for x in range(len(pops)) if pops[x] == z] for z in list(set(pops))
    }
    
    pop_freqs= {}
    
    for pop in pop_dict.keys():
        pop_gen= Window[pop_dict[pop],:]
        freqs= np.sum(pop_gen,axis= 0) / (2*pop_gen.shape[0])
        pop_freqs[pop]= freqs
        ## discount alleles outside freq range.

    return pop_dict, pop_freqs


def array_to_dict(vector):
    '''convert vector to set dictionary to save space'''

    set_dict= {
        z: [x for x in range(len(vector)) if vector[x] == z] for z in list(set(vector))
    }

    return set_dict


def pop_dict_SFS(Window, pop_dict, ploidy= 2):
    '''
    read allele frequencies.
    '''
    
    pop_freqs= {}
    
    for pop in pop_dict.keys():
        pop_gen= Window[pop_dict[pop],:]
        pop_gen= np.sum(pop_gen,axis= 0)
        freq_dict= array_to_dict(pop_gen)
        freq_dict= [(z,len(g)) for z,g in freq_dict.items()]

        pop_freqs[pop]= freq_dict
        
    return pop_freqs


