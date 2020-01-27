
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


fig_dir= 'Figures/MutVar'
os.makedirs(fig_dir, exist_ok=True)
fig_dir= fig_dir + '/'


mutlog= 'toMut.log'
min_size= 70
sampling= [5,100,10]
bases= 'ACGT'
collapsed= False
row= 48
col= 4
ksize= 3

sample_sim= 150

data, data_freqs = MC_sample_matrix_v1(min_size= min_size, samp= sampling, count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir,sim_dir= sims_dir,
                          muted_dir= muted_dir, diffs= diffs, row= row,bases= bases,
                       exclude= False,sample_sim= sample_sim,collapsed= collapsed)


print("data extracted: data= {} elements".format(len(data)))

from tools.fasta_utilities import (
    get_mutations, get_by_path, kmer_comp_index, kmer_mut_index,
    fasta_get_freq,kmer_dict_init
    )


def fasta_get_freq(seq,start= 0,end= 0,step= 1,ksize=3,bases= 'ATCG'):
    '''return count of kmer across fasta region'''
    kmer_dict= kmer_dict_init(ksize= ksize,bases=bases)
    if end == 0:
        end= len(seq) - ksize

    for ki in range(start,end,step):
        kmer= seq[ki:ki+ksize]
        if 'N' in kmer:
            continue
        get_by_path(kmer_dict, kmer[:-1])[kmer[-1]] += 1

    return kmer_dict


def kmer_freq_balance(kmer_dict, mutations, fasta_len= 10000, bases= 'ACGT',ksize= 3):
    '''return list of possible kmer mutations'''

    mutation_sum= []
    Nkmers= fasta_len - ksize

    for idx in range(len(mutations)):
        mut= mutations[idx]
        prop= get_by_path(kmer_dict,mut[0])
        prop= prop / Nkmers
        mutation_sum.append(prop)

    return np.array(mutation_sum).reshape(1,-1)


def get_fasta_prop(sim,sim_dir,mutations,ksize= 3,bases= 'ACGT'):
    
    chrom= sim.split('.')[0].split('C')[1]
    
    fasta_file= sim_dir + 'chr{}_{}.fa.gz'.format(chrom,sim)

    with gzip.open(fasta_file,'r') as f:
        lines= f.readlines()
        lines= [x.decode() for x in lines]

    refseq= lines[1].strip()

    kmer_dict= fasta_get_freq(refseq,start= 0,end= 0,step= 1,ksize=ksize,bases= bases)

    ref_kmer_prop =kmer_freq_balance(kmer_dict,mutations,fasta_len= len(refseq))
    
    return ref_kmer_prop


ksize= 3
#bases= 'ATCG'
bases= 'ACGT'

mutations= get_mutations(bases= bases,ksize= ksize)
kmers, kmer_idx= kmer_comp_index(mutations)

mut_lib= kmer_mut_index(mutations)


from tools.mcounter_tools import read_args

p_value= 1e-5
test_m= 'chi2'
individually= False
exclude= False
frequency_range= [0,1]
extract= 'pval'
Nbins= 100
tag_ref= '_ss'


### separate original and subset mutation count dicts.
#### this should be done in MC_sample_matrix_v1.

avail= list(data.keys())
ref_idx= [int(tag_ref in avail[x]) for x in range(len(avail) )]
categ= {
    z: [x for x in range(len(avail)) if ref_idx[x] == z] for z in [0,1]
}

print([len(categ[x]) for x in [0,1]])

######
###### extract reference mutation count proportions. 

### possible combinations per simulation.
### ref_mdict: store count proportions by mutation matrix (list of lists). PCA,(deprecatted distances )
ref_mdict= recursively_default_dict()
### ref_mats: mutation matrix file (label) for reference simulation. 
ref_mats= {}

### ref_mat_dict: mutation count by population by simulation by mut_matrix file (used for distances)
ref_mat_dict= recursively_default_dict()

### fasta_ref_dict: store kmer frequency by reference simulation (use for standardize count proportions)
fasta_ref_dict= {}
### store fasta frequencies by population for each mutation matrix label - PCA. 
fasta_pop= {}


for idx in categ[0]:
    ref= avail[idx]
    ref_dir= sims_dir + ref + '/'
    ref_args= read_args(ref,sim_dir=ref_dir)
    mut_matrix= ref_args['mut_file']
    
    fasta_kmer_prop= get_fasta_prop(ref,ref_dir,mutations,ksize= 3,bases= 'ACGT')
    fasta_ref_dict[ref]= fasta_kmer_prop
    
    ref_mats[ref]= mut_matrix
    
    batch= ref.split('C')[0]
    
    sizes= data[ref]['sizes']
    #
    
    chromosomes= [ref.split('.')[0].split('C')[1]]

    pop_counts= data[ref]['counts']

    pop_counts= {
    z: pop_counts[z] / np.sum(pop_counts[z]) for z in pop_counts.keys()
    }
    
    if mut_matrix not in ref_mdict.keys():
        ref_mdict[mut_matrix]= []
        fasta_pop[mut_matrix]= []
        
    
    for pop in pop_counts.keys():
        mlist= pop_counts[pop].reshape(1,np.prod(pop_counts[pop].shape))[0]
        mlist= mlist.reshape(1,-1)
        
        ## balance for kmer frequencies in fasta.
        #mlist= mlist * (1/fasta_kmer_prop.shape[1] / fasta_kmer_prop)
        mlist= mlist - (fasta_kmer_prop/3)
        
        fasta_pop[mut_matrix].append(list(fasta_kmer_prop[0]))
        ref_mdict[mut_matrix].append(list(mlist[0]))
        ref_mat_dict[mut_matrix][ref][pop]= mlist
    


## reference frequency and kmer counts as pop x kmer arrays. 
from sklearn import decomposition
from sklearn import preprocessing

mat_avail= list(ref_mdict.keys())

labels= np.repeat(mat_avail,[len(ref_mdict[x]) for x in mat_avail])
mat_dict= {
    z: [x for x in range(len(labels)) if labels[x] == z] for z in list(set(labels))
}
data_ref= list(it.chain(*[ref_mdict[x] for x in mat_avail]))
data_ref= np.array(data_ref)
data_freq= list(it.chain(*[fasta_pop[x] for x in mat_avail]))
data_freq=np.array(data_freq)


#################
#################  processing sub-sampled data sets

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
from sklearn.metrics import pairwise_distances
labels_ref= list(mat_dict.keys())

sub_pop= []
sub_prop= []
sub_label= []
sub_values= []

dub_diffs= []
euc_predict= []


for ref_sim in pop_asso.keys():
    #print(ref_sim)
    batch= ref_sim.split('C')[0]
    mut_matrix= ref_mats[ref_sim]
    # kmer_frequencies
    fasta_kmer_prop= fasta_ref_dict[ref]
    
    for pop in pop_asso[ref_sim].keys():
        for prop in pop_asso[ref_sim][pop].keys():
            
            for sub in pop_asso[ref_sim][pop][prop].keys():
                
                poppy= pop_asso[ref_sim][pop][prop][sub]
                
                pop_counts= data[sub]['counts']
                pop_counts= {z: g / np.sum(g) for z,g in pop_counts.items()}
                if not poppy:
                    continue
                mlist= pop_counts[poppy].reshape(1,np.prod(pop_counts[poppy].shape))[0]

                sub_pop.append(poppy)
                sub_values.append(list(mlist))
                sub_label.append(mut_matrix)
                sub_prop.append(prop)
                mlist=mlist.reshape(1,-1)
                ## balance for kmer frequency in respective fasta.
                mlist= mlist - (fasta_kmer_prop / 3)
                #mlist= mlist / np.sum(mlist)
                
                #ref_dists= [pairwise_distances(mlist,np.array(ref_mdict[z]),metric= 'euclidean') for z in labels_ref]
                #ref_dists= [pairwise_distances(mlist,np.mean(np.array(ref_mdict[z]),axis=0).reshape(1,-1)) for z in labels_ref]
                ref_dists= [[pairwise_distances(mlist,ref_mat_dict[z][r][p],metric= 'euclidean') for r in ref_mat_dict[z].keys() for p in ref_mat_dict[z][r].keys()  if [z,r,p] != [mut_matrix,ref,pop]] for z in labels_ref]
                ref_dists= [np.mean(x) for x in ref_dists]
                
                predict= np.argmin(np.array(ref_dists))
                predict= labels_ref[predict]
                euc_predict.append(predict)
                

print('Done with processing: {} observations.'.format(len(sub_values)))
sub_values= np.array(sub_values)


#####################
##################### Assignment to mutation grids
# Use mutation count arrays from full populations as reference label set for classification.
# implemented. classification by average proximity :
#             distance calculated between subset and ref count props, averaged by mut. matrix
#             classification label selected using minimum distance. subset original pop excluded. 
# available: sklearnn KNN classification. poor results at N=3.

#### KNN trainning
predictions= euc_predict


#### prop extraction.
############
############ compare predictions with known mutation matrix across populations.
############ Store by original mutation matrix, average by proportion bin. 

Nbins= 100
bins= np.linspace(0,1,Nbins)
bins= np.round(bins,4)
bins= [(bins[x-1],bins[x]) for x in range(1,len(bins))]


lab_dict= {
    z: [x for x in range(sub_values.shape[0]) if sub_label[x] == z] for z in list(set(sub_label))
}

## actual props
'''
lab_prop= {
    lab: {
        prop: [x for x in lab_dict[lab] if sub_prop[x] == prop] for prop in list(set(sub_prop))
    } for lab in lab_dict.keys()
}
'''

## bins
lab_prop= {
    lab: {
        sum(bi)/2: [x for x in lab_dict[lab] if sub_prop[x] >= bi[0] and sub_prop[x] < bi[1]] for bi in bins
    } for lab in lab_dict.keys()
}


comp_dict= {
    lab: {
        prop: [int(predictions[x] == sub_label[x]) for x in lab_prop[lab][prop]] for prop in lab_prop[lab].keys()
    } for lab in lab_prop.keys()
}


comp_std= {
    lab: {
        prop: np.std(comp_dict[lab][prop]) for prop in comp_dict[lab].keys()
    } for lab in comp_dict.keys()
}


comp_dict= {
    lab: {
        prop: np.mean(comp_dict[lab][prop]) for prop in comp_dict[lab].keys()
    } for lab in comp_dict.keys()
}



###############
############### number of kmer with mutation rate modifications by simulation.

mat_dir= ''
mat_lines= {}

for lab in comp_dict.keys():
    with open(mat_dir + lab,'r') as fp:
        lines= fp.readlines()
        mat_lines[lab]= len(lines)


###########################
############################ PLOTs
############### plots

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt


############ PCA and kmer frequency in fasta vs. kmer count proportion in full (reference) populations.


freq_array= data_freq.reshape(1,np.prod(data_freq.shape))[0] / 3
ref_array= data_ref.reshape(1,np.prod(data_ref.shape))[0]



plt.figure(figsize=(10, 10))

plt.scatter(freq_array,ref_array)
plt.xlabel("kmer freq. in fasta")
#plt.ylim(0,1.5)
plt.ylabel("kmer count")
plt.title('kmer frequency & count across data sets.')

plt.savefig(fig_dir + "frequency_count.png",bbox_inches='tight')
plt.close()


###############
############## plot accuracy

plt.figure(figsize=(10, 10))

for lab in comp_dict.keys():
    y= [comp_dict[lab][prop] for prop in sorted(comp_dict[lab].keys())]
    x= sorted(comp_dict[lab].keys())
    errory= [comp_std[lab][prop] for prop in sorted(comp_dict[lab].keys())]
    plt.errorbar(x,y,yerr=errory,label= lab.split('.')[0] + '_{}'.format(mat_lines[lab]))

plt.xlabel('relative sample size')
#plt.ylim(0,1.5)
plt.ylabel('accuracy')
plt.title('rate of accuracy. kmer freq controlled')

plt.legend()
plt.savefig(fig_dir + "Accuracy_by_mutMatrix.png",bbox_inches='tight')
plt.close()



####################
################## plot false positives

####################
################## FALSE POSITIVES

fp_dict= {
    lab: {
        prop: list(it.chain(*[[int(predictions[x] == lab and sub_label[x] != lab) for x in lab_prop[olab][prop]] for olab in lab_prop.keys() if olab != lab])) for prop in lab_prop[lab].keys()
    } for lab in lab_prop.keys()
}

fp_std= {
    lab: {
        prop: np.std(fp_dict[lab][prop]) for prop in fp_dict[lab].keys()
    } for lab in fp_dict.keys()
}


fp_dict= {
    lab: {
        prop: np.mean(fp_dict[lab][prop]) for prop in fp_dict[lab].keys()
    } for lab in fp_dict.keys()
}


plt.figure(figsize=(10, 10))

for lab in fp_dict.keys():
    props= sorted(fp_dict[lab].keys())
    print(props)

    y= [fp_dict[lab][prop] for prop in props],
    x= props,
    errory=  [fp_std[lab][prop] for prop in props][0],
    print(errory)
    plt.errorbar(x,y,yerr=errory, label= lab)

plt.xlabel('relative sample size')
#plt.ylim(0,1.5)
plt.ylabel('false positive rate')
plt.title('FP & sampling. kmer freq controlled')

plt.legend()
plt.savefig(fig_dir + "FalsePositives_by_mutMatrix.png",bbox_inches='tight')
plt.close()


