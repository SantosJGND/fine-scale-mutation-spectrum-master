import numpy as numpy
from itertools import product
import os
from functools import reduce  # forward compatibility for Python 3
import operator
import gzip
import tempfile
import collections

def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)



def reference_sequence(chromosome_number,reference= 'rheMac10',dir_launch='..'):

    file_path_template = dir_launch+'chr{}_{}.fa.gz'.format(chromosome_number,reference)
    
    with gzip.open(file_path_template) as infile:
        lines = infile.readlines()
    
    infile.close()
    
    f = tempfile.TemporaryFile()
    
    for v in range(1,len(lines)):
        processed_line= lines[v].upper().strip(b'\n')

        f.write(processed_line)
    
    f.seek(os.SEEK_SET)
    result= f.read()
    
    f.close()
    
    return result

####
####
def get_mutations(bases= 'ACGT',ksize= 3):
    '''return list of possible kmer mutations'''
    
    mutations=[]
    
    
    base_set= [bases]*ksize

    for trimer in product(*base_set):
        for base in bases:
            if trimer[int(ksize / 2)] != base:
                mutations.append((''.join(trimer), base))
    
    return mutations

####
def get_by_path(root, items):
    """Access a nested object in root by item sequence."""
    return reduce(operator.getitem, items, root)

def set_by_path(root, items, value):
    """Set a value in a nested object in root by item sequence."""
    get_by_path(root, items[:-1])[items[-1]] = value

####
####

def kmer_dict_init(ksize= 3,bases='ATCG'):
    '''produce nested dictionary of nucs for a particular kmer size'''
    mut_lib= recursively_default_dict()

    base_set= [bases]*ksize

    for trimer in product(*base_set):
        set_by_path(mut_lib,list(trimer),0)
    
    return mut_lib


def fasta_get_freq(seq,start= 0,end= 0,step= 1,ksize=3,bases= 'ATCG'):
    '''return count of kmer across fasta region'''
    kmer_dict= kmer_dict_init(ksize= ksize,bases=bases)
    if end == 0:
        end= len(seq) - ksize
    
    for ki in range(start,end-ksize,step):
        kmer= seq[ki:ki+ksize]
        if 'N' in kmer:
            continue
        get_by_path(kmer_dict, kmer[:-1])[kmer[-1]] += 1
    
    return kmer_dict



def get_complement(kmer):
    '''Return complement of a given kmer'''
    complements= {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    
    comp= [complements[x] for x in kmer][::-1]
    return comp


def complement_dicts(bases= 'ACGT',ksize= 3):
    '''return dict of comp kmers + index dict to parse with'''
    comp_dict= {}
    comp_index= {}
    d= 0
    base_set= [bases] * ksize
    mers= product(*base_set)
    
    for kmer in mers:
        kmer= ''.join(kmer)
        if kmer not in comp_dict.keys():
            comp_dict[kmer]= d
            
            comp= get_complement(kmer)
            comp= ''.join(comp)
            comp_dict[comp]= d
            
            comp_index[d]= (kmer,comp)
            d += 1
    
    return comp_dict,comp_index


def collapse_freqs(kmer_dict,comp_index):
    '''return vector of collapsed counts by kmer'''
    counts= []
    
    for kdx in comp_index.keys():
        total= [get_by_path(kmer_dict, list(comp)) for comp in comp_index[kdx]]
        total= sum(total)
        counts.append(total)
        
    return counts
            

