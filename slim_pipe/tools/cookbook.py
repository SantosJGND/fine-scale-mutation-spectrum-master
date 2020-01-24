import numpy as np
import os
from itertools import product
import collections
from functools import reduce 
import operator
def recursively_default_dict():
    return collections.defaultdict(recursively_default_dict)

from tools.SLiM_pipe_tools import (
    write_fastaEx, write_popIDs
    )


def write_args(args_dict,SIMname,SIM_dir):
    '''
    *args.txt file stores recipe arguments.
    '''
    with open(SIM_dir + '/' + SIMname + "_args.txt",'w') as f:
        for v,g in args_dict.items():
            if v != "other":
                f.write("\t".join([str(v),str(g)]) + '\n')


def cook_constants_v1(fasta_dict, s1= 216,s2=108,s3=206,dir_data= "./data/sims/", 
            dir_vcf= "vcf_data/", slim_dir= './', batch_name= ''):
    '''
    set up conditions.
    constants:
        - vcf_file;
        - fasta_file - writes fasta; mkdir fasta_dir
        - sampling: s1= 216;s2=108;s3=206
    '''
    cookID= 'v1'
    sim_store= {}
    
    for chrom in fasta_dict.keys():
        for start in fasta_dict[chrom].keys():
            fasta= fasta_dict[chrom][start]

            ### set up names and directories.
            SIMname= batch_name + 'C{}.{}'.format(chrom,str(start))
            SIM_dir= dir_data + SIMname + '/'
            
            os.makedirs(SIM_dir, exist_ok=True)
            vcf_file= SIM_dir + SIMname + "_chr{}.vcf".format(chrom)
            #ref_dir= SIM_dir + SIMname + '_reference'
               
            ### write fasta file for SLiM.
            fasta_file= write_fastaEx(fasta,chrom=chrom,start= start,
                          ID= SIMname,fasta_dir= SIM_dir)
            
            sim_store[SIMname]= {
                "vcf_file": vcf_file,
                "fasta_file": fasta_file,
                "s1": s1,
                "s2": s2,
                "s3": s3
            }
            
            ### write arguments to file
            write_args(sim_store[SIMname],SIMname,SIM_dir)
            ### population identifiers file
            sample_sizes= [sim_store[SIMname][x] for x in ["s1","s2","s3"]]
            write_popIDs(sample_sizes,file_dir= SIM_dir)
            
    
    return sim_store, cookID



def cook_constants_Gravel2sampleRange(fasta_dict, nrange= [.05,.5], step= 10,
            Nmax= 100, dir_data= "./data/", dir_vcf= "vcf_data/sims/", 
            slim_dir= './', batch_name= ''):
    '''
    set up conditions.
    constants:
        - vcf_file;
        - fasta_file - writes fasta; mkdir fasta_dir
        - sampling: 
            s1 (int) to vary in range=nrange as proportion of Nmax;
            s2= Nmax - s1;
            s3= 0
    '''
    cookID= 'Gravel2sampleRange'
    sim_store= {}
    
    s1range= np.linspace(nrange[0],nrange[1],step) * Nmax
    s1range= np.array(s1range,dtype=int)
    s2range= Nmax - s1range
    s3= 0
    
    d= 0
    
    for chrom in fasta_dict.keys():
        for start in fasta_dict[chrom].keys():
            fasta= fasta_dict[chrom][start]

            ### set up names and directories.
            SIMname= batch_name + 'C{}.{}'.format(chrom,str(start))
            SIM_dir= dir_data + SIMname
            
            os.makedirs(SIM_dir, exist_ok=True)
            
            vcf_file= SIM_dir + '/' + SIMname + "_chr{}.vcf".format(chrom)
            #ref_dir= SIM_dir + SIMname + '_reference'
               
            ### write fasta file for SLiM.
            fasta_file= write_fastaEx(fasta,chrom=chrom,start= start,
                          ID= SIMname,fasta_dir= SIM_dir)
            
            sim_store[SIMname]= {
                "vcf_file": vcf_file,
                "fasta_file": fasta_file,
                "s1": s1range[d],
                "s2": s2range[d],
                "s3": s3
            }
            
            ### write arguments to file
            write_args(sim_store[SIMname],SIMname,SIM_dir)
            ### population identifiers file
            sample_sizes= [sim_store[SIMname][x] for x in ["s1","s2","s3"]]
            write_popIDs(sample_sizes,file_dir= SIM_dir)
            
            d += 1
    
    return sim_store, cookID




def cook_constants_simple2split(fasta_dict, nrange= [.05,.5], step= 10,
            Nmax= 100, dir_data= "./data/", dir_vcf= "vcf_data/sims/", 
            slim_dir= './', batch_name= ''):
    '''
    set up conditions.
    constants:
        - vcf_file;
        - fasta_file - writes fasta; mkdir fasta_dir
        - sampling: 
            s1 (int) to vary in range=nrange as proportion of Nmax;
            s2= Nmax - s1;
            s3= 0
    '''
    cookID= 'simple2split'
    sim_store= {}
    
    s1range= np.linspace(nrange[0],nrange[1],step) * Nmax
    s1range= np.array(s1range,dtype=int)
    s2range= Nmax - s1range
    
    d= 0
    
    for chrom in fasta_dict.keys():
        for start in fasta_dict[chrom].keys():
            fasta= fasta_dict[chrom][start]

            ### set up names and directories.
            SIMname= batch_name + 'C{}.{}'.format(chrom,str(start))
            SIM_dir= dir_data + SIMname
            
            os.makedirs(SIM_dir, exist_ok=True)
            
            vcf_file= SIM_dir + '/' + SIMname + "_chr{}.vcf".format(chrom)
            #ref_dir= SIM_dir + SIMname + '_reference'
               
            ### write fasta file for SLiM.
            fasta_file= write_fastaEx(fasta,chrom=chrom,start= start,
                          ID= SIMname,fasta_dir= SIM_dir)
            
            sim_store[SIMname]= {
                "vcf_file": vcf_file,
                "fasta_file": fasta_file,
                "s1": s1range[d],
                "s2": s2range[d],
            }
            
            ### write arguments to file
            write_args(sim_store[SIMname],SIMname,SIM_dir)
            ### population identifiers file
            sample_sizes= [sim_store[SIMname][x] for x in ["s1","s2"]]
            write_popIDs(sample_sizes,file_dir= SIM_dir)
            
            d += 1
    
    return sim_store, cookID



def cook_constants_SizeChange(fasta_dict, s1= 1092, NeC= 2e5, Nef= 4e5, Grate= 1.03,
            dir_data= "./data/", dir_vcf= "vcf_data/sims/", 
            slim_dir= './', batch_name= ''):
    '''
    set up conditions.
    constants:
        - vcf_file;
        - fasta_file - writes fasta; mkdir fasta_dir
        - sampling: 
            s1 (int) to vary in range=nrange as proportion of Nmax;
        - NeC: initial population eff. size.
        - Nef: effective population size after change.
        - Grate: growth rate during change. 
    '''
    cookID= 'sizeChange'
    sim_store= {}

    for chrom in fasta_dict.keys():
        for start in fasta_dict[chrom].keys():
            fasta= fasta_dict[chrom][start]

            ### set up names and directories.
            SIMname= batch_name + 'C{}.{}'.format(chrom,str(start))
            SIM_dir= dir_data + SIMname
            
            os.makedirs(SIM_dir, exist_ok=True)
            
            vcf_file= SIM_dir + '/' + SIMname + "_chr{}.vcf".format(chrom)
            #ref_dir= SIM_dir + SIMname + '_reference'
               
            ### write fasta file for SLiM.
            fasta_file= write_fastaEx(fasta,chrom=chrom,start= start,
                          ID= SIMname,fasta_dir= SIM_dir)
            
            sim_store[SIMname]= {
                "vcf_file": vcf_file,
                "fasta_file": fasta_file,
                "s1": s1,
                "NeC": NeC,
                "Nef": Nef,
                "Grate": Grate,
            }

            ### write arguments to file
            write_args(sim_store[SIMname],SIMname,SIM_dir)            
            ### population identifiers file
            sample_sizes= [sim_store[SIMname][x] for x in ["s1"]]
            write_popIDs(sample_sizes,file_dir= SIM_dir)
    
    return sim_store, cookID




def cook_constants_Burnin(fasta_dict, bt= 50000, sdelay= 1000,s1= 1092, NeC= 2e5, Nef= 4e5, Grate= 1.03,
            dir_data= "./data/", dir_vcf= "vcf_data/sims/", 
            slim_dir= './', batch_name= ''):
    '''
    set up conditions.
    constants:
        - vcf_file;
        - fasta_file - writes fasta; mkdir fasta_dir
        - sampling: 
            s1 (int) to vary in range=nrange as proportion of Nmax;
        - NeC: initial population eff. size.
        - Nef: effective population size after change.
        - Grate: growth rate during change. 
    '''
    cookID= 'burnin'
    sim_store= {}

    possible= sum([len(fasta_dict[x]) for x in fasta_dict.keys()])
    burnin_list= np.linspace(1,bt,possible,dtype= int)
    d= 0

    for chrom in fasta_dict.keys():
        for start in fasta_dict[chrom].keys():
            fasta= fasta_dict[chrom][start]

            # burnin time
            evt= burnin_list[d]
            st= burnin_list[d] + sdelay

            ### set up names and directories.
            SIMname= batch_name + 'T{}'.format(evt) + 'C{}.{}'.format(chrom,str(start))
            SIM_dir= dir_data + SIMname
            
            os.makedirs(SIM_dir, exist_ok=True)
            
            vcf_file= SIM_dir + '/' + SIMname + "_chr{}.vcf".format(chrom)
            #ref_dir= SIM_dir + SIMname + '_reference'
            
            ### write fasta file for SLiM.
            fasta_file= write_fastaEx(fasta,chrom=chrom,start= start,
                          ID= SIMname,fasta_dir= SIM_dir)
            
            sim_store[SIMname]= {
                "vcf_file": vcf_file,
                "fasta_file": fasta_file,
                "s1": s1,
                "NeC": NeC,
                "Nef": Nef,
                "Grate": Grate,
                "evt": evt,
                "other": {
                    "//grow": "{}: ".format(evt) + "{\n",
                    "//sample": "{} late() ".format(st) + "{\n"
                }
            }
            
            ### write arguments to file
            write_args(sim_store[SIMname],SIMname,SIM_dir)
            ### population identifiers file
            sample_sizes= [sim_store[SIMname][x] for x in ["s1"]]
            write_popIDs(sample_sizes,file_dir= SIM_dir)
            
            d += 1
    
    return sim_store, cookID



##### mutation rate var.
#####
def cook_constants_rateVar(fasta_dict, mu= 1e-8, mut_file= 'mut_matrix_v0.txt', s1= 1092, NeC= 2e5, Nef= 4e5, Grate= 1.03,
            dir_data= "./data/", dir_vcf= "vcf_data/sims/", 
            slim_dir= './', batch_name= ''):
    '''
    set up conditions.
    constants:
        - vcf_file;
        - fasta_file - writes fasta; mkdir fasta_dir
        - sampling: 
            s1 (int) to vary in range=nrange as proportion of Nmax;
        - NeC: initial population eff. size.
        - Nef: effective population size after change.
        - Grate: growth rate during change. 
    '''
    cookID= 'rateVar'
    sim_store= {}

    for chrom in fasta_dict.keys():
        for start in fasta_dict[chrom].keys():
            fasta= fasta_dict[chrom][start]

            ### set up names and directories.
            SIMname= batch_name + 'C{}.{}'.format(chrom,str(start))
            SIM_dir= dir_data + SIMname
            
            os.makedirs(SIM_dir, exist_ok=True)
            
            vcf_file= SIM_dir + '/' + SIMname + "_chr{}.vcf".format(chrom)
            #ref_dir= SIM_dir + SIMname + '_reference'
               
            ### write fasta file for SLiM.
            fasta_file= write_fastaEx(fasta,chrom=chrom,start= start,
                          ID= SIMname,fasta_dir= SIM_dir)
            
            sim_store[SIMname]= {
                "vcf_file": vcf_file,
                "fasta_file": fasta_file,
                "s1": s1,
                "NeC": NeC,
                "Nef": Nef,
                "Grate": Grate,
                "mu": mu,
                "mut_file": mut_file
                #"other": {'//mut_file': '\tfile_mut= readFile("{}");\n'.format(mut_file)}  
            }

            ### write arguments to file
            write_args(sim_store[SIMname],SIMname,SIM_dir)            
            ### population identifiers file
            sample_sizes= [sim_store[SIMname][x] for x in ["s1"]]
            write_popIDs(sample_sizes,file_dir= SIM_dir)
    
    return sim_store, cookID


def get_by_path(root, items):
    """Access a nested object in root by item sequence."""
    return reduce(operator.getitem, items, root)


def mutation_dict_full(bases= 'ATCG',ksize= 3):
    mutations= []
    mut_lib= recursively_default_dict()
    mut_org= []
    
    base_set= [bases]*ksize

    for trimer in product(*base_set):
        mut_org.append(trimer)
        for base in bases:
            mutations.append((''.join(trimer), base))
            get_by_path(mut_lib, trimer[:-1])[trimer[-1]]= base
            
            
    
    return mut_lib,mutations,mut_org



def rate_mods(mut_org,rate_range= [1,5],rate_change= 10, bases= 'ACGT',mu= 1e-8):
    
    Nr_sim= np.random.randint(rate_range[0],high=rate_range[1],size=1)[0]

    which_muts= np.random.randint(0,high= len(mut_org) * 3,size= Nr_sim)

    mut_coords= []
    for cell in which_muts: 
        row= int(cell / len(bases)-1)
        col= cell % (len(bases)-1)

        mut_coords.append((row,col))

    rate_dict= {
        z: [mut_coords[x][1] for x in range(len(mut_coords)) if mut_coords[x][0] == z] for z in list(set([x[0] for x in mut_coords]))
    }

    mut_dict= {

    }

    for row in rate_dict.keys():
        mut= mut_org[row]
        idx_poss= [x for x in range(len(bases)) if x != mut[1]]
        rates= [[0,mu/3][int(bases[x] != mut[1])] for x in range(len(bases))]

        for col in rate_dict[row]:
            coord= (row,col)                                                                                                  
            idx_change= idx_poss[col]

            rates[idx_change]= rates[idx_change] * rate_change

        rates= np.array(rates) * (mu / sum(rates))

        mut_dict[''.join(mut)]= rates


    return mut_dict


def cook_constants_rateVarMat(fasta_dict, mu= 1e-8, bases= 'ACGT', rate_change= 10, rate_range= [1,5],Nmat= 5, s1= 2000, NeC= 2e5, 
            Nef= 4e5, Grate= 1.03, dir_data= "./data/", dir_vcf= "vcf_data/sims/", 
            slim_dir= './', batch_name= ''):
    '''
    set up conditions.
    constants:
        - vcf_file;
        - fasta_file - writes fasta; mkdir fasta_dir
        - sampling: 
            s1 (int) to vary in range=nrange as proportion of Nmax;
        - NeC: initial population eff. size.
        - Nef: effective population size after change.
        - Grate: growth rate during change. 
    '''
    cookID= 'rateVarII'
    sim_store= {}
    mutations_full_dict, mutations_full_list, mut_org= mutation_dict_full(bases= bases)
    var_store= {
        "M{}".format(x): rate_mods(mut_org,rate_range= rate_range,rate_change= rate_change, bases= bases,mu= mu) for x in range(1,Nmat+1)
    }
    var_store["M0"]= {}
    mat_names= {mat: batch_name + mat + '_grid.txt' for mat in var_store.keys()}

    for mat in var_store.keys():
        with open(mat_names[mat],'w') as fp:
            for mut in var_store[mat].keys():
                rates= var_store[mat][mut]
                rates= ','.join([str(x) for x in rates])
                fp.write('\t'.join([mut,rates]) + '\n')

    for chrom in fasta_dict.keys():
        for start in fasta_dict[chrom].keys():
            fasta= fasta_dict[chrom][start]

            for mat in var_store.keys():
            ### set up names and directories.
                SIMname= batch_name + mat + 'C{}.{}'.format(chrom,str(start))
                SIM_dir= dir_data + SIMname
                
                os.makedirs(SIM_dir, exist_ok=True)
                
                vcf_file= SIM_dir + '/' + SIMname + "_chr{}.vcf".format(chrom)
                #ref_dir= SIM_dir + SIMname + '_reference'
                   
                ### write fasta file for SLiM.
                fasta_file= write_fastaEx(fasta,chrom=chrom,start= start,
                              ID= SIMname,fasta_dir= SIM_dir)
                
                sim_store[SIMname]= {
                    "vcf_file": vcf_file,
                    "fasta_file": fasta_file,
                    "s1": s1,
                    "NeC": NeC,
                    "Nef": Nef,
                    "Grate": Grate,
                    "mu": mu,
                    "mut_file": mat_names[mat]
                    #"other": {'//mut_file': '\tfile_mut= readFile("{}");\n'.format(mut_file)}  
                }

                ### write arguments to file
                write_args(sim_store[SIMname],SIMname,SIM_dir)            
                ### population identifiers file
                sample_sizes= [sim_store[SIMname][x] for x in ["s1"]]
                write_popIDs(sample_sizes,file_dir= SIM_dir)
    
    return sim_store, cookID
