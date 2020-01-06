import numpy as np
import os

from tools.SLiM_pipe_tools import (
    write_fastaEx, write_popIDs
    )

def cook_constants_v1(fasta_dict, dir_data= "./data/sims/", 
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
                "s1": 216,
                "s2": 108,
                "s3": 206
            }
            
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
    cookID= 'burnin_test'
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
                    "//grow": "{}: ".format(evt) + "{",
                    "//sample": "{} late() ".format(st) + "{"
                }
            }

            with open(SIM_dir + '/' + SIMname + "_args.txt",'w') as f:
                for v,g in sim_store[SIMname].items():
                    if v != "other":
                        f.write("\t".join([str(v),str(g)]) + '\n')

            ### population identifiers file
            sample_sizes= [sim_store[SIMname][x] for x in ["s1"]]
            write_popIDs(sample_sizes,file_dir= SIM_dir)
            
            d += 1
    
    return sim_store, cookID


