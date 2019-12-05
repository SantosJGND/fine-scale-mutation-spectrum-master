import sys
import argparse
import numpy as np

from itertools import product
import itertools as it

import tempfile
import os
import gzip 
import subprocess

from datetime import datetime




def read_chrom_sizes(assembly,size_dir= 'chrom_sizes/',reject= ['M','Un','X','Y','_']):
    '''
    read chromosome size file. Store and return as {str(chrom): int(chrom_len)} dict.
    '''
    filename= size_dir + assembly + '.chrom.sizes'
    with open(filename,'r') as f:
        lines= f.readlines()
        lines= [line.strip().split() for line in lines]
        
        Sizes= {x[0].strip('chr'): int(x[1]) for x in lines if sum([int(y in x[0]) for y in reject]) == 0}
    
    return Sizes



def region_sample(L, chrom_sizes, N, fasta_file= ''):
    '''
    deprecated.
    extract sequences from fasta, by chrom.
    '''
    chroms= np.random.choice(list(chrom_sizes.keys()),N,replace= True)
    chroms= sorted(chroms)
    chroms= {z: chroms.count(z) for z in list(set(chroms))}
    
    rseqs= {z: {} for z in chroms.keys()}
    
    for Chr in chroms.keys():
        seq_dict= fasta_chrExtract(fasta_file,[Chr])
        
        starts= np.random.randint(0,high= chrom_sizes[Chr] - L,size= chroms[Chr])
        
        for s in starts:
            rseqs[Chr][s]= seq_dict[Chr][s:(s+L)]
    
    return rseqs



def region_samplev1(L, chrom_sizes, N, fasta_file= ''):
    '''
    prepare sequence dictionary: {chrom:list(start pos)}
    provide to fasta_Rextract.
    
    ## modify this function for more pertinent selections.
    '''
    chroms= np.random.choice(list(chrom_sizes.keys()),N,replace= True)
    chroms= sorted(chroms)
    chroms= {z: chroms.count(z) for z in list(set(chroms))}
    
    seqs_store= {
        z: np.random.randint(low= 0,high= chrom_sizes[z] - L,size= chroms[z]) for z in chroms.keys()
    }
    
    seqs= fasta_Rextract(fasta_file, seqs_store, L= L)
    
    
    return seqs



def fasta_Rextract(fasta,seq_store,L= 10000):
    ''' 
    Extract regions from fasta file.
    Takes dictionary {chrom: list(start points)};
    extracts sequences after reading each chromosome.
    '''
    refseqs= {x:{} for x in seq_store.keys()}
    
    Outfiles= {
        x: tempfile.TemporaryFile() for x in seq_store.keys()
    }
    
    d= 0
     
    with gzip.open(fasta,'rb') as fp:
        for line in fp:
            #line= line.decode()
            if line[0:1] == b'>':
                head= line.decode()
                head= head.strip()
                if d != 0:
                    d=0

                head= head[1:].split('\t')
                head= head[0].strip('chr')

                if head in seq_store.keys():
                    print('opening fasta chr: {}'.format(head))
                    d= head
                    outfile= Outfiles[d]
                    continue
            
            if d != 0:
                processed_line= line.upper().strip(b'\n')
                outfile.write(processed_line)
    
    for chrom in Outfiles.keys():
        f= Outfiles[chrom]
        f.seek(os.SEEK_SET)
        result= f.read().decode()
        f.close()
        
        starts= seq_store[chrom]
        
        for s in starts:
            refseqs[chrom][s]= result[s:(s+L)]
    
    return refseqs


def write_fastaEx(fasta,chrom= '1',start= 0, ID= 'SIM', fasta_dir= ''):
    ''' write extracted sequence as fasta file'''
    filename= ''.join([fasta_dir,'/chr',chrom,'_',ID,'.fa'])
    
    header= ' '.join(['>'+ID,chrom,str(start),str(start + len(fasta))]) + '\n'
    
    with open(filename,'w') as fp:
        fp.write(header)
        fp.write(fasta)
    
    return filename


def process_recipe(recipe,constant_dict, SIMname):
    '''add constant definitions to a SLiM recipe'''
    
    new_recipe= recipe.split('/')
    new_recipe[-1]= SIMname + '_' + new_recipe[-1]
    new_recipe= '/'.join(new_recipe)
    
    with open(recipe,'r') as f:
        lines= f.readlines()
    
    init_idx= [x for x in range(len(lines)) if 'initialize()' in lines[x]][0]
    defConst= '\tdefineConstant({},{});\n'
    
    for v,g in constant_dict.items():
        newline= defConst.format('"{}"'.format(v),'"{}"'.format(g))
        
        lines.insert(init_idx+1,newline)
    
    with open(new_recipe,'w') as f:
        f.write(''.join(lines))
    
    return new_recipe



def SLiM_dispenser(fasta_dict, sim_recipe, dir_data= "./data/", 
	dir_vcf= "vcf_data/", slim_dir= './', logSims= 'sims.log',
	mutlog= 'toMut.log',batch_name= ''):
    ''' execute SLiM program
    - simulation specific recipe:
    - recipe template is re-written to direct to new fasta.
    '''
    for chrom in fasta_dict.keys():
        for start in fasta_dict[chrom].keys():
            fasta= fasta_dict[chrom][start]
            
            ## SLiM takes only 'ACTG'.
            if len([x for x in fasta if x in 'ATCG']) < len(fasta):
                continue

            ### set up names and directories.
            SIMname= batch_name + 'C{}.{}'.format(chrom,str(start))
            vcf_file= dir_vcf + SIMname + "_chr{}.vcf".format(chrom)
            ref_dir= dir_data + SIMname + '_reference'
            os.makedirs(ref_dir, exist_ok=True)
            
            ### write fasta file for SLiM.
            fasta_file= write_fastaEx(fasta,chrom=chrom,start= start,
                          ID= SIMname,fasta_dir= ref_dir)

            command_line_constants= {
                "vcf_file": vcf_file,
                "fasta_file": fasta_file
            }
            
            ### generate modified slim recipe
            new_recipe= process_recipe(sim_recipe,command_line_constants, SIMname)
            
            ### Launch SLiM through shell.
            slim_soft= slim_dir + 'slim' 
            seed= np.random.randint(0,high= start,size= 1)[0]

            command_units= [slim_soft, '-m', '-s', str(seed), new_recipe]
            
            os.system(' '.join(command_units))
            
            os.system('gzip {}'.format(vcf_file))
            os.system('gzip {}'.format(fasta_file))
            os.system('rm {}'.format(new_recipe))
            
            now = datetime.now()
            tnow= now.strftime("%d/%m/%Y %H:%M:%S")
            
            with open(logSims,'a') as fp:
                fp.write(SIMname + '\t' + tnow + '\n')
            with open(mutlog,'a') as fp:
                fp.write(SIMname + '\n')



def mutation_counter_launch(logfile,count_dir= './count/', dir_launch= '..',main_dir= './'):
    '''
    launch mutation counter.
    - read mut.log to know which have not been yet processed.
    - launch process_chromosomes.py using simulation name. 
    '''
    with open(logfile,'r') as fp:
        sims= fp.readlines()
    
    print(sims)
    sims= [x.strip() for x in sims]
    chroms= [x.split('.')[0].split('C')[-1].strip('chr') for x in sims]
    
    job= 'python process_chromosomes.py -c {} -r {} -s {} -v {}_ -d {} -q vcf_data'
    
    sims= [job.format(chroms[x],*[sims[x]]*3,dir_launch) for x in range(len(sims))]
    
    os.chdir(count_dir)
    for sim in sims:
        
        os.system(sim)
    
    os.chdir(main_dir)

    open(logfile,'w').close()
    



#####
#####
def fasta_chrExtract(fasta,chromosomes):
    '''
    extract selected chromosomes all in one go. return as {chrom: fasta} dict.

    '''
    Outfiles= {
        x: tempfile.TemporaryFile() for x in chromosomes
    }
    seqs= {}
    
    d= 0
    
    with gzip.open(fasta,'rb') as fp:
        for line in fp:
            #line= line.decode()
            if line[0:1] == b'>':
                head= line.decode()
                head= head.strip()
                if d != 0:
                    d=0

                head= head[1:].split('\t')
                head= head[0].strip('chr')

                if head in chromosomes:
                    print('opening fasta chr: {}'.format(head))
                    d= head
                    outfile= Outfiles[d]
            
            if d != 0:
                processed_line= line.upper().strip(b'\n')
                outfile.write(processed_line)
    
    for chrom in Outfiles.keys():
        f= Outfiles[chrom]
        f.seek(os.SEEK_SET)
        result= f.read()
        f.close()
        seqs[chrom]= result
    
    return seqs

