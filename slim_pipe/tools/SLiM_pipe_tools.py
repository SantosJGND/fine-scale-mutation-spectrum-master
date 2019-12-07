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



def region_samplev2(L, chrom_sizes, N, fasta_file= ''):
    '''
    prepare sequence dictionary: {chrom:list(start pos)}
    provide to fasta_Rextract.
    
    ##
    '''
    
    chroms= list(chrom_sizes.keys())
    sizes= [int(chrom_sizes[x]) for x in chroms]
    sizes= np.array(sizes) / sum(sizes)
    
    choices= np.random.choice(chroms,size= N,p= sizes,replace= True)
    
    seqs_store= {
        z: len([x for x in range(N) if choices[x]==z]) for z in list(set(choices))
    }
    
    print(seqs_store)
    seqs= fasta_RextractUnif(fasta_file, seqs_store, L= L)
    
    return seqs



def fasta_RextractUnif(fasta,seq_store,L= 10000):
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
        
        chrom_seqs= return_seqs(result,size= seq_store[chrom],L= L)
        
        refseqs[chrom]= chrom_seqs
    
    return refseqs



def return_seqs(seq,size= 10,L= 1000,keep= ['A','T','G','C']):
    '''returns N=size segments of length L, unique.()=keep'''
    d= 0
    
    seqL= len(seq)
    seq_dict= {}
    
    while d < size:
        pos= np.random.randint(low= 0, high= seqL - L,size= 1)[0]
        given= seq[pos:(pos + L)]
        
        scrag= [x for x in given if x not in keep]
        
        if not scrag:
            seq_dict[pos]= given
            
            d += 1
    
    return seq_dict



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
        if isinstance(g,str):
            newline= defConst.format('"{}"'.format(v),'"{}"'.format(g))
        else:
            newline= defConst.format('"{}"'.format(v),str(g))
        
        lines.insert(init_idx+1,newline)
    
    with open(new_recipe,'w') as f:
        f.write(''.join(lines))
    
    return new_recipe




def SLiM_dispenserv1(sim_store, sim_recipe, cookID= 'ID', slim_dir= './', batch_name= '',
                    ID= 'v1',L= 10000, logSims= 'sims.log', mutlog= 'toMut.log'):
    ''' execute SLiM program
    - simulation specific recipe:
    - recipe template is re-written to direct to new fasta.
    '''
    nf= len(sim_store)
    for SIMname in sim_store.keys():
        
        command_line_constants= sim_store[SIMname]

        ### generate modified slim recipe
        new_recipe= process_recipe(sim_recipe,command_line_constants, SIMname)

        seed= np.random.randint(0,high= nf,size= 1)[0]
        ### Launch SLiM through shell.
        slim_soft= slim_dir + 'slim*' 

        command_units= [slim_soft, '-m', '-s', str(seed), new_recipe]

        os.system(' '.join(command_units))

        os.system('gzip {}'.format(sim_store[SIMname]["vcf_file"]))
        os.system('gzip {}'.format(sim_store[SIMname]["fasta_file"]))
        os.system('rm {}'.format(new_recipe))

        now = datetime.now()
        tnow= now.strftime("%d/%m/%Y %H:%M:%S")
        
        constant_str= ';'.join(['='.join([v,str(g)]) for v,g in command_line_constants.items()])

        with open(logSims,'a') as fp:
            INFO= [SIMname,tnow,sim_recipe,cookID,'L='+str(L),constant_str]
            fp.write('\t'.join(INFO) + '\n')
        
        with open(mutlog,'a') as fp:
            fp.write(SIMname + '\n')



def mutation_counter_launch(logfile,count_dir= './count/', 
                dir_launch= '..',main_dir= './', outlog= 'muted.log'):
    '''
    launch mutation counter.
    - read mut.log to know which have not been yet processed.
    - launch process_chromosomes.py using simulation name. 
    '''
    with open(logfile,'r') as fp:
        lines= fp.readlines()
    
    
    sims= [x.strip() for x in lines]
    chroms= [x.split('.')[0].split('C')[-1].strip('chr') for x in sims]
    
    job= 'python process_chromosomes.py -c {} -r {} -s {} -v {}_ -q {} -d {}'
    
    sims= [job.format(chroms[x],*[sims[x]]*4,dir_launch) for x in range(len(sims))]
    
    os.chdir(count_dir)
    for sim in sims:
        
        os.system(sim)
    
    os.chdir(main_dir)

    with open(outlog,'a') as fp:
        fp.write('\n' + ''.join(lines))

    open(logfile,'w').close()



def write_popIDs(sampleSizes,popSuff= 'pop',indSuff= 'i',
                file_dir= '/mnt/d/GitHub/fine-scale-mutation-spectrum-master/slim_pipe/mutation_counter/data/'):
    '''create identifiers file ind/pops'''
    
    popNames= [popSuff + str(x) for x in range(len(sampleSizes))]
    indIDs= [indSuff + str(x) for x in range(sum(sampleSizes))]
    
    popIDs= np.repeat(np.array(popNames),sampleSizes)
    data= np.array([
        indIDs,
        popIDs
    ]).T
    
    filename= file_dir + '/ind_assignments.txt'
    
    with open(filename,'w') as f:
        vector= ['\t'.join(x) for x in data]
        vector= '\n'.join(vector)
        f.write(vector)



####################################################################
####################################################################

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



############
############# deprecated

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


def return_seqs_deprecated(seq,size= 10,L= 1000,keep= ['A','T','G','C']):
    
    d= 0
    
    seqL= len(seq)
    seq_dict= {}
    
    while d <= size:
        pos= np.random.randint(low= 0, high= seqL - L,size= 1)[0]
        given= seq[pos:pos + L]
        scrag= [x for x in given if x not in keep]
        if not scrag:
            seq_dict[pos]= given
            d += 1



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


#############################################
#############################################
