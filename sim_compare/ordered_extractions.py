import itertools as it
import numpy as np
import tempfile
import os
import gzip



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




def region_samplev3(L, chrom_sizes, N, fasta_file= '', 
        fasta_utils= {}):
    '''
    prepare sequence dictionary: {chrom:list(start pos)}
    provide to fasta_Rextract.
        - v3: is provided functions that determine fasta retrieval behaviour.    
    ##
    '''
    
    chroms= list(chrom_sizes.keys())
    sizes= [int(chrom_sizes[x]) for x in chroms]
    
    sizes= np.array(sizes) / float(sum(sizes))
    
    choices= np.random.choice(chroms,size= N,p= sizes,replace= True)
    
    seqs_store= {
        z: len([x for x in range(N) if choices[x]==z]) for z in list(set(choices))
    }
    
    func= fasta_utils['p1']
    seqs= func(fasta_file, seqs_store, fasta_utils['p2'], L= L)
    
    return seqs



def fasta_RextractUnif_v1(fasta,seq_store,seqSeek,L= 10000):
    ''' 
    Extract regions from fasta file.
    Takes dictionary {chrom: list(start points)};
    extracts sequences after reading each chromosome.
    - v1: takes the seqSeek function, which can now be modified to return other output.
        Here we will want it to return start and end position alongside the sequence.
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
        
        chrom_seqs= seqSeek(result,size= seq_store[chrom],L= L)
        
        refseqs[chrom]= chrom_seqs
    
    return refseqs




def return_seqs_v1(seq,size= 10,L= 1000,keep= ['A','T','G','C']):
    '''
    returns N=size segments of length L, unique.()=keep
    - v1: returns start and end position as well as seq.
    '''
    d= 0
    
    seqL= len(seq)
    seq_dict= {}
    
    while d < size:
        pos= np.random.randint(low= 0, high= seqL - L,size= 1)[0]
        given= seq[pos:(pos + L)]
        
        scrag= [x for x in given if x not in keep]
        
        if len(scrag) == 0:
            seq_dict[str(pos)]= {
                'fasta':given,
                'start': str(pos),
                'end': str(pos + L)
                }
            
            d += 1
    
    return seq_dict



def line_subtract(line, start=0):
    '''
    Takes line of vcf file, subtracts start from 
    '''
    s= line.split(b'\t')
    pos= s[1].decode()
    pos= int(pos) - start
    pos= str(pos).encode()
    s[1]= pos
    s= b'\t'.join(s)

    return s


def write_fastaEx(fasta,chrom= '1',start= 0, ID= 'SIM', fasta_dir= ''):
    ''' write extracted sequence as fasta file'''
    filename= ''.join([fasta_dir,'/chr',chrom,'_',ID,'.fa'])
    
    header= ' '.join(['>'+ID,chrom,str(start),str(start + len(fasta))]) + '\n'
    
    with open(filename,'w') as fp:
        fp.write(header)
        fp.write(fasta)
    
    return filename




def vcf_subtract(vcf_file,start= 0, fasta= ''):
    '''
    take vcf, subtract start (int) from every position. 
    '''
    dir_vcf= vcf_file.split('/')
    dir_vcf= '/'.join(dir_vcf[:-1]) + '/'

    tempfile= dir_vcf + 'temp.vcf.gz'

    with gzip.open(tempfile,'w') as f:
        file_path = (vcf_file)
        infile = gzip.open(vcf_file,'r')

        line = infile.readline()
        f.write(line)
        while not line.startswith(b'#CHROM'):
            line = infile.readline()
            f.write(line)

        for line_counter, line in enumerate(infile):
            line= line_subtract(line, start)
            if fasta:
                s= line.split(b'\t')
                pos= int(s[1].decode())
                ref= s[3].decode()

                print('cc. {}'.format(pos))
                print('cc.{},{}'.format(fasta[pos-2:pos+1],ref))
            f.write(line)

        infile.close()

    os.system('rm {}'.format(vcf_file))
    os.system('mv {} {}'.format(tempfile,vcf_file))




def region_extract_launch(seq_dict,vcf_file,extract_bash= "region_extract.sh",popIDs= 'ind_assignments.txt',out_dir= "", 
                            batchName= '', logfile= 'regions.log'):
    '''
    takes dictionary : {CHROM: {POS: {PosI,PosF,Fasta}}},
    for each leaf of the dictionary, launches bash script that itself runs 
    vcf tools to extract SNPs in the region between PosI and PosF on CHROM.
    bash sctipt arguments:
    - vcf_dir=$1
    - vcf_file=$2
    - chrom=$3
    - start=$4
    - end=$5
    - out_dir=$6

    VCF and Fasta files are written to specific directory in *out_dir*, name dermined as:
    - _something_;
	
    VCF is reprocessed so every position = pos - PosI.

    '''
    dir_vcf= vcf_file.split('/')
    vcf_file= dir_vcf[-1]
    if len(dir_vcf) > 1:
    	dir_vcf= '/'.join(dir_vcf[:-1]) + '/'
    else: dir_vcf= './'

    for chrom in seq_dict.keys():
        for pos in seq_dict[chrom]:
            end= seq_dict[chrom][pos]['end']

            regionID= batchName + 'C{}.{}'.format(chrom,pos)
            region_dir= out_dir + regionID + '/'

            os.makedirs(region_dir, exist_ok=True)
            new_vcf= region_dir + regionID + '_chr{}.vcf.gz'.format(chrom)
            #command_lines= ['sbatch',extract_bash,dir_vcf,vcf_file,'chr'+chrom,pos,end,new_vcf]
            command_lines= ['vcftools', '--gzvcf', dir_vcf + vcf_file,'--chr', 'chr' + chrom,
            '--from-bp', pos,'--to-bp', end,
            '--recode', '--stdout', '|', 'bgzip', '>', new_vcf]

            os.system(' '.join(command_lines))
            
            fasta= seq_dict[chrom][pos]['fasta']
            print(fasta[:10])
            vcf_subtract(new_vcf,start= int(pos), fasta= fasta)

            L= str(len(fasta))
            fasta_file= write_fastaEx(fasta,chrom= chrom,start= int(pos),
                    ID= regionID, fasta_dir= region_dir)
            os.system('gzip {}'.format(fasta_file))

            ind_command= ['cp',popIDs,region_dir + 'ind_assignments']
            os.system(' '.join(ind_command))
            #write_popIDs(sample_sizes, file_dir= region_dir)

            INFO= [regionID,'L='+L]

            with open(logfile,'a') as f:
            	f.write('\n' + '\t'.join(INFO))



if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser()

    parser.add_argument('-v', '--vcf', type=str)

    parser.add_argument('-l', '--length', type=int, default= 1e5)

    parser.add_argument('-n', '--number', type=int, default= 10)

    parser.add_argument('-a', '--assembly', type=str, default= 'hg38')

    parser.add_argument('-b', '--batch', type=str,default= 'test')

    parser.add_argument('-o', '--out', type=str, default= '')

    args = parser.parse_args()


    ##
    ## Read chrom_sizes file to decide where to sample files from.
    size_dir= "/home/jgarc235/Rhesus/chrom_sizes/" 
    chrom_sizes= read_chrom_sizes(args.assembly, size_dir= size_dir)
    
    fastas_dir= "/home/jgarc235/Rhesus/Fastas/"
    fasta= fastas_dir + args.assembly + '.fa.gz'

    """
     region_samplev3(L, chrom_sizes, N, fasta_file= '', 
        fasta_utils= {}):
    """
    fasta_utils= {
        'p1': fasta_RextractUnif_v1,
        'p2': return_seqs_v1
    }
    rseqs= region_samplev3(args.length, chrom_sizes, args.number, fasta_file= fasta,
        fasta_utils= fasta_utils)  

    ##
    extractSH_dir= "/home/jgarc235/Rhesus/bash_commands/Extract/"
    bashExtract= "region_extract_mt.sh"

    assignments= '/home/jgarc235/Rhesus/Mutation_study/ind_assignments.txt'
    outdir= '/home/jgarc235/Sim_compare/data/'

    region_extract_launch(rseqs,args.vcf,extract_bash= bashExtract,popIDs= assignments,out_dir= args.out, 
                            batchName= args.batch, logfile= 'regions.log')








