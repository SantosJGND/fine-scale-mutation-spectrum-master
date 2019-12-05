import os
import errno
import argparse
import sys

from labels import populations
from common import open_infile, get_conserved, get_args
from mutation_counter import MutationCounter, IncludedRegion


def read_bed(chromosomes,bed_file):

    chrom_dict= {chrom:[] for chrom in chromosomes}

    with open(bed_file,'r') as f:
        for line in f:
            line= line.split()
            if line[0].strip('chr') in chromosomes:
                chrom_dict[chrom].append(line[2])

    return chrom_dict


def read_exclude(dirf='../data/bed_files/',filename='files_Regexclude.txt'):
    '''
    read files of regions to exclude in mutation analysis. 
    list file of files in "../bed_files/" to read.
    '''
    filename= dirf + filename

    with open(filename) as f:
        files= f.readlines()

    files= [x.strip() for x in files]

    return files


def get_finescale(mutation_counter,vcf_file,dir_launch,vcf_dir= 'vcf_data',bed_filter=[]):
    infile, line = open_infile(mutation_counter.chrom,vcf_file,vcf_dir=vcf_dir,dir_launch=dir_launch)
    
    print('configuring..')
    mutation_counter.configure(line,dir_launch)
    
    print('processing lines.')
    for line_counter, line in enumerate(infile):
        mutation_counter.process_line(line)
    
    print('writing.')
    mutation_counter.write_output()



if __name__ == '__main__':
    chromosomes, reference, short, vcf_file, dir_launch, vcf_dir, bed_filter, exc = get_args()
    
    bed_tag=['','_lb_'][int(len(bed_filter) > 0)]

    outfile_dir= dir_launch+'/'+ reference + bed_tag + '_finescale_mut_spectra_vcf.' + vcf_dir
    os.makedirs(outfile_dir, exist_ok=True)
    outfile_dir= outfile_dir +  '/'
    
    try:
        os.mkdir(outfile_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    if len(bed_filter):
        chrom_dict= read_bed(chromosomes,bed_filter)
    else:
        chrom_dict= {x:[] for x in chromosomes}

    if exc:
        exclude= read_exclude()
    else:
        exclude= []

    for chrom in chromosomes:
        chrom = str(chrom)
        included_regions = []

        output={population: 'Mut\n' for population in populations}
        outfile_path = outfile_dir + 'mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
        conserved = [[0, 1e12]]
        included_regions.append(IncludedRegion(chrom, output, outfile_path, conserved))

        for file in exclude:
            file_name= file.split('.')[0]
            output = {population: 'Ref Alt \n' for population in populations}
            outfile_path = outfile_dir + file_name + '_mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
            conserved = get_conserved(dir_launch+'/data/bed_files/' + file, chrom)
            included_regions.append(IncludedRegion(chrom, output, outfile_path, conserved))

            #output = {population: 'Ref Alt \n' for population in populations}
            #outfile_path = '../'+short+'_finescale_mut_spectra/phyloP_conserved_mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
            #conserved = get_conserved(dir_launch+'/data/bed_files/phastConsElements100way.txt.gz', chrom)
            #included_regions.append(IncludedRegion(chrom, output, outfile_path, conserved))

        mutation_counter = MutationCounter(chrom, included_regions,reference,short)
        get_finescale(mutation_counter,vcf_file,dir_launch,vcf_dir= vcf_dir,bed_filter=chrom_dict[chrom])

        print('finished chrom'), chrom


