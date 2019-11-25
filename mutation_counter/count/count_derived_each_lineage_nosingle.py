import gzip

from itertools import product
import numpy as np

from labels import (
    groups,
    population_to_group,
    sample_id_to_population,
)
from mutations import mutations, bases
from common import get_chromosomes_from_args, reference_sequence, get_human_chimp_differences


class BadDataQualityError(Exception):
    pass


def parse_line(line):
    ''' parse vcf line. alleles are returned as flattened vector.'''

    (
        chromosome_number,
        position,
        _,  # SNP id
        reference_allele,
        alternate_allele,
        _,  # quality score,
        filter_result,
        info,
        _,  # format (what ever this is???) ## good point
        haplotypes
    ) = line.split(None, 9)

    if (
        reference_allele not in list(bases) or
        alternate_allele not in list(bases) or
        filter_result != 'PASS'
    ):
        raise BadDataQualityError

    position = int(position)

    ## this doesn't do it.
    #alleles = haplotypes[::2]  # remove '\t' and '|' separators 

    ## this. simplifies allele addition later too.
    alleles= ''.join([x.split(':')[0] for x in haplotypes.split('\t')])
    alleles= [x for x in alleles if x.isdigit()]
    alleles= np.array(alleles,dtype= int)
    derived_count= sum(alleles)

    return reference_allele, alternate_allele, position, derived_count, alleles


def get_mutation(position, reference_allele, alternate_allele, refseq, chimp_alleles):
    '''compare vcf reference allele to reference in fasta. return binary conclusion and fasta kmer context'''
    context = refseq[position - 2 : position + 1]
    if 'N' in context:
        raise BadDataQualityError

    if chimp_alleles.get(position) == alternate_allele:
        derived_allele='0'
        this_mut=(context[0] + alternate_allele + context[2], reference_allele)
    else:
        derived_allele='1'
        this_mut=(context, alternate_allele)

    return derived_allele, this_mut


def process_line(line, refseq, chimp_alleles):
    ''' transit line between parse and get mutation functions '''
    (
        reference_allele,
        alternate_allele,
        position,
        derived_count,
        alleles
    ) = parse_line(line)
    derived_allele, this_mut = get_mutation(position, reference_allele,
        alternate_allele, refseq, chimp_alleles)
    return derived_allele, derived_count, this_mut, alleles



def update_counts(alleles, mutation_counts, derived_count, derived_allele, n_lineages, this_mut):
    ''' count mutations by individual. store in mutation_count dictionary by mutation type.'''
    if derived_count>1 and n_lineages-derived_count>1:
        if derived_allele == '0':
            derived_count=n_lineages-derived_count
            alleles= 1 - alleles

        der_observed=0

        for i, allele in enumerate(alleles):
            #if allele == derived_allele:
            mutation_counts[(this_mut, i)] += allele
            der_observed += allele


        if der_observed != derived_count:
            print('{}-{}'.format(derived_count,der_observed))
            print(alleles)

        assert der_observed == derived_count



def write_output(mutation_counts, sample_ids, mutations, n_lineages, chrom,outfile_dir= '..'):
    output='Mut_type '
    output += ' '.join(
        [sample_id_to_population[sample_id.encode()].decode() for sample_id in sample_ids]
    )
    output+='\n'

    for mut in mutations:
        output+=mut[0]+'_'+mut[1]
        for i in range(n_lineages):
            output+=' '+str(mutation_counts[(mut,i)])
        output+='\n'

    with open(outfile_dir + 'derived_each_lineage_chr'+chrom+'_nosingle.txt','w') as outfile:
        outfile.write(output)



def process_chromosome(chrom,ref,short,vcf_file,suff='chr',vcf_dir= 'vcfs',dir_launch='..',outfile_dir= '..'):
    '''
    Platform function
        read vcf file;
        read fasta;
        read reference to vcf differences; 

        dispatch each genotype line to process_line(); 
        update mutation / individual library counts;
        write
    '''
    # check args
    refseq=reference_sequence(chrom,ref,dir_launch=dir_launch)
    refseq= refseq.decode()

    # check args
    chimp_alleles=  get_human_chimp_differences(chrom, ref, short,dir_launch=dir_launch)

    ## mutations imported

    # check args: chrom,vcf_file,suff='chr',vcf_dir= 'vcfs',dir_launch='..'
    gzip_path= ''.join([dir_launch,'/data/',vcf_dir,'/',vcf_file,suff,chrom,'.vcf.gz'])
    # gzip_path = '../data/vcfs/ALL.chr'+chrom+'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
    with gzip.open(gzip_path) as infile:
        c= 0
        d= 0
        for line in infile:
            line= line.decode()
            c += 1

            if line.split()[0] == '#CHROM':
                
                sample_ids = line.split()[9:]
                
                n_lineages = 2 * len(sample_ids)

                mutation_counts = {
                    (mutation, haplotype_index): 0
                    for haplotype_index in range(n_lineages)
                    for mutation in mutations
                }
                d= 1

            if d == 1:

                try:
                    derived_allele, derived_count, this_mut, alleles = process_line(line, refseq, chimp_alleles)
                except BadDataQualityError:
                    continue

                update_counts(alleles, mutation_counts, derived_count, derived_allele, n_lineages, this_mut)
    
    write_output(mutation_counts, sample_ids, mutations, n_lineages, chrom, outfile_dir)




def get_args():
    valid_chromosomes = ['X'] + [str(i) for i in range(1, 21)]

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chromosomes', type=str, nargs='+',
                        default=valid_chromosomes)

    parser.add_argument('-r', '--ref', type=str,
                        default='rheMac8')

    parser.add_argument('-s', '--short', type=str,
                        default='panTrop5')

    parser.add_argument('-v', '--vcf', type=str,
                        default='')

    parser.add_argument('-d', '--dir', type=str,
                        default='..')

    parser.add_argument('-q', '--vcfdir', type=str,
                        default='..')

    parser.add_argument('-b', '--bedfilter', type=str,
                        default='')


    #parser.add_argument('-a', '--annotations', type=str, nargs='+', default= [])

    args = parser.parse_args()

    chromosomes = parser.parse_args(sys.argv[1:]).chromosomes
    for chrom in chromosomes:
        assert chrom in valid_chromosomes

    return chromosomes, args.ref, args.short, args.vcf, args.dir, args.vcfdir, args.bedfilter


    return chromosomes, args.ref, args.short, args.vcf, args.dir, args.vcfdir


if __name__ == '__main__':
    import argparse
    import sys

    chromosomes, reference, short, vcf_file, dir_launch, vcf_dir, bed_filter = get_args()

    bed_tag=['','_lb_'][int(len(bed_filter) > 0)]

    outfile_dir= dir_launch+'/'+ reference + bed_tag + '_finescale_mut_spectra_vcf.' + vcf_dir + '/'

    for chrom in chromosomes:
        chromosome = str(chrom)

        # def process_chromosome(chrom,ref,short,dir_launch,suff='chr',vcf_dir= 'vcfs',dir_launch='..'):
        process_chromosome(chromosome,reference,short,vcf_file,
            dir_launch= dir_launch,vcf_dir= vcf_dir, outfile_dir= outfile_dir)

        print('finished chrom {}'.format(chromosome))
