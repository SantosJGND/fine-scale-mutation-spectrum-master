import sys
import argparse
import gzip

from mutations import mutations

import tempfile
import os



def get_args():
    valid_chromosomes = ['X'] + [str(i) for i in range(1, 21)]

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--id', type=str,
                        default='')

    parser.add_argument('-c', '--chromosomes', type=str, nargs='+',
                        default=valid_chromosomes)

    parser.add_argument('-r', '--ref', type=str,
                        default='rheMac8')

    parser.add_argument('-s', '--short', type=str,
                        default='panTrop5')

    parser.add_argument('-v', '--vcf', type=str,
                        default='')

    parser.add_argument('-p', '--pops', type=str,
                        default='ind_assignments.txt')

    parser.add_argument('-d', '--dir', type=str,
                        default='..')

    parser.add_argument('-q', '--vcfdir', type=str,
                        default='..')

    parser.add_argument('-b', '--bedfilter', type=str,
                        default='')

    parser.add_argument('-e', '--exclude',
                        default=False, action='store_true')


    #parser.add_argument('-a', '--annotations', type=str, nargs='+', default= [])

    args = parser.parse_args()

    chromosomes = parser.parse_args(sys.argv[1:]).chromosomes

    return chromosomes, args.ref, args.short, args.vcf, args.dir, args.vcfdir, args.bedfilter, args.exclude, args.pops, args.id


def get_chromosomes_from_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chromosomes', type=str, nargs='+',
                        default=valid_chromosomes)

    chromosomes = parser.parse_args(sys.argv[1:]).chromosomes

    return chromosomes



def reference_sequence(chromosome_number,reference,dir_launch='..'):

    file_path_template = dir_launch+'/data/sims/{}/chr{}_{}.fa.gz'.format(reference,chromosome_number,reference)
    
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



def get_human_chimp_differences(chromosome_number,reference,short,dir_launch='..'):
    human_chimp_differences = {}

    file_path= ''.join([dir_launch,'/data/diffs/',reference,'_',short,'_diffs/',
        reference,'_',short,'_diffs_chr', chromosome_number,'.txt.gz'])
	
    file_path= (file_path)

    infile= gzip.open(file_path,'r')
    d= 0
    
    for line in infile:
        line=line.decode()
        if d > 0 and 'SNP' in line:
            position, _, _, chimp_allele = line.split()
            human_chimp_differences[int(position)] = chimp_allele
        
        d += 1

    return human_chimp_differences



def get_column_indices(column_labels, populations, sample_id_to_population):
    #column_labels= [x.strip(b'0_') for x in column_labels]

    population_to_column_indices = {
        population: [] for population in populations
    }

    sample_ids = column_labels[9:]
    
    for i, sample_id in enumerate(sample_ids):
        
        population_to_column_indices[
            sample_id_to_population[sample_id]
        ].append(i + 9)

    return population_to_column_indices


def get_column_index_to_population(column_labels, sample_id_to_population):
    #column_labels= [x.strip(b'0_') for x in column_labels]
    sample_ids = column_labels[9:]
    
    return {
        i + 9: sample_id_to_population[sample_id]
        for i, sample_id in enumerate(sample_ids)
    }



def initialize_mut_count(population_to_column_indices,populations):
    mut_count = {}
    for pop in populations:
        for mut in mutations:
            for i in range(1, 2*len(population_to_column_indices[pop])+1):
                mut_count[(mut, pop, i)] = 0

    return mut_count


def write_output(output, outfile_path, indices, mut_count, populations):
    for pop in populations:
        for mut in mutations:
            output[pop] += mut[0] + '_' + mut[1]
            for i in range(1, 2 * len(indices[pop]) + 1):
                output[pop] += ' ' + str(mut_count[(mut, pop, i)])
            output[pop] += '\n'
        filename= outfile_path % pop.decode()
        
        outfile = open(filename, 'w')
        outfile.write(output[pop])
        outfile.close()


def open_infile(chrom,vcf_file,suff='chr',vcf_dir= 'vcf_data',dir_launch='..'):

    filename= ''.join([dir_launch,'/data/sims/',vcf_dir,'/',vcf_file,suff,chrom,'.vcf.gz'])
    print('read from vcf: {}'.format(filename))

    file_path = (filename)
    infile = gzip.open(file_path)

    line = infile.readline()
    while not line.startswith(b'#CHROM'):
        line = infile.readline()

    return infile, line


def get_conserved(infile_path, chrom):
    print(infile_path)
    infile = gzip.open(infile_path)
    lines = infile.readlines()
    infile.close()

    chrom= str.encode(chrom)

    ind = 0
    print(lines[ind])
    s = lines[ind].split(b'\t')

    while not s[1] == b'chr' + chrom:
        ind += 1

        s = lines[ind].split(b'\t')

    conserved = [(int(s[2]), int(s[3]))]
    ind += 1
    s = lines[ind].split(b'\t')

    while ind < len(lines) - 1 and s[1] == b'chr' + chrom:
        if int(s[2]) == conserved[-1][-1] + 1:
            new_tup = (conserved[-1][0], int(s[3]))
            conserved.pop()
            conserved.append(new_tup)
        else:
            conserved.append((int(s[2]), int(s[3])))
        ind += 1
        line = lines[ind]
        s = line.strip(b'\n').split(b'\t')

    return conserved
