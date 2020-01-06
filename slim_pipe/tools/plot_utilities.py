
import sys
import numpy as np
from scipy.stats import chi2_contingency
from itertools import product



comp = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}
ypos, ylabel = [], []

mut_index = {}
row, col = 0, 0

for b2, d in [('A', 'T'), ('A', 'C'), ('A', 'G'),
              ('C', 'T'), ('C', 'G'), ('C', 'A')]:
    for b1 in 'ACGT':
        col = 0
        ypos.append(row+0.5)
        if b1 == 'T' and b2 == 'C' and d == 'A':
            ylabel.append('5\'-'+b1)
        elif b1 == 'C':
            ylabel.append(b2+r'$\to$'+d+r'  '+b1)
        else:
            ylabel.append(b1)
        for b3 in 'ACGT':
            mut_index[(b1+b2+b3, d)] = (row, col)
            mut_index[(comp[b3]+comp[b2]+comp[b1], comp[d])] = (row, col)
            col += 1
        row += 1



def Population(population):
    if population not in valid_populations:
        raise ValueError
    return population


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


def frequency_breakdown(path, chromosomes, frequency_range):
    
    count_array = np.zeros((row, col))
    for chrom in chromosomes:
        infile = open(path % str(chrom))
        lines = infile.readlines()
        infile.close()

        s = lines[1].strip('\n').split(' ')
        # TODO(mason) do these better.
        start_index = 2
        while float(start_index - 2) / (len(s) - 2) < frequency_range[0]:
            start_index += 1

        end_index = len(s) - 1
        while float(end_index) / (len(s)-2) > frequency_range[1]:
            end_index -= 1

        for line in lines[1:]:
            s = line.strip('\n').split(' ')
            for i in range(start_index, end_index):
                count_array[mut_index[(s[0][:3], s[0][4])]] += int(s[i])
    return count_array


def heatmap(chromosomes, population_pair, frequency_range, exclude, 
                p_value, short,muted_dir,tag= '',output= 'pval'):

    outdir= muted_dir + '{}{}_finescale_mut_spectra_vcf.{}/'.format(short,tag,short)

    if exclude:
        files= read_exclude()
    else:
        files= {}

    pop_counts = {}
    num_variants = {}

    for pop in population_pair:
        path = (outdir + 'mut_type_v_allele_freq_' + pop + '_chr%s_nosingle.txt')
        pop_counts[pop] = frequency_breakdown(path, chromosomes,
                                              frequency_range)
        if exclude:

            for file in files:
                file_name= file.split('.')[0]
                repeats_path = (outdir + file_name + '_mut_type_v_allele_freq_' +
                                pop + '_chr%s_nosingle.txt')
                pop_counts[pop] -= frequency_breakdown(repeats_path, chromosomes,
                                                       frequency_range)

                #conserved_path = ('../{}_finescale_mut_spectra/'.format(short) + 
                #                  'phyloP_conserved_mut_type_v_allele_freq_' +
                #                  pop + '_chr%s_nosingle.txt')
                #pop_counts[pop] -= frequency_breakdown(conserved_path, chromosomes,
                #                                       frequency_range)

        num_variants[pop] = pop_counts[pop].sum()


    refpop, pop = population_pair

    ratio_grid = np.zeros((row, col))
    sig_x, sig_y = [], []
    for i in range(row):
        for j in range(col):
            chi_array= np.array([
                    [pop_counts[pop][i][j], num_variants[pop]],
                    [pop_counts[refpop][i][j], num_variants[refpop]]
                ])

            chi_0= np.sum(chi_array,axis= 1)
            chi_1= np.sum(chi_array,axis= 0)

            if chi_0[0] == 0 or chi_0[1] == 0:
                ratio_grid[i][j] = np.nan
                sig_x.append(j+0.5)
                sig_y.append(i+0.5)
            
            elif chi_1[0] == 0 or chi_1[1] == 0:
                ratio_grid[i][j] = 1

            else:
                #print(chi_array)
                ##
                _, this_pval, _, _ = chi2_contingency(
                    chi_array
                )
                if output == 'pval':
                	ratio_grid[i][j] = this_pval
                else:
	                ratio_grid[i][j] = (pop_counts[pop][i][j] * num_variants[refpop] /
	                                    (num_variants[pop] * pop_counts[refpop][i][j]))
                if this_pval < p_value:
                    sig_x.append(j+0.5)
                    sig_y.append(i+0.5)

    return ratio_grid, (sig_x, sig_y)


def make_titles(chromosome_groups, population_pairs, frequency_range, exclude,
                p_value):
    title = ''

    distinct_population_pairs = len(set(population_pairs)) != 1
    if not distinct_population_pairs:
        title += '/'.join(population_pairs[0])

    distinct_chromosome_groups = (
        len(set(tuple(cg) for cg in chromosome_groups)) != 1
    )
    if not distinct_chromosome_groups:
        if set(chromosome_groups[0]) == set(str(c) for c in range(1, 23)):
            title += ' Autosomes'
        elif len(chromosome_groups[0]) == 1:
            title += ' Chromosome %s' % str(chromosome_groups[0][0])
        else:
            title += ' Chromosomes %s' % ' '.join(
                str(c) for c in chromosome_groups[0]
            )
    # if exclude:
    #     title += '\nExcluding Repeats and Conserved Regions'
    if frequency_range != [0, 1]:
        title += ' %.2f <= f <= %.2f' % tuple(frequency_range)
    title += ' p < %1.1e' % p_value

    column_titles = []
    for chromosomes, population_pair in product(chromosome_groups,
                                                population_pairs):
        column_title = ''
        if distinct_population_pairs:
            column_title += '/'.join(population_pair)
        if distinct_chromosome_groups:
            column_title += ' '
            column_title += ' '.join(str(c) for c in chromosomes)
        column_titles.append(column_title)

    return title, column_titles


