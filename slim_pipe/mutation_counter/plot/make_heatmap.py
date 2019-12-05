import sys
import argparse
import numpy as np
from scipy.stats import chi2_contingency
from itertools import product

import matplotlib
matplotlib.use('Agg')  # This prevents the plotting engine from starting up.
import matplotlib.pyplot as plt  # noqa E402

K= 5

valid_populations = ['littoralis', 'brevicaudus', 'tcheliensis', 'lasiotis', 'mulatta', 'CH']



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
    #print(frequency_range)
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


def heatmap(chromosomes, population_pair, frequency_range, exclude, p_value, short,vcf_data):
    pop_counts = {}
    num_variants = {}
    outdir= '../{}_finescale_mut_spectra_vcf.{}/'.format(short,vcf_data)

    for pop in population_pair:
        path = (outdir + 'mut_type_v_allele_freq_' +
                pop + '_chr%s_nosingle.txt')
        pop_counts[pop] = frequency_breakdown(path, chromosomes,
                                              frequency_range)
        if exclude:
            files= read_exclude()
            print("excluding regions in:")
            print(files)
            
            for file in files:
                file_name= file.split('.')[0]
                repeats_path = (outdir + file_name + '_mut_type_v_allele_freq_' +
                                pop + '_chr%s_nosingle.txt')
                pop_counts[pop] -= frequency_breakdown(repeats_path, chromosomes,
                                                       frequency_range)
                

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

            
            _, this_pval, _, _ = chi2_contingency(
                chi_array
            )
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


def make_plot(ratio_grids, significant_indices, plot_title, column_titles, pop_pair_names):
    combined_grids = np.hstack(ratio_grids)
    r = max(1-np.min(combined_grids), np.max(combined_grids)-1)
    plt.pcolor(np.hstack(ratio_grids), vmin=1-r, vmax=1+r, cmap='seismic')
    plt.colorbar()

    xtick_values = len(ratio_grids) * list('ACGT')
    xtick_values[0] = "3'-A"
    plt.xticks(np.arange(.5, 4 * len(ratio_grids)), xtick_values, rotation=90)
    plt.yticks(ypos, ylabel)

    for k in range(len(ratio_grids)):
        plt.axvline(x=4 * k, color='black')

    for k in range(0, len(ratio_grids[0]), 4):
        plt.axhline(y=k, color='black', linestyle='dashed')

    top_x_axis = plt.twiny()  # Create a twin Axes sharing the yaxis, lol
    top_x_axis.set_xticks(np.arange(2, 4 * len(ratio_grids), 4))
    top_x_axis.set_xticklabels(column_titles)

    combined_x_positions = []
    combined_y_positions = []
    for i, (x_positions, y_positions) in enumerate(significant_indices):
        for x, y in zip(x_positions, y_positions):
            combined_x_positions.append(x + 4 * i)
            combined_y_positions.append(y)
    plt.scatter(combined_x_positions, combined_y_positions, marker='.',
                color='white')

    plt.xlim(0, len(combined_grids[0]))
    plt.ylim(0, len(combined_grids))

    plt.title(plot_title + '\n')

    plt.gcf()
    pop_pair_names= '_'.join(pop_pair_names)
    plt.savefig('heatmap_{}.pdf'.format(pop_pair_names), format='pdf')
    plt.clf()



def Population(population):
    if population not in valid_populations:
        raise ValueError
    return population


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chromosomes', type=str, nargs='+',
                        default=[str(c) for c in range(1, 21)])
    parser.add_argument('-r', '--ref', type=str,
                        default='rheMac8')
    parser.add_argument('-v', '--vcf', type=str,
                        default='vcf_data')
    parser.add_argument('-i', '--individually', action='store_true')
    parser.add_argument('-p', '--pops', type=Population,
                        nargs='+')
    parser.add_argument('-f', '--frequency-range', type=float, nargs=2,
                        default=[0, 1])
    parser.add_argument('-e', '--exclude', action='store_true')
    parser.add_argument('--p-value', type=float, default=1e-5)

    args = parser.parse_args(sys.argv[1:])

    chromosomes = args.chromosomes

    for chromosome in chromosomes:
        assert chromosome == 'X' or int(chromosome) in range(1, 23)
    if args.individually:
        chromosome_groups = [[chromosome] for chromosome in chromosomes]
    else:
        chromosome_groups = [chromosomes]


    return (chromosome_groups, args.pops, args.frequency_range,
            args.exclude, args.p_value,args.ref, args.vcf)


if __name__ == '__main__':
    (chromosome_groups, population_pairs, frequency_range, exclude,
     p_value,reference,vcf_data) = parse_args()

    pop_pair_names= zip(population_pairs[::2],
                               population_pairs[1::2])

    pop_pair_names= ['-'.join(list(x)) for x in pop_pair_names]


    population_pairs = zip(population_pairs[::2],
                           population_pairs[1::2])

    heatmaps = [
        heatmap(
            chromosomes, population_pair, frequency_range, exclude, p_value,reference,vcf_data
        ) for chromosomes, population_pair in product(chromosome_groups,
                                                      population_pairs)
    ]

    ratio_grids, significant_indices = zip(*heatmaps)

    plot_title, column_titles = make_titles(
        chromosome_groups, population_pairs, frequency_range, exclude, p_value
    )

    make_plot(ratio_grids, significant_indices, plot_title, column_titles,pop_pair_names)
